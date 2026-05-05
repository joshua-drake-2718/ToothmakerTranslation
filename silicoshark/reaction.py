"""Reaction-diffusion step (eqs. 14-18 of the 2010 paper).

One forward-Euler update of Act/Inh/Sec. Pure Jacobi: all rates are
computed against the start-of-step concentrations, then applied at once.

Equations implemented (citations to Salazar-Ciudad and Jernvall 2010):

  d[Act]/dt = k_act * [Act] / (1 + k_inh * [Inh])
              - k_deg * [Act]
              + k_da  * lap[Act]                                 (eq. 14)

      ^ Note: paper prints `1 + k_inh * [Act]` in the denominator, but
      this is a typo — the FORTRAN computes `1 + k_inh * [Inh]`,
      which is the biologically meaningful self-saturating inhibition
      and the canonical form. PATH_B_DEFAULT follows the FORTRAN per
      the charter; `PAPER_LITERAL_2010` exposes the paper-as-written
      form via `Discretisation.eq14_denominator='act_typo'`.
      See docs/research/paper-review-2010-salazar-ciudad-jernvall.md
      §'Paper-vs-code findings'.

  d[Inh]/dt = -k_deg * [Inh] + k_di * lap[Inh]                   (eq. 15)
       (default for non-secreting cells)

  d[Inh]/dt =  [Act] - k_deg * [Inh] + k_di * lap[Inh]           (eq. 17)
       (cells where d_i >= k_int OR cell is a knot;
        FORTRAN variant uses (rate of [Act]) * d_i in place of [Act])

  d[Sec]/dt = -k_deg * [Sec] + k_ds * lap[Sec]                   (eq. 16)
       (default for non-secreting cells)

  d[Sec]/dt =  k_sec - k_deg * [Sec] + k_ds * lap[Sec]           (eq. 18)
       (cells where d_i >= k_set OR cell is a knot;
        FORTRAN variant uses k_sec * d_i in place of k_sec)

Knot detection rule (paper p. 586): "Differentiation is irreversible.
A cell becomes an enamel knot when [Act] >= 1." We apply this after the
Euler step against the new [Act] and store the boolean flag persistently
(once a knot, always a knot). Under the FORTRAN's
`knot_threshold_gate='first_border_cell'` rule, only cells whose index
is `>= state.first_border_cell` are eligible.

Path B v2 wires every paper-vs-FORTRAN choice through the
`Discretisation` dataclass (charter
`docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md`).
The default `PATH_B_DEFAULT` reproduces the v1 behaviour byte-for-byte
on the seal example.
"""
from __future__ import annotations

import numpy as np

from .params import Params
from .state import State
from .mesh import Mesh
from .discretisation import Discretisation, PATH_B_DEFAULT


def step_reaction_diffusion(
    state: State,
    params: Params,
    mesh: Mesh,
    dt: float,
    disc: Discretisation = PATH_B_DEFAULT,
) -> None:
    """Apply one forward-Euler step to act/inh/sec in-place.

    Pure Jacobi: rates computed against current state, applied once.

    Diffusion uses the length-weighted graph Laplacian on the triangular
    mesh (see `mesh.Mesh`). Reaction terms branch on
    `disc.eq14_denominator`, `disc.eq17_inh_source`,
    `disc.eq18_sec_source`, and `disc.knot_threshold_gate`. Each branch
    is a single short conditional with a one-line comment naming the
    field, so a reader scanning the function sees every option
    exercised.
    """
    act = state.act
    inh = state.inh
    sec = state.sec
    diff = state.diff
    knot = state.knot

    # --- Eq. 14 activator rate ---------------------------------------
    # disc.eq14_denominator: paper-as-written (`act_typo`) vs FORTRAN-
    # corrected (`inh_corrected`). The corrected form is the
    # biologically meaningful self-saturating inhibition by Inh.
    if disc.eq14_denominator == 'inh_corrected':
        rd_act = params.k_act * act / (1.0 + params.k_inh * inh) - params.k_deg * act
    elif disc.eq14_denominator == 'act_typo':
        rd_act = params.k_act * act / (1.0 + params.k_inh * act) - params.k_deg * act
    else:
        raise ValueError(f'unknown eq14_denominator: {disc.eq14_denominator!r}')

    # --- Eq. 17 inhibitor source -------------------------------------
    # disc.eq17_inh_source: paper Eq. 17 ([Act] concentration) vs
    # FORTRAN ((rate of [Act]) * d_i — the temp-variable form, used
    # only for cells where d_i > Int or the cell is a knot).
    rd_inh = -params.k_deg * inh
    secrete_inh = knot | (diff >= params.k_int)
    if disc.eq17_inh_source == 'act_concentration':
        rd_inh = np.where(secrete_inh, act - params.k_deg * inh, rd_inh)
    elif disc.eq17_inh_source == 'act_rate_times_di':
        rd_inh = np.where(secrete_inh, rd_act * diff - params.k_deg * inh, rd_inh)
    else:
        raise ValueError(f'unknown eq17_inh_source: {disc.eq17_inh_source!r}')

    # --- Eq. 18 secondary-signal source ------------------------------
    # disc.eq18_sec_source: paper (constant k_sec) vs FORTRAN (k_sec
    # ramped by d_i). The d_i factor smooths Sec onset; not in paper.
    rd_sec = -params.k_deg * sec
    secrete_sec = knot | (diff >= params.k_set)
    if disc.eq18_sec_source == 'constant_k_sec':
        rd_sec = np.where(secrete_sec, params.k_sec - params.k_deg * sec, rd_sec)
    elif disc.eq18_sec_source == 'k_sec_times_di':
        rd_sec = np.where(secrete_sec, params.k_sec * diff - params.k_deg * sec, rd_sec)
    else:
        raise ValueError(f'unknown eq18_sec_source: {disc.eq18_sec_source!r}')

    # Diffusion (computed once against current state).
    lap_act = mesh.laplacian(act)
    lap_inh = mesh.laplacian(inh)
    lap_sec = mesh.laplacian(sec)

    # Forward Euler update — Jacobi (single state copy).
    state.act = act + dt * (rd_act + params.k_da * lap_act)
    state.inh = inh + dt * (rd_inh + params.k_di * lap_inh)
    state.sec = sec + dt * (rd_sec + params.k_ds * lap_sec)

    # Concentrations are non-negative.
    np.maximum(state.act, 0.0, out=state.act)
    np.maximum(state.inh, 0.0, out=state.inh)
    np.maximum(state.sec, 0.0, out=state.sec)

    # --- Knot detection ----------------------------------------------
    # disc.knot_threshold_gate: paper (any cell with [Act] >= 1) vs
    # FORTRAN (only cells with index >= state.first_border_cell).
    # Knot status is irreversible: once True, always True.
    if disc.knot_threshold_gate == 'none':
        state.knot |= state.act >= 1.0
    elif disc.knot_threshold_gate == 'first_border_cell':
        idx = np.arange(state.act.shape[0])
        state.knot |= (state.act >= 1.0) & (idx >= state.first_border_cell)
    else:
        raise ValueError(f'unknown knot_threshold_gate: {disc.knot_threshold_gate!r}')
