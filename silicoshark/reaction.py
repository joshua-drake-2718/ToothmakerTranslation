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
        13.f90 literal variant uses [Act] * d_i in place of [Act];
        Path A's coreop2d.py rewrite uses (rate of [Act]) * d_i)

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
    # disc.eq17_inh_source: three forms in play, used only for cells
    # where d_i > Int or the cell is a knot.
    #   - 'act_concentration': paper Eq. 17 — source is plain [Act].
    #   - 'act_times_di': 13.f90:487 literal — source is
    #       q3d(i,1,1) * DiffState(i), i.e. [Act] concentration * d_i.
    #     Humppa's humppa_translate.f90:617 matches.
    #   - 'act_rate_times_di': Path A's coreop2d.py:647 rewrite —
    #       source is hq3d[i, 0, 0] * diff_state[i], i.e. the
    #       (post-eq.14) Act *rate* * d_i. NOT a faithful translation
    #       of 13.f90; the A8 audit identified this divergence.
    rd_inh = -params.k_deg * inh
    secrete_inh = knot | (diff >= params.k_int)
    if disc.eq17_inh_source == 'act_concentration':
        rd_inh = np.where(secrete_inh, act - params.k_deg * inh, rd_inh)
    elif disc.eq17_inh_source == 'act_times_di':
        rd_inh = np.where(secrete_inh, act * diff - params.k_deg * inh, rd_inh)
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

    # --- Vertical coupling to mesenchyme (eqs. 15, 16) ---------------
    # disc.mesenchyme: 'absent' (epithelium-only; no vertical coupling)
    # vs 'per_column_z_layers' (paper-faithful: each epithelial cell
    # couples vertically to mes_inh[:, 0] / mes_sec[:, 0] beneath it).
    # Act has NO vertical diffusion (paper: 'mesenchymal cells have no
    # reaction terms' and Act is not produced or transmitted there).
    # The vertical flux (mes[0] - epi) is added to the epithelial
    # Laplacian; the mesenchyme side of the flux (epi - mes[0]) is
    # applied in `step_mesenchyme_diffusion` so total mass is conserved.
    if disc.mesenchyme == 'per_column_z_layers':
        lap_inh = lap_inh + (state.mes_inh[:, 0] - inh)
        lap_sec = lap_sec + (state.mes_sec[:, 0] - sec)
    elif disc.mesenchyme != 'absent':
        raise ValueError(f'unknown mesenchyme: {disc.mesenchyme!r}')

    # Forward Euler update — Jacobi (single state copy).
    state.act = act + dt * (rd_act + params.k_da * lap_act)
    state.inh = inh + dt * (rd_inh + params.k_di * lap_inh)
    state.sec = sec + dt * (rd_sec + params.k_ds * lap_sec)

    # Concentrations are non-negative.
    np.maximum(state.act, 0.0, out=state.act)
    np.maximum(state.inh, 0.0, out=state.inh)
    np.maximum(state.sec, 0.0, out=state.sec)

    # --- Mesenchyme diffusion + degradation (eqs. 15, 16) ------------
    # No reaction terms in mesenchyme (paper p. 586). Each layer ℓ
    # accumulates horizontal Laplacian (across epithelial columns at
    # the same depth) plus vertical coupling to neighbouring layers
    # (top layer couples up to epi, bottom layer is reflective at the
    # substrate). Vertical flux on the layer-0 ↔ epi interface is
    # symmetric to the term added to `lap_inh` / `lap_sec` above
    # (using the SAME start-of-step `inh` / `sec` snapshots), so
    # total mass (epi + mesenchyme columns) is conserved under k_deg=0.
    step_mesenchyme_diffusion(
        state, params, mesh, dt, disc,
        pre_step_epi_inh=inh, pre_step_epi_sec=sec,
    )

    # --- Knot detection ----------------------------------------------
    # disc.knot_threshold_gate: paper (any cell with [Act] >= 1) vs
    # FORTRAN (only interior cells, identified by the topological
    # `first_border_cell` gate).
    #
    # In FORTRAN's reversed lattice the gate spells `i >= first_border_cell`
    # to select interior cells; silicoshark's lattice has interior cells
    # at LOW indices (centre + early rings), so the same semantic gate
    # is `idx < state.first_border_cell` here. Knots form in the
    # interior (the centre of the lattice), not on the border.
    # Knot status is irreversible: once True, always True.
    if disc.knot_threshold_gate == 'none':
        state.knot |= state.act >= 1.0
    elif disc.knot_threshold_gate == 'first_border_cell':
        idx = np.arange(state.act.shape[0])
        state.knot |= (state.act >= 1.0) & (idx < state.first_border_cell)
    else:
        raise ValueError(f'unknown knot_threshold_gate: {disc.knot_threshold_gate!r}')


def step_mesenchyme_diffusion(
    state: State,
    params: Params,
    mesh: Mesh,
    dt: float,
    disc: Discretisation = PATH_B_DEFAULT,
    *,
    pre_step_epi_inh: np.ndarray | None = None,
    pre_step_epi_sec: np.ndarray | None = None,
) -> None:
    """Forward-Euler update of `state.mes_inh` and `state.mes_sec` for
    each mesenchymal z-layer.

    Per the 2010 paper p. 586, mesenchymal cells have NO reaction
    terms; their dynamics are diffusion + first-order degradation
    only (eqs. 15, 16 with the Eq. 17 / 18 source terms zero). For
    each layer ℓ:

        d[mes_inh][:, ℓ]/dt = -k_deg * mes_inh[:, ℓ]
                              + k_di * (lap_h[ℓ] + lap_v[ℓ])
        d[mes_sec][:, ℓ]/dt = -k_deg * mes_sec[:, ℓ]
                              + k_ds * (lap_h[ℓ] + lap_v[ℓ])

    where `lap_h[ℓ]` is the lateral Laplacian on the (N,) array of
    layer-ℓ concentrations (same `mesh.laplacian` operator as the
    epithelial layer), and `lap_v[ℓ]` is the 1D vertical Laplacian
    along the column with reflective ghost layers at the top
    (epithelium boundary) and the bottom (substrate boundary). The
    epithelium ↔ layer-0 flux is symmetric to the term added to
    `lap_inh` / `lap_sec` in `step_reaction_diffusion`, so total
    mass (epi + all mesenchyme layers) is conserved under k_deg = 0
    when the same `pre_step_epi_inh` / `pre_step_epi_sec` snapshot
    is used in both halves.

    No-op when `disc.mesenchyme == 'absent'`.

    Vertical Laplacian form (zero-flux at substrate, coupling to
    epithelium at the top):

        layer 0 (top, against epi):
            lap_v[0] = (epi - mes_inh[:, 0]) + (mes_inh[:, 1] - mes_inh[:, 0])

        layer 1 .. n_mes-2 (interior; absent for N_MESENCHYME_LAYERS=2):
            lap_v[ℓ] = (mes_inh[:, ℓ-1] - mes_inh[:, ℓ])
                       + (mes_inh[:, ℓ+1] - mes_inh[:, ℓ])

        layer n_mes-1 (bottom, against substrate):
            lap_v[n-1] = (mes_inh[:, n-2] - mes_inh[:, n-1])
                         + 0  (zero-flux at substrate)

    Parameters
    ----------
    pre_step_epi_inh, pre_step_epi_sec : optional (N,) arrays
        Pre-update epithelial concentrations to use in the vertical
        flux into layer 0. When called from `step_reaction_diffusion`
        these are the start-of-step `inh` / `sec` snapshots so the
        epi/mesenchyme pair shares a consistent value (mass-conserving).
        When called standalone these default to the current
        `state.inh` / `state.sec` (small drift in mass conservation;
        acceptable for diagnostic / unit-test calls).
    """
    if disc.mesenchyme == 'absent':
        return
    if disc.mesenchyme != 'per_column_z_layers':
        raise ValueError(f'unknown mesenchyme: {disc.mesenchyme!r}')

    n_layers = state.mes_inh.shape[1]
    if n_layers < 1:
        return

    # Snapshot current mesenchyme state (Jacobi: rates against
    # start-of-step).
    mes_inh = state.mes_inh.copy()
    mes_sec = state.mes_sec.copy()
    epi_inh = pre_step_epi_inh if pre_step_epi_inh is not None else state.inh.copy()
    epi_sec = pre_step_epi_sec if pre_step_epi_sec is not None else state.sec.copy()

    # Horizontal Laplacian per layer (one mesh-laplacian call per
    # layer). The mesh's neighbour graph is the epithelial neighbour
    # graph; for the per-column layout the same graph applies to
    # every depth (mesenchymal cells are attached to their epithelial
    # column; lateral neighbours are the underlying mesenchyme of
    # adjacent epithelial columns).
    lap_h_inh = np.empty_like(mes_inh)
    lap_h_sec = np.empty_like(mes_sec)
    for ell in range(n_layers):
        lap_h_inh[:, ell] = mesh.laplacian(mes_inh[:, ell])
        lap_h_sec[:, ell] = mesh.laplacian(mes_sec[:, ell])

    # Vertical Laplacian per layer with zero-flux substrate boundary
    # and explicit coupling to the epithelium at layer 0.
    lap_v_inh = np.zeros_like(mes_inh)
    lap_v_sec = np.zeros_like(mes_sec)
    if n_layers == 1:
        # Single layer: couples up to epi and down to reflective
        # substrate (substrate term is 0).
        lap_v_inh[:, 0] = epi_inh - mes_inh[:, 0]
        lap_v_sec[:, 0] = epi_sec - mes_sec[:, 0]
    else:
        # Top layer (couples up to epi, down to layer 1).
        lap_v_inh[:, 0] = (epi_inh - mes_inh[:, 0]) + (mes_inh[:, 1] - mes_inh[:, 0])
        lap_v_sec[:, 0] = (epi_sec - mes_sec[:, 0]) + (mes_sec[:, 1] - mes_sec[:, 0])
        # Interior layers (none for N_MESENCHYME_LAYERS=2; loop kept
        # so future n_layers > 2 configurations are correct).
        for ell in range(1, n_layers - 1):
            lap_v_inh[:, ell] = (
                (mes_inh[:, ell - 1] - mes_inh[:, ell])
                + (mes_inh[:, ell + 1] - mes_inh[:, ell])
            )
            lap_v_sec[:, ell] = (
                (mes_sec[:, ell - 1] - mes_sec[:, ell])
                + (mes_sec[:, ell + 1] - mes_sec[:, ell])
            )
        # Bottom layer (couples up to layer n-2; zero flux at substrate).
        last = n_layers - 1
        lap_v_inh[:, last] = mes_inh[:, last - 1] - mes_inh[:, last]
        lap_v_sec[:, last] = mes_sec[:, last - 1] - mes_sec[:, last]

    # Forward-Euler update. Degradation acts on the layer's own value
    # (no production source).
    state.mes_inh = mes_inh + dt * (
        -params.k_deg * mes_inh + params.k_di * (lap_h_inh + lap_v_inh)
    )
    state.mes_sec = mes_sec + dt * (
        -params.k_deg * mes_sec + params.k_ds * (lap_h_sec + lap_v_sec)
    )

    np.maximum(state.mes_inh, 0.0, out=state.mes_inh)
    np.maximum(state.mes_sec, 0.0, out=state.mes_sec)
