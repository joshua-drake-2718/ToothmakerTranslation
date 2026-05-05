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
      and the canonical form. Path B follows the FORTRAN per charter.
      See docs/research/paper-review-2010-salazar-ciudad-jernvall.md
      §'Paper-vs-code findings'.

  d[Inh]/dt = -k_deg * [Inh] + k_di * lap[Inh]                   (eq. 15)
       (default for non-secreting cells)

  d[Inh]/dt =  [Act] - k_deg * [Inh] + k_di * lap[Inh]           (eq. 17)
       (cells where d_i >= k_int OR cell is a knot)

  d[Sec]/dt = -k_deg * [Sec] + k_ds * lap[Sec]                   (eq. 16)
       (default for non-secreting cells)

  d[Sec]/dt =  k_sec - k_deg * [Sec] + k_ds * lap[Sec]           (eq. 18)
       (cells where d_i >= k_set OR cell is a knot)

Knot detection rule (paper p. 586): "Differentiation is irreversible.
A cell becomes an enamel knot when [Act] >= 1." We apply this after the
Euler step against the new [Act] and store the boolean flag persistently
(once a knot, always a knot).

Path-B-vs-FORTRAN deviations recorded for posterity:

1. **Eq. 17 source term.** The FORTRAN computes Inh production as
   `(rate of change of Act) * d_i` for cells where d_i > Int — i.e. it
   uses the temporary Act-rate variable rather than [Act] itself. This
   is not the paper's eq. 17 form. Path B uses the paper's [Act].

2. **Eq. 18 d_i scaling.** The FORTRAN scales Sec production by d_i
   (`k_sec * d_i - k_deg * [Sec]`). The paper writes a constant rate
   (`k_sec - k_deg * [Sec]`). Path B uses the paper form.

3. **first_border_cell knot gate.** The FORTRAN restricts knot
   formation to cells with index >= first_border_cell (the outermost
   ring of the hexagonal lattice in `13.f90`). The paper has no such
   restriction. Path B allows knots at any cell where [Act] >= 1.

4. **Negative-rate clamping.** The FORTRAN clamps individual reaction
   rates to zero before integration. Path B does not pre-clamp; it
   clamps concentrations after the Euler step. Mathematically these
   are different operators, but for the parameters in seal.txt they
   converge to the same equilibrium (the reaction terms are positive
   in steady state).

If validation against the seal goldens fails because of these, each
deviation is a candidate for a `docs/findings/` entry and possibly a
charter amendment.
"""
from __future__ import annotations

import numpy as np

from .params import Params
from .state import State
from .mesh import Mesh


def step_reaction_diffusion(
    state: State,
    params: Params,
    mesh: Mesh,
    dt: float,
) -> None:
    """Apply one forward-Euler step to act/inh/sec in-place.

    Pure Jacobi: rates computed against current state, applied once.

    Diffusion uses the length-weighted graph Laplacian on the triangular
    mesh (see `mesh.Mesh`). Reaction terms follow the 2010 paper with
    the FORTRAN's denominator correction in eq. 14.
    """
    act = state.act
    inh = state.inh
    sec = state.sec
    diff = state.diff
    knot = state.knot

    # Reaction terms (state functions, not deltas).
    # Eq. 14 with FORTRAN denominator. The Hill-style numerator is
    # k_act * [Act] (linear); saturation is via the (1 + k_inh*[Inh])
    # denominator. This is autocatalysis with feedback inhibition by Inh.
    rd_act = params.k_act * act / (1.0 + params.k_inh * inh) - params.k_deg * act

    # Eq. 15 default + Eq. 17 override.
    rd_inh = -params.k_deg * inh
    secrete_inh = knot | (diff >= params.k_int)
    rd_inh = np.where(secrete_inh, act - params.k_deg * inh, rd_inh)

    # Eq. 16 default + Eq. 18 override.
    rd_sec = -params.k_deg * sec
    secrete_sec = knot | (diff >= params.k_set)
    rd_sec = np.where(secrete_sec, params.k_sec - params.k_deg * sec, rd_sec)

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

    # Knot detection: irreversible.
    state.knot |= state.act >= 1.0
