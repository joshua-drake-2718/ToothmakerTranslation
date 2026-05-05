"""Mechanical forces on cell centres (eqs. 1–13 of the 2010 paper).

All forces are computed against the start-of-iteration positions and
accumulated into a per-cell (N, 3) array of position deltas. Pure
Jacobi: a single state copy is read; a single update is applied per
step.

Equations implemented:

  Repulsion (eq. 1): f_ij = k_rep * (|d| - p0) * d
       Active when |d| < p0 (inward push). p0 = state.rest_length = 1.

  Adhesion (eq. 2): f_ij = k_adh * d / |d|
       Active when |d| > p0 (unit-vector pull).

  Nucleus traction (eq. 4):
       dp_i/dt = k_ntr * (1 - d_i) * (mean_neigh - p_i)
       Pulls cell toward the centroid of its mesh neighbours, gated
       by (1 - d_i) (knot/differentiated cells freeze).

  Epithelial growth (eq. 5):
       dp_i/dt = k_egr * (1 - d_i) * sum(u_ij) / |sum(u_ij)|
       u_ij = (p_j - p_i) / |p_j - p_i|.
       For a symmetric-neighbour cell this is zero; for a border
       cell it is non-zero in the direction of the bulk.

  Cervical-loop downgrowth (eqs. 7–9 with eqs. 10–12 mesenchyme):
       Border cells only. Outward horizontal direction
       o = (p_i - mean_neigh)_xy, normalised.
       Mesenchyme thickening: (|o| + k_mgr * [Sec] + k_umgr) / |o|.
       Vertical component: k_dgr (positive = into mesenchyme).
       Then scale by k_egr * (1 - d_i).

  Buoyancy (eq. 13):
       dp_i/dt -= k_boy * [Sec] * (1 - d_i) * apical_normal
       Apical = upward in the FORTRAN's z convention. For seal.txt
       [Sec] = 0 throughout the run, so this is a no-op for the
       seal validation, but implemented for completeness.

Documented Path-B-vs-FORTRAN deviations:

1. **Eq. 5 not gated by neighbour-z.** The FORTRAN restricts the
   sum-of-unit-vectors to neighbours that lie above in z (`b < -1e-4`).
   This is a 3D outward filter that the paper does not specify. Path
   B uses the symmetric form: every mesh neighbour contributes.

2. **Eq. 5 applied to all cells.** The FORTRAN restricts eq. 5 to
   inner cells (`i >= first_border_cell`); border cells receive only
   the cervical-loop force. The paper has no such restriction.
   Path B applies eq. 5 to all cells, in addition to the cervical-loop
   force on border cells.

3. **Adhesion (eq. 2) on all cells.** The FORTRAN restricts adhesion
   to inner cells. The paper has no such restriction.

4. **Nucleus traction (1 - d_i) gate on all components.** The FORTRAN
   gates only the z component. Path B follows the paper.

5. **Border bias multiplier semantics.** The FORTRAN's update_cell_position
   multiplies x-force by `Pbi`/`Abi` (depending on x sign) and z-force
   by `Bgr` for cells in the y-band |y| < Bwi (excluding x exactly = 0).
   For seal.txt this is non-trivial: Pbi = Abi = 0 zeroes the x force,
   and Bgr = 10 amplifies the z force on the y-axis-adjacent rim. Path
   B retains this multiplier behaviour because it is load-bearing for
   the 57-cell plateau.

If validation diverges by more than ±5%, each item is a candidate for
a `docs/findings/` entry.
"""
from __future__ import annotations

import numpy as np

from .params import Params
from .state import State
from .mesh import Mesh


_SAFE_DIV_EPS = 1e-12


def _safe_norm(v: np.ndarray) -> np.ndarray:
    """Norm along last axis, clipped at EPS to avoid divide-by-zero."""
    return np.maximum(np.linalg.norm(v, axis=-1), _SAFE_DIV_EPS)


def _edge_arrays(state: State, mesh: Mesh) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Build (rows, cols, d, dist) arrays from mesh CSR adjacency.

    rows[e], cols[e] are the two endpoint indices of edge e.
    d[e] = positions[cols[e]] - positions[rows[e]].
    dist[e] = |d[e]|.
    """
    n = state.num_active
    rows = np.repeat(
        np.arange(n, dtype=np.int64),
        np.diff(mesh.neigh_starts),
    )
    cols = mesh.neigh_idx
    d = state.positions[cols] - state.positions[rows]
    dist = np.linalg.norm(d, axis=1)
    return rows, cols, d, dist


def compute_forces(state: State, params: Params, mesh: Mesh) -> np.ndarray:
    """Total mechanical force on each cell, before Pbi/Abi/Bgr scaling.

    Returns (N, 3) array of dp/dt contributions to be multiplied by dt
    in `simulator.step` and applied to positions.

    The Pbi/Abi/Bgr multipliers from FORTRAN's update_cell_position are
    NOT applied here; they are applied in `apply_border_multipliers` as
    a separate step so the bias logic stays inspectable.
    """
    n = state.num_active
    pos = state.positions
    rows, cols, d, dist = _edge_arrays(state, mesh)

    forces = np.zeros((n, 3), dtype=np.float64)
    one_minus_d = (1.0 - state.diff)[:, None]

    # --- Repulsion (eq. 1) -------------------------------------------
    # Always-on Hookean spring per the paper: f_ij = k_rep * (|d| - p0) * d
    # smoothly transitions from repulsion (|d| < p0) to weak attraction
    # (|d| > p0). The FORTRAN's `c = min(Rep, 1)` clamp is preserved —
    # for seal.txt Rep = 1.5, so effective stiffness is 1.0. At |d| = p0
    # exactly the force is zero, which avoids the floating-point noise
    # near equilibrium that an `|d| < p0` gate would amplify.
    rest = state.rest_length
    rep_clip = min(params.k_rep, 1.0)
    rep_scalar = rep_clip * (dist - rest)
    f_rep = rep_scalar[:, None] * d
    np.add.at(forces, rows, f_rep)

    # --- Adhesion (eq. 2) --------------------------------------------
    # Additive on top of eq. 1, only when truly stretched. Use a small
    # tolerance band (1e-6 * rest) to avoid spurious activation at
    # equilibrium. Unit-vector pull, weighted by k_adh.
    stretch_tol = rest * 1e-6
    outside = dist > rest + stretch_tol
    adh_scalar = np.where(outside, params.k_adh / np.maximum(dist, _SAFE_DIV_EPS), 0.0)
    f_adh = adh_scalar[:, None] * d
    np.add.at(forces, rows, f_adh)

    # --- Nucleus traction (eq. 4) ------------------------------------
    counts = np.diff(mesh.neigh_starts).astype(np.float64)
    sum_neigh = np.zeros((n, 3), dtype=np.float64)
    np.add.at(sum_neigh, rows, pos[cols])
    safe_counts = np.maximum(counts, 1.0)[:, None]
    mean_neigh = sum_neigh / safe_counts
    forces += params.k_ntr * one_minus_d * (mean_neigh - pos)

    # --- Epithelial growth (eq. 5) -----------------------------------
    # f_i = k_egr (1-d_i) * sum_j(unit_ij) / |sum_j(unit_ij)|, summed
    # only over neighbours with strictly higher z than cell i. The
    # z-gate matches the FORTRAN's `b < -1e-4` condition in
    # ep_growth_border_force; without it, the symmetric-flat-z initial
    # condition produces near-zero sums whose normalised direction
    # amplifies floating-point noise into spurious inward forces.
    #
    # As cells curve under cervical-loop downgrowth, the inner cells'
    # neighbours sit higher in z and the gate begins to fire,
    # propagating apical motion through the bulk.
    z_diff = state.positions[cols, 2] - state.positions[rows, 2]
    z_gate = z_diff > 1e-4
    unit_ij = d / np.maximum(dist[:, None], _SAFE_DIV_EPS)
    gated = unit_ij * z_gate[:, None]
    sum_unit = np.zeros((n, 3), dtype=np.float64)
    np.add.at(sum_unit, rows, gated)
    sum_unit_norm = _safe_norm(sum_unit)
    egr_dir = sum_unit / sum_unit_norm[:, None]
    has_dir = sum_unit_norm > _SAFE_DIV_EPS
    forces[has_dir] += (
        params.k_egr * one_minus_d[has_dir] * egr_dir[has_dir]
    )

    # --- Cervical-loop downgrowth (eqs. 7–9 + 10–12) -----------------
    # Border cells get an additional outward+downward push from the
    # cervical-loop dynamics (paper Methods, p. 585). This is ADDITIVE
    # on top of eq. 5 — the paper's eqs. 7–9 are described as a
    # supplementary force on border cells, not as a replacement for the
    # general epithelial growth.
    #
    # Outward direction in xy: -(mean_neigh - p_i)_xy, normalised.
    # (Pointing away from the bulk of neighbours.)
    border = mesh.is_border
    if np.any(border):
        outward_xy = -(mean_neigh[border] - pos[border])[:, :2]
        out_h = np.maximum(np.linalg.norm(outward_xy, axis=1), _SAFE_DIV_EPS)
        # Mesenchyme thickening factor (eqs. 10–12): the underlying
        # mesenchyme depth grows with [Sec]; thicker mesenchyme means
        # more lateral extension before z descent.
        sec_b = state.sec[border]
        thicken = (out_h + params.k_mgr * sec_b + params.k_umgr) / out_h
        outward_xy = outward_xy * thicken[:, None]
        # Vertical: k_dgr units (positive z = into mesenchyme in FORTRAN
        # convention, preserved here).
        cc = np.full(border.sum(), params.k_dgr, dtype=np.float64)
        # Combine and normalise to k_egr * (1 - d_i).
        out3 = np.column_stack([outward_xy[:, 0], outward_xy[:, 1], cc])
        out_mag = _safe_norm(out3)
        scale = params.k_egr / out_mag
        gate = np.maximum(1.0 - state.diff[border], 0.0)
        f_cl = out3 * (scale * gate)[:, None]
        forces[border] += f_cl

    # --- Buoyancy (eq. 13) -------------------------------------------
    # Subtract apical-pointing component scaled by Sec * (1 - d_i).
    # Apical normal: simplest discrete approximation = unit z-axis
    # (positive z = into mesenchyme; apical = negative z). Apply only
    # where Sec > 0 (no-op for the seal example).
    sec_active = state.sec > 0.0
    if np.any(sec_active):
        # Force component currently in the (apical = -z) direction.
        # Subtract k_boy * [Sec] * (1 - d_i) of the apical-pointing
        # part. We approximate apical_normal = (0, 0, -1).
        boy_factor = params.k_boy * state.sec * np.maximum(1.0 - state.diff, 0.0)
        # Subtract from the z-component (push back toward mesenchyme).
        forces[sec_active, 2] -= -boy_factor[sec_active]
        # ^ Equivalent to forces[..., 2] += boy_factor — buoyancy resists
        # apical motion by pushing back into the mesenchyme.

    return forces


def apply_border_multipliers(
    forces: np.ndarray,
    state: State,
    params: Params,
    mesh: Mesh,
) -> None:
    """In-place application of FORTRAN's update_cell_position bias rules.

    For border cells whose y is in the band |y| < Bwi:
      x > 0  -> forces[i, 0] *= k_pbi;  forces[i, 2] *= k_bgr
      x < 0  -> forces[i, 0] *= k_abi;  forces[i, 2] *= k_bgr
      x == 0 (exact float zero) -> no change (FORTRAN quirk preserved
              for golden compatibility — see coreop2d.update_cell_position
              docstring; seal.txt's y-axis cells start at x = 0 exactly).

    Then any cell with negative z-force or knot status has its z-force
    zeroed (FORTRAN's "due to the pressure of the stelate" comment).
    """
    border = mesh.is_border
    pos = state.positions
    in_band = np.abs(pos[:, 1]) < params.k_bwi
    east = pos[:, 0] > 0
    west = pos[:, 0] < 0
    apply_pbi = border & in_band & east
    apply_abi = border & in_band & west
    apply_bgr = border & in_band & (east | west)

    forces[apply_pbi, 0] *= params.k_pbi
    forces[apply_abi, 0] *= params.k_abi
    forces[apply_bgr, 2] *= params.k_bgr

    # Z-clamp: knots and negative z-force become zero.
    z_clamp = (forces[:, 2] < 0) | state.knot
    forces[z_clamp, 2] = 0.0
