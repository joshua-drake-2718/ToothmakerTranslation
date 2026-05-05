"""Mechanical forces on cell centres (eqs. 1–13 of the 2010 paper).

All forces are computed against the start-of-iteration positions and
accumulated into a per-cell (N, 3) array of position deltas. Pure
Jacobi: a single state copy is read; a single update is applied per
step.

Equations implemented:

  Repulsion (eq. 1): f_ij = k_rep * (|d| - p0) * d
       Default form is `hookean_signed` — smoothly signed across rest
       length. `paper_gated` restricts the term to |d| < p0 (active
       only as repulsion; leaves the |d| > p0 case to adhesion alone).

  Adhesion (eq. 2): f_ij = k_adh * d / |d|
       Default form is `unit_vector` per the paper. `hookean_attraction`
       is the FORTRAN's else-branch: f = k_adh * d (no division by |d|),
       a Hookean attractive spring rather than a unit-vector pull. Both
       forms gate on |d| > rest + tol so they don't fire at exact-rest
       equilibrium.

  Nucleus traction (eq. 4):
       dp_i/dt = k_ntr * (1 - d_i) * (mean_neigh - p_i)

  Epithelial growth (eq. 5):
       dp_i/dt = k_egr * (1 - d_i) * sum(u_ij) / |sum(u_ij)|
       u_ij = (p_j - p_i) / |p_j - p_i|. Default: gated to neighbours
       strictly above in z (`eq5_z_gate=True`), and applied to all cells
       (`eq5_apply_to='all'`). The FORTRAN restricts the sum to
       above-z neighbours and the term to interior cells (index >=
       state.first_border_cell). Note that 'interior_only' uses the
       index gate (state.first_border_cell), NOT mesh.is_border —
       these are distinct semantics: the latter is geometric (computed
       each step from the triangulation), the former is FORTRAN's
       topological book-keeping (set at construction from the initial
       lattice).

  Cervical-loop downgrowth (eqs. 7–9 with eqs. 10–12 mesenchyme):
       Border cells only (mesh.is_border).

  Buoyancy (eq. 13):
       dp_i/dt -= k_boy * [Sec] * (1 - d_i) * apical_normal
       Apical = upward in the FORTRAN's z convention. For seal.txt
       [Sec] = 0 throughout the run, so this is a no-op for the
       seal validation, but implemented for completeness.

The `apply_border_multipliers` post-process applies FORTRAN's
update_cell_position bias rules. The `border_bias_x_zero_quirk`
field controls whether cells with `x` exactly == 0 are skipped
(FORTRAN behaviour) or receive the multiplier on the basis of
their x sign (the corrected behaviour).

Path B v2 wires every paper-vs-FORTRAN choice through the
`Discretisation` dataclass (charter
`docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md`).
The default `PATH_B_DEFAULT` reproduces the v1 behaviour byte-for-byte
on the seal example.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial import cKDTree

from .params import Params
from .state import State
from .mesh import Mesh
from .discretisation import Discretisation, PATH_B_DEFAULT


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


def compute_forces(
    state: State,
    params: Params,
    mesh: Mesh,
    disc: Discretisation = PATH_B_DEFAULT,
) -> np.ndarray:
    """Total mechanical force on each cell, before Pbi/Abi/Bgr scaling.

    Returns (N, 3) array of dp/dt contributions to be multiplied by dt
    in `simulator.step` and applied to positions.

    The Pbi/Abi/Bgr multipliers from FORTRAN's update_cell_position are
    NOT applied here; they are applied in `apply_border_multipliers` as
    a separate step so the bias logic stays inspectable.

    `disc` selects the Path B v2 implementer-choice configuration. Each
    decision point is a single short conditional with a one-line comment
    naming the field, so a reader scanning the function sees every
    option exercised.
    """
    n = state.num_active
    pos = state.positions
    rows, cols, d, dist = _edge_arrays(state, mesh)

    forces = np.zeros((n, 3), dtype=np.float64)
    one_minus_d = (1.0 - state.diff)[:, None]
    rest = state.rest_length
    rep_clip = min(params.k_rep, 1.0)

    # --- Repulsion (eq. 1) -------------------------------------------
    # disc.rep_form: paper-gated (active only when |d| < rest) vs
    # hookean_signed (smoothly signed across rest length, the v1 form).
    if disc.rep_form == 'hookean_signed':
        rep_scalar = rep_clip * (dist - rest)
    elif disc.rep_form == 'paper_gated':
        rep_scalar = np.where(dist < rest, rep_clip * (dist - rest), 0.0)
    else:
        raise ValueError(f'unknown rep_form: {disc.rep_form!r}')
    f_rep = rep_scalar[:, None] * d
    np.add.at(forces, rows, f_rep)

    # --- Adhesion (eq. 2) --------------------------------------------
    # disc.adh_form: unit_vector (paper, f = k_adh * d / |d|) vs
    # hookean_attraction (FORTRAN, f = k_adh * d). Both gate on
    # |d| > rest + tol to avoid firing at exact-rest equilibrium.
    stretch_tol = rest * 1e-6
    outside = dist > rest + stretch_tol
    if disc.adh_form == 'unit_vector':
        adh_scalar = np.where(outside, params.k_adh / np.maximum(dist, _SAFE_DIV_EPS), 0.0)
        f_adh = adh_scalar[:, None] * d
    elif disc.adh_form == 'hookean_attraction':
        f_adh = np.where(outside[:, None], params.k_adh * d, 0.0)
    else:
        raise ValueError(f'unknown adh_form: {disc.adh_form!r}')
    np.add.at(forces, rows, f_adh)

    # --- Repulsion against close non-mesh cells (FORTRAN pushingnovei) -
    # disc.rep_neighbour_set: 'mesh' (paper — eq. 1 sums only over
    # mesh neighbours) vs 'mesh_plus_all_close' (FORTRAN's extra
    # `pushingnovei` loop over any cell within ~1.4 of `i` that is
    # not already a mesh neighbour). The non-mesh form follows
    # FORTRAN's polynomial form (coreop2d.py:920-925):
    #   dd = 1 / (d + 1)**8;  f = -d_vec * (dd / d) * Rep
    # rather than the simpler Hookean: 13.f90 actually applies this
    # polynomial profile, and the strong short-range repulsion is the
    # mechanism that prevents mesh-tangling, which is plausibly load-
    # bearing for stability. The Hookean alternative would be the
    # same `(|d|-rest) * d_vec * rep_clip` we use for mesh edges; we
    # keep the FORTRAN form so `LEGACY_FORTRAN` reproduces the
    # FORTRAN's stability profile.
    if disc.rep_neighbour_set == 'mesh':
        pass
    elif disc.rep_neighbour_set == 'mesh_plus_all_close':
        # Build a set of mesh edge pairs (sorted endpoints) for filter.
        mesh_pairs = set()
        for e in range(rows.shape[0]):
            i = int(rows[e])
            j = int(cols[e])
            a, b = (i, j) if i < j else (j, i)
            mesh_pairs.add((a, b))
        tree = cKDTree(pos)
        # FORTRAN's gate is `d < 1.4` in the absolute coordinate sense.
        close_pairs = tree.query_pairs(r=1.4, output_type='set')
        for a, b in close_pairs:
            if (a, b) in mesh_pairs:
                continue
            # FORTRAN computes both directions independently in
            # repel_non_neigh: cell i's loop over ii reads
            # `ux = positions[ii] - positions[i]` and accumulates
            # `-ux * (dd/d) * Rep` on cell i. Symmetrically apply
            # to both endpoints.
            d_ab = pos[b] - pos[a]
            mag = float(np.linalg.norm(d_ab))
            if mag < 1e-15 or mag >= 1.4:
                continue
            # FORTRAN's per-axis 1.4 gates (`if ux > 1.4: continue` etc.)
            # are signed-only checks on the positive side; the magnitude
            # gate is the active criterion. We reproduce only the
            # magnitude gate here; the per-axis gates are vestigial in
            # FORTRAN since after the magnitude check they cannot fire.
            dd = 1.0 / (mag + 1.0) ** 8
            scale = dd / mag
            # Cell a: f_a += -(p_b - p_a) * (dd/mag) * Rep
            # Cell b: f_b += -(p_a - p_b) * (dd/mag) * Rep = +d_ab * (dd/mag) * Rep
            forces[a, 0] -= d_ab[0] * scale * params.k_rep
            forces[a, 1] -= d_ab[1] * scale * params.k_rep
            forces[a, 2] -= d_ab[2] * scale * params.k_rep
            forces[b, 0] += d_ab[0] * scale * params.k_rep
            forces[b, 1] += d_ab[1] * scale * params.k_rep
            forces[b, 2] += d_ab[2] * scale * params.k_rep
    else:
        raise ValueError(f'unknown rep_neighbour_set: {disc.rep_neighbour_set!r}')

    # --- Nucleus traction (eq. 4) ------------------------------------
    counts = np.diff(mesh.neigh_starts).astype(np.float64)
    sum_neigh = np.zeros((n, 3), dtype=np.float64)
    np.add.at(sum_neigh, rows, pos[cols])
    safe_counts = np.maximum(counts, 1.0)[:, None]
    mean_neigh = sum_neigh / safe_counts
    forces += params.k_ntr * one_minus_d * (mean_neigh - pos)

    # --- Epithelial growth (eq. 5) -----------------------------------
    # disc.eq5_z_gate: paper-symmetric (False — every neighbour
    # contributes) vs FORTRAN-gated (True — only neighbours strictly
    # above in z, `b < -1e-4` in coreop2d.py:702). Default True per
    # v1 evidence: without the gate, symmetric flat-z initial
    # conditions produce near-zero sums whose unit-vector
    # normalisation amplifies floating-point noise.
    if disc.eq5_z_gate:
        z_diff = state.positions[cols, 2] - state.positions[rows, 2]
        z_gate = z_diff > 1e-4
        unit_ij = d / np.maximum(dist[:, None], _SAFE_DIV_EPS)
        gated = unit_ij * z_gate[:, None]
    else:
        unit_ij = d / np.maximum(dist[:, None], _SAFE_DIV_EPS)
        gated = unit_ij
    sum_unit = np.zeros((n, 3), dtype=np.float64)
    np.add.at(sum_unit, rows, gated)
    sum_unit_norm = _safe_norm(sum_unit)
    egr_dir = sum_unit / sum_unit_norm[:, None]
    has_dir = sum_unit_norm > _SAFE_DIV_EPS

    # disc.eq5_apply_to: paper ('all' — every cell) vs FORTRAN
    # ('interior_only' — only cells with index >= first_border_cell).
    # The 'interior_only' branch uses the FORTRAN's index gate,
    # NOT mesh.is_border: in 13.f90 eq. 5 fires for cells in the
    # range [first_border_cell, num_active_cells), which is a
    # topological book-keeping criterion, not a geometric one.
    if disc.eq5_apply_to == 'all':
        eligible = has_dir
    elif disc.eq5_apply_to == 'interior_only':
        idx = np.arange(n)
        eligible = has_dir & (idx >= state.first_border_cell)
    else:
        raise ValueError(f'unknown eq5_apply_to: {disc.eq5_apply_to!r}')
    forces[eligible] += (
        params.k_egr * one_minus_d[eligible] * egr_dir[eligible]
    )

    # --- Cervical-loop downgrowth (eqs. 7–9 + 10–12) -----------------
    # Border cells (mesh.is_border — geometric) get an additional
    # outward+downward push from the cervical-loop dynamics.
    border = mesh.is_border
    if np.any(border):
        outward_xy = -(mean_neigh[border] - pos[border])[:, :2]
        out_h = np.maximum(np.linalg.norm(outward_xy, axis=1), _SAFE_DIV_EPS)
        # Mesenchyme thickening factor (eqs. 10–12).
        sec_b = state.sec[border]
        thicken = (out_h + params.k_mgr * sec_b + params.k_umgr) / out_h
        outward_xy = outward_xy * thicken[:, None]
        cc = np.full(border.sum(), params.k_dgr, dtype=np.float64)
        out3 = np.column_stack([outward_xy[:, 0], outward_xy[:, 1], cc])
        out_mag = _safe_norm(out3)
        scale = params.k_egr / out_mag
        gate = np.maximum(1.0 - state.diff[border], 0.0)
        f_cl = out3 * (scale * gate)[:, None]
        forces[border] += f_cl

    # --- Buoyancy (eq. 13) -------------------------------------------
    sec_active = state.sec > 0.0
    if np.any(sec_active):
        boy_factor = params.k_boy * state.sec * np.maximum(1.0 - state.diff, 0.0)
        forces[sec_active, 2] -= -boy_factor[sec_active]

    return forces


def apply_border_multipliers(
    forces: np.ndarray,
    state: State,
    params: Params,
    mesh: Mesh,
    disc: Discretisation = PATH_B_DEFAULT,
) -> None:
    """In-place application of FORTRAN's update_cell_position bias rules.

    For border cells whose y is in the band |y| < Bwi:
      x > 0  -> forces[i, 0] *= k_pbi;  forces[i, 2] *= k_bgr
      x < 0  -> forces[i, 0] *= k_abi;  forces[i, 2] *= k_bgr
      x == 0 -> behaviour depends on `disc.border_bias_x_zero_quirk`:
        True (default) — no change (FORTRAN quirk; see
            coreop2d.update_cell_position docstring; seal.txt's y-axis
            cells start at x = 0 exactly and the FORTRAN's `if x>0`
            / `elif x<0` chain skips them).
        False — choose Pbi for x>=0, Abi for x<0, so x==0 cells are
            treated as belonging to the east half-plane and receive
            the Bgr z multiplier.

    Then any cell with negative z-force or knot status has its z-force
    zeroed (FORTRAN's "due to the pressure of the stelate" comment).
    """
    border = mesh.is_border
    pos = state.positions
    in_band = np.abs(pos[:, 1]) < params.k_bwi

    # disc.border_bias_x_zero_quirk: True preserves the FORTRAN
    # `if x>0 / elif x<0` chain that omits cells at x == 0; False
    # promotes x == 0 to the east half-plane (>= 0) so those cells
    # receive Pbi and Bgr. This is documented as a load-bearing
    # FORTRAN accident in coreop2d.py:1011-1023.
    if disc.border_bias_x_zero_quirk:
        east = pos[:, 0] > 0
        west = pos[:, 0] < 0
    else:
        east = pos[:, 0] >= 0
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
