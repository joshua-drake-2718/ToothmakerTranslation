"""Main iteration loop for silicoshark.

One forward-Euler step per iteration. The exact ordering of sub-steps
depends on `Discretisation.update_order`:

  - `'jacobi'` (default): all rates and forces are evaluated against
    start-of-step state; positions are updated once at the end. This
    is the simplest reading of the paper's PDE form, and matches v1.
        1. Build mesh (and laplacian operator) from current positions.
        2. Reaction-diffusion update (Act, Inh, Sec).
        3. Differentiation update.
        4. Compute mechanical forces against start-of-step state.
        5. Apply border-bias multipliers; z clamp.
        6. positions += dt * forces.
        7. Cell division: split qualifying edges.

  - `'gauss_seidel_forces'` / `'mixed_jacobi_gs'` (alias): matches the
    FORTRAN's `iteration()` ordering at coreop2d.py:1632-1644. The
    load-bearing distinction relative to Jacobi is that nuclear
    traction's position update is applied BEFORE the remaining force
    components are computed, so the latter see the post-traction
    positions. Finer-grained per-force Gauss-Seidel ordering is a
    refinement deferred to A6/A7.

`Discretisation.topology` selects between two adjacency strategies:

  - `'delaunay_each_step'`: rebuild via `scipy.spatial.Delaunay` every
    iteration. Path B v1 default; preserved for comparison.
  - `'static_with_local_update'`: build once via Delaunay at t=0, then
    locally rewire on cell division (FORTRAN's `add_cell` strategy
    minus the topology-walk accidents). Default for v2.

`dt` matches FORTRAN's delta = 0.05 by default. Iteration counts are
specified by the caller; FORTRAN goldens use 100 iters/save, 5 saves.
"""
from __future__ import annotations

import numpy as np

from .params import Params
from .state import State, N_MESENCHYME_LAYERS
from .mesh import Mesh
from .topology import Topology
from .reaction import step_reaction_diffusion
from .forces import compute_forces, apply_border_multipliers
from .discretisation import Discretisation, PATH_B_DEFAULT


DEFAULT_DT: float = 0.05
DIVISION_FACTOR: float = 2.0  # edge length threshold = DIVISION_FACTOR * rest


def build_topology(state: State) -> Topology:
    """Build a Topology from the current state's positions.

    Used at the head of `run` (and by callers wanting to construct a
    topology to pass back into `step`) when
    `Discretisation.topology == 'static_with_local_update'`.
    """
    return Topology.from_positions(state.positions)


def step_differentiation(
    state: State,
    params: Params,
    dt: float,
    disc: Discretisation = PATH_B_DEFAULT,
) -> None:
    """Update differentiation state via the eq. 6 accumulator.

    d_i_{t+1} = clip(d_i_t + dt * k_dff * X, 0, 1)

    where X is `[Sec]` (FORTRAN; humppa line 659; default per the
    Path B v1 charter) or `[Act]` (paper Eq. 6). The branch is
    selected by `disc.diff_accumulator`.
    """
    # disc.diff_accumulator: paper Eq. 6 (`act`) vs FORTRAN (`sec`).
    if disc.diff_accumulator == 'sec':
        state.diff += dt * params.k_dff * state.sec
    elif disc.diff_accumulator == 'act':
        state.diff += dt * params.k_dff * state.act
    else:
        raise ValueError(f'unknown diff_accumulator: {disc.diff_accumulator!r}')
    np.clip(state.diff, 0.0, 1.0, out=state.diff)


def _build_mesh(
    state: State,
    disc: Discretisation,
    top: Topology | None,
) -> Mesh:
    """Build the per-step Mesh according to `disc.topology`."""
    laplacian_kind = disc.laplacian
    if disc.topology == 'delaunay_each_step':
        return Mesh.from_positions(state.positions, kind=laplacian_kind)
    if disc.topology == 'static_with_local_update':
        if top is None:
            raise ValueError(
                "topology='static_with_local_update' requires a non-None "
                'Topology argument; call simulator.build_topology() to '
                'construct one from the initial state.'
            )
        return Mesh.from_topology(top, state.positions, kind=laplacian_kind)
    raise ValueError(f'unknown topology: {disc.topology!r}')


def divide_cells(
    state: State,
    mesh: Mesh,
    division_dist: float,
    *,
    max_per_pass: int | None = None,
    max_edge_dist: float | None = None,
    disc: Discretisation = PATH_B_DEFAULT,
    top: Topology | None = None,
) -> bool:
    """One pass of cell division. Returns True if any cells were added.

    Splits qualifying edges in `mesh.neigh_idx`. Defaults to splitting
    every disjoint qualifying edge (vertex-disjoint to avoid duplicate
    daughters when several edges share an endpoint).

    Optional safeguards (active only under `topology='delaunay_each_step'`):
      - `max_per_pass`: cap on new cells emitted in this pass.
      - `max_edge_dist`: ignore edges longer than this multiple of rest.

    Daughter cell:
      - position = midpoint of the two mothers
      - concentrations = average of mothers
      - d_i: depends on `disc.knot_daughter_di`
        - `'zero_reset'` (default, paper): hard reset to 0 if either
          mother is a knot, else inherit average.
        - `'inherit_avg'` (FORTRAN): always inherit average.
      - knot = False (paper: 'the new cell is not a knot cell')
      - half-plane flags: inherited from the mother whose y has the
        same sign as the midpoint; lingual/buccal similarly on x.

    If `top` is given, the topology is mutated to insert each new
    daughter via `Topology.insert_daughter`. This keeps the static
    adjacency graph in sync with the state arrays.
    """
    n = state.num_active
    rows = np.repeat(
        np.arange(n, dtype=np.int64),
        np.diff(mesh.neigh_starts),
    )
    cols = mesh.neigh_idx
    # Each undirected edge appears twice (i,j) and (j,i); restrict to i<j.
    keep = rows < cols
    rows_u = rows[keep]
    cols_u = cols[keep]
    if rows_u.size == 0:
        return False
    d = state.positions[cols_u] - state.positions[rows_u]
    dist = np.linalg.norm(d, axis=1)
    mask = dist >= division_dist
    if max_edge_dist is not None:
        mask &= dist <= max_edge_dist
    if not np.any(mask):
        return False

    # Insert one cell per qualifying edge, in shortest-edge-first order
    # so the most physical splits happen first. To avoid creating
    # duplicate daughters from adjacent long edges sharing a vertex, we
    # deduplicate by vertex use: each cell can participate in at most
    # one split per pass.
    long_idx = np.flatnonzero(mask)
    long_idx = long_idx[np.argsort(dist[long_idx])]
    used = np.zeros(n, dtype=bool)
    new_rows: list[int] = []
    new_cols: list[int] = []
    for e in long_idx:
        a, b = int(rows_u[e]), int(cols_u[e])
        if used[a] or used[b]:
            continue
        used[a] = used[b] = True
        new_rows.append(a)
        new_cols.append(b)
        if max_per_pass is not None and len(new_rows) >= max_per_pass:
            break
    if not new_rows:
        return False

    a_idx = np.asarray(new_rows, dtype=np.int64)
    b_idx = np.asarray(new_cols, dtype=np.int64)

    midpoint = 0.5 * (state.positions[a_idx] + state.positions[b_idx])
    avg_act = 0.5 * (state.act[a_idx] + state.act[b_idx])
    avg_inh = 0.5 * (state.inh[a_idx] + state.inh[b_idx])
    avg_sec = 0.5 * (state.sec[a_idx] + state.sec[b_idx])
    avg_diff = 0.5 * (state.diff[a_idx] + state.diff[b_idx])
    knot_parent = state.knot[a_idx] | state.knot[b_idx]

    # disc.knot_daughter_di: 'zero_reset' (paper: 'the new cell is not
    # a knot cell' — Path B reads this as resetting d_i too) vs
    # 'inherit_avg' (FORTRAN: half the mothers' average; with two-knot
    # mothers a daughter's d_i can stay near 1).
    if disc.knot_daughter_di == 'zero_reset':
        daughter_diff = np.where(knot_parent, 0.0, avg_diff)
    elif disc.knot_daughter_di == 'inherit_avg':
        daughter_diff = avg_diff
    else:
        raise ValueError(
            f'unknown knot_daughter_di: {disc.knot_daughter_di!r}'
        )
    daughter_knot = np.zeros(a_idx.shape, dtype=bool)
    avg_mes_inh = 0.5 * (state.mes_inh[a_idx] + state.mes_inh[b_idx])
    avg_mes_sec = 0.5 * (state.mes_sec[a_idx] + state.mes_sec[b_idx])

    # Half-plane flags inherit from whichever mother sits on the same
    # side of the relevant axis as the daughter.
    def inherit_along_axis(flag: np.ndarray, axis: int) -> np.ndarray:
        sign_a = np.sign(state.positions[a_idx, axis])
        sign_b = np.sign(state.positions[b_idx, axis])
        sign_d = np.sign(midpoint[:, axis])
        choose_b = (sign_b == sign_d) & (sign_a != sign_d)
        return np.where(choose_b, flag[b_idx], flag[a_idx])

    new_anterior = inherit_along_axis(state.init_anterior, 1)
    new_posterior = inherit_along_axis(state.init_posterior, 1)
    new_lingual = inherit_along_axis(state.init_lingual, 0)
    new_buccal = inherit_along_axis(state.init_buccal, 0)

    state.positions = np.concatenate([state.positions, midpoint], axis=0)
    state.act = np.concatenate([state.act, avg_act])
    state.inh = np.concatenate([state.inh, avg_inh])
    state.sec = np.concatenate([state.sec, avg_sec])
    state.diff = np.concatenate([state.diff, daughter_diff])
    state.knot = np.concatenate([state.knot, daughter_knot])
    state.mes_inh = np.concatenate([state.mes_inh, avg_mes_inh], axis=0)
    state.mes_sec = np.concatenate([state.mes_sec, avg_mes_sec], axis=0)
    state.init_anterior = np.concatenate([state.init_anterior, new_anterior])
    state.init_posterior = np.concatenate([state.init_posterior, new_posterior])
    state.init_lingual = np.concatenate([state.init_lingual, new_lingual])
    state.init_buccal = np.concatenate([state.init_buccal, new_buccal])

    # Keep the static topology in sync with the appended state arrays.
    # The Topology's daughter index is `num_cells` at insert time, which
    # matches the position we appended at — so a, b above and the
    # topology's view of the new cell line up index-for-index.
    if top is not None:
        for a, b in zip(new_rows, new_cols):
            top.insert_daughter(int(a), int(b))

    return True


def _step_jacobi(
    state: State,
    params: Params,
    mesh: Mesh,
    dt: float,
    disc: Discretisation,
) -> None:
    """Pure-Jacobi step body: forces all evaluated against start-of-step
    state, single position update at the end. Matches the v1 step.
    """
    step_reaction_diffusion(state, params, mesh, dt, disc)
    step_differentiation(state, params, dt, disc)
    forces = compute_forces(state, params, mesh, disc)
    apply_border_multipliers(forces, state, params, mesh, disc)
    state.positions = state.positions + dt * forces


def _step_gauss_seidel(
    state: State,
    params: Params,
    mesh: Mesh,
    dt: float,
    disc: Discretisation,
) -> None:
    """Gauss-Seidel-flavoured step: split mechanical forces into two
    sub-steps so nuclear traction's position update is applied before
    the remaining force components are evaluated.

    The FORTRAN's full ordering at coreop2d.py:1632-1644 is finer-
    grained (each force component reads/writes positions in turn). The
    load-bearing distinction relative to Jacobi is that the nuclear-
    traction sub-step happens against the post-rd state, and the
    remaining mechanical forces happen against the post-traction state.
    Reproducing the full per-component ordering would require slicing
    `compute_forces` into per-equation primitives — A6/A7 territory.

    The implementation here:

      1. apply_diffusion + reaction (`step_reaction_diffusion`) and
         differentiation (`step_differentiation`) on start-of-step
         state — same as Jacobi for these.
      2. Compute the full force field on start-of-step positions.
      3. Apply nuclear traction's position-update slice now (eq. 4):
         positions += dt * k_ntr * (1 - d_i) * (mean_neigh - p_i).
         This is the FORTRAN's `apply_nuclear_traction` step.
      4. Re-compute forces against the post-traction positions; subtract
         the eq. 4 component already applied; apply border-bias and
         remaining forces.

    This is a pragmatic 'apply nuclear traction first' reduction;
    deeper Gauss-Seidel ordering is a refinement.
    """
    # 1. Reactions + differentiation (Jacobi for these — FORTRAN uses a
    # copy `hq3d` for reactions, which is equivalent).
    step_reaction_diffusion(state, params, mesh, dt, disc)
    step_differentiation(state, params, dt, disc)

    # 2. Full force field on the post-rd, pre-traction positions.
    pre_forces = compute_forces(state, params, mesh, disc)

    # 3. Apply nuclear traction (eq. 4) to positions immediately, so the
    # remaining force components see the updated geometry. eq. 4 is:
    #   f_ntr = k_ntr * (1 - d_i) * (mean_neigh - p_i)
    # Re-derive it here so we have the per-cell vector to step with.
    n = state.num_active
    counts = np.diff(mesh.neigh_starts).astype(np.float64)
    rows = np.repeat(
        np.arange(n, dtype=np.int64),
        np.diff(mesh.neigh_starts),
    )
    cols = mesh.neigh_idx
    sum_neigh = np.zeros((n, 3), dtype=np.float64)
    np.add.at(sum_neigh, rows, state.positions[cols])
    safe_counts = np.maximum(counts, 1.0)[:, None]
    mean_neigh = sum_neigh / safe_counts
    one_minus_d = (1.0 - state.diff)[:, None]
    f_ntr = params.k_ntr * one_minus_d * (mean_neigh - state.positions)
    state.positions = state.positions + dt * f_ntr

    # 4. Recompute the remaining forces against the updated positions
    # (this captures the Gauss-Seidel essence). Subtract the traction
    # component that has already been applied so we don't double-count.
    new_forces = compute_forces(state, params, mesh, disc)
    # The eq. 4 sub-step in `compute_forces` is recomputed against the
    # new geometry; subtract it because we already moved by f_ntr above
    # at the *old* geometry, and the post-traction f_ntr should be
    # absorbed into the next iteration rather than applied again here.
    counts = np.diff(mesh.neigh_starts).astype(np.float64)
    sum_neigh2 = np.zeros((n, 3), dtype=np.float64)
    np.add.at(sum_neigh2, rows, state.positions[cols])
    mean_neigh2 = sum_neigh2 / np.maximum(counts, 1.0)[:, None]
    f_ntr2 = params.k_ntr * one_minus_d * (mean_neigh2 - state.positions)
    remaining = new_forces - f_ntr2

    apply_border_multipliers(remaining, state, params, mesh, disc)
    state.positions = state.positions + dt * remaining


def step(
    state: State,
    params: Params,
    dt: float = DEFAULT_DT,
    disc: Discretisation = PATH_B_DEFAULT,
    top: Topology | None = None,
) -> Mesh:
    """Run one forward-Euler iteration. Returns the mesh used this step.

    Side effects: mutates `state` (positions, concentrations, diff, knot,
    iter_count). The caller can re-use the returned mesh for inspection
    or rendering between iterations.

    `disc` selects the Path B v2 implementer-choice configuration.
    `top` is required when `disc.topology == 'static_with_local_update'`;
    ignored otherwise.
    """
    mesh = _build_mesh(state, disc, top)

    # disc.update_order: 'jacobi' (paper-PDE form, v1 default) vs
    # 'gauss_seidel_forces' / 'mixed_jacobi_gs' (FORTRAN-style: nuclear
    # traction's position update is applied before the remaining force
    # components are evaluated). See _step_gauss_seidel for the
    # simplification rationale.
    if disc.update_order == 'jacobi':
        _step_jacobi(state, params, mesh, dt, disc)
    elif disc.update_order in ('gauss_seidel_forces', 'mixed_jacobi_gs'):
        _step_gauss_seidel(state, params, mesh, dt, disc)
    else:
        raise ValueError(f'unknown update_order: {disc.update_order!r}')

    # Cell division. Under delaunay_each_step we apply the v1 safeguards
    # (cap per-pass, reject very-long edges as topology artefacts);
    # under static_with_local_update neither safeguard is needed —
    # the topology cannot grow long Delaunay-hull edges because the
    # graph never re-triangulates. The disc fields override defaults.
    if disc.topology == 'delaunay_each_step':
        max_per_pass = (
            disc.division_per_step_cap
            if disc.division_per_step_cap is not None
            else 4
        )
        max_edge = (
            disc.division_max_edge
            if disc.division_max_edge is not None
            else DIVISION_FACTOR * state.rest_length * 1.5
        )
        divide_cells(
            state,
            mesh,
            division_dist=DIVISION_FACTOR * state.rest_length,
            max_per_pass=max_per_pass,
            max_edge_dist=max_edge,
            disc=disc,
            top=None,
        )
    else:
        divide_cells(
            state,
            mesh,
            division_dist=DIVISION_FACTOR * state.rest_length,
            max_per_pass=disc.division_per_step_cap,
            max_edge_dist=disc.division_max_edge,
            disc=disc,
            top=top,
        )
    state.iter_count += 1
    return mesh


def run(
    state: State,
    params: Params,
    n_iters: int,
    dt: float = DEFAULT_DT,
    disc: Discretisation = PATH_B_DEFAULT,
    top: Topology | None = None,
) -> Mesh:
    """Run `n_iters` iterations. Returns the final mesh.

    If `disc.topology == 'static_with_local_update'` and `top` is None,
    constructs a Topology from the initial state's positions before the
    loop. The same `top` is then mutated across all iterations; cell
    division updates it via `Topology.insert_daughter`.
    """
    if disc.topology == 'static_with_local_update' and top is None:
        top = build_topology(state)
    mesh: Mesh | None = None
    for _ in range(n_iters):
        mesh = step(state, params, dt, disc, top)
    if mesh is None:
        mesh = _build_mesh(state, disc, top)
    return mesh
