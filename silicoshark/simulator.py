"""Main iteration loop for silicoshark.

One forward-Euler step per iteration:

  1. Triangulate (Delaunay) and build the mesh adjacency.
  2. Reaction-diffusion update (Act, Inh, Sec). Pure Jacobi.
  3. Differentiation update: d_i += dt * k_dff * [Sec], clipped to [0, 1].
     (paper eq. 6 says [Act]; FORTRAN uses [Sec] — Path B follows the
      FORTRAN per the charter; see
      docs/findings/2026-05-04-differentiation-uses-sec-not-act.md.)
  4. Compute mechanical forces (Rep, Adh, Ntr, Egr, cervical-loop, Boy).
  5. Apply border-bias multipliers (Pbi, Abi, Bgr); z clamp.
  6. Update positions: positions += dt * forces.
  7. Cell division: any edge with |d| >= division_threshold (= 2 * rest)
     splits — insert a new cell at the midpoint, average concentrations,
     reset d_i if the daughter is non-knot but mothers are knots.

`dt` matches FORTRAN's delta = 0.05 by default. Iteration counts are
specified by the caller; FORTRAN goldens use 100 iters/save, 5 saves.
"""
from __future__ import annotations

import numpy as np

from .params import Params
from .state import State, N_MESENCHYME_LAYERS
from .mesh import Mesh
from .reaction import step_reaction_diffusion
from .forces import compute_forces, apply_border_multipliers
from .discretisation import Discretisation, PATH_B_DEFAULT


DEFAULT_DT: float = 0.05
DIVISION_FACTOR: float = 2.0  # edge length threshold = DIVISION_FACTOR * rest


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


def divide_cells(
    state: State,
    mesh: Mesh,
    division_dist: float,
    *,
    max_per_pass: int | None = None,
    max_edge_dist: float | None = None,
) -> bool:
    """One pass of cell division. Returns True if any cells were added.

    Splits qualifying edges in `mesh.neigh_idx`. Defaults to splitting
    every disjoint qualifying edge (vertex-disjoint to avoid duplicate
    daughters when several edges share an endpoint).

    Optional safeguards:
      - `max_per_pass`: cap on new cells emitted in this pass. Used to
        prevent cascade-style runaway when Delaunay re-triangulation
        produces many simultaneous long edges.
      - `max_edge_dist`: ignore edges longer than this. Long Delaunay
        edges (e.g. across the convex hull when border cells migrate
        outward) are topology artefacts, not biological neighbours,
        and should not trigger division.

    Daughter cell:
      - position = midpoint of the two mothers
      - concentrations = average of mothers
      - d_i = 0 if either mother is a knot (paper: "the new cell is
              not a knot cell" — Path B reads this as a hard reset)
      - knot = False
      - half-plane flags: inherited from the mother whose y has the
        same sign as the midpoint; lingual/buccal similarly on x.
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
    # Paper: "If one of the mother cells is an enamel knot cell, the
    # new cell is not a knot cell." Reset d_i for that case.
    daughter_diff = np.where(knot_parent, 0.0, avg_diff)
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
    return True


def step(
    state: State,
    params: Params,
    dt: float = DEFAULT_DT,
    disc: Discretisation = PATH_B_DEFAULT,
) -> Mesh:
    """Run one forward-Euler iteration. Returns the mesh used this step.

    Side effects: mutates `state` (positions, concentrations, diff, knot,
    iter_count). The caller can re-use the returned mesh for inspection
    or rendering between iterations.

    `disc` selects the Path B v2 implementer-choice configuration. The
    default `PATH_B_DEFAULT` reproduces the v1 behaviour byte-for-byte
    on the reaction-diffusion + differentiation path; force-side
    branching arrives in A4.
    """
    mesh = Mesh.from_positions(state.positions)
    step_reaction_diffusion(state, params, mesh, dt, disc)
    step_differentiation(state, params, dt, disc)
    forces = compute_forces(state, params, mesh, disc)
    apply_border_multipliers(forces, state, params, mesh, disc)
    state.positions = state.positions + dt * forces
    # Cap divisions per step and reject very-long Delaunay edges as
    # topology artefacts; otherwise the convex-hull-spanning edges of
    # a re-triangulated mesh trigger cascading splits.
    divide_cells(
        state,
        mesh,
        division_dist=DIVISION_FACTOR * state.rest_length,
        max_per_pass=4,
        max_edge_dist=DIVISION_FACTOR * state.rest_length * 1.5,
    )
    state.iter_count += 1
    return mesh


def run(
    state: State,
    params: Params,
    n_iters: int,
    dt: float = DEFAULT_DT,
    disc: Discretisation = PATH_B_DEFAULT,
) -> Mesh:
    """Run `n_iters` iterations. Returns the final mesh."""
    mesh: Mesh | None = None
    for _ in range(n_iters):
        mesh = step(state, params, dt, disc)
    if mesh is None:
        mesh = Mesh.from_positions(state.positions)
    return mesh
