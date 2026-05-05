"""Unit tests for Path B v2 mechanical-force branches.

Each test exercises one `Discretisation` field by constructing a
minimal `State` + `Params` and asserting that the resulting forces
array agrees with the analytical expectation under each option.

These tests are mathematical, not biological: they pin down the
arithmetic of a single force-evaluation on small (2- or 7-cell)
states. The biological-regression tests live in
`tests/test_simulator_smoke.py` and (eventually)
`tests/test_silicoshark_smoke.py`.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from silicoshark.discretisation import PATH_B_DEFAULT
from silicoshark.forces import compute_forces, apply_border_multipliers
from silicoshark.mesh import Mesh
from silicoshark.params import Params
from silicoshark.state import N_MESENCHYME_LAYERS, State, build_initial_state


# --- Fixture helpers -------------------------------------------------


def _zero_params(**overrides) -> Params:
    """Params with every coefficient zero. Caller overrides what they
    care about. Useful for isolating one force term at a time.
    """
    p = Params()
    for k, v in overrides.items():
        setattr(p, k, v)
    return p


def _state_from_positions(
    positions: np.ndarray,
    *,
    diff: np.ndarray | None = None,
    sec: np.ndarray | None = None,
    first_border_cell: int = 0,
) -> State:
    """Build a minimal State from a positions array. All concentrations
    zero unless overridden via kwargs.
    """
    n = positions.shape[0]
    return State(
        positions=positions.astype(np.float64),
        act=np.zeros(n, dtype=np.float64),
        inh=np.zeros(n, dtype=np.float64),
        sec=np.zeros(n, dtype=np.float64) if sec is None else sec.astype(np.float64),
        diff=np.zeros(n, dtype=np.float64) if diff is None else diff.astype(np.float64),
        knot=np.zeros(n, dtype=bool),
        mes_inh=np.zeros((n, N_MESENCHYME_LAYERS), dtype=np.float64),
        mes_sec=np.zeros((n, N_MESENCHYME_LAYERS), dtype=np.float64),
        init_anterior=np.zeros(n, dtype=bool),
        init_posterior=np.zeros(n, dtype=bool),
        init_lingual=np.zeros(n, dtype=bool),
        init_buccal=np.zeros(n, dtype=bool),
        first_border_cell=first_border_cell,
    )


# --- rep_form: hookean_signed vs paper_gated -------------------------


def test_rep_form_compressed_pair_matches():
    """Two cells at distance 0.8 (< rest=1.0). Both rep_form options
    must produce the SAME force on cell 0: the Hookean form is
    `k_rep_clip * (0.8 - 1.0) * d_vec` = `-0.2 * 0.8 * k_rep_clip` on
    the x-axis. Under paper_gated this branch is active because
    |d| < rest. Pick k_rep = 0.5 so rep_clip = min(0.5, 1.0) = 0.5.

    Use 4 cells in a tight cluster (all pairs < rest) so the test
    isolates the rep_form choice. The (0, 1), (0, 2), (0, 3) edges
    must all be < rest for hookean_signed to match paper_gated; we
    place cells so this is guaranteed.
    """
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.8, 0.0, 0.0],   # 0.8 < rest
            [0.4, 0.7, 0.0],   # |d| = sqrt(0.16+0.49) ≈ 0.806 < rest
            [-0.4, 0.7, 0.0],  # |d| = sqrt(0.16+0.49) ≈ 0.806 < rest
        ],
        dtype=np.float64,
    )
    state = _state_from_positions(positions)
    params = _zero_params(k_rep=0.5)  # rep_clip = 0.5
    mesh = Mesh.from_positions(state.positions)

    f_signed = compute_forces(
        state, params, mesh, replace(PATH_B_DEFAULT, rep_form='hookean_signed'),
    )
    f_gated = compute_forces(
        state, params, mesh, replace(PATH_B_DEFAULT, rep_form='paper_gated'),
    )
    # When all edges have |d| < rest the two forms agree exactly.
    np.testing.assert_array_equal(f_signed, f_gated)


def test_rep_form_stretched_pair_diverges():
    """Two cells along the x-axis at distance 1.2 (> rest=1.0). The
    other cells are placed far enough that they're also stretched
    edges, but we compare cells (0, 1)'s pair specifically through
    the difference.

    Tighter assertion: with k_rep=0.5 and rep_clip=0.5, hookean_signed
    contributes `0.5 * (1.2 - 1.0) * (1.2, 0, 0) = (0.12, 0, 0)` from
    the (0, 1) edge. paper_gated contributes 0 from that edge.

    The OTHER edges (0, 2) and (0, 3) ALSO differ between the two
    branches if their distances are also > rest. To isolate the (0, 1)
    edge we place cells 2, 3 within rest:
        cell 2 at (0.4, 0.7, 0) — distance from cell 0 = 0.806 < rest
        cell 3 at (0.8, 0.7, 0) — distance from cell 0 = 1.061 > rest

    Hmm, this is tricky. Simplest: place three cells where only (0, 1)
    is stretched and (0, 2), (0, 3) are below rest.
    """
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.2, 0.0, 0.0],   # 1.2 > rest
            [0.4, 0.7, 0.0],   # |d| ≈ 0.806 < rest
            [-0.4, 0.7, 0.0],  # |d| ≈ 0.806 < rest
        ],
        dtype=np.float64,
    )
    state = _state_from_positions(positions)
    params = _zero_params(k_rep=0.5)
    mesh = Mesh.from_positions(state.positions)

    # First check whether (0, 1) is in fact a Delaunay edge. If not,
    # the test is moot; force a check.
    rows = np.repeat(
        np.arange(state.num_active, dtype=np.int64),
        np.diff(mesh.neigh_starts),
    )
    cols = mesh.neigh_idx
    pairs = set()
    for k in range(rows.shape[0]):
        a, b = int(rows[k]), int(cols[k])
        pairs.add((min(a, b), max(a, b)))
    assert (0, 1) in pairs, 'expected (0, 1) to be a mesh edge'

    f_signed = compute_forces(
        state, params, mesh, replace(PATH_B_DEFAULT, rep_form='hookean_signed'),
    )
    f_gated = compute_forces(
        state, params, mesh, replace(PATH_B_DEFAULT, rep_form='paper_gated'),
    )
    # Compute the expected difference contributed by stretched edges
    # at cell 0. Find every neighbour of cell 0 with dist > rest and
    # sum the hookean_signed contribution they would give (which is
    # exactly what paper_gated zeroes out).
    rest = state.rest_length
    rep_clip = min(params.k_rep, 1.0)
    expected_delta = np.zeros(3)
    for j in mesh.neigh_idx[mesh.neigh_starts[0]:mesh.neigh_starts[0 + 1]]:
        d_vec = state.positions[j] - state.positions[0]
        d_mag = float(np.linalg.norm(d_vec))
        if d_mag > rest:
            expected_delta += rep_clip * (d_mag - rest) * d_vec
    diff = f_signed - f_gated
    np.testing.assert_allclose(diff[0], expected_delta, atol=1e-12)
    # And the difference is non-zero: cell 0 has at least one stretched edge.
    assert np.linalg.norm(diff[0]) > 1e-9


# --- adh_form: unit_vector vs hookean_attraction ---------------------


def test_adh_form_branches_differ_by_factor_d():
    """Two cells at distance 1.5 (> rest=1.0). With k_adh = 1.0 and
    k_rep = 0.0:
      unit_vector:        f_0 = 1.0 * (1.5, 0, 0) / 1.5 = (1.0, 0, 0)
      hookean_attraction: f_0 = 1.0 * (1.5, 0, 0)       = (1.5, 0, 0)

    Magnitudes differ by a factor of 1.5.
    """
    positions = np.array(
        [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.75, 5.0, 0.0], [-3.0, 5.0, 0.0]],
        dtype=np.float64,
    )
    state = _state_from_positions(positions)
    params = _zero_params(k_adh=1.0)
    mesh = Mesh.from_positions(state.positions)

    f_unit = compute_forces(
        state, params, mesh, replace(PATH_B_DEFAULT, adh_form='unit_vector'),
    )
    f_hook = compute_forces(
        state, params, mesh, replace(PATH_B_DEFAULT, adh_form='hookean_attraction'),
    )
    # Cell 0's adhesion contribution from (0, 1) edge:
    # unit_vector pulls by (1.0, 0, 0); hookean by (1.5, 0, 0).
    # The (0, 2), (0, 3) edges contribute equally under both forms iff
    # they are at the same distance — they aren't (cell 2 is at sqrt(0.75^2+25),
    # cell 3 is at sqrt(9+25)). So we just check the difference.
    delta_0 = f_hook[0] - f_unit[0]
    # The (0, 1) edge contributes (1.5 - 1.0) = 0.5 to delta on x.
    # The (0, 2), (0, 3) edges contribute differences too.
    # Cell 2 is at (0.75, 5, 0): d = (0.75, 5, 0), |d| = sqrt(25.5625) ≈ 5.0561
    # unit:  (0.75, 5, 0) / 5.0561 = (0.14834, 0.98893, 0)
    # hook:  (0.75, 5, 0)
    # delta: (0.75 - 0.14834, 5 - 0.98893, 0) = (0.60166, 4.01107, 0)
    # Cell 3 is at (-3, 5, 0): d = (-3, 5, 0), |d| = sqrt(34) ≈ 5.8310
    # unit:  (-3, 5, 0) / 5.8310 = (-0.51450, 0.85749, 0)
    # hook:  (-3, 5, 0)
    # delta: (-3 - -0.51450, 5 - 0.85749, 0) = (-2.48550, 4.14251, 0)
    # Total expected delta on cell 0:
    # x: 0.5 + 0.60166 + -2.48550 = -1.38384
    # y: 0 + 4.01107 + 4.14251 = 8.15358
    expected_x = 0.5 + (0.75 - 0.75 / np.sqrt(25.5625)) + (-3.0 - -3.0 / np.sqrt(34.0))
    expected_y = 0.0 + (5.0 - 5.0 / np.sqrt(25.5625)) + (5.0 - 5.0 / np.sqrt(34.0))
    assert delta_0[0] == pytest.approx(expected_x, abs=1e-10)
    assert delta_0[1] == pytest.approx(expected_y, abs=1e-10)
    # And the magnitudes are not equal.
    assert not np.allclose(f_unit[0], f_hook[0])


# --- eq5_z_gate: True vs False ---------------------------------------


def test_eq5_z_gate_branches_differ():
    """3-cell asymmetric setup: cell 0 at origin, neighbour 1 above
    (z=1.5), neighbour 2 below (z=0.5), neighbour 3 at same z (level).

    Under True (gated): only neighbour 1 contributes (z_diff > 1e-4).
        sum_unit = unit(p1 - p0)
    Under False (open): all three neighbours contribute.

    Other force terms (cervical-loop, repulsion, adhesion, traction) are
    invariant under the eq5_z_gate choice — both branches see the same
    contributions from those terms. So the DIFFERENCE
    `f_open[0] - f_gated[0]` equals the eq.5 contribution diff alone.
    """
    # Cell positions chosen so the (x, y) projection is non-degenerate
    # (Delaunay needs distinct (x, y) per cell). z varies per cell so
    # the gate fires asymmetrically.
    positions = np.array(
        [
            [0.0, 0.0, 1.0],    # cell 0 (focal)
            [1.0, 0.1, 1.5],    # cell 1 (above in z)
            [1.0, -0.1, 0.5],   # cell 2 (below in z)
            [-1.0, 0.0, 1.0],   # cell 3 (same z as cell 0)
        ],
        dtype=np.float64,
    )
    state = _state_from_positions(positions)
    params = _zero_params(k_egr=1.0)
    mesh = Mesh.from_positions(state.positions)

    # All four cells must be mesh neighbours of cell 0 for the test.
    cell0_neighs = mesh.neigh_idx[mesh.neigh_starts[0]:mesh.neigh_starts[0 + 1]]
    assert set(cell0_neighs.tolist()) == {1, 2, 3}, (
        f'cell 0 should neighbour {{1,2,3}}; got {set(cell0_neighs.tolist())}'
    )

    f_gated = compute_forces(
        state, params, mesh, replace(PATH_B_DEFAULT, eq5_z_gate=True),
    )
    f_open = compute_forces(
        state, params, mesh, replace(PATH_B_DEFAULT, eq5_z_gate=False),
    )

    # Compute expected eq.5 unit-vector sums on cell 0.
    def _unit(p_to, p_from):
        v = p_to - p_from
        return v / np.linalg.norm(v)

    u01 = _unit(positions[1], positions[0])  # z above
    u02 = _unit(positions[2], positions[0])  # z below
    u03 = _unit(positions[3], positions[0])  # z same

    sum_gated = u01.copy()             # only z-above neighbour
    egr_gated = sum_gated / np.linalg.norm(sum_gated)

    sum_open = u01 + u02 + u03         # all neighbours
    egr_open = sum_open / np.linalg.norm(sum_open)

    expected_eq5_diff = egr_open - egr_gated  # k_egr * (1 - d) factor = 1.0
    actual_diff = f_open[0] - f_gated[0]
    np.testing.assert_allclose(actual_diff, expected_eq5_diff, atol=1e-12)
    # And the difference is non-trivial.
    assert np.linalg.norm(actual_diff) > 1e-6


# --- eq5_apply_to: 'all' vs 'interior_only' --------------------------


def test_eq5_apply_to_branches_differ():
    """Build a seal initial state (rad=4 → 37 cells, first_border_cell=18).

    Zero out every coefficient except k_egr, then perturb the geometry so
    eq.5 has something to do (lift cell 0 up so its neighbours are below
    in z, then check force; we use k_egr=1.0 and check on a cell with
    index < first_border_cell).

    Under 'all': interior cells (index < 18) DO receive the eq.5
    contribution.
    Under 'interior_only': cells with index < 18 do NOT receive eq.5.

    For this test to fire the gate has to be open (z-asymmetric), so we
    perturb the z of one boundary so its neighbours are all above it.
    """
    params = Params.from_file('examples/seal.txt')
    state = build_initial_state(params)
    # Push cell 0 (centre) down so its 6 neighbours are above in z.
    state.positions[0, 2] = 0.5

    # Zero out every coefficient except k_egr; rest_length stays 1.0.
    test_params = _zero_params(k_egr=1.0)
    mesh = Mesh.from_positions(state.positions)

    f_all = compute_forces(
        state, test_params, mesh, replace(PATH_B_DEFAULT, eq5_apply_to='all'),
    )
    f_interior = compute_forces(
        state, test_params, mesh, replace(PATH_B_DEFAULT, eq5_apply_to='interior_only'),
    )

    # Cell 0 has index 0 < first_border_cell = 6 * (4 - 1) = 18. Under
    # 'all' it receives eq.5 (non-zero); under 'interior_only' it doesn't.
    fbc = state.first_border_cell
    assert fbc == 18
    assert np.linalg.norm(f_all[0]) > 1e-6, 'cell 0 eq.5 should be non-zero under all'
    assert np.allclose(f_interior[0], 0.0, atol=1e-12), (
        'cell 0 eq.5 should be zero under interior_only (index 0 < first_border_cell=18)'
    )

    # An interior-by-FORTRAN cell (index 18) should be equal under both
    # branches, because both grant it eq.5.
    if np.linalg.norm(f_all[fbc]) > 1e-6:
        np.testing.assert_allclose(f_all[fbc], f_interior[fbc], atol=1e-12)


# --- rep_neighbour_set: 'mesh' vs 'mesh_plus_all_close' --------------


def test_rep_neighbour_set_branches_differ():
    """4-cell square layout where one cell pair is mesh-adjacent and
    another non-adjacent pair is within 1.4 of each other.

    Place 4 cells in a 1.0 x 1.3 rectangle. The diagonal (length
    sqrt(1 + 1.69) ≈ 1.642) is too long, but Delaunay's diagonal
    choice will make one diagonal a neighbour and the other not.
    To get a clean non-mesh pair within 1.4, use a layout where the
    'across' diagonal is shorter.

    Layout:
        cell 0 at (0, 0)
        cell 1 at (1, 0)
        cell 2 at (0.5, 0.86)   (equilateral with 0,1)
        cell 3 at (0.5, -0.86)

    Distances: (0,1) = 1, (0,2) = 1, (0,3) = 1, (1,2) = 1, (1,3) = 1.
    The (2, 3) distance = 1.72 — outside 1.4. Bad layout.

    Better: a 4-cell rhombus where the 'short' diagonal is < 1.4 but
    not part of the Delaunay graph. Use:
        cell 0 at (0, 0)
        cell 1 at (1.2, 0)
        cell 2 at (0.6, 1.0)
        cell 3 at (0.6, -1.0)

    Pair (0, 1) distance 1.2 — adjacent (along x-axis).
    Pair (2, 3) distance 2.0 — too far.
    Pair (0, 2) = sqrt(0.36 + 1) = sqrt(1.36) ≈ 1.166 — adjacent.

    Try a 5th-cell setup. Actually simpler: place 4 cells in a square
    with side 1.0 plus one extra non-adjacent cell.

    Layout:
        cell 0 at (0, 0)
        cell 1 at (1, 0)
        cell 2 at (0, 1)
        cell 3 at (1, 1)
        cell 4 at (2, 0.5)

    Delaunay over (x, y): triangles will be (0,1,2), (1,2,3), (1,3,4) etc.
    Cell 0 and cell 4 are at distance sqrt(4 + 0.25) ≈ 2.06 — too far.

    Simpler: just check that adding the mesh_plus_all_close BRANCH
    introduces at least one extra repulsive contribution by setting
    up a non-mesh close pair explicitly.

    Easiest: 6-cell arrangement where Delaunay will draw an edge that
    isn't between two specific close cells. Use:
        cell 0 at (0.0, 0.0)
        cell 1 at (1.0, 0.0)
        cell 2 at (0.5, 0.866)
        cell 3 at (0.5, -0.866)

    The diagonal (2, 3) is at distance 1.732 — outside 1.4.

    Instead use cell 2 at (0.5, 0.6) and cell 3 at (0.5, -0.6) so
    the (2, 3) distance is 1.2 < 1.4 but Delaunay (in 2D) will make
    triangles (0,1,2), (0,1,3) — putting (2, 3) NOT as a mesh edge.

    Distances:
      (0, 1) = 1.0 (adjacent, mesh)
      (0, 2) = sqrt(0.25 + 0.36) = sqrt(0.61) ≈ 0.781 (adjacent)
      (0, 3) = sqrt(0.25 + 0.36) ≈ 0.781 (adjacent)
      (1, 2) = sqrt(0.25 + 0.36) ≈ 0.781 (adjacent)
      (1, 3) = sqrt(0.25 + 0.36) ≈ 0.781 (adjacent)
      (2, 3) = 1.2 (not adjacent in Delaunay because cells 0, 1
              partition them)

    Under 'mesh', cells 2 and 3 only receive forces from their mesh
    neighbours. Under 'mesh_plus_all_close', they additionally get a
    polynomial repulsion from each other (distance 1.2 < 1.4).
    """
    positions = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.5, 0.6, 0.0],
            [0.5, -0.6, 0.0],
        ],
        dtype=np.float64,
    )
    state = _state_from_positions(positions)
    # Set k_rep > 0 (and rest of params zero for clarity).
    params = _zero_params(k_rep=1.0)
    mesh = Mesh.from_positions(state.positions)

    # Verify (2, 3) is NOT a mesh edge.
    rows = np.repeat(
        np.arange(state.num_active, dtype=np.int64),
        np.diff(mesh.neigh_starts),
    )
    cols = mesh.neigh_idx
    pairs = set()
    for i in range(rows.shape[0]):
        a = int(rows[i])
        b = int(cols[i])
        pairs.add((min(a, b), max(a, b)))
    assert (2, 3) not in pairs, 'expected (2, 3) to be non-adjacent in Delaunay'

    f_mesh = compute_forces(
        state, params, mesh, replace(PATH_B_DEFAULT, rep_neighbour_set='mesh'),
    )
    f_plus = compute_forces(
        state, params, mesh, replace(PATH_B_DEFAULT, rep_neighbour_set='mesh_plus_all_close'),
    )

    delta_2 = f_plus[2] - f_mesh[2]
    delta_3 = f_plus[3] - f_mesh[3]

    # Both cells should have non-zero deltas (the (2, 3) pair contributes
    # a polynomial repulsion at distance 1.2).
    assert np.linalg.norm(delta_2) > 1e-12, 'cell 2 should gain non-mesh repulsion'
    assert np.linalg.norm(delta_3) > 1e-12, 'cell 3 should gain non-mesh repulsion'

    # And the deltas should be opposite in y (cell 2 pushed up by cell 3
    # below it).
    assert delta_2[1] > 0
    assert delta_3[1] < 0


# --- border_bias_x_zero_quirk: True vs False -------------------------


def test_border_bias_x_zero_quirk_branches_differ():
    """Build a state where a border cell sits at exactly x = 0, |y| < Bwi.
    Under True (default): no Pbi/Bgr applied. Under False: Pbi and Bgr
    applied.

    Use a 7-cell rad=2 hex (cells 1..6 are at unit distance from origin
    on a hex grid; index 4 is at (0, sqrt(3)/2 ≈ 0.866) which is NOT
    exactly x=0 in the rad=2 hex layout — let me check the hex_lattice
    output).

    Actually rad=2 produces 7 cells but none are at x=0 except the
    centre (cell 0, which is interior, not border). For this test we
    construct positions explicitly.
    """
    # 4 cells: a hex-ish arrangement. Cell 0 at (0, 0.5) — exactly x=0,
    # on the y-axis. Cells 1, 2, 3 around it.
    positions = np.array(
        [
            [0.0, 0.5, 1.0],     # focal cell — exactly x=0, |y|=0.5 < Bwi
            [1.0, 0.5, 1.0],     # east neighbour
            [-1.0, 0.5, 1.0],    # west neighbour
            [0.0, 1.5, 1.0],     # north neighbour
        ],
        dtype=np.float64,
    )
    state = _state_from_positions(positions)
    mesh = Mesh.from_positions(state.positions)

    # Cell 0 must be a border cell (Delaunay over 4 cells will give it
    # at most 3 neighbours — yes, cell 0 sits in the middle).
    # Actually with 4 cells in this layout cell 0 may have 3 neighbours
    # which makes it a border cell (< 6).
    assert mesh.is_border[0], 'cell 0 must be a border cell (< 6 neighbours)'

    # Hand-build a forces array with non-trivial values on cell 0.
    forces_quirk = np.array(
        [
            [1.0, 0.0, 1.0],   # cell 0: x-force=1, z-force=1
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        ],
        dtype=np.float64,
    )
    forces_no_quirk = forces_quirk.copy()

    params = _zero_params()
    params.k_pbi = 0.5  # << 1, so applying it will halve the x-force
    params.k_abi = 0.5
    params.k_bgr = 5.0  # >> 1, so applying it will multiply the z-force
    params.k_bwi = 1.0  # |y|=0.5 < 1.0 = Bwi → cell 0 is in band

    apply_border_multipliers(
        forces_quirk, state, params, mesh,
        replace(PATH_B_DEFAULT, border_bias_x_zero_quirk=True),
    )
    apply_border_multipliers(
        forces_no_quirk, state, params, mesh,
        replace(PATH_B_DEFAULT, border_bias_x_zero_quirk=False),
    )

    # Under quirk=True, cell 0 (x exactly 0) is skipped. So x-force
    # stays 1.0; z-force is 1.0 (then z-clamp leaves it alone since > 0).
    assert forces_quirk[0, 0] == pytest.approx(1.0, abs=1e-12)
    assert forces_quirk[0, 2] == pytest.approx(1.0, abs=1e-12)

    # Under quirk=False, cell 0 is treated as east (x >= 0): Pbi*x,
    # Bgr*z. So x-force = 0.5 * 1.0 = 0.5; z-force = 5.0 * 1.0 = 5.0.
    assert forces_no_quirk[0, 0] == pytest.approx(0.5, abs=1e-12)
    assert forces_no_quirk[0, 2] == pytest.approx(5.0, abs=1e-12)


# --- Default-preset equivalence to v1 (sanity check) -----------------


def test_default_disc_eq_omitted_disc_compute_forces():
    """`PATH_B_DEFAULT` must be the implicit default for compute_forces."""
    params = Params.from_file('examples/seal.txt')
    state = build_initial_state(params)
    mesh = Mesh.from_positions(state.positions)

    f_implicit = compute_forces(state, params, mesh)
    f_explicit = compute_forces(state, params, mesh, PATH_B_DEFAULT)
    np.testing.assert_array_equal(f_implicit, f_explicit)


def test_default_disc_eq_omitted_disc_apply_border():
    """`PATH_B_DEFAULT` must be the implicit default for apply_border_multipliers."""
    params = Params.from_file('examples/seal.txt')
    state = build_initial_state(params)
    mesh = Mesh.from_positions(state.positions)

    f_implicit = compute_forces(state, params, mesh, PATH_B_DEFAULT)
    f_explicit = f_implicit.copy()
    apply_border_multipliers(f_implicit, state, params, mesh)
    apply_border_multipliers(f_explicit, state, params, mesh, PATH_B_DEFAULT)
    np.testing.assert_array_equal(f_implicit, f_explicit)
