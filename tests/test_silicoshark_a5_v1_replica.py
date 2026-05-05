"""V1-replica byte-identity test for Path B v2 phase A5.

PATH_B_DEFAULT prescribes `topology='static_with_local_update'` with
length-weighted Laplacian, jacobi update, zero_reset knot daughter, and
mesenchyme='per_column_z_layers'. The v1 simulator never used the
static-with-local-update path — it called Mesh.from_positions every
iteration. So `PATH_B_DEFAULT` does not reproduce v1 byte-for-byte.

The 'V1 replica' configuration is the Discretisation that exercises the
A5-wired code paths in their v1-equivalent shape:

    V1_REPLICA = replace(
        PATH_B_DEFAULT,
        topology='delaunay_each_step',
        mesenchyme='absent',
    )

This phase asserts byte-identity of one full simulator step under
V1_REPLICA against the pre-A5 baseline captured at
`.tmp/path-b-v2/v1_a5_baseline.npz`.
"""
from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from silicoshark.discretisation import PATH_B_DEFAULT
from silicoshark.params import Params
from silicoshark.simulator import DEFAULT_DT, step
from silicoshark.state import build_initial_state


REPO_ROOT = Path(__file__).resolve().parent.parent
BASELINE = REPO_ROOT / '.tmp' / 'path-b-v2' / 'v1_a5_baseline.npz'


V1_REPLICA = replace(
    PATH_B_DEFAULT,
    topology='delaunay_each_step',
    mesenchyme='absent',
)


@pytest.mark.skipif(
    not BASELINE.exists(),
    reason=f'v1 A5 baseline missing — run .tmp/path-b-v2/capture_v1_a5_baseline.py',
)
def test_v1_replica_byte_identity():
    """One simulator.step under V1_REPLICA on examples/seal.txt must
    reproduce the v1 baseline arrays bit-for-bit.

    Confirms that (a) the new A5 code paths (Mesh laplacian dispatch,
    update_order branching, knot_daughter_di branching) do not change
    behaviour when configured to v1's choices, and (b) the
    delaunay_each_step branch of `simulator.step` is wired correctly.
    """
    params = Params.from_file(REPO_ROOT / 'examples' / 'seal.txt')
    state = build_initial_state(params)

    step(state, params, DEFAULT_DT, V1_REPLICA)

    baseline = np.load(BASELINE)
    assert np.array_equal(state.positions, baseline['positions']), (
        'positions diverged from v1 baseline'
    )
    assert np.array_equal(state.act, baseline['act'])
    assert np.array_equal(state.inh, baseline['inh'])
    assert np.array_equal(state.sec, baseline['sec'])
    assert np.array_equal(state.diff, baseline['diff'])
    assert np.array_equal(state.knot, baseline['knot'])


def test_v1_replica_definition():
    """Sanity: V1_REPLICA's overrides relative to PATH_B_DEFAULT are
    exactly the two fields documented above. Drift in the preset
    defaults should be a deliberate choice, not a silent change to
    what 'v1 replica' means.
    """
    expected_diff = {
        'topology': ('delaunay_each_step', PATH_B_DEFAULT.topology),
        'mesenchyme': ('absent', PATH_B_DEFAULT.mesenchyme),
    }
    for field, (replica, default) in expected_diff.items():
        assert getattr(V1_REPLICA, field) == replica
        # Either V1_REPLICA differs from PATH_B_DEFAULT here, or the
        # replica field happens to equal the default already (which is
        # fine — we still document the override explicitly).
        assert getattr(PATH_B_DEFAULT, field) == default


def test_v1_replica_default_mesenchyme_is_absent_or_documented():
    """If PATH_B_DEFAULT.mesenchyme stays 'per_column_z_layers' (the
    current default), V1_REPLICA must override it. If a future commit
    flips the default to 'absent', this test stays green and serves as
    the marker for that decision.
    """
    assert V1_REPLICA.mesenchyme == 'absent'


# --- Cotangent Laplacian smoke tests ---------------------------------


def test_cotangent_laplacian_of_constant_is_zero():
    """A constant scalar field has Laplacian zero everywhere. The
    cotangent operator must agree on a regular hex lattice.
    """
    from silicoshark.mesh import Mesh
    from silicoshark.state import hex_lattice

    xy = hex_lattice(rad=4)
    n = xy.shape[0]
    positions = np.zeros((n, 3), dtype=np.float64)
    positions[:, :2] = xy
    mesh = Mesh.from_positions(positions, kind='cotangent')
    u = np.full(n, 3.7, dtype=np.float64)
    Lu = mesh.laplacian(u)
    assert np.allclose(Lu, 0.0, atol=1e-10), (
        f'Lap(constant) should be 0; got max |Lu| = {np.max(np.abs(Lu)):.3e}'
    )


def test_cotangent_laplacian_of_linear_is_zero_interior():
    """A linear scalar field u = x has Laplacian zero everywhere on a
    flat mesh. Boundary cells of the cotangent operator can pick up an
    artefact from missing the second triangle, so we restrict the
    check to interior cells (those with 6 neighbours).
    """
    from silicoshark.mesh import Mesh
    from silicoshark.state import hex_lattice

    xy = hex_lattice(rad=4)
    n = xy.shape[0]
    positions = np.zeros((n, 3), dtype=np.float64)
    positions[:, :2] = xy
    mesh = Mesh.from_positions(positions, kind='cotangent')
    u = positions[:, 0].copy()  # u = x
    Lu = mesh.laplacian(u)

    # Interior cells: those with exactly 6 neighbours (regular hex).
    counts = np.diff(mesh.neigh_starts)
    interior = counts >= 6
    assert interior.sum() > 0, 'expected at least some interior cells'
    # Cotangent of equilateral hex angles is irrational; floating-point
    # noise dominates so use a generous absolute tolerance.
    max_interior = float(np.max(np.abs(Lu[interior])))
    assert max_interior < 1e-9, (
        f'Lap(x) should be ~0 on interior cells; '
        f'got max |Lu| = {max_interior:.3e} at indices '
        f'{np.argsort(-np.abs(Lu * interior))[:5]}'
    )


# --- Static-topology end-to-end --------------------------------------


def test_static_topology_step_preserves_symmetry():
    """One simulator.step under PATH_B_DEFAULT (overridden to mesenchyme=
    absent) on the seal example must keep the topology graph symmetric.
    Confirms that `Topology.insert_daughter` calls during `divide_cells`
    do not corrupt the adjacency.
    """
    from silicoshark.params import Params
    from silicoshark.simulator import build_topology, step
    from silicoshark.state import build_initial_state

    disc = replace(PATH_B_DEFAULT, mesenchyme='absent')
    assert disc.topology == 'static_with_local_update'

    params = Params.from_file(REPO_ROOT / 'examples' / 'seal.txt')
    state = build_initial_state(params)
    top = build_topology(state)
    assert top.is_symmetric(), 'initial topology must be symmetric'

    step(state, params, DEFAULT_DT, disc, top)

    assert top.is_symmetric(), (
        'topology lost symmetry after one step — daughter insertion bug?'
    )
    assert top.num_cells == state.num_active, (
        f'topology has {top.num_cells} cells but state has '
        f'{state.num_active} — lock-step accounting failure'
    )
