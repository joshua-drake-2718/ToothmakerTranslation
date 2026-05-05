"""Tests for `silicoshark.metrics`.

Each function gets at least one targeted test on hand-crafted input.
The metrics are simple by design (no PDE state, no integration), so
the tests are arithmetic checks rather than smoke checks against the
simulator.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from silicoshark.metrics import (
    cell_count_plateau,
    cusp_count,
    cusp_width,
    regime,
    vertex_envelope,
)
from silicoshark.mesh import Mesh
from silicoshark.state import State


def _make_state(positions: np.ndarray, knot: np.ndarray) -> State:
    """Construct a minimal State with the fields the metrics consume."""
    n = positions.shape[0]
    zeros_n = np.zeros(n, dtype=np.float64)
    zeros_b = np.zeros(n, dtype=bool)
    zeros_mes = np.zeros((n, 2), dtype=np.float64)
    return State(
        positions=positions.astype(np.float64),
        act=zeros_n.copy(),
        inh=zeros_n.copy(),
        sec=zeros_n.copy(),
        diff=zeros_n.copy(),
        knot=np.asarray(knot, dtype=bool),
        mes_inh=zeros_mes.copy(),
        mes_sec=zeros_mes.copy(),
        init_anterior=zeros_b.copy(),
        init_posterior=zeros_b.copy(),
        init_lingual=zeros_b.copy(),
        init_buccal=zeros_b.copy(),
    )


def test_cusp_count_counts_knot_cells():
    pos = np.zeros((5, 3))
    pos[:, 0] = np.arange(5)
    knot = np.array([False, True, True, False, True])
    state = _make_state(pos, knot)
    assert cusp_count(state) == 3


def test_cusp_count_zero_when_no_knots():
    pos = np.zeros((4, 3))
    state = _make_state(pos, np.zeros(4, dtype=bool))
    assert cusp_count(state) == 0


def test_cusp_width_zero_when_fewer_than_two_knots():
    # Build a mesh from a small triangle layout; knot count = 1.
    pos = np.array([
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.5, 1.0, 0.0],
        [1.5, 1.0, 0.0],
    ], dtype=np.float64)
    knot = np.array([False, True, False, False])
    state = _make_state(pos, knot)
    mesh = Mesh.from_positions(pos)
    assert cusp_width(state, mesh) == 0.0


def test_cusp_width_mean_pairwise_distance_between_knot_neighbours():
    # Triangle 0-1-2 plus a non-mesh-adjacent cell 3 far away.
    pos = np.array([
        [0.0, 0.0, 0.0],
        [3.0, 0.0, 0.0],
        [1.5, 2.0, 0.0],
        [10.0, 10.0, 0.0],
    ], dtype=np.float64)
    knot = np.array([True, True, True, False])
    state = _make_state(pos, knot)
    mesh = Mesh.from_positions(pos)
    # Mesh edges among knot cells are 0-1 (3.0), 1-2 (sqrt(1.5^2+2^2)=2.5),
    # 0-2 (sqrt(1.5^2+2^2)=2.5). Mean = (3.0 + 2.5 + 2.5) / 3 = 2.6666...
    expected = (3.0 + 2.5 + 2.5) / 3.0
    assert cusp_width(state, mesh) == pytest.approx(expected)


def test_cell_count_plateau_short_history_falls_back_to_final():
    assert cell_count_plateau([1, 2, 3, 4]) == 4.0


def test_cell_count_plateau_long_history_uses_last_tenth():
    history = [10] * 90 + [50] * 10
    # 100 entries, last 10% (10 entries) are all 50.
    assert cell_count_plateau(history) == 50.0


def test_cell_count_plateau_empty_returns_nan():
    assert math.isnan(cell_count_plateau([]))


def test_vertex_envelope_returns_min_max_per_axis():
    positions = np.array([
        [-1.0, 2.0, 3.5],
        [4.0, -2.0, 0.5],
        [0.0, 1.0, 7.5],
    ])
    env = vertex_envelope(positions)
    assert env['x_min'] == -1.0 and env['x_max'] == 4.0
    assert env['y_min'] == -2.0 and env['y_max'] == 2.0
    assert env['z_min'] == 0.5 and env['z_max'] == 7.5


def test_vertex_envelope_empty_returns_nan_dict():
    env = vertex_envelope(np.zeros((0, 3)))
    assert all(math.isnan(v) for v in env.values())


def test_regime_nan_for_nonfinite_history():
    assert regime([10, float('nan'), 20]) == 'NaN'


def test_regime_collapse_when_final_below_peak_ratio():
    # Peak = 100, final = 30, ratio = 0.5 → collapse.
    assert regime([10, 50, 100, 80, 30]) == 'collapse'


def test_regime_oscillate_for_alternating_diffs():
    # Last half is [50, 60, 50, 60, 50]; diffs alternate sign repeatedly.
    history = [10, 20, 30, 40, 50, 60, 50, 60, 50]
    assert regime(history) == 'oscillate'


def test_regime_plateau_for_steady_tail():
    # Last 20% of 50 entries = last 10 entries; all equal -> std=0.
    history = list(range(40)) + [40] * 10
    assert regime(history) == 'plateau'


def test_regime_monotone_grow_for_steady_increase():
    history = list(range(1, 21))  # strictly monotone, no plateau in tail
    assert regime(history) == 'monotone-grow'


def test_regime_empty_history_returns_nan():
    assert regime([]) == 'NaN'
