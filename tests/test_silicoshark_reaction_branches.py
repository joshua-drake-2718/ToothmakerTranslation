"""Unit tests for Path B v2 reaction-term and differentiation branches.

Each test exercises one `Discretisation` field by constructing a
minimal `State` + `Params`, running one step of
`step_reaction_diffusion` (or `step_differentiation`) under each
option of the field, and asserting that the resulting state arrays
agree with the analytical expectation.

These tests are mathematical, not biological: they pin down the
forward-Euler arithmetic of a single step on small (1- or 7-cell)
states. The biological-regression tests live in
`tests/test_simulator_smoke.py` and (eventually)
`tests/test_silicoshark_smoke.py`.
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from silicoshark.discretisation import PATH_B_DEFAULT
from silicoshark.mesh import Mesh
from silicoshark.params import Params
from silicoshark.reaction import step_reaction_diffusion
from silicoshark.simulator import step_differentiation
from silicoshark.state import N_MESENCHYME_LAYERS, State


# --- Fixture builders -------------------------------------------------


def _single_cell_state(
    *,
    act: float = 0.0,
    inh: float = 0.0,
    sec: float = 0.0,
    diff: float = 0.0,
    knot: bool = False,
    first_border_cell: int = 0,
) -> State:
    """One isolated cell at the origin. Mesh built from three colinear
    helpers so Delaunay can triangulate, but the helpers are far
    enough away that the focal cell's diffusion contribution from
    them is irrelevant when their concentrations are zero.

    For the reaction-side tests we drive `k_da = k_di = k_ds = 0` so
    diffusion is inert and only the local reaction terms matter. The
    mesh is still required because `step_reaction_diffusion` calls
    `mesh.laplacian` on each species.
    """
    # 4 cells: focal at origin, three at (10, 0), (0, 10), (10, 10).
    # Far enough apart that they don't crash the triangulation.
    n = 4
    positions = np.array(
        [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [10.0, 10.0, 0.0]],
        dtype=np.float64,
    )
    act_arr = np.zeros(n, dtype=np.float64)
    inh_arr = np.zeros(n, dtype=np.float64)
    sec_arr = np.zeros(n, dtype=np.float64)
    diff_arr = np.zeros(n, dtype=np.float64)
    knot_arr = np.zeros(n, dtype=bool)
    act_arr[0] = act
    inh_arr[0] = inh
    sec_arr[0] = sec
    diff_arr[0] = diff
    knot_arr[0] = knot
    return State(
        positions=positions,
        act=act_arr,
        inh=inh_arr,
        sec=sec_arr,
        diff=diff_arr,
        knot=knot_arr,
        mes_inh=np.zeros((n, N_MESENCHYME_LAYERS)),
        mes_sec=np.zeros((n, N_MESENCHYME_LAYERS)),
        init_anterior=np.zeros(n, dtype=bool),
        init_posterior=np.zeros(n, dtype=bool),
        init_lingual=np.zeros(n, dtype=bool),
        init_buccal=np.zeros(n, dtype=bool),
        first_border_cell=first_border_cell,
    )


def _diffusion_off_params(**overrides) -> Params:
    """Params with all diffusion coefficients zero; only reaction
    terms matter. Caller overrides the reaction parameters they care
    about.
    """
    p = Params(
        k_act=0.0,
        k_inh=0.0,
        k_sec=0.0,
        k_da=0.0,
        k_di=0.0,
        k_ds=0.0,
        k_deg=0.0,
        k_dff=0.0,
        k_int=0.0,
        k_set=0.0,
    )
    for k, v in overrides.items():
        setattr(p, k, v)
    return p


# --- Eq. 14 denominator ----------------------------------------------


def test_eq14_denominator_branches_differ():
    """With Act=0.3, Inh=0.5, k_act=1.0, k_inh=2.0, k_deg=0.1:

      inh_corrected: 1.0 * 0.3 / (1 + 2.0 * 0.5) - 0.1 * 0.3
                   = 0.3 / 2.0 - 0.03 = 0.15 - 0.03 = 0.12
      act_typo:      1.0 * 0.3 / (1 + 2.0 * 0.3) - 0.1 * 0.3
                   = 0.3 / 1.6 - 0.03 = 0.1875 - 0.03 = 0.1575

    Expected new Act with dt=1.0: 0.3 + rd_act, so 0.42 vs 0.4575.
    """
    params = _diffusion_off_params(k_act=1.0, k_inh=2.0, k_deg=0.1)
    dt = 1.0

    state_corr = _single_cell_state(act=0.3, inh=0.5)
    mesh = Mesh.from_positions(state_corr.positions)
    step_reaction_diffusion(
        state_corr, params, mesh, dt,
        replace(PATH_B_DEFAULT, eq14_denominator='inh_corrected'),
    )
    assert state_corr.act[0] == pytest.approx(0.3 + 0.12, abs=1e-12)

    state_typo = _single_cell_state(act=0.3, inh=0.5)
    mesh = Mesh.from_positions(state_typo.positions)
    step_reaction_diffusion(
        state_typo, params, mesh, dt,
        replace(PATH_B_DEFAULT, eq14_denominator='act_typo'),
    )
    assert state_typo.act[0] == pytest.approx(0.3 + 0.1575, abs=1e-12)

    assert state_corr.act[0] != state_typo.act[0]


# --- Eq. 17 inhibitor source -----------------------------------------


def test_eq17_inh_source_branches_differ():
    """1-cell state with knot=False, diff=0.5 (> Int=0.4), Act=0.7.

      act_concentration: rd_inh = 0.7 - Deg * Inh
      act_rate_times_di: rd_inh = (rate-of-Act) * 0.5 - Deg * Inh

    Where rate-of-Act = k_act * Act / (1 + k_inh * Inh) - k_deg * Act.
    With k_act=1.0, k_inh=2.0, k_deg=0.1, Act=0.7, Inh=0.2:
      rate-of-Act = 1.0 * 0.7 / (1 + 2.0 * 0.2) - 0.1 * 0.7
                  = 0.7 / 1.4 - 0.07 = 0.5 - 0.07 = 0.43
      act_concentration: rd_inh = 0.7 - 0.1 * 0.2 = 0.68
      act_rate_times_di: rd_inh = 0.43 * 0.5 - 0.1 * 0.2 = 0.215 - 0.02 = 0.195

    With dt=1.0: new Inh = 0.2 + rd_inh, so 0.88 vs 0.395.
    """
    params = _diffusion_off_params(
        k_act=1.0, k_inh=2.0, k_deg=0.1,
        k_int=0.4, k_set=2.0,  # k_set high so eq.18 doesn't fire
    )
    dt = 1.0

    state_paper = _single_cell_state(act=0.7, inh=0.2, diff=0.5, knot=False)
    mesh = Mesh.from_positions(state_paper.positions)
    step_reaction_diffusion(
        state_paper, params, mesh, dt,
        replace(PATH_B_DEFAULT, eq17_inh_source='act_concentration'),
    )
    assert state_paper.inh[0] == pytest.approx(0.2 + 0.68, abs=1e-12)

    state_fortran = _single_cell_state(act=0.7, inh=0.2, diff=0.5, knot=False)
    mesh = Mesh.from_positions(state_fortran.positions)
    step_reaction_diffusion(
        state_fortran, params, mesh, dt,
        replace(PATH_B_DEFAULT, eq17_inh_source='act_rate_times_di'),
    )
    assert state_fortran.inh[0] == pytest.approx(0.2 + 0.195, abs=1e-12)

    assert state_paper.inh[0] != state_fortran.inh[0]


# --- Eq. 18 secondary-signal source ----------------------------------


def test_eq18_sec_source_branches_differ():
    """1-cell state with diff=0.6 (> Set=0.5), Sec=0.1.

      constant_k_sec: rd_sec = k_sec - Deg * Sec
      k_sec_times_di: rd_sec = k_sec * 0.6 - Deg * Sec

    With k_sec=0.5, k_deg=0.1, Sec=0.1, dt=1.0:
      constant: rd_sec = 0.5 - 0.01 = 0.49 -> new Sec = 0.59
      times_di: rd_sec = 0.5 * 0.6 - 0.01 = 0.29 -> new Sec = 0.39
    """
    params = _diffusion_off_params(
        k_sec=0.5, k_deg=0.1,
        k_int=2.0, k_set=0.5,  # k_int high so eq.17 doesn't fire
    )
    dt = 1.0

    state_paper = _single_cell_state(sec=0.1, diff=0.6, knot=False)
    mesh = Mesh.from_positions(state_paper.positions)
    step_reaction_diffusion(
        state_paper, params, mesh, dt,
        replace(PATH_B_DEFAULT, eq18_sec_source='constant_k_sec'),
    )
    assert state_paper.sec[0] == pytest.approx(0.59, abs=1e-12)

    state_fortran = _single_cell_state(sec=0.1, diff=0.6, knot=False)
    mesh = Mesh.from_positions(state_fortran.positions)
    step_reaction_diffusion(
        state_fortran, params, mesh, dt,
        replace(PATH_B_DEFAULT, eq18_sec_source='k_sec_times_di'),
    )
    assert state_fortran.sec[0] == pytest.approx(0.39, abs=1e-12)

    assert state_paper.sec[0] != state_fortran.sec[0]


# --- Differentiation accumulator -------------------------------------


def test_diff_accumulator_branches_differ():
    """State with Act=0.4, Sec=0.1, k_dff=0.5, dt=1.0.

      'act' branch: diff_after = 0 + 1.0 * 0.5 * 0.4 = 0.2
      'sec' branch: diff_after = 0 + 1.0 * 0.5 * 0.1 = 0.05
    """
    params = _diffusion_off_params(k_dff=0.5)
    dt = 1.0

    state_act = _single_cell_state(act=0.4, sec=0.1, diff=0.0)
    step_differentiation(
        state_act, params, dt,
        replace(PATH_B_DEFAULT, diff_accumulator='act'),
    )
    assert state_act.diff[0] == pytest.approx(0.2, abs=1e-12)

    state_sec = _single_cell_state(act=0.4, sec=0.1, diff=0.0)
    step_differentiation(
        state_sec, params, dt,
        replace(PATH_B_DEFAULT, diff_accumulator='sec'),
    )
    assert state_sec.diff[0] == pytest.approx(0.05, abs=1e-12)

    assert state_act.diff[0] != state_sec.diff[0]


# --- Knot threshold gate ---------------------------------------------


def _seven_cell_state(*, act_cell0: float, first_border_cell: int) -> State:
    """7-cell rad=2 state. Positions are the standard hex lattice;
    all concentrations zero except `act_cell0` at index 0.
    """
    # Use the actual hex lattice (rad=2 → 7 cells) so the geometry
    # is realistic; the test only inspects `state.knot` after one
    # step, so the positions don't materially matter beyond letting
    # Delaunay triangulate.
    from silicoshark.state import hex_lattice
    xy = hex_lattice(2)
    n = xy.shape[0]
    assert n == 7
    positions = np.zeros((n, 3), dtype=np.float64)
    positions[:, :2] = xy
    positions[:, 2] = 1.0
    act = np.zeros(n, dtype=np.float64)
    act[0] = act_cell0
    return State(
        positions=positions,
        act=act,
        inh=np.zeros(n),
        sec=np.zeros(n),
        diff=np.zeros(n),
        knot=np.zeros(n, dtype=bool),
        mes_inh=np.zeros((n, N_MESENCHYME_LAYERS)),
        mes_sec=np.zeros((n, N_MESENCHYME_LAYERS)),
        init_anterior=np.zeros(n, dtype=bool),
        init_posterior=np.zeros(n, dtype=bool),
        init_lingual=np.zeros(n, dtype=bool),
        init_buccal=np.zeros(n, dtype=bool),
        first_border_cell=first_border_cell,
    )


def test_knot_threshold_gate_branches_differ():
    """7-cell rad=2 state, Act=2.0 in cell 0 only, first_border_cell=1
    (= 1 + 3*(rad-1)*(rad-2) for rad=2; rad=2 has only the centre as
    interior, and the 6 surrounding cells as the outer ring).

      'none': cell 0 (Act >= 1) becomes a knot.
      'first_border_cell': cell 0 (index 0 < first_border_cell=1) IS
        an interior cell, so it DOES become a knot. The gate excludes
        cells at indices >= first_border_cell (the outer ring) from
        knot formation. With first_border_cell=1, only cell 0 (the
        centre) is interior; the 6 ring-1 cells are border.

    Cell 5 is in the outer ring (index 5 >= 1) so we additionally
    seed it with Act=2.0 and assert that under 'first_border_cell'
    cell 5 does NOT become a knot.

    All reaction parameters set to zero so the only effect of the
    step is the knot-detection clause and a degraded copy of Act
    (here, with k_deg=0, Act stays at 2.0).
    """
    params = _diffusion_off_params()  # all zero
    dt = 1.0
    fbc = 1 + 3 * (2 - 1) * (2 - 2)  # = 1 (centre only is interior)

    state_none = _seven_cell_state(act_cell0=2.0, first_border_cell=fbc)
    state_none.act[5] = 2.0  # additionally seed an outer-ring cell
    mesh = Mesh.from_positions(state_none.positions)
    step_reaction_diffusion(
        state_none, params, mesh, dt,
        replace(PATH_B_DEFAULT, knot_threshold_gate='none'),
    )
    assert state_none.knot[0] == True, 'paper rule: cell 0 should be a knot'
    assert state_none.knot[5] == True, 'paper rule: cell 5 should be a knot'

    state_gate = _seven_cell_state(act_cell0=2.0, first_border_cell=fbc)
    state_gate.act[5] = 2.0
    mesh = Mesh.from_positions(state_gate.positions)
    step_reaction_diffusion(
        state_gate, params, mesh, dt,
        replace(PATH_B_DEFAULT, knot_threshold_gate='first_border_cell'),
    )
    assert state_gate.knot[0] == True, (
        'FORTRAN rule: cell 0 (index 0 < first_border_cell=1) IS '
        'interior, so it should become a knot.'
    )
    assert state_gate.knot[5] == False, (
        'FORTRAN rule: cell 5 (index 5 >= first_border_cell=1) is on '
        'the outer ring, so the knot gate should exclude it.'
    )


# --- Default-preset equivalence to v1 (sanity check) -----------------


def test_default_disc_eq_omitted_disc():
    """`PATH_B_DEFAULT` must be the implicit default. Calling
    `step_reaction_diffusion` without a `disc` argument must produce
    the same result as passing `PATH_B_DEFAULT` explicitly.
    """
    params = _diffusion_off_params(k_act=1.0, k_inh=2.0, k_deg=0.1)
    dt = 1.0

    state_implicit = _single_cell_state(act=0.3, inh=0.5)
    mesh = Mesh.from_positions(state_implicit.positions)
    step_reaction_diffusion(state_implicit, params, mesh, dt)

    state_explicit = _single_cell_state(act=0.3, inh=0.5)
    mesh = Mesh.from_positions(state_explicit.positions)
    step_reaction_diffusion(state_explicit, params, mesh, dt, PATH_B_DEFAULT)

    np.testing.assert_array_equal(state_implicit.act, state_explicit.act)
    np.testing.assert_array_equal(state_implicit.inh, state_explicit.inh)
    np.testing.assert_array_equal(state_implicit.sec, state_explicit.sec)
