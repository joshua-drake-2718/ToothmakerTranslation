"""Unit tests for Path B v2 B3: per-column-z-layers mesenchyme model.

The mesenchyme branch couples the epithelial Inh and Sec fields to a
column of mesenchymal layers underneath each epithelial cell. Per the
2010 paper p. 586:

  - Mesenchymal cells have NO reaction terms (only diffusion +
    degradation).
  - Inh and Sec diffuse vertically between epithelium and the top
    mesenchymal layer, and between adjacent mesenchymal layers.
  - Act has NO mesenchymal counterpart (no vertical Act diffusion).
  - Cervical-loop thickening (eqs. 10–12) and buoyancy (eq. 13) read
    Sec from the top mesenchymal layer, not the epithelium.

These tests pin down the arithmetic of the implementation in
`silicoshark/reaction.py` (`step_mesenchyme_diffusion` plus the
vertical-coupling term added to `step_reaction_diffusion`) and the
mesenchyme-aware reads in `silicoshark/forces.py` (cervical-loop
thickening and buoyancy).
"""
from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from silicoshark.discretisation import (
    LEGACY_FORTRAN, PATH_B_DEFAULT,
)
from silicoshark.forces import compute_forces
from silicoshark.mesh import Mesh
from silicoshark.params import Params
from silicoshark.reaction import step_reaction_diffusion
from silicoshark.state import N_MESENCHYME_LAYERS, State, build_initial_state


# --- Fixture builders -------------------------------------------------


def _zero_params(**overrides) -> Params:
    """Params with every coefficient zero. Caller overrides the
    diffusion / reaction coefficients they care about.
    """
    p = Params(
        k_act=0.0, k_inh=0.0, k_sec=0.0,
        k_da=0.0, k_di=0.0, k_ds=0.0,
        k_deg=0.0, k_dff=0.0,
        k_int=2.0, k_set=2.0,  # high so eq. 17/18 sources never fire
    )
    for k, v in overrides.items():
        setattr(p, k, v)
    return p


def _seven_cell_state() -> State:
    """rad=2 hex lattice. Returns a State with all concentrations zero
    and the topology / first_border_cell consistent with build_initial_state.
    """
    p = Params(rad=2)
    return build_initial_state(p)


# --- Test 1: 'absent' is no-op ---------------------------------------


def test_absent_keeps_mesenchyme_zero_after_many_steps():
    """Under disc.mesenchyme='absent', mes_inh and mes_sec must stay
    bit-for-bit zero across many iterations of step_reaction_diffusion,
    and state.inh / state.sec must be unaffected by the mesenchyme
    branch.
    """
    params = Params.from_file('examples/seal.txt')
    state = build_initial_state(params)
    mesh = Mesh.from_positions(state.positions)
    disc = replace(PATH_B_DEFAULT, mesenchyme='absent')

    # Seed the epithelial Inh on one cell so there's something to track.
    state.inh[0] = 1.0

    for _ in range(100):
        step_reaction_diffusion(state, params, mesh, 0.05, disc)

    np.testing.assert_array_equal(
        state.mes_inh, np.zeros_like(state.mes_inh),
        err_msg='mes_inh must stay zero under mesenchyme=absent',
    )
    np.testing.assert_array_equal(
        state.mes_sec, np.zeros_like(state.mes_sec),
        err_msg='mes_sec must stay zero under mesenchyme=absent',
    )


# --- Test 2: Mass conservation under per_column_z_layers -------------


def test_mass_conservation_inh_under_per_column_z_layers():
    """Total Inh mass (sum of state.inh + sum of state.mes_inh) must
    be conserved within 1e-6 over 100 iterations under
    disc.mesenchyme='per_column_z_layers' with k_deg=0 and no reaction
    sources (no knots, no diff>=Int). Diffusion alone shuffles mass
    between cells and layers; total stays constant.
    """
    state = _seven_cell_state()
    # Seed Inh on the centre cell only (epithelial layer); other cells
    # and all mesenchyme layers start at zero.
    state.inh[0] = 1.0
    initial_total = state.inh.sum() + state.mes_inh.sum()

    params = _zero_params(k_di=0.5, k_ds=0.5)  # diffusion only, no degradation, no reaction
    mesh = Mesh.from_positions(state.positions)
    disc = replace(PATH_B_DEFAULT, mesenchyme='per_column_z_layers')

    for _ in range(100):
        step_reaction_diffusion(state, params, mesh, 0.01, disc)

    final_total = state.inh.sum() + state.mes_inh.sum()
    assert abs(final_total - initial_total) < 1e-6, (
        f'Inh mass not conserved: initial={initial_total:.10f}, '
        f'final={final_total:.10f}, drift={final_total - initial_total:.3e}'
    )


def test_mass_conservation_sec_under_per_column_z_layers():
    """Same as above but for Sec (state.sec + state.mes_sec)."""
    state = _seven_cell_state()
    state.sec[0] = 1.0
    initial_total = state.sec.sum() + state.mes_sec.sum()

    params = _zero_params(k_di=0.5, k_ds=0.5)
    mesh = Mesh.from_positions(state.positions)
    disc = replace(PATH_B_DEFAULT, mesenchyme='per_column_z_layers')

    for _ in range(100):
        step_reaction_diffusion(state, params, mesh, 0.01, disc)

    final_total = state.sec.sum() + state.mes_sec.sum()
    assert abs(final_total - initial_total) < 1e-6, (
        f'Sec mass not conserved: initial={initial_total:.10f}, '
        f'final={final_total:.10f}, drift={final_total - initial_total:.3e}'
    )


# --- Test 3: Diffusion to / from mesenchyme --------------------------


def test_inh_leaks_from_epi_to_mesenchyme_in_one_step():
    """With state.inh[0] = 1.0 and all other concentrations zero,
    after one step under per_column_z_layers:
      - mes_inh[0, 0] should be > 0 (some Inh leaked down into the
        top mesenchyme layer).
      - state.inh[0] should have decreased (it lost mass to mes layer 0
        and to its lateral neighbours via mesh.laplacian, and gained
        nothing from the still-zero mes_inh[0, 0]).
    """
    state = _seven_cell_state()
    state.inh[0] = 1.0
    starting_inh_0 = state.inh[0]

    params = _zero_params(k_di=0.5, k_ds=0.5)
    mesh = Mesh.from_positions(state.positions)
    disc = replace(PATH_B_DEFAULT, mesenchyme='per_column_z_layers')

    step_reaction_diffusion(state, params, mesh, 0.01, disc)

    assert state.mes_inh[0, 0] > 0, (
        'mes_inh[0, 0] should be > 0 after one step (Inh leaks down)'
    )
    assert state.inh[0] < starting_inh_0, (
        f'state.inh[0] should decrease from {starting_inh_0}, '
        f'got {state.inh[0]}'
    )


def test_inh_leaks_from_mesenchyme_to_epi_in_one_step():
    """Symmetric: with mes_inh[0, 0] = 1.0 and all other concentrations
    zero, after one step under per_column_z_layers:
      - state.inh[0] should be > 0 (Inh leaked up from layer 0 to epi).
      - mes_inh[0, 0] should have decreased (it lost mass upward to
        epi and downward to mes_inh[0, 1]).
    """
    state = _seven_cell_state()
    state.mes_inh[0, 0] = 1.0
    starting_mes_inh_0_0 = state.mes_inh[0, 0]

    params = _zero_params(k_di=0.5, k_ds=0.5)
    mesh = Mesh.from_positions(state.positions)
    disc = replace(PATH_B_DEFAULT, mesenchyme='per_column_z_layers')

    step_reaction_diffusion(state, params, mesh, 0.01, disc)

    assert state.inh[0] > 0, (
        'state.inh[0] should be > 0 after one step (Inh leaks up from '
        'layer 0 mesenchyme)'
    )
    assert state.mes_inh[0, 0] < starting_mes_inh_0_0, (
        f'mes_inh[0, 0] should decrease from {starting_mes_inh_0_0}, '
        f'got {state.mes_inh[0, 0]}'
    )


def test_inh_propagates_from_mes_top_to_mes_bottom():
    """With mes_inh[0, 0] = 1.0 (top layer) and all else zero, after
    one step the bottom layer mes_inh[0, 1] should be > 0 (mass has
    moved down between layers via vertical Laplacian).
    """
    state = _seven_cell_state()
    state.mes_inh[0, 0] = 1.0

    params = _zero_params(k_di=0.5, k_ds=0.5)
    mesh = Mesh.from_positions(state.positions)
    disc = replace(PATH_B_DEFAULT, mesenchyme='per_column_z_layers')

    step_reaction_diffusion(state, params, mesh, 0.01, disc)

    assert state.mes_inh[0, 1] > 0, (
        'mes_inh[0, 1] should be > 0 after one step (Inh propagates '
        'downward between mesenchyme layers)'
    )


# --- Test 4: Cervical-loop reads mesenchyme under per_column_z_layers


def test_cervical_loop_reads_mesenchyme_sec_under_per_column():
    """Build a state where state.sec is zero but mes_sec[border, 0] is
    non-zero. Forces under disc.mesenchyme='per_column_z_layers' should
    reflect the mes_sec value in the cervical-loop thickening factor;
    forces under 'absent' should not.

    The seal initial state has zero Sec everywhere; we keep state.sec
    at zero and seed mes_sec on border cells. Under absent, the
    thickening factor uses state.sec=0 so thicken = 1 + k_umgr/out_h.
    Under per_column_z_layers, the factor uses mes_sec=non-zero so
    thicken differs.
    """
    # Use LEGACY_FORTRAN as the disc base so we have border cells with
    # non-trivial cervical-loop forces (and `border_definition='topological_descendants'`
    # so the border set is well-defined).
    params = Params.from_file('examples/seal.txt')
    state = build_initial_state(params, lattice_orientation='fortran')
    # k_mgr must be non-zero for the Sec source in thicken to matter;
    # peek at the seal params.
    if params.k_mgr == 0:
        params.k_mgr = 1.0  # paper-typical: ensure the test exercises the path

    mesh = Mesh.from_positions(state.positions)
    border = mesh.is_border  # geometric border for the test
    assert border.any(), 'expected at least one border cell'

    # Seed mes_sec on border cells only, layer 0.
    state.mes_sec[border, 0] = 0.5

    disc_absent = replace(LEGACY_FORTRAN, mesenchyme='absent', laplacian='length_weighted')
    disc_per_col = replace(LEGACY_FORTRAN, mesenchyme='per_column_z_layers', laplacian='length_weighted')

    f_absent = compute_forces(state, params, mesh, disc_absent)
    f_per_col = compute_forces(state, params, mesh, disc_per_col)

    # The xy magnitudes on border cells must differ between branches:
    # under per_column_z_layers, thicken is amplified by mes_sec; under
    # absent, thicken sees state.sec = 0.
    border_xy_absent = np.linalg.norm(f_absent[border, :2], axis=1)
    border_xy_per_col = np.linalg.norm(f_per_col[border, :2], axis=1)
    assert not np.allclose(border_xy_absent, border_xy_per_col), (
        'cervical-loop xy magnitudes on border cells must differ between '
        'mesenchyme=absent and mesenchyme=per_column_z_layers when '
        'mes_sec[border, 0] is non-zero'
    )


def test_buoyancy_reads_mesenchyme_sec_under_per_column():
    """Buoyancy (eq. 13) under 'absent' reads state.sec; under
    'per_column_z_layers' it reads mes_sec[:, 0]. With state.sec=0 but
    mes_sec[:, 0] = non-zero, the per_column_z_layers branch produces
    a non-zero z-force component on cells with mes_sec, while 'absent'
    produces zero buoyancy.
    """
    params = Params.from_file('examples/seal.txt')
    state = build_initial_state(params)
    if params.k_boy == 0:
        params.k_boy = 1.0  # ensure the buoyancy term has a coefficient

    mesh = Mesh.from_positions(state.positions)
    # Seed mes_sec on every cell, layer 0.
    state.mes_sec[:, 0] = 0.3

    disc_absent = replace(PATH_B_DEFAULT, mesenchyme='absent')
    disc_per_col = replace(PATH_B_DEFAULT, mesenchyme='per_column_z_layers')

    f_absent = compute_forces(state, params, mesh, disc_absent)
    f_per_col = compute_forces(state, params, mesh, disc_per_col)

    # The buoyancy term contributes only to z-force. Under absent
    # (state.sec=0) it is zero; under per_column_z_layers it is
    # non-zero (because mes_sec[:, 0] = 0.3).
    z_diff = f_per_col[:, 2] - f_absent[:, 2]
    assert np.any(np.abs(z_diff) > 1e-9), (
        'buoyancy z-force must differ between branches when state.sec=0 '
        'but mes_sec[:, 0] is non-zero'
    )
