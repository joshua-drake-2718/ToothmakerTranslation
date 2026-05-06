"""Unit tests for `Mesh.fortran_margins_laplacian`.

The operator is the FORTRAN-style per-edge `pes`-weighted Laplacian
(`silicoshark/mesh.py:fortran_margins_laplacian`), introduced in Path B
v2 B3 to close rebuttal §B.3's structural-no-op gap. The tests exercise:

  - constant-field annihilation (any Laplacian's most basic invariant);
  - the analytical pes/area_p values on a regular 7-cell hex lattice;
  - sign and symmetry on a linear field gradient;
  - degenerate cells with < 3 neighbours contribute no flux;
  - the positions=None error path;
  - end-to-end CLI smoke on the cusp-forming preset.

The full FORTRAN apply_diffusion (substrate sink, vertical-z flux,
substrate-edge layer) is NOT byte-matched here — that scope is
deferred per the May 2026 architectural decision documented in
`docs/plans/2026-05-05-fortran-margins-laplacian.md`.
"""
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from silicoshark.mesh import Mesh
from silicoshark.topology import Topology

REPO_ROOT = Path(__file__).resolve().parent.parent


def _hex_7_positions() -> np.ndarray:
    """Centre + 6 unit-distance hexagonal ring, in the xy plane at z=0."""
    pts = [(0.0, 0.0, 0.0)]
    for k in range(6):
        angle = math.radians(60.0 * k)
        pts.append((math.cos(angle), math.sin(angle), 0.0))
    return np.asarray(pts, dtype=np.float64)


def _hex_7_mesh() -> Mesh:
    pos = _hex_7_positions()
    top = Topology.from_positions(pos)
    return Mesh.from_topology(top, pos, kind='fortran_margins')


def test_constant_field_annihilates():
    mesh = _hex_7_mesh()
    u = np.full(7, 3.14159, dtype=np.float64)
    out = mesh.fortran_margins_laplacian(u)
    assert np.allclose(out, 0.0, atol=1e-14), out


def test_centre_cell_pes_normed_matches_analytical():
    """For the centre of a regular hex lattice with 6 unit-distance
    neighbours, the pes weight applied to each ring neighbour can be
    derived analytically. With margin midpoints at radius 0.5 spaced
    60° apart:

        pes_per_seg = 2 * 0.5 * sin(30°) = 0.5
        area_p_per  = 0.5 * |(0.5, 0) x (0.25, 0.25*sqrt(3))| = sqrt(3)/16
        area_bottom = 6 * sqrt(3)/16 = 3*sqrt(3)/8
        sum_a       = 6 * 0.5 + 2 * 3*sqrt(3)/8 = 3 + 3*sqrt(3)/4
        pes_normed  = 0.5 / sum_a

    The Laplacian at the centre with u[ring]=1, u[centre]=0 is then
    6 * pes_normed.
    """
    mesh = _hex_7_mesh()
    u = np.zeros(7, dtype=np.float64)
    u[1:] = 1.0
    out = mesh.fortran_margins_laplacian(u)

    sum_a = 3.0 + 3.0 * math.sqrt(3.0) / 4.0
    pes_normed = 0.5 / sum_a
    expected_centre = 6.0 * pes_normed
    assert math.isclose(out[0], expected_centre, rel_tol=1e-12), (
        out[0], expected_centre,
    )


def test_centre_cell_sign_flips_with_field():
    mesh = _hex_7_mesh()
    u_up = np.zeros(7, dtype=np.float64)
    u_up[1:] = 1.0
    u_down = -u_up
    out_up = mesh.fortran_margins_laplacian(u_up)
    out_down = mesh.fortran_margins_laplacian(u_down)
    assert math.isclose(out_up[0], -out_down[0], rel_tol=1e-14)
    assert out_up[0] > 0.0
    assert out_down[0] < 0.0


def test_linear_field_centre_is_zero():
    """A field linear in x has zero Laplacian at the centre by symmetry:
    the contributions from the +x neighbour and the -x neighbour cancel,
    and the off-axis ring members pair up symmetrically. The pes weights
    are equal across all 6 ring slots on the regular hex, so the sum
    `sum_k pes_normed * (u[k] - u[centre])` reduces to
    `pes_normed * sum_k u[k]`, which is zero for u(x,y) = x.
    """
    mesh = _hex_7_mesh()
    pos = _hex_7_positions()
    u = pos[:, 0].copy()
    out = mesh.fortran_margins_laplacian(u)
    assert math.isclose(out[0], 0.0, abs_tol=1e-13), out[0]


def test_degenerate_cell_under_3_neighbours_contributes_zero():
    """A cell with < 3 neighbours has no closed polygon and the operator
    returns zero flux there. We construct a Mesh by hand to exercise
    this branch (Topology.from_positions over a 3-cell line wouldn't
    triangulate cleanly; we go straight to Mesh).
    """
    n = 3
    pos = np.array(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        dtype=np.float64,
    )
    # Build a degenerate CSR by hand: cell 0 has neighbour 1; cell 1
    # has neighbours 0 and 2; cell 2 has neighbour 1.
    neigh_starts = np.array([0, 1, 3, 4], dtype=np.int64)
    neigh_idx = np.array([1, 0, 2, 1], dtype=np.int32)
    edge_w = np.zeros(neigh_idx.shape[0], dtype=np.float64)
    diag_w = np.zeros(n, dtype=np.float64)
    is_border = np.array([True, True, True])
    triangles = np.zeros((0, 3), dtype=np.int32)
    mesh = Mesh(
        triangles=triangles,
        neigh_starts=neigh_starts,
        neigh_idx=neigh_idx,
        edge_w=edge_w,
        diag_w=diag_w,
        is_border=is_border,
        kind='fortran_margins',
        positions=pos,
    )
    u = np.array([1.0, 2.0, 3.0], dtype=np.float64)
    out = mesh.fortran_margins_laplacian(u)
    assert np.allclose(out, 0.0), out


def test_positions_none_raises():
    mesh = _hex_7_mesh()
    mesh.positions = None
    with pytest.raises(ValueError, match='requires positions'):
        mesh.laplacian(np.zeros(7))


def test_cli_runs_with_fortran_margins():
    """End-to-end smoke: the CLI accepts laplacian=fortran_margins and
    produces 5 OFF files with finite vertex coordinates.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / 'out'
        out.mkdir()
        args = [
            sys.executable, '-m', 'silicoshark',
            str(REPO_ROOT / 'examples' / 'wt-tribosphenic-2014-paper-ina.txt'),
            str(out), 'run', '100', '5',
            '--preset', 'PAPER_2010',
            '--override', 'mesenchyme=absent',
            '--override', 'laplacian=fortran_margins',
        ]
        proc = subprocess.run(args, cwd=REPO_ROOT, capture_output=True, text=True)
        assert proc.returncode == 0, (
            f'silicoshark exited {proc.returncode}\n'
            f'stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}'
        )
        final = out / '500_run_.off'
        assert final.exists(), f'no final OFF at {final}'
        # Parse header and first vertex line; check finiteness.
        with final.open() as f:
            f.readline()  # COFF
            n_v, n_f, _ = f.readline().split()
            n_v = int(n_v)
            for _ in range(n_v):
                line = f.readline()
                xyz = [float(x) for x in line.split()[:3]]
                assert all(math.isfinite(v) for v in xyz), line
        assert n_v == 19, f'expected 19 vertices, got {n_v}'


def test_legacy_fortran_byte_identical_under_either_laplacian():
    """Regression: LEGACY_FORTRAN on `wt-tribosphenic-2014.txt` must
    produce byte-identical OFF outputs under `laplacian=length_weighted`
    and `laplacian=fortran_margins`. This locks in the empirical null
    on the knock-down direction of the disentanglement matrix and
    catches any future change to the in-plane operator that
    accidentally drifts LEGACY_FORTRAN's saturated-regime output.

    See findings doc
    `docs/findings/2026-05-05-path-b-v2-fortran-margins-implementation.md`
    and the §B.3 rebuttal entry for the rationale.

    Iteration budget kept short (100×3 = 300 iters) for fast-suite
    inclusion. The 2026-05-05 byte-identicality confirmation was on
    500×5 = 2500 iters; if that ever drifts, this 300-iter test will
    pick it up first.
    """
    import hashlib
    import tempfile

    def _run_and_hash(out_dir: Path, laplacian: str) -> dict[str, str]:
        args = [
            sys.executable, '-m', 'silicoshark',
            str(REPO_ROOT / 'examples' / 'wt-tribosphenic-2014.txt'),
            str(out_dir), 'run', '100', '3',
            '--preset', 'LEGACY_FORTRAN',
            '--override', f'laplacian={laplacian}',
            '--override', 'mesenchyme=absent',
        ]
        proc = subprocess.run(args, cwd=REPO_ROOT, capture_output=True, text=True)
        assert proc.returncode == 0, (
            f'silicoshark with laplacian={laplacian} exited {proc.returncode}\n'
            f'stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}'
        )
        hashes = {}
        for save in (100, 200, 300):
            off = out_dir / f'{save}_run_.off'
            assert off.exists(), f'missing OFF at {off}'
            hashes[save] = hashlib.sha256(off.read_bytes()).hexdigest()
        return hashes

    with tempfile.TemporaryDirectory() as tmp:
        out_lw = Path(tmp) / 'lw'
        out_fm = Path(tmp) / 'fm'
        out_lw.mkdir()
        out_fm.mkdir()
        h_lw = _run_and_hash(out_lw, 'length_weighted')
        h_fm = _run_and_hash(out_fm, 'fortran_margins')
    assert h_lw == h_fm, (
        f'LEGACY_FORTRAN OFF outputs differ between length_weighted and '
        f'fortran_margins:\n  length_weighted: {h_lw}\n  fortran_margins: {h_fm}'
    )
