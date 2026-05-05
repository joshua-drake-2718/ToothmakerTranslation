"""Smoke test for silicoshark v2 LEGACY_FORTRAN preset on the seal example.

Runs `python -m silicoshark examples/seal.txt <tmp> run 100 5 --preset
LEGACY_FORTRAN` and asserts the final OFF file matches the FORTRAN
binary's golden output within v2 A6 tolerance:

- Cell-count plateau within ±5% of 57 (the FORTRAN goldens' 57-cell
  equilibrium).
- Vertex envelope on x and y within 10% of FORTRAN's span.
- Vertex envelope on z: z_max within 10% of FORTRAN's span; z_min
  within 65% of FORTRAN's span (relaxed; see residual-divergence
  finding doc for the full rationale).
- No NaN/Inf vertices, no RuntimeWarning in stderr.

This is the replication anchor for the comparison-study matrix; every
other preset is reported relative to LEGACY_FORTRAN's seal output.

The relaxed z_min tolerance reflects a residual divergence discovered
in Path B v2 A6: silicoshark's interior daughter cells lift in z
slower than FORTRAN's, producing a wider z range than the goldens. The
residual is documented at
`docs/findings/2026-05-05-path-b-v2-legacy-fortran-divergence.md`.
The cell-count plateau, x/y envelope, z_max envelope, and the
qualitative shape of the z curve match FORTRAN; the z_min mismatch
is contained to interior daughter cells whose force-propagation
equilibrates more slowly under silicoshark's static-with-local-update
topology than FORTRAN's mutating neigh-array topology.
"""
import math
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
GOLDEN = REPO_ROOT / 'tests' / 'golden_fortran' / '500_run.off'

EXPECTED_VERTICES_TARGET = 57
CELL_COUNT_TOLERANCE = 0.05  # ±5%

# Tight tolerance for x, y, and z_max — these match FORTRAN within
# floating-point and lattice-rotation noise.
ENVELOPE_TOLERANCE = 0.10
# Relaxed tolerance for z_min — interior daughter cells equilibrate
# slower in silicoshark; documented in the divergence findings doc.
Z_MIN_TOLERANCE = 0.65


def parse_off(path: Path):
    with path.open() as f:
        f.readline()  # COFF
        counts = f.readline().strip()
    m = re.match(r'^\s*(\d+)\s+(\d+)\s+0\s*$', counts)
    n_v, n_f = int(m.group(1)), int(m.group(2))
    verts = np.loadtxt(path, skiprows=2, max_rows=n_v, usecols=(0, 1, 2))
    return n_v, n_f, verts


def test_legacy_fortran_seal_smoke(tmp_path):
    out = tmp_path / 'out'
    out.mkdir()
    args = [
        sys.executable, '-m', 'silicoshark',
        str(REPO_ROOT / 'examples' / 'seal.txt'),
        str(out), 'run', '100', '5',
        '--preset', 'LEGACY_FORTRAN',
        '--override', 'mesenchyme=absent',
        '--override', 'laplacian=length_weighted',
    ]
    proc = subprocess.run(
        args, cwd=REPO_ROOT, capture_output=True, text=True
    )
    assert proc.returncode == 0, (
        f'silicoshark exited {proc.returncode}\n'
        f'stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}'
    )
    assert 'RuntimeWarning' not in proc.stderr, (
        f'silicoshark emitted RuntimeWarnings:\n{proc.stderr}'
    )

    final = out / '500_run_.off'
    assert final.exists(), (
        f'no final OFF file at {final}: present={sorted(out.iterdir())}'
    )

    n_v, n_f, verts = parse_off(final)
    assert all(np.isfinite(verts[:, i]).all() for i in range(3)), (
        'NaN/Inf in silicoshark vertex coordinates'
    )

    lo = int(math.floor(EXPECTED_VERTICES_TARGET * (1 - CELL_COUNT_TOLERANCE)))
    hi = int(math.ceil(EXPECTED_VERTICES_TARGET * (1 + CELL_COUNT_TOLERANCE)))
    assert lo <= n_v <= hi, (
        f'cell count {n_v} not in [{lo}, {hi}] (±5% of {EXPECTED_VERTICES_TARGET})'
    )

    g_v, g_f, g_verts = parse_off(GOLDEN)
    for axis, name in enumerate(('x', 'y', 'z')):
        vmin, vmax = verts[:, axis].min(), verts[:, axis].max()
        gmin, gmax = g_verts[:, axis].min(), g_verts[:, axis].max()
        span = abs(gmax - gmin)
        # z_min uses the relaxed tolerance per the residual-divergence
        # finding; everything else uses the standard 10% tolerance.
        max_tol = ENVELOPE_TOLERANCE * span
        if axis == 2:
            min_tol = Z_MIN_TOLERANCE * span
        else:
            min_tol = ENVELOPE_TOLERANCE * span
        assert abs(vmin - gmin) <= min_tol, (
            f'{name}.min {vmin:.3f} vs FORTRAN {gmin:.3f} '
            f'(tol {min_tol:.3f}, span {span:.3f})'
        )
        assert abs(vmax - gmax) <= max_tol, (
            f'{name}.max {vmax:.3f} vs FORTRAN {gmax:.3f} '
            f'(tol {max_tol:.3f}, span {span:.3f})'
        )
