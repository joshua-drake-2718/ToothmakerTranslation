"""Regression smoke for `examples/wt-tribosphenic-2014.txt`.

The 2014 wild-type tribosphenic mouse parameter set (Harjunmaa et al.,
*Nature*, Ext Data Fig. 2) is the canonical published parameter set
that exercises the activator–inhibitor coupling in the silicoshark
code path. Unlike `examples/seal.txt` — whose `ina = 0.15` initial Act
never crosses the knot threshold under the silicoshark
`mesenchyme='absent'` workaround, leaving Inh permanently at zero —
this set seeds Act above threshold and produces non-trivial knot
formation within the first save block.

This is a regression test for *cusp formation*, not byte-equality.
The headline assertion is that LEGACY_FORTRAN forms at least one knot
cell within 100 iterations on the wt-tribosphenic set. If that
property breaks, the eq.14 typo question (PAPER_2010 vs
PAPER_LITERAL_2010) becomes untestable on this parameter file, and
the comparison-study results in
`experiments/discretisation-study/cusp-forming/` would no longer
reflect the differentiable biology the parameter file is meant to
exercise.

Knot cells are identified by the cyan RGBA `(0, 1, 1, 0)` colour
written by `silicoshark.io.cell_colour`.
"""
import re
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
PARAMS = REPO_ROOT / 'examples' / 'wt-tribosphenic-2014.txt'

# Match the runner's tolerance.
KNOT_RGBA = (0.0, 1.0, 1.0, 0.0)
KNOT_TOL = 1e-9


def _parse_off(path: Path) -> tuple[int, np.ndarray, np.ndarray]:
    with path.open() as f:
        f.readline()  # COFF
        header = f.readline().strip()
    m = re.match(r'^\s*(\d+)\s+(\d+)\s+0\s*$', header)
    assert m is not None, f'malformed OFF header: {header!r}'
    n_v = int(m.group(1))
    if n_v == 0:
        return 0, np.zeros((0, 3)), np.zeros((0, 4))
    arr = np.loadtxt(path, skiprows=2, max_rows=n_v, usecols=(0, 1, 2, 3, 4, 5, 6))
    if arr.ndim == 1:
        arr = arr[None, :]
    return n_v, arr[:, :3], arr[:, 3:7]


def _knot_count(colours: np.ndarray) -> int:
    if colours.shape[0] == 0:
        return 0
    target = np.asarray(KNOT_RGBA, dtype=np.float64)
    mask = np.all(np.abs(colours - target) <= KNOT_TOL, axis=1)
    return int(mask.sum())


def test_wt_tribosphenic_2014_forms_knots_under_legacy_fortran(tmp_path):
    """Run wt-tribosphenic-2014 for 100 iters under LEGACY_FORTRAN; assert
    at least one knot cell forms.

    The cusp-forming property is what makes this parameter file useful
    for the comparison study. If a future code change suppresses knot
    formation here, the eq.14 typo's biological footprint becomes
    invisible to the study.
    """
    out = tmp_path / 'out'
    out.mkdir()
    args = [
        sys.executable, '-m', 'silicoshark',
        str(PARAMS),
        str(out), 'run', '100', '1',
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

    final = out / '100_run_.off'
    assert final.exists(), (
        f'no final OFF file at {final}: present={sorted(out.iterdir())}'
    )
    n_v, positions, colours = _parse_off(final)
    assert n_v >= 1
    knots = _knot_count(colours)
    assert knots >= 1, (
        f'wt-tribosphenic-2014 produced {knots} knot cells under '
        'LEGACY_FORTRAN at iter 100; expected >= 1. The cusp-forming '
        'property the comparison study depends on has regressed.'
    )
