"""End-to-end test for `scripts/run-discretisation-study.py
--single-field-mode`.

The single-field-mode flag is the Path B v2 B1 disentanglement entry
point. This test exercises it on a deliberately tiny configuration
(50 iters x 5 saves, knock-down direction only) and asserts the three
contracts the report needs:

1. The output `results-table.md` contains the LEGACY_FORTRAN anchor
   row (the knock-down baseline).
2. At least one perturbation row appears in the table (the anchor on
   its own would be useless; the table only earns its place if it
   shows the perturbation matrix).
3. `results.json` is valid JSON and keys at least one row whose entry
   carries the headline metric fields the table renders from.

50 iters x 5 saves on the seal example takes ~2-3 seconds per run; the
full knock-down sweep (1 anchor + 13 perturbations = 14 runs) fits in
under a minute on a developer machine.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


def test_single_field_mode_knock_down_smoke(tmp_path):
    out_dir = tmp_path / 'single-field-knock-down'
    cmd = [
        sys.executable,
        str(REPO_ROOT / 'scripts' / 'run-discretisation-study.py'),
        '--single-field-mode', 'knock-down',
        '--iters', '50',
        '--saves', '5',
        '--out', str(out_dir),
        '--override', 'mesenchyme=absent',
        '--override', 'laplacian=length_weighted',
        # Tight per-run timeout: B1 found that knocking out
        # division_total_cap, knot_threshold_gate, and rep_form
        # produces runaway division. The runner records each as
        # regime=NaN and continues; without the timeout this test
        # would hang for hours.
        '--per-run-timeout', '20',
    ]
    proc = subprocess.run(
        cmd, cwd=REPO_ROOT, capture_output=True, text=True,
        timeout=600,
    )
    assert proc.returncode == 0, (
        f'runner exited {proc.returncode}\n'
        f'stdout:\n{proc.stdout}\nstderr:\n{proc.stderr}'
    )

    # (a) The output table contains the anchor row.
    table_path = out_dir / 'results-table.md'
    assert table_path.exists(), (
        f'no results-table.md at {table_path}; '
        f'out_dir contents: {sorted(p.name for p in out_dir.iterdir())}'
    )
    table_text = table_path.read_text()
    assert '`LEGACY_FORTRAN` (anchor)' in table_text, (
        f'LEGACY_FORTRAN anchor row missing from results-table.md:\n'
        f'{table_text}'
    )

    # (b) At least one perturbation row.
    assert 'LEGACY_FORTRAN_minus_' in table_text, (
        f'no perturbation rows in results-table.md:\n{table_text}'
    )

    # (c) results.json valid + has at least one row with headline metrics.
    results_path = out_dir / 'results.json'
    assert results_path.exists(), (
        f'no results.json at {results_path}'
    )
    payload = json.loads(results_path.read_text())  # raises on bad JSON
    assert isinstance(payload, dict) and payload, (
        f'results.json is not a non-empty dict: {payload!r}'
    )
    # The anchor row must be in there with at least one params_basename.
    assert 'LEGACY_FORTRAN' in payload, (
        f'results.json missing LEGACY_FORTRAN anchor: keys={list(payload)}'
    )
    anchor_runs = payload['LEGACY_FORTRAN']
    assert anchor_runs, (
        f'LEGACY_FORTRAN entry is empty: {anchor_runs!r}'
    )
    # Pick any params_basename and check the metric fields exist.
    any_run = next(iter(anchor_runs.values()))
    for field in ('cell_count_plateau', 'cusp_count', 'regime',
                  'vertex_envelope', 'cell_count_history'):
        assert field in any_run, (
            f'metric {field!r} missing from results.json entry: '
            f'keys={list(any_run)}'
        )
