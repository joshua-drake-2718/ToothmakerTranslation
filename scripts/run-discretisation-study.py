#!/usr/bin/env python3
"""Run the Path B v2 comparison study across all named Discretisation
presets and one or more parameter files.

For each (preset, params_file) combination the runner:

1. Resolves the effective Discretisation (preset + study-wide
   overrides + auto-overrides for the not-yet-implemented `mesenchyme`
   / `laplacian='fortran_margins'` fields). The auto-overrides are
   conditional: a preset that already specifies `mesenchyme='absent'`
   and a non-`fortran_margins` Laplacian gets no auto-overrides; the
   currently-shipping presets all need at least one of them.
2. Invokes silicoshark via subprocess (the simulator currently lives
   behind a CLI; a thin in-process wrapper would have to duplicate
   the CLI's preset-and-override plumbing). Captures stdout/stderr to
   the run directory.
3. Reads the run's OFF files back to reconstruct (a) the cell-count
   history (one int per save block), and (b) the final state's
   positions plus knot inference (knot cells are tagged with cyan
   RGBA `(0, 1, 1, 0)` in `silicoshark.io.cell_colour`).
4. Computes the metrics defined in `silicoshark.metrics` and writes
   them to `<run_dir>/metrics.json`. The OFF files, parameter-file
   snapshot, and preset snapshot also live in the run directory.

After all runs complete, the runner writes:

- `<out>/results.json`: nested {preset → params_basename → metrics
  + regime + cell_count_history}.
- `<out>/results-table.md`: one markdown table per parameter file,
  rows = presets, columns = the headline metrics.

ProgressReporter integration (per CLAUDE.md §Batch runner
instrumentation): the runner emits a batch-level progress.json with
experiment name `batch_runner_discretisation_study`. exptop picks
this up as the WORK line; each silicoshark subprocess writes its
own progress.json (TASK line). On non-zero subprocess exit the
runner overwrites the per-run progress.json status to `'failed'` so
exptop doesn't show a dead run as still 'running' forever (the
subprocess died without ProgressReporter's __exit__ firing).
"""
from __future__ import annotations

import argparse
import dataclasses
import json
import os
import re
import shutil
import subprocess
import sys
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))

from silicoshark.discretisation import ALL_PRESETS, Discretisation
from silicoshark.metrics import (
    cell_count_plateau,
    cusp_count,
    cusp_width,
    regime,
    vertex_envelope,
)
from silicoshark.utils import ProgressReporter


# Knot cells are written with the cyan RGBA `(0, 1, 1, 0)` colour by
# `silicoshark.io.cell_colour`. The OFF reader compares to this with
# fp tolerance (the writer emits the literal `0` and `1` so equality
# is bit-exact, but we keep a small tolerance for safety).
KNOT_RGBA = (0.0, 1.0, 1.0, 0.0)
KNOT_TOL = 1e-9


# Discretisation fields whose options are not implemented in the
# current code path; if a preset uses one of these, the runner adds an
# `--override` to fall back to the implemented option. This mirrors
# `silicoshark.__main__._check_implemented`.
DEFERRED_FIELD_OVERRIDES: dict[str, tuple[str, str]] = {
    # field name -> (forbidden_value, replacement_value)
    'mesenchyme': ('per_column_z_layers', 'absent'),
    'laplacian': ('fortran_margins', 'length_weighted'),
}


# ------------------------------------------------------------ helpers


def _coerce(value: str):
    """Coerce a string to int -> float -> bool -> str. Mirrors the CLI."""
    low = value.lower()
    if low in ('true', 'false'):
        return low == 'true'
    if low == 'none':
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def _parse_overrides(overrides: list[str]) -> dict[str, object]:
    """Parse `key=value` strings into a dict of coerced values."""
    out: dict[str, object] = {}
    for ov in overrides:
        if '=' not in ov:
            raise SystemExit(f'malformed --override (expected key=value): {ov!r}')
        k, _, v = ov.partition('=')
        out[k.strip()] = _coerce(v.strip())
    return out


def _auto_overrides_for(disc: Discretisation) -> dict[str, str]:
    """Return any `--override` entries needed to make `disc` runnable
    under the current silicoshark code path.

    Currently: `mesenchyme=absent` if the preset has
    `mesenchyme='per_column_z_layers'`; `laplacian=length_weighted`
    if the preset has `laplacian='fortran_margins'`.
    """
    out: dict[str, str] = {}
    for field_name, (forbidden, replacement) in DEFERRED_FIELD_OVERRIDES.items():
        if getattr(disc, field_name) == forbidden:
            out[field_name] = replacement
    return out


def _effective_disc(
    preset: Discretisation,
    user_overrides: dict[str, object],
) -> Discretisation:
    """Apply user overrides (and only those) to `preset`. Used for the
    snapshot saved alongside the run; the subprocess receives the
    same overrides via its own CLI parsing.
    """
    if not user_overrides:
        return preset
    valid = {f.name for f in dataclasses.fields(preset)}
    bad = set(user_overrides) - valid
    if bad:
        raise SystemExit(f'unknown Discretisation fields: {sorted(bad)}')
    return dataclasses.replace(preset, **user_overrides)


# ------------------------------------------------------------ OFF reader


_OFF_HEADER_RE = re.compile(r'^\s*(\d+)\s+(\d+)\s+0\s*$')


def parse_off(path: Path) -> tuple[int, np.ndarray, np.ndarray]:
    """Read a COFF file. Returns (n_v, positions, colours).

    `positions` is (n_v, 3) float64; `colours` is (n_v, 4) float64
    (RGBA, in the writer's order). Faces are skipped — the metrics
    only need vertex data.
    """
    with path.open() as f:
        f.readline()  # 'COFF'
        header = f.readline().strip()
    m = _OFF_HEADER_RE.match(header)
    if m is None:
        raise ValueError(f'malformed OFF header in {path}: {header!r}')
    n_v = int(m.group(1))
    if n_v == 0:
        return 0, np.zeros((0, 3)), np.zeros((0, 4))
    arr = np.loadtxt(path, skiprows=2, max_rows=n_v, usecols=(0, 1, 2, 3, 4, 5, 6))
    if arr.ndim == 1:
        arr = arr[None, :]
    positions = arr[:, :3].astype(np.float64)
    colours = arr[:, 3:7].astype(np.float64)
    return n_v, positions, colours


def _knot_mask_from_colours(colours: np.ndarray) -> np.ndarray:
    """Boolean mask over cells: True where the colour matches knot RGBA."""
    if colours.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    target = np.asarray(KNOT_RGBA, dtype=np.float64)
    return np.all(np.abs(colours - target) <= KNOT_TOL, axis=1)


# ------------------------------------------------------------ subprocess runner


def _build_subprocess_args(
    params_file: Path,
    run_dir: Path,
    out_name: str,
    iters: int,
    saves: int,
    preset_name: str,
    overrides: dict[str, object],
) -> list[str]:
    """Construct the silicoshark CLI argv for one run."""
    args = [
        sys.executable, '-m', 'silicoshark',
        str(params_file),
        str(run_dir),
        out_name,
        str(iters),
        str(saves),
        '--preset', preset_name,
    ]
    for k, v in overrides.items():
        # The CLI's _coerce reverses int/float/bool/str; we must emit
        # a string the CLI can re-coerce back to the same value.
        args.extend(['--override', f'{k}={v}'])
    return args


def _mark_progress_failed(progress_path: Path, reason: str) -> None:
    """Overwrite a per-run progress.json's status to 'failed'.

    Required by CLAUDE.md §Stale progress cleanup: when a subprocess
    dies (OOM, exception during init, kill -9) the ProgressReporter
    context manager's __exit__ does not fire, leaving status='running'
    permanently. exptop then displays dead runs as live. The batch
    runner detects non-zero exit codes and marks the file failed.
    """
    if not progress_path.exists():
        return
    try:
        payload = json.loads(progress_path.read_text())
    except (OSError, json.JSONDecodeError):
        return
    payload['status'] = 'failed'
    payload['message'] = reason
    payload['updated_at'] = datetime.now(timezone.utc).isoformat()
    progress_path.write_text(json.dumps(payload, indent=2) + '\n')


# ------------------------------------------------------------ per-run analysis


def _final_off_path(run_dir: Path, iters: int, saves: int, out_name: str) -> Path:
    """Final-block OFF file; matches __main__.py's naming convention."""
    iter_label = str(iters * saves)
    return run_dir / f'{iter_label}_{out_name}_.off'


def _cell_count_history(run_dir: Path, iters: int, saves: int, out_name: str) -> list[int]:
    """Read each save block's OFF and return its cell count.

    Returns a list of length `saves` (one entry per save block); a
    block whose OFF file is missing contributes a `np.nan` entry,
    which the regime classifier picks up as a divergence.
    """
    history: list[int] = []
    for block in range(1, saves + 1):
        iter_label = str(block * iters)
        path = run_dir / f'{iter_label}_{out_name}_.off'
        if not path.exists():
            history.append(int(-1))  # placeholder; converted to nan in regime
            continue
        try:
            n_v, _, _ = parse_off(path)
        except (ValueError, OSError):
            history.append(int(-1))
            continue
        history.append(n_v)
    return history


def _analyse_run(
    run_dir: Path,
    iters: int,
    saves: int,
    out_name: str,
) -> dict:
    """Compute metrics for a completed run from its on-disk artefacts.

    Returns a dict of metrics. `cell_count_history` is included; the
    aggregator decides whether to record it in the master results.
    """
    raw_history = _cell_count_history(run_dir, iters, saves, out_name)
    # Convert -1 placeholder to nan for regime classification (regime
    # treats any non-finite as 'NaN').
    history_for_regime: list[float] = [
        float('nan') if h < 0 else float(h) for h in raw_history
    ]

    final = _final_off_path(run_dir, iters, saves, out_name)
    if not final.exists():
        return {
            'cell_count_history': raw_history,
            'cell_count_plateau': float('nan'),
            'cell_count_final': None,
            'cusp_count': None,
            'cusp_width': None,
            'vertex_envelope': None,
            'regime': 'NaN',
            'note': 'no final OFF file',
        }

    n_v, positions, colours = parse_off(final)
    knot = _knot_mask_from_colours(colours)
    finite_pos = np.isfinite(positions).all() if positions.size else False

    plateau = cell_count_plateau([h for h in raw_history if h >= 0])
    if any(h < 0 for h in raw_history) or not finite_pos:
        run_regime = 'NaN'
    else:
        run_regime = regime(history_for_regime)

    # cusp_width needs a Mesh built from the final positions. We
    # construct it directly here (the run's State is no longer in
    # memory; positions + knot are sufficient to call the metric).
    if knot.sum() < 2:
        width = 0.0
    elif finite_pos and positions.shape[0] >= 3:
        from silicoshark.mesh import Mesh
        from silicoshark.state import State

        n = positions.shape[0]
        zeros_n = np.zeros(n, dtype=np.float64)
        zeros_b = np.zeros(n, dtype=bool)
        zeros_mes = np.zeros((n, 2), dtype=np.float64)
        state = State(
            positions=positions,
            act=zeros_n.copy(),
            inh=zeros_n.copy(),
            sec=zeros_n.copy(),
            diff=zeros_n.copy(),
            knot=knot,
            mes_inh=zeros_mes.copy(),
            mes_sec=zeros_mes.copy(),
            init_anterior=zeros_b.copy(),
            init_posterior=zeros_b.copy(),
            init_lingual=zeros_b.copy(),
            init_buccal=zeros_b.copy(),
        )
        try:
            mesh = Mesh.from_positions(positions)
            width = cusp_width(state, mesh)
        except Exception:  # pragma: no cover — degenerate triangulations
            width = float('nan')
    else:
        width = float('nan')

    envelope = vertex_envelope(positions) if finite_pos else None

    return {
        'cell_count_history': raw_history,
        'cell_count_plateau': plateau,
        'cell_count_final': int(n_v),
        'cusp_count': int(knot.sum()),
        'cusp_width': float(width),
        'vertex_envelope': envelope,
        'regime': run_regime,
    }


# ------------------------------------------------------------ aggregation


def _format_envelope(env: dict[str, float] | None) -> str:
    if env is None:
        return '—'
    return (
        f'x [{env["x_min"]:.2f}, {env["x_max"]:.2f}] · '
        f'y [{env["y_min"]:.2f}, {env["y_max"]:.2f}] · '
        f'z [{env["z_min"]:.2f}, {env["z_max"]:.2f}]'
    )


def _format_float(x: float | None, fmt: str = '.2f') -> str:
    if x is None:
        return '—'
    try:
        if not np.isfinite(x):
            return 'NaN'
    except TypeError:
        return str(x)
    return format(x, fmt)


def _write_results_table(
    out_dir: Path,
    results: dict[str, dict[str, dict]],
    presets: list[str],
    params_basenames: list[str],
) -> None:
    """Write one markdown table per parameter file.

    Columns: preset, plateau, final, cusps, cusp width, envelope,
    regime. Rows: presets, in CLI order.
    """
    lines: list[str] = []
    lines.append('# Path B v2 discretisation comparison study\n')
    lines.append('Generated: '
                 + datetime.now(timezone.utc).isoformat() + '\n')
    for params_basename in params_basenames:
        lines.append(f'\n## `{params_basename}`\n')
        lines.append('| Preset | Plateau | Final | Cusps | Cusp width | Envelope | Regime |')
        lines.append('|---|---:|---:|---:|---:|---|---|')
        for preset_name in presets:
            entry = results.get(preset_name, {}).get(params_basename)
            if entry is None:
                lines.append(f'| `{preset_name}` | — | — | — | — | — | (no run) |')
                continue
            plateau = _format_float(entry.get('cell_count_plateau'))
            final = entry.get('cell_count_final')
            final_s = '—' if final is None else str(final)
            cusps = entry.get('cusp_count')
            cusps_s = '—' if cusps is None else str(cusps)
            width = _format_float(entry.get('cusp_width'), '.3f')
            env_s = _format_envelope(entry.get('vertex_envelope'))
            run_regime = entry.get('regime', '—')
            lines.append(
                f'| `{preset_name}` | {plateau} | {final_s} | {cusps_s} | '
                f'{width} | {env_s} | {run_regime} |'
            )
    (out_dir / 'results-table.md').write_text('\n'.join(lines) + '\n')


# ------------------------------------------------------------ main


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog='run-discretisation-study',
        description=(
            'Run the Path B v2 comparison study across Discretisation '
            'presets and parameter files. Writes results.json and '
            'results-table.md to --out.'
        ),
    )
    parser.add_argument(
        '--params', nargs='+', default=[str(REPO_ROOT / 'examples' / 'seal.txt')],
        help='Parameter files (default: examples/seal.txt).',
    )
    parser.add_argument(
        '--presets', nargs='+', default=sorted(ALL_PRESETS.keys()),
        help='Presets to include (default: every entry of ALL_PRESETS).',
    )
    parser.add_argument(
        '--iters', type=int, default=500,
        help='Iterations per save block (default: 500).',
    )
    parser.add_argument(
        '--saves', type=int, default=5,
        help='Number of save blocks (default: 5).',
    )
    parser.add_argument(
        '--out', default=str(REPO_ROOT / 'experiments' / 'discretisation-study' / 'results'),
        help='Output directory.',
    )
    parser.add_argument(
        '--override', action='append', default=[],
        metavar='KEY=VALUE',
        help='Override a Discretisation field on every preset. '
             'May be repeated. Auto-overrides for the deferred '
             'mesenchyme / fortran_margins fields are applied '
             'automatically and need not be passed here.',
    )
    args = parser.parse_args(argv)

    params_files = [Path(p).resolve() for p in args.params]
    for p in params_files:
        if not p.exists():
            raise SystemExit(f'params file not found: {p}')
    presets = list(args.presets)
    for name in presets:
        if name not in ALL_PRESETS:
            raise SystemExit(
                f'unknown preset: {name!r}. Valid: {sorted(ALL_PRESETS)}'
            )

    user_overrides = _parse_overrides(args.override)

    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Batch-level ProgressReporter (WORK line for exptop). Total phases
    # = one per (preset × params_file) combination.
    total_phases = len(presets) * len(params_files)
    batch_progress_path = out_dir / 'progress.json'

    results: dict[str, dict[str, dict]] = {}
    params_basenames = [p.stem for p in params_files]

    # `flush=True` per CLAUDE.md §Flushed output under nohup so that
    # status messages between subprocess launches reach the log
    # synchronously, not when the buffer fills.
    print(f'Running discretisation study: {len(presets)} presets × '
          f'{len(params_files)} params = {total_phases} runs',
          flush=True)
    print(f'Output: {out_dir}', flush=True)

    with ProgressReporter(
        batch_progress_path,
        experiment='batch_runner_discretisation_study',
        total_phases=total_phases,
    ) as batch_progress:
        for preset_name in presets:
            results.setdefault(preset_name, {})
            preset_disc = ALL_PRESETS[preset_name]

            for params_file in params_files:
                params_basename = params_file.stem
                batch_progress.begin_phase(
                    name=f'{preset_name}/{params_basename}',
                    description='running silicoshark',
                    total_steps=args.iters,
                )

                run_dir = out_dir / preset_name / params_basename
                # Wipe stale OFF files so we don't read previous-run
                # artefacts if the new run terminates early.
                if run_dir.exists():
                    shutil.rmtree(run_dir)
                run_dir.mkdir(parents=True, exist_ok=True)

                # Effective overrides = user overrides + auto-overrides
                # for any deferred field options the preset uses.
                auto_ov = _auto_overrides_for(preset_disc)
                effective_ov: dict[str, object] = dict(auto_ov)
                effective_ov.update(user_overrides)
                # Snapshot the resolved Discretisation alongside.
                effective_disc = _effective_disc(preset_disc, effective_ov)

                # Snapshots for reproducibility.
                shutil.copy(params_file, run_dir / 'params-snapshot.txt')
                (run_dir / 'preset-snapshot.json').write_text(
                    json.dumps({
                        'preset_name': preset_name,
                        'overrides': {k: str(v) for k, v in effective_ov.items()},
                        'effective': dataclasses.asdict(effective_disc),
                    }, indent=2) + '\n'
                )

                out_name = 'run'
                cmd = _build_subprocess_args(
                    params_file=params_file,
                    run_dir=run_dir,
                    out_name=out_name,
                    iters=args.iters,
                    saves=args.saves,
                    preset_name=preset_name,
                    overrides=effective_ov,
                )

                stdout_path = run_dir / 'silicoshark.stdout'
                stderr_path = run_dir / 'silicoshark.stderr'
                # Fresh env: ensure REPO_ROOT is on PYTHONPATH for the
                # subprocess too (matters when this script is launched
                # from outside the repo via an absolute path).
                env = os.environ.copy()
                pp = env.get('PYTHONPATH', '')
                env['PYTHONPATH'] = (
                    f'{REPO_ROOT}{os.pathsep}{pp}' if pp else str(REPO_ROOT)
                )

                print(f'  [{preset_name}/{params_basename}] launching '
                      f'{args.iters}×{args.saves} iters',
                      flush=True)
                with stdout_path.open('w') as so, stderr_path.open('w') as se:
                    proc = subprocess.run(
                        cmd, cwd=REPO_ROOT, stdout=so, stderr=se,
                        env=env, text=True,
                    )

                # Per-run progress.json status correction (CLAUDE.md
                # §Stale progress cleanup). On non-zero exit we mark the
                # per-run progress.json as failed; on zero exit the CLI's
                # ProgressReporter context manager has already marked it
                # 'completed' on its way out.
                run_progress = run_dir / 'progress.json'
                if proc.returncode != 0:
                    err_tail = stderr_path.read_text().strip().splitlines()[-5:]
                    msg = (
                        f'silicoshark exited {proc.returncode}; '
                        f'last stderr: {" | ".join(err_tail)}'
                    )
                    _mark_progress_failed(run_progress, msg)
                    print(f'    FAILED: {msg}', flush=True)

                # Analyse whatever the run produced (even on failure
                # the partial OFFs are useful).
                try:
                    metrics = _analyse_run(
                        run_dir, args.iters, args.saves, out_name
                    )
                except Exception as exc:  # noqa: BLE001 — runner shouldn't crash
                    metrics = {
                        'cell_count_history': [],
                        'regime': 'NaN',
                        'note': f'analysis failed: {exc!r}',
                        'traceback': traceback.format_exc(),
                    }

                metrics['preset'] = preset_name
                metrics['params_file'] = str(params_file)
                metrics['returncode'] = proc.returncode
                (run_dir / 'metrics.json').write_text(
                    json.dumps(metrics, indent=2) + '\n'
                )

                results[preset_name][params_basename] = metrics

                cells = metrics.get('cell_count_final')
                cusps = metrics.get('cusp_count')
                run_regime = metrics.get('regime')
                print(f'    done: cells={cells} cusps={cusps} regime={run_regime}',
                      flush=True)

                # Mark phase step at end.
                batch_progress.step(args.iters, 'finished')

    # Aggregate.
    (out_dir / 'results.json').write_text(
        json.dumps(results, indent=2) + '\n'
    )
    _write_results_table(out_dir, results, presets, params_basenames)
    print(f'Wrote {out_dir / "results.json"} and {out_dir / "results-table.md"}',
          flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
