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

from silicoshark.discretisation import (
    ALL_PRESETS,
    LEGACY_FORTRAN,
    PATH_B_DEFAULT,
    Discretisation,
)
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


# A single row in the results table. In the default `--presets` mode each
# preset becomes one row with `base_preset == row_label`; in single-field
# mode each row is a perturbation of one of the two anchor presets,
# expressed as `base_preset` plus a one-element `perturbation` dict.
@dataclasses.dataclass(frozen=True)
class RowSpec:
    """Recipe for one row in the results table.

    `row_label` is the table row and run-directory name. `base_preset`
    is the named preset passed to silicoshark via `--preset`.
    `perturbation` is the additional `{field: value}` patch on top of the
    base preset, applied via `--override` in addition to user overrides
    and auto-overrides.
    """
    row_label: str
    base_preset: str
    perturbation: dict[str, object]


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


# ------------------------------------------------------------ single-field mode


def _differing_fields(a: Discretisation, b: Discretisation) -> list[str]:
    """Return the names of Discretisation fields where `a` and `b`
    disagree. Drives the single-field-mode perturbation matrix so the
    suite picks up new fields automatically as they are added.
    """
    out: list[str] = []
    for f in dataclasses.fields(a):
        if getattr(a, f.name) != getattr(b, f.name):
            out.append(f.name)
    return out


def _build_single_field_rows(mode: str) -> list[RowSpec]:
    """Build the row specs for single-field-mode.

    The two anchor presets (`LEGACY_FORTRAN` and `PATH_B_DEFAULT`) are
    each emitted once; for each of the fields where they disagree, a
    knock-down row (LEGACY_FORTRAN with that field replaced by
    PATH_B_DEFAULT's value) and/or a knock-up row (the symmetric one)
    is generated, depending on `mode`.

    The returned list is ordered: anchors first in their canonical
    pairing direction, then the perturbations (knock-down then
    knock-up), with field names sorted so the table is reproducible.
    """
    if mode not in ('knock-down', 'knock-up', 'both'):
        raise SystemExit(
            f'unknown --single-field-mode: {mode!r}. '
            "Valid: 'knock-down', 'knock-up', 'both'."
        )
    diff_fields = sorted(_differing_fields(LEGACY_FORTRAN, PATH_B_DEFAULT))
    if not diff_fields:
        raise SystemExit(
            'LEGACY_FORTRAN and PATH_B_DEFAULT carry identical field '
            'values; nothing to disentangle.'
        )

    rows: list[RowSpec] = []
    seen: set[str] = set()

    def _push(row: RowSpec) -> None:
        if row.row_label in seen:
            return
        seen.add(row.row_label)
        rows.append(row)

    # Anchor for the knock-down half (LEGACY_FORTRAN baseline).
    if mode in ('knock-down', 'both'):
        _push(RowSpec(
            row_label='LEGACY_FORTRAN',
            base_preset='LEGACY_FORTRAN',
            perturbation={},
        ))
        for fname in diff_fields:
            _push(RowSpec(
                row_label=f'LEGACY_FORTRAN_minus_{fname}',
                base_preset='LEGACY_FORTRAN',
                perturbation={fname: getattr(PATH_B_DEFAULT, fname)},
            ))

    # Anchor for the knock-up half (PATH_B_DEFAULT baseline).
    if mode in ('knock-up', 'both'):
        _push(RowSpec(
            row_label='PATH_B_DEFAULT',
            base_preset='PATH_B_DEFAULT',
            perturbation={},
        ))
        for fname in diff_fields:
            _push(RowSpec(
                row_label=f'PATH_B_DEFAULT_plus_{fname}',
                base_preset='PATH_B_DEFAULT',
                perturbation={fname: getattr(LEGACY_FORTRAN, fname)},
            ))

    return rows


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
        # The CLI's _coerce reverses int/float/bool/str/None; we must
        # emit a string the CLI can re-coerce back to the same value.
        # Python None must be passed as the literal `null` (JSON
        # convention); the string `'none'` is a valid Discretisation
        # field value (knot_threshold_gate='none') and must NOT be
        # coerced to None on the receiving side.
        rendered = 'null' if v is None else str(v)
        args.extend(['--override', f'{k}={rendered}'])
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


def _format_axis_range(env: dict[str, float] | None, axis: str) -> str:
    """Format one axis of an envelope dict as `[min, max]`.

    `axis` is `'x'`, `'y'`, or `'z'`. Returns `'—'` if env is None.
    """
    if env is None:
        return '—'
    lo = env.get(f'{axis}_min')
    hi = env.get(f'{axis}_max')
    if lo is None or hi is None:
        return '—'
    return f'[{lo:.2f}, {hi:.2f}]'


def _write_single_field_table(
    out_dir: Path,
    rows: list[RowSpec],
    results: dict[str, dict[str, dict]],
    params_basenames: list[str],
    mode: str,
) -> None:
    """Write the single-field-mode results table.

    One section per direction (knock-down, knock-up); within each
    section, one markdown table per parameter file. Columns:
    `Preset, Plateau, Cusps, x range, y range, z range, Regime`.

    The interpretation paragraph below each table calls out the
    plateau delta induced by each perturbation, classifying each
    field as strongly / moderately / not consequential.
    """
    lines: list[str] = []
    lines.append('# Path B v2 single-field disentanglement\n')
    lines.append(
        'Generated: ' + datetime.now(timezone.utc).isoformat() + '\n'
    )
    lines.append(
        '\nFor each field where `LEGACY_FORTRAN` and `PATH_B_DEFAULT` '
        'differ, this study runs a single-field perturbation:\n'
        '\n- *Knock-down*: `LEGACY_FORTRAN` with that field reset to '
        'its `PATH_B_DEFAULT` value. If the plateau collapses toward '
        '`PATH_B_DEFAULT`, the field is doing the work in the '
        'FORTRAN-flavoured bundle.\n'
        '- *Knock-up*: `PATH_B_DEFAULT` with that field set to its '
        '`LEGACY_FORTRAN` value. If the plateau lifts toward '
        '`LEGACY_FORTRAN`, that field alone is sufficient.\n'
    )

    # Partition rows by direction (the row_label prefix encodes the
    # base anchor and direction).
    knock_down = [
        r for r in rows
        if r.base_preset == 'LEGACY_FORTRAN'
    ]
    knock_up = [
        r for r in rows
        if r.base_preset == 'PATH_B_DEFAULT'
    ]

    sections: list[tuple[str, list[RowSpec], str]] = []
    if mode in ('knock-down', 'both') and knock_down:
        sections.append(('Knock-down: LEGACY_FORTRAN minus each field',
                         knock_down, 'LEGACY_FORTRAN'))
    if mode in ('knock-up', 'both') and knock_up:
        sections.append(('Knock-up: PATH_B_DEFAULT plus each field',
                         knock_up, 'PATH_B_DEFAULT'))

    for section_title, section_rows, anchor_label in sections:
        lines.append(f'\n## {section_title}\n')
        for params_basename in params_basenames:
            lines.append(f'\n### `{params_basename}`\n')
            lines.append(
                '| Preset | Plateau | Cusps | x range | y range | '
                'z range | Regime |'
            )
            lines.append(
                '|---|---:|---:|---|---|---|---|'
            )
            anchor_plateau: float | None = None
            for row in section_rows:
                entry = results.get(row.row_label, {}).get(params_basename)
                if entry is None:
                    label = row.row_label
                    suffix = ' (anchor)' if row.row_label == anchor_label else ''
                    lines.append(
                        f'| `{label}`{suffix} | — | — | — | — | — | (no run) |'
                    )
                    continue
                plateau_val = entry.get('cell_count_plateau')
                plateau = _format_float(plateau_val)
                cusps = entry.get('cusp_count')
                cusps_s = '—' if cusps is None else str(cusps)
                env = entry.get('vertex_envelope')
                x_range = _format_axis_range(env, 'x')
                y_range = _format_axis_range(env, 'y')
                z_range = _format_axis_range(env, 'z')
                run_regime = entry.get('regime', '—')

                # Anchor row gets a label suffix; capture its plateau
                # so the interpretation paragraph can compute deltas.
                if row.row_label == anchor_label:
                    label_disp = f'`{row.row_label}` (anchor)'
                    if isinstance(plateau_val, (int, float)) and \
                            plateau_val == plateau_val:  # not nan
                        anchor_plateau = float(plateau_val)
                else:
                    label_disp = f'`{row.row_label}`'

                lines.append(
                    f'| {label_disp} | {plateau} | {cusps_s} | '
                    f'{x_range} | {y_range} | {z_range} | {run_regime} |'
                )

            # Interpretation paragraph.
            lines.append('')
            if anchor_plateau is None:
                lines.append(
                    '_Anchor plateau unavailable; skipping '
                    'interpretation paragraph._\n'
                )
                continue
            # Determine the comparison baseline (the other anchor's
            # plateau) so we know what 'collapse toward' / 'lift toward'
            # means quantitatively.
            other_anchor = (
                'PATH_B_DEFAULT' if anchor_label == 'LEGACY_FORTRAN'
                else 'LEGACY_FORTRAN'
            )
            other_entry = results.get(other_anchor, {}).get(params_basename)
            other_plateau: float | None = None
            if other_entry is not None:
                op = other_entry.get('cell_count_plateau')
                if isinstance(op, (int, float)) and op == op:  # not nan
                    other_plateau = float(op)
            span = (
                abs(anchor_plateau - other_plateau)
                if other_plateau is not None else 0.0
            )

            verdicts: list[str] = []
            for row in section_rows:
                if row.row_label == anchor_label:
                    continue
                # The single perturbed field is the one entry in the
                # perturbation dict. Anchor rows (perturbation={})
                # short-circuited above.
                if not row.perturbation:
                    continue
                fname = next(iter(row.perturbation))
                entry = results.get(row.row_label, {}).get(params_basename)
                if entry is None:
                    continue
                pv = entry.get('cell_count_plateau')
                if not isinstance(pv, (int, float)) or pv != pv:  # nan
                    verdicts.append(
                        f'- Field `{fname}`: run failed or produced '
                        'NaN plateau; effect cannot be assessed.'
                    )
                    continue
                pv = float(pv)
                delta = pv - anchor_plateau
                # 'Strong' = perturbation moves at least 50% of the way
                # toward the other anchor; 'moderate' = 15-50%; 'not'
                # = < 15%. If the two anchors are equal (span=0), every
                # field is reported as 'not consequential'.
                if span > 0:
                    fraction = abs(delta) / span
                else:
                    fraction = 0.0
                if fraction >= 0.5:
                    verdict = 'strongly'
                elif fraction >= 0.15:
                    verdict = 'moderately'
                else:
                    verdict = 'not'
                if anchor_label == 'LEGACY_FORTRAN':
                    direction = 'Removing'
                else:
                    direction = 'Adding'
                verdicts.append(
                    f'- {direction} field `{fname}` '
                    f'{"from" if anchor_label == "LEGACY_FORTRAN" else "to"} '
                    f'`{anchor_label}` shifts the plateau from '
                    f'{anchor_plateau:.2f} to {pv:.2f} '
                    f'(delta {delta:+.2f}; '
                    f'{fraction*100:.0f}% of the {span:.2f}-cell span '
                    f'between anchors) — field `{fname}` is therefore '
                    f'**{verdict} consequential**.'
                )
            if verdicts:
                lines.append('Interpretation:\n')
                lines.extend(verdicts)
                lines.append('')

    (out_dir / 'results-table.md').write_text('\n'.join(lines) + '\n')


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
    parser.add_argument(
        '--single-field-mode',
        choices=('knock-down', 'knock-up', 'both'),
        default=None,
        help='Run a single-field perturbation matrix instead of '
             'iterating over `--presets`. `knock-down` derives one row '
             'per field that differs between LEGACY_FORTRAN and '
             'PATH_B_DEFAULT, replacing that field on LEGACY_FORTRAN '
             "with PATH_B_DEFAULT's value; `knock-up` is the symmetric "
             'replacement on PATH_B_DEFAULT; `both` runs the union. '
             'The two anchor presets are always included as rows.',
    )
    parser.add_argument(
        '--per-run-timeout', type=float, default=300.0,
        help='Maximum wall-clock seconds for any single silicoshark '
             'subprocess. On timeout the run is killed, its progress.json '
             'is marked failed, and the metrics row records regime=NaN. '
             'Default 300s. Set to 0 to disable. The B1 disentanglement '
             'matrix discovered that knocking out stability fields '
             '(division_total_cap, knot_threshold_gate; rep_form on the '
             'knock-up direction) produces runaway division — the timeout '
             'gives the same protection as a manual watchdog.',
    )
    args = parser.parse_args(argv)

    params_files = [Path(p).resolve() for p in args.params]
    for p in params_files:
        if not p.exists():
            raise SystemExit(f'params file not found: {p}')

    user_overrides = _parse_overrides(args.override)

    out_dir = Path(args.out).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build the row list. In default mode each `--presets` entry
    # becomes one row whose base preset is itself with no perturbation;
    # in single-field mode the anchors plus per-field perturbations are
    # generated programmatically.
    if args.single_field_mode:
        rows = _build_single_field_rows(args.single_field_mode)
        mode_desc = f'single-field-mode={args.single_field_mode}'
    else:
        presets = list(args.presets)
        for name in presets:
            if name not in ALL_PRESETS:
                raise SystemExit(
                    f'unknown preset: {name!r}. Valid: {sorted(ALL_PRESETS)}'
                )
        rows = [
            RowSpec(row_label=p, base_preset=p, perturbation={})
            for p in presets
        ]
        mode_desc = f'{len(rows)} presets'

    # Batch-level ProgressReporter (WORK line for exptop). Total phases
    # = one per (row × params_file) combination.
    total_phases = len(rows) * len(params_files)
    batch_progress_path = out_dir / 'progress.json'

    results: dict[str, dict[str, dict]] = {}
    params_basenames = [p.stem for p in params_files]

    # `flush=True` per CLAUDE.md §Flushed output under nohup so that
    # status messages between subprocess launches reach the log
    # synchronously, not when the buffer fills.
    print(f'Running discretisation study: {mode_desc} × '
          f'{len(params_files)} params = {total_phases} runs',
          flush=True)
    print(f'Output: {out_dir}', flush=True)

    with ProgressReporter(
        batch_progress_path,
        experiment='batch_runner_discretisation_study',
        total_phases=total_phases,
    ) as batch_progress:
        for row in rows:
            results.setdefault(row.row_label, {})
            base_disc = ALL_PRESETS[row.base_preset]

            for params_file in params_files:
                params_basename = params_file.stem
                batch_progress.begin_phase(
                    name=f'{row.row_label}/{params_basename}',
                    description='running silicoshark',
                    total_steps=args.iters,
                )

                run_dir = out_dir / row.row_label / params_basename
                # Wipe stale OFF files so we don't read previous-run
                # artefacts if the new run terminates early.
                if run_dir.exists():
                    shutil.rmtree(run_dir)
                run_dir.mkdir(parents=True, exist_ok=True)

                # Effective overrides = auto-overrides + the row's
                # single-field perturbation + user overrides. User
                # overrides take precedence over the perturbation,
                # which takes precedence over auto-overrides.
                auto_ov = _auto_overrides_for(base_disc)
                effective_ov: dict[str, object] = dict(auto_ov)
                effective_ov.update(row.perturbation)
                effective_ov.update(user_overrides)
                # If the perturbation moves a field into a deferred
                # value (e.g. mesenchyme=per_column_z_layers, which
                # is unimplemented), we have to add the auto-override
                # for the perturbed field too. Recompute auto-overrides
                # against the post-perturbation Discretisation for safety.
                post_perturb = _effective_disc(base_disc, dict(row.perturbation))
                for fname, (forbidden, replacement) in \
                        DEFERRED_FIELD_OVERRIDES.items():
                    # User overrides win; only auto-fill if user didn't.
                    if fname in user_overrides:
                        continue
                    if getattr(post_perturb, fname) == forbidden:
                        effective_ov[fname] = replacement

                # Snapshot the resolved Discretisation alongside.
                effective_disc = _effective_disc(base_disc, effective_ov)

                # Snapshots for reproducibility.
                shutil.copy(params_file, run_dir / 'params-snapshot.txt')
                (run_dir / 'preset-snapshot.json').write_text(
                    json.dumps({
                        'row_label': row.row_label,
                        'base_preset': row.base_preset,
                        'perturbation': {
                            k: str(v) for k, v in row.perturbation.items()
                        },
                        'overrides': {
                            k: str(v) for k, v in effective_ov.items()
                        },
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
                    preset_name=row.base_preset,
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

                print(f'  [{row.row_label}/{params_basename}] launching '
                      f'{args.iters}×{args.saves} iters',
                      flush=True)
                timeout = args.per_run_timeout if args.per_run_timeout > 0 else None
                timed_out = False
                with stdout_path.open('w') as so, stderr_path.open('w') as se:
                    try:
                        proc = subprocess.run(
                            cmd, cwd=REPO_ROOT, stdout=so, stderr=se,
                            env=env, text=True, timeout=timeout,
                        )
                    except subprocess.TimeoutExpired:
                        timed_out = True
                        proc = None

                # Per-run progress.json status correction (CLAUDE.md
                # §Stale progress cleanup). On non-zero exit OR timeout
                # we mark the per-run progress.json as failed; on zero
                # exit the CLI's ProgressReporter context manager has
                # already marked it 'completed' on its way out.
                run_progress = run_dir / 'progress.json'
                if timed_out:
                    msg = f'silicoshark timed out after {timeout}s (runaway protection)'
                    _mark_progress_failed(run_progress, msg)
                    print(f'    TIMEOUT: {msg}', flush=True)
                elif proc.returncode != 0:
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

                metrics['row_label'] = row.row_label
                metrics['base_preset'] = row.base_preset
                metrics['perturbation'] = {
                    k: str(v) for k, v in row.perturbation.items()
                }
                metrics['params_file'] = str(params_file)
                metrics['returncode'] = -1 if timed_out else proc.returncode
                metrics['timed_out'] = timed_out
                (run_dir / 'metrics.json').write_text(
                    json.dumps(metrics, indent=2) + '\n'
                )

                results[row.row_label][params_basename] = metrics

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
    if args.single_field_mode:
        _write_single_field_table(
            out_dir, rows, results, params_basenames,
            args.single_field_mode,
        )
    else:
        _write_results_table(
            out_dir, results,
            [r.row_label for r in rows], params_basenames,
        )
    print(f'Wrote {out_dir / "results.json"} and {out_dir / "results-table.md"}',
          flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
