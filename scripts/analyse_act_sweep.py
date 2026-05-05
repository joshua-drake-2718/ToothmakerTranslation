"""Aggregate the Act-sweep results into a tidy CSV + cusp-vs-Act plot.

Reads `experiments/discretisation-study/act-sweep-2014/results/results.json`
(produced by `scripts/run-discretisation-study.py` over the 16 Act
variants × 5 presets), pivots it into a long-form table indexed by
(preset, Act), and writes:

- `experiments/discretisation-study/act-sweep-2014/cusp-vs-act.csv` —
  one row per (preset, Act) with cusps, cell-count plateau, regime.
- `docs/figures/fig5-act-sweep/fig5-act-sweep.{pdf,png}` — cusp count
  vs Act under each of the five named presets.

Run from the repo root: `.venv/bin/python scripts/analyse_act_sweep.py`.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parent.parent
RESULTS = REPO_ROOT / 'experiments' / 'discretisation-study' / 'act-sweep-2014' / 'results' / 'results.json'
CSV_OUT = REPO_ROOT / 'experiments' / 'discretisation-study' / 'act-sweep-2014' / 'cusp-vs-act.csv'
FIG_DIR = REPO_ROOT / 'docs' / 'figures' / 'fig5-act-sweep'


PRESETS = [
    'HUMPPA_LITERAL',
    'LEGACY_FORTRAN',
    'PAPER_2010',
    'PAPER_LITERAL_2010',
    'PATH_B_DEFAULT',
]
PRESET_COLOURS = {
    'HUMPPA_LITERAL': '#d62728',
    'LEGACY_FORTRAN': '#ff7f0e',
    'PAPER_2010': '#2ca02c',
    'PAPER_LITERAL_2010': '#1f77b4',
    'PATH_B_DEFAULT': '#9467bd',
}
PRESET_MARKERS = {
    'HUMPPA_LITERAL': 'o',
    'LEGACY_FORTRAN': 's',
    'PAPER_2010': '^',
    'PAPER_LITERAL_2010': 'D',
    'PATH_B_DEFAULT': 'v',
}


def parse_act_from_basename(basename: str) -> float:
    return float(basename.removeprefix('Act_'))


def main() -> None:
    raw = json.loads(RESULTS.read_text())

    rows: list[tuple[str, float, int | None, float | None, str]] = []
    for preset in PRESETS:
        per_preset = raw.get(preset, {})
        for params_basename, entry in per_preset.items():
            act = parse_act_from_basename(params_basename)
            cusps = entry.get('cusp_count')
            plateau = entry.get('cell_count_plateau')
            regime = entry.get('regime', '')
            rows.append((preset, act, cusps, plateau, regime))
    rows.sort(key=lambda r: (PRESETS.index(r[0]), r[1]))

    CSV_OUT.parent.mkdir(parents=True, exist_ok=True)
    with CSV_OUT.open('w') as f:
        f.write('preset,Act,cusp_count,cell_count_plateau,regime\n')
        for preset, act, cusps, plateau, regime in rows:
            cusps_s = '' if cusps is None else str(cusps)
            plateau_s = '' if plateau is None or (isinstance(plateau, float) and np.isnan(plateau)) else f'{plateau:g}'
            f.write(f'{preset},{act:g},{cusps_s},{plateau_s},{regime}\n')
    print(f'Wrote {CSV_OUT}')

    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=120)
    for preset in PRESETS:
        xs: list[float] = []
        ys: list[float] = []
        for r in rows:
            if r[0] != preset:
                continue
            if r[2] is None:
                continue
            xs.append(r[1])
            ys.append(float(r[2]))
        ax.plot(
            xs, ys,
            label=preset,
            color=PRESET_COLOURS[preset],
            marker=PRESET_MARKERS[preset],
            linewidth=1.5,
            markersize=6,
        )
    ax.set_xlabel('k_act (Act parameter)')
    ax.set_ylabel('Cusp count at iteration 2500')
    ax.set_title(
        'Cusp count vs Act on examples/wt-tribosphenic-2014.txt (ina=0.5)\n'
        '14,000 iterations not used here; this is the comparison-study '
        '2,500-iter cadence for cross-comparability'
    )
    ax.set_xticks(np.arange(0.1, 1.7, 0.1))
    ax.grid(True, alpha=0.3)
    ax.legend(loc='best', fontsize=8)
    fig.tight_layout()
    for ext in ('pdf', 'png'):
        out = FIG_DIR / f'fig5-act-sweep.{ext}'
        fig.savefig(out)
        print(f'Wrote {out}')
    plt.close(fig)


if __name__ == '__main__':
    main()
