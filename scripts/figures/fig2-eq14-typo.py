"""Figure 2: PAPER_2010 vs PAPER_LITERAL_2010 cusp comparison on the
2014 wild-type tribosphenic mouse parameter set.

Reads the iteration-2500 OFF file from each preset's results directory,
parses cell positions and per-cell colour (cyan = knot, grey = non-knot
in the silicoshark writer), and renders a 2x2 grid: top-down (x, y) on
the top row, side-on (x, z) on the bottom row, one preset per column.

Run from the repo root:

    .venv/bin/python scripts/figures/fig2-eq14-typo.py
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
STUDY_DIR = REPO_ROOT / 'experiments' / 'discretisation-study' / 'cusp-forming'
OUT_DIR = REPO_ROOT / 'docs' / 'figures'

PRESETS = (
    (
        'PAPER_2010',
        'eq.14 with corrected denominator:\n1 + k_inh × [Inh] (FORTRAN)',
    ),
    (
        'PAPER_LITERAL_2010',
        'eq.14 with paper denominator:\n1 + k_inh × [Act] (typo)',
    ),
)


def parse_off(path: Path):
    """Return (positions (n, 3), rgba (n, 4)) for the cell list in the OFF.

    The silicoshark writer emits 'COFF' on line 1, '{nv} {nf} 0' on line 2,
    then `nv` lines of `x y z r g b a` floats.
    """
    with path.open() as f:
        f.readline()  # COFF
        counts = f.readline().strip()
    m = re.match(r'^\s*(\d+)\s+(\d+)\s+0\s*$', counts)
    if m is None:
        raise ValueError(f'unexpected OFF counts line in {path}: {counts!r}')
    n_v = int(m.group(1))
    data = np.loadtxt(path, skiprows=2, max_rows=n_v, usecols=(0, 1, 2, 3, 4, 5, 6))
    return data[:, :3], data[:, 3:]


def is_knot(rgba: np.ndarray) -> np.ndarray:
    """Knot cells are cyan (r=0, g=1, b=1) per silicoshark.io.cell_colour."""
    return (rgba[:, 0] == 0.0) & (rgba[:, 1] == 1.0) & (rgba[:, 2] == 1.0)


def load_metrics(preset: str) -> dict:
    path = STUDY_DIR / preset / 'wt-tribosphenic-2014' / 'metrics.json'
    with path.open() as f:
        return json.load(f)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    runs = []
    for preset, annotation in PRESETS:
        off_path = STUDY_DIR / preset / 'wt-tribosphenic-2014' / '2500_run_.off'
        positions, rgba = parse_off(off_path)
        metrics = load_metrics(preset)
        runs.append({
            'preset': preset,
            'annotation': annotation,
            'positions': positions,
            'knot_mask': is_knot(rgba),
            'metrics': metrics,
        })

    # Common axis ranges so the panels are visually comparable.
    all_pos = np.vstack([r['positions'] for r in runs])
    pad = 0.15
    x_lo, x_hi = all_pos[:, 0].min() - pad, all_pos[:, 0].max() + pad
    y_lo, y_hi = all_pos[:, 1].min() - pad, all_pos[:, 1].max() + pad
    z_lo, z_hi = all_pos[:, 2].min() - 0.02, all_pos[:, 2].max() + 0.02

    fig, axes = plt.subplots(2, 2, figsize=(8, 6))

    knot_colour = '#00bcd4'
    nonknot_colour = '#9e9e9e'

    for col, run in enumerate(runs):
        positions = run['positions']
        knot_mask = run['knot_mask']
        cusp_count = int(run['metrics'].get('cusp_count') or 0)

        # Top row: top-down (x, y).
        ax_top = axes[0, col]
        ax_top.scatter(
            positions[~knot_mask, 0], positions[~knot_mask, 1],
            s=42, c=nonknot_colour, edgecolors='#555555', linewidths=0.4,
            label='non-knot',
        )
        ax_top.scatter(
            positions[knot_mask, 0], positions[knot_mask, 1],
            s=42, c=knot_colour, edgecolors='#00838f', linewidths=0.6,
            label='knot',
        )
        ax_top.set_xlim(x_lo, x_hi)
        ax_top.set_ylim(y_lo, y_hi)
        ax_top.set_aspect('equal')
        ax_top.set_xlabel('x')
        ax_top.set_ylabel('y')
        ax_top.set_title(f'{run["preset"]} — top-down', fontsize=10)
        ax_top.text(
            0.03, 0.97, f'cusps: {cusp_count}',
            transform=ax_top.transAxes,
            fontsize=9, fontweight='bold',
            va='top', ha='left',
            bbox=dict(facecolor='white', edgecolor='#888888', pad=2.5),
        )
        ax_top.text(
            0.03, 0.03, run['annotation'],
            transform=ax_top.transAxes,
            fontsize=7.5,
            va='bottom', ha='left',
            bbox=dict(facecolor='#fff8d0', edgecolor='#caa600', pad=2.5),
        )
        if col == 1:
            ax_top.legend(loc='upper right', fontsize=7, framealpha=0.9)

        # Bottom row: side-on (x, z).
        ax_bot = axes[1, col]
        ax_bot.scatter(
            positions[~knot_mask, 0], positions[~knot_mask, 2],
            s=42, c=nonknot_colour, edgecolors='#555555', linewidths=0.4,
        )
        ax_bot.scatter(
            positions[knot_mask, 0], positions[knot_mask, 2],
            s=42, c=knot_colour, edgecolors='#00838f', linewidths=0.6,
        )
        ax_bot.set_xlim(x_lo, x_hi)
        ax_bot.set_ylim(z_lo, z_hi)
        ax_bot.set_xlabel('x')
        ax_bot.set_ylabel('z')
        ax_bot.set_title(f'{run["preset"]} — side-on', fontsize=10)

    fig.suptitle(
        'eq.14 denominator typo: cusp formation on wt-tribosphenic-2014',
        fontsize=11,
    )
    plt.tight_layout(rect=(0, 0, 1, 0.96))

    pdf_path = OUT_DIR / 'fig2-eq14-typo.pdf'
    png_path = OUT_DIR / 'fig2-eq14-typo.png'
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, bbox_inches='tight', dpi=180)
    plt.close(fig)
    print(f'wrote {pdf_path}')
    print(f'wrote {png_path}')


if __name__ == '__main__':
    main()
