"""Figure 3: single-field disentanglement heatmap.

Reads the single-field disentanglement results.json. Builds a (fields ×
directions) matrix where rows are field names (sorted alphabetically),
columns are knock-down (LEGACY_FORTRAN minus this field) and knock-up
(PATH_B_DEFAULT plus this field). Cells carry the absolute cell-count
delta from the corresponding anchor (LEGACY_FORTRAN's plateau = 60 for
knock-down; PATH_B_DEFAULT's plateau = 37 for knock-up). NaN / failed
runs are marked with red hatching.

Run from the repo root:

    .venv/bin/python scripts/figures/fig3-disentanglement.py
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

REPO_ROOT = Path(__file__).resolve().parents[2]
RESULTS_JSON = (
    REPO_ROOT
    / 'experiments'
    / 'discretisation-study'
    / 'single-field-disentanglement'
    / 'results.json'
)
OUT_DIR = REPO_ROOT / 'docs' / 'figures'

KNOCKDOWN_RE = re.compile(r'^LEGACY_FORTRAN_minus_(.+)$')
KNOCKUP_RE = re.compile(r'^PATH_B_DEFAULT_plus_(.+)$')

ANCHOR_KD = 60.0  # LEGACY_FORTRAN plateau on the seal example
ANCHOR_KU = 37.0  # PATH_B_DEFAULT plateau on the seal example


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with RESULTS_JSON.open() as f:
        # The results-table contains literal `NaN` tokens; tolerate them.
        text = f.read()
    text = text.replace('NaN', 'null')
    data = json.loads(text)

    knockdown: dict[str, float | None] = {}
    knockup: dict[str, float | None] = {}
    for key, by_params in data.items():
        # Each row is keyed by params-file basename; the seal example is
        # the only one in this file.
        seal = by_params.get('seal')
        if seal is None:
            continue
        plateau = seal.get('cell_count_plateau')
        m_kd = KNOCKDOWN_RE.match(key)
        m_ku = KNOCKUP_RE.match(key)
        if m_kd is not None:
            knockdown[m_kd.group(1)] = plateau
        elif m_ku is not None:
            knockup[m_ku.group(1)] = plateau

    fields = sorted(set(knockdown) | set(knockup))
    n_rows = len(fields)
    matrix = np.full((n_rows, 2), np.nan)
    failure_mask = np.zeros((n_rows, 2), dtype=bool)
    annotations: list[list[str]] = [['', ''] for _ in range(n_rows)]

    for r, field in enumerate(fields):
        for c, (anchor, source) in enumerate(((ANCHOR_KD, knockdown), (ANCHOR_KU, knockup))):
            plateau = source.get(field)
            if plateau is None or (isinstance(plateau, float) and math.isnan(plateau)):
                failure_mask[r, c] = True
                annotations[r][c] = 'NaN'
            else:
                delta = abs(float(plateau) - anchor)
                matrix[r, c] = delta
                # Show '0' as '0' (still informative) and integer deltas.
                annotations[r][c] = f'{delta:.0f}' if delta == round(delta) else f'{delta:.1f}'

    # Mask the failure cells so the colourmap does not see them.
    masked = np.ma.array(matrix, mask=failure_mask)

    fig, ax = plt.subplots(figsize=(8, 6))
    cmap = plt.colormaps['YlGnBu'].copy()
    vmax = float(np.nanmax(matrix)) if np.isfinite(np.nanmax(matrix)) else 1.0
    norm = Normalize(vmin=0.0, vmax=vmax)
    im = ax.imshow(masked, cmap=cmap, norm=norm, aspect='auto')

    # Overlay red hatching on failure cells.
    for r in range(n_rows):
        for c in range(2):
            if failure_mask[r, c]:
                rect = Rectangle(
                    (c - 0.5, r - 0.5), 1.0, 1.0,
                    facecolor='#ffd6d6',
                    edgecolor='#aa0000',
                    hatch='///',
                    linewidth=0.8,
                    zorder=2,
                )
                ax.add_patch(rect)

    # Annotate each cell with its delta value or 'NaN'.
    for r in range(n_rows):
        for c in range(2):
            text = annotations[r][c]
            if failure_mask[r, c]:
                colour = '#990000'
                weight = 'bold'
            else:
                # Choose dark or light text based on the cell's colour-bar value.
                value = matrix[r, c]
                rel = value / vmax if vmax > 0 else 0.0
                colour = 'white' if rel > 0.55 else '#1a1a1a'
                weight = 'bold' if value >= max(1.0, vmax * 0.5) else 'normal'
            ax.text(
                c, r, text,
                ha='center', va='center',
                fontsize=9, color=colour, fontweight=weight,
                zorder=3,
            )

    ax.set_xticks([0, 1])
    ax.set_xticklabels(
        ['knock-down\n(LEGACY_FORTRAN\nminus field)',
         'knock-up\n(PATH_B_DEFAULT\nplus field)'],
        fontsize=9,
    )
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([f'`{field}`' for field in fields], fontsize=9)
    ax.set_xlabel('')
    ax.set_ylabel('Discretisation field')
    ax.set_title(
        'Single-field disentanglement on examples/seal.txt\n'
        '(absolute cell-count delta from anchor)',
        fontsize=10, pad=10,
    )

    # Move tick marks off the edges for readability.
    ax.tick_params(axis='x', length=0)
    ax.tick_params(axis='y', length=0)

    cbar = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.04)
    cbar.set_label('|plateau − anchor| (cells)', fontsize=9)

    plt.tight_layout()

    pdf_path = OUT_DIR / 'fig3-disentanglement.pdf'
    png_path = OUT_DIR / 'fig3-disentanglement.png'
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, bbox_inches='tight', dpi=180)
    plt.close(fig)
    print(f'wrote {pdf_path}')
    print(f'wrote {png_path}')


if __name__ == '__main__':
    main()
