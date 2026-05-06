"""Figure 1: code-base lineage of the tooth-morphogenesis model.

Renders an 11-node, 10-edge tree showing how the 2010 Salazar-Ciudad and
Jernvall *Nature* model's FORTRAN implementation has been translated and
forked over fifteen years, with this work (silicoshark v2) at one of the
leaf positions. Nodes are placed manually (no graphviz dependency); each
edge is annotated with the kind of transition.

Run from the repo root:

    .venv/bin/python scripts/figures/fig1-lineage.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = REPO_ROOT / 'docs' / 'figures'


# Node (x, y, label, kind) — kind controls fill/border styling.
# Coordinate system: x in [0, 12], y in [0, 14] (top is high y).
NODES = {
    'humppa': (
        6.0, 13.0,
        'humppa_translate.f90\n(Catalan, ~2010)',
        'origin',
    ),
    'jernvall': (
        1.5, 10.5,
        'jernvall-lab\nToothMaker\n(C++ port)',
        'fork',
    ),
    'tgrohens': (
        1.5, 8.0,
        'tgrohens fork',
        'fork',
    ),
    '13f90': (
        6.0, 10.5,
        '13.f90\n(English rename;\ntranslation drift)',
        'fork',
    ),
    'silicoshark_f': (
        10.5, 10.5,
        'silicoshark FORTRAN\n(shark extension)',
        'fork',
    ),
    'zimm': (
        10.5, 8.0,
        'Zimm et al. 2023\nPNAS',
        'fork',
    ),
    'coreop': (
        6.0, 8.0,
        'coreop2d.py (Path A)\n(Python translation;\nsilent rewrite of eq.17;\nFMA-fluke discovery)',
        'fork',
    ),
    'sv1': (
        6.0, 5.0,
        'silicoshark v1\n(numpy re-impl;\nmesh-degeneration finding)',
        'fork',
    ),
    'sv2': (
        6.0, 2.0,
        'silicoshark v2\n(configurable\ncomparative study;\nTHIS WORK)',
        'this_work',
    ),
}


# (src, dst, label) for each edge.
EDGES = [
    ('humppa', 'jernvall', 'C++ port'),
    ('jernvall', 'tgrohens', 'fork'),
    ('humppa', '13f90', 'rename'),
    ('humppa', 'silicoshark_f', 'shark extension'),
    ('silicoshark_f', 'zimm', 'publication'),
    ('13f90', 'coreop', 'Python translation'),
    ('coreop', 'sv1', 'numpy re-implementation'),
    ('sv1', 'sv2', 'configurable comparative'),
]


STYLE = {
    'origin': {
        'facecolor': '#fdf2d0',
        'edgecolor': '#7a5c00',
        'linewidth': 1.4,
    },
    'fork': {
        'facecolor': '#eaf2fb',
        'edgecolor': '#214872',
        'linewidth': 1.0,
    },
    'this_work': {
        'facecolor': '#d8efd8',
        'edgecolor': '#1f5d1f',
        'linewidth': 2.4,
    },
}


def draw_node(ax, x: float, y: float, label: str, kind: str) -> None:
    style = STYLE[kind]
    # Box dimensions are sized to the label content.
    n_lines = label.count('\n') + 1
    width = 3.0
    height = 0.45 * n_lines + 0.5
    box = FancyBboxPatch(
        (x - width / 2, y - height / 2),
        width,
        height,
        boxstyle='round,pad=0.05,rounding_size=0.15',
        facecolor=style['facecolor'],
        edgecolor=style['edgecolor'],
        linewidth=style['linewidth'],
        zorder=2,
    )
    ax.add_patch(box)
    weight = 'bold' if kind == 'this_work' else 'normal'
    ax.text(
        x, y, label,
        ha='center', va='center',
        fontsize=8.5,
        fontweight=weight,
        zorder=3,
    )


def draw_edge(ax, src: tuple[float, float], dst: tuple[float, float], label: str) -> None:
    arrow = FancyArrowPatch(
        src, dst,
        arrowstyle='->',
        mutation_scale=12,
        color='#444444',
        linewidth=0.9,
        connectionstyle='arc3,rad=0.0',
        zorder=1,
    )
    ax.add_patch(arrow)
    # Label at the midpoint, with a slight perpendicular offset.
    mx = (src[0] + dst[0]) / 2
    my = (src[1] + dst[1]) / 2
    dx = dst[0] - src[0]
    dy = dst[1] - src[1]
    length = (dx ** 2 + dy ** 2) ** 0.5 or 1.0
    nx, ny = -dy / length, dx / length  # unit normal
    offset = 0.18
    ax.text(
        mx + nx * offset,
        my + ny * offset,
        label,
        ha='center', va='center',
        fontsize=7.5,
        style='italic',
        color='#333333',
        bbox=dict(facecolor='white', edgecolor='none', pad=1.0, alpha=0.85),
        zorder=4,
    )


def edge_anchor(node: str, towards: tuple[float, float]) -> tuple[float, float]:
    """Return the point on the node box edge closest to `towards`.

    A simple heuristic that stops arrows from running into node interiors:
    if the target is mostly above/below, anchor on the top/bottom edge;
    otherwise anchor on the left/right edge.
    """
    x, y, label, _ = NODES[node]
    n_lines = label.count('\n') + 1
    half_w = 1.5
    half_h = 0.5 * (0.45 * n_lines + 0.5)
    dx = towards[0] - x
    dy = towards[1] - y
    if abs(dy) * half_w >= abs(dx) * half_h:
        # vertical dominates
        return (x + dx * half_h / max(abs(dy), 1e-6), y + (half_h if dy > 0 else -half_h))
    return (x + (half_w if dx > 0 else -half_w), y + dy * half_w / max(abs(dx), 1e-6))


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 7))
    ax.set_xlim(-0.5, 12.5)
    ax.set_ylim(0.0, 14.5)
    ax.set_aspect('equal')
    ax.axis('off')

    for src_key, dst_key, label in EDGES:
        sx, sy, _, _ = NODES[src_key]
        dx, dy, _, _ = NODES[dst_key]
        src = edge_anchor(src_key, (dx, dy))
        dst = edge_anchor(dst_key, (sx, sy))
        draw_edge(ax, src, dst, label)

    for key, (x, y, label, kind) in NODES.items():
        draw_node(ax, x, y, label, kind)

    ax.set_title(
        'Code-base lineage of the 2010 tooth-morphogenesis model',
        fontsize=11, pad=8,
    )

    plt.tight_layout()
    pdf_path = OUT_DIR / 'fig1-lineage.pdf'
    png_path = OUT_DIR / 'fig1-lineage.png'
    fig.savefig(pdf_path, bbox_inches='tight')
    fig.savefig(png_path, bbox_inches='tight', dpi=180)
    plt.close(fig)
    print(f'wrote {pdf_path}')
    print(f'wrote {png_path}')


if __name__ == '__main__':
    main()
