"""Figure 4: eq. 17 inhibitor source clause across the code-base lineage.

Visualises how a single equation morphs across three implementations
(paper -> 13.f90 -> coreop2d.py) and is then captured as three named
silicoshark Discretisation options bound to named presets. The
silent-rewrite edge between 13.f90 and coreop2d.py is highlighted to
mark the translation drift discussed in the briefing's section 6.

Run from the repo root:

    .venv/bin/python scripts/figures/fig4-eq17-lineage.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = REPO_ROOT / 'docs' / 'figures'


# Node layout. Coordinate system: x in [0, 12], y in [0, 28] (top is high y).
# Each node has a heading (sans-serif, bold), a code/equation block
# (monospace), and a short trailing label (italic). The render function
# below positions these three elements within the box.
NODES = {
    'paper': {
        'x': 6.0, 'y': 25.0,
        'heading': 'Paper Eq. 17  (Salazar-Ciudad and Jernvall 2010, p. 586)',
        'code': (
            r'$\partial$[Inh]/$\partial$t  =  [Act]'
            r'  $-$  k$_{deg}$[Inh]  +  k$_{di}$ $\nabla^{2}$[Inh]'
            '\n'
            r'(applies when knot OR d$_i$ $\geq$ k$_{int}$)'
        ),
        'note': 'Source: activator concentration',
        'kind': 'paper',
    },
    'fortran': {
        'x': 6.0, 'y': 18.0,
        'heading': '13.f90:487  (humppa_translate.f90:617 byte-identical)',
        'code': (
            'hq3d(i,1,2) = q3d(i,1,1)*DiffState(i)\n'
            '              - mu * q3d(i,1,2)'
        ),
        'note': 'Source: [Act] x d_i  (concentration; smooth d_i ramp added)',
        'kind': 'fortran',
    },
    'python': {
        'x': 6.0, 'y': 11.0,
        'heading': 'coreop2d.py:647  (Path A Python translation)',
        'code': (
            'hq3d[i, 0, 1] = hq3d[i, 0, 0] * cls.diff_state[i] \\\n'
            '                - cls.Deg * cls.q3d[i, 0, 1]'
        ),
        'note': 'Source: (rate of [Act]) x d_i  -- rate variable, NOT concentration',
        'kind': 'python',
    },
    'silicoshark': {
        'x': 6.0, 'y': 4.0,
        'heading': 'silicoshark v2 Discretisation.eq17_inh_source  (THIS WORK)',
        'code': (
            "  'act_concentration'   ->  PAPER_2010 / PATH_B_DEFAULT\n"
            "  'act_times_di'        ->  LEGACY_FORTRAN / HUMPPA_LITERAL\n"
            "  'act_rate_times_di'   ->  PATH_A_REWRITE"
        ),
        'note': 'All three forms named, cited, and exercisable in the comparison study',
        'kind': 'this_work',
    },
}


# Edges: (src, dst, label, kind). 'silent' triggers the highlighted style.
EDGES = [
    (
        'paper', 'fortran',
        'FORTRAN translation\n(d_i ramp added; cosmetic)',
        'normal',
    ),
    (
        'fortran', 'python',
        'Python translation\n*** silent variable substitution ***\nq3d(i,1,1)  ->  hq3d[i,0,0]',
        'silent',
    ),
    (
        'python', 'silicoshark',
        'configurable comparative\n(captured as three named options)',
        'normal',
    ),
]


# Box fills colour-coded by kind. Borders chosen to remain
# distinguishable in greyscale (light vs medium vs heavy linewidth).
STYLE = {
    'paper': {
        'facecolor': '#fdf2d0',
        'edgecolor': '#7a5c00',
        'linewidth': 1.2,
    },
    'fortran': {
        'facecolor': '#eaf2fb',
        'edgecolor': '#214872',
        'linewidth': 1.2,
    },
    'python': {
        'facecolor': '#fde2e2',
        'edgecolor': '#8a1a1a',
        'linewidth': 1.4,
    },
    'this_work': {
        'facecolor': '#d8efd8',
        'edgecolor': '#1f5d1f',
        'linewidth': 2.6,
    },
}


# Edge styles. The 'silent' edge uses a warm orange-red plus a heavier
# arrow so it remains distinguishable from the two black edges in
# greyscale (different linewidth and dash pattern).
EDGE_STYLE = {
    'normal': {
        'color': '#444444',
        'linewidth': 1.0,
        'linestyle': '-',
        'mutation_scale': 14,
        'label_colour': '#333333',
        'label_face': 'white',
        'label_edge': 'none',
    },
    'silent': {
        'color': '#c43d1d',
        'linewidth': 2.4,
        'linestyle': '-',
        'mutation_scale': 20,
        'label_colour': '#7a1f0a',
        'label_face': '#fff1ec',
        'label_edge': '#c43d1d',
    },
}


# Box dimensions. Height is computed from the number of code lines so
# all four boxes share the same width but vary in height as needed.
BOX_HALF_WIDTH = 5.0


def box_height(node: dict) -> float:
    """Total box height in data units; tuned to leave generous margins
    around the heading + code + note triplet."""
    code_lines = node['code'].count('\n') + 1
    # base + per-code-line + heading + note rows
    return 1.6 + 0.55 * code_lines + 1.1


def draw_node(ax, key: str) -> None:
    node = NODES[key]
    x, y = node['x'], node['y']
    style = STYLE[node['kind']]
    h = box_height(node)
    w = 2 * BOX_HALF_WIDTH
    box = FancyBboxPatch(
        (x - BOX_HALF_WIDTH, y - h / 2),
        w, h,
        boxstyle='round,pad=0.06,rounding_size=0.18',
        facecolor=style['facecolor'],
        edgecolor=style['edgecolor'],
        linewidth=style['linewidth'],
        zorder=2,
    )
    ax.add_patch(box)

    # Heading: top of the box.
    ax.text(
        x, y + h / 2 - 0.45,
        node['heading'],
        ha='center', va='top',
        fontsize=9.5,
        fontweight='bold',
        family='sans-serif',
        zorder=3,
    )

    # Code/equation block: monospace, centred vertically below heading.
    ax.text(
        x, y + 0.05,
        node['code'],
        ha='center', va='center',
        fontsize=9.0,
        family='monospace',
        zorder=3,
    )

    # Trailing note: italic, near the bottom of the box.
    ax.text(
        x, y - h / 2 + 0.45,
        node['note'],
        ha='center', va='bottom',
        fontsize=8.5,
        style='italic',
        color='#1a1a1a',
        zorder=3,
    )


def edge_anchor(key: str, towards_y: float) -> tuple[float, float]:
    """Return the top or bottom centre of a node, depending on
    whether the target is above or below it. The lineage is purely
    vertical so we never anchor on the side edges."""
    node = NODES[key]
    h = box_height(node)
    if towards_y > node['y']:
        return (node['x'], node['y'] + h / 2)
    return (node['x'], node['y'] - h / 2)


def draw_edge(ax, src_key: str, dst_key: str, label: str, kind: str) -> None:
    src = edge_anchor(src_key, NODES[dst_key]['y'])
    dst = edge_anchor(dst_key, NODES[src_key]['y'])
    style = EDGE_STYLE[kind]
    arrow = FancyArrowPatch(
        src, dst,
        arrowstyle='-|>',
        mutation_scale=style['mutation_scale'],
        color=style['color'],
        linewidth=style['linewidth'],
        linestyle=style['linestyle'],
        connectionstyle='arc3,rad=0.0',
        zorder=1,
    )
    ax.add_patch(arrow)

    # Label sits at the midpoint of the edge with a small horizontal
    # offset so it doesn't overlap the arrow shaft.
    mx = (src[0] + dst[0]) / 2 + 0.4
    my = (src[1] + dst[1]) / 2
    ax.text(
        mx, my,
        label,
        ha='left', va='center',
        fontsize=8.5,
        style='italic',
        color=style['label_colour'],
        bbox=dict(
            facecolor=style['label_face'],
            edgecolor=style['label_edge'],
            pad=2.5,
            boxstyle='round,pad=0.18',
        ),
        zorder=4,
    )


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # A4 portrait: 8 x 11 inches.
    fig, ax = plt.subplots(figsize=(8, 11))
    fig.patch.set_facecolor('#fbfbf7')  # off-white background
    ax.set_xlim(-0.5, 12.5)
    ax.set_ylim(0.0, 28.5)
    ax.set_aspect('equal')
    ax.axis('off')

    # Edges first so the box fills cover the arrow tails cleanly.
    for src_key, dst_key, label, kind in EDGES:
        draw_edge(ax, src_key, dst_key, label, kind)

    for key in NODES:
        draw_node(ax, key)

    ax.set_title(
        'Figure 4. The eq. 17 inhibitor source clause across the lineage',
        fontsize=11, pad=10,
    )

    plt.tight_layout()
    pdf_path = OUT_DIR / 'fig4-eq17-lineage.pdf'
    png_path = OUT_DIR / 'fig4-eq17-lineage.png'
    fig.savefig(pdf_path, bbox_inches='tight', facecolor=fig.get_facecolor())
    fig.savefig(png_path, bbox_inches='tight', dpi=150, facecolor=fig.get_facecolor())
    plt.close(fig)
    print(f'wrote {pdf_path}')
    print(f'wrote {png_path}')


if __name__ == '__main__':
    main()
