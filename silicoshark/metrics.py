"""Metrics for the Path B v2 comparison study.

Each function takes a `State` (or arrays) and returns a scalar (or a
small dict in the case of the envelope). Functions are deliberately
free-standing rather than methods so that the comparison-study runner
can call them on states rebuilt from disk (where re-instantiating a
full State is heavier than passing an `np.ndarray`).

Metrics catalogue, per the v2 charter §Validation strategy:

- **`cusp_count(state)`** — topological cusp count: the number of knot
  cells (`state.knot.sum()`) at the time of measurement.
- **`cusp_width(state, mesh)`** — geometric cusp width: mean pairwise
  Euclidean distance between knot cells that are mesh-adjacent. Zero
  if fewer than two knot cells exist.
- **`cell_count_plateau(history)`** — mean of the last 10% of the
  cell-count history. Falls back to the final value if the history is
  too short for a 10% slice (< 10 entries).
- **`vertex_envelope(positions)`** — bounding-box of the vertex set:
  `{x_min, x_max, y_min, y_max, z_min, z_max}`.
- **`regime(history)`** — qualitative classification of the cell-count
  trajectory: `'NaN'`, `'collapse'`, `'oscillate'`, `'plateau'`, or
  `'monotone-grow'`.

These are the comparison axes from the v2 charter §Validation
strategy: each preset × parameter file combination produces one set
of these metrics, which the markdown table renders for review.
"""
from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover — type-only imports
    from .mesh import Mesh
    from .state import State


def cusp_count(state: 'State') -> int:
    """Number of knot cells at the time of measurement.

    In the 2010 paper an enamel knot is a cell with `[Act] >= 1`; the
    `state.knot` boolean array carries that flag, set irreversibly in
    `step_reaction_diffusion`.
    """
    return int(state.knot.sum())


def cusp_width(state: 'State', mesh: 'Mesh') -> float:
    """Mean pairwise distance between mesh-adjacent knot cells.

    Returns 0.0 if fewer than two knot cells exist (no meaningful
    pairwise distance). Uses `mesh.neigh_idx` / `mesh.neigh_starts`
    (CSR neighbour graph) so only cell pairs that are direct mesh
    neighbours contribute. This is closer to the paper's 'distance
    between adjacent cusps' than an all-pairs mean would be.
    """
    knot = np.asarray(state.knot, dtype=bool)
    if int(knot.sum()) < 2:
        return 0.0
    n = state.num_active
    rows = np.repeat(
        np.arange(n, dtype=np.int64),
        np.diff(mesh.neigh_starts),
    )
    cols = mesh.neigh_idx
    # Restrict to undirected edges (i < j) so each pair counts once.
    keep = rows < cols
    rows_u = rows[keep]
    cols_u = cols[keep]
    # Among undirected edges, keep those with both endpoints knot cells.
    edge_knot = knot[rows_u] & knot[cols_u]
    if not edge_knot.any():
        return 0.0
    a = state.positions[rows_u[edge_knot]]
    b = state.positions[cols_u[edge_knot]]
    dist = np.linalg.norm(a - b, axis=1)
    return float(dist.mean())


def cell_count_plateau(history: list[int]) -> float:
    """Mean of the last 10% of the cell-count history.

    If the history is short (fewer than 10 entries), returns the final
    value as a fallback so callers always get a meaningful number.
    Returned as `float` so subsequent averaging across runs is clean
    (an `int` mean would force re-casting at every aggregate).
    """
    if not history:
        return float('nan')
    n = len(history)
    if n < 10:
        return float(history[-1])
    tail_n = max(1, n // 10)
    tail = history[-tail_n:]
    return float(sum(tail) / len(tail))


def vertex_envelope(positions: np.ndarray) -> dict[str, float]:
    """Bounding-box of the vertex set.

    Parameters
    ----------
    positions : (N, 3) float array of cell centres.

    Returns
    -------
    dict with keys `x_min`, `x_max`, `y_min`, `y_max`, `z_min`, `z_max`.
    """
    if positions.size == 0:
        return {
            'x_min': float('nan'), 'x_max': float('nan'),
            'y_min': float('nan'), 'y_max': float('nan'),
            'z_min': float('nan'), 'z_max': float('nan'),
        }
    return {
        'x_min': float(positions[:, 0].min()),
        'x_max': float(positions[:, 0].max()),
        'y_min': float(positions[:, 1].min()),
        'y_max': float(positions[:, 1].max()),
        'z_min': float(positions[:, 2].min()),
        'z_max': float(positions[:, 2].max()),
    }


def regime(history: list[int], peak_collapse_ratio: float = 0.5) -> str:
    """Classify a cell-count trajectory.

    Categories (in evaluation order):

    - `'NaN'` if any entry is non-finite. Callers signal a divergent
      run by passing `np.nan` somewhere in the history.
    - `'collapse'` if the final value is less than `peak_collapse_ratio`
      times the historical maximum. Captures runs that grew then died.
    - `'oscillate'` if the first-difference of the last half of the
      history changes sign more than once. Captures bouncing runs.
    - `'plateau'` if the last 20% of the history has standard deviation
      below 5% of its mean. Captures the seal-style equilibrium.
    - `'monotone-grow'` otherwise. Catches both steady growth and
      runs that haven't yet equilibrated.

    Empty history returns `'NaN'` (no data is treated as failure).
    """
    if not history:
        return 'NaN'
    arr = np.asarray(history, dtype=np.float64)
    if not np.isfinite(arr).all():
        return 'NaN'
    peak = arr.max()
    final = arr[-1]
    if peak > 0 and final < peak * peak_collapse_ratio:
        return 'collapse'
    n = arr.shape[0]
    if n >= 4:
        # Sign changes in the first-difference of the last half.
        half = arr[n // 2:]
        diffs = np.diff(half)
        # Drop exact-zero diffs (flat plateau steps shouldn't be
        # counted as sign changes).
        nz = diffs[diffs != 0]
        if nz.size >= 2:
            sign_changes = int(np.sum(np.sign(nz[:-1]) != np.sign(nz[1:])))
            if sign_changes > 1:
                return 'oscillate'
    if n >= 5:
        tail_n = max(1, n // 5)
        tail = arr[-tail_n:]
        mean = float(tail.mean())
        std = float(tail.std())
        if mean > 0 and std < 0.05 * mean:
            return 'plateau'
    return 'monotone-grow'
