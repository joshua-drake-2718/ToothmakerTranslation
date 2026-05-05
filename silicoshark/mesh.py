"""Mesh — Delaunay triangulation, neighbour lists, edge-weighted Laplacian.

Path B's mesh handling. The 2010 paper specifies a triangular mesh with
"flux proportional to contact area" between neighbouring cells, but does
not give the discrete formula. We use:

- `scipy.spatial.Delaunay` over the (x, y) projection to triangulate the
  cell centres each iteration. This eliminates the FORTRAN's topology
  walk and its associated cell-division accidents (charter §What NOT to
  reproduce).

- A length-weighted graph Laplacian for diffusion:
      L u[i] = sum_{j in N(i)}  (u[j] - u[i]) / |p_j - p_i|
  This is the simplest Fick's-law discretisation that preserves
  conservation between neighbouring cells. It approaches the continuous
  Laplacian as the mesh refines.

The full cotangent Laplace–Beltrami operator (charter recommendation) is
a future upgrade; this v1 form is sufficient for the seal example's
±5% cell-count tolerance.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import Delaunay


@dataclass
class Mesh:
    """Triangulation and adjacency derived from a cell layout.

    Built once per iteration via `Mesh.from_positions(positions)`.

    Attributes
    ----------
    triangles : (T, 3) int array
        Vertex indices of each triangle (Delaunay over (x, y)).
    neigh_starts : (N+1,) int array
        CSR-style row pointer into `neigh_idx` for each cell's
        neighbour block.
    neigh_idx : (E,) int array
        Concatenated sorted neighbour indices. `neigh_idx[neigh_starts[i]
        : neigh_starts[i+1]]` is cell `i`'s neighbour list.
    edge_w : (E,) float array
        Off-diagonal weight for each (i, j) edge: 1 / |p_j - p_i|.
        Indexed identically to `neigh_idx`.
    diag_w : (N,) float array
        Diagonal of the Laplacian: -sum of edge weights for each cell.
    is_border : (N,) bool array
        True if cell `i` has fewer than 6 neighbours. The 2010 paper
        defines border cells as those with fewer epithelial neighbours
        than the interior count (6 in the hex lattice). Used downstream
        for cervical-loop and bias logic.
    """

    triangles: np.ndarray
    neigh_starts: np.ndarray
    neigh_idx: np.ndarray
    edge_w: np.ndarray
    diag_w: np.ndarray
    is_border: np.ndarray

    @classmethod
    def from_positions(cls, positions: np.ndarray) -> 'Mesh':
        n = positions.shape[0]
        if n < 3:
            raise ValueError(f'need at least 3 cells to triangulate, got {n}')
        triangles = Delaunay(positions[:, :2]).simplices.astype(np.int32)

        # Build undirected edge set (i, j) with i < j from triangles, then
        # symmetrise into a CSR-style neighbour list. Numpy unique gives
        # both deduplication and sorting in one pass.
        edges = np.vstack([
            triangles[:, [0, 1]],
            triangles[:, [1, 2]],
            triangles[:, [2, 0]],
        ])
        edges = np.sort(edges, axis=1)
        edges = np.unique(edges, axis=0)

        # Symmetrise: for each (i, j) emit both (i, j) and (j, i) so each
        # cell sees all its neighbours.
        sym_rows = np.concatenate([edges[:, 0], edges[:, 1]])
        sym_cols = np.concatenate([edges[:, 1], edges[:, 0]])

        # Sort by row, then by column within each row, to get CSR layout.
        order = np.lexsort((sym_cols, sym_rows))
        rows_sorted = sym_rows[order]
        cols_sorted = sym_cols[order]

        # CSR row pointer: counts of each row's neighbours, cumsum'd.
        counts = np.bincount(rows_sorted, minlength=n)
        neigh_starts = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(counts, out=neigh_starts[1:])
        neigh_idx = cols_sorted.astype(np.int32)

        # Edge weights: 1 / |p_j - p_i| for each edge in the CSR list.
        delta = positions[neigh_idx] - positions[rows_sorted]
        edge_len = np.linalg.norm(delta, axis=1)
        # Avoid division by zero for degenerate co-located cells; clamp.
        edge_w = np.where(edge_len > 1e-12, 1.0 / np.maximum(edge_len, 1e-12), 0.0)

        # Diagonal: -sum of off-diagonal weights for each cell.
        diag_w = np.zeros(n, dtype=np.float64)
        np.add.at(diag_w, rows_sorted, -edge_w)

        is_border = counts < 6

        return cls(
            triangles=triangles,
            neigh_starts=neigh_starts,
            neigh_idx=neigh_idx,
            edge_w=edge_w,
            diag_w=diag_w,
            is_border=is_border,
        )

    def laplacian(self, u: np.ndarray) -> np.ndarray:
        """Apply the Laplacian to a per-cell scalar field.

        L u[i] = sum_{j in N(i)} (u[j] - u[i]) / |p_j - p_i|
              = (sum_j w_ij u[j]) - (sum_j w_ij) u[i]
              = (off-diag part) + diag * u[i]

        where diag = -sum of off-diag weights.
        """
        n = u.shape[0]
        out = self.diag_w * u
        # Scatter-add: out[row] += w_ij * u[col]. Each edge (row, col) in
        # the CSR contributes w * u[col] to out[row].
        rows = np.repeat(
            np.arange(n, dtype=np.int64),
            np.diff(self.neigh_starts),
        )
        np.add.at(out, rows, self.edge_w * u[self.neigh_idx])
        return out
