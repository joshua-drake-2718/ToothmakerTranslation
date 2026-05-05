"""Mesh — triangulation, neighbour lists, Laplacian operators.

Path B's mesh handling. The 2010 paper specifies a triangular mesh with
"flux proportional to contact area" between neighbouring cells, but does
not give the discrete formula. Two construction paths are supported:

- `Mesh.from_positions(positions)`: re-triangulate every iteration via
  `scipy.spatial.Delaunay` over the (x, y) projection. Used by the
  `topology='delaunay_each_step'` discretisation option (Path B v1
  default; preserved as a comparison anchor).

- `Mesh.from_topology(top, positions)`: build from an already-maintained
  `Topology` (CSR adjacency graph) plus current positions. The CSR
  defines neighbours; positions only contribute the edge-length weights
  and triangle-area calculations. Used by the
  `topology='static_with_local_update'` option, which avoids the
  v1 mesh-degeneration failure mode (charter §Mesh and topology).

Two Laplacian forms are exposed, selected by `Discretisation.laplacian`
at the call site (typically `simulator.step` or `reaction.step`):

- `laplacian(u)` — length-weighted graph Laplacian, Fick's-law form:
      L u[i] = sum_{j in N(i)}  (u[j] - u[i]) / |p_j - p_i|
  Simple, conservation-preserving, used as the default.

- `cotangent_laplacian(u)` — textbook discrete Laplace–Beltrami:
      L u[i] = (1 / 2 A_i) sum_j (cot α_ij + cot β_ij) (u[j] - u[i])
  α, β are the angles opposite edge (i, j) in its two adjacent
  triangles. A_i is the lumped (barycentric one-third) area at i.

The FORTRAN's per-edge `pes` margins / `area_p` z-layer weighting (the
`fortran_margins` discretisation option) is not implemented in this
phase; calling it raises `NotImplementedError`.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial import Delaunay

from .topology import Topology


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
    kind: str = 'length_weighted'
    """Which Laplacian `mesh.laplacian(u)` evaluates.

    Selected by the simulator from `Discretisation.laplacian` before
    calling reaction code. Allowed values: `'length_weighted'` (the
    default Fick's-law form), `'cotangent'` (textbook Laplace-Beltrami;
    requires `positions` to be set), `'fortran_margins'` (deferred —
    raises NotImplementedError). Reaction code remains unaware of the
    branch: it always calls `mesh.laplacian(u)`.
    """
    positions: np.ndarray | None = None
    """Per-cell positions snapshot, retained for the cotangent
    Laplacian's triangle-area calculations. `None` for length-weighted
    meshes, set by `from_topology`/`from_positions` for cotangent.
    """

    @classmethod
    def from_positions(
        cls,
        positions: np.ndarray,
        *,
        kind: str = 'length_weighted',
    ) -> 'Mesh':
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
            kind=kind,
            positions=positions if kind == 'cotangent' else None,
        )

    @classmethod
    def from_topology(
        cls,
        top: 'Topology',
        positions: np.ndarray,
        *,
        kind: str = 'length_weighted',
    ) -> 'Mesh':
        """Build a Mesh from a `Topology` (CSR neighbour graph) plus
        per-cell positions.

        The CSR defines adjacency; the positions only contribute the
        per-edge length weights, the triangle list (derived as 3-cycles
        in the CSR; see `Topology.triangles`), and the border flag.
        Used by the `topology='static_with_local_update'` option.
        """
        n = positions.shape[0]
        if top.num_cells != n:
            raise ValueError(
                f'topology has {top.num_cells} cells but positions has {n}'
            )

        neigh_starts, neigh_idx = top.to_csr()
        triangles = top.triangles()

        # Symmetric (rows, cols) for CSR-style scatter operations.
        rows = np.repeat(
            np.arange(n, dtype=np.int64),
            np.diff(neigh_starts),
        )

        delta = positions[neigh_idx] - positions[rows]
        edge_len = np.linalg.norm(delta, axis=1)
        edge_w = np.where(
            edge_len > 1e-12, 1.0 / np.maximum(edge_len, 1e-12), 0.0
        )

        diag_w = np.zeros(n, dtype=np.float64)
        np.add.at(diag_w, rows, -edge_w)

        is_border = top.is_border()

        return cls(
            triangles=triangles,
            neigh_starts=neigh_starts,
            neigh_idx=neigh_idx,
            edge_w=edge_w,
            diag_w=diag_w,
            is_border=is_border,
            kind=kind,
            positions=positions if kind == 'cotangent' else None,
        )

    def laplacian(self, u: np.ndarray) -> np.ndarray:
        """Apply the discrete Laplacian selected by `self.kind`.

        Reaction code is unaware of the branch: it calls `mesh.laplacian(u)`
        and the simulator has set `self.kind` from `Discretisation.laplacian`
        before constructing the mesh.

        - `'length_weighted'` — Fick's-law graph Laplacian (default).
        - `'cotangent'` — discrete Laplace-Beltrami; requires positions
          to have been retained on the mesh.
        - `'fortran_margins'` — raises NotImplementedError; the FORTRAN
          per-edge margin scheme is documented as future work.
        """
        if self.kind == 'length_weighted':
            return self._length_weighted_laplacian(u)
        if self.kind == 'cotangent':
            if self.positions is None:
                raise ValueError(
                    "Mesh.kind='cotangent' requires positions to be set "
                    'on the mesh; got positions=None.'
                )
            return self.cotangent_laplacian(u, self.positions)
        if self.kind == 'fortran_margins':
            return self.fortran_margins_laplacian(u)
        raise ValueError(f'unknown Mesh.kind: {self.kind!r}')

    def _length_weighted_laplacian(self, u: np.ndarray) -> np.ndarray:
        """Length-weighted graph Laplacian.

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

    def cotangent_laplacian(
        self, u: np.ndarray, positions: np.ndarray
    ) -> np.ndarray:
        """Apply the discrete cotangent Laplace–Beltrami operator.

        L u[i] = (1 / 2 A_i) sum_{j in N(i)} (cot α_ij + cot β_ij) (u[j] - u[i])

        where α_ij and β_ij are the angles opposite edge (i, j) in the
        two triangles sharing that edge. For boundary edges with one
        adjacent triangle, only one cotangent contributes. A_i is the
        lumped barycentric area at i (one-third of each incident
        triangle's area).

        Parameters
        ----------
        u : (N,) float array
            Per-cell scalar field.
        positions : (N, 3) float array
            Cell centres. Triangles are taken from `self.triangles`.

        Returns
        -------
        Lu : (N,) float array
            Discrete Laplacian of u at each cell.
        """
        n = u.shape[0]
        tris = self.triangles
        if tris.shape[0] == 0:
            return np.zeros_like(u)

        # Per-triangle vertex positions and edge vectors.
        i_idx = tris[:, 0]
        j_idx = tris[:, 1]
        k_idx = tris[:, 2]
        pi = positions[i_idx]
        pj = positions[j_idx]
        pk = positions[k_idx]

        # Edge vectors. e_ij = p_j - p_i (the edge opposite vertex k);
        # cot of the angle at k applies to edge (i, j).
        e_ij = pj - pi
        e_jk = pk - pj
        e_ki = pi - pk

        # Cotangent of the angle at vertex v (in radians) opposite edge e:
        #   cot(theta) = (a . b) / |a x b|
        # where a, b are the two edges meeting at v. Use the standard
        # identity that for triangle vertices i, j, k with edges
        # u = p_j - p_i and v = p_k - p_i meeting at i:
        #   cot(angle at i) = (u . v) / |u x v|
        # The angle at i is opposite edge (j, k).
        def _cot(a: np.ndarray, b: np.ndarray) -> np.ndarray:
            cross = np.cross(a, b)
            denom = np.linalg.norm(cross, axis=1)
            return np.where(
                denom > 1e-15, np.einsum('ij,ij->i', a, b) / np.maximum(denom, 1e-15), 0.0
            )

        # angle at i is opposite edge (j, k); contributes to weight of (j, k)
        cot_at_i = _cot(pj - pi, pk - pi)
        # angle at j is opposite edge (i, k); contributes to weight of (i, k)
        cot_at_j = _cot(pi - pj, pk - pj)
        # angle at k is opposite edge (i, j); contributes to weight of (i, j)
        cot_at_k = _cot(pi - pk, pj - pk)

        # Triangle area for the lumped vertex area: A_tri = 0.5 |e_ij x e_ik|
        area_tri = 0.5 * np.linalg.norm(np.cross(e_ij, pk - pi), axis=1)

        # Accumulate per-vertex barycentric area: A_i += A_tri / 3.
        vertex_area = np.zeros(n, dtype=np.float64)
        np.add.at(vertex_area, i_idx, area_tri / 3.0)
        np.add.at(vertex_area, j_idx, area_tri / 3.0)
        np.add.at(vertex_area, k_idx, area_tri / 3.0)

        # Cotangent weights per (undirected) edge, accumulated across
        # both adjacent triangles. We accumulate symmetric edge
        # contributions to a sparse representation indexed by the CSR.
        # Off-diagonal edge weights: build a dict-like accumulator on
        # (i, j) ordered pairs, then materialise into the CSR layout.
        # For vectorisation, we accumulate into a flat (N,N)-equivalent
        # via three np.add.at calls per directed edge.
        cot_sum = np.zeros(self.neigh_idx.shape[0], dtype=np.float64)

        # Build a (rows -> position-in-CSR-block) lookup so we can map
        # (i, j) to the slot in `cot_sum` that holds edge (i, j). The
        # CSR rows are sorted by row, and within each row neigh_idx is
        # sorted ascending, so we can binary-search within each row.
        def _edge_slot(i_arr: np.ndarray, j_arr: np.ndarray) -> np.ndarray:
            """Return the CSR slot for each (i, j) directed edge."""
            slots = np.empty(i_arr.shape[0], dtype=np.int64)
            for e in range(i_arr.shape[0]):
                row_i = int(i_arr[e])
                row_j = int(j_arr[e])
                start = self.neigh_starts[row_i]
                end = self.neigh_starts[row_i + 1]
                block = self.neigh_idx[start:end]
                # neigh_idx blocks are sorted ascending.
                pos = int(np.searchsorted(block, row_j))
                slots[e] = start + pos
            return slots

        # Edge (j, k) ↔ angle at i.
        slots_jk = _edge_slot(j_idx, k_idx)
        slots_kj = _edge_slot(k_idx, j_idx)
        np.add.at(cot_sum, slots_jk, cot_at_i)
        np.add.at(cot_sum, slots_kj, cot_at_i)

        # Edge (i, k) ↔ angle at j.
        slots_ik = _edge_slot(i_idx, k_idx)
        slots_ki = _edge_slot(k_idx, i_idx)
        np.add.at(cot_sum, slots_ik, cot_at_j)
        np.add.at(cot_sum, slots_ki, cot_at_j)

        # Edge (i, j) ↔ angle at k.
        slots_ij = _edge_slot(i_idx, j_idx)
        slots_ji = _edge_slot(j_idx, i_idx)
        np.add.at(cot_sum, slots_ij, cot_at_k)
        np.add.at(cot_sum, slots_ji, cot_at_k)

        # cot_sum[e] is the sum of cotangents over both adjacent
        # triangles for the directed edge at slot e (or one cotangent
        # for boundary edges). The discrete Laplacian:
        #   L u[i] = (1 / 2 A_i) sum_j w_ij (u[j] - u[i])
        # where w_ij is the un-doubled (cot α + cot β). Our cot_sum
        # is exactly (cot α + cot β) per directed edge.
        rows = np.repeat(
            np.arange(n, dtype=np.int64),
            np.diff(self.neigh_starts),
        )

        out = np.zeros(n, dtype=np.float64)
        # off-diagonal: sum_j w_ij u[j]
        np.add.at(out, rows, cot_sum * u[self.neigh_idx])
        # diagonal: u[i] * (-sum_j w_ij)
        diag = np.zeros(n, dtype=np.float64)
        np.add.at(diag, rows, cot_sum)
        out -= diag * u

        # Lumped area normalisation: divide by (2 A_i).
        safe_area = np.maximum(vertex_area, 1e-15)
        out = out / (2.0 * safe_area)
        # Cells with zero area (no incident triangles) get zero.
        out[vertex_area <= 1e-15] = 0.0
        return out

    def fortran_margins_laplacian(self, u: np.ndarray) -> np.ndarray:
        """The FORTRAN's per-edge `pes`-margin Laplacian. Not implemented
        in this phase; documented in the discretisation audit as its own
        sub-project (the FORTRAN scheme is interleaved with z-layer
        mesenchyme handling and is not a drop-in textbook discretisation).
        """
        raise NotImplementedError(
            "Discretisation.laplacian='fortran_margins' is not implemented "
            "in this phase. The FORTRAN's per-edge margin scheme is "
            'documented as a future-work item in the discretisation audit.'
        )
