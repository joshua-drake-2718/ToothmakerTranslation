"""Topology — static neighbour graph with local update on cell division.

The v1 Path B used `scipy.spatial.Delaunay` re-triangulation each
iteration. v1 evidence
(`docs/findings/2026-05-04-path-b-v1-mesh-degeneration.md`) showed
this is unstable on the seal example: convex-hull edges grow long as
border cells migrate outward, triggering inappropriate divisions and
mesh degeneration.

This module implements the alternative: build the neighbour graph
once at t = 0 from `scipy.spatial.Delaunay`, and update it locally
each time a cell divides. New daughter inherits a subset of each
parent's neighbours and replaces the parent-parent direct edge. No
full re-triangulation. This matches the FORTRAN's `add_cell` strategy
without the topology-walk accidents (txungu else-branch, iiii
preservation, panic-and-return) that 13.f90 carries.

Storage uses CSR-style int arrays for vectorised access in
`mesh.py`'s Laplacian and `forces.py`'s edge ops. The graph is
mutated in place as cells are added; this module owns the rebalance.

Charter: `docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md`
§Mesh and topology.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
from scipy.spatial import Delaunay


@dataclass
class Topology:
    """Mutable adjacency graph stored as a list of neighbour-lists.

    The internal representation during incremental updates is a
    Python list of `set[int]` — divides involve neighbour-set
    surgery and the set semantics make duplicate-prevention free.
    Use `to_csr()` to materialise the read-only CSR arrays consumed
    by `mesh.Mesh` and the force / diffusion code.

    Attributes
    ----------
    neighbours : list[set[int]]
        `neighbours[i]` is the set of cell indices adjacent to cell i.
        Symmetric: if `j in neighbours[i]` then `i in neighbours[j]`.
    """

    neighbours: list[set[int]] = field(default_factory=list)

    @property
    def num_cells(self) -> int:
        return len(self.neighbours)

    @classmethod
    def from_positions(cls, positions: np.ndarray) -> 'Topology':
        """Build the initial graph from a Delaunay triangulation of
        the (x, y) projection.

        For the typical hexagonal initial condition this gives every
        interior cell 6 neighbours and every boundary cell 3 or 4.
        """
        n = positions.shape[0]
        if n < 3:
            raise ValueError(f'need >= 3 cells to triangulate, got {n}')
        triangles = Delaunay(positions[:, :2]).simplices
        neighbours: list[set[int]] = [set() for _ in range(n)]
        for tri in triangles:
            a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
            neighbours[a].update((b, c))
            neighbours[b].update((a, c))
            neighbours[c].update((a, b))
        return cls(neighbours=neighbours)

    # --- Local update on cell division ------------------------------

    def insert_daughter(
        self,
        parent_a: int,
        parent_b: int,
    ) -> int:
        """Insert a new cell on the (parent_a, parent_b) edge.

        Topology rewiring (paper SI fig. 1, FORTRAN add_cell):
        - The new daughter's neighbours are {parent_a, parent_b} plus
          the cells in N(parent_a) ∩ N(parent_b) — the two cells
          (or one, on a boundary edge) that share a triangle with
          the (a, b) edge.
        - parent_a and parent_b lose each other as direct neighbours
          and gain the daughter.
        - The cells in N(a) ∩ N(b) gain the daughter as a neighbour
          but do not lose anything else: their original triangle
          (a, b, k) becomes two triangles (a, daughter, k) and
          (daughter, b, k).

        Returns the daughter's cell index, which is `num_cells` at
        call time.
        """
        if parent_a == parent_b:
            raise ValueError(f'parent_a == parent_b == {parent_a}')
        if parent_b not in self.neighbours[parent_a]:
            raise ValueError(
                f'parents {parent_a}, {parent_b} are not direct neighbours'
            )

        shared = self.neighbours[parent_a] & self.neighbours[parent_b]
        # The shared set is the cells that border the (a, b) edge in
        # both triangles. Discard parent endpoints in case of any
        # bookkeeping pathology — they should not appear here in a
        # well-formed graph.
        shared.discard(parent_a)
        shared.discard(parent_b)

        daughter = self.num_cells
        # Daughter's initial neighbour set.
        self.neighbours.append({parent_a, parent_b} | shared)

        # Parents lose the direct (a, b) edge and gain the daughter.
        self.neighbours[parent_a].discard(parent_b)
        self.neighbours[parent_a].add(daughter)
        self.neighbours[parent_b].discard(parent_a)
        self.neighbours[parent_b].add(daughter)

        # Shared cells gain the daughter; they keep their existing
        # links to both parents. (Their triangles are split, not
        # collapsed.)
        for k in shared:
            self.neighbours[k].add(daughter)

        return daughter

    def insert_daughters(
        self,
        parent_pairs: Sequence[tuple[int, int]],
    ) -> list[int]:
        """Insert a batch of daughters in one call.

        Pairs are inserted in order; each daughter's index reflects
        cells added before it. Behaviour is well-defined even when
        consecutive pairs share an endpoint: the topology after
        insertion N reflects insertions 0..N-1, so the (N+1)th
        pair is interpreted against the updated graph.

        Returns the list of new cell indices in insertion order.
        """
        return [self.insert_daughter(a, b) for a, b in parent_pairs]

    # --- Read interfaces --------------------------------------------

    def neighbour_array(self, i: int) -> np.ndarray:
        """Return cell i's neighbours as a sorted int32 array."""
        return np.fromiter(sorted(self.neighbours[i]), dtype=np.int32)

    def to_csr(self) -> tuple[np.ndarray, np.ndarray]:
        """Materialise the CSR pair (`neigh_starts`, `neigh_idx`)
        consumed by `mesh.Mesh`-style code.

        Each row is sorted ascending, so an external `i < j` filter
        for undirected-edge enumeration is correct.

        Returns
        -------
        neigh_starts : (N+1,) int64
            Row pointer; cell i's neighbours occupy
            `neigh_idx[neigh_starts[i] : neigh_starts[i+1]]`.
        neigh_idx : (E,) int32
            Concatenated sorted neighbour indices.
        """
        n = self.num_cells
        counts = np.fromiter((len(s) for s in self.neighbours), dtype=np.int64, count=n)
        neigh_starts = np.zeros(n + 1, dtype=np.int64)
        np.cumsum(counts, out=neigh_starts[1:])
        neigh_idx = np.empty(int(neigh_starts[-1]), dtype=np.int32)
        for i, s in enumerate(self.neighbours):
            block = neigh_idx[neigh_starts[i]:neigh_starts[i + 1]]
            block[:] = sorted(s)
        return neigh_starts, neigh_idx

    # --- Diagnostics ------------------------------------------------

    def is_symmetric(self) -> bool:
        """True iff `j in N(i) <=> i in N(j)` for all i, j."""
        for i, ns in enumerate(self.neighbours):
            for j in ns:
                if i not in self.neighbours[j]:
                    return False
        return True

    def is_border(self, threshold: int = 6) -> np.ndarray:
        """Bool array, True for cells with fewer than `threshold`
        neighbours. Default threshold = 6 = interior count of a
        regular hex lattice (paper main p. 583 fig. 1 caption).
        """
        return np.fromiter(
            (len(s) < threshold for s in self.neighbours),
            dtype=bool,
            count=self.num_cells,
        )
