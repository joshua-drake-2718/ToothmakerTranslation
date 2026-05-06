"""OFF-format writer matching the FORTRAN binary's output structure.

The file format used by the FORTRAN binary (and reproduced by `esclec.py`):

    COFF
    V F 0
    x y z r g b a    (V lines)
    n i1 i2 ... iN   (F lines, n = vertex count of face)

Vertices are epithelial cell centres (mesenchymal cells are not exported).
Faces are triangles or quads found by walking the neighbour graph.

The reference Python writer (esclec.guardaveinsoff_2) walks the neighbour
table for triangle/quad detection. Path B uses scipy.spatial.Delaunay and
takes the resulting triangle list directly — no quads are emitted (the
FORTRAN's quad emission corresponds to mesh defects in the topology walk
which we avoid by re-triangulating). This means face counts will not
match byte-for-byte, but the cell-count plateau and vertex envelope
(the validation criteria from the charter) are the comparable metrics.
"""
from __future__ import annotations

from typing import IO

import numpy as np

from .state import State


def cell_colour(state: State, i: int) -> tuple[float, float, float, float]:
    """RGB+A colour for cell i, mirroring esclec.get_rainbow_knot semantics.

    Knot cells -> cyan (0, 1, 1, 0). Differentiating cells span a rainbow
    based on differentiation state. The FORTRAN's `mat()` builds an `ma`
    array (1.0 for knots, 1.0 if diff > Set, 0.1 if diff > Int else 0.0)
    and rainbow-maps that. We reproduce the same pattern here.
    """
    if state.knot[i]:
        return (0.0, 1.0, 1.0, 0.0)
    # Without knowledge of Int/Set thresholds in this writer we use a
    # plain rainbow over the differentiation state. The colour is purely
    # cosmetic and not part of validation.
    f = float(state.diff[i])
    if f < 0.07:
        return (0.6, 0.6, 0.6, 0.8)
    if f < 0.2:
        return (1.0, f, 0.0, 0.5)
    if f < 1.0:
        return (1.0, f * 3, 0.0, 1.0)
    return (1.0, 1.0, 0.0, 1.0)


def write_off(
    state: State,
    triangles: np.ndarray,
    f: IO[str],
) -> None:
    """Emit a COFF file for `state` with face list `triangles`.

    Parameters
    ----------
    state : State
    triangles : (T, 3) int array of vertex indices, or empty array.
    f : open writable text file.
    """
    n_v = state.num_active
    n_f = int(triangles.shape[0]) if triangles.size else 0

    f.write('COFF\n')
    f.write(f'{n_v} {n_f} 0\n')
    for i in range(n_v):
        x, y, z = state.positions[i]
        r, g, b, a = cell_colour(state, i)
        f.write(f'{x} {y} {z} {r} {g} {b} {a}\n')
    for tri in triangles:
        f.write(f'3 {int(tri[0])} {int(tri[1])} {int(tri[2])}\n')
