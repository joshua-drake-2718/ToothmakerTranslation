"""silicoshark — paper-faithful re-implementation of the tooth-morphogenesis
model from Salazar-Ciudad and Jernvall, *Nature* 464, 583 (2010), extended
with the tribosphenic parameters from Harjunmaa et al., *Nature* 512, 44
(2014).

This is Path B in the project plan. It is intentionally separate from the
FORTRAN-faithful translation in `coreop2d.py`/`esclec.py`/`main.py`, which
remains the cross-validation oracle.

See `docs/plans/2026-05-04-path-b-charter.md` for scope, paper-vs-FORTRAN
precedence rules, and validation strategy.
"""
from .params import Params
from .state import State, build_initial_state, hex_lattice
from .mesh import Mesh
from .io import write_off
from .simulator import run, step

__all__ = [
    'Params',
    'State',
    'Mesh',
    'build_initial_state',
    'hex_lattice',
    'run',
    'step',
    'write_off',
]
