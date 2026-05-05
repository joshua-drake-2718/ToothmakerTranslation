"""Parameter loading for Path B (silicoshark).

Reads the FORTRAN-compatible parameter file format used by `examples/seal.txt`
(value name, one pair per line, FORTRAN parap order). Maps each FORTRAN
symbol to the paper symbol on the dataclass, so equations downstream can read
`p.k_egr` rather than `p.Egr` and stay close to the 2010 paper's notation.

References:
- 2010 paper, SI Table 1 (parameter inventory):
  docs/research/paper-review-2010-salazar-ciudad-jernvall.md
- 13.f90 setparams() ordering: 13.f90:setparams (preserved in esclec.py)
- 13.f90 parap[6] (Swi) is silently dropped — see
  docs/research/13f90-vs-humppa-divergences.md note 7. We accept Swi here
  but currently route it nowhere; document and revisit if needed.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


# FORTRAN parap-slot ordering of seal.txt-style parameter files. Index = line
# number (1-based) - 1; the value here is the (FORTRAN-name, attr-name) pair
# we expect on that line. Lines marked None are slots present in 13.f90's
# parap[] array but not stored in the file (parap[5]/parap[9]/parap[11]/
# parap[15] in the 0-indexed view).
#
# The list order mirrors esclec.set_params (FORTRAN parap[2..31]).
_PARAM_ORDER: list[tuple[str, str]] = [
    ('Egr',     'k_egr'),     # eq. 5  — epithelial growth
    ('Mgr',     'k_mgr'),     # eq. 12 — mesenchymal growth
    ('Rep',     'k_rep'),     # eq. 1  — repulsion stiffness
    ('unused6', '_unused6'),  # parap[5] in FORTRAN — Swi/tadi border-band
                              # threshold; 13.f90 setparams drops it.
    ('Adh',     'k_adh'),     # eq. 2  — adhesion
    ('Act',     'k_act'),     # eq. 14 — activator autoregulation strength
    ('Inh',     'k_inh'),     # eq. 14 — inhibitor strength on Act
    ('unused10', '_unused10'),
    ('Sec',     'k_sec'),     # eq. 18 — Sec secretion rate at knot/Set cells
    ('unused12', '_unused12'),
    ('Da',      'k_da'),      # eq. 14 — Act diffusion
    ('Di',      'k_di'),      # eq. 15 — Inh diffusion
    ('Ds',      'k_ds'),      # eq. 16 — Sec diffusion
    ('unused16', '_unused16'),
    ('Int',     'k_int'),     # initial differentiation threshold for Inh secretion
    ('Set',     'k_set'),     # secondary differentiation threshold for Sec secretion
    ('Boy',     'k_boy'),     # eq. 13 — buoyancy
    ('Dff',     'k_dff'),     # eq. 6  — differentiation efficiency
    ('Bgr',     'k_bgr'),     # eqs. 7-9 — border (cervical-loop) growth
    ('Abi',     'k_abi'),     # anterior bias multiplier
    ('Pbi',     'k_pbi'),     # posterior bias multiplier
    ('Bbi',     'k_bbi'),     # buccal bias multiplier
    ('Lbi',     'k_lbi'),     # lingual bias multiplier
    ('Rad',     'rad'),       # initial hex radius (integer)
    ('Deg',     'k_deg'),     # eqs. 14-16 — generic degradation
    ('Dgr',     'k_dgr'),     # cervical-loop downgrowth multiplier
    ('Ntr',     'k_ntr'),     # eq. 4 — nucleus traction
    ('Bwi',     'k_bwi'),     # 2014 border-band width — applied in apply_border_bias
    ('ina',     'init_act'),  # FORTRAN initact: initial Act concentration on
                              # all epithelial cells. Paper: "an arbitrary
                              # constant source of Act in the border cells
                              # initially". Path B follows the FORTRAN here
                              # so seal.txt golden is reproducible.
    ('umgr',    'k_umgr'),    # mesenchyme-growth multiplier (under-cell)
]


@dataclass
class Params:
    """Tooth-morphogenesis parameters, paper-symbol-named.

    Construct via `Params.from_file(path)` to load a FORTRAN-compatible file.
    """

    # Mechanical
    k_rep: float = 0.0
    k_adh: float = 0.0
    k_ntr: float = 0.0
    k_egr: float = 0.0
    k_bgr: float = 0.0
    k_mgr: float = 0.0
    k_boy: float = 0.0
    k_dgr: float = 0.0

    # Reaction-diffusion
    k_act: float = 0.0
    k_inh: float = 0.0
    k_sec: float = 0.0
    k_da: float = 0.0
    k_di: float = 0.0
    k_ds: float = 0.0
    k_deg: float = 0.0

    # Differentiation
    k_dff: float = 0.0
    k_int: float = 0.0
    k_set: float = 0.0

    # Border biases
    k_abi: float = 0.0
    k_pbi: float = 0.0
    k_bbi: float = 0.0
    k_lbi: float = 0.0
    k_bwi: float = 0.0
    k_umgr: float = 0.0

    # Geometry
    rad: int = 0
    init_act: float = 0.0

    # Unused parap slots (kept for round-trip readability of files)
    _unused6: float = 0.0
    _unused10: float = 0.0
    _unused12: float = 0.0
    _unused16: float = 0.0

    # Source provenance
    source_path: str = ''

    @classmethod
    def from_file(cls, path: str | Path) -> 'Params':
        path = Path(path)
        out = cls(source_path=str(path))
        with path.open() as f:
            lines = [ln for ln in f.readlines() if ln.strip()]
        if len(lines) < len(_PARAM_ORDER):
            raise ValueError(
                f'{path}: expected at least {len(_PARAM_ORDER)} non-blank '
                f'lines, got {len(lines)}'
            )
        for i, (fortran_name, attr) in enumerate(_PARAM_ORDER):
            parts = lines[i].split()
            if len(parts) < 2:
                raise ValueError(
                    f'{path}:{i+1}: expected "value name", got {lines[i]!r}'
                )
            value = float(parts[0])
            if attr == 'rad':
                value = int(value)
            setattr(out, attr, value)
        return out
