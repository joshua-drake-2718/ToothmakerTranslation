---
title: 'Parameter set: 2014 wild-type tribosphenic mouse (Harjunmaa et al.)'
date: 2026-05-05
---

## Source

`docs/research/paper-review-2014-harjunmaa.md`, table under
'## New parameters introduced (or first published values)'.

The values come from Harjunmaa et al. 2014, *Nature*, Extended Data
Fig. 2 — the ToothMaker GUI screenshot for the 'Tribosphenic tooth'
model variant set to wild-type mouse.

## Expected behaviour

Per Ext Data Fig. 4a, the wild-type set with `Act = 1.6` asymptotes
to ~5 cusps as iterations progress. Ext Data Table 1 reports an
in-vivo wild-type cusp count of ~5–6 with talonid height ~1.0. The
silicoshark replication of this set under `LEGACY_FORTRAN` is
expected to form at least two knot cells (Act >= 1) within 500
iterations and to accumulate non-trivial inhibitor on the
epithelium (max[Inh] > 0.1 at some point in the run).

## Deviations from the literal Ext Data Fig. 2 values

`ina = 0.5` (the paper caption notes 'Ina is not used in the model'
for wild-type, i.e. `ina = 0.0`). The paper's mechanism for seeding
cusp formation under `ina = 0` relies on the lingual/buccal border
floors (`Bbi = 1.34`, `Lbi = 1.31`) lifting border-cell `[Act]` above
threshold over the published 14,000-iteration sweep. silicoshark's
`mesenchyme='absent'` workaround does not exercise those floors as
strongly within 500 iterations under the static-with-local-update
topology, so I seed Act uniformly with `ina = 0.5` to bring
threshold-crossing forward into the 500-iter window. This is the
same convention `examples/seal.txt` uses (`ina = 0.15`) and is the
smallest deviation from the paper that produces non-trivial Inh
accumulation under the silicoshark code path. The 2010 paper itself
uses uniform `ina` in its example simulations.

`umgr = 0.0`. Not in the 2014 paper at all. A post-2014 addition to
`13.f90` with no GUI control in Ext Data Fig. 2. Set to zero to
match its absence in the published parameter set.

## File format

FORTRAN `parap`-slot order, one `(value name)` pair per line. The
silicoshark loader at `silicoshark/params.py:_PARAM_ORDER` does not
support `#`-prefixed comment lines; this sidecar carries the
provenance instead.

## Validation status

- Knot formation: confirmed under `LEGACY_FORTRAN` within 500
  iterations (see `experiments/discretisation-study/cusp-forming/`).
- Non-trivial Inh: confirmed (see same).
