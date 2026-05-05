---
title: 'Path B v2: cusp-forming disentanglement and Act sweep (B.1 strengthening)'
date: 2026-05-05
related-rebuttal: docs/paper-rebuttal-preemption-2026-05-05.md §B.1, §B.4
related-briefing: docs/paper-briefing-2026-05-05.md §5B, §5C
---

## Question

Two questions for the rebuttal §B.1 strengthening pass:

1. Does the single-field disentanglement, when run on the cusp-forming
   dataset (`examples/wt-tribosphenic-2014.txt`), reveal a load-bearing
   structure consistent with — or different from — the seal example's
   `rep_form` / `adh_form` finding?
2. Does the eq. 14 typo divergence (PAPER_2010 = 19 cusps,
   PAPER_LITERAL_2010 = 0) survive a sweep of `k_act` across the full
   published range Ext Data Fig. 4a of Harjunmaa et al. specifies
   (Act ∈ [0.1, 1.6] in 0.1 increments)?

## Method

**Disentanglement** — `scripts/run-discretisation-study.py
--params examples/wt-tribosphenic-2014.txt --single-field-mode both
--iters 500 --saves 5`. Output:
`experiments/discretisation-study/single-field-cusp-forming/`. 28 rows
(2 anchors + 14 knock-down + 14 knock-up). 2,500 iterations per row to
match the existing committed cusp-forming evidence cadence.

**Act sweep** — 16 parameter files
(`experiments/discretisation-study/act-sweep-2014/params/Act_{value}.txt`)
identical to `examples/wt-tribosphenic-2014.txt` except for the `Act`
parameter. All 5 named presets × 16 Act values = 80 runs at 2,500
iterations. Output:
`experiments/discretisation-study/act-sweep-2014/results/`. Aggregator:
`scripts/analyse_act_sweep.py` (writes `cusp-vs-act.csv` and
`docs/figures/fig5-act-sweep/`).

## Result A: cusp-forming disentanglement

| Anchor | Cells | Cusps | Knock direction |
|---|---:|---:|---|
| `LEGACY_FORTRAN` | 19 | 7 | knock-down direction |
| `PATH_B_DEFAULT` | 19 | 19 | knock-up direction |

Of the 14 fields that differ between the two anchors:

- `knot_threshold_gate`: the single field where the choice moves the
  cusp count between the two regimes.
  - `LEGACY_FORTRAN_minus_knot_threshold_gate`: 7 → 19 cusps
    (knock-down to PATH_B_DEFAULT's `'none'`).
  - `PATH_B_DEFAULT_plus_knot_threshold_gate`: 19 → 7 cusps
    (knock-up to LEGACY_FORTRAN's `'first_border_cell'`).
  - The field is symmetric: it carries 100 % of the cusp-count span on
    both directions.
- `update_order`: dominant for cell count.
  - `LEGACY_FORTRAN_minus_update_order`: cell count 19 → 60 (knocking
    out the FORTRAN's Gauss–Seidel order in favour of Jacobi triggers
    division beyond the 19-cell plateau).
  - Cusp count unchanged at 7.
- `rep_form`: runaway-on-knock-up.
  - `PATH_B_DEFAULT_plus_rep_form`: NaN (timed out at 600 s; cell
    division runaway). Same failure mode as on the seal example,
    confirming the field's instability is not a parameter-set artefact.
- 11 other fields: dormant in both directions (no observable change in
  cell count or cusp count).

The auto-generated `results-table.md` interpretation reads each row as
'not consequential' because its 'percent of span' formula divides by
the 0.0-cell span between anchors (both anchors have plateau 19); the
auto-text is misleading on this dataset and the cusp-count signal must
be read off the table directly.

## Result A — methodological reading

The shape of the load-bearing structure on the cusp-forming dataset is
similar to the seal example — a small set of dominant fields plus a
long tail of dormant ones — but the *identity* of the dominant fields
is parameter-set-specific:

| Dataset | Dominant for cusps | Dominant for cells | Runaway-on-knock-up | Dormant count |
|---|---|---|---|---|
| seal.txt | not exercised (no cusp formation) | `rep_form` (100 %), `adh_form` (87 %) | `rep_form`, `division_total_cap`, `knot_threshold_gate` | 11 |
| wt-tribosphenic-2014.txt | `knot_threshold_gate` (100 %) | `update_order` (delta 41 cells) | `rep_form` | 11 |

`rep_form` carries instability across both datasets. The other dominant
fields differ. The dormant count is the same in both cases — eleven of
the fourteen differing fields are biologically inert under the chosen
metric and parameter set. This is the methodologically interesting
result: the catalogue's *structure* (small dominant subset, long
dormant tail) is parameter-set-robust, but its *contents* are
parameter-set-specific. Without the catalogue and the perturbation
matrix, neither the structural regularity nor the parametric
specificity would be visible.

The eq. 14 typo finding (PAPER_2010 = 19 cusps, PAPER_LITERAL_2010 = 0)
is unaffected by this disentanglement — the typo divergence is between
two paper-flavoured presets, not between LEGACY_FORTRAN and
PATH_B_DEFAULT. The disentanglement quantifies what *else* is going on
inside the FORTRAN bundle on this dataset; it does not bear on the
typo finding directly.

## Result B: Act sweep

| Preset | Cusps at every Act ∈ [0.1, 1.6] |
|---|---|
| `HUMPPA_LITERAL` | 7 |
| `LEGACY_FORTRAN` | 7 |
| `PAPER_2010` | 19 |
| `PAPER_LITERAL_2010` | 0 |
| `PATH_B_DEFAULT` | 19 |

The cusp count is **flat across the full published Act range** for
every preset. Cell counts are uniformly 19. Vertex envelopes match
within floating-point noise. See
`experiments/discretisation-study/act-sweep-2014/cusp-vs-act.csv` and
`docs/figures/fig5-act-sweep/fig5-act-sweep.{pdf,png}` for the full
data and the cusp-vs-Act plot.

## Result B — methodological reading

Two distinct findings come out of the Act sweep, and the paper should
keep them separate.

**Finding B1: the eq. 14 typo divergence is robust across the full
published Act range.** PAPER_2010 produces 19 cusps and
PAPER_LITERAL_2010 produces 0 at every Act value Harjunmaa et al.'s
Ext Data Fig. 4a sweeps over. This is a much stronger statement than
the single-Act-point comparison the briefing originally reported: the
typo's catastrophic biology is invariant under the published
activator-strength range. Combined with the ina = 0 14,000-iteration
run (which confirmed mathematical dormancy in the Act = 0 limit), the
typo finding now rests on 17 parameter-set variants of cusp-forming
plus the analytical Act = 0 reduction.

**Finding B2: the silicoshark replication is in a saturated regime
where Act does not move biology.** Harjunmaa et al.'s Ext Data Fig. 4a
shows a monotonic 1 → 5 cusp progression as Act sweeps from 0.1 to
1.6; silicoshark's `mesenchyme='absent'` / `topology='static_with_local_update'`
configuration produces no such progression on any preset. The
FORTRAN-flavoured presets cap at 7 cusps (the
`knot_threshold_gate='first_border_cell'` allowance, which restricts
knot eligibility to a fixed 7-cell band irrespective of Act); the
paper-flavoured presets cap at 19 cusps (every cell crosses the knot
threshold and the absent-mesenchyme configuration provides no
inhibitory bath to suppress neighbour-cell knot formation); the
typo'd preset stays at 0 across the entire range.

The silicoshark configuration is therefore not a quantitative
replication of the FORTRAN's published Act-sweep curve. Reproducing
the curve quantitatively requires implementing
`mesenchyme='per_column_z_layers'` with parameter values that allow
inhibitor accumulation between knot loci, and/or implementing the
`fortran_margins` laplacian (§Future work item 3 of the briefing).
Until that work is done, silicoshark can answer 'does the typo
matter?' but not 'how many cusps does the FORTRAN form at Act = X?'.

These two questions are methodologically separate and the rebuttal
should treat them as such.

## Rebuttal status updates

**B.1** ('two parameter sets is not enough'): WEAK → ARGUABLE. The
strengthening: 17 parameter-set variants of cusp-forming all
producing the same 19-vs-0 cusp divergence under the typo / corrected
denominator pair, plus the disentanglement on a second dataset
showing the catalogue's structural regularity. Residual: all 17
variants are perturbations of one biological problem (the 2014
wild-type tribosphenic mouse), so 'breadth' in the strict sense is
still bounded. A reviewer can still legitimately request a second,
independent biological problem (e.g. a different species' parameter
set, a knock-out variant); §Future work names that as the next step.

**B.4** ('the rep_form result might be parameter-tuning artefact'):
ARGUABLE → ARGUABLE (no change). The disentanglement on cusp-forming
shows that `rep_form` carries its instability signature across two
parameter sets (NaN runaway on PATH_B_DEFAULT in both cases), which
is reassuring; but the seal-example's quantitative claim
(`rep_form` carries 100 % of cell-count span) does not test on
cusp-forming because both anchors' cell-count plateaus are 19.
Still ARGUABLE because the dataset-specific quantitative claim is
the part B.4 questioned.

## Files added

- `experiments/discretisation-study/single-field-cusp-forming/`
  (28-row results, OFFs, metrics)
- `experiments/discretisation-study/act-sweep-2014/`
  (16 parameter files + 80-row results + cusp-vs-act.csv)
- `docs/figures/fig5-act-sweep/fig5-act-sweep.{pdf,png}`
- `scripts/analyse_act_sweep.py` (figure generator)

## Bug fix

The disentanglement uncovered a CLI bug: `silicoshark/__main__.py`'s
`_coerce` mapped the literal string `'none'` to Python `None`, which
collided with `Discretisation.knot_threshold_gate`'s legitimate
string-valued option `'none'`. The first disentanglement run reported
`unknown knot_threshold_gate: None` for the
`LEGACY_FORTRAN_minus_knot_threshold_gate` row. The fix changes the
None sentinel from `'none'` to `'null'` (JSON convention); the runner
in `scripts/run-discretisation-study.py:_build_subprocess_args` now
emits `'null'` when an override value is Python None. New unit tests
in `tests/test_silicoshark_cli_coerce.py` pin the round-trip behaviour.
