---
title: 'Path B v2 charter — configurable discretisation as methodological contribution'
author: Lyndon Drake (with Claude Code)
date: 2026-05-05
supersedes: docs/plans/2026-05-04-path-b-charter.md (relevant decision points; v1 charter remains the historical record)
---

## Aim and methodological framing

I aim to do two things at once.

The first is to re-implement the tooth-morphogenesis model from
Salazar-Ciudad and Jernvall, *Nature* 464, 583 (2010), and the
tribosphenic extension in Harjunmaa et al., *Nature* 512, 44 (2014),
in idiomatic vectorised numpy. This was the v1 aim and remains in
scope. The FORTRAN-faithful Path A in `coreop2d.py` is the
cross-validation oracle.

The second is to use this re-implementation as a vehicle for a
methodological study. Path B v1 surfaced a striking finding: where
the paper's prose under-specifies the model, the FORTRAN's choice is
neither uniquely-correct nor obviously-wrong. Many of the FORTRAN's
choices are defensible numerical decisions inside an ambiguity
window left open by the paper. A few are paper-vs-code tensions
where one side is more credibly 'the intended model'. A handful
are FORTRAN accidents (load-bearing or otherwise). The v1 work
showed that LLM-assisted analysis of the paper, the FORTRAN
source, and the broader code-base lineage (humppa, C++ port,
tgrohens, silicoshark) can identify, classify, and justify each
of these choices in a way that human-only review struggles with at
this depth.

I argue that this matters beyond developmental modelling. The
problem of knowing whether a piece of code matches its informal
human-language description is a long-standing one in software
engineering. The formal-specification literature has tools to
relate code to mathematically-precise specifications, but those
tools are not usable for most real problems, including most of
science. The pragmatic gap between informal paper-prose and
executable code is bridged today, in practice, by trust, by
expert reading, and by the implementer's intuition. AI-assisted
analysis offers a third option: a partial bridge that does not
formalise the specification but does make the gap *visible and
auditable*. This work is a case study in that proposition, with
the tooth-morphogenesis model as a non-trivial worked example.

The methodological contribution is the disciplined enumeration
of implementer-choice points, each tied to a paper passage, a
FORTRAN line, and a numpy-idiomatic alternative. The scientific
contribution is what the comparison study reveals about which
choices materially change biological output (cusp count, cusp
shape, regime identification) and which are numerically
irrelevant.

This double aim is reflected in two outputs:

1. The `silicoshark/` package, configurable across each named
   choice point, with named presets (`LEGACY_FORTRAN`,
   `HUMPPA_LITERAL`, `PAPER_2010`, `PAPER_LITERAL_2010`,
   `PATH_B_DEFAULT`).
2. A comparison study running the presets across the seal example
   plus the 2014 paper's qualitative validation targets (Ext Data
   Figs. 2, 4, 5), reporting cusp count, cusp width, vertex
   envelope, cell-count plateau, and qualitative regime.

A follow-up paper draws on (2) to argue (a) which choices in the
2010 model are load-bearing, (b) which are noise, and (c) what the
methodology of LLM-assisted comparative analysis adds beyond
human-only paper-vs-code review. The seal-golden replication is
no longer the validation target — it is one cell in the comparison
matrix, occupied by the `LEGACY_FORTRAN` preset.

## Scope (unchanged from v1)

Two papers, in priority order:

1. Salazar-Ciudad and Jernvall, *Nature* 464, 583 (2010). DOI
   [10.1038/nature08838](https://doi.org/10.1038/nature08838).
   Equation source of truth where unambiguous.
2. Harjunmaa et al., *Nature* 512, 44 (2014). DOI
   [10.1038/nature13613](https://doi.org/10.1038/nature13613).
   Tribosphenic extension and parameter set.

Out of scope: Zimm et al., *PNAS* 120(15) e2216959120 (2023). The
shark-specific extensions are not modelled. Reading the SI for
parameter naming conventions is fine.

## Architecture

### `Discretisation` dataclass

A frozen dataclass carrying every implementer-choice point. Travels
alongside `Params` through every simulator function. Each field is a
`Literal[...]` of named options, with a default that is the most
defensible numpy-idiomatic choice. Module-level presets construct
canonical configurations.

```python
@dataclass(frozen=True)
class Discretisation:
    # Spatial discretisation
    laplacian: Literal['length_weighted', 'cotangent', 'fortran_margins'] = 'length_weighted'

    # Temporal discretisation
    update_order: Literal['jacobi', 'gauss_seidel_forces', 'mixed_jacobi_gs'] = 'jacobi'

    # Mesh / topology
    topology: Literal['delaunay_each_step', 'static_with_local_update'] = 'static_with_local_update'
    division_per_step_cap: int | None = None
    division_max_edge: float | None = None

    # Eq. 5 epithelial growth
    eq5_z_gate: bool = True
    eq5_apply_to: Literal['all', 'interior_only'] = 'all'

    # Eqs. 1-2 cell-cell mechanics
    rep_form: Literal['hookean_signed', 'paper_gated'] = 'hookean_signed'
    adh_form: Literal['unit_vector', 'hookean_attraction'] = 'unit_vector'
    rep_neighbour_set: Literal['mesh', 'mesh_plus_all_close'] = 'mesh'

    # Reaction terms
    eq14_denominator: Literal['act_typo', 'inh_corrected'] = 'inh_corrected'
    eq17_inh_source: Literal['act_concentration', 'act_rate_times_di'] = 'act_concentration'
    eq18_sec_source: Literal['constant_k_sec', 'k_sec_times_di'] = 'constant_k_sec'

    # Differentiation
    diff_accumulator: Literal['act', 'sec'] = 'sec'
    knot_threshold_gate: Literal['none', 'first_border_cell'] = 'none'
    knot_daughter_di: Literal['zero_reset', 'inherit_avg'] = 'zero_reset'

    # Mesenchyme
    mesenchyme: Literal['absent', 'per_column_z_layers'] = 'per_column_z_layers'
    n_mesenchyme_layers: int = 2

    # Border identification and biases
    border_definition: Literal['neighbour_count', 'topological_descendants'] = 'neighbour_count'
    border_bias_x_zero_quirk: bool = True
```

### Decision-point branching

Each branch is a single short conditional at the call site. No
abstract base classes; no strategy pattern. The branching is
deliberately visible and grep-able:

```python
def step_reaction_diffusion(state, params, mesh, dt, disc):
    ...
    if disc.eq14_denominator == 'inh_corrected':
        rd_act = k_act * act / (1.0 + k_inh * inh) - k_deg * act
    elif disc.eq14_denominator == 'act_typo':
        rd_act = k_act * act / (1.0 + k_inh * act) - k_deg * act
    ...
```

A reader scanning the function sees every option exercised and the
charter row that justifies it (via a sibling docstring referencing
`docs/research/discretisation-audit.md`).

### Named presets

| Preset | Aim | Key non-default fields |
|---|---|---|
| `LEGACY_FORTRAN` | Reproduce 13.f90 byte-for-byte semantically, including its accidents | `laplacian='fortran_margins'`, `update_order='gauss_seidel_forces'`, `topology='static_with_local_update'`, `eq5_z_gate=True`, `eq5_apply_to='interior_only'`, `rep_form='paper_gated'`, `adh_form='hookean_attraction'`, `rep_neighbour_set='mesh_plus_all_close'`, `eq17_inh_source='act_rate_times_di'`, `eq18_sec_source='k_sec_times_di'`, `diff_accumulator='sec'`, `knot_threshold_gate='first_border_cell'`, `knot_daughter_di='inherit_avg'`, `border_bias_x_zero_quirk=True` |
| `HUMPPA_LITERAL` | The upstream Catalan FORTRAN (`humppa_translate.f90`) — distinguishes 13.f90's renaming-introduced divergences from genuine model behaviour | Same as `LEGACY_FORTRAN` minus the txungu/Swi/parap[6] quirks |
| `PAPER_2010` | The 2010 paper as written, with the eq. 14 typo corrected and choices made by best-faith reading where the paper is silent | All defaults |
| `PAPER_LITERAL_2010` | The 2010 paper as written, with the eq. 14 typo *not* corrected — for testing whether the typo materially changes biology | `eq14_denominator='act_typo'` |
| `PATH_B_DEFAULT` | The cleanest numpy-idiomatic choices for production use; identical to `PAPER_2010` initially, may diverge as we learn | All defaults |

Each preset has a docstring naming what the comparison study uses
it to test for.

## Validation strategy

Replaces the v1 single-tolerance approach with a comparison
matrix.

### Replication anchor

The `LEGACY_FORTRAN` preset must reproduce the FORTRAN binary's
output on `examples/seal.txt` to within the v1 tolerance:

- Cell count plateau within ±5% of 57.
- Vertex envelope within 10% of `tests/golden_fortran/500_run.off`.
- No NaN/Inf.
- No RuntimeWarning.

If `LEGACY_FORTRAN` cannot meet this, the presets-and-options
decomposition has a bug that needs fixing before the comparison
study is meaningful.

### Comparison matrix

The 2014 paper provides four qualitative validation targets that
exercise the model across parameter space (charter v1 §Phase 2):

1. Wild-type tribosphenic mouse (Ext Data Fig. 2, 25 parameters).
2. Act sweep (Ext Data Fig. 4, Act ∈ [0.1, 1.6] in 0.1 increments).
3. Inh / Di perturbation lattice (Ext Data Fig. 5, 5 (Inh, Di) points).
4. Parameter sweep map (Ext Data Fig. 3, 9 parameter directions).

For each preset × parameter file combination, capture:

- **Cusp count** (topological — number of `[Act] >= 1` knot
  cells at the plateau).
- **Cusp width** (geometric — paper's Fig. 5 metric).
- **Cell-count plateau** (final cell count).
- **Vertex envelope** (min/max x, y, z).
- **Plateau trajectory** (cell-count vs iteration).
- **Qualitative regime** (monotone-grow / oscillate / collapse / NaN).

Run `LEGACY_FORTRAN` on the same inputs as the cross-reference. The
comparison matrix is the deliverable for the methodological paper.

### Audit-trail validation

Every option of every Discretisation field must be traceable to:

(a) A passage in the 2010 or 2014 paper that frames the question
    (or explicitly notes its absence).
(b) A line in `13.f90`, `humppa_translate.f90`, or another
    code-base in the lineage that resolves the question one way.
(c) Any third source (C++ port, tgrohens, silicoshark FORTRAN)
    that resolves it differently, if it exists.

`docs/research/discretisation-audit.md` is the canonical record.
A pre-commit check (or just visual review) ensures no option is
added without a citation.

## Implementer-choice catalogue

Each row below corresponds to one Discretisation field. The full
audit lives in `docs/research/discretisation-audit.md`; this charter
section names the question, sketches the options, and points at
where the paper / FORTRAN evidence lives.

### Summary table — kind of disagreement per field

The 'kind' column classifies what each option pair (or trio) really
is. There are three categories:

- **PaperAmbiguity**: the paper is silent or sketches the answer in
  prose without pinning a discrete formula. Multiple options are
  defensible readings; none is uniquely right; this is exactly the
  ambiguity the methodological study is meant to expose.
- **PaperVsCodeTension**: the paper says one thing and the FORTRAN
  does another, and one side is more credibly the intended model.
  These are the rows where v1's charter took a position; the
  presets `PAPER_*` and `LEGACY_FORTRAN` exercise both sides.
- **FortranAccident**: a 13.f90-only artefact (renaming bug, fp
  fluke, dead branch) that is preserved in `LEGACY_FORTRAN` only
  for golden compatibility, not because it represents a real
  modelling choice.

| Field | Kind | Paper says | FORTRAN does | Notes |
|---|---|---|---|---|
| `laplacian` | PaperAmbiguity | 'finite-volume on triangular mesh, flux ∝ contact area' (no formula) | per-edge `pes` margins with `area_p` z-weighting | Cotangent and length-weighted are textbook alternatives |
| `update_order` | PaperVsCodeTension | Prose strongly suggests Jacobi | Mixed: reactions Jacobi, mechanical Gauss–Seidel | Paper-prose is reasonable but ambiguous |
| `topology` | PaperAmbiguity | SI fig. 1 picture only — no algorithm | Static graph + local update via topology walk | Delaunay re-triangulation is a v1 numpy-idiomatic alternative; v1 evidence shows it is unstable on the seal example |
| `division_per_step_cap` | PaperAmbiguity | Silent | No cap (all qualifying edges split) | Safeguard knob for `delaunay_each_step`; inert under static |
| `division_max_edge` | PaperAmbiguity | Silent | No cap | Same as above |
| `eq5_z_gate` | PaperAmbiguity | 'outward unit vector' along apical direction | Gates on `b < -1e-4` (neighbours strictly above in z) | Gate is FORTRAN's reading of 'outward' for a curved sheet |
| `eq5_apply_to` | PaperVsCodeTension | All cells (no restriction) | Inner cells only (`i >= first_border_cell`) | The paper formula has no border / interior split |
| `rep_form` | PaperAmbiguity | Eq. 1 written as Hookean spring; prose treats rep + adh as separate | Hookean (always-on) when shifted by sign | Both are defensible readings of eq. 1 |
| `adh_form` | PaperVsCodeTension | Unit-vector pull (eq. 2) | Hookean (`Adh * d_vec`, not unit-normed) | FORTRAN form is stronger at distance |
| `rep_neighbour_set` | PaperVsCodeTension | Mesh neighbours only | Mesh + extra `pushingnovei` loop over close non-mesh cells | Paper formula has no extra loop; FORTRAN's is plausibly load-bearing for stability |
| `eq14_denominator` | PaperVsCodeTension | `1 + k_inh * [Act]` (typo) | `1 + Inh * [Inh]` (biologically right) | Charter follows FORTRAN; PAPER_LITERAL_2010 tests the typo's biological cost |
| `eq17_inh_source` | PaperVsCodeTension | `[Act]` (concentration) | `(rate of [Act]) * d_i` (uses temp variable) | FORTRAN form is suspicious as a faithful translation |
| `eq18_sec_source` | PaperVsCodeTension | `k_sec` (constant) | `k_sec * d_i` (smooth ramp) | FORTRAN's `d_i` modulation is not in paper |
| `diff_accumulator` | PaperVsCodeTension | Eq. 6: `[Act]` | `[Sec]` (humppa line 659) | v1 charter followed FORTRAN; this preset matrix is what justifies revisiting |
| `knot_threshold_gate` | PaperVsCodeTension | Any cell with `[Act] >= 1` | Only cells with `i >= first_border_cell` | FORTRAN excludes innermost cells from knot formation; biologically odd |
| `knot_daughter_di` | PaperAmbiguity | 'New cell is not a knot cell' (silent on `d_i`) | Inherit avg of mothers' `d_i` | Paper-aligned reading is to reset to 0 |
| `mesenchyme` | PaperVsCodeTension | 'Two layers of mesenchymal cells underneath' | 4 z-layers per epithelial column | Path B v1 used `absent` as a shortcut; not paper-faithful |
| `n_mesenchyme_layers` | PaperVsCodeTension | 2 | 2 + substrate-edge boundary layer | FORTRAN's 4th layer is a boundary condition, not a real layer |
| `border_definition` | PaperAmbiguity | Conflates topological (descendants of initial border) and geometric (low neighbour count) | Mostly topological (`first_border_cell`) with some geometric promotion via `update_border_cells` | Paper conflates; FORTRAN is also mixed |
| `border_bias_x_zero_quirk` | FortranAccident | n/a (paper does not describe the multiplier scheme at this level of detail) | Cells with x exactly == 0 don't receive Pbi/Abi/Bgr (because of `if x>0` / `elif x<0` chain) | Documented as load-bearing for the seal 57-cell plateau in coreop2d.py docstring |

The classification is itself an output of LLM-assisted analysis,
informed by reading the 2010 paper, the 2014 paper, all four
FORTRAN code-bases (humppa, 13.f90, tgrohens, silicoshark) and the
C++ port. A goal of the comparison study is to test how well this
classification predicts which fields actually move biology.

### Spatial discretisation

**`laplacian`** — How is `nabla^2 c` discretised on the triangular
mesh? Paper says 'finite-volume with flux ∝ contact area', no
formula. Options: `length_weighted` (1/|edge_ij|, simplest),
`cotangent` (textbook discrete Laplace–Beltrami), `fortran_margins`
(reproduces 13.f90's `apply_diffusion`). Citations: paper main p.
585; coreop2d.py:506 (apply_diffusion).

### Temporal discretisation

**`update_order`** — Are forces and reactions computed against the
start-of-step state ('Jacobi') or against in-place updated state
('Gauss–Seidel')? Paper prose suggests Jacobi; FORTRAN is mixed
(reactions Jacobi, mechanical forces Gauss–Seidel). Options:
`jacobi`, `gauss_seidel_forces`, `mixed_jacobi_gs`. Citations:
paper main p. 583; coreop2d.iteration().

### Mesh and topology

**`topology`** — Do we re-triangulate each step, or maintain a
static graph updated locally on division? Paper does not address
this. v1 evidence (`docs/findings/2026-05-04-path-b-v1-mesh-degeneration.md`)
shows Delaunay-each-step is unstable on the seal example. Options:
`delaunay_each_step`, `static_with_local_update`. Citations:
v1 finding; FORTRAN `add_cell` for the static-with-local pattern.

**`division_per_step_cap`** and **`division_max_edge`** — Per-step
division-rate safeguards for the Delaunay-each-step path.

### Mechanical forces

**`eq5_z_gate`** — Restrict the eq. 5 `sum(unit_ij)` to neighbours
strictly above in z (`b < -1e-4`)? Paper says 'outward unit vector'
without specifying. FORTRAN gates. Options: `True`, `False`.

**`eq5_apply_to`** — Apply eq. 5 to all cells (paper) or interior
only (FORTRAN)? Options: `'all'`, `'interior_only'`.

**`rep_form`** — Eq. 1 always-on Hookean (smoothly signed) or paper-
gated (active only when `|d| < p_0`)? Both are defensible readings of
the paper's eq. 1. Options: `hookean_signed`, `paper_gated`.

**`adh_form`** — Eq. 2 unit-vector pull (paper) or Hookean attraction
(FORTRAN's `repulse_neighbour` does this)? Options: `unit_vector`,
`hookean_attraction`.

**`rep_neighbour_set`** — Repulsion only over mesh neighbours, or
also against non-mesh cells too close? FORTRAN has both
(`pushing` + `pushingnovei`); paper has only the first. Options:
`mesh`, `mesh_plus_all_close`.

### Reaction terms

**`eq14_denominator`** — Paper writes `1 + k_Inh [Act]`; FORTRAN
computes `1 + Inh * [Inh]`. The latter is biologically meaningful;
the former is the paper's text. Options: `inh_corrected`, `act_typo`.

**`eq17_inh_source`** — Paper Eq. 17: `[Inh]' = [Act] - Deg [Inh]`.
FORTRAN for `d_i > Int`: `[Inh]' = (rate of Act) * d_i - Deg [Inh]`.
The FORTRAN form is suspicious — it uses a temporary rate variable
where the paper uses the concentration. Options: `act_concentration`
(paper), `act_rate_times_di` (FORTRAN).

**`eq18_sec_source`** — Paper: `[Sec]' = k_sec - Deg [Sec]`. FORTRAN
for `d_i > Set`: `[Sec]' = k_sec * d_i - Deg [Sec]`. The d_i factor
makes the production smoothly ramp up; not in paper. Options:
`constant_k_sec` (paper), `k_sec_times_di` (FORTRAN).

### Differentiation

**`diff_accumulator`** — Paper Eq. 6: `d_i' = k_dff [Act]`. FORTRAN:
`d_i' = k_dff [Sec]`. Path B v1 followed FORTRAN. Options:
`act` (paper), `sec` (FORTRAN). Citations:
`docs/findings/2026-05-04-differentiation-uses-sec-not-act.md`.

**`knot_threshold_gate`** — Paper: any cell with `[Act] >= 1` becomes
a knot. FORTRAN: only cells with `i >= first_border_cell` (the inner
2/3 of the lattice). Options: `none`, `first_border_cell`.

**`knot_daughter_di`** — Paper: 'the new cell is not a knot cell'
when the mother is. Path B v1 reset d_i to 0; FORTRAN inherits the
average. Options: `zero_reset` (paper-aligned), `inherit_avg`
(FORTRAN).

### Mesenchyme

**`mesenchyme`** — Model the underlying mesenchymal layer
explicitly (FORTRAN-style, q3d z-layers) or absent (Path B v1
shortcut)? Paper: 'two layers of mesenchymal cells underneath'.
Options: `absent`, `per_column_z_layers`.

**`n_mesenchyme_layers`** — How many z-layers? FORTRAN: 4
(epithelium + 2 mesenchyme + substrate edge). Paper: 2.

### Border and biases

**`border_definition`** — Topological (descendants of initial border
cells) or geometric (neighbour-count threshold)? Paper conflates
both. FORTRAN uses topological (via `first_border_cell`). Options:
`neighbour_count`, `topological_descendants`.

**`border_bias_x_zero_quirk`** — Preserve FORTRAN's
update_cell_position behaviour where cells with x exactly == 0 don't
receive Pbi/Abi/Bgr multipliers? This is a coreop2d.py:1011
docstring-documented FORTRAN accident. Options: `True`, `False`.

## Out of scope for v2

- Modifying `coreop2d.py`, `esclec.py`, `main.py`, `vector.py`, or
  `13.f90`. Path A remains the cross-validation oracle.
- Modifying `tests/test_simulator_smoke.py` or
  `tests/golden_fortran/*.off`.
- Implementing the Zimm 2023 shark variant.
- Reimplementing the GUI or pre-2010 model versions.

## Deliverables

1. `silicoshark/` package with the `Discretisation` dataclass, named
   presets, and call-site branching.
2. `tests/test_silicoshark_smoke.py` — `LEGACY_FORTRAN` reproduces the
   FORTRAN seal goldens within v1 tolerance.
3. `scripts/run-discretisation-study.py` — preset × parameter-file
   matrix runner with per-cell metrics.
4. `docs/research/discretisation-audit.md` — citation trail per
   Discretisation field option.
5. (Eventually) a briefing document to support a human-written manuscript draft drawing on (3) and (4).
