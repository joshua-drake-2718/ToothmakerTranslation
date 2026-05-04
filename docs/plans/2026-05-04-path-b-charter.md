---
title: 'Path B charter — paper-faithful re-implementation'
author: Lyndon Drake (with Claude Code)
date: 2026-05-04
---

## Aim

Implement the tooth-morphogenesis model directly from the published
equations, in idiomatic NumPy, in a new module that does not modify
the existing FORTRAN-faithful `coreop2d.py`. The result should read
like the papers, run faster than the FORTRAN-faithful translation
(which is ~2× the FORTRAN binary on the seal example), and be
maintainable by a domain scientist who has never seen the FORTRAN.

The FORTRAN-faithful Python remains the cross-validation oracle.

## Paper scope

Two papers, in priority order:

1. Salazar-Ciudad and Jernvall, 'A computational model of teeth and
   the developmental origins of morphological variation', *Nature*
   464, 583–586 (2010). DOI [10.1038/nature08838](https://doi.org/10.1038/nature08838).
   The canonical model. Equation source of truth.
2. Harjunmaa, Seidel, Häkkinen et al., 'Replaying evolutionary
   transitions from the dental fossil record', *Nature* 512, 44–48
   (2014). DOI [10.1038/nature13613](https://doi.org/10.1038/nature13613).
   The tribosphenic extension that produced `humppa_translate.f90`,
   itself the ancestor of `13.f90`. Use for parameters specific to
   tribosphenic morphology and for any equations not present in the
   2010 paper.

**Out of scope:** Zimm, Berio, Debiais-Thibaud, Goudemand, *PNAS*
120(15) e2216959120 (2023). This is a *shark* extension that
branched off the same code base in 2019 (`silicoshark/toothmodel.f90`).
Reading the SI is fine for parameter naming conventions, but it is
not the model `13.f90` implements. See
`docs/research/path-b-research.md` §Paper for the full provenance
chain.

## Reference materials, in priority order

1. The two Nature papers and their SIs (`.tmp/papers/`):
   - `nature08838.pdf` + `41586_2010_BFnature08838_MOESM283_ESM.pdf`
     — Salazar-Ciudad and Jernvall 2010, the canonical model. The SI
     has an 18-equation listing and the parameter table (SI Table 1).
   - `nature13613.pdf` + `41586_2014_BFnature13613_MOESM78_ESM.pdf`
     — Harjunmaa et al. 2014, the tribosphenic extension. Note: the
     SI has *no* equations and *no* parameter table; the 2014 paper
     uses the 2010 model as a black box. The wild-type tribosphenic
     parameter set is in Ext Data Fig. 2 of the main article.
   - The full reviews are in
     `docs/research/paper-review-2010-salazar-ciudad-jernvall.md`
     and `docs/research/paper-review-2014-harjunmaa.md`. Read those
     before reading the PDFs themselves; they extract every equation,
     parameter, and known paper-vs-code gap.

2. `examples/seal.txt` — the parameter file the existing FORTRAN-
   faithful Python is validated against. Use this as the default-
   values backstop for any parameter the papers don't pin down.

3. The 2014 paper's Ext Data Fig. 2 — fully-specified wild-type
   tribosphenic mouse parameter set (25 named parameters). Best
   single authoritative parameter listing for tribosphenic teeth.
   Captured into a table in `paper-review-2014-harjunmaa.md`.

4. `cpp_README.md` from jernvall-lab's C++ port — variable and
   parameter naming glossary only. Do **not** read the C++ source as
   a structural template. Located at
   `.tmp/jernvall-toothmaker/cpp_README.md`.

5. `humppa_translate.f90` — the original Catalan FORTRAN. Use as
   line-by-line ground truth where the papers are ambiguous, and to
   distinguish 'what the model does' from 'what `13.f90` happens to
   do'. Located at `.tmp/jernvall-toothmaker/humppa_translate.f90`.

6. `tgrohens/coreop2d.f90` — minimal, clean FORTRAN fork. Easier to
   diff than humppa's monolith, no semantic differences. Located at
   `.tmp/tgrohens-toothmaker/coreop2d.f90`.

7. `13.f90` and the existing `coreop2d.py` — last resort, with full
   awareness of the renaming-introduced divergences from humppa
   catalogued in `docs/research/13f90-vs-humppa-divergences.md`.

**Use cautiously:** `silicoshark/multicusp.ini` was previously listed
as a parameter glossary, but the 2023 review
(`docs/research/paper-review-2023-zimm.md`) established that the file
is not a clean 2010 baseline — it appends shark-specific 'Further
additions' (`Egr(i)`, `A(i)`, Shh threshold) that are not in the
2010+2014 mammalian model. Reading it for symbol conventions is fine;
do not use it as the parameter source-of-truth. `examples/seal.txt`
plus the 2014 Ext Data Fig. 2 are better.

## Paper-vs-FORTRAN precedence

Two known cases where the paper and the FORTRAN disagree, with the
disagreement large enough to need an explicit Path B decision:

1. **2010 paper Eq. (14) denominator typo.** The paper prints the
   inhibitor saturation term as `1 + k_Inh [Act]`. The FORTRAN
   computes `1 + Inh · [Inh]`, which is the biologically meaningful
   form (an inhibitor saturating its own response). Path B follows
   the FORTRAN. See
   `docs/research/paper-review-2010-salazar-ciudad-jernvall.md`
   §'Paper-vs-code findings' for the evidence.

2. **Differentiation accumulator.** The 2010 paper describes
   differentiation as accumulating activator (Act) over time. The
   FORTRAN accumulates secretion / growth-factor (Sec) instead
   (humppa line 659; preserved in `13.f90` and the existing
   `coreop2d.py`). This is undocumented in either the 2010 or 2014
   paper. Path B should explicitly choose: follow the paper (Act
   accumulation) or follow the FORTRAN (Sec accumulation).
   Recommendation: follow the FORTRAN, with the choice documented
   in the docstring and the rationale recorded as a `docs/findings/`
   entry. See
   `docs/findings/2026-05-04-differentiation-uses-sec-not-act.md`.

For all other paper-vs-FORTRAN gaps, prefer the paper.

## Underspecified — implementer's choice

A consolidated list of points where the 2010 and 2014 papers leave
the implementation choice open. Path B is free to make the choice
that produces clean, vectorisable numpy. The relevant reviews
enumerate these in detail; key headline freedoms:

- **Diffusion Laplacian discretisation.** The papers describe
  reaction-diffusion in continuous form; the spatial discretisation
  (per-edge margin walk vs cotangent weights vs Voronoi-area
  weights) is not specified. Path B can use any choice that gives
  the right qualitative behaviour. Cotangent weights from
  `scipy.spatial.Delaunay` triangulation are an idiomatic numpy
  choice.

- **Update order (Jacobi vs Gauss–Seidel).** The papers' prose is
  consistent with pure Jacobi (compute all deltas first, then apply
  globally). The FORTRAN mixes: some terms read in-place updated
  state. Pure Jacobi is simpler and not load-bearing for the
  equilibrium.

- **Cell-division mesh repair.** The 2010 SI gives the mesh-repair
  step as a *picture* (SI Fig. 1), not as code or pseudocode. Path
  B can use `scipy.spatial.Delaunay` re-triangulation around the
  inserted node, eliminating the entire `13.f90` topology-walk
  chain (the txungu else-branch, iiii preservation,
  panic-and-return). This is the single biggest win.

- **Border-cell definition.** Two definitions exist in the 2010
  paper — topological (descended from the initial 12 outer cells)
  and geometric (cells with fewer than 6 neighbours). The paper
  conflates them. Path B can use either, with the choice
  documented.

- **Initial cell layout.** The 2010 paper specifies a 7-cell
  hexagon; `13.f90` uses 19 cells (Rad=4). The 2014 wild-type
  parameter set also uses Rad=4. The 19-cell layout is a 2014
  choice not specified in either paper. Path B should accept Rad
  as a parameter and not hardcode.

- **Knot-daughter `d_i` reset, knot threshold semantics, 'tenth
  border' boundary handling, `pushingnovei` non-neighbour
  repulsion details** — all enumerated in the 2010 review's
  Path-B-implications list. Path B has freedom on each.

For the full enumeration, see
`docs/research/paper-review-2010-salazar-ciudad-jernvall.md`
(15-item P/I list) and
`docs/research/paper-review-2014-harjunmaa.md`
(8-item decision-points list).

## What NOT to reproduce

`13.f90` and our FORTRAN-faithful Python preserve a chain of
accidents (some introduced by the Catalan→English renaming, some
present in the canonical FORTRAN). Path B should reproduce none of
them:

- **The txungu else-branch's single-neighbour cells.** This is a
  `13.f90`-only renaming inconsistency (see chain step 1 of the
  `coreop2d.py` docstring). humppa does not have it. Path B should
  not generate single-neighbour cells in the first place.
- **The FMA-divide-by-zero dodge in `apply_diffusion`.** With no
  single-neighbour cells, the divide is never zero. No FMA-limit
  guard is needed; the whole `sum_a == 0` branch can disappear.
- **The label-77 fall-through in `calculate_margins`.** Empty
  neighbour slots should be empty, not filled with the previous
  slot's value.
- **The silent-panic-and-return in the topology walk.** The
  topology walk should either converge on a valid placement or
  fail loudly. Better: if vectorised numpy lets us avoid topology
  walks entirely, do that.
- **The C++ port's `if (suma > 0)` skip in apply_diffusion.** It
  has no validation evidence and produces different deltas for
  degenerate cells than either the FORTRAN binary or our limit
  substitution. See `docs/research/cpp-port-review.md`.
- **The C++ port's file split.** Designed for C++ compilation
  units; in NumPy several of the files become trivially small or
  artificially separated.

One latent `13.f90`-only quirk to be aware of (per
`docs/research/13f90-vs-humppa-divergences.md` entry 7): `13.f90`'s
`setparams` silently drops `Swi` (parap slot 6) — the slot is marked
`!Nothing`. The seal example has `Swi = 0` so this hasn't bitten the
existing translation, but Path B's parameter loader should accept
`Swi` as a real parameter and route it into the border-band logic in
`apply_border_bias` (where humppa uses it). If you encounter a
parameter file expecting `Swi != 0` to take effect, this is the
reason the FORTRAN-faithful Python ignores it.

## Architecture goals

1. **Vectorised arrays, not per-cell loops.** The diffusion,
   repulsion, and growth-pushing operations all admit vectorised
   formulations against `(num_cells, ...)` arrays. The FORTRAN's
   per-cell-per-neighbour scalar loops are a translation artefact
   of FORTRAN-77 idioms; the equations themselves are matrix
   operations.
2. **Named variables matching the papers.** Use `cell_positions`,
   `neighbours`, `cell_margins`, `concentrations` (or domain-
   appropriate analogues like `act`, `inh`, `fgf`, `ect`),
   `growth_rate`, etc. Do not use `malla`, `vei`, `marge`, `q3d`,
   `tacre`, `acac`. Map the FORTRAN names to the paper names in
   one place (a glossary module or a comment block) for those who
   need to cross-reference.
3. **Equation numbers in docstrings.** Each function that
   implements a paper equation should cite it: '(Eq. 1, 2010 paper)'
   or 'p. S2 of the 2014 SI'.
4. **Module name:** `silicoshark/`. (We checked; it's not in use
   elsewhere in this repo.)
5. **No FORTRAN-isms.** No `goto`-equivalents, no implicit
   variable lifetimes carrying state across iterations, no
   one-letter variable names except `i` for indices.

## Validation strategy

The FORTRAN-faithful Python's golden outputs (`tests/golden_fortran/
*.off`, captured from the FORTRAN binary on `examples/seal.txt`)
remain the cross-validation oracle. Tolerance:

1. **Cell count plateau within ±5%** of the FORTRAN's value (57 on
   the seal example). 'Plateau' means the value at the last save
   block (block 5).
2. **Vertex envelope match.** Min/max x, y, z over all vertices
   should be within 10% of the FORTRAN's, accounting for sign of
   the asymmetry (left side, right side, top, bottom).
3. **No NaN or Inf** anywhere in the output.
4. **Qualitative cell-count trajectory.** The growth curve should
   be monotonically increasing then plateau; not exponentially
   diverging, not collapsing.

Byte equality is **not** a goal. The FORTRAN's specific numerical
choices (rounding to 9 decimals in distance checks, FMA noise on
degenerate cross products, particular normalisation orders) are
implementation accidents, not paper equations.

A second optional golden: tgrohens' `nosea_10261__0.txt` parameter
file (`.tmp/tgrohens-toothmaker/nosea_10261__0.txt`) exercises a
non-seal regime with `Rad=2` (7-cell initial condition), Egr 0.0135,
Mgr 204, Inh 26. Useful for testing generalisation. Capture a
golden run from the FORTRAN binary (or from our FORTRAN-faithful
Python) on this input and add to `tests/golden_fortran/` as
`nosea_*.off`. Defer this until Path B passes the seal goldens.

### Phase 2 validation goldens (from the 2014 paper)

Once the seal golden passes, the 2014 paper review identifies several
qualitative validation targets that exercise the model across
parameter space and provide stronger evidence than a single example:

- **Wild-type tribosphenic mouse** (Ext Data Fig. 2 of the 2014
  paper) — 25-parameter set producing a ~5-cusp tribosphenic
  morphology. Best single authoritative target.
- **Act sweep** (Ext Data Fig. 4) — Act values from 0.1 to 1.6 in
  0.1 increments produce a monotonic 1→5 cusp progression. The
  cleanest single-axis validation Path B can run.
- **Inh / Di perturbation lattice** (Ext Data Fig. 5) — five
  concrete (Inh, Di) points showing the Eda-null baseline plus
  four multi-cusp recoveries.
- **Parameter sweep map** (Ext Data Fig. 3) — qualitative
  directions (per of nine parameters) for cusp-count reduction.
  Useful as a 'sanity-check' validator: each direction should
  reduce cusp count if Path B is correct.

Capture goldens for these from the existing FORTRAN-faithful Python
(no need to rebuild the FORTRAN binary), with relaxed tolerances
where the original paper used qualitative comparisons.

The new test module is `tests/test_silicoshark_smoke.py` (not
`test_simulator_smoke.py`, which is reserved for the FORTRAN-
faithful translation). Both pytest modules should pass.

## Investigation method (carry-over from Path A)

Per-stage stderr tracing of both implementations, diffed at
progressively finer granularity until each divergence has a
single root cause. See the 'Investigation method' section of
`docs/findings/2026-05-04-path-a-fma-and-over-division.md`.

When Path B disagrees with the FORTRAN goldens beyond tolerance,
investigate deliberately: is the difference (a) a paper-vs-FORTRAN
discrepancy worth recording in `docs/findings/`, (b) a Path B bug,
or (c) a benign floating-point difference?

## Out of scope for Path B

- Modifying `coreop2d.py`, `esclec.py`, `main.py`, `vector.py`, or
  `13.f90`.
- Modifying `tests/test_simulator_smoke.py`.
- Modifying the existing FORTRAN goldens.
- Re-implementing the GUI, the parameter-file format, or any
  pre-2010 versions of the model.
- Implementing the Zimm 2023 shark variant.
