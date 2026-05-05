---
title: 'Path B v2 — Discretisation field audit trail'
author: Lyndon Drake (with Claude Code)
date: 2026-05-05
---

## Purpose

This document is the per-field, per-option audit of every choice
point exposed by the `silicoshark.discretisation.Discretisation`
dataclass. Where the v2 charter
(`docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md`)
sets policy and sketches each field at the level needed to design
the dataclass, this document supplies the granular evidence: for
each option of each field, the paper passage that frames the
question, the FORTRAN line(s) that resolve it one way, any
alternative code-base in the lineage that resolves it differently,
and the comparison-study question the option exists to answer. It
is the audit trail the charter promised but did not contain.

The methodological role is to make the implementer-choice surface
of an under-specified scientific code visible and traceable. The
2010 paper has the equations but routinely under-specifies their
discrete realisation; the FORTRAN ancestors (humppa, `13.f90`,
tgrohens, silicoshark, the C++ port) make concrete choices that the
paper does not constrain and rarely even names; and a reader who
asks 'where did this number come from?' has, until now, had to
read all four code-bases to answer. The follow-up paper draws on
this audit to argue that LLM-assisted analysis can identify and
ground each implementer-choice point in published text and in
alternative code-bases at a level of resolution that human-only
review struggles to sustain, and that the resulting catalogue is
itself a methodological output worth publishing alongside the
science.

## How to read this document

Each field is a section. Each section opens with the question being
decided, the default (the field's value in `PATH_B_DEFAULT`), and
the 'kind' of disagreement (`PaperAmbiguity`, `PaperVsCodeTension`,
or `FortranAccident`, following the charter taxonomy). Each option
is a sub-section with bullets covering paper evidence, FORTRAN
evidence (line numbers in `13.f90` where I have them, otherwise
`coreop2d.py` line numbers — Path A maps onto `13.f90`
line-by-line), alternative-code-base evidence where any such
evidence exists, and the comparison-study question the option
isolates. Where a finding has been written up separately, I cite
the `docs/findings/` file rather than reproduce its content.

Citations follow the convention 'main p. N' for the 2010 *Nature*
letter and Methods, 'SI p. N' for the supplementary PDF, and 'EDF'
for the 2014 Extended Data figures. Equation numbers follow the
2010 paper.

## Fields

### `laplacian` — kind: PaperAmbiguity

**Question:** how is `nabla^2 c` discretised on the triangular mesh?

**Default:** `length_weighted` — the simplest defensible Fick's-law
form, recommended by the v1 charter as the cleanest numpy-idiomatic
choice. The cotangent form is the textbook discrete
Laplace–Beltrami operator and is on the shortlist for promotion
once the comparison study reports.

#### Option `length_weighted`

- **Paper evidence:** main p. 585 prescribes 'finite volume on the
  triangular mesh, with flux proportional to contact area' but
  gives no formula. The length-weighted form
  `L u[i] = sum_j (u[j] - u[i]) / |p_j - p_i|` is one defensible
  Fick's-law reading: it is dimensionally consistent with the prose
  description if 'contact area' is read as 'inverse distance'
  on a uniform-edge mesh, and it matches the most common graph-
  Laplacian convention in numerical practice.
- **FORTRAN evidence:** not used. Neither `humppa_translate.f90`
  nor `13.f90` computes a length-weighted Laplacian; both go via
  the per-edge `pes` margins (see `fortran_margins`).
- **Alternative code-bases:** none in the lineage compute a
  length-weighted graph Laplacian. This is a numpy-idiomatic
  choice introduced by Path B v1.
- **Comparison-study question:** does this simpler form change
  cusp count, cusp width, or qualitative regime relative to the
  cotangent operator or the `fortran_margins` form?

#### Option `cotangent`

- **Paper evidence:** no direct support. The paper's 'flux ∝
  contact area' prose is more sympathetic to a contact-area-based
  weight, of which the cotangent operator is the standard discrete
  realisation on triangulated surfaces.
- **FORTRAN evidence:** not in any FORTRAN. Neither
  `humppa_translate.f90`, `13.f90`, nor tgrohens' `coreop2d.f90`
  implements a cotangent operator. The C++ port likewise does not.
- **Alternative code-bases:** the cotangent discrete
  Laplace–Beltrami operator is standard textbook material for
  curved-surface diffusion (see e.g. Meyer, Desbrun, Schröder,
  Barr, 'Discrete differential-geometry operators for triangulated
  2-manifolds', *Visualization and Mathematics III*, 2003). It is
  the form most defensible from the paper's prose if the
  triangulated mesh is read as a discrete curved surface rather
  than a flat graph.
- **Comparison-study question:** does the cotangent operator's
  curvature awareness change behaviour on parameter sets where
  z-displacement is large (i.e. the seal example, where
  cervical-loop downgrowth produces ~125-unit z-spread)?

#### Option `fortran_margins` (deferred — `NotImplementedError`)

- **Paper evidence:** no direct support. The paper does not
  describe per-edge margins or area_p z-weighting at any level of
  detail.
- **FORTRAN evidence:** `13.f90:392-466`
  (`coreop2d.py:506-666`). The FORTRAN computes a per-edge `pes`
  margin from the neighbour `border` array and a per-cell
  `area_p` z-weighting that drives vertical diffusion between
  z-layers. The hand-tuned `0.044D1` factor at the substrate
  boundary (`13.f90:431, 444, 459`) is a load-bearing magic
  number — see findings on FMA and over-division for context. The
  approach is described in the C++ port glossary (`cpp_README.md`)
  as 'internode positions' but the inner indexing
  (`cellMargins[i][j][k]` with k=0..7 carrying margin x/y/z plus
  distance components) is not paper-specified.
- **Alternative code-bases:** the C++ port at
  `tooth_model_diffusion.cpp` reproduces the FORTRAN margins
  (with the `if (suma > 0)` guard documented in
  `docs/research/cpp-port-review.md` §2(a) producing a third
  behaviour distinct from both FORTRAN and Path A's Python).
  tgrohens' `coreop2d.f90:401-571` is byte-identical to humppa.
- **Comparison-study question:** how much of the FORTRAN's
  seal-example trajectory is driven by the per-edge margin
  geometry versus the simpler graph-Laplacian alternatives? This
  option is deferred behind a `NotImplementedError` because
  reproducing the FORTRAN margins inside the v2 vectorised
  pipeline is its own sub-project; the residual z-min
  divergence in
  `docs/findings/2026-05-05-path-b-v2-legacy-fortran-divergence.md`
  may or may not close once `fortran_margins` is implemented.

---

### `update_order` — kind: PaperVsCodeTension

**Question:** are forces and reactions evaluated against the
start-of-step state ('Jacobi') or against in-place updated state
('Gauss–Seidel')?

**Default:** `jacobi` — the paper's prose strongly suggests Jacobi
(each `partial / partial t` is a function of the current state),
and pure Jacobi is naturally vectorisable.

#### Option `jacobi`

- **Paper evidence:** main p. 583 (fig. 1 caption) writes 'all
  equations are integrated using the Euler method', with no
  mention of in-place updates. The paper's eqs. (1)–(13) are each
  written as `partial p_i / partial t = ...` against a single
  state, which reads as Jacobi to anyone trained in numerical
  PDE.
- **FORTRAN evidence:** not used for forces. The reaction–
  diffusion pass at `13.f90:392-518`
  (`coreop2d.py:506-665`) is Jacobi (uses copy `hq3d`); the
  mechanical-force passes are not.
- **Alternative code-bases:** none.
- **Comparison-study question:** how much of the FORTRAN's
  trajectory is driven by the Gauss–Seidel ordering of mechanical
  forces, as opposed to the choice of force formulae?

#### Option `gauss_seidel_forces`

- **Paper evidence:** none. The paper does not describe a
  Gauss–Seidel update.
- **FORTRAN evidence:** `13.f90:1632-1645`
  (`coreop2d.iteration`). The iteration order is
  `apply_diffusion → apply_border_bias → apply_differentiation →
  ep_growth_border_force → boy_force → repel_non_neigh →
  repulse_neighbour → apply_nuclear_traction →
  update_cell_position`. Each force pass reads from `positions`
  and writes to `forces`; `update_cell_position` then commits
  forces to positions. But `apply_nuclear_traction`
  (`coreop2d.py:955-1005`) writes positions directly, not forces,
  so any later force pass sees post-traction positions.
- **Alternative code-bases:** the C++ port follows the FORTRAN
  ordering (`cpp_README.md` and the partially-mirrored
  `tooth_model_iteration.cpp`); tgrohens is byte-identical to
  humppa.
- **Comparison-study question:** under FORTRAN-faithful update
  order, how does the seal-example trajectory differ from pure
  Jacobi?

#### Option `mixed_jacobi_gs`

- **Paper evidence:** none.
- **FORTRAN evidence:** alias for `gauss_seidel_forces` in the
  current implementation (see `silicoshark/simulator.py` A5
  notes). Reserved for future variants where Gauss–Seidel is
  applied at finer per-component granularity (e.g. boy / repel /
  repulse each writing positions in turn rather than all
  reading the same post-traction state).
- **Alternative code-bases:** none.
- **Comparison-study question:** placeholder for variants the
  matrix may need; not currently distinguishable from
  `gauss_seidel_forces`.

---

### `topology` — kind: PaperAmbiguity

**Question:** how is the cell-cell adjacency graph maintained
across cell division?

**Default:** `static_with_local_update` — v1 evidence
(`docs/findings/2026-05-04-path-b-v1-mesh-degeneration.md`) shows
that re-triangulating each step is unstable on the seal example.

#### Option `delaunay_each_step`

- **Paper evidence:** main p. 583 (SI fig. 1) shows the
  triangulated mesh in a static picture; no algorithm is
  described. SI p. ~1 says 'a new cell adds a new cell centre and
  connections to the original triangular mesh', without giving
  the connection-update rule.
- **FORTRAN evidence:** not used. The FORTRAN's `add_cell`
  (`coreop2d.py:1034-1554`, `13.f90:1034-1554`) walks the
  topology incrementally rather than re-triangulating.
- **Alternative code-bases:** the v1 silicoshark prototype used
  this approach via `scipy.spatial.Delaunay` per iteration. The
  C++ port keeps the FORTRAN incremental walk.
- **Comparison-study question:** does global re-triangulation
  produce qualitatively different cusp morphology from
  static-with-local-update on parameter sets where it does not
  outright fail (e.g. early iterations of any parameter set)?
  v1 evidence shows it fails on the seal example by iter ~270
  via mesh collapse and `QHullError`.

#### Option `static_with_local_update`

- **Paper evidence:** indirect. The paper describes the initial
  hex lattice's neighbour count
  ('each epithelial cell has six epithelial neighbours, except
  for cells on the borders that have three or four', main
  p. 586) and the SI-fig.-1 connection update on division. The
  static-graph reading is consistent with both, so long as the
  initial graph is built from the hex topology and locally
  updated on division.
- **FORTRAN evidence:** `13.f90:1034-1554`
  (`coreop2d.add_cell`) walks the topology to insert the
  daughter and rewires neighbour lists in-place. The walk has
  three `13.f90`-only quirks documented in
  `docs/research/13f90-vs-humppa-divergences.md` §1, §6 (cosmetic);
  the walk's panic-and-return on degenerate cells is
  load-bearing for the seal-example plateau (see
  `docs/findings/2026-05-04-path-a-fma-and-over-division.md`).
- **Alternative code-bases:** humppa and tgrohens are
  byte-identical to each other; the C++ port follows the FORTRAN
  walk with one off-by-one in the panic threshold (cpp-port-review
  §2(c)).
- **Comparison-study question:** what qualitative behaviour
  changes when the topology-walk panic is replaced by
  `division_total_cap`? See `division_total_cap` for the
  empirical caveat.

---

### `division_per_step_cap` — kind: PaperAmbiguity

**Question:** maximum new cells that may be emitted in a single
iteration.

**Default:** `None` (unlimited). Inert under
`static_with_local_update`.

#### Option `None`

- **Paper evidence:** main p. 585 paragraph following eq. (6):
  cells divide when the connection between two centres reaches
  two units, with no cap stated.
- **FORTRAN evidence:** `13.f90:1034-1554` does not cap per-step
  divisions; it caps via the topology-walk panic
  (`coreop2d.py:1050` documents this behaviour).
- **Alternative code-bases:** none cap.
- **Comparison-study question:** baseline.

#### Option `int` (any positive integer)

- **Paper evidence:** none. Implementer-supplied safeguard.
- **FORTRAN evidence:** none. The FORTRAN gets bounded growth
  for free from the panic.
- **Alternative code-bases:** v1 silicoshark used a per-step
  cap (4 per pass per the v1 finding) to dampen
  Delaunay-each-step cascade behaviour.
- **Comparison-study question:** under
  `topology='delaunay_each_step'`, does a small per-step cap
  rescue stability without distorting the equilibrium? Inert
  under `static_with_local_update`.

---

### `division_max_edge` — kind: PaperAmbiguity

**Question:** reject cell-division on edges longer than this
multiple of the rest length.

**Default:** `None` (no cap). Inert under
`static_with_local_update`.

#### Option `None`

- **Paper evidence:** main p. 585: divide when the connection
  reaches two units. No upper edge-length cut-off.
- **FORTRAN evidence:** none. The FORTRAN's edges are
  intrinsically local because the topology is local, so this
  question does not arise.
- **Alternative code-bases:** none.
- **Comparison-study question:** baseline under
  `static_with_local_update`.

#### Option `float` (any positive float)

- **Paper evidence:** none. Implementer-supplied safeguard.
- **FORTRAN evidence:** none.
- **Alternative code-bases:** v1 silicoshark used a `max_edge_dist`
  filter to suppress long Delaunay convex-hull edges from triggering
  spurious divisions in the v1 mesh-degeneration scenario.
- **Comparison-study question:** does an edge-length cap rescue
  Delaunay-each-step stability on parameter sets v1 found
  unstable? Inert under `static_with_local_update`.

---

### `division_total_cap` — kind: FortranAccident (phenomenological match)

**Question:** maximum total cell count after division.

**Default:** `None` (unlimited). Set to 60 in `LEGACY_FORTRAN` and
`HUMPPA_LITERAL` to match the FORTRAN's seal-example plateau.

#### Option `None`

- **Paper evidence:** none. The paper describes no cap.
- **FORTRAN evidence:** none directly; but the FORTRAN's
  `add_cell` topology-walk panic at
  `13.f90:1248, 1257, 1291, 1299` (`coreop2d.add_cell`) aborts
  division when the walk fails to find a path, which on the seal
  example fires at iteration ~14 and pins the cell count at 57
  for the remainder of the run.
- **Alternative code-bases:** the C++ port panics one step
  earlier (`cj >= MAX_NEIGHBORS` rather than `cj > nvmax`) per
  cpp-port-review §2(c). humppa and tgrohens are byte-identical
  to each other in the panic logic.
- **Comparison-study question:** baseline; default behaviour for
  `PATH_B_DEFAULT` and `PAPER_2010` where the topology-walk
  panic is not modelled.

#### Option `60` (in `LEGACY_FORTRAN`)

- **Paper evidence:** none.
- **FORTRAN evidence:** the value 60 is ~5% above the seal
  example's observed 57-cell plateau, chosen as a tolerance band
  rather than an exact match. The FORTRAN's plateau is a
  phenomenological consequence of the panic, not a designed
  value.
- **Alternative code-bases:** none model a total cap explicitly.
- **Comparison-study question:** the cap exists because
  silicoshark's `static_with_local_update` topology has no analog
  of FORTRAN's `add_cell` topology-walk panic; without it the
  seal example divides exponentially. The cap is documented as a
  phenomenological match in the field's docstring; reproducing
  the panic faithfully is documented in
  `docs/findings/2026-05-05-path-b-v2-legacy-fortran-divergence.md`
  as future work.

---

### `eq5_z_gate` — kind: PaperAmbiguity

**Question:** does eq. 5's `sum(unit_ij)` include neighbours
strictly above in z, or all neighbours?

**Default:** `True` — v1 evidence shows the gate stabilises
flat-z initial conditions where `sum_unit ≈ 0` would otherwise be
amplified by unit-vector normalisation.

#### Option `True`

- **Paper evidence:** main p. 585, eq. (5):
  `partial p_i / partial t = k_egr (1 - d_i) sum_j unit_ij /
  |sum_j unit_ij|` along 'the average outward unit vector'. The
  paper does not define 'outward' for an initially-flat
  epithelial sheet; on a curved sheet the natural reading is
  'apically pointing'. The z-gate (`b < -1e-4` in `coreop2d.py`)
  is a discrete realisation of 'apical' that works for the
  initial flat-z hex.
- **FORTRAN evidence:** `coreop2d.py:702`
  (`13.f90:541` — `if (b<-1e-4)` inside
  `EpGrowthBorderForce`). The FORTRAN gates on `b < -1e-4` where
  `b = z_i - z_k`, so cells include only those neighbours
  strictly higher than themselves.
- **Alternative code-bases:** humppa and tgrohens are byte-
  identical here; the C++ port's mechanics file is missing from
  the mirror so I cannot confirm but the README does not flag a
  divergence.
- **Comparison-study question:** baseline. The v1 finding
  (`docs/findings/2026-05-04-path-b-v1-mesh-degeneration.md` §1)
  shows the gate is load-bearing for stability under flat-z
  initial conditions.

#### Option `False`

- **Paper evidence:** the paper's eq. (5) as written has no
  z-gate — `sum_j unit_ij` is over all mesh neighbours. The
  ungated form is the literal reading.
- **FORTRAN evidence:** not used.
- **Alternative code-bases:** none.
- **Comparison-study question:** does ungating break the
  initial-step force-free property? v1 evidence says yes for the
  seal example; the comparison study should test other parameter
  sets to see if the gate is necessary in regimes where
  cervical-loop curvature develops earlier.

---

### `eq5_apply_to` — kind: PaperVsCodeTension

**Question:** does eq. 5 apply to all cells (paper) or interior
cells only (FORTRAN)?

**Default:** `'all'` — the paper's eq. (5) has no
interior/border distinction.

#### Option `'all'`

- **Paper evidence:** main p. 585, eq. (5) is written without
  any restriction on cell index. The paper's prose elsewhere
  ('Growth biases', main p. 585) splits anterior and posterior
  border cells but does not exclude interior cells from eq. (5).
- **FORTRAN evidence:** not used.
- **Alternative code-bases:** none.
- **Comparison-study question:** paper-faithful baseline.

#### Option `'interior_only'`

- **Paper evidence:** none directly. The paper splits the
  cervical-loop downgrowth (eqs. 7–9) onto border cells and the
  growth (eq. 5) onto all cells; the FORTRAN's interior-only
  split is one possible reading but not the literal one.
- **FORTRAN evidence:** `coreop2d.py:692`
  (`13.f90:536` — `do i=first_border_cell,num_active_cells`).
  The FORTRAN's `EpGrowthBorderForce` runs eq. 5 only over cells
  `i >= first_border_cell` (the inner ~2/3 of the lattice).
  Border cells (`i < first_border_cell`) receive the
  cervical-loop downgrowth force only, with eq. (5) excluded by
  loop bounds.
- **Alternative code-bases:** humppa and tgrohens preserve this;
  the C++ port's mechanics file is missing from the mirror.
- **Comparison-study question:** how much of the FORTRAN's
  cervical-loop dynamics is driven by removing eq. (5) on
  border cells? On the seal example this is one of the gates A6
  identified as load-bearing for the z=125 envelope.

---

### `rep_form` — kind: PaperAmbiguity

**Question:** is eq. 1 always-on Hookean (smoothly signed across
the rest length) or paper-gated (active only when `|d| < p_0`)?

**Default:** `hookean_signed` — v1 evidence shows the gated form
amplifies floating-point noise at exact-rest-length initial edges.

#### Option `hookean_signed`

- **Paper evidence:** main p. 585, eq. (1):
  `f_ij = k_rep (|d| - |p_0|)(p_j - p_i)`. The formula as
  written is a Hookean spring centred on `|d| = p_0`. Reading
  the paper's eq. 1 in isolation, without the eq. 2 framing as
  'attraction' that follows, the natural reading is a smoothly
  signed force that transitions from repulsion to weak
  attraction at the rest length.
- **FORTRAN evidence:** not exactly. The FORTRAN's
  `repulse_neighbour` (`coreop2d.py:841-882`,
  `13.f90:671-719`) tests `dr < rd` and applies the Hookean
  formula only on that branch.
- **Alternative code-bases:** none in the lineage use a
  smoothly-signed Hookean form. This is a numpy-idiomatic choice
  introduced by Path B v1 to avoid floating-point edge cases at
  exact equilibrium.
- **Comparison-study question:** does the smooth form cause
  weak-attraction artefacts at the rest length under non-seal
  parameter sets?

#### Option `paper_gated`

- **Paper evidence:** main p. 585. The paper's prose treats
  eq. 1 (repulsion) and eq. 2 (adhesion) as separate forces, with
  eq. 1 'a constant attraction force resulting from the adhesion
  of cells' kicking in 'when the distance between two cell
  centres is larger than `|p_0|`'. The gating reading takes this
  prose seriously: eq. 1 is repulsion, active only when
  `|d| < p_0`.
- **FORTRAN evidence:** `coreop2d.py:863-872`
  (`13.f90:691-700`). The FORTRAN gates on `dr < rd` (where `rd`
  is the per-edge rest length, also held at `|p_0| = 1.0` for
  initial edges) before applying the Hookean formula. Else
  branch handles attraction via Adh (see `adh_form`).
- **Alternative code-bases:** humppa and tgrohens preserve the
  gating; the C++ port follows.
- **Comparison-study question:** under FORTRAN-faithful gating,
  how much does the floating-point noise at unit edges
  contribute to the trajectory? v1 evidence says it amplifies
  noticeably; the comparison study quantifies under the
  static-topology code path.

---

### `adh_form` — kind: PaperVsCodeTension

**Question:** is eq. 2 a unit-vector pull (paper) or a Hookean
attraction (FORTRAN)?

**Default:** `unit_vector` per the paper.

#### Option `unit_vector`

- **Paper evidence:** main p. 585, eq. (2):
  `f_ij = k_adh (p_j - p_i) / |p_j - p_i|`. This is a unit-vector
  pull with magnitude `k_adh`, active only when
  `|p_j - p_i| > |p_0|`. The paper's Methods note explicitly
  flags `|p_0|` as 'the initial distance between cells', not a
  parameter.
- **FORTRAN evidence:** not used.
- **Alternative code-bases:** none.
- **Comparison-study question:** paper-faithful baseline.

#### Option `hookean_attraction`

- **Paper evidence:** none directly. The paper's eq. 2 as
  written is unit-vector; the Hookean form is what the FORTRAN
  computes regardless.
- **FORTRAN evidence:** `coreop2d.py:870`
  (`13.f90:698`). In `repulse_neighbour`'s else-branch (entered
  when `dr >= rd`, i.e. for cells beyond rest length), the
  FORTRAN sets `persu[j, :] = ux/uy/uz * Adh`, multiplying the
  un-normalised displacement vector by `Adh`. The result is a
  Hookean attraction (force grows linearly with separation),
  not a unit-vector pull (force is constant in magnitude).
- **Alternative code-bases:** humppa and tgrohens preserve the
  Hookean attraction.
- **Comparison-study question:** the FORTRAN form is stronger
  at distance, which may matter for late-stage cell separation;
  the paper form is constant-magnitude and stops growing once
  beyond rest length. Which produces more biologically
  plausible cusp shapes?

---

### `rep_neighbour_set` — kind: PaperVsCodeTension

**Question:** does repulsion act only on mesh neighbours, or also
on non-mesh cells whose centres are close?

**Default:** `mesh` per the paper.

#### Option `mesh`

- **Paper evidence:** main p. 585, eq. (3):
  `partial p_i / partial t = sum_{j in N(i)} f_ij` over the
  paper's neighbour set `N(i)`. The neighbour set is not
  defined explicitly, but the rest of the paper consistently
  treats `N(i)` as the triangulated-mesh neighbour set.
- **FORTRAN evidence:** `coreop2d.py:841-882`
  (`13.f90:671-719`, `repulse_neighbour`). This is the
  paper-faithful loop over mesh neighbours.
- **Alternative code-bases:** none deviate.
- **Comparison-study question:** paper-faithful baseline.

#### Option `mesh_plus_all_close`

- **Paper evidence:** none. This is implementer-added.
- **FORTRAN evidence:** `coreop2d.py:886-934`
  (`13.f90:723-790`, `repel_non_neigh` /
  Catalan `pushingnovei`). After the mesh-neighbour repulsion
  pass, the FORTRAN runs an additional loop over all cells whose
  centre is within ~1.4 units, regardless of mesh adjacency, and
  applies a polynomial repulsion `dd = 1/(d+1)**8; f = -d_vec *
  (dd/d) * Rep`. The 1.4-unit cut-off and the polynomial form
  are not in the paper.
- **Alternative code-bases:** humppa and tgrohens preserve the
  loop; the C++ port has it (declared in
  `tooth_model_iteration.cpp:24` even though
  `tooth_model_mechanics.cpp` is missing from the mirror).
- **Comparison-study question:** is this loop load-bearing for
  mesh-tangling prevention? The 2010 paper review
  (§'Discrepancies between paper and `13.f90`', point 3) notes
  the loop is 'plausibly load-bearing for stability'. The
  comparison study should test whether removing it produces
  mesh tangling on the seal example or any other parameter set.

---

### `eq14_denominator` — kind: PaperVsCodeTension

**Question:** does the eq. 14 inhibitor saturation term divide by
`(1 + k_inh * [Act])` (paper text) or `(1 + k_inh * [Inh])`
(FORTRAN — the biologically meaningful form)?

**Default:** `inh_corrected` — the paper's prose elsewhere
consistently treats Inh as inhibiting Act, so the FORTRAN form is
biologically right and the paper text is a typo.

#### Option `inh_corrected`

- **Paper evidence:** the paper's prose at main p. 585–586
  describes Inh as 'an enamel knot-secreted inhibitor of enamel
  knot formation'; the eq. 14 denominator should saturate the
  Act-inhibition by Inh. The 2010 paper review document
  (`docs/research/paper-review-2010-salazar-ciudad-jernvall.md`
  §'Activator (eq. 14)') concludes the paper's `[Act]` in the
  denominator is a typo and the FORTRAN's `[Inh]` is the
  authoritative form.
- **FORTRAN evidence:** `13.f90:485`
  (`coreop2d.py:645`):
  `hq3d(i,1,1) = a/(1 + Inh*q3d(i,1,2)) - Deg*q3d(i,1,1)`. The
  third index `2` selects Inh.
- **Alternative code-bases:** humppa, tgrohens and the C++ port
  all use the Inh-corrected form. There is no published
  defender of the literal paper form.
- **Comparison-study question:** baseline.

#### Option `act_typo`

- **Paper evidence:** main p. 586, eq. (14):
  `partial[Act] / partial t = k_act [Act] / (1 + k_inh [Act]) -
  k_deg [Act] + k_da nabla^2 [Act]`. Read literally, the
  denominator depends on Act, making eq. (14) a saturating
  self-inhibition of the activator with no Inh involvement.
- **FORTRAN evidence:** none.
- **Alternative code-bases:** none.
- **Comparison-study question:** does the typo materially change
  biology on the seal example or on the 2014 EDF Act-sweep
  validation set? The `PAPER_LITERAL_2010` preset isolates this
  one field at the typo'd value to answer this question.
  Indistinguishable behaviour means the typo is harmless;
  divergent behaviour means the typo would have produced
  different published results had the paper text been used as
  the implementation.

---

### `eq17_inh_source` — kind: PaperVsCodeTension

**Question:** for cells with `d_i >= Int`, is the eq. 17 Inh
source `[Act]` (paper) or `(rate of [Act]) * d_i` (FORTRAN
temp-variable form)?

**Default:** `act_concentration` per the paper.

#### Option `act_concentration`

- **Paper evidence:** main p. 586, eq. (17):
  `partial[Inh] / partial t = [Act] - k_deg [Inh] +
  k_di nabla^2 [Inh]` for cells with `d_i >= k_int`. The source
  is the Act *concentration*.
- **FORTRAN evidence:** not used.
- **Alternative code-bases:** none.
- **Comparison-study question:** paper-faithful baseline.

#### Option `act_rate_times_di`

- **Paper evidence:** none.
- **FORTRAN evidence:** the form differs between
  `13.f90` and Path A's `coreop2d.py`. `13.f90:487` reads
  `hq3d(i,1,2) = q3d(i,1,1)*DiffState(i) - Deg*q3d(i,1,2)`,
  which is `[Act] * d_i - Deg * [Inh]` (Act concentration
  times d_i). `coreop2d.py:647` reads `hq3d[i, 0, 1] =
  hq3d[i, 0, 0] * cls.diff_state[i] - cls.Deg * cls.q3d[i, 0, 1]`,
  where `hq3d[i, 0, 0]` is the eq. 14 *rate* (post-numerator,
  post-denominator, post-degradation), not the Act
  concentration. The Path A `coreop2d.py` form uses the
  temporary variable that 13.f90 *would* have used had
  someone thought the d_i factor through, so silicoshark's
  `act_rate_times_di` follows the Path A reading. Either way
  the paper's plain `[Act]` source is replaced by a d_i-modulated
  form.
- **Alternative code-bases:** humppa and tgrohens preserve
  13.f90's `q3d(i,1,1)*DiffState(i)` form. The C++ port at
  `tooth_model_diffusion.cpp:218` matches the FORTRAN.
- **Comparison-study question:** the d_i ramp smooths the Inh
  source onset; the paper's plain `[Act]` is binary on/off at
  threshold. Does the smoother onset materially change the
  knot dynamics? Note also the loose match between `13.f90`'s
  literal form and Path A's temp-variable form is itself a
  paper-vs-code finding worth reporting.

---

### `eq18_sec_source` — kind: PaperVsCodeTension

**Question:** for cells with `d_i >= Set`, is the eq. 18 Sec
source the constant `k_sec` (paper) or `k_sec * d_i` (FORTRAN
ramp)?

**Default:** `constant_k_sec` per the paper.

#### Option `constant_k_sec`

- **Paper evidence:** main p. 586, eq. (18):
  `partial[Sec] / partial t = k_sec - k_deg [Sec] +
  k_ds nabla^2 [Sec]` for cells with `d_i >= k_set`. The
  source is the constant `k_sec`.
- **FORTRAN evidence:** not used.
- **Alternative code-bases:** none.
- **Comparison-study question:** paper-faithful baseline.

#### Option `k_sec_times_di`

- **Paper evidence:** none.
- **FORTRAN evidence:** `13.f90:494`
  (`coreop2d.py:652`):
  `a = Sec*DiffState(i) - Deg*q3d(i,1,3)`. The Sec source for
  cells with `d_i >= Set` is multiplied by `d_i`, ramping
  smoothly from `k_sec * Set` at threshold up to `k_sec` at
  `d_i = 1`.
- **Alternative code-bases:** humppa, tgrohens, C++ port all
  preserve.
- **Comparison-study question:** as for eq. 17 — does the
  smoother onset change the knot dynamics?

---

### `diff_accumulator` — kind: PaperVsCodeTension

**Question:** does eq. 6 differentiation accumulate Act (paper) or
Sec (FORTRAN)?

**Default:** `sec` per the v1 charter — this is the one
paper-vs-FORTRAN gap where v1 followed the FORTRAN; the comparison
study is what justifies revisiting.

#### Option `act`

- **Paper evidence:** main p. 585, eq. (6):
  `partial d_i / partial t = k_dff [Act]`. The accumulator is
  Act. This reading is consistent with the paper's prose
  description of Act as 'the activator that defines knot cells'.
- **FORTRAN evidence:** not used.
- **Alternative code-bases:** none.
- **Comparison-study question:** paper-faithful baseline. Worth
  flagging: under `act`, knot cells self-reinforce because Act
  peaks at knots; under `sec`, non-knot cells in the knot's
  neighbourhood accumulate towards differentiation because Sec
  is secreted by knots and diffuses outward. The two produce
  different spatial patterns of differentiation.

#### Option `sec`

- **Paper evidence:** none. The paper's eq. (6) writes `[Act]`,
  not `[Sec]`.
- **FORTRAN evidence:** `13.f90:524`
  (`coreop2d.py:671`):
  `DiffState(i) = DiffState(i) + Dff*(q3d(i,1,3))`. The third
  index `3` selects Sec, not Act. See finding
  `docs/findings/2026-05-04-differentiation-uses-sec-not-act.md`
  for the full reasoning. The 2014 paper does not specify which
  species drives `d_i`, so the FORTRAN's choice is unconstrained
  by the published text but was used to produce the 2014
  Extended Data Fig. 2 results.
- **Alternative code-bases:** humppa and tgrohens preserve. The
  C++ port matches (`tooth_model_diffusion.cpp:218` per the C++
  port review).
- **Comparison-study question:** is the paper-faithful Act form
  observably different on the seal example? The findings doc
  notes this is one of the rare paper-vs-FORTRAN gaps where v1
  followed FORTRAN without testing the alternative; the
  comparison study reports both.

---

### `knot_threshold_gate` — kind: PaperVsCodeTension

**Question:** does any cell with `[Act] >= 1` become a knot
(paper) or only cells with `i >= first_border_cell` (FORTRAN)?

**Default:** `none` per the paper.

#### Option `none`

- **Paper evidence:** main p. 586, 'Differentiation is
  irreversible': 'a cell becomes an enamel knot once `[Act] >= 1`
  in that cell'. The threshold value 1.0 is hard-coded and
  applies to all cells without geometric restriction.
- **FORTRAN evidence:** not used.
- **Alternative code-bases:** none.
- **Comparison-study question:** paper-faithful baseline.

#### Option `first_border_cell`

- **Paper evidence:** none. The paper has no interior/border
  split for knot formation.
- **FORTRAN evidence:** `13.f90:480`
  (`coreop2d.py:640`):
  `if (q3d(i,1,1)>1) then ; if (i>=first_border_cell) knots(i)=1`.
  Only cells with index `>= first_border_cell` (the inner ~2/3
  of the lattice in the FORTRAN's lattice numbering;
  the outer ~1/3 in silicoshark's centre-out numbering — see
  the A6 `first_border_cell` correction in the v2 progress
  log) become knots. The innermost cells are biologically odd
  to exclude from knot formation; this looks like an
  optimisation against premature knotting at the lattice
  centre, but the paper has no such restriction.
- **Alternative code-bases:** humppa, tgrohens, C++ port all
  preserve.
- **Comparison-study question:** how much of the FORTRAN's
  knot-pattern is driven by this gate?

---

### `knot_daughter_di` — kind: PaperAmbiguity

**Question:** when a knot mother divides, what `d_i` does the
(non-knot) daughter inherit?

**Default:** `zero_reset` per the paper.

#### Option `zero_reset`

- **Paper evidence:** main p. 585: 'If one of the mother cells
  is an enamel knot cell, the new cell is not a knot cell.'
  Reading 'is not a knot cell' as a hard reset of `d_i` to 0
  is one defensible reading of this prose.
- **FORTRAN evidence:** not used.
- **Alternative code-bases:** none. v1 silicoshark adopted this
  reading.
- **Comparison-study question:** paper-aligned baseline. Note
  the paper does not explicitly say what `d_i` value the
  daughter inherits — only that it is not a knot cell. This is
  one of the few PaperAmbiguity rows where the paper's silence
  is unequivocal: it does not specify, and either reset-to-zero
  or inherit-average is defensible.

#### Option `inherit_avg`

- **Paper evidence:** the paper says species concentrations of
  the new cell are 'the average of the two mother cells', and
  is silent on `d_i`. Reading the average rule as extending to
  `d_i` is one way to fill the silence.
- **FORTRAN evidence:** `coreop2d.add_cell` sets the daughter's
  `DiffState` to `0.5 * (DiffState[a] + DiffState[b])` per the
  add_cell logic in `13.f90:1034-1554`. A daughter of two knots
  inherits `d_i ≈ 1`, which means the daughter is frozen even
  though the paper says it is not a knot. The 2010 paper review
  flags this as a paper-vs-code gap (§'What the paper does NOT
  specify', point 10).
- **Alternative code-bases:** humppa and tgrohens preserve.
- **Comparison-study question:** does the inherited `d_i`
  change knot growth dynamics relative to the zero-reset?

---

### `mesenchyme` — kind: PaperVsCodeTension

**Question:** is the underlying mesenchymal layer modelled
explicitly?

**Default:** `per_column_z_layers` per the paper.

#### Option `absent`

- **Paper evidence:** none. The paper requires mesenchyme.
- **FORTRAN evidence:** none. The FORTRAN models mesenchyme.
- **Alternative code-bases:** v1 silicoshark used this as a
  fast-prototype shortcut. Acceptable for parameter sets where
  Sec stays 0 (e.g. seal.txt — Path A confirmed Sec stays at 0
  throughout the seal run); broken for parameter sets where
  Sec becomes non-zero.
- **Comparison-study question:** v1-style fast-prototype
  baseline. Should be exercised on the seal example only as a
  sanity check that v1 reproduces v2 in this regime.

#### Option `per_column_z_layers`

- **Paper evidence:** main p. 586: 'two layers of mesenchymal
  cells underneath' each epithelial cell. Inh and Sec diffuse
  vertically between layers; Act is epithelial-only.
- **FORTRAN evidence:** the `q3d(i, kk, k)` array in
  `13.f90:392-466` (`coreop2d.py:506-666`) carries
  per-cell-per-z-layer concentrations. The vertical diffusion
  between layers `kk-1, kk, kk+1` is computed via the `area_p`
  weighting; the substrate boundary layer at `kk =
  max_z_layers - 1` has the `0.044D1` sink coefficient.
- **Alternative code-bases:** humppa, tgrohens, the C++ port
  all model mesenchyme.
- **Comparison-study question:** baseline for any parameter set
  where mesenchymal Sec dynamics matter (most non-seal sets).

---

### `n_mesenchyme_layers` — kind: PaperVsCodeTension

**Question:** how many mesenchymal z-layers per epithelial
column?

**Default:** 2 per the paper.

#### Option `2`

- **Paper evidence:** main p. 586: 'two layers of mesenchymal
  cells underneath'.
- **FORTRAN evidence:** the FORTRAN allocates `max_z_layers = 4`
  (epithelium + 2 mesenchyme + substrate edge). The substrate
  edge is a boundary condition (sink with `0.044D1` coefficient
  per `13.f90:431, 444, 459`), not a real layer. So the FORTRAN
  agrees with the paper on the count of *modelled* mesenchyme
  layers (2), but adds a boundary-edge layer the paper does
  not.
- **Alternative code-bases:** humppa, tgrohens, C++ port all
  use `max_z_layers = 4`.
- **Comparison-study question:** baseline.

#### Option `int` (any positive integer)

- **Paper evidence:** none beyond the value 2.
- **FORTRAN evidence:** the FORTRAN's array dimension is
  parametric on `max_z_layers` but the value is hard-coded to 4
  in `allocate_initial_state` (`13.f90:135-266`). Increasing it
  would require re-tuning the substrate-sink coefficient.
- **Alternative code-bases:** none vary the count.
- **Comparison-study question:** does the count materially
  change Inh / Sec equilibration depth on parameter sets where
  Sec is non-zero? Inert when `mesenchyme = 'absent'`.

---

### `border_definition` — kind: PaperAmbiguity

**Question:** are border cells defined topologically (descendants
of initial border) or geometrically (low neighbour count)?

**Default:** `neighbour_count` because it composes cleanly with
the static-topology code path and matches the simpler reading of
the paper text.

#### Option `neighbour_count`

- **Paper evidence:** main p. 583, fig. 1 caption: 'cells with
  fewer epithelial neighbours than the interior count'. This
  is the geometric reading.
- **FORTRAN evidence:** indirect. The FORTRAN's
  `update_border_cells` (`coreop2d.py:1556-1628`,
  `13.f90:1556-1628`) promotes geometrically-emerging cells to
  border status, so the FORTRAN is partly geometric — but
  growth biases use the topological reading via
  `first_border_cell`.
- **Alternative code-bases:** none use a pure geometric reading.
- **Comparison-study question:** baseline for the
  `static_with_local_update` topology.

#### Option `topological_descendants`

- **Paper evidence:** main p. 585, 'Growth biases': 'The
  anterior border is defined as the border cells that in the
  initial condition are anterior to the centre (positive y
  values), and all of their descendants that lay in the tenth
  border. The posterior border is defined in the same way for
  negative values of y.' This is the topological reading and is
  load-bearing for the bias logic.
- **FORTRAN evidence:** the FORTRAN's bias logic at
  `coreop2d.py:938-946` (`13.f90:780-790`,
  `apply_border_bias`) uses `first_border_cell` (the index of
  the first non-interior cell in the lattice ordering) as the
  topological border marker. The `add_cell` logic
  (`coreop2d.py:1066-1070`) marks new cells as
  `new_cell_is_external` if both their parents were external,
  inheriting the topological border status across division.
- **Alternative code-bases:** humppa and tgrohens preserve. The
  C++ port at `tooth_model_division.cpp:749, 795` has the
  hard-coded `if (i == 2 || i == 5) continue` exclusions
  (faithful to FORTRAN, broken for any radius other than 3 per
  cpp-port-review §4).
- **Comparison-study question:** the A6 finding
  (`docs/findings/2026-05-05-path-b-v2-legacy-fortran-divergence.md`)
  shows this field is load-bearing for the seal example: under
  `neighbour_count`, the cervical-loop spreads to
  geometrically-emerging daughters and produces exponential
  over-division; under `topological_descendants`, the border
  set stays bounded (18 cells throughout) and the FORTRAN's
  trajectory is matched. This was the discovery that prompted
  wiring the field through to `forces.py` in A6, even though
  the field had existed since A1.

---

### `lattice_orientation` — kind: FortranAccident (load-bearing)

**Question:** is the initial hex lattice in axial layout (v1
default) or rotated 30° to match the FORTRAN's lattice
(`LEGACY_FORTRAN`)?

**Default:** `axial` to preserve the v1-byte-identity baseline.

#### Option `axial`

- **Paper evidence:** main p. 586, 'Initial conditions':
  'Seven epithelial cells in a hexagon (six surrounding one
  central cell), all at distance `|p_0| = 1`'. The paper does
  not specify the orientation of the hexagon. Either axial or
  rotated layout is consistent with the prose.
- **FORTRAN evidence:** not used.
- **Alternative code-bases:** v1 silicoshark used this layout.
  Outer-ring cells include `(rad, 0)` and `(-rad, 0)`; none at
  `(0, rad)`.
- **Comparison-study question:** baseline for v1-byte-identity
  and `PATH_B_DEFAULT`.

#### Option `fortran`

- **Paper evidence:** none. The paper does not specify
  orientation.
- **FORTRAN evidence:** `13.f90:135-266`
  (`coreop2d.allocate_initial_state`, FORTRAN
  `posarrad`-equivalent). The FORTRAN's hex lattice is rotated
  30° counter-clockwise relative to the axial convention; the
  outer ring includes cells at exactly `(0, rad)` and
  `(0, -rad)`. This was identified during the LEGACY_FORTRAN
  replication run (see A6 progress notes, residual divergence
  finding §'Cause' point 1).
- **Alternative code-bases:** humppa, tgrohens, the C++ port
  all use the FORTRAN orientation.
- **Comparison-study question:** required for
  `border_bias_x_zero_quirk` to have any effect on the seal
  example: the quirk excludes cells at `x exactly == 0` from
  the Bgr z-multiplier, anchoring those cells low in z and
  creating the gradient that drives the FORTRAN's 57-cell
  plateau and z=125 envelope. Under axial orientation no cells
  sit at x=0, the quirk is inert, and the lattice cannot
  reproduce the FORTRAN's z dynamics. Field added in A6 once
  this lattice-rotation dependency was identified.

---

### `border_bias_x_zero_quirk` — kind: FortranAccident (load-bearing)

**Question:** preserve the FORTRAN's `update_cell_position`
behaviour where cells with x exactly `== 0` don't receive
Pbi/Abi/Bgr multipliers?

**Default:** `True` to preserve compatibility with the FORTRAN
goldens. The corrected behaviour (`False`) is biologically
cleaner but produces different numerics on the seal example.

#### Option `True`

- **Paper evidence:** none. The paper does not describe the
  multiplier scheme at this level of detail. Main p. 585
  ('Growth biases') describes the Pbi/Abi multipliers
  conceptually but does not specify the half-plane chain.
- **FORTRAN evidence:** `coreop2d.py:1009-1031`
  (`13.f90:1009-1031`, `update_cell_position`). The
  per-half-plane bias logic uses `if x > 0` and `elif x < 0`,
  which omits cells at exactly `x == 0`. For the seal initial
  condition with y-axis cells at `x == 0` exactly (under
  `lattice_orientation = 'fortran'`), this means those cells
  don't receive the Bgr z-amplification. The chain-of-accidents
  docstring at the top of `coreop2d.py` documents this as
  load-bearing for the 57-cell plateau and z=125 envelope.
- **Alternative code-bases:** humppa and tgrohens preserve.
  The C++ port's `updatePositions` is in the missing
  `tooth_model_mechanics.cpp` per cpp-port-review §2(b), so I
  cannot confirm; the C++ review notes the editor's spot-check
  of `tooth_model_mechanics.cpp:537-545` confirms the nesting
  matches FORTRAN.
- **Comparison-study question:** the FORTRAN-faithful baseline
  for `LEGACY_FORTRAN`.

#### Option `False`

- **Paper evidence:** none.
- **FORTRAN evidence:** none — this option corrects the
  half-plane logic to include `x == 0` cells in one of the
  half-planes (currently `x >= 0` per the A4 progress notes).
- **Alternative code-bases:** none correct.
- **Comparison-study question:** does the corrected behaviour
  change the seal-example trajectory? On the seal example with
  `lattice_orientation = 'fortran'`, removing the quirk
  uniformly lifts the y-axis cells in z, eliminating the
  gradient that produces the 57-cell plateau. This is one of
  the few rows where I expect a large biological change, with
  the corrected behaviour being arguably the right answer
  scientifically but losing FORTRAN compatibility. The
  comparison study quantifies the difference.

## Cross-field interactions

A handful of field combinations are notable enough to flag.

The combination `update_order = 'jacobi'` ×
`topology = 'static_with_local_update'` is the cleanest
numpy-idiomatic preset and corresponds to `PATH_B_DEFAULT`. Pure
Jacobi composes cleanly with the static graph because there is no
intermediate-state ordering to maintain across daughter insertions.

The combination
`lattice_orientation = 'fortran'` ×
`border_bias_x_zero_quirk = True` is required for
LEGACY_FORTRAN to reproduce the seal-example z dynamics. Either
field alone is insufficient: without the rotation, no cells sit at
x=0 and the quirk has no targets; without the quirk, the y-axis
cells receive the same Bgr multiplier as everyone else and the
lattice lifts uniformly. The interaction is documented in the
`lattice_orientation` field's docstring and in the residual
divergence finding.

`LEGACY_FORTRAN` bundles every FORTRAN-side choice plus the
phenomenological `division_total_cap = 60`. Each FORTRAN-flavoured
option has been A/B-tested individually during A6 (see progress
notes); the bundle reproduces the FORTRAN goldens within v2
tolerance modulo the residual z_min divergence.

`PAPER_LITERAL_2010` differs from `PAPER_2010` in
`eq14_denominator` only. Comparing the two on the seal example
isolates the biological cost of the eq. 14 typo. If the
trajectories are indistinguishable, the typo is harmless; if they
diverge, the typo would have produced different published results
had the paper text been used as the implementation.

`HUMPPA_LITERAL` is currently identical to `LEGACY_FORTRAN`
because all the 13.f90-only divergences catalogued in
`docs/research/13f90-vs-humppa-divergences.md` live in
machinery v2 either replaces or has not yet wired through.
The Swi/parap[6] handling that genuinely distinguishes humppa
from 13.f90 will need a new Discretisation field once the band
becomes observable.

## What the comparison study tests

| Preset | Question answered | Key fields differing from `PATH_B_DEFAULT` |
|---|---|---|
| `PATH_B_DEFAULT` | What does the cleanest numpy-idiomatic implementation produce? | (defaults) |
| `PAPER_2010` | What does the 2010 paper as written (with eq. 14 typo corrected) produce? | (defaults — same as `PATH_B_DEFAULT` initially) |
| `PAPER_LITERAL_2010` | Does the eq. 14 typo materially change biology? | `eq14_denominator='act_typo'` |
| `LEGACY_FORTRAN` | Does the v2 framework reproduce the FORTRAN goldens? | every FORTRAN-flavoured option, plus `division_total_cap=60`, `lattice_orientation='fortran'`, `border_definition='topological_descendants'`, `border_bias_x_zero_quirk=True` |
| `HUMPPA_LITERAL` | Are the 13.f90-vs-humppa divergences observable in the v2 framework? | currently identical to `LEGACY_FORTRAN`; will diverge once Swi-band-style fields are added |

The comparison-matrix runner in
`scripts/run-discretisation-study.py` exercises each preset on the
seal example plus the 2014 EDF Act-sweep and Inh/Di-lattice
validation targets, capturing cusp count, cusp width, cell-count
plateau, vertex envelope, plateau trajectory, and qualitative
regime per cell of the matrix. The audit trail above grounds each
cell's reported behaviour in a paper passage, a FORTRAN line, and
where applicable an alternative-code-base citation, so a reader
who asks 'why does cell X look like that?' can trace the answer
back to a specific implementer-choice point.
