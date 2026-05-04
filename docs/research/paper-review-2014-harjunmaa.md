---
title: 'Paper review — Harjunmaa et al. 2014 (tribosphenic extension)'
author: Lyndon Drake (with Claude Code)
date: 2026-05-04
---

## Aim

Review the 2014 *Nature* paper by Harjunmaa, Seidel, Häkkinen et al.,
'Replaying evolutionary transitions from the dental fossil record', and
its Supplementary Information, focusing on what is *new* relative to the
2010 Salazar-Ciudad and Jernvall paper. The 2014 paper is the
tribosphenic extension that produced `humppa_translate.f90`, the
ancestor of `13.f90`. For Path B I want to know which model decisions
this paper introduces, which it leaves unspecified, and which parameter
sets the paper provides as candidate validation targets.

## Headline finding

The 2014 paper introduces no new equations and defines no new model
mechanisms in its published text. It is an *experimental* paper that
uses the 2010 model as a black box to interpret EDA dosage and SHH
inhibition results. The Methods section names the three signalling
species, points the reader to ref. 21 (the 2010 paper), and otherwise
describes only laboratory technique. The Supplementary Information is
two pages of experimental cusp tabulation (Supplementary Table 1) with
no equations and no parameter table at all (SI pp. 1–2).

This is consequential for Path B: every equation Path B implements
must come from the 2010 paper. The 2014 paper contributes (a) a single
fully-specified parameter set for the wild-type tribosphenic mouse
(Extended Data Fig. 2), (b) a parameter sweep map (Extended Data
Fig. 3) that names which parameters move which morphological features,
and (c) candidate Eda-null and SHH-inhibited parameter combinations
(Extended Data Figs 4–5) which can serve as validation goldens.

The detailed model decisions that Path B has to make — initial cell
layout, asymmetry handling, knot placement rules, what 'Tribosphenic
tooth' means as a model variant in the GUI — are present in
`humppa_translate.f90` but not in the paper. They are implementation
choices that the FORTRAN authors filled in by judgement.

## What is new in the model (text claims)

### Three diffusible signals named explicitly

The Methods section ('Modelling') states the model has three
diffusible signals: 'an activator inducing enamel knots, an enamel
knot-secreted inhibitor of enamel knot formation, and a growth factor
regulating growth of the epithelium and the mesenchyme' (Methods,
p. S1, paragraph 'Modelling'). This matches the 2010 Nature paper's
Act / Inh / Fgf species.

The paper does not state what 'Sec' (the FORTRAN's growth-factor
species `q3d(i,1,3)`) corresponds to in molecular terms beyond
'growth factor'. It does not introduce a fourth or fifth signalling
species. The FORTRAN's `q3d` array has five slots (`act, inh, fgf, ect, p`
per `cpp_README.md` line 56), so two slots — `ect` (ectodin) and `p`
(an unnamed quantity) — are not equations the 2014 paper documents.
The 2010 paper or the FORTRAN itself must be the source for those.

### EDA signalling treated as parameter modulation, not a fourth species

EDA is the *experimental* lever. The paper says: 'EDA signalling is
known to interact with multiple signalling pathways, and accordingly
our focus was to model the effects of EDA through the genetic
parameters' (Methods, p. S1). EDA does not enter the model as a
distinct diffusible quantity. Instead, varying EDA dose in the lab is
mapped onto varying *one or more* of the existing 2010 parameters in
silico. This is consistent with the FORTRAN having no `eda`
concentration in its `q3d` array.

### SHH inhibition treated similarly

SHH inhibition (the second experimental lever, used to recover
basal-rodent characters) is also not modelled as a distinct chemical
species. The paper explicitly equates it with reducing 'cusp spacing'
(p. 47, 'Retrieving ancestral character states') and operationalises
this in the model as 'reducing activator inhibition by inhibitor
(Inh) or diffusion of inhibitor (Di)' (Extended Data Fig. 5 caption,
p. ED-5). So SHH-low corresponds to lower `Inh` and/or lower `Di` in
the 2010 model's Eq. 2 inhibitor dynamics.

### Bifurcation-like irreversible knot differentiation

Methods states: 'Because the enamel knot differentiation is
irreversible, the model represents an irreversible bifurcation-like
model' (p. S1). This is an interpretive label rather than an equation.
The irreversibility is implemented in the FORTRAN by the `q2d(i,1)==1`
gate in `diferenciacio` and by `knots(i)` flags being set without
reset; the 2014 paper does not redefine the differentiation criterion.

### A new parameter named in the GUI: `Pbi` for posterior bias

Extended Data Fig. 2 shows the ToothMaker GUI 'Tribosphenic tooth'
model. The bias panel contains four directional bias parameters —
`Abi`, `Pbi`, `Lbi`, `Bbi` (anterior, posterior, lingual, buccal) —
plus a `Bgr` ('bias factor') and a `Bwi` ('bias width', actually the
AP-bias centre radius per `cpp_README.md` lines 87–88, 102–104). The
2010 paper's bias system used `BA` for buccal in symmetric mouse
geometry; the four-direction split is needed to break tribosphenic
mesiodistal asymmetry (anteroconid/protoconid vs. talonid). The
paper does not describe the four-bias mechanism in its text. The
mechanism lives in `humppa_translate.f90:1044–1056`: when a border
cell's `|y| < radibii`, its x-velocity is multiplied by `bia` (if x>0)
or `bip` (if x<0), and its z-velocity is multiplied by `fac` (= Bgr).
Path B will have to reproduce this bias geometry from the FORTRAN
source; the paper does not specify it.

## New parameters introduced (or first published values)

The wild-type tribosphenic mouse parameter values appear visually in
Extended Data Fig. 2 (the GUI screenshot). I list them here with
FORTRAN-name and `cpp_README.md` mapping. Values are read from the
GUI controls in the screenshot; some are partially obscured.

| GUI label | Fortran | Role | Wild-type value (Ext Data Fig. 2) |
|-----------|---------|------|-----------------------------------|
| Egr  | `tacre`  | Epithelial proliferation rate | 0.013 |
| Mgr  | `tahor`  | Mesenchymal proliferation rate | 16000 |
| Rep  | `elas`   | Repulsion (Young's modulus) | 1.5 |
| Adh  | `crema`  | Neighbour traction | 0.005 |
| Act  | `acac`   | Activator auto-activation | 1.6 |
| Inh  | `ihac`   | Inhibition of activator | 800.0 |
| Sec  | `ih`     | Growth factor secretion | 0.14 |
| Da   | (Act diffusion)  | Activator diffusion | 0.2 |
| Di   | (Inh diffusion)  | Inhibitor diffusion | 0.2 |
| Ds   | (Sec diffusion)  | Growth factor diffusion | 1.0 |
| Int  | `us`     | Inhibitor production threshold | 0.19 |
| Set  | `ud`     | Growth factor production threshold | 0.05 |
| Boy  | `difq2d(2)` | Buoyancy / stellate strength | 0.17 |
| Dff  | `tadif`  | Differentiation rate | 0.00013 |
| Bgr  | `fac`    | Bias factor (z-component on AP bias) | 3.4 |
| Abi  | `bia`    | Anterior bias (x-multiplier when x>0) | 3.2 |
| Pbi  | `bip`    | Posterior bias (x-multiplier when x<0) | 7.0 |
| Bbi  | `bib`    | Buccal bias (BMP4 floor at buccal border) | 1.34 |
| Lbi  | `bil`    | Lingual bias (BMP4 floor at lingual border) | 1.31 |
| Rad  | `radi`   | Initial grid radius | 3.0 |
| Deg  | `mu`     | Protein degradation | 0.076 |
| Dgr  | `tazmax` | Sharpness maximum | 10500.0 |
| Ntr  | `radibi` | Border-to-nucleus traction | 3.5 |
| Bwi  | `radibii`| AP-bias centre radius | 1.32 |
| Ina  | `ina`    | Initial activator | 0.0 (wild-type — see caption) |

The Extended Data Fig. 2 caption explicitly notes 'Initial activator
concentration (Ina) is not used in the model' for the published
tribosphenic wild type. This is at odds with `13.f90` examples like
`seal.txt` which set `ina = 0.5`. Read the paper's claim narrowly:
*for this particular wild-type-mouse parameter set*, Ina is zero. The
mechanism still exists in the model.

`Pbi` (posterior bias = 7.0) being more than twice `Abi` (3.2)
encodes the mesiodistal asymmetry of the mouse first lower molar:
the tooth grows faster posteriorly, producing the talonid shelf.
This is the parameter setting that distinguishes the 'Tribosphenic
tooth' GUI variant from the 2010 'Mouse tooth' variant. The paper
does not say this in words; it is recoverable only by reading the
GUI image.

### Two parameters absent from Extended Data Fig. 2 but present in `13.f90`

- `Swi` (`tadi`): border definition distance. The GUI image shows a
  `Swi` field but the value digit is partially cut off in the
  screenshot. The seal example file has `Swi = 0`. Per the Path B
  charter, `13.f90`'s setparams silently drops Swi anyway.
- `umgr` (`uMgr`): basal mesenchymal rate (Sec-independent). Not in
  the GUI image. Present in seal example as `umgr = 0`. Likely
  introduced post-2014.

## What the paper says about reproducing fossil tooth morphologies

This is the most useful contribution for Path B validation. The paper
provides three layers of parameter→morphology mapping.

### Layer 1: single-parameter sweeps (Extended Data Fig. 3)

Methods 'Modelling' (p. S1, last paragraph) and Extended Data Fig. 3
report a 14,000-iteration scan of each of nine parameters individually,
varied at 10% intervals up to 90% departure from the wild-type mouse,
to map cusp count as a function of each parameter:

- `Act-` (decreasing) — reduces cusps; only single parameter that on
  its own can produce a single-cusped morphology.
- `Da-`, `Int-`, `Inh+`, `Di-`, `Set-`, `Sec-`, `Ds+`, `Dff+`, `Egr-`,
  `Mgr-` — all decrease cusp number when moved in the noted direction.

Extended Data Fig. 3b shows the resulting tooth shape for each
parameter at its cusp-decreasing extreme. These are not full goldens
(no parameter file is published) but they are *qualitative*
predictions a Path B implementation should reproduce: e.g., reducing
`Act` from 1.6 to ~0.5 should give a single-cusped tooth.

### Layer 2: Act sweep at 0.1 increments (Extended Data Fig. 4)

Extended Data Fig. 4a shows simulated shapes at `Act` ∈ {0.1, 0.2,
0.3, …, 1.6} (sixteen values, all other parameters at wild-type). The
caption gives a clean monotonic relationship: low Act → single tip,
intermediate → 2 cusps near 0.5, then 3, then 4, asymptoting to 5 at
Act = 1.6. Extended Data Fig. 4b plots inhibitor-domain size against
cusp count and shows a near-linear trend up to ~3,000 'arbitrary
units' for cusp count 5.

This is the cleanest single-axis validation for Path B: hold seal
example or wild-type-mouse parameters fixed, sweep `Act`, count
cusps, and verify the same monotone progression.

### Layer 3: Eda-null and SHH-inhibited parameter combinations (Extended Data Fig. 5)

Extended Data Fig. 5 gives concrete parameter sets that the paper
claims produce specific tribosphenic recoveries from an Eda-null
baseline. The figure shows four panels with fully-specified
two-parameter changes:

| Modification | Inh | Di | Phenotype claim |
|--------------|-----|----|----|
| Baseline (Eda null simulated) | 800 | 0.2 | Single-cusped Eda-null phenotype |
| Reduce inhibition strength | 20 | 0.2 | Multi-cusp formation |
| Reduce inhibition strength further | 10 | 0.2 | Multi-cusp formation |
| Reduce inhibitor diffusion strongly | 800 | 0.01 | Multi-cusp formation |
| Reduce inhibitor diffusion further | 800 | 0.0001 | Multi-cusp formation |

These five points define a low-dimensional validation lattice. If
Path B reproduces (a) Eda-null single cusp at the baseline parameters,
and (b) multi-cusp at any of the four perturbation points, that is
strong evidence the inhibitor dynamics are correctly implemented.

### Layer 4: cusp-number tabulation in Extended Data Table 1

Extended Data Table 1 reports anterior cusps, posterior cusps, and
talonid height for *Eda* null at EDA doses 0, 2.5, 10, 25, 50, 100,
500, 1000 ng/ml plus wild type and *Tribosphenomys minutus*. These
numbers are experimental, not modelled outputs, but they map onto the
in-silico Act sweep (Eda dose 0 ≈ Act 0.1; wild type ≈ Act 1.0–1.6).
For Path B this is a downstream validation: if the Act sweep gives
the expected cusp counts at each ng/ml dose, the model and the
experiment agree. Wild-type values: 3.20 anterior, 2.50 posterior,
talonid height 0.97. Eda-null: 1.54 anterior, 1.00 posterior,
talonid height 0.43.

### What no parameter set looks like

The paper does *not* provide explicit parameter files for the seal,
lemur, vole, or non-mouse morphologies. The validation set is
constrained to: (a) the wild-type mouse Tribosphenic parameter set
(Ext Data Fig. 2), (b) the Act sweep across 16 values, (c) the
Inh/Di lattice in Ext Data Fig. 5. The seal example in `13.f90`
parameters comes from a different lineage (likely Zimm 2023 / shark
fork or hand-tuned by users) and is not a paper golden.

## Implementation details specific to tribosphenic teeth

The paper does not describe these. They are recoverable only from
`humppa_translate.f90`. The most consequential are:

1. **Initial cell layout: hexagonal grid with `Rad = 3`.** This
   produces 19 cells (1 + 6 + 12) in the wild-type tribosphenic GUI
   image. `posarrad` in the FORTRAN constructs concentric rings of six
   cells each.
2. **Mesiodistal asymmetry via two-sided x-bias.** `humppa_translate.f90:
   1044–1056` (subroutine `actualitza`): for any cell with `|y| <
   radibii` (i.e. inside the AP-bias band), if `x > 0` multiply
   x-velocity by `bia` (Abi); if `x < 0` multiply by `bip` (Pbi). z-
   velocity is multiplied by `fac` (Bgr) regardless of sign. With the
   wild-type values `Abi = 3.2`, `Pbi = 7.0`, the posterior side grows
   ~2× faster in x than the anterior side, producing the talonid.
3. **Lingual–buccal floor on activator concentration.** `biaixbl`
   (line 966 onwards): for border cells (which the routine identifies
   by some criterion the paper does not state), if `q3d(i,1,1) < bil`
   (lingual floor) or `< bib` (buccal floor), reset to that floor.
   This is the BMP4-as-bias mechanism; it is described in the 2010
   paper conceptually but the floor-replacement form (rather than
   additive forcing) is a FORTRAN choice.
4. **Sec-independent basal mesenchymal rate (`umgr`).** Not in
   tribosphenic wild type per Extended Data Fig. 2 but present in
   `13.f90`. A post-2014 addition.
5. **`Tribosphenic tooth` as a named GUI model variant.** The
   ToothMaker GUI dropdown ('Model: Tribosphenic tooth' in Ext Data
   Fig. 2) implies there is a separate 'Mouse tooth' or 2010 model
   variant in the same binary. The two variants likely differ only in
   parameter defaults and possibly in whether the four-direction bias
   panel is exposed. The paper does not document the variant
   distinction; the FORTRAN source has not been audited for whether
   any code paths are actually variant-gated.

## What the paper does NOT specify that `13.f90` implements

I list these as Path B decision points. Each is a place where the
FORTRAN authors chose an implementation that the paper does not
constrain.

1. **Triangulation method.** The paper shows triangulated meshes in
   every figure but never says how cells are triangulated. The
   FORTRAN's `calculmarges` and topology walk implement a particular
   neighbour-cycling convention. Path B can choose a different (e.g.
   Voronoi-based or Delaunay) approach without violating the paper.
2. **Cell-division geometry.** When does a cell divide, and where do
   the two daughters land? Paper says only that growth happens via
   proliferation rates Egr and Mgr. The FORTRAN's `afegircel` uses
   `dmax = 2.0` as the division-trigger distance and a specific
   placement rule using neighbour topology. The paper does not state
   either.
3. **Diffusion discretisation.** The paper writes the reaction-
   diffusion as continuous equations (in the 2010 paper). The FORTRAN
   uses a per-cell sum-over-neighbours discretisation with the
   `0.044D1` boundary-sink coefficient (per `cpp_README.md` line 252).
   The 0.44 magic number is not in the paper.
4. **Differentiation criterion exact form.** Paper says enamel-knot
   differentiation is 'irreversible' and triggered when activator
   exceeds threshold; FORTRAN sets `q2d(i,1) = q2d(i,1) + tadif *
   q3d(i,1,3)` (line 659) accumulating 'differentiation' as an
   integral of the *third* `q3d` slot (Sec/growth factor), not the
   activator itself, and treats `q2d(i,1) >= 1` as the irreversibility
   gate. This is paper-counter-intuitive — knots accumulate based on
   growth factor, not directly on the activator that defines them.
   Path B has to decide whether to follow the FORTRAN or the natural
   reading of the paper.
5. **Stellate-reticulum / buoyancy mechanism.** GUI exposes `Boy`;
   the FORTRAN implements it in `stelate`. The paper does not mention
   stellate reticulum at all. Path B can omit `Boy` (set to wild-type
   0.17 or just zero) for the seal example without violating the
   paper.
6. **'AP-bias centre radius' (`Bwi`/`radibii`).** The bias only
   applies to cells with `|y| < radibii`. Wild-type value 1.32, with
   `la = 1.0` inter-cell distance and `Rad = 3` initial layout, means
   a band one cell wide at y=0. Why 1.32 and not 1.0? Not in the
   paper.
7. **Knot placement rules in initial conditions.** Paper says knots
   form where activator exceeds threshold. FORTRAN's `initact` simply
   sets `q3d(i,1,1) = ina` for all cells uniformly. With `Ina = 0` for
   wild-type tribosphenic, the initial activator field is zero, and
   knots only emerge from the lingual/buccal BMP4 floor seeding the
   border cells. This is not stated in the paper.
8. **Time-step size and total iteration count.** Paper says the
   simulation 'accounts for mouse molar development up to day 16–19'
   but specifies only that all parameter scans run for 14,000
   iterations (Ext Data Fig. 3 caption). The FORTRAN's `delta = 0.05`
   is not in the paper. The 14,000-iteration figure is a choice, not
   a derivation.

## Path B implications

### Decisions Path B faces that the 2010 paper alone does not answer

The 2014 paper is helpful here because it confirms which parameters
are exposed by the GUI of a Jernvall-lab-blessed implementation. Path
B should expose at least the same set:

1. **Adopt the four-direction bias parameter set.** Path B's
   parameter loader needs `Abi`, `Pbi`, `Bbi`, `Lbi`, `Bgr`, `Bwi`.
   The FORTRAN's two-sided x-multiplier scheme (lines 1044–1056) is
   the implementation pattern.
2. **Adopt `Tribosphenic tooth` as the only model variant.** Do not
   try to support a 'Mouse tooth' 2010 variant simultaneously. The
   parameter file format is shared; the variant difference lives in
   the parameter values.
3. **Decide whether to adopt the 'differentiation accumulates Sec,
   not Act' choice.** This is documented in `humppa_translate.f90:659`
   and is a paper-vs-FORTRAN gap worth a `docs/findings/` entry. My
   recommendation: implement the FORTRAN form, mark it explicitly in
   the docstring, and note the gap.
4. **Decide on triangulation approach.** Path B can use a Delaunay or
   Voronoi triangulation of cell positions per timestep, computed
   in NumPy via SciPy's `Delaunay`, instead of maintaining the
   FORTRAN's running neighbour list. This is the place where
   vectorisation buys most.
5. **Decide on the irreversibility gate.** Once a cell has `q2d(i,1)
   >= 1`, its contribution to position updates is suppressed (knots
   sink, no growth). The form of the gate is in
   `humppa_translate.f90:1062–1064`. Path B should match.

### Validation goldens this paper supports

In addition to the existing `tests/golden_fortran/seal.off` golden,
Path B can use the following from the 2014 paper as validation
targets. None of these are byte-equality goldens; all are
*qualitative* or *low-dimensional* checks.

1. **Wild-type tribosphenic mouse parameter set** (Ext Data Fig. 2):
   the values listed in the table above. Expected output: 5–6 cusps
   on a tribosphenic profile, talonid height ≈ 1.0 relative to
   trigonid (Ext Data Table 1, Mus musculus row: 4 anterior, 3
   posterior, talonid 1.00). Capture from `13.f90` once and add as
   `tests/golden_fortran/wt_mouse.off`.
2. **Eda-null parameter set**: same as wild-type but `Act = 0.1`
   (per Ext Data Fig. 4a's lowest value). Expected output: 1 cusp.
   Capture as `tests/golden_fortran/eda_null.off`.
3. **Act sweep**: 16 parameter sets at `Act` ∈ {0.1, 0.2, …, 1.6}.
   Expected: monotonic increase in cusp count from 1 (at 0.1) to 5
   (at 1.6), matching Ext Data Fig. 4b. Implement as a
   parameterised pytest, not 16 OFF goldens.
4. **Inh/Di lattice from Ext Data Fig. 5**: five (Inh, Di) points
   atop the Eda-null baseline. Expected: single-cusp baseline,
   multi-cusp at the four perturbations. Implement as a
   parameterised pytest with cusp-count assertions.
5. **Single-parameter sweep cusp predictions** from Ext Data Fig. 3:
   nine parameters, each scanned ±10–90% of wild-type, expected
   monotone direction of cusp-count change. Implement as a small
   regression test on the cusp-count function.
6. **Cusp number ↔ primary-knot-area scaling** (Fig. 1b, Extended
   Data Fig. 1). Reduced major axis regression slope 0.0533, r² =
   0.613 across 46 simulated teeth. A statistical sanity check Path
   B can run by sweeping any single parameter that affects cusp count.

### What this paper rules out as validation source

- A seal-specific parameter set published by the authors. There
  isn't one.
- Per-fossil parameter sets (lemur, vole, *Tribosphenomys*
  recovered). The paper *experiments* with SHH inhibition on real
  mouse molars to retrieve *Tribosphenomys*-like morphology in vivo
  but never gives the in-silico parameter set used to model that
  recovery. Ext Data Fig. 5's lattice is the closest the paper comes,
  and it is for cusp number not full morphology.
- Equations not in the 2010 paper. There are none new.

### Summary table of validation goldens grouped by source paper

| Source | Identifier | Validation target |
|--------|------------|-------------------|
| `seal.txt` (project) | `seal.off` | Existing FORTRAN-Python cross-validation |
| 2014 Ext Data Fig. 2 | `wt_mouse.off` (new) | 5–6 cusps, tribosphenic profile |
| 2014 Ext Data Fig. 4a | Act sweep | Monotone cusp-count progression |
| 2014 Ext Data Fig. 5 | Inh/Di lattice | Eda-null + four perturbations |
| 2014 Ext Data Fig. 3 | Single-param sweeps | Direction of cusp-count change |
| `nosea_10261__0.txt` (charter §Validation) | `nosea_*.off` | Non-seal regime, Rad=2 |

## Closing observation

That a 2014 *Nature* paper offers no new model equations is itself a
finding worth recording. The 2014 contribution is *experimental and
parametric*, not *theoretical*. Path B should implement the 2010
equations, parameterise them per Extended Data Fig. 2's wild-type
tribosphenic mouse values, and treat the 2014 paper as a source of
qualitative validation targets and as the canonical specification of
which parameters the published GUI exposes. The model decisions that
distinguish 'Tribosphenic tooth' from the bare 2010 model live in
`humppa_translate.f90`, not in the paper, and Path B will need to
reproduce or replace them by reading the FORTRAN source directly.
