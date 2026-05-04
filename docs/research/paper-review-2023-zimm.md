---
title: 'Paper review — Zimm et al. 2023 PNAS (shark ToothMaker extension)'
author: Lyndon Drake (with Claude Code)
date: 2026-05-04
---

## Scope and aim

The 2023 Zimm et al. PNAS paper [@zimm2023; DOI 10.1073/pnas.2216959120]
extends the Salazar-Ciudad and Jernvall 2010 ToothMaker model to
*Scyliorhinus canicula* (small-spotted catshark) odontogenesis. It is
**out of scope for Path B**, which targets the 2010 mammalian model
plus the 2014 tribosphenic extension. I review it here for
transferable presentation conventions, not for biology. The companion
code at `.tmp/silicoshark/toothmodel.f90` is a sibling of `13.f90`'s
ancestor, not the same model.

I read the SI in full (21 pp.) and pages 1–3 of the main article. I
sampled but did not exhaustively read main-article pages 4–10
(Results/Discussion, mostly biology).

The headline finding for Path B is that the SI is **less useful than
hoped** for clean conventions, and that several of its most
prominent presentation choices reflect shark-specific changes rather
than the canonical 2010 model. A handful of small things are worth
borrowing; most are not.

## Parameter naming and the SI parameter table

The SI provides Table S1 ('List of developmental parameters', p. 9 of
21) and Table S2 ('Limits of developmental parameter values', p. 10 of
21). Table S1 uses a four-column layout:

> `Parameter | Type | Symbol | Description`

Twenty parameters are listed. The **Type** column groups them as
`diffusion rate`, `genetic`, `other genetic`, `threshold`,
`growth rate`, `mechanical`, and `biases`. This grouping is genuinely
useful: it maps onto the same partition Path B will need internally
(diffusion arrays, GRN reaction terms, mechanical forces, boundary
biases). I would borrow the type-column structure, but **not** the
parameter list itself, because Zimm has substituted shark-specific
names (`S_1`, `S_2`, `S_3`) for the mammalian morphogen names (Fgf,
Bmp, Shh).

The naming substitution is explicit in the SI footnote to Table S1
(p. 9, lines under the table):

> 'The placeholder variables i,j stand for the morphogen kernel
> members S_1 (analogous to Fgf), S_2 (analogous to Bmp), S_3
> (analogous to Shh) and Sec (only relevant for diffusion).'

The SI text repeats this substitution at lines 53–55 of p. 2 of 21
('we use S_1, S_2 and S_3 to denote the signals in our model whose
dynamics are compared, respectively, to Fgf, Bmp and Shh'). Path B
uses the mammalian model and should keep `Fgf`, `Bmp`, `Shh` (or
abstract `act`, `inh`, `sec` — see comparison with `multicusp.ini`
below). The S-subscript notation actively obscures the
correspondence with `13.f90`'s field arrays and would harm
maintainability.

The shorter parameter symbols **are** worth borrowing where they map
1:1 onto the 2010 names already used in `coreop2d.py` and the
`multicusp.ini` glossary:

| 2023 SI symbol (Table S1) | `multicusp.ini` label | 2010 paper use |
|---------------------------|-----------------------|---------------------|
| `D_i` | 'Diffusion rates' | diffusion rate per morphogen |
| `k_{ij}` | (row labels FgfFgf, FgfBmp, BmpFgf etc.) | reaction-coupling matrix entries |
| `a_i` | (no direct label; `A(i)` row 29) | upstream activation rate |
| `Sec` | 'Secretion' | secondary activator production |
| `Diff` | 'Differentiation' | differentiation rate |
| `Deg` (μ) | 'Degradation' | degradation rate |
| `T(S_2)`, `T(Sec)`, `T(S_3)` | 'Bmp threshold', 'Sec threshold', 'Shh threshold' | activation thresholds |
| `Egr`, `Mgr`, `Egr_i`, `Dgr` | 'Epithelial and mesenchymal growth', 'Downgrowth' | epithelial, mesenchymal, morphogen-induced, downgrowth rates |
| `Boy`, `Adh`, `Rep`, `Ntr` | 'Buoyancy', 'Adh traction', 'Repulsion', 'Nuclear traction' | mechanical parameters |
| `abia`, `pbia`, `zbia`, `rbia` | 'AP biases', 'z force bias', 'bias radius' | biases |

The lower-case mechanical symbols (`Egr`, `Mgr`, `Boy`, `Adh`, `Rep`,
`Ntr`) are stable across the 2010 paper, the 2014 paper, the
`multicusp.ini` file, and the 2023 SI. These are the names Path B
should adopt verbatim. The Greek `μ` for degradation (Table S1 row 6,
'Deg; μ') and the Greek `θ` (`ϑ`) for thresholds (eqs. 1–4) are also
canonical.

**Useful from Table S2.** Table S2 gives min/max for each parameter
in a simple two-column layout. It is a sanity check on plausible
ranges and is one place worth borrowing the layout — but the
**values** are shark-fitted. Path B's parameter loader can record
expected ranges in similar form, but the numbers should come from
2010/2014 supplementary material (or from `multicusp.ini` for a
seal/no-sea baseline), not from Table S2.

## Cross-check against `multicusp.ini`

The `.tmp/silicoshark/multicusp.ini` file matches the SI's parameter
groupings reasonably well, with a few naming mismatches that confirm
the SI is the right reference for symbols but not for the row order.

The `.ini` file's order — diffusion rates, reaction-coupling rates
(FgfFgf, FgfBmp, BmpFgf), Degradation, Differentiation, Secretion,
'Interactions Shh' (4 numbers), Bmp threshold, Sec threshold, then
mechanical block, then biases, then 'Further additions' (z force
bias, bias radius, Egr(i), A(i), Shh threshold) — is not the SI Table
S1 order. The 'Further additions' block is the giveaway: parameters
that were tacked on for the shark extension (`Egr(i)`, the
per-morphogen growth-induction rates; `A(i)`, the upstream activation
vector; the Shh threshold) are appended rather than slotted into the
table. This confirms two things:

1. The SI's Table S1 is the right reference for parameter symbols
   and groupings.
2. The companion `multicusp.ini` is a hybrid: 2010 baseline + 2023
   shark additions. Path B should not use it as a parameter
   template; it has shark biases (`zbia=1.005`, `rbia=0.8`) and the
   `Egr(i)` and `A(i)` vectors that do not exist in the 2010 model.

## Equation typesetting

Equations 1–4 (SI p. 3 of 21) give the morphogen dynamics:

- Eq. 1: dS_1/dt = (Hill self-activation) + a_1 − μ[S_1] − D_{S_1}∇²[S_1]
- Eq. 2: dS_2/dt = piecewise on [S_1] threshold ϑ_{S_1}
- Eq. 3: dS_3/dt = piecewise on [S_2] threshold ϑ_{S_2}
- Eq. 4: dSec/dt = piecewise on [S_1] threshold ϑ_{Sec}

The typesetting is conventional: bracket notation `[S_i]` for
concentrations, partial-derivative `δ/δt` notation, Laplacian `∇²` for
diffusion, piecewise `\begin{cases}` for threshold-triggered
behaviour. Two observations:

1. **Cleaner than the 2010 paper for the diffusion term.** The 2023
   SI consistently writes `−D_{X}∇²[X]` whereas the 2010 paper text
   describes diffusion verbally and the 2010 SI uses a different
   notation. Path B docstrings should adopt `−D∇²` for the diffusion
   sign convention regardless of which paper is being cited.
2. **The piecewise threshold structure is genuinely clearer than
   what `13.f90` looks like.** Equations 2–4 each have two cases —
   a 'baseline' case below threshold and a 'modified' case at or
   above. The FORTRAN encodes these as `if (q3d(...) > th) then`
   blocks scattered through `diffusion`/`reaccio`. Path B can keep
   the piecewise shape in code by constructing boolean masks per
   threshold and using `np.where` or array slicing — which both
   reads like the equations and vectorises naturally.

The Hill term `k[S_1]/(1+k_{S_2,S_1}[S_2])` in eq. 1 is the
auto-activation of S_1 under inhibition by S_2. This has the same
algebraic shape as a Hill-1 saturation. Path B can lift the form
directly; the parameter `k_{S_2,S_1}` corresponds to one of the
`BmpFgf`/`FgfBmp` rows in `multicusp.ini`.

**Caveat.** Equation 1's auto-activation term is **not** in the 2010
mammalian model in this form. The 2010 model uses a simpler
threshold-driven activator dynamic. The 2023 paper notes this
explicitly (main article p. 2, right column): the original kernel
'consists of one molecular activator and its inhibitor' and Zimm 'we
changed the original signalling kernel to a three-agent signalling
network'. Path B should **not** adopt eqs. 1–4 verbatim; the
relevant equations are in the 2010 paper's SI.

## Differences from the 2010 model

The SI is admirably explicit about what was changed. From SI p. 2,
'Mathematical modeling' section:

> 'Overall, we chose a conservative approach by keeping most parts of
> the original Toothmaker. We only changed or added parts about which
> we either had direct or indirect evidence for specific differences
> between the two classes' (i.e. mammals vs sharks).

The substantive changes (paraphrasing the SI) are:

1. **Three-component signalling kernel instead of activator–inhibitor
   pair.** Main article p. 2 right col., and SI p. 2 lines 42–55.
   Three signals `S_1` (Fgf-like), `S_2` (Bmp-like), `S_3`
   (Shh-like), each with independent dynamics. The 2010 model uses
   one activator–inhibitor pair plus secondary factor Sec.
   **Path B avoids this entirely.**
2. **Substrate-depletion vs activator-inhibitor topology.** Fig. S3
   illustrates the distinction: Zimm's network can toggle between
   the two regimes via parameters. The 2010 model is purely
   activator-inhibitor. **Path B should not implement the
   substrate-depletion alternative.**
3. **Per-morphogen growth induction `Egr_i`.** SI p. 3 lines 99–104.
   In Zimm's model, epithelial growth is augmented by a sum of
   per-morphogen products (rate × concentration). This is the
   `Egr(i)` row in `multicusp.ini`. The 2010 model has constitutive
   epithelial growth with a single rate. **Path B uses the 2010
   single-rate Egr only.**
4. **Conveyer-belt downgrowth.** SI p. 3 lines 105 onwards.
   Elasmobranch teeth form on a conveyer-belt-like dental lamina, so
   Zimm adds downgrowth pressure (Dgr) and a buoyancy parameter
   (Boy) to model lamina mechanics. `13.f90` already has both
   `downgrowth` and `Bgr` (buoyancy). The 2010 mammalian model has
   downgrowth but not buoyancy in this form. The 2014 tribosphenic
   model is closer to the FORTRAN; Path B should follow the 2014
   formulation, **not** the 2023 conveyer-belt formulation.
5. **Force-multiplier biases instead of GRN-activation biases.**
   SI p. 3 lines 110+ (very explicit):
   > 'Unlike in the original model, where biases described
   > differences in the initial activation of the GRN, we use biases
   > to describe force differences that locally affect how teeth
   > grow in silico.'

   Equation 5 (SI p. 4) gives the force-multiplier form:
   `F_j(i) = F_j(i) · b_j` where `j ∈ {x, y, z}`. This is not the
   2010 model. **Path B should follow the 2010 GRN-activation bias
   semantics**, which is what `13.f90`'s `apply_border_bias` does
   (it modulates Fgf/Bmp activation in border cells, not force
   vectors).
6. **Auto-activation Hill term in eq. 1.** Already noted above; not
   in the 2010 mammalian dynamics.
7. **Conditional inhibition** (k_{S_2,S_3+} vs k_{S_2,S_3-}, see
   `multicusp.ini` 'Interactions Shh' row of four numbers and Table
   S2 entries). Two coupling rates that switch sign based on a
   threshold (Zhang et al. 2000, cited at SI ref. 24). Path B should
   not implement this — the 2010 model uses simple unsigned
   coupling.

## What the SI clarifies about the 2010 model

A few points where the 2023 SI re-states 2010 mechanics in a way
clearer than the 2010 paper itself:

1. **Differentiation as a continuous variable.** SI p. 3 lines 91–93:
   > 'In the model, this is implemented by the variable D(i),
   > where D denotes the differentiation stage of the cell i,
   > increasingly slowing down mitotic activity and, eventually,
   > leading to its arrest as differentiation reaches a value of 1.'

   This matches `13.f90`'s `marge` (cell-margin / differentiation
   state) variable. The framing as 'D(i) ∈ [0, 1] saturating' is
   cleaner than what the 2010 paper provides and worth using as a
   docstring.
2. **Threshold-as-developmental-stage interpretation.** SI p. 3
   lines 85–89:
   > 'These thresholds account for the observation that some
   > processes (such as the induction of tissue differentiation) are
   > only seen ensuing sustained morphogen expression, effectively
   > emulating developmental delays between different developmental
   > stages of odontogenesis.'

   This is a clearer biological gloss on the threshold
   parameterisation than the 2010 paper offers. Worth recording in
   Path B's threshold-parameter docstring.
3. **Bias radius semantics (rbia).** SI p. 4 lines 1–4:
   > 'Within the radius (in x,y direction) at which the effects of
   > biases apply, we multiply the sum of all force vector
   > components ⃗F_j resulting from morphogenetic processes in cell i
   > with their respective bias parameters b_j.'

   Even though the 2023 force-multiplier semantics differ from the
   2010 GRN semantics, the **radius** concept (`rbia` in
   `multicusp.ini`, `radibias` in `13.f90`) is consistent across
   both: a circular zone around the tooth centre within which biases
   apply. This is a clearer statement of what `radibias` does than
   the 2010 paper provides. The role of biases ('what is multiplied')
   differs between the two papers; the radius semantics do not.

## Validation methodology

Two strands. Both are interesting but neither transfers to Path B
unchanged.

1. **Morphological-distance metrics (SI 'Morphological distances'
   p. 5).** Two metrics: discrete cosine Fourier-coefficient
   Euclidean distance (eq. 6) on apical tooth outlines, and an
   average-outline-distance metric (eq. 7) inspired by Procrustes
   methods. Both quantify shape differences between in-silico tooth
   outlines and target outlines. Useful intellectually, but Path B's
   validation is against FORTRAN goldens (cell count, vertex
   envelope) per the charter, not against real-tooth outlines. If
   Path B ever needs to compare against real teeth, the
   Fourier-coefficient approach is well-described and citation-ready.
2. **In vivo bead-implantation experiments (SI 'Bead experiments'
   p. 7).** Affi-Gel beads soaked in Fgf3 or SU5402 (FgfR
   inhibitor), implanted next to developing tooth buds in
   *S. canicula* embryos, then microCT-scanned to compare against
   `a_1`-perturbed in silico variants. This is genuinely impressive
   experimental validation but is shark-specific and not relevant
   to Path B's mammalian-FORTRAN-faithful target.

The validation logic — perturb a parameter in silico, perturb the
homologous biological process in vivo, compare phenotype — is a
template for any future biological validation of the Path B
implementation, but is out of scope for the initial work.

## Other SI features

- **Phenotypic-region degeneracy analysis (SI eqs. 8–12)**:
  parameter-enrichment, Hamming-distance, anisotropy. All on top of
  large variant ensembles (5000–30000 simulations per round, SI
  p. 4). Not relevant to Path B but worth noting that the
  `multicusp.ini` is one *parameterisation*, and the model produces
  a high-dimensional phenotypic landscape under perturbation.
- **R packages (SI p. 8 line 274)**: PCAtools, ICSPN (HotellingT2),
  RVAideMemoire, plotting via gnuplot 5.2. Irrelevant for Path B
  but worth noting if anyone ever rebuilds the analysis pipeline.

## Path-B implications

Honest answer: **mostly nothing**. The 2023 paper is more
shark-divergent from the 2010 mammalian model than the SI's
'conservative approach' framing suggests. Five substantive changes
(three-component kernel, substrate-depletion topology,
per-morphogen Egr_i, conveyer-belt downgrowth, force-multiplier
biases) are exactly the bits Path B must avoid.

Things to **adopt** from the 2023 SI:

1. **Type-column grouping in any internal parameter glossary.** The
   four-class partition (diffusion / GRN / mechanical / bias / growth
   / threshold) is sound and matches the natural module structure
   of `silicoshark/`. SI Table S1, p. 9.
2. **Stable parameter symbols where the 2010 and 2023 names agree:
   `Egr`, `Mgr`, `Dgr`, `Boy`, `Adh`, `Rep`, `Ntr`, `μ`, `D_i`,
   `k_{ij}`, `a_i`, `θ` (or `ϑ`) for thresholds.** SI Table S1,
   p. 9. These are the docstring symbols Path B should cite when a
   2010 equation is implemented.
3. **The `−D∇²` diffusion-term sign convention** in equation
   docstrings. SI eqs. 1–4, p. 3.
4. **Piecewise-threshold equation shape with boolean masks in
   NumPy.** SI eqs. 2–4 are written as `\begin{cases}` per
   threshold; the natural NumPy idiom is `np.where(mask, branch_a,
   branch_b)` or boolean fancy-indexing. This reads like the
   equations and vectorises across cells.
5. **The differentiation-as-saturating-`D(i)` framing** for the
   `marge` docstring. SI p. 3 lines 91–93.
6. **The bias-radius `rbia` semantics docstring** (radius within
   which biases apply, regardless of what they multiply). SI p. 4
   lines 1–4.

Things to **avoid** from the 2023 paper:

1. The `S_1`/`S_2`/`S_3` placeholder names. Path B uses Fgf, Bmp,
   Shh (or abstract `act`, `inh`, `sec`) per the 2010 paper. SI
   Table S1 footnote, p. 9.
2. Equations 1–4 of the 2023 SI as model equations. They describe
   the three-component shark kernel, not the 2010 mammalian
   activator–inhibitor pair. SI p. 3.
3. Equation 5 of the 2023 SI (force-multiplier biases). The 2010
   model uses GRN-activation biases. SI p. 4.
4. Table S2's parameter values. They are shark-fitted. SI p. 10.
5. The `multicusp.ini` file as a parameter template — it is a
   2010-baseline-plus-shark-extensions hybrid. Path B should derive
   defaults from `examples/seal.txt` (which the FORTRAN-faithful
   translation already validates against) and from the 2010 and
   2014 paper SI tables.
6. The `Egr(i)` per-morphogen growth-induction vector. Not in the
   2010 model. SI p. 3 lines 99–104.
7. The conditional sign-flipping coupling
   (`k(S_2,S_3)+`, `k(S_2,S_3)-`) implied by the 'Interactions
   Shh' row in `multicusp.ini` and Table S2. Not in the 2010 model.

## Summary

The 2023 SI is competent, conservative-sounding, but in fact
substantially redefines the model's biology: morphogen kernel,
biases, growth induction, downgrowth dynamics. For Path B's
purpose — a clean re-implementation of the 2010 mammalian model
plus 2014 tribosphenic extension — the 2023 paper is reference
material for **symbol conventions** (Egr, Mgr, μ, θ, D_i, k_{ij},
a_i) and **typesetting style** (piecewise-threshold equations,
`−D∇²` diffusion sign), and a useful negative example for
everything else. The 2010 paper and its SI remain the canonical
source.
