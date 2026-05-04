---
title: 'Paper review — Salazar-Ciudad and Jernvall, *Nature* 464 (2010), and SI'
author: Lyndon Drake (with Claude Code)
date: 2026-05-04
---

## Purpose

This review extracts, from the 2010 paper and its SI, the material Path B
needs to re-implement the model: equations, parameter table, geometric
primitives, discretisation choices, validation evidence, and — most
useful for planning — the design decisions the paper leaves to the
implementer. The biology narrative (Figs 1–4 of the main article) is not
reviewed.

Citations follow the form 'main p. X' for the *Nature* letter and Methods,
and 'SI p. X' for the supplementary PDF (`41586_2010_BFnature08838_MOESM283_ESM.pdf`).
Equation numbers follow the paper.

## Model equations

The model has three coupled parts: a reaction–diffusion gene network
(eqs 14–18), mechanical forces on cell centres (eqs 1–9, 13), and a cell
division and differentiation rule (no equation; described in prose). All
parameters are dimensionless (main p. 586, last paragraph of 'Initial
conditions').

### Mechanical forces

Forces are summed per cell and integrated by forward Euler with a single
global step (main p. 583, fig. 1 caption: 'all equations are integrated
using the Euler method').

**Repulsion (eq. 1, main p. 585):**

$$
\mathbf{f}_{ij} = k_{\text{Rep}}\,\bigl(\lVert\mathbf{p}_j - \mathbf{p}_i\rVert - \lVert\mathbf{p}_0\rVert\bigr)\,(\mathbf{p}_j - \mathbf{p}_i)
$$

— $k_\text{Rep}$ is a stiffness, $\lVert\mathbf{p}_0\rVert$ is the
*initial* distance between cells $i$ and $j$ (Methods explicitly notes:
'$\lVert\mathbf{p}_0\rVert$ is *not* a parameter of the model … When the
distance between two cell centres is larger than $\lVert\mathbf{p}_0\rVert$
there is a constant attraction force resulting from the adhesion of cells').
Maps to FORTRAN `Rep`.

**Adhesion (eq. 2, main p. 585):**

$$
\mathbf{f}_{ij} = k_{\text{Adh}}\,\frac{\mathbf{p}_j - \mathbf{p}_i}{\lVert\mathbf{p}_j - \mathbf{p}_i\rVert}
$$

— a unit-vector pull. Active only when $\lVert\mathbf{p}_j - \mathbf{p}_i\rVert > \lVert\mathbf{p}_0\rVert$. Maps to FORTRAN `Adh`.

**Net positional update from (1)–(2) (eq. 3, main p. 585):**

$$
\frac{\partial \mathbf{p}_i}{\partial t} = \sum_{j \in \mathcal{N}(i)} \mathbf{f}_{ij}
$$

— neighbour set $\mathcal{N}(i)$ (note the paper does *not* define which
neighbours; see §What the paper does NOT specify).

**Nucleus traction (eq. 4, main p. 585):**

$$
\frac{\partial \mathbf{p}_i}{\partial t} = k_{\text{Ntr}}\,(1 - d_i)\,\Bigl(-\mathbf{p}_i + \tfrac{1}{n_i}\sum_{j \in \mathcal{N}(i)} \mathbf{p}_j\Bigr)
$$

— $n_i$ = neighbour count; $d_i$ = differentiation state (clamped to
$[0,1]$, see eq. 6); the $(1 - d_i)$ factor freezes differentiated
(knot) cells. Maps to FORTRAN `Ntr`.

**Epithelial growth (eq. 5, main p. 585):**

$$
\frac{\partial \mathbf{p}_i}{\partial t} = k_{\text{Egr}}\,(1 - d_i)\,\frac{\sum_j \hat{\mathbf{u}}_{ij}}{\bigl\lVert \sum_j \hat{\mathbf{u}}_{ij} \bigr\rVert}
$$

— each cell shoves its neighbours away along the average outward unit
vector. Magnitude $k_\text{Egr}(1 - d_i)$. Maps to FORTRAN `Egr`. Crucial:
this is the term modulated by Diff (the differentiation field) via
$1 - d_i$, *not* a separate growth rule.

**Differentiation evolution (eq. 6, main p. 585):**

$$
\frac{\partial d_i}{\partial t} = k_{\text{Dff}}\,[\text{Sec}]
$$

— Sec drives differentiation linearly. $d_i$ ranges $[0, 1]$ (stated in
prose). Maps to FORTRAN `Dff`.

**Cervical-loop downgrowth (eqs. 7–9, main p. 585), border cells only:**

$$
\frac{\partial x_{\text{border},i}}{\partial t} = \frac{d_m\,k_{\text{Egr}}\,(1 - d_i)}{\sqrt{d_x^2 + d_y^2 + k_{\text{Bgr}}^2}},\qquad
\frac{\partial y_{\text{border},i}}{\partial t} = \frac{d_m\,k_{\text{Egr}}\,(1 - d_i)}{\sqrt{d_x^2 + d_y^2 + k_{\text{Bgr}}^2}}
$$

$$
\frac{\partial z_{\text{border},i}}{\partial t} = \frac{-d_m\,k_{\text{Egr}}\,(1 - d_i)}{\sqrt{d_x^2 + d_y^2 + k_{\text{Bgr}}^2}}
$$

— note the *minus* sign on $z$: border cells are pushed downwards.
$d_m$ is the local mesenchymal-tissue depth (defined in eqs 10–12).
$k_{\text{Bgr}}$ ('border growth') controls the tendency of border cells
to grow downward when there is little underlying mesenchyme. Maps to
FORTRAN `Bgr`.

**Mesenchyme depth (eqs. 10–12, main p. 586):**

$$
d_x = d_{xo}\,c, \qquad d_y = d_{yo}\,c, \qquad c = 1 + \frac{k_{\text{Mgr}}\,[\text{Sec}]}{\sqrt{d_{xo}^2 + d_{yo}^2}}
$$

— $d_{xo}, d_{yo}$ are the $x, y$ components of eq. (5)'s growth vector
*before* application of the Sec-driven scaling. $k_{\text{Mgr}}$
('mesenchymal growth') = paper symbol for FORTRAN `Mgr`. Sec causes
mesenchyme to thicken under signalling regions, which then resists
invagination via the buoyancy term.

**Mesenchymal buoyancy (eq. 13, main p. 586):**

$$
\frac{\partial \mathbf{p}_i}{\partial t} = k_{\text{Boy}}\,(1 - d_i)\,[\text{Sec}]\,\hat{\mathbf{n}}
$$

— $\hat{\mathbf{n}}$ is the outward unit normal at cell $i$ (apical
direction); buoyancy is *only* applied to cells whose averaged
neighbour-vector-product points apically (Methods note). Maps to FORTRAN
`Boy`.

### Reaction–diffusion gene network

Three molecular species are modelled in the epithelium: Act
(activator), Inh (inhibitor), Sec (secondary signal). Mesenchyme cells
hold only Sec and Inh. All three species diffuse via the discretised
Fick's law (described in prose, main p. 585; SI fig. 1 shows the
mechanical grid this reuses).

**Activator (eq. 14, main p. 586):**

$$
\frac{\partial[\text{Act}]}{\partial t} = \frac{k_{\text{Act}}\,[\text{Act}]}{1 + k_{\text{Inh}}\,[\text{Act}]} - k_{\text{Deg}}\,[\text{Act}] + k_{\text{Da}}\,\nabla^2[\text{Act}]
$$

— the only term with non-trivial gene-network behaviour. The
$k_{\text{Inh}}\,[\text{Act}]$ in the denominator is, *per the paper*,
the inhibition of Act by Inh; but the paper writes $[\text{Act}]$ where
one expects $[\text{Inh}]$. Reading the paper literally, the activator
self-inhibits via a saturating Hill-style term in its own concentration.
Cross-checking against the FORTRAN: humppa's `reaccio_difusio` (and
hence `13.f90`'s `apply_diffusion`) divides by $1 + \text{Inh} \cdot [\text{Inh}]$,
i.e. the denominator depends on Inh, not Act. The paper's $[\text{Act}]$
in the denominator is therefore a **typo**; the FORTRAN is the
authoritative form. This is the single most important paper-vs-FORTRAN
discrepancy and Path B should follow the FORTRAN. Maps: $k_{\text{Act}}\!=\!$ FORTRAN
`Act`; $k_{\text{Inh}}\!=\!$ FORTRAN `Inh`; $k_{\text{Deg}}\!=\!$ FORTRAN
`Deg`; $k_{\text{Da}}\!=\!$ FORTRAN `Da`.

**Inhibitor (eq. 15, main p. 586) — applies only to enamel-knot cells (see
eq. 17):**

$$
\frac{\partial[\text{Inh}]}{\partial t} = -k_{\text{Deg}}\,[\text{Inh}] + k_{\text{Di}}\,\nabla^2[\text{Inh}]
$$

— Inh has only degradation and diffusion in non-knot cells. $k_{\text{Di}}\!=\!$
FORTRAN `Di`.

**Secondary signal (eq. 16, main p. 586):**

$$
\frac{\partial[\text{Sec}]}{\partial t} = -k_{\text{Deg}}\,[\text{Sec}] + k_{\text{Ds}}\,\nabla^2[\text{Sec}]
$$

— same form as Inh. $k_{\text{Ds}}\!=\!$ FORTRAN `Ds`.

**Knot-secretion overrides (eqs. 17–18, main p. 586):** a cell becomes
an enamel knot once $[\text{Act}] \ge 1$ (an arbitrary threshold,
Methods 'Differentiation is irreversible' paragraph). Knot cells
override eq. (15) with

$$
\frac{\partial[\text{Inh}]}{\partial t} = [\text{Act}] - k_{\text{Deg}}\,[\text{Inh}] + k_{\text{Di}}\,\nabla^2[\text{Inh}].
$$

Cells whose differentiation $d_i \ge k_{\text{Int}}$ ('initial threshold')
also begin secreting Inh; cells whose $d_i \ge k_{\text{Set}}$
('secondary threshold') begin secreting Sec under

$$
\frac{\partial[\text{Sec}]}{\partial t} = k_{\text{Sec}} - k_{\text{Deg}}\,[\text{Sec}] + k_{\text{Ds}}\,\nabla^2[\text{Sec}].
$$

Mesenchymal cells: same equations as epithelium for Inh and Sec; *no*
Act reaction terms (only diffusion). Stated in prose, main p. 586:
'mesenchymal cells have no reaction terms, and thus only epithelial cells
secrete the three growth factors'.

### Cell division

No equation. The rule (Methods, main p. 585, paragraph following eq. 6):
'Cells divide after the connection between two cell centres is equal to
or exceeds two units of space (initial distance between cells equals one
unit) and is implemented by placing a new cell at the midpoint between
the two mother cells.' SI fig. 1 caption clarifies that a new cell adds
both a centre and connections to the original triangular mesh. The
concentrations of all species in the new cell are 'the average of the
two mother cells'. Knots: 'If one of the mother cells is an enamel knot
cell, the new cell is not a knot cell.'

### Border (cervical-loop) bias

No standalone equation. Anterior border = the initial cells anterior to
the centre with $y > 0$, plus all of their descendants by division.
Posterior border = same with $y < 0$ (main p. 585, 'Growth biases'
paragraph). The anterior bias parameter $k_{\text{Abi}}$ (FORTRAN `Abi`)
multiplies the eq. (5) growth magnitude on anterior border cells;
$k_{\text{Pbi}}$ (FORTRAN `Pbi`) does the same for posterior border. The
paper notes 'Because of the symmetric nature of these biases, we only
modulated $k_{\text{Pbi}}$ in the variation analysis' — i.e. for the
seal runs, anterior and posterior biases are kept equal in spirit but
$k_{\text{Pbi}}$ is the one tuned. There is no equation in the paper for
left/right (lingual/buccal) bias, but FORTRAN has `Lbi`/`Rbi` parameters,
which SI Table 1 confirms are kept at 1.0 for the seal teeth.

## Parameter table

The full parameter inventory comes from SI Table 1 (SI p. 6) for the
seal runs. Units are dimensionless throughout (main p. 586). I have
added FORTRAN/`silicoshark/multicusp.ini` mappings against the
divergence catalogue (`docs/research/13f90-vs-humppa-divergences.md`,
note 7).

| Paper symbol | Role (paper) | Ringed-seal value | Grey-seal value | FORTRAN name | Notes |
|---|---|---|---|---|---|
| Di | Inh diffusion rate ($k_{\text{Di}}$) | 0.2 | 0.2 | `Di` | eq. 15 |
| Act | Act autoregulation strength ($k_{\text{Act}}$) | 1.1 | 0.29 | `Act` | eq. 14 |
| Inh | Inh strength on Act ($k_{\text{Inh}}$) | 26 | 26 | `Inh` | eq. 14 (paper writes $[\text{Act}]$ in denominator; FORTRAN uses $[\text{Inh}]$) |
| Set | Secondary differentiation threshold ($k_{\text{Set}}$) | 0.95 | 0.95 | `Set` | eq. 18 onset |
| Da | Act diffusion rate ($k_{\text{Da}}$) | 0.2 | 0.2 | `Da` | eq. 14 |
| Bgr | Border growth ($k_{\text{Bgr}}$) | 1 | 1 | `Bgr` | eqs. 7–9 |
| Egr | Epithelial growth ($k_{\text{Egr}}$) | 0.017 | 0.0146 | `Egr` | eq. 5 |
| Abi | Anterior bias ($k_{\text{Abi}}$) | 15 | 15 | `Abi` | not in paper as equation |
| Pbi | Posterior bias ($k_{\text{Pbi}}$) | 18 | 16 | `Pbi` | not in paper as equation |
| Adh | Adhesion ($k_{\text{Adh}}$) | 0.001 | 0.001 | `Adh` | eq. 2 |
| Ntr | Nuclear traction ($k_{\text{Ntr}}$) | 0.00001 | 0.00001 | `Ntr` | eq. 4 |
| Mgr | Mesenchymal growth ($k_{\text{Mgr}}$) | 200 | 200 | `Mgr` | eq. 12 |
| Deg | Generic degradation rate ($k_{\text{Deg}}$) | 0.1 | 0.1 | `Deg` | eqs. 14, 15, 16 |
| Ds | Sec diffusion rate ($k_{\text{Ds}}$) | 0.2 | 0.2 | `Ds` | eq. 16 |
| Dff | Differentiation efficiency ($k_{\text{Dff}}$) | 0.0002 | 0.0002 | `Dff` | eq. 6 |
| Int | Initial differentiation threshold ($k_{\text{Int}}$) | 0.19 | 0.19 | `Int` | onset of Inh secretion |
| Rep | Repulsion stiffness ($k_{\text{Rep}}$) | 1 | 1 | `Rep` | eq. 1 |
| Dgr | Cervical-loop downgrowth | 10000 | 10000 | `Dgr` | not directly in main equations; Methods 'Cervical loop' paragraph |
| Boy | Buoyancy ($k_{\text{Boy}}$) | 0.1 | 0.1 | `Boy` | eq. 13 |
| Sec | Sec secretion rate ($k_{\text{Sec}}$) | 0.03 | 0.03 | `Sec` | eq. 18 |

P2–P5 ringed and grey seal runs (the tooth-row results) use additional
italicised values for `Egr`, `Act`, `Pbi` per SI Table 1.

The 'Seal 1–9' alternative-parameter runs (SI Table 5, p. 9–10) use a
much wider range — `Mgr` 2000–14000, `Dgr` 10000–20000, `Egr`
0.018–0.03, `Inh` 10–50 — and additionally tune a `Time`
(iteration-count) parameter that varies from 5000–10000 between virtual
teeth. This is the closest the paper comes to publishing the iteration
budget; the seal 'wild-type' (Fig. 1c) uses 10000 iterations explicitly
labelled in that figure.

**Parameters present in FORTRAN but not in this paper.** `Bwi` (border
band width, humppa `Bwi`/`bwi`) and `Swi` (humppa `tadi`) appear in the
2014 paper / FORTRAN parameter file but not in 2010. The 2010 SI never
mentions a 'border band' — the cervical-loop logic is described as
applying to all border cells without any band thickness. This is a 2014
addition. (See note 1 of the divergences catalogue: `13.f90`'s
`apply_border_bias` in fact ignores `Swi` because `setparams` skips
`parap(6)`. Path B can omit `Swi` from the 2010-only path.)

**One unmapped paper parameter.** Symbol $\lVert\mathbf{p}_0\rVert$ in
eqs. (1)–(2): the *initial* equilibrium distance, fixed to one unit by
the initial-conditions paragraph ('initial distance between cells
equals one unit'). Methods explicitly state it 'is not a parameter of
the model'. Path B should treat it as the constant 1.0.

## Geometric primitives

**Cell layout (main p. 583, fig. 1 caption; main p. 585 'Cell repulsion';
SI fig. 1):**

- A cell is a *centre point* with neighbour links — it is *not* a sphere
  or a disc node. It is described as enclosing 'all the points in space
  that are closer to its centre than to the centre of any other cell',
  i.e. a Voronoi region implicitly defined by the centre set.
- A second 'triangular mesh' is built whose corners are the centres of
  three mutually-adjacent cells — these corners are 'Voronoi nodes',
  equidistant from the three cells. SI fig. 1c–f shows this triangular
  mesh.
- Each epithelial cell is 'a three-dimensional volume including the
  cell itself and its immediate extracellular surrounding space'. The
  basal corners (mesenchymal side) are placed directly below the apical
  corners by a fixed z-unit — main p. 586 'each epithelial cell has six
  epithelial neighbours, except for cells on the borders that have three
  or four depending on the location'. (The exact six-or-three-or-four
  rule is specific to the *initial* hexagonal lattice; the paper does
  not describe how neighbour counts evolve under cell division.)

**Initial condition (main p. 586 'Initial conditions'; SI fig. 1a, c):**

- Seven epithelial cells in a hexagon (six surrounding one central
  cell), all at distance $\lVert\mathbf{p}_0\rVert = 1$.
- Two layers of mesenchymal cells underneath. Each epithelial cell has
  one underlying mesenchymal neighbour.
- An 'arbitrary constant source of Act' exists in the border epithelial
  cells initially (main p. 586). All other species concentrations and
  $d_i$ values are zero in all cells initially.

**Border cells.** Defined twice in the paper:

1. *Topologically* (Methods, main p. 585, 'Growth biases'): 'The
   anterior border is defined as the border cells that in the initial
   condition are anterior to the centre (positive $y$ values), and all
   of their descendants that lay in the tenth border. The posterior
   border is defined in the same way for negative values of $y$.'
2. *Geometrically* (main p. 583, fig. 1 caption): cells with fewer
   epithelial neighbours than the interior count.

The descendant rule is the load-bearing one for biases; the geometric
definition is what the topology walk in `13.f90` actually computes.

**Knots ('signalling centres').** A cell becomes an enamel knot when
$[\text{Act}] \ge 1$ in that cell (main p. 586, 'Differentiation is
irreversible'). The threshold value 1.0 is hard-coded. Knot cells (i)
override eq. (15) with eq. (17) Inh-secretion, (ii) freeze growth
(via $1 - d_i$ saturating to 0), (iii) do not produce knot daughters
when they divide ('the new cell is not a knot cell'). The paper does
*not* say what value of $d_i$ a knot cell has; the FORTRAN forces it to
1.0.

## Discretisation choices

- **Time integration:** forward Euler ('Euler method', main p. 583
  fig. 1 caption; main p. 585 'All equations are integrated using the
  Euler numerical method'). One scalar time-step, applied to all
  per-cell variables in lockstep. Iteration count is part of the
  'Time' parameter (SI Table 5).
- **Spatial discretisation of the diffusion Laplacian:** finite-volume
  on the triangular mesh ('finite volume method', main p. 585): 'making
  the molecular flux between two cells proportional to the area in
  which they are in contact or near each other'. The flux through a
  Voronoi face replaces $\nabla^2 c$. The paper does *not* give the
  formula for that area, only the verbal rationale; it is left to the
  implementer. The FORTRAN computes it from the per-neighbour `marge`
  array (the 'cell margins') seeded by `calculate_margins`.
- **Per-cell update order:** the paper writes each force as a separate
  $\partial \mathbf{p}_i / \partial t$, with the implicit understanding
  that they sum into a single Euler increment per iteration. The
  paper does not specify the order in which the forces are evaluated,
  whether all are computed against the same start-of-iteration positions
  ('Jacobi') or whether earlier forces update positions read by later
  forces ('Gauss–Seidel'). The FORTRAN is Jacobi for the gene network
  and Gauss–Seidel for the mechanical forces — this is one of the
  strongest implementer-choice points (see §What the paper does NOT
  specify).
- **Stability constraints:** the paper does not state a CFL-like
  bound. Empirically, the seal parameters give a stable run at 10000
  iterations with $k_{\text{Egr}} = 0.017$ and a unit cell spacing,
  which sets a baseline for what 'stable' means.
- **Knot freeze:** $d_i$ ranges in $[0, 1]$ (stated in prose). The
  $(1 - d_i)$ multiplier in eqs. (4), (5), (7)–(9), (13) is the *only*
  mechanism by which differentiated cells stop moving. Once $d_i = 1$,
  the only force a knot cell can experience is repulsion (eq. 1) and
  adhesion (eq. 2), neither of which carries a $(1 - d_i)$ factor.
- **Cell division side effects:** when a new cell is added at the
  midpoint, the triangular mesh is updated (SI fig. 1 reference). The
  paper does not give the algorithm for this update; the FORTRAN's
  `add_cell` walks the topology and re-stitches neighbour lists.

## Validation in the paper

The paper validates against four reference points:

1. **Real seal teeth ($n = 70$, *Phoca hispida ladogensis*).** Cusp tip
   $(x, y)$ coordinates measured from photographs of P4s (SI Table 2,
   pp. 6–8) and the top-cusp angle. Geometric morphometrics
   (Procrustes) compare four-cusp seal teeth to virtual teeth with
   four cusps.
2. **Heterodont grey seal teeth (*Halichoerus grypus*).** Tooth-row
   variation produced *in silico* by changing two parameters
   (`Egr`, `Act`) is compared visually to real grey-seal P2–P5
   morphology (Fig. 4).
3. **Cusp-width pattern.** The cusp-c-narrows-when-cusp-d-present
   pattern (SI Fig. 5) reproduces in real teeth ($P = 0.009$
   Mann–Whitney) and in *in silico* teeth.
4. **Population-level PCA of cusp shape variation.** The first two PCs
   of real seal cusp positions account for 46% and 30% of variance;
   *in silico* runs perturbing each parameter ±2% are mapped onto the
   same PC axes (Fig. 2c, Fig. 3).

The validation outputs are: cusp count, top-cusp angle, cusp-tip
$(x, y)$, and PCA loadings on cusp position. **Vertex-level
mesh equality is never used.** Path B's golden-tolerance regime
(±5% cell count, ±10% vertex envelope, no NaN/Inf) sits in the same
spirit as the paper's own evaluation: cusp count is a topological
invariant rather than a per-vertex match. This makes the existing
golden tolerance scientifically defensible.

## What the paper does NOT specify

These are the load-bearing implementer choices that Path B will have to
make. Listed in roughly decreasing order of likely numerical impact.

1. **The exact form of the diffusion Laplacian.** Paper says 'finite
   volume on the triangular mesh, with flux ∝ contact area'. It does
   *not* give the formula. The FORTRAN computes a per-edge weight via
   `calculate_margins`. Path B must choose between (a) reproducing
   FORTRAN-style margins, (b) the standard cotangent-weight Laplacian
   (which is the textbook discrete Laplace–Beltrami operator on a
   triangular mesh), or (c) an equal-weight neighbour average. These
   are not equivalent; cotangent weights match curvature behaviour
   most cleanly for curved-surface triangulations.
2. **Update order: Jacobi vs Gauss–Seidel.** The paper's prose strongly
   suggests Jacobi (each $\partial / \partial t$ is a function of the
   *current* state). The FORTRAN is mixed: gene network is Jacobi (uses
   a copy `hq3d`), but mechanical forces in `iteration` are
   Gauss–Seidel — repulsion runs first and updates positions in place,
   then adhesion, then nucleus traction, etc. Path B should default to
   pure Jacobi and treat any divergence from the FORTRAN as
   FORTRAN-implementation-specific.
3. **Border-cell identification each iteration.** Paper says border
   cells are those with fewer-than-interior neighbours. With cell
   division, this needs recomputing each iteration. The FORTRAN has the
   topology walk in `update_border_cells` / `update_cell_position`.
   Path B can either (a) compute border-status from neighbour-count
   thresholds each step, vectorised (b) maintain it incrementally on
   division. Option (a) is cleanest and avoids the panic-and-return
   silent bug.
4. **Neighbour-set update on cell division.** Paper says SI fig. 1 and
   no more. The FORTRAN uses a topology walk that the divergence
   catalogue documents (see `13f90-vs-humppa-divergences.md` and
   recent fix `coreop2d.py:add_cell` — the `iiii` preservation bug).
   This is the largest source of accidental implementation choice in
   the FORTRAN. A vectorised Path B re-triangulation (e.g.
   `scipy.spatial.Delaunay` over the projected positions) could
   replace the walk entirely.
5. **The 'tenth border' definition for biases.** Methods say bias
   applies to 'all of their descendants that lay in the tenth border'.
   The 'tenth border' is undefined elsewhere in the paper. The
   FORTRAN's `apply_border_bias` reads `parap(?)` for the bias
   half-plane condition (`y > 0` or `y < 0`) and applies the multiplier
   to all border cells in that half-plane. Path B should default to
   the FORTRAN reading (any border cell, not just descendants of the
   initial-condition border cells).
6. **The neighbour set in eq. (3) repulsion sum.** Paper says
   '$\sum_{j \in \mathcal{N}(i)}$' but does not define $\mathcal{N}(i)$
   for repulsion — direct mesh neighbours, all cells within some cutoff,
   or all cells full-stop? FORTRAN distinguishes `pushing`
   (mesh-neighbour repulsion, eq. 1) from `pushingnovei` (non-neighbour
   repulsion against everyone too close). Both are present. The paper
   describes only the first; the FORTRAN's `pushingnovei` is therefore
   an implementer addition, but plausibly necessary to avoid mesh
   tangling.
7. **The 'apical normal' $\hat{\mathbf{n}}$ in eq. (13).** Paper says
   'a unit vector normal to the epithelial cell surface', and
   parenthetically 'this assumes that the mesenchyme behaves like a
   compressed fluid'. It does *not* give the discrete formula. FORTRAN
   computes per-cell normals from neighbour cross products
   (`promig`/`apply_nucleus_traction` and similar code paths).
   Standard discrete options are area-weighted face-normal averages or
   angle-weighted normals.
8. **Buoyancy gating: 'apically pointing'.** Methods say buoyancy is
   applied 'only the resulting vectors pointing apically are taken for
   each product'. This is opaque. The FORTRAN gates by sign of the
   z-component of the averaged neighbour cross-product. Path B should
   document its choice explicitly.
9. **Cell-division mesh repair.** Paper says 'a new cell adds a new
   cell centre and connections to the original triangular mesh
   (Supplementary Fig. 1)'. SI fig. 1 is a *picture*, not an algorithm.
   Path B has full freedom here — local edge-flip / Delaunay restore /
   global re-triangulation are all defensible.
10. **Knot daughter-cell concentrations.** Paper says daughters get the
    average of the two mother cells' species concentrations, but does
    *not* say what $d_i$ they inherit. The FORTRAN sets it to the
    average of the mothers' $d_i$ and *not* zero, which means a
    daughter of two knot cells (which the paper forbids being a knot)
    can still inherit $d_i$ near 1, freezing it.
11. **Repulsion 'rest length' on division.** Eq. (1) uses the
    *initial* spacing as rest length. After division, what is the rest
    length of the new edges? Paper is silent. FORTRAN treats it as the
    same constant `1.0` for all edges, including new ones.
12. **The Inh-by-Act denominator.** As noted under eq. (14): paper
    writes $1 + k_{\text{Inh}} [\text{Act}]$, FORTRAN computes
    $1 + \text{Inh} \cdot [\text{Inh}]$. The FORTRAN form is the
    biologically meaningful one (Inh inhibits Act, as the paper's
    prose says everywhere else); Path B should treat the paper as a
    typo and follow the FORTRAN.

## Discrepancies between paper and `13.f90`

1. **Eq. (14) denominator typo.** Paper $[\text{Act}]$ vs FORTRAN
   $[\text{Inh}]$. Already covered.
2. **`Swi` parameter (humppa `tadi`) absent from paper.** Paper does
   not mention a 'border band' anywhere; the parameter exists in
   humppa and was carried into `13.f90` then silently dropped. No
   conflict — `Swi = 0` is consistent with the paper's silence.
3. **`pushingnovei` (non-neighbour repulsion).** Not described. The
   FORTRAN runs an additional repulsion loop against all cells whose
   distance is below threshold, regardless of mesh adjacency. This is
   plausibly load-bearing for stability.
4. **Cell-division topology walk.** Paper offers SI fig. 1 only; the
   FORTRAN's walk has at least three documented `13.f90`-only quirks
   (the `iiii` preservation, the silent panic-and-return, the
   txungu else-branch) — none of which are in the paper. Path B
   should not reproduce any of them.
5. **FMA-divide-by-zero guard in `apply_diffusion`.** This is purely
   a FORTRAN-binary numerical artefact. Eq. (14)–(16) are silent on
   degenerate cells.
6. **Label-77 fall-through in `calculate_margins`.** Pure FORTRAN
   accident. Margins should come from a mesh, not from previous
   slot's value.
7. **Border-band threshold (`Bwi` in `apply_border_bias`).** 2014
   addition; not in 2010. The 2010 paper applies the bias to the
   whole border half-plane unconditionally.
8. **`Dgr` (cervical-loop downgrowth, value 10000).** Paper Methods
   describe cervical-loop downgrowth in prose ('Epithelial cells at
   the border of the tooth (cervical loop) grow to engulf the
   underlying dental mesenchyme'), but the parameter `Dgr` is in
   SI Table 1 without a corresponding equation. The FORTRAN treats
   `Dgr` as a multiplier on the eq. (7)–(9) downward force when the
   underlying mesenchyme is thin. This is an implementer-supplied
   formulation.

## Path B implications

Numbered design decisions Path B will face. 'P' = paper specifies; 'I' =
implementer's choice (paper silent or under-specified).

1. **(P)** The full set of 20 parameters in SI Table 1 with the
   ringed-seal values is the seal `examples/seal.txt` ground truth.
   Path B's parameter loader should read these by paper symbol;
   `Bwi`/`Swi` are 2014 additions that can be initialised to defaults
   (0 and 0 respectively) for 2010-faithful runs.
2. **(P)** Reaction–diffusion for Act, Inh, Sec via eqs. (14)–(18),
   with the FORTRAN denominator $1 + \text{Inh} \cdot [\text{Inh}]$
   in eq. (14) (paper typo).
3. **(P)** Forward Euler integration with a fixed iteration count.
   The seal example uses 10 000 iterations.
4. **(I)** Discrete Laplacian: choose between FORTRAN-style edge
   margins, cotangent weights, or simple neighbour averaging. Default
   recommendation: cotangent weights for the curved surface, mass-lumped
   on the dual cells. Document the choice in a docstring with the
   paper-equation citation.
5. **(I)** Update ordering: pure Jacobi for everything (gene network
   *and* forces), against a single start-of-iteration state copy. This
   diverges from the FORTRAN's mixed Jacobi/Gauss–Seidel scheme but
   matches the paper's prose and is naturally vectorisable.
6. **(I)** Cell-division mesh repair: replace the topology walk with
   either a global `scipy.spatial.Delaunay` re-triangulation on the
   $(x, y)$-projected coordinates each iteration (cheap; ~milliseconds
   for 60–200 cells), or local edge-flip restore. Either avoids all of
   the FORTRAN topology-walk accidents.
7. **(I)** Border-cell identification: vectorised neighbour-count
   threshold, recomputed each iteration. No half-plane dead band
   (`Swi` defaults to 0).
8. **(I)** Repulsion neighbour set: include both mesh neighbours
   *and* non-mesh-but-too-close cells (i.e. retain
   `pushingnovei` semantics). Justify in a comment as 'mesh stability,
   not in paper'.
9. **(I)** Apical normal $\hat{\mathbf{n}}$: angle-weighted normal
   from the one-ring of triangular faces. Document deviation from
   FORTRAN.
10. **(I)** Buoyancy gating: explicit sign-of-z-component check on
    $\hat{\mathbf{n}}$, vectorised.
11. **(I)** New-cell concentrations: average of mother cells (paper
    specifies). $d_i$ of a newly-born non-knot daughter of two knots
    should be set to **0**, not the average — the paper says 'the new
    cell is not a knot cell', which Path B should read as a hard
    reset of $d_i$ for that case (departure from FORTRAN, justified by
    paper text).
12. **(I)** Repulsion rest length on new edges: 1.0 (paper-implied
    constant; same as FORTRAN).
13. **(P)** Initial condition: 7-cell hexagon plus 14 mesenchymal
    cells (2 layers × 7), Act source at the six border epithelial
    cells, all other concentrations and $d_i$ = 0.
14. **(P)** Validation against `tests/golden_fortran/*.off` from the
    seal run, with the existing tolerance (±5% cell count, ±10%
    vertex envelope, no NaN/Inf). The paper's own validation regime
    is at the same level of rigour (cusp count, top-cusp angle, cusp
    tip positions), so equality of tolerance regimes is defensible.
15. **(I)** Knot daughter $d_i$: explicit zero-reset on knot-mother
    division (see item 11).

The two strongest cleanups that Path B can execute relative to `13.f90`,
without departing from the paper, are item 6 (replace topology walk
with global re-triangulation) and item 5 (pure Jacobi). Both directly
eliminate categories of FORTRAN-only accidents documented in
`docs/research/13f90-vs-humppa-divergences.md` and
`docs/findings/2026-05-04-path-a-fma-and-over-division.md`.
