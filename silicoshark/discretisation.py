"""Discretisation — the configurable choice points of silicoshark.

The 2010 paper under-specifies the model in many places; the FORTRAN
resolves each ambiguity one way; alternative numpy-idiomatic choices
exist for several. Path B v2 makes each choice point a named field
in this dataclass, with paper-citation, FORTRAN-citation, and
rationale recorded in `docs/research/discretisation-audit.md`.

Five named presets cover the canonical configurations:

- `LEGACY_FORTRAN`: reproduces 13.f90 semantically (modulo
  paper-vs-FORTRAN tensions where Path B has chosen the paper —
  e.g. the eq. 14 typo).
- `HUMPPA_LITERAL`: the upstream Catalan FORTRAN. Distinguishes
  13.f90's renaming-introduced divergences (txungu, Swi/parap[6])
  from genuine model behaviour. Currently identical to
  `LEGACY_FORTRAN`; differences will be teased out as the
  comparison study runs and the divergences become observable.
- `PAPER_2010`: the 2010 paper as written, with the eq. 14
  denominator typo corrected and choices made by best-faith
  reading where the paper is silent.
- `PAPER_LITERAL_2010`: the 2010 paper exactly as written,
  including the eq. 14 typo. For testing whether the typo
  materially changes biology.
- `PATH_B_DEFAULT`: the cleanest numpy-idiomatic choices.
  Currently identical to `PAPER_2010`. May diverge as we learn.

This module has no behaviour. It defines the dataclass and the
preset constants only. Wiring through `reaction.py`, `forces.py`,
`mesh.py`, `simulator.py` follows in subsequent phases.

Charter: `docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


# Type aliases for option literals — keeps the dataclass field
# annotations readable, and gives the type checker something to lean
# on at call sites.

LaplacianOpt = Literal['length_weighted', 'cotangent', 'fortran_margins']
UpdateOrderOpt = Literal['jacobi', 'gauss_seidel_forces', 'mixed_jacobi_gs']
TopologyOpt = Literal['delaunay_each_step', 'static_with_local_update']

Eq5ApplyOpt = Literal['all', 'interior_only']
RepFormOpt = Literal['hookean_signed', 'paper_gated']
AdhFormOpt = Literal['unit_vector', 'hookean_attraction']
RepNeighOpt = Literal['mesh', 'mesh_plus_all_close']

Eq14DenomOpt = Literal['act_typo', 'inh_corrected']
Eq17SourceOpt = Literal['act_concentration', 'act_rate_times_di']
Eq18SourceOpt = Literal['constant_k_sec', 'k_sec_times_di']

DiffAccumOpt = Literal['act', 'sec']
KnotGateOpt = Literal['none', 'first_border_cell']
KnotDaughterOpt = Literal['zero_reset', 'inherit_avg']

MesenchymeOpt = Literal['absent', 'per_column_z_layers']
BorderDefOpt = Literal['neighbour_count', 'topological_descendants']
LatticeOrientOpt = Literal['axial', 'fortran']


@dataclass(frozen=True)
class Discretisation:
    """All implementer-choice points carried through silicoshark.

    Frozen dataclass so a single instance can be passed by reference
    through the simulator without risk of mid-iteration mutation. To
    construct a study-specific configuration, start from a named
    preset and use `dataclasses.replace`:

        from silicoshark.discretisation import LEGACY_FORTRAN, replace
        cfg = replace(LEGACY_FORTRAN, eq14_denominator='act_typo')

    Each field is documented inline with:
      - Question being decided.
      - Citations for paper passage(s) and FORTRAN line(s).
      - The default (most defensible numpy-idiomatic choice) and
        why.

    Field defaults form the `PATH_B_DEFAULT` preset.
    """

    # --- Spatial discretisation -------------------------------------

    laplacian: LaplacianOpt = 'length_weighted'
    """Discrete form of `nabla^2 c` on the triangular mesh.

    Paper main p. 585 says 'finite volume on the triangular mesh'
    with 'flux proportional to contact area'. The formula is not
    given. Three options:

      - `length_weighted`: `L u[i] = sum_{j in N(i)} (u[j] - u[i]) / |p_j - p_i|`.
        Simplest Fick's-law discretisation; Path B v1 default.
      - `cotangent`: textbook discrete Laplace–Beltrami operator
        `L u[i] = (1/2A_i) sum_j (cot α_ij + cot β_ij) (u[j] - u[i])`.
        Paper-faithful for curved-surface diffusion.
      - `fortran_margins`: reproduces 13.f90's `apply_diffusion`
        with per-edge `pes` margins and z-layer `area_p` weighting.

    Citations: paper main p. 585; coreop2d.py:506.
    """

    # --- Temporal discretisation ------------------------------------

    update_order: UpdateOrderOpt = 'jacobi'
    """Jacobi vs Gauss–Seidel ordering of force / reaction terms.

    Paper prose suggests Jacobi (each `partial / partial t` is a
    function of the current state). FORTRAN is mixed: gene-network
    reactions are Jacobi (uses a copy `hq3d`); mechanical forces
    are Gauss–Seidel (repulsion, adhesion, traction each read and
    write `positions` in turn). Three options:

      - `jacobi`: pure Jacobi for everything.
      - `gauss_seidel_forces`: Jacobi for reactions, Gauss–Seidel
        for forces. Reproduces FORTRAN.
      - `mixed_jacobi_gs`: same as `gauss_seidel_forces` for
        symmetry; reserved for future variants.

    Citations: paper main p. 583 (fig. 1 caption);
    coreop2d.iteration().
    """

    # --- Mesh and topology ------------------------------------------

    topology: TopologyOpt = 'static_with_local_update'
    """How is the cell-cell adjacency graph maintained?

      - `delaunay_each_step`: re-triangulate via
        `scipy.spatial.Delaunay` every iteration. Eliminates the
        FORTRAN's topology-walk accidents but creates long-range
        convex-hull edges that drove the v1 mesh degeneration.
      - `static_with_local_update`: build the initial neighbour
        graph from Delaunay at t = 0; on cell division, locally
        rewire so the daughter inherits a subset of each parent's
        neighbours. Equivalent to FORTRAN's `add_cell` minus the
        topology-walk accidents.

    Citations: charter v2 §Architecture; v1 finding
    `docs/findings/2026-05-04-path-b-v1-mesh-degeneration.md`.
    Default is `static_with_local_update` because v1 evidence
    shows `delaunay_each_step` fails on the seal example.
    """

    division_per_step_cap: int | None = None
    """Maximum new cells emitted per iteration. `None` = unlimited.

    Used as a safeguard for the `delaunay_each_step` topology where
    re-triangulation can produce many simultaneous long edges.
    Inert under `static_with_local_update`.
    """

    division_max_edge: float | None = None
    """Reject division on edges longer than this multiple of rest
    length. `None` = no cap.

    Long Delaunay edges (e.g. across the convex hull when border
    cells have migrated outward) are topology artefacts, not
    biological neighbours. Inert under `static_with_local_update`.
    """

    division_total_cap: int | None = None
    """Maximum total cell count after division. `None` = unlimited.

    A pragmatic stand-in for FORTRAN 13.f90's `add_cell` topology-walk
    panic. The FORTRAN's `add_cell` aborts (panic; no cells added)
    when its neighbour-walk encounters a degenerate cell (one with
    too many or too few neighbours). On the seal example this fires
    after ~14 iterations and prevents further division beyond ~57
    cells, even though long edges continue to form. silicoshark's
    `static_with_local_update` topology has no such panic — its
    daughter-insertion logic always succeeds — so without this cap
    the seal example divides exponentially.

    `LEGACY_FORTRAN` sets this to 60 (~5% above the FORTRAN's 57-cell
    plateau) so the LEGACY_FORTRAN preset reproduces the FORTRAN's
    cell-count plateau within the seal smoke-test tolerance. The cap
    is a phenomenological match, not a faithful translation of
    FORTRAN's add_cell logic. Reproducing the topology-walk panic
    exactly is documented as future work in the discretisation audit.
    """

    # --- Eq. 5 epithelial growth ------------------------------------

    eq5_z_gate: bool = True
    """Restrict the eq. 5 `sum(unit_ij)` to neighbours strictly
    above in z (`z_j > z_i + 1e-4`)?

    The FORTRAN's `ep_growth_border_force` includes only neighbours
    with `b < -1e-4` (i.e. higher z). With flat-z initial conditions
    this means eq. 5 is force-free until cervical-loop curvature
    develops. Without the gate, perfect-symmetry cells produce
    near-zero `sum_unit` whose unit-vector normalisation amplifies
    floating-point noise into spurious forces (v1 finding).

    Default `True` because v1 evidence shows it stabilises the
    initial condition. The cleaner long-term answer may be a
    curvature-based outward-normal computation that does not need
    a z-gate.

    Citations: coreop2d.py:702; v1 finding.
    """

    eq5_apply_to: Eq5ApplyOpt = 'all'
    """Apply eq. 5 to all cells (paper) or interior cells only
    (FORTRAN)?

    The FORTRAN restricts eq. 5 to `i >= first_border_cell` (the
    inner ~2/3 of the lattice); border cells receive only the
    cervical-loop force. The 2010 paper has no such restriction.

    Default `'all'` per the paper.

    Citations: paper main p. 585; coreop2d.py:692.
    """

    # --- Eqs. 1-2 cell-cell mechanics -------------------------------

    rep_form: RepFormOpt = 'hookean_signed'
    """Eq. 1 form: always-on Hookean (signed) or paper-gated?

      - `hookean_signed`: `f = k_rep * (|d| - p_0) * d` for all
        edges. Smoothly transitions from repulsion (|d| < p_0)
        to weak attraction (|d| > p_0) with f = 0 at |d| = p_0.
        Avoids the floating-point edge-case at exact equilibrium.
      - `paper_gated`: same formula, active only when |d| < p_0.

    The paper's eq. 1 is written as the Hookean form; the gating
    interpretation comes from the paper's separate framing of
    eq. 1 (repulsion) and eq. 2 (adhesion) as 'the' two forces.

    Default `hookean_signed` because v1 evidence shows the gated
    form amplifies fp noise at exact-rest-length initial edges.

    Citations: paper main p. 585.
    """

    adh_form: AdhFormOpt = 'unit_vector'
    """Eq. 2 form: unit-vector pull (paper) or Hookean attraction
    (FORTRAN)?

      - `unit_vector`: `f = k_adh * d / |d|`, paper eq. 2 as written.
      - `hookean_attraction`: `f = k_adh * d`, FORTRAN's
        `repulse_neighbour` else-branch (Hookean, not unit-vector).

    Default `unit_vector` per the paper.

    Citations: paper main p. 585; coreop2d.py:870.
    """

    rep_neighbour_set: RepNeighOpt = 'mesh'
    """Repulsion only over mesh neighbours, or also against any
    cell whose centre is closer than ~1.4 units?

      - `mesh`: paper-only — eq. 1 sums over `j in N(i)` where
        `N` is the triangular-mesh neighbour set.
      - `mesh_plus_all_close`: FORTRAN runs an additional
        `pushingnovei` loop over all cells with distance < 1.4,
        regardless of mesh adjacency. Plausibly load-bearing for
        mesh-tangling prevention.

    Default `mesh` per the paper. Comparison study to test
    whether the extra repulsion is needed for stability.

    Citations: paper main p. 585; coreop2d.py:886 (repel_non_neigh).
    """

    # --- Reaction terms ---------------------------------------------

    eq14_denominator: Eq14DenomOpt = 'inh_corrected'
    """Eq. 14 inhibitor saturation term: paper text or FORTRAN
    correction?

      - `act_typo`: `1 + k_inh * [Act]` — paper as written. Reads
        as activator self-inhibition.
      - `inh_corrected`: `1 + k_inh * [Inh]` — FORTRAN's form.
        Inhibitor saturating its own response, which is the
        biologically meaningful reading the paper's prose
        elsewhere implies.

    Default `inh_corrected` per the charter (FORTRAN-on-this-one).

    Citations: paper main p. 586; coreop2d.py:645;
    docs/research/paper-review-2010-salazar-ciudad-jernvall.md
    §'Paper-vs-code findings'.
    """

    eq17_inh_source: Eq17SourceOpt = 'act_concentration'
    """Eq. 17 Inh production source for cells with `d_i >= Int`:

      - `act_concentration`: paper Eq. 17 — `[Inh]' = [Act] - Deg [Inh]`.
      - `act_rate_times_di`: FORTRAN — `[Inh]' = (rate of Act) * d_i - Deg [Inh]`.
        Uses the temporary Act-rate variable, not the Act
        concentration itself. Suspicious as a faithful translation
        of the paper.

    Default `act_concentration` per the paper.

    Citations: paper main p. 586 (eq. 17); coreop2d.py:647.
    """

    eq18_sec_source: Eq18SourceOpt = 'constant_k_sec'
    """Eq. 18 Sec production source for cells with `d_i >= Set`:

      - `constant_k_sec`: paper Eq. 18 — `[Sec]' = k_sec - Deg [Sec]`.
      - `k_sec_times_di`: FORTRAN — `[Sec]' = k_sec * d_i - Deg [Sec]`.
        Smooths production via `d_i` ramp.

    Default `constant_k_sec` per the paper.

    Citations: paper main p. 586 (eq. 18); coreop2d.py:652.
    """

    # --- Differentiation --------------------------------------------

    diff_accumulator: DiffAccumOpt = 'sec'
    """Eq. 6 differentiation accumulator: paper or FORTRAN?

      - `act`: paper Eq. 6 — `d_i' = k_dff [Act]`.
      - `sec`: FORTRAN — `d_i' = k_dff [Sec]`. Path B v1 default
        per charter v1.

    Default `sec` per the v1 charter (this is the one paper-vs-
    FORTRAN gap where v1 followed the FORTRAN; the comparison
    study is what justifies revisiting this).

    Citations: paper main p. 585 (eq. 6); coreop2d.apply_differentiation;
    docs/findings/2026-05-04-differentiation-uses-sec-not-act.md.
    """

    knot_threshold_gate: KnotGateOpt = 'none'
    """Restrict knot formation to cells with index above some
    threshold?

      - `none`: any cell with `[Act] >= 1` becomes a knot (paper).
      - `first_border_cell`: only cells with `i >= first_border_cell`
        become knots. FORTRAN's behaviour. Excludes the
        innermost cells from knot formation, which is biologically
        odd but observable in 13.f90.

    Default `none` per the paper.

    Citations: paper main p. 586; coreop2d.py:640.
    """

    knot_daughter_di: KnotDaughterOpt = 'zero_reset'
    """When a knot mother divides, what `d_i` does the (non-knot)
    daughter inherit?

      - `zero_reset`: hard reset to 0. Reads the paper's 'the new
        cell is not a knot cell' as also resetting differentiation.
      - `inherit_avg`: 0.5 * (d_i_a + d_i_b). FORTRAN's behaviour;
        a daughter of two knots can keep `d_i` near 1, freezing it.

    Default `zero_reset` per the paper.

    Citations: paper main p. 585; coreop2d.add_cell.
    """

    # --- Mesenchyme -------------------------------------------------

    mesenchyme: MesenchymeOpt = 'per_column_z_layers'
    """Model the underlying mesenchyme explicitly?

      - `absent`: epithelium only. Path B v1 shortcut. Inh and Sec
        diffuse only laterally. Buoyancy and mesenchyme-thickening
        terms have no Sec source. Acceptable for parameter sets
        where Sec stays 0 (e.g. seal.txt), broken for those where
        it does not.
      - `per_column_z_layers`: each epithelial cell carries a
        column of mesenchymal Inh/Sec values. Vertical diffusion
        between layers. Matches the paper's 'two layers underneath
        each epithelial cell'. FORTRAN-compatible.

    Default `per_column_z_layers`. Pinning to `absent` is a
    Path B v1-style fast-prototype option.

    Citations: paper main p. 586; coreop2d.apply_diffusion z-layer
    handling.
    """

    n_mesenchyme_layers: int = 2
    """Number of mesenchymal layers per epithelial column.

    Paper specifies 2. FORTRAN allocates 4 (epi + 2 mesenchyme +
    substrate edge); the substrate edge is a boundary condition,
    not a real layer. Inert when `mesenchyme == 'absent'`.

    Citations: paper main p. 586; coreop2d.max_z_layers.
    """

    # --- Border identification and biases ---------------------------

    border_definition: BorderDefOpt = 'neighbour_count'
    """How are border cells identified each iteration?

      - `neighbour_count`: cells with fewer than 6 mesh neighbours.
        Geometric, recomputed each step from the current
        triangulation. Vectorised, cheap.
      - `topological_descendants`: cells whose ancestry traces to
        the initial-condition border cells. Maintained
        incrementally on division. Paper's 'Growth biases'
        definition.

    The 2010 paper conflates these definitions. FORTRAN's biases
    use the topological reading (`first_border_cell`), but
    `update_border_cells` then promotes geometrically-emerging
    cells to border status — so the FORTRAN is also mixed.

    Default `neighbour_count` because it composes cleanly with
    the static-topology code path and matches the simpler reading
    of the paper text.

    Citations: paper main p. 583 (fig. 1 caption), p. 585 ('Growth
    biases'); coreop2d.update_border_cells.
    """

    lattice_orientation: LatticeOrientOpt = 'axial'
    """Initial hex lattice orientation.

      - `'axial'`: conventional axial layout. Ring 1 cells at angles
        0°, 60°, 120°, 180°, 240°, 300°. Outer-ring cells include
        (rad, 0), (-rad, 0); none at (0, rad). The v1 simulator was
        developed against this orientation; the v1-byte-identity
        baseline depends on it.
      - `'fortran'`: rotated 30° counter-clockwise so the outer ring
        includes cells at exactly (0, rad), (0, -rad). Matches
        `13.f90`'s `allocate_initial_state`. Required for
        LEGACY_FORTRAN's `border_bias_x_zero_quirk` to have any
        effect on the seal example: that quirk excludes cells at
        x exactly == 0 from the Bgr z-multiplier, anchoring those
        cells low in z and creating the gradient that drives the
        FORTRAN's 57-cell plateau and z=125 envelope. Without the
        rotation, no cells sit at x=0 and the quirk is inert,
        producing a flat border-z curve and ~50% z-envelope error.

    This field was added in Path B v2 A6 (2026-05-05) when the
    lattice-rotation dependency on the seal example's z dynamics
    was identified during the LEGACY_FORTRAN replication run.
    """

    border_bias_x_zero_quirk: bool = True
    """Preserve FORTRAN's update_cell_position quirk where cells
    with x exactly == 0 don't receive Pbi/Abi/Bgr multipliers?

    The FORTRAN's per-half-plane bias logic uses `if x > 0` and
    `elif x < 0`, which omits cells at x == 0. For the seal
    initial condition with y-axis cells at x = 0 exactly, this
    means those cells don't get the Bgr z-amplification. The
    chain-of-accidents docstring at the top of `coreop2d.py`
    documents this as load-bearing for the 57-cell plateau.

    Default `True` to preserve compatibility with the FORTRAN
    goldens. The corrected behaviour (`False`) is biologically
    cleaner but produces different numerics on the seal example.

    Citations: coreop2d.update_cell_position;
    docs/research/13f90-vs-humppa-divergences.md.
    """


# --- Named presets ---------------------------------------------------
#
# Each preset is what the comparison study uses to test a specific
# question. The aim string is the headline; the docstring spells out
# what the preset isolates relative to its neighbours in the matrix.


PATH_B_DEFAULT = Discretisation()
"""The cleanest numpy-idiomatic choices for production use.

Currently identical to `PAPER_2010` (the field defaults of
`Discretisation` reflect the paper's choices where it speaks, and
the most defensible numpy-idiomatic choice where it is silent).
The preset will diverge from `PAPER_2010` if and only if the
comparison study reveals a numpy-idiomatic choice that is
preferable to the paper's literal reading on biological grounds.
"""


PAPER_2010 = Discretisation(
    # Identical to PATH_B_DEFAULT until/unless they diverge. We keep
    # this preset distinct so a later edit to one does not silently
    # alter the other.
)
"""The 2010 paper as written, with the eq. 14 denominator typo
corrected and choices made by best-faith reading where the paper
is silent. Cusp-formation behaviour reported here is what the
paper claims its model produces.
"""


PAPER_LITERAL_2010 = Discretisation(
    eq14_denominator='act_typo',
)
"""The 2010 paper exactly as written, including the eq. 14 typo.
The comparison with PAPER_2010 quantifies whether the typo
materially changes biology; if the dynamics are
indistinguishable, the typo is harmless. If they differ, the
typo would have produced different published results had the
paper's text been used as the implementation.
"""


LEGACY_FORTRAN = Discretisation(
    laplacian='fortran_margins',
    update_order='gauss_seidel_forces',
    topology='static_with_local_update',
    eq5_z_gate=True,
    eq5_apply_to='interior_only',
    rep_form='paper_gated',
    adh_form='hookean_attraction',
    rep_neighbour_set='mesh_plus_all_close',
    eq14_denominator='inh_corrected',
    eq17_inh_source='act_rate_times_di',
    eq18_sec_source='k_sec_times_di',
    diff_accumulator='sec',
    knot_threshold_gate='first_border_cell',
    knot_daughter_di='inherit_avg',
    border_definition='topological_descendants',
    border_bias_x_zero_quirk=True,
    lattice_orientation='fortran',
    division_total_cap=60,
)
"""Reproduces 13.f90 semantically, including the FORTRAN's
implementer choices and (where they are load-bearing) its
accidents. This is the replication anchor for the comparison
study: every other preset is reported relative to this one's
seal-example output.

Caveats:
- The eq. 14 denominator stays at `inh_corrected` (not the
  paper's `act_typo`) because the FORTRAN itself uses the
  corrected form; this is one place where 13.f90 and the paper
  text disagree and 13.f90 is biologically right.
- 13.f90's renaming-introduced quirks (txungu, Swi/parap[6])
  are *not* preserved here because they live in the topology
  walk that v2 replaces with static-update. Set `topology` to
  `delaunay_each_step` to test what happens without those
  quirks under the v2 code path.
- `border_definition='topological_descendants'` was added in
  Path B v2 A6 (2026-05-05) when the cervical-loop's spread to
  geometric `<6 neighbours` daughters was found to drive
  exponential over-division on the seal example. The FORTRAN
  uses topological border (`first_border_cell` block + new
  cells whose `i, k < first_border_cell`), and LEGACY_FORTRAN
  must match. This is the only post-A5 modification to the
  preset; it is a bug fix (the field existed in A5 but was not
  wired through to forces.py).
"""


HUMPPA_LITERAL = Discretisation(
    laplacian='fortran_margins',
    update_order='gauss_seidel_forces',
    topology='static_with_local_update',
    eq5_z_gate=True,
    eq5_apply_to='interior_only',
    rep_form='paper_gated',
    adh_form='hookean_attraction',
    rep_neighbour_set='mesh_plus_all_close',
    eq14_denominator='inh_corrected',
    eq17_inh_source='act_rate_times_di',
    eq18_sec_source='k_sec_times_di',
    diff_accumulator='sec',
    knot_threshold_gate='first_border_cell',
    knot_daughter_di='inherit_avg',
    border_definition='topological_descendants',
    border_bias_x_zero_quirk=True,
    lattice_orientation='fortran',
    division_total_cap=60,
)
"""The upstream Catalan FORTRAN (`humppa_translate.f90`).

Currently identical to `LEGACY_FORTRAN` because all the
13.f90-only divergences catalogued in
`docs/research/13f90-vs-humppa-divergences.md` live in the
topology-walk machinery that v2 replaces with static-update
(rendering the divergences inert) — except for the `Swi` /
parap[6] handling, which v1 already accepts as a parameter
and routes nowhere yet. As humppa-specific behaviour becomes
observable through future Discretisation fields, this preset
will diverge from `LEGACY_FORTRAN`.
"""


ALL_PRESETS: dict[str, Discretisation] = {
    'PATH_B_DEFAULT': PATH_B_DEFAULT,
    'PAPER_2010': PAPER_2010,
    'PAPER_LITERAL_2010': PAPER_LITERAL_2010,
    'LEGACY_FORTRAN': LEGACY_FORTRAN,
    'HUMPPA_LITERAL': HUMPPA_LITERAL,
}
"""Registry of named presets. Used by the CLI's `--preset` flag and
by the comparison-study runner to enumerate configurations.
"""
