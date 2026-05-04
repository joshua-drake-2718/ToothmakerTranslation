import sys
import math
import vector
import numpy as np


"""
Translation note — chain of accidental behaviours in the FORTRAN reference

This file is a Python translation of `13.f90`, which is itself an
English-renamed copy of `humppa_translate.f90` from `jernvall-lab/
ToothMaker` (see `docs/research/path-b-research.md` for the full
provenance chain). Where this docstring refers to 'the FORTRAN', it
means `13.f90` specifically; the upstream Catalan-named ancestor
(humppa) and its smaller fork (`tgrohens/toothmaker`) differ in the
ways noted below.

`13.f90` reaches a stable numerical result on the seal example only
through a chain of compile-environment- and renaming-process accidents
that compose. The Python translation preserves these accidents verbatim
to match the captured goldens. The chain, in order of execution:

  1. `add_cell`'s txungu else-branch executes a `rtt:` loop whose
     body never runs because of a renaming-introduced inconsistency at
     `13.f90:1156-1158`: the loop assigns to `jjj` (`jjj=pillats(kkk)`)
     but its conditionals test `jjjj` (a loop-invariant value set
     before the loop). The branch fires ~100 times per save block and
     always produces `temp_new_neigh[i, :] = [ini, 0, 0, ...]` —
     single-neighbour cells.

     Provenance: this inconsistency is *not* in the upstream FORTRAN.
     `humppa_translate.f90:1336-1345` and `tgrohens/coreop2d.f90:
     1276-1285` both write `jjjj=pillats(kkk)` and test `jjjj` —
     internally consistent, the loop runs. The mismatch was
     introduced when `13.f90` was renamed from Catalan to English:
     the assignment was renamed `jjjj` → `jjj` but the conditionals
     were left untouched. In humppa, the `rtt:` loop is live code
     and produces multi-neighbour cells in this branch.

     The `goto 899` immediately after the loop (`13.f90:1166`,
     humppa `:1346`, tgrohens `:1286`) is identical in all three
     files and skips a large disabled block to the continuation at
     label 899. That `goto` is unrelated to the renaming typo.

  2. `calculate_margins` has label 77 *outside* the `if (neigh /= 0)`
     block (FORTRAN 13.f90:393), so for every empty neighbour slot it
     re-runs `border(i, j, 1:3) = a/cont` using the LAST computed values
     of `a, b, c, cont` from the previous iteration. Empty slots get
     filled with the same border position as the previous non-empty
     slot. Present in humppa and tgrohens identically; not a
     translation artefact.

  3. `apply_diffusion`'s fallback branch then computes a cross product of
     two bit-equal vectors (border[i, j, :] - positions[i, :]) ×
     (border[i, 1, :] - positions[i, :]). Mathematically zero.

  4. With FMA emission enabled (default on ARM at -O2; Intel needs
     Haswell+ and `-march=native` or `-mfma`), the cross product
     `uy*dz - uz*dy` is computed as `fma(uy, dz, -uz*dy)` with one
     rounding step instead of two, leaving an ulp-level non-zero
     residual (~1e-18). Verified in PR #14: rebuilding the FORTRAN with
     `-ffp-contract=off` produces the same NaN cascade as Python.

  5. That ε noise makes `sum_a` tiny but positive, so the normalisation
     `area_bottom /= sum_a` is finite. Vertical diffusion proceeds at
     `area_bottom = 0.5` (first norm) and `area_bottom = 1.0` (second norm
     after the renormalisation block) — the same values the mathematical
     limit gives.

If any link 1-4 breaks in `13.f90` (txungu inconsistency "fixed", label 77
moved, fall-through patched, FMA disabled), `13.f90` also NaN-cascades
— confirmed by PR #14's `-ffp-contract=off` rebuild. Pure Python breaks
link 4 by default (`*` and `-` don't use FMA). The Python translation
handles link 4 by substituting the explicit limit at link 5
(`if sum_a == 0: area_bottom = 0.5/1.0`); see the FMA-LIMIT GUARD
comment in `apply_diffusion`.

humppa and tgrohens, run on equivalent inputs, would not necessarily
reproduce `13.f90`'s goldens — link 1 is `13.f90`-specific, and the
txungu else-branch's behaviour differs. Path B (paper-faithful
re-implementation) targets the upstream model behaviour, not
`13.f90`'s renaming-typo accidents.

The original handover claimed that fixing link 5 alone exposed a
single coupled "over-division" bug (cells multiplying exponentially:
37 → 49 → 79 → 139 → 199 → 229 → 289 → 349 in 70 iterations vs
FORTRAN's plateau at 57). PR #16 traced that to TWO unrelated bugs
that the NaN cascade had been masking. Both turned out to be
artefacts of *our* Python migration (and not present in `13.f90` or
its upstream); the fixes restore the FORTRAN structure verbatim:

  - `update_cell_position`: an earlier Python migration had moved
    `forces[:, 2] *= Bgr` outside the x-sign branches. `13.f90:
    863-872`, `humppa:1044-1070` and `tgrohens:984-1010` all nest
    it inside each branch, so cells with x == 0 don't receive the
    boost. PR #16 restored the FORTRAN nesting.
  - `add_cell`: an earlier Python migration had added an explicit
    `iiii = kkk = 0` reset that the FORTRAN doesn't have — the
    FORTRAN relies on local-variable lifetime to carry `iiii`
    forward across topology-walk iterations, which is the
    mechanism that makes degenerate walks panic instead of
    adding bogus cells (`13.f90:1045-1054, 1072-1081`,
    `humppa:1226-1252`). PR #16 removed the spurious reset.

These two were independent of FMA and independent of the link-1-to-4
chain. The chain documented in steps 1-4 remains real and load-bearing
in `13.f90` at exactly one site (the apply_diffusion divide); it is not
a global equilibrium effect.

Each step in the chain has a `KNOWN ISSUE` (or, for step 4-5, an
`FMA-LIMIT GUARD`) comment in this file at the relevant code site.
See PRs #11-#16 for the full investigation; see
`docs/research/cpp-port-review.md` and `docs/research/tgrohens-review.md`
for the cross-fork comparisons that established the link-1 provenance.
"""


class Coreop2d():
    positions:  np.ndarray
    border:     np.ndarray
    neigh:      np.ndarray
    knots:      np.ndarray
    num_neigh:  np.ndarray
    q3d:        np.ndarray
    diff_state: np.ndarray
    forces:     np.ndarray
    force_snapshot: np.ndarray
    px: np.ndarray
    py: np.ndarray
    pz: np.ndarray

    num_active_cells: int
    num_all_cells: int
    max_z_layers: int
    nv_max: int = 30
    Rad: int
    num_species_in_q3d: int = 3
    temps: int
    npas: int
    la: float = 1

    # true parameters
    Egr: float      # Egr (epithelial proliferation rate, "tacre")
    Mgr: float      # Mgr (mesenchymal proliferation rate, "tahor")
    Rep: float      # Rep (Young's modulus, stiffnes, "elas")

    Adh: float      # Adh (traction between neigh, "Adh")
    Act: float      # Act (activator auto-activation, "acac")
    Inh: float      # Inh (inhibition of activator, "iMac")

    Sec: float      # Sec (growth factor secretion rate, "ih")

    Da: float       # Activator (BMP4) diffusion rate
    Di: float       # Inhibitor (SHH) diffusion rate
    Ds: float       # FGF (Sec) diffusion rate

    Int: float      # Int (differentiation threshold, "us")
    Set: float      # Set (differentiation threshold, "ud")

    Boy: float      # mesenchymal buoyancy (Eq. 13). In subroutine BoyForce. Paper Eq. 13. Mech. resistance of mesench. to invagination of epithelium . 
    Dff: float      # differentiation rate (Eq. 6). Formerly tadif). Equation 6. ! Originally difq2d(3). 
    Bgr: float      # border growth to invagination coupling (Used In Subroutine "updatecellposition". Originally difq2d(4).)
    dmax: float = 2 # maximum distance before making a new node
    Abi: float      # bias posterior, anterior, lingual, buccal pabl (NOTE: Abi="bip", Pbi="bia", Bbi="bib", Lbi="bil")
    Pbi: float
    Bbi: float
    Lbi: float
    Deg: float      # Deg (protein degradation rate, "mu")
    Dgr: float      # Dgr (sharpness maxima, indicates how the epithelium pulls down if we do not have pressure of the mesenchyme, "tazmax")
                    # it doesn't affect much unless it's very low, 100 is a good value
    Ntr: float      # Ntr (mechanical traction from the borders to the nucleus. "radibi")
    Bwi: float      # radius of the center where we apply the bias ap
    ina: float      # initial activator concentration
    umgr: float     # basal mesenchymal prolif. rate (independent of Sec level)

    # semi-constants
    umelas: float
    maxcels: int

    # implementation detail
    delta: float = 0.05
    v_min: float = 0.15
    nca: int
    i_centre: int
    centre: int
    first_border_cell: int
    x: float
    xx: float
    y: float
    yy: float
    csu: float
    ssu: float
    csd: float
    ssd: float
    cst: float
    sst: float
    csq: float
    ssq: float
    csc: float
    ssc: float
    css: float
    sss: float
    num_new_cells: int
    num_border_cells: int
    n_map: int
    border_cell_list: np.ndarray
    m_map: np.ndarray

    # for visualisation
    vlinies: int
    vrender: int
    showborders: int
    vvec: int
    vvecx: int
    vveck: int
    vex: int
    vn: int
    nc: int
    pin: int
    pina: int
    nivell: int
    kko: int
    iti: int
    iteration_total: int
    nfpro: int


    @classmethod
    def initial_conditions(cls):
        cls.max_z_layers = 4
        cls.temps = 0
        cls.npas = 1
        # for visualization
        cls.vlinies = 0
        cls.vrender = 1
        cls.showborders = 0
        cls.vvec = 0
        cls.vvecx = 0
        cls.vveck = 0
        cls.vex = 0
        cls.vn = 0
        cls.pin = 1
        cls.pina = 1
        cls.nivell = 1


    @classmethod
    def allocate_initial_state(cls):
        cls.umelas = 1 - cls.Rep

        # Variable initial values
        Rad = max(0, cls.Rad)   # if j is zero then num_all_cells=1 to begin with 
        cls.num_all_cells =    3 * Rad * (Rad + 1) + 1
        cls.num_active_cells = 3 * (Rad - 1) * Rad + 1
        s60 = math.sqrt(3) / 2
        cls.csu = 0.0;  cls.ssu = 1.0
        cls.csd = s60;  cls.ssd = 0.5
        cls.cst = s60;  cls.sst =-0.5
        cls.csq = 0.0;  cls.ssq =-1.0
        cls.csc =-s60;  cls.ssc =-0.5
        cls.css =-s60;  cls.sss = 0.5

        # allocations
        cls.positions =     np.zeros((cls.num_all_cells,3))
        cls.neigh =         np.zeros((cls.num_all_cells, cls.nv_max), dtype=np.int32)
        cls.border =        np.zeros((cls.num_all_cells, cls.nv_max,8))
        cls.knots =         np.zeros((cls.num_all_cells), dtype=np.int32)
        cls.num_neigh =     np.zeros((cls.num_all_cells), dtype=np.int32)
        cls.diff_state =    np.zeros((cls.num_all_cells))
        cls.q3d =           np.zeros((cls.num_all_cells, cls.max_z_layers, cls.num_species_in_q3d))
        cls.m_map =         np.zeros((cls.Rad), dtype=np.int32)
        cls.forces =        np.zeros((cls.num_all_cells,3))
        cls.force_snapshot =    np.zeros((cls.num_all_cells,3))
        cls.border_cell_list =  np.zeros((cls.Rad), dtype=np.int32)

        # visualisation matrices
        cls.px = np.zeros((cls.num_all_cells))
        cls.py = np.zeros((cls.num_all_cells))
        cls.pz = np.zeros((cls.num_all_cells))

        # zero values (np.zeros already zeros arrays; sentinel for neigh is -1)
        cls.neigh[:] = -1

        # initial values
        cls.positions[0, 0] = 0
        cls.positions[0, 1] = 0
        cls.positions[0, 2] = 1
        cls.nca = 1
        cls.num_neigh[:] = 6

        # initial mesh values
        for i_centre in range(cls.num_active_cells): # al
            x = cls.positions[i_centre, 0]
            y = cls.positions[i_centre, 1]

            cls.initialise_cell_positions(i_centre, x+cls.csu*cls.la, y+cls.ssu*cls.la, 0, 3)
            cls.initialise_cell_positions(i_centre, x+cls.csd*cls.la, y+cls.ssd*cls.la, 1, 4)
            cls.initialise_cell_positions(i_centre, x+cls.cst*cls.la, y+cls.sst*cls.la, 2, 5)
            cls.initialise_cell_positions(i_centre, x+cls.csq*cls.la, y+cls.ssq*cls.la, 3, 0)
            cls.initialise_cell_positions(i_centre, x+cls.csc*cls.la, y+cls.ssc*cls.la, 4, 1)
            cls.initialise_cell_positions(i_centre, x+cls.css*cls.la, y+cls.sss*cls.la, 5, 2)
        for i in range(1, cls.num_active_cells):
            for j in range(cls.nv_max):
                if cls.neigh[i, j] >= cls.num_active_cells:
                    cls.neigh[i, j] = cls.num_all_cells
        for k in range(3):
            for i in range(1, cls.num_active_cells):
                for j in range(cls.nv_max-1):
                    if cls.neigh[i, j] == cls.num_all_cells and cls.neigh[i, j+1] == cls.num_all_cells:
                        for jj in range(j, cls.nv_max-1):
                            cls.neigh[i, jj] = cls.neigh[i, jj+1]
        for i in range(1, cls.num_active_cells):
            k = False
            for j in range(cls.nv_max):
                if cls.neigh[i, j] == cls.num_all_cells:
                    if k:
                        cls.neigh[i, j] = -1
                        break
                    k = True
        cls.positions = cls.positions.round(decimals=14)
        for i in range(cls.num_active_cells):
            for j in range(3):
                if abs(cls.positions[i, j]) < 1e-14:
                    cls.positions[i, j] = 0
        
        # original distance calculation between nodes

        # investment so that the former are on the margins

        cv = cls.neigh.copy()
        temp_positions = cls.positions.copy()
        for i in range(cls.num_active_cells):
            ii = cls.num_active_cells - 1 - i
            cls.neigh[i, :] = cv[ii, :]
            cls.positions[i, :] = temp_positions[ii, :]

        cv = cls.neigh.copy()
        for i in range(cls.num_active_cells):
            ii = cls.num_active_cells - 1 - i
            for j in range(cls.num_active_cells):
                for jj in range(cls.nv_max):
                    if cv[j, jj] == i:
                        cls.neigh[j, jj] = ii

        cls.calculate_margins()
        cls.num_neigh[:] = 3
        cls.num_neigh[0] = 6
        cls.border[:, :, 3:5] = cls.la
        cls.centre = cls.num_active_cells
        cls.first_border_cell = 6 * (cls.Rad - 1)
        cls.border_cell_list[:] = 0
        cls.m_map[:] = 0
        # margin values ATTENTION NOT SCALABLE WITH Rad, only for Rad=2
        for i in range(cls.Rad):
            cls.m_map[i] = i
        for i in range(cls.Rad):
            cls.border_cell_list[i] = cls.first_border_cell // 2 + i
        cls.num_border_cells = cls.Rad
        cls.n_map = Rad
        cls.calculate_margins()
        cls.q3d[:] = 0

    @classmethod
    def initialise_cell_positions(cls, i_centre: int, x: float, y: float, j: int, jj: int):
        for i in range(cls.nca):
            if i == i_centre: continue
            x_approx = math.isclose(x, cls.positions[i, 0], abs_tol=1e-6)
            y_approx = math.isclose(y, cls.positions[i, 1], abs_tol=1e-6)
            if x_approx and y_approx:
                for ii in range(cls.nv_max):
                    if cls.neigh[i_centre, ii] == i: return
                cls.neigh[i_centre, j] = i
                cls.neigh[i, jj] = i_centre
                cls.num_neigh[i] += 1
                cls.num_neigh[i_centre] += 1
                return
        cls.neigh[i_centre, j] = cls.nca
        cls.neigh[cls.nca, jj] = i_centre
        cls.positions[cls.nca, 0] = x
        cls.positions[cls.nca, 1] = y
        cls.positions[cls.nca, 2] = 1
        cls.num_neigh[i_centre] += 1
        cls.num_neigh[cls.nca] += 1
        cls.nca += 1

    @classmethod
    def calculate_margins(cls):
        count: float
        kl: int
        ii: int
        iii: int
        a: float
        b: float
        c: float

        # always executes goto 77, but sometimes early
        def calculate_margins_a():
            nonlocal count, kl, ii, iii, a, b, c
            for jj in range(j-2, -1, -1):  # FORTRAN: do jj=j-1,1,-1 (1-based) → j-2..0
                if cls.neigh[i, jj] != -1:
                    if cls.neigh[i, jj] < cls.num_active_cells:
                        ii = cls.neigh[i, jj]
                        a += cls.positions[ii, 0]
                        b += cls.positions[ii, 1]
                        c += cls.positions[ii, 2]
                        count += 1
                    return # goto 77
            for jj in range(cls.nv_max-1, j-1, -1):  # FORTRAN: do jj=nvmax,j+1,-1 → nv_max-1..j
                if cls.neigh[i, jj] != -1:
                    if cls.neigh[i, jj] < cls.num_active_cells:
                        ii = cls.neigh[i, jj]
                        a += cls.positions[ii, 0]
                        b += cls.positions[ii, 1]
                        c += cls.positions[ii, 2]
                        count += 1
                    return # goto 77

        # returns early if it executes goto 77
        def calculate_margins_b():
            nonlocal count, kl, ii, iii, a, b, c
            looping = True
            while looping:
                # 66
                if ii == i:
                    return # goto 77
                kl += 1
                if kl > 100:
                    calculate_margins_a()
                    return # goto 77
                a += cls.positions[ii, 0]
                b += cls.positions[ii, 1]
                c += cls.positions[ii, 2]
                count += 1
                looping = calculate_margins_c() # goto 66 or 77

        # returns true if it executes goto 66, false if goto 77
        def calculate_margins_c() -> bool:
            nonlocal count, kl, ii, iii, a, b, c
            jjj = 0
            for jj in range(cls.nv_max):
                if cls.neigh[ii, jj] == iii:
                    jjj = jj
                    break
            for jj in range(jjj+1, cls.nv_max):
                if cls.neigh[ii, jj] != -1:
                    if cls.neigh[ii, jj] >= cls.num_active_cells:
                        return False # goto 77
                    iii = ii
                    ii = cls.neigh[iii, jj]
                    return True # goto 66
            for jj in range(0, jjj-1):  # FORTRAN: do jj=1,jjj-1 (forward loop, 1-based) → 0..jjj-2
                if cls.neigh[ii, jj] != -1:
                    if cls.neigh[ii, jj] >= cls.num_active_cells:
                        return False # goto 77
                    iii = ii
                    ii = cls.neigh[iii, jj]
                    return True # goto 66
            return False # goto 77

        cls.border[:, :, 0:3] = 0
        for i in range(cls.num_active_cells):
            kl = 0
            for j in range(cls.nv_max):
                if cls.neigh[i, j] != -1:
                    a = cls.positions[i, 0]
                    b = cls.positions[i, 1]
                    c = cls.positions[i, 2]
                    count = 1
                    ii = cls.neigh[i, j]
                    iii = i
                    if ii >= cls.num_active_cells:
                        calculate_margins_a()
                    else:
                        calculate_margins_b()
                    # 77 — KNOWN ISSUE: FORTRAN 13.f90:393 has this assignment
                    # OUTSIDE the `if (neigh /= 0)` block, so for empty
                    # neighbour slots (`neigh(i, j) == 0`) the assignment runs
                    # using the LAST computed values of `a, b, c, count` from
                    # the previous non-empty j. Effect: a 1-neighbour cell
                    # (which `add_cell`'s txungu else-branch deliberately
                    # produces) gets `border[i, j, 0:3]` filled with the same
                    # position for every j, rather than zeros.
                    #
                    # Verified by FORTRAN trace (PR #15): for cell 42 with
                    # neigh = [2, 0, 0, ...], FORTRAN's border(42, 1, 1:3) =
                    # border(42, 5, 1:3) = border(42, 10, 1:3) — all the same
                    # real position. Python's structure (assignment INSIDE
                    # the if-block, as the migration moved it) sets only
                    # border[42, 0, 0:3] and leaves the rest at zero.
                    #
                    # This difference is then load-bearing for `apply_diffusion`
                    # — see the KNOWN ISSUE comment near `area_bottom /= sum_a`.
                    # If border[i, j, 0:3] were also filled here for empty
                    # slots, the apply_diffusion fallback would have non-zero
                    # input vectors and the cross product would be small but
                    # non-zero (without needing FMA). That alone wouldn't fix
                    # the simulator (the over-division bug is downstream of
                    # apply_diffusion), but it is one of the chain's load-
                    # bearing accidents — see the module-level docstring.
                    cls.border[i, j, 0] = a / count
                    cls.border[i, j, 1] = b / count
                    cls.border[i, j, 2] = c / count


    @classmethod
    def apply_diffusion(cls):
        pes = np.zeros((cls.num_all_cells, cls.nv_max))  # area of contact between i and neighbor (i, j)
        area_p = np.zeros((cls.num_all_cells, cls.nv_max))   # fraction of a cell’s surface area that is in contact with the “bottom” (the z-direction neighbors / substrate) and is used to weight vertical diffusion.areabottom
        sum_a: float
        area_bottom: float
        hq3d = np.zeros((cls.num_all_cells, cls.max_z_layers, cls.num_species_in_q3d))   # num_all_cells = number of cells (within “real radius”), max_z_layers = z depth, num_species_in_q3d = 3 (constant)
        ux: float
        uy: float
        uz: float
        dx: float
        dy: float
        dz: float
        ua: float
        ub: float
        uc: float

        for i in range(cls.num_active_cells):
            pes[i, :] = 0
            area_p[i, :] = 0
            for j in range(cls.nv_max): # ui
                if cls.neigh[i, j] != -1:
                    ua = cls.positions[i, 0]
                    ub = cls.positions[i, 1]
                    uc = cls.positions[i, 2]
                    for jj in range(j+1, cls.nv_max):
                        if cls.neigh[i, jj] != -1:
                            pes[i, j] = vector.a_distance_between(cls.border[i, j, :], cls.border[i, jj, :])
                            ux = cls.border[i, j, 0] - ua
                            uy = cls.border[i, j, 1] - ub
                            uz = cls.border[i, j, 2] - uc
                            dx = cls.border[i, jj, 0] - ua
                            dy = cls.border[i, jj, 1] - ub
                            dz = cls.border[i, jj, 2] - uc
                            area_p[i, j] = 0.5 * vector.a_magnitude(vector.cross_product(ux, uy, uz, dx, dy, dz))
                            break # cycle ui (skip fallback)
                    else:
                        # FORTRAN fallback: no valid jj > j was found, use border[i, 0]
                        pes[i, j] = vector.a_distance_between(cls.border[i, j, :], cls.border[i, 0, :])
                        ux = cls.border[i, j, 0] - ua
                        uy = cls.border[i, j, 1] - ub
                        uz = cls.border[i, j, 2] - uc
                        dx = cls.border[i, 0, 0] - ua
                        dy = cls.border[i, 0, 1] - ub
                        dz = cls.border[i, 0, 2] - uc
                        area_p[i, j] = 0.5 * vector.a_magnitude(vector.cross_product(ux, uy, uz, dx, dy, dz))
            area_bottom = area_p[i, :].sum()
            sum_a = pes[i, :].sum() + 2 * area_bottom
            # FMA-LIMIT GUARD: when sum_a == 0 (degenerate cell with
            # collapsed border geometry — produced by add_cell's txungu
            # else-branch, whose `rtt:` loop is dead code due to a `jjjj`
            # vs `jjj` typo at FORTRAN 13.f90:1158 — see the matching
            # KNOWN ISSUE near `temp_new_neigh[i, :] = -1` later in this
            # file), this would divide exact zero by exact zero and NaN
            # would propagate to every cell via diffusion coupling.
            #
            # FORTRAN avoids the NaN by accident: with FMA emission
            # enabled, the upstream cross product (uy*dz − uz*dy) on
            # bit-equal vectors compiles to fma(uy, dz, −uz*dy) with one
            # rounding step instead of two, leaving an ulp-level residual
            # (~1e-18) instead of an exact zero. That ε makes sum_a tiny
            # but positive, so the divide is finite. Verified in PR #14
            # by rebuilding the FORTRAN with `-ffp-contract=off`: the
            # FORTRAN binary then produces the same NaN cascade as Python.
            #
            # FMA emission depends on the compile environment: ARMv8 has
            # FMA in the base ISA so any -O2 build uses it; Intel needs
            # Haswell+ (FMA3, 2013) plus `-march=native` or `-mfma`
            # (gcc/gfortran default `-march=x86-64` does NOT enable FMA
            # on x86). Pure Python's `*` and `-` never use FMA, so the
            # cross product is exactly zero.
            #
            # Substitute the mathematical limit explicitly. As areap_sum → 0+
            # with pes_sum = 0:
            #   first  norm: area_bottom = areap_sum / (2·areap_sum) = 0.5,  pes = 0
            #   second norm: area_bottom = areap_sum / areap_sum = 1.0,      pes = 0
            sum_a_is_zero = (sum_a == 0)
            if sum_a_is_zero:
                area_bottom = 0.5
                pes[i, :] = 0
            else:
                area_bottom /= sum_a
                pes[i, :] /= sum_a
            for k in range(cls.num_species_in_q3d):
                for kk in range(1, cls.max_z_layers-1):  # FORTRAN: do kk=2,max_z_layers-1 → 0-based 1..max_z_layers-2
                    hq3d[i, kk, k] += area_bottom * (cls.q3d[i, kk-1, k] - cls.q3d[i, kk, k])
                    hq3d[i, kk, k] += area_bottom * (cls.q3d[i, kk+1, k] - cls.q3d[i, kk, k])
                    for j in range(cls.nv_max):
                        if cls.neigh[i, j] != -1:
                            ii = cls.neigh[i, j]
                            if ii == cls.num_all_cells:
                                hq3d[i, kk, k] -= 0.44 * pes[i, j] * cls.q3d[i, kk, k]
                            else:
                                hq3d[i, kk, k] += pes[i, j] * (cls.q3d[ii, kk, k] - cls.q3d[i, kk, k])
                top = cls.max_z_layers - 1  # 1-based max_z_layers → 0-based max_z_layers-1
                hq3d[i, top, k] = -0.44 * area_bottom * cls.q3d[i, top, k]
                hq3d[i, top, k] += area_bottom * (cls.q3d[i, top-1, k] - cls.q3d[i, top, k])
                for j in range(cls.nv_max):
                    if cls.neigh[i, j] != -1:
                        ii = cls.neigh[i, j]
                        if ii == cls.num_all_cells:
                            hq3d[i, top, k] -= 0.44 * pes[i, j] * cls.q3d[i, top, k] # sink
                        else:
                            hq3d[i, top, k] += pes[i, j] * (cls.q3d[ii, top, k] - cls.q3d[i, top, k])
            if sum_a_is_zero:
                area_bottom = 1.0
                # pes stays 0
            else:
                pes[i, :] *= sum_a
                area_bottom *= sum_a    # FORTRAN 13.f90:451 — restore absolute, then subtract
                sum_a -= area_bottom
                pes[i, :] /= sum_a
                area_bottom /= sum_a    # FORTRAN 13.f90:452 — re-normalise to new sum
            for k in range(cls.num_species_in_q3d): # ATTENTION
                hq3d[i, 0, k] = area_bottom * (cls.q3d[i, 1, k] - cls.q3d[i, 0, k])
                for j in range(cls.nv_max):
                    if cls.neigh[i, j] != -1:
                        ii = cls.neigh[i, j]
                        if ii == cls.num_all_cells:
                            hq3d[i, 0, k] -= 0.44 * pes[i, j] * cls.q3d[i, 0, k]
                        else:
                            hq3d[i, 0, k] += pes[i, j] * (cls.q3d[ii, 0, k] - cls.q3d[i, 0, k])

        # Activator diffusion
        cls.q3d[:, :, 0] += cls.delta * cls.Da * hq3d[:, :, 0]

        # Inhibitor diffusion
        cls.q3d[:, :, 1] += cls.delta * cls.Di * hq3d[:, :, 1]

        # FGF diffusion
        cls.q3d[:, :, 2] += cls.delta * cls.Ds * hq3d[:, :, 2]

        # REACTION
        hq3d[:] = 0
        for i in range(cls.num_active_cells):
            if cls.q3d[i, 0, 0] > 1 and i >= cls.first_border_cell:
                cls.knots[i] = 1

            a = cls.Act * cls.q3d[i, 0, 0]  # Act = Activator (activator auto-activation)
            if a < 0: a = 0
            hq3d[i, 0, 0] = a / (1 + cls.Inh * cls.q3d[i, 0, 1]) - cls.Deg * cls.q3d[i, 0, 0]   # Eq. (14) sans diffusion
            if cls.diff_state[i] > cls.Int: # Int (initial inhibitor threshold)
                hq3d[i, 0, 1] = hq3d[i, 0, 0] * cls.diff_state[i] - cls.Deg * cls.q3d[i, 0, 1]  # Eq. (17). NOTE: DiffState <= 1.0
            else:
                if cls.knots[i] == 1:
                    hq3d[i, 0, 1] = cls.q3d[i, 0, 0] - cls.Deg * cls.q3d[i, 0, 1]   # Eq. (17)
            if cls.diff_state[i] > cls.Set:     # Set (growth factor threshold)
                a = cls.Sec * cls.diff_state[i] - cls.Deg * cls.q3d[i, 0, 2]    # Eq. (18). NOTE: DiffState <= 1.0
                if a < 0: a = 0
                hq3d[i, 0, 2] = a
            elif cls.knots[i] == 1:
                a = cls.Sec - cls.Deg * cls.q3d[i, 0, 2]  # Eq. (18). Sec = Sec (growth factor secretion rate)
                if a < 0: a = 0
                hq3d[i, 0, 2] = a

        max_abs = abs(hq3d[:, 0, 0:2]).max()
        if max_abs > 1e100:
            sys.exit()
        cls.q3d[:, 0, 0:3] += cls.delta * hq3d[:, 0, 0:3]   # explicit Euler method

        cls.q3d = cls.q3d.clip(min=0)

    
    @classmethod
    def apply_differentiation(cls):
        for i in range(cls.num_active_cells):
            cls.diff_state[i] += cls.Dff * cls.q3d[i, 0, 2]     # paper Eq. 6
            if cls.diff_state[i] > 1: cls.diff_state[i] = 1
    
    
    @classmethod
    def ep_growth_border_force(cls):
        uux: float
        uuy: float
        uuz: float
        ua: float
        ub: float
        uc: float
        uaa: float
        ubb: float
        uuux: float
        uuuy: float
        duux: float
        duuy: float

        cls.force_snapshot[:] = 0
        cls.forces[:] = 0
        for i in range(cls.first_border_cell, cls.num_active_cells):
            if cls.knots[i] == 1: continue
            ua = cls.positions[i, 0]
            ub = cls.positions[i, 1]
            uc = cls.positions[i, 2]
            aa = bb = cc = 0.0
            for j in range(cls.nv_max):
                k = cls.neigh[i, j]
                if k == -1 or k >= cls.num_active_cells: continue
                b = uc - cls.positions[k, 2]
                if b < -1e-4:
                    uux = ua-cls.positions[k, 0]
                    uuy = ub-cls.positions[k, 1]
                    uuz = uc-cls.positions[k, 2]
                    d = 1 / vector.magnitude(uux, uuy, uuz)
                    aa -= uux * d
                    bb -= uuy * d
                    cc -= uuz * d
            d = vector.magnitude(aa, bb, cc)
            if d > 0:
                d = cls.Egr / d
                a = 1 - cls.diff_state[i]
                if a < 0: a = 0
                d *= a
                cls.forces[i, 0] = aa * d
                cls.forces[i, 1] = bb * d
                cls.forces[i, 2] = cc * d

        for i in range(cls.first_border_cell):
            aa = bb = 0.0
            a = -0.3
            b = 0.0
            c = 0.0
            ua = cls.positions[i, 0]
            ub = cls.positions[i, 1]
            for j in range(cls.nv_max):
                k = cls.neigh[i, j]
                if k < 0 or k >= cls.num_active_cells: continue
                if k >= cls.first_border_cell:
                    uux = ua - cls.positions[k, 0]
                    uuy = ub - cls.positions[k, 1]
                    d = vector.magnitude(uux, uuy, 0)
                    if d > 0:
                        c = math.acos(uux / d)
                        if uuy < 0: c = 2 * math.pi - c
                else:
                    uux = ua - cls.positions[k, 0]
                    uuy = ub - cls.positions[k, 1]
                    d = vector.magnitude(uux, uuy, 0)
                    if d > 0:
                        if a == -0.3:
                            a = math.acos(uux / d)
                            if uuy < 0: a = 2 * math.pi - a
                            dd = 1 / d
                            uuux = -uuy * dd
                            uuuy = uux * dd
                            uaa = math.acos(uuux)
                            if uuuy < 0: uaa = 2 * math.pi - uaa
                        else:
                            b = math.acos(uux / d)
                            if uuy < 0: b = 2 * math.pi - b
                            dd = 1 / d
                            duux = -uuy * dd
                            duuy = uux * dd
                            ubb = math.acos(duux)
                            if duuy < 0: ubb = 2 * math.pi - ubb

            if a < b: a, b = b, a
            if c < a and c > b:
                if uaa < a and uaa > b:
                    # is on the inside side and then we have to invert it
                    uuux *= -1
                    uuuy *= -1
                if ubb < a and ubb > b:
                    duux *= -1
                    duuy *= -1
            else:
                if uaa > a or uaa < b:
                    # is on the inside side and then we have to invert it
                    uuux *= -1
                    uuuy *= -1
                if ubb > a or ubb < b:
                    duux *= -1
                    duuy *= -1
            aa = -uuux-duux
            bb = -uuuy-duuy

            # now let's see if it's outward from the tooth to the shabby
            a = ua + aa
            b = ub + bb
            c = ua - aa
            d = ub - bb
            dd = vector.magnitude(a, b, 0)
            ddd = vector.magnitude(c, d, 0)
            if ddd > dd:
                aa *= -1
                bb *= -1

            # We also have the downward traction due to the adhesion to the mesenchyme.
            d = vector.magnitude(aa, bb, 0)
            if d > 0:
                d = (d + cls.Mgr * cls.q3d[i, 0, 2] + cls.umgr) / d     # Eq. (12) + basal Mgr
                aa *= d     # Eq. 10
                bb *= d     # Eq. 11
            cc = cls.Dgr
            d = vector.magnitude(aa, bb, cc)
            if d > 0:
                d = cls.Egr / d
                a = 1 - cls.diff_state[i]
                if a < 0: a = 0     # DiffState(i) = differentiation state
                d *= a
                cls.forces[i, 0] = aa * d
                cls.forces[i, 1] = bb * d
                cls.forces[i, 2] = cc * d

        cls.force_snapshot = cls.forces.copy()


    @classmethod
    def boy_force(cls):
        ax: float
        ay: float

        # does it push perpendicularly? This is equation 13. 

        for i in range(cls.num_active_cells):
            ax = cls.forces[i, 0]
            ay = cls.forces[i, 1]
            d = vector.magnitude(ax, ay, 0)
            if d != 0:
                c = cls.forces[i, 2]
                a = vector.magnitude(ax, ay, c)
                a = -c / a
                ax *= a
                ay *= a
                dd = vector.magnitude(ax, ay, d)
                dd = cls.Boy * cls.q3d[i, 0, 2] / dd    # I believe this is Eq. 13. Boy corresponds to k_Boy
                if dd > 0:                              # q3d(i,0,2) = [Sec], DiffState(i) = d_i
                    a = 1 - cls.diff_state[i]           # this is the differentiation gate (DiffState(i) is DIFF state variable)
                    if a < 0: a = 0
                    ax *= dd * a
                    ay *= dd * a
                    d *= dd * a
                    cls.forces[i, 0] -= ax
                    cls.forces[i, 1] -= ay
                    cls.forces[i, 2] -= d
        
    
    @classmethod
    def repulse_neighbour(cls):
        persu = np.zeros((cls.nv_max, 3))

        # finite element roll
        for i in range(cls.num_active_cells):
            ua = cls.positions[i, 0]
            ub = cls.positions[i, 1]
            uc = cls.positions[i, 2]
            persu[:] = 0
            for j in range(cls.nv_max):
                k = cls.neigh[i, j]
                if k >= 0 and k < cls.num_active_cells:
                    ux = cls.positions[k, 0] - ua
                    uy = cls.positions[k, 1] - ub
                    uz = cls.positions[k, 2] - uc   # eq1: r_ij = p_j - p_i
                    if abs(ux) < 1e-15: ux = 0
                    if abs(uy) < 1e-15: uy = 0
                    if abs(uz) < 1e-15: uz = 0
                    dr = vector.magnitude(ux, uy, uz)   # eq1: d_r = ||p_j - p_i||
                    rd = cls.border[i, j, 4]
                    if dr < 1e-8: dr = 0
                    if rd < 1e-8: rd = 0
                    if (cls.knots[i] == 1 and cls.knots[k] == 1) or dr < rd:
                        d = dr - rd     # eq1: d= (||p_j - p_i|| -||p_ij^0||) (||p_j - p_i||)
                        dr = d / dr     # introduce normalisation
                        persu[j, 0] = ux * dr   # persu = instantenaous "force" contribution vector. will be updated & stored as Force
                        persu[j, 1] = uy * dr   # normalized eq1: d= (||p_j - p_i|| -||p_ij^0||) (||p_j - p_i||)
                        persu[j, 2] = uz * dr
                    elif i >= cls.first_border_cell:    # the cells from "0" to "first_border_cell-1" are interior epithelial cells, so this is saying "if i IS a border cell!"
                        persu[j, 0] = ux * cls.Adh
                        persu[j, 1] = uy * cls.Adh
                        persu[j, 2] = uz * cls.Adh

            # fast version without sorting (possible biases for floats)
            c = cls.Rep
            if c > 1: c = 1
            a = sum(persu[j, 0] for j in range(cls.nv_max))    # x direction forces
            cls.forces[i, 0] += a * c   # eq1 = kRep * (||p_j - p_i|| - ||p_o||)(p_j - p_i)
            a = sum(persu[j, 1] for j in range(cls.nv_max))    # y direction forces
            cls.forces[i, 1] += a * c
            a = sum(persu[j, 2] for j in range(cls.nv_max))    # z direction forces
            cls.forces[i, 2] += a * c
    

    @classmethod
    def repel_non_neigh(cls):
        conta = 0
        espai = 20
        persu = np.zeros((espai, 3))

        def repel_non_neigh_a(i: int, ii: int) -> bool:
            for j in range(cls.nv_max):
                if cls.neigh[i, j] == ii:
                    return True
            return False

        for i in range(cls.num_active_cells):
            ua = cls.positions[i, 0]
            ub = cls.positions[i, 1]
            uc = cls.positions[i, 2]
            persu[:] = 0
            conta = 0
            for ii in range(cls.num_active_cells):
                if i == ii: continue
                if repel_non_neigh_a(i, ii): continue
                ux = cls.positions[ii, 0] - ua
                if ux > 1.4: continue
                uy = cls.positions[ii, 1] - ub
                if uy > 1.4: continue
                uz = cls.positions[ii, 2] - uc
                if uz > 1.4: continue
                if abs(ux) < 1e-15: ux = 0
                if abs(uy) < 1e-15: uy = 0
                if abs(uz) < 1e-15: uz = 0
                d = vector.magnitude(ux, uy, uz)
                if d < 1.4:
                    if conta >= espai:    # FORTRAN 13.f90:747 — grow before writing
                        espai += 20
                        persu = np.pad(persu, ((0, 20), (0, 0)))
                    dd = 1 / (d + 1)**8
                    d = dd / d
                    d = np.floor(d*1e8)*1e-8
                    persu[conta, 0] = -ux * d
                    persu[conta, 1] = -uy * d
                    persu[conta, 2] = -uz * d
                    conta += 1

            # quick unsorted version (possible biases for floats)
            a = sum(persu[j, 0] for j in range(espai))
            cls.forces[i, 0] += a * cls.Rep
            a = sum(persu[j, 1] for j in range(espai))
            cls.forces[i, 1] += a * cls.Rep
            a = sum(persu[j, 2] for j in range(espai))
            cls.forces[i, 2] += a * cls.Rep
    

    @classmethod
    def apply_border_bias(cls):
        for i in range(cls.first_border_cell):
            if cls.positions[i, 1] < 0:
                if cls.q3d[i, 0, 0] < cls.Lbi:  # Won't reset activator concentration to value of bias; to make Ina work in borders.
                    cls.q3d[i, 0, 0] = cls.Lbi
            elif cls.positions[i, 1] > 0:
                if cls.q3d[i, 0, 0] < cls.Bbi:  # Same as above
                    cls.q3d[i, 0, 0] = cls.Bbi


    @classmethod
    def initact(cls):   # Adds initial activator concentration in each cell. Defined by Ina parameter.
        for i in range(cls.num_active_cells):
            cls.q3d[i, 0, 0] = cls.ina
    

    @classmethod
    def apply_nuclear_traction(cls):    # creates bias
        positions_after_traction = cls.positions.copy()
        for i in range(cls.first_border_cell, cls.num_active_cells):
            if cls.diff_state[i] == 1: continue
            a = b = c = 0.0
            n = 0
            for j in range(cls.nv_max):
                k = cls.neigh[i, j]
                if k != -1 and k < cls.num_active_cells:
                    a += cls.positions[k, 0]
                    b += cls.positions[k, 1]
                    c += cls.positions[k, 2]
                    n += 1
            a /= n
            b /= n
            c /= n
            a -= cls.positions[i, 0]
            b -= cls.positions[i, 1]
            c -= cls.positions[i, 2]
            positions_after_traction[i, 0] += cls.delta * cls.Ntr * a
            positions_after_traction[i, 1] += cls.delta * cls.Ntr * b
            if cls.knots[i] == 0:
                a = 1 - cls.diff_state[i]
                if a < 0: a = 0
                positions_after_traction[i, 2] += cls.delta * cls.Ntr * c * a

        # for the margins
        for i in range(cls.first_border_cell):
            if cls.diff_state[i] == 1: continue
            a = b = c = 0.0
            n = 0
            for j in range(cls.nv_max):
                k = cls.neigh[i, j]
                if k >= 0 and k < cls.first_border_cell and k < cls.num_active_cells:
                    a += cls.positions[k, 0]
                    b += cls.positions[k, 1]
                    c += cls.positions[k, 2]
                    n += 1
            a /= n
            b /= n
            c /= n
            a -= cls.positions[i, 0]
            b -= cls.positions[i, 1]
            c -= cls.positions[i, 2]
            positions_after_traction[i, 0] += cls.delta * cls.Ntr * a
            positions_after_traction[i, 1] += cls.delta * cls.Ntr * b
            if cls.knots[i] == 0:
                a = 1 - cls.diff_state[i]
                if a < 0: a = 0
                positions_after_traction[i, 2] += cls.delta * cls.Ntr * c * a
        cls.positions = positions_after_traction
    

    @classmethod
    def update_cell_position(cls):
        # we determine the extremes
        # FORTRAN 13.f90:863-872 — the forces[:,2] *= Bgr is INSIDE each
        # x-sign branch. Cells with x exactly == 0 (initial y-axis cells)
        # don't receive the Bgr boost. Applying Bgr unconditionally when
        # |y| < Bwi over-pushes y-axis cells in z by Bgr×, breaking the
        # equilibrium and causing exponential over-division downstream.
        for i in range(cls.first_border_cell):
            if abs(cls.positions[i, 1]) < cls.Bwi:
                if cls.positions[i, 0] > 0:
                    cls.forces[i, 0] *= cls.Pbi
                    cls.forces[i, 2] *= cls.Bgr
                elif cls.positions[i, 0] < 0:
                    cls.forces[i, 0] *= cls.Abi
                    cls.forces[i, 2] *= cls.Bgr

        for i in range(cls.num_active_cells):
            if cls.forces[i, 2] < 0 or cls.knots[i] == 1:
                cls.forces[i, 2] = 0    # it is due to the pressure of the stelate

        for i in range(cls.num_active_cells):
            cls.positions[i, :] += cls.delta * cls.forces[i, :]     # Eq3 kinda! This is the final position update (in the case where Swi = 0)

    
    @classmethod
    def add_cell(cls):
        new_cell_pairs = np.zeros((cls.num_active_cells * cls.nv_max, 2), dtype=np.int32)
        new_cell_is_external = np.zeros((cls.num_active_cells * cls.nv_max), dtype=bool)
        pillats = np.zeros((cls.nv_max), dtype=np.int32)
        num_new_cells = 0
        # FORTRAN 13.f90:1045-1054, 1072-1081: iiii is preserved across the
        # neighbour search when no new neighbour is found (only kkk is reset,
        # not iiii). FORTRAN's local-variable lifetime carries the previous
        # iiii forward; the topology walk relies on this to get stuck and
        # panic on degenerate cells (e.g. cell with only one neighbour),
        # which is the safety mechanism that prevents bogus cell additions.
        # Python must explicitly preserve iiii across both pair iterations
        # and within-loop searches. A naive `iiii = kkk = 0` reset jumps
        # the walk to cell index 0 (a real cell in 0-based) instead of
        # staying stuck, so the walk continues incorrectly and adds cells
        # FORTRAN refuses to add.
        iiii = 0

        # first we identify and name the new nodes and rescale the mesh matrix and see
        for i in range(cls.num_active_cells):
            for j in range(cls.nv_max):
                k = cls.neigh[i, j]
                ua = cls.positions[i, 0]
                ub = cls.positions[i, 1]
                uc = cls.positions[i, 2]
                if k != -1 and k > i and k < cls.num_active_cells:
                    ux = cls.positions[k, 0] - ua
                    uy = cls.positions[k, 1] - ub
                    uz = cls.positions[k, 2] - uc
                    a = vector.magnitude(ux, uy, uz)
                    a = round(a, 9)
                    if a > cls.dmax:    # we add new node
                        new_cell_pairs[num_new_cells, 0] = i
                        new_cell_pairs[num_new_cells, 1] = k
                        if i < cls.first_border_cell and k < cls.first_border_cell:
                            new_cell_is_external[num_new_cells] = True
                        num_new_cells += 1

        if num_new_cells > 0:
            for i in range(cls.num_all_cells):
                for j in range(cls.nv_max):
                    if cls.neigh[i, j] == -1:
                        for jj in range(j, cls.nv_max-1):
                            cls.neigh[i, jj] = cls.neigh[i, jj+1]

            new_num_active_cells = cls.num_all_cells + num_new_cells

            temp_positions =    np.zeros((new_num_active_cells, 3))
            temp_neigh =        np.zeros((new_num_active_cells, cls.nv_max), dtype=np.int32)
            temp_num_neigh =    np.zeros((new_num_active_cells), dtype=np.int32)
            c_diff_state =      np.zeros((new_num_active_cells))
            c_q3d =             np.zeros((new_num_active_cells, cls.max_z_layers, cls.num_species_in_q3d))
            temp_new_neigh =    np.zeros((num_new_cells, cls.nv_max), dtype=np.int32)
            c_knots =           np.zeros((new_num_active_cells), dtype=np.int32)
            temp_border =       np.zeros((new_num_active_cells, cls.nv_max, 8))
            temp_neigh[:] = -1
            temp_new_neigh[:] = -1

            for i in range(cls.num_active_cells + num_new_cells, new_num_active_cells):
                temp_positions[i, :] =  cls.positions[i - num_new_cells, :]
                temp_neigh[i, :] =      cls.neigh[i - num_new_cells, :]   # FORTRAN 13.f90:963
                temp_num_neigh[i] =     cls.num_neigh[i - num_new_cells]
                c_diff_state[i] =       cls.diff_state[i - num_new_cells]
                c_q3d[i, :, :] =        cls.q3d[i - num_new_cells]
                c_knots[i] =            cls.knots[i - num_new_cells]
                temp_border[i, :, 3:8] = cls.border[i - num_new_cells, :, 3:8]

            for i in range(cls.num_active_cells):
                temp_positions[i, :] =  cls.positions[i, :]
                temp_neigh[i, :] =      cls.neigh[i, :]
                temp_num_neigh[i] =     cls.num_neigh[i]
                c_diff_state[i] =       cls.diff_state[i]
                c_q3d[i, :, :] =        cls.q3d[i, :, :]
                c_knots[i] =            cls.knots[i]
                temp_border[i, :, :] =  cls.border[i, :, :]

            for i in range(cls.num_active_cells + num_new_cells):
                for j in range(cls.nv_max):
                    if temp_neigh[i, j] >= cls.num_active_cells:
                        temp_neigh[i, j] = new_num_active_cells

            for i in range(num_new_cells):
                ii = new_cell_pairs[i, 0]
                kk = new_cell_pairs[i, 1]
                jj = cls.num_active_cells + i
                temp_neigh[jj, :] = -1
                temp_neigh[jj, 0] = ii
                temp_neigh[jj, 1] = kk

                a = temp_positions[ii, 0] + temp_positions[kk, 0]
                b = temp_positions[ii, 1] + temp_positions[kk, 1]
                d = vector.magnitude(a, b, 0)
                a /= d
                b /= d
                d /= 2
                d = round(d, 10)
                a *= d
                b *= d
                temp_positions[jj, 0] = a
                temp_positions[jj, 1] = b

                temp_positions[jj, 2] = 0.5 * (cls.positions[ii, 2] + cls.positions[kk, 2])
                c_q3d[jj, :, :] =       0.5 * (cls.q3d[ii, :, :] + c_q3d[kk, :, :])
                c_diff_state[jj] =      0.5 * (cls.diff_state[ii] + cls.diff_state[kk])
                # The 0.333 instead of 0.25 is a trick because before we multiplied the parents by 3/4

            for i in range(num_new_cells):
                # FORTRAN 13.f90:999-1004: replace each parent's reference to
                # the other parent with the new cell, so the parents are now
                # separated by the newly-inserted cell instead of directly
                # adjacent.
                ii = new_cell_pairs[i, 0]
                kk = new_cell_pairs[i, 1]
                jj = cls.num_active_cells + i
                for j in range(cls.nv_max):
                    if temp_neigh[ii, j] == kk:
                        temp_neigh[ii, j] = jj
                        break
                for j in range(cls.nv_max):
                    if temp_neigh[kk, j] == ii:
                        temp_neigh[kk, j] = jj
                        break

            for i in range(num_new_cells):
                pillats[:] = -1
                ii = new_cell_pairs[i, 0]
                kk = new_cell_pairs[i, 1]
                jj = cls.num_active_cells + i
                # now I have to look at which parent is to follow it (it will be called ini) it is the one that has no external nodes towards j+1
                jjj = 0
                for j in range(cls.nv_max):
                    if cls.neigh[ii, j] == kk:
                        jjj = j
                        break
                ini = fi = kkk = 0
                for jjjj in range(jjj+1, cls.nv_max):
                    if cls.neigh[ii, jjjj] >= 0:
                        kkk = 1
                        if cls.neigh[ii, jjjj] < cls.num_active_cells:
                            ini = ii
                            fi = kk
                        else:
                            ini = kk
                            fi = ii
                        break
                if kkk == 0:
                    for jjjj in range(jjj):
                        if cls.neigh[ii, jjjj] >= 0:
                            if cls.neigh[ii, jjjj] < cls.num_active_cells:
                                ini = ii
                                fi = kk
                            else:
                                ini = kk
                                fi = ii
                            break
                iii = ini
                cj = 0
                pillats[cj] = iii
                # now we look for the j in which iii has jj (the cell whose neigh we are looking for)
                for j in range(cls.nv_max):
                    if temp_neigh[iii, j] == jj:
                        jjj = j
                        break
                # side j+1(right); we follow the neighbor towards the j+1 side of iii
                kkk = 0  # iiii preserved (see comment in add_cell head)
                for j in range(jjj+1, cls.nv_max):
                    jji = temp_neigh[iii, j]
                    if jji != -1 and jji < cls.num_active_cells + num_new_cells:
                        iiii = jji
                        kkk = 1
                        break
                if kkk == 0:    # I couldn't find the neighbor and I need to go over it again.
                    for j in range(jjj):
                        jji = temp_neigh[iii, j]
                        if jji != -1 and jji < cls.num_active_cells + num_new_cells:
                            iiii = jji
                            break

                cj += 1
                if cj >= cls.nv_max:
                    return  # FORTRAN: panic=1; return
                pillats[cj] = iiii
                for j in range(cls.nv_max):
                    if temp_neigh[iiii, j] == iii:
                        jjjj = j
                        break

                # False if goto 88 executed
                looping = True
                cycling = False
                while looping:
                    iii = iiii
                    jjj = jjjj

                    kkk = 0  # iiii preserved (see comment in add_cell head)
                    for j in range(jjj+1, cls.nv_max):
                        jji = temp_neigh[iii, j]
                        if jji != -1 and jji < cls.num_active_cells + num_new_cells:
                            iiii = jji
                            kkk = 1
                            break
                    if kkk == 0:    # I couldn't find the neighbor and I need to go over it again.
                        for j in range(jjj):
                            jji = temp_neigh[iii, j]
                            if jji != -1 and jji < cls.num_active_cells + num_new_cells:
                                iiii = jji
                                break

                    cj += 1
                    if cj >= cls.nv_max:
                        return  # FORTRAN: panic=1; return
                    pillats[cj] = iiii
                    for j in range(cls.nv_max):
                        if temp_neigh[iiii, j] == iii:
                            jjjj = j
                            break
                    if iiii == fi:  # equinox
                        kkk = sjj = 0
                        for kkkk in range(jjjj+1, cls.nv_max):
                            if temp_neigh[iiii, kkkk] != -1 and kkk == 1:
                                kkk = 2
                                if temp_neigh[iiii, kkkk] >= cls.num_active_cells + num_new_cells:
                                    iiii = ini
                                    cj += 1
                                else:
                                    sjj = kkkk
                                break
                            if temp_neigh[iiii, kkkk] != -1 and kkk == 0:
                                kkk = 1
                                sjj = kkkk
                        if kkk < 2:
                            for kkkk in range(jjjj):
                                if temp_neigh[iiii, kkkk] != -1 and kkk == 1:
                                    kkk = 2
                                    if temp_neigh[iiii, kkkk] >= cls.num_active_cells + num_new_cells:
                                        iiii = ini
                                        cj += 1
                                    else:
                                        sjj = kkkk
                                    break
                                if temp_neigh[iiii, kkkk] != -1 and kkk == 0:
                                    kkk = 1
                                    sjj = kkkk
                        jjjj = sjj - 1

                    # we have gone all the way around
                    if iiii == ini:
                        if cj >= cls.nv_max:
                            return  # FORTRAN: panic=1; return
                        c_pillats = np.full((cls.nv_max), -1, dtype=np.int32)
                        pillats[cj] = -1
                        cj -= 1
                        for jjj in range(cj + 1):
                            c_pillats[cj - jjj] = pillats[jjj]
                        pillats = c_pillats
                        # now let's see what nodes can actually be
                        jjj = 0
                        for kkk in range(cj + 1):
                            kkkk = pillats[kkk]
                            if kkkk >= cls.num_active_cells and kkkk < cls.num_active_cells + num_new_cells:
                                jjj += 1
                        if jjj == 0:
                            temp_new_neigh[i, :] = pillats.copy()   # we don't have new nodes on the sides then the crossing is impossible
                            for j in range(cls.nv_max):
                                k = temp_new_neigh[i, j]
                                if k != -1:
                                    ii = new_cell_pairs[i, 0]
                                    iiii = new_cell_pairs[i, 1]
                                    if k != ii and k != iiii and k < cls.num_active_cells + num_new_cells: # It's one of those that I have to connect with.
                                        kkkk = 0
                                        for kk in range(cls.nv_max):
                                            for kkk in range(cls.nv_max):
                                                if temp_neigh[k, kk] != -1 and temp_neigh[k, kk] == pillats[kkk]:
                                                    if kkkk == 1:
                                                        ji = kk
                                                    else:
                                                        ij = kk
                                                    kkkk += 1
                                                    break
                                            if kkkk == 2: break
                                        if ji - ij == 1:
                                            for kk in range(cls.nv_max-1, ji, -1):
                                                temp_neigh[k, kk] = temp_neigh[k, kk-1]
                                                temp_border[k, kk, 3:8] = temp_border[k, kk-1, 3:8]
                                            temp_neigh[k, ji] = jj
                                        else:
                                            temp_neigh[k, ji+1] = jj
                        else:   # the txungu will be in sinomes we have a new one because then there will only be 3 neigh
                            temp_new_neigh[i, :] = -1
                            temp_new_neigh[i, 0] = ini
                            kkkk = 0
                            if cj >= cls.nv_max:
                                return  # FORTRAN: panic=1; return
                            for kkk in range(cj + 1):
                                jjj = pillats[kkk]
                                # KNOWN ISSUE — `13.f90` renaming inconsistency,
                                # preserved verbatim because it is load-bearing
                                # for `13.f90`'s goldens.
                                #
                                # `13.f90:1156-1163` assigns `jjj=pillats(kkk)`
                                # but tests `jjjj` (a loop-invariant value set
                                # at `13.f90:1096` to `sjj-1` — a *slot index*
                                # in temp_neigh, runs 1..nv_max). The two checks
                                # `if (jjjj == fi)` and `if (jjjj > num_active_cells)`
                                # compare a slot index to a cell index —
                                # semantically nonsense. The loop body never
                                # fires; the loop is dead code in `13.f90`.
                                #
                                # Provenance: this inconsistency is *not* in
                                # the upstream Catalan FORTRAN.
                                # `humppa_translate.f90:1336-1345` and
                                # `tgrohens/coreop2d.f90:1276-1285` both write
                                # `jjjj=pillats(kkk)` and test `jjjj` —
                                # internally consistent, the loop runs. The
                                # Catalan→English renaming changed the
                                # assignment to `jjj` but missed the
                                # conditional tests. See the module-level
                                # docstring for the full chain of accidents
                                # this triggers.
                                #
                                # Trace evidence (PR #15) from gfortran -O2 on
                                # the seal example: the txungu else-branch
                                # fires 506 times across 5 save blocks, and
                                # EVERY firing produces
                                # temp_new_neigh[i, :] = [ini, 0, 0, ...].
                                # The conditional checks never fire because
                                # jjjj is always 1..5 while fi is ~20+ and
                                # num_active_cells is 37+.
                                #
                                # The dead-loop output (single-neighbour cells)
                                # is what the rest of `13.f90` expects:
                                # calculate_margins's label-77 fall-through
                                # fills empty border slots with the
                                # last-computed value; apply_diffusion's
                                # fallback then computes a cross product of
                                # bit-equal vectors that FMA turns into an
                                # ulp-noise non-zero. Break this chain at any
                                # point — fix this inconsistency, fix
                                # label-77, fix the fall-through, disable FMA
                                # — and `13.f90` also produces a NaN cascade.
                                #
                                # Naively changing `jjjj` to `jjj` here (or
                                # equivalently renaming the test in `13.f90`
                                # to match) produces 629 cells/block (target
                                # 57) — confirms that the rest of `13.f90`
                                # has been tuned to the dead-loop output.
                                # See PR #12, PR #14, PR #15 for the full
                                # investigation chain. Path B should NOT
                                # reproduce this branch's behaviour; humppa
                                # does not.
                                if jjjj == fi:
                                    kkkk += 1
                                    temp_new_neigh[i, kkkk] = fi
                                elif jjjj >= cls.num_active_cells:
                                    kkkk += 1
                                    temp_new_neigh[i, kkkk] = jjjj

                        # now we need to add the connections to external nodes
                        ii = new_cell_pairs[i, 0]
                        kk = new_cell_pairs[i, 1]
                        kkk = 0
                        for j in range(cls.nv_max):
                            if temp_neigh[ii, j] >= cls.num_active_cells + num_new_cells:
                                kkk = 1
                                break
                        for j in range(cls.nv_max):
                            if temp_neigh[kk, j] >= cls.num_active_cells + num_new_cells:
                                kkk += 1
                                break
                        if kkk == 2:
                            ij = jjj = 0
                            for j in range(cls.nv_max):
                                if temp_new_neigh[i, j] == ii:
                                    ij = j
                                    break
                            for j in range(cls.nv_max):
                                if temp_new_neigh[i, j] == kk:
                                    jjj = j
                                    break
                            if ij > jjj:
                                ji = ij
                                ij = jjj
                            else:
                                ji = jjj
                            # we have to connect between ij and ji
                            if ji - ij == 1:
                                for kk in range(cls.nv_max-1, ji, -1):
                                    temp_new_neigh[i, kk] = temp_new_neigh[i, kk-1]
                                    temp_border[i, kk, 3:8] = temp_border[i, kk-1, 3:8]
                                temp_new_neigh[i, ji] = new_num_active_cells
                            else:
                                temp_new_neigh[i, ji+1] = new_num_active_cells # temp_border(i,ji+1,4:5)=0.
                        looping = False # cycle statement ends while loop
                    # goto 88 just continues while loop

            # now we need to add the external connections

            # now we replace
            temp_neigh[cls.num_active_cells:cls.num_active_cells+num_new_cells, :] = temp_new_neigh[0:num_new_cells, :]
            cls.positions = temp_positions

            # we calculate the new basal distances of the new cells
            for i in range(cls.num_active_cells, cls.num_active_cells + num_new_cells):
                ua = cls.positions[i, 0]
                ub = cls.positions[i, 1]
                uc = cls.positions[i, 2]
                for j in range(cls.nv_max):
                    ii = temp_neigh[i, j]
                    if ii >= 0 and ii < cls.num_active_cells + num_new_cells:
                        ux = cls.positions[ii, 0] - ua
                        uy = cls.positions[ii, 1] - ub
                        uz = cls.positions[ii, 2] - uc
                        if abs(ux) < 1e-13: ux = 0
                        if abs(uy) < 1e-13: uy = 0
                        if abs(uz) < 1e-13: uz = 0
                        d = vector.magnitude(ux, uy, 0)
                        temp_border[i, j, 4] = d
                        d = vector.magnitude(ux, uy, uz)
                        temp_border[i, j, 3] = d
                        temp_border[i, j, 5] = ux
                        temp_border[i, j, 6] = uy
                        temp_border[i, j, 7] = uz

            # we calculate the basal distances of the new connections between old and new cells
            for i in range(cls.num_active_cells):
                ua = cls.positions[i, 0]
                ub = cls.positions[i, 1]
                uc = cls.positions[i, 2]
                for j in range(cls.nv_max):
                    ii = temp_neigh[i, j]
                    if ii >= cls.num_active_cells and ii < cls.num_active_cells + num_new_cells:
                        ux = cls.positions[ii, 0] - ua
                        uy = cls.positions[ii, 1] - ub
                        uz = cls.positions[ii, 2] - uc
                        if abs(ux) < 1e-13: ux = 0
                        if abs(uy) < 1e-13: uy = 0
                        if abs(uz) < 1e-13: uz = 0
                        d = vector.magnitude(ux, uy, 0)
                        temp_border[i, j, 4] = d
                        d = vector.magnitude(ux, uy, uz)
                        temp_border[i, j, 3] = d
                        temp_border[i, j, 5] = ux
                        temp_border[i, j, 6] = uy
                        temp_border[i, j, 7] = uz

            cls.num_all_cells = new_num_active_cells
            prev_num_active_cells = cls.num_active_cells
            cls.num_active_cells += num_new_cells

            cls.forces = np.zeros((new_num_active_cells, 3))
            cls.force_snapshot = np.zeros((new_num_active_cells, 3))
            cls.px = np.zeros(new_num_active_cells)
            cls.py = np.zeros(new_num_active_cells)
            cls.pz = np.zeros(new_num_active_cells)

            cls.neigh = temp_neigh
            cls.num_neigh = temp_num_neigh
            cls.diff_state = c_diff_state
            cls.q3d = c_q3d
            cls.knots = c_knots
            cls.border = temp_border
            cls.forces[:] = 0

            # we remove the sentinel -1 values (compact neigh)
            for i in range(cls.num_all_cells):
                for j in range(cls.nv_max):
                    if cls.neigh[i, j] == -1:
                        for jj in range(j, cls.nv_max-1):
                            cls.neigh[i, jj] = cls.neigh[i, jj+1]
            for i in range(cls.num_all_cells):
                ii = 0
                for j in range(cls.nv_max):
                    if cls.neigh[i, j] >= 0:
                        ii += 1
                cls.num_neigh[i] = ii
        # END IF WE HAVE NEW CELLS

        for iii in range(num_new_cells):
            if new_cell_is_external[iii]: # we have a new external cell
                ii = prev_num_active_cells + iii
                if cls.first_border_cell == cls.centre:
                    cls.centre = ii
                
                snapshot_neigh = cls.neigh.copy()
                snapshot_num_neigh = cls.num_neigh.copy()
                post_division_temp_positions = cls.positions.copy()
                sc_diff_state = cls.diff_state.copy()
                sc_q3d = cls.q3d.copy()
                temp_border_of_external_cells = cls.border.copy()
                sc_knots = cls.knots.copy()

                cls.neigh[ii, :] = snapshot_neigh[cls.first_border_cell, :]
                cls.neigh[cls.first_border_cell, :] = snapshot_neigh[ii, :]
                cls.num_neigh[ii] = snapshot_num_neigh[cls.first_border_cell]
                cls.num_neigh[cls.first_border_cell] = snapshot_num_neigh[ii]
                cls.positions[ii, :] = post_division_temp_positions[cls.first_border_cell, :]
                cls.positions[cls.first_border_cell, :] = post_division_temp_positions[ii, :]
                cls.diff_state[ii] = sc_diff_state[cls.first_border_cell]
                cls.diff_state[cls.first_border_cell] = sc_diff_state[ii]
                cls.q3d[ii, :] = sc_q3d[cls.first_border_cell, :]
                cls.q3d[cls.first_border_cell, :] = sc_q3d[ii, :]
                cls.border[ii, :] = temp_border_of_external_cells[cls.first_border_cell, :]
                cls.border[cls.first_border_cell, :] = temp_border_of_external_cells[ii, :]
                cls.knots[ii] = sc_knots[cls.first_border_cell]
                cls.knots[cls.first_border_cell] = sc_knots[ii]
                snapshot_neigh = cls.neigh.copy()

                for i in range(cls.num_active_cells):
                    for j in range(cls.nv_max):
                        if snapshot_neigh[i, j] == ii:
                            cls.neigh[i, j] = cls.first_border_cell
                for i in range(cls.num_active_cells):
                    for j in range(cls.nv_max):
                        if snapshot_neigh[i, j] == cls.first_border_cell:
                            cls.neigh[i, j] = ii
                cls.first_border_cell += 1     # FORTRAN 13.f90:1384: expand border-cell block by 1


        if num_new_cells > 0:
            cls.update_border_cells()

    @classmethod
    def update_border_cells(cls):
        # makes the new cells that are on the margin between two extremes for the bias extreme
        new_border_cells = np.zeros(cls.first_border_cell, dtype=np.int32)
        num_new_border_cells = 0

        # now let's correct the ends
        # code inside loop labelled er
        def er(i: int):
            nonlocal num_new_border_cells
            for ii in range(cls.num_border_cells):
                iii = cls.border_cell_list[ii]
                if iii == i: return # cycle er
            kk = False
            for j in range(cls.nv_max): # err
                k = cls.neigh[i, j]
                if k < cls.first_border_cell:
                    for ii in range(cls.num_border_cells):
                        iii = cls.border_cell_list[ii]
                        if k == iii:
                            if kk:
                                new_border_cells[num_new_border_cells] = i
                                num_new_border_cells += 1
                                return # cycle era
                            else:
                                kk = True
                                break # cycle err
        for i in range(cls.first_border_cell):
            if i == 2 or i == 5: continue  # 0-based equivalents of 1-based 3 and 6 (ATTENTION NOT SCALABLE WITH Rad)
            er(i)

        if num_new_border_cells > 0:
            old_num_border_cells = cls.num_border_cells
            old_border_cell_list = cls.border_cell_list
            cls.num_border_cells += num_new_border_cells
            cls.border_cell_list = np.zeros(cls.num_border_cells, dtype=np.int32)
            cls.border_cell_list[0:cls.num_border_cells-num_new_border_cells] = old_border_cell_list[:]
            for i in range(num_new_border_cells):
                cls.border_cell_list[old_num_border_cells+i] = new_border_cells[i]

        new_border_cells[:] = -1
        num_new_border_cells = 0
        # code inside loop labelled era
        def era(i: int):
            nonlocal num_new_border_cells
            for ii in range(cls.n_map):
                iii = cls.m_map[ii]
                if iii == i: return # cycle era
            kk = False
            for j in range(cls.nv_max): # erra
                k = cls.neigh[i, j]
                if k < cls.first_border_cell:
                    for ii in range(cls.n_map):
                        iii = cls.m_map[ii]
                        if k == iii:
                            if kk:
                                new_border_cells[num_new_border_cells] = i
                                num_new_border_cells += 1
                                return # cycle era
                            else:
                                kk = True
                                break # cycle erra
        for i in range(cls.first_border_cell):
            if i == 2 or i == 5: continue  # 0-based equivalents of 1-based 3 and 6 (ATTENTION NOT SCALABLE WITH Rad)
            era(i)

        if num_new_border_cells > 0:
            old_border_cells = cls.n_map
            old_border_cell_list = cls.m_map
            cls.n_map += num_new_border_cells
            cls.m_map = np.zeros(cls.n_map, dtype=np.int32)
            cls.m_map[0:cls.n_map-num_new_border_cells] = old_border_cell_list
            for i in range(num_new_border_cells):
                cls.m_map[old_border_cells+i] = new_border_cells[i]
    

    @classmethod
    def iteration(cls, tbu: int):
        for _ in range(tbu):
            cls.apply_diffusion()
            cls.apply_border_bias()
            cls.apply_differentiation()
            cls.ep_growth_border_force()
            cls.boy_force()
            cls.repel_non_neigh()
            cls.repulse_neighbour()
            cls.apply_nuclear_traction()
            cls.update_cell_position()
            cls.add_cell()
            cls.calculate_margins()
            cls.temps += 1