---
title: 'Code review — C++ port of jernvall-toothmaker (Path B planning)'
author: Lyndon Drake (with Claude Code)
date: 2026-05-04
---

> **Editor's note (added after the review).** At review time the local
> mirror was missing three of the nine C++ files
> (`tooth_model_mechanics.cpp`, `file_io.cpp`, `main.cpp`) because the
> initial fetch used the wrong default-branch name (`main` rather than
> the repo's actual `master`). All nine files were re-fetched into
> `.tmp/jernvall-toothmaker/` after the review completed and now match
> the upstream. The 'cannot link / unverifiable' caveats below were
> accurate at review time. Spot-checks against the now-complete mirror:
> Bug 2 (Bgr branch placement) is handled correctly by the C++ —
> `tooth_model_mechanics.cpp:537-545` nests `positionDeltas[i][2] *=
> biasFactor` inside each x-sign branch, matching FORTRAN. The other
> review findings (Bug 1 'if suma > 0' skip semantic, Bug 3 partial
> cross-pair regression, off-by-one panic threshold, no validation
> harness, not vectorised, README glossary accurate) stand.

## Scope and material on disk

The local mirror at `/home/lyndon/repo/project/ToothmakerTranslation/.tmp/jernvall-toothmaker/` contains six source files, the README, and the FORTRAN reference. Three files named in `cpp_README.md` (and in the review brief) are absent from the mirror:

- `tooth_model_mechanics.cpp` (~475 lines, claimed)
- `file_io.cpp` (~663 lines, claimed)
- `main.cpp` (~104 lines, claimed)

A `find` over the entire `.tmp/` tree returns no matches for these names. The C++ port as supplied therefore cannot link: `tooth_model_iteration.cpp:24-29` calls `calculateGrowthPushing`, `calculateBuoyancy`, `checkNonNeighborRepulsion`, `calculateNeighborRepulsion`, `calculateNucleusTraction`, and `updatePositions`, all declared in the header but defined nowhere in the present sources. This is the dominant fact about the review and constrains every other answer.

## Test/validation status

There is no test infrastructure in the mirror. `cpp_README.md` lines 494-514 describe a `make test` target that runs three tests including a Fortran-vs-C++ cross-validation with `tolerance 0.001`, but no `Makefile`, `test/`, or `reference/` directory exists in `.tmp/jernvall-toothmaker/` or its parent. There is no `tests/` directory in the C++ port itself, no golden output file, no comment of the form 'matches FORTRAN at iteration N', and no fixture under `.tmp/` referenced from the C++ source.

The README's claim of cross-validation is unbacked by artefacts in the mirror. Either the test directory was sibling to the C++ port and not copied across when the mirror was made (plausible — README references `cd ..; make test`), or it was never produced. Either way, **the reviewer cannot verify that the C++ port has ever been compared bit-for-bit, or even within tolerance, against a FORTRAN run**, and the port cannot even be built locally to check. For Path B planning, this should be treated as 'no validation' until the missing pieces are sourced.

## The three Path-A bugs

### 2(a) `apply_diffusion` divide-by-zero when sum_a == 0

The C++ at `tooth_model_diffusion.cpp:96-101` and `:148-153` reads:

```cpp
if (suma > 0) {
    areaBelow /= suma;
    for (int j = 0; j < MAX_NEIGHBORS; j++) {
        rdWeight[i][j] /= suma;
    }
}
```

This is a **third behaviour, distinct from both the FORTRAN and our Python**. FORTRAN at `humppa_translate.f90:516` and `:550-551` divides unconditionally and survives only via FMA-induced ulp noise in the cross-product. Our Python's PR #16 fix substitutes the analytic limit (`area_bottom = 0.5` first norm, `1.0` second norm; weights = 0). The C++ guards by skipping the division entirely, leaving `rdWeight[i][j]` and `areaBelow` at their unnormalised raw-area values. For a degenerate cell with `sum_a == 0`, the C++ then enters the diffusion update loop with full unnormalised areas — propagating quantitatively different (and probably much larger in magnitude) deltas than either the FORTRAN or our limit-based Python. Whether this matters for the seal parameter set is empirical, but it is not the FORTRAN behaviour and not our Python behaviour.

### 2(b) `update_cell_position` Bgr placement

`updatePositions` is declared in `tooth_model.hpp:196` and called from `tooth_model_iteration.cpp:29`, but its definition lives in the missing `tooth_model_mechanics.cpp`. **The reviewer cannot check whether the C++ preserves the FORTRAN nesting `if (malla(i,1)>0) ... hmalla(i,3)=hmalla(i,3)*fac` (`humppa_translate.f90:1049-1054`) or unconditionally applies `fac` to all y-axis cells.** This is precisely the bug we just fixed in our Python. It is the most economically valuable piece of the review and it is unverifiable on what is on disk.

### 2(c) `add_cell` iiii preservation in topology walk

This one is checkable, and the C++ does **not** match FORTRAN here.

In `tooth_model_division.cpp:295`, `int iiii = -1;` is declared inside the `for (int i = 0; i < numNewCells; i++)` per-pair loop (line 220). Each new-cell pair iteration starts with `iiii = -1`. FORTRAN's `iiii` is declared at subroutine scope and persists across pairs within one call to `afegircel`, which our `coreop2d.py:1010` (`iiii = 0` outside the outer pair loop) explicitly preserves with a long comment.

The within-pair carry-forward (when both forward and wrap searches fail and `iiii` retains its prior value, FORTRAN `humppa_translate.f90:1241-1247`) **is** correctly preserved by the C++: lines 340-359 only assign `iiii` inside the if-found branches, and `pillats[cj] = iiii + 1` at line 367 records the carried-forward value. So one of the two iiii-preservation semantics is correct; the across-pair one is broken.

The same site has one further per-call divergence: at lines 317, 362, 429, the C++ panics on `cj >= MAX_NEIGHBORS` (i.e. `cj >= 30`). FORTRAN at `humppa_translate.f90:1248`, `:1257`, `:1291`, `:1299` panics on `cj > nvmax` (i.e. `cj > 30`, or `cj >= 31`). The C++ panics one step earlier than FORTRAN. For chains of length 30, FORTRAN would proceed; C++ would panic-and-return. The Path-A finding flagged the panic-and-return as the safety net that prevents bogus cell additions; tightening the threshold by one means more iterations panic where FORTRAN would not, dropping more cell-addition requests.

Net answer to question 2: of the three Path-A bugs, **(a) is handled differently from both FORTRAN and our fix, (b) is unverifiable on present files, (c) has one half right and one half wrong, plus a related off-by-one in the panic threshold**.

## Vectorisation status (question 3)

Not vectorised at all. No Eigen, no `std::valarray`, no SIMD pragmas, no `#pragma omp` (the README mentions OpenMP only as a future improvement, `cpp_README.md:397-401`). Storage is `std::vector<std::array<double, 3>> cellPositions`, `std::vector<std::vector<std::array<double, 8>>> cellMargins`, etc. — array-of-structs with one outer level per cell. Inner loops are explicit per-cell-per-neighbour scalar loops. The shape is essentially a 1:1 transliteration of the FORTRAN's `do i=1,ncels ; do j=1,nvmax` patterns. Re-translating to NumPy would mean rewriting nearly every loop from scratch as broadcasted operations on `(N_cells, MAX_NEIGHBORS)` arrays. The C++ files give you the *equations* in cleaner form than the FORTRAN, but they do not give you a vectorisable structure for free.

The one place the C++ has done structural rework is `getTriangles()` in `tooth_model_geometry.cpp:151-240`, which uses sorted neighbour lists and `std::set_intersection` to find triangles — this is decent C++, but it is mesh-extraction-only, not part of the simulation hot loop, and the algorithm is incompatible with the original FORTRAN `guardaveinsoff` (the README notes this at lines 258-264 and that triangles are written twice with both orientations to dodge winding-order recovery).

## Risky translation choices (question 4)

These are real concerns I would flag in a code review of this port even if it built:

1. **`placeCell()` first-empty-slot fallback** (`tooth_model_core.cpp:323-329`): the comment is explicit — 'For now, find first empty slot'. FORTRAN's `posar` (`humppa_translate.f90:218-223`) is called with explicit direction arguments `j` and `jj` per call; the C++ caller (`allocateAndInit`) does not set them, and `placeCell` improvises. The hexagonal slot ordering in `neighbors[i][j]` therefore does not correspond to a stable compass direction. Subsequent uses of `(j + 3) % 6` as 'opposite direction' (lines 330, 348) inherit this ambiguity. This may happen to produce a topologically correct grid for radius-3 starts because all cells are filled in the same order, but it is structurally fragile for any later code that assumes slot index = compass direction.

2. **`numBorderCells` semantic shift** (`tooth_model_core.cpp:228`, README:118): the C++ uses `(radius - 1) * 6`, dropping the FORTRAN's `+1`, claiming this is an intentional 0-based-indexing adjustment. The change is sound *if* every downstream consumer is also adjusted. The README says results are identical but with no test artefact to verify, this is an assertion not a measurement. Several boundary-condition checks in the diffusion code use `i >= numBorderCells` (e.g. `tooth_model_diffusion.cpp:191`) where one off-by-one in the count gives a 'border' cell that should be 'interior', or vice versa. Worth auditing carefully.

3. **`if (knotMarkers[i] > diffThresholdSet)`** (`tooth_model_diffusion.cpp:218`): integer 0/1 compared to a double parameter typically in [0, 1]. Faithful to FORTRAN `humppa_translate.f90:628` `if (knots(i)>ud)` so this is a faithful port of a FORTRAN bug. The inhibitor branch above (`tooth_model_diffusion.cpp:207`) also does `if (knotMarkers[i] == 1)` which has the integer interpretation. Both deserve a comment in any port forward; the C++ silently transliterates without comment.

4. **Markers in `markBorderCells()` hard-coded indices** (`tooth_model_division.cpp:749, :795`): `if (i == 2 || i == 5) continue;`. The comment 'Skip certain indices (not scalable with radius - keeping original logic)' admits this is wrong for any radius other than 3. Faithful to FORTRAN, but it's a known-broken corner that any radius-parametric work has to deal with.

5. **`std::round(1000000 * position)` for cell-position equality** (`tooth_model_core.cpp:313-314`): tolerance ~1e-6 for cells that get rounded to 1e-14 elsewhere (line 184). The README itself flags this as a moderate issue. Anything in a NumPy port should use vectorised `np.isclose` with a documented tolerance, not bit-rounded comparison.

6. **`addNewCells` complexity is staggering** (`tooth_model_division.cpp:6-739`, ~730 lines for one logical operation). The C++ comments are diligent about pointing back to FORTRAN line numbers, which helps, but the function reads as a transliteration with `goto`-eliminated control flow rather than as a re-thought implementation. Six separate panic exits, three different in-place neighbour-array swaps, and a 'temporary copy of full neighbour array' allocation per external-cell swap (line 712, inside an O(numNewCells) loop with O(numCells * MAX_NEIGHBORS) updates per iteration). This is the function that needs the most rethinking before a NumPy port; reading the C++ here will not save much over reading the FORTRAN directly.

## File-split appropriateness for Python (question 5)

The C++ split is shaped by C++'s declaration/definition separation and by some functions being large enough to warrant a file each. Mapping to Python:

- `tooth_model.hpp` collapses into module-level data structures and class definitions in Python — no separate header needed.
- `tooth_model_iteration.cpp` is 36 lines of orchestration; in Python it would be one method on the model class, not a file.
- `tooth_model_core.cpp` (constructor + initialisation) is also small; in Python it is reasonably part of the model class itself.
- `tooth_model_geometry.cpp`, `tooth_model_diffusion.cpp`, `tooth_model_division.cpp` are the three substantive computational areas. Each is large enough to deserve its own module in Python.
- The missing `tooth_model_mechanics.cpp` (475 lines per README) is similarly substantive and would warrant its own module.
- `file_io.cpp` (claimed 663 lines) is bloated by humppa-format reading, ToothMaker-format reading, binary I/O, OFF mesh export, parameter store/load, and visualisation colour mapping. In Python this is reasonably split into `io_humppa.py`, `io_toothmaker.py`, `mesh_export.py`, and `viz.py`.
- `main.cpp` is the CLI entry point; in Python it is `__main__.py` or `cli.py`.

A reasonable Python layout: a single `tooth_model/` package with `__init__.py` containing the model class and the iteration orchestrator, plus `geometry.py`, `diffusion.py`, `mechanics.py`, `division.py`, `io.py`, `mesh.py`, `cli.py`. The C++ split gives a useful indication of the natural seam lines. The `iteration` and `core` files do not survive as separate modules.

## Glossary fidelity (question 6)

The README's variable table is largely accurate where it can be checked against the present sources. Spot checks:

- `quantities3D` — used as `quantities3D[i][kk][k]` throughout `tooth_model_diffusion.cpp` ✓ matches the README's `[cell][z][type]` description.
- `cellMargins` — declared `std::vector<std::vector<std::array<double, 8>>>` (`tooth_model.hpp:58`); the README does not document the inner `[8]` (margin x/y/z plus distance components 3-7). Slight under-documentation; the README says 'internode positions' but indices 3-7 carry distances and direction vectors (`tooth_model_division.cpp:613-617`). Worth expanding.
- `borderWidth` is flagged in the README as 'actually differentiation rate' (line 99), and the C++ source carries the same warning at `tooth_model.hpp:124` and `tooth_model_diffusion.cpp:267` — this is one of the rare places where the README, the C++ comment, and the actual usage all agree on the misnomer.
- `numBorderCells` semantic difference is noted at README:118 — but only there, not in the source comment near the off-by-one (`tooth_model_core.cpp:226-228`), where it is mentioned but not as load-bearing as the README suggests.
- `diffusionCoeffs2D[1] == Boy (buoyancy, not a diffusion coeff)` — flagged in README, in `tooth_model.hpp:70`, but I cannot verify the actual usage because that's in `calculateBuoyancy` which is in the missing mechanics file.

The README is internally consistent with the present source files and explicitly calls out the FORTRAN's parameter-name confusions. For glossary purposes it is reliable; for verifying that the *semantics* match it is incomplete because it leans on functions that were not delivered.

## Anything else load-bearing (question 7)

The 'Removed Unused Code' list (README:317-352) and 'Improvement Suggestions' (README:355-415) reveal that this port has been through a cleanup pass — at least 20 GUI-only members were stripped, several debug helpers removed, parameter-store/load was kept despite an admitted index-overlap bug. This is good hygiene and indicates the porter understood what they were doing. But the same pass added `if (suma > 0)` guards that change semantics relative to FORTRAN (point 2(a) above), introduced a 'first empty slot' shortcut that changes structural assumptions about hexagonal slot ordering, and produced an off-by-one in the panic threshold for the topology walk. The cleanup is mixed with silent semantic edits, and there is no test artefact with which to distinguish 'intended cleanup' from 'inadvertent semantic drift'.

## Verdict

**Use partially, with restraint.** The C++ port is best treated as a heavily-annotated reading aid for the FORTRAN — its footnoted line-number references back to `humppa_translate.f90` are genuinely useful, the variable-rename table in the README accelerates comprehension of what each FORTRAN identifier means, and the iteration order and broad function decomposition translate directly to a Python module layout. None of that requires the port to be runnable. Beyond that, this port is not a sound structural template: three of its files are absent in the mirror so half the simulator (mechanics + I/O + entrypoint) cannot be inspected, the present files contain at least one Path-A-equivalent bug (cross-pair `iiii` reset in `addNewCells`) and one different-from-FORTRAN-and-different-from-our-Python semantic in `apply_diffusion`'s degenerate path, the README's cross-validation claim is unbacked by any artefact in the mirror, and the code is unvectorised per-cell scalar loops that would have to be rewritten almost in their entirety for NumPy. For Path B I would read it the way one reads an interlinear gloss — useful for comprehension, not authoritative — and would not import any of its design choices wholesale. The FORTRAN itself, plus our existing `coreop2d.py` (which carries the three Path-A fixes verified against the FORTRAN binary), is the better source of truth.
