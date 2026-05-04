---
title: 'Divergences between `13.f90` and `humppa_translate.f90`'
author: Lyndon Drake
date: 2026-05-04
---

## Purpose and method

`13.f90` (project root) is a renamed-to-English copy of `humppa_translate.f90` (mirrored at `/home/lyndon/repo/project/ToothmakerTranslation/.tmp/jernvall-toothmaker/humppa_translate.f90`). I have catalogued every divergence I could identify, classified as cosmetic (a), structural (b), or semantic (c). The renaming map used is the variable and subroutine table from `/home/lyndon/repo/project/ToothmakerTranslation/.tmp/jernvall-toothmaker/cpp_README.md`. I applied the renames to humppa via `sed` (script saved at `.tmp/divergence-analysis/rename.sed`) and then took per-subroutine unified diffs against `13.f90` (script `.tmp/divergence-analysis/diff_all.sh`, output `.tmp/divergence-analysis/all_diffs.txt`). Coverage: every subroutine in `module coreop2d` was diffed; the `program tresdac` and `module esclec` sections were inspected directly. Lines cited below are checked file:line locations.

## Cosmetic divergences (a)

These do not affect behaviour:

1. **Variable and subroutine renames** throughout, e.g. `malla` → `positions`, `vei` → `neigh`, `marge` → `border`, `q2d(:,1)` → `DiffState(:)`, `ncels` → `num_active_cells`, `ncals` → `num_all_cells`, `ncils` → `first_border_cell`, `ng=5/ngg=4` → `num_species_in_q3d=3`, plus all parameter renames per `cpp_README.md`. Subroutine renames: `ciinicial` → `initialconditions`, `dime` → `allocateinitialstate`, `posar` → `initializecellpositions`, `calculmarges` → `calculatemargins`, `reaccio_difusio` → `applydiffusion`, `diferenciacio` → `applydifferentiation`, `empu` → `EpGrowthBorderForce`, `stelate` → `BoyForce`, `pushing` → `repulseneighbor`, `pushingnovei` → `repelnonneigh`, `biaixbl` → `applyborderbias`, `promig` → `applynucleartraction`, `actualitza` → `updatecellposition`, `afegircel` → `addcell`, `perextrems` → `updatebordercells`, `iteracio` → `iteration`. The `setparams` subroutine in `13.f90` corresponds to `posarparap` in humppa — see semantic note 7 below.
2. **Catalan-to-English comment translations** throughout (e.g. humppa `!em de passar de veinatge a faces` → 13.f90 `!I have to go from neighboring to faces`).
3. **Indentation and whitespace adjustments**, plus splitting some single-line `;`-separated statements into multi-line blocks (`13.f90:284-388` `calculatemargins`; `13.f90:671-719` `repulseneighbor`).
4. **Removed in-line commented-out code** in many subroutines (e.g. dead `! q3d(:,:,i)=q3d(:,:,i)*0.075D1` in `addcell`, dead OpenGL constant in module preamble at `humppa:107` vs `13.f90:107`).
5. **Identifier overlap reduced**: humppa's `do(nvmax)` is gone in `13.f90` (`do` is a Fortran keyword used as identifier in humppa, removed in `13.f90:888` ff).

## Structural divergences (b)

Code present in one file and not the other:

1. **OpenGL hooks removed.** `humppa:29 use opengl_gl` → commented out. The visualisation member `submenuid` (`humppa:101`) and the unused real `gldouble pi` (`humppa:107`) are gone. The visualisation flag `vmarges` is renamed `showborders` and retained.
2. **Subroutines deleted from `13.f90` (present in humppa).** From `module coreop2d`: `ci` (`humppa:136-156`), `redime` (`humppa:332-369`), `referci` (`humppa:371-377`), `refercid` (`humppa:379-383`), `ordenarepe` (`humppa:1072-1084`), `controlz` (`humppa:1649-1663`), `vbia` (`humppa:1665-1686`). From `module esclec`: `guardapara`, `guardaforma`, `guardacon`, `guardaex`, `guardaknots`, `guardaveins`, `guardamarges`, `escriuparatxt`, `llegirpara`, `llegirforma`, `llegirex`, `llegirknots`, `llegirveins`, `agafarparap`, `posarparap`, `llegir`, `llegirinicial` — i.e. every binary `.dad` save/load and the parameter-array round-trip helpers (`humppa:1754-2029`).
3. **Subroutines added to `13.f90` (not in humppa).** `read_param_file` (`13.f90:1541-1552`, equivalent to humppa `llegirparatxt`), `setparams(imap)` (`13.f90:1556-1590`, equivalent to humppa `posarparap` — see semantic note 7), `initialize_from_parameter_file` (`13.f90:1594-1600`), `guardaveinsoff_2` (`13.f90:1604-1723`, replaces humppa `guardaveinsoff`), `get_rainbow_knot` (`13.f90:1780-1819`, used only by `guardaveinsoff_2`).
4. **Initial-condition setup simplifications.** `13.f90:135-266` `allocateinitialstate` drops humppa's commented-out alternate hexagonal layout (`humppa:215-220`), the `condme` parameter, the `ampl=Rad*0.75` line, and several large blocks of commented-out adjustment loops (`humppa:274-321`).
5. **Iteration loop file-logging removed.** `humppa:1707-1709` writes a progress-bar line per step; `13.f90:1492-1511` `iteration` does not.
6. **Diffusion array reduced from five species to three.** `q3d` third dimension is `ng=5` in humppa, `num_species_in_q3d=3` in `13.f90`. Inner loops `do k=1,4` (`humppa:517,552`) become `do k=1,num_species_in_q3d` (`13.f90:423,453`). Final update `do i=1,4 ; q3d(:,:,i)=q3d(:,:,i)+delta*difq3d(i)*hq3d(:,:,i)` (`humppa:588-590`) becomes three explicit assignments using the scalar parameters `Da`, `Di`, `Ds` (`13.f90:467-474`). The `q2d` array is removed entirely; cell-level state migrates to the 1-D array `DiffState(:)`. The `difq3d`/`difq2d` array machinery is gone; the relevant coefficients become scalars (`Da`, `Di`, `Ds`, `Boy`, `Dff`, `Bgr`).
7. **`program tresdac` rewritten.** `13.f90:1832-1920` takes five command-line arguments (input, output folder, output basename, iterations, save-block count); `humppa:2251-2396` takes four. `13.f90` writes only the `.off` mesh file via `guardaveinsoff_2`; humppa writes `.dad`, `.off`, and `.txt`. `13.f90` calls `setparams(1)`; humppa's equivalent flow is `agafarparap(1)` then `posarparap(1)` (`humppa:2368,2389-2390`). The trailing `print *,"fora"` (`humppa:2395`) is gone.
8. **`applyborderbias` no longer references the `Swi` band.** See semantic note 1 below — the explicit threshold is replaced with a hard-coded zero, dropping the `Swi` (humppa `tadi`) parameter from this subroutine. Since `13.f90`'s `setparams` (`13.f90:1556-1590`) skips `parap(6)` (commented `!Nothing`), the `Swi` parameter is never loaded anywhere in `13.f90` — only the `Swi=0` comment at `13.f90:1` and the dead-cite at `13.f90:882`. So this is a structural removal, not just a per-subroutine simplification.

## Semantic divergences (c)

These are actual behavioural differences. Each entry gives (i) file:line locations in both files, (ii) what differs, (iii) effect on simulation output, (iv) confidence.

### 1. `addcell` `txungu` else-branch typo (CONFIRMED, the bug already known)

- **Location**: `13.f90:1156-1157` (`rtt: do kkk=1,cj ; jjj=pillats(kkk)`); `humppa_translate.f90:1336-1337` (`rtt: do kkk=1,cj ; jjjj=pillats(kkk)`).
- **What differs**: humppa assigns `jjjj=pillats(kkk)` and then tests `if (jjjj==fi)` and `if (jjjj>ncels)` against the freshly read pillat. `13.f90` writes `jjj=pillats(kkk)` (three j's) but the subsequent tests at `13.f90:1158, 1161` read `jjjj` (four j's), so they read whatever `jjjj` contained from before entering the loop (the value left over from `13.f90:1115-1118` or earlier).
- **Effect**: in the `txungu` ('thorny') branch of cell-division topology resolution — the branch entered when `jjj>0`, i.e. when at least one of the surrounding nodes was itself a newly-added node — the loop fails to inspect each pillat individually and instead tests the same stale value `cj` times. The `temp_new_neigh(i,:)` neighbour list constructed for the new cell `i` will therefore omit some neighbours and may double-count others, depending on the value of `jjjj` carried in from outside. This produces incorrect connectivity in clustered division events, which then propagates through `calculatemargins` and the subsequent diffusion pass.
- **Confidence**: high. Verified by direct inspection of the two files at the cited lines and by `grep "pillats(kkk)"` across both.

### 2. `applyborderbias` threshold band collapse

- **Location**: `13.f90:780, 785` (`if (positions(i,2) < 0.0d0)` / `if (positions(i,2) > 0.0d0)`); `humppa_translate.f90:968, 973` (`if (malla(i,2)<-tadi)` / `if (malla(i,2)>tadi)`, where `tadi`=`Swi`).
- **What differs**: humppa applies the buccal/lingual bias only outside a `[-Swi, +Swi]` band around y=0; `13.f90` applies it in both half-planes with no dead band.
- **Effect**: when `Swi=0` (the case stated by `13.f90:1`), `< -0` is identical to `< 0` and `> 0` is identical to `> 0`, so the subroutines are equivalent. When `Swi != 0` they diverge. Because `13.f90`'s `setparams` does not load `Swi` from `parap` at all (see structural note 8), `Swi` is permanently 0 in any `13.f90` run — so no observable simulation effect under the current parameter pipeline.
- **Confidence**: high on the textual difference; high that there is no observable effect under the current pipeline because `Swi` is never written.

### 3. `applydiffusion` reaction term: ectodin coupling dropped

- **Location**: `13.f90:483` (`a = Act*q3d(i,1,1)`); `humppa_translate.f90:613` (`a = acac*q3d(i,1,1) - q3d(i,1,4)`).
- **What differs**: humppa subtracts `q3d(i,1,4)` (the species-4 ectodin/Not2 concentration) from the activator-production rate before the `a/(1+Inh*q3d(i,1,2))` Michaelis–Menten term. `13.f90` drops this term, and the entire species-4 update block at `humppa:635-637` (`a = acec*q3d(i,1,1) - mu*q3d(i,1,4) - acaca*q3d(i,1,3)` followed by `hq3d(i,1,4)=a`) is also gone.
- **Effect**: when `acec=0` and `acaca=0` (as `13.f90:1` declares: 'Ect = 0, Not3 = 0'), `q3d(:,:,4)` is initialised to zero and never produced, so the missing `- q3d(i,1,4)` term is identically zero in any `13.f90` run. This is the consistent dead-code removal that follows from species 4 having been eliminated. The coupling matters only if Not2/ectodin is enabled, which `13.f90` cannot do.
- **Confidence**: high on the textual difference; high that there is no observable effect under the current pipeline.

### 4. `applydiffusion` array-update collapse

- **Location**: `13.f90:467-474` (three explicit `q3d(:,:,k) = q3d(:,:,k) + delta * {Da,Di,Ds} * hq3d(:,:,k)` assignments); `humppa_translate.f90:587-590` (`do i=1,4 ; q3d(:,:,i)=q3d(:,:,i)+delta*difq3d(i)*hq3d(:,:,i) ; end do`).
- **What differs**: humppa updates four species (1=Activator, 2=Inhibitor, 3=Sec, 4=Ectodin) using array `difq3d`. `13.f90` updates only species 1–3 using scalars `Da`, `Di`, `Ds`. The mapping is `Da ↔ difq3d(1)`, `Di ↔ difq3d(2)`, `Ds ↔ difq3d(3)`, all loaded from the same `parap(13:15)` slots (per cpp_README and `13.f90:1569-1571` setparams; humppa `posarparap` `humppa:1930` `do j=1,ng ; difq3d(j)=parap(12+j,imap)`).
- **Effect**: numerical update of species 1–3 is byte-identical given the same parameter file. The dropped species-4 update is dead under the constraint `acec=0, acaca=0` (see note 3). No observable effect under the current pipeline.
- **Confidence**: high.

### 5. `applydiffusion` boundary-sink condition simplified

- **Location**: `13.f90:455, 491, 502` (boundary-sink branch always fires, gated only by `if (ii==num_all_cells)`); `humppa_translate.f90:570, 588, 605` (commented-out outer guard `! if (k/=1) then ... ! end if` left dangling around the same `if (ii==ncals)` guard).
- **What differs**: humppa has commented-out code that, if active, would have prevented the sink from being applied to species 1 (Activator). The comments in humppa say `! ACHTUNG aixo es perque als marges de la dent tenim molt activador`. In `13.f90` these comments are removed and the (already-active) sink for all species is preserved.
- **Effect**: none. Humppa's behaviour is identical to `13.f90`'s — the guards were already commented out in humppa. The textual divergence is removal of dormant comments, not a logic change.
- **Confidence**: high.

### 6. `addcell` post-division border-cell swap: `first_border_cell` increment guard

- **Location**: `13.f90:1384` (`first_border_cell=first_border_cell+1                              ! wxpand the contiguous border-cell block by 1first_border_cell`); `humppa_translate.f90:1544` (`ncils=ncils+1`).
- **What differs**: textual only — a typo'd English comment ('wxpand', 'block by 1first_border_cell') has been appended in `13.f90`. The increment itself is identical.
- **Effect**: none.
- **Confidence**: high.

### 7. `setparams` parameter-array index overlaps removed

- **Location**: `13.f90:1556-1590` `setparams(imap)`; `humppa_translate.f90:1925-1939` `posarparap(imap)`.
- **What differs**: humppa loads parameters into both scalar variables and the `difq3d`/`difq2d` arrays from overlapping slots: `parap(17,imap)` is read into both `us` (=`Int`) and `difq3d(5)`; `parap(18,imap)` is read into both `ud` (=`Set`) and `difq2d(1)`; `parap(19,imap)` into both the (unused) overlap and `difq2d(2)` (=`Boy`); `parap(20,imap)` into `difq2d(3)` (=`Dff`) and via `tadif=parap(14+ng+1,imap)` into `tadif`; `parap(21,imap)` into `difq2d(4)` and via `fac=parap(15+ng+1,imap)` into `fac` (=`Bgr`). The `cpp_README` describes this as 'Parameter Array Index Overlap'. `13.f90` removes the overlap by dropping the array machinery: each scalar reads from a single slot — `Set=parap(18)`, `Boy=parap(19)`, `Dff=parap(20)`, `Bgr=parap(21)`.
- **Effect**: for the same parameter file, the values loaded into `Da, Di, Ds, Set, Boy, Dff, Bgr, Int, Egr, Mgr, Rep, Adh, Act, Inh, Sec, Abi, Pbi, Bbi, Lbi, Rad, Deg, Dgr, Ntr, Bwi, ina, umgr` are identical in both files. The slots that humppa loads into species 4–5 of `difq3d` (`parap(16:17)`) and into the unused `difq2d` overlaps are simply not read in `13.f90`. The `Swi` slot at `parap(6,imap)` is read into `tadi` in humppa but skipped (commented `!Nothing`) in `13.f90`. So `Swi` is the one parameter genuinely lost — it would always have value 0 in any `13.f90` run regardless of what is in the input file. See semantic note 2 for the downstream consequence.
- **Confidence**: high on the slot-by-slot mapping (verified by `cpp_README` and direct inspection); medium on whether any production input file ever sets `Swi != 0` (not investigated within scope).

### 8. `guardaveinsoff_2` colour and OFF-header changes

- **Location**: `13.f90:1604-1723` `guardaveinsoff_2`; `humppa_translate.f90:2031-2191` `guardaveinsoff`.
- **What differs**: (i) the per-vertex colour call changes from `get_rainbow(ma(i),vamin,vamax,c)` to `get_rainbow_knot(ma(i),vamin,vamax,c,i)`, which paints knot cells cyan instead of using the rainbow gradient; (ii) the OFF header writes `0` for the edge count instead of `num_active_cells`; (iii) format specifiers widen from `4I4`/`5I4` to `4I10`/`5I10`; (iv) extensive dead deduplication code (`ja`, `nfacre`, `pasos`, `npasos`, `nfa`) is removed.
- **Effect**: only the visualisation `.off` output differs. The simulation state (positions, neigh, q3d, DiffState, knots, border) is unchanged.
- **Confidence**: high.

### 9. `iteration` file-logging side-effect removed

- **Location**: `13.f90:1492-1511`; `humppa_translate.f90:1688-1732`.
- **What differs**: humppa writes a per-step counter to `nfpro` via `open(1,file=nfpro,access='append') ; write (1,*) iteraciototal*(iti-1)+temps ; close(1)` (`humppa:1707-1709`). `13.f90` does not.
- **Effect**: no effect on simulation state. Side-effect on the filesystem only.
- **Confidence**: high.

### 10. `allocateinitialstate` removes redundant trailing `addcell` / `calculatemargins` calls

- **Location**: `13.f90:135-266`; `humppa_translate.f90:158-330`.
- **What differs**: humppa has commented-out blocks (`humppa:255-280`, beginning `!do iit=1,int(radibi*2)`) that, if active, would loop over `addcell`/`calculmarges`. In `13.f90` these comments are gone. The active code paths are identical (`call calculatemargins ; q3d=0` in both at the corresponding locations).
- **Effect**: none — the removed code was already inactive in humppa.
- **Confidence**: high.

## Coverage and unverified ground

- All 18 subroutines listed in the `module coreop2d` header comment of both files were diff'd via the rename script and inspected. The diffs are saved at `/home/lyndon/repo/project/ToothmakerTranslation/.tmp/divergence-analysis/all_diffs.txt` (2,008 lines).
- The `program tresdac` block, the `module esclec` headers (variable declarations, `read_param_file`/`llegirparatxt`, `setparams`/`posarparap`, `guardaveinsoff_2`/`guardaveinsoff`, `mat`, `get_rainbow`, `get_rainbow_knot`) were inspected directly.
- The deleted humppa subroutines (`ci`, `redime`, `referci`, `refercid`, `ordenarepe`, `controlz`, `vbia`, all `guarda*`/`llegir*`/`agafarparap`/`posarparap`/`llegir`/`llegirinicial`/`escriuparatxt`) were not inspected for further internal divergences from anything in `13.f90` because `13.f90` does not contain them.
- Sed-based renaming with `\b` word boundaries can in principle leave residual textual differences where an identifier appears inside a string literal or a comment; I have not exhaustively audited every comment and string. Any such residue would manifest as cosmetic-only diff lines and should not change behaviour.
- I did not search for divergences in subroutines added to `13.f90` (`read_param_file`, `setparams`, `initialize_from_parameter_file`, `guardaveinsoff_2`, `get_rainbow_knot`) against humppa's distantly-related counterparts beyond the comparisons noted in entries 7 and 8 above.

If a future audit needs to confirm absence of further bugs, the scaffolding is in place: `.tmp/divergence-analysis/rename.sed`, `.tmp/divergence-analysis/pairs.txt`, `.tmp/divergence-analysis/diff_all.sh`, and `.tmp/divergence-analysis/all_diffs.txt`. Re-running `diff_all.sh` regenerates the per-subroutine unified diffs against a fresh humppa copy.
