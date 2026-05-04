---
title: 'tgrohens/toothmaker — code review for Path B planning'
author: Lyndon Drake (review notes)
date: 2026-05-04
---

## Scope

Read-only comparison of `tgrohens/toothmaker`
(`/home/lyndon/repo/project/ToothmakerTranslation/.tmp/tgrohens-toothmaker/`)
against the Catalan-named ancestor
`/home/lyndon/repo/project/ToothmakerTranslation/.tmp/jernvall-toothmaker/humppa_translate.f90`
and our English-renamed reference
`/home/lyndon/repo/project/ToothmakerTranslation/13.f90`.
The question is whether tgrohens fixes any of the three Path A bugs.

## Path A bug status in tgrohens

### Bug 1 — `reaccio_difusio` divide-by-zero

The body of `reaccio_difusio` is byte-identical between humppa and tgrohens
(`diff` of `humppa_translate.f90:483-653` vs `coreop2d.f90:401-571` returns
empty). The unguarded division at the head of the cell loop is preserved:

- humppa `humppa_translate.f90:516`:
  `suma=sum(pes(i,:))+2*areasota ; areasota=areasota/suma ; pes(i,:)=pes(i,:)/suma`
- tgrohens `coreop2d.f90:434`: identical text.

There is no `if (suma > 0)` guard, no early-return for degenerate
single-neighbour cells, and the cross-product formulation at lines
`coreop2d.f90:423`/`humppa_translate.f90:505` is unchanged. The Path A
finding — that gfortran's default FMA contract is what saves the cross
product from being bit-identical to zero — applies verbatim to tgrohens.

**Status: bug present, no fix.**

### Bug 2 — `actualitza` Bgr-equivalent placement

`actualitza` (Catalan name; renamed `update_cell_position` in `13.f90`)
is also byte-identical between humppa (`humppa_translate.f90:1044-1070`)
and tgrohens (`coreop2d.f90:984-1010`). The Bgr-equivalent
parameter is `fac` in both Catalan files (`coreop2d.f90:68`,
`humppa_translate.f90:74`); `13.f90:1578` confirms the rename
`Bgr=parap(21,imap)` with the in-source comment 'formerly fac'.

The `fac` multiplications are placed correctly inside each x-sign branch,
matching `13.f90:866` and `:869`:

- tgrohens `coreop2d.f90:989-994`:
  `if (malla(i,1)>0) then ; hmalla(i,1)=hmalla(i,1)*bia ; hmalla(i,3)=hmalla(i,3)*fac ; end if`
  `if (malla(i,1)<0) then ; hmalla(i,1)=hmalla(i,1)*bip ; hmalla(i,3)=hmalla(i,3)*fac ; end if`

So Bug 2 is a translation artefact our Python migration introduced when
moving from FORTRAN's nested control flow to a NumPy boolean-mask
formulation. The Catalan FORTRAN — both jernvall's monolith and tgrohens'
fork — has the correct structure. tgrohens offers no help on Bug 2 because
neither FORTRAN file ever had the bug.

**Status: not a FORTRAN bug; tgrohens preserves the correct structure.**

### Bug 3 — `afegircel` `iiii` topology-walk preservation

`afegircel` is byte-identical between humppa (`humppa_translate.f90:1086-1553`,
468 lines) and tgrohens (`coreop2d.f90:1026-1493`, 468 lines): `diff`
returns zero. The critical line at the `88` label, which Path A's
investigation found relies on `iiii` retaining its previous-iteration
value, reads:

- humppa `humppa_translate.f90:1250`: `88    iii=iiii ; jjj=jjjj ; kkk=0`
- tgrohens `coreop2d.f90:1190`: `88    iii=iiii ; jjj=jjjj ; kkk=0`
- 13.f90 `13.f90:1054`: `88  iii=iiii ; jjj=jjjj` (translator dropped the
  `kkk=0` from the same line; harmless because `kkk=0` is set again at line 1056)

Neither FORTRAN ever resets `iiii` at the loop head. As with Bug 2,
this is a Python translation artefact we introduced; the FORTRAN
preservation of `iiii`-across-iterations was implicit in the lack of an
explicit reset, and our migration to a Python control-flow rewrite added
one. tgrohens contributes no fix because there was never a FORTRAN bug.

**Status: not a FORTRAN bug; tgrohens preserves the correct structure.**

## The `jjjj` typo at the `rtt:` loop

The chain you traced uses 'humppa line ~1158' as a landmark, but the
actual `rtt:` loop with the dead-code symptom is at humppa
`humppa_translate.f90:1336-1345` and tgrohens `coreop2d.f90:1276-1285`:

```
rtt:      do kkk=1,cj
            jjjj=pillats(kkk)
            if (jjjj==fi) then
              kkkk=kkkk+1 ; ccvei(i,kkkk)=fi
            else
              if (jjjj>ncels) then
                kkkk=kkkk+1 ; ccvei(i,kkkk)=jjjj
              end if
            end if
          end do rtt
goto 899
```

Humppa and tgrohens both write `jjjj=pillats(kkk)` and then test `jjjj`
— internally consistent. The dead-code aspect is not a within-loop typo
but the unconditional `goto 899` immediately afterwards (humppa
`:1346`, tgrohens `:1286`, 13.f90 `:1166`), which jumps over a large
disabled block to the continuation at 899 (humppa `:1369`, tgrohens
`:1309`, 13.f90 `:1203`). Our `13.f90:1157` does have a `jjj=pillats(kkk)`
followed by `if (jjjj==fi)` — that is a translation typo introduced when
porting the loop body, not a FORTRAN typo carried over from upstream.

**Status: tgrohens matches humppa exactly; the load-bearing 'dead code'
is the `goto 899` skip, not the loop body.**

## What tgrohens has removed

humppa is a single 2,401-line monolithic file containing module
`coreop2d` (`:25-1734`), module `esclec` (`:1739-2242`), and
program `tresdac` (`:2251-2396`). tgrohens splits these into three
files (`coreop2d.f90`, `esclec.f90`, `toothmaker.f90`). The line-count
delta (humppa 2,401 → tgrohens 1,635 + 392 + 141 = 2,168, a saving of
~233 lines, or ~10%) comes from removing dead/legacy subroutines, not
from semantic changes.

Removed from `coreop2d`:

- `redime` (humppa `:332-369`, marked 'not in use' in the header at line 5)
- `referci` (`:371-377`), `refercid` (`:379-383`) — not referenced
  elsewhere
- `llegirparam` (`:397-418`) — superseded by `llegirparatxt` in the
  esclec module
- `controlz` (`:1649-1663`), `vbia` (`:1665-1686`) — neither subroutine
  is referenced from `iteracio` or anywhere reachable

Removed from `esclec`:

- `llegirpara` (humppa `:1845-1862`), `llegirforma` (`:1864-1870`),
  `llegirex` (`:1872-1878`), `llegirknots` (`:1880-1890`),
  `llegirveins` (`:1892-1902`), `llegir` (`:1943-2021`) — binary-format
  load functions, replaced by the simpler text-only `llegirinicial` →
  `llegirparatxt` flow
- `guardaveinsoff` (`:2031-2191`, 161 lines) — replaced by the much
  shorter `guardaveinsoff_2` (tgrohens `esclec.f90:149-298`, 150 lines)
- Various OpenGL hooks (commented out in both files; the `use opengl_gl`
  line is commented in both)

Nothing load-bearing for the `iteracio`-driven simulation loop has
been cut. The diferenciacio (3 lines), iteracio, empu, stelate,
biaixbl, pushing, pushingnovei, promig, calculmarges, ci, ciinicial,
posar, dime, ordenarepe, perextrems and afegircel subroutines are
either byte-identical or differ only in trailing whitespace. Spot
checks: `diff humppa:666-790 vs tgrohens:584-708` for `empu` returns
empty; `:794-821 / :712-739` for `stelate` returns empty; `:966-980 /
:906-920` for `biaixbl` returns empty; the only diff in `dime` is one
whitespace line.

## What tgrohens has modified

Almost nothing. The module-public list at `coreop2d.f90:25` drops
`llegirparam` (since the routine was removed); otherwise the public API
is identical to humppa `:31`. No Fortran-2003 features, no `deallocate`
discipline added, no parameter additions, no renaming of variables —
the file still uses Catalan-named globals (`malla`, `marge`, `vei`,
`hmalla`, `nveins`, `q3d`, `q2d`, `bip/bia/bil/bib`, `fac`, `tadi`,
`crema`, `acac`, etc.). Our `13.f90` is the only one of the three that
renames to English. The example input file uses the same 30-line
parameter layout as `examples/seal.txt`, with one cosmetic label
mismatch on line 28 (`Rcb_Radius_center_bias` vs our `Bwi`); both load
into `parap(28)` via `llegirparatxt` (`esclec.f90:83-91` reading
file lines 1..30 into `parap(3..32)`), and `posarparap` then assigns
`tadif=parap(20)` (Bwi) and `radibii=parap(30)` (Rcb), so the labels
are actually orthogonal — `seal.txt` happens to label the slot at
file-line 28 as Bwi while tgrohens' example labels file-line 28 as
Rcb, but both routines parse the values blindly by position.

## Compilation flags

`/home/lyndon/repo/project/ToothmakerTranslation/.tmp/tgrohens-toothmaker/Makefile:5`:

```
gfortran -O3 $(SOURCES) -o runt.e
```

This is **worse** for our purposes than the documented Path A baseline
of `-O2`. Critical points:

1. `-O3` enables `-ffp-contract=fast` (the gfortran default for any
   non-`-Og` optimisation level), so FMA contraction is on. The Bug 1
   workaround — having `(uy*dz-uz*dy)`, `(uz*dx-ux*dz)` and
   `(ux*dy-uy*dx)` accumulate FMA noise that prevents `areap(i,j)=0`
   for collinear single-neighbour cases — is therefore active in
   tgrohens just as in our Path A `-O2` baseline. Building tgrohens
   with `-ffp-contract=off` would expose the same divide-by-zero we
   reproduced for `13.f90`.
2. There is no `-march=native`, no `-fno-finite-math-only`, no
   `-ffloat-store`, and no explicit FMA control. `-O3` over `-O2` is
   primarily about enabling vectorisation and inlining; on this code's
   loop structure (heavy `nvmax=30`-bounded scalar control flow with
   tight panic-and-return paths) the difference will be small.
3. The debug target adds `-g -pg` (line 8) but no array-bounds checks,
   no `-ffpe-trap=invalid,zero,overflow`. There is no instrumentation
   for catching the divide-by-zero even in 'debug' mode.

## The `nosea_10261__0.txt` example input

It is the same physical 30-line text format as `examples/seal.txt` —
parameter values one per line, file-line 1 mapping to `parap(3)`
(tacre/Egr) via `llegirparatxt` (`esclec.f90:83-91`). The values are
visibly different from `seal.txt`:

- Egr 0.0135 (vs seal 0.5 — much slower epithelial growth)
- Mgr 204.6 (vs 16000 — far weaker mesenchyme)
- Act 1.10 (vs 0.15)
- Inh 26.0 (vs 1.0 — much stronger inhibition)
- Sec 0.026 (vs 3.0)
- Da 0.20 (vs 0.008 — 25× faster activator diffusion)
- Bgr 1.0 (vs 10.0)
- Rad 2 (vs 4 — *7-cell* initial condition vs the seal's 19)

The `Rad=2` (7 cells) is explicit in the parameter label
('2_for_having_7_cells'). The filename `nosea_` (literally 'no-seal-'
in colloquial Catalan) plus values like Inh=26 and the very small Egr
suggest this is a non-seal parameter set deliberately tuned for slower,
finer-grained diffusion-driven dynamics — possibly a tribosphenic
mammal exemplar from the Salazar-Ciudad / Jernvall parameter sweep
work, but without the README or a paper reference I cannot identify the
target organism. It is not a seal-equivalent input; it would not load
or evolve like our `seal.txt` golden.

## Genuinely surprising observation

tgrohens' README is the single line 'This repo holds my work on the
core simulation code of the [ToothMaker repo]' — yet the diff to
the upstream Catalan core is essentially zero except for dead-code
removal and a file split. The fork has not begun the modernisation
work the framing implies; it is a clean substrate for future work
rather than an evolved branch. For our purposes this is good news —
tgrohens cannot have introduced silent semantic divergences we would
need to chase.

## Verdict

tgrohens is a useful negative cross-check but does not change any Path
A finding. It confirms that all three Path A bugs are properties of the
Catalan-named upstream that survive into both `humppa_translate.f90`
and tgrohens unchanged: Bug 1 (the FMA-shielded divide-by-zero in
`reaccio_difusio`) is preserved verbatim and would manifest in tgrohens
under `-ffp-contract=off`; Bug 2 (Bgr branch placement) and Bug 3
(`iiii` topology-walk preservation) were not FORTRAN bugs at all but
artefacts of our Python migration, and tgrohens preserves the
correct FORTRAN structure as one would expect. tgrohens offers no
parameter-handling improvements, no English renaming, no `deallocate`
discipline, and no FMA-disabling compile flags — so for Path B planning
it adds two things only: (i) a clean three-file split that is easier to
diff than the 2,401-line monolith, and (ii) one extra example
parameter file (`nosea_10261__0.txt`) that exercises a non-seal regime
with `Rad=2` (7-cell initial condition) and very different growth/
diffusion balance, which could be useful as a second golden if we
care about generalisation beyond the seal case.
