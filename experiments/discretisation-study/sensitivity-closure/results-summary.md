# Sensitivity-closure runs — closing the §5E unknowns

Generated: 2026-05-05.

This directory holds the targeted runs that close two cells of the
§5E parameter-set sensitivity table previously marked `unknown`:

1. `mesenchyme` under FORTRAN-flavoured presets on
   `examples/wt-tribosphenic-2014.txt`.
2. `n_mesenchyme_layers` under paper-flavoured presets on the same.

Both runs use `--iters 100 --saves 5` for a 500-iteration total
budget per configuration, the same cadence as the existing
comparison study.

## (i) `mesenchyme = absent` vs `per_column_z_layers` under LEGACY_FORTRAN

```
Subdir                Cells  Cusps  z range
legacy-mes-absent     19     7      [12.9614, 13.8233]
legacy-mes-percol     19     7      [12.9614, 13.8233]
```

**Result: bit-for-bit identical.** Cells, cusp count, and per-vertex
positions all agree to 1e-9. The mesenchyme branch *executes* every
iteration, but Inh and Sec remain at zero throughout the FORTRAN-
flavoured run (a consequence of `eq18_sec_source='k_sec_times_di'`
combined with `diff_accumulator='sec'` locking d_i at 0; documented
in §5 of the briefing). With nothing to diffuse vertically, the
choice between `absent` and `per_column_z_layers` is observationally
inert on this parameter set.

## (ii) `n_mesenchyme_layers` ∈ {1, 2, 3} under PATH_B_DEFAULT

```
Subdir              Cells  Cusps  z range
paper-mes-1layer    19     19     [1.000617, 1.013941]
paper-mes-2layer    19     19     [1.000617, 1.013941]
paper-mes-3layer    19     19     [1.000617, 1.013941]
```

**Result: bit-for-bit identical across the three layer counts.**
Per-vertex positions agree to numerical precision. Under PATH_B_DEFAULT
the eq.17 secrete clause does fire (every cell becomes a knot by
iteration ~100), so Inh and Sec accumulate in the epithelial layer
and propagate to mesenchyme via the vertical Laplacian. But the per-
iteration vertical flux is small (k_di and k_ds are 0.2 and 0.1
respectively in the wt-tribosphenic parameter set; the dt is 0.05),
so over 2500 iterations the additional substrate layers' presence
makes no observable difference at the cell-count or cusp-count
resolution. The choice of layer count is observationally inert on
this parameter set within this iteration budget.

## Implication for §5E

Both cells previously marked `unknown` resolve to **'does not fire
observably'** — the branch executes mechanically but produces no
state change visible at the resolution of cell count, cusp count,
or vertex envelope. This is distinct from a branch that genuinely
does not execute (e.g. `eq14_denominator` on the seal example,
where Inh = 0 makes both denominator forms reduce to 1). Both are
dormant for biological consequence purposes, but the audit
classification differs.

A briefing footnote in §5E now distinguishes:

- **fires** — branch executes AND output observably depends on the choice.
- **does not fire** — branch does not execute (the regime is never reached).
- **inert** — branch executes but output does not depend on the choice
  (the regime is reached but the choice is mathematically immaterial,
  or the time budget is too short to amplify the difference).

Both `mesenchyme` and `n_mesenchyme_layers` on
`wt-tribosphenic-2014.txt` are **inert** under the relevant preset
families. A longer iteration budget or a parameter set with larger
k_di / k_ds is the natural follow-up to test whether they would
become observable.

### Update (2026-05-05): re-test under `laplacian=fortran_margins`

Path B v2 B3 (2026-05-05) wired the in-plane fortran_margins
operator. Re-running case (ii) under `laplacian=fortran_margins`
(`experiments/discretisation-study/sensitivity-closure-fortran-margins-2026-05-05/`)
shows:

```
Subdir                          1/2/3 layer hashes
paper-mes-{1,2,3}layer-fortran  5ff49087634d (all three identical)
```

`n_mesenchyme_layers` remains **inert** under fortran_margins —
the three layer counts produce byte-identical OFF outputs. The
absolute output differs from the length_weighted runs (consistent
with the −57 % cusp_width shift on `PATH_B_DEFAULT_plus_laplacian`
documented in
`docs/findings/2026-05-05-path-b-v2-fortran-margins-implementation.md`),
but the layer-count sensitivity is unchanged.

Reading: the operator's choice shifts where PATH_B_DEFAULT
equilibrates, but does not unlock the substrate-edge layers'
contribution to the dynamics. The §5E inert classification for
`n_mesenchyme_layers` therefore survives the operator change.

## Reproduction

```bash
TMPDIR=.tmp .venv/bin/python -m silicoshark examples/wt-tribosphenic-2014.txt \
    experiments/discretisation-study/sensitivity-closure/legacy-mes-absent run 100 5 \
    --preset LEGACY_FORTRAN \
    --override mesenchyme=absent --override laplacian=length_weighted

TMPDIR=.tmp .venv/bin/python -m silicoshark examples/wt-tribosphenic-2014.txt \
    experiments/discretisation-study/sensitivity-closure/legacy-mes-percol run 100 5 \
    --preset LEGACY_FORTRAN \
    --override laplacian=length_weighted

for n in 1 2 3; do
  TMPDIR=.tmp .venv/bin/python -m silicoshark examples/wt-tribosphenic-2014.txt \
      experiments/discretisation-study/sensitivity-closure/paper-mes-${n}layer run 100 5 \
      --preset PATH_B_DEFAULT \
      --override n_mesenchyme_layers=$n --override laplacian=length_weighted
done

# 2026-05-05 fortran_margins re-test
for n in 1 2 3; do
  TMPDIR=.tmp .venv/bin/python -m silicoshark examples/wt-tribosphenic-2014.txt \
      experiments/discretisation-study/sensitivity-closure-fortran-margins-2026-05-05/paper-mes-${n}layer-fortran \
      run 100 5 --preset PATH_B_DEFAULT \
      --override n_mesenchyme_layers=$n --override laplacian=fortran_margins
done
```
