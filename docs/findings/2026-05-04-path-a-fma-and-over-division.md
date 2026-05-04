---
title: 'Path A findings — FMA framing rebutted, over-division root-caused'
author: Lyndon Drake (with Claude Code)
date: 2026-05-04
---

## Summary

Path A is closed (PR #16). Three coupled bugs in `coreop2d.py` were fixed,
the smoke test was tightened to assert FORTRAN-matching cell count and finite
coordinates, and the Python simulator now reaches the same equilibrium as the
FORTRAN reference on `examples/seal.txt` (57 cells, 236 faces, no NaN
z-coords). What follows is what was learnt about the system that *doesn't*
fit on a commit message.

## The FMA framing was partially wrong

The handover (HANDOVER.md, now removed) and the docstring on `Coreop2d`
described a 'chain of accidents': `add_cell`'s txungu branch deliberately
produces single-neighbour cells; `calculate_margins`'s label-77 fall-through
fills empty border slots with the previous slot's value; `apply_diffusion`'s
fallback then computes a cross product of bit-equal vectors that FMA turns
into ulp-noise instead of zero; that ε noise makes `sum_a` finite, division
proceeds, vertical diffusion happens. The framing implied that the FORTRAN's
equilibrium *depends* on this entire chain, and that breaking any link in
isolation would either NaN-cascade or exponentially over-divide.

The investigation rebuts the strong form of this claim.

- **FMA-as-divide-by-zero-dodge: confirmed.** When `sum_a == 0` in
  `apply_diffusion` (genuinely degenerate single-neighbour cells from the
  txungu branch), the FORTRAN survives only because FMA in the cross
  product leaves an ulp-level non-zero residual. Disabling FMA via
  `gfortran -ffp-contract=off` reproduces the NaN cascade. Pure Python's
  multiplication and subtraction never use FMA, so the cross product
  is exactly zero and the division NaNs.
- **FMA-as-equilibrium-stabiliser-everywhere: rebutted.** Once the
  FMA-limit guard is in place at that one site (the math limit is
  `area_bottom = 0.5` first norm, `1.0` second norm; `pes = 0`), Python's
  positions match the FORTRAN binary bit-for-bit through iteration 1,
  and within ulp tolerance for the entire 500-iteration run. The
  exponential over-division observed previously was not an FMA-chain
  symptom — it was an honest translation bug in `update_cell_position`,
  multiplying y-axis cells' z-force by `Bgr` (=10×) when FORTRAN does
  not. After fixing the Bgr placement, Python's equilibrium becomes
  insensitive to FMA at every other cross product in the simulator.

The 'must be fixed together' caveat was therefore correct in observation
(applying the FMA guard alone *did* expose the over-division) but wrong
in mechanism. The FMA guard didn't *cause* the over-division; the over-
division was already latent. The NaN cascade had been masking it by
killing the simulation's force vectors before they could move y-axis
cells too far.

## What the 'load-bearing accidents' actually are

The handover treated three FORTRAN-side phenomena as a co-evolved
system. With more visibility, they're better described as one
undefined-behaviour reliance plus several bugs whose only saving grace
is a panic check:

- **`apply_diffusion`'s sum_a == 0 path** is a genuine FMA dependency,
  but only as a way to dodge division by zero. The FORTRAN code never
  intended to compute diffusion for degenerate cells; it gets a 'small
  but finite' contribution that happens not to ruin anything because
  the cell is degenerate enough to dominate nothing else.
- **`add_cell`'s txungu branch** is a probable typo (FORTRAN
  13.f90:1158-1163 tests `jjjj` where the author almost certainly meant
  `jjj`) that produces 1-neighbour cells the simulator can't handle.
  These cells exist because of the bug, not because the model needs
  them.
- **`add_cell`'s topology walk** doesn't *handle* degenerate cells, it
  *gives up on them*. When the walk gets stuck on a one-neighbour cell,
  `cj` exceeds `nvmax` and the routine returns with `panic=1`, leaving
  the simulation state untouched for that iteration. The Python
  translation looked smarter (kept walking past the stuck cell) but
  was wrong, because the panic IS the only thing preventing bogus cell
  additions in this corner.
- **`calculate_margins`'s label-77 fall-through** fills empty border
  slots with the previous slot's value because of the way `goto 77`
  cascades in the FORTRAN. The Python translation zeroes empty slots,
  which differs in principle, but in practice doesn't break parity for
  the seal parameter set — and the divergence is at ulp level after
  iteration 1 with the three Path-A fixes in place. So this site is
  documented but not actively load-bearing for the equilibrium.

## Practical implication for Path B

A from-scratch implementation that follows the Zimm et al. PNAS 2023
equations directly will not need any of these accidents:

- No reason to ever produce a single-neighbour cell, so no FMA-limit
  guard is needed. The whole `apply_diffusion` `sum_a == 0` branch can
  go away.
- The Bgr-on-y-axis behaviour is a per-cell decision that should be
  made deliberately — either 'all border cells in the central band get
  Bgr' or 'only off-axis ones do' — not as a side-effect of x-sign
  branching that excludes the y-axis by accident.
- The topology walk should either converge or fail loudly. The
  FORTRAN's silent-panic-and-return behaviour preserved stability at
  the cost of dropping requested cell additions; a clean rewrite
  should either not need to walk topologies (vectorised numpy) or
  should walk them in a way that always terminates.

## Verification of byte-equality across machines

The committed `tests/golden_fortran/*.off` files reproduce byte-for-byte
from a fresh `gfortran -O2 -ffree-line-length-none -std=legacy` build
on this Linux x86_64 system. Per the handover, they also reproduce on
Apple Silicon (ARMv8). So at least for `-O2` on these two ISAs, the FMA
residual at the one degenerate-cell site evidently produces the same
downstream rounding choices on both, and the goldens are stable. A
future Path-B implementation could reasonably assert
structural-equivalence-plus-ulp-tolerance against these goldens without
worrying about platform drift.

## Investigation method (for the record)

The technique that landed all three fixes was the same: per-stage
stderr tracing on both implementations, diffed at progressively finer
granularity until each divergence had a single root cause.

1. **Per-iteration cell counts** — found that Python diverges from
   FORTRAN at iteration 8 (FORTRAN adds 4 cells, Python adds 0), then
   at iteration 13 (Python adds 16, FORTRAN adds 0).
2. **Per-iteration positions** — found cell 6 z differs by 0.225 in
   iteration 1 alone, before any cell additions.
3. **Per-stage positions within iteration 1** — narrowed the iter-1
   divergence to `update_cell_position` exactly. All earlier stages
   matched bit-for-bit.
4. **Per-stage forces** — confirmed forces match through
   `repulse_neighbour`, isolating the bug to the post-force-computation
   stages, and trivially to `update_cell_position` since
   `apply_nuclear_traction` did not change positions for cell 6.
5. **Reading FORTRAN side-by-side with the Python** — found the Bgr
   placement bug.
6. **Re-running with fix 2 applied** — Python now matches FORTRAN
   through iteration 12 then diverges at iteration 13's `add_cell`.
7. **Per-pair `cj` progression in the topology walk** — found that
   FORTRAN gets stuck (iiii=42 forever, eventual panic) while Python
   walks the topology successfully but adds wrong cells. Reading the
   FORTRAN found the iiii-preservation semantic.

The full debug instrumentation was on a throwaway `debug/fortran-trace`
branch (not pushed). Diagnostic Python scripts in `.tmp/` (gitignored)
diffed the traces. Both are removed at the end of the session — they
were scaffolding, not artifacts. The pattern is in the handover's
'Instrumentation patterns to copy' section if needed again.
