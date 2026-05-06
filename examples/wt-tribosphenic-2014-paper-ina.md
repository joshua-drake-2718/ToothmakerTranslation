---
title: 'Parameter set: 2014 wild-type tribosphenic mouse, paper-literal ina = 0'
date: 2026-05-05
---

## Source

Identical to `examples/wt-tribosphenic-2014.txt` except `ina = 0`,
matching the Harjunmaa et al. 2014 caption ('Ina is not used in the
model' for wild-type — i.e. `ina = 0.0`). All other 29 parameter values
are unchanged from the wild-type tribosphenic set in Ext Data Fig. 2.

## Purpose

This file exists to support the rebuttal §B.5 strengthening run in
`docs/paper-rebuttal-preemption-2026-05-05.md`. The headline cusp-forming
finding (PAPER_2010 = 19 cusps; PAPER_LITERAL_2010 = 0) currently uses
`ina = 0.5` to seed Act above the knot threshold within 500 iterations —
a deviation from the paper's literal `ina = 0`. The eq. 14 typo's
saturating effect depends on Act being non-zero, so a paper-literal
`ina = 0` run with extended iterations tests whether the typo's
catastrophic biology survives or vanishes when no Act seed is present.

Three outcomes are possible:

1. **Both PAPER presets produce zero cusps over the extended run** —
   confirms that Act-seeding (mesenchymal Sec back-diffusion or
   border-cell lift via `Bbi`/`Lbi`) does not bring Act above the knot
   threshold under silicoshark's mesenchyme-absent / static-with-local-
   update default within the iteration budget. The headline becomes
   'on parameter sets that drive Act above zero, the typo is
   catastrophic'.
2. **PAPER_2010 forms cusps (via border-bias seeding) but
   PAPER_LITERAL_2010 does not** — converts the headline finding from
   'on a parameter set with one tuning deviation' to 'on the paper-
   literal parameter set, given enough iterations'.
3. **Both PAPER presets form similar cusp counts** — would mean the
   typo is dormant on the paper-literal parameter set and the
   silicoshark replication is missing the mechanism that drives Act
   above zero in `13.f90`'s long-iteration runs. This would prompt
   further investigation of the residual seed mechanism.

## Deviations from the literal Ext Data Fig. 2 values

`umgr = 0.0`. Not in the 2014 paper at all. A post-2014 addition to
`13.f90` with no GUI control in Ext Data Fig. 2. Set to zero to
match its absence in the published parameter set, identically to
`wt-tribosphenic-2014.txt`.

There are no other deviations from Ext Data Fig. 2 in this file.

## File format

FORTRAN `parap`-slot order, one `(value name)` pair per line. Identical
in format to `wt-tribosphenic-2014.txt`.
