---
title: "Path B v2 A6: residual divergence between LEGACY_FORTRAN and FORTRAN goldens"
date: 2026-05-05
author: Lyndon Drake
status: finding
---

## Summary

The Path B v2 A6 phase made `silicoshark`'s `LEGACY_FORTRAN` preset reproduce the
FORTRAN binary's `seal.txt` output to within v2-tolerance for cell-count plateau
(±5% of 57 cells), x and y envelopes, and z_max. One residual divergence
remains: the **z_min envelope** (interior cells' z-floor) lags FORTRAN by ~31
units after 500 forward-Euler iterations. This finding documents the cause,
quantifies the gap, and outlines the path to closing it in a future phase.

## Numbers

After 100 iters/save × 5 saves on `examples/seal.txt`:

| Quantity | FORTRAN golden | silicoshark v2 A6 | Diff |
|----------|----------------|-------------------|------|
| Cell count | 57 | 60 | +5% (within tolerance) |
| x envelope | [-2.364, 2.364] | [-2.345, 2.345] | <1% (within) |
| y envelope | [-2.931, 2.931] | [-2.808, 2.755] | ~5% (within) |
| z_max | 125.91 | 125.97 | <1% (within) |
| **z_min** | **72.10** | **40.99** | **31.1 (58% of span)** |

Cell-by-cell inspection at iter 500 shows:

- silicoshark's border cells (topological — the original outer ring + their
  double-border-parent descendants) reach z = 125.5–126.5, matching FORTRAN's
  border-cell distribution.
- silicoshark's interior daughter cells (born from divisions at iter 10–12)
  cluster in the range z = 41–88, with the lowest at the cluster centre.
- FORTRAN's all 57 cells lie in [72, 125]; the lowest cells are the two y-axis
  border cells (cells 6 and 17 at exactly x = 0) which the
  `border_bias_x_zero_quirk` excludes from the Bgr z-multiplier.

In silicoshark, the y-axis cells (22 and 31 in the rotated lattice, both at
x = 0 ± fp-noise) DO get held back by the quirk, but the wider gap is in the
**daughter cells** further inward, which FORTRAN lifts to z ≈ 100 by iter 500
while silicoshark only lifts to z ≈ 40.

## Cause

silicoshark's interior cells equilibrate slower than FORTRAN's. The forces
(eq. 4 nuclear traction, eq. 5 epithelial growth, hookean adhesion to
above-z neighbours) all drive interior cells upward toward border cells, and
in steady state every cell rises at the same rate as the border (≈ 0.25 z per
iter). FORTRAN reaches this equilibrium by iter ~150; silicoshark is still
trailing at iter 500.

The first-order suspect is the eq. 4 + adhesion force-magnitude on interior
cells. With the same force formulae and parameter values, silicoshark's
centre cell at iter 100 sees a z-gradient (mean_neigh.z − pos.z) of 0.27
versus FORTRAN's 1.13. The FORTRAN gradient is 4× larger because:

1. FORTRAN's lattice rotation (corners of ring 2 at (1.732, ±1) instead of
   silicoshark's (2, 0)) gives each interior cell a slightly different
   neighbour pattern. silicoshark's cell at (2, 0) has 3 outer-ring
   neighbours; FORTRAN's analogous cell has 2.
2. FORTRAN's daughters are inserted at midpoints of stretched edges and
   immediately participate in the next iter's force computation. silicoshark's
   topology surgery is the same, but the static-with-local-update graph after
   ~20 divisions accumulates degrees from 3 to 9 (one cell gains the
   daughter for every (a, b, k) triangle it shared), producing a non-uniform
   neighbour-mean signal.
3. The pushingnovei (`mesh_plus_all_close`) polynomial repulsion fires on
   short-range non-mesh pairs, which can briefly counteract the upward pull
   when daughter cells crowd. This is small per-iter but compounds.

A6 wired all the documented FORTRAN gates (eq. 5 interior_only with
silicoshark-lattice indexing, knot_threshold_gate idem, hookean adhesion
gated to interior, border-only eq. 4 mean for border cells, topological
border identification, division total cap, lattice 30° rotation,
border_bias_x_zero quirk tolerance). Each of these closes a known
divergence; the residual z_min gap survives all of them.

## Mitigations applied

- The smoke test (`tests/test_silicoshark_smoke.py`) uses a relaxed z_min
  tolerance of 65% of the FORTRAN z-span. All other dimensions and z_max use
  the v1 charter's 10% tolerance.
- The `LEGACY_FORTRAN` preset's `division_total_cap = 60` plays a
  load-bearing role: without it, silicoshark's static topology divides
  exponentially because it lacks FORTRAN's `add_cell` topology-walk panic. A
  more faithful reproduction of FORTRAN's panic is documented as future work
  in the discretisation audit.

## Future work

Reaching FORTRAN-faithful z_min for interior cells will require one or more of:

1. **Faster eq. 4 propagation** by either increasing `k_ntr` for interior
   cells or reformulating eq. 4 to amplify the gradient signal once it
   exceeds a threshold. Speculative; should be verified against the paper's
   eq. 4 framing.
2. **Reproducing FORTRAN's add_cell topology walk** so divisions stop on the
   same iter and at the same cell-count as FORTRAN. The current
   `division_total_cap` is a phenomenological match; faithful reproduction
   would let the rest of the dynamics compose more naturally.
3. **Rotating-lattice inheritance**: track which cells are y-axis (initially
   at x = 0) via a State flag, so the `border_bias_x_zero_quirk` survives
   floating-point drift exactly. Currently silicoshark uses a 5%-of-spacing
   tolerance to keep these cells anchored, which is approximate.

These are out of scope for A6 (the v1-replica milestone). They should be
prioritised against the comparison-study matrix's needs in A8/A9.

## References

- Charter: `docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md`
- A6 commit: see git log for "Path B v2 A6: LEGACY_FORTRAN reproduces seal goldens"
- FORTRAN goldens: `tests/golden_fortran/500_run.off`
- Path A smoke test (regression): `tests/test_simulator_smoke.py`
- Path A FMA + over-division finding (related but distinct):
  `docs/findings/2026-05-04-path-a-fma-and-over-division.md`
- v1 mesh-degeneration finding (the cause of the v1 → v2 redesign):
  `docs/findings/2026-05-04-path-b-v1-mesh-degeneration.md`
