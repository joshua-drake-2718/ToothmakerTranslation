---
title: 'Path B v1 attempt — Delaunay-driven mesh degeneration on the seal example'
author: Lyndon Drake (with Claude Code)
date: 2026-05-04
---

## Context

First-pass `silicoshark/` implementation. Phases 1 (scaffold + I/O), 2
(triangulation + reaction-diffusion), and most of 3 (mechanical
forces) and 4 (cell division) are written; Phases 5 (border biases)
and 6 (validation) are stubbed. The simulator runs end-to-end on
`examples/seal.txt` for a few hundred iterations, but degenerates
before reaching the 500-iter golden checkpoint.

## What works

- Hex-lattice initial conditions, parameterised by `rad` (matches
  `13.f90`'s `num_active_cells = 3*(rad-1)*rad + 1`).
- Length-weighted graph Laplacian on `scipy.spatial.Delaunay`-derived
  edges. Lap(1) = 0, Lap(x) = 0 on interior cells (planar
  harmonic), Lap(bump) negative at peak. Stable for `dt = 0.05`.
- Reaction-diffusion (eqs. 14–18) with paper-faithful denominators
  and knot/threshold overrides. Does the right thing on the seal
  parameters: Act decays from 0.5 toward 0, no knots form (matches
  the FORTRAN-faithful behaviour on this input).
- Repulsion (eq. 1, Hookean), adhesion (eq. 2, additive
  unit-vector pull), nucleus traction (eq. 4), epithelial growth
  (eq. 5, z-gated to avoid floating-point amplification at flat
  initial conditions), cervical-loop downgrowth (eqs. 7–9 with
  mesenchyme-thickening eqs. 10–12).
- Pbi/Abi/Bgr border-bias multipliers and z-clamp from the
  FORTRAN's update_cell_position, including the x = 0 quirk.
- Cell division by edge-length threshold, with knot-daughter
  d_i hard reset per the paper, and per-step caps to dampen
  cascade behaviour from Delaunay re-triangulation.

## What does not work

The seal example's mechanical equilibrium is unstable in the v1
implementation. Cells contract in y over the first ~100 iterations
(initial y-range ±3 collapses to ±0.7 by iter 200), and after another
~70 iterations the (x, y) projection becomes near-coplanar enough that
`scipy.spatial.Delaunay` raises `QHullError` on the next
re-triangulation. Forces and positions briefly explode to ~1e88 just
before the failure.

## Root causes

1. **Eq. 5 cannot be unconditionally applied to the symmetric initial
   hex.** Without the FORTRAN's z-gate (`b < -1e-4`), perfectly
   symmetric neighbour configurations have `sum_unit ≈ 0` in
   floating-point but not exactly zero; the unit-vector normalisation
   amplifies this to a phantom inward force of ~k_egr per step. With
   the z-gate (now in `forces.py`), the initial step is force-free,
   but once cervical-loop curvature develops, eq. 5 fires
   asymmetrically and pulls cells inward in y.

2. **Adhesion fires from floating-point edge-length noise at
   equilibrium.** Initial unit edges have lengths in
   `{1.0 ± 1e-16}` due to floating-point arithmetic; a fraction
   register as `dist > rest`, triggering eq. 2's unit-vector pull.
   Fixed in v1 by making eq. 1 a fully Hookean spring (smoothly
   signed across the rest length) and gating eq. 2 with a tolerance
   band, but the underlying issue — cells starting at exactly the
   FORTRAN-pinned rest length — is fragile.

3. **Delaunay re-triangulation produces convex-hull edges that drive
   inappropriate divisions.** When border cells migrate outward and
   inner cells stay put, Delaunay's outer edges can grow to length
   ~3, triggering eq. 1 attraction and cell division. The new
   daughter then sits between cells that the FORTRAN's
   static-topology view would never have called neighbours,
   creating new long edges and a cascade. Capped per-step
   division (currently 4 per pass) and a `max_edge_dist` filter
   slow this down but do not fundamentally fix it.

## What the FORTRAN does differently

`13.f90` uses a static neighbour graph: the hex topology is built
once at initial conditions and updated only locally on cell division
(`add_cell` rewires the immediate parent-daughter links and the
nearby topology walk). Long-distance edges to migrated border cells
do not appear in the FORTRAN's neighbour table because the topology
is intrinsically local — cells that were never neighbours in the
initial hex never become neighbours later. This is what holds the
seal example at the 57-cell plateau over 500+ iterations: the
topology only sees ~6 neighbours per cell, division splits only those
6 edges, and the geometry stays bounded.

Path B's Delaunay-each-iteration approach trades static-topology
robustness for vectorisation simplicity. The charter
(`docs/plans/2026-05-04-path-b-charter.md`) explicitly endorses
Delaunay re-triangulation as a way to eliminate the FORTRAN's
topology-walk accidents (`txungu` else-branch, `iiii` preservation,
panic-and-return). The current v1 evidence suggests Delaunay
re-triangulation needs to be paired with either a tighter division
criterion or a static-topology fallback to be stable.

## Recommended next steps for Path B v2

In rough priority order:

1. **Static topology with local update on division (FORTRAN-style
   without the walk).** Build the neighbour graph once from
   `scipy.spatial.Delaunay` at t = 0; on cell division, locally
   rewire so that the new daughter inherits a subset of each parent's
   neighbours and the parents' direct edge is replaced by two edges
   through the daughter. No full re-triangulation. This eliminates
   the cascade source of v1.

2. **Investigate FORTRAN's epithelial-growth z-gate.** The gate's
   intent is unclear from the paper. Either match it (current v1)
   or replace with a curvature-based outward direction (e.g.,
   discrete mean-curvature normal at each vertex).

3. **Re-examine the eq. 5 / cervical-loop relationship.** The
   FORTRAN treats them as exclusive (one or the other on a cell);
   the paper treats them as additive. Decide and document.

4. **Reaction-diffusion mesenchyme.** Currently absent. The seal
   example does not exercise it (Sec stays 0), but other parameter
   sets will need vertical Inh/Sec diffusion between epithelial and
   mesenchymal layers (paper p. 586).

5. **Cotangent-weight Laplacian.** The current length-weighted
   graph Laplacian works for the v1 paper-vs-FORTRAN gap, but the
   charter recommendation is the cotangent operator for proper
   curvature handling.

## Files added in v1

- `silicoshark/__init__.py`
- `silicoshark/params.py` — paper-symbol-named parameter dataclass
- `silicoshark/state.py` — vectorised state arrays + hex_lattice +
  build_initial_state
- `silicoshark/mesh.py` — Delaunay triangulation, CSR neighbour list,
  length-weighted Laplacian
- `silicoshark/reaction.py` — eqs. 14–18 forward Euler
- `silicoshark/forces.py` — eqs. 1–13 mechanical forces with
  apply_border_multipliers
- `silicoshark/simulator.py` — `step()`, `run()`, `divide_cells()`
- `silicoshark/io.py` — OFF writer

The FORTRAN-faithful Path A test (`tests/test_simulator_smoke.py`)
still passes (1.0s, 57 vertices/236 faces).

No `tests/test_silicoshark_smoke.py` yet — Phase 6 is blocked on the
mesh-stability fix above.
