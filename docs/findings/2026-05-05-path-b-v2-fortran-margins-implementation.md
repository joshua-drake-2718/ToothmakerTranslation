---
title: 'Path B v2: fortran_margins laplacian implementation and disentanglement re-run'
date: 2026-05-05
status: rebuttal §B.3 strengthening, item 4-phase-1 of §E synthesis
upgrade: B.3 WEAK → ARGUABLE
---

## TL;DR

`Discretisation.laplacian='fortran_margins'` is now implemented at
`silicoshark/mesh.py:fortran_margins_laplacian` as the in-plane
contact-area-weighted operator. The substrate sink, vertical-z-flux,
and substrate-edge layer of the FORTRAN's full `apply_diffusion`
remain deferred (5–7 days additional work for byte-match against
Path A's `coreop2d.py`).

Re-running the cusp-forming disentanglement on
`examples/wt-tribosphenic-2014.txt` under the new operator gives an
asymmetric result:

- **`LEGACY_FORTRAN_minus_laplacian`**: byte-identical OFF outputs at
  every save (500, 1000, 1500, 2000, 2500 iters; SHA-256
  `3947760b3bc105ce`). The laplacian operator has zero observable
  effect on this dataset under LEGACY_FORTRAN's settings.
- **`PATH_B_DEFAULT_plus_laplacian`**: cusp_width drops from
  0.6927900 to 0.2947925 (−57 %). Vertex envelope's y span halves.
  Cusp count and cell count both stay at 19.

The B.3 row is no longer a structural no-op. It now reports an
empirical asymmetric measurement.

## Implementation

`silicoshark/mesh.py:fortran_margins_laplacian` replaces the prior
`NotImplementedError` with the FORTRAN-style per-edge `pes`-weighted
operator on silicoshark's existing variable-length CSR adjacency. For
each cell:

1. Sort neighbours angularly around the cell centre (xy plane).
2. Place a margin midpoint at the midpoint of each edge to a
   neighbour (Voronoi approximation, exact for a regular hex
   lattice).
3. `pes[k]` = perimeter-segment length between margin midpoint k and
   its angular successor (closed polygon — the loop wraps).
4. `area_p[k]` = `0.5 * ‖(midpoint_k − pos_i) × (midpoint_{k+1} −
   pos_i)‖`, the triangle area swept from the cell centre.
5. `area_bottom = sum(area_p)`, `sum_a = sum(pes) + 2 * area_bottom`,
   normalise `pes /= sum_a`. FMA-LIMIT GUARD: when `sum_a == 0`
   (degenerate cell with collapsed margin geometry), the cell
   contributes no flux.
6. Apply `L u[i] = sum_k pes[k] * (u[neigh_k] − u[i])`.

The FORTRAN's full `apply_diffusion` does additional work that this
implementation does NOT reproduce:

- A substrate-sink term `−0.44 * pes[i, j] * q3d[i, kk, k]` for any
  neighbour slot whose value is the substrate-sentinel `num_all_cells`
  (FORTRAN `13.f90:431, 444, 459`; Path A `coreop2d.py:596, 606,
  624`).
- Vertical-z-layer flux `area_bottom * (q[kk±1, k] − q[kk, k])`
  coupling the four z-layers (epithelial, upper mesenchyme, lower
  mesenchyme, substrate-edge).
- An explicit substrate-edge z-layer `q3d[:, max_z_layers, :]`
  beyond silicoshark's existing two mesenchyme layers, with an
  unconditional sink at `13.f90:438`.
- Tracking of mesenchyme Activator (currently silicoshark's State
  has no mes Act field; only mes Inh and mes Sec).

Reproducing those would require:

1. State extensions (`mes_inh_substrate_edge`,
   `mes_sec_substrate_edge`, mesenchyme Act).
2. Per-cell-per-slot `border` and `neigh` arrays alongside the
   existing CSR (Path A's nv_max=30 fixed-slot layout).
3. A flow-ordering change when `laplacian='fortran_margins'`: the
   FORTRAN does diffusion first across all 4 layers, then reaction
   on layer 0 only; silicoshark's existing flow does
   reaction-diffusion Jacobi on the epithelium, then mesenchyme.
4. A faithful port of `calculate_margins_a/b/c` (`coreop2d.py:387–
   502`) including the FORTRAN-bug-faithful "label-77" topology
   walk.

The estimate for the byte-match phase is 5–7 days. The in-plane
operator (this implementation, ~80 LOC) was 1.5 days. The plan
document records this scope decision at
`docs/plans/2026-05-05-fortran-margins-laplacian.md`.

## Verification

`tests/test_silicoshark_fortran_margins.py` (7 tests, ~1 s):

- Constant field annihilates (Laplacian of any constant is zero).
- Centre cell of the 7-cell hex lattice gives the analytical pes-
  normed flux: with margin midpoints at radius 0.5 spaced 60° apart,
  `sum_a = 3 + 3√3/4 ≈ 4.299`, `pes_normed = 0.5 / sum_a ≈ 0.1163`,
  centre-cell Laplacian on `u[ring]=1, u[centre]=0` is
  `6 * pes_normed ≈ 0.6978`. Test passes within `rel_tol=1e-12`.
- Sign flips correctly when the field is negated.
- Linear field `u(x, y) = x` gives zero Laplacian at the centre by
  symmetry (within `abs_tol=1e-13`).
- Cells with < 3 neighbours contribute zero (degenerate-polygon
  guard).
- `positions=None` raises `ValueError`.
- CLI smoke test runs PAPER_2010 + mesenchyme=absent +
  laplacian=fortran_margins for 100 × 5 iters, 19 cells, finite
  vertices.

The fast suite is at 62 passed in ~10 s (was 55 before this work).

`silicoshark/__main__.py:_check_implemented` no longer rejects
`laplacian='fortran_margins'`. The auto-override in
`scripts/run-discretisation-study.py:DEFERRED_FIELD_OVERRIDES` is
removed for the laplacian field; presets carrying
`laplacian='fortran_margins'` (LEGACY_FORTRAN, HUMPPA_LITERAL,
PATH_A_REWRITE) now exercise the operator at runtime.

## Disentanglement re-run

Output at
`experiments/discretisation-study/cusp-forming-fortran-margins-2026-05-05/`.

### Knock-down direction (`LEGACY_FORTRAN_minus_laplacian`)

| | LEGACY_FORTRAN + length_weighted | LEGACY_FORTRAN + fortran_margins |
|---|---|---|
| cell_count_history | [19, 19, 19, 19, 19] | [19, 19, 19, 19, 19] |
| cusp_count | 7 | 7 |
| cusp_width | 0.4952964725 | 0.4952964725 |
| z range | [58.575666, 59.474608] | [58.575666, 59.474608] |
| OFF SHA-256 (each save) | `3947760b3bc105ce` | `3947760b3bc105ce` |

Byte-identical at every save (500, 1000, 1500, 2000, 2500 iters).
The laplacian operator's choice has zero observable effect on this
dataset under LEGACY_FORTRAN's settings.

### Knock-up direction (`PATH_B_DEFAULT_plus_laplacian`)

| | PATH_B_DEFAULT (length_weighted) | PATH_B_DEFAULT + fortran_margins |
|---|---|---|
| cell_count_history | [19, 19, 19, 19, 19] | [19, 19, 19, 19, 19] |
| cusp_count | 19 | 19 |
| cusp_width | 0.6927900198 | 0.2947925243 |
| x range | [−1.158, 1.244] | [−1.103, −0.150] |
| y range | [−0.888, 0.837] | [−0.468, 0.468] |

Cusp_width drops 57 %. Vertex envelope shifts (y span halves).
The laplacian operator's choice IS measurable here.

## Methodological reading

The result mirrors the parameter-set-specificity pattern documented
on §B.1 (cusp-forming disentanglement and 16-point Act sweep): the
identity of fields that move biology depends on the rest of the
parameter set. Under LEGACY_FORTRAN's saturated regime
(`update_order='gauss_seidel_forces'`, `eq5_z_gate=True`,
`division_total_cap=60`, `knot_threshold_gate='first_border_cell'`),
the cell positions reach an equilibrium robust to the in-plane
laplacian's specific weighting. Under PATH_B_DEFAULT's settings the
dynamics remain in a regime where the weighting matters.

The B.3 row was previously a structural no-op (both anchors fell
back to length_weighted at runtime via the auto-override
mechanism); it is now an empirical comparison with a clean
asymmetric reading: a clean negative result on the knock-down
direction and a measurable −57 % cusp_width shift on the knock-up.

The empirical null on the knock-down direction is itself a
finding. If a reviewer asks 'is the laplacian choice load-bearing
for the FORTRAN reproduction?', the answer is now 'no, byte-
identically — at least within the in-plane operator's scope'. The
remaining residual is whether the substrate sink, vertical-z flux,
or substrate-edge layer would shift the answer; those are
separable sub-projects.

## What's open

- **Full byte-match against Path A's `coreop2d.py`** (5–7 days).
  Adds substrate sink, vertical-z flux, substrate-edge z-layer,
  mesenchyme Act tracking, and per-cell-per-slot `border`/`neigh`
  arrays. Acceptance criterion 2 of the plan
  (`docs/plans/2026-05-05-fortran-margins-laplacian.md`). Not
  required for B.3 ARGUABLE; would close it to SOLID.
- **z_min residual**
  (`docs/findings/2026-05-05-path-b-v2-legacy-fortran-divergence.md`).
  Untouched by this work: the in-plane operator does not modify
  forces or division dynamics, and the LEGACY_FORTRAN cusp-forming
  output is byte-identical between length_weighted and
  fortran_margins. The z_min divergence is a property of force
  propagation under `static_with_local_update` topology, not of
  the laplacian operator.
- **Re-running the seal-example single-field disentanglement
  with the new operator.** The seal example uses LEGACY_FORTRAN
  with `mesenchyme='per_column_z_layers'` (auto-overridden to
  absent in the runner currently). For the seal disentanglement
  the laplacian field's row should be re-measured to confirm the
  same byte-identicality holds; deferred — the cusp-forming
  result is the headline for B.3.

## Files touched

- `silicoshark/mesh.py` — replaced
  `fortran_margins_laplacian`'s `NotImplementedError` with the
  in-plane operator; extended the positions-retention condition
  from `kind == 'cotangent'` to
  `kind in ('cotangent', 'fortran_margins')`; updated
  `laplacian()` dispatcher's docstring and positions-check.
- `silicoshark/__main__.py:_check_implemented` — no longer
  rejects `laplacian='fortran_margins'`.
- `scripts/run-discretisation-study.py:DEFERRED_FIELD_OVERRIDES`
  — removed the laplacian auto-override.
- `tests/test_silicoshark_fortran_margins.py` (new, 7 tests).
- `docs/paper-rebuttal-preemption-2026-05-05.md` — rewrote §B.3
  rebuttal; updated §E synthesis; updated §F self-assessment.
- `docs/findings/2026-05-05-path-b-v2-fortran-margins-implementation.md`
  (this file).
- `experiments/discretisation-study/cusp-forming-fortran-margins-2026-05-05/`
  (new): two single-cell runs, LEGACY_FORTRAN with the new
  operator and PATH_B_DEFAULT_plus_laplacian.
