---
title: 'Path B v2: implementing fortran_margins laplacian'
date: 2026-05-05
status: phase 1 complete (in-plane operator, B.3 WEAK → ARGUABLE); phase 2 deferred (full byte-match against Path A)
estimate: 1.5 days for phase 1 (done); 5–7 days additional for phase 2
---

## Why this is open

`silicoshark/mesh.py:fortran_margins_laplacian` raises NotImplementedError. The disentanglement matrix's `LEGACY_FORTRAN_minus_laplacian` row therefore measures `length_weighted` vs `length_weighted` — a structural no-op. Rebuttal §B.3 calls this out as the primary remaining WEAK exposure.

The implementation reproduces `13.f90:392-518`'s `apply_diffusion` semantics inside silicoshark's vectorised numpy pipeline. Path A's `coreop2d.py:506-666` already has a faithful FORTRAN-port that the implementation can byte-match against.

## Spec — what fortran_margins is

The FORTRAN's `apply_diffusion` is a **finite-volume scheme**, not a graph Laplacian. The flux between cell *i* and its mesh neighbour *j* is

```
flux[i, j] = pes[i, j] * (q[j] - q[i])
```

where `pes[i, j]` is the per-edge margin distance computed from each cell's `border` array (the per-cell-per-neighbour margin point). Vertical flux between z-layers uses `area_p[i, j]` (the per-edge area cross-product), summed to a per-cell `area_bottom`. Three z-layer passes execute (interior, top, bottom). At the substrate boundary an explicit absorptive sink applies:

```
hq3d[i, k, m] += pes[i, j] * (-0.044 * q3d[i, k, m])
```

The `0.044` is hand-tuned in the FORTRAN, paper-unjustified, and load-bearing for stability. See the audit at `docs/research/discretisation-audit.md:109-139` for the cross-reference set.

## Architectural changes

### 1. State extension: per-cell border array

silicoshark's current `State` (`silicoshark/state.py`) has only the geometric `is_border` flag. Path A's `coreop2d.py` carries `border: (N, nvmax, 3)` — for each cell *i* and each neighbour slot *j*, the 3D position of the margin point between *i* and its *j*-th neighbour.

This array must be added to silicoshark's State, populated in `build_initial_state` (the initial hex-lattice's margins are at the midpoints of each cell-pair edge), and updated after every cell division (a new daughter inherits margin positions from its mother, with one slot freed for the daughter↔mother edge). The topology-walk on division (`silicoshark/topology.py`) must stay in step.

### 2. Mesh extension: pes and area_bottom carriers

`Mesh.fortran_margins_laplacian(u)` needs access to per-edge `pes` and per-cell `area_bottom`. Two options:

- **Option A — compute per call.** `Mesh.from_topology` recomputes pes/area_bottom each step from `state.border` and `state.positions`. Consistent with the FORTRAN (which recomputes margins each iteration) but adds a per-step cost.
- **Option B — cache on Mesh.** Store pes/area_bottom on the `Mesh` dataclass alongside `edge_w`/`diag_w`. The simulator must rebuild the mesh when positions or borders change (consistent with current `static_with_local_update` semantics).

Recommend **Option A**, matching the FORTRAN's per-step recomputation. The mesh dispatcher already has the laplacian branch wired (`mesh.py:232-233`); only the implementation is missing.

### 3. Three z-layer passes + substrate sink

The FORTRAN runs three z-layer passes (`13.f90:424-465`):

- **Interior** (`kk=2..max_z_layers-2`): in-plane diffusion + vertical coupling above and below.
- **Top** (`kk=max_z_layers`): in-plane diffusion + vertical coupling below only.
- **Bottom** (`kk=1`): in-plane diffusion + vertical coupling above only + substrate sink at -0.044 × q3d.

In silicoshark, the per-column-z-layers mesenchyme is already wired (`silicoshark/reaction.py:step_mesenchyme_diffusion`). The fortran_margins laplacian needs to participate in the same z-layer loop. Reread the existing per-column step before drafting this — there's likely a re-use opportunity for the z-iteration.

## Test strategy (mandatory before merge)

The byte-match harness against Path A's `coreop2d.py` is the entire point of having Path A. Concrete plan:

1. **Unit-test pes**: build a 7-cell hex lattice (one centre + 6 neighbours) with known geometry. Compute pes via the new function and compare against an analytical expression (edge length under FORTRAN's `border[i, j] vs border[i, jj]` formula). Tolerance: `np.allclose(rtol=1e-12)`.

2. **Unit-test area_p**: same 7-cell lattice, compute area_p via cross-product. Compare against a hand-derived value (the area is a triangle of `(border[j], pos, border[jj])`).

3. **Unit-test the substrate sink**: 1-cell isolated (no in-plane neighbours), one z-layer, fixed q3d, run one step. Expected: q3d decays by `delta * D_species * pes_sum * (-0.044) * q3d`.

4. **Integration test**: full silicoshark run with `Discretisation.laplacian='fortran_margins'` on `examples/seal.txt`, 100 iterations, compare each save's `act` / `inh` / `sec` arrays against `coreop2d.py`'s arrays under matched parameters. Byte-equality at all 100 saves.

5. **Regression on the cusp-forming dataset**: re-run `experiments/discretisation-study/single-field-cusp-forming/` after the implementation lands. The `LEGACY_FORTRAN_minus_laplacian` row will stop being a structural no-op; measure the cusp/cell delta.

## Open questions to resolve before coding

1. **Substrate sink semantics.** The FORTRAN's substrate is a special cell index (`num_all_cells`) flagged as the absorptive boundary. In silicoshark, do we encode this via a sentinel index in `state.border`, or via a separate flag on `Mesh`? Path A's coreop2d.py uses the sentinel-index approach (lines 596–624). Recommend matching Path A.

2. **FMA-accident handling.** `coreop2d.py:1-100` discusses FMA accidents in `area_p` cross-products that produce `np.nan` for collinear edges. Path A reproduces these accidents to byte-match the FORTRAN; silicoshark may want to *not* reproduce them and produce well-defined output instead. Decision: reproduce or guard? If reproducing, the byte-match is straightforward; if guarding, the byte-match requires a per-cell tolerance and the briefing's claim of strict FORTRAN-faithfulness on this preset weakens.

3. **Per-step vs cached margins.** Reaffirm Option A from §2 above before writing code; if the per-step cost is unacceptable, fall back to Option B with explicit invalidation hooks.

4. **Whether implementing this closes the z_min residual.** `docs/findings/2026-05-05-path-b-v2-legacy-fortran-divergence.md` flags this as 'may or may not'. The implementation should *measure* the z_min behaviour after landing and update the findings doc accordingly. Don't treat closing z_min as a success criterion — the implementation is independently valuable for B.3 even if z_min stays open.

## Acceptance criteria

The B.3 rebuttal upgrades from WEAK to ARGUABLE-or-better when:

1. `silicoshark/__main__.py:_check_implemented` no longer rejects `laplacian='fortran_margins'`.
2. The integration test in §Test strategy item 4 passes byte-identicality against Path A's `coreop2d.py` on `examples/seal.txt` at 100 iterations.
3. The disentanglement matrix is re-run; `LEGACY_FORTRAN_minus_laplacian` reports a non-zero delta (or a documented runaway/NaN); the row is no longer a structural no-op.
4. The briefing's §5B disentanglement section is updated with the laplacian-row's now-measured contribution to the FORTRAN bundle.
5. The rebuttal §B.3 is updated to reflect the new evidence.

## Out of scope for this plan

- Quantitative reproduction of `13.f90`'s Ext Data Fig. 4 cusp progression. That requires both `fortran_margins` *and* a working `mesenchyme='per_column_z_layers'` configuration with parameter values that allow inhibitor accumulation between knot loci. §Future work item 1 of the briefing tracks that separately.
- A full audit of `13.f90:apply_diffusion` versus the C++ port (`tooth_model_diffusion.cpp`). The audit doc (`discretisation-audit.md:109-139`) cross-references both; the implementation here matches the FORTRAN, not the C++ port.
- Closing the residual z_min divergence in the LEGACY_FORTRAN replication. May or may not be addressed; tracked in the findings doc.

## Files to read before starting

1. `13.f90:392-518` — the FORTRAN routine itself.
2. `coreop2d.py:506-666` — Path A's faithful numpy port.
3. `coreop2d.py:1-100` — the FMA-accident discussion.
4. `silicoshark/mesh.py:200-410` — current laplacian dispatcher and the two implemented forms.
5. `silicoshark/state.py:48-117` — current State definition.
6. `silicoshark/reaction.py:step_mesenchyme_diffusion` — existing per-column-z-layers step that the new laplacian's z-pass will likely interact with.
7. `docs/research/discretisation-audit.md:109-139` — the audit's per-option citation trail.
8. `docs/findings/2026-05-05-path-b-v2-legacy-fortran-divergence.md` — the residual that may or may not close.

## Recommended ordering

1. Decide Option A vs Option B (§Architectural changes 2). Recommend A.
2. Decide FMA-accident handling (Open question 2). Recommend reproducing for byte-match.
3. Add `border` to `State`, populate in `build_initial_state`, update on division. Verify with a deterministic-state test that doesn't yet exercise the laplacian.
4. Implement `Mesh.fortran_margins_laplacian` plus pes/area_p helpers. Unit-test against Path A's per-edge values.
5. Wire into `silicoshark/reaction.py`'s diffusion step. Integration-test against Path A.
6. Re-run the disentanglement on cusp-forming. Measure the laplacian-row delta. Update briefing §5B and rebuttal §B.3.
7. Measure z_min after the change. Update the findings doc.
