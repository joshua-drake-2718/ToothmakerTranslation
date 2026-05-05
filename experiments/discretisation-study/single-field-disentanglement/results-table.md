# Path B v2 single-field disentanglement

Generated: 2026-05-05T11:57:12.671500+00:00


For each field where `LEGACY_FORTRAN` and `PATH_B_DEFAULT` differ, this study runs a single-field perturbation:

- *Knock-down*: `LEGACY_FORTRAN` with that field reset to its `PATH_B_DEFAULT` value. If the plateau collapses toward `PATH_B_DEFAULT`, the field is doing the work in the FORTRAN-flavoured bundle.
- *Knock-up*: `PATH_B_DEFAULT` with that field set to its `LEGACY_FORTRAN` value. If the plateau lifts toward `LEGACY_FORTRAN`, that field alone is sufficient.


## Knock-down: LEGACY_FORTRAN minus each field


### `seal`

| Preset | Plateau | Cusps | x range | y range | z range | Regime |
|---|---:|---:|---|---|---|---|
| `LEGACY_FORTRAN` (anchor) | 60.00 | 0 | [-1.49, 1.73] | [-2.13, 1.93] | [455.31, 609.25] | plateau |
| `LEGACY_FORTRAN_minus_adh_form` | 40.00 | 0 | [-1.56, 1.28] | [-1.86, 2.03] | [115.21, 118.74] | plateau |
| `LEGACY_FORTRAN_minus_border_definition` | 60.00 | 0 | [-0.79, 0.83] | [-1.34, 1.28] | [384.90, 494.38] | plateau |
| `LEGACY_FORTRAN_minus_division_total_cap` | NaN | — | — | — | — | NaN |
| `LEGACY_FORTRAN_minus_eq17_inh_source` | 60.00 | 0 | [-1.49, 1.73] | [-2.13, 1.93] | [455.31, 609.25] | plateau |
| `LEGACY_FORTRAN_minus_eq18_sec_source` | 60.00 | 0 | [-1.49, 1.73] | [-2.13, 1.93] | [455.31, 609.25] | plateau |
| `LEGACY_FORTRAN_minus_eq5_apply_to` | 60.00 | 0 | [-1.41, 1.40] | [-1.32, 0.38] | [849.48, 1084.14] | plateau |
| `LEGACY_FORTRAN_minus_knot_daughter_di` | 60.00 | 0 | [-1.49, 1.73] | [-2.13, 1.93] | [455.31, 609.25] | plateau |
| `LEGACY_FORTRAN_minus_knot_threshold_gate` | NaN | — | — | — | — | NaN |
| `LEGACY_FORTRAN_minus_laplacian` | 60.00 | 0 | [-1.49, 1.73] | [-2.13, 1.93] | [455.31, 609.25] | plateau |
| `LEGACY_FORTRAN_minus_lattice_orientation` | 60.00 | 0 | [-1.30, 1.30] | [-1.61, 1.63] | [474.32, 627.42] | plateau |
| `LEGACY_FORTRAN_minus_rep_form` | 37.00 | 0 | [-1.28, 1.28] | [-2.55, 2.55] | [160.14, 162.96] | plateau |
| `LEGACY_FORTRAN_minus_rep_neighbour_set` | 60.00 | 0 | [-1.49, 1.75] | [-2.13, 1.92] | [455.19, 609.26] | plateau |
| `LEGACY_FORTRAN_minus_update_order` | 60.00 | 0 | [-2.60, 2.60] | [-1.95, 1.95] | [292.16, 387.35] | plateau |

Interpretation:

- Removing field `adh_form` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 40.00 (delta -20.00; 87% of the 23.00-cell span between anchors) — field `adh_form` is therefore **strongly consequential**.
- Removing field `border_definition` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 60.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `border_definition` is therefore **not consequential**.
- Field `division_total_cap`: run failed or produced NaN plateau; effect cannot be assessed.
- Removing field `eq17_inh_source` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 60.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `eq17_inh_source` is therefore **not consequential**.
- Removing field `eq18_sec_source` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 60.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `eq18_sec_source` is therefore **not consequential**.
- Removing field `eq5_apply_to` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 60.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `eq5_apply_to` is therefore **not consequential**.
- Removing field `knot_daughter_di` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 60.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `knot_daughter_di` is therefore **not consequential**.
- Field `knot_threshold_gate`: run failed or produced NaN plateau; effect cannot be assessed.
- Removing field `laplacian` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 60.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `laplacian` is therefore **not consequential**.
- Removing field `lattice_orientation` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 60.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `lattice_orientation` is therefore **not consequential**.
- Removing field `rep_form` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 37.00 (delta -23.00; 100% of the 23.00-cell span between anchors) — field `rep_form` is therefore **strongly consequential**.
- Removing field `rep_neighbour_set` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 60.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `rep_neighbour_set` is therefore **not consequential**.
- Removing field `update_order` from `LEGACY_FORTRAN` shifts the plateau from 60.00 to 60.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `update_order` is therefore **not consequential**.


## Knock-up: PATH_B_DEFAULT plus each field


### `seal`

| Preset | Plateau | Cusps | x range | y range | z range | Regime |
|---|---:|---:|---|---|---|---|
| `PATH_B_DEFAULT` (anchor) | 37.00 | 0 | [-3.00, 3.00] | [-0.32, 0.32] | [138.85, 141.38] | plateau |
| `PATH_B_DEFAULT_plus_adh_form` | 37.00 | 0 | [-3.00, 3.00] | [-1.11, 1.11] | [223.13, 225.88] | plateau |
| `PATH_B_DEFAULT_plus_border_definition` | 37.00 | 0 | [-3.00, 3.00] | [-0.39, 0.39] | [137.79, 140.32] | plateau |
| `PATH_B_DEFAULT_plus_division_total_cap` | 37.00 | 0 | [-3.00, 3.00] | [-0.32, 0.32] | [138.85, 141.38] | plateau |
| `PATH_B_DEFAULT_plus_eq17_inh_source` | 37.00 | 0 | [-3.00, 3.00] | [-0.32, 0.32] | [138.85, 141.38] | plateau |
| `PATH_B_DEFAULT_plus_eq18_sec_source` | 37.00 | 0 | [-3.00, 3.00] | [-0.32, 0.32] | [138.85, 141.38] | plateau |
| `PATH_B_DEFAULT_plus_eq5_apply_to` | 37.00 | 0 | [-3.00, 3.00] | [-0.26, 0.26] | [128.55, 131.13] | plateau |
| `PATH_B_DEFAULT_plus_knot_daughter_di` | 37.00 | 0 | [-3.00, 3.00] | [-0.32, 0.32] | [138.85, 141.38] | plateau |
| `PATH_B_DEFAULT_plus_knot_threshold_gate` | 37.00 | 0 | [-3.00, 3.00] | [-0.32, 0.32] | [138.85, 141.38] | plateau |
| `PATH_B_DEFAULT_plus_laplacian` | 37.00 | 0 | [-3.00, 3.00] | [-0.32, 0.32] | [138.85, 141.38] | plateau |
| `PATH_B_DEFAULT_plus_lattice_orientation` | 37.00 | 0 | [-2.60, 2.60] | [-0.86, 0.86] | [120.14, 122.24] | plateau |
| `PATH_B_DEFAULT_plus_rep_form` | NaN | — | — | — | — | NaN |
| `PATH_B_DEFAULT_plus_rep_neighbour_set` | 37.00 | 0 | [-3.00, 3.00] | [-0.54, 0.54] | [138.52, 140.96] | plateau |
| `PATH_B_DEFAULT_plus_update_order` | 37.00 | 0 | [-1.10, 1.18] | [-1.50, 1.46] | [149.56, 152.25] | plateau |

Interpretation:

- Adding field `adh_form` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `adh_form` is therefore **not consequential**.
- Adding field `border_definition` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `border_definition` is therefore **not consequential**.
- Adding field `division_total_cap` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `division_total_cap` is therefore **not consequential**.
- Adding field `eq17_inh_source` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `eq17_inh_source` is therefore **not consequential**.
- Adding field `eq18_sec_source` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `eq18_sec_source` is therefore **not consequential**.
- Adding field `eq5_apply_to` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `eq5_apply_to` is therefore **not consequential**.
- Adding field `knot_daughter_di` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `knot_daughter_di` is therefore **not consequential**.
- Adding field `knot_threshold_gate` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `knot_threshold_gate` is therefore **not consequential**.
- Adding field `laplacian` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `laplacian` is therefore **not consequential**.
- Adding field `lattice_orientation` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `lattice_orientation` is therefore **not consequential**.
- Field `rep_form`: run failed or produced NaN plateau; effect cannot be assessed.
- Adding field `rep_neighbour_set` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `rep_neighbour_set` is therefore **not consequential**.
- Adding field `update_order` to `PATH_B_DEFAULT` shifts the plateau from 37.00 to 37.00 (delta +0.00; 0% of the 23.00-cell span between anchors) — field `update_order` is therefore **not consequential**.

