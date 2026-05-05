# Path B v2 single-field disentanglement

Generated: 2026-05-05T21:40:55.926945+00:00


For each field where `LEGACY_FORTRAN` and `PATH_B_DEFAULT` differ, this study runs a single-field perturbation:

- *Knock-down*: `LEGACY_FORTRAN` with that field reset to its `PATH_B_DEFAULT` value. If the plateau collapses toward `PATH_B_DEFAULT`, the field is doing the work in the FORTRAN-flavoured bundle.
- *Knock-up*: `PATH_B_DEFAULT` with that field set to its `LEGACY_FORTRAN` value. If the plateau lifts toward `LEGACY_FORTRAN`, that field alone is sufficient.


## Knock-down: LEGACY_FORTRAN minus each field


### `wt-tribosphenic-2014`

| Preset | Plateau | Cusps | x range | y range | z range | Regime |
|---|---:|---:|---|---|---|---|
| `LEGACY_FORTRAN` (anchor) | 19.00 | 7 | [2.55, 4.08] | [7.89, 9.02] | [58.58, 59.47] | plateau |
| `LEGACY_FORTRAN_minus_adh_form` | 19.00 | 7 | [2.46, 3.99] | [9.09, 10.23] | [57.19, 58.09] | plateau |
| `LEGACY_FORTRAN_minus_border_definition` | 19.00 | 7 | [1.92, 6.36] | [-0.11, -0.11] | [30.47, 31.66] | plateau |
| `LEGACY_FORTRAN_minus_division_total_cap` | 19.00 | 7 | [2.55, 4.08] | [7.89, 9.02] | [58.58, 59.47] | plateau |
| `LEGACY_FORTRAN_minus_eq17_inh_source` | 19.00 | 7 | [2.55, 4.08] | [7.89, 9.02] | [58.58, 59.47] | plateau |
| `LEGACY_FORTRAN_minus_eq18_sec_source` | 19.00 | 7 | [0.80, 5.75] | [-0.25, 0.20] | [104.72, 105.85] | plateau |
| `LEGACY_FORTRAN_minus_eq5_apply_to` | 19.00 | 7 | [2.36, 3.85] | [6.86, 7.98] | [61.44, 62.35] | plateau |
| `LEGACY_FORTRAN_minus_knot_daughter_di` | 19.00 | 7 | [2.55, 4.08] | [7.89, 9.02] | [58.58, 59.47] | plateau |
| `LEGACY_FORTRAN_minus_knot_threshold_gate` | 19.00 | 19 | [-0.08, 5.11] | [-0.60, 0.60] | [1.01, 1.01] | plateau |
| `LEGACY_FORTRAN_minus_laplacian` | 19.00 | 7 | [2.55, 4.08] | [7.89, 9.02] | [58.58, 59.47] | plateau |
| `LEGACY_FORTRAN_minus_lattice_orientation` | 19.00 | 7 | [2.15, 6.97] | [0.35, 0.35] | [44.44, 46.31] | plateau |
| `LEGACY_FORTRAN_minus_rep_form` | 19.00 | 7 | [0.82, 5.38] | [-0.29, 0.26] | [76.57, 77.43] | plateau |
| `LEGACY_FORTRAN_minus_rep_neighbour_set` | 19.00 | 7 | [0.71, 5.90] | [0.00, 0.00] | [70.11, 70.93] | plateau |
| `LEGACY_FORTRAN_minus_update_order` | 60.00 | 7 | [-0.30, 0.66] | [-0.51, 0.70] | [1.00, 9.18] | plateau |

Interpretation:

- Removing field `adh_form` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `adh_form` is therefore **not consequential**.
- Removing field `border_definition` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `border_definition` is therefore **not consequential**.
- Removing field `division_total_cap` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `division_total_cap` is therefore **not consequential**.
- Removing field `eq17_inh_source` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `eq17_inh_source` is therefore **not consequential**.
- Removing field `eq18_sec_source` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `eq18_sec_source` is therefore **not consequential**.
- Removing field `eq5_apply_to` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `eq5_apply_to` is therefore **not consequential**.
- Removing field `knot_daughter_di` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `knot_daughter_di` is therefore **not consequential**.
- Removing field `knot_threshold_gate` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `knot_threshold_gate` is therefore **not consequential**.
- Removing field `laplacian` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `laplacian` is therefore **not consequential**.
- Removing field `lattice_orientation` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `lattice_orientation` is therefore **not consequential**.
- Removing field `rep_form` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `rep_form` is therefore **not consequential**.
- Removing field `rep_neighbour_set` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `rep_neighbour_set` is therefore **not consequential**.
- Removing field `update_order` from `LEGACY_FORTRAN` shifts the plateau from 19.00 to 60.00 (delta +41.00; 0% of the 0.00-cell span between anchors) — field `update_order` is therefore **not consequential**.


## Knock-up: PATH_B_DEFAULT plus each field


### `wt-tribosphenic-2014`

| Preset | Plateau | Cusps | x range | y range | z range | Regime |
|---|---:|---:|---|---|---|---|
| `PATH_B_DEFAULT` (anchor) | 19.00 | 19 | [-1.16, 1.24] | [-0.89, 0.84] | [1.00, 1.01] | plateau |
| `PATH_B_DEFAULT_plus_adh_form` | 19.00 | 19 | [-1.00, 1.41] | [-0.85, 0.89] | [1.00, 1.01] | plateau |
| `PATH_B_DEFAULT_plus_border_definition` | 19.00 | 19 | [-5.06, -3.41] | [-1.00, 0.99] | [1.00, 1.02] | plateau |
| `PATH_B_DEFAULT_plus_division_total_cap` | 19.00 | 19 | [-1.16, 1.24] | [-0.89, 0.84] | [1.00, 1.01] | plateau |
| `PATH_B_DEFAULT_plus_eq17_inh_source` | 19.00 | 19 | [-1.16, 1.24] | [-0.89, 0.84] | [1.00, 1.01] | plateau |
| `PATH_B_DEFAULT_plus_eq18_sec_source` | 19.00 | 19 | [-1.01, -0.04] | [-0.44, 0.44] | [1.00, 1.01] | plateau |
| `PATH_B_DEFAULT_plus_eq5_apply_to` | 19.00 | 19 | [-1.27, 1.13] | [-0.89, 0.83] | [1.00, 1.01] | plateau |
| `PATH_B_DEFAULT_plus_knot_daughter_di` | 19.00 | 19 | [-1.16, 1.24] | [-0.89, 0.84] | [1.00, 1.01] | plateau |
| `PATH_B_DEFAULT_plus_knot_threshold_gate` | 19.00 | 7 | [-0.90, 0.00] | [-0.38, 0.38] | [1.00, 1.28] | plateau |
| `PATH_B_DEFAULT_plus_laplacian` | 19.00 | 19 | [-1.10, -0.15] | [-0.47, 0.47] | [1.00, 1.01] | plateau |
| `PATH_B_DEFAULT_plus_lattice_orientation` | 19.00 | 19 | [-1.05, 0.94] | [-0.92, 0.93] | [1.00, 1.01] | plateau |
| `PATH_B_DEFAULT_plus_rep_form` | NaN | — | — | — | — | NaN |
| `PATH_B_DEFAULT_plus_rep_neighbour_set` | 19.00 | 19 | [-1.00, 0.34] | [-0.63, 0.63] | [1.00, 1.01] | plateau |
| `PATH_B_DEFAULT_plus_update_order` | 19.00 | 19 | [-0.10, 4.50] | [-0.00, 0.00] | [1.01, 1.01] | plateau |

Interpretation:

- Adding field `adh_form` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `adh_form` is therefore **not consequential**.
- Adding field `border_definition` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `border_definition` is therefore **not consequential**.
- Adding field `division_total_cap` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `division_total_cap` is therefore **not consequential**.
- Adding field `eq17_inh_source` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `eq17_inh_source` is therefore **not consequential**.
- Adding field `eq18_sec_source` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `eq18_sec_source` is therefore **not consequential**.
- Adding field `eq5_apply_to` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `eq5_apply_to` is therefore **not consequential**.
- Adding field `knot_daughter_di` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `knot_daughter_di` is therefore **not consequential**.
- Adding field `knot_threshold_gate` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `knot_threshold_gate` is therefore **not consequential**.
- Adding field `laplacian` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `laplacian` is therefore **not consequential**.
- Adding field `lattice_orientation` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `lattice_orientation` is therefore **not consequential**.
- Field `rep_form`: run failed or produced NaN plateau; effect cannot be assessed.
- Adding field `rep_neighbour_set` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `rep_neighbour_set` is therefore **not consequential**.
- Adding field `update_order` to `PATH_B_DEFAULT` shifts the plateau from 19.00 to 19.00 (delta +0.00; 0% of the 0.00-cell span between anchors) — field `update_order` is therefore **not consequential**.

