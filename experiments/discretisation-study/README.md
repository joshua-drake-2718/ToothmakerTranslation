---
title: 'Path B v2 discretisation comparison study'
date: 2026-05-05
---

## Purpose

This study runs each named `Discretisation` preset across one or
more parameter files and captures cusp count, cusp width, cell-count
plateau, vertex envelope, and qualitative regime. It is the
empirical evidence base for the methodological paper described in
the v2 charter (`docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md`,
§Validation strategy): the comparison matrix is what shows which
implementer-choice points materially change biological output and
which are numerically inert.

The replication anchor is `LEGACY_FORTRAN`; every other preset's
results are read against that one. The seal-baseline cell of the
matrix is the v1 single-tolerance test repurposed: a presets-and-
options decomposition that fails to reproduce `LEGACY_FORTRAN`'s
seal output indicates a wiring bug in the catalogue, not a finding.

## Directory layout

```
experiments/discretisation-study/
├── README.md                 (this file)
├── seal-baseline/            (committed seal.txt results from a full pass)
│   ├── results.json
│   ├── results-table.md
│   └── <preset>/seal/        (per-run artefacts)
│       ├── *.off              ← OFF files (one per save block)
│       ├── metrics.json
│       ├── params-snapshot.txt
│       ├── preset-snapshot.json
│       ├── progress.json      (silicoshark CLI's per-run progress)
│       ├── silicoshark.stdout
│       └── silicoshark.stderr
└── results/                  (uncommitted; default --out target)
```

`results.json` is a nested dict
`{preset: {params_basename: {metrics + regime + cell_count_history}}}`
covering every (preset × params) combination.

`results-table.md` renders one markdown table per parameter file —
rows are presets, columns are the headline metrics.

`progress.json` (per-run, in each `<preset>/<params>/` directory)
is the silicoshark CLI's `ProgressReporter` output. The batch-level
`progress.json` (in `--out`) carries the experiment name
`batch_runner_discretisation_study` so `exptop` displays it as the
WORK line; per-run `progress.json` files are picked up as the TASK
line.

## How to run

The full pass on `seal.txt` with all five presets:

```bash
TMPDIR=.tmp .venv/bin/python scripts/run-discretisation-study.py \
    --iters 500 --saves 5 \
    --out experiments/discretisation-study/results \
    --override mesenchyme=absent --override laplacian=length_weighted
```

Smoke-run on a single preset, fewer iterations:

```bash
TMPDIR=.tmp .venv/bin/python scripts/run-discretisation-study.py \
    --presets PATH_B_DEFAULT --iters 50 --saves 5 \
    --out .tmp/path-b-v2/study-smoke \
    --override mesenchyme=absent --override laplacian=length_weighted
```

The two `--override` flags are the deferred-field workarounds for
the not-yet-implemented `mesenchyme='per_column_z_layers'` and
`laplacian='fortran_margins'`. The runner detects whether each
preset needs them and adds them automatically as well; passing
them on the command line is harmless (the user values win), and is
useful as documentation that the run is operating with a known
caveat.

The currently-shipping presets all need both auto-overrides:

| Preset | mesenchyme | laplacian | Effect |
|---|---|---|---|
| `PATH_B_DEFAULT` | `per_column_z_layers` | `length_weighted` | mes auto-overridden |
| `PAPER_2010` | `per_column_z_layers` | `length_weighted` | mes auto-overridden |
| `PAPER_LITERAL_2010` | `per_column_z_layers` | `length_weighted` | mes auto-overridden |
| `LEGACY_FORTRAN` | `per_column_z_layers` | `fortran_margins` | both auto-overridden |
| `HUMPPA_LITERAL` | `per_column_z_layers` | `fortran_margins` | both auto-overridden |

Closing these gaps (per-column z-layer mesenchyme;
`fortran_margins` Laplacian) is on the v2 backlog
(see `docs/findings/2026-05-05-path-b-v2-legacy-fortran-divergence.md`).

## Currently supported parameter files

`examples/seal.txt` is the single supported parameter file at A7.
The 2014 paper's qualitative validation targets are deferred to
their own phase:

- **Wild-type tribosphenic mouse** (Ext Data Fig. 2): the 25
  parameters need transcribing into a `wt-trib.txt`-style file.
  Source: `docs/research/paper-review-2014-harjunmaa.md` §Ext Data
  Fig. 2.
- **Act sweep** (Ext Data Fig. 4): a sweep-level wrapper around the
  wild-type parameter file with `Act ∈ [0.1, 1.6]` in 0.1
  increments. Source: same review doc, §Fig. 4.
- **Inh / Di perturbation lattice** (Ext Data Fig. 5): five (Inh,
  Di) points perturbing the wild-type baseline.

When those parameter files are added, run them via the same
`--params` flag (it accepts multiple paths):

```bash
TMPDIR=.tmp .venv/bin/python scripts/run-discretisation-study.py \
    --params examples/seal.txt examples/wt-trib.txt \
    --iters 500 --saves 5 \
    --out experiments/discretisation-study/results
```

## Metrics catalogue

Each per-run `metrics.json` carries:

- `cell_count_history`: `[n_cells_at_block_1, ..., n_cells_at_block_S]`.
- `cell_count_plateau`: mean of the last 10% of the history (final
  value if the history is shorter than 10 entries).
- `cell_count_final`: integer final cell count.
- `cusp_count`: number of knot cells (`[Act] >= 1`) at the final
  save block. Inferred from the OFF colour codes (knot cells are
  cyan RGBA `(0, 1, 1, 0)` per `silicoshark.io.cell_colour`).
- `cusp_width`: mean pairwise distance between mesh-adjacent knot
  cells. `0.0` when fewer than two knots exist.
- `vertex_envelope`: `{x_min, x_max, y_min, y_max, z_min, z_max}`.
- `regime`: one of `'NaN'`, `'collapse'`, `'oscillate'`,
  `'plateau'`, `'monotone-grow'` (see `silicoshark/metrics.py`).

`returncode` and `note` (if any) flag failed runs. The runner does
not crash on failures: a NaN or exception in one preset produces a
`'NaN'` regime entry; the remaining presets continue.
