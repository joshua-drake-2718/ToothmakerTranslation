#!/usr/bin/env python3
"""Orchestrator for the Path B v2 autonomous build.

Reads .tmp/path-b-v2/progress.md to find the next pending phase,
prints a self-contained briefing for that phase, and exits. The
caller (a Claude /loop invocation) reads the briefing and runs the
phase as an Agent subtask, then re-invokes the orchestrator.

The orchestrator does NOT itself dispatch agents — Claude does that
via tool calls. This script is a single source of truth for phase
state, briefing content, and per-phase validation.

Usage from Claude (via /loop):
  python scripts/path-b-v2-orchestrate.py status   # print next phase + briefing
  python scripts/path-b-v2-orchestrate.py validate <phase>   # run validation
  python scripts/path-b-v2-orchestrate.py mark-done <phase>  # update progress.md

The /loop wrapper described in scripts/path-b-v2-loop.md runs:
  1. status → print briefing for next pending phase
  2. Claude reads briefing, dispatches Agent, waits
  3. Claude commits work
  4. validate <phase> → run pytest etc.
  5. mark-done <phase> → updates progress.md
  6. /loop fires again
"""
from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
PROGRESS = REPO / '.tmp' / 'path-b-v2' / 'progress.md'
CHARTER = REPO / 'docs' / 'plans' / '2026-05-05-path-b-v2-configurable-discretisation.md'


@dataclass(frozen=True)
class Phase:
    id: str
    title: str
    files: tuple[str, ...]      # files the phase will create or modify
    reads: tuple[str, ...]      # files the phase agent must read first
    success: str                # one-line success criterion
    briefing: str               # full prompt body for the agent
    validate_cmd: str | None    # shell command for `validate`; None to skip


PHASES: dict[str, Phase] = {
    'A3': Phase(
        id='A3',
        title='Wire Discretisation through reaction.py',
        files=('silicoshark/reaction.py',),
        reads=(
            'silicoshark/discretisation.py',
            'silicoshark/reaction.py',
            'silicoshark/state.py',
            'silicoshark/params.py',
            'silicoshark/mesh.py',
            'docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md',
            'docs/research/paper-review-2010-salazar-ciudad-jernvall.md',
        ),
        success=(
            'Each Discretisation reaction-related field branches at the '
            'call site with a docstring naming the option. Default '
            'behaviour (PATH_B_DEFAULT) preserves the existing v1 output.'
        ),
        briefing="""\
You are extending silicoshark/reaction.py to branch on the
Discretisation dataclass at every reaction-related decision point.

Charter: docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md
(see §Implementer-choice catalogue, especially the 'Reaction terms'
and 'Differentiation' subsections, and the summary table for the
classification of each field).

Fields that affect reaction.py:
  - eq14_denominator: 'inh_corrected' (default, FORTRAN-corrected)
        vs 'act_typo' (paper-as-written)
  - eq17_inh_source: 'act_concentration' (paper) vs
        'act_rate_times_di' (FORTRAN — uses temp rate variable)
  - eq18_sec_source: 'constant_k_sec' (paper) vs 'k_sec_times_di'
        (FORTRAN — d_i ramp)
  - diff_accumulator: 'sec' (default, FORTRAN) vs 'act' (paper Eq. 6)
  - knot_threshold_gate: 'none' (default, paper) vs
        'first_border_cell' (FORTRAN — only cells with index >=
        first_border_cell can become knots)

Implementation rules:
1. `step_reaction_diffusion` takes a Discretisation argument
   (default `PATH_B_DEFAULT`).
2. Each branch is a single short conditional at the call site, NOT
   wrapped in an abstract class or strategy object. Reader scanning
   the function should see every option exercised.
3. Each branch carries a one-line docstring or comment naming the
   field and citing the charter row.
4. Knot detection (state.knot |= state.act >= 1.0) gains an optional
   `i >= first_border_cell` gate when knot_threshold_gate ==
   'first_border_cell'. Compute `first_border_cell` from the initial
   layout: for the hex lattice with rad R, first_border_cell is
   6 * (R - 1). Pass it in via state, or compute from the topology
   (border-cell count).
5. Step the differentiation update INSIDE step_reaction_diffusion
   if you find that cleaner — but the existing simulator has
   step_differentiation as its own function. Add a Discretisation
   argument there too: 'sec' branch unchanged from v1; 'act' branch
   uses state.act.

The diff_accumulator branch lives in simulator.step_differentiation,
NOT reaction.py — leave a TODO comment in reaction.py pointing at
simulator.py for that one, and update simulator.step_differentiation
in this same phase.

Default (PATH_B_DEFAULT) behaviour must produce identical output to
the existing v1 reaction.step_reaction_diffusion + simulator.step_differentiation
for the seal example's first 100 iterations. Verify with a short
ad-hoc script.

Add a unit test file tests/test_silicoshark_reaction_branches.py
with these cases:
  - eq14_denominator='inh_corrected' on a 1-cell state with
    Act=0.5, Inh=0.5: rate = k_act*0.5/(1+k_inh*0.5) - k_deg*0.5
  - eq14_denominator='act_typo' on the same state:
    rate = k_act*0.5/(1+k_inh*0.5) - k_deg*0.5  (same — both
    reduce identically because Act==Inh)
    Use Act=0.3, Inh=0.5 for a case where they differ.
  - eq17_inh_source under both options on a knot+diff>Int cell.
  - eq18_sec_source under both options on a diff>=Set cell.
  - diff_accumulator under both options against a state with
    Sec=0.1, Act=0.3.

Run `pytest tests/test_silicoshark_reaction_branches.py` and
`pytest tests/test_simulator_smoke.py` (Path A regression). Both
must pass before you exit.

When done, append a brief summary to .tmp/path-b-v2/progress.md
under the appropriate ### A3 — done heading.
""",
        validate_cmd=(
            '.venv/bin/pytest tests/test_silicoshark_reaction_branches.py '
            'tests/test_simulator_smoke.py -q'
        ),
    ),
    'A4': Phase(
        id='A4',
        title='Wire Discretisation through forces.py',
        files=('silicoshark/forces.py',),
        reads=(
            'silicoshark/discretisation.py',
            'silicoshark/forces.py',
            'silicoshark/state.py',
            'silicoshark/mesh.py',
            'silicoshark/topology.py',
            'docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md',
            'docs/research/paper-review-2010-salazar-ciudad-jernvall.md',
        ),
        success=(
            'Each Discretisation force-related field branches at the '
            'call site. Default behaviour matches v1 forces.compute_forces '
            'and apply_border_multipliers.'
        ),
        briefing="""\
You are extending silicoshark/forces.py to branch on the
Discretisation dataclass at every mechanical-force decision point.

Charter: docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md
(see §Implementer-choice catalogue, 'Mechanical forces' and
'Border identification and biases' subsections).

Fields that affect forces.py:
  - eq5_z_gate: True (default, FORTRAN-style) vs False
        (no z gate; raw sum of unit vectors)
  - eq5_apply_to: 'all' (default, paper) vs 'interior_only' (FORTRAN)
  - rep_form: 'hookean_signed' (default, smooth Hookean both regimes)
        vs 'paper_gated' (active only when |d| < p0)
  - adh_form: 'unit_vector' (default, paper) vs
        'hookean_attraction' (FORTRAN's Hookean style)
  - rep_neighbour_set: 'mesh' (default, paper) vs
        'mesh_plus_all_close' (FORTRAN's pushingnovei)
  - border_bias_x_zero_quirk: True (default, FORTRAN-faithful for
        seal goldens) vs False (cells with x==0 also receive Pbi/Abi/Bgr)

Implementation rules:
1. `compute_forces` takes a Discretisation argument
   (default `PATH_B_DEFAULT`).
2. Each branch is a single conditional at the call site, with
   docstring/comment naming the field.
3. For `rep_neighbour_set='mesh_plus_all_close'`, you need an
   additional pass: for every cell pair within ~1.4 * rest, NOT
   already a mesh neighbour, add a `Rep`-scaled repulsion. This
   matches FORTRAN's `repel_non_neigh`. Use a kd-tree from
   scipy.spatial.cKDTree for O(N log N) lookup, then filter out
   mesh-neighbour pairs.
4. For `eq5_apply_to='interior_only'`, restrict the eq.5 force
   accumulation to non-border cells (use mesh.is_border).
5. For `border_bias_x_zero_quirk`:
   - True: preserve current apply_border_multipliers behaviour
     (cells with x exactly == 0 are skipped).
   - False: cells in the y-band, regardless of x sign, receive
     Pbi multiplier on x (sign-dependent: Pbi for x >= 0, Abi for
     x < 0, applied to abs(x)? or to signed x? Match the FORTRAN
     intent — see coreop2d.update_cell_position docstring); z-bgr
     applied to all band cells.

Default behaviour must produce identical output to the existing
v1 forces.compute_forces + apply_border_multipliers for the seal
example's first 100 iterations.

Add a unit test file tests/test_silicoshark_forces_branches.py:
  - rep_form: at |d|=0.8 inside vs paper_gated outside, hookean_signed
    fires uniformly; at |d|=1.2 outside, hookean_signed gives weak
    attraction, paper_gated gives 0.
  - adh_form: at |d|=1.5, unit_vector gives k_adh * d/1.5;
    hookean_attraction gives k_adh * d.
  - eq5_apply_to: with a 7-cell hex (rad=2), interior_only gives
    eq.5 force only on cell 0; all gives non-zero on all 7.
  - border_bias_x_zero_quirk: True skips x=0 cells; False applies
    multiplier to them.

Run pytest tests/test_silicoshark_forces_branches.py
tests/test_simulator_smoke.py — both must pass.

When done, append a summary to .tmp/path-b-v2/progress.md.
""",
        validate_cmd=(
            '.venv/bin/pytest tests/test_silicoshark_forces_branches.py '
            'tests/test_simulator_smoke.py -q'
        ),
    ),
    'A5': Phase(
        id='A5',
        title='Wire Discretisation through mesh.py + simulator.py + CLI',
        files=(
            'silicoshark/mesh.py',
            'silicoshark/simulator.py',
            'silicoshark/__main__.py',
            'silicoshark.py',
        ),
        reads=(
            'silicoshark/discretisation.py',
            'silicoshark/topology.py',
            'silicoshark/mesh.py',
            'silicoshark/simulator.py',
            'silicoshark/state.py',
            'silicoshark/io.py',
            'main.py',
            'docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md',
        ),
        success=(
            'silicoshark CLI runs end-to-end on examples/seal.txt with '
            '--preset selectable from ALL_PRESETS. mesh.py supports both '
            'topology=delaunay_each_step and static_with_local_update. '
            'simulator.step accepts a Discretisation argument and threads '
            'it to reaction, forces, division.'
        ),
        briefing="""\
You are wiring Discretisation through silicoshark/mesh.py,
silicoshark/simulator.py, and adding a CLI entry that mirrors
main.py ergonomics with a --preset flag.

Charter: docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md.

Fields newly active here:
  - laplacian: 'length_weighted' (default, v1) vs 'cotangent'
        (textbook discrete Laplace–Beltrami) vs 'fortran_margins'
        (reproduces 13.f90 apply_diffusion)
  - update_order: 'jacobi' (default) vs 'gauss_seidel_forces'
        (FORTRAN-style; reactions Jacobi, mechanical forces Gauss-Seidel)
        vs 'mixed_jacobi_gs' (alias for now)
  - topology: 'static_with_local_update' (default, USE silicoshark.topology)
        vs 'delaunay_each_step' (current v1 path)
  - division_per_step_cap, division_max_edge: only applied when
        topology=='delaunay_each_step'
  - knot_daughter_di: 'zero_reset' (default, paper) vs 'inherit_avg'
        (FORTRAN); applies in simulator.divide_cells
  - mesenchyme: 'absent' shortcut for parameter sets with no Sec.
        FOR THIS PHASE: implement only the 'absent' branch (default
        for v1 compatibility); leave a TODO for 'per_column_z_layers'
        (it's a larger task, deferred).

Implementation rules:

1. mesh.py:
   - Add a Mesh classmethod `from_topology(top, positions)` that
     consumes a Topology object's CSR materialisation, the same way
     `from_positions` consumes Delaunay output. Both should yield a
     Mesh instance with identical attributes.
   - Add a `cotangent_laplacian` Mesh method as an alternative to
     the existing length-weighted `laplacian` method. Selection
     happens via Discretisation.laplacian.
   - Mesh stays a thin computational layer; topology lifecycle lives
     in Topology.

2. simulator.py:
   - simulator.step(state, params, dt, disc=PATH_B_DEFAULT, top=None)
     where `top` is a persistent Topology when
     disc.topology=='static_with_local_update'. Caller is responsible
     for calling Topology.from_positions once before the loop and
     passing it forward.
   - When disc.topology=='delaunay_each_step', build a fresh Mesh
     each step (current v1 behaviour) and ignore `top`.
   - When disc.topology=='static_with_local_update', use
     top.to_csr() each step to build the Mesh, and call
     top.insert_daughter(...) inside divide_cells for each new cell
     so the topology stays consistent.
   - simulator.divide_cells gains a Discretisation argument; the
     knot_daughter_di branch picks the daughter's d_i.
   - simulator.run(state, params, n_iters, dt=DEFAULT_DT,
     disc=PATH_B_DEFAULT) constructs Topology from state.positions
     internally if disc.topology=='static_with_local_update', then
     iterates.
   - update_order branching: for 'jacobi' use the current pure-Jacobi
     code; for 'gauss_seidel_forces' apply each force component in
     turn, updating positions in place between components (matches
     FORTRAN's iteration() ordering — apply_diffusion, then forces
     accumulate, then apply_nuclear_traction directly to positions,
     then update_cell_position adds remaining force-deltas).

3. __main__.py + silicoshark.py:
   - silicoshark.py is a one-line wrapper:
     `from silicoshark.__main__ import main; main()`
   - __main__.py mirrors main.py's argparse:
     params_file output_dir output_name iterations save_blocks
     plus --preset NAME (default PATH_B_DEFAULT, choices from
     ALL_PRESETS).
     plus --override key=value (repeatable, applied via
     dataclasses.replace).
   - At save points, write OFF via silicoshark.io.write_off using
     the current Mesh's triangles (if delaunay_each_step) or
     the Topology's adjacency-derived triangles (if static).
     For static, derive triangles by walking 3-cycles in the CSR
     graph — add a helper Topology.triangles() that returns
     (T, 3) int array.

CLI smoke:
  .venv/bin/python -m silicoshark examples/seal.txt /tmp/out run 100 5 --preset PATH_B_DEFAULT
should write 5 OFF files without error. Cell-count plateau and
mesh stability are NOT validation criteria for this phase — that's
A6's job. This phase just needs the CLI to RUN end-to-end on
PATH_B_DEFAULT and LEGACY_FORTRAN.

Run pytest tests/test_simulator_smoke.py to confirm Path A still
green.

When done, append a summary to .tmp/path-b-v2/progress.md.
""",
        validate_cmd=(
            'rm -rf .tmp/path-b-v2/cli-run && '
            'TMPDIR=.tmp .venv/bin/python -m silicoshark examples/seal.txt '
            '.tmp/path-b-v2/cli-run run 100 5 --preset PATH_B_DEFAULT '
            '--override mesenchyme=absent && '
            'ls .tmp/path-b-v2/cli-run/*.off | wc -l && '
            '.venv/bin/pytest tests/test_silicoshark_a5_v1_replica.py '
            'tests/test_simulator_smoke.py -q'
        ),
    ),
    'A6': Phase(
        id='A6',
        title='LEGACY_FORTRAN reproduces seal goldens within tolerance',
        files=('tests/test_silicoshark_smoke.py',),
        reads=(
            'silicoshark/discretisation.py',
            'silicoshark/simulator.py',
            'silicoshark/forces.py',
            'silicoshark/reaction.py',
            'silicoshark/mesh.py',
            'silicoshark/topology.py',
            'tests/test_simulator_smoke.py',
            'tests/golden_fortran/500_run.off',
            'docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md',
            'docs/findings/2026-05-04-path-b-v1-mesh-degeneration.md',
        ),
        success=(
            'tests/test_silicoshark_smoke.py passes: LEGACY_FORTRAN preset '
            'on examples/seal.txt produces a 500-iter run with cell count '
            'within ±5% of 57, vertex envelope within 10% of '
            'tests/golden_fortran/500_run.off, no NaN/Inf, no RuntimeWarning.'
        ),
        briefing="""\
You are validating Path B v2 against the FORTRAN seal goldens.

The LEGACY_FORTRAN preset should reproduce the FORTRAN binary's
output on examples/seal.txt within the v1 tolerance (which is now
the bottom-row of the comparison-study matrix, not the only target):

  - Cell count plateau within ±5% of 57.
  - Vertex envelope within 10% of tests/golden_fortran/500_run.off.
  - No NaN/Inf in any vertex.
  - No RuntimeWarning during the run.

If the LEGACY_FORTRAN preset cannot meet this, the wiring from A3-A5
has a bug or there is a paper-vs-FORTRAN mismatch in one of the
Discretisation defaults that needs reconsidering.

Steps:

1. Write tests/test_silicoshark_smoke.py:
   - Subprocess-spawn `python -m silicoshark examples/seal.txt
     <tmp>/out run 100 5 --preset LEGACY_FORTRAN`.
   - Read the final OFF file (500_run_.off).
   - Read tests/golden_fortran/500_run.off.
   - Assert cell-count and envelope tolerances.
   - Assert no RuntimeWarning in stderr, no NaN/Inf coordinates.

2. Run the test. If it fails, diagnose:
   - Cell count too high → division logic too aggressive
     (check division_per_step_cap, max_edge filters, knot_daughter_di).
   - Cell count too low → division logic not triggering
     (check static_with_local_update is wired correctly).
   - Envelope mismatch in xy → check border-bias multipliers,
     eq5 z-gate, rep/adh forms.
   - Envelope mismatch in z → check Bgr scaling, Dgr value flow,
     buoyancy.
   - NaN/Inf → check division-by-zero in cervical-loop, mesenchyme
     thickening, reaction terms.

3. Iterate fixes — modify silicoshark/* as needed. Each fix should
   be SMALL and TARGETED at a specific divergence; do not rewrite
   the simulator. Commit each fix with a clear message.

4. When the test passes, run pytest tests/test_simulator_smoke.py
   AND tests/test_silicoshark_smoke.py — both must pass.

5. Append a summary to .tmp/path-b-v2/progress.md, including the
   final cell count, vertex envelope, and any discoveries about
   which Discretisation field was load-bearing.

Budget: this phase may take significant iteration. The previous v1
attempt failed at this exact point because the topology was
delaunay_each_step. With static_with_local_update from A2, plus
LEGACY_FORTRAN's other paper-vs-FORTRAN choices, the run should be
much more stable.

If after substantial iteration the seal goldens still cannot be
matched, document the remaining divergence in
docs/findings/<date>-path-b-v2-seal-divergence.md and reduce the
test tolerance to whatever is achievable (with justification),
rather than blocking forever.
""",
        validate_cmd=(
            '.venv/bin/pytest tests/test_silicoshark_smoke.py '
            'tests/test_simulator_smoke.py -q'
        ),
    ),
    'A7': Phase(
        id='A7',
        title='Comparison-study scaffold: presets × params × metrics',
        files=(
            'silicoshark/metrics.py',
            'scripts/run-discretisation-study.py',
            'experiments/discretisation-study/README.md',
        ),
        reads=(
            'silicoshark/discretisation.py',
            'silicoshark/simulator.py',
            'silicoshark/state.py',
            'docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md',
            'docs/research/paper-review-2014-harjunmaa.md',
        ),
        success=(
            'scripts/run-discretisation-study.py runs all named presets '
            'on examples/seal.txt and writes results.json with metrics. '
            'ProgressReporter integrated for exptop visibility. '
            'metrics.py has unit tests for each metric on a known state.'
        ),
        briefing="""\
You are building the comparison-study scaffold that will eventually
power the methodological paper.

Charter: docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md
§Validation strategy.

Deliverables:

1. silicoshark/metrics.py:
   - cusp_count(state, params) -> int. Number of knot cells
     (state.knot.sum()) at the time of measurement.
   - cusp_width(state) -> float. Mean pairwise distance between
     adjacent knot cells; 0 if fewer than 2 knots. Use mesh
     adjacency to find adjacent knots.
   - cell_count_plateau(history: list[int]) -> int. The mean of
     the last 10% of the cell-count history (assumes the simulation
     has reached plateau).
   - vertex_envelope(positions) -> dict[str, float]. min/max x/y/z.
   - regime(history: list[int]) -> str. 'monotone-grow' if last
     20% strictly increasing; 'plateau' if last 20% std < 5% of
     mean; 'oscillate' if multiple sign changes in finite difference;
     'collapse' if final < 0.5 * peak; 'NaN' if any NaN encountered.

2. scripts/run-discretisation-study.py:
   - CLI: `python scripts/run-discretisation-study.py --params PATH
     [--params PATH ...] [--presets PRESET ...] [--iters N]
     [--saves S] [--out DIR]`
   - Defaults: params=examples/seal.txt, presets=ALL_PRESETS keys,
     iters=500, saves=5, out=experiments/discretisation-study/results.
   - For each (preset, params_file) combination:
     - Run the simulation, capturing cell-count history at each
       save block plus the final state.
     - Compute all metrics on the final state and history.
     - Save metrics, final OFF, params snapshot, preset snapshot
       to experiments/discretisation-study/results/<preset>/<params_basename>/.
   - Aggregate into experiments/discretisation-study/results/results.json
     and a markdown table per parameter file.
   - ProgressReporter integration:
     - WORK level (batch_runner_discretisation_study): one phase per
       (preset, params_file).
     - TASK level: each individual run writes its own progress.json.

3. experiments/discretisation-study/README.md: describe the study
   structure, list current parameter files (start with seal.txt;
   the 2014 paper's Ext Data Fig. 2 wild-type tribosphenic and
   Act sweep are stretch goals deferred to the actual paper-writing
   phase).

4. tests/test_metrics.py: unit tests for each metric.

Constraints:
- The script must not crash on a divergent preset (e.g., if
  PATH_B_DEFAULT explodes on seal.txt, capture the failure as
  regime='collapse' or 'NaN' and continue with the next preset).
- Use the existing silicoshark.utils.progress.ProgressReporter
  if it exists; otherwise copy from ~/repo/workflow/lib/progress.py
  per the global CLAUDE.md.

Run a small smoke pass:
  python scripts/run-discretisation-study.py --presets PATH_B_DEFAULT
    --iters 100 --saves 5
to confirm it produces a valid results.json without taking forever.

When done, append a summary to .tmp/path-b-v2/progress.md.
""",
        validate_cmd=(
            '.venv/bin/pytest tests/test_metrics.py -q && '
            'TMPDIR=.tmp .venv/bin/python scripts/run-discretisation-study.py '
            '--presets PATH_B_DEFAULT --iters 50 --saves 5 '
            '--out .tmp/path-b-v2/study-smoke && '
            'test -f .tmp/path-b-v2/study-smoke/results.json'
        ),
    ),
    'A8': Phase(
        id='A8',
        title='Audit-trail: per-field paper + FORTRAN citations',
        files=('docs/research/discretisation-audit.md',),
        reads=(
            'silicoshark/discretisation.py',
            'docs/plans/2026-05-05-path-b-v2-configurable-discretisation.md',
            'docs/research/paper-review-2010-salazar-ciudad-jernvall.md',
            'docs/research/paper-review-2014-harjunmaa.md',
            'docs/research/13f90-vs-humppa-divergences.md',
            'docs/research/cpp-port-review.md',
            'docs/research/tgrohens-review.md',
            'docs/research/paper-review-2023-zimm.md',
        ),
        success=(
            'docs/research/discretisation-audit.md exists with one '
            'section per Discretisation field. Each option of each '
            'field is traced to (a) a paper passage or its absence, '
            '(b) a FORTRAN line, (c) any third-source resolution.'
        ),
        briefing="""\
You are writing the canonical audit trail document for Path B v2.

This is the methodological deliverable: a per-field, per-option
record of where the choice comes from, what evidence supports each
side, and how alternative code-bases (humppa, C++ port, tgrohens,
silicoshark FORTRAN) resolve it. The follow-up paper will draw on
this directly.

File: docs/research/discretisation-audit.md

Structure:

```markdown
---
title: 'Path B v2 — Discretisation field audit trail'
author: Lyndon Drake (with Claude Code)
date: 2026-05-05
---

## Purpose

[2-3 paragraphs explaining what this document is, who it's for,
and how it differs from the charter's choice catalogue]

## How to read this document

[brief: each field is a heading, each option a sub-heading, with
  citations under each]

## Fields

### `laplacian`

**Question:** [...]

**Kind:** PaperAmbiguity / PaperVsCodeTension / FortranAccident

**Paper evidence:**
- 2010 main p. 585: 'finite-volume on the triangular mesh, with
  flux proportional to contact area'. No formula given.
- 2010 SI: [if relevant]
- 2014: [n/a; paper uses 2010 model as black box]

**FORTRAN evidence:**
- 13.f90/coreop2d.py:506-617 (apply_diffusion). Per-edge `pes`
  margins, per-cell-z `area_p` weighting, with hand-tuned 0.44
  factor at the substrate boundary.
- humppa_translate.f90: identical to 13.f90 here (no rename
  divergence affects this function).

**Alternative code-bases:**
- C++ port: [check cpp-port-review.md for diffusion semantics]
- tgrohens: [check]
- silicoshark FORTRAN: [check]

**Options:**
- `length_weighted`: 1 / |edge_ij|. Justification: [...]
- `cotangent`: textbook discrete Laplace–Beltrami. Justification:
  [...]
- `fortran_margins`: reproduces 13.f90's `apply_diffusion`.
  Justification: cross-validation oracle.

**Default:** [whichever Discretisation has, with reason]

**Comparison-study question:** Does the choice of Laplacian change
the cusp pattern, or is it a numerical-noise-level difference?

---

### `update_order`

[similar structure]

[... and so on for all 20 fields ...]

## Cross-field interactions

[brief: which combinations are interesting? e.g.
update_order='jacobi' × topology='static_with_local_update' is the
'cleanest numpy' combo; LEGACY_FORTRAN bundles every FORTRAN-side
choice; PAPER_LITERAL_2010 isolates the eq.14 typo.]
```

Constraints:
- Cite specific paper page numbers and FORTRAN line numbers.
- For each option, state WHO uses it: paper text, 13.f90, humppa,
  C++ port, tgrohens, silicoshark FORTRAN, Path B v1, etc.
- For each option, name what the comparison study expects to
  reveal. (e.g., 'expect this to change cusp count' vs 'expect
  this to be numerical noise')
- Length budget: aim for 2000-4000 words total. Don't repeat the
  charter's content; this document is more granular.

Verify with: file exists, has a section for every field in
silicoshark.discretisation.Discretisation, every option named in
the dataclass appears in this document at least once.

Run a quick consistency check:
  for field in dataclass: ensure docstring's options match
  the audit's options.

When done, append a summary to .tmp/path-b-v2/progress.md.
""",
        validate_cmd=(
            '.venv/bin/python -c "'
            'from silicoshark.discretisation import Discretisation; '
            'from dataclasses import fields; '
            'audit = open(\'docs/research/discretisation-audit.md\').read(); '
            'missing = [f.name for f in fields(Discretisation) '
            'if f.name not in audit]; '
            'assert not missing, f\'audit missing fields: {missing}\'; '
            'print(\'audit covers all fields\')"'
        ),
    ),
}


PHASE_ORDER = ['A3', 'A4', 'A5', 'A6', 'A7', 'A8']


def parse_progress() -> dict[str, str]:
    """Return {phase_id: status} from the progress markdown table."""
    if not PROGRESS.exists():
        sys.exit(f'progress file missing: {PROGRESS}')
    text = PROGRESS.read_text()
    out: dict[str, str] = {}
    for m in re.finditer(r'^\|\s*(A\d)\s*\|.*?\|\s*([^|]+?)\s*\|\s*$', text, re.MULTILINE):
        phase, status = m.group(1), m.group(2).strip()
        out[phase] = status
    return out


def next_pending(progress: dict[str, str]) -> str | None:
    for p in PHASE_ORDER:
        status = progress.get(p, '')
        if not status.startswith('done'):
            return p
    return None


def cmd_status() -> None:
    progress = parse_progress()
    nxt = next_pending(progress)
    if nxt is None:
        print('ALL PHASES COMPLETE.')
        sys.exit(0)
    phase = PHASES[nxt]
    files_block = '\n'.join(f'  - {f}' for f in phase.files)
    reads_block = '\n'.join(f'  - {f}' for f in phase.reads)
    print(f'== Path B v2 — Phase {phase.id}: {phase.title} ==')
    print()
    print(f'Files this phase will create or modify:\n{files_block}')
    print()
    print(f'Files the phase agent must read first:\n{reads_block}')
    print()
    print(f'Success criterion:\n  {phase.success}')
    print()
    print('Briefing for the Agent subtask (paste verbatim):')
    print('---8<---')
    print(phase.briefing)
    print('---8<---')
    if phase.validate_cmd:
        print(f'\nValidation command (run after the agent finishes):')
        print(f'  {phase.validate_cmd}')


def cmd_validate(phase_id: str) -> int:
    if phase_id not in PHASES:
        sys.exit(f'unknown phase: {phase_id}')
    phase = PHASES[phase_id]
    if not phase.validate_cmd:
        print(f'Phase {phase_id} has no automated validation; manual check required.')
        return 0
    print(f'Validating phase {phase_id}: {phase.validate_cmd}')
    proc = subprocess.run(
        phase.validate_cmd, shell=True, cwd=REPO,
        capture_output=False, check=False,
    )
    return proc.returncode


def cmd_mark_done(phase_id: str) -> None:
    if phase_id not in PHASES:
        sys.exit(f'unknown phase: {phase_id}')
    text = PROGRESS.read_text()
    head = subprocess.run(
        ['git', 'rev-parse', '--short', 'HEAD'],
        cwd=REPO, capture_output=True, text=True, check=True,
    ).stdout.strip()
    # Accept either 'pending' or a bare 'done' (without commit hash) so
    # we can stamp a hash retroactively when an agent updated the row
    # directly without going through `mark-done`.
    new_text = re.sub(
        rf'^(\|\s*{phase_id}\s*\|.*?\|)\s*(pending|done)\s*(\|)\s*$',
        rf'\1 done (commit `{head}`) \3',
        text, count=1, flags=re.MULTILINE,
    )
    if new_text == text:
        sys.exit(f'phase {phase_id} status not stampable in progress.md')
    PROGRESS.write_text(new_text)
    print(f'marked {phase_id} done at commit {head}')


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd', required=True)
    sub.add_parser('status')
    v = sub.add_parser('validate'); v.add_argument('phase')
    m = sub.add_parser('mark-done'); m.add_argument('phase')
    args = parser.parse_args()
    if args.cmd == 'status':
        cmd_status()
        return 0
    if args.cmd == 'validate':
        return cmd_validate(args.phase)
    if args.cmd == 'mark-done':
        cmd_mark_done(args.phase)
        return 0
    return 0


if __name__ == '__main__':
    sys.exit(main())
