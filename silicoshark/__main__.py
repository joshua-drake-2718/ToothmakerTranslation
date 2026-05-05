"""silicoshark CLI entry point.

Mirrors `main.py`'s argparse signature so existing run scripts work with
either Path A (`python main.py …`) or Path B (`python -m silicoshark …`):

    python -m silicoshark PARAMS_FILE OUT_FOLDER OUT_NAME ITERATIONS SAVE_BLOCKS
                          [--preset NAME] [--override key=value ...]

Each save block runs ITERATIONS iterations, then writes one OFF file
named `<OUT_FOLDER>/<iter*block>_<OUT_NAME>_.off` mirroring the FORTRAN
binary's emission cadence.

`--preset NAME` selects one of the named `Discretisation` presets
(`PATH_B_DEFAULT`, `PAPER_2010`, `PAPER_LITERAL_2010`, `LEGACY_FORTRAN`,
`HUMPPA_LITERAL`). `--override key=value` patches individual fields on
the chosen preset via `dataclasses.replace`; values are coerced
int -> float -> bool -> str, and unknown keys are rejected.

Example:

    python -m silicoshark examples/seal.txt out/ run 100 5 \\
        --preset PATH_B_DEFAULT
"""
from __future__ import annotations

import argparse
import os
import sys
from dataclasses import fields, replace
from pathlib import Path

import numpy as np

from .discretisation import ALL_PRESETS, Discretisation
from .io import write_off
from .params import Params
from .simulator import DEFAULT_DT, build_topology, run
from .state import build_initial_state
from .utils import ProgressReporter


def _coerce(value: str):
    """Coerce a string to None -> int -> float -> bool -> str.

    The literal `null` (JSON convention) coerces to Python `None`. The
    string `none` is *not* coerced to `None` — `knot_threshold_gate`
    accepts the literal string value `'none'` as one of its valid
    options, and that case must round-trip through the CLI without
    collapsing to Python's None singleton.
    """
    low = value.lower()
    # Bool first — `'1'`/`'0'` should be int, not bool, so we check
    # bool only against the explicit literals.
    if low in ('true', 'false'):
        return low == 'true'
    if low == 'null':
        return None
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        pass
    return value


def _apply_overrides(disc: Discretisation, overrides: list[str]) -> Discretisation:
    """Apply a list of `key=value` strings to `disc` via dataclasses.replace.

    Unknown keys raise SystemExit with a clear error.
    """
    if not overrides:
        return disc
    valid = {f.name for f in fields(disc)}
    patches: dict[str, object] = {}
    for ov in overrides:
        if '=' not in ov:
            raise SystemExit(f'malformed override (expected key=value): {ov!r}')
        key, _, raw = ov.partition('=')
        key = key.strip()
        if key not in valid:
            raise SystemExit(
                f'unknown Discretisation field: {key!r}. '
                f'Valid: {sorted(valid)}'
            )
        patches[key] = _coerce(raw.strip())
    return replace(disc, **patches)


def _check_implemented(disc: Discretisation) -> None:
    """Reject preset/override combinations whose code path is not yet
    implemented. Currently:
      - laplacian == 'fortran_margins' (deferred to its own sub-project).

    `mesenchyme='per_column_z_layers'` was unimplemented in A5 and is
    now wired through `silicoshark/reaction.step_mesenchyme_diffusion`
    plus the cervical-loop / buoyancy reads in `silicoshark/forces.py`
    (Path B v2 B3, 2026-05-05).
    """
    if disc.laplacian == 'fortran_margins':
        raise SystemExit(
            "discretisation field laplacian='fortran_margins' is not "
            'implemented in this phase (the FORTRAN per-edge margin '
            'scheme is documented as a separate sub-project). Use '
            '--override laplacian=length_weighted (or =cotangent).'
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog='silicoshark',
        description=(
            'silicoshark — Path B re-implementation of the 2010 tooth-'
            'morphogenesis model. CLI mirrors main.py with --preset and '
            '--override flags for Discretisation configuration.'
        ),
    )
    parser.add_argument('params_file', help='Parameter file path')
    parser.add_argument('output_folder', help='Folder for output files')
    parser.add_argument('output_name', help='Base name for output files')
    parser.add_argument(
        'iterations', type=int,
        help='Number of forward-Euler iterations per save block',
    )
    parser.add_argument(
        'save_blocks', type=int,
        help='Number of save blocks (one OFF file per block)',
    )
    parser.add_argument(
        '--preset', default='PATH_B_DEFAULT',
        choices=sorted(ALL_PRESETS.keys()),
        help='Discretisation preset (default: PATH_B_DEFAULT)',
    )
    parser.add_argument(
        '--override', action='append', default=[],
        metavar='KEY=VALUE',
        help='Override a Discretisation field. May be repeated. '
        'Values are coerced int->float->bool->str.',
    )
    args = parser.parse_args(argv)

    disc = ALL_PRESETS[args.preset]
    disc = _apply_overrides(disc, args.override)
    _check_implemented(disc)

    params = Params.from_file(args.params_file)
    state = build_initial_state(params, lattice_orientation=disc.lattice_orientation)
    top = build_topology(state) if disc.topology == 'static_with_local_update' else None

    out_folder = Path(args.output_folder)
    out_folder.mkdir(parents=True, exist_ok=True)

    sstep = abs(args.save_blocks)

    # ProgressReporter integration (per CLAUDE.md §Experiment monitoring).
    # One phase per save block so the TASK line in exptop shows
    # something like 'block 3/5'. Experiment name embeds the preset and
    # output name so the parent batch runner's WORK line plus this TASK
    # line together identify the run.
    progress_path = out_folder / 'progress.json'
    experiment_name = f'silicoshark_{args.preset}_{args.output_name}'
    with ProgressReporter(
        progress_path, experiment=experiment_name, total_phases=sstep
    ) as progress:
        for block in range(1, sstep + 1):
            progress.begin_phase(
                name=f'block_{block}',
                description=f'Save block {block}/{sstep}',
                total_steps=args.iterations,
            )
            mesh = run(state, params, args.iterations, DEFAULT_DT, disc, top)

            iter_label = str(block * args.iterations)
            nff = out_folder / (iter_label + '_' + args.output_name)
            # Mirror main.py's output-name normalisation (spaces -> underscores).
            nff = out_folder / nff.name.replace(' ', '_')
            nfoff = nff.with_name(nff.name + '_.off')

            # Triangles come from whichever Mesh the run last built.
            with open(nfoff, 'w') as f:
                write_off(state, mesh.triangles, f)

            # Mark step at end of block so step_pct reflects completion.
            progress.step(args.iterations, f'wrote {nfoff.name}')
            print(f'Block {block}/{sstep} complete: {nfoff}')

    return 0


if __name__ == '__main__':
    sys.exit(main())
