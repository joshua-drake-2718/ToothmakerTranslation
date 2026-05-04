"""Structural smoke test for the simulator end-to-end run.

Runs the simulator on examples/seal.txt and asserts that the OFF outputs
match the FORTRAN reference structurally: 57 cells / 236 faces (the
equilibrium plateau for this parameter set), no NaN/Inf coordinates, and
no stderr noise from the diffusion code.

Catches regressions of the three Path-A fixes:
  - apply_diffusion FMA-limit guard (NaN z-coords if missing)
  - update_cell_position Bgr placement (cell count != 57 if missing)
  - add_cell iiii preservation (cell count >> 57 if missing)
"""
import math
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

EXPECTED_VERTICES = 57
EXPECTED_FACES = 236


def test_seal_run_smoke(tmp_path):
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    result = subprocess.run(
        [
            sys.executable, "main.py",
            str(REPO_ROOT / "examples" / "seal.txt"),
            str(out_dir),
            "run", "100", "5",
        ],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        f"simulator exited {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    assert "RuntimeWarning" not in result.stderr, (
        "diffusion code emitted RuntimeWarnings (likely divide-by-zero):\n"
        f"{result.stderr}"
    )

    off_files = sorted(out_dir.glob("*.off"))
    assert len(off_files) == 5, (
        f"expected 5 OFF files, got {len(off_files)}: "
        f"{[p.name for p in off_files]}"
    )

    for off in off_files:
        with off.open() as f:
            header = f.readline().strip()
            counts = f.readline().strip()
            vertex_lines = [f.readline() for _ in range(EXPECTED_VERTICES)]
        assert header == "COFF", f"{off.name}: expected 'COFF' header, got {header!r}"
        m = re.match(r"^\s*(\d+)\s+(\d+)\s+0\s*$", counts)
        assert m, f"{off.name}: expected 'V F 0' counts line, got {counts!r}"
        n_vertices, n_faces = int(m.group(1)), int(m.group(2))
        assert n_vertices == EXPECTED_VERTICES, (
            f"{off.name}: expected {EXPECTED_VERTICES} vertices, got {n_vertices}"
        )
        assert n_faces == EXPECTED_FACES, (
            f"{off.name}: expected {EXPECTED_FACES} faces, got {n_faces}"
        )

        # Each vertex line: x y z r g b a — three coordinates then four colour
        # channels. None should be NaN/Inf.
        for line_no, line in enumerate(vertex_lines):
            x, y, z = (float(v) for v in line.split()[:3])
            assert all(math.isfinite(v) for v in (x, y, z)), (
                f"{off.name} vertex {line_no}: non-finite coordinate "
                f"x={x} y={y} z={z}"
            )
