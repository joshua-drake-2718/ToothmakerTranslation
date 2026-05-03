"""Structural smoke test for the simulator end-to-end run.

Asserts: the simulator completes a 5-block run on examples/seal.txt without
crashing, writes 5 OFF files, and each file has a valid COFF header with
sensible vertex/face counts.

Does NOT byte-compare to the FORTRAN goldens — the Python iteration code has
known bugs (NaN z-coords, frozen cell count) that are out of scope for this
migration. See Task 9 for the FORTRAN-comparison follow-up.
"""
import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent


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

    off_files = sorted(out_dir.glob("*.off"))
    assert len(off_files) == 5, (
        f"expected 5 OFF files, got {len(off_files)}: "
        f"{[p.name for p in off_files]}"
    )

    for off in off_files:
        with off.open() as f:
            header = f.readline().strip()
            counts = f.readline().strip()
        assert header == "COFF", f"{off.name}: expected 'COFF' header, got {header!r}"
        m = re.match(r"^\s*(\d+)\s+(\d+)\s+0\s*$", counts)
        assert m, f"{off.name}: expected 'V F 0' counts line, got {counts!r}"
        n_vertices, n_faces = int(m.group(1)), int(m.group(2))
        assert n_vertices > 0, f"{off.name}: zero vertices"
        assert n_faces > 0, f"{off.name}: zero faces"
