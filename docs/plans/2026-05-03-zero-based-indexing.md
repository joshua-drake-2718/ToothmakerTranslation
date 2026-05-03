# Zero-Based Indexing Migration Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the 1-based `indexing.py` wrappers (a transient FORTRAN-translation aid) and convert the entire Python codebase to native 0-based numpy indexing, while preserving simulator output bit-for-bit.

**Architecture:** Capture a reference output from the current 1-based code first ("golden" OFF files). Then convert each source file to use plain numpy arrays and the builtin `range`. The biggest semantic risk is the sentinel value `0` in `neigh[]` (meaning "no neighbour") — under 0-based indexing, `0` becomes a valid cell index, so the sentinel must become `-1`. After migration, the same simulator run on `examples/seal.txt` must produce byte-identical OFF output.

**Tech Stack:** Python 3.11+, NumPy, pytest.

---

## Background

`indexing.py` provides four wrapper classes (`int_array`, `float_array`, `bool_array`, `str_array`) that subtract 1 from every index, plus an overridden `range` whose semantics are FORTRAN-style (`range(N)` iterates `1..N` inclusive; `range(a, b)` iterates `a..b` inclusive). All four source files (`coreop2d.py`, `esclec.py`, `main.py`, `vector.py`) `from indexing import *` and use 1-based access throughout.

The wrappers were a translation convenience so the Python could mirror the FORTRAN line-for-line. Now that the translation is largely complete, the wrappers add overhead, hide bugs, and make the code unidiomatic. Removing them lets us use numpy slicing, broadcasting, and any third-party Python tool that expects ordinary arrays.

## Translation cheat sheet

These patterns recur. Apply them mechanically except where noted.

| 1-based pattern (current) | 0-based pattern (target) | Notes |
|---|---|---|
| `from indexing import *` | (delete; import only what's needed from numpy) | Final state. |
| `int_array((N, M))` | `np.zeros((N, M), dtype=np.int32)` | Already pre-zeroed; remove subsequent `arr[:] = 0`. |
| `float_array((N,))` | `np.zeros(N)` | dtype defaults to float64. |
| `bool_array((N,))` | `np.zeros(N, dtype=bool)` |  |
| `str_array(N)` | `np.empty(N, dtype=object)` | Strings of arbitrary length. |
| `int_array(buffer=[1,2,3])` | `np.array([1,2,3], dtype=np.int32)` |  |
| `range(N)` (overridden, ⇒ 1..N) | `range(N)` (builtin, ⇒ 0..N-1) | **Loop body indices shift by -1.** |
| `range(a, b)` (overridden, ⇒ a..b) | `range(a-1, b)` (builtin, ⇒ a-1..b-1) | Same N-1 iterations, all shifted by -1. |
| `range(a, b, -1)` (overridden, ⇒ a..b decreasing) | `range(a-1, b-2, -1)` | Both endpoints shift by -1; stop is exclusive in builtin. |
| `arr[1, 1]` (literal indices) | `arr[0, 0]` | Subtract 1 from each numeric literal. |
| `arr[1:5]` | `arr[0:4]` | Both endpoints shift by -1; stop is exclusive in builtin. |
| `arr[i, j]` where `i, j` are loop vars | `arr[i, j]` (no change if loop bounds were also shifted) | Only change literals; shifted-loop variables are correct as-is. |
| `arr.copy()` (wrapper) | `arr.copy()` (numpy) | Identical method name. |
| `arr.inner` | `arr` | The wrapper's underlying numpy array; no longer needed. |
| `c = float_array(4); c[1]..c[4]` | `c = np.zeros(4); c[0]..c[3]` | Allocation size unchanged; literal access shifts. |
| `cls.neigh[i, j] == 0` (sentinel) | `cls.neigh[i, j] == -1` | **Critical: see Sentinel section.** |
| `cls.neigh[i, j] = 0` (sentinel) | `cls.neigh[i, j] = -1` | **Critical.** |
| `cls.neigh[:] = 0` (init) | `cls.neigh[:] = -1` | **Critical.** |
| `cls.neigh[i, j] > cls.num_active_cells` | `cls.neigh[i, j] >= cls.num_active_cells` | Boundary check shifts: in 1-based, indices > N are out-of-active; in 0-based, indices >= N are out-of-active. |
| `cls.neigh[i, j] == cls.num_all_cells` | `cls.neigh[i, j] == cls.num_all_cells` (no change) OR `>= cls.num_all_cells` | The "boundary cell marker" value `num_all_cells` is itself a sentinel — see notes in Task 7. |
| `f.write(f"{i - 1}")` (OFF output adjustment) | `f.write(f"{i}")` | The `-1` was converting 1-based Python indices to 0-based file-format indices; now Python is already 0-based. |
| `for i in range(2, N): ... arr[i, j]` (FORTRAN `do i=2,N`) | `for i in range(1, N): ... arr[i, j]` | The "skip first cell" idiom from FORTRAN. |
| `arr[1] = X` then loop `for j in range(N): arr[j] = ...` | `arr[0] = X` then `for j in range(N): arr[j] = ...` | "Special first row, then fill the rest" pattern. |

### Sentinel value change (highest-risk item)

In the current code, `neigh[i, j] == 0` means "slot j of cell i has no neighbour." Cell indices run `1..num_all_cells`, so `0` is safely outside that range. After converting to 0-based, valid cell indices are `0..num_all_cells - 1`, and `0` is a real cell. We must change the sentinel to `-1` everywhere `neigh` is initialised, assigned-to-clear, or compared-against-zero.

**All four sentinel sites** (verified by grep):
- `coreop2d.py:172` — `cls.neigh[:] = 0`  ⇒  `cls.neigh[:] = -1`
- `coreop2d.py:214` — `cls.neigh[i, j] = 0`  ⇒  `cls.neigh[i, j] = -1`
- `coreop2d.py:883` — `if cls.neigh[i, j] == 0`  ⇒  `if cls.neigh[i, j] == -1`
- `coreop2d.py:924` — `temp_neigh[jj, :] = 0`  ⇒  `temp_neigh[jj, :] = -1`
- `coreop2d.py:1122` — `temp_new_neigh[i, :] = 0`  ⇒  `temp_new_neigh[i, :] = -1`
- `coreop2d.py:1244` — `if cls.neigh[i, j] == 0`  ⇒  `if cls.neigh[i, j] == -1`
- `esclec.py` — every `if iii == 0` / `if ii == 0` / `if iiii == 0` inside `guardaveinsoff_2` (the "neighbour walk" in the OFF face-emission loops) refers to a value pulled from `core.neigh`. Every such comparison must become `== -1`.

`knots[]` and `num_neigh[]` use `0` as a real value (boolean false / count zero), **not** as a sentinel. Do not change those.

### `nca` semantics

`cls.nca` is the count-of-cells-allocated-so-far, used in `initialise_cell_positions` as both a count and a 1-based index. The current pattern is:

```python
cls.nca += 1
cls.neigh[i_centre, j] = cls.nca
cls.positions[cls.nca, 1] = x
```

After 0-based migration, change to:

```python
cls.neigh[i_centre, j] = cls.nca
cls.positions[cls.nca, 0] = x
cls.nca += 1
```

(i.e., write into the next free slot, then increment). Update the initialisation `cls.nca = 1` to `cls.nca = 1` if cell 0 is pre-populated (positions[0, :]=[0,0,1]) — `nca` then means "next free index." Verify all uses of `nca` (loops over `range(cls.nca)`, etc.) after the change.

### Parameter file lines

`esclec.read_param_file` reads 30 lines from the file and stores them at `parap[3..32, map]`. `set_params` reads `parap[1..32, map]`. The parameter-file format is 1-based by tradition (line 1 = first parameter). After migration:

- Allocation `parap = float_array((32, ma_max))` → `parap = np.zeros((32, ma_max))` (same size — 32 slots).
- Loop `for i in range(3, 33): parap[i, cls.map] = ...` → `for i in range(2, 32): parap[i, cls.map] = ...` (still reads 30 lines into 30 slots, now indexed 2..31).
- `set_params` indices `[1, imap], [2, imap], ..., [32, imap]` → `[0, imap], [1, imap], ..., [31, imap]`.
- `cls.map = 1` → `cls.map = 0`.
- Add a comment at the top of `read_param_file` explaining: "Parameter file lines 1..30 map to `parap[2..31]` (positions 0 and 1 are reserved for `temps` and `num_active_cells`, set elsewhere)." This makes the 1-based file format explicit.

### Slicing edge cases

- `cls.border[:, :, 4:5]` (1-based slice 4..5 inclusive ⇒ memory `[3:5]`) → `cls.border[:, :, 3:5]`.
- Fancy slices like `cv[ii, :]` don't change — `:` means "all" in both worlds.

---

## File structure

No new files. Final state:

- `indexing.py` — **deleted**.
- `vector.py` — uses plain `numpy.ndarray`; no `from indexing` import.
- `main.py` — uses builtin `range`; no `from indexing` import.
- `esclec.py` — uses `numpy.ndarray` and `numpy.zeros`; documented parameter-file mapping comment.
- `coreop2d.py` — uses `numpy.ndarray`; `neigh` sentinel is `-1`; `nca` semantics documented.
- `tests/test_simulator_smoke.py` — **new**, end-to-end regression test.
- `tests/golden/` — **new**, captured OFF reference output.

---

## Tasks

### Task 1: Pin Python and NumPy environment

**Files:**
- Modify: `pyproject.toml` (create if absent)

- [ ] **Step 1: Confirm Python version**

Run: `python --version`
Expected: Python 3.11.x or later (the code uses `from typing import Self`, which is 3.11+).

- [ ] **Step 2: Confirm numpy is available and pin version for reproducibility**

Run: `python -c "import numpy; print(numpy.__version__)"`
Expected: a version string (e.g. `2.0.x` or `1.26.x`). Note the version.

- [ ] **Step 3: Create `pyproject.toml` recording the pinned versions**

If `pyproject.toml` does not already exist, create it:

```toml
[project]
name = "toothmaker-translation"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = [
    "numpy>=1.26,<3.0",
]

[project.optional-dependencies]
dev = ["pytest>=7"]

[tool.pytest.ini_options]
testpaths = ["tests"]
```

- [ ] **Step 4: Commit**

```bash
git add pyproject.toml
git commit -m "chore: add pyproject.toml pinning Python 3.11 and numpy"
```

---

### Task 2: Verify the current simulator runs and capture golden output

**Files:**
- Create: `tests/golden/run/100_seal_.off` (and friends — 5 files)
- Create: `tests/golden/RUN_COMMAND.md`

- [ ] **Step 1: Run the simulator with the seal example**

Run:

```bash
mkdir -p .tmp/run-baseline
python main.py examples/seal.txt .tmp/run-baseline run 100 5
```

Expected: prints `Block 1/5 complete: ...` through `Block 5/5 complete: ...`, and produces `.tmp/run-baseline/100_run_.off`, `200_run_.off`, `300_run_.off`, `400_run_.off`, `500_run_.off`.

- [ ] **Step 2: If the run fails, stop and report**

If Python raises an exception, stop here and surface the error in the task report. Do not attempt fixes — the conversion plan assumes a working baseline. The user will decide whether to fix the bug before proceeding or to switch to a "no golden output, manual inspection" mode.

- [ ] **Step 3: Capture the OFF files as the golden reference**

Run:

```bash
mkdir -p tests/golden
cp .tmp/run-baseline/*.off tests/golden/
```

- [ ] **Step 4: Document the exact run command**

Write `tests/golden/RUN_COMMAND.md`:

```markdown
# Golden output

Captured from the pre-0-based-migration codebase by running:

    python main.py examples/seal.txt .tmp/run-baseline run 100 5

Output: 5 OFF files, one per save block (`100_run_.off` ... `500_run_.off`).

The files in this directory are the exact bytes the simulator produced before
the 0-based indexing migration. The post-migration simulator must reproduce them
byte-for-byte.
```

- [ ] **Step 5: Commit**

```bash
git add tests/golden/
git commit -m "test: capture pre-migration golden OFF output"
```

---

### Task 3: Add a regression smoke test

**Files:**
- Create: `tests/test_simulator_smoke.py`
- Create: `tests/__init__.py` (empty, if needed)

- [ ] **Step 1: Write the test**

Create `tests/test_simulator_smoke.py`:

```python
"""End-to-end regression test: simulator output must match captured golden files."""
import filecmp
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
GOLDEN = REPO_ROOT / "tests" / "golden"


def test_seal_run_matches_golden(tmp_path):
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

    expected_files = sorted(p.name for p in GOLDEN.glob("*.off"))
    actual_files = sorted(p.name for p in out_dir.glob("*.off"))
    assert actual_files == expected_files, (
        f"file list mismatch: golden={expected_files}, actual={actual_files}"
    )

    for name in expected_files:
        assert filecmp.cmp(GOLDEN / name, out_dir / name, shallow=False), (
            f"{name} differs from golden output"
        )
```

- [ ] **Step 2: Run the test against the current (pre-migration) code**

Run: `pytest tests/test_simulator_smoke.py -v`
Expected: PASS (the simulator currently produces exactly the bytes we just captured).

- [ ] **Step 3: Commit**

```bash
git add tests/test_simulator_smoke.py
git commit -m "test: add end-to-end smoke test against golden output"
```

---

### Task 4: Convert `vector.py`

`vector.py` is 39 lines and only uses `float_array` for the four-component (1-based-indexed) shim — not for storage. After migration it imports nothing from `indexing` and accepts plain numpy arrays.

**Files:**
- Modify: `vector.py`

- [ ] **Step 1: Rewrite `vector.py`**

Replace the entire file contents with:

```python
import math
import numpy as np

def magnitude(x: float, y: float, z: float) -> float:
    return math.sqrt(x*x + y*y + z*z)

def a_magnitude(v) -> float:
    if isinstance(v, np.ndarray):
        return magnitude(v[0], v[1], v[2])
    return magnitude(v[0], v[1], v[2])

def distance_between(
    x1: float, y1: float, z1: float,
    x2: float, y2: float, z2: float,
) -> float:
    return magnitude(x1-x2, y1-y2, z1-z2)

def a_distance_between(u, v) -> float:
    return distance_between(u[0], u[1], u[2], v[0], v[1], v[2])

def cross_product(
    x1: float, y1: float, z1: float,
    x2: float, y2: float, z2: float,
):
    return np.array([
        y1 * z2 - z1 * y2,
        z1 * x2 - x1 * z2,
        x1 * y2 - y1 * x2,
    ])

def a_cross_product(u, v) -> np.ndarray:
    return cross_product(u[0], u[1], u[2], v[0], v[1], v[2])
```

- [ ] **Step 2: Do not run the smoke test yet**

The other files still expect 1-based wrappers; the smoke test will fail until the whole migration is done. This is expected. Continue.

- [ ] **Step 3: Commit**

```bash
git add vector.py
git commit -m "refactor: convert vector.py to plain numpy (0-based)"
```

---

### Task 5: Convert `main.py`

**Files:**
- Modify: `main.py`

- [ ] **Step 1: Update imports and the iteration loop**

Replace `main.py` with:

```python
import sys
import argparse
import os
from coreop2d import Coreop2d
from esclec import Esclec


def main():
    parser = argparse.ArgumentParser(
        description='ToothMaker tooth morphogenesis simulator (Python translation)'
    )
    parser.add_argument('input_file', help='Parameter file path')
    parser.add_argument('output_folder', help='Folder for output files')
    parser.add_argument('output_name', help='Base name for output files')
    parser.add_argument('iterations', type=int,
                        help='Number of iterations per save block')
    parser.add_argument('save_blocks', type=int,
                        help='Number of save blocks')
    args = parser.parse_args()

    cac = args.input_file
    caufolder = args.output_folder
    cau = args.output_name
    iteration_total = args.iterations
    sstep = args.save_blocks

    core = Coreop2d
    io = Esclec

    core.max_z_layers = 4
    core.initial_conditions()
    io.initialize_from_parameter_file(core, cac)

    prev_num_active_cells = core.num_active_cells
    core.allocate_initial_state()
    io.set_params(core, 0)
    core.num_active_cells = prev_num_active_cells

    core.initact()

    os.makedirs(caufolder, exist_ok=True)

    for iti in range(1, abs(sstep) + 1):
        iter_label = str(iti * iteration_total)

        nff = os.path.join(caufolder, iter_label + '_' + cau)
        nff = nff.replace(' ', '_')

        nfoff = nff + '_.off'
        nfes = nff + '_.txt'

        core.temps = 0
        nt = iteration_total
        core.iteration(nt)

        with open(nfoff, 'w') as f:
            io.guardaveinsoff_2(core, core.neigh, f)

        print(f'Block {iti}/{abs(sstep)} complete: {nfoff}')


if __name__ == '__main__':
    main()
```

Key changes:
- Removed `from indexing import *`.
- `io.set_params(core, 1)` → `io.set_params(core, 0)` (parameter set index becomes 0-based).
- The save-block `for iti in range(1, abs(sstep) + 1)` is **kept as-is** (this loop runs `1..sstep` for human-friendly labels — it has nothing to do with array indexing).

- [ ] **Step 2: Commit**

```bash
git add main.py
git commit -m "refactor: convert main.py to builtin range and 0-based set_params"
```

---

### Task 6: Convert `esclec.py`

**Files:**
- Modify: `esclec.py`

- [ ] **Step 1: Replace imports**

Change line 2 from `from indexing import *` to `import numpy as np`.

- [ ] **Step 2: Convert class-level allocations**

Lines 8-11 currently:

```python
positionsp = float_array((1000, 3, ma_max))
parap = float_array((32, ma_max))
param_names = str_array(32)
knotsp = int_array((1000, ma_max))
```

Replace with:

```python
positionsp = np.zeros((1000, 3, ma_max))
parap = np.zeros((32, ma_max))
param_names = np.empty(32, dtype=object)
knotsp = np.zeros((1000, ma_max), dtype=np.int32)
```

Also change the type annotation `ma: float_array` to `ma: np.ndarray`.

- [ ] **Step 3: Convert `read_param_file`**

Replace the body with:

```python
@classmethod
def read_param_file(cls, filepath):
    # Parameter file format (1-based for human readability):
    #   lines 1..30 map to parap[2..31, map].
    #   parap[0..1, map] are reserved for temps and num_active_cells,
    #   set elsewhere (set_params reads them but they default to 0).
    cls.read_failed = 0
    try:
        with open(filepath, 'r') as f:
            for i in range(2, 32):  # 30 slots: indices 2..31
                line = f.readline()
                if not line:
                    cls.read_failed = 1
                    return
                parts = line.split()
                if len(parts) < 2:
                    cls.read_failed = 1
                    return
                cls.parap[i, cls.map] = float(parts[0])
                cls.param_names[i] = parts[1]
    except (IOError, ValueError):
        print("reading error for", filepath)
        cls.read_failed = 1
```

- [ ] **Step 4: Convert `set_params`**

Subtract 1 from every numeric index in `parap[N, imap]`:

```python
@classmethod
def set_params(cls, core: Coreop2d, imap: int):
    core.temps =            cls.parap[0, imap]
    core.num_active_cells = cls.parap[1, imap]
    core.Egr = cls.parap[2, imap]
    core.Mgr = cls.parap[3, imap]
    core.Rep = cls.parap[4, imap]
    # parap[5, imap] reserved (file slot 6 unused)
    core.Adh = cls.parap[6, imap]
    core.Act = cls.parap[7, imap]
    core.Inh = cls.parap[8, imap]
    # parap[9, imap] reserved
    core.Sec = cls.parap[10, imap]
    # parap[11, imap] reserved
    core.Da  = cls.parap[12, imap]
    core.Di  = cls.parap[13, imap]
    core.Ds  = cls.parap[14, imap]
    # parap[15, imap] reserved
    core.Int = cls.parap[16, imap]
    core.Set = cls.parap[17, imap]
    core.Boy = cls.parap[18, imap]
    core.Dff = cls.parap[19, imap]
    core.Bgr = cls.parap[20, imap]
    core.Abi = cls.parap[21, imap]
    core.Pbi = cls.parap[22, imap]
    core.Bbi = cls.parap[23, imap]
    core.Lbi = cls.parap[24, imap]
    core.Rad = cls.parap[25, imap]
    core.Deg = cls.parap[26, imap]
    core.Dgr = cls.parap[27, imap]
    core.Ntr = cls.parap[28, imap]
    core.Bwi = cls.parap[29, imap]
    core.ina = cls.parap[30, imap]
    core.umgr = cls.parap[31, imap]
```

- [ ] **Step 5: Convert `initialize_from_parameter_file`**

Change `cls.map = 1` to `cls.map = 0`.

- [ ] **Step 6: Convert `guardaveinsoff_2`**

Apply the cheat sheet to the entire method:

- `c = float_array(4); mic = float_array(4)` → `c = np.zeros(4); mic = np.zeros(4)`.
- `face = int_array((core.num_active_cells * 20, 5))` → `face = np.zeros((core.num_active_cells * 20, 5), dtype=np.int32)`.
- `cls.ma = float_array(core.num_active_cells)` → `cls.ma = np.zeros(core.num_active_cells)`.
- Every `iii == 0`, `ii == 0`, `iiii == 0` (sentinel checks pulled from `core.neigh`) → `== -1`.
- Every `iii > core.num_active_cells` → `>= core.num_active_cells` (boundary shifts).
- The vertex-write loop `for i in range(core.num_active_cells): ... x = core.positions[i, 1]; y = core.positions[i, 2]; z = core.positions[i, 3]; f.write(... c[1]..c[4] ...)` → loop variable `i` now natively iterates `0..N-1`; rewrite as `x = core.positions[i, 0]; y = core.positions[i, 1]; z = core.positions[i, 2]; f.write(f"{x} {y} {z} {c[0]} {c[1]} {c[2]} {c[3]}\n")`.
- The face-write lines `f.write(f"3 {i - 1} {ii - 1} {iii - 1}\n")` → `f.write(f"3 {i} {ii} {iii}\n")` (Python is now 0-based, matching the OFF format).
- Same for the quad-face line: `f.write(f"4 {i - 1} {ii - 1} {iii2 - 1} {iiii2 - 1}\n")` → `f.write(f"4 {i} {ii} {iii2} {iiii2}\n")`.
- Inner loops over `range(core.nv_max)` and `range(core.num_active_cells)` are correct as-is once we are using the builtin `range` (now they iterate `0..N-1`).

- [ ] **Step 7: Convert `mat`, `get_rainbow`, `get_rainbow_knot`**

- `mat`: loop `for i in range(core.num_active_cells)` is correct under builtin range; `cls.ma[:] = 0` is fine; the `va_max` and `va_min` comprehensions work as-is.
- `get_rainbow` and `get_rainbow_knot`: every `c[1], c[2], c[3], c[4] = ...` becomes `c[0], c[1], c[2], c[3] = ...`.

- [ ] **Step 8: Commit**

```bash
git add esclec.py
git commit -m "refactor: convert esclec.py to 0-based indexing and plain numpy"
```

---

### Task 7: Convert `coreop2d.py`

This is the largest file (1385 lines). The conversion is mechanical — apply the cheat sheet — but the volume makes mistakes likely. Work through it method by method, committing after each logical block. After every commit, re-grep for stragglers (see Step 0).

**Files:**
- Modify: `coreop2d.py`

- [ ] **Step 0: Establish a per-method work-list and verification grep**

List the methods in order:

```bash
/usr/bin/grep -n "    def \|    @classmethod" coreop2d.py
```

After each per-method commit, run this verification grep — it must return zero hits when the file is fully converted:

```bash
/usr/bin/grep -nE "from indexing|float_array\(|int_array\(|bool_array\(|str_array\(" coreop2d.py
```

(While the conversion is in progress this will return hits in unconverted methods; that is fine.)

- [ ] **Step 1: Replace the import line**

Line 5: `from indexing import *` → `import numpy as np`.

Do **not** commit yet — the file will not run until later steps are also done. The smoke test will be re-run only at the end of Task 8.

- [ ] **Step 2: Convert type annotations (top of class)**

Lines 9-21 currently use `float_array` / `int_array` annotations. Replace each with `np.ndarray`:

```python
positions:      np.ndarray
border:         np.ndarray
neigh:          np.ndarray
knots:          np.ndarray
num_neigh:      np.ndarray
q3d:            np.ndarray
diff_state:     np.ndarray
forces:         np.ndarray
force_snapshot: np.ndarray
px: np.ndarray
py: np.ndarray
pz: np.ndarray
```

Also `border_cell_list: int_array` and `m_map: int_array` lower in the class → `np.ndarray`.

- [ ] **Step 3: Convert `initial_conditions`**

No array indexing in this method. No changes.

- [ ] **Step 4: Convert `allocate_initial_state` (lines ~137-258)**

Apply the cheat sheet. Specific landmark changes:

- All `float_array((...))` / `int_array((...))` → `np.zeros((...), dtype=...)`.
- `cls.neigh[:] = 0` → `cls.neigh[:] = -1` **(sentinel change).**
- `cls.positions[1, 1] = 0`, `cls.positions[1, 2] = 0`, `cls.positions[1, 3] = 1` → `cls.positions[0, 0] = 0`, `cls.positions[0, 1] = 0`, `cls.positions[0, 2] = 1`.
- `cls.nca = 1` → `cls.nca = 1` (still means "next free index"; cell 0 is now the populated centre).
- The mesh-init loop `for i_centre in range(cls.num_active_cells)` reads `cls.positions[i_centre, 1]`, `[i_centre, 2]` → `[i_centre, 0]`, `[i_centre, 1]`. The loop bounds (`range(N)`) are correct under builtin range.
- The `range(2, cls.num_active_cells)` loops → `range(1, cls.num_active_cells)`.
- The `range(2, cls.num_active_cells)` zero-out loop body `cls.neigh[i, j] = 0` → `cls.neigh[i, j] = -1` **(sentinel)**.
- The compaction loop `for jj in range(j, cls.nv_max-1): cls.neigh[i, jj] = cls.neigh[i, jj+1]` — `j` is itself a shifted loop variable, so the body is unchanged; just verify the `range(j, cls.nv_max-1)` becomes `range(j, cls.nv_max-1)` (already correct because both endpoints are shifted variables, not literals).
- `cls.positions.inner = cls.positions.inner.round(...)` → `cls.positions = cls.positions.round(decimals=14)`.
- The `for i in range(cls.num_active_cells): for j in range(3): if abs(cls.positions[i, j]) < 1e-14:` loop is correct under builtin range (no literal indices to shift).
- The `cls.num_active_cells - i + 1` mirror calculation: under 1-based, valid indices are `1..N`, so `N - i + 1` is the mirrored 1-based index. Under 0-based, valid indices are `0..N-1`, so the mirror should be `N - 1 - i`. Update both occurrences (lines ~230 and ~236).
- `cls.num_neigh[1] = 6` → `cls.num_neigh[0] = 6`.
- `cls.border[:, :, 4:5] = cls.la` → `cls.border[:, :, 3:5] = cls.la` (1-based slice 4..5 inclusive ⇒ 0-based slice 3..5 exclusive).
- `cls.first_border_cell = 6 * (cls.Rad - 1) + 1` → `cls.first_border_cell = 6 * (cls.Rad - 1)` (this stores a 1-based cell index; under 0-based it loses the +1).
- The two trailing `for i in range(cls.Rad)` loops: `cls.m_map[i] = i` and `cls.border_cell_list[i] = cls.first_border_cell // 2 + i` work as-is under builtin range — `i` runs `0..Rad-1` and the assigned values shift correspondingly.

- [ ] **Step 5: Convert `initialise_cell_positions` (lines ~260-281)**

Refactor the increment-then-write pattern:

```python
@classmethod
def initialise_cell_positions(cls, i_centre: int, x: float, y: float, j: int, jj: int):
    for i in range(cls.nca):
        if i == i_centre: continue
        x_approx = math.isclose(x, cls.positions[i, 0], abs_tol=1e-6)
        y_approx = math.isclose(y, cls.positions[i, 1], abs_tol=1e-6)
        if x_approx and y_approx:
            for ii in range(cls.nv_max):
                if cls.neigh[i_centre, ii] == i: return
            cls.neigh[i_centre, j] = i
            cls.neigh[i, jj] = i_centre
            cls.num_neigh[i] += 1
            cls.num_neigh[i_centre] += 1
            return
    cls.neigh[i_centre, j] = cls.nca
    cls.neigh[cls.nca, jj] = i_centre
    cls.positions[cls.nca, 0] = x
    cls.positions[cls.nca, 1] = y
    cls.positions[cls.nca, 2] = 1
    cls.num_neigh[i_centre] += 1
    cls.num_neigh[cls.nca] += 1
    cls.nca += 1
```

Note `j` and `jj` here are **call-site arguments** (1, 2, 3, 4, 5, 6 in `allocate_initial_state`). After migration, the call sites must pass `0, 1, 2, 3, 4, 5` instead. **Update those six call sites** (lines 193-198 in `allocate_initial_state`):

```python
cls.initialise_cell_positions(i_centre, x+cls.csu*cls.la, y+cls.ssu*cls.la, 0, 3)
cls.initialise_cell_positions(i_centre, x+cls.csd*cls.la, y+cls.ssd*cls.la, 1, 4)
cls.initialise_cell_positions(i_centre, x+cls.cst*cls.la, y+cls.sst*cls.la, 2, 5)
cls.initialise_cell_positions(i_centre, x+cls.csq*cls.la, y+cls.ssq*cls.la, 3, 0)
cls.initialise_cell_positions(i_centre, x+cls.csc*cls.la, y+cls.ssc*cls.la, 4, 1)
cls.initialise_cell_positions(i_centre, x+cls.css*cls.la, y+cls.sss*cls.la, 5, 2)
```

- [ ] **Step 6: Commit progress so far**

```bash
git add coreop2d.py
git commit -m "refactor: convert coreop2d setup methods to 0-based"
```

- [ ] **Step 7: Convert `calculate_margins` (lines ~283-394)**

The two reverse-direction inner loops are the trickiest:

- `for jj in range(j-1, 1, -1)` (overridden ⇒ `j-1, j-2, ..., 1`) → `for jj in range(j-1, -1, -1)` (builtin ⇒ `j-1, j-2, ..., 0`). Note the new stop `-1` not `0` because builtin `range`'s stop is exclusive, and we want `0` included.

  Wait — that gives **one extra iteration** compared to the original. Recompute: the original yielded values `j-1` down to `1` (inclusive), which is `j-1` values. After the cell-index shift by -1, those values should be `j-2` down to `0`, which is also `j-1` values: `range(j-2, -1, -1)`. Use that.

- `for jj in range(cls.nv_max, j+1, -1)` (overridden ⇒ `nv_max, nv_max-1, ..., j+1`) → `range(cls.nv_max-1, j-1, -1)` (builtin ⇒ `nv_max-1, ..., j`).

- `for jj in range(1, jjj, -1)` (overridden ⇒ `1` only, since stop>start with negative step yields nothing useful — confirm by re-reading the surrounding context whether this loop is intended to iterate at all). After conversion: `range(0, jjj-1, -1)` (builtin) yields nothing for `jjj >= 1`. **Verify this loop's intent against the FORTRAN source `13.f90` before converting.**

- All `cls.positions[ii, 1]`, `[ii, 2]`, `[ii, 3]` → `[ii, 0]`, `[ii, 1]`, `[ii, 2]`.
- `cls.neigh[i, jj] != 0` checks: these are checking the sentinel. → `!= -1`.
  **Wait — re-check.** The current code uses `0` as the "no neighbour" sentinel only at the four sites listed in the Sentinel section. Verify: at line ~297 `if cls.neigh[i, jj] != 0` — this is checking "if the slot is occupied", so yes, sentinel check. → `!= -1`.

- [ ] **Step 8: Continue method by method through the rest of the file**

For each remaining method (`initact`, `iteration`, the diffusion/reaction routines, the mechanical-force routines, the cell-division routines, etc.):

1. Identify all `range(...)` loops and adjust per the cheat sheet.
2. Identify all numeric-literal index expressions (`arr[1]`, `arr[1, 1]`, `arr[i, 1]`, etc.) and shift literals by -1.
3. Identify all `arr.inner` accesses and replace with `arr` directly.
4. Identify any `cls.neigh[X] == 0`, `cls.neigh[X] = 0`, `cls.neigh[X] != 0` — change `0` to `-1` (sentinel).
5. Identify any `cls.neigh[X] > cls.num_active_cells` boundary checks — change `>` to `>=`.
6. Leave `cls.knots[X] == 0`, `cls.knots[X] == 1`, `cls.num_neigh[X] = 0` etc. **alone** — these are real values, not sentinels.

Commit after every two or three methods to keep diffs reviewable:

```bash
git add coreop2d.py
git commit -m "refactor: convert <method names> to 0-based"
```

- [ ] **Step 9: Final verification grep for stragglers**

Run:

```bash
/usr/bin/grep -nE "from indexing|float_array\(|int_array\(|bool_array\(|str_array\(|\.inner" coreop2d.py
```

Expected: **no output**. If anything matches, fix it before continuing.

Also:

```bash
/usr/bin/grep -nE "neigh\[.*\]\s*==\s*0|neigh\[.*\]\s*=\s*0|neigh\[:\]\s*=\s*0" coreop2d.py
```

Expected: **no output**. Every neigh sentinel must now be `-1`.

- [ ] **Step 10: Commit any final cleanups**

```bash
git add coreop2d.py
git commit -m "refactor: complete coreop2d.py 0-based conversion"
```

---

### Task 8: Run the regression smoke test

**Files:**
- (none modified — verification only)

- [ ] **Step 1: Run the smoke test**

Run: `pytest tests/test_simulator_smoke.py -v`

Expected: PASS — the post-migration simulator produces byte-identical OFF output to the golden files.

- [ ] **Step 2: If the test fails, diff and diagnose**

If files differ, run:

```bash
diff <(head -50 tests/golden/100_run_.off) <(head -50 .tmp/run-baseline/100_run_.off)
```

(Recapture the live output first by re-running the simulator manually.)

The most likely failure modes, in order:
1. **Sentinel miss** — somewhere `neigh[...] == 0` was not changed to `== -1`. Symptoms: missing or extra OFF faces, garbled face indices.
2. **Off-by-one in `range` conversion** — a loop bound was shifted incorrectly. Symptoms: missing first or last cell, extra cell.
3. **Boundary check still uses `>` not `>=`** — a cell on the boundary is misclassified. Symptoms: faces along the rim differ.
4. **Mirror calculation `N - i + 1` not converted to `N - 1 - i`** — geometry mirrored to the wrong row. Symptoms: positions mirrored inconsistently.
5. **`positions[i, 1/2/3]` not all converted to `[i, 0/1/2]`** — one coordinate wrong. Symptoms: vertex coordinates differ in one column.

Do not commit any "fix" until the test passes.

- [ ] **Step 3: When green, commit any fixes**

```bash
git add -p
git commit -m "fix: <specific issue> after 0-based migration"
```

(Use `gitops.sh add` per your wrappers if `git add -p` is awkward.)

---

### Task 9: Delete `indexing.py`

**Files:**
- Delete: `indexing.py`

- [ ] **Step 1: Confirm nothing imports from indexing**

Run:

```bash
/usr/bin/grep -rn "from indexing\|import indexing" .
```

Expected: no matches in any `.py` file.

- [ ] **Step 2: Delete the file**

Run: `git rm indexing.py`

- [ ] **Step 3: Re-run the smoke test**

Run: `pytest tests/test_simulator_smoke.py -v`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git commit -m "refactor: remove obsolete 1-based indexing wrappers"
```

---

### Task 10: Open the PR

**Files:**
- (none — git/gh workflow only)

- [ ] **Step 1: Push the branch**

Branch should be named `refactor/zero-based-indexing` (created at the start of this work).

```bash
gitops.sh push-new
```

- [ ] **Step 2: Open a PR with `main` as the base**

```bash
ghops.sh pr create --base main --title "refactor: convert codebase to 0-based indexing" --body "$(cat <<'EOF'
## Summary
- Removes the 1-based `indexing.py` wrappers (an aid for translating from FORTRAN).
- Converts `vector.py`, `main.py`, `esclec.py`, and `coreop2d.py` to native numpy 0-based indexing.
- Changes the `neigh[]` sentinel from `0` to `-1` (since `0` is now a valid cell index).
- Adds an end-to-end smoke test that runs the simulator on `examples/seal.txt` and compares output against captured golden OFF files.

## Test plan
- [ ] `pytest tests/test_simulator_smoke.py -v` passes
- [ ] OFF output for `examples/seal.txt` is byte-identical to the pre-migration baseline

🤖 Generated with [Claude Code](https://claude.com/claude-code)
EOF
)"
```

- [ ] **Step 3: Report the PR URL**

---

## Self-review

**Spec coverage:** the goal was to remove `indexing.py` and switch to native 0-based indexing while preserving simulator output. Tasks 1-3 establish a regression baseline; Tasks 4-7 convert the four source files; Task 8 verifies; Task 9 deletes the wrapper; Task 10 opens the PR. Covered.

**Placeholder scan:** the cheat sheet, the sentinel section, and each per-file step state explicit pattern transformations and exact code where the engineer needs it. Task 7 Step 8 ("continue method by method") is the one place I rely on the engineer to apply the cheat sheet rather than spelling out every method — but `coreop2d.py` is 1385 lines and listing every method explicitly would inflate the plan beyond usefulness. The cheat sheet, the verification greps in Step 9, and the smoke test together prevent this from being a blank check.

**Type consistency:** `np.ndarray` is used consistently for type annotations; `np.zeros` for allocations (matching the wrappers' pre-zeroed behaviour); `-1` consistently as the new sentinel. The `nca` semantics change (write-then-increment vs increment-then-write) is described once in the cheat sheet and once in Task 7 Step 5 with a worked example, and the call-site updates are spelled out.

**Known risks the plan flags:**
- The `range(1, jjj, -1)` loop in `calculate_margins` is suspicious in the original and the plan says to verify against the FORTRAN before blindly converting.
- The "no working baseline" branch (Task 2 Step 2) — if the pre-migration code does not run, the engineer must stop and report rather than silently continue without a regression check.
