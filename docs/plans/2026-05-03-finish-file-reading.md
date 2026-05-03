# Finish File Reading and Main Loop Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Get `python main.py examples/seal.txt out run 100 5` to run end-to-end and produce 5 OFF files, while staying on the existing 1-based `indexing.py` wrappers. Capture the resulting OFF files as the golden output that will be used to validate the subsequent 0-based migration.

**Architecture:** Two known bugs block the end-to-end run today: (1) `esclec.read_param_file` iterates 31 times over a 30-line file, and (2) `main.main` saves `prev_num_active_cells` before `allocate_initial_state` populates it. Fix each, retry, and resolve any further faults that surface. Stop at the first byte-stable OFF output; do not over-engineer or refactor.

**Tech Stack:** Python 3.12, NumPy 2.x, pytest. Project venv at `.venv/` (already created with numpy and pytest installed).

---

## Background

The current state of `feat/off-output` (now merged to `main`) implements parameter-file reading, the main simulation loop, and OFF mesh output. The user has confirmed the file reading "is not yet working." Investigation today confirms:

- `examples/seal.txt` contains 30 data lines (Egr through umgr), aligned with FORTRAN `do i=3,32` reading into `parap(3:32, map)`.
- `esclec.read_param_file` uses the overridden `range(3, 33)` which iterates `3..33` inclusive — **31 iterations**, one more than the FORTRAN. The 31st `readline()` returns empty/short and trips `read_failed = 1`.
- Because `read_failed == 1`, `set_params` is skipped, so `Rad`, `Egr`, etc. never get written.
- `main.py` then crashes on `core.num_active_cells` (never set).

Even after fixing the off-by-one, the next problem is in `main.py`'s sequencing: it saves `prev_num_active_cells = core.num_active_cells` *before* calling `allocate_initial_state`, which is the only thing that derives `num_active_cells` from `Rad`. The intended FORTRAN sequence (lines 1830-1837 in `13.f90`) is:

```fortran
call initialconditions
call initialize_from_parameter_file        ! sets Rad (and clobbers num_active_cells = parap(2) = 0)
call allocateinitialstate                  ! sets num_active_cells = 3*(Rad-1)*Rad + 1

prev_num_active_cells = num_active_cells   ! save the Rad-computed value
call setparams(1)                          ! clobbers num_active_cells back to 0
num_active_cells = prev_num_active_cells   ! restore

call initact
```

The Python translation got the save/allocate steps in the wrong order.

These two bugs are in scope. **Other latent bugs may surface once the run gets further** — fix them as they appear, but only enough to unblock the run. Do not refactor for cleanliness; that is the next milestone (the 0-based migration).

---

## File structure

No new files in the source tree. Final state:

- `esclec.py` — read loop fixed (30 iterations, not 31).
- `main.py` — initialisation order fixed; matches FORTRAN.
- (any other source file touched only if a third blocking bug surfaces).
- `tests/test_simulator_smoke.py` — **new**, end-to-end test asserting the simulator exits 0 and emits 5 OFF files matching the captured golden bytes.
- `tests/golden/100_run_.off` ... `500_run_.off` — **new**, captured reference output.
- `tests/golden/RUN_COMMAND.md` — **new**, documents how the goldens were captured.
- `pyproject.toml` — **new**, pins Python and numpy versions.

---

## Tasks

### Task 1: Pin the environment

**Files:**
- Create: `pyproject.toml`

- [ ] **Step 1: Confirm the venv already exists**

Run: `ls .venv/bin/python && .venv/bin/python -c "import numpy, pytest; print(numpy.__version__, pytest.__version__)"`
Expected: prints `2.4.4 <pytest-version>` (or whatever versions are installed).

If the venv is missing, run `python -m venv .venv && .venv/bin/pip install numpy pytest` first.

- [ ] **Step 2: Create `pyproject.toml`**

Write `pyproject.toml`:

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

- [ ] **Step 3: Commit**

```bash
gitops.sh add pyproject.toml
gitops.sh commit "chore: add pyproject.toml pinning Python 3.11 and numpy"
```

---

### Task 2: Fix the off-by-one in `read_param_file`

**Files:**
- Modify: `esclec.py:28`

- [ ] **Step 1: Read the current loop**

Run: `sed -n '23,42p' esclec.py`

Expected: shows the `read_param_file` body with `for i in range(3, 33):` at line 28.

- [ ] **Step 2: Fix the loop bound**

Edit `esclec.py:28` — change `for i in range(3, 33):  # indices 3 to 32 inclusive` to:

```python
            for i in range(3, 32):  # indices 3..32 inclusive (30 lines, mirroring FORTRAN do i=3,32)
```

(The overridden `range` is end-inclusive, so `range(3, 32)` iterates 30 times: 3..32.)

- [ ] **Step 3: Verify the param file now loads**

Run:

```bash
.venv/bin/python -c "
from coreop2d import Coreop2d
from esclec import Esclec
core = Coreop2d
core.max_z_layers = 4
core.initial_conditions()
Esclec.initialize_from_parameter_file(core, 'examples/seal.txt')
print('read_failed:', Esclec.read_failed)
print('Egr =', core.Egr, '(expected 0.5)')
print('Rad =', core.Rad, '(expected 4.0)')
print('umgr =', core.umgr, '(expected 0.0)')
"
```

Expected:

```
read_failed: 0
Egr = 0.5 (expected 0.5)
Rad = 4.0 (expected 4.0)
umgr = 0.0 (expected 0.0)
```

If `read_failed` is still 1, stop and inspect — there may be a second issue (file encoding, trailing-blank-line handling, or `parts < 2`).

- [ ] **Step 4: Commit**

```bash
gitops.sh add esclec.py
gitops.sh commit "fix(esclec): read 30 param lines, not 31 (off-by-one in read_param_file)"
```

---

### Task 3: Fix the initialisation order in `main.main`

**Files:**
- Modify: `main.py:31-40`

- [ ] **Step 1: Read the current sequence**

Run: `sed -n '28,42p' main.py`

Expected: shows

```python
    core.max_z_layers = 4
    core.initial_conditions()
    io.initialize_from_parameter_file(core, cac)

    prev_num_active_cells = core.num_active_cells
    core.allocate_initial_state()
    io.set_params(core, 1)
    core.num_active_cells = prev_num_active_cells
```

- [ ] **Step 2: Reorder to match FORTRAN**

Edit `main.py` — swap the order so `allocate_initial_state` runs *before* the save:

```python
    core.max_z_layers = 4
    core.initial_conditions()
    io.initialize_from_parameter_file(core, cac)

    core.allocate_initial_state()
    prev_num_active_cells = core.num_active_cells
    io.set_params(core, 1)
    core.num_active_cells = prev_num_active_cells
```

(The two changed lines: `prev_num_active_cells = core.num_active_cells` moves *after* `core.allocate_initial_state()`.)

- [ ] **Step 3: Verify the sequence runs without an `AttributeError`**

Run:

```bash
.venv/bin/python -c "
from coreop2d import Coreop2d
from esclec import Esclec
core = Coreop2d
core.max_z_layers = 4
core.initial_conditions()
Esclec.initialize_from_parameter_file(core, 'examples/seal.txt')
core.allocate_initial_state()
prev_num_active_cells = core.num_active_cells
Esclec.set_params(core, 1)
core.num_active_cells = prev_num_active_cells
print('num_active_cells =', core.num_active_cells)
print('num_all_cells =', core.num_all_cells)
print('Rad =', core.Rad)
"
```

Expected: prints sensible cell counts. With `Rad = 4.0`, `num_all_cells = 3*4*5 + 1 = 61` and `num_active_cells = 3*3*4 + 1 = 37`.

- [ ] **Step 4: Commit**

```bash
gitops.sh add main.py
gitops.sh commit "fix(main): allocate state before saving prev_num_active_cells (mirror FORTRAN)"
```

---

### Task 4: Run the simulator end-to-end and triage further failures

**Files:**
- (no commits in this task unless additional bugs surface)

- [ ] **Step 1: Try the run**

Run:

```bash
mkdir -p .tmp/run-baseline
.venv/bin/python main.py examples/seal.txt .tmp/run-baseline run 100 5 2>&1 | tail -30
```

Expected (best case): prints `Block 1/5 complete: ...` through `Block 5/5 complete: ...` and produces five `.off` files in `.tmp/run-baseline/`.

- [ ] **Step 2: If the run completes, skip to Task 5**

Confirm: `ls -la .tmp/run-baseline/*.off` shows exactly 5 files. Skip to Task 5.

- [ ] **Step 3: If the run fails, triage the failure**

If Python raises an exception, capture the traceback. The most likely failure points, with diagnostics:

**3a. NumPy deprecation/error in `indexing.py`:** numpy 2.x removed several aliases. If you see `AttributeError: module 'numpy' has no attribute 'bool'` or similar, fix the wrapper in `indexing.py` (e.g., `np.bool_` is still valid; some others need replacement). Commit as `fix(indexing): numpy 2.x compatibility`.

**3b. Index-out-of-bounds in `coreop2d`:** the wrappers translate every index by -1, so an off-by-one in a literal slice (e.g. `arr[N:M]` where N or M was meant to be 1-based) shows up as a numpy `IndexError`. Read the offending line, compare to FORTRAN `13.f90` at the corresponding subroutine, and adjust. Commit each fix individually.

**3c. Type mismatch (float vs int):** `cls.parap` is float-typed but several values become array indices or loop bounds. The pattern is `int(cls.parap[N, imap])` or `int(core.Rad)`. Add the explicit `int()` cast where needed. Commit as `fix(<file>): cast float param to int for indexing`.

**3d. `set_params` writes integer params as float:** `num_active_cells = parap[2, imap]` becomes a numpy float64. In numpy 2.x this can no longer be used as an array index without explicit casting. If this is the failure, change `set_params` to cast `temps`, `num_active_cells`, and `Rad` to `int`:

```python
core.temps = int(cls.parap[1, imap])
core.num_active_cells = int(cls.parap[2, imap])
...
core.Rad = int(cls.parap[26, imap])
```

(Other params stay float.)

**3e. `int_array(buffer=...)` with a list of floats:** if the buffer pathway raises, inspect the call site.

**3f. Crash in `coreop2d.iteration` past the first iteration:** this means the math is running but a stateful bug surfaces later. Fix only if it blocks the 5-block run; otherwise note and continue.

For each fix:
1. Make the minimal change.
2. Re-run the simulator command from Step 1.
3. If it gets further, iterate. If not, re-triage.
4. Commit each fix individually with a descriptive message.

- [ ] **Step 4: Loop until 5 OFF files are produced**

Repeat Step 3 until `.tmp/run-baseline/*.off` shows 5 files. Do not move on without that.

- [ ] **Step 5: Sanity-check the OFF output**

Run:

```bash
head -3 .tmp/run-baseline/100_run_.off
wc -l .tmp/run-baseline/*.off
```

Expected: each file starts with `COFF`, then a header line `<num_vertices> <num_faces> 0`, then has `num_vertices + num_faces + 2` total lines. Sizes should be roughly consistent block-to-block (cell count grows over time as cells divide; `500_run_.off` should have more vertices than `100_run_.off`).

If output looks malformed (zero faces, garbage numbers), stop and triage — the smoke test would lock in broken output as "golden."

---

### Task 5: Capture golden output

**Files:**
- Create: `tests/golden/100_run_.off`
- Create: `tests/golden/200_run_.off`
- Create: `tests/golden/300_run_.off`
- Create: `tests/golden/400_run_.off`
- Create: `tests/golden/500_run_.off`
- Create: `tests/golden/RUN_COMMAND.md`

- [ ] **Step 1: Copy the OFF files into `tests/golden/`**

Run:

```bash
mkdir -p tests/golden
cp .tmp/run-baseline/*.off tests/golden/
ls tests/golden/
```

Expected: lists the 5 `.off` files.

- [ ] **Step 2: Document the capture command**

Write `tests/golden/RUN_COMMAND.md`:

```markdown
# Golden simulator output

Captured from the post-bugfix, pre-0-based-migration codebase by running:

    .venv/bin/python main.py examples/seal.txt .tmp/run-baseline run 100 5

Output: 5 OFF files, one per save block (`100_run_.off` ... `500_run_.off`).

These files are the bytes the simulator produces with the 1-based `indexing.py`
wrappers in place. The smoke test (`tests/test_simulator_smoke.py`) asserts the
simulator continues to produce identical bytes — first to verify ongoing 1-based
work, and later to lock down the 0-based migration.
```

- [ ] **Step 3: Commit**

```bash
gitops.sh add tests/golden/
gitops.sh commit "test: capture golden OFF output from working 1-based simulator"
```

---

### Task 6: Add the regression smoke test

**Files:**
- Create: `tests/test_simulator_smoke.py`
- Create: `tests/__init__.py` (empty)

- [ ] **Step 1: Create the empty package marker**

Run: `touch tests/__init__.py`

- [ ] **Step 2: Write the smoke test**

Write `tests/test_simulator_smoke.py`:

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

- [ ] **Step 3: Run the test**

Run: `.venv/bin/pytest tests/test_simulator_smoke.py -v`

Expected: PASS — the simulator produces byte-identical output to the goldens captured one minute ago.

If it fails (most likely cause: nondeterminism — e.g., dict iteration order, set iteration, hash randomisation), stop and investigate. The migration plan depends on byte stability; nondeterminism is a blocker.

- [ ] **Step 4: Commit**

```bash
gitops.sh add tests/__init__.py tests/test_simulator_smoke.py
gitops.sh commit "test: add end-to-end smoke test against golden output"
```

---

### Task 7: Open the PR

**Files:**
- (none — git/gh workflow only)

- [ ] **Step 1: Confirm branch and create one if needed**

Run: `git branch --show-current`

If the current branch is `main`, create a new branch first:

```bash
git checkout -b fix/finish-file-reading
```

If it's already a feature branch (e.g. `fix/finish-file-reading`), continue.

- [ ] **Step 2: Push and open the PR**

Run:

```bash
gitops.sh push-new
ghops.sh pr create --base main --title "fix: finish file reading and capture golden output" --body "$(cat <<'EOF'
## Summary
- Fixes the off-by-one in `esclec.read_param_file` so the parameter file actually loads.
- Fixes the initialisation order in `main.main` so `num_active_cells` survives the parameter-file clobber-and-restore dance.
- Captures golden OFF output from a working end-to-end run on `examples/seal.txt`.
- Adds a pytest smoke test that compares simulator output to the goldens byte-for-byte.

This unblocks the next milestone (converting the codebase from 1-based to 0-based indexing — see `docs/plans/2026-05-03-zero-based-indexing.md`) by giving us a regression check.

## Test plan
- [ ] `.venv/bin/pytest tests/test_simulator_smoke.py -v` passes
- [ ] `.venv/bin/python main.py examples/seal.txt .tmp/out run 100 5` produces 5 OFF files
EOF
)"
```

- [ ] **Step 3: Report the PR URL**

---

## Self-review

**Spec coverage:** the goal was to get the file-reading + main loop working end-to-end on `examples/seal.txt`, on the existing 1-based wrappers, and capture a golden baseline for the upcoming migration. Tasks 1-3 fix the two bugs identified during investigation; Task 4 absorbs any further failures with explicit triage notes; Tasks 5-6 capture and lock in the baseline; Task 7 lands it.

**Placeholder scan:** Task 4 Step 3 is the loosest part — it lists six likely failure modes (3a-3f) with diagnostic guidance rather than exact code. This is intentional: until the simulator runs further than it does now, the next bug is not fully knowable. The exhaustive bug catalogue would be padding. The Step structure (try → triage → fix → re-run → loop) gives the engineer a clear stopping rule.

**Type consistency:** `Esclec.set_params(core, 1)` keeps the existing 1-based parameter-set index — no change. `core.num_active_cells` is treated as int throughout. `core.Rad` may need an explicit `int()` cast (flagged in Task 4 Step 3c/3d) but that is a fix-as-needed rather than a precondition.

**Known risks the plan flags:**
- Task 4 may surface bugs not yet known. Plan budgets time for triage rather than predicting them.
- Task 6 Step 3 calls out nondeterminism as a possible blocker — if found, the entire migration regression strategy needs rethinking.
