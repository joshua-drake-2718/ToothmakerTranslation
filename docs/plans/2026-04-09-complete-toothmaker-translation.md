# Complete ToothMaker FORTRAN-to-Python Translation

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete Joshua's translation of the ToothMaker tooth morphogenesis simulator from FORTRAN 90 to Python, fixing all identified bugs and implementing the missing file I/O and main loop.

**Architecture:** The Python translation mirrors the FORTRAN structure: `coreop2d.py` (simulation engine), `esclec.py` (file I/O and parameter management), `main.py` (entry point), `indexing.py` (1-based array wrappers), `vector.py` (vector maths). The simulation models tooth crown development via reaction-diffusion signalling coupled with mechanical forces on a hexagonal cell mesh.

**Tech Stack:** Python 3.11+, NumPy, argparse (stdlib)

---

## Context: The File Reading Issue (Joshua's Blocker)

The FORTRAN program `tresdac` (lines 1858-1862 of `13.f90`) reads 5 command-line arguments:

```fortran
call getarg(1,cac)       ! input parameter file path
call getarg(2,caufolder) ! output folder
call getarg(3,cau)       ! output file name base
call getarg(4,cad)       ! iterations per save block (string -> int)
call getarg(5,cass)      ! number of save blocks (string -> int)
```

The **parameter file** (`cac`) contains 30 lines (indices 3-32), each with a float value and a parameter name:

```
0.013   Egr
16000   Mgr
1.5     Rep
0.0     unused6
0.5     Adh
...
```

The output file name is constructed as:
```
nff = "{caufolder}/{iterations_total * block_number}_{cau}"
nfoff = nff + "_.off"    (3D mesh output)
nfes  = nff + "_.txt"    (parameter snapshot)
```

Example invocation: `./runt.e ./P2.txt . P2.txt 9000 2` runs 2 blocks of 9000 iterations, saving `9000_P2.txt_.off` and `18000_P2.txt_.off`.

None of this is implemented in the Python version. `main.py` has only variable declarations; `esclec.py` has empty stubs for `read_param_file()` and `initialize_from_parameter_file()`.

---

## File Structure

| File | Responsibility | Status |
|------|---------------|--------|
| `main.py` | Entry point, CLI args, simulation loop, output | **Mostly missing** |
| `coreop2d.py` | Simulation engine (all physics subroutines) | **Translated but has bugs** |
| `esclec.py` | Parameter file I/O, mesh output (OFF format) | **Partially translated, 5 stubs** |
| `indexing.py` | 1-based array wrappers for NumPy | **Has 3 bugs** |
| `vector.py` | Vector magnitude, distance, cross product | **Complete** |
| `.gitignore` | Git ignore rules | **Missing** |

---

## Complete Bug Inventory

### indexing.py

| Line | Bug | Fix |
|------|-----|-----|
| 17 | `__len__` can return negative for empty ranges | Wrap in `max(0, ...)` |
| 34 | `tuple(transform_index(index for index in indices))` passes generator as single arg | `tuple(transform_index(index) for index in indices)` |
| 116 | `np.bool` deprecated in NumPy 1.24+ | Use `np.bool_` |

### coreop2d.py

| Line | Bug | Fix |
|------|-----|-----|
| 466 | Inhibitor diffusion uses `cls.Da` | Use `cls.Di` |
| 469 | FGF diffusion uses `cls.Da` | Use `cls.Ds` |
| 497 | Reaction update reads `hq3d[:, 0, 1:3]` (z-layer 0) | Use `hq3d[:, 1, 1:3]` (z-layer 1) |
| 489 | FGF knot condition: `cls.knots[i] > cls.Set` | Should be `cls.knots[i] == 1` (FORTRAN line 508) |
| 836 | `abs(cls.positions[i, 2] < cls.Bwi)` -- comparison inside abs() | `abs(cls.positions[i, 2]) < cls.Bwi` |
| 909 | `cls.num_neigh[i, :]` -- num_neigh is 1D | `cls.num_neigh[i]` |
| 910 | `cls.diff_state[i, :]` -- diff_state is 1D | `cls.diff_state[i]` |
| 918 | `new_num_active_cells - 1` in boundary remap | `new_num_active_cells` (FORTRAN line 968) |
| 1212-1214 | `abs(ux < 1e-13)` -- comparison inside abs() | `abs(ux) < 1e-13` (same for uy, uz) |
| 1275-1276 | `cls.diff_state[ii, :]` -- diff_state is 1D | `cls.diff_state[ii]` etc. |
| 1281-1282 | `cls.knots[ii, :]` -- knots is 1D | `cls.knots[ii]` etc. |
| 1292 | `cls.neigh[i, j] == ii` -- comparison not assignment | `cls.neigh[i, j] = ii` |
| 1338 | `num_new_border_cells[:] = 0` -- int not array | `num_new_border_cells = 0` |

### esclec.py

| Line | Bug | Fix |
|------|-----|-----|
| 143 | `core.neigh[i, k]` uses wrong variable `k` from inner scope | Should be `core.neigh[i, j]` |

---

## PR Plan

### Task 1: Project hygiene -- .gitignore and remove tracked artefacts

**Files:**
- Create: `.gitignore`

This PR adds a `.gitignore` and removes `.DS_Store` and `__pycache__/` from git tracking.

- [ ] **Step 1: Create .gitignore**

```
__pycache__/
*.pyc
.DS_Store
.tmp/
*.egg-info/
.venv/
```

- [ ] **Step 2: Remove tracked artefacts from git index**

```bash
git rm --cached .DS_Store
git rm -r --cached __pycache__/
```

- [ ] **Step 3: Commit and create PR**

```bash
git checkout -b chore/gitignore
git add .gitignore
git commit -m "chore: add .gitignore, remove tracked artefacts"
```

---

### Task 2: Fix indexing.py bugs

**Files:**
- Modify: `indexing.py:17` (`__len__`), `indexing.py:34` (`transform_indices`), `indexing.py:116` (`np.bool`)

These are foundational bugs that affect every array operation in the simulation.

- [ ] **Step 1: Fix `range.__len__` negative return (line 17)**

```python
def __len__(self):
    return max(0, (self.stop - self.start + abs(self.step)) // self.step)
```

- [ ] **Step 2: Fix `transform_indices` generator bug (line 34)**

The current code passes a generator expression as a single argument to `transform_index` instead of mapping over the tuple elements.

```python
def transform_indices(indices):
    if isinstance(indices, tuple):
        return tuple(transform_index(index) for index in indices)
    else:
        return transform_index(indices)
```

- [ ] **Step 3: Fix `np.bool` deprecation (line 116)**

Replace `np.bool` with `np.bool_` in both occurrences in the `bool_array` class (lines 116 and 118).

```python
class bool_array():
    def __init__(self, shape=None, buffer=None):
        if not shape is None:
            self.inner = np.ndarray(shape, buffer=buffer, dtype=np.bool_)
        elif not buffer is None:
            self.inner = np.array(buffer, dtype=np.bool_)
```

- [ ] **Step 4: Commit and create PR**

```bash
git checkout -b fix/indexing-bugs
git add indexing.py
git commit -m "fix: three bugs in indexing.py (range len, transform_indices, np.bool)"
```

---

### Task 3: Fix bugs in coreop2d.py -- diffusion and reaction

**Files:**
- Modify: `coreop2d.py:466,469,489,497`

These bugs affect the core physics. The diffusion coefficient bug means inhibitor and FGF diffuse at the wrong rate; the reaction update reads from the wrong z-layer.

- [ ] **Step 1: Fix diffusion coefficients (lines 466, 469)**

```python
# Inhibitor diffusion
cls.q3d[:, :, 2] += cls.delta * cls.Di * hq3d[:, :, 2]

# FGF diffusion
cls.q3d[:, :, 3] += cls.delta * cls.Ds * hq3d[:, :, 3]
```

- [ ] **Step 2: Fix FGF knot reaction condition (line 489)**

The FORTRAN (line 508) checks `if(knots(i)==1)` but the Python has `elif cls.knots[i] > cls.Set`:

```python
elif cls.knots[i] == 1:
    a = cls.Sec - cls.Deg * cls.q3d[i, 1, 3]
    if a < 0: a = 0
    hq3d[i, 1, 3] = a
```

- [ ] **Step 3: Fix reaction Euler update z-layer (line 497)**

The update reads from `hq3d[:, 0, ...]` (z-layer 0) but all reaction terms were computed at z-layer 1:

```python
cls.q3d[:, 1, 1:3] += cls.delta * hq3d[:, 1, 1:3]
```

- [ ] **Step 4: Commit and create PR**

```bash
git checkout -b fix/diffusion-reaction-bugs
git add coreop2d.py
git commit -m "fix: diffusion coefficients (Di/Ds), knot condition, Euler z-layer"
```

---

### Task 4: Fix bugs in coreop2d.py -- mechanics and cell division

**Files:**
- Modify: `coreop2d.py:836,909,910,918,1212-1214,1275-1276,1281-1282,1292,1338`

These are bugs in `update_cell_position`, `add_cell`, and `update_border_cells`.

- [ ] **Step 1: Fix abs() parenthesis in update_cell_position (line 836)**

```python
if abs(cls.positions[i, 2]) < cls.Bwi:
```

- [ ] **Step 2: Fix 1D array access in add_cell (lines 909-910)**

`num_neigh` and `diff_state` are 1D arrays, not 2D:

```python
temp_num_neigh[i] =     cls.num_neigh[i]
c_diff_state[i] =       cls.diff_state[i]
```

- [ ] **Step 3: Fix off-by-one in boundary remap (line 918)**

The FORTRAN (line 968) maps out-of-range neighbours to `new_num_active_cells`, not `new_num_active_cells - 1`:

```python
if temp_neigh[i, j] > cls.num_active_cells:
    temp_neigh[i, j] = new_num_active_cells
```

- [ ] **Step 4: Fix abs() parenthesis in distance calculation (lines 1212-1214)**

```python
if abs(ux) < 1e-13: ux = 0
if abs(uy) < 1e-13: uy = 0
if abs(uz) < 1e-13: uz = 0
```

- [ ] **Step 5: Fix 1D array access in border cell swap (lines 1275-1282)**

```python
cls.diff_state[ii] = sc_diff_state[cls.first_border_cell]
cls.diff_state[cls.first_border_cell] = sc_diff_state[ii]
# ...
cls.knots[ii] = sc_knots[cls.first_border_cell]
cls.knots[cls.first_border_cell] = sc_knots[ii]
```

- [ ] **Step 6: Fix comparison-as-assignment bug (line 1292)**

```python
cls.neigh[i, j] = ii  # was == (comparison, not assignment)
```

- [ ] **Step 7: Fix int-as-array bug in update_border_cells (line 1338)**

```python
num_new_border_cells = 0  # was num_new_border_cells[:] = 0
```

- [ ] **Step 8: Commit and create PR**

```bash
git checkout -b fix/mechanics-division-bugs
git add coreop2d.py
git commit -m "fix: parenthesis errors, 1D array access, assignment bug in cell division"
```

---

### Task 5: Implement file reading and main simulation loop (Joshua's blocker)

**Files:**
- Modify: `main.py` (complete rewrite), `esclec.py:23-25,62-64`

This is the core missing functionality. It implements:
1. Command-line argument parsing (replacing FORTRAN `getarg`)
2. Parameter file reading (replacing FORTRAN `read_param_file`)
3. Parameter initialisation (replacing FORTRAN `initialize_from_parameter_file`)
4. The main simulation loop with output file generation
5. A sample parameter file for testing

- [ ] **Step 1: Implement `read_param_file` in esclec.py (line 23-25)**

The FORTRAN reads 30 lines from the opened file, each containing a float value and a parameter name string:

```python
@classmethod
def read_param_file(cls, filepath):
    cls.read_failed = 0
    try:
        with open(filepath, 'r') as f:
            for i in range(3, 33):  # indices 3 to 32 inclusive
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

- [ ] **Step 2: Implement `initialize_from_parameter_file` in esclec.py (line 62-64)**

```python
@classmethod
def initialize_from_parameter_file(cls, core, filepath):
    cls.map = 1
    cls.read_param_file(filepath)
    if cls.read_failed == 0:
        cls.set_params(core, cls.map)
```

- [ ] **Step 3: Rewrite main.py with argparse and simulation loop**

```python
import sys
import argparse
import os
from indexing import *
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
    io.set_params(core, 1)
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

- [ ] **Step 4: Create sample parameter file**

Create `examples/seal.txt` with 30 lines (indices 3-32) based on the default values from `initialconditions` in the FORTRAN and the original ToothMaker examples:

```
0.5     Egr
16000   Mgr
1.5     Rep
0.0     unused6
0.5     Adh
0.15    Act
1.0     Inh
0.0     unused10
3.0     Sec
0.0     unused12
0.008   Da
0.2     Di
0.1     Ds
0.0     unused16
0.5     Int
0.2     Set
0.005   Boy
0.03    Dff
10.0    Bgr
0.0     Abi
0.0     Pbi
0.0     Bbi
0.0     Lbi
4.0     Rad
0.2     Deg
100.0   Dgr
0.1     Ntr
3.0     Bwi
0.5     ina
0.0     umgr
```

- [ ] **Step 5: Commit and create PR**

```bash
git checkout -b feat/file-reading-main-loop
git add main.py esclec.py examples/seal.txt
git commit -m "feat: implement parameter file reading and main simulation loop

Resolves the file naming/reading issue. The parameter file is a plain
text file with 30 lines, each containing a float value and parameter
name. The program takes 5 CLI arguments: input file, output folder,
output name, iterations per block, and number of blocks."
```

---

### Task 6: Implement output file writing (OFF format mesh)

**Files:**
- Modify: `esclec.py:67-170`

This completes the mesh output by implementing the stub methods and adding the file-writing statements to `guardaveinsoff_2`.

- [ ] **Step 1: Implement `mat()` (line 161-162)**

Computes per-cell material values for colouring based on knot/differentiation state:

```python
@classmethod
def mat(cls, core):
    cls.ma[:] = 0
    for i in range(core.num_active_cells):
        if core.knots[i] == 1:
            cls.ma[i] = 1.0
        else:
            if core.diff_state[i] > core.Int:
                cls.ma[i] = 0.1
            if core.diff_state[i] > core.Set:
                cls.ma[i] = 1.0
    cls.va_max = max(cls.ma[i] for i in range(core.num_active_cells))
    cls.va_min = min(cls.ma[i] for i in range(core.num_active_cells))
```

- [ ] **Step 2: Implement `get_rainbow()` (line 165-166)**

Maps a value to RGBA colour based on normalised position between min and max:

```python
@classmethod
def get_rainbow(cls, val, min_val, max_val, c):
    if max_val > min_val:
        f = (val - min_val) / (max_val - min_val)
    else:
        f = 0.5

    if f < 0.07:
        c[1], c[2], c[3], c[4] = 0.6, 0.6, 0.6, 0.8
    elif f < 0.2:
        c[1], c[2], c[3], c[4] = 1.0, f, 0.0, 0.5
    elif f < 1.0:
        c[1], c[2], c[3], c[4] = 1.0, f * 3, 0.0, 1.0
    else:
        c[1], c[2], c[3], c[4] = 1.0, 1.0, 0.0, 1.0
```

- [ ] **Step 3: Implement `get_rainbow_knot()` (line 169-170)**

Same as `get_rainbow` but knot cells are coloured cyan:

```python
@classmethod
def get_rainbow_knot(cls, core, val, min_val, max_val, c, i):
    if max_val > min_val:
        f = (val - min_val) / (max_val - min_val)
    else:
        f = 0.5

    if core.knots[i] == 1:
        c[1], c[2], c[3], c[4] = 0, 1, 1, 0
    elif f < 0.07:
        c[1], c[2], c[3], c[4] = 0.6, 0.6, 0.6, 0.8
    elif f < 0.2:
        c[1], c[2], c[3], c[4] = 1.0, f, 0.0, 0.5
    elif f < 1.0:
        c[1], c[2], c[3], c[4] = 1.0, f * 3, 0.0, 1.0
    else:
        c[1], c[2], c[3], c[4] = 1.0, 1.0, 0.0, 1.0
```

- [ ] **Step 4: Complete `guardaveinsoff_2()` with file-writing statements**

The existing face-counting logic needs the `k` -> `j` variable fix on line 143, and the method needs to accept a file handle and write the COFF header, vertices, and faces. The FORTRAN writes:

1. `"COFF"` header
2. `num_active_cells nfaces 0`
3. Per vertex: `x y z r g b a`
4. Per triangle face: `3 v0 v1 v2` (0-based indices)
5. Per quad face: `4 v0 v1 v2 v3`

Update the method signature to accept a file handle parameter, fix the variable bug, and add write statements matching the FORTRAN output format (lines 1658-1713 of `13.f90`).

- [ ] **Step 5: Commit and create PR**

```bash
git checkout -b feat/off-output
git add esclec.py
git commit -m "feat: implement OFF mesh output (mat, get_rainbow, guardaveinsoff_2)"
```

---

## PR Dependency Order

```
PR 1 (gitignore) -----> independent, merge first
PR 2 (indexing) -------> independent, merge early
PR 3 (diffusion bugs) -> depends on PR 2 (correct indexing needed)
PR 4 (mechanics bugs) -> depends on PR 2
PR 5 (file reading) ---> depends on PR 2 (correct indexing for param arrays)
PR 6 (OFF output) -----> depends on PR 4 and PR 5
```

PRs 1 and 2 can be created and merged independently. PRs 3, 4, and 5 depend on PR 2. PR 6 depends on 4 and 5.

---

## Verification Strategy

After all PRs are merged, the simulation should be runnable as:

```bash
python main.py examples/seal.txt output seal 1000 1
```

This should:
1. Read the parameter file without errors
2. Run 1000 iterations of the simulation
3. Write `output/1000_seal_.off` containing a valid COFF mesh

The `.off` file can be viewed in MeshLab or any COFF-compatible viewer to visually confirm tooth-like morphology.
