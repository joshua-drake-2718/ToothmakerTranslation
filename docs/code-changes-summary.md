---
title: "ToothMaker Translation: Code Changes Summary"
author: Lyndon Drake
date: 2026-04-09
---

# ToothMaker Translation: Code Changes Summary

This document walks you through every change we made to your
Python translation. There are 6 pull requests on GitHub. Merge them in
order (1 through 6) and you will have a working simulator. I have
grouped the changes by theme so you can see the pattern behind each
bug, not just the fix.

---

## 1. Project hygiene (PR #1)

A `.gitignore` file tells git which files to leave out of version
control. Your repo was tracking `.DS_Store` (a macOS metadata file)
and `__pycache__/` (compiled Python bytecodes). These change every
time you open the folder or run the code, so they create noise in your
commit history.

The `.gitignore` we added covers the common Python artefacts:

```
__pycache__/
*.pyc
.DS_Store
.tmp/
*.egg-info/
.venv/
```

We also ran `git rm --cached` to untrack the files already committed.
They stay on your disk; git just stops watching them.

---

## 2. Three bugs in `indexing.py` (PR #2)

Your `indexing.py` provides 1-based array wrappers so the Python code
can use the same index numbering as the FORTRAN. Every array access in
the entire simulation goes through these wrappers, so bugs here affect
everything.

### 2.1 `range.__len__` could return a negative number

```python
# Before (line 17)
def __len__(self):
    return (self.stop - self.start + abs(self.step)) // self.step

# After
def __len__(self):
    return max(0, (self.stop - self.start + abs(self.step)) // self.step)
```

If `stop < start` (an empty range), the old formula produced a
negative length. Wrapping in `max(0, ...)` ensures it returns 0
instead.

### 2.2 Generator passed as a single argument

```python
# Before (line 34)
return tuple(transform_index(index for index in indices))

# After
return tuple(transform_index(index) for index in indices)
```

Look carefully at where the closing parenthesis sits. The old code
passed the *entire generator* `index for index in indices` as a single
argument to `transform_index`. The fix moves the closing paren so that
`transform_index` is called on *each* index individually, and the
results are collected into a tuple.

### 2.3 `np.bool` removed in NumPy 1.24

```python
# Before (lines 116, 118)
dtype=np.bool

# After
dtype=np.bool_
```

NumPy removed the plain `np.bool` alias. The underscore version
`np.bool_` is the correct NumPy scalar type.

---

## 3. Four bugs in the reaction-diffusion physics (PR #3)

These are in the `diffusion` method of `coreop2d.py` -- the core
physics of the simulator.

### 3.1 Wrong diffusion coefficients (lines 466, 469)

```python
# Before
cls.q3d[:, :, 2] += cls.delta * cls.Da * hq3d[:, :, 2]   # Inhibitor
cls.q3d[:, :, 3] += cls.delta * cls.Da * hq3d[:, :, 3]   # FGF

# After
cls.q3d[:, :, 2] += cls.delta * cls.Di * hq3d[:, :, 2]   # Inhibitor
cls.q3d[:, :, 3] += cls.delta * cls.Ds * hq3d[:, :, 3]   # FGF
```

The activator uses `Da`, the inhibitor uses `Di`, and the FGF
(secondary signal) uses `Ds`. All three lines had `Da`, which meant
the inhibitor and FGF were diffusing at the activator's rate. These
are different parameters in the FORTRAN (lines 497-499 of `13.f90`),
so the three signals need their own coefficients.

### 3.2 FGF knot condition (line 489)

```python
# Before
elif cls.knots[i] > cls.Set:

# After
elif cls.knots[i] == 1:
```

The FORTRAN (line 508) checks `if(knots(i)==1)` -- a simple binary
test for whether the cell is a knot cell. The Python had a threshold
comparison against `Set`, which is a completely different test. Knot
status is a flag (0 or 1), not a continuous value.

### 3.3 Euler update reads from the wrong z-layer (line 497)

```python
# Before
cls.q3d[:, 1, 1:3] += cls.delta * hq3d[:, 0, 1:3]

# After
cls.q3d[:, 1, 1:3] += cls.delta * hq3d[:, 1, 1:3]
```

All the reaction terms above are computed at z-layer 1
(`hq3d[i, 1, ...]`). But the Euler update was reading from z-layer 0,
which is all zeroes. This effectively disabled all reaction dynamics.
Changing `0` to `1` means the update actually uses the values that
were just computed.

---

## 4. Nine bugs in mechanics and cell division (PR #4)

These are spread across `update_cell_position`, `add_cell`, and
`update_border_cells` in `coreop2d.py`. They fall into four patterns.

### Pattern A: comparison inside `abs()` (lines 836, 1212-1214)

```python
# Before
if abs(cls.positions[i, 2] < cls.Bwi):    # line 836
if abs(ux < 1e-13): ux = 0                 # line 1212

# After
if abs(cls.positions[i, 2]) < cls.Bwi:
if abs(ux) < 1e-13: ux = 0
```

This is a subtle precedence bug. `abs(x < val)` first evaluates
`x < val` (which gives `True` or `False`), then takes the absolute
value of that boolean (which is always 0 or 1). What we actually want
is `abs(x) < val` -- take the absolute value of `x`, then compare.
The parenthesis was in the wrong place in four locations.

### Pattern B: 1D arrays accessed as 2D (lines 909-910, 1275-1276, 1281-1282)

```python
# Before
temp_num_neigh[i] = cls.num_neigh[i, :]    # line 909
cls.diff_state[ii, :] = ...                 # line 1275
cls.knots[ii, :] = ...                      # line 1281

# After
temp_num_neigh[i] = cls.num_neigh[i]
cls.diff_state[ii] = ...
cls.knots[ii] = ...
```

`num_neigh`, `diff_state`, and `knots` are all 1D arrays (one value
per cell). The `, :` slice notation is for 2D arrays. This would raise
a runtime error. Six occurrences fixed across the three arrays.

### Pattern C: comparison instead of assignment (line 1292)

```python
# Before
cls.neigh[i, j] == ii

# After
cls.neigh[i, j] = ii
```

`==` is a comparison that returns `True` or `False` and throws the
result away. `=` is an assignment that actually stores the value. This
meant that after swapping border cells, the neighbour references were
never updated -- the mesh topology would be wrong.

### Pattern D: int treated as array (line 1338)

```python
# Before
num_new_border_cells[:] = 0

# After
num_new_border_cells = 0
```

`num_new_border_cells` is a plain integer, not an array. The `[:]`
slice only works on arrays. Simple assignment is what we need.

### Off-by-one in boundary remap (line 918)

```python
# Before
temp_neigh[i, j] = new_num_active_cells - 1

# After
temp_neigh[i, j] = new_num_active_cells
```

The FORTRAN (line 968) maps out-of-range neighbours to
`new_num_active_cells`, not `new_num_active_cells - 1`. Because the
Python uses 1-based indexing (via your `indexing.py` wrappers), the
FORTRAN value translates directly without adjustment.

---

## 5. Parameter file reading and main loop (PR #5)

This was the feature you were blocked on. The FORTRAN program
(`tresdac`, lines 1858-1862 of `13.f90`) reads 5 command-line
arguments:

```
./runt.e ./P2.txt . P2.txt 9000 2
         ↑        ↑ ↑      ↑    ↑
         input    | output  iters blocks
         file     | name
                  output
                  folder
```

### `read_param_file` (esclec.py)

The parameter file is a plain text file with 30 lines. Each line has a
float value and a parameter name:

```
0.5     Egr
16000   Mgr
1.5     Rep
...
```

The new method reads each line, splits it into value and name, and
stores the value in `cls.parap[i, cls.map]` for indices 3 through 32.
These indices correspond to the parameter slots that `set_params`
already reads from -- so once the file is loaded, `set_params` picks
up the values automatically.

### `initialize_from_parameter_file` (esclec.py)

A small wrapper that sets the map index to 1, calls `read_param_file`,
and if reading succeeded, calls `set_params` to apply the values to
the simulation engine.

### `main.py`

Replaced your variable declarations with a complete program:

1. Parse 5 CLI arguments with `argparse`
2. Initialise the simulation engine (`initial_conditions`,
   `allocate_initial_state`, `initact`)
3. Loop over save blocks, running `iteration(nt)` each time
4. Write `.off` output file after each block

### `examples/seal.txt`

A sample parameter file with 30 lines of default seal tooth
parameters, so you can test the simulator immediately.

---

## 6. OFF mesh output (PR #6)

The `.off` file format (specifically COFF, for Colour OFF) stores a 3D
mesh. You can open these in MeshLab to see the tooth shape. The format
is:

```
COFF
num_vertices num_faces 0
x y z r g b a          ← one line per vertex
x y z r g b a
...
3 v0 v1 v2             ← triangle face (0-based vertex indices)
4 v0 v1 v2 v3          ← quad face
...
```

### `mat(core)`

Computes a material value for each cell. Knot cells get 1.0.
Differentiated cells get 0.1 or 1.0 depending on their
differentiation state relative to the `Int` and `Set` thresholds.
These values drive the colour mapping.

### `get_rainbow` and `get_rainbow_knot`

Map a scalar value to RGBA colour. Knot cells are coloured cyan;
everything else follows a yellow-orange-grey ramp based on the
normalised material value.

### `guardaveinsoff_2`

The existing face-counting logic was already translated but had no
file-writing statements. We added:

- COFF header and vertex count line
- Per-vertex output (position + colour from `get_rainbow_knot`)
- Per-face output (triangle and quad faces with 0-based indices)

Two bugs were also fixed:

- Line 162: `core.neigh[i, k]` used the wrong loop variable `k`
  (from an inner scope) instead of `j` (from the outer loop)
- Line 155: `cls.ma(i)` used function-call parentheses instead of
  array-subscript brackets `cls.ma[i]`

---

## How to merge

The PRs are stacked: each branch is based on the previous one. Merge
them in order on GitHub:

1. PR #1 (`chore/gitignore`) -- merge into `main`
2. PR #2 (`fix/indexing-bugs`) -- retarget to `main`, then merge
3. PR #3 (`fix/diffusion-reaction-bugs`) -- retarget to `main`, then merge
4. PR #4 (`fix/mechanics-division-bugs`) -- retarget to `main`, then merge
5. PR #5 (`feat/file-reading-main-loop`) -- retarget to `main`, then merge
6. PR #6 (`feat/off-output`) -- retarget to `main`, then merge

After each merge, GitHub should let you retarget the next PR to `main`
(it may prompt you automatically).

## How to run

After all PRs are merged:

```bash
python main.py examples/seal.txt output seal 1000 1
```

This reads the seal parameters, runs 1000 iterations, and writes
`output/1000_seal_.off`. Open the `.off` file in MeshLab to see the
tooth shape.
