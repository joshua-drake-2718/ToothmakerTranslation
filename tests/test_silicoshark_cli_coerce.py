"""Tests for the silicoshark CLI's `_coerce` and override pass.

`Discretisation` accepts both Python `None` (e.g. `division_total_cap`)
and the literal string `'none'` (e.g. `knot_threshold_gate='none'`) as
valid field values. The CLI must distinguish these:

- `--override knot_threshold_gate=none`     → string `'none'`
- `--override division_total_cap=null`      → Python `None`

The first form previously coerced to Python `None` (the CLI's earlier
sentinel was `'none'`), which broke when `_apply_overrides` then asked
the simulator to use `knot_threshold_gate=None` — `step_reaction_diffusion`
raises `ValueError: unknown knot_threshold_gate: None`. The fix moved
the None sentinel to `'null'` (JSON convention) and these tests pin
the new behaviour.
"""
from __future__ import annotations

from silicoshark.__main__ import _coerce


def test_coerce_null_returns_python_none():
    assert _coerce('null') is None


def test_coerce_uppercase_null_returns_python_none():
    assert _coerce('NULL') is None


def test_coerce_none_returns_string_none():
    """The literal string 'none' must NOT collapse to Python None.

    `Discretisation.knot_threshold_gate` accepts the string `'none'`
    as one of its two valid values; preserving that round-trip is the
    point of the fix.
    """
    assert _coerce('none') == 'none'
    assert _coerce('none') is not None


def test_coerce_true_false_returns_bool():
    assert _coerce('true') is True
    assert _coerce('false') is False


def test_coerce_int_string_returns_int():
    assert _coerce('60') == 60
    assert isinstance(_coerce('60'), int)


def test_coerce_float_string_returns_float():
    assert _coerce('1.5') == 1.5
    assert isinstance(_coerce('1.5'), float)


def test_coerce_general_string_passes_through():
    assert _coerce('hookean_signed') == 'hookean_signed'
    assert _coerce('first_border_cell') == 'first_border_cell'


def test_coerce_zero_one_are_int_not_bool():
    """Documented edge case: '1' and '0' coerce to int, not bool."""
    assert _coerce('1') == 1
    assert _coerce('1') is not True
    assert _coerce('0') == 0
    assert _coerce('0') is not False
