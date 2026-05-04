---
title: 'Differentiation accumulator: FORTRAN uses Sec, paper says Act'
author: Lyndon Drake (with Claude Code)
date: 2026-05-04
---

## Summary

The 2010 paper describes cell differentiation as a function of
*activator* (Act) concentration accumulating over time. The FORTRAN
ancestor (`humppa_translate.f90` line 659, preserved verbatim in
`13.f90` and the existing `coreop2d.py`) accumulates *secretion /
growth-factor* (Sec) instead. The 2014 paper does not specify which
species drives differentiation; it does not contradict the 2010
paper, but it does not endorse the FORTRAN's choice either. This is
a paper-vs-code divergence, not a bug.

## Evidence

**2010 paper (Salazar-Ciudad and Jernvall, *Nature* 464:583).** SI
§'Differentiation' describes `d_i` (cell-i differentiation state) as
incrementing in proportion to local activator concentration when
above a threshold, saturating at `d_i = 1` (the 'knot' state). The
paper's prose names Act consistently throughout this section.

**FORTRAN code.** `humppa_translate.f90:659` (in subroutine
`diferenciacio` / `applydifferentiation`):

```fortran
diff_state(i) = diff_state(i) + delta * q3d(i, 1, 2) * (1d0 - diff_state(i))
```

The third index `2` selects the second species in `q3d`, which is
*secretion* (Sec), not activator (`q3d(:, :, 1)` = Act, `q3d(:, :, 2)`
= Sec, per the cpp_README glossary). This is preserved in `13.f90`
and `coreop2d.py`; both call this slot Sec / `growth_factor` in their
docstrings.

**2014 paper (Harjunmaa et al., *Nature* 512:44).** The 2014 paper
treats the 2010 model as a black box. Its main text and SI do not
specify which species drives `d_i` accumulation. The paper does
discuss SHH-pathway perturbations affecting cusp number, which is
consistent with either interpretation (SHH ≈ Sec in the model
naming, but the model's signalling species don't map one-to-one
onto specific real molecules).

## Why it matters

Path B implements the model from the papers. If Path B follows the
paper literally, differentiation accumulates Act and the simulator
will produce different equilibria from the FORTRAN goldens.
Specifically: Act peaks at knots, so accumulating Act means knots
self-reinforce; Sec is secreted *by* knots and diffuses outward, so
accumulating Sec means non-knot cells in the knot's neighbourhood
also accumulate towards differentiation (a different spatial
pattern).

The FORTRAN's choice (Sec) is therefore not equivalent to the
paper's prose (Act). One of the two has to be wrong, or the paper
is using 'activator' in a looser sense than the model's `q3d(:, :, 1)`.

## Recommended Path B decision

**Follow the FORTRAN.** Reasons:

1. The FORTRAN goldens are the validation oracle. Following the
   paper would make `seal.txt` produce a different equilibrium and
   blow the ±5% tolerance.
2. The FORTRAN's choice has produced biologically plausible tooth
   morphologies in the published 2014 simulations (Ext Data Fig. 2),
   so it's empirically validated even if the paper text is loose.
3. The 2010 paper's Act-driven prose is consistent with the
   FORTRAN's behaviour if 'activator' is read loosely as 'the local
   signalling concentration that triggers knot-cell behaviour' —
   plausible given that Sec is itself produced by knots and
   correlates strongly with Act in space.

Document the choice in the Path B docstring with an explicit
reference to this finding. Do NOT silently follow the paper without
flagging the divergence.

## Alternative: implement both and compare

A cleaner long-term path would be to implement Path B with a
parameterised choice (`differentiation_species: {'act', 'sec'}`),
default to `'sec'` (FORTRAN-compatible), and run a one-off
comparison to see how much the equilibrium changes. If the
difference is small, the 'paper-faithful' interpretation can be
restored as the default. If the difference is large, the FORTRAN's
choice deserves explicit empirical justification in a published
note.

## Provenance of the finding

This was identified by the 2014 paper review subagent during Path B
planning research. See
`docs/research/paper-review-2014-harjunmaa.md` §'What `13.f90`
implements but the paper does not specify' for the original
identification, and `docs/research/paper-review-2010-salazar-
ciudad-jernvall.md` for the 2010 paper's Act-based prose.
