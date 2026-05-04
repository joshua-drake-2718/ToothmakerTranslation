---
title: 'Research notes — index'
author: Lyndon Drake (with Claude Code)
date: 2026-05-04
---

This directory holds background research compiled while planning Path B
(a paper-faithful re-implementation of the tooth-morphogenesis model).
Each document was produced by a focused background agent and reviewed
manually; cite them where their findings inform code or other docs.

## Contents

### Code provenance and forks

| File | Subject | Use when |
|------|---------|----------|
| [`path-b-research.md`](path-b-research.md) | The Zimm et al. 2023 PNAS paper, the wider Salazar-Ciudad and Jernvall paper lineage, and a survey of related GitHub repositories (jernvall-lab/ToothMaker, RolandZimm/silicoshark, tgrohens/toothmaker, jupander/BITES, others). Includes the subroutine and variable rename tables and the GUI parameter-label glossary. | Understanding which paper to cite; finding a particular FORTRAN function in the upstream sources; mapping cryptic FORTRAN names to paper symbols. |
| [`cpp-port-review.md`](cpp-port-review.md) | Code review of jernvall-lab's C++ port of `humppa_translate.f90` (Feb 2026, described as 'translated with Claude Opus 4.5'). Verdict: *use partially* — its README glossary is reusable, but its file split, `if (suma > 0)` divide-by-zero handling, partial iiii-preservation regression, and per-cell loop shape should not be adopted. | Sanity-checking how the C++ port handles a particular FORTRAN feature; deciding whether to adopt a structural choice it makes. |
| [`tgrohens-review.md`](tgrohens-review.md) | Code review of `tgrohens/toothmaker`, a 2025 minimal FORTRAN fork. Verdict: *useful negative cross-check* — confirms upstream FORTRAN matches humppa; offers a clean three-file split easier to diff, plus one extra example input (`nosea_10261__0.txt`) usable as a second-regime golden. | Confirming whether a particular FORTRAN behaviour is canonical or `13.f90`-specific; finding non-seal example inputs. |
| [`13f90-vs-humppa-divergences.md`](13f90-vs-humppa-divergences.md) | Catalogue of every identified divergence between `13.f90` (English-renamed) and `humppa_translate.f90` (Catalan original), classified as cosmetic, structural, or semantic. The txungu typo is the only confirmed semantic divergence at time of writing; this document is the place to record any further ones. | Verifying that a `13.f90` behaviour is also in humppa before treating it as canonical; ensuring Path B does not re-create `13.f90`-only accidents. |

### Paper reviews

| File | Subject | Use when |
|------|---------|----------|
| [`paper-review-2010-salazar-ciudad-jernvall.md`](paper-review-2010-salazar-ciudad-jernvall.md) | Review of the 2010 *Nature* canonical paper and its SI. Extracts all 18 numbered equations, the SI Table 1 parameter table, the geometric primitives, the discretisation choices, and a 15-item Path-B-implications list distinguishing 'paper specifies' from 'implementer's choice'. Identifies one paper typo (Eq. 14) and several paper-vs-code gaps. | Implementing or cross-checking any specific equation; deciding what's bound by the paper vs. left to Path B's discretion. |
| [`paper-review-2014-harjunmaa.md`](paper-review-2014-harjunmaa.md) | Review of the 2014 *Nature* tribosphenic-extension paper and its (mostly experimental) SI. The paper introduces *no new equations*; the value is the wild-type tribosphenic parameter set in Ext Data Fig. 2, the Act / Inh / Di sweeps as qualitative validation goldens, and an 8-item list of `13.f90`-implements-but-paper-doesn't-specify decision points. Identifies the differentiation-accumulator paper-vs-code finding. | Validating Path B against more than one parameter regime; understanding what's tribosphenic-specific vs. inherited from the 2010 model. |
| [`paper-review-2023-zimm.md`](paper-review-2023-zimm.md) | Review of the 2023 *PNAS* shark paper and its SI. The paper is *out of scope* for Path B's primary implementation but useful for parameter-naming conventions and for distinguishing what's core to the 2010 model vs. shark-specific extensions. Verdict: small transferable findings (Table S1 layout, sign convention, piecewise-threshold typesetting); five substantive 2023-only changes to *avoid*; `multicusp.ini` is not a clean 2010 baseline. | Confirming what's NOT in the 2010 model; borrowing parameter-listing conventions; sanity-checking that `multicusp.ini` is shark-flavoured. |

## Related documents elsewhere

- `docs/findings/2026-05-04-path-a-fma-and-over-division.md` — Path A's
  investigation findings. Contains the 'load-bearing accidents' analysis
  and the verification-method record. Updated after this research with
  inline corrections; see its editor's note for the diff.
- `docs/findings/2026-05-04-differentiation-uses-sec-not-act.md` —
  paper-vs-code finding identified during the 2014 paper review:
  `13.f90` accumulates Sec for differentiation, the 2010 paper's
  prose says Act. Path B has to choose; recommendation is to follow
  the FORTRAN with the choice documented.
- `docs/plans/2026-05-04-path-b-charter.md` — Path B charter, drawing
  on all seven research documents above to define scope, references,
  validation strategy, and what NOT to reproduce.
- The chain-of-accidents docstring at the top of `coreop2d.py` —
  references this directory for the cross-fork provenance evidence.
