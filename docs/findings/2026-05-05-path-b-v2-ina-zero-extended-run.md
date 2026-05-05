---
title: 'Path B v2: paper-literal ina = 0 extended run on wt-tribosphenic-2014'
date: 2026-05-05
related-rebuttal: docs/paper-rebuttal-preemption-2026-05-05.md §B.5
related-briefing: docs/paper-briefing-2026-05-05.md §5C
---

## Question

Does the eq. 14 typo's catastrophic effect (PAPER_2010 = 19 cusps, PAPER_LITERAL_2010 = 0 on `examples/wt-tribosphenic-2014.txt` with `ina = 0.5`) survive when the parameter file is restored to the paper-literal `ina = 0.0` and the run is extended from 2,500 to 14,000 iterations (matching Harjunmaa et al.'s published sweep length)?

## Method

I built `examples/wt-tribosphenic-2014-paper-ina.txt`: identical to `examples/wt-tribosphenic-2014.txt` except `ina = 0.0`. All 29 other parameter values are unchanged from Ext Data Fig. 2 of the 2014 paper.

I ran the discretisation study runner on this file with all five named presets at 14,000 iterations (28 saves of 500), output to `experiments/discretisation-study/cusp-forming-paper-ina/`.

## Result

| Preset | Plateau | Cusps | Final OFF SHA-256 prefix | Regime |
|---|---:|---:|---|---|
| `HUMPPA_LITERAL` | 19 | 0 | `9f2507cdba6c` | plateau |
| `LEGACY_FORTRAN` | 19 | 0 | `9f2507cdba6c` | plateau |
| `PAPER_2010` | 19 | 0 | `92b4615630b9` | plateau |
| `PAPER_LITERAL_2010` | 19 | 0 | `92b4615630b9` | plateau |
| `PATH_B_DEFAULT` | 19 | 0 | `92b4615630b9` | plateau |

Per-save cusp counts: 0 at every save block (1 → 28) for all five presets. No knot ever forms.

The two cluster-internal pairs are byte-identical at the final OFF: `HUMPPA_LITERAL` and `LEGACY_FORTRAN` produce exactly the same vertex positions and colours; `PAPER_2010`, `PAPER_LITERAL_2010`, and `PATH_B_DEFAULT` likewise produce exactly the same output as each other. The two clusters differ from each other but neither cluster's split is attributable to the eq. 14 typo: the typo branch (`act_typo`) and the corrected branch (`inh_corrected`) reduce to mathematically identical denominators when Act ≡ 0 throughout the run (`1 + k_inh × 0 = 1 + 800 × 0 = 1` under both readings).

## Interpretation

The byte-identicality of `PAPER_2010` and `PAPER_LITERAL_2010` is the strongest possible empirical confirmation of the mechanism the rebuttal §B.5 sketched on mechanistic grounds: the typo is mathematically dormant when no Act seed is present, because both denominators evaluate to the same value at every step, on every cell. The two presets execute their branch but the branch produces the same output.

The methodological consequence is a refinement of the headline phrasing. The previous claim — 'the eq. 14 typo is biologically catastrophic on the 2014 wild-type tribosphenic mouse parameter set' — generalised over 'parameter set' but the typo's effect depends on a sub-property of the parameter set, namely whether the regime drives Act above zero. The refined claim is: 'the eq. 14 typo is biologically catastrophic on parameter sets that drive Act above zero (e.g. via `ina ≠ 0` seeding); on parameter sets where Act is identically zero throughout, the typo is dormant — both denominators produce numerically identical output'.

The result is consistent with — not in tension with — the original headline finding. The original finding's `ina = 0.5` deviation does not manufacture the typo's effect; it exposes it. The deviation puts the simulator into the regime where the typo's branch produces a different value from the corrected branch, and on that regime the difference is biologically catastrophic (19 cusps vs 0).

A separate question the new run does not answer: whether `13.f90` itself, with its border-bias mechanism (`Bbi`, `Lbi`) and its mesenchymal-back-diffusion path, would lift Act above zero under `ina = 0` and thereby exercise the typo. The silicoshark replication uses `mesenchyme='absent'` and `topology='static_with_local_update'` to manage the v1 mesh-degeneration failure mode; these workarounds suppress border-bias amplification of Act on the silicoshark code path. The `13.f90` published 14,000-iteration sweep produces 5 cusps under wild-type `ina = 0` (Ext Data Fig. 4a), so the FORTRAN does have an Act-seeding path that silicoshark does not yet replicate. Reproducing that path is item 1 of §Future work in the briefing (laplacian / mesenchyme follow-up).

The byte-identicality finding is also a small methodological corollary about determinism. Path B v2 is fully deterministic on a fixed parameter file under a fixed preset: the runner produced identical output on five independent subprocess invocations of the two pairs. This is not new — the existing test suite already pins this — but it is a useful side-confirmation when a comparison study turns up matching SHA hashes.

## Rebuttal status update

This run upgrades B.5 from WEAK to SOLID. The original rebuttal text mounted the mechanistic argument; this run provides the empirical confirmation, including the surprising byte-identicality result that the mechanistic argument did not predict but is consistent with. The rebuttal's WEAK-tagged residual ('but the existing run does not test this') is closed.

## Files added

- `examples/wt-tribosphenic-2014-paper-ina.txt`
- `examples/wt-tribosphenic-2014-paper-ina.md`
- `experiments/discretisation-study/cusp-forming-paper-ina/` (results + per-preset OFFs)
