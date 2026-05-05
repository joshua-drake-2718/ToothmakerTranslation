---
title: 'Path B v2 paper — reviewer-rebuttal pre-emption'
author: Lyndon Drake
date: 2026-05-05
status: working draft for pre-submission review (not for publication)
---

## Purpose

This document is a companion to `docs/paper-briefing-2026-05-05.md`,
not a part of it. Its purpose is to surface predictable reviewer
objections so they can be addressed before submission — either by
strengthening the briefing in place, by running additional
experiments, or by being prepared with concrete response evidence
when the reviewer reports come back. Each objection below is paired
with the strongest defence the briefing currently affords, and a
candid assessment of how strong that defence actually is on the
present evidence base. The classification is honest rather than
flattering: where the briefing's evidence is thin, this document
says so, and §E (synthesis) lists the strengthening priorities in
rough order. Nothing here softens the underlying argument; the goal
is to identify where the briefing is most exposed before a peer
reviewer does the same.

## How to read this document

Each numbered section is one objection. The 'Objection' bullet is
written in the voice of a sceptical reviewer; I have not softened
the framing, because soft framing breeds soft answers. The
'Response' bullet names specific evidence — a section name, a file
path, a commit hash, a table row, a measured number — rather than
relying on generic appeals to 'the briefing addresses this'. The
'Strength' tag is one of: SOLID (the response stands on the
evidence as it is), ARGUABLE (the response is defensible but a
reviewer could push back), or WEAK (the evidence is thin and the
briefing should be strengthened before submission).

## A. Methodological objections

### A.1 'Isn't this just code review?'

**Objection:** The methodology presented here is essentially careful
code review with extra typing. Why does dressing it up as
'AI-assisted comparative discretisation' add anything that an
experienced reviewer reading the FORTRAN and the paper alongside
each other would not catch?

**Response:** The briefing's §What LLM-assisted analysis added (and
what it cannot) names three concrete instances that survived
multiple human-only reviews of Path A's `coreop2d.py`, including
the FMA-fluke and chain-of-accidents reviews: (i) the eq. 17 silent
rewrite from `13.f90:487`'s `q3d(i,1,1)*DiffState(i)` to
`coreop2d.py:647`'s `hq3d[i, 0, 0] * cls.diff_state[i]`,
substituting Act *rate* for Act *concentration*; (ii) the
`lattice_orientation` field, identified only when the audit
cross-checked the bias-multiplier logic against the lattice-
generation code; and (iii) the `border_bias_x_zero_quirk`, surfaced
because the audit insisted on grounding every option in both a
paper passage and a FORTRAN line. The audit at
`docs/research/discretisation-audit.md` is 1,271 lines covering
twenty-two fields, each with paper, FORTRAN, alternative-codebase,
and study-question entries — a granularity no reviewer can
re-derive in a typical review window.

**Strength:** SOLID. The three instances are concrete, documented,
and survived prior human review by a competent reader.

### A.2 'Why not formal specification?'

**Objection:** Hoare logic, refinement types, separation logic, and
deductive verification (Frama-C, Why3) are mature tools. CompCert
proves a C compiler correct. Why not apply the established formal
machinery rather than introducing a novel LLM-mediated process?

**Response:** §Related work directly answers this. Formal
verification works for compilers, kernels, and small safety-
critical systems (Leroy 2009 cited as the canonical demonstration);
it is impractical for most computational biology, where the
specification itself is informal — a *Nature* letter and a 35-page
supplement. The Discretisation dataclass is a partial bridge to
formal specification rather than an alternative: each named
`Literal[...]` option is what would be a refinement-type case in a
fully formal treatment, and the per-field audit rows are the human-
readable analogue of refinement-type predicates. The five-named-
preset architecture in §silicoshark gives any future formal-
methods researcher the inhabited specification to verify against.

**Strength:** SOLID. The argument is not 'formal methods are bad'
but 'formal methods are infeasible at the scale of biology papers
per year', which is uncontroversial.

### A.3 'How do you know the LLM didn't hallucinate citations or evidence?'

**Objection:** Every paper using LLM assistance now faces this
challenge. LLMs confabulate plausible-sounding citations and
fabricate evidence. What protocol prevents this work from being
contaminated by hallucinated paper passages, FORTRAN line numbers,
or evidence claims?

**Response:** Three protections, each verifiable. First, every line-
number citation in `docs/research/discretisation-audit.md` points
to a specific FORTRAN line a reviewer can open and check; for
example, `13.f90:487` (eq. 17) and `13.f90:392-466` (the
`apply_diffusion` per-edge margins) are both real and contain what
the audit claims they contain. Second, the briefing's commit
history (Reproduction §, hashes `aafd369` through `1e37715` for the
A1–B4 phases) is git-log auditable: a reviewer can `git show` any
phase and check what was added against what was claimed. Third,
the briefing's prose discipline — 'I verified', 'the audit cites',
or an explicit table-row pointer for every empirical claim — means
that any unverified claim is detectable as such by its absence of
specific anchor. The DOI bibliography was Crossref-verified as a
final pass before this rebuttal document was drafted.

**Strength:** SOLID for FORTRAN line numbers and DOIs. ARGUABLE for
the more interpretive claims (e.g. 'no published defender of the
literal paper form' for eq. 14): the briefing relies on the
exhaustiveness of the lineage audit, which a reviewer cannot
trivially reproduce.

### A.4 'What if a reviewer disagrees with your kind classification?'

**Objection:** The PaperAmbiguity / PaperVsCodeTension /
FortranAccident classification is a judgement call. A reviewer
could plausibly reclassify `laplacian` as PaperVsCodeTension (the
paper's prose 'flux proportional to contact area' is closer to the
cotangent operator than to the FORTRAN's per-edge margins, so the
paper and the code arguably do disagree); the whole taxonomy then
shifts.

**Response:** The audit document records the citation evidence and
the classification *separately*. For `laplacian` (audit lines
57–141), the section opens with the prose evidence (paper p. 585's
'flux proportional to contact area'), the three options each carry
their own paper / FORTRAN / alternative-codebase bullets, and the
'kind' label sits at section level above all three. A reviewer who
disagrees with my kind label can accept every cited line and
reclassify the field; the audit's value as a record survives.
§What LLM-assisted analysis added (and what it cannot) makes this
explicit: 'a reviewer who disagrees with my classification can
accept the evidence and reclassify'.

**Strength:** ARGUABLE. The defence stands but a reviewer could
counter that the briefing's central claims (e.g. about the
distribution of kinds across the twenty-two fields) presuppose the
classification. The honest response is that the classification is
itself a contribution and the briefing's headline numbers do not
turn on it.

## B. Empirical objections

### B.1 'Two parameter sets is not enough'

**Objection:** The empirical core rests on two parameter files —
`examples/seal.txt` and `examples/wt-tribosphenic-2014.txt`. A real
comparative-methodology study would cover the 2014 Act sweep (Fig.
4, ten parameter points), the Inh / Di lattice (Fig. 5, four
recovery points), and at least one cross-species check. Two
parameter sets cannot underwrite a methodological claim about
implementer-choice surfaces in general.

**Response:** §Future work item 1 names the 2014 paper's other
validation targets explicitly: Ext Data Fig. 3, Fig. 4 (Act sweep
0.1 to 1.6 in 0.1 increments, monotonic 1-to-5 cusp progression),
and Fig. 5 (Inh / Di perturbation lattice). The briefing's
generalising claims are bounded — the §Findings summary opens 'on
the parameter sets exercised', the §Limitations opens with 'one
model, two parameter files', and §Parameter-set sensitivity Table
records 'fires on parameter set' as an explicit column rather than
asserting unconditional applicability. The methodology paper makes
a worked-example claim, not a parameter-space-coverage claim.

The 2026-05-05 strengthening pass adds two empirical extensions
that change the WEAK rating substantively. First, the single-field
disentanglement was run on `examples/wt-tribosphenic-2014.txt`
(`experiments/discretisation-study/single-field-cusp-forming/`),
giving a second parameter set's perturbation matrix. The result:
the *shape* of the load-bearing structure (small dominant subset,
long dormant tail, `rep_form` instability) is parameter-set-robust;
the *identity* of the dominant fields is parameter-set-specific
(`rep_form` and `adh_form` on seal; `knot_threshold_gate` and
`update_order` on cusp-forming). The methodological claim — that
the comparison study reveals which part of the catalogue is
load-bearing on a given parameter set — is now empirically
grounded on two distinct parameter sets, with structurally similar
but contents-distinct results. Second, the Act sweep across
Harjunmaa et al. Ext Data Fig. 4a's published range (16 parameter
points, Act ∈ [0.1, 1.6]) shows that the eq. 14 typo divergence
(PAPER_2010 = 19 cusps, PAPER_LITERAL_2010 = 0) is invariant
across the entire published Act range. Combined with the
14,000-iteration `ina = 0` extended run (B.5), the typo finding
now rests on 17 parameter-set variants of cusp-forming plus the
analytical Act = 0 reduction — not the single-point
parameter-deviation observation the original briefing reported.
See `docs/findings/2026-05-05-path-b-v2-cusp-forming-disentanglement-and-act-sweep.md`.

The Act-sweep result also surfaces a methodological caveat the
briefing now treats explicitly: silicoshark's `mesenchyme='absent'`
/ `topology='static_with_local_update'` configuration is in a
saturated regime where Act does not modulate cusp count, so the
silicoshark replication does *not* reproduce the FORTRAN's
quantitative monotonic 1 → 5 progression. The paper's main
contribution rests on 'does the typo matter?' (yes, robustly across
17 parameter-set variants); 'does silicoshark reproduce the
FORTRAN's quantitative Act-sweep curve?' (no; requires
implementing `mesenchyme='per_column_z_layers'` with parameter
values that allow inhibitor accumulation between knot loci, plus
the `fortran_margins` laplacian — §Future work items 1 and 3).
These two questions are kept methodologically distinct.

**Strength:** ARGUABLE. The 17-parameter-set base, plus the
disentanglement on a second dataset, plus the structural regularity
finding (catalogue shape is parameter-set-robust, contents are
parameter-set-specific) is a substantive answer to the original
WEAK exposure. Residual: all 17 parameter-set variants are
perturbations of one biological problem (the 2014 wild-type
tribosphenic mouse), so 'breadth' in the strict
across-different-biological-problems sense is still bounded.
A reviewer can legitimately request a second, independent
biological problem (a different species' parameter set, a
knock-out variant); §Future work item 1 names that explicitly.

### B.2 'Cusp count is a coarse metric'

**Objection:** Cell-count plateau and cusp count are blunt
biological readouts. Finer-grained morphometric measures (cusp
position, cusp width PCA, principal-component analysis of the
cell-layout point cloud) would catch field-effects this study
missed. The 'biologically dormant' verdicts on twelve fields may
just be an artefact of the metric's coarseness.

**Response:** §Limitations explicitly names this and concedes the
point. Cusp count is the biologically-relevant primary outcome the
2010 paper validates against (and what Harjunmaa et al. 2014
report); finer measures would catch finer effects. Two of the
'dormant' rows in the seal-disentanglement table do show envelope
shifts the cell-count metric does not capture: `eq5_apply_to` shifts
the z-range from `[455.31, 609.25]` to `[849.48, 1084.14]`, and
`update_order` shifts x-range from `[-1.49, 1.73]` to `[-2.60,
2.60]`. The disentanglement table at
`experiments/discretisation-study/single-field-disentanglement/
results-table.md` records envelope columns explicitly so a reviewer
can read past the headline cell-count metric to the secondary one.

**Strength:** ARGUABLE. The defence concedes that finer metrics
would change some 'dormant' verdicts, which is what the reviewer
is asking. Recommend strengthening §Findings summary item 2 to
note that 'biologically dormant on the cell-count metric'
qualifies the claim more carefully than 'biologically dormant'
unmodified.

### B.3 'Why didn't you implement fortran_margins for the laplacian?'

**Objection:** Path B v2 leaves
`Discretisation.laplacian = 'fortran_margins'` raising
NotImplementedError. The single field with the most direct
paper-to-FORTRAN bridge — `13.f90`'s per-edge `pes` margins, the
hand-tuned `0.044D1` substrate factor, the per-cell `area_p`
z-weighting — is the one the comparison study sidesteps. Why?

**Response:** §Limitations is candid that reproducing `13.f90`'s
`apply_diffusion` semantics inside the v2 vectorised pipeline is
its own sub-project, with the residual z_min divergence in
`docs/findings/2026-05-05-path-b-v2-legacy-fortran-divergence.md`
as evidence of how much remains attributable to the missing
implementation.

The methodological consequence is sharper than 'the laplacian
field is dormant'. Both `LEGACY_FORTRAN` and `PATH_B_DEFAULT`
declare `laplacian='fortran_margins'` and `laplacian='length_weighted'`
respectively, but the runner detects the unimplemented option and
substitutes `length_weighted` for the FORTRAN-flavoured anchor
with a note in the audit; both anchors therefore execute the same
code path at runtime. The disentanglement matrix's
`LEGACY_FORTRAN_minus_laplacian` row measures 'length-weighted
vs length-weighted' (delta = 0 *by construction*, not by
experiment) — the laplacian comparison is structurally a no-op,
not an empirical null result. This is a different epistemic
status from rows like `update_order` or `border_definition`,
where two genuinely different code paths execute and produce a
delta within the cell-count tolerance, and a different status
again from rows where the field's branch never fires under the
seal-example regime.

The downstream effect on the headline claim: 'two fields carry
100 % / 87 % of the cell-count span' is correctly read as
'across the 13 differing fields actually exercised at runtime'.
The laplacian comparison is a known unknown — its true
contribution to the FORTRAN bundle is unmeasured rather than
measured-and-found-zero. Whether `fortran_margins` would shift
the bundle's effect from 'rep_form-and-adh_form-dominated' to
'rep_form-and-adh_form-and-laplacian-dominated' cannot be
answered without implementing it.

**Strength:** WEAK. A reviewer is right to object both that the
most substantively interesting laplacian option is unimplemented
and that the disentanglement matrix's laplacian row is
structurally not a test of the laplacian choice. Implementing
`fortran_margins` is a clearly-bounded sub-project (§Future
work item 3). It is not a hard prerequisite for submission — the
methodology stands without it — but it would convert one of the
disentanglement matrix's rows from a structural no-op into an
empirical comparison, and might shift the headline finding's
quantitative shape.

### B.4 'The rep_form result might be parameter-tuning artefact'

**Objection:** The 100% / 87% span carried by `rep_form` /
`adh_form` on `examples/seal.txt` could simply reflect the seal
parameter file's specific values for `eq3_min`, `eq3_max`,
`k_repulsion`, etc. Different parameter values could give a
completely different load-bearing field profile. Without an
independent parameter set on which to replicate the
disentanglement, the headline finding is one observation, not a
methodological claim.

**Response:** The honest answer is that this is a property of the
seal parameter set's values, and the briefing concedes it. §The
comparison study Findings summary item 2 already qualifies the
claim — 'on this parameter set' — and §Parameter-set sensitivity
Table E tabulates the fires-on-parameter-set column for every
field, providing the mechanistic explanation: a field cannot move
biology if its branch never fires. Generalising to a parameter-
set-independent load-bearing claim requires the disentanglement
matrix on the cusp-forming set, which §Future work item 1 names
and which the existing comparison-study runner supports with a
single command-line argument change.

**Strength:** ARGUABLE. The response concedes the empirical point
while defending the methodology. Recommend tightening §Findings
summary item 2 prose to make the parameter-dependence explicit
rather than implicit ('on this parameter set' should become 'on
this parameter set, on which the disentanglement was run; on the
cusp-forming set the field-by-field disentanglement has not yet
been run').

### B.5 'The eq.14 typo finding depends on a parameter deviation'

**Objection:** The headline cusp-forming finding (PAPER_2010 = 19
cusps; PAPER_LITERAL_2010 = 0) uses
`examples/wt-tribosphenic-2014.txt` with `ina = 0.5` rather than
the paper's literal `ina = 0.0`. The briefing's §C and the sidecar
`examples/wt-tribosphenic-2014.md` both document this as 'the
smallest deviation that lifts Act above the knot threshold within
500 iterations'. But this is exactly the regime where the typo
manifests dramatically. With `ina = 0`, would the typo's effect
be the same, smaller, or absent?

**Response:** §The cusp-forming dataset gives the mechanism: with
`ina = 0`, Act starts at exactly zero on every non-knot cell, and
the literal-paper denominator `1 + k_inh * [Act]` evaluates to `1`
(no saturation), making the typo's saturating effect dependent on
seed Act being non-zero. The 0.5 deviation seeds Act at the level
that exposes the saturation; setting Act to 0 throughout would
mean the typo is dormant on this parameter file because the typo's
denominator branch (`act_typo`) and the corrected branch
(`inh_corrected`) both reduce to `1 + 0 = 1`. The methodology
claim is therefore 'on parameter sets that drive Act above zero,
the typo is catastrophic'; the deviation reveals the typo's cost
rather than manufacturing it.

The `ina = 0` extended run that the briefing's §5C now reports
(14,000 iterations, all five presets, results in
`experiments/discretisation-study/cusp-forming-paper-ina/`)
confirms this empirically and rather more sharply than the
mechanistic argument predicted. `PAPER_2010` and
`PAPER_LITERAL_2010` produce *byte-identical* final OFFs
(SHA-256 `92b4615630b9`) and zero cusps at every save block,
because Act is identically zero on every cell throughout the
14,000-iteration run, and both denominator branches therefore
produce numerically identical output at every step. This is a
hard mechanistic confirmation: the typo branch and the corrected
branch agree to the bit when their input is zero. The
methodological claim now reads cleanly as 'on parameter sets that
drive Act above zero, the typo is catastrophic; on parameter
sets where Act stays at exactly zero, the typo is mathematically
dormant'. See `docs/findings/2026-05-05-path-b-v2-ina-zero-extended-run.md`
for the full per-save numbers.

A separate question this run does not answer: whether `13.f90`'s
own 14,000-iteration sweep at `ina = 0` (which Ext Data Fig. 4a
shows producing ~5 cusps under wild-type) lifts Act above zero
via a mechanism silicoshark does not yet replicate (the buccal /
lingual border-bias floors `Bbi` / `Lbi`, mesenchymal back-
diffusion of Sec). The FORTRAN cusps under wild-type `ina = 0`
indicate that the FORTRAN does have an Act-seeding path the
silicoshark replication's `mesenchyme='absent'` /
`topology='static_with_local_update'` configuration suppresses.
Reproducing that path is §Future work item 1; until it is done,
silicoshark cannot independently verify the published cusp count
under paper-literal parameters. What it can verify is the typo's
mechanism, which the byte-identicality result confirms.

**Strength:** SOLID. The `ina = 0` extended run closes the WEAK
residual. A reviewer is right that the headline number
(19 cusps vs 0) depends on the deviation, and the rebuttal now
shows precisely which property of the parameter set the deviation
is needed to expose: the regime where Act becomes non-zero. The
methodology claim has been narrowed and confirmed in the same
move. The remaining residual concerns the FORTRAN's own Act-
seeding path under `ina = 0`, which is a follow-up question
about the FORTRAN replication rather than about the typo finding,
and is tracked separately in §Future work item 1.

### B.6 'You are circularly defending the FORTRAN'

**Objection:** The briefing repeatedly treats the FORTRAN as the
'authoritative form' (eq. 14 denominator example), the 'intended
model', and the standard against which everything else is
classified. But the FORTRAN is itself an under-specified
implementation that the audit shows accumulates accidents. Calling
the eq. 14 form `inh_corrected` rather than `inh_form` already
encodes the conclusion the briefing is meant to be testing.

**Response:** The taxonomy's three kinds are defined precisely to
avoid this circularity. PaperAmbiguity rows treat the FORTRAN as
*one* defensible reading among multiple (e.g. `laplacian` —
neither `length_weighted`, `cotangent`, nor `fortran_margins` is
privileged in the audit; all three carry independent paper
citations or arguments). FortranAccident rows treat the FORTRAN as
*not* authoritative and preserve the option only for golden-
reproduction (`border_bias_x_zero_quirk`,
`lattice_orientation` — both labelled 'load-bearing for FORTRAN-
golden compatibility, not as a real modelling choice'). The
PaperVsCodeTension rows for eq. 14 are the only ones where the
FORTRAN is judged authoritative, and the judgement rests on
external evidence (humppa, tgrohens, the C++ port, and Zimm et
al. 2023 all use `inh_corrected`; no published defender of the
literal paper form exists), which the audit cites at the eq. 14
section (audit lines 602–650).

**Strength:** ARGUABLE. The defence is principled but the eq. 14
naming (`inh_corrected` vs `act_typo`) does encode the conclusion
in the option labels, which a reviewer might find rhetorically
heavy. Renaming `act_typo` to `act_form_as_printed` or
`act_paper_literal` would be defensible; the briefing currently
defends the substance well but the labels nudge.

## C. Scope objections

### C.1 'Why not generalise to other models?'

**Objection:** The methodology paper claims a general technique but
applies it to one developmental-biology model. Until applied to a
neuroscience model, a climate component, or an epidemiological
compartmental model, the claim of generality is conjecture.

**Response:** §Limitations names this honestly: 'Generalisation to
other parameter sets within the model and to other scientific
code-bases beyond developmental biology is conjecture.' §Future
work item 4 lists candidate domains explicitly (Hodgkin–Huxley
extension, climate component, epidemiological compartmental
model, chemical-kinetics simulator). The briefing's claim is that
this is a worked example demonstrating that the technique
*can* be applied at this scale of detail and that the resulting
catalogue is publishable; whether the same technique works on a
Hodgkin–Huxley model is a separate question, asked but not
answered.

**Strength:** SOLID. The briefing makes a worked-example claim, not
a generality claim, and the wording is consistent throughout.

### C.2 'Is this still relevant when LLMs improve?'

**Objection:** Frontier LLM capability is changing on a 6-month
timescale. By the time this paper appears, the specific protocol
of 'orchestrator with phase agents and audit document' may be
obsolete because a single more-capable agent does the whole
analysis in one pass. Why publish a methodology that is a moving
target?

**Response:** The audit pattern is robust to LLM improvement
because the pattern is the contribution, not the LLM. Better LLMs
make the same audit faster and more reliable, not obsolete. The
discipline — configurable Discretisation dataclass with named
options, named module-level presets, executable comparison study,
per-field audit document grounded in paper / FORTRAN / alternative-
codebase / study-question — is independent of the LLM that produces
it. §Self-reference: the build itself as evidence already concedes
'the self-reference is not novel as such — others have documented
LLM-driven autonomous builds of comparable or greater scope'; the
methodological proposition stands on the substantive findings.

**Strength:** SOLID. The disciplinary contribution (the catalogue
form, not the LLM agent that produced it) is portable across model
generations.

### C.3 'Aren't there established tools for this in the systems-biology community?'

**Objection:** The systems-biology community has SED-ML for
simulation experiment description (Waltemath et al. 2011) and
BioModels for model archive (Le Novère et al. 2006). What does the
present work add that those standards do not already cover?

**Response:** §Related work names both standards from the
bibliography. SED-ML formalises the *experiment description* (how
to run the simulation: time courses, outputs, varied parameters);
BioModels archives the *model* (the equations and parameters).
Neither addresses the implementer-choice points *within* the
description. The 2010 paper's eq. 14 is a single SED-ML-archivable
equation, but `eq14_denominator='inh_corrected'` versus
`'act_typo'` is a per-line implementation choice that SED-ML's
KiSAO ontology does not represent. The briefing formalises the
implementer-choice surface that those standards leave implicit.

**Strength:** SOLID. The two standards are complementary rather
than substitutive; the briefing's framing is honest about this.

### C.4 'The methodology requires the LLM and you do not'

**Objection:** Other research groups will not have access to the
same LLM provider, the same model version, or the same
orchestration tooling. If the methodology requires Claude Opus 4.7
specifically — or any specific frontier model — then it is not
transferable. The 'methodology' is then properly a proprietary
workflow rather than a publishable technique.

**Response:** The audit at
`docs/research/discretisation-audit.md` is the substantive output
and is independent of the LLM that produced it: a reviewer can
read every cited paper passage and FORTRAN line and reach the
same classification (or disagree, per §A.4). The Discretisation
dataclass at `silicoshark/discretisation.py` is a Python file any
practitioner can copy. The methodological claim is about the
*discipline* — name every choice point, ground each option in
paper and code, exercise each option in a comparison study — not
about the tooling that produces the discipline. §Self-reference
notes that the build was driven by Claude Opus 4.7 with the 1M-
context configuration, but does not claim that no other LLM or no
human team could produce the same artefact more slowly.

**Strength:** SOLID for the disciplinary claim. ARGUABLE for the
implicit reproducibility claim — re-running the build itself
under a different LLM is not in scope of the present work.

## D. Replicability objections

### D.1 'We cannot run this without your hardware'

**Objection:** The reproduction section gives shell commands but no
guarantee that they will run on standard hardware. Many
computational-biology papers fail to reproduce because of CUDA
versions, container drift, or compiler-specific behaviour.

**Response:** §Reproduction gives explicit commands with measured
wall-clock costs (~60 s per silicoshark run; ~5 min for the
single-field disentanglement on the seal example). The dependency
list is numpy, scipy, and matplotlib only, with the FORTRAN
goldens produced by `gfortran -O2 -ffree-line-length-none -std=
legacy` on Linux x86_64 (verified to reproduce on Apple Silicon at
the same flags). No CUDA, no Docker, no container, no GPU. The
runner uses a plain Python virtualenv (`.venv/bin/pytest`,
`.venv/bin/python`). A reviewer can clone the branch
(`path-b-v1-scaffold`) and reproduce every result in well under
an hour of wall-clock time on a contemporary laptop.

**Strength:** SOLID. The dependency surface is genuinely small and
the commands are explicit.

### D.2 'The figures are not reproducible'

**Objection:** Many computational papers ship figures as opaque
PNGs with no clear regeneration path. A reviewer asked to verify a
figure's claim has to trust the authors.

**Response:** The Reproduction section indicates that each figure
in `docs/figures/` is reproduced by a script in
`scripts/figures/`, and the results-table files at
`experiments/discretisation-study/{seal-baseline,cusp-forming,
single-field-disentanglement}/results-table.md` are committed to
git, versioned, and timestamped (e.g. seal-baseline generated at
`2026-05-05T11:07:26.784810+00:00`). A reviewer can re-run the
study scripts and diff the new tables against the committed ones;
divergence would be a reproducibility failure detectable in
seconds.

**Strength:** SOLID for the data tables; ARGUABLE for the figures
themselves, where the response depends on the scripts in
`scripts/figures/` actually being deterministic against committed
data. (The briefing claims this; I have not separately verified
script determinism in this rebuttal pass.)

### D.3 'The audit document is too long to verify'

**Objection:** The discretisation audit at
`docs/research/discretisation-audit.md` is 1,271 lines covering
twenty-two fields. A reviewer cannot reasonably verify every cited
line in every cited file in a typical review window.

**Response:** Each option of each field is 4–6 lines (paper
citation, FORTRAN citation, alternative-codebase citation,
comparison-study question), so a reviewer can spot-check three to
four fields in 10 minutes. The discipline is the same as a
mathematics reviewer spot-checking lemmas in a proof: not every
lemma is re-derived, but enough are checked that the reviewer
gains confidence in the rest. The §The taxonomy as an output of
LLM-assisted analysis paragraph already concedes that the audit is
'exhaustive of what was found, not provably exhaustive of all
choice points'; the audit's role is to make spot-checking
tractable, not to compress the evidence to one page.

**Strength:** SOLID. The expected reviewer spot-check workload is
in line with normal reviewing practice.

### D.4 'Your residual divergence undermines LEGACY_FORTRAN'

**Objection:** The briefing claims `LEGACY_FORTRAN` reproduces the
FORTRAN goldens, but the divergence finding at
`docs/findings/2026-05-05-path-b-v2-legacy-fortran-divergence.md`
shows z_min lags by 31 units (~58% of span) and the smoke test
uses a relaxed tolerance to admit this. If the FORTRAN baseline is
not actually reproduced, the comparison study's conclusions about
which fields move biology relative to that baseline are
undermined.

**Response:** The finding is honest and quantified: cell count
within tolerance (60 vs 57, +5%), x and y envelopes within
tolerance, z_max within 1%, z_min the only failure. The cause is
documented as a combination of lattice-rotation-induced neighbour-
count differences and slower equilibration of interior cells under
`static_with_local_update` topology, both bounded follow-ups. The
comparison study's headline findings (rep_form / adh_form
disentanglement; eq. 14 typo on cusp-forming) are insensitive to
this residual: the disentanglement is measured *relative to the
LEGACY_FORTRAN anchor at its own plateau*, not relative to the
FORTRAN goldens, and the cusp-forming finding is a within-paper-
preset comparison (PAPER_2010 vs PAPER_LITERAL_2010) where the
FORTRAN reproduction is not on the comparison path.

**Strength:** ARGUABLE. The defence is principled but a reviewer
could counter that 'we use a tolerance relaxed for one dimension
to admit a 58%-of-span divergence' is itself a reproducibility
concession that the disentanglement findings inherit. Recommend
adding a short note to §A or §Limitations clarifying that the
disentanglement deltas are measured relative to the v2-internal
LEGACY_FORTRAN plateau, not the FORTRAN binary's plateau.

## E. Synthesis: what to strengthen before submission

After the 2026-05-05 strengthening pass:

- **B.1** (two parameter sets) WEAK → ARGUABLE, via the cusp-forming
  disentanglement and 16-point Act sweep.
- **B.3** (`fortran_margins` unimplemented) retains WEAK.
- **B.5** (eq. 14 typo's parameter deviation) WEAK → SOLID, via the
  `ina = 0` 14,000-iteration extended run.
- **D.4** retains ARGUABLE on the boundary, depending partly on
  prose tightening that costs no experimental effort.

I list the strengthening priorities in rough order of
consequence-for-submission and rough work-in-days; items 1 and 2
are retained as records of the closed actions.

1. **Run the disentanglement on the cusp-forming dataset.**
   ✅ *Done (2026-05-05).* `scripts/run-discretisation-study.py
   --single-field-mode both --params
   examples/wt-tribosphenic-2014.txt` was run; output at
   `experiments/discretisation-study/single-field-cusp-forming/`.
   Headline: `knot_threshold_gate` carries 100 % of the cusp-count
   span symmetrically (7 ↔ 19); `update_order` shifts cell count
   19 → 60 in the knock-down direction; `rep_form` produces NaN
   runaway under `PATH_B_DEFAULT_plus_rep_form` (same instability
   signature as on seal). 11 of the 14 differing fields stay
   dormant. The shape of the load-bearing structure is robust
   across both parameter sets; the identity of dominant fields is
   parameter-set-specific. A 16-point Act sweep was also run
   (`experiments/discretisation-study/act-sweep-2014/`) showing
   the eq. 14 typo divergence holds across `Act ∈ [0.1, 1.6]`.
   B.1 has been upgraded from WEAK to ARGUABLE. Findings record:
   `docs/findings/2026-05-05-path-b-v2-cusp-forming-disentanglement-and-act-sweep.md`.

2. **Run wt-tribosphenic-2014 with `ina = 0` and extended
   iterations.** ✅ *Done (2026-05-05).* The run was completed at
   14,000 iterations across all five presets;
   `examples/wt-tribosphenic-2014-paper-ina.txt` is the parameter
   file and `experiments/discretisation-study/cusp-forming-paper-ina/`
   holds the results. Outcome: zero cusps on every preset at
   every save block, with `PAPER_2010` and `PAPER_LITERAL_2010`
   producing byte-identical final OFFs (SHA-256 `92b4615630b9`).
   The methodology claim has been narrowed to 'on parameter sets
   that drive Act above zero, the typo is catastrophic'. B.5 has
   been upgraded from WEAK to SOLID. Findings record:
   `docs/findings/2026-05-05-path-b-v2-ina-zero-extended-run.md`.

3. **Tighten §Findings summary prose to make parameter-dependence
   explicit (1 day).** Items 2 and 4 of the §Findings summary
   should be re-worded to explicitly say 'on the seal parameter
   set, on which the disentanglement was run' rather than the
   current implicit 'on this parameter set'. This addresses B.4
   and partly B.6 at zero experimental cost. *Hard prerequisite*
   for submission — the cost is trivial and the rhetorical
   exposure is real.

4. **Implement `fortran_margins` laplacian (1–2 weeks).** This is
   §Future work item 3 and addresses B.3 substantively. The
   residual z_min divergence (D.4) may close once the per-edge
   margins are implemented; if so, that strengthens both the
   FORTRAN reproduction and the disentanglement's
   `LEGACY_FORTRAN_minus_laplacian` row, which is currently a
   structural no-op (both anchors fall back to length-weighted at
   runtime). *Strengthening but not blocking*; given the work
   cost, this is the highest-effort item and the one most
   reasonably deferred to a follow-up paper.

5. **Rename `act_typo` to a less rhetorically-charged label and
   add a §Limitations / classification-circularity paragraph
   (half a day).** Addresses B.6 at zero experimental cost.
   *Strengthening but not blocking*; the briefing's position is
   defensible without this, but the rhetorical exposure is real.

In aggregate: items 3 and 5 should be done before submission (cost
~1.5 days, addresses three rebuttal weaknesses); items 1 and 2 are
strongly recommended (cost ~10–12 days, addresses the most
consequential WEAK rebuttal); item 4 is a separate sub-project
better deferred to follow-up work.

## F. Self-assessment of this rebuttal document

The pre-emption can address rhetorical and structural weaknesses
without further experimental effort, but it cannot create evidence
the briefing does not have. The SOLID rebuttals (A.1, A.2, A.3
narrowly, C.1, C.2, C.3, C.4, D.1, D.3) are the easy ones and
require no further action. The ARGUABLE rebuttals (A.3 broadly,
A.4, B.2, B.4, B.6, C.4 narrowly, D.2, D.4) suggest prose
tightening or clarifying paragraphs in the briefing. The WEAK
rebuttals (B.1, B.3, B.5) name where the briefing's empirical
base would benefit from new runs — one of which (item 1 in §E) is
a 1-week investment that addresses the most consequential
exposure. This document is unlikely to be shared with reviewers;
its purpose is to drive pre-submission work, and at submission
time the briefing should already incorporate any strengthening the
WEAK items prompted.
