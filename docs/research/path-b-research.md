---
title: 'Path B research: Zimm 2023 PNAS paper and ToothMaker forks'
author: Background research subagent
date: 2026-05-04
---

## Paper

### Citation

Zimm R, Berio F, Debiais-Thibaud M, Goudemand N. 'A shark-inspired general model of tooth morphogenesis unveils developmental asymmetries in phenotype transitions.' *Proceedings of the National Academy of Sciences USA* 120 (15) e2216959120 (April 2023). DOI [10.1073/pnas.2216959120](https://doi.org/10.1073/pnas.2216959120). PMID 37027430. PMC PMC10104537. License CC BY-NC-ND 4.0.

The article is 10 pages, with 5 figures.

### Files saved

- `/home/lyndon/repo/project/ToothmakerTranslation/.tmp/papers/zimm-2023-pnas.pdf` — main article, 13 MB, 10 pages, complete with figures. Sourced from Europe PMC at `https://europepmc.org/articles/PMC10104537?pdf=render`.

### Supplementary information — could not fetch

I was unable to download the SI PDF (`pnas.2216959120.sapp.pdf`). PNAS, PMC and HAL all sit behind anti-bot challenges (Cloudflare or Anubis) that block scripted fetches. The known canonical URLs are:

- `https://www.pnas.org/doi/suppl/10.1073/pnas.2216959120/suppl_file/pnas.2216959120.sapp.pdf` (PNAS, returns Cloudflare 403)
- `https://europepmc.org/articles/instance/10104537/bin/pnas.2216959120.sapp.pdf` (Europe PMC, redirects to backend that closes the stream)
- `https://www.ncbi.nlm.nih.gov/pmc/articles/instance/10104537/bin/pnas.2216959120.sapp.pdf` (NCBI legacy, returns generic redirect page)

The SI is not embedded in the main 10-page PDF. The user can grab the SI manually from any of these URLs in a browser session.

A French preprint deposit exists at HAL (`hal-04894916`) under <https://nantes-universite.hal.science/IGFL/hal-04894916v1>, also Anubis-protected.

### What the paper covers (abstract; relevance to Path B)

The paper extends the Salazar-Ciudad and Jernvall (2010, Nature) tooth-development reaction–diffusion–mechanics model to shark teeth, replacing the ToothMaker signalling kernel with a three-component network for Fgf, Bmp and Shh. The companion code is the FORTRAN file `toothmodel.f90` in `RolandZimm/silicoshark` (see Forks section), which is a sibling of our `13.f90`, sharing the same `coreop2d` / `esclec` module structure and 90 % of the subroutine bodies.

This means the 2023 paper is **not** the strict reference for `13.f90` itself — `13.f90` is the older 2014 mammalian tribosphenic version, and the Zimm 2023 work is the shark variant that branched off the same code in 2019. For Path B (paper-faithful re-implementation), the 2010 Nature paper by Salazar-Ciudad and Jernvall is the more authoritative equation reference, with the 2023 SI useful for the modern parameter table conventions.

### Related primary references worth tracking down

- Salazar-Ciudad I, Jernvall J. 'A computational model of teeth and the developmental origins of morphological variation.' *Nature* 464, 583–586 (2010). DOI [10.1038/nature08838](https://doi.org/10.1038/nature08838). The ancestral model.
- Harjunmaa E, Seidel K, Häkkinen T, et al. 'Replaying evolutionary transitions from the dental fossil record.' *Nature* 512, 44–48 (2014). DOI [10.1038/nature13613](https://doi.org/10.1038/nature13613). The 2014 tribosphenic extension that produced `humppa_translate.f90` / `13.f90`.

## Forks

A code search across GitHub found seven repositories descending from the Salazar-Ciudad and Jernvall FORTRAN model. They fall into three groups: the canonical line, vectorisation/translation efforts, and downstream user repositories.

### Canonical line

#### jernvall-lab/ToothMaker

<https://github.com/jernvall-lab/ToothMaker>. Stars 3, forks 0, last push 28 Feb 2026, MIT, primary language C++ (with FORTRAN sources). The GUI repository for the Salazar-Ciudad and Jernvall tooth model. Contains:

- `models/triconodont/src/humppa.f90` — the 2010 triconodont FORTRAN model (2,361 lines).
- `models/tribosphenic/src/fortran/humppa_translate.f90` — **the direct ancestor of our `13.f90`** (2,401 lines). Module names `coreop2d` and `esclec` match exactly; subroutine bodies are essentially identical with Catalan/Spanish names instead of English.
- `models/tribosphenic/src/cpp/` — a complete C++ port of `humppa_translate.f90` published February 2026. The README states it was 'Translated with Claude Opus 4.5'. Files: `tooth_model.hpp` (278 lines), `tooth_model_core.cpp` (361), `tooth_model_geometry.cpp` (236), `tooth_model_diffusion.cpp` (296), `tooth_model_mechanics.cpp` (475), `tooth_model_division.cpp` (843), `tooth_model_iteration.cpp` (36), `file_io.cpp` (663), `main.cpp` (104). The README contains a complete variable-name translation table — `malla` → `cellPositions`, `vei` → `neighbors`, `marge` → `cellMargins`, `q3d` → `quantities3D`, etc., plus a parameter table mapping FORTRAN names to GUI labels (`tacre` → `epithelialGrowthRate` / Egr, `acac` → `activatorAutoActivation` / Act, etc.).

**Path B relevance**: extremely high. This is a paper-faithful, freshly produced (Feb 2026) C++ translation of the same FORTRAN source family as our `13.f90`. The variable-name mapping table is a ready-made glossary; the file split is a structural template; the C++ code itself is small enough (under 4 K lines) to be re-translated to NumPy. This C++ port should be the first reference to consult before re-implementing equations from the paper.

The repository's variable mapping is saved locally at `/home/lyndon/repo/project/ToothmakerTranslation/.tmp/jernvall-toothmaker/cpp_README.md`, alongside copies of `humppa.f90`, `humppa_translate.f90`, and the seven C++ files.

### Vectorisation and translation efforts

#### tgrohens/toothmaker

<https://github.com/tgrohens/toothmaker>. Stars 0, forks 0, last push 16 Apr 2025, language Fortran. Description: 'Fork of the core functionality of the ToothMaker project.' Files: `coreop2d.f90` (1,635 lines), `esclec.f90` (392), `toothmaker.f90` (141), `Makefile`, `run.sh`, plus an example input `nosea_10261__0.txt`.

**This is the cleanest, smallest, most self-contained derivative** of the Salazar-Ciudad and Jernvall FORTRAN model. The combined `coreop2d.f90` + `esclec.f90` + `toothmaker.f90` total 2,168 lines — slightly more than `13.f90` (1,923) but with the same module structure. Subroutine names are still Catalan (not yet renamed). Useful as a minimal compilable reference and for cross-checking Path A behaviour.

**Path B relevance**: useful as the smallest verifiable FORTRAN reference. Its small size (only the 'core functionality') makes it the best starting point if you want to trace exactly which sections Path A's Python translation chose to include and which to drop.

Saved locally at `/home/lyndon/repo/project/ToothmakerTranslation/.tmp/tgrohens-toothmaker/`.

#### tgrohens/ToothEvol

<https://github.com/tgrohens/ToothEvol>. Stars 0, last push 7 Feb 2025, size 0 KB. **Empty placeholder** — repository description claims 'Run evolutionary simulations with ToothMaker' but no code is published. Path B relevance: none currently.

### Companion paper code

#### RolandZimm/silicoshark

<https://github.com/RolandZimm/silicoshark>. Stars 1, forks 0, single release v2.0 on 29 Nov 2022, language Fortran. The companion code for the Zimm 2023 PNAS shark paper. Files: `toothmodel.f90` (2,401 lines), `make_mesh_for_gnuplot.f90` (43 lines), `multicusp.ini` (parameter set), `tric_target.txt` and `pentac_target.txt` (CT-scanned shark tooth target outlines). README references the 2010 Salazar-Ciudad and Jernvall Nature paper as the source.

**Path B relevance**: medium. The FORTRAN is a sibling of `13.f90` (also 2,401 lines, like `humppa_translate.f90`), so structural diffs help calibrate which differences in our Python are intentional shark-extensions versus translation artefacts. The `multicusp.ini` is a self-documenting parameter file that makes the Egr / Mgr / Act / Inh / Sec / Deg / Bgr / etc. names concrete with values. This is currently the best human-readable parameter glossary outside the SI.

Saved locally at `/home/lyndon/repo/project/ToothmakerTranslation/.tmp/silicoshark/`.

#### RolandZimm/shark_heterodonty

<https://github.com/RolandZimm/shark_heterodonty>. Stars 0, last push 27 Nov 2025, MIT, no language detected (data only). Companion to a 2024 bioRxiv preprint (Zimm, Tobias-Santos, Goudemand) on multi-level dental diversity across sharks. Contains tooth half-outlines from 51 species, ecological trait data, and CCA correlation files. **No simulation code.** Path B relevance: none — this is shape data only.

### Downstream users (parameter optimisation, evolutionary studies)

#### jupander/BITES

<https://github.com/jupander/BITES>. Stars 0, last push 13 Mar 2026, language Python, Di-Poï Lab (University of Helsinki). 'Bayesian Inference for Tooth Emergence Simulation.' A 2026 GUI tool that wraps the FORTRAN tooth model in an Optuna-based Bayesian parameter optimiser, comparing simulated outputs against 3D tooth scans via Chamfer distance. Currently Windows-only.

Of particular note: `src/humppa_2_array_input_v8.f90` (101 KB, 2,868 lines) is a heavily modified FORTRAN that takes parameters as a delimited string argument rather than reading from a file, and has version comments documenting eight rounds of modification (`v1` … `v8`). Subroutines added include `parse_string_to_real_array`, `parse_string_to_char_array`, `assign_input_to_parap_and_noms`. This is useful if Path B wants to drive simulation runs programmatically without writing parameter files.

Cited paper: Razmadze et al. 2026, DOI 10.64898/2026.03.12.711251.

**Path B relevance**: medium. Not a Python re-implementation of the model itself, but an authoritative source for parameter names, value ranges, and shape-comparison metrics. The custom FORTRAN demonstrates how real users wrap this model for optimisation.

Saved locally at `/home/lyndon/repo/project/ToothmakerTranslation/.tmp/bites/humppa_2_array_input_v8.f90`.

#### sam-osullivan/tooth-evo-devo

<https://github.com/sam-osullivan/tooth-evo-devo>. Stars 1, last push 18 Jun 2025, language Python. Code for evolutionary studies of tooth complexity using the ToothMaker model. Files include `new_humppa_translate2.f90` (2,515 lines, another ToothMaker FORTRAN variant) plus Python scripts for OPC (orientation patch count), DNE (Dirichlet normal energy), RFI (relief index), Latin Hypercube Sampling, mutant generation, and `.off` to `.ply` conversion. Workflow involves cluster job submission via `addqueue` and per-batch processing of 19,000 simulated teeth.

**Path B relevance**: low for the model itself, high for downstream tooth-shape metrics. If Path B wants to validate against shape complexity measures (cusps, OPC, RFI), this is the authoritative open-source pipeline.

Saved locally at `/home/lyndon/repo/project/ToothmakerTranslation/.tmp/sam-osullivan/`.

#### sam-osullivan/Complexity-distributions-in-GP-maps

<https://github.com/sam-osullivan/Complexity-distributions-in-GP-maps>. Stars 0, last push 25 Mar 2026, language Jupyter Notebook. A research repository on genotype-phenotype maps including a `teeth/` subdirectory that calls back into `tooth-evo-devo`. Path B relevance: low.

#### jupander/Toothpic-browser

<https://github.com/jupander/Toothpic-browser>. Stars 0, language Python. 'Rudimentary tool for browsing screenshots of tooth shapes generated by ToothMaker.' Path B relevance: none — visualisation only.

### Unrelated namespace collisions

- `CUHK-AIM-Group/ToothMaker` (stars 6) — a 2024 deep-learning project for generating panoramic dental radiographs. Unrelated to the Salazar-Ciudad and Jernvall morphogenesis model.
- `msemon/ToothTranscriptomeAnalyses` — R code for transcriptome analysis. Unrelated.

### Summary table

| Repo | Stars | Last commit | Language | Path B value |
|------|-------|-------------|----------|--------------|
| jernvall-lab/ToothMaker | 3 | 2026-02-28 | C++ + FORTRAN | **High** — ground-truth FORTRAN ancestor and a fresh paper-faithful C++ port with full variable glossary |
| tgrohens/toothmaker | 0 | 2025-04-16 | FORTRAN | **High** — minimal self-contained FORTRAN, smallest verifiable reference |
| RolandZimm/silicoshark | 1 | 2022-11-29 | FORTRAN | **Medium** — Zimm 2023 companion code, `multicusp.ini` parameter glossary |
| jupander/BITES | 0 | 2026-03-13 | Python | **Medium** — modern wrapper, parameter-as-string FORTRAN, shape-comparison metrics |
| sam-osullivan/tooth-evo-devo | 1 | 2025-06-18 | Python | **Medium** — tooth complexity metrics (OPC, DNE, RFI) |
| RolandZimm/shark_heterodonty | 0 | 2025-11-27 | data | None — shape data only |
| tgrohens/ToothEvol | 0 | 2025-02-07 | empty | None — placeholder |
| sam-osullivan/Complexity-distributions-in-GP-maps | 0 | 2026-03-25 | Jupyter | Low — wraps tooth-evo-devo |
| jupander/Toothpic-browser | 0 | — | Python | None — viewer only |

### Subroutine name mapping for `13.f90`

Confirmed by line-by-line comment correspondence with `humppa_translate.f90` (2014 tribosphenic FORTRAN in `jernvall-lab/ToothMaker`):

| `13.f90` (English) | Catalan/Spanish original |
|---|---|
| `initialconditions` | `ciinicial` |
| `allocateinitialstate` | `dime` |
| `initializecellpositions` | `posar` |
| `calculatemargins` | `calculmarges` |
| `applydiffusion` | `reaccio_difusio` |
| `applydifferentiation` | `diferenciacio` |
| `EpGrowthBorderForce` | `empu` |
| `BoyForce` | `stelate` |
| `repulseneighbor` | `pushing` |
| `repelnonneigh` | `pushingnovei` |
| `applyborderbias` | `biaixbl` |
| `applynucleartraction` | `promig` |
| `updatecellposition` | `actualitza` |
| `addcell` | `afegircel` |
| `updatebordercells` | `perextrems` |
| `iteration` | `iteracio` |

Variable-name glossary (from the C++ port README):

| FORTRAN | C++ name | Description |
|---|---|---|
| `malla` | `cellPositions` | Cell node positions (x, y, z) |
| `vei` | `neighbors` | Neighbour indices per cell |
| `marge` | `cellMargins` | Internode margin positions |
| `knots` | `knotMarkers` | Knot markers (1 = knot, 0 = not) |
| `nveins` | `neighborCount` | Number of neighbours per cell |
| `hmalla` | `positionDeltas` | Position changes per iteration |
| `q3d` | `quantities3D` | 3D quantities `[cell][z][type]` for Act, Inh, Fgf, Ect, p |
| `q2d` | `quantities2D` | 2D quantities `[cell][type]` |
| `difq3d` | `diffusionCoeffs3D` | Diffusion coefficients (3D) |
| `difq2d` | `diffusionCoeffs2D` | Diffusion coefficients (2D) |

GUI parameter labels (from the C++ port README), mapping the cryptic FORTRAN names to the human-readable labels seen in `multicusp.ini`:

| FORTRAN | GUI label | Description |
|---|---|---|
| `tacre` | Egr | Epithelial proliferation rate |
| `tahor` | Mgr | Mesenchymal proliferation rate |
| `umgr` | uMgr | Basal Mgr (Sec-independent) |
| `acac` | Act | Activator auto-activation |
| `acec` | Not2 | Ectodin rate (hidden) |
| `acaca` | Not3 | Not3 rate (hidden) |
| `ihac` | Inh | Inhibition of activator |
| `ih` | Sec | Growth-factor secretion rate |
| `mu` | Deg | Protein degradation rate |
| `ina` | Ina | Initial activator concentration |
| `elas` | Rep | Young's modulus / stiffness |
| `crema` | Adh | Traction between neighbours |
| `radibi` | Ntr | Border-to-nucleus traction |
| `tazmax` | Dgr | Sharpness maxima |
| `tadi` | Swi | Border definition distance |
| `ud` | Set | Growth-factor threshold |
| `us` | Int | Initial inhibitor threshold |

### Recommendations for Path B

1. Treat `jernvall-lab/ToothMaker/models/tribosphenic/src/cpp/` as the primary structural reference. It is freshly authored, paper-faithful, and split along the same conceptual boundaries Path B will need (`core`, `geometry`, `diffusion`, `mechanics`, `division`, `iteration`, `file_io`).
2. Use `humppa_translate.f90` as the line-by-line ground truth where the C++ is ambiguous, and `tgrohens/toothmaker/coreop2d.f90` as a minimal cross-check to confirm a behaviour is in the canonical core rather than a 2014 tribosphenic addition.
3. Use `silicoshark/multicusp.ini` as the human-readable parameter file format reference. The labels there match the BITES / 2023 PNAS conventions.
4. Try once more for the SI PDF from a real browser session — `https://www.pnas.org/doi/suppl/10.1073/pnas.2216959120/suppl_file/pnas.2216959120.sapp.pdf` should serve directly when the request is not coming from a headless agent.
