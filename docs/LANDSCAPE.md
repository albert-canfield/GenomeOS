# Landscape (September 2026)

What already exists, what each piece covers, and what GenomeOS adds. Compiled
from a research pass on 2026-09-10; URLs are the primary sources.

## 1. Composition and whole-cell runtimes

| Project | What it is | Verdict for GenomeOS |
|---|---|---|
| [process-bigraph](https://github.com/vivarium-collective/process-bigraph) (Vivarium 2.0, Apache-2.0, Python) | Typed hierarchical state, processes with ports, multi-timescale scheduling, structural rewrites. vivarium-core 1.x is deprecated in its favour. | **Strongest candidate for the BioVM scheduler.** Expose GenomeOS engines as processes. |
| [vEcoli](https://github.com/CovertLab/vEcoli) (MIT) | Covert lab E. coli whole-cell model on Vivarium processes. | Proof that the composition approach scales to a whole cell. Reference architecture. |
| [Minimal_Cell_4DWCM](https://github.com/Luthey-Schulten-Lab/Minimal_Cell_4DWCM) (Cell, Apr 2026) | Spatial 4D simulation of JCVI-syn3A's full cell cycle: RDME + CME + ODE + chromosome mechanics. | The high-water mark for a cell. GPU-bound, research code. Target to compare against, not embed. |
| [Lattice Microbes](https://github.com/Luthey-Schulten-Lab/Lattice_Microbes) | GPU stochastic reaction-diffusion solver. | Later, for spatial molecular events. Needs NVIDIA GPU. |

## 2. Network simulators (molecule to cell)

| Project | Regime | Python | Verdict |
|---|---|---|---|
| [libRoadRunner / Tellurium](https://github.com/sys-bio/roadrunner) v2.10 (Aug 2026, Apache-2.0) | ODE, SSA, steady state; JIT-compiled SBML; Antimony language | native | **Adopt** as the ODE/SSA engine for imported SBML and physiology modules. |
| [COPASI / basico](https://github.com/copasi/basico) 4.47 | ODE, SSA, parameter estimation | native | Use for parameter fitting. |
| [PySB](https://github.com/pysb/pysb) 1.17 / [BioNetGen](https://github.com/RuleWorld/bionetgen) | rule-based molecular events, combinatorial complexes | native | **Adopt** for the stochastic molecular regime. Its Python DSL is the closest cousin of BioLang. |
| [MaBoSS](https://github.com/sysbio-curie/MaBoSS) 2.6 | continuous-time Boolean GRNs | pyMaBoSS | **Adopt** for large per-cell regulatory logic. Already co-simulates inside PhysiCell. |
| libSBML 5.21, libCellML 0.7 | format libraries | yes | Import paths for the compiler. |

## 3. Spatial and tissue

| Project | Regime | Python | Verdict |
|---|---|---|---|
| [CompuCell3D](https://github.com/CompuCell3D/CompuCell3D) 4.9 (MIT) | Cellular Potts + PDE fields + per-cell SBML; GPU solvers in 4.10 | native | Leading candidate for the development engine. |
| [Tissue Forge](https://github.com/tissue-forge/tissue-forge) 0.2 (LGPL-3) | particle-based, molecular to multicellular | native, Jupyter | Lighter alternative; good for prototypes. |
| [PhysiCell](https://github.com/MathCancer/PhysiCell) 1.14 (BSD-3) | off-lattice agent-based + diffusion; CSV rules grammar | no core bindings | Scientifically strongest ABM; its "rules" grammar is worth studying for BioLang. |
| Morpheus 2.4, Chaste 2026.1, BioDynaMo | CPM/ODE/PDE; cardiac + cell-based FE; HPC ABM | limited/none | Reference points, not foundations. |

## 4. Whole organism

[OpenWorm](https://github.com/openworm/OpenWorm) 0.9.8: c302 (NeuroML connectome), Sibernetic (SPH body), Geppetto. Low velocity, stitched by scripts; no unified organism runtime exists anywhere. C. elegans (959 somatic cells) remains the right intermediate target for a genome-to-organism run.

## 5. Design languages and compilers

- **SBOL 3** ([pySBOL3](https://github.com/SynBioDex/pySBOL3)): formal data model for genetic designs. Import/export target for BioLang designs.
- **Cello**: compiles Verilog logic into genetic circuits. The existence proof that a biological compiler works, at circuit scale.
- **PhysiCell rules CSV** and **Antimony**: two small, successful biology DSLs. BioLang should stay at least that simple.

## 6. AI models (the "predicted" evidence kind)

| Model | Input → output | Access | Licence |
|---|---|---|---|
| [AlphaGenome](https://github.com/google-deepmind/alphagenome) + [Atlas](https://alphagenome.google/atlas) (Sep 2026) | ≤1 Mb DNA → expression, splicing, chromatin, contacts; Atlas: all ~9 billion SNVs scored | hosted API (`pip install alphagenome`), weights on Kaggle/HF | non-commercial for now |
| [Evo 2](https://github.com/ArcInstitute/evo2) (Nature, Mar 2026) | DNA ≤1 Mb → likelihoods, embeddings, generation | HF weights 1B/7B/40B; NVIDIA NIM | Apache-2.0 |
| [Arc STATE](https://github.com/ArcInstitute/state) | single-cell expression + perturbation → post-perturbation expression | `pip install arc-state`, HF weights | non-commercial |
| Virtual Cell Challenge 2026 | zero-shot prediction in six unperturbed cell lines; deadline 5 Nov 2026 | — | — |
| Geneformer, scGPT, C2S-Scale, Tahoe-x1, TranscriptFormer | cell-level foundation models | HF | permissive |
| Enformer, Borzoi, Caduceus, HyenaDNA, NTv3 | sequence-to-function | HF | mostly permissive; NTv3 non-commercial |
| [biolearn](https://github.com/bio-learn/biolearn), [pyaging](https://github.com/rsinghlab/pyaging) | methylation → epigenetic age (Horvath, GrimAge, DunedinPACE, ...) | pip, CPU | open |
| Telomerecat, TelomereHunter2 | WGS BAM → telomere length | pip, CPU | open |

## 7. What is missing, and therefore what GenomeOS is

No system ties a **specific genome** to **reusable executable biology** with
**evidence and uncertainty as first-class values**, runs it **upward through
cells and development**, and lets you **fork and compare**. The pieces above
cover every layer individually. GenomeOS's contribution is the language
(BioLang), the representation (BioIR, with evidence and UNKNOWN), and the
integration contract that lets those engines run one organism together.
