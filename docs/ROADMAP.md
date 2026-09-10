# Roadmap

Status after the 2026-09-10 build: ✅ done, ◐ partial. Details per task in [PROGRESS.md](PROGRESS.md).

Build upward. Each level must run, be tested against known biology, and
report its confidence before the next begins.

| Version | Target | Proof it works |
|---|---|---|
| **0.1** | Genome engine, BioIR, BioLang v0.1, central dogma, network runtime, cell ageing, library catalogue | mitochondrial genome translates to the known proteins; repressilator oscillates; ageing clocks match cohort attrition; 24 tests |
| 0.2 ✅ | **Compiler from public data**: GENCODE → genes/transcripts, VCF → an individual's variants applied to the reference, GO/Reactome → library membership and rules. Diploid `Genome`. | HG002's genome loads; every protein-coding gene has a `Gene` entity; `genomeos libs` counts come from data, not the hand catalogue |
| 0.3 ✅ (in-house engines) | **Cell runtime**: `CellType`, differentiation and division events, context-gated rule sets per cell type, import of SBML models via libRoadRunner, MaBoSS for large GRNs | a hematopoietic differentiation module reproduces known lineage choices; a Reactome pathway runs unchanged |
| 0.4 ✅ spike | **Composition**: BioVM engines as process-bigraph processes; one organism = a composite; variable resolution per tissue | the ageing runtime, a GRN and an SBML metabolism model run in one composite with a shared clock |
| 0.5 ✅ (1-D/2-D in-house) | **Multicellular and spatial**: CompuCell3D or Tissue Forge adapter; morphogen gradients; segmentation clock forms segments | a 2D somite-like pattern emerges from HES7/Notch rules plus a gradient |
| 0.6 ◐ early embryo | **Minimal organism**: a C. elegans-scale developmental run from a real genome plus bootstrap state | lineage tree matches known early divisions within stated confidence |
| 0.7 ✅ | **BioTwin**: a person's genome + measured state (methylation age, telomere length from reads) + environment; fork and compare | twin A versus twin B differ in the predicted direction for a known variant |
| 0.8 ✅ debugger | **Debugger and IDE**: step, breakpoints on biological conditions, evidence trace | a differentiation event can be explained back to its rules and sources |
| 1.0 ◐ search only | **BioForge**: design under constraints, in-silico experiments, confidence-scored predictions | an external group reproduces a published perturbation result from a BioLang module |

| 1.x | **BioLang as an engine**: the language, IR, VM and standard library packaged on their own (`bio run | check | compile | test | repl`), embeddable, with GenomeOS as the first application | a `.bio` program runs without GenomeOS installed |

Not on the roadmap: wet-lab synthesis or any translation of designs into
real organisms. GenomeOS stops at simulation and analysis.
