# Action plan: to an operative and functional GenomeOS

Written 2026-09-10 after the research pass (see LANDSCAPE.md and DATA.md).
The plan is ordered so that every phase produces something runnable and a
test that proves it against known biology.

## Where we are (updated after execution, 2026-09-10)

All five phases have been executed to a first working level; see PROGRESS.md for the per-task evidence. The core remains dependency-free; process-bigraph, alphagenome and biolearn are optional extras. What is still open: a real 30x BAM for telomere validation, an AlphaGenome API key for live predictions, the C. elegans lineage beyond the early embryo, and replacing in-house engines with libRoadRunner / MaBoSS / CompuCell3D once they ship Python 3.14 wheels.

At the start of execution the platform had zero external dependencies and 24 passing tests:

- Genome engine: any FASTA, real human chromosomes, `fetch(locus)`.
- BioIR v0.1 with evidence, confidence and `UNKNOWN`; JSON round-trip.
- BioLang v0.1: declarative modules compile with reference checking.
- BioVM: central dogma (validated on the human mitochondrial genome), Hill-function network runtime (repressilator oscillates), cell ageing runtime (four clocks, cited parameters).
- BioLib catalogue: 45 libraries across core / blueprint / timer / systems / parts.

## Guiding decisions

1. **Python + uv, permissive licences only in the core.** Every adopted engine below is Apache/MIT/BSD. Non-commercial AI weights (AlphaGenome, STATE) are used only through optional adapters.
2. **Do not build a scheduler.** Evaluate process-bigraph (Vivarium 2.0) in Phase 2 and adopt it unless it fails the spike. BioIR remains ours because it carries evidence and UNKNOWN.
3. **HG002 is the test human.** Open consent, best-characterised genome, T2T diploid assembly, benchmark variants.
4. **Every rule has a source.** No PR merges a rule with `evidence: none` outside a clearly labelled placeholder.
5. **Build upward.** No human-scale simulation until a cell-scale one is validated.

## Phase 1 (weeks 1–6): the compiler from public data

Goal: an individual's real genome becomes a BioIR module set automatically.

| # | Task | Output | Test |
|---|---|---|---|
| 1.1 | GFF3 reader (GENCODE 50) → `Gene`, `Transcript`, exons, CDS | `genomeos compile-annotation` | every MT- and chr21 protein-coding gene has a CDS that translates without internal stops |
| 1.2 | Random-access FASTA (write a faidx-style index; optional `pyfaidx`) | O(1) `fetch(locus)` on hg38 | fetch 10,000 random loci in under a second |
| 1.3 | VCF reader (GIAB HG002 benchmark) → apply SNVs/indels to reference → diploid `Genome` (maternal/paternal from phased GT) | `genomeos twin build` | HG002 chr21 haplotypes differ from reference at the expected variant count |
| 1.4 | Alternative start codons and mitochondrial gene table | correct MT-ND1/ND2/ND5 ORFs | 13 mtDNA proteins found at their published coordinates |
| 1.5 | GO and Reactome import → library membership tables replacing the hand catalogue | `genomeos libs` driven by data | catalogue genes are all members of their claimed GO terms |
| 1.6 | Variant effect on CDS (synonymous / missense / nonsense / frameshift) | `genomeos twin variants --gene X` | ClinVar pathogenic nonsense variants are classified nonsense |

Dependencies allowed in Phase 1: `pysam` or `pyfaidx`, `gffutils` or a hand-written GFF3 parser, `cyvcf2`. Keep the zero-dependency core importable without them.

## Phase 2 (weeks 6–12): cell runtime and composition spike

| # | Task | Output | Test |
|---|---|---|---|
| 2.1 | BioLang v0.2: `cell_type`, `event` (divide, differentiate, die) blocks, nested `transcript`, imports across files, a `bio.std` prelude | grammar + parser | existing v0.1 files still compile |
| 2.2 | `CellType` entities from Cell Ontology; context-gated rule sets per type | `genomeos check --context cell_type=X` | rule counts per cell type are reported |
| 2.3 | libRoadRunner adapter: run an SBML model from BioModels as a BioVM engine | `genomeos run model.xml` | BioModels reference trajectories reproduced |
| 2.4 | MaBoSS adapter for Boolean GRNs | large per-cell logic | a published Boolean model (e.g. T-cell differentiation) reproduces its attractors |
| 2.5 | **Spike:** wrap `NetworkRuntime` and `CellRuntime` as process-bigraph processes; run both in one composite with a shared clock | go / no-go decision documented in ARCHITECTURE.md | the composite reproduces standalone results |
| 2.6 | Ligand–receptor table (CellPhoneDB) → `Rule`s with `when: tissue` | the tissue dependency graph | graph reproduces known signalling axes (e.g. liver ↔ adipose) |

## Phase 3 (months 3–6): twin and timers on real data

| # | Task | Output | Test |
|---|---|---|---|
| 3.1 | Epigenetic age from a methylation matrix via `biolearn` or `pyaging`; store in `CellState` | `genomeos twin age --methylation file` | public GEO samples yield published clock ages |
| 3.2 | Telomere length from WGS reads via Telomerecat / TelomereHunter2 | measured `telomere_bp` for HG002 | within the published range for a male in his 40s at sampling |
| 3.3 | BioTwin = genome + measured state + environment; `fork`, `run`, `diff` | `genomeos twin fork/run/diff` | a variant in a `timer.*` gene shifts the simulated trajectory in the documented direction |
| 3.4 | AlphaGenome adapter (hosted API) producing `predicted` rules with model name and score | regulatory rules for UNKNOWN regions | predicted expression direction agrees with GTEx eQTLs for known variants |
| 3.5 | Uncertainty report per level (molecular / cellular / tissue / organism) in every run | part of `genomeos run` output | report is present and monotonic with evidence kinds |

## Phase 4 (months 6–12): space and development

| # | Task | Output | Test |
|---|---|---|---|
| 4.1 | Spatial adapter (CompuCell3D first; Tissue Forge fallback) | cells with position, diffusion fields | a gradient-reading rule produces a French-flag pattern |
| 4.2 | Segmentation clock module in BioLang: HES7/LFNG/Notch + moving FGF/WNT gradient | segments form in simulation | segment count and period within the human-cell range |
| 4.3 | Gastrulation and germ-layer module with master factors | three-population outcome | proportions match literature within stated confidence |
| 4.4 | Debugger: step, breakpoints on conditions, evidence trace | `genomeos debug` | a differentiation event is explained back to its rules |

## Phase 5 (year 2): organism scale and BioForge

Minimal organism run (C. elegans scale, real genome + bootstrap state);
resolution levels per tissue; BioForge design-under-constraints; external
validation by reproducing published perturbation experiments.

## UI track (runs alongside every phase)

v0.1 ships `genomeos serve`: a localhost web UI with four views (Genome,
Program, Ageing, Libraries), stdlib server, one HTML file, inline SVG charts.
Each phase adds a view rather than a separate tool:

| Phase | View | What it shows |
|---|---|---|
| 1 | **Genes** | GENCODE gene track over a chromosome; click a gene to see transcripts, CDS translation, and the individual's variants applied (HG002) |
| 1 | **Libraries (data-driven)** | membership counts from GO/Reactome; click a gene to see which libraries it belongs to and which cell types express it |
| 2 | **Cell** | a cell type's active rule set, its regulatory graph drawn as nodes and edges, run and watch expression |
| 2 | **Evidence explorer** | filter every rule by evidence kind and confidence; export the weak ones as a review list |
| 3 | **Twin** | genome + measured state; fork, edit a variant or the environment, run both, show the diff per clock |
| 4 | **Space** | 2D cell field with morphogen gradients and segment formation |
| 4 | **Debugger** | step, breakpoints on conditions, and the trace from an event back to its rules and sources |

## Team and infrastructure

- One developer can execute Phases 1–2 on a laptop. Phase 3.4 needs an API key; 3.2 needs a 30x BAM (~100 GB) and CPU hours. Phase 4 benefits from a GPU but does not require one.
- A biologist reviewer from Phase 2 onward: every library module gets a "known phenotypes it must reproduce" test list written by someone who knows the tissue.
- CI: `uv run pytest` on every push; a nightly job that runs the real-data tests against `data/reference`.

## Risks

| Risk | Mitigation |
|---|---|
| process-bigraph API instability | pin a version; keep adapters thin; BioIR is the stable contract |
| Non-commercial AI licences | adapters optional; core never depends on them |
| Over-claiming biology | evidence + confidence on every rule; the uncertainty report is mandatory output |
| Scope creep toward "simulate a human" | roadmap gates: no level starts before the previous is validated |
| Hand-curated gene lists drift from data | Phase 1.5 replaces them with computed membership; the catalogue becomes a test, not a source |

## First five concrete tickets

1. GFF3 parser for GENCODE 50 → `Gene`/`Transcript` entities (1.1).
2. FASTA index for random access on hg38 (1.2).
3. VCF apply → diploid HG002 chr21 (1.3).
4. Alternative start codons; all 13 mitochondrial proteins (1.4).
5. process-bigraph spike with the two existing runtimes (2.5, can run in parallel).
