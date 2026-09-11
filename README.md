# GenomeOS

**An executable model of biology.**

**Ambition.** GenomeOS is built to become the software geneticists reach
for first: to work, research, investigate and build. A real bio-engineering
tool that understands the biological code rather than a viewer of it. Every
feature is judged by that bar: does it help a geneticist decode, model or
build, and does it say how sure it is.

GenomeOS treats a genome as source material and biology as a runtime. It is not a
DNA editor and not a genome browser. It is the layer that sits *between* sequence
and simulation: a compiler that turns a real genome plus scientific knowledge into
an executable model, and a virtual machine that runs that model through time.

```
Genome source          Biological knowledge
(FASTA / VCF / GFF3)   (GO, Reactome, Uberon, SBML, CellML, AlphaGenome, ...)
        │                        │
        └──────────┬─────────────┘
                   ▼
            Biological Compiler
                   │
                   ▼
              BioIR  (intermediate representation)
                   │
                   ▼
              BioVM  (runtime)
                   │
                   ▼
   cell state ──▶ tissue ──▶ organ ──▶ organism   over time
```

## The central equation

```
S(t + Δt) = Runtime(G, S(t), E, Δt)
```

| Symbol | Meaning | Analogy |
|---|---|---|
| `G` | genome (mostly constant over a lifetime) | the program |
| `S` | biological state: epigenome, telomeres, proteins, damage, cell types | the save file |
| `E` | environment: nutrition, signals, toxins, temperature | input / IO |
| `Runtime` | biological semantics (compiled rules, evidence, uncertainty) | the VM |

A genome alone cannot run. A zygote inherits both information (`G`) and state
(`S0`: ribosomes, mitochondria, maternal RNAs, spatial organisation). GenomeOS
therefore always executes `Genome + Bootstrap State`, never DNA in isolation.

## Project family

| Component | Role | Status (2026-09-11) |
|---|---|---|
| **BioLang** | biological programming language / DSL | v0.3: 18 block kinds from `gene` to `organism` and `experiment`; evidence on every block; programs test themselves |
| **BioIR** | common intermediate representation joining genomics, pathways, cells, physiology | 0.3: 20 types with JSON round-trip, `UNKNOWN` as a value |
| **BioVM** | runtime: discrete events, stochastic, continuous, spatial engines | eleven engines (central dogma, network, Boolean, SBML, ageing, spatial, segmentation, gastrulation, Body, debugger, uncertainty); process-bigraph composition |
| **`bio`** | the toolchain on its own | `bio check | compile | run | test | repl` |
| **BioLib** | reusable biological modules (`cell.core`, `animal.development`, `organs.heart`, ...) | 45 libraries with data-computed membership from GO and Reactome |
| **BioTwin** | a specific genome + specific biological state | HG002 with measured state; fork, run, diff; genome-wide build in progress |
| **BioForge** | in-silico design and experiment system | design search under constraints; experiments as input next |

Where each stands in detail, what is missing and what comes next:
[docs/ROADMAP.md](docs/ROADMAP.md) §2.1; the document index is
[docs/README.md](docs/README.md).

## What runs today (v0.1)

- **Genome engine.** Stream any FASTA (plain or gzip), including full human
  chromosomes, into typed `Chromosome` / `Sequence` / `Locus` objects.
- **Central dogma runtime.** Transcription with strand awareness, splicing from
  exon loci, translation with the standard and vertebrate-mitochondrial codon
  tables, open-reading-frame discovery.
- **Cell ageing runtime.** A stochastic per-cell model with telomere attrition,
  somatic mutation accumulation, epigenetic-clock drift and senescence. Every
  parameter carries its evidence and confidence, so the output reports how much
  of it is measured versus guessed.
- **BioIR v0.1 types.** Entities, loci, rules with context, evidence levels,
  confidence, and an explicit `UNKNOWN` role. See [docs/BIOIR-v0.1.md](docs/BIOIR-v0.1.md).
- **BioLang v0.1.** A declarative language that compiles to BioIR with
  compile-time reference checks. The repressilator is the first program.
- **BioLib catalogue.** 45 biological "libraries" found in the genome, in five
  layers: core, blueprint, timer, systems, parts. `genomeos libs -v`.
- **Therapeutic target reasoning.** From a tumour VCF: protein localisation and
  topology, whether a circulating binder can reach the altered protein, healthy
  versus tumour expression, trafficking after binding, the mutant peptide/HLA
  route for intracellular proteins, and mechanism compatibility across antibody,
  engager, conjugate, radionuclide, TCR and experimental delivery mechanisms.
  Every conclusion carries its evidence level; everything unestablished stays
  unknown. See [docs/THERAPEUTICS.md](docs/THERAPEUTICS.md).

Validated on real human DNA: the mitochondrial genome translates to the known
proteins (MT-CO1 513 aa, MT-CO2 227 aa, MT-ATP6 226 aa) and chromosome 21
loads in half a second.

## Status (2026-09-10)

The action plan in [docs/ACTION-PLAN.md](docs/ACTION-PLAN.md) has been executed
through all five phases to a first working level; [docs/PROGRESS.md](docs/PROGRESS.md)
records the evidence per task. Highlights:

- a real person's genome (GIAB HG002) compiles to BioIR: 55,210 chromosome-21 variants applied to both haplotypes, every coding transcript translated, ClinVar consequences reproduced at >90%
- 45 genome "libraries" verified against Gene Ontology and Reactome by data (95.5% agreement, the gap documented)
- SBML models from BioModels and Boolean models run unchanged; the Fauré cell cycle gives its published attractors
- epigenetic clocks validated on real blood methylation (Horvath r > 0.8); telomere estimator; digital twins that fork, run and diff
- development: French flag, segmentation clock at the human 5-hour period, germ layers along a NODAL gradient
- an organism from one cell: BioLang v0.3 `organism` programs run by the Body runtime (`genomeos grow`); C. elegans grows from the zygote to the 959-cell adult with every cell, fate and programmed death matching Sulston's lineage, founders decided by maternal factors, Wnt and Notch, timing within a median 12 minutes ([docs/BIOLANG-v0.3.md](docs/BIOLANG-v0.3.md), [docs/ORGANISM-FROM-ONE-CELL.md](docs/ORGANISM-FROM-ONE-CELL.md))
- knockouts as experiments: `experiment` blocks remove a maternal factor or a signal and run against the wild type; the six classic C. elegans mutants (pop-1, skn-1, pie-1, apx-1, glp-1, pal-1) reproduce their published founder phenotypes from the cited rules (`genomeos grow data/organisms/celegans/mutants.bio --experiments`)
- the same runtime at population resolution grows a human body: Carnegie-stage timers, germ layers by share, sixteen tissue populations that reach the adult counts of Sender & Milo 2021 and turn over at their measured rates (2.8e13 cells at 20 years, 3.2e11 replaced per day), every number cited and the whole labelled low confidence because it is counts, not mechanism (`genomeos grow data/organisms/human/body.bio --until "20 yr"`)
- a debugger with biological breakpoints and evidence traces, and BioForge design search whose outputs are labelled predicted
- the DNA → RNA → protein relationship as one traceable object: every base mapped to its mRNA position, codon and residue and back, with the consequence of a single-base change (Flow tab, `genomeos flow`)
- a federated protein compiler: Ensembl, UniProt, InterPro, PDB, AlphaFold, Reactome, STRING and the Human Protein Atlas compiled into one definition per protein keyed by UniProt accession, evidence and confidence per section, predicted never treated as observed (`genomeos protein X --compile`, docs/PROTEIN.md)
- nodes above genes: domains between CTCF boundaries with the genes and enhancers inside (`genomeos domains`)
- therapeutic target reasoning: for every tumour alteration, where the protein sits and whether a binder can physically reach it, how the tumour differs from healthy tissue, what happens after binding, the peptide/HLA route, and which of 15 therapeutic mechanisms the biology supports, each with its evidence level and its missing data (`genomeos therapeutic`, docs/THERAPEUTICS.md)
- a local web UI with a dozen views, including a Progress tab with live background jobs and the project's task log

## Quick start

```bash
uv sync
uv run genomeos info data/demo/demo.fa
uv run genomeos orfs data/demo/demo.fa --min-aa 30
uv run genomeos age --cell-type fibroblast --years 90 --cells 500 --evidence
uv run genomeos libs --layer timer -v
uv run genomeos grow data/organisms/celegans/embryo.bio --until 6000 --compare   # one cell to the adult worm
uv run genomeos grow data/organisms/human/body.bio --until "20 yr"               # a human body as populations
uv run genomeos serve --open           # light web UI at http://127.0.0.1:8765
uv run bio test data/demo              # the BioLang toolchain on its own: check, compile, run, test, repl
uv run genomeos therapeutic --tumour data/demo/cancer_tumour.vcf --report
uv run pytest                          # ~80 tests; real-data tests skip until data is fetched
uv run ruff check .
```

More commands: `annotate`, `gene`, `index`, `variant`, `twin build|new|fork|run`,
`clock`, `telomere`, `cells`, `lr`, `debug`, `develop`, `forge`, `organism`, `grow`,
`flow`, `regulation`, `domains`, `repeats`, `rna`, `protein --compile|--bio`, `proteome`,
`graph`, `pathway`, `lookup`, `report`, `verify`, `unknown`, `cancer compare|tumour|expression`, `therapeutic`,
`design`, `data fetch --chrom C [--individual]`; the `bio` toolchain: `bio check|compile|run|test|repl`.
Optional extras: `uv sync --extra compose` (process-bigraph), `--extra predict`
(AlphaGenome client, needs `ALPHAGENOME_API_KEY`), `--extra clocks` (biolearn, needs torch).

The web UI is a single HTML page served by the standard library: Start,
Genome, Anatomy, Blocks (2-D block map with domains, genes, UNKNOWN classes
and ENCODE elements), Flow (DNA → RNA → protein, linked lanes), Program,
Ageing, Twin, Space, Debugger, Molecules (protein definition and AlphaFold
structure), Cancer, Targets (therapeutic candidates, mechanisms and the
evidence behind each), Progress and Libraries.

To work with a real human chromosome:

```bash
scripts/fetch_reference.sh chrM      # human mitochondrial genome, 6 KB
scripts/fetch_reference.sh chr21     # hg38 chromosome 21, 13 MB
uv run genomeos orfs data/reference/chrM.fa.gz --table mito --min-aa 200
uv run genomeos info data/reference/chr21.fa.gz
scripts/fetch_reference.sh hg002     # a real individual's diploid genome (open consent), ~1 GB
```
See [docs/DATA.md](docs/DATA.md) for every open dataset and its verified URL.

## Storage principle

Stream, distil, discard: raw data is streamed or downloaded once, every real-data run leaves a small committed summary under `data/results/`, and the raw file is deleted (`genomeos data status|distil|clean`). See [docs/STORAGE.md](docs/STORAGE.md).

## Documents

- [docs/GENOME-AS-CODE.md](docs/GENOME-AS-CODE.md) — the reverse-engineering guide: timers, blueprint, parts list, system integration
- [docs/ACTION-PLAN.md](docs/ACTION-PLAN.md) — phased plan to an operative platform, with tests per task
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — layers, engines, state model, uncertainty as a type
- [docs/BIOIR-v0.1.md](docs/BIOIR-v0.1.md) — the intermediate representation
- [docs/BIOLANG-v0.1.md](docs/BIOLANG-v0.1.md) — the language
- [docs/LANDSCAPE.md](docs/LANDSCAPE.md) — existing systems (research pass, Sept 2026) and what GenomeOS adds
- [docs/DATA.md](docs/DATA.md) — open genome data for testing, URLs verified
- [docs/STORAGE.md](docs/STORAGE.md) — stream, distil, discard: keeping the project small on disk
- [docs/LESSONS.md](docs/LESSONS.md) — what the core has learned from data, and where it is embedded
- [docs/SEQUENCE-GRAMMAR.md](docs/SEQUENCE-GRAMMAR.md) — the genome's own delimiters, measured, and how BioLang builds on them
- [docs/GENOME-ANATOMY.md](docs/GENOME-ANATOMY.md) — counting the blocks of three genomes, and the budgets for building one from zero
- [docs/FLOW.md](docs/FLOW.md) — the upward flow (DNA → regulation → RNA → protein → cell → tissue → organism) and the types that carry it
- [docs/NODES-READER-WRITER.md](docs/NODES-READER-WRITER.md) — nodes, reader, writer, executor: the conceptual model against the biology
- [docs/UNKNOWN.md](docs/UNKNOWN.md) — classifying the space between genes from sequence alone, with a pattern file you can extend
- [docs/CANCER.md](docs/CANCER.md) — healthy versus tumour: somatic differences, drivers from cBioPortal, targets, agent packet
- [docs/THERAPEUTICS.md](docs/THERAPEUTICS.md) — therapeutic targeting: accessibility, selectivity, trafficking, neoantigens, mechanism compatibility, and the design dataset
- [docs/DESIGN-MINIMAL-CELL.md](docs/DESIGN-MINIMAL-CELL.md) — how much genome a neuron needs, and why it still is not a small extract
- [docs/ROADMAP.md](docs/ROADMAP.md) — from a molecular cell to a digital twin

## Goals

1. **Decode.** Take real genomes and let their blocks, elements and rules
   emerge from the data (anatomy, signals, libraries), with evidence and
   uncertainty on everything.
2. **Model.** Run biology forward through time, from a specific genome and a
   measured state: cells, ageing, development, twins.
3. **Build.** Design under constraints in silico, and eventually produce the
   budgets and layouts a genome built from zero would need.
4. **BioLang as an engine.** Today BioLang is one part of GenomeOS, developed
   here. The long-term goal is for BioLang to stand on its own the way Node
   stands for JavaScript: an engine and toolchain people use to write, script,
   compile and run biology, with GenomeOS as its first application. How to get
   there is open; the architecture keeps the boundaries that make it possible
   (see [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md), "BioLang as a separable engine").

## Principles

1. **Lossless genome.** The whole sequence is kept, not a gene list. Internally
   the design target is a pangenome graph, not one linear reference.
2. **Genome and state are separate.** Age, methylation, telomeres and damage are
   state, not sequence.
3. **Unknown is a type.** Regions and rules without evidence are marked
   `UNKNOWN`, never invented.
4. **Evidence is executable.** Every rule carries provenance and confidence, and
   simulation output reports confidence per level (molecular, cellular, tissue, organism).
5. **Import, don't reinvent.** FASTA, VCF, GFF3, SBML, CellML, GO, Uberon,
   Reactome are inputs to the compiler.
6. **Build upward.** Molecule, cell, minimal organism, C. elegans-scale, vertebrate
   subsystems, then human digital twins. Never "simulate a human" as step one.

## Scope note

GenomeOS is a computational platform: design, simulate, analyse. Turning designs
into real organisms is a separate, safety-sensitive discipline outside this project.
