# GenomeOS documents

GenomeOS decodes the full genome (its blocks, their function, their
structure) and mimics its behaviour through BioLang, BioIR and BioVM. Open
source (Apache 2.0 for the BioLang engine, AGPL-3.0-or-later for the
application), author Albert Canfield. The documents are layered from
wide to narrow; read them top down.

## Wide: what and why

| Document | Holds |
|---|---|
| [../README.md](../README.md) | the ambition, the central equation, what runs today, quick start, command list |
| [ROADMAP.md](ROADMAP.md) | the one plan: scope, architecture split, where each component stands, ten areas (A to J) with goals, requirements, gaps and next steps, data jobs, the ordered next steps (§5; item 14 is the current order), milestones, ideas pool, working rules |
| [ARCHITECTURE.md](ARCHITECTURE.md) | the stack, the state model, execution regimes, uncertainty as a type, context gating, resolution levels, the separable engine |
| [DECISIONS.md](DECISIONS.md) | every architecture decision, numbered and dated, with its consequence |
| [LANDSCAPE.md](LANDSCAPE.md) | what exists elsewhere and what GenomeOS adds |
| [GENOME-AS-CODE.md](GENOME-AS-CODE.md) | the reverse-engineering view: timers, blueprint, parts, systems |
| [../CONTRIBUTING.md](../CONTRIBUTING.md) | branches, the every-step cycle, `commit_own.sh` and `push_own.sh`, the push gate |
| [CODEBASE_MAP.md](CODEBASE_MAP.md) | the code, package by package: what each module does and where to start reading |
| [../LICENSING.md](../LICENSING.md) | the two licences, which paths each covers and why; data terms; dual licensing |
| [../ACKNOWLEDGEMENTS.md](../ACKNOWLEDGEMENTS.md) | every library, model and data source with its licence, and the terms of use |

## Specifications: the engine

| Document | Holds |
|---|---|
| [BIOLANG-GRAMMAR.md](BIOLANG-GRAMMAR.md) | the language and the IR as they are: the v0.4 grammar (26 block kinds) and the BioIR 0.4 types (25), generated from the parser and the dataclasses; a test fails when it is out of date |
| [BIOLANG-v0.1.md](BIOLANG-v0.1.md), [BIOLANG-v0.2.md](BIOLANG-v0.2.md), [BIOLANG-v0.3.md](BIOLANG-v0.3.md), [BIOLANG-v0.4-ECONOMY.md](BIOLANG-v0.4-ECONOMY.md) | the language design, one document per version (v0.3 adds the organism layer and experiments; v0.4 structure, economy and control: stages 1, 2 and 4 implemented, stage 3, energy and mitochondria, not) |
| [BIOIR-v0.1.md](BIOIR-v0.1.md) | the intermediate representation as first designed; the current types are in BIOLANG-GRAMMAR.md |
| [BIO-TOOLCHAIN.md](BIO-TOOLCHAIN.md) | `bio check | compile | run | test | repl`, self-testing programs, the REPL |
| [SEQUENCE-GRAMMAR.md](SEQUENCE-GRAMMAR.md) | the genome's own delimiters as a grammar |

## Designs: one per area

| Area | Document |
|---|---|
| Genome decoding | [GENOME-ANATOMY.md](GENOME-ANATOMY.md), [UNKNOWN.md](UNKNOWN.md), [NODES-READER-WRITER.md](NODES-READER-WRITER.md) |
| The 98% (attribution of the non-coding genome) | [ATTRIBUTION.md](ATTRIBUTION.md); the known-locus benchmark [LOCI-BENCHMARK.md](LOCI-BENCHMARK.md); the one independent prediction result [CRISPRI-RESULT.md](CRISPRI-RESULT.md) |
| The grammar by comparison (constraint on two axes, origin per gene, duplication, operators, the decompiled locus) | [GRAMMAR-BY-COMPARISON.md](GRAMMAR-BY-COMPARISON.md) |
| Molecules and flow | [PROTEIN.md](PROTEIN.md), [FLOW.md](FLOW.md) |
| Prediction (optional models) | [ALPHAGENOME.md](ALPHAGENOME.md) |
| Organism | [ORGANISM-FROM-ONE-CELL.md](ORGANISM-FROM-ONE-CELL.md) |
| Design | [DESIGN-MINIMAL-CELL.md](DESIGN-MINIMAL-CELL.md), [BIOFORGE-CONFIDENCE.md](BIOFORGE-CONFIDENCE.md) (why a design answer carries no stated probability) |
| Cancer and therapeutics | [CANCER.md](CANCER.md), [THERAPEUTICS.md](THERAPEUTICS.md) |
| Data, storage and jobs | [DATA.md](DATA.md), [STORAGE.md](STORAGE.md), [JOBS.md](JOBS.md) |

## Records

| Document | Holds |
|---|---|
| [PROGRESS.md](PROGRESS.md) | append-only log of finished work with its evidence; rendered in the Progress tab |
| [LESSONS.md](LESSONS.md) | what the code learned from data, and where it is embedded |
| [FALSE-RED-VERDICTS.md](FALSE-RED-VERDICTS.md) | the four check verdicts of 2026-10-02 that were statements about no tree, and what changed |
| [ACTION-PLAN.md](ACTION-PLAN.md) | the original phased plan of 2026-09-10, executed; kept as history |

Rules: a new area gets one design document here and one row in ROADMAP.md
§3; a decision that changes how things are built gets a row in
DECISIONS.md; a finished task gets a PROGRESS.md entry; a language change
gets a versioned spec.
