# Architecture decisions

The decisions that shape GenomeOS, numbered, dated, with their consequence.
A decision stays until a later entry supersedes it. New decisions are
appended; a superseded one keeps its number and gains a "superseded by"
line. Wide-level decisions come first; narrower ones follow by area.

## Project

| # | Date | Decision | Consequence |
|---|---|---|---|
| D1 | 2026-09-10 | GenomeOS is open source under the MIT licence, authored by Albert Canfield; the core has zero dependencies and admits permissive licences only (Apache, MIT, BSD). Non-commercial AI weights (AlphaGenome, STATE) enter only through optional adapters | `pyproject` has no runtime dependencies; process-bigraph, biolearn and alphagenome are optional extras |
| D2 | 2026-09-10 | The name is GenomeOS; GenCode was rejected because GENCODE is the annotation consortium | |
| D3 | 2026-09-10 | The ambition: decode the full genome (its blocks, their function, their structure) and mimic its behaviour through BioLang, BioIR and BioVM; the organism is a functional system and DNA is its seed; timers, blueprint, parts and system connections stay first-class | every feature is judged by whether it helps a geneticist decode, model or build, and says how sure it is |
| D4 | 2026-09-10 | Every rule and parameter carries evidence and a confidence; `UNKNOWN` is a legal value; the uncertainty report per level is mandatory output; nothing is filled in | the `Evidence` type in BioIR; `genomeos check`; the report printed by every run |
| D5 | 2026-09-10 | Genome and state are separate: a run is always Genome + bootstrap state, never DNA alone | `S(t + Δt) = Runtime(G, S(t), E, Δt)`; twins hold measured state apart from the genome |
| D6 | 2026-09-10 | HG002 (GIAB, open consent) is the test human | all individual-level tests use HG002 |
| D7 | 2026-09-10 | Stream, distil, discard: raw downloads never enter the repository; one real-data run keeps one summary in `data/results`; the engine embeds what was learned | `genomeos data status / distil / clean`; the CI fetches 60 MB and reads summaries for the rest |
| D8 | 2026-09-10 | No wet-lab synthesis, no clinical advice, and (2026-09-11) no therapeutic sequences, constructs, vectors, formulations, dosing or protocols; GenomeOS stops at describing what recognition is required | the therapeutic dataset is checked to contain no nucleotide run |
| D9 | 2026-09-10 | Build upward: no level starts before the previous one is validated against known biology | roadmap gates; mtDNA proteins before chr21, chr21 before the genome, worm before human |
| D10 | 2026-09-10 | All work on `dev`; Albert opens the pull requests to `main`; commits carry no tool attribution | CONTRIBUTING.md; sessions never run `scripts/promote.sh` |
| D11 | 2026-09-11 | Several sessions share one checkout: commit through a private index, own files only, format only own files, announce hunks in shared files | LESSONS.md "Several sessions in one checkout"; ROADMAP.md §7 |

## Language, IR and engine

| # | Date | Decision | Consequence |
|---|---|---|---|
| D12 | 2026-09-10 | BioLang, BioIR, BioVM and the standard library are a separable engine; GenomeOS is its first application. Boundaries: BioIR JSON is the contract, engines take a Module and return a typed result, the standard library is written in BioLang, one thin `bio` entry point, no dependencies in the core, engines import nothing from genome, knowledge or web | ARCHITECTURE.md §10; ROADMAP.md §2 lists the two files that still break the last rule |
| D13 | 2026-09-10 | Do not build a scheduler: process-bigraph (Vivarium 2.0) is the composition layer, adopted after the spike passed | `runtime/compose.py`; updates are additive deltas registered with `core.register_link` |
| D14 | 2026-09-10 | In-house SBML, Boolean and spatial engines stand in until libRoadRunner, MaBoSS and CompuCell3D ship Python 3.14 wheels; they keep the interfaces so adapters can replace them | `runtime/sbml.py`, `boolean.py`, `spatial.py` |
| D15 | 2026-09-10 | The genome works in nodes: TADs bounded by CTCF are the nodes, the epigenome is the reader, replication / epigenetic drift / germline are the writers, the runtime is the executor | NODES-READER-WRITER.md; `Domain` block in BioIR 0.3; `genomeos domains`, `genomeos reader` |
| D16 | 2026-09-11 | Mechanism where it is known, the observed program elsewhere, both cited; the uncertainty report says which is which per level; timers run with their measured spread by default | `founders.bio` versus `lineage_*.bio`; the Body runtime |
| D17 | 2026-09-11 | Programs test themselves: `# test:` lines and `assert:` clauses are the program's claims, run by `bio test` in CI | BIO-TOOLCHAIN.md; `.github/workflows/ci.yml` |
| D18 | 2026-09-11 | Experiments are a language construct: an `experiment` block perturbs (knockouts, absent factors) and compares against the wild type; BioForge will take experiments as input | `genomeos grow FILE --experiments`; `mutants.bio` |

## Genome decoding

| # | Date | Decision | Consequence |
|---|---|---|---|
| D19 | 2026-09-10 | The genome's block delimiters are statistical signals, not tokens; segments are parsed by a grammar with context and never asserted from one signal | SEQUENCE-GRAMMAR.md; `genomeos signals` reports relative scores |
| D20 | 2026-09-11 | UNKNOWN blocks are classified in a fixed order, curated repeats (RepeatMasker) and chromatin evidence (ENCODE) before sequence-only classes such as long ORFs; sequence signatures are the fallback for a chromosome not yet distilled | UNKNOWN.md; the `long_orf` class shrank 17-fold on chr21 once repeats came first |
| D21 | 2026-09-11 | Enhancer targets are inferred within the CTCF domain (nearest coding genes first) and labelled inferred until Hi-C or a predictive model measures them | `genomeos regulation`; confidences 0.4 / 0.25 |
| D22 | 2026-09-11 | A fetched chromosome arrives analysed: sequence, models, elements, repeats, the curated UNKNOWN pass and domains in one job | `genomeos data fetch --chrom C --analyse`; the `fetch_<chrom>` jobs |

## Molecules

| # | Date | Decision | Consequence |
|---|---|---|---|
| D23 | 2026-09-10 | Proteins are compiled by a federated compiler from public sources, never kept as an own database; the UniProt accession is the protein id; the model is Gene → Transcript(s) → Protein isoform(s) → modified states, never Gene → Protein; definition (what it is) is separate from state (where, when, how much); a structure carries its source and method and predicted never equals experimental; a STRING association is not a proven interaction; knowledge is imported and versioned locally, never fetched per simulation tick | PROTEIN.md; `molecules/compiler.py` |
| D24 | 2026-09-10 | Reactome pathways execute as reachability graphs with knockouts because the export carries no rate laws; kinetics come from SBML where BioModels has the model | `genomeos pathway` |
| D25 | 2026-09-11 | Translation is verified against the curators chromosome by chromosome; a disagreement is a bug until shown otherwise (this found the synonym-match bug and the selenocysteine case) | `genomeos verify`; run by the proteome job |

## Cancer and therapeutics

| # | Date | Decision | Consequence |
|---|---|---|---|
| D26 | 2026-09-10 | cBioPortal's public REST API is used, not its MCP server (needs auth, not reproducible); Ensembl VEP REST annotates any chromosome; driver knowledge is distilled to one committed summary | `cancer/cbioportal.py`; `data/results/cancer_msk_impact_2017.json` |
| D27 | 2026-09-10 | The internal `Variant.pos` is 0-based; VEP, VCF and every report are 1-based at the boundary | caught by the demo disagreeing with the local classifier |
| D28 | 2026-09-11 | Three target spaces are searched in parallel (direct surface, peptide/HLA, indirect pathway-induced) and a nuclear mutation never gets "target = no"; target, binder, mechanism, cargo and effector are separate entities; the negative set is as important as the positive one; a conclusion is never read above its patient data level | THERAPEUTICS.md; `TherapeuticDesignDataset` v1.0 |
| D29 | 2026-09-11 | External accessibility is claimed only with positive evidence of an outward-facing part; a membrane annotation alone is "side not established" | the SRC case in PROGRESS.md |

## Product

| # | Date | Decision | Consequence |
|---|---|---|---|
| D30 | 2026-09-10 | One HTML page, inline SVG, no build step; each phase adds a view, never a separate tool | `web/static/index.html` and per-tab scripts |
| D31 | 2026-09-11 | Long jobs run once, through the registry, with the pid recorded on disk; a restarted server recognises a live job | `jobs.py` liveness check after five copies of the proteome job ran at once |
| D32 | 2026-09-11 | Governing documents are layered: README and ROADMAP (wide), ARCHITECTURE and DECISIONS (how and why), one design document per area, specs per language version, PROGRESS and LESSONS as records; docs/README.md is the index | docs/README.md |
| D33 | 2026-09-11 | The cycle is finish → `scripts/check.sh` (lint, format, tests, `bio test`) → commit to `dev` → push once per finished piece; a `pre-push` hook in the shared checkout refuses a red push; CI cancels a run superseded by a newer push; pull requests to `main` are opened by Albert only, never by a session or a script | CONTRIBUTING.md; `.github/workflows/ci.yml` concurrency; the CI failure notifications of 2026-09-11 came from intermediate pushes |
