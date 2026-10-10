# GenomeOS: scope, architecture, areas, plan

Reviewed and reorganised on 2026-09-17 against the code, the tests, the
result files and the running jobs (first written 2026-09-11). Reviewed again on
2026-10-10: the header, §1, §2, §2.1 and §8 were brought current, and the order of
work from that date is the execution plan in §5. This is the one document that says what
GenomeOS is, how it is built, what each area is for, what is done, what is
missing and what comes next. Evidence per finished task stays in
[PROGRESS.md](PROGRESS.md); the original phased plan is kept as history in
[ACTION-PLAN.md](ACTION-PLAN.md); lessons the code carries are in
[LESSONS.md](LESSONS.md).
The roadmap as it stood before this cut, every row and correction verbatim, is
[ROADMAP-HISTORY.md](ROADMAP-HISTORY.md), frozen at `c00250f` and read-only.

Measured state on the morning of 2026-10-10: **1,556 commits**, 63 `genomeos` commands plus the
`bio` toolchain, **5,784 tests** collected (5,571 passed in CI on `06d0761`, and every
daily CI run since has been green), 44 `.bio` programs at 227/227 `bio test` checks,
lint clean, every human chromosome fetched and analysed, the whole human proteome
compiled and verified against UniProt, a genome-wide knowledge graph, the whole
worm grown from one cell, a human body grown as populations, haematopoiesis as
mechanism, and the engine packaged and installable on its own. Milestones 0.9, 1.0,
1.1, 1.2 and 2.0 are reached; 1.3 (the 98%) is partly met (§6). No commit had landed
from 2026-10-03 until that morning, and the paid AstroREG-2 run 2 stopped at 363 of its
1,232 answers (§5 item 5). Work resumed the same day under the execution plan in §5: the
`reader` construct landed in the language, and the first two registered tests of what
the decoded genome predicts when run returned "no difference detected" (§5 items 1 and 3).
On 2026-09-17 the same line read 536 commits and 1,050 tests.

Both genome-wide sweeps are finished: **961,227 enhancer elements deleted** in
AlphaGenome across 24 chromosomes, and the **HPRC panel of 90 human genomes read on
all 24** (1,925,587 storage units, 187,966 executor-ready). What those two sweeps
established is smaller than what they were expected to: the node beats random
boundaries — by +2.9 points against a count-matched control and **+5.9 on 661 measured CRISPRi
pairs**, the size depending on how the control is matched (+1.2 to +6.6; re-audited 2026-09-27) — no unknown tier can be ordered against another once the
comparison is standardised, and **0.45% of the unknown space has ever been measured
by any assay this project holds**. The experiment that would change that is designed
and unrun: 284,001 oligos, controls matched in advance (§5 item 3b).

The four corrections of 2026-09-16 and 2026-09-17 are in the record rather than
edited out, and the arithmetic behind them now lives in `genomeos/compare.py` so the
next one is caught by a test instead of by a reader.

---

## 1. Main scope

The ambition, in Albert's words: decode the full genome on its blocks,
their function and their structure, and mimic its behaviour through
BioLang, BioIR, BioVM and the rest of the family. GenomeOS is therefore an
executable model of biology for geneticists: a compiler that turns a real
genome plus public scientific knowledge into an executable model, and a
virtual machine that runs that model through time, from one cell to an
organism, with evidence and uncertainty on every fact. Open source: Apache 2.0
for the BioLang engine and AGPL-3.0-or-later for the application
([LICENSING.md](../LICENSING.md), DECISIONS D39). The governing documents are indexed in [README.md](README.md) and
the decisions behind them in [DECISIONS.md](DECISIONS.md).

The organism is treated as a functional system and DNA as its seed. Four
things stay first-class: timers (construction and maintenance cycles),
blueprint (structure), the parts list (organs, tissues, cell types) and how
systems connect. `UNKNOWN` is a legal value; nothing is filled in.

Ambition: the software a geneticist reaches for first to work, research,
investigate and build. Every feature is judged by that bar: does it help
decode, model or build, and does it say how sure it is.

Non-goals, permanent: no wet-lab synthesis, no therapeutic sequences,
constructs, vectors, formulations, dosing or protocols, no clinical advice,
no claim that a therapy works. Storage principle: stream, distil, discard;
raw downloads never enter the repository, summaries do.

---

## 2. Architecture

Two products share one repository today and separate later without a
rewrite: **BioLang** (the language, IR, VM and standard library, with the
`bio` toolchain) and **GenomeOS** (the application: genome decoding,
molecules, twins, organisms, cancer, therapeutics, web UI). The split is a
move of packages, not a redesign, as long as the boundaries below hold.

```
BioLang (engine, standalone later)         GenomeOS (application on the engine)
─────────────────────────────────          ────────────────────────────────────
genomeos/lang      parser → BioIR          genomeos/genome      sequence, annotation, blocks, unknown,
genomeos/ir        BioIR types, JSON                            repeats, regulation, domains, signals,
genomeos/runtime   engines (BioVM):                             anatomy, fetch, lookup, telomere
                   central dogma, GRN,     genomeos/molecules   protein compiler, RNA, proteome, verify,
                   Boolean, SBML, cell,                         graph, Reactome pathways
                   spatial, segmentation,  genomeos/flow        DNA → RNA → protein trace
                   gastrulation, body,     genomeos/twin        clocks, twin fork / run / diff
                   debugger, uncertainty,  genomeos/organism    reference lineages, diff, human, experiments
                   compose (bigraph)       genomeos/cancer      cBioPortal, tumour-only pipeline, compare
genomeos/std       bio.std.* prelude       genomeos/therapeutics targets, mechanisms, design dataset
genomeos/bio.py    bio check|compile|run|  genomeos/knowledge   GO, Reactome, CellPhoneDB, Cell Ontology
                   test|repl               genomeos/lib         library catalogue and membership
genomeos/coords    Locus (shared)          genomeos/predict     AlphaGenome adapter (optional)
                                           genomeos/design.py   minimal-cell design
                                           genomeos/forge       BioForge search
                                           genomeos/web         server + one-page UI
                                           genomeos/jobs        background jobs registry
                                           genomeos/storage     stream / distil / clean
                                           genomeos/results     result summaries (data/results)
                                           genomeos/attribution the 98%: tiers, attributions, compiled programs
                                           genomeos/benchmark   known loci, CRISPRi, held-out scores
                                           genomeos/provenance  what each hand-written number rests on
                                           genomeos/cli.py      63 commands
```

Boundary rules (from ARCHITECTURE.md §10) and where they are broken today:

| Rule | State |
|---|---|
| BioIR JSON is the contract; the language layer reads no GenomeOS data paths | holds |
| Engines take a `Module` and return a trajectory or a typed result | holds |
| Standard library written in BioLang, not Python | holds (ageing, cell_types, development) |
| One thin entry point (`bio`) | holds (`genomeos/bio.py`, 481 lines on 2026-10-10) |
| No dependencies in the language core | holds |
| Engines import nothing from genome, knowledge or web | **holds since 2026-09-11** (b162919 and b90e741): a grep over `lang`, `ir`, `runtime`, `coords.py`, `version.py` and `bio.py` on dev finds no import of the application. `bio` runs on an engine-side `lang/tools.py` instead of the CLI, and the version string moved to `genomeos/version.py` so the Apache half never reads the AGPL package root |
| `lang` and `ir` depend only on themselves | `Locus` lives in `genomeos/coords.py`; move it under `ir` when the packages split |

### 2.1 Where each component stands (2026-09-11; the BioLang, BioIR, `bio`, BioForge, cancer and web rows brought current 2026-10-10)

"Functional" means a geneticist can use it today on real data with evidence
attached. "Missing for complete" is what separates it from the ambition.

| Component | Version | Functional today | Missing for complete |
|---|---|---|---|
| **BioLang** (language) | v0.4 + v0.5 `reader` | 27 block kinds plus the `module` and `import` directives (v0.3's 18 counted the directives; v0.5's `reader` landed 2026-10-10, `5682e6c`): v0.3's `gene` to `experiment`, then `field`, `design`, `competence`, `commitment`, `order` and the v0.4 structure and economy blocks (`compartment`, `transport`, `pool`, `allocation`, `regime`); evidence and confidence on every block; `bio.std` prelude; programs test themselves; the grammar is generated from the parser (BIOLANG-GRAMMAR.md, a test fails when it is stale) and unknown block keys are refused | the `reader` construct, so a rule's context comes from chromatin (§5, M2); `population` and `writer` constructs; a registry for shared libraries; std modules for signalling and human development; v0.4 stage 3's gate as §6 words it (built 2026-10-10 with 10 of 10 registered clauses met; three asks not askable) |
| **BioIR** (representation) | 0.4 | 26 types with JSON round-trip, listed in BIOLANG-GRAMMAR.md (v0.3's 20 plus Compartment, Transport, Regime, Field, Design and, since 2026-10-10, Reader; the generated list omits six dataclasses, §5 item 13); `UNKNOWN` as a value and `UNSTATED` for a confidence nobody stated | ProteinState and Population as types; `Locus` moved under `ir` |
| **BioVM** (engines) | 0.3 | twelve engines: central dogma (exact on mtDNA and verified against UniProt), gene network (Hill ODE), Boolean (attractors), SBML (MathML + RK4), cell ageing, 2-D spatial fields, segmentation clock, gastrulation, Body (discrete events, one cell to the whole worm and a human as populations), debugger, uncertainty; composition through process-bigraph | one composite that runs Body + network + SBML on a shared clock; division and movement in space; external engines (libRoadRunner, MaBoSS, CompuCell3D) behind the same interfaces |
| **`bio` toolchain** | 0.1 | `bio check | compile | run | test | repl`; `bio test` runs 44 programs at 227/227 checks in CI; packaged with the engine as the separately installable `biolang` package with its own test suite (milestone 2.0, 2026-09-27) | `bio fmt`; language server; publication to a package index (the owner's) |
| **BioLib** (libraries) | 45 libraries + the proteome | data-computed membership from GO and Reactome (95.5% agreement with the curated catalogue); five layers (core, blueprint, timer, systems, parts); ligand-receptor protocol; haematopoiesis as a runnable mechanism module; **the whole human proteome packaged** (19,283 proteins in 2.3 MB, `genomeos protein X --lib`, answers offline) | `cancer.*` layer; more human development modules as runnable BioLang; "phenotypes it must reproduce" lists per library |
| **BioTwin** (individual) | 0.9 | any GRCh38 genome imported locally and used by carrier lookup, the gene report, the twin and the gene-by-gene walk; HG002 across all 22 autosomes with measured telomere and epigenetic age, the telomere also read from the 601 GB BAM by range without a download; fork, run, diff; predicted effect of a person's own regulatory variants; ClinVar carrier screen, genome-wide truncating-variant scan, coding inventory by consequence with AlphaMissense scores kept off the repository, a one-page dossier; the GIAB trio read for inheritance, 96.4% of HG002's calls inherited | true de novo variants out of the trio (phasing waits for a phased import); a population distribution for the polygenic scores; the effect per regulatory variant (area I's model) |
| **BioForge** (design) | search, experiments as input | design search under constraints with predicted labels; experiments as input (milestone 1.0); a `design` block over organism programs; a design answer carries its effect and score and no stated probability, because no calibration record exists for any outcome (R4a, `16a0434`); minimal-cell design estimate | a population of design answers with published outcomes to calibrate against; a constraint language |
| **Genome decoding** (GenomeOS) | 0.9 reached | every chromosome fetched and analysed (2026-09-11): inventory, UNKNOWN classified genome-wide with curated repeats (98.5%), domains on 24, curated repeats and ENCODE elements on 25, reader on 25 with eleven cell types; the node model tested genome-wide (90.2% of enhancers act inside their node) and on two mouse chromosomes; the segment parser on chr21 at 92.7% gene precision once RNA over the exons counts as evidence | the parser on every chromosome with measured RNA; enhancer targets and boundaries from Hi-C |
| **Molecules** (GenomeOS) | 0.9 reached | the whole human proteome compiled from seven public databases (19,478 coding genes, 25 chromosomes) and packaged as a 2.3 MB offline library, translation verified on 19,249 of them, genome-wide knowledge graph, RNA layer with GTEx, pathways as reachability and as kinetics where a curated ODE model exists; 96,362 modifiable sites with their 352 writers as `modifies` edges in the graph; the dominant isoform per tissue from GTEx (APP695 in the cerebellum) | measured modification state; isoform expression per cell type or stage; an engine that handles 100-species models |
| **The 98%** (GenomeOS) | genome budgeted | every UNKNOWN block of the 24 chromosomes tiered with Zoonomia constraint read per base (1,009 Mb): structural 23%, fossil 33%, regulatory 34%, neutral 7%, constrained-unknown 3.2% (1,098 blocks, 32 Mb, 1% of the genome); `genomeos budget`, the Progress tab card; the attributions compiled to a BioLang program **on all 24 chromosomes** (2026-09-17, 940,803 stated facts, 93.6% predicted, none experimental); scored against VISTA, eQTL, lentiMPRA, GWAS and ClinVar; the organiser reads copies first and leaves a real unknown of 882 blocks, 69 of them constrained on both axes; **0.45% of the unknown space has ever been measured by any assay held here, and the experiment that would change it is designed and unrun (284,001 oligos)** | measurement over the real unknown, which is the binding constraint; closure at cell and organism level |
| **Cancer and therapeutics** (GenomeOS) | 1.2 reached | tumour-only pipeline, cohort expression, altered-protein reconstruction, CNA and SV tables, 15 mechanisms, design dataset with negative set, the retrospective benchmark (milestone 1.2, 2026-09-27: ten tumour cases, the tenth hand-written from cohort medians; target recovered, route and top mechanism scored apart) | surface scoring that reads the alteration rather than the gene's annotation; a scored structural-variant case; a matched normal instead of the population consensus |
| **Web UI and CLI** (GenomeOS) | 19 views, 63 commands | every layer visible with evidence pills; jobs with progress; gene dossier; the Cells and Cell views; the Progress tab shows the work board, the ordered next steps of §5 and the task log; the Evidence explorer reads the hand-written programs and the 24 compiled chromosomes; daily CI on `dev` | a release tag (the owner's); the AlphaGenome live test needs a repository secret (the owner's) |

Data architecture: `data/reference` and `data/knowledge` are caches (git
ignored, rebuilt by `genomeos data fetch` and the compilers); `data/results`
holds the committed summaries, one per real-data run; `data/organisms` and
`data/demo` hold BioLang programs; `data/jobs` holds job metadata and logs.

---

## 3. Areas

Each area lists its goal, the code that implements it, the design documents, the results it has
committed, what it requires, what is built (one dated bullet per finding, with its figures and the
commit that established it), what is missing, and its open next steps in order. Brought current on
2026-10-10 against the code and the §5 rows: a step that was finished has left the Next list and
appears among the built bullets. "Owner" names the live lane that holds the area's files, if any;
lanes are assigned in §5. Each area ends with the lines of ROADMAP-HISTORY.md that hold its full
record, every row and correction as it was written.

### A. BioLang and BioVM (the engine)

- **Goal.** A language a bio-engineer writes biology in, an IR that carries
  evidence, confidence and UNKNOWN, engines that run it, and a toolchain that
  stands alone (`bio run | check | compile | test | repl`).
- **Code.** `lang/parser.py`, `lang/grammar.py`, `lang/tools.py`, `ir/model.py`, `runtime/*`
  (`grn.py`, `economy.py`, `located.py`, `division.py`, `methylation.py`, `body.py`), `std/*.bio`,
  `bio.py`; the engine package is built by `scripts/package_engine.py`. Versions: BioLang v0.4 with
  the v0.5 `reader` construct, BioIR 0.4 (`bioir_version`), the `biolang` package 0.1.0.
- **Design.** BIOLANG-v0.1/0.2/0.3.md, BIOLANG-v0.4-ECONOMY.md, BIOLANG-v0.5-READER.md,
  BIOLANG-GRAMMAR.md (generated from the parser), BIOIR-v0.1.md, BIO-TOOLCHAIN.md, ARCHITECTURE.md
  §10, SEQUENCE-GRAMMAR.md.
- **Data.** `data/demo/*.bio` and `genomeos/std/*.bio`; `bio test data/demo genomeos/std
  data/organisms` runs 45 programs at 246/246 checks since the reader construct landed on 2026-10-10 (44 at
  227/227 before it).
  Results: `located_mitochondrion`, `plasticity_gate`, `burden_gate`, `abundance_gate`,
  `division_gate`, `engine_package`.
- **Requirements.** Every construct has evidence and confidence; programs
  test themselves (`# test:` and `assert:`); a `.bio` file runs without the
  rest of GenomeOS; the same program runs on any engine that fits it. State on 2026-10-10: the
  language, the IR, the runtimes and the toolchain are built and packaged on their own (milestone
  2.0); v0.4 stages 1, 2 and 4 are built and stage 3 is not; the `reader` construct's engine half
  landed on 2026-10-10.
- **BioLang v0.4 stage 1, structure (2026-09-14, `b40959f`).** `compartment`, `signals`, `transport`
  and `regime`, run by `runtime/located.py`: 13 of 13 chrM proteins stay inside the mitochondrion, all
  1,123 nuclear-encoded MitoCarta3.0 proteins reach it only through a transport, and the rho0
  phenotype comes out right (King & Attardi 1989).
- **Fate: commitment and competence, gated by a published series (2026-09-14, `126e65b`).**
  `plasticity.bio` runs Fukushige & Krause 2005 and Yuzyuk et al. 2009 under `bio test` (forced at 60
  min, 610 of 610 cells become muscle); removing either construct breaks an arm.
- **`share:` for a partition (2026-09-17, `e362ef2`).** The sixteen germ-layer splits of
  `tissues.bio` state the published share of the layer, and its first use found an error: Glia,
  LungEpithelium and Myocyte, each the last split of its layer, had held no cells, and **Glia is 24.5%
  of the ectoderm**. Two committed numbers were corrected with their reason.
- **Stage 2 of v0.4, the economy: built, and gate (b) failed (2026-09-21, `d334b48`, `eaa603d`, `6929a40`, `26b76bc`, `e6d4e4f`, `b57de43`, `a9f566f`).** `pool`, `cost` and `allocation` parse and
  refuse; `runtime/economy.py` charges a cost and scales each entity by its tightest pool. **Gate (a),
  burden: passed on three clauses of four**; the order-of-magnitude clause against Ceroni 2015 and
  Frei 2020 is registered **not askable**. **Gate (b), absolute abundance: failed**: on 15,159 genes,
  log10 MAE **0.952005** for the baseline and **0.952005** for the allocation arm, improvement
  **0.000000000**, because `_share` gives every demander of a pool the same factor and a fitted
  baseline absorbs it. **The pool layer stays optional.** `competitive` is as gene-blind as
  `proportional`, which closed §10 decision 5 with a negative.
- **The engine packaged on its own, milestone 2.0 (2026-09-27, `5720c7c`, `14f9ada`).** `biolang`
  0.1.0 runs in an interpreter where `import genomeos` fails, with its grammar, specs, `bio repl` and
  an `import protein:` that stops with a named error. The committed `engine_package.json` (2026-09-29)
  reads 12 of 12 checks, the engine's own pytest suite among them. Not done and not a clause:
  publication to PyPI, which is the owner's.
- **Unknown block keys refused (2026-09-27, `61c5838` census, `b41ddd6` parser).** A key a block does
  not read is an error naming the key, the block, the closest match and the accepted keys.
- **v0.4 stage 4, partitioning division (2026-09-28, `4ff3cd4` registered, `73af2d5` built).**
  `partition: duplicate | contents | binomial` on an `event`; `runtime/division.py` conserves every
  species exactly. **8 of 8 registered clauses met**, and the spec's own dilution gate for the stage
  is refuted, as derived before the instrument existed: partitioning is gene-blind and contributes 0.0
  of the 9.392136 decades between a 240 h and a 6 h protein.
- **The `reader` construct, the engine half (2026-10-10, `5682e6c`).** `reader <CellType> { open;
  poised; closed; unknown }` and `through: E` on a rule: only an open element runs the rule, and every
  other state, or no reading at all, holds it at strength 0 rather than dropping it. It says
  accessibility, not function, and it does not read `context_evidence`; the data half is area B's.
- **Missing.** A per-demander saturable allocation policy: `runtime/economy.py` offers
  `proportional`, `priority` and `competitive`, and the two saturating ones are gene-blind, so gate
  (b) cannot be answered. v0.4 stage 3 (core metabolism, mitochondrial copies, heteroplasmy) is not
  built, and until it supplies ATP over time a declared division cost refuses every division. Imports
  from a remote registry (the resolver hook exists). The control layer of v0.4 (`homeostat` with set
  points and gains, `role`) is specified and not built. Gate (a)'s magnitude clause and decision 2's
  volume arm are named and unrun.
- **Next.** 1. A per-demander saturable allocation policy, registered against gate (b) before it is
  run. 2. v0.4 stage 3: core metabolism and mitochondrial copies, with its falsifying measurement:
  **done 2026-10-10** (lane-stage3; `e09edf7` registered as spec §9.5, `ac89b41`, `fa5a4f1`, `4668ee2`,
  `f99d634` measured as §9.6): `heteroplasmy:` on a chrM compartment, `bio.std.core_metabolism`,
  `runtime/metabolism.py`, a per-pool capacity in `Economy`; **10 of 10 registered clauses met, and the
  gate as §6 words it is not passed**: three of its asks are not askable here, the budget hierarchy is
  reproduced as an input rather than predicted, and the Pasteur shift holds by construction; decision 4
  (whose cost accounting) stays open. 3. Imports from a remote registry. 4. Splitting the repository, which is a release decision rather
  than an engineering one now that the package builds. 5. The control layer of v0.4 (`homeostat`,
  `role`), specified with its falsifying measurements, which waits for Albert's steer.
- **Owner.** Unowned: lane-readerlang built the reader's engine half on 2026-10-10 and finished;
  lanes are assigned in §5.

History: ROADMAP-HISTORY.md §3 A, lines 160–332.

### B. Genome decoding (reverse-engineering the sequence)

- **Goal.** Every base of a real genome accounted for: gene, transcript,
  exon, regulatory element, repeat, domain, or an UNKNOWN with a class and
  a confidence; the delimiters the cell reads learned from the sequence.
- **Code.** `genome/*` (annotation, blocks, unknown, repeats, regulation, domains, signals, anatomy,
  fetch, lookup, segments, reader, epigenome, hic, mouse), `genome/reader.py` above all;
  `predict/rna_tracks.py`, `predict/contact_maps.py`, `predict/enhancer_target.py`; scripts
  `reader_genome_wide.py`, `node_containment_audit.py`, `crispri_direction.py`.
- **Design.** GENOME-AS-CODE.md, GENOME-ANATOMY.md, SEQUENCE-GRAMMAR.md, UNKNOWN.md,
  NODES-READER-WRITER.md, BIOLANG-v0.5-READER.md.
- **Data.** `anatomy_hg38_by_chromosome`, `unknown_chr*` (25), `rmsk_chr*` (25), `ccres_chr*` (25),
  `domains_chr*` (24), `segments_chr21`, `mouse_mm10_chr*` (2), `epigenome_genome_wide`,
  `domains_oriented_comparison`, `reader_genome_wide`, `reader_normalised`,
  `node_containment_measured`, `crispri_direction_both`, `silenceragree`.
- **Requirements.** Any chromosome on demand (`genomeos data fetch --chrom`);
  classification order fixed (curated repeats before ORFs); every class
  names its evidence layer; the block map shows it. State on 2026-10-10: every chromosome is
  fetched, classified and cut into CTCF nodes; the reader runs on thirteen biosamples; the node is
  held against measured CRISPRi pairs; enhancer-to-gene is read by a model and, for 661 pairs, by
  measurement.
- **Enhancer to gene, genome-wide (2026-09-11).** Of 4,800 distal enhancers deleted in AlphaGenome,
  2,994 move a gene and 2,291 name a coding gene, 90.2% of those inside the element's CTCF node;
  qualified on 2026-09-14 to 81.7% against 79.1% at random over 113,399 elements. The node figure
  now stands on the audit below. 1,157 of the calls are silencer-like: expression rises on deletion,
  a direction and not a function.
- **The segment parser, three runs and one conclusion (2026-09-12, `abfbd4c`).** On chr21, predicted
  and measured together, 69 of 72 candidates are genes (95.8% precision at 54.3% sensitivity); on
  chr22 the same configuration reads 95.1% at 80.1%. The grammar finds the genes, splice-site
  prediction draws their exons, and transcription evidence says which are real; the series is
  closed, every run kept side by side.
- **The second mammal (2026-09-12, `48c48bc`).** 379 tested mouse nodes sit at 93 to 94% in the same
  human neighbourhood; Compara and MGI agree on 97.3% (chr19) and 96.3% (chr11), so the node finding
  does not depend on the orthology source.
- **Measured Hi-C boundaries for the node edges (2026-09-12).** On chr21 only, GM12878 supports 58%
  of the 227 CTCF-only edges against 42% at random; qualified on 2026-09-14, because its dense calls
  put a random position within 20 kb of a boundary 54% of the time. Still a proxy; edges stay at 0.4.
- **The epigenome layer (2026-09-14, `46cb5d2`).** Five histone marks in eleven biosamples and WGBS in
  eight; the marks move direction AUC 0.53 to 0.62 against an AlphaGenome outcome, so nothing here is
  measured. The orientation-aware caller is kept as an option beside the CTCF-only default.
- **The stricter site call, closed as a negative (2026-09-21, lane-sites).**
  Item 1 of the Next list below is answered, both halves. Registered in
  `ad07dfb` before scoring and run in `56e2c50`: best-hit strand plus a 0.95
  motif threshold, judged on the same four measurements as the 2026-09-14
  caller, with the split decision fixed in advance. **Hi-C passes** — enrichment
  1.75, 2.36, 2.27, 1.60 in H1, K562, HepG2 and IMR-90 against 1.74, 1.96, 1.99,
  1.45 — and it is a precision bought by keeping 2,171 edges where the old
  caller keeps 17,971, so reach falls from 41–48% of measured boundaries to
  5–7%. **Node content fails**: excess over its own random control is **-0.3
  points against the default's +2.9**, the raw 94.5% being pure resolution — and
  that raw number being the highest in the table is exactly what fixing the
  judged statistic in advance was for. **Mouse synteny fails** on chr11, 85.3%
  against 92.4%, below even the orientation caller's 88.0%; chr19 tests 73 nodes,
  under the registered minimum of 80, and is not judged. **HOXD fails**, nearest
  edge 1,097 kb from the published interval against the orientation caller's
  16 kb. One pass and three failures is the split case the registration decided
  in advance: **the default stays CTCF-only**, the stricter call is a named
  option, and the item's second half is answered **yes** — Hi-C questions and
  enhancer-to-gene questions need two node sets.
  **The attribution matters more than the verdict, and it is a clean negative on
  the half that was blamed.** Best-hit strand does recover the 11,696 elements
  (9.5%) dropped for a weak opposite-strand hit, which is the mechanism the
  2026-09-14 note guessed — and it moves **all four measurements inside their
  registered noise floors**. Those elements are not where the boundaries are.
  The whole effect is the threshold, and what the threshold does is delete: 0.95
  keeps **8,616 of the 123,166** elements that have a site at 0.85, and 8,613 of
  the survivors are already single-strand, which makes the strand rule moot.
  A 0.90 call added **after** the registered run and labelled post-hoc,
  registered as unable to decide anything, is the best Hi-C caller of the seven
  and still fails node content, chr11 synteny and HOXD — so 0.95 was not simply
  the wrong number. **No further site-call variant is worth running against these
  four measurements**; what the node comparison needs next is a measured
  enhancer-gene set in place of the model's own reading, which is item 4 of the
  list below. Node content was re-measured on 440,377 archive elements rather
  than 2026-09-14's 113,399, every caller re-run in the same job, and the
  default's excess survives the 3.9× growth (+2.6 → +2.9).
  NODES-READER-WRITER.md, "A stricter site call, judged on all four measurements
  at once".

  *(2026-10-10: the bullet above is kept verbatim because the Progress tab quotes it; the "Next
  list below" it names is the one in ROADMAP-HISTORY.md, and its item 4, the measured set, is the
  node containment bullet below.)*
- **Biosamples twelve and thirteen, and the per-biosample family banded against depth (2026-09-22, `ca8aa0f`, `1680ebc`, `41bc7d6`, `86d668f`).** Testis and ovary joined the eleven. `nodes_open`
  turned out to read 9,988 to 10,016 of 20,002 in all thirteen and was withdrawn as a
  between-biosample reading; the normalised readings that replaced it reached milestone 1.1 on
  2026-09-27 (`63929e1`, `b13143e`).
- **Direction against measured perturbations (2026-09-22, `5c842ca`, `d492834`).** The downward arm
  passes, 41 of 44 (0.9318, 95% [0.8177, 0.9765]); **the upward arm fails, 15 of 36 (0.4167
  [0.2714, 0.5780])**, and the combined balanced accuracy is 0.6221 [0.5354, 0.7088]: undecidable.
  The `represses` half of `action` is measured and unsupported.
- **The node containment claim, audited and held against measurement (2026-09-27, `05b71ae`, `2b77663`).** Against a control matched on boundary count only, **+2.90, 95% CI +2.03 to +3.81**,
  and +1.21 to +6.58 across four defensible baselines. On the CRISPRi-measured regulated pairs (n =
  594 of 661 scored), **+5.89, 95% CI +3.18 to +8.46**, with 89.1% of the pairs from K562. The
  supported wording: nodes hold together enhancer-gene links at distances where randomly placed
  boundaries would separate them.
- **`silencer_like` checked against assay-validated silencers (2026-10-03, `c3ebf70` registered, `a5a43b3` result).** 420 of 1,181 (35.563%) against a matched control's 37.707%: NO_SIGNAL. The
  key is reported as a direction; it licenses nothing about function.
- **Missing.** No curated silencer source in the pipeline: the derived `silencer_like` is the only
  call, and the reachable curated source covers half the scored cell lines (0 validated silencers of
  23,894 in GM12878). Measured Hi-C node edges exist on chr21 only. The read-by-marks call is shown
  nowhere. A neural or embryonic biosample beyond SK-N-SH and H1 is not added, and §5 item 10's stop
  list says no more biosamples until a reading needs one. The `reader` construct has its engine half
  and no generated reader blocks yet.
- **Next.** 1. The `reader` construct's data half (§5, M2b): reader blocks generated per cell
  type from the reader and epigenome results, so a compiled rule's context is gated by chromatin.
  The engine half landed on 2026-10-10 (`5682e6c`); its three open language decisions (UNKNOWN held
  or refused, poised, the version number) are Albert's, in BIOLANG-v0.5-READER.md. 2. A curated
  silencer source joined in the pipeline beside the derived call. 3. Measured Hi-C node edges beyond
  chr21. 4. The read-by-marks call shown in the Blocks lane.
- **Owner.** lane-readercompile works on Next 1 from 2026-10-10 (§5, M2b: the compiler emits
  reader gating), holding `attribution/compile.py` and `attribution/readergate.py`; the area's own
  files are unowned since 2026-10-03, and lanes are assigned in §5.

History: ROADMAP-HISTORY.md §3 B, lines 333–955.

### C. Molecules (RNA, proteins, pathways, the knowledge graph)

- **Goal.** From gene to transcript to protein isoform to modified state,
  compiled from public sources into one definition per protein, verified against the curators, and
  executable as pathways.
- **Code.** `molecules/compiler.py`, `molecules/proteome.py`, `molecules/verify.py`,
  `molecules/rna.py`, `molecules/graph.py`, `molecules/reactome.py`, `molecules/uniprot.py`,
  `molecules/alphafold.py`, `molecules/ptm.py`, `molecules/biomodels.py`; `flow/trace.py`;
  `runtime/sbml.py` (the in-house kinetic engine); `lib/biolang.py` and the packaged
  `lib/data/proteome.json.gz`; scripts `phosphosite_observation.py`, `cptac_phospho.py`.
- **Design.** PROTEIN.md, FLOW.md.
- **Data.** `proteome_chr*` (25 of 25) plus `proteome_genome_wide`; `translation_vs_uniprot_chr*`
  (25 of 25), with no genome-wide file of its own; `graph_chr21`, `graph_genome`, `ptm_genome_wide`,
  `phosphosite_observation`, `phosphosite_tumour_abundance`; compiled definitions in
  `data/knowledge/proteins` (19,453 files).
- **Requirements.** UniProt accession is the id; Gene → Transcript →
  Protein isoform, never Gene → Protein; predicted never equals
  experimental; STRING association is not interaction; never call an API
  per simulation tick. State on 2026-10-10: the whole proteome is compiled, verified and packaged;
  modification sites carry observation but not occupancy; kinetics run where BioModels has a curated
  model and one recorded model does not reproduce.
- **The whole human proteome (2026-09-11, `3a64298`).** 19,478 coding genes over 25 chromosomes:
  sequence 99.1%, domains 99.0%, function 86.3%, interactions 81.5%, pathways 58.2%, experimental
  structure 45.2%, predicted structure 98.2%, disease 25.5%. Translation verified against UniProt on
  19,249 genes: 90.9% identical to the canonical entry and 97.8% exact for some isoform. The graph
  genome-wide reads 42,045 nodes and 330,018 edges, largest component 15,165 of 19,283 compiled
  proteins. Packaged as a library, `genomeos protein X --lib` answers with no network at all.
- **Kinetics where a curated model exists, and a recorded failure (2026-09-11).** RAF1 out of
  Hornberg2005 erases the downstream transients (peaks 0.201 and 0.562 to zero); MEK out of
  Kholodenko2000 stops the oscillation. Schoeberl2002 (100 species, 125 reactions) does not reproduce
  on the in-house engine: the cascade species stay flat, and the result is committed beside the two
  that work.
- **Post-translational state and isoform expression (2026-09-11, `96f624b`).** 96,362 sites on
  13,083 proteins with 352 named writers; a site that can be modified is not a measurement that it is.
  Isoform expression reads GTEx per tissue.
- **Modification observed (2026-09-27, `58f84fb` registered, `869158d` result).** 29,279 of the
  41,661 curated phosphosites (70.3%) carry the cell types or tissues they were observed in, from
  Ochoa et al. 2020's reanalysis of 6,801 PRIDE runs at 1% site FDR; curation and reference share
  source studies, so the 70% partly agrees with itself. "Observed modified in N experiments" is not
  deliverable from either source.
- **Tumour abundance joined (2026-09-28, `cae8d7d` registered, `5938b5d` result).** 19,104 curated
  phospho sites on 5,105 proteins (45.9%) carry a relative abundance in lung adenocarcinoma,
  glioblastoma and paediatric brain tumours, at 0.02% residue mismatch: never occupancy, never per
  patient. Coverage missed high on both registered ranges, reported as a miss, and the seven
  gene-keyed profiles stay uncommitted at 14.3% residue mismatch until their isoforms are known.
- **Missing.** The person's own expression: GM24385 has no ENCODE experiment and the GIAB direct-RNA
  runs are raw reads, so `individual protein` states GTEx's median and says why. Measured
  modification occupancy stays blocked. The libRoadRunner adapter does not exist
  (`runtime/sbml.py` is the stand-in), and Schoeberl2002 still fails *(2026-10-10: on the stand-in only; it passes 33 of 33 on the libRoadRunner adapter, `bdb30d3`, and the stand-in's failure was a divergence that now raises)*. Isoform expression is GTEx's
  adult tissues, not a cell type or a stage.
- **Next.** 1. The libRoadRunner adapter, scored first on Schoeberl2002 (§5, M6): **done 2026-10-10** (`1dcf998`, `bdb30d3`; Schoeberl2002 33 of 33 on libRoadRunner, 0 of 33 on the stand-in). 2. The
  seven gene-keyed CPTAC profiles, once their isoforms are known. 3. The person's own expression in
  place of GTEx's median, which waits for a processed track of GM24385. 4. Measured modification
  occupancy, which waits for a source that publishes it.
- **Owner.** Unowned since 2026-10-03; lanes are assigned in §5.

History: ROADMAP-HISTORY.md §3 C, lines 956–1122.

### D. The individual (BioTwin)

- **Goal.** A person's genome plus measured state plus environment, forked, run and compared.
- **Code.** `twin/twin.py`, `twin/clocks.py`, `twin/calibrate.py`, `twin/methylation.py`,
  `genome/telomere.py`, `genome/variants.py`, `genome/individuals.py`, `genome/missense.py`,
  `genome/polygenic.py`, `genome/diff.py`, `genome/regdiff.py`, `genome/clinvar.py`,
  `genome/bam_range.py`, `genome/bgzf.py`, `runtime/variant_effect.py`; scripts
  `trio_q100_phase.py`, `trio_q100_population.py`, `trio_vaf.py`.
- **Design.** ACTION-PLAN.md Phase 3, STORAGE.md, DATA.md, LESSONS.md (timers).
- **Data.** `hg002_chr21`, `hg002_telomere_stream`, `hg002_methylation_stream`, `clock_GSE41169`,
  `clinvar_chr21_agreement`, `calibration_hematopoietic_stem_attrition`, `telomere_range_HG002`,
  `trio_q100_phase`, `trio_q100_population`, `trio_vaf`. Per-person output stays under the git-ignored
  `data/individuals` and `data/twins`.
- **Requirements.** HG002 is the test human; measured state is labelled
  with its caveats (cell-line signature); every inferred parameter carries
  its confidence. State on 2026-10-10: HG002's genome is imported and read for what a person learns
  from it (screen, knockouts, coding inventory, missense effect, polygenic scores, telomere, trio);
  `genomeos/twin` itself is untouched since 2026-09-10 (`33c3629`), and every later result is in
  `genomeos/genome/`.
- **What a person learns from their own file (2026-09-11).** ClinVar streamed once (4.47 M rows,
  346,571 pathogenic and likely pathogenic kept locally); HG002's four million variants give two hits
  in seven seconds, the F11 carrier allele among them, and 88 knockouts in 84 genes. Research
  annotation, not a clinical report, and the page says so.
- **The coding inventory and the missense effect (2026-09-11 to 2026-09-12).** HG002 on 22
  autosomes: 21,019 coding SNVs, 9,830 missense, 61 nonsense; 5,438 genes carry a protein-changing
  variant. AlphaMissense, streamed once (71.7 million rows), scores 8,933 of the 9,830: 188 likely
  pathogenic, 254 ambiguous, 8,659 likely benign, confidence capped at 0.7; its licence keeps the
  scores out of every committed result. HG002 matches the reference in 5 of the 30 hg38 truncations
  on file.
- **Telomere from streamed reads (2026-09-12, GC index 2026-09-13, `9f55511`).** On GIAB's 300x BAM
  (601 GB), 218 MB in 48 range requests give about 2,676 bp per end at confidence 0.3, a number to
  compare between people read the same way, not to quote. A GC-corrected dimensionless index,
  `telomeric_per_gc_comparable_read` 0.001915, is reported beside it and deliberately not converted
  to base pairs.
- **The trio, normalised (2026-09-12).** With `normalise_variant`, de novo candidates fall from
  12,296 to 1,430 (1,173 SNVs, 257 indels) and Mendelian errors from 2,930 to 7, inherited 96.4% to
  96.9%.
- **Polygenic scores (2026-09-12).** Four PGS Catalog scores at 90 to 99% coverage, each a raw score
  with coverage, not a percentile, because the catalog publishes no population distribution; the
  child sits between the parents on three of four, the only check available without a reference
  population.
- **Phasing the 1,430 on the Q100 assembly (2026-09-27, `457702c` registered, `4d63910` result).**
  Of 1,411 assessable candidates, 1,382 sit on one parental chromosome, 703 paternal to 679 maternal
  (share 0.509). **P1–P3 failed as registered** (they expected 40–150 germline de novos at a ~0.8
  paternal share); P4 holds. 3.0% of the candidate SNVs are seen in 88 HPRC haplotypes against 92.8%
  of inherited calls: the candidates are real and new, most likely post-zygotic or cell-line
  mutations, which is an interpretation, not a measurement.
- **Read allele fractions (2026-09-28, `bf789a3` registered, `0125926` result).** The 1,375
  parent-assigned candidates sit at a median read fraction of 0.442 against 0.495 for 1.28 million
  inherited heterozygous SNVs: one shifted group, read as one clone's mutations carried in about 88%
  of cells, not germline. HG002's DNA is from a cell line.
- **Missing.** A population distribution for the polygenic scores. The trace follows a gene's
  canonical transcript, not the tissue's (`flow/trace.py`). A per-candidate germline call, and
  telling early-embryonic from line-founder mutations, both need DNA that never went through culture.
  The note `trio()` writes into every trio result (`genome/individuals.py`, around line 1111) still
  says a count above 60 to 100 is "representation differences and residual no-calls", which DATA.md
  records as wrong after normalisation.
- **Next.** 1. A population distribution for the polygenic scores (§5, M6). 2. Correct the
  `trio()` note to what the phasing and the read fractions found. 3. The tissue's own transcript in
  the trace. 4. A per-candidate germline call, which waits for DNA that never went through culture.
- **Owner.** Unowned since 2026-10-03; lanes are assigned in §5.

History: ROADMAP-HISTORY.md §3 D, lines 1123–1344.

### E. From one cell to an organism

- **Goal.** Grow an organism from one cell with mechanism where it is known
  and the observed program where it is not, timed by the genome's own
  timers, and say at each level how much is mechanism.
- **Code.** `runtime/body.py`, `organism/*` (`celegans.py`, `fate_rules.py`, `commitment.py`,
  `lateral.py`, `terminality.py`, `haematopoiesis.py`, `human.py`, `digital_development.py`,
  `packer.py`), `runtime/cell.py`, `runtime/grn.py`, `runtime/boolean.py`, `runtime/sbml.py`,
  `runtime/spatial.py`, `runtime/segmentation.py`, `runtime/gastrulation.py`, `runtime/compose.py`;
  scripts `celegans_fate_reads.py`, `celegans_terminality.py`, `mesoderm_diagnosis.py`,
  `sox17_repressor_search.py`.
- **Design.** ORGANISM-FROM-ONE-CELL.md, AREA-E-NEXT.md, BIOLANG-v0.3.md.
- **Data.** `celegans_lineage_cells` (2,183 cells), `celegans_packer2019`, `human_cell_turnover`,
  `celegans_fate_rules`, `celegans_terminality`, `gastrulation_census_comparison`,
  `mesoderm_diagnosis_single`, `sox17_repressor_single`, `number_provenance_census`,
  `rule_number_sources_census`; programs in `data/organisms/celegans` and `data/organisms/human`.
- **Requirements.** Every cell by name against Sulston; deaths and fates
  exact; timers with their measured spread; the human body reproduced by
  growth and turnover, not asserted; organism-level uncertainty populated. State on 2026-10-10: the
  worm grows from one cell with factor and exposure rules deciding the fates they can, and a ceiling
  measured on what they decide without the lineage; the human body grows as populations, with
  haematopoiesis its only mechanism module; gastrulation is a demo whose numbers cite no source.
- **Built by 2026-09-11.** Haematopoiesis, the Body as a process-bigraph process gated by a network,
  space (the French flag grown from one cell), the measured reader, Packer 2019 at 84% on
  single-tissue ids, and Digital Development as the published knockout set.
- **Factor fates, contacts, AC/VU, glia and PAR (2026-09-14).** Factors decide 496 of 555 embryonic
  terminal fates; 8 of 8 founder identities and 8 of 8 published knockouts from a curated contact
  table; AC/VU from Z1.ppp in 96 of 199 seeds (48.2%) and the Kimble ablations 100 of 100; sheath glia
  6/22 → 11/22; PAR polarity takes Digital Development from 7 of 11 to 11 of 11. `commitment` and
  `competence` earn their place by ablation (`126e65b`).
- **Fate rules on exposure reads (2026-09-15, `8042a57`).** 501 → 522 of 555, and 476 → 483 of 555
  held out (85.8% → 87.0%), with three negatives recorded; corrected on 2026-09-29 (`b587bb1`): the
  factors still decide the fate of 5, and without the commitment the score falls 522 → 412.
- **`share:` landed (2026-09-15, `f6ca4a8`) and was first used (2026-09-17, `e362ef2`).** The human
  tissue splits state published shares; three populations came back from nothing. Body: 2.83e13
  cells at 20 years.
- **The ceiling on what the worm's factors decide (2026-09-22, `dddda7b`, `d5215bf`, `de7a135`, `9d266a3`).** 522 with the fallback, 381 honest, 136 free, **0 decided without the lookup having
  written the fate first**. The terminality gate's two registered clauses fail (factors 0.8079
  against per-founder depth 0.8773 held out), and generation 7 is at chance.
- **Gastrulation, measured honestly (2026-09-28 to 2026-10-02).** Against the one measured human
  gastrula (Tyser et al. 2021, CS7) the proportions are falsified on all three layers (`d55cb19`). A
  sourced repressor of SOX17 restores a mesoderm band in a registered search (`430b510`, `53d9936`:
  MESP1 and HAND1 tail-min 0.5333 against 0.0000) while deleting endoderm (0.367 → 0.000), no
  comparison with the CS7 range is made, and the repressor is not in the shipped program. 0 of the
  module's 42 number slots are sourced (`00e7028`).
- **The next mechanism chosen in writing (2026-10-10, `460c46f`).** The small-intestinal crypt, the
  largest-turnover population of the counts-only half of the body program, with a closed eight-factor
  sourced list and two phenotypes declared partial in advance, registered before a number was fitted
  (docs/AREA-E-NEXT.md). The worm's next number was set aside because it needs runtime machinery or a
  data source the project does not have.
- **Missing.** The human program is counts, not mechanism, outside haematopoiesis. The worm's next
  real number needs contact and geometry at division, or a measured cell-cycle regulator: the
  contact table is curated, not measured, and no public dataset resolves terminal division. A
  delayed `differentiate` is not built (`recheck: crossings`, `9d42485`, gives a re-decision). The
  external engines (libRoadRunner, MaBoSS) have no adapter: `runtime/sbml.py` and
  `runtime/boolean.py` are stand-ins. `haematopoiesis.to_bio` still emits `fraction:` rather than
  `share:`.
- **Next.** 1. Contact and geometry at division, or a measured cell-cycle regulator, for the worm's
  next real number; set aside on 2026-10-10 for want of both (AREA-E-NEXT.md). 2. The second human mechanism module after haematopoiesis:
  the small-intestinal crypt as registered in `460c46f`, then its run and its mutants. 3. A measured contact table for the embryo over time, to replace the curated one.
  4. libRoadRunner and MaBoSS adapters tested in CI on Python 3.12. 5. A `differentiate` that may be
  read after a stated delay. 6. `share:` in `haematopoiesis.to_bio`.
- **Owner.** lane-organism holds `data/organisms/human/` and docs/AREA-E-NEXT.md from 2026-10-10
  (§5, M4); lanes are assigned in §5.

History: ROADMAP-HISTORY.md §3 E, lines 1345–1655; gastrulation rows in §5 item 10, lines 4312–4320.

### F. Cancer and therapeutics

- **Goal.** From a tumour's alterations to what physically distinguishes
  those cells, what a therapy could reach and which mechanism the biology
  supports, with evidence per claim and a design dataset a binder-design
  engine can consume. It stops at describing the recognition required.
- **Code.** `cancer/*`, `therapeutics/*` (pipeline, mechanisms, scoring,
  evidence, design, providers, expression, neoantigen, structure,
  trafficking, localisation, logic, report, scan, binding, depmap, synthetic_lethality,
  selectivity); scripts `therapeutic_benchmark.py`, `normal_tissue.py`, `synthetic_lethality.py`.
- **Design.** CANCER.md, THERAPEUTICS.md.
- **Data.** `cancer_msk_impact_2017`, `expression_brca_tcga_pan_can_atlas_2018`,
  `therapeutic_benchmark`, `therapeutic_benchmark_v2`, `normal_tissue_atlas`, `synthetic_lethality`,
  `marker_selectivity`; `data/demo/therapeutics/*`, caches under `data/knowledge/therapeutics`.
- **Requirements.** Three target spaces searched in parallel; target,
  binder, mechanism, cargo and effector kept apart; negative set as
  important as positive; patient data levels reported; nothing above its
  data. State on 2026-10-10: the pipeline runs end to end on public tumours and the retrospective
  benchmark reads 10 of 10 cases on all three questions (milestone 1.2, reached 2026-09-27).
- **Benchmarked against known answers (2026-09-11).** Six tumours whose
  driver and approved therapy are public run through the pipeline
  (`scripts/therapeutic_benchmark.py`, result committed, tests read it). All
  six targets are recovered. Both approved antibody targets, EGFR and ERBB2,
  are found as surface targets; the four whose approved drug is a
  small-molecule inhibitor are correctly *not* called surface targets, and
  each still offers the peptide/HLA route, which is the pipeline's
  distinctive claim. The third question, whether the mechanism ranked first is
  defensible, was 4 of 6 and is **6 of 6 since 2026-09-15**; no known defect
  remains, and the pin in `tests/test_therapeutic_benchmark.py` still reads 2
  and may be lowered.

  *(2026-10-10: kept verbatim because the Progress tab quotes it. The pin reads 0, and the benchmark
  holds ten cases.)*
- **Four fixes the benchmark prompted.** (1) A small-molecule precedent
  counted as evidence that a mechanism needing an extracellular epitope could
  reach the target, so a kinase inhibitor "supported" a radioligand against
  cytoplasmic BRAF; such precedent now requires positive evidence of
  reachability. (2) An agonist antibody was ranked first against an activating
  driver, which would push the pathway the tumour already over-drives;
  agonism now needs positive evidence that triggering the target is the
  intent, and is refused on ignorance rather than offered on it. (3) A
  mechanism whose hard requirement merely went unanswered could outrank one
  whose requirements were established; such a mechanism is now capped at 0.25
  and carries the reason. (4) Capping it was not enough: a capped mechanism
  still headed the list, because for a cytoplasmic driver nothing else scored
  at all, so BRAF's preferred mechanism was a blocking antibody and PIK3CA's
  was ADCP. A gate has three answers and not two, so a mechanism whose hard
  requirement is *unanswered* is now provisional and can never be a
  candidate's preferred mechanism; it stays scored, listed and named with its
  open requirement, and `best_mechanism` returns none instead. Both rows
  closed, 4 of 6 to 6 of 6, with no target, rank, class or verdict moved. The
  control that says this is a fix and not a silencing: the same BRAF case
  given an HLA genotype gets a preferred mechanism, `tcr_based` at 0.44, the
  peptide/HLA route, which is the right answer for a protein no binder
  reaches; and EGFR and ERBB2 keep `blocking_antibody` unchanged, their
  surface requirement being answered rather than merely unrefused. All four
  are the same mistake in different clothes: treating an unanswered question
  as permission.
- **Libraries and a healthy atlas (2026-09-15, `9c0e6ce`).** The `cancer.*` libraries and a healthy
  reference of 5,588 genes over 20 tissues, against which the expression scan names FCRL5 and CD19.
- **DepMap synthetic lethality (2026-09-16, `2305bee`).** 44 hits at q ≤ 0.10 with five of seven
  declared controls recovered; **marker selectivity FAILS its claim** (median percentile 25.1 against
  51.5).
- **Every route scored, 9/9 (2026-09-21, `e21d34a`, `3aee02f`, `deee0a2`).** CD19 by an expression
  call (7/7), ERBB2 by a copy number (8/8), and EML4-ALK as a structural variant: ALK is the 3'
  partner, so `intracellular_only` at accessibility 0.0. With the tenth case, the driver calls are six
  point mutations, two copy numbers, one structural variant and one expression call.
- **The evidence tier, as a gate (2026-09-27, `5013406` registered, `bdd5058` built).** The
  registered weighted form fired its own falsifier; as a gate, the twelve rows where a guess
  outranked the target are zero and `BURIED_SURFACE_TARGETS` is 0. Milestone 1.2 reached.
- **The tenth case and the gates that closed nothing (2026-09-28 to 2026-09-30, `028085c`, `2805552`, `aa7f92c`).** Five of the nine targets are `no_established_mechanism`; **the gate closed
  no case**, and on the tenth case (MYC with ERBB2, 0.246 for MYC, last of six) the tiebreak was never
  asked. The committed result reads 10/10/10; ten benchmark cases, not ten validated therapies.
- **The fusion partner's end (2026-10-03, `3fe7c27`).** For a 5' partner or an unoriented fusion the
  surface score returns none instead of reading the full-length annotation; no published value moved.
- **Missing.** A surface score still reads the gene's curated annotation, not the alteration; magnitude
  is a tiebreak and the tier a gate, so the ranking consequence is closed, and whether a surface score
  should read the alteration is a decision. What a deep deletion is actually worth (the dependency it
  creates) is not wired into the pipeline, though `therapeutics/depmap.py` and
  `therapeutics/synthetic_lethality.py` exist. The demo runs at data level 1 only. The scan's healthy
  reference is 20 tissues of population consensus, not a matched normal.
- **Next.** 1. Wire the deep-deletion dependency (DepMap) into the pipeline's evidence, registered
  before it is scored: **done 2026-10-10** (lane-depmapwire; `c2f1fa9` registered, `1181068`, `064ad9e`,
  `e26ebf9`, `52edaac` measured). `therapeutics/dependency.py` turns a deep deletion in the patient's
  copy-number table into dependency candidates from the 2026-09-16 sweep's 44 frozen pairs, each with its
  sweep row and a test on DepMap lines of the stated lineage, in a list kept apart from the targets and
  labelled cell-line or computational, never patient-derived. **The ten benchmark rows are identical in
  every field**; no case carries a deep deletion. The CDKN2A glioma test case yields WDR77 and PRMT5,
  unknown without a lineage and supported in CNS/Brain lines (47 deleted against 28 kept); exploratory,
  after the registered result: in those lines the dependency rides on MTAP. A cell-line dependency is a
  hypothesis for the tumour, not a measurement in it. 2. Decide whether a surface score reads the alteration, then build it, since
  it moves every case's numbers: **done 2026-10-10** (lane-surface; `7c0f234` registered, `1350f51`,
  `7242ea8`, `9602c6d` measured). `therapeutics/product.py` reads the product's extent from the
  reconstructed altered sequence or from `hgvsp` and holds it against the curated topology: a stop
  before the only curated transmembrane segment answers 0.0, an in-frame deletion of the whole
  anchor 0.0, a deletion inside the outward face scores the residues that survive; a frameshift, a
  cut inside the anchor, a splice-site change, a lost start and an indel with no span answer
  unknown, and decay is not predicted and is named. **The ten benchmark rows are identical in every
  field**, and no overall score rose anywhere. 3. A demo case at a data level above 1. 4. Cohort expression for
  the other TCGA studies (§4).
- **Owner.** Unowned since 2026-10-03; lanes are assigned in §5.

History: ROADMAP-HISTORY.md §3 F, lines 1656–1908; the 1.2 notes in §6, lines 5128–5173.

### G. Product: web UI, CLI, packaging, CI

- **Goal.** A geneticist installs it, opens one page and gets every layer with its evidence; a
  developer trusts CI.
- **Code.** `web/server.py`, `web/static/*`, `cli.py`, `jobs.py`, `storage.py`, `results.py`,
  `manifest.py`, `roadmap.py`, `work.py`, `evidence.py`, `.github/workflows/ci.yml`,
  `pyproject.toml`; scripts `check.sh`, `pre-push.sh`, `fast_prepush.py`, `promote_main.sh`,
  `roadmap_truth_pass.py`.
- **Design.** README.md, DATA.md, STORAGE.md, ACTION-PLAN.md UI track.
- **Data.** `data/jobs/*` (git-ignored), `docs/PROGRESS.md` (rendered in the Progress tab),
  `evidence_compiled_genome`.
- **Requirements.** One HTML page, inline SVG, no build step; every control
  explained; jobs visible with progress; the Progress tab reads this
  project's own log. State on 2026-10-10: version 1.0.0, 63 `genomeos` commands and 19 web views;
  the daily CI was green on `06d0761` every day from 2026-10-04 to 2026-10-09; `main` sits at the
  work to 2026-10-01 and `dev` is not promoted since.
- **Project progress as current work and done (2026-09-12).** `genomeos/roadmap.py` reads this file
  into areas, planned steps and finished items, and `genomeos/work.py` keeps the board; both are
  served at `/api/roadmap` and `/api/work`. The bullets of §3 are its parsing contract.
- **The Evidence explorer (2026-09-13, compiled programs 2026-09-17, `cc81fd2`).** One row per stated
  fact; re-run 2026-09-27, 26,845 facts over 41 programs, mean confidence 0.617; the compiled genome
  adds 940,803 facts and not one of them is experimental. Parse failures are reported since
  2026-09-30 (`91b085c`, `854c2a4`).
- **Nightly real-data CI and the Cells view (2026-09-17).** A second daily job runs the readers that
  reach Zoonomia, gnomAD Gnocchi, UniProt and AlphaFold, as `continue-on-error` so someone else's
  downtime does not train everyone to ignore it. The Cells view reads `reader_genome_wide` per cell
  type: the same genome reads 70.3% of its coding genes in hepatocyte and 53.2% in GM12878.
- **The Cell view (2026-09-21, `7abb450`).** `/api/cell`, `/api/cell/programs`, `/api/cell/run`;
  every view promised in the UI track is built.
- **Three status measures kept apart (2026-09-28, `9f18960`).** `/api/state` returns software
  delivered, biology independently validated and the release checks, never summed, with the
  investigations closed on a negative listed apart; it quotes its rows from this file.
- **Missing.** Every view promised in the UI track is now built. The last, **Cell**, shipped on
  2026-09-21 as `7abb450` with `/api/cell`, `/api/cell/programs`, `/api/cell/run` and
  `tests/test_cell_view.py`; this line said it was unbuilt until 2026-09-27. Its "no nightly
  real-data CI" clause was stale too — Next item 4 records it done on 2026-09-17. **What is genuinely
  outstanding: no git tag and no PyPI package** (the version reads 1.0.0 since 2026-09-27), and no
  "known phenotypes it must reproduce" list per library from a biologist, which is a person this
  project does not have rather than a task. *(2026-10-10: "Next item 4" is the history's; the live
  AlphaGenome test stays out of CI because it needs a paid key as a repository secret, and the fast
  pre-push gate ruled on 2026-10-03 is not wired, §5, E1.)*
- **Next.** 1. Back up `dev` to `main` and cut the 1.0.0 tag, both Albert's (§5, A2). 2. A
  PyPI package, outward-facing and the owner's. 3. The AlphaGenome live test in CI, which waits on a
  paid key as a repository secret, Albert's decision. 4. A "known phenotypes it must reproduce" list
  per library, which needs a biologist. 5. The fast pre-push gate wired as ruled on 2026-10-03 (E1).
- **Owner.** Unowned since 2026-10-03; lanes are assigned in §5.

History: ROADMAP-HISTORY.md §3 G, lines 1909–1985.

### H. BioForge (design under constraints)

- **Goal.** In-silico experiments and design search whose outputs are labelled predicted with a
  confidence, validated by reproducing published perturbation results.
- **Code.** `forge/design.py`, `forge/calibration.py`, `forge/network_experiment.py`,
  `organism/forge.py`, `organism/experiment.py`, `certainty.py`, `runtime/debugger.py`.
- **Design.** DESIGN-MINIMAL-CELL.md, BIOFORGE-CONFIDENCE.md, ACTION-PLAN.md Phase 5.
- **Data.** `design_neuron`, `grn_clamp_census_after`; programs
  `data/organisms/human/haematopoiesis_designs.bio` and the worm's design blocks.
- **Requirements.** None stated in this area beyond its goal. State on 2026-10-10: `design` searches knockouts, additions and knobs toward targets under
  constraints and a budget, and network knockouts are ranked against matched random ones; an answer
  carries effect, score and loss, and its probability stays unavailable because no evaluation
  population exists.
- **The `design` block and the published mutant sets (2026-09-11).** BioForge takes experiments as
  input and finds POP-1, PIE-1, GFI1 and IKZF1 as answers; ten worm mutants and nine blood mutants run
  as `experiment` programs in CI.
- **The predicted confidence and its verdict (2026-09-22, `e13df9e`, `eea443e`; corrected 2026-09-28, `3b29509`).** The emitted confidence was a hand-set 0.3 on every answer, and the census found 2
  designs whose published outcome is cited and 0 fixtures pairing a design with an outcome. The
  2026-09-22 verdict `UNCHECKABLE_BY_CONSTRUCTION` was corrected: a constant equal to the base rate is
  calibrated with zero error and AUC 0.5, so the value was unsupported and no evaluation population
  exists (`UNSUPPORTED_VALUE_NO_EVALUATION_POPULATION`, kept beside the first verdict).
- **The Certainty record (2026-09-28, `38087a1` tests, `16a0434` build).** Evidence category, effect
  with its unit, measurement uncertainty, model score and probability are kept apart; a probability
  cannot exist without a calibration naming outcome, population and method. The fixed 0.3 is retired:
  an answer carries effect = wild-type loss − answer loss and score = loss, and the emitted
  `experiment` states no confidence. The six `min(x, 0.3)` caps remain as evidence-quality caps, not
  probabilities.
- **Budget constraints in `design` (2026-09-28, `faf75c4` registered, `6c3c66d`).** A run under a
  budget never answers over the limit and names what the limit excluded: two_intestinal_founders
  answers −POP-1 at limit 1 and the wild type at limit 0. A hard limit hides the fix in the excluded
  list rather than pricing it.
- **Network knockouts in `experiment` (2026-09-28, `8e705f7`; clamp fixed `51c10b8`).** A node held at
  zero, or one edge read as zero, is ranked against matched random knockouts: CycD off turns the Fauré
  2006 cycle into the G1 fixed point, consequence 0.556, p = 0.1 over all nine matched knockouts. The
  rank is coarse on small networks.
- **Missing.** A population of design answers with published outcomes to calibrate against: no test
  asserts the answers of `haematopoiesis_designs.bio` against a published outcome, although CI runs
  its designs. Network knockouts have no CLI or grammar surface.
- **Next.** 1. A population of design answers with published outcomes, to calibrate against (§5 item
  14, M6). 2. A CLI and grammar surface for network knockouts.
- **Owner.** Unowned since 2026-10-03; lanes are assigned in §5.

History: ROADMAP-HISTORY.md §3 H, lines 1986–2057.

### I. The 98% (attribution of the non-coding genome)

- **Goal.** Every UNKNOWN block with a best guess of what it is for the
  organism and the evidence behind it: from "98% no clue" to "98% a very good
  guess", tested by quantities (the house budget), by measured ground truth,
  and by the organism run forward. Albert's framing of 2026-09-11, the DNA as
  the blueprint and the living organism as the built house, is the opening of
  ATTRIBUTION.md.
- **Code.** `attribution/*`, among them `budget.py`, `constraint.py`, `bigwig.py`, `organise.py`,
  `human_panel.py`, `closure.py`, `closure_cell.py`, `crispri.py`, `measured.py`, `targets.py`,
  `compile.py`, `bridge.py`, `runrules.py` and `astrorun.py`; `runtime/grn.py`, `response_map.py`,
  `response_map2.py`; scripts `budget_genome_wide.py`, `re2g_likeforlike.py`, `astroreg2_send.py`.
- **Design.** ATTRIBUTION.md, CRISPRI-RESULT.md, NODES-READER-WRITER.md, ALPHAGENOME.md.
- **Data.** `budget_chr*` (24 of 24), `budget_genome_wide`, `organised_genome_wide`,
  `human_panel_genome_wide`, `enhancer_targets_all_chr*` (24), `vista_genome_wide`, `eqtl_targets`,
  `mpra_genome_wide`, `closure_chr21`, `closure_chr22`, `unknown_coverage`, `crispri_published`,
  `re2g_likeforlike`, `target_prevalence`, `node_independence_audit`, `clause2_matched_control`,
  `gene_rates`, `astroreg2_registration`, `runrules_result`, `closure_cell_result`; the conserved elements per chromosome
  are a cache under `data/knowledge/constraint`.
- **Requirements.** Constraint is read per base and never stored (the Zoonomia
  track is 9.6 GB; chromosome 21 costs 33 MB of ranges); a tier is a best
  guess with a confidence, never a verdict; every threshold is a named constant
  repeated in the result; "neutral" is a claim with evidence (unconstrained, no
  element), not the absence of one. State on 2026-10-10: the genome is budgeted and every enhancer
  element scored by deletion; the attribution is compiled into BioLang programs of 440,589 rules, of
  which 22,580 can be simulated, every rate borrowed from mouse, and one registered run of them reads
  "no difference detected" against measured RNA; the one result scored on held-out
  measurement is the CRISPRi deletion gain, in one cell type; milestone 1.3's clause 2 is not met.
- **The composition budget, genome-wide (2026-09-12).** 1,009 Mb of UNKNOWN blocks in 3,088 Mb of
  genome, 780 Mb of them holding sequence, of which 15.0 Mb are constrained (1.93%, against 10.7% for
  the genome as a whole); the constrained-unknown tier is 1,098 blocks and 32 Mb. After copies are
  set aside, the real unknown is 882 blocks and 30.6 Mb.
- **Measured ground truth: VISTA, GTEx and lentiMPRA (2026-09-12).** VISTA: among elements the deletion
  names a gene for, a strong effect in 52% of positives (422 of 806) against 39% of negatives (271 of
  695); the deletion effect stays a target finder first and an activity call second. GTEx: the deletion's target is an eGene 71.7% of the time
  (n = 1,851) against 66.2% for the nearest TSS (n = 3,725). lentiMPRA: rank correlation 0.443 in
  K562 and 0.323 in HepG2; the model knows the cell, not the amount.
- **Closure at the gene level (2026-09-12 to 2026-09-13).** On a gene's whole regulatory input, the
  cell's own DNase-gated score names the most-expressed of four cell types 41% of the time against
  25% (p 0.001, 141 genes); chr22 32% (p 0.003, 375 genes). The effect is modest (mean rank
  correlation 0.19), and 72 attributions are rejected by name. The chr19 replication was never run.
- **Many human genomes at once (2026-09-14 to 2026-09-16, `17a76d8`).** The HPRC panel of 90
  assemblies on 24 of 24 chromosomes: 1,925,587 storage units, 187,966 executor-ready. The
  constrained-unknown depletion ends at 0.829 of the neutral tier on 19 of 23 (p 0.0013).
- **Every enhancer element deleted, and the executor claim (2026-09-14 to 2026-09-17, `cbff5f6`, `2a07c80`).** 961,227 elements inside their node for 778,780 AlphaGenome requests (the node figure is re-audited in area B). The discovery run read
  +0.268 (p 0.0021) and its hold-out +0.206; the replication read **E2R +0.091 (p 0.0043), verdict
  `negative`**, and **E3R −0.030, verdict `futility`**; E1 genome-wide +0.014 (p 0.38) decided
  nothing. The registered extension returned +0.075, one-sided p 0.00037: the WEAK band by its own
  bar, not the declared success, carried by GM12878 with Jurkat flat. The ENCODE objection is
  narrowed, not closed. Widening further is a no-go at 0 requests (2026-10-02).
- **The tier reading, withdrawn (2026-09-16 to 2026-09-17, `e9245d9`, `4a2c199`, `66a45ac`).** Every
  non-coding tier reads above the matched background, which carries +6 to +11 points of its own; once
  standardised on length, GC and promoter distance, no tier differs from the neutral tier, and what
  holds is only that elements inside UNKNOWN blocks act less than the genome's scored elements.
- **Measurement coverage of the unknown (2026-09-17, `3124bb4`).** 4.554 Mb of the 1,008.8 Mb unknown
  space has been measured by any assay held here, 0.45%; the real unknown is 30.6 Mb with 0.52%
  measured, in 202 of its 882 blocks. The 284,001-oligo library is a design, not an order (D8).
- **The CRISPRi deletion result (2026-09-22 to 2026-10-02, `df6c5f4`, `bda8a83`, `9dedf52`, `bdba855`).** Held-out K562 AUPRC 0.633 → 0.691 once structural zeros were repaired; rebuilt under
  its own name, the K562 gain is +0.1361, 95% [0.0734, 0.2270] over 22 clusters. Against a
  reconstructed ENCODE-rE2G on the same 1,918 K562 pairs (118 positives), the like-for-like
  difference is **+0.0593, 95% [−0.001, 0.1485], "no difference detected"**. GM12878 falls with its
  interval crossing zero. The second cell type, HCT116 (`a39073d`): "passes, not replicated", +0.0222
  with no interval read, from 5 clusters. README's wording and its qualifications are in
  CRISPRI-RESULT.md.
- **Calibration and re-banding (2026-09-22, `5881d45`, `9a44b65`).** The 2026-09-17 calibration
  verdict of failed stands, and re-banding takes the band from 559,607 of 593,765 targets (94.2%).
- **The prevalence term, registered then fitted — it does not transport, and the genome-wide band table is withdrawn as a quotable confidence (2026-09-22, lane-prevalence, `1a1a852` registration, `e4914ef` fit, 0 requests).** The census that had to be written first settles the question on
  its own: **not one screen the calibration is READ on is a screen it was FITTED
  on** — the two sets are disjoint — so a per-screen intercept is a thing that
  exists only after a screen has been run, and says nothing about a target in no
  screen at all, which is almost the whole genome. That was registered as the
  expected answer before the fit and the fit did not move it.
  **Two things the fit found are worse for the published table than the failure
  it was sent to repair.** The single +0.6793 shift that restored 9 of 10 bins
  is **a mean over per-screen shifts that run in both directions**: −1.303 for
  K562_DC_TAP, which is 63% of the held-out population, then +0.945, +2.050 and
  +4.228. **It has the wrong sign for the largest screen**, and the pooled table
  looks repaired only because 38 pairs needing +4.2 are averaged against 1,084
  needing −1.3. And **the genome-wide band table is one screen's**: under
  Gasperini2019's own base rate the top band reads 56,166 against the published
  **55,957 — a difference of 0.37%** — because Gasperini2019 is 207 of the 285
  on-gate pairs. The spread *inside* each screen set is far larger than the gap
  between them: 3.5× across the fitted screens and 99× across the read ones,
  against the ×1.31 the whole failure was named for.
  **The same 593,765 targets band from 27,514 to 198,475 in 0.9–1 depending on
  whose base rate is used — a factor of 7.2.** Under one screen's rate 61.8% of
  targets fall below 0.25; under another, none do. So the prevalence term does
  not unlock that table: **it is the reason the table was never quotable.**
  **Decision, 2026-09-22: the 593,765-target band table is retired rather than
  corrected.** What stands is the re-banding's strata and the union axis —
  **5.8% and 9.70% coverage** — with "no band" said plainly for the rest, and
  **the shape rather than the level**: leave-one-screen-out moves AUPRC by at
  most 0.022 and by under 0.005 on four of six screens, so the ordering is what
  survives. **No level claim is registered again until a screen exists on both
  sides of a split.** If a band is wanted for a planned experiment it is quoted
  as a function of that experiment's expected base rate with the range shown,
  never as one number; `genome_under_each_screen` is the machinery for it.
  The offset itself failed leave-one-screen-out on 3 of 4 clauses (3 of 6
  screens improving against a needed 4, weighted ECE falling 0.0034 against
  0.005, AUPRC falling on three). **The clause that passed is the informative
  one**: median absolute residual **0.3235**, less than half of +0.6793, so a
  screen's own base rate removes rather more than half the level error and
  leaves the rest — and it fails worst on the screen that is 54% of the pooled
  pairs, where offset and fitted intercept double-count.
  **A hazard for every lane, found mid-run and fixed beside rather than in
  place**: `crispri.logistic_fit` is an undamped Newton fit and blows up on a
  near-separable problem — weights to 1e12 and every prediction 0. The first run
  of this lane would have reported "0 of 10 bins, failed" for the opposite
  reason. A backtracking line search fixes it, reproduces the shared fit to 1e-6
  with every offset zero, and that equivalence is now a test. **The shared
  `logistic_fit` is still undamped.**
- **The real unknown against chance (2026-09-27 to 2026-09-28, `d717b28`, `7f7c8c1`, `0c8b82d`).**
  The real unknown's blocks name a target in 62.3% of 531 against 86.0% of 44,100 random windows; on
  chr21 the unknown scoring reads 49.1% against 65.3%; matched on length and local gene density, the
  gap is **−27.25 points** (95% −30.91 to −23.58, n = 531). Milestone 1.3's clause 2 is not met (§6).
- **How independent that measured confirmation is, component by component (2026-09-28, lane-nodeindep, `d3f51f0` registered, `366db9a` repair, `f752e5e` result, 0 requests).** Eleven components audited:
  **nine CLEAN, two EVALUATION-ONLY, one EXPOSED.** The default caller imported in a fresh interpreter
  pulls in 17 modules and **none can reach a CRISPRi source**; the 50 kb floor, the 5 kb merge, the
  CTCF-only filter and the node confidence were written **2026-09-10**, the control seed and draws
  **2026-09-14**, against the benchmark arriving **2026-09-16** and the first containment figure on it
  **2026-09-27**; stage 2 recomputes identically on all 23 covered chromosomes with the AlphaGenome
  archive never opened. The pairs are **471 training and 190 held-out**, so 71% of the claim rests on
  the training split, K562 only. **The one exposure is the choice of caller**: the default was kept over
  six orientation callers on four criteria, one of which is this same statistic on model output, and on
  the 661 pairs **0 of the 6 rejected callers clear zero** (+1.86 down to −14.98, a spread of 20.87
  points). By the pre-registered reading rule **the direction survives and the size carries a selection
  premium**: +5.89 is the maximum of a selected family, not an unbiased effect. So the claim may be
  called *confirmed on a measurement the caller never read*; it may **not** be called independent,
  out-of-sample or held-out. **No independent arm exists at 0 requests**: 358 pairs are needed for 80%
  power at the claim's own effect, and the IGVF MHC screen (19, all in one 4 Mb locus) plus the DC-TAP
  K562 remainder (36) bring 55 between them; both are left unread so a pooled arm stays possible.

  *(2026-10-10, read from `node_independence_audit.json`: the result lists 10 components, 7 CLEAN, 2
  EVALUATION-ONLY and 1 EXPOSED. The counts in the sentence above are kept as written; the result
  file's counts are the ones to quote.)* *(Reconciled 2026-10-10, lane-docfix: the result file is right.
  Its ten statuses are the ones registration `d3f51f0` fixed; `366db9a` changed none, and the file was
  written once at `f752e5e`. "Eleven" counts the rows of the ATTRIBUTION table added at `f752e5e`, whose
  eleventh row ("Which split: 471 training / 190 held-out") has no status; "nine CLEAN" comes from the
  `f752e5e` commit message, whose nine, two and one add to twelve. The reading does not change: the one
  EXPOSED component and its selection-premium reading are the same in every record.)*
- **Compiled programs, the bridge and the rates (2026-09-28 to 2026-09-29).** Every compiled rule
  carries its cell as an executable `when` (R1, `8dcdade`, `6a4a8a9`): 440,589 rules, 212 of them
  experimental. The bridge audit found 0 of 440,589 rules simulable as compiled (R3, `06b7557`); with
  declared rates, 22,580 of 440,589 compiled rules can be simulated, and every rate is from mouse NIH
  3T3 fibroblasts (`ff2062a`, `e921fae`), because no human rate exists to borrow (`4c3a225`). No fixed
  activator-combination rule survives the measured enhancer pairs (`93caf61`). Runnable is not
  validated (S5, `0ea0a90`).
- **Joint inference and independence (2026-09-28 to 2026-09-29).** The R8 pretest retired two score
  transformations (`4098753`) and the coherence pilot stopped on its rule (`197c560`); both are in §5's
  kept rows. No independent target benchmark exists: the IGVF MHC screen is NOT READABLE
  (`042a4aa`). Without AlphaGenome, 93 of 440,377 predicted links survive (C4, `b7e4bf0`).
- **The response map (2026-09-30 to 2026-10-02, `ed734fe`, `5fefcf1`).** Increment 2 covers 110 of
  473 loci with 527 assertions, reproducible on this machine only: a map of where the assays were
  pointed, not of the genome. A fourth increment is a no-go.
- **AstroREG-2 run 2 stopped (2026-10-03).** Registered (`0f4c372`) and sent under Albert's sign-off
  (`f029d47`), it stopped at 363 of its 1,232 answers, 5 of 23 chromosome clusters, below the
  10-cluster floor: the 363 paid answers decide nothing on their own (§5). The resume code
  (`f0f3cbb`) and a post-run runner (`393ccbc`) landed on 2026-10-10.
- **M1, a run of the compiled rules against measured RNA (2026-10-10, `8420637` census, `008f5c6` registered, `747b207` result).** The census first: 0 of 440,589 compiled rules has a gene as its
  source, so a run of the compiled programs is one independent one-gene system per regulated gene.
  The one registered run reads **"no difference detected"**: running the compiled rules did not change
  the squared error of measured log2(TPM+1) against the no-rule baseline by a detectable amount (Delta
  MSE −0.01565 [95% CI −0.04284, +0.01612] over 13,375 gene-tissue pairs, 8,451 genes, 54 GTEx
  tissues). The rules were gated from AlphaGenome tracks trained on GTEx-sourced RNA, so the outcome is
  not independent of their source, and no reading here is evidence of a mechanism.
- **M3, closure at the cell level (2026-10-10, `ffe4919` registered, `80bf3db` result).** The
  registered reading: "No difference was detected at the cell level": over 35 elements whose
  attributed gene is one K562 depends on and whose nearest coding gene is not, the median log2 relative
  guide abundance in the CRISPRi screen against the control is −0.012, against +0.006 over the 5,427
  elements where neither gene is one K562 depends on (one-sided Mann-Whitney p 0.1024). The positive
  control separates in the same population (757 promoter elements, p 5.0e-62), so this is a negative
  result about the attribution, not about the assay. No difference detected is not no difference.
- **Missing.** Every target claim rests on the reused ENCODE CRISPRi file and one cell type. Only
  0.52% of the real unknown is measured, in 202 of its 882 blocks; deletion requests on them are
  stopped (§5 item 10's stop list) and the oligo library is not ordered (D8), so clause 2 moves only on
  public measurement. Closure has run at the gene level on two chromosomes and once at the cell level, on 35
  elements of one screen, and not at the organism level. No simpler genome has been read.
- **Next.** 1. Compiled rules that run (§5, M1): after the registered run's "no difference
  detected" (`747b207`), the next step, raising the simulable share or the strength-dependent
  combination form R3 left open, is registered before it is run. 2. AstroREG-2
  run 2, resumed or closed (§5, A1), which waits on Albert's own words: the resume code
  landed in `f0f3cbb`. 3. The model's second targets as a shortlist (the 52 off-gate pairs above a
  drop of 0.1, refused a band at 26 and 24 pairs), held against an independent perturbation; no level
  claim is registered again until a screen exists on both sides of a split. 4. Milestone 1.3 clause 2
  from public measurement (M5, lane-measuredunknown): a public MPRA or CRISPR screen whose elements
  overlap the constrained-unknown blocks, searched for first and registered only if found. 5. Closure beyond the
  gene level (M3): the registered cell-level run detected no difference on 35 elements (`80bf3db`);
  what follows it, a larger population or another screen, and then the organism level, is registered
  before it is run. 6. The remaining
  one-target consumers (§5 item 10, W2): nine moved on 2026-10-03 (`7b1db55`, `601261d`, `6faddae`);
  `crispri_contact` and the five then held by other lanes (compile, organise, argmaxcell,
  context_evidence, cli) remain; `genome/motifs.py` and `genome/regulation.py` stay frozen inside the
  paid study's closure, and `report.py` and `knowledge/across.py` are deliberately not moved. 7. Measured RNA over the 69 syntax
  candidates. 8. Predicted chromatin tracks and in-silico mutagenesis over the regulatory blocks, which
  need a self-hosted AlphaGenome (§7). 9. Simpler genomes, where closure can be complete: **first test done 2026-10-10** (lane-simplergenome; `72ed222`, `595d7c2` registered, `eede7d9` result; docs/SIMPLER-GENOMES.md): on 5,814 C. elegans putative enhancers accessible in one of five tissues (Serizay et al. 2020), the domain rule's gene is expressed highest in the element's tissue for 54.6% against 51.7% for nearest TSS, +2.89 points [+1.83, +3.91], above a matched-distance control; accessibility is not function (§5 item 13, I). 10. A
  repeat-aware VNTR reading of the panel.
- **Owner.** Live lanes on 2026-10-10: lane-readercompile holds `attribution/compile.py` and
  `attribution/readergate.py` (M2b), lane-measuredunknown holds `attribution/measured_unknown.py` (M5),
  and lane-resumeprep holds `attribution/astrorun.py`, `attribution/astropost.py` and
  `scripts/astroreg2_send.py` (A1). lane-runrules (M1) and lane-closurecell (M3) finished that day. The
  rest of the area is unowned since 2026-10-03; lanes are assigned in §5.

History: ROADMAP-HISTORY.md §3 I, lines 2058–3293; the area's §5 rows are items 1–3b, 6, 7, 10–14.

### J. The grammar by comparison (species × people, origin, duplication, operators)

- **Goal.** Infer the genome's grammar the way an undocumented executable is
  reverse-engineered from its many versions: compare across species and
  across people, find what is conserved and what varies, and turn the
  comparison into syntax (invariant), value slots (variable inside a
  constrained frame), operators (recurring motif logic), libraries with an
  age (the deepest clade that has them) and cross-references (libraries
  gained and lost together). Albert's framing of 2026-09-12 (a language has
  syntax, operators and values; eye colour keeps the syntax and changes the
  value) and the two DNA-as-code conversations of the same day are the
  brief; GRAMMAR-BY-COMPARISON.md holds them against the code.
- **Code.** `attribution/variation.py` and `genomeos variation --chrom C`, `attribution/syntax.py`
  and `genomeos syntax --gene G --chrom C`, `genome/duplications.py`, `genome/motifs.py`,
  `knowledge/homology.py`, `knowledge/profiling.py`, `knowledge/across.py`, `knowledge/satmut.py`,
  `decompile.py`, `attribution/motif_grammar.py`, `attribution/motif_transfer.py`,
  `attribution/satmut_vista.py`, `molecules/graph.py`; scripts `origin_genome_wide.py`,
  `satmut_grammar.py`, `profiling_genes.py`. It builds on `attribution/bigwig.py`,
  `attribution/constraint.py`, `genome/mouse.py` and `attribution/compile.py`.
- **Design.** GRAMMAR-BY-COMPARISON.md; GENOME-AS-CODE.md §7 (the mapping
  table with the construct per row); SEQUENCE-GRAMMAR.md §2 and §5 (the
  motif scan, `requires:`).
- **Data.** `variation_chr*` (24), `variation_genome_wide`, `variation_vista_genome_wide`,
  `origin_genome_wide`, `duplication_chr*` (24), `duplication_genome_wide`, `motifs_chr*` (24),
  `motifs_genome_wide`, `promoter_composition_genome_wide`, `profiling_genes`, `across_ZRS`,
  `across_HERC2_OCA2`, `motif_grammar`, `motif_transfer`, `satmut_vista`.
- **Requirements.** Every comparison enters as evidence with a confidence
  (D4, D42, D43); within-human constraint is never read as proof of
  neutrality; curated duplication before heuristics (D20's order); the
  species axis and the human axis are reported side by side in the same
  table, never merged into one score; each step is scored on chr21 first and
  run genome-wide second, as the budget was. State on 2026-10-10: all seven steps of the area
  are built and genome-wide; of six falsifiable claims, one survives (§5 item 10's count), and the
  last substantive section of GRAMMAR-BY-COMPARISON.md is from 2026-09-17.
- **The human axis, copies and origin (2026-09-12).** Elements by case name a gene 73.9% (syntax,
  825), 56.0% (relaxed), 50.0% (recent), 45.3% (tolerant); the constrained-unknown tier is 492
  relaxed, 346 tolerant, 77 syntax and 30 recent of 948 measured, chrY unmeasured, not zero. 216 of
  1,098 constrained-unknown blocks (19.7%) are mostly copies. Origin over 19,392 coding genes, with
  64,050 `paralogue_of` edges.
- **Operators, searched (2026-09-12).** No TFClass family pair recurs beyond its matched expectation
  in any of the 42 libraries; REST stays enriched in the nervous system at 5.7x. On the ZRS and a
  148-element VISTA panel, held-site density does not separate.
- **Measured readout and the registered grammar test (2026-09-13 to 2026-09-16).** Against lentiMPRA,
  no factor family survives correction, so §12's reading is unsupported by measurement. **The
  registered claim fails in all three cell lines**: arrangement adds +0.0021 [−0.0049, +0.0085] in
  K562, −0.0056 [−0.0120, +0.0014] in HepG2 and +0.0022 [−0.0050, +0.0097] in WTC11; what worked is
  the counts, +0.234 Spearman in K562 and +0.235 in HepG2.
- **The count positive transfers (2026-09-17, `b7e1405`).** Conservation plus counts 0.655 against
  conservation alone 0.596, +0.060 [+0.022, +0.099] on 643 held-out elements (351 positive), weak and
  on an unmatched comparison; over 152,556 windows of the unknown blocks the model cannot tell real
  sequence from its own shuffle.
- **Gene-level profiling and the in-silico mutagenesis panel, both negatives (2026-09-17, `9d216e3`, `08bcf14`).** The profiling is a negative about its own method: 14 of 42 libraries sit at or
  below a null made of real genes, and nothing is called a dependency. The ISM panel
  reads +0.0453 at one-sided p 0.1741, below the registered weak bar of +0.05, which is failure and
  was the registered expectation.
- **Missing.** Enhancer pair logic read outside promoters (GATA plus T-box in heart has never been
  read in enhancers); a species-tree null, the only thing that would let a pair be called a
  dependency; and a model that reads the sequence itself, scored on §17's own held-out elements, the
  one readout §17 names as able to overturn its negative.
- **Next.** 1. A species-tree null for library pairs (§5, M6): **done 2026-10-10** (`eb2967d`; pairs: no reading, the null is not calibrated; libraries: 26 of 42 clear, all 21 of §19's plus 5); a calibrated pair null follows (lane-pairnull). 2. Enhancer pair logic outside
  promoters, starting with GATA plus T-box in heart. 3. A sequence model scored on
  GRAMMAR-BY-COMPARISON §17's held-out elements.
- **Owner.** Unowned since 2026-10-03; lanes are assigned in §5.

History: ROADMAP-HISTORY.md §3 J, lines 3294–3562.

## 4. Data jobs to be done

Brought current on 2026-10-10. Every finished job is condensed into one row per kind with the date
it finished; the rows as they were written, with their notes and corrections, are in
ROADMAP-HISTORY.md §4 (lines 3563–3603). Results land in `data/results` and are committed.

| Layer | Done | Job | Notes |
|---|---|---|---|
| Reference layers per chromosome: sequence, GENCODE rows, ENCODE elements, RepeatMasker, the curated UNKNOWN pass, CTCF domains, the anatomy inventory and the UNKNOWN classification (`fetch_<chrom>`) | 25 of 25 (2026-09-11) | done | sequences stay local in `data/reference`; CTCF domains on 24 (chrM has none); 98.5% of the UNKNOWN space classified genome-wide |
| Reader (open nodes per cell type) | 25 of 25 (2026-09-11) | done | eleven cell types on that date; a twelfth and thirteenth biosample added 2026-09-22 (area B) |
| Proteome compiled, translation verified against UniProt, knowledge graph, modifiable sites | 25 of 25 (2026-09-11) | done | 19,478 coding genes compiled; 90.9% of 19,249 canonical translations identical; 42,045 nodes and 330,018 edges; 96,362 sites on 13,083 proteins (area C) |
| Enhancer targets predicted (AlphaGenome) and all-elements deletion scoring (every ENCODE enhancer inside a node) | 24 of 24 (2026-09-16) | done | 961,227 elements for 778,780 requests; summaries committed, element tables local and packed (area I) |
| Human panel (HPRC release 1, 90 assemblies) | 24 of 24 (chr2 on 2026-09-16, `17a76d8`) | done | `human_panel_genome_wide` (area I) |
| HG002 twin | 22 of 22 autosomes (2026-09-11) | done | chrX and chrY are not phased in the GIAB benchmark, so out of scope rather than pending (area D) |
| Composition budget: constraint per UNKNOWN block, a tier per block | 24 of 24 (2026-09-12) | done | constrained-unknown 32 Mb (area I) |
| Segment parser scored against GENCODE | chr21 and chr22, series closed 2026-09-12 | done | every run kept side by side (area B) |

| Job | Needs | Owner |
|---|---|---|
| One-off jobs, each finished on its date: telomere length from HG002 reads (2026-09-12), measured Hi-C boundaries against the node edges (2026-09-12), Packer 2019 disagreements resolved to one labelling depth (2026-09-11), the mouse chromosome as the second mammal (2026-09-11), Sender & Milo turnover as maintenance timers in `bio.std` (2026-09-11), the retrospective therapeutic benchmark (six cases on 2026-09-11, ten since 2026-09-30) | **done**, each on its date | — |
| AlphaGenome live predictions for UNKNOWN regulatory blocks | `ALPHAGENOME_API_KEY` and Albert's grant. The stop list of §5 item 10 (2026-09-27) says no AlphaGenome requests on the 882 blocks, and the execution plan in §5 allows no model requests unless Albert grants them, so this row waits on him and is argued against by the plan | Albert |
| Cohort expression for the other TCGA studies (`genomeos cancer expression --study`) | 30 s each; pick the studies | unowned since 2026-10-03 |

---

## 5. Next steps, consolidated: the execution plan from 2026-10-10

**How this section works.** This is the order of work from the review of 2026-10-10 (item 14 of the
previous list, kept in [ROADMAP-HISTORY.md](ROADMAP-HISTORY.md)), approved by Albert the same day
("proceed"). Each numbered item names its owner (a lane on the work board), the files it holds, what it
waits on, and what counts as done. The Progress tab reads the items below.

Four rules govern the whole plan:
- **Model first.** At least two of every three lane slots go to model work, meaning work the model can be held to. One slot at most goes to engineering, and only when that work unblocks a model item or a decision.
- **Register before reading.** Every test commits its registration, with its frozen analysis code, before any outcome is read. A negative result is a result and is reported in its registered words.
- **One paid run at a time.** A model request is sent only under Albert's own words, given in his own session, and recorded verbatim.
- **The coordinator dispatches and reviews.** Lanes write code and commit on `dev`; the coordinator pushes between pieces, folds the rows in here, and keeps at least three lanes working while open items remain.

**Wave 0 — started 2026-10-10.**

1. **Done 2026-10-10 — M1, compiled rules that run** (lane-runrules; `8420637` census, `008f5c6` registered, `747b207` result; 0 model requests). Census: 0 of 440,589 compiled rules has a gene as its source, so a run has no network effects and each regulated gene is its own one-gene system. Every resolved pair runs "on" at its borrowed rate in every cell (0 of 160,392 has no basal rate), so a run cannot predict "off". A run's only cross-tissue prediction is element present over element absent, which equals its input observation; the planned on/off test was therefore replaced by a level test, on GTEx v8 because the reader's cell lines have no cached RNA values. Registered reading: "no difference detected: running the compiled rules did not change the squared error of measured log2(TPM+1) against the no-rule baseline by a detectable amount (Delta MSE -0.01565 [95% CI -0.04284, +0.01612] over 13,375 gene-tissue pairs, 8,451 genes, 54 GTEx tissues; control p = 0.000999)". Genes were higher in their rule's tissue in 41.7% of pairs. The rules come from AlphaGenome tracks trained on GTEx, so the outcome is not independent of them, and no reading is evidence of a mechanism. Open, as its own registration: the promoter-competition form (no pair in this population has two mechanisms). Next: item 7b.

2. **Done 2026-10-10 — M2, the `reader` construct, engine half** (lane-readerlang, `5682e6c`, `56f0e65`). `reader <CellType> { open/poised/closed/unknown: elements }` with evidence and confidence, and `through: E` on an activates/inhibits rule. The gate is at rule selection (`Module.active_rules`): a rule runs only where the cell's reader states E open. Poised, closed and every kind of unknown hold it at strength 0, reported with the reason and never run. It is held, not dropped, because a gene with no activator runs at full max_rate: on the demo, closed and held 0.05, dropped 60.05. A reader states accessibility, not function. The language has 27 block kinds and BioIR 26 types; `bio test` 45 files at 246/246 (the 44 existing unchanged); 34 new tests. Spec in docs/BIOLANG-v0.5-READER.md, demo in data/demo/reader_gating.bio. Open for Albert, each built on the most conservative choice: whether unknown holds or refuses the run, whether poised ever lets a rule run, and the version number (v0.4 §10 decision 1). Next: item 7.

3. **Done 2026-10-10 — M3, closure at the cell level** (lane-closurecell; `ff04f9b` census, `ffe4919` registered, `80bf3db` result; 0 model requests). `attribution/closure_cell.py` chains element → gene → cell: the sweep's K562-track prediction, K562's DepMap 24Q4 dependency on the attributed gene, and a readout the project had never seen, ENCODE's genome-wide CRISPRi proliferation screen in K562 (ENCSR162HMH) against its own control (ENCSR625TAF); 1.39 GB was streamed with every md5 matching ENCODE's, then discarded. Of 440,252 attributed elements, 47,912 carry at least three library guides. The arms cross the attributed gene's dependency with the nearest coding gene's, which leaves 35 elements where the attribution stands alone, against a floor of 5,427, with 757 promoter elements as the positive control. The detectable separation, 0.622, was stated in advance. Registered reading: "No difference was detected at the cell level: over 35 elements whose attributed gene is one K562 depends on and whose nearest coding gene is not, the median log2 relative guide abundance in the CRISPRi screen against the control is -0.012, against +0.006 over the 5427 elements where neither gene is one K562 depends on (one-sided Mann-Whitney p 0.1024) … The positive control does separate in the same population (757 promoter elements of genes K562 depends on, median -0.047, p 5.0e-62), so the readout is responsive and this is a negative result about the attribution, not about the assay. No difference detected is not no difference." Beside it, deciding nothing: depletion tracks the attributed gene where both genes are essential (145, p 1.1e-05), not proximity (147, p 0.874). Open: the files carry no non-targeting guide row, so the assay's zero is not established; mediation through the gene cannot be established from a pooled screen.

4. **Done 2026-10-10 — A1 prep, AstroREG-2 run 2 made resumable on Albert's word; nothing sent** (lane-resumeprep; `f0f3cbb`, `393ccbc`; no key read, no network, ledger and cache untouched). `--resume-authorisation WORDS` has no default and is refused when blank, or for a fresh or completed ledger; on `--send` the words go verbatim into the ledger before the first new charge. Only this ledger's 363 paired answers are exempt from the cache rule; an uncharged cached element, a request-364-shaped element, and an element answered in run 1's ledger are still refused. A dry check writes nothing (both ledgers byte-identical afterwards). At the real numbers a resume buys exactly 868 in plan order and never request 364, and the pilot does not re-fire. Seven mutations were caught; 281 targeted tests pass. `scripts/astroreg2_analyse.py --run 2` refuses today and computes item (h) as registered. **The primary claim's assembly is refused:** the recording path writes `by_cell` only for K562, HepG2, GM12878 and IMR-90, so the registered `crispri` path finds no astrocyte value on any paid answer (the values exist only in the full vectors), and the screen pairs carry no distance. A blind amendment of the primary claim's assembly therefore comes before any answer is read, whether or not the run resumes. Comparison (b), the 667-column split (adapter v2 records 371 columns) and the secondary are not computed, each with its reason. The lane left an orphan `.claude/worktrees/resumeprep-check`, which waits on the removal tool's orphan path.

5. **A1 — Resume run 2 — Albert's decision.** Recommendation: resume. 868 requests remain, exactly the cap's remainder, and the 363 already paid for decide nothing alone. Two things are needed from him:
   - his own sentence authorising a resume of run 2 on its ledger;
   - who signs as supervisor, because the code names `genomeos-fe` and that session is gone.

   Then, in order:
   1. one commit changes `SUPERVISOR_SESSION`, quoting his words;
   2. the closure files are held fixed (the closure differs from the signed one in five files today: the two A1-prep files, and `ir/model.py`, `lang/grammar.py` and `lang/parser.py` from the `reader` construct);
   3. the named signer reviews the dry run and commits its own sign-off line (`by` equal to `SUPERVISOR_SESSION`, `run_id` 2, "dry run reviewed", the digest, and the signed closure equal to the sender's closure at that commit);
   4. CI passes on that sha;
   5. from a worktree at that sha, on mains power, the operator runs `caffeinate -ims uv run python scripts/astroreg2_send.py --run 2 --resume-authorisation "<Albert's words, verbatim>"` as the dry check, then the same with `--send`;
   6. a read over all answers;
   7. the blind amendment of the primary claim's assembly (item 4), registered before any answer is read;
   8. the registered analyses, through item 4's runner;
   9. the ledger committed as its own commit.

   To close instead: commit the ledger and record the run as stopped incomplete.

6. **A2 — Back up `dev` to `main`, then tag 1.0.0 — Albert's step.** On 2026-10-10 CI was green on `c00250f` (run 38040550668: test, test-bare and live-readers all success), and `enforce_admins` was on. `scripts/promote_main.sh --dry-run c00250f` refused: "c00250f does not contain origin/main (a1ce607); moving main there is not a fast-forward". `main` carries PR #13's merge commits (`054eb46`, `a1ce607`), which `dev` does not; `main`'s files are identical to `747b7ff`'s, which `dev` contains. How to proceed is Albert's decision. After the merge: `scripts/promote_main.sh --verify <sha>`, then the 1.0.0 tag on the merged commit (milestone 1.0 was reached on 2026-09-22, and §6 records the tag as the owner's). After that, promote daily, at most every two days, on a green CI sha.

**Wave 1 — as wave-0 lanes finish, refilling slots so three lanes run and two are model.**

7. **Done 2026-10-10 — M2b, the `reader` construct, application half** (lane-readercompile; `98a1bcb` emitter, `1a6205a` census and RNA sources, `4c84b48` registered, `ff35b09` result; 0 model requests). `compile.gate_program` gives the 24 compiled programs one `cell_type` and one `reader` block per reader biosample (13, with the narrowPeak accession as evidence; 459,449 elements stated in every cell; open from 38,537 in GM12878 to 138,063 in cardiac muscle) and puts `through:` on all 440,589 rules. `not_assessable` maps to `unknown`; reader v1's "no peak" maps to `closed`, and the program says so. In K562 the gate holds 333,754 rules closed and 62 unknown, and runs 106,752. The gated programs are built in memory, so no tracked program and no census pin moved. RNA: 27 ENCODE polyA gene quantifications (294.7 MB, md5-checked, git-ignored). Registered reading: "worse: the reader-gated run raised the squared error of measured log2(TPM+1) across the thirteen reader cell types above the no-rule baseline's (Delta MSE -0.06154 [95% CI -0.08665, -0.03963] over 6,773 gene-cell pairs, 521 genes; control p = 0.00799) and was further from measured RNA than the ungated run as compiled (MSE difference -0.05389 [95% CI -0.0753, -0.03508])". Post hoc, not registered: the information term (+0.0089) beats 992 of 1,000 shuffles, but the cost of the shifts is about eight times as large. Accessibility is not function. Open: genes with several elements cannot run until a combination form is registered (only 738 of 18,973 genes have one rule); no element-level poised call exists; the bridge drops `through` on the rules it rebuilds. Next: a multi-element combination registration, together with 7b.

7b. **Done 2026-10-10 — M1b, compiled rules with gene sources** (lane-regulators; `c267869` census, `c04cf5e` registered, `a951f52` result; 0 model requests). Census: 46,889 of 440,589 compiled rules (10.6%) gain a factor source in their own cell, all from measured ENCODE ChIP occupancy and none from function (30,571 also carry the factor's JASPAR motif). Only 71 of 317 contexts have any ChIP data. A motif alone matches a median of 361 factors per element, so a motif cannot choose a source. Direction from curated GO: 122,590 pairs activate, 16,042 inhibit, 472,589 unknown. Registered test on 34 ENCODE knockdown units in K562 and HepG2: "the compiled regulator assignment did not anticipate a factor's knockdown better than another factor's assignment did: mean stratified AUROC 0.5022 for the factor's own assigned genes against a median of 0.5016 for another factor's, control p = 0.4026 over 34 factor-cell units; the paired difference against genes near the factor's peaks, -0.0150 [95% CI -0.0228, -0.0071], is reported and decides nothing". Emitting these rules waits on `compile.py` and two engine questions: a rule with no stated strength runs at 1.0, so `strength: unknown` is needed (item 13, A); and unknown-direction rows need the grammar token in Albert's queue (item 16). The specification is `regulators.COMPILE_SPEC`. Replogle Perturb-seq was not read, because its perturbed rows are N1b step 2, the owner's decision.

7c. **Done 2026-10-10 — M1c, where the compiled layer loses its signal** (lane-signaltrace; `9637584` registered with the 1,918 pair ids pinned by sha256, `1103719` result; 0 model requests). The 1,918 held-out K562 CRISPRi pairs (118 regulated) were traced through each compiled stage with the anchor's own metric (AP gain over activity and distance). The anchor reproduced exactly: +0.1361 [0.0749, 0.2258]. Target choice kept 39 pairs (35 regulated), gain +0.0792 [0.0282, 0.1562], below the anchor by −0.0569 [−0.0970, −0.0248]: 1,209 pairs sit on elements whose compiled rule names another gene, and 670 have no compiled element. The bridge in K562 keeps none: every kept gene has more than one compiled rule, so none is identifiable, and the reader, regulator and run stages receive no pairs; as compiled, the bridge keeps 6 (+0.0260 [0.0027, 0.0711]). Registered reading: "the measured signal is lost at s4_bridge (the strength the bridge fits in K562): none kept: the stage keeps none of the 1,918 held-out K562 pairs (0 of the 9,237 fit rows), so no compiled prediction reaches the measurement at this stage; gain over activity and distance +0.0000 [95% CI +0.0000, +0.0000], against the anchor -0.1361 [95% CI -0.2255, -0.0745]; it was first weakened at s2_target (target choice)". One defect, reported beside the readings: `save_result` quarantined the first write because the run read its own registration without declaring it; `d1df762` declares it and writes the payload byte-identical, with no re-run. Consequence for the plan: the multi-element combination (lane-combine2) is the stage that lacks identifiability, and a second registration re-running stage s4 on the same pinned pairs with that form would be the decisive test.

7d. **Done 2026-10-10 — M2c, the multi-element combination** (lane-combine2; `bc9fc7a` code, `e79302d` census, `917ce25` registered, `948cfc9` result; 0 model requests). Under the reader's hold a held rule sits at strength 0, so the mean is linear and a joint fit runs genes with several elements: 4,937 runnable under the mean and 2,687 under promoter competition, of 17,983 fittable (18,229 of 18,973 genes have 2+ rules). The competition form A=W/(K+W) is Bothma 2015 (eLife 4:e07956, Methods eqs. 7-9), with one transcriptional efficiency per gene and THETA=0.5 declared; nothing was fitted. It is opt-in in `runtime/grn.py`, and the default is unchanged. Registered reading: "the runtime's mean: worse; the reader-gated run of genes with two or more compiled elements raised the squared error of measured log2(TPM+1) across the thirteen reader cell types above the no-rule baseline's (Delta MSE -0.4435 [97.5% CI -0.4843, -0.4085] over 61,425 gene-cell pairs, 4,725 genes; control p = 0.3); the promoter-competition form (Bothma 2015, declared half-occupied promoter): worse; ... (Delta MSE -0.8102 [97.5% CI -0.8994, -0.7252] over 32,513 gene-cell pairs, 2,501 genes; control p = 0.734). Head-to-head over the 2,405 genes in both populations, the promoter-competition form was further from measured RNA than the mean (MSE(mean) - MSE(competition) -0.6944 [95% CI -0.7613, -0.6288])." Neither form beats its shuffled control. The `grn.py` edit moves the AstroREG-2 closure digest, so item 5's re-sign covers one more file. Open: a super-additive form (it needs an engine floor), the roughly 13,000 genes no form can fit, and item 7e.

7e. **Done 2026-10-10 — M2d, the bridge stage re-run with the joint fit** (lane-s4joint; `9f199d6` frozen, `a9055e0` registered, `86ddf73` result; 0 model requests). Same 1,918 held-out K562 pairs, metric, bootstrap and decision as 7c; every check before a joint stage was read passed (the anchor exact, s2 and s3 equal to 7c's, every run within 3.6e-15 log2 of its closed form). Registered readings: under the mean, "the measured signal does not survive the bridge: at s4_joint_mean (…), lost: no gain over activity and distance detected, +0.0051 [95% CI -0.0001, +0.0232] on the 1,918 held-out K562 pairs (1 kept by the stage, 1 of them regulated), and below the anchor, -0.1310 [95% CI -0.2146, -0.0724]"; under promoter competition, "none kept: the stage keeps none of the 1,918 held-out K562 pairs (2 of the 9,237 fit rows)". Chain, both forms: LOST_AT s4_joint, first weakened at s2_target. **Why:** 38 of the 39 target-choice pairs (mean) and all 39 (competition) are dropped as `basal_not_positive`: the gene's single-element removal effects add up to more than its whole expression. The pairs reach 210 genes with a median of 28 rules each, and the joint fit runs 10 (mean) and 2 (competition). An additive decomposition of a gene's expression into its elements' single-deletion effects is inconsistent with those effects on these genes.

7g. **Done 2026-10-10 — multi-target compiled choice** (lane-targetchoice; `3d85122` training census, `79aa1e2` frozen code, `b746a23` registered, `715d830` result; 0 model requests). The parameters were fixed on the training split only (9,237 fit rows; every compiled target is the first coding gene at the compiler's own bar, 0 disagreements): primary, every protein-coding gene with `MIN_EFFECT` 0.1 log2 or more, compiled target first, nothing tuned; secondary, top-k with k = 6 by a rule written beforehand, deciding nothing. Disclosed in the registration: a dry run on a stand-in population of training pairs read NO_BETTER before registering, and no parameter changed after it. Registered reading: "better than the one-target choice and below the anchor: keeping every coding gene at the compiler's bar ranks the 1,918 held-out K562 pairs better than keeping one target per element, +0.0318 [95% CI +0.0027, +0.0721] in average precision (67 pairs kept, 57 of them regulated, against 39 and 35); its gain over activity and distance is +0.1110 [95% CI +0.0582, +0.1859], the one-target choice's +0.0792 [95% CI +0.0282, +0.1562]; it stays below the anchor, -0.0250 [95% CI -0.0447, -0.0067]". Of the 1,209 pairs on an element naming another gene, 28 are added and 1,181 stay unlisted because their gene is below the bar; 670 have no compiled element. Only target choice changes here; the bridge still loses every pair (7c, 7e). Next: the compiler emits every coding target at its bar (a `compile.py` change, registered against the compiled programs' census tests), once 7f says whether any composition carries the signal through the bridge.

7f. **M2e — A saturating composition, registered** (lane-saturate, running since 2026-10-10). 7e's cause names the one structural hypothesis left on the bridge: single-deletion effects that sum past the whole are what a saturating (sub-additive) composition produces, e.g. a gene's output as a saturating function of its elements' summed input, so each single deletion removes more than its additive share. The test fixes the form and its parameter source in advance (a cited saturating form, or one parameter set from the training split only), and is registered before it is run: first how many of the 210 genes become identifiable, then 7c's stage s4 on the same pinned pairs, read in 7c's words. If it fails too, the compiled execution path is recorded as unable to carry the deletion signal, and the plan moves model effort to the parts with positive readings (item 13, I's worm test; item 8's mechanism modules).

8. **Done 2026-10-10 — M4, area E's second human mechanism module: the small-intestinal crypt** (lane-organism; `460c46f` choice and registration, `7646bd7` module, `ad8c662` result; 0 model requests). The choice was written before anything was built (docs/AREA-E-NEXT.md): the worm's next number needs runtime machinery or a data source nobody has, so mechanism went where the body program is counts. The closed list is eight gating factors, each with a fetched loss-of-function citation (TCF7L2, ASCL2, RBPJ, ATOH1, NEUROG3, GFI1, SPDEF, SOX9), and five excluded with the reason. **All eight registered knockouts reproduce the registered pattern, 32 of 32 lineage checks, and nothing about a mutant was fitted.** One number was fitted, the villus residence at 4.70 d, which closes the module to the published 1.2e11 cells and 2.2e10 per day; that is the fit, not a result. Two registered readings about the formalism were wrong and are pinned in tests: conservation does produce Gfi1-/- supernumerary enteroendocrine cells (6.25×; direction published, magnitude unmeasured) and inflates enterocytes 1.64× without ATOH1, where the source reports them unaffected. `bio test` 246/246 → 288/288. Open: colon and stomach are still counts; the tuft lineage waits on a sourced share or residence; gates are presence or absence, so no hypomorph; body.bio is deliberately not wired, because tissues.bio counts the same cells.

9. **Done 2026-10-10 — M5, milestone 1.3 clause 2 from public measurement: inventory only, no test registered** (lane-measuredunknown; `0d05f97` code, `7906142` result; 0 model requests). Pooled public element-perturbation screens with a gene readout tested an element against a coding gene in **13 of 882** real-unknown blocks (14 at most); Gasperini 2019 gives 12 and the CRISPRi benchmark 3. The smallest calibrated design needs 75, so the shortfall is **62 blocks** (518 for the half design). Reporter assays reach up to 385 blocks and growth screens 153, but name no gene. GTEx fine-mapped eQTL reaches 879 blocks by a common-variant proxy, but has no calibrated sizing and no chosen control, and its track was not read. Clause 2 stays not met.

9b. **Done 2026-10-10 — M5b, an eQTL arm for milestone 1.3 clause 2: calibrated, and no fair test exists, so none was registered** (lane-eqtlcal; `68a273b` frozen code, `6033113` calibration; 0 model requests). Calibrated outside the 882 blocks only (blocks ±1 kb excluded, 32.3 Mb): on CRISPRi positives, a fine-mapped GTEx eQTL (PIP ≥0.5) names the measured gene at 3.61% (19 of 526) against 2.59% at windows matched on length, common-variant count and TSS distance (excess +1.03, −0.84 to +3.14); tested negatives score 0.115% (15 of 13,031); at block length, 6.65% against 5.26%. For 80% detection, all attributions right needs 2,000 to 3,000 blocks and half right more than 5,000; the blocks offer 179 to 284, which give 0.07 to 0.09 power at half. An eQTL here adds little beyond distance to the gene, and it is association, not mechanism. **Clause 2 stays not met, and no public route to it is known today** (item 9: 13 of 882 blocks by perturbation screens; this item: underpowered by an order of magnitude by eQTL). Disclosed: a test-only partial decode of chr21 records inside blocks (nothing printed or counted) was narrowed to a region between blocks in `16a415f`. The raw dbSNP file (1.8 GB) was deleted after distillation, with its URL and sha256 kept.

10. **Done 2026-10-10 — E1, the fast gate is wired** (lane-fastgate, `41b14b7`; `scripts/pre-push.sh`, `scripts/fast_prepush.py`, `scripts/last_green_ci.py`, CONTRIBUTING steps 2 and 4). A `dev` push runs `ruff check` and `ruff format --check`, then the always-run tests and the test files changed since the last green CI, in both the normal and store-free legs. Its verdict is `fast-<tree>.json`, scope "fast", which `genomeos/verdict.py` never reads as a project verdict. The full suite returns when the daily CI is red or unread (`last_green_ci.py --daily`, older than 26 h counts as unread), on `main` and `promote-<sha>` branches, and under `GENOMEOS_FULL_GATE=1`, the switch for a paid-path sha. The trade, as ruled: a red in an unchanged test file does not block a `dev` push. 21 planted tests pass, 13 of which fail against the unwired hook. The full `scripts/check.sh` was green on tree `85eb696` (`tree_begin == tree_end ==` the commit tree; 5,909 passed, 0 failed, 15 skipped). Until a daily CI reading is recorded, every push stays on the full gate; the last-green record is `747b7ff`, so the first fast window would be 156 changed test files. After a fast push, `push_own.sh`'s retry and `GENOMEOS_SKIP_CHECK=1` refuse, because they need a project verdict.

11. **Done 2026-10-10 — E2, the engine package writer** (lane-enginepkg; `48c65e3` writer, `050a9b2` result). `scripts/package_engine.py` records the package path relative to the repository ("build/biolang/biolang"; a build outside it is recorded relative to `--out`, with a flag), and only the default build writes `data/results/engine_package.json`; any other `--out` or `--name` is a scratch run that writes beside its build and says data/results was not touched. The canonical result was rebuilt in a clean worktree of `48c65e3`: 12 of 12 checks, `dirty: false`, engine suite "366 passed, 3 skipped, 3 xfailed in 34.45s". It replaces both the stale result of 2026-09-29 and the uncommitted scratch rebuild of 2026-10-02, which was not committed.

**Wave 2 — from the middle of October, in this order unless a result changes it.**

12. **AstroREG-2 analyses**, after item 5 completes: the registered primary analysis and item (h), each from its committed registration only, reported verbatim. Then Albert's decision on the second cell type (725 requests), with the run-2 outcome in hand.

13. **M6 — The dormant areas, one per free model slot** (from 2026-10-10 under Albert's "start working unstoppable on each lane", each area has a lane: A lane-engdebts, C lane-roadrunner, D lane-pgs, H lane-forgecal, J lane-speciesnull; and area F's Next 1, surface scoring that reads the alteration, is lane-surface).
    - A, the engine's debts: **done 2026-10-10 except two (lane-engdebts; `8df8270`, `be2d00a`, `a39af38`, `0ada532`, `5794d8b`).** The grammar builds its IR list from the dataclasses (32, the hand list had 26), and `check()` fails on an unexported one. `strength: unknown` is accepted on activates and inhibits rules and held at strength 0 at rule selection and in the located runtime, named with its reason by `bio check` and `bio run`; a rule with no strength line still runs at 1.0, because about 32 inline programs in tests and scripts rely on it while every such rule in the 47 committed programs states one. Event rate units must be per time and `param` units come from a closed set; `mrna_half_life` and `protein_half_life` refuse non-hour units, because the runtimes read hours without converting (`30 min` ran as 30 h). `bio check` lists the 98 blocks stated above v0.4 §2's 0.3 cap on ungated constructs (92 are worm fate decisions on integrated reads) without changing any. The labels stay 0.4, held equal by a test. `bio test` 47 files, 288/288 before and after.
      - open, Albert's (item 16): whether v0.4 §2's cap is amended to govern values a construct computes rather than an author's stated confidence (the lane's recommendation: amend; enforced, it would turn 28 cited experimental fate rules into 0.3), and the version number;
      - **v0.4 stage 3 built 2026-10-10 (lane-stage3; `e09edf7` registered, `f99d634` measured):** 10 of 10 registered clauses met (32 and 28 ATP per glucose, 6 O₂, anoxia ×16.0 and ×14.0 when uptake is not the limit, the red cell 1.62 at every oxygen level, the heteroplasmy threshold between 0.80 and 0.82 against a registered 0.8125, stage 4's division paid from ATP regenerated over 24 h, 13 of 13 refusals); the gate as §6 words it is not passed (the budget's magnitude and split and the red cell's ranges are not askable; the hierarchy is an input; Pasteur by construction), so decision 4 stays open. A pool's `regenerates` rate now keeps its unit (`/min` ran as `/h`; no committed number changes). `bio test` 54 files, 354/354.
      - open: `policy: optimise` compiles and the runtime refuses it, and `regenerates:` is stored and read by nothing;
      - **done 2026-10-10 except one file (lane-strengthread, `f06f95b`):** the Cell view's rule table lists a rule of unknown strength as `held`, its strength shown as "unknown" and the reason stated (before, the page crashed on it), and the bridge reports `strength_unknown` instead of fitting the rule's note or reading 0 as no effect; 7 planted tests, 5 of which fail on the old code; 270 web and bridge tests pass. `forge/design.py` (`Knob.get`) still reads it as a number: its two-line change is specified and waits on the forge licence-header commit in Albert's queue (item 16).
    - I, simpler genomes (from 2026-10-10, area I's Next 9): **the worm's first test done** (lane-simplergenome; `72ed222` code and note, `595d7c2` inventory, census and registration, `eede7d9` result). Without AlphaGenome, the human attribution reduces to its node model and distance; the worm has no CTCF, so the node was redefined from gene structure (host gene, else the gene the element sits upstream of). Registered reading: "On 5,814 C. elegans putative enhancers accessible in exactly one of five tissues (Serizay et al. 2020, young adults), the gene the domain rule attributes is expressed at 5 TPM or more and highest in the element's own tissue for 3,174 of 5,814 (54.6%), against 3,006 of 5,814 (51.7%) for the gene with the nearest 5' end; the rules differ on 1,248 elements; difference +2.89 points, 95% +1.83 to +3.91 over 1-Mb clusters. On this endpoint the domain rule is ahead of nearest TSS. The domain rule's genes match the element's tissue more often than genes at a matched distance (24.1%, one-sided p 0.001), so position carries tissue information here. Accessibility is not function: this says the attributed gene is expressed where the element is open, not that the element acts on it." The positive control responded (coding promoters 70.1% against 31.6% at matched distance); transcription through gene bodies and outrons is not separated. Next: nematode phyloP tiers by range read, then Drosophila (lane-simpler2).
    - D: a population distribution for the polygenic scores (today they are raw sums).
    - J: **done 2026-10-10 (lane-speciesnull; `5a87e2f` frozen code, `b37b198` registered, `eb2967d` result).** The tree is Ensembl Compara's pan tree (r63) with the vertebrate tree (r116) grafted in, all 355 species placed. Library pairs (§14): "No reading: the species-tree null is not calibrated on this data. Of 860 pairs of random gene sets drawn to the tested pairs' origin-stratum make-up, 205 (23.8%) have a tree-null p at or below 0.05, above the registered ceiling of 10%, so the null finds random genes alike more often than it should and a library pair's p under it is not evidence. §14's 35 of 860 pairs stand as written, with their caveat." Library coherence (§19): "26 of 42 libraries' members co-occur beyond what the species tree and real genes predict", under a stochastic Dollo model against 200 birth-node-matched random gene sets, which is all 21 of §19's clear libraries plus 5 (about 2.1 expected with no effect); "it names no pair, calls nothing a dependency, and an orthologue's presence is not function". Then the calibrated pair null, **done 2026-10-10** (lane-pairnull; `5b614d0`, `3abde37` registered, `c61a335` result; GRAMMAR-BY-COMPARISON §23): Pagel's lambda, binomially scaled, REML per side, calibrated on random pairs before any pair was read (44 and 66 of 860 at p ≤ 0.05, none at q; the Brownian null 205); "154 of 860 testable library pairs are more alike than age and shared ancestry predict … 6 of the 35 pairs §14 called more alike than age, and 148 that §14 did not call. The other 29 of §14's 35 do not survive this null … No library pair is called a dependency." The calibration did not cover two large library deviations, so a mixed calibration comes before the 154 are read further (lane-pairmix).
    - C: **done 2026-10-10 (lane-roadrunner; `1dcf998` adapter, `b933129` references and tolerance committed before any scored run, `bdb30d3` result, `c0c03f2` divergence now raises):** `runtime/roadrunner.py` sits behind the `roadrunner` extra and is chosen by name; the stand-in stays the default and an absent extra raises, never falls back. Against BioModels curation figures, with the tolerance fixed before the runs (a level or time within max(10%, 3 figure pixels), a period within 3%), **Schoeberl2002 passes 33 of 33 checks at three EGF doses on libRoadRunner and fails all 33 on the stand-in**; Kholodenko2000, Elowitz2000 and Hornberg2005 pass on both engines, which agree within 6e-5 of each species' range. **Correction to the record:** the "flat cascade" of 2026-09-11 was RK4 diverging within 7 steps, with NaN clamped to zero by `max(0.0, nan)`; a diverged run now raises `FloatingPointError`. Open: MaBoSS. *(Done 2026-10-10 by lane-ciextras: CI's `test` job installs the extra, so the 6 libRoadRunner tests run there, `e4239ac`; the packaged `biolang` declares `roadrunner`, with a test holding every engine third-party import to an extra declared as the application declares it, `cba4077`.)*
    - H: **done 2026-10-10 (lane-forgecal; `c35ee24` author's confidence kept, `5df76f2` minimum set before counting, `953c61d` design answers pinned in CI, `04dbc35` census).** `to_experiment()` now carries the design author's own stated confidence, labelled as the author's, and BioForge's probability stays None. The calibration minimum, set before counting, is 38 independent (design question, answer, published outcome) triples, the smallest n at which the exact 95% interval is narrower than 1/3 at every count. Census: 4 design triples and 27 inverse-design triples from experiment blocks, **all circular** (same source, the model's own acceptance check, or a hand-written candidate list): 0 independent, shortfall 38, so no calibration was registered and no rate computed. An answer that reproduces its model's acceptance phenotype is not a prediction. The `no_lymphocytes_at_all` design ties IKZF1 with TCF3, broken by label order. Next for H: held-out publications set aside before a module is written.

14. **The known-locus benchmark** (previous item 4): the frame drawn for the question it asks, and the re-score plan registered by lane-readerfix, which awaits the benchmark owner's decision. Recommendation: Albert names the coordinator as owner, and the frame is registered as its own lane.

15. **Engineering backlog, one slot at most** (the documentation findings below are lane-docfix, and the engine package writer is lane-enginepkg, item 11):
    - the remaining one-target consumers (`syntax_tiling` and the two deliberately not moved);
    - the CPTAC seven profiles, uncommitted at 14.3% mismatch;
    - the removal tool's orphan path, for `.claude/worktrees/lane-quarantine`;
    - disk: 92% used and 37 GiB free in the morning of 2026-10-10, then 93% and 33 GiB free by midday (`dd84686`);
    - from the documentation pass (`dd84686`); **done 2026-10-10 by lane-docfix** (`d58343c`, `7843655`, `c8f5f22`; LICENSING's engine-test count is now 33, 32 headed, after `bc9fc7a`):
      - a comment in `tests/test_engine_licence_headers.py` still gives the old engine-test count (31/30; it is 32 tests, 31 Apache-headed);
      - DESIGN-MINIMAL-CELL does not yet carry the SOX17-repressor mesoderm result (registration `430b510`, result `53d9936`);
      - the ACKNOWLEDGEMENTS tables lag the code (DepMap, CPTAC and Ochoa, 4DN, AlphaMissense, PGS Catalog and lentiMPRA are used but not tabled);
    - figures that disagree with their result files (`74a4a0c`); **done 2026-10-10 by lane-docfix** (`7ffcc31`, `67971c5`: each figure now names the population it describes; ACKNOWLEDGEMENTS records that AlphaMissense states two different terms, and the stricter one governs):
      - CRISPRI-RESULT's coverage row "every held-out K562 pair" gives the weighted rebuild [0.0734, 0.2270], while that row's own unweighted estimator rebuilds to [0.0749, 0.2258] (`crispri_published_v2.json`, `coverage_arms_k562_heldout.arm1_all_pairs`; README's figure is correct);
      - ALPHAGENOME's constrained-versus-uniform table says "of 2,995" where `enhancer_targets_genome_wide.json` gives 2,994, and its 70.7% and 86.0% (`constrained_targets_genome_wide.json`, `uniform_samples`) differ from the 71.2% and 90.2% in the table above it;
      - **reconciled 2026-10-10 (lane-docfix; the result file is right, see area I):** area I's node-independence bullet says "Eleven components … nine CLEAN, two EVALUATION-ONLY, one EXPOSED", while `node_independence_audit.json` lists 10 components: 7 CLEAN, 2 EVALUATION-ONLY, 1 EXPOSED (kept verbatim with a dated note in area I, not reconciled);
    - **Done 2026-10-10 (`7e15a22`, lane-enginepkg):** `genomeos/verdict.py`'s `working_tree_hash` built its index with `git add -A`, so a local ignore rule that matches a TRACKED file leaves that file out, and the clean worktree's hash stops equalling the commit tree. On 2026-10-10 a `msg-*.txt` line the coordinator added to `.git/info/exclude` at 10:22 matched the tracked `msg-lane-silenceragree-manifest-1791050698.txt` (from `d352da7`), and every push gate refused every sha with "no green verdict exists for the tree being pushed" until the line was removed. lane-fastgate found it. It now adds every tracked path an ignore rule hid with `update-index`; the planted test fails on the old code and passes on the new, and a clean worktree of `050a9b2` hashes to its commit tree;
    - from the codebase-map refresh (`fc21c29`):
      - `benchmark/loci_gene_input_plan.py` is imported by nothing;
      - **done 2026-10-10 (lane-headers; `4c6e2d5`, `ab1d9d4`, `544ffe9`, `a60f30e`, test `fefe807`):** 94 of 249 application modules and the 7 `std/*.bio` engine files carried no licence header, and no test checked presence. 81 files gained the AGPL-3.0-or-later line as a one-line insertion, and `tests/test_licence_header_presence.py` requires one licence per `.py` file under `genomeos/` and `scripts/`, matched to its side. Still pinned: 30 files in the AstroREG-2 sender's signed closure (a header changes the closure digest, so they go with Albert's run-2 decision), 2 forge files held by lane-forgecal, and the 7 `std/*.bio` files (no SPDX convention for `.bio` yet);
      - `genome/` imports six `attribution` modules, and `predict/alphagenome_adapter.py` imports the paid-run guard `attribution.astrorun`, so the application's layers run both ways;
      - the `parser.py` docstring says v0.3 while the grammar says 0.4 and the parser accepts v0.5's `reader`;
      - the map itself is 62 KB.

**Standing.**

16. **Albert's decision queue**, each with its recommendation; the previous item rows hold the detail.
    - Resume of run 2 and its signer: item 5.
    - Promotion merge and the 1.0.0 tag: item 6.
    - The second cell type and E2/E3: hold both until item 12.
    - S8 floor values.
    - A grammar token for an unresolved direction, and `dirv2` in programs.
    - Re-cutting the published chr21 programme to v3.
    - The cell-provenance proposal.
    - N1b step 2.
    - The hook rule for `push_own.sh`.
    - `.gitignore` lines for the caches.
    - The `.env` backup.
    - The `manifest_headlines.json` writer. Recommendation: a new-named registry through a committed writer, with the old bytes kept and marked superseded.
    - The `reader` construct's three choices, each built on the most conservative option (item 2): does an unknown chromatin state hold a rule or refuse the run (built: hold); may poised ever let a rule run (built: never); the language version number. Recommendation: keep hold and never; number the language and the BioIR file format separately, because a BioIR file can now carry an `UNKNOWN` strength that older readers take as a string (item 13, A).
    - Whether v0.4 §2's 0.3 confidence cap is amended to govern values an ungated construct computes rather than an author's stated confidence (item 13, A). Recommendation: amend; enforced as written it would turn 28 cited experimental worm fate rules into 0.3.
    - With the run-2 decision (item 5): 30 files in the AstroREG-2 sender's signed closure still lack their licence header, because adding it changes the closure digest; they can be headed in the same commit as the re-sign.
    - A forced commit for the two forge licence headers: lane-ciextras inserted `# SPDX-License-Identifier: AGPL-3.0-or-later` at line 1 of `genomeos/forge/__init__.py` and `design.py` and emptied `DEFERRED_HELD` in `tests/test_licence_header_presence.py` (licence tests 37 passed); `commit_own.sh`'s stale-base guard refuses because emptying the pin removes 4 lines of `fefe807`, and the permission classifier denied the lane's `--force`. The change is left uncommitted in the tree, draft message `msg-lane-ciextras-forge-headers-1791636617.txt`. Recommendation: allow it; the removed lines are the pending entries the headers resolve. After it, `forge/design.py`'s `Knob.get` takes the two-line UNKNOWN check lane-strengthread specified (`f06f95b`'s report).

17. **Keeping the record current.** At each lane's end:
    - the coordinator folds its row into the item above and marks it with `**done YYYY-MM-DD**`;
    - PROGRESS.md gets the day's entry;
    - the area's Next list in §3 changes in the same commit, because the Progress tab reads it;
    - the board is retired weekly (`genomeos work retire`).

### The review of 2026-10-10: the state it measured and what it found

*Albert asked for a full review on 2026-10-10; the order of work it set is the plan above, which supersedes the A, B and C lists the review first carried. The full item as first committed is in `c00250f`, and ROADMAP-HISTORY.md keeps the roadmap as it stood before the restructure.*

**State, measured 2026-10-10.**

| Fact | Value |
|---|---|
| commits on `dev` | 1,556; the last at 2026-10-03 22:50 (`06d0761`) |
| `main` | `a1ce607` (PR #13, the work to 2026-10-01); `dev` not promoted since |
| tests | 5,784 collected; 5,571 passed in CI on `06d0761`; six daily CI runs green since; lint and format clean |
| `bio test` | 44 programs, 227/227 checks |
| product | 63 `genomeos` commands, 19 web views, version 1.0.0, no git tag |
| milestones | 0.9, 1.0, 1.1, 1.2, 2.0 reached; 1.3 partly met, clause 2 not met (§6) |
| jobs | 91 in the registry, none running or stalled |
| disk | 92% used, 37 GiB free |

**What the review found.**

   1. **AstroREG-2 run 2 stopped at 363 of 1,232 answers, and no row recorded the launch until this one.** Albert's words, first-hand on 2026-10-03 (`f029d47`): "Record the sign-off and start run 2 under my approval: at most 1,232 requests, its own ledger, pilot check on the first 10 answers." The coordinator sent it from `06d0761` at 23:49:31; the pilot passed (no `pilot_stop`, all ten answers usable, the supervisor's independent read agreeing). `data/ledgers/astroreg2_run2.jsonl` holds 364 requests and 363 answers, every answer usable; its last line is request 364 (element `EH38E3663949`) at 23:57:59, with no answer, no `error` line and no stop event. The sender writes an `error` line on any Python exception, Ctrl-C included, so the process most likely ended by a signal or the host sleeping (not established). Request 364 counts as charged by design: the ledger line is written before the request. The ledger is untracked, and `data/results/astroreg2_send.json` was never written. The answers cover chr1 to chr4 and 14 of the 27 elements on chr5: 5 of the 23 chromosome clusters. Under the registered rules that is below the 10-cluster floor, so item (h) would read "inconclusive" and the primary analysis has no ci95. **The 363 paid answers decide nothing on their own.**
   2. **Model work is about a third of the commits, against two thirds.** Of 684 commits from 2026-09-28 to 2026-10-03, a random sample of 60 subjects (seed 7) reads as 20 model work, 7 records of model results and 33 verification, tooling or process. Commits touching science code were 14% to 26% a day, and the share did not rise after the 2026-10-02 rule.
   3. **The model work went to one area.** Area I took almost all of it after 2026-09-21. D (the `twin` package untouched since 2026-09-10), J (the last substantive section of GRAMMAR-BY-COMPARISON is from 2026-09-17) and E (last model result 2026-09-22, the human body program untouched since 2026-09-17, no owner) are dormant; C and H have no commit since 2026-09-28.
   4. **The largest gap against §1 is execution, not decoding.** 22,580 of 440,589 compiled rules can be simulated, and every rate they use is borrowed from mouse NIH 3T3 fibroblasts (R3, `e921fae`). No fixed activator-combination rule survives the measured enhancer pairs (R3, `93caf61`), and the strength-dependent form is open as its own registration. `context_evidence` gates nothing, because the `reader` construct is unbuilt (area B, Next 1). Outside haematopoiesis, the human body program is counts, not mechanism.
   5. **The plan's text has fallen behind its rows.** Several areas' Next lists name work that was done and recorded only in §5 rows: H's three items (`16a0434` and the two R4 rows), J's whole list, A's refusal of unknown block keys, F's neoantigen provider, and E's `share:`. The Progress tab reads those lists, so it overstates what is left. The header, §1 (which still said MIT), §2, §2.1 and §8 were stale and are corrected in the same commit as this item.
   6. **Engineering loose ends.**
   - The fast gate ruled on 2026-10-03 is not wired: `scripts/pre-push.sh` still runs the full `scripts/check.sh`, and `1b2e768` was the ruling's row only.
   - `data/results/engine_package.json` carries an uncommitted scratch rebuild (a temporary path, a dirty tree, the engine suite at "1 failed, 349 passed"), because `scripts/package_engine.py` saves the canonical result on every run, a `--out` scratch run included.
   - 395 commit-message files and one already-applied patch sit untracked at the repository root.
   - 121 work-board entries belong to sessions that are gone.
   - Three push worktrees are prunable, and `.claude/worktrees/lane-quarantine` holds an orphan data copy whose removal waits for the removal tool's orphan path.
   - This file is 849 KB and 5,192 lines, too long for a person to read and costly for a session.

**Items 1–13 of §5 as they ended**, with where their full rows are kept. Nothing here is current
work unless the row says open.

| Item | As it ended | History (ROADMAP-HISTORY.md) |
|---|---|---|
| 1. The executor claim's external endpoint | **Done 2026-09-17** (`2a07c80`): the WEAK band by its own bar; the ENCODE objection is narrowed, not closed | lines 3655–3673 |
| 2. The panel's chr2 and the 24-chromosome rollup | **Done 2026-09-16** (`17a76d8`): chr2 read with every control passing, and it moved nothing | lines 3674–3677 |
| 3. Whether the neutral tier is neutral | **Answered 2026-09-16** (`e9245d9`): every non-coding tier reads above the matched background, so tier statements are read against the neutral tier | lines 3678–3684 |
| 3a. The 69 syntax blocks | Registered (`9b2a4ea`), then cancelled unrun and the observation withdrawn (`4a2c199`, `66a45ac`); 2,400 requests not spent, nothing claimed | lines 3685–3690 |
| 3b. Measurement over the real unknown | Coverage measured (`3124bb4`, 0.45%); the 284,001-oligo library designed and not ordered (D8, 2026-09-27) | lines 3691–3767 |
| 4. The known-locus benchmark with more loci | **Open.** No row since 2026-09-22; lane-readerfix registered a re-score plan on 2026-10-03 (`1f6d4ab`), and running it is the benchmark owner's decision | lines 3768–4010 |
| 5. Albert's language decisions, then economy stage 2 | Stage 2 built 2026-09-21: gate (a) passed on three clauses of four, gate (b) **FAILED**; the leftovers and stage 3 are in area A | lines 4011–4060 |
| 6. Four lanes opened 2026-09-16 | Done: all four landed on 2026-09-16 (`d3c5c52`, `5a31c39`, `2305bee` and area J's), three of them negatives | lines 4061–4121 |
| 7. Two lanes opened 2026-09-17 | Done: calibration landed with its reliability claim FAILED on the held-out population (`bdc2364`); area J's transfer landed (`b7e1405`) | lines 4122–4160 |
| 8. Disk | **Open**: 92% used, 37 GiB free on 2026-10-10 (§5) | lines 4161–4170 |
| 9. A coordinating session, from 2026-09-21 | A standing process rather than an item; §8 says how the work is organised | lines 4171–4233 |
| 10. The execution plan of 2026-09-27 | Waves 0–3 as below; its stop list and lane rules are kept in the history | lines 4234–4387 |
| 10 · W0 node containment, `unknown_scoring`, repository size | **Done** 2026-09-27 (`05b71ae`, `2b77663`; `77de3bb`, `7f7c8c1`; `753edb8`) | lines 4339–4342 |
| 10 · W0 records truth pass | Ran 2026-10-03 (lane-truthpass, checker `7d21174`); the audit row still reads UNVERIFIED because it was written before | lines 4341, 4325 |
| 10 · W1 deletion result publishable | **Partial**: like-for-like against ENCODE-rE2G +0.0593 [−0.001, 0.1485], "no difference detected" (`bdba855`); the second cell type is Albert's, spend withheld | line 4343 |
| 10 · W1 reader family | **Done** 2026-09-27 (`63929e1`, `b13143e`); milestone 1.1 reached by it | line 4344 |
| 10 · W1 therapeutic evidence tier | **Done** (`5013406`, `bdd5058`); milestone 1.2 reached by it; `BURIED_SURFACE_TARGETS` 0 since 2026-09-27 | lines 4345–4346 |
| 10 · W1/W2 BioLang 2.0 parts one and two | **Done** 2026-09-27 (`5720c7c`, `14f9ada`); milestone 2.0 reached | lines 4347–4348 |
| 10 · W2 remaining one-target consumers | **Open, unheld**: 18 live consumers counted, 11 free, 3 inside a frozen closure (`24adf33`); see area I | line 4349 |
| 10 · W2 D: GIAB Q100 phasing | **Done** 2026-09-27 (`457702c`, `4d63910`) | line 4350 |
| 10 · W2 C: Ochoa 2020 and CPTAC | **Done** 2026-09-27/28 (`58f84fb`, `869158d`, `cae8d7d`, `5938b5d`), target wording withdrawn | line 4351 |
| 10 · W3 E2/E3 widened | **Closed as a registered no-go** 2026-10-02: E2R `negative`, E3R `futility` | line 4352 |
| 10 · W3 A: `none` facts and v0.4 stage 4 | **Done** (`61c5838`, `b41ddd6`; `4ff3cd4`, `73af2d5`): the stage's own dilution gate refuted, as registered | lines 4353, 4307 |
| 10 · W3 H: budget in `design`, network knockouts | **Done** 2026-09-28 (`faf75c4`, `6c3c66d`, `8e705f7`, `66675d8`) | line 4354 |
| 10 · C4 hide evidence, predict it | **Done** 2026-09-28 (`cc824f4`, `b7e4bf0`) | lines 4355, 4655–4657 |
| 10 · C5 natural variation as experiments | **Partial**: probe done (`3156e5b`), the build gated on a phase that does not follow | lines 4356, 4658–4659 |
| 10 · C1–C3 the coherence pilot | **Discontinued** 2026-09-29 (`a9660ff` registered, `197c560` result) | lines 4357–4359, 4765–4768 |
| 10 · C6 hidden redundancy | **Blocked, dead premise**: it depends on the discontinued pilot | line 4360 |
| 11. The external review of 2026-09-28 (R1–R9) | Every item done or closed; R8 closed on its pretest, reopened as S3 and then discontinued with the pilot. Its table is kept below because the Progress tab reads it | lines 4388–4458 |
| 12. The second external review (S1–S8) | S1, S2, S4, S5, S6 and S7 done; S3 discontinued (the pilot); S8 built (`bee8f20`), its floor values the owner's. Its table is kept below | lines 4459–4643 |
| 13. The coherence programme (C1–C6 and what followed) | The pilot discontinued; the rows that followed it, each closed: | lines 4644–4784 |
| 13 · C4 and C5 | C4 done (`b7e4bf0`); C5 probe done (`3156e5b`) | lines 4655–4659 |
| 13 · audits A, B and C | A done (accepted at `7ae0681`); B parked; C no-go (`e7f4c3e`) | lines 4660–4668, 4762–4764 |
| 13 · the response map | Increments 1 and 2 done (`ed734fe`, `d2043c3`; `5fefcf1`); MAP-4 a no-go at 0 requests | lines 4669–4673, 4689–4690, 4747 |
| 13 · N1 and N1b | N1 closed 2026-10-01 under its stop rule: completed, inconclusive. N1b found T = X·√n calibrated on the non-targeting rows; its step 2 is the owner's decision | lines 4674–4682, 4709, 4713 |
| 13 · wiring and paired enhancers | Wiring closed 2026-10-01 as infeasible on the K562 CRISPRi benchmark (`31e624e`, `f58a2a3`); paired enhancers a no-go accepted by the owner | lines 4683–4687 |
| 13 · E2/E3 | A no-go at 0 requests (2026-10-02) | line 4688 |
| 13 · P1 | Done 2026-09-29: replicated, not activity (`d19ccac`, `a48e9e8`) | line 4769 |
| 13 · the second cell type | 725 requests approved on conditions, then recommended against; the owner's | lines 4758–4760 |
| 13 · AstroREG-2 run 2 | The launch checklist written 2026-10-03; run 2 then started and stopped at 363 of 1,232 answers (§5, A1) | lines 4724–4725 |

**The rows below are kept as history in this file because the Progress tab quotes them.** They are
the two external reviews' tables (items 11 and 12), item 13's opening and the rows that closed
investigations on a negative, as they ended. A long detail cell is cut to a pointer into
ROADMAP-HISTORY.md; every lead, code cell and state cell is unchanged. None of it is current work.

*From the 2026-09-21 preamble of §5 (ROADMAP-HISTORY.md lines 3630–3634):*

Two more that were never wrong, only unsaid: **0.45% of the unknown space has ever been measured by
any assay held here**, and the compiled genome states **940,803 facts of which none is experimental**.
The experiment that would change the first is designed and unrun (284,001 oligos with mappability at
four read lengths). Area J's gene-level profiling closed as a negative about its own method: 14 of 42
libraries do not beat a null made of real genes.

*From item 5, as it ended on 2026-09-21 (ROADMAP-HISTORY.md lines 4011–4060):*

   **Three things the gate bought, none of them about allocation.** At a real transcriptome's
   demand §5.1's capacities **are not even binding**; the pools had to be driven 1e4 times harder
   before anything went short, so at proteome scale the layer currently does nothing. **`competitive`
   is as gene-blind as `proportional`** — it saturates on the pool's total demand rather than per
   demander — which **closes decision 5 with a negative**, since gate (b) was the named instrument
   for choosing between policies that turn out to make identical predictions. And the one unfitted
   number lands: total protein per cell from the declared ribosome capacity reads **1.36e10 against
   Milo 2013's 4e9 to 1.2e10**, when it could have been wrong by decades. The newly registered
   quantity, the slope of log10 protein on log10 transcript at **0.515 [0.499, 0.531]**, sizes what
   a future per-demander policy would have to supply and allocation supplies none of.
   **Still open and named rather than dropped:** gate (a)'s magnitude clause, which needs a
   published curve digitised with its axis stated, and decision 2's proposed volume arm, which was
   never what resolved decision 2 and has still not been run. BIOLANG-v0.4-ECONOMY.md §9.2,
   `data/results/abundance_gate.json`.

*From item 10, 2026-09-27 (ROADMAP-HISTORY.md lines 4239–4245 and 4307):*

   **Where the project stands, on two denominators that must be read together.** Against the plan as
   written it is **about 70% built** — 99 of 129 area items (77%), 20 of 24 data jobs (83%), milestones
   at 64% counting a partial as half. Against the goal it exists for it is **0.52% measured**: 160,447
   of the real unknown's 30,602,182 bases have ever been read by any assay this project holds. Per area,
   each on its own stated denominator: A 88% of requirements, B 78%, C 63% of goal clauses, D 78%,
   E 64% of its agenda with **0 fates decided independently of the lineage**, F 9/9 routes, G 100% of
   requirements with **release readiness 0 of 3**, H 50% of its goal, I 89–93% of its ledger, J **100%
   built and 17% surviving** (one of six falsifiable claims).
   *(2026-09-28, named so the figure can be checked: the three are the ones area G's **Missing** bullet lists —
   **no git tag** (1.0.0 is untagged), **no PyPI package**, and **no "known phenotypes it must reproduce" list per
   library from a biologist**, which area G notes is a person this project does not have rather than a task.)*

   | wave | lane | why here | lane-sessions |
   |---|---|---|---|
   | 3 | ↳ **v0.4 stage 4 done 2026-09-28 (lane-v04, `4ff3cd4` registered, `73af2d5` built)** | `lang/grammar.py`, `ir/model.py`, `runtime/division.py` | One new key, `partition: duplicate | contents | binomial` on an `event`, with the checkpoint expressed as the event's own `when` through `ir.model.matches` rather than a private mechanism. `runtime/division.py` conserves every species exactly, splits binomially below the regime threshold under `units: copies` and by exact halves above, and refuses half a molecule. The default stays `duplicate`, so no committed program changes. **The registered negative is the finding: the spec's own gate for this stage is refuted, not passed, and the refutation was derived from the arithmetic and committed before the instrument existed.** §8 asked that dilution separate a stable protein from a short-lived one; partitioning is gene-blind, so it contributes **0.0** of the 9.392136 decades between a 240 h and a 6 h protein (identical whether contents are partitioned or duplicated), and in a steady state dividing **compresses** the separation 8.80-fold, 40.0000 → 4.5455. What passes is the stage's own claim: a protein given no half-life falls by exactly 2^−8 partitioned and not at all duplicated; 2,000 of 2,000 binomial splits of 20 copies conserve exactly (mean 9.8895 against 10, variance 4.8993 against 5); variance is exactly 0 above the threshold; an unmet checkpoint fires 0 divisions. **8 of 8 registered clauses met**, 24 tests added, check.sh 221/221, engine package 12/12 (36 → 38 files, 177 → 192 packaged test functions). **Reported and excluded from the verdict** because the registration did not name it: the demo declares a division cost of 1e10 ATP against a 3e9 pool, so the resource checkpoint refuses every division — a pool is a standing stock and a cost is a draw over a cycle, and nothing converts between them until stage 3 supplies ATP over time. It cannot claim any measured half-life (none is held), nor that real partitioning is binomial (no partitioning-error data held) | the engine result file is the coordinator's to rebuild | done |

11. **The external review of 2026-09-28, adopted as the next order of work.** Nine items, R1–R9,
   each closed: done, or for R8 closed on its pretest and reopened as item 12's S3, which the pilot
   then discontinued. The full item is ROADMAP-HISTORY.md lines 4388–4458.

   | # | Item (review's priority) | Code it names | Acceptance (the review's) | Lane, order |
   | --- | --- | --- | --- | --- |
   | R1 | Context in executable rules (P1) | `attribution/compile.py:173`, `measured.py:970` | a K562-specific rule is inactive in HepG2; conflicting results from two cell types survive compilation and round-trip serialisation | after lane-split (same files) |
   | ↳ R1 | **Done 2026-09-28 (lane-context, `8dcdade` registered, `6a4a8a9` build, `25f81c6` chr21)** | `compile.py`, `measured.py`, `ir/model.py`, `lang/grammar.py` | detail: ROADMAP-HISTORY.md line 4414 | done |
   | ↳ R1 | **Rule conditions fixed 2026-09-28 (lane-when, `49aca5d` census, `1d8713a` fix)** | `ir/model.py` | detail: ROADMAP-HISTORY.md line 4415 | done |
   | ↳ R1 | **`!=` refused 2026-09-28 (lane-ne, `ff0d767`)** | `lang/parser.py` | detail: ROADMAP-HISTORY.md line 4416 | done |
   | R2 | CRISPRi negatives read as outcomes (P1) | `measured.py:362` | significant decrease, significant increase, non-significant, missing and invalid kept apart with units, power, study and context; a significant increase never becomes "no effect"; a low-powered null cannot strongly reject; an increase does not establish a silencer | taken into lane-split (it holds `measured.py` and already carries power) |
   | R3 | A validated bridge from annotation to dynamics (P1) | `compile.py:295`, `runtime/grn.py:109` | an end-to-end fixture shows a supported expression response; missing regulator state gives an explicit unresolved-model diagnostic; another citation for the same mechanism does not change its simulated strength; predicted and measured forms of one relation are not counted twice | after R1 |
   | ↳ R3 | **Done 2026-09-28 (lane-bridge, `06b7557` registered, `b986238` amended for A8, `eda2728` build)** | `runtime/grn.py`, `attribution/bridge.py`, `compile.py` | detail: ROADMAP-HISTORY.md line 4419 | done |
   | ↳ R3 | **Activator rule decided on evidence 2026-09-28 (lane-combine, `a6fc5c4` census, `7f3e597` registered, `93caf61` result): mean retained, undetermined** | `runtime/grn.py` (`combine_activators`) | detail: ROADMAP-HISTORY.md line 4420 | done |
   | ↳ R3 | **Rates declared 2026-09-28 (lane-rates, `ff2062a` registered, `e921fae` result)** | `attribution/bridge.py`, `runtime/grn.py` | detail: ROADMAP-HISTORY.md line 4421 | done |
   | ↳ R3 | **No human rate to borrow 2026-09-28 (lane-rates2, `4c3a225` registered, `8a0ef9d` result)** | `attribution/bridge.py`, `build_gene_rates.py` | **The only absolute per-cell rate published is mouse** (Shao 2022, mouse ESC, "cell⁻¹ min⁻¹", not distributed per gene). (detail: ROADMAP-HISTORY.md line 4422) | done |
   | ↳ R3 | *Correction, 2026-09-29, lane-s5: the "1.56 h against 2.73 h" quoted in the two rows above is Schofield's own pair of cells (human K562, mouse MEF), not the half-lives the rate table carries, which are Schwanhäusser's NIH 3T3 values (median 9.96 h). Per gene, the K562 half-life is a median 0.23 times the table's value, and **most of that gap is measurement method, not species**. The factor the rows attribute to species is therefore mostly a method difference.* | — | — |
   | R4 | Confidence, effect size and evidence quality separated (P1) | `compile.py:271`, `organism/forge.py:170`, BIOFORGE-CONFIDENCE.md | effect magnitude alone never raises certainty; every probability names its outcome, calibration population and method, and an unavailable probability stays unavailable; the calibration text corrected (a constant can be calibrated to a base rate and still not discriminate); tests that only preserve the old constant updated | after lane-h releases `forge.py`; compiler half after R1 |
   | ↳ R4a | **BioForge half done 2026-09-28 (lane-conf, `38087a1` tests, `16a0434` build, `3b29509` docs)** | `genomeos/certainty.py`, `organism/forge.py` | detail: ROADMAP-HISTORY.md line 4425 | R4a done |
   | ↳ R4b | **Attribution half outside compile.py done 2026-09-28 (lane-conf2, `ef31028` census and strict-xfail tests, `3c5e122` build)** | `human_panel.py`, `variation.py`, `budget.py` | detail: ROADMAP-HISTORY.md line 4426 | done outside compile.py |
   | ↳ R4c | **Compiler half done 2026-09-28 (lane-bridge, `96bc5e5`)** | `attribution/compile.py` | detail: ROADMAP-HISTORY.md line 4427 | **R4 done** |
   | ↳ R4d | **Follow-ups done 2026-09-28 (lane-evid, `8d57243` calibration, `262ea4e` explorer)** | `evidence.py`, `ir/model.py`, `lang/parser.py` (two lines), `attribution/confidence_calibration.py` | detail: ROADMAP-HISTORY.md line 4428 | done |
   | ↳ R4e | **Leftovers done 2026-09-28 (lane-evid2, `a378182` census, `2085279` links, `37172b6` weak counts)** | `predict/enhancer_target.py`, two census scripts | detail: ROADMAP-HISTORY.md line 4429 | done |
   | ↳ R4f | **Last leftovers done 2026-09-28 (lane-conf3, `1219917` census, `7485df0` compiler, `6e833c3` adapter and missense)** | `compile.py` (`element_certainty` only), `predict/alphagenome_adapter.py`, `genome/missense.py` | detail: ROADMAP-HISTORY.md line 4430 | **R4 done** |
   | R5 | Evaluation independence through provenance (P1) | `measured.py:199` | immutable study, assay, source-file and partition identifiers on every record; validation outcomes kept out of features, candidate selection, starting labels and search objectives; automated overlap checks including related intervals and variants; reused benchmarks reported as reused | lane-split, running (its audit is the first half) |
   | ↳ R5, R2 | **done 2026-09-28 (lane-split: `0af7e1b` audit, `36fa298` split and provenance, `1cff24a` outcomes)** | `attribution/measured.py`, `compile.py` | detail: ROADMAP-HISTORY.md line 4432 | done |
   | R6 | Assay observations kept before aggregation (P2) | `measured.py:390` | one strong reporter tile cannot set a region's label; maxima only as labelled descriptive statistics; conflicting observations inspectable; lentiMPRA described as integrated, not episomal | after R2 |
   | ↳ R6 | **Done 2026-09-28 (lane-assay, `e0c911a` census, `bea3478` registered, `180ad41` + `f063eff` build)** | `measured.py`, `mpra.py` | detail: ROADMAP-HISTORY.md line 4434 | done |
   | R7 | A compositional ontology (P2) | `compile.py:283` | origin, role, activity, target and evidence status separated; a repeat-derived regulatory element keeps both; a predicted increase does not force a silencer label; alternatives survive serialisation; constraint read as selection, not mechanism | after R3 |
   | ↳ R7 | **Done 2026-09-28 (lane-ontology, `2ba726b` census, `8812e1a` registered, `329c2f9` build, `d2ab0ae` header fix)** | `compile.py`, `lang/grammar.py`, `lang/parser.py` | detail: ROADMAP-HISTORY.md line 4436 | **R7 done** |
   | ↳ R7 | **Follow-up done 2026-09-28 (lane-budget, `b6e6d17` census, `bbce3d8` registered, `308484a` build)** | `genome/unknown.py`, `attribution/budget.py` | detail: ROADMAP-HISTORY.md line 4437 | done |
   | ↳ R7 | **Consumers moved 2026-09-28 (lane-consumers, `775d4be` registered, `690a71c` build)** | `budget.read_axes` and 18 direct readers | detail: ROADMAP-HISTORY.md line 4438 | done |
   | R8 | The joint inference engine, bounded (P2) | `attribution/closure.py`; genomeos-8a's brief | pretest first (coupling against independent scoring on validation folds; a negative closes it); then fixed-boundary search solving a synthetic two-change case, and gains on independent evidence over unchanged annotation and independent-block baselines; convergence never reported as validation | after R1–R5 |
   | ↳ indep | **Independent benchmark probed 2026-09-28 (lane-indep, `042a4aa` registered, `0924a5c` result): NOT READABLE** | `attribution/indep.py` | Coverage was fine (68% within reach), but **only 13 covered distal decreases against a floor of 20**: 74 of the screen's 93 decreases are promoters, and its distal hits are 181 increases against 19 decreases, which anyone reusing it needs the authors' reading of. (detail: ROADMAP-HISTORY.md line 4440) | open |
   | ↳ indep | **Second probe 2026-09-28 (lane-indep2, `12ad319`): nothing scorable at 0 requests** | `scripts/indep_igvf_probe.py` | detail: ROADMAP-HISTORY.md line 4441 | open |
   | ↳ indep | **HCT116 arm run 2026-09-28 (lane-hct116, `6257b95` protocol, `a39073d` result; the 705 requests approved by the owner): passes, NOT replicated** | `scripts/crispri_hct116.py` | It passes the registered rule, but the interval includes zero, so **the CRISPRi result stays replicated in one cell type**. (detail: ROADMAP-HISTORY.md line 4442) | done |
   | ↳ 1.3 | **Clause 2 re-tested against a matched control 2026-09-28 (lane-clause2, `e1dbcf3` registered, `0c8b82d` result): still not met** | `scripts/clause2_matched_control.py` | detail: ROADMAP-HISTORY.md line 4443 | clause 2 not met |
   | ↳ R8 | **Pretest done 2026-09-28 (lane-joint, `dcc277b` registered, `4098753` result): NEGATIVE, R8 closes** | `attribution/joint_pretest.py` | **No joint search engine is built**; the synthetic two-change case and the search are cancelled. (detail: ROADMAP-HISTORY.md line 4444) | closed |
   | ↳ R8 | *Correction, 2026-09-28, second review (item 12, S3): the row above overstates its negative. The pretest rejected **two particular score transformations** (element competition and gene budget) on the available K562 pairs; it did not show that joint inference cannot help. "No joint search engine is built" and "cancelled" are withdrawn as conclusions; the two transformations are retired, the negative stands for them, and a bounded pilot reopens the objective.* | — | reopened as S3 |
   | R9 | Result manifests and public claims (P2) | `genomeos/results.py:23`, README | every result carries source accessions and versions, checksums, assembly, coordinate convention, code revision, parameters, exclusions and partitions; a second environment rebuilds one representative result or names what is missing; public summaries separate sequence, annotation and assay coverage, software correctness and independent prediction | now; README corrected 2026-09-28 by the coordinator |
   | ↳ R9 | **manifests done 2026-09-28 (lane-manifest, `fe0880a`, `4b8e8b9`; evidence fix `c22c8c6`)** | `genomeos/manifest.py`, under `result_manifest` | detail: ROADMAP-HISTORY.md line 4447 | done, one writer |
   | ↳ R9 | **Headlines rebuilt 2026-09-28 (lane-manifest2, `42e786f`, `190a144`, `65d2d23`, `f48b909`)** | the five README headline writers | detail: ROADMAP-HISTORY.md line 4448 | done, five of six |
   | ↳ R9 | **Sixth headline rebuilt 2026-09-28 (lane-manifest3, `0126346`, `86e3c52`)** | `crispri_published`, `manifest.py`, `manifest_rebuild.py` | detail: ROADMAP-HISTORY.md line 4449 | done, six of six |
   | ↳ R9 | **Model dependency recorded 2026-09-28 (lane-pin, `837d5da`)** | `manifest.py`, five headline writers | detail: ROADMAP-HISTORY.md line 4450 | done |
   | ↳ R9 | **Five manifests re-recorded 2026-10-02 (lane-rerecord; `8f7d215`, the writer only)** | the five results whose manifests predate `members` | detail: ROADMAP-HISTORY.md line 4451 | done, nothing published |
   | ↳ R9 | **Model version pinned 2026-09-28 (lane-pin2, `981b5b0`)** | `predict/alphagenome_adapter.py`, `enhancer_target.py`, `manifest.py` | detail: ROADMAP-HISTORY.md line 4452 | done |

12. **The second external review, adopted 2026-09-28: the original objective is unfinished, and the
   status measures used so far do not measure it.** Eight items, S1–S8: S1, S2, S4, S5, S6 and S7
   done, S3 discontinued with the pilot, S8 built with its floor values left to the owner. The full
   item is ROADMAP-HISTORY.md lines 4459–4643.

   | # | Item | Code it names | Acceptance (the review's) | Order |
   | --- | --- | --- | --- | --- |
   | S1 | Uncertainty survives serialisation | `ir/model.py:43` | unavailable confidence is explicit and preserved through parsing, JSON, exports and summaries; the test that accepts the loss is replaced by a preservation test | wave 1 |
   | ↳ S1 | **Done 2026-09-28 (lane-s1, `1307fee` census, `92654f5` build)** | `ir/model.py`, `lang/parser.py`, `lang/tools.py`, `bio.py`, `runtime/{uncertainty,located,debugger,methylation}.py`, `cli.py`, `web/server.py`, `index.html`, `organism/{experiment,forge}.py` | detail: ROADMAP-HISTORY.md line 4505 | done |
   | S2 | Calibrate the power model before any sample size is used | `scripts/clause2_design_power.py:460` | the simulation reproduces its own anchor's quantity (the review reports that it takes a 13.45% window-level anchor and produces 16.80% element detection and 56.96% positive windows; reproduce first); sensitivity and false positives used consistently; observed and latent correlation distinguished; the planned analysis simulated under the null and alternatives, with shared controls and clustering; designs compared at equal total assay cost; sizes above the eligible block population flagged. **Until then 20–5,000 is exploratory, not an experimental specification** | wave 1 |
   | ↳ S2 | **Done 2026-09-28 (lane-s2, `b903e1e` reproduction, `019a9ab` registered, `cab70d9` result)** | `scripts/clause2_design_power.py`, `clause2_design_power_calibrated.json` | detail: ROADMAP-HISTORY.md line 4507 | done |
   | ↳ S2 | **Follow-ups done 2026-09-29 (lane-rule, `1cb7559` registered, `1ef8990` build, `d85fc5d` rebuild, `4f094bf` result)** | `scripts/clause2_design_power.py`, `scripts/clause2_measured_arm.py`, `clause2_reading_rule.json` | detail: ROADMAP-HISTORY.md line 4508 | done |
   | S6 | The evidence contract enforced | `results.py:52`, `manifest.py:83` | a failed new result goes to a separate location and only validated results enter the registry; an explicit legacy allowlist; untracked executable source counted in the revision stamp | wave 1 |
   | ↳ S6 | **Done 2026-09-28 (lane-s6, `0e50036` census, `e1464b2` registered, `f0dbfd0` build)** | `results.py:save_result`, `manifest.py:code_revision`, `data/results_legacy.txt` | detail: ROADMAP-HISTORY.md line 4510 | done |
   | ↳ S6 | **Every writer under the contract, 2026-09-29 (lane-contract, `d6a77dd` census, `07caa6c` clause 2 corrections, `8934a72` build)** | `scripts/results_writer_census.py`, `tests/test_results_writers_guard.py`, `organism/provenance.py` | detail: ROADMAP-HISTORY.md line 4511 | done |
   | ↳ S6 | **Two worm results reproduce again, 2026-09-29 (lane-repro, `b587bb1` code, `236fa2b` results)** | `organism/commitment.py`, `scripts/celegans_fate_reads.py`, `tests/test_celegans_determinism.py`, `tests/test_fate_rules.py` | detail: ROADMAP-HISTORY.md line 4512 | done |
   | ↳ S6 | **Backup gap found 2026-10-01 (coordinator, read-only)** | git-ignored data | detail: ROADMAP-HISTORY.md line 4513 | open (owner's choice) |
   | ↳ S6 | **Backup made and verified 2026-10-01 12:16 (coordinator, on the owner's authorisation)** | `/Volumes/Albert Canfield/GenomeOs/GenomeOS-snapshot-2026-10-01T1216` | detail: ROADMAP-HISTORY.md line 4514 | done |
   | S7 | Release checks follow the revision being promoted | `.github/workflows/ci.yml:8` | CI runs on later pushes to a pull request; the scheduled check exercises the branch being worked on; checks are required on the revision promoted to `main`; local hooks stay, but not as the only enforcement | wave 1 |
   | ↳ S7 | **Done 2026-09-28 (lane-s7, `b50c873`)** | `.github/workflows/ci.yml`, `scripts/promote_main.sh` | detail: ROADMAP-HISTORY.md line 4516 | done |
   | ↳ S7 | **Gate found insufficient 2026-09-29 (coordinator)** | `scripts/promote_main.sh`, `CONTRIBUTING.md`, the checkout guard | detail: ROADMAP-HISTORY.md line 4517 | open |
   | ↳ S7 | **Cause narrowed 2026-09-29 (the external reviewer's correction, checked against the source by the coordinator)** | — | detail: ROADMAP-HISTORY.md line 4518 | open |
   | ↳ S7 | **Gate made honest 2026-09-30 (lane-gate, night shift; `3a893fc`)** | `scripts/promote_main.sh`, `scripts/promote.sh`, CONTRIBUTING.md | detail: ROADMAP-HISTORY.md line 4519 | done (clean-up `8df8799`, 2026-10-02) |
   | ↳ S7 | **Promoted 2026-10-02: `main` = `a1ce607` (PR #13, opened by the supervisor under Albert's delegation of 2026-10-01, merged by Albert; verified independently by the coordinator)** | — | detail: ROADMAP-HISTORY.md line 4520 | done |
   | ↳ S7 | **Qualifying event, not eligibility 2026-09-30 (lane-gate; `03497e7`, additions only)** | `scripts/promote_main.sh`, CONTRIBUTING.md, tests | detail: ROADMAP-HISTORY.md line 4521 | done (`8df8799`, 2026-10-02) |
   | ↳ S7 | **Gate finding, 2026-09-30 (narrowed by the external reviewer)** | `scripts/promote_main.sh` | detail: ROADMAP-HISTORY.md line 4522 | open (finding) |
   | ↳ S7 | **Promoted 2026-09-30: main = `f7016ff` (PR #12, merged 19:33 UTC)** | — | detail: ROADMAP-HISTORY.md line 4523 | done |
   | ↳ S7 | **Promoted 2026-09-29: main = `7ff4e37` (PR #11, merged by the owner)** | — | detail: ROADMAP-HISTORY.md line 4524 | done |
   | S3 | Reopen the coherence objective with one bounded pilot | `attribution/joint_pretest.py:70` | R8's negative is kept and its two transformations retired, but its code and this plan stop claiming more than it showed; a small model with explicit roles, target links and cell context solves a synthetic case that needs two simultaneous corrections, then is scored on untouched evidence against unchanged labels and independent-block scoring; it returns ranked alternatives, supporting and conflicting evidence, unresolved alternatives and the next measurement that would separate them; **success is better prediction on untouched evidence at useful coverage, never a higher internal coherence score**; expansion stops if it adds nothing | wave 2, after wave 1 |
   | ↳ S3 | **Discontinued 2026-09-29 (lane-pilot; `a9660ff` registered, `1a492c8` gate 1, `a89bf5c` gate 2 registered, `197c560` result, `6bbf1d6` + `40b66e9` merge recount)** | `attribution/pilot.py`, `attribution/pilot_bio.py`, `attribution/joint_pretest.py` (R8's reading corrected beside the original; element competition and gene budget RETIRED) | **Gate 2 failed, `beats_unchanged_only`:** on chr1/6/8/9/10/11/19 the pilot beat the unchanged labels (Gasperini2019 +0.165 [0.074, 0.266], Morris +0.431, Xie +0.277, also at their coverage) but **passed 0 of 4 CRISPRi decrease endpoints**; minus independent-block, Gasperini2019 **−0.065 [−0.125, −0.012]**, the rest span 0; minus distance, above 0 on Xie only (+0.105 [0.015, 0.261]). (detail: ROADMAP-HISTORY.md line 4526) | **discontinued investigation (stop rule fired)** |
   | S4 | What a correct attribution means | — | origin, function, activity, target and context kept separate; not naming a coding gene is not absence of function; reporter activity, perturbation, contact and conservation never substitute for one another; "the observation model is inadequate" kept as an explanation; target accuracy reported apart from role accuracy and from coverage | wave 3 |
   | ↳ S4 | **Done 2026-09-29 (lane-s4, `b87f59e` registered, `ab62d99` result)** | `attribution/correctness.py`, `scripts/s4_correctness_run.py` | detail: ROADMAP-HISTORY.md line 4528 | done |
   | ↳ S4 | **Traced 2026-09-29 (lane-repress; `5d14d37` classes registered, `54c47f1` code, `0762073` tests, `3cc59bf` doc, `641909e` result)** | `scripts/repression_trace.py` | detail: ROADMAP-HISTORY.md line 4529 | done |
   | ↳ S4 | **Judge v2 done 2026-09-29 (lane-judge2; `5931609` registered, `172572c` code, `4567219` result, `8f3689e` doc)** | `attribution/correctness.py` (`RULE_V1`, `RULE_V2`), `scripts/s4_correctness_run.py --rule v2` | detail: ROADMAP-HISTORY.md line 4530 | done |
   | ↳ S4 | **The owner's rulings, 2026-09-29** | — | detail: ROADMAP-HISTORY.md line 4531 | running |
   | ↳ S4 | **Judge v3 done 2026-09-29 (lane-judge3; `2c95f66` registered, `b7cb814` code, `d1e09ae` result, `d9c9ea9` doc)** | `attribution/correctness.py` (`RULE_V3`, now the default), `scripts/s4_correctness_run.py --rule v3` | detail: ROADMAP-HISTORY.md line 4532 | done |
   | S5 | Runnable simulation separated from validated human kinetics | `attribution/bridge.py:174` | 22,580 of 440,589 rules (5.1%) simulable on mouse NIH 3T3 rates supports assumption-dependent simulation only; transferred rates labelled; parameter ranges propagated and the predictions stable across them identified; one human cell context with matched measurements preferred over more rules run on defaults | wave 3 |
   | ↳ S5 | **Done 2026-09-29 (lane-s5, `0ea0a90` registered, `83013e6` label and pre-run amendment, `a54ea12` result, `74e6a2f` audit)** | `attribution/bridge.py`, `runtime/grn.py`, `scripts/s5_rate_ranges.py`, `scripts/bridge_audit.py` | detail: ROADMAP-HISTORY.md line 4534 | done |
   | S8 | The next experiment's information value | — | a representative arm (how common a function is) and a disagreement-selected arm (which model is right) kept apart, selection probabilities recorded; independent loci, assay feasibility and distinguishable hypotheses before volume; HCT116 stays inconclusive, and re-reading held-out data is not fresh validation | wave 3 |
   | ↳ S8 | **Built 2026-10-01, pushed `bee8f20` (lane-s8): the gate that asks first, checked by the coordinator against the files** | `genomeos/attribution/experiment_value.py`, `tests/test_experiment_value.py` | detail: ROADMAP-HISTORY.md line 4536 | done |
   | ↳ S8 | **Four things S8 deliberately left to a decision, 2026-10-01 (lane-s8, recorded by the coordinator; none chosen by the lane)** | (1) **The floor values.** `Floors` ships no defaults: a real design must fix `min_units`, `min_positive_units`, `min_label_prevalence`, `min_independent_loci`, `min_distinguishing_candidates`, `target_loci` and `units_per_locus` before any count. The lane invented no number, which is right; **the numbers are the owner's**. (2) **One shared locus constant.** `experiment_value.py` and `attribution/cell2.py` registered the same 1 Mb chained convention on the same day (`7e45429`), independently and with a comment in each saying the two must not diverge. **Coordinator's decision: they become one constant, after lane-cell2 finishes**, because that module is live and a shared file must not be edited under a working lane; until then the comment is the guard. (3) **Where an incomplete design fails.** An undeclared readout comparability or an unstated detection limit refuses the whole design, while an unavailable assay drops only the candidate; the lane's reason is in the module and the alternative is a two-line change plus a test. (4) **Volume is an allocation, not power.** A detectable-effect calculation would need a variance assumption registered per assay, with data and money implications, so the lane did not start it. Also left undone: no design is instantiated — no script, CLI, API surface or result file — and `independent_loci` treats positions as points, not intervals, so wide elements need checking before they are passed to it | detail: ROADMAP-HISTORY.md line 4537 | owner decides (1); coordinator holds (2) |

13. **The coherence programme, 2026-09-28: six approaches, organised as the design of item 12's pilot (S3)
   and what follows it.** Brought by the owner the same day as item 12. The organising idea: stop asking
   "is this block X or Y" one block at a time, and ask which **smallest set of assumptions cannot explain
   the observations together**, which repairs remove the contradiction, and which repairs the evidence
   cannot tell apart. **The metric is fixed before anything is built:** validated corrections per
   compute-hour, with errors and abstentions reported beside it, where a correction counts as validated
   only if it improves prediction of evidence it did not see. Internal coherence is never the score.
   Scale is decided on prediction, not on how many consistent stories are found.

   | # | Approach | What it does | Depends on | Phase, lane |
   | --- | --- | --- | --- | --- |
   | ↳ C1–C3 | **Discontinued 2026-09-29 (lane-pilot)** | one pilot, as item 12's S3 row records | The stop rule fired (gate 2 `beats_unchanged_only`): **the joint debugger never beat its own blocks revised alone, and beat distance to the gene on 1 of 4 CRISPRi decrease endpoints**; 6 validated and 6 wrong of 261 committed corrections, 118.34 of each per CPU-hour. Recorded as a discontinued investigation, apart from implemented capabilities. **Phase C (C6, C5's build, S8's design) does not follow from it** | — | discontinued |

   **Sequence.** Phase A runs now beside item 12's wave 1. Phase B is one pilot combining C1, C2 and C3 on
   several chromosomes, with chromosomes as validation partitions and supported connections kept across
   neighbourhood boundaries; it starts when C4's harness exists and S1 and S6 have landed. Its first gate
   is item 12's: a synthetic case that needs two simultaneous corrections. Its second is the metric
   above on C4's withheld evidence, against unchanged labels and independent-block scoring. **If the
   pilot does not improve prediction, it stops and is recorded as a discontinued investigation**, apart
   from implemented capabilities. Phase C (C6, C5's build, S8's experiment design) follows only a pilot
   that helped. Ordinary CPUs and cached data throughout; no model requests unless the owner approves.

The rest of item 13 is ROADMAP-HISTORY.md lines 4652–4784.

---

## 6. Milestones

Each row as it stands on 2026-10-10. A Proof cell longer than a paragraph is cut to its verdict,
in the row's own words, and a pointer; the full rows and the notes that stood below them are in
ROADMAP-HISTORY.md §6 (lines 4923–5176).

| Version | Milestone | Proof |
|---|---|---|
| 0.1–0.8 | engine, compiler from public data, cell runtime, composition spike, spatial, whole worm, twin, debugger | done, see PROGRESS.md |
| **0.9 whole genome** ✅ | every chromosome fetched, classified with curated repeats, domains found, proteome compiled and verified, HG002 twin genome-wide, graph genome-wide | all reached 2026-09-11; the proteome also ships as a packaged offline library. Version is 0.9.0 and the README states what that means; the git tag is cut when the pull request to `main` is merged |
| **1.0 experiments** ✅ | `experiment` block; C. elegans mutants reproduced; three published perturbations as tests; BioForge takes experiments as input; Evidence explorer | **Reached, and the row said "planned" until 2026-09-22 while every clause was already met.** `bio test` ran 41 programs at 220/220 checks in CI on that date (45 at 246/246 since `5682e6c`, 2026-10-10), including ten worm mutants and nine blood mutants as `experiment` programs. **One caveat the version does not carry: BioForge's own emitted confidence is a hand-set `0.3` and is `UNCHECKABLE_BY_CONSTRUCTION` (docs/BIOFORGE-CONFIDENCE.md), so "takes experiments as input" is true and "prices its answer" is not.** *(2026-09-28: since `16a0434` an answer no longer carries the fixed 0.3, its probability stays unavailable with its reason, and `3b29509` corrected "could never be calibrated" to "unsupported"; "prices its answer" is still not true. Area H.)* Detail: ROADMAP-HISTORY.md line 4929 |
| **1.1 human mechanism** ✅ | haematopoiesis as a mechanism module inside the human body program; reader v1 (open nodes per cell type) | **Reached 2026-09-27 (`63929e1` registered, `b13143e` scored), on a criterion fixed before the run.** Haematopoiesis landed 2026-09-11 with its mutants in CI, and reader v1 runs on thirteen biosamples across every chromosome. Three per-biosample readings each ship a normalised counterpart beside the raw count, and all four clear the registered bars. **"Open nodes per cell type" has a between-cell-type reading for the first time**: nodes at one fixed 4.79 peaks per 100 kb, residualised on DNase depth. **Two caveats travel with the version.** The flagship "K562 reads 74%, HepG2 65%" is about 70% DNase depth: at matched depth K562 still reads more, by **2.7 points rather than 9.0**; and the residuals remove DNase and H3K27me3 depth, not FRiP, fragment length or library complexity, none of which is on disk. Detail: ROADMAP-HISTORY.md line 4930 |
| **1.2 therapeutics benchmark** ✅ | approved targets recovered from public tumours; CNA and SV; `cancer.*` libraries | **Reached 2026-09-27 (`bdd5058`).** Clause by clause. *Approved targets recovered from public tumours*: nine cases, **9/9 recovered, 9/9 routes called correctly, 9/9 top mechanisms defensible**, and since 2026-09-27 every target also ranks first among the candidates measured in its own tumour (CD19 second, behind a candidate from the same RNA call), with no case where a hypothesis measured nowhere outranks it (`outranked_by_hypotheses` 12 rows → 0, `BURIED_SURFACE_TARGETS` 2 → 0, `AMPLIFIED_TARGET_RANK` 4 → 1, every score byte-identical). *CNA and SV*: scored cases driven by an amplification (`3aee02f`) and a fusion (`deee0a2`). *`cancer.*` libraries*: landed 2026-09-14 (`9c0e6ce`). One precision: the evidence tier passed on an **amended** design; the registered one (`5013406`) fired its own falsifier and is left in the record. The benchmark holds ten cases at 10/10/10 since 2026-09-30 (area F). Detail: ROADMAP-HISTORY.md line 4931 |
| **1.3 the 98%** ◑ | every UNKNOWN block with a tier and a confidence; the constrained-unknown blocks attributed to a gene and a tissue; the attributions scored against measured elements. **The third clause cannot be met as written and the milestone says so since 2026-09-17: only 0.45% of the unknown space has ever been measured by any assay this project holds, so "scored against measured elements" can be satisfied for a sliver and not for the space. The exit criterion is therefore split: the attribution and its scoring where measurement exists, and a designed experiment for the rest, which is built (**284,001 oligos** over 878 blocks, 41,037 of the test arm from blocks no assay has touched; this said 312,129 until 2026-09-27, the first build's count before the orderability filter dropped 4,264 untouched windows, and 45,301 − 4,264 = 41,037 exactly) and unrun.** | the genome budgeted 2026-09-12: every block tiered with a confidence of 0.5 or above. The attribution read on 2026-09-16 was withdrawn on 2026-09-17: standardised on length, GC and promoter distance, no tier differs from the neutral tier. Clause 2 is held not met (the notes below); clause 3 is split as the goal cell says, its experiment designed and not ordered (D8). Detail: ROADMAP-HISTORY.md line 4932 |
| **2.0 BioLang standalone** ✅ | `biolang` package: lang, ir, runtime, std, `bio`; GenomeOS depends on it | a `.bio` program runs with GenomeOS uninstalled; two test suites. **Reached 2026-09-27 (`5720c7c`), on the second clause as restated that morning in section 5 item 10, before the work began:** "GenomeOS depends on it" was a means, and is read as "the engine is separately installable and independently tested", because the repository is deliberately not split while nothing outside it consumes the engine. The package is `biolang` **0.1.0**, **9 of 9 checks** in an interpreter where `import genomeos` fails, the ninth being its own pytest suite of **24 files and 157 tests — 154 passed, 3 skipped**. *(The committed `engine_package.json` of 2026-09-29 reads 12 of 12 checks; area A.)* **Not done, and not a clause: publication to PyPI**, which is outward-facing and the owner's. Detail: ROADMAP-HISTORY.md line 4933 |
| **1.4 the decoded genome runs** (proposed 2026-10-10) | the compiled programme runs in the BioVM per cell type with its unresolved parts reported; chromatin gates its rules through the language (`reader`), generated from the decoded data; compiled rules carry the regulator that binds each element, so a run has network effects; one registered test at gene level and one at cell level of what a run predicts, against measurement the model never saw, each reported in its registered words whatever the outcome | proposed by the review of 2026-10-10, for the owner to adopt or change. Clause 1 holds in part since 2026-10-10 (M1, `747b207`: every resolved pair runs, as one-gene systems); clause 2's engine half landed (`5682e6c`) and its data half too (`98a1bcb`), with the gated run registered as "worse" against measured RNA (`ff35b09`); clause 3 is §5 item 7b; clause 4's first two tests were run (M1 `747b207` and M3 `80bf3db`, both "no difference detected"), and a test whose run has network effects waits on clause 3 |
| **1.5 human mechanism beyond blood** (proposed 2026-10-10) | a second human mechanism module in the body program, written from a closed sourced list, its published mutants as `experiment` blocks in CI, every fitted number named as fitted | proposed by the review of 2026-10-10, for the owner to adopt or change; §5 item 8 (M4, the small-intestinal crypt: registered at `460c46f`, module `7646bd7`, result `ad8c662`, eight of eight registered knockouts reproduced with one number fitted and named) |

**1.3's clauses, condensed** (the notes of 2026-09-27 to 2026-09-29 are ROADMAP-HISTORY.md lines
4935–5110). Clause 1 is met. Clause 2, against four matched controls (`0c8b82d`, `faeb0da`, `a7f207f`,
`e8aa571`): 9.18% of the 3,280 elements inside the real unknown name a coding gene against 44.04% of the
456,573 in their windows, a difference in target-naming frequency by the model, not in prediction
accuracy; the measured arm cannot decide (59 of 882 tier blocks hold a qualifying measured element); the
experiment's size is exploratory until the model is calibrated (`cab70d9`, `4f094bf`).

**The strongest defensible conclusion:** target naming is substantially less frequent in the analysed
constrained-unknown elements, and the difference persists under several adjustments. The available
measurements cannot distinguish biological differences from model limitations. Experimental
requirements remain provisional. Clause 2 stays not met.

The other notes that stood below this table (the version brought to 1.0.0 on 2026-09-27 with no git
tag, and 1.2's moves to 8/8 and 9/9 on 2026-09-21) are in ROADMAP-HISTORY.md, lines 5118–5173.

---

## 7. Improvement and creativity pool

Ideas judged worth keeping, not yet scheduled. Each would be judged by the
ambition in §1 before it is started.

- **Cell2Sentence, as a comparator and never as evidence.** C2S-Scale (van
  Dijk lab; Apache 2.0; Pythia 160M-1B and Gemma-2 2B/27B on HuggingFace since
  2025-10-15; trained on 57M+ CellxGene and HCA cells) turns expression
  profiles into ranked gene sentences for LLMs: annotation, generation,
  perturbation prediction. **It reads expression, not DNA, so it cannot
  attribute function to a non-coding element.** *Trigger:* only when a
  perturbation-response test (N1 type) clears its label-coverage floors. Then
  C2S-Scale enters as a **pre-registered comparator, never as evidence**.
  *Conditions:* (1) a **contamination check first** -- is any benchmark
  dataset, or a study it derives from, in CellxGene or HCA? If yes, it is
  development evidence only; (2) **2B or smaller locally** (27B does not fit
  18 GiB RAM); no 27B without Albert's approval, since it means paid compute;
  (3) the same floors, bootstrap and reading rules as every comparator.
  Logged 2026-10-02 by the supervisor, on Albert's question; **no lane opened,
  and the figures in this row are the supervisor's and are not verified
  here.**

- **Programs that read the genome directly.** A `gene` block whose sequence
  is taken from a locus at compile time, so a BioLang program can be
  compiled against HG002 and against the reference and diffed.
- **The reader as a first-class object.** `reader cell_type X { open: [domains] }`
  distilled from ENCODE DNase, so context gating comes from chromatin, not
  from a hand list of expressed genes.
- **Writers.** Replication errors (somatic mutation rate per division),
  epigenetic drift (clock CpGs), germline recombination; each a `writer`
  block with its rate and evidence, run by the Body runtime.
- **Genome diff as code review.** Two individuals, or an individual against
  the reference, as a diff of BioIR modules with consequences and evidence,
  not a VCF. Done on 2026-09-12 at both levels (`genomeos individual diff`
  and `--regulatory`, area D); the effect per regulatory variant waits on the
  self-hosted model (area I).
- **Time-lapse view.** The Organism tab replayed with the block map: which
  nodes are being read at each stage.
- **`bio` package registry.** `import bio.std.*` from a versioned registry
  with licences, so libraries such as `human.haematopoiesis` can be shared.
- **Question-driven UI.** A single box ("what does APP do in the brain", "what
  breaks if TP53 is lost") routed to the commands that answer it, each answer
  carrying its evidence pills.
- **Uncertainty budgets.** A program declares the confidence it needs
  (`require: confidence >= 0.6 at cellular`) and `bio test` fails when the
  evidence does not reach it.
- **Cross-species compilation.** The same `organism` program compiled
  against worm and mouse genomes where orthologues exist, reporting what has
  no counterpart.
- **Agent packet as a stable API.** The cancer and therapeutic packets
  versioned and documented as the contract other tools consume.
- **GenomeOS as an agent skill.** DeepMind's `science-skills` collection
  (Apache 2.0 / CC-BY) packages one scientific database per skill so an agent
  can use it. GenomeOS already implements two thirds of that list natively and
  with evidence, so the value is the other direction: publish a skill that
  drives `genomeos` and `bio`, and let an agent reach every layer through one
  tool instead of thirty. Their list is also a source checklist: JASPAR
  (the promoter motif scan SEQUENCE-GRAMMAR.md already plans), UniBind,
  Foldseek and the EBI Ontology Lookup Service are sources we do not use yet.
- **The two-axis constraint map as a track.** Once `variation_chr*` exists
  genome-wide, a Blocks-tab lane that colours every block by its case
  (syntax, slot, recent, relaxed) beside the tier, so the "syntax versus
  value" reading is a picture of the chromosome, not a table.
- **Recurring operators as a BioLang library.** If step 4 of area J finds
  motif combinations that recur across a library's promoters, write them as
  `signal` and `rule` blocks in `bio.std.operators` with their enrichment as
  evidence, so a program can `import` the heart-field or blood operator the
  way it imports a pathway.
- **Genome diff across species as code review.** The same locus decompiled
  in human, mouse and zebrafish, diffed as BioIR: which elements, targets and
  rules are shared, which are lineage-specific; the first output of area J's
  step 7 in the form the "genome diff as code review" idea above proposes for
  individuals.
- **A self-hosted AlphaGenome.** The weights are published; self-deployment
  removes the per-request quota and turns enhancer-to-gene, the predicted
  reader and in-silico mutagenesis from per-chromosome jobs into genome-wide
  ones (docs/ALPHAGENOME.md). Since 2026-09-11 it is the prerequisite of area
  I's attribution step: a million blocks cannot go through the hosted quota.

---

## 8. How the work is organised (several sessions, one checkout)

- **At least 2 of every 3 lane slots are MODEL lanes** (raised from 2 of every 4 after Albert's
  direction of 2026-10-03 that git and verification had taken the project's focus; on 2026-10-10
  the measured share was about a third, §5) (correct, or complete) —
  work that adds or repairs something the model can be held to. A VERIFICATION
  lane takes a slot only when it blocks a model result or a promotion, or when
  it fixes a defect found in published figures. Tooling for its own sake waits
  in §7. Set 2026-10-02, after a night in which almost every lane was
  verification: a rebuild that passed without hashing the largest thing it read,
  17 published results reading an undeclared table, 10 tests passing on
  hand-installed packages, and a verdict that could read green about a tree that
  no longer existed. That work was worth doing and the proportion still was not
  right, so the ratio is written down rather than left to judgement.
- **The order of work from 2026-10-10 is the execution plan in §5** (items 1 to 17). As it stood on
  2026-10-02, kept for the record: AstroREG once its
  calibration and the owner's answer are in; then a registered test of the
  increase-derived links, with the links answerable IN THEIR OWN CELL counted
  before any direction is read (the direction test fell from 48 measured to 21
  answerable, and only noticed because it counted); then response-map
  increment 4 on the next assay its count names.
- All work on `dev`; Albert alone opens pull requests to `main`, by hand;
  no session or script opens or merges one. No attribution trailers.
- The cycle: finish the piece, its targeted tests and lint on its own files
  green, commit to `dev` with `scripts/commit_own.sh`, push once per finished
  piece with `scripts/push_own.sh`. The `pre-push` hook checks the pushed sha
  in a clean worktree and refuses a red push: since 2026-10-10 (`41b14b7`, §5
  item 10) a `dev` push takes the fast gate (lint, format, the always-run tests
  and the test files changed since the last green CI), and the full
  `scripts/check.sh` returns when the daily CI is red or unread, on `main` and
  promotion branches, and under `GENOMEOS_FULL_GATE=1`. A lane runs the full
  suite itself only for a paid path or a verification tool (since 2026-10-02).
  CONTRIBUTING.md has the commands.
- Each session owns files (listed in §3) and commits through a private
  index; shared files (`cli.py`, `server.py`, `index.html`, `README.md`,
  `PROGRESS.md`) take isolated hunks only, announced to the peers. The
  shared `.git/index` is stale and must not be committed from.
- Lint and format only your own files; the full suite must pass at the
  push (5,784 tests collected on 2026-10-10).
- Long jobs run once, through the registry, with their pid recorded; a
  result per chromosome is committed as it lands.
- PROGRESS.md is append-only and is what the Progress tab shows; this file
  is edited by the coordinating session named on the work board (session
  names change on restart), or, when none is live, by the session Albert
  asks; other sessions send finished rows. It is updated at each milestone.
  Documents are layered as docs/README.md describes: wide
  (README, ROADMAP, ARCHITECTURE, DECISIONS), specifications per language
  version, one design document per area, records (PROGRESS, LESSONS).
- docs/ROADMAP-HISTORY.md is the roadmap as it stood at `c00250f`, frozen on 2026-10-10 and
  read-only: nobody edits it, each area here ends with the lines there that hold its record, and a
  correction to anything in it is made additively in this file.
