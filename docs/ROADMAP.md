# GenomeOS: scope, architecture, areas, plan

Reviewed and reorganised on 2026-09-17 against the code, the tests, the
result files and the running jobs (first written 2026-09-11). This is the one document that says what
GenomeOS is, how it is built, what each area is for, what is done, what is
missing and what comes next. Evidence per finished task stays in
[PROGRESS.md](PROGRESS.md); the original phased plan is kept as history in
[ACTION-PLAN.md](ACTION-PLAN.md); lessons the code carries are in
[LESSONS.md](LESSONS.md).

Measured state on 2026-09-17: **536 commits**, 63 `genomeos` commands plus the
`bio` toolchain, **1,050 tests** collected and green with 212/212 `.bio` checks,
lint clean, every human chromosome fetched and analysed, the whole human proteome
compiled and verified against UniProt, a genome-wide knowledge graph, the whole
worm grown from one cell, a human body grown as populations, haematopoiesis as
mechanism, and the engine separable from the application under its own licence.

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
organism, with evidence and uncertainty on every fact. Open source, MIT
licence. The governing documents are indexed in [README.md](README.md) and
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
                                           genomeos/design      minimal-cell design
                                           genomeos/forge       BioForge search
                                           genomeos/web         server + one-page UI
                                           genomeos/jobs        background jobs registry
                                           genomeos/storage     stream / distil / clean
                                           genomeos/results     result summaries (data/results)
                                           genomeos/cli.py      42 commands
```

Boundary rules (from ARCHITECTURE.md §10) and where they are broken today:

| Rule | State |
|---|---|
| BioIR JSON is the contract; the language layer reads no GenomeOS data paths | holds |
| Engines take a `Module` and return a trajectory or a typed result | holds |
| Standard library written in BioLang, not Python | holds (ageing, cell_types, development) |
| One thin entry point (`bio`) | holds (`genomeos/bio.py`, 333 lines) |
| No dependencies in the language core | holds |
| Engines import nothing from genome, knowledge or web | **holds since 2026-09-11** (b162919 and b90e741): a grep over `lang`, `ir`, `runtime`, `coords.py`, `version.py` and `bio.py` on dev finds no import of the application. `bio` runs on an engine-side `lang/tools.py` instead of the CLI, and the version string moved to `genomeos/version.py` so the Apache half never reads the AGPL package root |
| `lang` and `ir` depend only on themselves | `Locus` lives in `genomeos/coords.py`; move it under `ir` when the packages split |

### 2.1 Where each component stands (2026-09-11)

"Functional" means a geneticist can use it today on real data with evidence
attached. "Missing for complete" is what separates it from the ambition.

| Component | Version | Functional today | Missing for complete |
|---|---|---|---|
| **BioLang** (language) | v0.3 | 18 block kinds: `module`, `import`, `gene`, `transcript`, `protein`, `region`, `element`, `domain`, `rule`, `param`, `cell_type`, `event`, `organism`, `stage`, `timer`, `signal`, `decision`, `experiment`; evidence and confidence on every block; `bio.std` prelude; programs test themselves | a written grammar and versioned spec; `population` and `reader`/`writer` constructs; a registry for shared libraries; std modules for signalling and human development |
| **BioIR** (representation) | 0.3 | 20 types with JSON round-trip: Evidence, Entity, Region, RegulatoryElement, Gene, Transcript, Protein, CellType, Effect, Event, Parameter, Rule, Domain, Signal, Timer, Stage, Decision, Experiment, Organism, Module; `UNKNOWN` as a value | a BIOIR-v0.3 spec (the doc is v0.1); ProteinState and Population as types; `Locus` moved under `ir` |
| **BioVM** (engines) | 0.3 | twelve engines: central dogma (exact on mtDNA and verified against UniProt), gene network (Hill ODE), Boolean (attractors), SBML (MathML + RK4), cell ageing, 2-D spatial fields, segmentation clock, gastrulation, Body (discrete events, one cell to the whole worm and a human as populations), debugger, uncertainty; composition through process-bigraph | one composite that runs Body + network + SBML on a shared clock; division and movement in space; external engines (libRoadRunner, MaBoSS, CompuCell3D) behind the same interfaces |
| **`bio` toolchain** | 0.1 | `bio check | compile | run | test | repl`; `bio test` runs the demo, std and organism programs in CI; depends on the engine alone, so it can be packaged as it stands | its own package and test suite; `bio fmt`; language server |
| **BioLib** (libraries) | 45 libraries + the proteome | data-computed membership from GO and Reactome (95.5% agreement with the curated catalogue); five layers (core, blueprint, timer, systems, parts); ligand-receptor protocol; haematopoiesis as a runnable mechanism module; **the whole human proteome packaged** (19,283 proteins in 2.3 MB, `genomeos protein X --lib`, answers offline) | `cancer.*` layer; more human development modules as runnable BioLang; "phenotypes it must reproduce" lists per library |
| **BioTwin** (individual) | 0.9 | any GRCh38 genome imported locally and used by carrier lookup, the gene report, the twin and the gene-by-gene walk; HG002 across all 22 autosomes with measured telomere and epigenetic age, the telomere also read from the 601 GB BAM by range without a download; fork, run, diff; predicted effect of a person's own regulatory variants; ClinVar carrier screen, genome-wide truncating-variant scan, coding inventory by consequence with AlphaMissense scores kept off the repository, a one-page dossier; the GIAB trio read for inheritance, 96.4% of HG002's calls inherited | true de novo variants out of the trio (phasing waits for a phased import); a population distribution for the polygenic scores; the effect per regulatory variant (area I's model) |
| **BioForge** (design) | search only | random-restart search under constraints with predicted labels; minimal-cell design estimate | experiments as input; published perturbations reproduced; constraint language |
| **Genome decoding** (GenomeOS) | 0.9 reached | every chromosome fetched and analysed (2026-09-11): inventory, UNKNOWN classified genome-wide with curated repeats (98.5%), domains on 24, curated repeats and ENCODE elements on 25, reader on 25 with eleven cell types; the node model tested genome-wide (90.2% of enhancers act inside their node) and on two mouse chromosomes; the segment parser on chr21 at 92.7% gene precision once RNA over the exons counts as evidence | the parser on every chromosome with measured RNA; enhancer targets and boundaries from Hi-C |
| **Molecules** (GenomeOS) | 0.9 reached | the whole human proteome compiled from seven public databases (19,478 coding genes, 25 chromosomes) and packaged as a 2.3 MB offline library, translation verified on 19,249 of them, genome-wide knowledge graph, RNA layer with GTEx, pathways as reachability and as kinetics where a curated ODE model exists; 96,362 modifiable sites with their 352 writers as `modifies` edges in the graph; the dominant isoform per tissue from GTEx (APP695 in the cerebellum) | measured modification state; isoform expression per cell type or stage; an engine that handles 100-species models |
| **The 98%** (GenomeOS) | genome budgeted | every UNKNOWN block of the 24 chromosomes tiered with Zoonomia constraint read per base (1,009 Mb): structural 23%, fossil 33%, regulatory 34%, neutral 7%, constrained-unknown 3.2% (1,098 blocks, 32 Mb, 1% of the genome); `genomeos budget`, the Progress tab card; the attributions compiled to a BioLang program **on all 24 chromosomes** (2026-09-17, 940,803 stated facts, 93.6% predicted, none experimental); scored against VISTA, eQTL, lentiMPRA, GWAS and ClinVar; the organiser reads copies first and leaves a real unknown of 882 blocks, 69 of them constrained on both axes; **0.45% of the unknown space has ever been measured by any assay held here, and the experiment that would change it is designed and unrun (284,001 oligos)** | measurement over the real unknown, which is the binding constraint; closure at cell and organism level |
| **Cancer and therapeutics** (GenomeOS) | 1.0 of the design dataset | tumour-only pipeline, cohort expression, altered-protein reconstruction, 15 mechanisms, design dataset with negative set | CNA/SV profiles; benchmark against approved targets; peptide/HLA predictor |
| **Web UI and CLI** (GenomeOS) | 17 views, 63 commands | every layer visible with evidence pills; jobs with progress; gene dossier; the Progress tab says what is going on, what is planned and, since 2026-09-17, the ordered next steps of §5; the Evidence explorer reads the project's own confidence (26,623 facts in the hand-written programs, and 967,426 with the 24 compiled chromosomes switched on, where 93.6% is predicted and none is experimental) | Cell view; release tags; nightly CI |

Data architecture: `data/reference` and `data/knowledge` are caches (git
ignored, rebuilt by `genomeos data fetch` and the compilers); `data/results`
holds the committed summaries, one per real-data run; `data/organisms` and
`data/demo` hold BioLang programs; `data/jobs` holds job metadata and logs.

---

## 3. Areas

Each area lists its goal, the code that implements it, the design document,
the data it produces, what is missing as of the review, and its next steps
in order. "Owner" is the session that holds the files today (see §7).

### A. BioLang and BioVM (the engine)

- **Goal.** A language a bio-engineer writes biology in, an IR that carries
  evidence, confidence and UNKNOWN, engines that run it, and a toolchain that
  stands alone (`bio run | check | compile | test | repl`).
- **Code.** `lang/parser.py`, `ir/model.py`, `runtime/*`, `std/*.bio`,
  `bio.py`. Versions: BioLang v0.3 (organism layer), BioIR 0.3.
- **Design.** BIOLANG-v0.1/0.2/0.3.md, BIOIR-v0.1.md, BIO-TOOLCHAIN.md,
  ARCHITECTURE.md §10, SEQUENCE-GRAMMAR.md.
- **Data.** `data/demo/*.bio`, `genomeos/std/*.bio`; `bio test` runs eight
  files in CI.
- **Requirements.** Every construct has evidence and confidence; programs
  test themselves (`# test:` and `assert:`); a `.bio` file runs without the
  rest of GenomeOS; the same program runs on any engine that fits it.
- **Stage 2 of v0.4: the economy, built and gated, and one of its two gates
  is a negative (2026-09-21).** `pool`, `cost` and `allocation` parse and
  refuse (`lang/parser.py`, `ir/model.py`); `runtime/economy.py` charges a
  cost, compares demand with capacity and scales each entity by its tightest
  pool. **Gate (a), burden:** the unrelated gene's factor falls 1.0 → 0.05
  over seven pre-registered demand levels and returns to 1.0 for all of them
  when the pools are multiplied by 1e6, which is the clause a bug would have
  passed; the order-of-magnitude clause against Ceroni 2015 and Frei 2020 is
  registered **not askable**, because neither curve is held here, and a test
  stops it contributing to a pass. **Gate (b), absolute abundance: failed.**
  PaxDb's integrated human set against GTEx v8 fibroblast median TPM, 15,159
  genes: log10 MAE **0.952005** for the one-to-one baseline and **0.952005**
  for the allocation arm, improvement **0.000000000**, bootstrap [0.0, 0.0],
  bar 0.05. The reason was derived from the runtime and committed before the
  transcript arm was fetched (`b57de43`): `_share` returns `capacity / wanted`
  to every demander, so the pool layer multiplies every gene's prediction by
  one constant and a fitted baseline absorbs it exactly. The run was done with
  the factor at 1.754e-4, moving every prediction 3.8 decades, and the error
  did not change in the sixth decimal. §5.3 had already written the
  consequence: **the pool layer stays optional.** It buys a burden, it refuses
  six modelling errors, and it improves no abundance prediction. Findings
  beside the verdict: at a real transcriptome's demand §5.1's capacities never
  bind; **`competitive` is as gene-blind as `proportional`**, saturating on the
  pool's total demand rather than per demander, which closes §10 decision 5
  with a negative; the unfitted secondary lands at **1.36e10 proteins per cell**
  against Milo 2013's 4e9 to 1.2e10; and the newly registered slope of log10
  protein on log10 transcript, **0.515 [0.499, 0.531]**, says a future
  per-demander policy would have to supply 0.485 decades per decade of
  compression. 20 tests across `tests/test_pools.py`, `test_economy.py`,
  `test_burden_gate.py` and `test_abundance_gate.py`, three of the eight new
  ones being contrasts that would catch machinery unable to see a difference at
  all. Spec BIOLANG-v0.4-ECONOMY.md §5 and §9.2; results
  `data/results/burden_gate.json` and `abundance_gate.json`, from
  `scripts/burden_gate.py` and `scripts/abundance_gate.py`. **Missing after
  this:** a per-demander saturable policy, the only change that could make gate
  (b) answerable, and stage 3.
- **2.0 part two, done (2026-09-27, lane-engine2, `14f9ada`).** The `biolang` package describes its own
  language: the grammar (regenerated from the packaged parser and re-rendered in isolation to match)
  and the v0.1–v0.4 specs ship in `biolang/docs/` inside the wheel, under CC BY 4.0. `bio repl` is
  driven through stdin in the isolated check. `import protein:` has an honest default: the engine
  ships no proteome and stops with an error naming the missing resolver (never an empty module, never
  the network); GenomeOS supplies it as the `protein` entry point of `biolang_import_resolvers`.
  12 of 12 isolated checks; engine suite 154 passed, 3 skipped. **Negatives first in the lane's
  report:** the 2.0 `bio run` smoke check had passed on an all-zero run (`basal_rate:` is accepted and
  ignored by the parser; now `basal:` with a level-above-zero test), and `bio check` crashed on
  `import protein:` even with GenomeOS installed. **Missing after this:** the parser rejecting unknown
  block keys; the engine package ships copies of `tests/*.py` that LICENSING.md lists as AGPL.
- **Missing.** A `biolang` package with its own tests; imports from a
  remote registry (the resolver hook exists: `import protein:TP53` reads the
  packaged proteome through `IMPORT_RESOLVERS`). Done 2026-09-11: the grammar
  and the IR type list are generated from the parser and the dataclasses
  (BIOLANG-GRAMMAR.md, held current by a test); `bio.std.human_stages`
  (Carnegie), `bio.std.signalling` (148 pairs from CellPhoneDB) and
  `bio.std.human_turnover` (Sender & Milo lifespans); `express` (reader
  state), `field`, `design`, space and populations in the language.
- **BioLang v0.4: structure, economy and control (2026-09-14).** The
  specification is [BIOLANG-v0.4-ECONOMY.md](BIOLANG-v0.4-ECONOMY.md): a cell's
  economy and structure belong in the language, not outside it. Stage 1 is
  built and has passed its gates — `compartment` (a containment tree with
  volumes, genomes and ribosomes), a location and targeting `signals` on every
  gene and protein, `transport` between adjacent compartments with one shared
  saturable capacity, and a `regime` block that declares and records how a run
  was executed (continuous or stochastic per species, update scheme, allocation
  policy, seed). `runtime/located.py` runs it: 13 of 13 chrM proteins are made
  and kept inside the mitochondrion with no transport, all 1,123 nuclear-encoded
  MitoCarta3.0 proteins reach it only through a transport (0 with the routes
  closed), the rho0 phenotype comes out right (complexes I, III, IV, V dead,
  complex II intact: King & Attardi 1989), removing one subunit's presequence
  strands it and kills its complex alone, and a red blood cell with no genome at
  all still runs. Six modelling errors are compile errors instead of silent
  answers. Stages 2 to 4 (pools, costs and a named allocation policy; core
  metabolism and heteroplasmy; partitioning division) and the control layer
  (`homeostat` with set points and gains, `role` for phenomenological
  behaviour, specialisation as state rather than inheritance) are specified with
  their falsifying measurements and wait for Albert's steer. Programs:
  `data/organisms/human/oxphos.bio` (generated, self-testing) and
  `data/organisms/human/erythrocyte.bio`; gate result
  `data/results/located_mitochondrion.json` from
  `scripts/located_mitochondrion.py`.
- **Fate: commitment and competence, gated by a published series (2026-09-14).**
  §7.2a of the v0.4 specification had both constructs inferred at 0.3, with our
  own atlas test negative and in the wrong direction. The falsifier it names is
  a perturbation, and it is now a program: `data/organisms/celegans/plasticity.bio`
  runs Fukushige & Krause 2005 and Yuzyuk et al. 2009 under `bio test` — forced
  at 60 min, 610 of 610 cells become muscle; forced at 350 min the embryo keeps
  every fate it had; without MES-2 the window never closes and 594 convert;
  without PHA-4 (their own control) nothing changes; and with the window held
  open, no cell that had already differentiated converts. Removing either
  construct breaks an arm (`scripts/plasticity_gate.py`,
  `data/results/plasticity_gate.json`), so neither is decoration. The runtime
  gained a perturbation that arrives at a time (`add: F at T`), and `bio test`
  now runs every organism program's `experiment` blocks — the eight worm and
  nine haematopoietic knockouts included.
- **Four defects and one language gap, from area E using it (2026-09-15).** A
  `cell_network` cadence let a rule that had lost the precedence contest win it
  later, costing the worm 17 terminal fates (484 of 555 against 501) for no
  reason but the cadence; two decisions sharing an id silently made one
  unreachable; a birth made its neighbours re-read their contacts but a death did
  not; and `asymmetric: X -> D` created X in the keeper when the mother had none.
  All four are fixed, counted (`overruled_fates`, `duplicate_decision_ids`) and
  regression-tested, with every committed organism program identical before and
  after. The **integrated reads** of §7.2a (`F.exposure(lineage)`,
  `F.mean(cell)`) are now a language read the Body computes, in absolute units,
  so area E can write their fate rules as declared windows instead of a
  precomputed lookup; a read that does not name its window is a compile error.
- **Order, and the engine packaged on its own (2026-09-15).** `order` (§7.6) is
  the benchmark's untestable sentence written down: HOXD colinearity as a claim
  the compiler refuses when a sequence contradicts its own coordinates, and a
  sequence along time reported against a run with its inversions. It drives
  nothing, which is what keeps it from being a second `stage`. And the engine now
  **builds and runs as a package of its own** (`scripts/package_engine.py`,
  `docs/ARCHITECTURE.md` §10.1, milestone 2.0): **36 files on 2026-09-27, 22 with imports rewritten, 8 of 8 checks** (31 files and 20 rewritten when first built on 2026-09-15, which this row quoted until today), Apache-2.0, no
  dependencies, run in a Python with no site-packages where `import genomeos`
  fails outright — all four verbs, every module, the standard library and an
  organism, with the application absent. *(Later on 2026-09-27, `5720c7c`: **9 of 9 checks**, the ninth being the engine's own pytest suite — 24 test files and 157 tests selected by their imports, 21 fixtures followed through `.bio` imports, 154 passed and 3 skipped for the optional `process-bigraph` extra — and its own version, `biolang` **0.1.0**. Milestone 2.0 is reached; see section 6.)*
- **A `population` type: judged, and it does not earn its place (2026-09-15).**
  Judged by the rule `order` was judged by — does it let a program *say* anything
  it cannot say now? A counted `Cell` already carries the count, the factors, the
  levels, the type and all five actions with `fraction`; a `Population` type would
  reorganise fifteen branches in the runtime and leave every program able to say
  exactly what it says today. The structure a real population has and a count does
  not — an age distribution, a size spread — is asked for by no committed program
  and by no measurement the project holds. It is refactoring, so it is not built.
  **What judging it did find** is in the body program's numbers. Splits at one
  decision point run in sequence, so a splitting `fraction` is a share *of what is
  left*, and `genomeos/organism/human.py` converts each population's published
  adult count into one. The run is correct, but every number after the first in a
  germ layer is not the share its own `evidence` line claims it is: ectoderm's
  glia read `fraction: 1.0000` for a 24.5% share, and mesoderm's myocytes read
  `1.0000` for a share of 0.0000071 — **140,000× apart**, with the evidence text
  saying "share of the layer's adult cell count" in both cases. The value also
  depends on the order the decisions are emitted in, which is the order-dependence
  §7.3 removed for fates. Two fixes, neither of them mine to make alone: the
  generator should say what its number is, and the language wants a `share:` that
  is absolute and normalised across the splits at one decision point, so a program
  can state the published number rather than a derived one. Specified here rather
  than built, because a construct no program uses is what this area has spent the
  day refusing to ship — it should land with the generator that will use it.
- **Next.** 1. `share:` for a partition, so the body program's numbers mean what their
  evidence says: **done 2026-09-17, and its first real use found an error.** The sixteen
  germ-layer splits of `tissues.bio` now state the published share of the layer instead of a
  fraction of what earlier splits left. Three populations came back from nothing — Glia,
  LungEpithelium and Myocyte, each the LAST split of its layer, which under `fraction` took
  1.0000 of a remainder that the earlier splits, printed at four decimals, had already
  exhausted. **Glia is 24.5% of the ectoderm**, so a quarter of a germ layer was being
  absorbed by its siblings and the population held no cells in a twenty-year run. Nothing
  caught it because the layer totals were conserved: the 20-year count moves by -1.8e-6 and
  every assert passed before and after. Two committed numbers were wrong and are corrected
  with their reason (`body.bio` cells, and the same figure in `test_human_body`). 2. Splitting the repository, which is a
  release decision rather than an engineering one now that the package builds.
  *Closed since this list was written:* the precursor-level `commitment` (measured
  closed — 0 of 12, §7.2a decision 3), `bio` packaged (above), and PAR polarity,
  which landed with area E.
- **Owner.** genomeos-c2 (the engine: language, IR, runtime, toolchain, v0.4).

### B. Genome decoding (reverse-engineering the sequence)

- **Goal.** Every base of a real genome accounted for: gene, transcript,
  exon, regulatory element, repeat, domain, or an UNKNOWN with a class and
  a confidence; the delimiters the cell reads learned from the sequence.
- **Code.** `genome/annotation.py`, `blocks.py`, `unknown.py`, `repeats.py`,
  `regulation.py`, `domains.py`, `signals.py`, `anatomy.py`, `fetch.py`,
  `lookup.py`, `segments.py`; `predict/rna_tracks.py`.
- **Design.** GENOME-AS-CODE.md, GENOME-ANATOMY.md, SEQUENCE-GRAMMAR.md,
  UNKNOWN.md, NODES-READER-WRITER.md.
- **Data.** `anatomy_hg38_by_chromosome`, `unknown_chr*` (25),
  `rmsk_chr*` (18), `ccres_chr*` (24), `domains_chr*` (2),
  `signals_chr21`, `segments_chr21`, `mouse_mm10_chr*` (2),
  `gencode_v50_chr*` rows (18).
- **Requirements.** Any chromosome on demand (`genomeos data fetch --chrom`);
  classification order fixed (curated repeats before ORFs); every class
  names its evidence layer; the block map shows it.
- **Missing.** Silencers have no source; enhancer → gene is
  inferred from the CTCF domain, never measured (Hi-C or a predictive model);
  the sequence grammar has signals but no segment parser (Viterbi over the
  grammar) and no JASPAR promoter scan; the **reader** of
  NODES-READER-WRITER.md (which nodes are open in which cell type, from
  DNase/ATAC/methylation) has no code; no second mammal for node comparison.
- **Done since the review (2026-09-11).** A fetched chromosome arrives
  analysed (`data fetch --analyse`: curated UNKNOWN pass and domains);
  reader v1 (`genomeos reader`: ENCODE DNase peaks per cell type over
  nodes, promoters and enhancers; K562 and HepG2 on chr21, chr22 analysed).
- **Enhancer to gene, genome-wide (2026-09-11).** Deleting an element in its
  1 Mb window and reading which gene moves turns "nearest coding gene in the
  CTCF domain, inferred 0.4" into a named gene with a tissue and a magnitude.
  Across all 24 chromosomes, 4,800 distal enhancers deleted one at a time:
  2,994 move a gene, 2,291 name a coding gene, **90.2% of those inside the
  element's CTCF node** (79.8% to 91.8% per chromosome) and 71.2% exactly the
  nearest TSS. **Qualified on 2026-09-14**, and the qualification matters more
  than the headline: those 4,800 were a sample of distal enhancers, and over
  the whole deletion archive (113,399 elements, chr15 to chr22 and chrY
  complete) the same caller keeps a coding target inside the element's node for
  81.7% of elements against 79.1% for the same number of boundaries placed at
  random. The node's advantage over a random partition of the same resolution
  is 2.6 points, so "the boundary predicts where an element acts" is a weak
  effect measured on a model's reading, not the nine-times-out-of-ten the
  sample suggested (NODES-READER-WRITER.md, "The orientation-aware caller"). Two findings worth as much as the headline: the nearest-gene
  heuristic names the wrong gene almost a third of the time, and 1,157
  elements the registry calls enhancer-like behave as silencers, rising when
  deleted. The registry names a class of element, not a direction.
- **Segment parser: three runs, one conclusion (2026-09-11).** The parser
  runs Viterbi over the grammar and is scored against every coding transcript
  on chr21. All three results are kept side by side so the comparison stays
  checkable.

  | signals used | candidates | exon precision | gene precision | gene sensitivity | canonical exons found |
  |---|---|---|---|---|---|
  | learned donor/acceptor matrices | 556 | 5.2% | 21.9% | 89.1% | 8.7% |
  | + AlphaGenome splice sites (feature c) | 1,003 | 46.2% | 20.0% | 84.2% | **74.8%** |
  | + starts anchored at ENCODE promoters | 158 | 61.9% | 53.8% | 59.3% | 56.3% |

  The last column matters more than it looks. Scored against every coding
  transcript, the parser finds 40% of CDS segments; scored against the
  *canonical* transcript of each gene it finds 74.8%. Most of what looked like
  failure was the parser not reproducing isoform-specific exons, which is a
  different and much smaller problem than not finding the gene's exons at all.

  Better splice sites multiply exon and site accuracy eight-fold and leave
  gene finding untouched. Anchoring the start at a curated promoter is the
  first thing that moves gene precision, 2.5-fold, but it buys that precision
  with sensitivity: it can only find genes the registry already marks, so a
  third of the genes drop out. Two independent signals, neither of which found
  a gene the grammar was missing. **The coding model is the limit**, and that
  is a measured claim now rather than a suspicion.
- **Reader lane in the block map (2026-09-11)**: nodes fill with their open
  fraction, silent ones take a red edge, and each coding gene carries a read
  or silent dot for the chosen cell type. The same chromosome through a
  different cell is a different map.
- **The second mammal (2026-09-11).** `genomeos mouse --chrom chr19` brings
  mouse mm10 chromosome 19 through the same three doors as a human chromosome
  (UCSC sequence, GENCODE vM25 models, ENCODE mouse elements) and infers its
  nodes with the *unchanged* human code: 371 nodes, 256 with coding genes,
  median 109 kb, against 20,002 human nodes indexed for comparison. Of the 91
  mouse nodes holding two or more genes whose symbol matches a human gene, 49
  land inside a single human node and 39 more split across *adjacent* human
  nodes, leaving 3 scattered. Synteny is well known, so the interesting number
  is not that the genes stay together; it is that when a mouse node does not
  map to one human node, the pieces are next to each other 39 times out of 42.
  The node behaves as a unit in both genomes, and the CTCF-only boundary is
  the resolution limit rather than the biology. That first run used gene-symbol
  identity for orthology and was written up as a lower bound; re-run with MGI's
  curated homology (48c48bc) the tested set grows from 91 nodes to 109 and the
  shape holds: 53% in one human node, 94% in the same neighbourhood, adjacent
  45 times out of 51. The stronger method moved the count and left the
  conclusion, which is what a lower bound is supposed to do. Both columns are
  kept in the result.
- **The second mouse chromosome (2026-09-11).** mm10 chr11 through the same
  code: 836 nodes, 270 tested, 59% in one human node, 93% in the same human
  neighbourhood, 20 scattered. With chr19 that is 379 tested mouse nodes at
  93 to 94% in the same human neighbourhood: two chromosomes, one shape.
- **The coding model was not the limit either (2026-09-11).** The claim
  above was tested as it asked to be. A 3-periodic fifth-order Markov model of
  coding sequence (`segments --coding markov`) in place of the codon table
  gains one to two points of precision in every configuration and no
  sensitivity on chr21. Both models are kept. What separates a real gene from
  a plausible reading frame is not a better score of the sequence but
  evidence that the sequence is transcribed.
- **Transcription evidence, two ways (2026-09-11).** Confining gene starts to
  the reader's DNase peaks (eleven cell types, 38% of chr21) cuts candidates
  1,003 to 614 and lifts exon precision 46% to 53% at 76% gene sensitivity, a
  milder trade than the promoter class. RNA over the exons is the lever:
  AlphaGenome's predicted RNA-seq coverage for eight tissues, cached per 1 Mb
  window (`genomeos/predict/rna_tracks.py`), applied after the parse as a
  filter (`segments --rna-filter`).

  | signals used | candidates | exon precision | gene precision | gene sensitivity |
  |---|---|---|---|---|
  | AlphaGenome splice sites | 1,003 | 46.2% | 20.0% | 84.2% |
  | + starts confined to open chromatin | 614 | 53% | not recorded | 76% |
  | + RNA over the exons, eight tissues | **82** | **68.6%** | **92.7%** | 57.9% |
  | + measured RNA, ENCODE K562 total RNA-seq | 55 | 72.1% | 89.1% | 44.8% |
  | + measured RNA, six ENCODE cell lines pooled | 127 | 66.3% | 72.4% | 63.8% |
  | predicted panel ∪ measured lines | 137 | 64.4% | 72.3% | 67.0% |
  | predicted panel ∩ measured lines | 72 | 71.2% | **95.8%** | 54.3% |

  Nine candidates in ten were sequence that is never transcribed in the
  panel. The panel is eight tissues rather than the body, and that is where
  the lost sensitivity went. SEQUENCE-GRAMMAR.md item 1 closes the series:
  the grammar finds the genes, transcription evidence says which are real.
- **Model and measurement agree (2026-09-11).** `segments --rna-measured
  K562` reads ENCODE's strand-specific total RNA-seq bigWigs (ENCSR860DWK)
  over the candidate exons only, through the range reader of
  `attribution/bigwig.py` (0.7 MB moved for 3,248 exons), and keeps the
  candidates with signal on their strand (`genome/rna_measured.py`). Same
  exon and gene precision as the predicted panel, fewer genes found, because
  one cell line expresses fewer of chr21's genes than eight tissues do; the
  signal threshold was swept (1.0 leaves 12 candidates, 0.05 the 55) and the
  sweep sits next to the table in SEQUENCE-GRAMMAR.md. Pooling six lines
  (K562, HepG2, GM12878, IMR-90, A549, MCF-7; HeLa-S3 and SK-N-SH have no
  released strand-specific track, recorded as missing) buys sensitivity, 44.8%
  to 63.8% of genes, and pays for it in precision, 89.1% to 72.4%, because
  measured transcription includes the non-coding transcription the model's
  gene-level tracks exclude. A model that knows genes plus a measurement that
  knows the cell is the honest filter. Results are
  `segments_chr21_predicted_sites_measured_<cell|panelN>`. A measurement
  lesson found by the closure test (area I) and fixed here on 2026-09-12:
  ENCODE labels a strand track by the read, and IMR-90's total RNA-seq reads
  antisense, so `MeasuredRna` now probes sixty first exons per strand in both
  orientations and swaps when the swapped one carries twice the signal (K562
  26.4 against 1.8, HepG2 19.7 against 1.4, GM12878 36.3 against 6.4, A549 29.8
  against 5.6, MCF-7 36.2 against 2.2, IMR-90 0.85 against 22.4, swapped).
  Rerun oriented (abfbd4c), the tables above carry the six-line numbers: the
  inverted line had not merely contributed nothing, read backwards it passed
  antisense candidates over real genes, which is why the measured panel's gene
  precision read 63% and 59%; oriented it is 72% and 74% at the same
  sensitivity, and the series' ordering does not move. SK-N-SH and HeLa-S3
  have no total RNA-seq bigWig on ENCODE.
- **The two filters combined (2026-09-12).** `--rna-combine union|intersection`
  applies when both the predicted panel and the measured lines are given.
  The union is the most sensitive call of the series, 67.0% of genes; the
  intersection is the strictest, 69 of 72 candidates are genes (95.8%
  precision at 54.3% sensitivity). The model supplies the gene-level
  judgement and the measurement the cell; asking for both is how a call
  becomes actionable. Results
  `segments_chr21_predicted_sites_rna_measured_panel6_{union,intersection}`.
- **The second chromosome (2026-09-12).** chr22 (447 coding genes, 8,669
  CDS segments) through the same five configurations, nothing tuned:

  | signals used, chr22 | candidates | exon precision | gene precision | gene sensitivity | canonical exons |
  |---|---|---|---|---|---|
  | learned matrices | 297 | 8.3% | 30.0% | 72.3% | 5.9% |
  | AlphaGenome splice sites | 1,293 | 49.4% | 18.7% | 97.3% | 70.4% |
  | + predicted RNA, eight tissues | 151 | 63.1% | 92.7% | 81.2% | 62.8% |
  | + measured RNA, six lines | 225 | 61.0% | 74.2% | 88.8% | 65.2% |
  | predicted ∩ measured | 142 | 63.7% | **95.1%** | 80.1% | 62.0% |

  The same ordering as chr21 in every row, with higher sensitivity because
  the eight-tissue panel covers chr22's genes better; chr21 carries more
  tissue-restricted and keratin-associated genes. Two chromosomes, ten runs,
  the series holds, and the parser thread of this area closes: the grammar
  finds the genes, splice-site prediction draws their exons, transcription
  evidence says which are real, and asking the model and the measurement
  together is the actionable call. The Progress tab shows the whole series
  ("The segment parser, by evidence", one table per chromosome from the
  committed results), next to the composition budget of area I.
- **Node edges against a predicted contact map (2026-09-12).** `genomeos
  domains --chrom chr21 --contact-map` takes the insulation minima of
  AlphaGenome's predicted contact map (45 windows, 2 kb bins, 28 4DN Micro-C
  and Hi-C cell types averaged) as predicted boundaries and holds the 227
  CTCF-only boundaries against them (`predict/contact_maps.py`,
  `domains_chr21_contact`). 83 of 227 inferred boundaries sit within 20 kb of
  one of the 353 predicted minima (37%), against 28% for the same number of
  boundaries placed at random; 24% of predicted minima have an inferred
  boundary; median distance 34 kb; the depth threshold swept 0.02 to 0.3
  never helps beyond chance. Agreement a little above chance, recorded as
  such: node content has the evidence (enhancer deletions 90% inside, mouse
  93 to 94% same neighbourhood), node edges stay a proxy at 0.4 until
  measured Hi-C or a better insulation predictor.
- **Orthology from Compara (2026-09-12).** `genome/mouse.py` streams the
  mouse-side Ensembl Compara 116 dump once (112 MB) and keeps the mouse to
  human orthologue pairs for the mouse chromosomes fetched (2,176 pairs for
  chr19 and chr11, `compara_mouse_human_orthology.tsv.gz`); `genomeos mouse`
  compares nodes under MGI and under Compara and reports their agreement.
  Same human neighbourhood: chr19 94.5% under MGI and 95.5% under Compara,
  chr11 92.6% and 91.8%; in one node 53.2% and 55.9%, 58.9% and 56.1%. Where
  both sources name orthologues for a gene they name the same human genes
  97.3% (chr19) and 96.3% (chr11) of the time; Compara covers 65 and 86 genes
  MGI lacks, MGI 34 and 42 Compara lacks. The node finding does not depend on
  the orthology source. A lesson recorded in LESSONS.md: Ensembl's human
  homology dump has Mus caroli, Mus spretus and the rat but not Mus musculus;
  the mouse to human pairs live in the mouse dump, and the first pass wrote an
  empty file and fell back to MGI silently, so the header now counts what it
  kept.
- **Next.** 1. The `reader` construct in
  BioLang, so context gating comes from
  chromatin (the parser half belongs to the language owner). 2. Done
  2026-09-12: measured Hi-C boundaries for the node edges (`genomeos domains
  --chrom C --hic BIOSOURCE`, `genome/hic.py`; the portal's Authorization
  header was the 403 with a valid key and is dropped on the redirect to S3).
  chr21's 227 CTCF-only boundaries against 4DN boundary calls within 20 kb,
  the same random control as the predicted map:

  | biosource, assay | measured boundaries | CTCF-only edges supported | at random | of the measured covered |
  |---|---|---|---|---|
  | GM12878, in situ Hi-C | 557 | **58%** | 42% | 27% |
  | H1-hESC, Micro-C | 149 | 24% | 13% | 36% |
  | K562 | 93 | 14% | 8% | 34% |
  | HepG2 | 93 | 15% | 8% | 37% |
  | IMR-90, dilution Hi-C | 126 | 12% | 11% | 22% |

  Clearer than the predicted contact map (37% against 28%) and in the
  direction the node model needs: the deepest measurement finds a boundary
  near six in ten CTCF-only edges, and the edges find one in four of its
  boundaries. Still a proxy, edges stay at 0.4; any 4DN biosource is one
  command away (`domains_chr21_hic_<biosource>`). The L1 ORF2
  and Alu sequence signatures are superseded by the RepeatMasker pass and
  stay as the fallback for chromosomes not yet distilled.
- **The epigenome layer (2026-09-14, genomeos-e1).** Openness said a node is
  open, not what it carries. `genome/epigenome.py` reads the chromatin itself
  from ENCODE for the reader's eleven biosamples: H3K4me3, H3K27ac, H3K4me1,
  H3K27me3 and H3K9me3 in all eleven, from one biosample each; GRCh38 WGBS in
  eight (cardiac muscle cell, keratinocyte and astrocyte have none: RRBS on
  hg19 only, or none). One record per locus per cell type carries openness,
  each mark (peak, fold change), methylation with coverage and an inferred
  state, with evidence and confidence per field and UNKNOWN with the reason
  (`epigenome_at`, `genomeos epigenome coverage|at|summary`,
  `epigenome_chr*` for all 24 chromosomes, `epigenome_genome_wide`). It is a
  data layer; no BioIR construct reads it. Three measured answers
  (docs/NODES-READER-WRITER.md):
  - **Mark identity beats registry class.** On 83,814 element-cell pairs of
    chr21 and chr22 scored per line, deletion direction goes from AUC 0.53
    (class) to 0.59 (with DNase) to 0.62 (with the line's marks), and
    magnitude from Spearman 0.18 to 0.22 to 0.33. Another line's marks and
    shuffled marks both fall short. Polycomb, bivalent and poised elements
    rise on deletion 44 to 46% of the time, active ones 24 to 25%.
    Caveat: the outcome is AlphaGenome, trained on these tracks.
  - **H3K27me3 explains the HOX clusters.** It sits on 56 to 100% of HOX
    promoters against 12 to 31% of GC- and CpG-matched promoters, and the
    exceptions are the expressed clusters. In H1 the DNase reader calls every
    HOXA promoter read and the marks call ten of eleven bivalent. It does not
    explain silent promoters in general: methylation covers 88 to 90% of them
    where the methylome is intact.
  - **Methylation does not explain the fossil tier.** Matched for GC and CpG
    density, it is methylated like unique unconstrained sequence: −1.5 to +0.8
    points in H1, hepatocyte and monocytes, split in hypomethylated lines, no
    H3K9me3 excess. All 328.5 Mb were read, 85 to 93% of windows measured in
    six biosamples (IMR-90 8%, monocytes 24%). **Negative.** Alu remains keep
    15 to 18 points more methylation in hypomethylated lines.
- **CTCF orientation at the node edges (2026-09-14, genomeos-e1).** Every
  CTCF-only element was scanned with JASPAR MA0139 and the 19,930 domain edges
  held against 4DN boundaries within 20 kb (`scripts/ctcf_orientation.py`,
  `domains_ctcf_orientation`).
  - **Positive control.** Orientation is real at measured boundaries: sites
    upstream are 65 to 67% reverse and downstream 65 to 66% forward in H1, K562
    and HepG2.
  - **Splitting our edges by orientation: negative.** Convergent pairs score at
    their strand-shuffle medians (p 0.17 to 0.69), divergent pairs no worse; on
    chr21 the 33 convergent edges score 55% and the 14 divergent 86% against
    GM12878. GM12878's dense calls put a random position within 20 kb of a
    boundary 54% of the time, so its 58% was mostly density.
  - **What limits the edges.** Having a site at all separates edges (1.5 to 2
    times the support), and 58% of edges have none. Used to place boundaries,
    orientation works: reverse-to-forward strand flips beat random 1.4 to 1.6
    times and forward-to-reverse flips fall below random.
  - **Decision pending.** The orientation-aware caller would move every node,
    so it waits for a decision.
- **The reader uses the marks (2026-09-14, genomeos-e1).** Not read: an open
  promoter with H3K27me3 and no H3K27ac (poised; H1 3,151, monocytes 2,434).
  Read: a closed promoter with H3K4me3 and H3K27ac (keratinocyte 2,836,
  GM12878 1,096). `silent_genes` is no longer cut at 200, which had shown every
  later silent gene as read. Against measured RNA in four lines, poised genes
  are expressed 6.5 to 16.5% of the time, like closed ones. Precision of "read"
  rises 56 to 65% (K562), 55 to 62% (HepG2) and 47 to 52% (IMR-90), recall
  falls under one point, and GM12878's recall rises 75 to 84%
  (`reader_poised_check`).
- **The Alu lead, closed (2026-09-14).** Alu minus L1 methylation in fossil
  blocks falls from 15 to 19 points raw to 4.6 to 12.7 once GC, CpG density,
  age and solo-WCGW share are matched. The matched L1s are atypical CpG-rich
  outliers, and methylation falls with solo-WCGW share inside both families, so
  the excess is context and region; a family programme is not shown
  (`epigenome_fossil_alu`).
- **Fold-change profiles, all 24 chromosomes (2026-09-14, genomeos-e2).** The
  broad job was deferred once as too slow; the handshake was the reason. One
  profile is 25 to 27 range requests over **two** connections (chr21: 233,550
  bins, 18.5 MB, 22.0 s; chr20: 322,221 bins, 21.8 MB, 34.3 s), and before the
  range reader kept connections open (46cb5d2) each of those requests paid 0.5
  to 17 s of TLS, 12 to 425 s a profile. Retried through that reader it is the
  bytes that dominate, at 0.6 to 0.8 MB/s a stream and six to ten profiles a
  minute across eight workers. All 1,320 profiles (eleven biosamples × five
  marks × 24 chromosomes) are now read: 2.2 GB kept from roughly 55 GB
  streamed, the last 593 in 59 minutes. Methylation is complete on all 24 for
  the eight biosamples with GRCh38 WGBS; three have none at all. `genomeos
  epigenome coverage` reports it, and the layer no longer answers UNKNOWN with
  "this chromosome not read".
- **The marks' gain transfers, and is a rounding error beside the model's own
  lines (2026-09-14, genomeos-e2).** The chr21/chr22 comparison was refitted on
  chr21 and chr22 alone and scored on chr14 to chr20: 564,824 element-line
  units over 141,209 elements, 115,096 acting
  (`scripts/epigenome.py direction-transfer`, `epigenome_direction_transfer`).
  chr1 (1.6% swept), chr13 (being written) and chrY (two of the four lines are
  female) were refused as held out, with the reason recorded.
  - **It transfers.** Direction goes 0.601 to 0.626 AUC off the training
    chromosomes, against 0.594 to 0.616 in sample, positive on all seven
    separately (+0.019 to +0.038). Shuffled marks score the fit without them.
  - **Half of it is not the cell's own chromatin.** Another line's marks reach
    0.614 of the 0.626; the line-specific part is +0.010 to +0.014 AUC.
  - **The model's own direction in the other three lines scores 0.907 alone**,
    against 0.626 for the whole marks model (+0.282 to +0.291). Adding the
    line's marks to it moves direction by +0.0027 to +0.0035 and acts by
    +0.0003 to +0.0022, and makes it *worse* for acts on chr15 and chr19 and
    for magnitude on chr15 and chr19. **Report the marks at that size.**
  - **Nothing here is measured.** Outcome and competing feature are both
    AlphaGenome, trained on these lines' ENCODE tracks.
- **The orientation-aware caller, as an alternative (2026-09-14, genomeos-e1).**
  `infer_domains(..., orientation=strands)` places boundaries at reverse-to-forward
  CTCF motif flips among elements with CTCF ChIP support. The default output is
  byte-identical, and every committed node is unchanged. Four measurements on
  the same chromosomes, CTCF-only against oriented (`domains_oriented_comparison`,
  NODES-READER-WRITER.md):
  - **Hi-C.** Oriented wins. Enrichment over random is 1.74, 1.96, 1.99 and
    1.45 in H1, K562, HepG2 and IMR-90, against 1.38, 1.36, 1.30 and 1.16.
  - **Node content.** Oriented loses. Over 113,399 archive elements naming a
    coding gene, 63.2% keep the gene inside the node, against 81.7%. Random
    boundaries give 77.7% and 79.1%, and at matched resolution the shares are
    67.4% against 82.2%. The old caller's own excess over random is only 2.6
    points on the full archive.
  - **Mouse synteny.** Oriented loses. The same human neighbourhood holds for
    91.3% and 88.0% of mouse nodes (chr19, chr11), against 95.9% and 92.4%.
  - **HOXD.** Neither caller recovers the boundary. Oriented splits the cluster
    between HOXD8 and HOXD4. Its rule does place a flip inside the published
    interval, which the inherited 50 kb node floor then removes.
  - **Verdict.** Keep both; the default stays CTCF-only. Insulation boundaries
    and the model's enhancer reach disagree about where regulation is bounded.
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
- **A twelfth and thirteenth biosample, and the reading they were fetched to
  test turns out to be assay depth (2026-09-22, lane-biosample, `ca8aa0f`
  inventory and registration, `1680ebc` result).** The inventory came first and
  had never been written down, and two of its facts decided the choice: **not
  one of the eleven biosamples was a tissue** — six immortalised or cancer
  lines, two directed-differentiation products of a stem line, three cultured
  primaries — and **neural was already occupied twice**, by SK-N-SH and
  astrocyte. So the committed rule (unoccupied lineage first, then whether the
  direction predicts a reading outside the eleven's range for a stated
  mechanism, availability only as a tiebreak) chose **gonadal**: testis, with
  ovary as the control that separates "germline" from "bulk tissue", since
  testis moves both at once.
  **The result that matters is a negative about a reading this project already
  ships. `enhancers_active` is predicted out of sample by DNase peak count
  alone**: fitted on the eleven, `enhancers = 57,968 + 0.3910 × peaks`, testis
  lands at **z = −0.11** and ovary at −1.29. A lineage and a material the layer
  had never seen have their enhancer count predicted by assay depth. **So
  `enhancers_active` must stop being read as a count of enhancers.** The one
  reading that escapes it is `genes_read`, where testis reads **15,078 of
  20,094, above the eleven's entire range**, on mid-range depth — the expected
  direction, since testis has the broadest transcriptome of any human tissue.
  **The germline methylation prediction was refuted and the test it was scored
  against shown non-discriminating**, which is the more useful half: the somatic
  range is a union of two non-overlapping populations, transformed lines at
  0.074–0.604 and non-transformed at 0.805–0.852 with **nothing between**, so no
  value could have failed it. The methylation layer separates transformed from
  non-transformed by a margin larger than any lineage effect in it.
  **And a prediction made from portal byte sizes before downloading anything was
  confirmed exactly**: testis yields **321 H3K9me3 peaks genome-wide and 1,902
  H3K27me3** against floors of 1,770 and 8,540, so those two marks are
  **under-called and not measured absence** — its 0.0209 silent-promoter share
  is a peak-call artefact. Ovary's H3K27ac sits 4% under the floor, unregistered
  and reported as such.
  All **264 of 264** existing (chromosome, biosample) blocks reproduced
  byte-identically. 1.91 GB over the wire, 113 MB on disk, 0 model requests.
  Signal profiles for the two were deliberately not fetched, so the direction,
  fossil and transfer analyses still exclude them.
  **What area B's Missing should now say**: the set is thirteen and holds its
  first bulk tissue and first gonadal lineage, and it remains **without a
  germline-resolved sample** — bulk testis dilutes spermatogenic cells with
  Sertoli, Leydig, peritubular and blood cells, so that question needs a sorted
  or single-cell source rather than another bulk tissue. "No second mammal for
  node comparison" stands unchanged.
  **Operational, and it is Albert's**: disk is at **95%, 21 GiB free**, worse
  than the 92% recorded on 2026-09-17. This lane's 113 MB is negligible and the
  trend is not.
- **The whole per-biosample family banded against depth, and two more defects
  under it (2026-09-22, lane-depth, `41bc7d6` registration, `86d668f` result).**
  *(Correction to the bullet above and to LESSONS.md: `enhancers_active` is the
  **DNase reader**, `genomeos/genome/reader.py:203`, not the epigenome layer.)*
  The confound is a property of the pipeline, so all nine readings that can
  carry a claim were checked against the covariate that actually drives each.
  **`nodes_open` turns out never to have measured anything**: `reader.py:197`
  takes the median peak density over that biosample's own nodes and calls a node
  open at or above it, so a median split returns half the nodes by construction
  — **9,988 to 10,016 of 20,002 for all thirteen biosamples, a span of 1.00
  against a depth span of 6.81.** Withdrawn as a between-biosample reading. And
  **DNase depth was the wrong covariate for the mark-derived readings**:
  `genes_poised` reads 0.38 against DNase and looks clean, 0.49 against
  H3K27me3, the mark that actually calls it — so a biosample with an
  under-called broad mark reads as un-poised, which is exactly what testis does.
  **The repair was registered with three ways to be wrong before it was built,
  and all three passed.** The depth residual replicates across disjoint halves
  of the genome at **rho 0.9835** (p 0.0001 over 10,000 shuffles), is not peak
  width (−0.033) and is not curvature (a log fit raises replication to 0.9945
  and drops the residual's own correlation with depth to 0.1374). The rate reads
  −0.522 against peak count, which clears the failure bar and **misses its own
  target band**, reported as the partial result it is.
  **The registered bar was "must move a published conclusion or it is
  decoration", and it was cleared by a reversal**: K562 reads 193,255 active
  enhancers against HepG2's 124,478, and **per 100k peaks K562 reads 37,272
  against HepG2's 70,472** — the published comparison points the other way.
  `enhancers_active_per_100k_peaks` now ships beside the count, and the CLI, the
  Cells view and the Progress tables carry the rate instead of the count.
  **Four published conclusions softened in the open**, including one from a
  commit six hours old: testis's 284 poised genes had been offered as evidence
  for the tissue-mixture argument two paragraphs after the same section recorded
  testis's H3K27me3 as under-called — reading one artefact twice — so testis is
  withdrawn from that claim and ovary carries it alone.
  **The limit is stated rather than glossed**: replicating across the genome and
  surviving peak count, width and curvature makes the residual a stable property
  of the biosample, and does not make it biology. FRiP, fragment length,
  crosslinking and library complexity are all biosample-stable assay properties
  and none is on disk. **Closed on 2026-09-22** (`1e57e3e` registration, `c6fc229` merge, `f21365f` the field and its
views); this row said "Still open" until 2026-09-27. There is now one definition,
`node_open_threshold()` in `genome/reader.py`, called from both sites. Against its registration:
**0 disagreements over 260,026 node calls** in 312 files, and **69 of 69** published candidates
identical on both fields. The verdict is the negative one: **a maintenance hazard and not a defect**,
because both fields it feeds are write-only, and because the construction that makes `nodes_open`
useless as a *count* is depth-invariant as a *membership*.
- **The node containment claim, audited and then held against a measurement it had never met
  (2026-09-27, lane-node, `05b71ae` audit and registration, `2b77663` measured arm, `f853fc6`,
  `698d10d`, 0 requests).** The "+2.88 points over random boundaries" was quoted seven times across
  three documents as settled shorthand, with no interval and **a baseline defined in no document**. It
  is as many uniform positions as the caller has post-merge edges, matched on boundary count per
  chromosome and on nothing else — **not on the 50 kb floor every real node passes**, so the control
  contains short nodes the real set cannot, and is biased low by its own geometry. On the same 440,377
  pairs, four defensible baselines give **+1.21, +2.90, +5.72 and +6.58**: a five-fold range, the
  published figure the least matched, the sign positive under all four. The interval, by a paired
  bootstrap over the 24 chromosomes because 440,377 pairs sit on about 20,000 nodes and are not
  440,377 observations: **+2.90, 95% CI +2.03 to +3.81**, confirmed to a tenth of a point by a
  node-cluster resample.
  **Four of the claim's digits could not be reproduced from the files it cites.** +2.88 is a 3-dp
  figure from one file minus a 4-dp figure from another; one of those files still carries a superseded
  control (0.791, from the 113,399-element archive), against which the node *loses* by 3.7 points.
  "18 of 24" counts chromosomes above **0.700, a threshold in no file, script or document**, and "sign
  test one-sided p 0.011" is exactly that count's binomial; against the 0.725 the same sentence names,
  it is **16 of 24, p 0.076 — not significant**. And the "87% of matched random windows" quoted beside
  the 63.7% is 52 of 60 windows in the twelve-locus panel, scored by a different layer, and was never a
  control for that rate. The denominators, by contrast, reconcile exactly: 961,227 rows → 612,323 name
  a gene (the 63.7%) → 593,765 in the annotation → 440,377 protein-coding, the same set as
  ATTRIBUTION.md's "compiled elements that state a confidence".
  **Then held against the 661 CRISPRi-measured regulated pairs already on disk, pre-registered with its
  power arithmetic**: at n = 594 the modelled +2.90 is 1.58 sigma, the minimum detectable excess 4.56
  points, and "undecidable" was registered as the most likely verdict. **Measured excess: +5.89, 95% CI
  +3.18 to +8.46**, larger than modelled under all four baselines. The registration's likeliest outcome
  was wrong, in the direction it had registered as plausible. Not clean, and recorded as such: 89.1% of
  the pairs are K562, and an unpaired node-cluster resample (+0.90 to +10.49) clears zero but not +2.9.
  **Where the increase comes from is the finding.** The measured share (0.7579) sits within 0.4 points
  of the modelled one; the whole difference is in the *control* (0.699 against 0.7252), because
  CRISPRi pairs are long-range and a random boundary cuts a long link far more often than a short one.
  So the claim the data supports is stronger than the one that was quoted, and differently worded:
  **nodes hold together enhancer-gene links at distances where randomly placed boundaries would
  separate them.** About 1,500 measured pairs would give 80% power at +2.9 points. NODES-READER-WRITER.md,
  the random-boundary audit (2026-09-27), which is the first time that document mentions the reading
  at all.
- **Direction, held against measured perturbations, and it passes (2026-09-22,
  lane-crispri).** Every enhancer-gene call this project makes carries an
  `action`, `activates` or `represses`, and that word is only the sign of a
  predicted log2 fold change. It had never been measured — §21 of the locus
  benchmark registered it **undecidable at n = 2**. The measurement was already
  on disk and unused: the ENCODE CRISPRi benchmark carries a **signed**
  `EffectSize`, and `attribution/crispri.py` reads `Regulated`, which is defined
  as `Significant AND EffectSize < 0`, so **the entire upward half of the
  measured signal was discarded before that module saw it.** Read two-sided
  against the finished genome-wide sweep for **0 requests** (`c115bbf` the
  registration, `5c842ca` the result): on the held-out arm the sign agrees
  **41 of 44 times, 0.9318, 95% [0.8177, 0.9765]** — K562 33/36, GM12878 8/8.
  **The registration expected this to be vacuous and said so in advance**, which
  is what makes the pass worth anything: every answerable pair carries a measured
  decrease, so a model that always says "down" would score 1.000, and the
  registered baseline was therefore not 0.5 but the sweep's own marginal
  down-rate, uncomputed at the time. It is **0.5361** — 325,989 of 610,034 K562
  elements — so the sweep is not a constant-down caller and both pre-registered
  margins are met. **A positive was then given a negative's scrutiny**: the 44
  are selected for being the top predicted target, which selects large effects,
  and the sweep's down-rate climbs 0.51 to 0.91 with magnitude, so the flat
  baseline was too easy; re-weighted to the scored pairs' own magnitudes it is
  **0.6537**, still cleared by 0.278 with the lower bound above it. Labelled
  post-hoc and reported beside the registered verdict, never in place of it.
  **All three errors sit below ‖log2fc‖ 0.055, and above 0.1 the model is 28 of
  28** — so the direction is worth reading in proportion to the magnitude it is
  read from, which is a sharper claim than the rate.
  **Two limits carry as much as the result.** Zero upward measured pairs are
  answerable, so **the `represses` half is untested**, and the only three up-calls
  the model made were all wrong. And the answerable set is the close half: the 111
  covered pairs that could not be asked about have a median element-to-TSS
  distance of **102,522 bp against 22,553**, and 37 of them are the upward pairs.
  **Closing that arm costs 30 requests, 35 with GM12878** — costed in advance,
  and a partial draw was refused on the ground that an arm sampled for its sign
  and scored in part is not the registered set. Authorised by the coordinator on
  2026-09-22 as the next spend, because it is the arm that could falsify a
  passing result. NODES-READER-WRITER.md, "Measured perturbations, and the
  direction question"; `data/results/crispri_direction.json`.
- **The upward arm, and the direction pass does not survive it (2026-09-22,
  lane-updir).** `c2a16d8` the registration, `d492834` the result, **0 requests
  of the 35 authorised** — because the 30/35 costing was reproduced exactly and
  then turned out to buy something already on disk. **"The sweep stores one
  target per element" is true of a derived table, not of the measurement**: the
  compact `all_elements/<chrom>.json` keeps one gene, while the per-element
  response cache written by the same run at `threshold=0.0` carries **every gene
  in the scorer's 1 Mb window with a signed per-cell log2fc, uncensored**. The
  answerable set goes from **44 one-signed pairs to 152 two-signed** ones, 116
  down and 36 up.
  **With both signs present the primary had to change, and the registration said
  so before scoring**: the set is not balanced, so a constant-down caller scores
  0.763 on raw agreement and quoting that against 0.5 would have manufactured a
  success one step further along. The registered primary is therefore **balanced
  accuracy**, where chance is 0.5 for any class mix and a constant-sign caller
  scores exactly 0.5 by construction.
  **THE UPWARD ARM FAILS: 15 of 36, 0.4167 [0.2714, 0.5780]** — not clearing
  chance, and below the model's own unconditional up-rate of 0.4643, so on pairs
  where an increase was measured the layer does no better than its own habit.
  Downward 96/116. **Combined balanced accuracy 0.6221 [0.5354, 0.7088]:
  undecidable, as the registration predicted to within a decimal.**
  **The sharpest number is raw agreement 0.7303 against a majority-class rate of
  0.7632 — the layer is worse at calling direction than a caller that says "down"
  every time and never looks at the element.** Its marginals are almost exactly
  right (predicted 117 down / 35 up against measured 116 / 36) and its pairings
  are wrong. **This is precisely what the first half could not have found**: on an
  all-down set a caller worse than constant-down still scores 0.93.
  **So the 41/44 pass of `5c842ca` is correct for what it measured and is now a
  statement about downward effects only**, annotated in place rather than
  rewritten, and no subset was searched for where the upward arm holds. All three
  registered falsifiers fired, including top-target selection: 42/44 on pairs the
  compact table could see against **69/108, balanced 0.5833**, on pairs only the
  cache can answer. **"All errors below ‖0.055‖" is falsified** by a 0.3172 miss
  on HBG1 in K562 against a measured -0.7356; what survives and sharpens is that
  in the weak band the sensitivities are 0.7069 down and 0.3939 up — **balanced
  accuracy 0.5504, chance** — so all of the layer's direction signal lives in
  large-magnitude downward calls. The registered magnitude-matched sub-test has
  only 2 upward pairs above 0.1 and is reported as **unable to settle** whether
  the upward failure is about sign or magnitude, rather than leaned on.
  **One registered check did not pass and was traced rather than waived**: one of
  the 44 changes sign, and it is a candidate-set difference rather than a reader
  disagreement — three scored elements overlap that perturbation and the wider
  set resolves the largest-magnitude rule elsewhere. It is recorded as a flaw in
  the lane's own registration, and reported because it moves the other way: it
  would turn one of the first half's three errors into a hit. **The general form
  is carried forward: "largest absolute change among overlapping matches" is not
  well defined unless the candidate set is pinned too.**
  **The `represses` half of `action` is now measured and unsupported, which is a
  better place to stand than measured nowhere.** NODES-READER-WRITER.md;
  `data/results/crispri_direction_both.json`.
- **Next (epigenome and nodes).** 1. Done as an alternative; the next test is a
  stricter site call (best-hit strand, stronger motifs) judged on all four
  measurements at once, and whether Hi-C questions and enhancer-to-gene
  questions need two node sets. *(Closed 2026-09-21 — see the bullet above.)*
  2. The poised call in the Blocks lane:
  **done 2026-09-17**. Poised genes arrive inside `silent_genes` on purpose, because every
  consumer reads "absent from that list" as read and leaving them out would report a held
  gene as expressed; the lane therefore drew all of them as silent. It now reads the
  narrower `poised_genes` first and draws three states: read (filled green), poised (amber)
  and silent (hollow). In one 4 Mb window of chr21 that is 23 genes shown as shut when they
  are held ready. `decompile` and `report` carried the same defect and are fixed with it: a gene's
  dossier now lists `poised_in` beside `read_in` and `silent_in`, and the decompiled line
  says "poised in K562" where it previously dropped that cell type from BOTH groups, which
  reads as "never measured" rather than as a claim. The read-by-marks call is still shown
  nowhere. 3. A neural, gonadal
  or embryonic biosample beyond SK-N-SH and H1. 4. Measured perturbations
  (CRISPRi in K562) for the direction finding.
- **Owner.** genomeos-e1 holds areas B, C and D since 2026-09-14, when the
  previous owner's session ended (genomeos-8e; genomeos-fe before the restart
  of 2026-09-12).

### C. Molecules (RNA, proteins, pathways, the knowledge graph)

- **Goal.** From gene to transcript to protein isoform to modified state,
  compiled from public sources into one definition per protein, verified
  against the curators, and executable as pathways.
- **Code.** `molecules/compiler.py`, `proteome.py`, `verify.py`, `rna.py`,
  `graph.py`, `reactome.py`, `uniprot.py`, `alphafold.py`, `ptm.py`;
  `flow/trace.py`.
- **Design.** PROTEIN.md, FLOW.md.
- **Data.** Counted on disk 2026-09-27. This row read "`proteome_chr*` (4 of 25: chrM, 21, 22, Y),
  `translation_vs_uniprot_chr*` (3) … `data/knowledge/proteins` (18 MB)" and understated every
  figure in it for sixteen days. `proteome_chr*` **25 of 25** plus `proteome_genome_wide`;
  `translation_vs_uniprot_chr*` **25 of 25**, with no genome-wide file of its own; `graph_chr21`,
  `graph_genome`, `ptm_genome_wide`; compiled definitions in `data/knowledge/proteins`, **19,453
  files and 291 MB**. A row that understates its own area is read as a to-do list by the next lane,
  which is how the same work gets scheduled twice.
- **Requirements.** UniProt accession is the id; Gene → Transcript →
  Protein isoform, never Gene → Protein; predicted never equals
  experimental; STRING association is not interaction; never call an API
  per simulation tick.
- **Done (2026-09-11).** The whole human proteome: 19,478 coding genes over
  25 chromosomes, sequence 99.1%, domains 99.0%, function 86.3%, interactions
  81.5%, pathways 58.2%, experimental structure 45.2%, predicted 98.2%,
  disease 25.5%. Translation verified against UniProt on 19,249 genes: 90.9%
  identical to the canonical entry, 97.8% exact for some isoform, 99 real
  disagreements left to explain. Knowledge graph genome-wide: 41,982 nodes,
  327,024 edges, largest component 15,165 of 19,283 compiled proteins.
- **Packaged as a library (2026-09-11).** The compiled proteome is distilled
  into `genomeos/lib/data/proteome.json.gz`: 19,283 proteins in 2.3 MB,
  shipped inside the package, with function on 87.0%, pathways 58.7%,
  partners 63.6%, experimental structure 45.6% and disease 25.7%, and the 168
  symbol mismatches marked rather than hidden. `genomeos protein X --lib`
  answers with no network at all, which is what makes the protein layer usable
  on a laptop and in CI.
- **Kinetics where a curated model exists (2026-09-11).** `genomeos pathway
  ID --kinetic` finds a curated ODE model by the Reactome pathway's name, runs
  it on the in-house SBML engine and compares a knockout against the baseline
  on final levels and peaks; a gene symbol reaches species named `x1`/`x2`
  through the model's MIRIAM annotations. Two worked results: taking RAF1 out
  of Hornberg2005 erases the downstream transients entirely (peaks 0.201 and
  0.562 to zero), and taking MEK out of Kholodenko2000 stops the ERK
  oscillation (Erk2-PP peak 298.8 to 10.0, final 37.0 to 0). The engine gained
  function definitions, rate rules, species names and annotations, and
  amount-versus-concentration semantics.
- **A recorded failure, kept on purpose.** Schoeberl2002 (100 species, 125
  reactions) does not reproduce on the in-house engine: the cascade species
  stay flat. The result file is committed next to the two that work. That is
  what keeps the libRoadRunner adapter (§5 item 12) an honest open item rather
  than an aspiration, and it marks the size of model where our engine stops.
- **Post-translational state (2026-09-11).** The layer the protein model
  names (`ProteinState`) is populated from UniProt's modified-residue
  features: `genomeos ptm X` lists a protein's modifiable sites with class
  and writer, `genomeos ptm --writer PKA` lists a writer's substrates, and the
  genome-wide summary counts 96,362 sites on 13,083 proteins with 352 named
  writers. The 2,985 writer to substrate pairs are `modifies` edges in the
  knowledge graph, rebuilt to 330,018 edges. A site that can be modified is
  not a measurement that it is, and the result says which it holds.
- **Isoform-level expression (2026-09-11).** Expression had been recorded per
  gene; `genomeos rna X --isoforms` reads GTEx v8 median TPM per transcript
  per tissue (paged from the portal, cached under
  `data/knowledge/expression`), labels each transcript with the compiled
  definition's name, canonical flag and protein length, and names the
  dominant isoform per tissue against the canonical one. TP53-201 dominates
  in all 54 tissues, so there the canonical choice is the expression. APP-201
  (770 aa) dominates in 33 tissues, APP-204 (751 aa) in 14, and APP-202, the
  695-residue neuronal isoform, in the cerebellum at 79% of the gene's TPM,
  which is the textbook. The layer says which protein a tissue makes, which
  is what `ProteinState.isoform` was for. The Molecules protein card shows
  the modification classes and the writers UniProt names.
- **Modification observation, partly unblocked (2026-09-27, lane-phospho, `58f84fb` registered, `869158d` result).** Occupancy stays
  blocked, but observation now exists: **29,279 of the 41,661 curated phosphosites (70.3%**, registered
  45–70%, a miss high by 0.3 points, range not moved; 30.4% of all 96,362 curated sites) carry
  `observed_in_cell_types_or_tissues` (1–83, median 12) from Ochoa et al. 2020's reanalysis of 6,801
  PRIDE runs at 1% site FDR, via the authors' funscoR data (LGPL; its PhosphoSitePlus-derived table
  excluded). No per-site experiment count is published, so the field is not "observed in N
  experiments". Checks pass (NPM1 S125 present, TP53 M1 absent, nothing off S/T/Y, 0.07% residue
  disagreement); curation and reference share source studies, so the 70% partly agrees with itself.
  CPTAC through cBioPortal is reachable (12 phosphoprotein profiles, ODbL, log2 ratios to a pooled
  reference) and usable only as `relative_abundance_in_tumours`, not built: it is keyed by RefSeq
  position and needs a RefSeq-to-UniProt residue mapping first. docs/PROTEIN.md has the detail.
- **Tumour abundance joined (2026-09-28, lane-cptac, `cae8d7d` registered, `5938b5d` result).** CPTAC
  phosphoproteomics through cBioPortal, mapped RefSeq → UniProt (10,644 of 12,025 proteins; 114,236
  single localised sites carried by identical sequence or a unique 15-residue window, residue letter
  checked, 0.02% mismatch) and joined to the curated sites as `relative_abundance_in_tumours` (per
  study: tumours with a value and the median log2 ratio to the pooled reference; never occupancy, never
  per patient; the result is an ODbL 1.0 derived database with its attribution inside). **19,104
  curated phospho sites on 5,105 proteins (45.9%)** from lung adenocarcinoma, glioblastoma and paediatric
  brain; EGFR Y1092 present, TP53 M1 absent, nothing off S/T/Y. **Negatives:** the pancreatic 2021
  profile holds log2 intensities, not ratios (median 18.8), and is excluded — the registered pooled
  check passed with it inside, so that check was too weak and a per-profile check caught it; coverage
  missed high on both registered ranges (55.3% against 35–55%, 68.6% against 45–65%), reported as a miss;
  the seven gene-keyed profiles stay uncommitted at 14.3% residue mismatch until their isoforms are
  known. Occupancy stays blocked. docs/PROTEIN.md, "Relative abundance in tumours (2026-09-28)".
- **Missing, and blocked on data (checked 2026-09-12).** Modification state
  is possibility, not occupancy: measured phosphoproteomics has no open
  per-site, per-tissue table the project could stream and distil as it does
  ClinVar (PhosphoSitePlus needs registration, PRIDE holds raw or per-study
  tables), so it is recorded as blocked rather than attempted. The person's
  own expression: HG002's cell line GM24385 has no ENCODE experiment, and the
  GIAB direct-RNA nanopore runs on ENA are raw reads with no aligner in the
  zero-dependency core, so `individual protein` states GTEx's median and says
  why; the day a processed signal track for the person exists,
  `rna_measured.py` reads it as it reads ENCODE's. Isoform expression is
  GTEx's adult tissues, not a cell type or a stage; kinetics only where
  BioModels has a curated model, which is a small fraction of Reactome.
- **The disagreements, read with the isoform the body makes (2026-09-11).**
  `scripts/disagreements_isoforms.py` took GTEx's dominant transcript for each
  of the 96 disagreement genes, translated it locally and compared it with
  UniProt (`translation_disagreements_isoforms`, 151 s). Not one is a
  canonical-choice artefact: 38 make the canonical transcript and still differ
  (reference alleles), 18 make another isoform that differs too, 16 have a
  non-coding dominant transcript in the local models, and 24 have no expressed
  isoform in GTEx or are absent from it, olfactory receptors mostly. The
  disagreements are reference alleles and gene-model differences, not isoform
  choice, cross-tabulated against the mechanism triage in the file. The item
  closes.
- **Writers as rules (2026-09-11).** `genomeos protein X --bio` and the
  packaged `--lib --bio` block end with one `rule WRITER modifies GENE`
  per UniProt-named writer (strength 1.0, evidence curated "UniProt: n
  modified residues written by WRITER", confidence 0.8), the writers declared
  first as bare proteins because a rule may only name a declared entity; a
  program that imports the writer itself drops the stub. The packaged
  proteome carries a writers field per protein (2.33 MB). A kinase knockout
  now reaches its substrates through the rules rather than the neighbourhood.
- **The disagreements held against three real genotypes (2026-09-11).**
  `scripts/disagreements_genotype.py` applies the SNVs HG002, HG003 and HG004
  carry inside each disagreement gene's canonical CDS, translates and compares
  with UniProt again (`translation_disagreements_genotype`). HG002: 40 genes
  with no variant inside the CDS, 52 whose variants leave the disagreement as
  it is or worsen it, 1 with only indels (not applied), 0 restored; HG003
  46/46/1/0, HG004 38/54/1/0. The 96 are not common alleles this trio carries
  the other way. PROTEIN.md now bounds them from three sides: not isoform
  choice, not this trio's alleles, mechanism by triage; population allele
  frequency at the exact site is what would settle each one.
- **The twin's isoform (2026-09-11).** `genomeos individual protein --name ME
  --gene G --chrom C` traces GTEx's dominant transcript in each tissue and
  reads the person's coding variants on that transcript rather than the
  canonical one: which protein each tissue makes in this person. HG002 and
  APP: APP-201 (770 aa) in 33 tissues, APP-204 (751 aa) in 14, APP-202 (695
  aa, neuronal) in six brain regions, APP-203 in the cortex, no coding variant
  on any. HG002 and FUT2: FUT2-201 in 49 tissues carrying p.Trp154Ter
  homozygous, the truncated non-secretor protein everywhere the gene is made.
  The output states that the isoform per tissue is GTEx's population median,
  not the person's own expression; that half stays open. DATA.md "the twin's
  isoform".
- **The disagreements settled by allele frequency (2026-09-11).**
  `scripts/disagreements_frequency.py` takes the 38 same-length disagreements,
  derives for each differing residue the single-base change that would give
  UniProt's residue, and looks that allele up in gnomAD and dbSNP through VEP
  (`translation_disagreements_frequency`, 112 s). Two genes where the curated
  allele is the common one, so hg38 carries a minor allele (OR9H1 rs7555046 at
  0.99; FCGBP, four sites at 1.0); six known polymorphisms where UniProt chose
  one haplotype (HLA-DQA1, MICA, GSTT2, SERPINA2, PRB4, SAMD1); ten in dbSNP
  but rare or unmeasured; sixteen unknown to both, so gene-model or entry
  differences; four with no residue-level difference. PROTEIN.md closes the
  series: the 96 read from four sides are 31 hg38 frameshift or nonsense
  alleles, 8 residue-level population variants, 7 different products, and the
  rest gene-model differences. The item closes.
- **Next.** 1. The person's own expression in place of GTEx's median, blocked
  until a processed track for the person exists (see Missing). 2. Measured
  modification state, blocked on an open per-site table (see Missing). 3. Done
  2026-09-12 in area D: the 31 hg38 frameshift and nonsense alleles surfaced
  in the coding inventory and the dossier, a carrier of the reference allele
  being a carrier of a truncation.
- **Owner.** genomeos-8e (genomeos-fe before the restart of 2026-09-12).

### D. The individual (BioTwin)

- **Goal.** A person's genome plus measured state plus environment, forked,
  run and compared.
- **Code.** `twin/twin.py`, `twin/clocks.py`, `genome/telomere.py`,
  `genome/variants.py`, `genome/individuals.py`, `runtime/variant_effect.py`,
  `calibrate`.
- **Design.** ACTION-PLAN.md Phase 3, STORAGE.md, LESSONS.md (timers).
- **Data.** `hg002_chr21`, `hg002_telomere_stream`,
  `hg002_methylation_stream`, `clock_GSE41169`, `clinvar_chr21_agreement`,
  `calibration_hematopoietic_stem_attrition`.
- **Requirements.** HG002 is the test human; measured state is labelled
  with its caveats (cell-line signature); every inferred parameter carries
  its confidence.
- **Done (2026-09-11).** The twin covers all 22 autosomes, and a person's own
  genome is a first-class input: `genomeos individual import FILE.vcf` splits
  any GRCh38 VCF into per-chromosome PASS files under a git-ignored directory,
  and carrier lookup, the gene report, twin build and the gene-by-gene walk
  all see it. Nothing leaves the machine, which is the only acceptable design
  for somebody's own genome.
- **All four AlphaGenome features are built** (a–d): variant effect per
  tissue, enhancer to gene, splice sites in the segment parser, and the
  predicted effect of a person's own regulatory variants on one gene, summed
  per haplotype where the genotype is phased.
- **Privacy is enforced by the repository, not by intent.** A person's
  imported genome (`data/individuals`), their per-variant predictions
  (`data/knowledge/alphagenome`) and now their twin (`data/twins`, which holds
  chronological age, telomere length and epigenetic age) are all git-ignored;
  only the open-consent public test subject stays tracked. Per-person results
  are never written under `data/results`.
- **What a person learns from their own file (2026-09-11).** Three commands,
  all reading the imported genome and writing under the person's own
  directory, never `data/results`. `genomeos individual screen` streams
  ClinVar once (4.47 M rows, the 346,571 pathogenic and likely pathogenic
  rows kept locally) and intersects the person's files: HG002's four million
  variants in seven seconds give two hits, the F11 factor XI deficiency
  carrier allele among them, each with genotype, zygosity, review stars and
  conditions. `genomeos individual knockouts` lists nonsense, start-lost and
  stop-lost SNVs on canonical transcripts across every chromosome, homozygous
  first, with the fraction of protein lost: HG002 has 88 in 84 genes, the
  FUT2 non-secretor allele homozygous among them; calls past a stop the
  reference itself carries are left out and counted. `genomeos individual
  report` writes the dossier as one Markdown page from whatever has been
  computed (import, reference checks, screen, truncating variants), and a
  section not run says so. Each has a Twin tab button. Research annotation,
  not a clinical report, and the page says that too.
- **The coding inventory (2026-09-11).** `genomeos individual coding` counts
  every coding SNV on canonical transcripts by consequence and by gene,
  protein-changing and homozygous first, as a dossier section and a Twin tab
  button. HG002 on 22 autosomes: 21,019 coding SNVs, of which 11,101
  synonymous, 9,830 missense, 61 nonsense, 16 stop lost and 11 start lost;
  5,438 genes carry a protein-changing variant and 2,529 a homozygous one,
  MUC16 and the HLA genes on top, as length and polymorphism predict. Calls
  past a stop the reference carries are left out, as in the knockout scan.
  DATA.md "Your own genome" describes it.
- **Missense ranked by the residue (2026-09-11).** The inventory ranks a
  person's missense variants by what UniProt records at the residue: an
  annotated site first (active or binding site, modified residue,
  glycosylation, the two cysteines of a disulfide bridge, a motif), then inside
  a domain, then nothing, homozygous before heterozygous; bridges and
  cross-links count only at their two residues. HG002: 46 of 9,830 missense
  variants sit on an annotated site and 5,374 more inside a domain; CD52
  p.Asn40Ser homozygous removes an N-glycosylation site, GALNTL5 p.Cys124Arg
  and OR56B1 p.Cys106Arg each lose a bridge cysteine. In the CLI table, the
  dossier and the person's coding.json. A rank from annotation, not a
  predicted effect.
- **The first pedigree (2026-09-11).** `individual import` streams from a URL,
  so HG003 and HG004 arrive as local individuals from the GIAB v4.2.1
  benchmarks; `individual regions --name X BED` keeps the regions a person's
  calls are trusted in; `individual trio --name HG002 --father HG003 --mother
  HG004` reads inheritance per chromosome. 4,096,122 child variants, 96.4%
  inherited, 2,336,470 from both parents. Without regions, 100,114 calls are
  absent from both parents and 48,848 are Mendelian errors; with the three
  trusted region sets that falls to 12,296 and 2,930, with 133,736 child calls
  counted as untrusted rather than as events. DATA.md says plainly that the
  12,296 are mostly representation differences between three call sets, since
  true de novo variants number 60 to 100 per child: the trio measures caller
  agreement first. Stored under the child's directory.
- **The reference's own truncations (2026-09-12).** For the 31 genes where
  hg38 carries a frameshift or nonsense allele against the curated protein
  (the verified disagreements by mechanism), `individual coding` and the
  dossier say whether the person differs from the reference anywhere inside
  the gene; no variant means the person carries hg38's truncation, which no
  caller lists. HG002 matches the reference in 5 of the 30 on file (OR2T7,
  OR4C45, OR4K3, OR1P1, SCYGR10) and differs inside the other 25, none of
  which restores the curated protein. `reference_alleles_carried()` in
  `genome/individuals.py`; DATA.md paragraph.
- **The missense effect predicted (2026-09-12).** `genomeos individual coding
  --name N --predict` streams DeepMind's AlphaMissense hg38 table once (643
  MB, 71.7 million rows, 37 s from the public bucket) and keeps only that
  person's missense variants under the person's git-ignored directory
  (`genome/missense.py`; a second call streams nothing). Each variant gets the
  score on the person's transcript where the table has it, else the highest
  across transcripts, the model's class (likely pathogenic at 0.564 or more,
  likely benign at 0.34 or less), confidence capped at 0.7, evidence
  `predicted`. HG002: 9,830 missense SNVs on canonical transcripts, 8,933
  scored, 188 likely pathogenic (37 homozygous), 254 ambiguous, 8,659 likely
  benign, 729 unscored (indels, transcripts the table lacks). The inventory,
  the dossier and the Twin card list them strongest first beside UniProt's
  site annotation, which stays the ranking where no score exists. The table is
  CC BY-NC-SA 4.0, so the scores never enter a committed result; the web
  endpoint reads the person's cache only and the stream is a CLI step.
- **Telomere from streamed reads (2026-09-12).** CRAM cannot be decoded
  without htslib, so the streamed range reads an indexed BAM instead:
  `genome/bam_range.py` reads a coordinate-sorted BAM behind a URL with the
  standard library (BGZF blocks, the .bai's linear index, its pseudo-bins for
  the mapped and unmapped counts, records), and `genomeos telomere <bam-url>
  --bai <index> --save NAME` runs TelSeq's arithmetic without downloading the
  file: the index gives the read totals, a 100 MB sample of the unmapped tail
  (where the pure TTAGGG reads sit) is scaled by the index's unmapped count,
  and the first and last 10 kb of assembled sequence of each chromosome, past
  the terminal N runs, are read whole for the boundary reads. HG002 on GIAB's
  300x GRCh38 BAM (601 GB): 6,182,841,975 reads in the index (372,318,760
  unmapped); 910,358 unmapped reads sampled, 1,130 with seven TTAGGG repeats;
  44 of 46 chromosome ends carry reads at mean coverage 203x (chr5 p and chrX q
  are telomeric repeat in the assembly itself); telomere about 2,676 bp per
  end, from 218 MB in 48 range requests, inferred at 0.3 with no GC
  normalisation, so a number to compare between people read the same way, not
  to quote (`telomere_range_HG002`, GIAB open consent).
- **Missing.** Polygenic scores are raw sums without a population
  distribution to place them on; the trace still uses the canonical
  transcript rather than the one the tissue makes; the trio's 1,430
  normalised candidates are not yet phased by parent or reduced to the 60 to
  100 a child really carries; the telomere estimate has no GC normalisation.
- **The trio normalised (2026-09-12).** The comparison normalises every
  allele before matching (common suffix and prefix trimmed to one anchor base,
  indels left-aligned along the local reference; `normalise_variant` in
  `genome/individuals.py`), because the three GIAB files write one insertion
  three ways. HG002 against HG003 and HG004 over 22 autosomes with all three
  region sets applied: de novo candidates 12,296 to 1,430 (1,173 SNVs, 257
  indels), Mendelian errors 2,930 to 7, inherited 96.4% to 96.9%, 124,881
  child calls outside a parent's region. What remains is child PASS calls
  inside both parents' high-confidence regions with no parent row within a
  base (97% of a 300-candidate sample): the true de novos, 60 to 100 expected,
  plus the benchmark files' own disagreements. `individual trio` prints the
  SNV and indel split and the representation used. The three items area D's
  Next held this morning are closed; ARCHITECTURE.md §7 now describes the two
  standard-library range readers (`bigwig.py`, `bam_range.py`) and the rule
  that only the distilled answer is kept.
- **Genome diff as code review (2026-09-12, from the pool).** `genomeos
  individual diff --name A --against B|reference` (`genome/diff.py`) renders
  two people's protein-changing variants gene by gene as a review: alleles
  normalised, `+` only A, `-` only B, the shared count per gene, each line
  with its consequence on the canonical transcript, zygosity, AlphaMissense
  class where the person is scored and the ClinVar row where screened; genes
  ordered ClinVar first, then truncations, then predicted class. The coding
  inventory keeps every protein-changing variant under the person for it, the
  Twin card gains a "compare with" control, and the Markdown goes to the
  terminal or a file, never under `data/results`. HG002 against HG003 over 22
  autosomes: 9,745 and 9,649 protein-changing variants, 7,513 shared, 2,857
  genes differ; only HG002 2,232 (23 truncating, 75 predicted likely
  pathogenic), only HG003 2,136 (30 truncating). HG002 against the reference:
  5,344 genes carry a change. This is the pool bullet's diff at the
  protein-changing level; the regulatory level (a person's elements and their
  predicted effects) is the next layer and not done.
- **Polygenic scores (2026-09-12).** `genomeos individual pgs --name N
  [--score ID ...]` (`genome/polygenic.py`) streams a PGS Catalog score's
  harmonised GRCh38 weight table once (cached git-ignored under
  `data/knowledge/pgs`, metadata from the catalog's REST: trait, publication,
  ancestries, licence) and sums effect-allele dosages over the person's
  genotypes: called variants as called; a no-call inside the person's trusted
  regions read as homozygous reference against the local FASTA; positions
  outside the regions missing and counted; allele mismatches counted. Weights
  `curated`, the sum `derived` at 0.3; a raw score with coverage, not a
  percentile, since the catalog publishes no population distribution, placed
  against the other local people with a within-group z and the ancestry caveat
  in every result. Defaults: breast cancer PGS000004 (313 variants), LDL
  PGS000115 (223), height PGS000297 (3,290), coronary artery disease
  PGS000018 (1.7 million). Results under the person, never under
  `data/results`; `/api/individual/pgs` and a Twin card button read them. The
  trio at 90 to 99% coverage: LDL HG003 +0.93, HG002 +0.77, HG004 +0.61;
  height 55.0, 54.5, 53.3; coronary disease -1.14, -0.82, -0.48; breast cancer
  +0.83, +0.85, +0.75. The child sits between the parents on three of four,
  which is what a correct sum should do and the only check available without
  a reference population.
- **The diff at the regulatory level (2026-09-12).** `genomeos individual
  diff --name A --against B|reference --regulatory --chrom C [--constraint]`
  (`genome/regdiff.py`) places both people's normalised variants inside the
  chromosome's ENCODE elements and writes the difference per target gene: the
  element's class, the CTCF node's inferred target with basis and distance,
  the AlphaGenome deletion target where the element has been scored, and with
  `--constraint` the base's own phyloP read by range; constrained and
  predicted lines first. Served by `/api/individual/diff` with a "regulatory
  (one chromosome)" checkbox on the Twin card. HG002 against HG003 on chr21:
  5,730 and 5,579 variants inside the chromosome's 13,356 elements, 4,352
  shared; only HG002 1,378 (14 in an element the deletion job scored, 27 on a
  base with phyloP at 2.27 or more), only HG003 1,227 (18 scored, 19
  constrained). The page says what it is: a variant in an element is a
  candidate, not an effect; the effect per variant is the self-hosted model's
  job (area I).
- **Next.** 1. Phasing by parent on the 1,430 candidates, which waits for a
  phased import (the GIAB files carry unphased genotypes); until then the
  1,430 stand as the trio's disagreement set.
  **Done 2026-09-27 (lane-q100, `457702c` registered, `4d63910` result), and not by the route
  expected.** GIAB v5.0q ships the Q100 v1.1 assembly's own GRCh38 dipcall VCF (40 MB, NIST public
  domain, acknowledgement only), labelled paternal|maternal and checked on the file: every chrX non-PAR
  call `.|1`, every chrY call `1|.`. Controls pass (99.9% found, 99.999% on the right parent). Of 1,411
  assessable candidates **1,382 sit on one parental chromosome, 703 paternal to 679 maternal (share
  0.509)**; 29 disagree with the assembly. **P1–P3 failed as registered** (they expected 40–150 germline
  de novos at a ~0.8 paternal share); P4 holds. The rival reading (missed parent calls) is contradicted
  too, unregistered: 3.0% of the candidate SNVs are seen in 88 HPRC haplotypes against 92.8% of
  inherited calls, and 69% lie in NIST's de novo and mosaic regions. So the candidates are real and
  new, most likely post-zygotic or cell-line mutations — an interpretation, not a measurement. A parent
  label says which chromosome carries an allele, not whether it was inherited. Result files hold counts
  only. 2. Germline versus post-zygotic per candidate, which needs read allele fractions. 3. The note in
  `trio()` still calls an excess "representation differences"; docs/DATA.md records that this is wrong
  after normalisation.
  **Read allele fractions, done 2026-09-28 (lane-vaf, `bf789a3` registered, `0125926` result).** GIAB
  v4.2.1's own pooled read depths (ADALL; 156 MB, NIST public domain, no BAM needed). **The 1,375
  parent-assigned candidates sit at a median read fraction of 0.442 against 0.495 for 1.28 million
  inherited heterozygous SNVs** (all checks pass; median depth 329). They form one shifted group, not a
  germline group at 0.5 plus a subclonal tail: the registered classes, which expected two groups,
  count 11.3% as low, mostly the lower edge of that one group. Read as one clone's mutations carried in
  about 88% of cells — the cell that founded the line, or its lineage in the donor — not germline.
  Exploratory: at most about 208 can be germline, and the 89 at or above 0.5 are 65% paternal (p ≈
  0.006), the excess germline de novos would show, while the rest split evenly. Limits: the benchmark
  excluded calls below 0.2 or above 0.8, and HG002 DNA is from a cell line. Still open: a per-candidate
  germline call, and telling early-embryonic from line-founder mutations, both of which need DNA that
  never went through culture.
- **Owner.** genomeos-8e (genomeos-fe before the restart of 2026-09-12).

### E. From one cell to an organism

- **Goal.** Grow an organism from one cell with mechanism where it is known
  and the observed program where it is not, timed by the genome's own
  timers, and say at each level how much is mechanism.
- **Code.** `runtime/body.py`, `organism/*`, `runtime/cell.py`, `grn.py`,
  `boolean.py`, `sbml.py`, `spatial.py`, `segmentation.py`,
  `gastrulation.py`, `compose.py`.
- **Design.** ORGANISM-FROM-ONE-CELL.md, BIOLANG-v0.3.md.
- **Data.** `celegans_lineage_cells` (2,183 cells), `celegans_packer2019`,
  `human_cell_turnover`; programs in `data/organisms/celegans` and
  `data/organisms/human`.
- **Requirements.** Every cell by name against Sulston; deaths and fates
  exact; timers with their measured spread; the human body reproduced by
  growth and turnover, not asserted; organism-level uncertainty populated.
- **Missing.** The human program is counts, not mechanism, except
  haematopoiesis; the worm's terminal fates below the founders were the
  observed lineage program, although every cell carries its measured
  transcription factors (Ma 2021 atlas, `express`) and factor rules now
  decide the embryonic fates they can (below); the PAR polarity rules are
  written since 2026-09-14 and Digital Development is 11 of 11;
  `commitment` and `competence` are in the language since 2026-09-14, and the
  worm's own program declares the first of them only since 2026-09-15, when a
  fate read from a growing integral turned out not to be stable without it;
  the external engines (libRoadRunner, MaBoSS, CompuCell3D)
  are deferred for lack of Python 3.14 wheels, although CI runs 3.12 and
  could test them as optional extras. Done 2026-09-11: haematopoiesis,
  the Body as a process-bigraph process gated by a network, space (sites,
  fields, division in place, contact inhibition, migration; the French flag
  grown from one cell), the measured reader, Packer 2019 at 84% on
  single-tissue ids (the rest is labelling depth), Digital Development as
  the published knockout set (7 of 11 modelled transformations then; 11 of 11
  since the PAR rules of 2026-09-14, below).
- **Done 2026-09-14.** Terminal fates from measured factors
  (`embryo_factors.bio`, `celegans_fate_rules`): with factor rules taking
  precedence over the lookup, 496 of 555 embryonic terminal fates (lookup
  555); textbook rules read as exposure integrated along the lineage 93.2%
  with the lookup as fallback, the instantaneous threshold 89.0% and worst
  at every threshold; held out by sublineage, factors decide 449 cells and
  are right 365 times where the sublineage majority is right 321; glia and
  coelomocytes are where factors do worse. Commitment, competence and
  lateral inhibition tested in the atlas (`celegans_commitment`): no
  sharpening of states or programme exclusion with time against shuffled
  time and curveball nulls, competence not separable from lineage history,
  sister divergence below the two-reporter measurement floor; all negative,
  recorded in ORGANISM-FROM-ONE-CELL.md with the construct proposal
  (`commitment`, `competence`, exposure reads, contacts) sent to genomeos-c1.
  Fate precedence explicit since genomeos-c1's `regime { fates: first }`:
  `fates.bio` carries a `priority` per rule, all 1,439 fates to 800 min are
  unchanged, 45 order-dependent decision points become 0 and 620 fewer
  decisions fire. **The founders decided by contact** (`founders_contacts.bio`,
  `contacts_embryo.tsv`, `embryo_contacts.bio`, `celegans_contacts`): the three
  founder signals name no sender and no receiver, each receiver reading the
  ligand summed over the cells touching it. 8 of 8 founder identities, 8 of 8
  published knockouts, and against Sulston **496 of 555 terminal fates, exactly
  the number before** — 1,438 of 1,439 cells identical in fate, terminal name and
  birth time, the odd one being EMS's own type. What contacts buy is not score
  but falsifiability: swapping ABa and ABp in the contact table alone swaps their
  fates, which is the blastomere rearrangement of Priess & Thomson 1987 and
  cannot be written with named senders. The table is curated geometry, not
  measured contact areas; no measured table is held locally.
  **AC/VU as an equivalence group** (`acvu.bio`, `contacts_acvu.tsv`,
  `celegans_acvu`): Z1.ppp and Z4.aaa identical, touching, and told apart only
  by seeded noise. Over 200 seeds, 199 give exactly one anchor cell and
  **Z1.ppp wins 96 of them, 48.2%**; with noise off 0 of 200 diverge and the two
  cells' Delta is identical to the last digit; creating each division's
  daughters in the opposite order gives the same winner in 200 of 200 seeds.
  Against the reference, which records one animal's outcome, the program is
  right 96 of 200 times, and **48% is the ceiling for an honest program**: the
  worm program names Z1.ppp and takes 100%, so this is the one terminal cell
  where the lineage score is the wrong instrument. Two ways the flip was
  silently loaded, both measured: a daughter inherits its mother's network
  state (list Z1.ppp as touching Z4.aa and Z1.ppp wins 100 of 100), and an
  instantaneous level read at an arbitrary time is not a fate (§7.2a's
  sustained reads are what is missing). Kimble 1981's ablations are run too:
  remove either cell and the survivor is the anchor cell in **100 of 100**
  seeds, remove both and there is **no anchor cell in any** of 100, remove a
  flanking cell and the split is untouched.
  **Glia from factors** (`celegans_glia`, `CITED_GLIA_RULES`): sheath glia
  partly, socket glia not at all. All 40 glia and all 40 sisters are in the
  atlas and no glial cell shares a factor set with a non-glial one, so the
  question is answerable. PROS-1/Prospero (Wallace et al. 2016) alone has
  precision 0.27; with the class factors NHR-25 and SOX-2 it claims 5 cells and
  all 5 are sheath glia, and the rule is now in `fates.bio`: **sheath 6/22 →
  11/22** with nothing lost elsewhere, **496 → 501 of 555** and 902 → 907 of
  961 to the adult. Socket glia get no rule and that is the measurement: the
  best socket rule that exists (SOX-2 ∧ PAG-3) takes 12 neurons to win 9 glia
  and costs 4 fates net. Exhaustively, over all 250 factors, 25,867 pairs and
  every triple, the in-sample ceiling is **+6 fates for sheath and +3 for
  socket**, the best sheath rules all contain PROS-1 and the best socket rules
  contain nothing anyone has named; rules fitted at in-sample precision 1.0
  score 0.21 (sheath) and 0.29 (socket) on a held-out sublineage against a 7%
  base rate. A glial cell differs from its sister by 39 factors, two sisters of
  the same tissue by 35, and the atlas's own two-reporter floor is 20.7%
  disagreement: the separation is below the measurement. What is missing is not
  another transcription factor but the ligand a glial cell receives and the
  terminal genes after the bean stage where the atlas stops.
  **PAR polarity** (`scripts/celegans_par.py`, `celegans_digital_development`):
  the first two divisions segregate the maternal factors only while the PAR
  domain that does the segregating is there — `div_P0` conditional on PAR-3,
  `div_P1` on PAR-2, the EMS fate rule reading `cell = AB|EMS|P2`, and MOM-2
  presented by P2 only while it has PIE-1. Against Digital Development
  **7 of 11 → 11 of 11 modelled transformations, 0 missed**, in both founder
  layers: par-3 gives AB adopting EMS, and in par-2 one lost identity (P2 keeps
  no PIE-1) produces all three of its observed changes. The wild type does not
  move — 1,439 cells, 501 of 555 fates, deaths 110/110 — because PAR-2 and
  PAR-3 are maternal factors of the zygote, and `mutants.bio` gains the two
  experiments so the claims are asserts. One extra prediction recorded: the
  same mechanism predicts E adopting MS in pie-1, which the table does not list.
- **Done 2026-09-15. The fate rules rewritten on the read the runtime computes**
  (`scripts/celegans_fate_reads.py`, `celegans_fate_reads`). The engine landed the
  integrated reads of §7.2a (8042a57), so `fates.bio` no longer tests
  `ELT-2_integrated = present` against a generated per-cell lookup whose threshold
  was 20% of a factor's largest path exposure **over the finished run** — a
  statistic of the answer. `exposure.bio` is deleted and the rules read
  `ELT-2.exposure(lineage) >= 15`, which the Body integrates from the reader's own
  presence calls. **15 min is one AB cell cycle** (Sulston 1983), so the rule says
  the path carried the factor across at least one division, and it is a plateau
  rather than a knob: 1 to 60 min give the identical 164 cells and 25 errors on the
  cited rules, 1 to 30 min the identical program, and above the plateau the score
  rises only by claiming fewer cells. Scored through the Body and the existing diff,
  with the baseline regenerated and run in the same eight folds:
  **501 → 522 of 555 embryonic terminal fates, 907 → 928 of 961 to the adult, and
  held out by founder sublineage 476 → 483 of 555 (85.8% → 87.0%)**; 1,439 cells,
  0 parent mismatches, deaths 110/110. Muscle 113 → 121 of 122, sheath glia
  11 → 20 of 22, socket glia 9 → 16 of 18, neurons 212 → 208 of 226, coelomocytes
  still 0 of 4. **The three negatives.** The held-out gain is a third of the
  in-sample gain (+7 against +21): the rewrite overfits harder, its in-to-held-out
  gap 39 fates against the lookup's 25. It claims 105 fewer cells, and scored
  without the lookup as fallback (factors, else the sublineage majority) it is
  **worse** — 377 against 414 — because the precomputed statistic normalised each
  factor by its own maximum over the run, which was carrying real information that
  no cell has. And the cited sheath-glia rule stopped being load-bearing: 522 with
  it and 522 without it, sheath 20/22 either way, kept and recorded as worth
  nothing. **The ordering survives a real runtime and the cell window wins
  nothing.** On the cited rules with nothing fitted, read at each cell's own
  decision point, the instantaneous read claims 225 cells and gets 61 wrong and
  `exposure(lineage)` claims 164 and gets 25 wrong; through the whole program 483
  against 522 in sample and 467 against 483 held out. But every one of the 555
  embryonic terminal cells gets **exactly one decision point, at its birth**, so
  `F.exposure(cell)` and `F.mean(cell)` are zero for every factor and every cell
  and decide nothing at all. "The mean over the cell's own life", the best read in
  the earlier measurement at 29 errors, is a summary of what the cell went on to
  carry *after* it chose — a second statistic of the future, larger than the first.
  **A fate written on a growing integral is only a fate if the cell cannot take it
  back.** The integral grows, so a rule whose threshold was not met when the fate
  was settled meets it later for no reason but the clock: with a 6-minute
  `cell_network` cadence the worm fell from 522 to 386 where the time-invariant
  lookup did not move. The obvious engine patch is not the fix and that was
  measured — dropping `d.applies(c.fate_ctx)` from the precedence refusal recovers
  only 31 of the 136, because most revisions come from rules of higher priority
  whose guard was genuinely false at birth. The fix is §7.2a's own construct, and
  `embryo_factors.bio` and `embryo_contacts.bio` now declare
  `commitment terminal_fate` over every terminal type — the first use of
  `commitment` outside the plasticity series. With it: **522 at cadence 0 and 522
  at cadence 6**, and at cadence 0 the program is identical to the cell with the
  block and without it, so its ablation is the whole of its justification. The
  cadence invariant is an invariant **for time-invariant guards**, and any program
  whose fate guards are integrated reads needs a `commitment` to have a stable fate
  at all; that belongs in §7.2a beside the two limits it already states.
  Two consequences of the rewrite elsewhere, both measured rather than assumed.
  **The precursor-commitment hole of §7.2a decision 3 is closed**: of the 12 cells
  born after the 700 min induction in the plasticity series' terminal arm, **0 now
  take the forced fate**, where the note recorded 12 of 610 taking it. What closed
  it is the fate rules, not the new lock — the count is 0 with `terminal_fate`
  stripped as well — because the rewrite makes those parents reach a terminal fate
  before they divide, so `plasticity.bio`'s own `commitment` with
  `inherit: daughters` already protects the daughters. **The human tissue
  generator's split fractions said what they are not.** Splits at one decision
  point run in sequence, so a splitting `fraction` is a share of what is left
  (docs/BIOLANG-v0.3.md); `human.py` converted each published adult count as though
  it were a share of the layer and every line carried the evidence "share of the
  layer's adult cell count". After the first split that sentence is false and the
  gap is enormous: ectoderm's glia read `fraction: 1.0000` for a 24.5% share and
  mesoderm's myocytes read `fraction: 1.0000` for 0.00071%. Every emitted line now
  carries both numbers — the cells and the published share of the layer, then the
  fraction with the split's position and what the earlier splits left — and the
  arithmetic is untouched: all 80 decisions in `tissues.bio` keep byte-identical
  fractions and `body.bio` still reaches 2.83e13 cells at 20 years. Checked and
  **not** changed: `haematopoiesis.py` has the same shape and is not the same
  defect, because its solve already encodes the sequential semantics
  (`retained = prod(1 - f_i)`), so each compartment's total outflow is right; what
  is off there is the split between two siblings, by 1.7% (MPP → CLP) to 7.9%
  (GMP → Monoblast) of the solved share. Compensating the emitted fractions by
  hand was tried and reverted: it moves a solved steady state the 39 mutant
  asserts are tuned to, taking erythrocytes at 3 years from 2.5e13 to 8.8e12.
  *(2026-09-29, lane-repro, `b587bb1` + `236fa2b`: the sentence above describes the runtime before
  `9d42485` (§7.5, 2026-09-15). Since then a cell decides again when a read it waits on crosses its
  threshold, so the `F.exposure(cell) >= 1` arm **fires on 170 cells**, 165 of them re-asserting the type
  the lookup already wrote 1.0 min after birth; **the factors still decide the fate of 5**, and no crossing
  changes a type (555 of 555 identical with rechecks off). `F.mean(cell)` read at birth is now the
  instantaneous read (487 of 555 in sample, not "nothing"). With the commitment removed the score now
  falls at cadence 0 as well (522 → 412). `celegans_fate_reads.json` is regenerated with both counts side
  by side; the conclusion about the commitment stands.)*
- **Next.** 1. `commitment` and `competence` landed on 2026-09-14
  (genomeos-c2, 126e65b) and earn their place by ablation: strip competence and
  both late arms of the published plasticity series fail, strip commitment and the
  terminal arm fails while the early arm collapses from 610 converts to 174. The
  exposure and mean reads landed on 2026-09-15 and the rules are rewritten on
  them (above), which also put the first `commitment` into the worm's own
  program. 2. PAR polarity rules:
  **done 2026-09-14** (above), and Digital Development is now 11 of 11.
  3. Contacts and AC/VU: **done 2026-09-14** (above). 4. Glia from factors:
  **done 2026-09-14** (above); what would move socket glia is a measurement the
  atlas does not contain, so this stays closed until there is one. 5. Three
  runtime changes area E asked for and did not make, written up with their
  measurements in ORGANISM-FROM-ONE-CELL.md: **all three fixed 2026-09-15** by
  genomeos-c2 (8042a57), along with `asymmetric` conjuring a factor the mother
  lacked. One request stands in their place, and it is a line of the specification
  rather than a patch: the cadence invariant holds only for time-invariant guards
  (above). 6. A measured time-resolved contact table for the embryo, to replace
  the curated one. 7. libRoadRunner and MaBoSS adapters tested in CI on 3.12.
  8. A terminal cell decides once, at birth, so nothing in the worm can read its
  own window; a `differentiate` that may be read after a stated delay, or a
  re-decision that is not a network cadence, is what would let the `cell` window
  be tested at all. 9. **`share:` is wanted** (area A offered
  to build it only with a generator that uses it): an absolute split normalised
  across the decisions at one decision point, so a partition line states the
  published number instead of a ratio of leftovers. Two generators would use it
  the day it lands — `human.to_bio_tissues`, whose lines would carry Sender &
  Milo's share directly and be checkable against the paper without replaying a
  sequence, and `haematopoiesis.to_bio`, whose solve would no longer have to
  encode the runtime's application order in `retained` and whose sibling splits
  would be exactly the solved ones. It also removes an order-dependence of the
  same kind §7.3 removed for fates. Until it lands both generators state both
  numbers on every line, which is honest but is a sentence standing in for a
  construct.
- **The two things `542d613` left open, settled and the shipped read left alone
  (2026-09-22, lane-worm2, `dddda7b` registration, `d5215bf` result).**
  **The 332 → 414 move in cells claimed is §7.5 `recheck: crossings`**
  (`9d42485`) and nothing in area E's program, which is byte-identical in the
  measuring arm: `recheck="none"` reproduces **332 exactly**, 82 cells join, all
  82 claimed-and-right and **not one wrong**, and the score is 522 either way.
  The mechanism is the lock rather than luck — the lookup writes the terminal
  type at the cell's birth, `commitment terminal_fate` is established at the end
  of that same decision point, and `_pick_fate` thereafter refuses any decision
  whose target differs. **250 of the 414 claimed cells fire only at a later
  crossing, every one already committed to the very type the rule then names.**
  So the added cells were decided before the factors were asked, and the score
  cannot see them because there is nothing there to see.
  **`mean(lineage)` wins every number this area has published and the shipped
  read does not change.** It takes honest credit 392/383 against 381/374,
  precision 0.947/0.905 against 0.920/0.839, and the fallback score 533/515
  against 522/483 — and it **loses free credit 56 to 136**, because 354 of its
  414 claims are made after the fate is locked, against 164 that fire at the
  cell's own birth under `exposure`. A fraction of a growing window crosses its
  threshold later than an absolute total crosses 15 minutes. The generous
  reading of "free" orders them the same way, 263 to 299.
  **The registration named that trap before the run, in its own words**: no
  factor rule has ever decided a fate the lineage had not already written, 414
  of 414, so a read that scores higher while that count stays at **0** has not
  made the factors decide one more thing — it has agreed with the answer sheet
  more often on cells already filled in. Both reads decide 0, and **no read can
  move that**, because every rule is guarded on a terminal `cell_type` that only
  the lookup sets. All six registered predictions held and nothing in the
  program was rewritten; both numbers are published side by side.
  **Four numbers now belong together** and the doc prints them that way: 522
  with the fallback, 381 honest, 136 free, **0 decided without the lookup having
  written the fate first**. The metric had removed the lookup as a fallback and
  not as a lock, and the second borrowing is the larger one.
  **Next in this area, and it is not a read choice**: the terminality gate — 44
  of 555 when the lookup's fates are deleted and the guard opened, with more
  non-terminal cells claimed than terminal ones. A rule set that can say *when*
  a cell is done dividing, not only what it becomes, is the next real
  measurement here.
- **The terminality gate, asked and answered: no (2026-09-22, lane-terminal,
  `de7a135` registration, `9d266a3` result, 0 requests).** On **1,326**
  embryonic cells — 555 terminal, 771 dividing — with the lookup's terminality
  deleted and **every cell called**, so that abstaining cannot help, the
  measured factors reach **0.8079** balanced accuracy held out by founder
  sublineage. That is **ahead of a single depth threshold (0.7815)** and
  **behind a per-founder depth threshold: 0.8998 in sample, 0.8773 held out** —
  and giving the factors depth as well makes them **worse**, 0.8148, an
  increment of **−0.085**. Both registered gate clauses fail.
  **So the 0 of 414 is no longer only an artefact of the `cell_type` guard: it
  now has a measurement behind it.** There is nothing behind that guard to lift
  it to. A rule set calling terminality at 0.65–0.81 where the lineage calls it
  at 0.88–0.90 would not be a better program, it would be a worse one that owed
  the lineage less. **522 / 381 / 136 / 0 stands as the ceiling**, and the last
  of those four is now measured rather than structural.
  **The circularity audit came first and excluded the signal that would have
  looked most like a measurement**: tracked lifetime ends at the cell's division
  *or at the last frame*, and the movie stops at 400.0 minutes where **454 cells
  end — 322 terminal against 78 dividing**, so a short lifetime is a division
  that was seen. Also excluded: the cell's own cycle length (defined only for a
  cell that divides), atlas coverage (455 of 555 terminal named against 717 of
  771 dividing, eleven points of signal about which cells the imaging resolved),
  `cell_type`, and the lineage name itself, which spells out its divisions.
  **The lane then corrected its own registered statistic, which is the part to
  keep.** Clause (c)'s pooled within-band accuracy is inflated by the stratum
  prior: a caller naming each stratum's own majority scores **0.7886** while
  scoring **exactly 0.500 inside every stratum**. The committed number stays,
  the macro average is published beside it, and a test pins both — including
  that the naive pooling exceeds 0.70, so the artefact cannot come back
  unnoticed. On the corrected statistic the factors do know something real:
  **0.6541** with depth and founder both fixed and sisters split across folds,
  against a floor of exactly 0.500 and a shuffle null of 0.5003 ± 0.0087. They
  are not silent, only quieter than the lineage's own count of how many times
  each founder divides. **Generation 7 is at chance** — a whole band where the
  atlas says nothing.
  **Four of six registered predictions broke, every one in the factors' favour**,
  which is the right direction for a lane reporting a negative: they beat the
  held-out depth threshold, they read factor identities and not only the clock
  (0.8079 against 0.6721 for a count with identities erased), and the shortfall
  is not the rule language.
  **Area E's next real number is not in Ma 2021.** It needs contact and geometry
  at the moment of division, or a measured cell-cycle regulator.
- **Owner.** genomeos-d3 (from 2026-09-15; genomeos-d2, genomeos-d1 and
  genomeos-73 before); **ownerless since that session ended, and lane-worm2 took
  the two open items on 2026-09-22.**

### F. Cancer and therapeutics

- **Goal.** From a tumour's alterations to what physically distinguishes
  those cells, what a therapy could reach and which mechanism the biology
  supports, with evidence per claim and a design dataset a binder-design
  engine can consume. It stops at describing the recognition required.
- **Code.** `cancer/*`, `therapeutics/*` (pipeline, mechanisms, scoring,
  evidence, design, providers, expression, neoantigen, structure,
  trafficking, localisation, logic, report).
- **Design.** CANCER.md, THERAPEUTICS.md.
- **Data.** `cancer_msk_impact_2017`, `expression_brca_tcga_pan_can_atlas_2018`,
  `data/demo/therapeutics/*`, caches under `data/knowledge/therapeutics`.
- **Requirements.** Three target spaces searched in parallel; target,
  binder, mechanism, cargo and effector kept apart; negative set as
  important as positive; patient data levels reported; nothing above its
  data.
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
- **Copy number and structural variants, done (2026-09-15).** The same
  cBioPortal client reads the study's `_cna` and `_structural_variants`
  profiles and distils them beside the mutation table
  (`cancer_alterations_msk_impact_2017`, 110 genes, 62 s, 139 KB). The
  mutation table alone calls CCND1 and MYC passengers at under 1%; both are
  amplified in over 4%. CDKN2A is deep-deleted in 7.6% of tumours and 32.6%
  of gliomas, ERBB2 amplified in 4.0% overall and 14.1% of breast
  carcinomas. A gene now reaches the tumour comparison and the target
  ranking through a copy-number or structural call exactly as through a
  variant: same origin record, same cohort grading, same score scale
  (`cancer/alterations.py`). Two refusals are the substance of it. A
  homozygously deleted gene is classed `unsuitable`, dropped from the
  surface targets, and fails a new hard requirement
  (`gene_product_present`) on every mechanism that must recognise a
  product — the tumour makes none of it. And the two copy-number
  conventions disagree about the number 2, so the format is decided by a
  stated rule and an ambiguous table is read the way that invents no
  amplification.
- **What the benchmark says about it: nothing, and that is the finding.**
  Before and after are identical in all six rows — same targets, same
  ranks, same classes, same mechanisms, same scores, 6/6 recovered, 6/6
  verdicts, 4/6 mechanisms defensible. The benchmark cannot see this work,
  because its one copy-number case carries a VCF passenger *inside* ERBB2
  (chr17:39700064), so ERBB2 was always recovered through the variant and
  never through the amplification, against that fixture's own stated
  intent. The control that does see it is in
  `tests/test_cancer_alterations.py` with fixtures in
  `data/demo/alterations/`: run at 8a4fe6e, a 12-copy ERBB2 with no ERBB2
  variant yields one candidate and it is not ERBB2; run after, ERBB2 is
  rank 1, direct surface, blocking antibody 0.521. A case whose driver is
  only a copy-number or structural call belongs in the benchmark;
  genomeos-f7 owns it.
- **The `cancer.*` library layer, done (2026-09-15).** Five libraries whose
  membership is computed from the two committed tables rather than written
  down (`cancer/libraries.py`, `genomeos cancer libraries`): `cancer.mutated`
  87 genes, `cancer.hotspots` 8, `cancer.amplified` 12, `cancer.deleted` 3,
  `cancer.rearranged` 25, each member carrying the frequency that put it
  there. `--gene CDKN2A` answers the question the layer exists for: deleted in
  7.6% of tumours and 32.6% of gliomas, mutated in 4.3% with truncating
  hotspots R80* and R58*, rearranged in 0.23% with MTAP among the partners.
  Three refusals. A hotspot gets no per-cancer-type breakdown, because the
  distillation counts recurrent changes study-wide and printing the gene's
  *mutation* distribution beside a hotspot frequency would read as the
  hotspot's. A gene off the panel is reported as never looked at, not as
  unbroken: `breaks_in("CD19")` says so in the sentence. And membership means
  selected, not driving, and not treatable. The layer does not join the shared
  `LIBRARIES` catalogue on import, since `genomeos/lib` is area C's file;
  `register(LIBRARIES, LAYERS)` wires it in one call whenever that area wants
  it, and a test asserts that importing `genomeos.cancer` mutates nothing.
- **Healthy tissue shipped rather than fetched (2026-09-15).** The Human
  Protein Atlas answers a whole protein class in one request, so its 5,573
  predicted membrane proteins and 384 CD markers are fetched once, kept to
  the twenty tissues GenomeOS weighs, and shipped gzipped in the package
  (5,588 genes, 396 KB, two requests, 7 s; `scripts/normal_tissue.py`, record
  in `normal_tissue_atlas`). Normal-tissue safety is now answered offline for
  every surface gene instead of being unknown and capping the score, and a
  scan over the membrane universe becomes possible, which is the missing
  input for expression-driven targets. **Two of the twenty tissues had never
  been read**: the Atlas returns `t_RNA_skin_1` and `t_RNA_stomach_1` under
  the titles "skin 1" and "stomach 1", the lookup asked for "skin" and
  "stomach", and both were silently recorded as absent on every gene ever
  scored. Skin is where an EGFR antibody's classic toxicity shows, and EGFR
  carries 48.4 nTPM there. Column titles are derived from the field ids now
  and a test asserts all twenty are read. The benchmark does not move: both
  surface cases already sit on the poor-safety cap (EGFR raw 0.572, held at
  0.5), which is the reason a benchmark is not the only control a change
  needs.
- **Expression-driven targets, done (2026-09-15).** CD19 and BCMA are not
  mutated, amplified or rearranged, so a pipeline that reaches a gene only
  through an alteration could not propose either: a hole in the middle of the
  target space. `--scan N` closes it from the patient's own RNA against the
  packaged atlas (`therapeutics/scan.py`). On the demo B-cell lymphoma it
  proposes FCRL5 at 5.9x and CD19 at 4.3x, both direct surface targets ranking
  above the TP53 driver and neither reachable by any other route here. The
  rejections are the more instructive half and are reported rather than
  dropped: CD20 at 2.2x, CD79A at 1.4x, CD22 at 1.0x — the targets of the most
  successful antibodies in oncology, each failing a selectivity screen because
  healthy spleen is full of the same lineage (CD20 sits at 247 nTPM there). A
  threshold here is a sort order and never a verdict. Three refusals hold: it
  runs only on patient RNA and never on a cohort, raised transcript is never
  called protein on the surface, and every candidate carries what the CD19
  story costs — B-cell aplasia is not a side effect of CD19 therapy but the
  same event from the other side. The benchmark does not move: the scan is off
  unless asked for, and none of its six cases supply RNA.
- **Evidence tier, done (2026-09-27, `5013406` registered, `bdd5058` built).** Each candidate carries
  how its alteration is known — observed in this tumour (1.0), measured in this patient (0.6), or an
  association hypothesis measured nowhere (0.2) — and ranking is tier-major, outside the weighted mean
  (`DIAGNOSTIC`). The registered shape put the tier *inside* the mean at weight 1.1 and fired its
  falsifier: PIK3CA rose over ERBB2 in the HER2-amplified case, because the mean runs over available
  dimensions and a poor-safety cap bites unevenly. Amended to a gate: all nine scores unchanged,
  ERBB2-amp 4 → 1, KRAS 5 → 1, PIK3CA 5 → 1, IDH1 and ALK 2 → 1, CD19 held at 2; the twelve rows where
  a guess outranked the target are zero. Left open by it: no rule asks whether a candidate has any
  mechanism at all, and inside the top tier twelve copies and one missense read the same.
- **Mechanism gate and magnitude tiebreak, done (2026-09-28, lane-mech, `028085c` registered,
  `2805552` built).** The two items the evidence tier left open are closed on the order, not in the
  mean. `mechanism_reach` partitions candidates by whether any modelled modality reaches them with
  its hard requirements answered — read after the tier, never before it, provisional-only in the
  same class as nothing — and `alteration_magnitude` publishes what the patient's data measure
  (copy count against the diploid 2, allele fraction, hotspot status), stating absences and breaking
  ties only between candidates equal on tier, gate and score. No falsifier fired: all nine rows are
  identical field by field; new pin `outranked_by_unreachable` 0 for all nine. **Five of the nine
  targets are `no_established_mechanism`** — the small-molecule cases — which the ranking now says
  rather than ordering them as though a modality existed. **Negative: the gate closed no case**, so
  its value is untested against a tumour where an unreachable candidate outranks a target, and the
  tiebreak has never fired end to end because the demo VCFs carry no allele fraction.
- **The gate's case (2026-09-28, lane-gatecase, `eb18c40` registered, `64f0afe` run, `204ce6c` tests).** A
  tenth benchmark case, chosen and registered before it ran: HER2-positive gastroesophageal
  adenocarcinoma with a co-amplified MYC, from TCGA stomach adenocarcinoma (ERBB2 amplified in 58 of 440,
  MYC in 53, both in 20), at that subgroup's median copy numbers, with TP53 R175H at the cohort's median
  allele fraction 0.49 — the benchmark's first allele fraction; cohort summaries only. Approved targets:
  trastuzumab (FDA 2010) and trastuzumab deruxtecan (FDA 2021). **Negative: the gate closed nothing
  again** — the rank without the gate equals the rank in all ten cases; MYC scores 0.246, last of six,
  because the annotations that make a gene unreachable are the ones that score it low. The configuration
  the gate prevents needs an unreachable surface receptor beside an amplified approved target, two
  independent drivers in one tumour: in MSK-IMPACT 2017, of 42 tumours whose structural variants make a
  surface receptor the 3′ partner, one also has an amplified ERBB2. **Second negative: the tiebreak was
  never asked** — every score tie is between hypotheses, which carry no measured quantity. A registered
  prediction was falsified by the registration (a row publishes the target's magnitude, so the 0.49 sits
  on TP53); amended additively. ERBB2 rank 1; the nine earlier rows unchanged. **Open:** the regenerated
  result is not committed — it is a headline result, so landing the tenth case is one commit that
  rebuilds it cleanly and moves the (9, 9, 9) pins and README's "9/9" to ten; its four tests skip until
  then.
- **Missing.** The demo runs at data level 1 only; what a deep deletion is
  actually worth (the dependency it creates) is not modelled; the scan's
  healthy reference is 20 tissues of population consensus, not a matched
  normal.
- **Next.** 1. A `NeoantigenProvider`
  behind a licence-checked optional extra. 2. ~~A benchmark case whose driver is
  a copy-number, structural or expression call, and lowering the pinned defect
  count from 2 to 0.~~ **Done 2026-09-21 (`e21d34a`), and both halves were
  smaller and larger than the row said.** The pin was payable: both defects
  closed on 2026-09-15 with `6b45bfe`, when a preferred mechanism was required
  to be *established* rather than merely unrefused, and the constant read 2
  against code that gave 0 for six days because the file belonged to another
  lane. Lowering it alone would have been vacuous — the zero is reached by
  `best_mechanism` being `None`, and a `None` cannot be indefensible, so a
  silent pipeline would have scored as sane. Each row therefore now records
  whether its preferred mechanism was established, the nearest provisional one
  and the requirement left open, and a test asserts the property the two
  defects violated instead of the number they produced. BRAF still reads
  `blocking_antibody` with `surface_accessible` unanswered and PIK3CA `adcp`
  with the same: demoted with their reason, not deleted.
  The seventh case is **CD19 in a B-cell lymphoma, reached by an expression
  call and by no DNA event at all** — the tumour's only DNA driver is TP53 — so
  it is the first scored case that is not a point mutation, and it tests the
  route rather than only the recovery, because every approved CD19 drug is
  antibody-like. **7/7 recovered, 7/7 verdicts, 7/7 defensible**, CD19 at rank 2
  with `adcc` 0.655 established, which is tafasitamab's mechanism.
  **Probing the other routes first turned up a third defect of the same class,
  and it was the day's finding.** Given EML4-ALK the pipeline offered a blocking
  antibody at 0.75 against what is, in that tumour, a cytoplasmic kinase: ALK is
  the 3' partner and its ectodomain is not in the product, so there is no
  approved antibody against it and there could not be. Same mistake as BRAF's
  one layer down — a curated compartment describing the full-length protein,
  asserted for a product that is not it. A fusion-only origin now leaves the
  surface requirement unanswered. **Half of it is pinned rather than argued
  away**: the candidate keeps `direct_surface` and accessibility 1.0, because
  class and score come from the gene and not the product, and closing that needs
  a measurement this project does not hold — the junction, the 5'/3' orientation
  or transcript evidence for the retained domains. A test asserts today's wrong
  answer on purpose so it fails when the measurement arrives. Copy number and
  structural variants still have no scored case, which `cases_by_driver_call`
  now states in the result instead of leaving it inferable from the case list.
  **The copy-number half of that gap closed the same night (`3aee02f`), and it
  is a pass with a negative inside it.** The eighth case is **ERBB2 a second
  time, amplified and not mutated** — the repetition is the design, since the
  point-mutation case reaches the same gene through a coding change, so the
  route is the only thing that differs, where CD19 varies target and route
  together and cannot separate them; it is also the call the clinic makes,
  amplification being the companion diagnostic trastuzumab is prescribed on.
  **8/8 recovered, 8/8 verdicts, 8/8 defensible**, and `cases_by_driver_call`
  reads six point mutations, one copy number, one expression.
  **But the three questions could not see what was wrong with the answer.**
  ERBB2 is recovered as `direct_surface` with an established blocking antibody
  and **ranks fourth, behind KDR, EGFR and PDGFRB — three surface proteins with
  no alteration in this tumour at all**, named only for being STRING partners
  of the mutated PIK3CA. **Recovering a target and preferring it are different
  results, and every verdict read as a pass.** A fourth measurement per row,
  `outranked_by_hypotheses`, now records it, and it shows the same pattern in
  two of the four approved-antibody cases and three of the four intracellular
  ones — so it was not invented for this case.
  **Under it, a duplicate and the reason it was invisible.** ERBB2 was proposed
  twice, once carrying its twelve copies and once as a hypothesis about PIK3CA,
  because `pathway_induced` seeded its exclusion set from the disrupted drivers
  alone; the existing control could never have caught it, since it passes
  `indirect=False`. **Both entries scored 0.494, identically**: for a surface
  target the score reads the gene's curated annotation and not the alteration,
  so twelve copies and a guess about a neighbour are worth the same. Fixing it
  freed a slot, KDR entered at 0.575 above everything, and ERBB2 went from
  third to fourth — **the fix made the ranking look worse and is still right.**
  The duplicate is fixed and pinned (`BURIED_SURFACE_TARGETS` at 2, may fall
  and never rise); **the scoring is recorded and not fixed**, because teaching
  the score to read the alteration changes every case's numbers and is a
  decision rather than a side effect of adding a case. It is the next thing
  this area should be asked for.
  3. What a deep deletion is worth: the dependency the loss creates, which
  nothing models.
- **Owner.** genomeos-f2 since 2026-09-14; the benchmark and the 2026-09-21
  defect work ran as lane-cancerF under the coordinating session.

### G. Product: web UI, CLI, packaging, CI

- **Goal.** A geneticist installs it, opens one page and gets every layer
  with its evidence; a developer trusts CI.
- **Code.** `web/server.py`, `web/static/*`, `cli.py`, `jobs.py`,
  `storage.py`, `results.py`, `.github/workflows/ci.yml`, `pyproject.toml`.
- **Design.** README.md, DATA.md, STORAGE.md, ACTION-PLAN.md UI track.
- **Data.** `data/jobs/*`, `docs/PROGRESS.md` (rendered in the Progress
  tab).
- **Requirements.** One HTML page, inline SVG, no build step; every control
  explained; jobs visible with progress; the Progress tab reads this
  project's own log.
- **Project progress as current work and done (2026-09-12).** The Progress
  tab's Project progress card reads the roadmap and a work board as data:
  `genomeos/roadmap.py` parses ROADMAP.md §3, §4 and §6 into areas, planned
  steps and finished items; `genomeos/work.py` is the board, one git-ignored
  `data/work/<session>.json` per session written by `genomeos work
  start|update|done --who SESSION`; served by `/api/work` and
  `/api/roadmap`. The panel shows who is working on what, jobs running,
  changes not yet committed, what landed in the last day and what is planned
  per area (genomeos-77). It reads this file's formats as they are: a bullet
  whose bold lead is not one of Goal, Code, Design, Data, Requirements,
  Missing, Next or Owner is a done item; a "Done YYYY-MM-DD:" sentence inside
  Missing is done; Next steps are numbered "1. … 2. …", and a step that
  starts with "Done" is done, or partial if it also says "next".
- **The Evidence explorer (2026-09-13).** The view the UI track promised in
  2026-09-10 and never built. `genomeos/evidence.py` reads every BioLang
  program in the project (the demos, the organism programs, the `bio.std`
  prelude) into one row per stated fact: the block it comes from, its
  evidence kind, source, note and confidence. A program is credited only with
  the facts its imports do not already declare, so `bio.std` is counted once
  rather than once per importer. The Evidence tab filters by evidence kind,
  by a confidence ceiling, by program and by text, shows the mix as a bar and
  one line per program, and exports the selection as a CSV review list;
  `genomeos evidence [--kind K] [--max-confidence X] [--by-program] [--csv F]`
  answers the same from the command line, and `/api/evidence` serves it.
  The first reading of the project's own model: 22,419 facts over 27
  programs, mean confidence 0.57, 48% experimental, 46% predicted, 4%
  curated, 2% inferred and 2 facts with no evidence at all; 9,974 facts sit
  at or below 0.5, and 9,782 of those are the compiled chr21 non-coding
  program, whose mean confidence is 0.27. **Re-run 2026-09-27: 26,845 facts over 41 programs, mean confidence 0.617, 10,075 at or below 0.5 — experimental 14,655 (54.6%), predicted 10,352 (38.6%), curated 1,343 (5.0%), inferred 492 (1.8%), none 3.** The first reading above is kept beside it; mean confidence rose because what arrived is mostly measured, and the weak half did not shrink (10,075 against 9,974). Measured and predicted are almost
  equal in number, which is what a project that compiles public data and then
  predicts on it should look like; the weak half is concentrated in exactly
  one place, area I's attributions, rather than spread through the
  hand-written biology.
- **Missing.** Every view promised in the UI track is now built. The last, **Cell**, shipped on
  2026-09-21 as `7abb450` with `/api/cell`, `/api/cell/programs`, `/api/cell/run` and
  `tests/test_cell_view.py`; this line said it was unbuilt until 2026-09-27. Its "no nightly
  real-data CI" clause was stale too — Next item 4 records it done on 2026-09-17. **What is genuinely
  outstanding: no git tag and no PyPI package** (the version reads 1.0.0 since 2026-09-27), and no
  "known phenotypes it must reproduce" list per library from a biologist, which is a person this
  project does not have rather than a task.
- **Next.** 1. Done 2026-09-13: the Evidence explorer over every program,
  with the CSV review list; next, the same reading over the compiled
  chromosome programs of area I, which is where the weak half sits.
  2. Done 2026-09-11: README refresh and the version bump to 0.9.0; the git
  tag follows the merge to `main`, which only Albert opens. 3. Done
  2026-09-11: `.gitignore` for `data/jobs`, keeping the registry in code.
  4. Nightly CI job that runs the real-data tests against cached reference
  data: **done 2026-09-17, and it was mostly already there.** The daily run has cached
  chr21, chrM, GO, GOA and CL since the cadence change and runs the whole suite against
  them, so "real data" was covered; what ran nowhere were the readers that reach the
  outside world, gated behind `GENOMEOS_LIVE` and `GENOMEOS_NET`. A second job in the
  same daily workflow sets both and runs the three that read Zoonomia and gnomAD Gnocchi
  by HTTP range and UniProt/AlphaFold by API — 26 assertions, 3.6 s locally. It is
  `continue-on-error`: an outage at UCSC is not a defect here, and a scheduled job that
  reddens on someone else's downtime trains everyone to ignore it. The AlphaGenome live
  test stays out because it needs a paid key as a repository secret, which is Albert's
  decision. Still one scheduled run a day. 5. Cell view: **done 2026-09-17**. The reader has run on **thirteen** cell types over (this said "eleven" until 2026-09-27; `reader_genome_wide.json` has named thirteen since 2026-09-14, the last two being testis and ovary)
  all 24 chromosomes since 2026-09-14 and its numbers lived pooled inside the Progress
  tab's genome-wide card, where cells could not be compared. `/api/cells` and the Cells
  view read `reader_genome_wide` per cell: genes read, promoter open, poised, read by
  marks, enhancers active and nodes silent, with the per-chromosome spread behind each
  row. The same genome reads 70.3% of its coding genes in hepatocyte and 53.2% in
  GM12878, which is the claim the reader exists to make checkable.
- **Owner.** genomeos-9c (views, packaging, CI, the Project progress card, the work board and the Evidence explorer; genomeos-77, genomeos-c6 and genomeos-f7 before the restarts); genomeos-f3 (jobs).

### H. BioForge (design under constraints)

- **Goal.** In-silico experiments and design search whose outputs are
  labelled predicted with a confidence, validated by reproducing published
  perturbation results.
- **Code.** `forge/design.py`, `organism/forge.py`, `forge/calibration.py`,
  `runtime/debugger.py`. *(This row said `forge.py` until 2026-09-22; there is
  no `genomeos/forge.py`.)*
- **Design.** DESIGN-MINIMAL-CELL.md, ACTION-PLAN.md Phase 5.
- **Data.** `design_neuron`.
- **The predicted confidence cannot be calibrated, and the obstruction is not n
  (2026-09-22, lane-forge, `e13df9e` registration, `eea443e` result, 0
  requests).** It is **a hand-set literal `0.3`** written onto every experiment
  a `design` block emits (`organism/forge.py:143`; *`:170` since 2026-09-27, `66675d8`*), plus six `min(x, 0.3)` caps
  — one of them in **area I, for an unrelated quantity**, which is the tell: a
  number that means the same thing for a worm-embryo design search and a phyloP
  budget is not a probability about either. It is a provenance marker reading
  *predicted, not validated*. It ignores loss, feasibility, evaluations and the
  count of rival candidates — across the four shipped designs the feasible
  counts run 2, 6, 2, 22 and all four answers are stamped 0.3 — and it
  overwrites the design's own stated 0.5 and 0.6. **A single-valued predictor
  has one reliability bin at every sample size, so no curve exists at any n**,
  which is read off the source rather than the data.
  **Census, counted before anything was designed: 2 designs whose published
  outcome is cited in the file, 4 counting two recorded only in prose, 0
  fixtures and 0 tests pairing a design with an outcome** — and all four solve
  at loss 0.0, so the outcome column is degenerate too. The 31 `experiment`
  blocks cannot substitute: an experiment *states* a perturbation and checks the
  simulator, while a design must *find* it.
  **The rate was computed and is withheld, which is the part worth keeping.**
  At 4 of 4 the exact 95% interval is [0.3976, 1.0000] and **excludes 0.3**,
  so the arithmetic invites raising the number. It is not published because the
  answer is planted in the population: `founders.bio:47` carries the same Lin
  1995 clause as the design's own answer key, at confidence 0.9, and every
  candidate list contains the published factor because the author put it there.
  BioForge is not predicting POP-1; it is evaluating a rule authored from POP-1
  over a list containing POP-1, so 4 of 4 measures that the simulator reads its
  own rules correctly — which existing tests already establish.
  **Five obstructions are named with n deliberately last**: a negative case, a
  blind candidate list, a design whose answer is not already in the module, a
  varying confidence, and then n. n is the smallest of the five and the only one
  more of the same work would fix, so "collect more designs" would have been the
  wrong response. Verdict **`UNCHECKABLE_BY_CONSTRUCTION`**, in
  docs/BIOFORGE-CONFIDENCE.md and `genomeos/forge/calibration.py`, with tests
  pinning the two claims a later edit could falsify. **Two description-only
  defects are named and not fixed**: the CLI prints the constant directly
  beneath a truthful loss-based verdict it bears no relation to, and the parser
  stores a stated confidence the emitter discards — one of those two behaviours
  is wrong.
- **Missing.** ~~No calibration of the predicted confidence against outcomes~~
  *(resolved 2026-09-22, above: it is not calibratable and not for want of n)*;
  knockouts of network species (the `experiment` block perturbs organisms;
  network perturbations are covered by tests: the repressilator's `# test:`
  claims and the Fauré cell-cycle knockout in `test_boolean.py`). Done
  2026-09-11: the `design` block (BioForge takes experiments as input;
  knockouts, additions and knobs searched toward targets under constraints,
  the answer emitted as a predicted experiment; POP-1, PIE-1, GFI1 and IKZF1
  found), eight worm mutants and nine blood mutants as experiment programs
  in CI, Digital Development as the published set. *(The worm mutants are
  **ten** since `par2` and `par3` were added; the nine blood mutants are right.
  And `haematopoiesis_designs.bio` is referenced only from prose — `bio test`
  runs `experiment` blocks, not `design` blocks, so no CI path executes it.)*
- **Next.** 1. Budget constraints in `design` (genes, bases, perturbations
  already there). 2. Network knockouts in `experiment` (a species held at
  zero). 3. Calibrate predicted confidence: how often a design's answer
  matches the published outcome.
- **Owner.** genomeos-73 after the experiment block.

---

### I. The 98% (attribution of the non-coding genome)

- **Goal.** Every UNKNOWN block with a best guess of what it is for the
  organism and the evidence behind it: from "98% no clue" to "98% a very good
  guess", tested by quantities (the house budget), by measured ground truth,
  and by the organism run forward. Albert's framing of 2026-09-11, the DNA as
  the blueprint and the living organism as the built house, is the opening of
  ATTRIBUTION.md.
- **Code.** `attribution/bigwig.py` (a bigWig reader over HTTP range requests,
  standard library only), `attribution/constraint.py` (Zoonomia phyloP over
  241 mammals, the 100-vertebrate conserved elements), `attribution/budget.py`
  (the tiers and the budget), `scripts/budget_genome_wide.py`, `genomeos budget`.
- **Design.** ATTRIBUTION.md.
- **Data.** `budget_chr*` (24 of 24), `budget_genome_wide`;
  `data/knowledge/constraint`, the conserved elements per chromosome (cache).
- **Requirements.** Constraint is read per base and never stored (the Zoonomia
  track is 9.6 GB; chromosome 21 costs 33 MB of ranges); a tier is a best
  guess with a confidence, never a verdict; every threshold is a named constant
  repeated in the result; "neutral" is a claim with evidence (unconstrained, no
  element), not the absence of one.
- **The first budget (2026-09-11).** Chromosome 21's 446 UNKNOWN blocks, 20.8
  Mb, with Zoonomia phyloP over every sequence-bearing base and the
  100-vertebrate elements over each block:

  | tier | blocks | Mb | of UNKNOWN | of chromosome | constrained kb |
  |---|---|---|---|---|---|
  | structural | 27 | 9.75 | 46.9% | 20.9% | 0.8 |
  | fossil | 104 | 3.91 | 18.8% | 8.4% | 39.9 |
  | regulatory | 194 | 3.93 | 18.9% | 8.4% | 50.6 |
  | constrained_unknown | 20 | 0.29 | 1.4% | 0.6% | 18.0 |
  | neutral | 101 | 2.93 | 14.1% | 6.3% | 39.8 |

  Constrained bases are 1.3% of the UNKNOWN space against 10.7% of the
  genome: constraint lives in and around genes, where the UNKNOWN pass does
  not look. Twenty blocks, 287 kb, are the real unknown of this chromosome;
  the most constrained, 47 kb at 16.82 Mb, is unique sequence with 179
  conserved elements, which is what an unannotated regulatory region or gene
  looks like. The registry's regulatory class is only 1.2% constrained: a cCRE
  is a chromatin state, mostly lineage-specific, not a conserved element.
- **The genome (2026-09-12).** The `budget_genome_wide` job ran the 24
  chromosomes in eight hours, smallest first, 2,252 MB of range requests over
  a 9.6 GB track that was never downloaded. 1,009 Mb of UNKNOWN blocks in
  3,088 Mb of genome, 780 Mb of them holding sequence, of which 15.0 Mb are
  constrained (1.93%, against 10.7% for the genome as a whole):

  | tier | blocks | Mb | of UNKNOWN | of genome | constrained Mb |
  |---|---|---|---|---|---|
  | structural | 238 | 229.7 | 22.8% | 7.4% | 0.05 |
  | fossil | 7,302 | 328.5 | 32.6% | 10.6% | 4.13 |
  | regulatory | 15,536 | 345.2 | 34.2% | 11.2% | 8.10 |
  | constrained_unknown | 1,098 | 32.0 | 3.2% | 1.0% | 1.66 |
  | neutral | 2,632 | 73.3 | 7.3% | 2.4% | 1.07 |

  Every block carries a guess at confidence 0.5 or above. The real unknown of
  the UNKNOWN space is 1,098 blocks and 32 Mb, one percent of the genome,
  holding 1.66 Mb of constrained bases; the regulatory tier holds five times
  more constrained sequence (8.1 Mb) with its target genes still to be named.
  The chromosomes sort themselves: chrY and chrX are the controls (0.66% and
  0.88% constrained, fossils 27% and 73%), the gene-dense chr17, chr19 and
  chr20 put half their space in the regulatory tier, the gene-poor chr4,
  chr13 and chr18 carry the most neutral and constrained-unknown sequence.
  Constraint per class ranges from 0.6% (centromere) to 2.5% (unique
  intergenic); no class of the UNKNOWN space approaches the genome's average,
  which is the quantitative form of "constraint lives in and around genes".
- **The constrained elements first (2026-09-12, genomeos-fe).** The first
  attribution of targets follows the budget's pointer: `scripts/
  constrained_targets.py` ranks a chromosome's distal enhancers by Zoonomia
  constraint through `attribution/constraint.py` (chr21: 6,618 elements, 37
  MB of ranges, 130 s; 412 have 20% or more constrained bases) and deletes
  the 100 most constrained in AlphaGenome (`constrained_targets_chr21`). 73%
  name a gene against 63.5% of the uniform sample of 200; 28 are
  silencer-like; the coding target is the nearest TSS in 60.9% (uniform
  67.8%) and inside the CTCF node in 78.1% (uniform 87.4%), so the constrained
  elements reach beyond the CTCF-only boundary twice as often. Strongest: a
  46%-constrained element moving LINC00945 by 1.2 log2 in testis, KCNE1 by
  0.93, SIM2 by 0.84. ALPHAGENOME.md feature b, "The constrained ones first".
  Genome-wide the same day (`constrained_targets_genome_wide`, a resumable
  job): 2,305 elements, the 100 most constrained per chromosome (59,415 of
  513,972 distal enhancers have 20% or more constrained bases).

  | distal enhancers deleted | name a gene | strong (≥ 0.3 log2) | silencer-like | target = nearest TSS | inside the node |
  |---|---|---|---|---|---|
  | 4,800 uniform | 62.4% | 24.5% | 38.6% | 70.7% | 86.0% |
  | 2,305 most constrained | **75.7%** | **31.9%** | 36.8% | 71.9% | 89.8% |

  Constraint predicts function on 22 of 24 chromosomes (chr1 65% to chr19
  92%; chr9 and chr15 the exceptions) and the strongest effects sit there:
  C1QTNF7-AS1 by 4.0 log2 (chr4, 76% constrained), MEF2C-AS2 by 3.2, ARX by
  3.0 in Purkinje cells (chrX, 92% constrained). Chromosome 21's "reach
  beyond the node twice as often" did not generalise: genome-wide the
  constrained elements sit inside their node as often as the rest, which is
  the node model holding on the elements that matter most. Area I's step 2
  has begun: a constrained regulatory block now comes with a predicted target,
  tissue and magnitude wherever AlphaGenome has been asked.
- **Against measured enhancers (2026-09-12, genomeos-fe; step 3).** VISTA's
  2,442 human elements tested in transgenic mouse embryos at e11.5, positive
  with tissues or negative, read blind through the registry, the CTCF node,
  Zoonomia constraint and the AlphaGenome deletion of the whole element
  (`attribution/vista.py`, job `vista_genome_wide`, `vista_chr*`,
  `vista_genome_wide`; the loci cached under `data/knowledge/vista`). Over 23
  chromosomes in five hours, 2,223 non-overlapping elements, 1,133 positive
  and 1,090 negative:

  | layer | positives | negatives | chromosomes leaning this way |
  |---|---|---|---|
  | an ENCODE enhancer-like element sits there | 84.5% | 66.1% | 22 of 23 |
  | constrained (≥ 20% phyloP bases) | 65.0% | 57.7% | 18 of 23 |
  | the deletion moves some gene | 72.7% | 64.9% | 17 of 23 |
  | strong effect (≥ 0.3 log2) among those named | 52% (422 of 806) | 39% (271 of 695) | |
  | predicted tissue in the measured tissue group | 68% (195 of 286 judged) | 35.5% by shuffling the labels | |

  Every layer leans the right way over the genome. The registry separates
  best, constraint least because VISTA chose its elements for conservation in
  the first place, and the two-chromosome reversal on the deletion (75% of
  positives against 91% of negatives on chr21 and chr22) was fifty elements'
  noise: over 2,200 the deletion names a gene more often for positives and
  the difference is in the strong effects. The tissue remains the part of the
  prediction the assay confirms most. A VISTA negative is conserved sequence
  next to a gene that did not drive expression on one embryonic day, while
  the model reads adult and cell-line tracks; the deletion effect stays a
  target finder first and an activity call second.
- **Against measured targets (2026-09-12, genomeos-fe).** GTEx v8 eQTLs as
  the ground truth for the target gene itself: the 1.56 GB archive streamed
  member by member over HTTP ranges in 135 s into one TSV per tissue under
  `data/knowledge/gtex` (34 MB), nothing else stored (`attribution/eqtl.py`,
  job `eqtl_targets`, `eqtl_targets`). A match is an eQTL variant inside the
  element or within 500 bp whose eGene is the named target; the control is
  the nearest coding TSS in the node.

  | elements | carry an eQTL | deletion target is an eGene (coding) | nearest TSS is an eGene | where they disagree |
  |---|---|---|---|---|
  | 4,800 uniform | 3,725 (78%) | **71.7%** (n = 1,851) | 66.2% (n = 3,725) | model 364, node 309 |
  | 2,305 most constrained | 1,394 (60%) | 56.0% | 54.9% | 148 against 154, level |

  Held to the same universe as the node (coding genes), the model beats the
  heuristic on the measured targets as it did on the enhancer deletions. Asked
  for any gene the model scores 52.9%, the drop being lncRNAs and pseudogenes
  GTEx has little power on, not evidence they are wrong. Constrained elements
  carry fewer common variants, so eQTLs are rarer there, which is why the
  deletion model was pointed at that tier in the first place. Where the
  predicted track is a GTEx tissue (101 elements) the eQTL for the same gene
  is significant in that tissue for 40%.
- **Against consequence (2026-09-12, genomeos-fe).** GWAS Catalog lead
  variants within 1 kb of an element, read against the same elements shifted
  100 kb as chance (`attribution/gwas.py`, job `consequence_targets`,
  `consequence_targets`): uniform elements hold one in 33.2% against 29.4%
  shifted (1.13 times), the most constrained 37.4% against 30.8% (1.21), VISTA
  47.5% against 42.9% (1.11). Within VISTA a named target raises it (50.6%
  against 41.0% with none), a strong effect more (55.8%), and constrained
  elements hold fewer (42.4% against 55.6%), the eQTL lesson again: selection
  removes the common variants a GWAS needs. The traits are the polygenic ones
  (height, BMI, blood pressure, type 2 diabetes, insomnia). The catalog's
  mapped gene equals the node's nearest TSS for 63% of hits and the model's
  target for 39%, which measures the catalog's own nearest-gene heuristic, not
  truth. ClinVar pathogenic variants with a non-coding consequence sit in 20,
  30 and 24 elements of the three sets; ClinVar's gene is the model's target
  for 11, 20 and 19 of them and the node's for 13, 17 and 12. Reading:
  consequence is the level where the regulatory attribution is least testable
  today, barely above chance and moved a few points the right way by a named
  target and a strong effect; the measured layers (VISTA, eQTL, lentiMPRA, the
  reader) are where it is tested.
- **The attributions as a program (2026-09-12).** `genomeos budget --chrom
  chr21 --bio` compiles what the lane has computed into BioLang
  (`attribution/compile.py`, `data/organisms/human/noncoding_chr21.bio`, 244
  kB): every UNKNOWN block as a `region` whose role is its tier and label with
  constraint as evidence, the constrained_unknown tier keeping `role: unknown`
  so the program's own count of unknowns is the chromosome's real unknown;
  every element with a predicted coding target as an `element` with its
  domain, target and basis plus one `rule` (activates or inhibits) carrying
  the predicted magnitude as strength, the targets declared once; the 63 CTCF
  nodes those elements sit in, each with the reader's view as a comment (open
  in which cell types, silent in which). chr21 from the two samples: 747
  entities (446 regions, 155 elements, 83 genes, 63 domains), 155 rules, 20
  unknowns, three `# test:` lines that `bio test` passes in 0.13 s (with every
  element scored, 5,972 entities and 5,174 rules, see the closure below); predicted evidence capped at 0.7,
  mean confidence 0.62 for regions and 0.25 for elements, which is what the
  evidence supports. The non-coding space of a chromosome is now something
  the engine reads, not a table. Since 848abc6 (genomeos-bb's hunk) a region's
  evidence note also carries the human axis where the chromosome has been read
  on it: "people N% of kilobases constrained, case syntax, relaxed, recent or
  tolerant"; the program's counts and checks do not move, since the axis is a
  note and not a tier.
- **Closure at the gene level (2026-09-12).** `genomeos closure --chrom chr21`
  (`attribution/closure.py`, `closure_chr21`) reads, for each of the 221
  coding genes and each of four ENCODE cell lines with both a DNase reader and
  a total RNA-seq track (K562, HepG2, GM12878, IMR-90), three things
  independently: the promoter's openness (the reader), the regulatory input
  from the attributed elements (an element counts as active where a DNase
  peak overlaps it, its predicted magnitude signed by action), and the
  measured expression (fraction of canonical-exon bases with RNA-seq signal on
  the gene's strand; 2.9 MB of ranges). The reader's claim holds in every
  cell: genes with an open promoter are expressed in 64% to 84% of cases
  against 2% to 22% with a closed one. The elements add nothing at this
  resolution: within a cell, open-promoter genes with activating input are
  expressed no more often than those with no element (57% to 100% on 3 to 15
  genes against 65% to 85%); across cells, for the 39 genes whose input and
  expression both vary, the cell where the elements are most active is the
  most expressed cell 15% of the time, against 25% by chance and 33% under
  shuffled cells. Twelve attributions are rejected by name (MAP3K7CL in HepG2
  with input 0.77 and no expression, SIM2 and EVA1C in K562, OLIG1 in IMR-90)
  and 43 expressed genes have a closed promoter and no active element, most
  in GM12878, whose reader has the fewest peaks. Reading: the target gene may
  be right and still not reproduce the cell, because the deletion's magnitude
  was scored in whichever tissue moved most and "active" is a DNase overlap;
  the closure asks for the deletion scored in the cell's own track. One
  measurement lesson on the way: IMR-90's ENCODE RNA-seq carries its strands
  inverted, so the module probes each cell's orientation and swaps when needed.
  **Rerun with the deletion scored on each cell's own track** (genomeos-fe
  back-filled `predicted_coding_by_cell` for the 200 named chr21 elements the
  same day): across cells no better. Ungated, over 71 genes, the most active
  cell is the most expressed 24% of the time, the shuffled rate exactly; gated
  by a DNase peak, 15% over 39 genes. Within a cell the sign carries a little:
  open-promoter genes with a repressing input are expressed less often than
  those with an activating one in all four cells (K562 56% against 63%, HepG2
  53% against 74%, GM12878 85% against 91%, IMR-90 61% against 81%), which is
  the model's direction being right where its magnitude is too small to rank
  cells. The negative stands and is sharper: 83 of 221 genes have any
  attributed element and most have one, with effects of 0.1 to 0.2 log2,
  while a gene's expression in a cell is set by its promoter and by the
  elements the two sampled runs never scored. An attribution of "this
  element, this gene" is not yet an attribution of "this cell", and the
  elements that would carry it are the ones not yet asked.
- **Syntax against values (2026-09-12).** Albert's framing, that a language has
  syntax, operators and values and a variable trait should show the same
  syntax in everyone with only the value changing, as a command: `genomeos
  syntax --gene G --chrom C` reads Zoonomia constraint per base (the syntax),
  every imported person's variants (the values) and the gene's features, and
  classes each variable position as a value or a value in syntax
  (`attribution/syntax.py`, `syntax_<GENE>`). HERC2: 211 kb, 5.6% of bases
  constrained, 160 variable positions across the GIAB trio, 3 on syntax
  (1.9% against 5.6% if variation ignored syntax); OCA2: 345 kb, 0.9%
  constrained, 544 variable, 3 on syntax (0.6% against 0.9%). Variation avoids
  syntax in both, and the second most constrained variable position of HERC2
  (phyloP 3.41, intronic) is rs12913832, the enhancer base that sets OCA2
  expression and eye colour, HG002 and HG003 carrying the alternative twice
  and HG004 once. The tool did not know the answer; the literature's value fell
  out of the reading.
- **Against measured activity (2026-09-12, genomeos-fe; step 3 closes).**
  ENCODE4's joint lentiMPRA library, 51,376 non-overlapping 200 bp elements
  assayed in K562, HepG2 and WTC11 (active at log2 RNA/DNA of 1 or more:
  7.8%, 6.6%, 8.2%), read through the registry, constraint, the measured
  reader and AlphaGenome's predicted DNase on the cell's own track
  (`attribution/mpra.py`, `predict/chromatin_tracks.py`, job
  `mpra_genome_wide`, `mpra_chr*`, `mpra_genome_wide`; 2,727 model windows).

  | layer, same cell | K562 | HepG2 |
  |---|---|---|
  | active when promoter-like / distal enhancer-like / no cCRE | 27.5% / 7.6% / 6.9% | 18.5% / 7.1% / 4.9% |
  | active when constrained / not | 8.2% / 7.7% | 7.7% / 6.5% |
  | measured reader: precision (base rate), recall | 11.0% (7.8%), 75.5% | 8.5% (6.6%), 73.8% |
  | predicted reader: rank correlation with activity | 0.443 | 0.323 |
  | predicted reader, the other cell's track | 0.030 | 0.119 |
  | top against bottom quartile active | 15.7% / 2.0% | 12.0% / 3.0% |

  Of 3,913 elements active in exactly one of K562 and HepG2, the predicted
  DNase is higher in the active line for 84.6%; the measured reader is open in
  the active line only for 50.8%. Reading: the registry orders activity as it
  should, but a distal enhancer-like element is barely more active than no
  class at all, because the reporter measures a sequence out of its node;
  constraint adds a point or two; the cell is what the two readers know, the
  predicted one more sharply than the measured one. The model knows the cell,
  not the amount: a rank correlation of 0.4 leaves most of the variance to
  promoter-like strength, which openness does not measure. Use 5 of
  ALPHAGENOME.md (the predicted reader) is built by this. The scoring step of
  area I closes: tissue (VISTA), target (eQTL), activity (lentiMPRA), openness
  (the readers) and consequence (GWAS, ClinVar) have each been asked.
- **Copies first (2026-09-12, genomeos-bb).** Curated segmental duplications
  (UCSC genomicSuperDups, 1 kb or more at 90% identity or more) read per
  UNKNOWN block (`genomeos duplications --chrom C`, `duplication_<chrom>`, a
  minute per chromosome): on chr21 the constrained_unknown tier is 60.6%
  duplicated, 15 of its 20 blocks mostly so, against 7.8% for the regulatory
  tier and 21.5% for neutral; the most constrained block of the chromosome
  (chr21:5,502,603-5,553,601, 51 kb) is 100% duplicated with 8 partners. On
  chr15 the tier is 11.2% duplicated (23 of 44 blocks mostly so) and neutral
  63%. Reading: a good part of what reads as deep mammalian constraint with no
  human score in intergenic space is a copy, where the alignment is paralogous
  and variant mapping fails; the two axes disagreeing had a cause. **With the
  genome in** (`duplication_genome_wide`, 24 chromosomes, 64,216 curated
  pairs, 166.9 Mb, 5.4% of the genome) the chr21 number is the acrocentric
  short arm, not the tier: genome-wide the constrained_unknown tier is 4.7%
  duplicated by bases, concentrated on chrY (72%), chr21 (61%) and chr22
  (42%), every other chromosome at 11% or below. What holds everywhere: 216 of
  the tier's 1,098 blocks (19.7%) are mostly copies, so "copies first,
  attribution second" applies to one block in five of the real unknown; the
  neutral tier is the most duplicated by bases (10.3%, 26% of its blocks) and
  the regulatory tier the least (2.8%, 4.4%); 2.0% of the 7,063 scored
  elements sit inside a duplication. The organiser reads `duplicated_fraction`
  at 0.5 or more as the copy flag. The classifier's own `similar_to` heuristic
  fired once genome-wide, and that pair is supported; the curated pairs replace
  it.
- **The human axis (2026-09-12, genomeos-bb).** gnomAD's Gnocchi constraint,
  the depletion of variation among 76,156 people at 1 kb, read per UNKNOWN
  block and per attributed element beside Zoonomia (`genomeos variation
  --chrom C`, `attribution/variation.py`, `variation_chr21`,
  `variation_chr15`, `variation_vista_genome_wide`; area J's step 1). The two
  axes agree where they should and add where it matters: canonical coding
  segments are 40% human-constrained against 2% for neutral blocks; on chr21's
  231 measured elements, those constrained on both axes name a gene 75% of the
  time against 60% (mammals only), 49% (people only) and 39% (neither). Two
  cautions the table carries: the constrained_unknown tier reads 1.7%
  human-constrained, no more than neutral, so deep mammalian constraint in
  intergenic blocks is not matched by human depletion at 1 kb; and VISTA
  positives are 20.2% human-constrained against negatives 16.9% over 2,223
  elements, a lean rather than a separation. Gnocchi counts whole kilobases,
  so sub-kb elements read their window. **Genome-wide** (`variation_genome_wide`,
  24 chromosomes, 137 MB of ranges): kilobases human-constrained are 38.1% of
  coding segments, 18.9% of introns, 11.9% of the regulatory tier, 3.0% fossil,
  2.7% neutral, 2.4% constrained_unknown, 0.4% structural.

  | element case (mammals × people) | elements | name a coding gene |
  |---|---|---|
  | syntax (constrained on both) | 825 | **73.9%** |
  | relaxed (mammals only) | 1,337 | 56.0% |
  | recent (people only) | 782 | 50.0% |
  | tolerant (neither) | 2,235 | 45.3% |

  The pair beats either axis alone on every count and on every chromosome
  (51% on chr1 to 87% on chr6 for the syntax elements). The
  constrained_unknown tier's 948 measured blocks split 492 relaxed, 346
  tolerant, 77 syntax and 30 recent: half the real unknown is held across
  mammals and variable among people, which is where the copies sit. chrY has no
  Gnocchi coverage and every block there is recorded as unmeasured, never zero.
- **The organiser (2026-09-12).** `genomeos organise --chrom C`
  (`attribution/organise.py`, `organised_chr*`, `organised_genome_wide`) joins
  the budget, the human axis and the copy flag per block, recomputes nothing,
  and re-reads every block copies first: a block half or more duplicated is a
  copy before anything else is said of it; a constrained_unknown block that is
  not a copy is read by the case the two axes make. Genome-wide:

  | tier | blocks | copies | after copies | syntax | relaxed | recent | tolerant | unmeasured |
  |---|---|---|---|---|---|---|---|---|
  | structural | 238 | 56 | 182 | 0 | 0 | 1 | 37 | 144 |
  | fossil | 7,302 | 1,240 | 6,062 | 0 | 0 | 732 | 5,130 | 200 |
  | regulatory | 15,536 | 689 | 14,847 | 612 | 1,011 | 3,831 | 8,932 | 461 |
  | constrained_unknown | 1,098 | 216 | 882 | **69** | 437 | 29 | 329 | 18 |
  | neutral | 2,632 | 693 | 1,939 | 0 | 0 | 306 | 1,537 | 96 |

  The real unknown after copies is 882 blocks and 30.6 Mb, and it is not one
  thing: 69 blocks (387 kb) are constrained on both axes and are the sharpest
  candidates for something unannotated; 437 (13.4 Mb) are held across mammals
  yet variable among people, a frame whose value varies or a function lost;
  329 (16.4 Mb) are free on both scales at kilobase resolution and want their
  alignment checked before anything is attributed. chr21's twenty blocks
  become five after copies, none of them syntax; chr2 has nine syntax
  candidates, chr1 and chr6 seven each, chr7 seven. The compiled program
  carries the copy flag too: a copy's region reads "copy, ..." in its role and
  the duplicated fraction and partner count in its note (116 of chr21's 446
  regions), with the checks and the unknown count unchanged.
- **The closure passes on a gene's whole regulatory input (2026-09-13).**
  genomeos-8e's job scored every one of chr21's 12,139 enhancer elements in
  AlphaGenome on the four cell lines' own tracks (102 minutes, no quota
  answer; `enhancer_targets_all_chr21`, the table kept local, 5,174 elements
  naming a coding gene), and the closure reran on the same gene by cell
  table; `attribution/targets.py` follows the committed summary to that local
  table, so no reader falls back to a sample. With 193 of the 221 genes now
  carrying a median of 11 elements, the verdict turns. The judge was made
  fair first: a tie among cells is broken at random rather than by cell
  order, and the null permutes each gene's inputs across cells 1,000 times.

  | element input across the four cells | genes | most active cell is the most expressed | null | p |
  |---|---|---|---|---|
  | tissue-agnostic magnitude, DNase-gated | 138 | 30% | 25% | 0.07 |
  | the deletion on the cell's own track | 146 | **36%** | 25% | 0.002 |
  | the cell's own score, DNase-gated | 141 | **41%** | 25% | 0.001 |

  Within a cell the whole input now carries too: open-promoter genes with
  activating input are expressed more often than those with no element in
  K562 (68% against 55%), HepG2 (80% against 67%) and IMR-90 (79% against
  57%), not in GM12878 (81% against 85%), and the unexplained genes fall from
  43 to 9, all but one in GM12878. What passes is specific: the deletion
  scored in the cell's own track, read where the cell's chromatin is open.
  The tissue-agnostic magnitude still fails, which is the earlier negative
  read correctly; a sample of one element per gene could not carry a cell,
  the whole input can. The effect is modest (mean rank correlation 0.19) and
  72 attributions are rejected by name (S100B in K562 with input 7.0 and no
  expression, TFF1 in HepG2 and IMR-90, RIPK4 in IMR-90), which is the list
  to read next. The compiled chr21 program now holds every attributed
  element: 5,972 entities (446 regions, 5,174 elements, 193 genes, 159
  domains), 5,174 rules, 20 unknowns, 2.7 MB, checked by `bio test` in 0.26 s.
- **Replicated on chr22 (2026-09-13).** The same test on the second
  chromosome, 19,708 elements scored per cell (genomeos-79's chain), 432 of
  447 genes carrying a median of 13:

  | element input across the four cells | chr21 genes | chr21 | chr22 genes | chr22 | null |
  |---|---|---|---|---|---|
  | tissue-agnostic magnitude, DNase-gated | 138 | 30% (p 0.07) | 373 | 31% (p 0.003) | 25% |
  | the deletion on the cell's own track | 146 | 36% (p 0.002) | 386 | 30% (p 0.013) | 25% |
  | the cell's own score, DNase-gated | 141 | **41%** (p 0.001) | 375 | **32%** (p 0.003) | 25% |

  The across-cell pass replicates, smaller: every input beats the null on
  chr22, the cell's own DNase-gated score best on both chromosomes. Within a
  cell the two chromosomes differ, and the difference is informative: on chr22
  repressing input lowers the expressed fraction of open-promoter genes in
  K562 (53% against 68% with activating input), HepG2 (59% against 68%) and
  IMR-90 (55% against 72%), but in K562 and HepG2 activating input does not
  raise it above genes with no element (68% against 91% on 21 genes, 68%
  against 81% on 41), so the model's repressors carry a cell better than its
  activators. In the tissue-agnostic reading 235 attributions are rejected by
  name (VPREB1, UPK3A and VPREB3 in K562 first, B-cell and urothelial genes
  the model activates in an erythroid line; 220 in the cell's own reading) and 59
  expressed genes are unexplained, nearly all in GM12878, whose reader has the
  fewest peaks. A measurement note: `MeasuredRna` now orients itself, and
  the closure reports its decision rather than probing a second time, so IMR-90
  reads "swapped" on both chromosomes.
- **The rejected attributions as a set (2026-09-13).** genomeos-8e read
  VPREB1, UPK3A and VPREB3 from the scorer's side: their tissue-agnostic
  magnitude comes from Jurkat and skeletal muscle and their K562 values sit
  near zero, so "activated in an erythroid line" was that input mode reading
  another tissue's effect as K562's. The closure's headline is now the cell's
  own score counted where a DNase peak sits on the element whenever per-cell
  scores exist (`primary_mode`), the tissue-agnostic reading kept as a
  comparison. In the headline reading the rejected set splits in two:

  | chromosome | rejected pairs | genes no line of the panel makes | genes another line makes |
  |---|---|---|---|
  | chr21 | 73 | 59% (OLIG1, OLIG2, GRIK1, S100B, KRTAPs) | 22 genes (COL6A1, COL6A2, ITGB2, JAM2) |
  | chr22 | 220 | 42% (CRYBB1, CRYBB2, CYP2D6, CLDN5) | 96 genes (APOL1, ADORA2A, BIK, ARSA) |

  Half the rejections are not wrong targets at all: neural, keratin, lens and
  liver genes with an open promoter and active elements that no line of a
  panel of erythroid, liver, lymphoid and fibroblast cells makes, because the
  tissue's own factors are missing. The rest are the genuine wrong-cell cases,
  and they share a shape: against the attributions the cell honours they carry
  more elements (57 against 26 on chr21, 33 against 18 on chr22), a higher
  summed input, and promoters open in fewer of the four lines (2.6 against 3.6,
  2.8 against 3.6). The closure is summing many weak elements into a large
  input for genes whose promoters are only marginally available, which is the
  model's systematic error and the next thing to correct.
- **The corrections, tested and rejected (2026-09-13).** Four other inputs
  on the same gene by cell table, with a promoter-only control:

  | input across the four cells | chr21 | chr22 |
  |---|---|---|
  | sum of the cell's own element scores, DNase-gated (headline) | **41%**, p 0.001 | **32%**, p 0.003 |
  | strongest open element | 36%, p 0.001 | 30%, p 0.005 |
  | mean over open elements | 37%, p 0.002 | 31%, p 0.008 |
  | sum weighted by the promoter's DNase percentile | 29%, p 0.14 | 29%, p 0.05 |
  | the promoter's DNase percentile alone, no elements | 27%, p 0.27 | 28%, p 0.12 |

  None improves on the sum, the rejected count barely moves (73 to 71, 220
  to 214), and sum and strongest element disagree in sign on 3% of pairs. The
  control is the finding: the promoter's own openness does not pick the cell
  that makes a gene, while the elements' summed input does, so the many-element
  genes are not a summing artefact to compress away and the elements carry
  which-cell information the promoter does not. The wrong-cell cases stay
  open, and genomeos-8e's per-element reading suggests why: some (ITGB2, JAM2)
  have one strongly cell-specific element that a per-element rule could
  separate, others (COL6A1) are genuine multi-line genes whose expression
  threshold rather than their input decides.
- **The 69 syntax candidates read one by one (2026-09-13, genomeos-i1).**
  `attribution/candidates.py` (`scripts/syntax_candidates.py`,
  `syntax_candidates_genome_wide`, six minutes, no model call) reads each
  block through the flanking GENCODE genes and the CTCF node, the reader and
  registry on its conserved bases, an open-frame codon score against the
  chromosome's canonical exons, the segment parser, the local proteome and
  the measured ground truth, every layer also on up to forty same-length
  intergenic control windows nearby.

  | reading | blocks | kb | mean confidence |
  |---|---|---|---|
  | regulatory element (registry on the conserved bases, or the reader above controls) | 23 | 193.1 | 0.34 |
  | promoter-like (PLS or DNase-H3K4me3 on the conserved bases) | 8 | 20.6 | 0.39 |
  | coding exon candidate (open frames beyond chance, or parser joins a neighbour) | 5 | 21.5 | 0.36 |
  | extended 3' end of a coding gene, unmeasured | 9 | 23.5 | 0.19 |
  | unexplained | 24 | 127.8 | 0 |

  31 plausible regulatory, 5 unannotated coding and 9 possible unannotated
  RNA, 24 unexplained; no copy of a known protein, no structured-RNA copy, no
  primate-repeat artefact. Against the controls the regulatory count is the
  weak one: a registry element sits in 54% of control windows (candidates
  70%, p 0.049), openness is not above chance (2.78 against 2.33 cells, p
  0.15), and the candidates' conserved bases are less often a registry
  element (7.6% against 19.8%) and less often open (2.2% against 3.7% per
  cell) than conserved bases nearby. The coding signal stands: 13.2% of 227
  conserved segments are coding-like against 5.6% of 2,711, and a parser exon
  lands on conserved sequence in 17.4% of blocks against 4.3% (p 0.024; the
  parser recovers 13.7% of known internal exons). Strongest: promoter-like
  at FAM237A with the one active lentiMPRA element (K562), enhancer-like
  chr7:25.26 Mb and chr22:38.57 Mb open in 11 and 10 of 11 cells, and the
  long_orf block chr7:55.59 Mb with 9 of 11 segments exon-like (p 3e-10) and
  a seven-exon parser structure and no known protein behind it. Nine blocks
  may borrow their human axis from a neighbour's exons. ATTRIBUTION.md "The
  syntax candidates read one by one".
- **What the whole-chromosome scoring bought the 98% (2026-09-13, genomeos-i1).**
  Decided on chr21 from the tables on disk (`attribution/unknown_scoring.py`,
  `unknown_scoring_chr21`, no model call). 2,334 of 12,139 scored elements
  overlap an UNKNOWN block (1,155 move a gene; only 17 over 5 of the 20
  constrained_unknown blocks). The ablation: the closure rebuilt from its
  committed table reproduces 40.7% (p 0.001); without the 665 unknown-block
  elements of its input it keeps 40.5% (random removals of the same size
  39.5%, rho 0.191 to 0.213, rejections 73 to 67); with only them it falls to
  28.7% on 61 genes, p 0.26 (random subsets 35.8%). Inside 480 strata of
  length, GC and distance to the nearest coding TSS they move a gene less
  often (49.5% against 68.0%), name a coding gene less often (28.6% against
  43.1%) and act in fewer cell lines (0.45 against 0.80), p 0.0005 each; the
  magnitude barely differs (0.309 against 0.338, p 0.048) and the silencer
  share not at all, and the rest is 99.9% intragenic, a confound no stratum
  removes. chr21 has no syntax block; mammal-constrained blocks show nothing
  that survives a block-level null; of the 69 candidates the two on complete
  chromosomes gained one weak target (DMC1, 0.108 log2) and no cell. The
  model agrees with measurement where they meet (signed effect against
  lentiMPRA activity, rho 0.22 K562 and 0.24 HepG2, n 329, p 0.0005; 59 of
  them over unknown blocks, too few to split). Coverage: unknown blocks with a
  named coding target went from 11 to 140 of 446 (90 with a cell), but the
  naming elements cover 185 kb, 0.9% of the unknown space. Finishing the
  sweep is about 674,000 requests at the measured 0.76 per element; scoring
  the 882 non-copy constrained-unknown blocks and the 69 candidates directly
  is a few thousand. **Verdict:** for area I the sweep does not pay; stop it
  as an instrument for the 98% after chr20 and spend the next quota on the
  constrained unknown and the candidates directly. ATTRIBUTION.md "What the
  whole-chromosome scoring bought the 98%".
- **The lexicon, chr21 and chr22 (2026-09-13, genomeos-i1).**
  - **What it is.** One index of the genome's units at six levels, built per chromosome from
    the sequence and the results on disk: k-mers per context, recurring 16-mer seeds and
    their families, curated units, blocks, nodes and libraries. Every count carries a GC by
    repeat-share stratified expectation, Poisson tails on occurrence-equivalents,
    Benjamini-Hochberg and the number of tests.
  - **Where it lives.** `attribution/lexicon.py` and `scripts/lexicon.py`; summaries in
    `lexicon_chr21` and `lexicon_chr22`; the index under `data/knowledge/lexicon`.
  - **Tests.** 22,076 and 22,795 context tests, 2,172 and 2,279 passing.
  - **Seeds.** The vocabulary found from sequence is the repeat library: 1,965 and 1,980 of
    the 2,000 seeds are RepeatMasker and the rest are tandem words, with no recurring
    segment that is neither. **Negative.**
  - **k-mers.** The top words are CpG words in every context, the heterogeneity an order-2
    chain cannot absorb. **Negative as a lexicon.**
  - **Fossils within a node.** They are alike against a label permutation (p 0.044 and
    0.001) but not against the same partition shifted (p 0.47 and 0.11, combined 0.21):
    proximity, not the CTCF node. **Negative.**
  - **Fossils and the node genes' libraries.** 4 and 5 pairs pass against 7 to 16 and 8 to
    23 for random gene sets. **Negative.**
  - **Syntax against value slots.**
    - 21 and 22 units are syntax by thresholds, and 12 and 12 survive the context-matched
      null.
    - Every JASPAR factor among them fails against its own element (0 of 22 pass; pooled,
      0.335 against 0.308).
    - What survives is known: promoter cCREs, one rRNA copy and one (GATG)n tandem.
    - Value slots: 1 and 0. **Negative:** the axes lack the resolution. Conservation is
      element membership, not per-base phyloP, and the human axis is block-level.
  - **Genome-wide cost.** About 1.4 s per Mb (71 minutes for the genome on one core), a
    5 GB peak for chr1, about 29 MB of index and under 1 MB of summaries. ATTRIBUTION.md
    "The lexicon".
- **Many human genomes at once (2026-09-14, genomeos-h1).** Albert's question, what is common,
  what varies among a few values and what is dispensable, asked of 89 human genomes rather than
  of an aggregate score. `attribution/human_panel.py` (`scripts/human_panel.py`,
  `human_panel_chr21`, no model call) streams the HPRC release-1 Cactus alignment of 90 human
  assemblies (UCSC `hprc90way`; hg38, T2T-CHM13 and 88 haplotypes of 44 people) by the block
  offsets the UCSC API returns and distils it per block into presence and status of every
  assembly, every variable column with each assembly's allele, and inserted bases: chr21's
  3.28 GB in 236 s into an 8.4 MB store under `data/knowledge/human_panel`, 1,054,284 variable
  columns, 424,828 events. Beside it: the HPRC arrangement and v2.1 SV tracks, gnomAD v4.1.1
  frequencies, the TRExplorer repeat histograms, GTEx DAP-G eQTLs and MPRAVarDB allele pairs,
  and, as the mutation-rate control, ENCODE Repli-seq replication timing for 11 lines lifted
  from hg19 per kilobase (the project's first replication-timing data).
  - **Controls pass, and bound the claim.** Canonical CDS carries 0.41 of the recurring
    variation a GC- and timing-matched background predicts, the neutral tier 0.91; 200-bp
    coding units are fixed 70.3% against a matched expectation of 44.1%, neutral units 48.2%
    against 47.6%; 58% of callable genes are core against 8% of neutral blocks and 14% of
    length-, GC- and timing-matched control windows. Not "overwhelmingly": 64 of 153 callable
    genes are variable, because synonymous variation is real.
  - **Replication timing survives the control.** Late DNA carries up to 1.5 times the
    background rate within a GC stratum, but 82% of the coding-against-neutral gap and 98% of
    the coding units' excess of fixed windows survive matching on timing, and unit classes are
    spread evenly over timing tertiles. What is late is what cannot be placed (74% of unplaced
    units). The regulatory tier's apparent core blocks were mostly timing (9.3% to 3.5%).
  - **The three catalogues of the unknown space** (200-bp units, 20.8 Mb): fixed 24,272 units
    (4.85 Mb), storage 26,969 (5.40 Mb), cannot place 22,055 (4.41 Mb, of which 774
    hypervariable). Storage domains are small (11,673 two-value, 9,024 three-value, median
    effective 1.68) and 81.9% of their events have a gnomAD frequency (panel minor share
    against gnomAD MAF, Spearman 0.69). 2,216 storage units carry a GTEx eQTL or an MPRA
    allele pair, 1,053 of them early-replicating with at most three values: the shortlist for
    the executor test. **But storage is the default**: half of all placed units on the
    chromosome hold a recurring alternative and no unknown tier holds more than its matched
    background.
  - **Calibration.** rs12913832 (HERC2/OCA2) and rs4988235 (LCT) come out as two-value
    columns, the DRD4 and SLC6A3 VNTRs as copy-number columns with the textbook alleles, ABO
    as a two-value column by length and ten-value on substitutions, HLA-A as a
    twenty-value slot; the INS class I/III VNTR fails, because the alignment breaks a 2-kb
    expansion into replaced blocks.
  - **Against Gnocchi** (area J's axis): weak agreement, Spearman -0.08 over 34,392
    kilobases, a Gnocchi-constrained kilobase depleted in the panel 25% of the time against
    18% for the rest. Gnocchi's constrained kilobases are early and GC-rich; the panel's are
    late and GC-poor, where genealogy dominates (counts are 15.6 times overdispersed, so 6,374
    depleted kilobases stand against 1,091 expected and the matched windows, not the Poisson
    tail, are the test). Where Gnocchi is silent and the panel reads depleted, 58% of the
    kilobases are segmental duplications with 45% of assembly-bases missing. At block scale
    the two human axes disagree: of 13 core blocks with a case, 3 are human-constrained.
  - **The 69 candidates** carry 1.21 times their own flanks' recurring variation (median
    1.03), 6 core and 59 variable: blocks held across mammals and by Gnocchi are not invariant
    among people at base resolution. The lexicon's negative, sharpened.
  - **Cost.** The genome is 266.6 GB of MAF, 5.3 hours of streaming at the measured 13.9 MB/s,
    about 1 GB of stores, 5 to 6 GB of memory for chr1, and 16 GB of gnomAD ranges for the
    storage units. The chr21 run took 8.6 hours end to end, almost all of it in the unprofiled
    regional calibration and candidate stage, which is what to fix before a genome run.
    ATTRIBUTION.md "Many human genomes"; the source row in GRAMMAR-BY-COMPARISON.md is now
    read rather than pending.
  - **Slow stage fixed, chr22 replicates (2026-09-14).** The 8.5 hours were TLS handshakes,
    one per range request (0.5 s quiet, 4 to 17 s on UCSC's loaded host), not analysis (0.1 s
    per locus). Kept-alive connections with a mirror fallback, readers opened once, and local
    Repli-seq copies bring chr21's whole run from 31,099 s to 1,392 s with every committed
    number identical; calibration is 91 s, the 69 candidates 2 s, and gnomAD transfer is what
    remains. chr22 (`human_panel_chr22`) reproduces the headline results. Coding exons read
    0.423 of the matched rate (chr21 0.411); 71% of the coding-neutral gap survives timing
    (82%); half of all units hold a recurring alternative; no unknown tier is more fixed than
    matched background, and chr22's blocks are core less often than matched windows in every
    tier. The two human axes barely agree (Spearman -0.02) and counts are 26-fold
    overdispersed. Three things do not replicate: on chr22 variable units lean late (37%
    against 31%), Gnocchi's constrained kilobases are not strongly early, and core blocks agree
    with Gnocchi at block scale (45% against 28%). chr22's 10 non-copy constrained-unknown
    blocks are all variable against their flanks; its syntax candidate chr22:38,571,213 is the
    most depleted (ratio 0.38, tail 0.0101). The view passed to area J: trust Gnocchi for
    constraint among people, phyloP per base, and the panel for structure, alleles and value
    domains, pooled or against matched windows only. ATTRIBUTION.md "replicated on chr22".
- **The lexicon re-asked at base resolution, chr21 and chr22 (2026-09-14, genomeos-i1).**
  - **What changed.** `attribution/lexicon_axes.py` (`scripts/lexicon_axes.py`,
    `lexicon_axes_chr21`, `lexicon_axes_chr22`, no model call) replaces all three soft
    axes:
    - Zoonomia phyloP per base, streamed once (97 MB and 4 minutes per chromosome) into a
      23 MB byte track;
    - the HPRC panel's recurring substitutions per base, read from genomeos-h1's stores
      with a trinucleotide by context by GC expectation (Gnocchi is joined but not
      trusted);
    - JASPAR at 0.95 over the whole chromosome (13.5 and 11.6 million merged sites), with
      column-permuted decoys scanned the same way.
  - **Nulls and calibration.** Every (unit, context) row is tested against its strata,
    against the phyloP bin and against its own flanks. Variances are built from the
    unit's own occurrences, BH is applied per family, and each family is calibrated on
    shifted copies (0 of 7,088 and 2 of 7,183 pass the flank tests among people).
  - **Mammals see motifs.** Sites beat their own decoys in every context: 549 of 743 and
    571 of 749 factors, median site-over-flank difference +0.15. Against the rest of
    their own cCRE, 334 of 690 and 405 of 723 factors pass, against 37 and 69 decoys.
  - **People barely do.** The paired median difference is -0.018 and -0.024. Pooled
    against their own element, sites give 0.986 and 0.975 of expected variation, against
    0.999 and 1.006 for decoys. Among rows held across mammals against flanks, the
    people-depleted share is 9.0% for motifs against 9.2% for decoys, and 4.5% against
    7.9%.
  - **Question 2.** No syntax class separates.
    - Decoys pass the people test nearly as often as motifs (192 against 233, 189 against
      193).
    - The one value-slot row that replicates is the dELS class (more conserved than its
      flanks, z 12.5 and 12.7; more variable among people, z 6.7 and 7.9). Open chromatin's
      mutation rate is not excluded.
    - **Negative, and stronger than the first.**
  - **The fossil tier.** Matched on replication timing, fossil copies of a subfamily are
    as variable as the same subfamily's intronic copies (0.994 on chr21, 1.086 on chr22)
    and no less variable than regulatory-tier copies (0.946 and 1.001). The chr21 hint of
    motif sites less variable inside fossils (-0.05) reverses on chr22 (+0.038).
    **Negative:** a fourth independent failure on that tier, after node libraries,
    methylation and h1's variability.
  - **The two human measures do not agree at unit scale.** Rank correlation of Gnocchi
    with panel depletion is -0.067 on chr21 and 0.013 on chr22.
  - **Candidate.** chr22's one candidate carries 4 recurring substitutions against 19
    expected.
  - **Decision.** The other chromosomes are not built. ATTRIBUTION.md "The lexicon at base
    resolution"; LESSONS.md "The UNKNOWN space".
- **Compression, chr21, chr22 and chr18 (2026-09-14, genomeos-m1).** Albert's probe of which way
  of classifying the genome buys the most per unit of effort, done the information-theoretic way.
  Each piece of knowledge becomes a model of the next base, and its worth is the bits it saves net
  of its annotation. The tool is `attribution/compress.py` (`scripts/compress.py`,
  `compress_chr21`, `compress_chr22`, `compress_chr18`), with no model call and one CpG-island
  request per chromosome.
  - **Honesty.** Static models are fitted on another chromosome; adaptive models pay as they
    learn; annotations are charged; a layer is inactive outside its claim. A control primes the
    generic models with the duplication partners' sequence.
  - **Baselines, bits per base, chr21 / chr22 / chr18.** xz 1.706 / 1.673 / 1.700. Naive Markov
    mixture 1.692 / 1.633 / 1.745. Generic with blind copy finding 1.590 / 1.545 / 1.650. Every
    layer with its annotation paid 1.532 / 1.499 / 1.665. The best net is primed generic plus
    duplications, 1.500 / 1.446 / 1.630.
  - **What pays.** Segmental duplications: +85 / +91 / +18 kb per Mb, two thirds of it from
    access to the partner's sequence, the alignment +0.22 to +0.33 bits per claimed base. GC and
    CpG: +7 / +7 / +6 kb per Mb.
  - **What does not. Negative.**
    - Repeats, over blind copy finding: -13 / -20 / -27 kb per Mb, 31 bits per labelled copy
      against 0.04 to 0.07 bits a base saved. This one is set by the training: fitted on two
      chromosomes, chr21's repeats turn to +3.8 kb per Mb, and still lose 12.5 left out of
      the full stack.
    - Coding phase: +0.03 to 0.05 bits per coding base against 0.14 to point to it.
    - cCREs: 0.01 bits per base or less against 0.08.
    - Motif sites: 0.5 bits per site base against 1.6 to address them.
    - Tiers: nothing over the repeat-aware stack.
  - **The unknown space. Negative, loudly.** Unique sequence costs 1.89 to 1.94 bits per base
    under every model in every tier. chr18's constrained_unknown (2.87 Mb, 43 blocks, 1.5%
    duplicated) is 0.0105 bits per base harder than neutral (0.0073 to 0.0146, 38 against 82
    blocks): the predicted sign at half a percent of the alphabet. Fossil repeat bases still cost
    1.55 bits per base. chr21's tier compresses easier than neutral only because 61% of it is
    copies.
  - **Cost.** 28 to 29 s per Mb and 7.5 to 10 GB peak. A genome run is about 23 h as run, or
    6 to 7 h for a production pass, with chr1 and chr2 needing chunked hashing on this machine.
    ATTRIBUTION.md "Compression".
- **Compression genome-wide, 22 of 24 chromosomes (2026-09-14, genomeos-m2).** The probe's blocker
  was memory, and the fix is committed (a097bc3): context hashing by doubling and counting in hash
  groups, one mixture for all 33 stacks instead of 47, and the chromosome mixed 32 Mb at a time.
  chr21 and chr22 reproduce their committed runs to four decimals; chr18 moves by at most 0.0028
  bits per base because the pass fixes the mixing window at W = 8 where chr18's grid had chosen 16,
  and no verdict changes. Then the reduced pass over the genome (`scripts/compress_genome_wide.py`,
  `compress_pass_<chrom>`, `compress_genome_wide`): 2,467 Mb, 6.8 s per Mb, 10.2 GB peak, 4 h 50 m,
  no model call.
  - **The tier gap holds, and is still tiny. The result the run was for.** Constrained_unknown
    against neutral on unique bases, over 901 blocks on 22 chromosomes instead of chr18's 43:
    +0.0117 bits per base (0.0095 to 0.0145) under the full stack, against chr18's +0.0099 (0.0064
    to 0.0136). Seventeen times the blocks moved the estimate by 0.002 and halved the interval.
    1.9299 against 1.9182 bits per base is one part in 170, half a percent of a two-bit alphabet.
    The 32 Mb the project calls most interesting is not information-rich; the sign is right and the
    size is a rounding error.
  - **The tiers stop paying at genome scale. Negative, and new.** +0.019 bits per claimed base over
    generic on chr21 became -0.0003 genome-wide, -0.0024 over the repeat-aware stack and -0.0038
    left out of the full stack. What made them look positive was the acrocentrics.
  - **Fossils are level with neutral genome-wide. Negative.** +0.0013 (-0.0013 to 0.0044) over
    5,350 blocks.
  - **Still paying:** duplications +0.541 bits per claimed base (+21.6 kb per Mb, a third of
    chr21's +85 because the genome is 5.7% duplicated, not 12.6%), GC and CpG +4.4 kb per Mb.
    **Still not:** repeats -32.7 kb per Mb (worse at scale), cCREs -7.3, coding -1.1, motifs -0.6.
  - **chr1 and chr2 did not run. Negative, and the machine's fault.** The disk floor stopped them
    four times, and the cost is measured: counting one adaptive order over chr2 consumed 17.5 GB of
    free space (35.1 down to 17.6), of which only 5.8 GB is the count rows, so chr2 needs about
    38 GB free to start and this 460 GB disk at 94% offers 30 to 35. Spending memory instead
    (--rows-memory-gb) does not help, because 19 GB of RAM is already 11 GB spoken for. They hold
    197 of the 1,098 blocks and 5.3 of the 32 Mb, which cannot overturn the interval. The chain
    resumes on them when there is room. ATTRIBUTION.md "Compression genome-wide".
- **The all-elements sweep is complete, and the node claim is now its final size
  (2026-09-16).** Every chromosome scored: 961,227 enhancer elements deleted inside
  their node in AlphaGenome, 778,780 requests, one scorer at a time behind the key
  lock. Genome-wide the node beats as many randomly placed boundaries: **0.754 against
  0.725, +2.88 points over 440,377 coding-target elements, ahead on 18 of 24
  chromosomes (sign test one-sided p 0.011)**. *(Re-audited 2026-09-27, NODES-READER-WRITER.md: the
  figure is **+2.90, 95% CI +2.03 to +3.81**, on a control matched on boundary count alone; across
  four defensible baselines it runs **+1.2 to +6.6**, the published one the least matched. "18 of 24"
  counts chromosomes above 0.700, a threshold in no file; against the 0.725 this sentence names it is
  **16 of 24, p 0.076**, and 19 of 24 against each chromosome's own control. Held against 661
  measured CRISPRi pairs the excess is **+5.89, CI +3.18 to +8.46** — the sign stands on measured data;
  the magnitude above was never identified.)* Real and small; the headline 90.2% of
  the 4,800-element sample was mostly what random boundaries give too. The large
  chromosomes pulled every rate down (inside the domain 0.772 at 18 chromosomes to
  0.754; agreement with nearest TSS 0.620 to 0.611).
  - **Two rules about where the excess lives were fixed before the data and tested
    on it.** The sign rule (the node wins above 5.8 boundaries per Mb) ends 10 right,
    1 wrong, 1 refused; its miss, chrX, carries an excess of -0.07 points, which is
    none. The measurability rule fixed at fold 18 (an excess is distinguishable from
    random iff the caller cuts at 6.55 per Mb or finer) predicted all six unseen
    chromosomes measurable and got **4 right, 2 wrong: chr7 (7.10 per Mb, p 0.24) and
    chr5 (6.85, p 0.07). FAILED as written.** Both misses are positive excesses that did
    not reach p 0.05. The per-chromosome scatter does not support a third rule; the
    genome-wide number is the claim.
  - **The panel at 23 of 24.** chr1, chr3 and chr4 landed. The claims that held at
    twenty still hold: coding exons at 0.387 of the matched rate (median), genes core
    58.6%. The constrained-unknown depletion weakened and did not break: below the
    matched rate on 16 of 22 (p 0.026), below the neutral tier on 18 of 22 (p 0.002),
    and 13 of 16 either way on the chromosomes that pass every control (p 0.011);
    chr4 (1.127) reads against it. **Seven chromosomes fail a control**, and on all six
    autosomes and X the failing check is the neutral tier missing its GC- and
    timing-matched rate (chr1, chr8, chr16, chr17, chr20, chrX); chrY fails the coding
    checks, as a gene-poor chromosome should. The checks were written after chr21's
    first read, so they describe an ordering rather than predict one. 174,210
    executor-ready units. chr2 is unread: the session reading it ended at the weekly
    usage limit, not on the disk.
- **Next.** 0. Done 2026-09-14 (the bullet above): the lexicon's axes at base resolution
  answered Question 2 negatively, so the genome-wide lexicon stays unbuilt. What could still
  sort sites into syntax and slots is measurement, not diversity: saturation mutagenesis and
  MPRA allele effects at motif sites (the satmut and MPRAVarDB layers), and per-base gnomAD
  frequencies for population depth. 1. What remains of the candidate reading: measured RNA (ENCODE
  total RNA-seq, as the closure reads it) over the 69 blocks and their
  controls, to test the 5 coding and 9 3'-extension readings; then, with the
  next quota (the recommendation above), the deletion of the 882 non-copy
  constrained-unknown blocks and the 69 candidates' conserved segments with
  matched control windows. The closure's replication on chr19 needs no quota. 3. The
  remaining attribution: AlphaGenome's predicted chromatin tracks for activity
  and tissue on the 15,536 regulatory blocks, in-silico mutagenesis for the
  bases that matter; the self-hosted model gates the move from chromosome 21
  to the genome. 4. Scoring against measured ground truth the way the parser
  was scored against GENCODE: VISTA done genome-wide, GTEx eQTLs done for the
  target gene, GWAS and ClinVar done as enrichment tests with a shifted
  control, weak by construction (all above); ENCODE4 lentiMPRA (K562, HepG2,
  WTC11) for activity built and running, same cell against the model's own
  track and cross-cell reported apart (chr21: rank correlation 0.39 same-cell
  against 0.07 cross-cell for K562, 0.35 against 0.13 for HepG2; predicted
  DNase higher in the active cell for 74% of 47 cell-specific elements) and
  done genome-wide (above). Step 4 is closed. 5. Closure at the cell level (a deleted element
  moves its gene and the cell program; IMPC, DepMap) and the organism level
  (the Body runtime still produces a human's cell counts over decades); the
  gene level ran once (above) and returns after step 1. 6. Simpler genomes as
  the comparison, C. elegans first, then Drosophila enhancers, fugu as the
  compact vertebrate, mouse for transfer. 7. The human panel's next step
  (genomeos-h1): the regional stage is fixed and chr22 replicates (above); the executor test is done
  (2026-09-14, genomeos-x1 after genomeos-h1): swapping a stored value that fine-mapping
  calls causal moves the model's read-out of the measured gene in the measured direction
  0.679 of the time, against 0.411 at matched units storing a value nobody has measured
  (+0.268, p 0.0021); a value that is only linked to a cause moves it +0.062, and the
  pre-registered hold-out gap of +0.206 (one-sided lower bound +0.058) says the agreement
  follows causality rather than the neighbourhood, which is the reading the pre-registration
  named as the one that would have mattered more. Weak points stated: E3's own difference is
  small and positive rather than zero, everything rests on one model whose ENCODE training is
  kin to GTEx, and E1, the MPRA endpoint that would answer that objection, has 20 allele pairs
  on chr21 and chr22 and needs more chromosomes. Next: E1 widened to every chromosome, which it no longer waits for (23 of 24 panel
  chromosomes and all 24 swept, 2026-09-16), and the
  882 non-copy constrained-unknown blocks and their flanks read across the panel on the
  chromosomes that carry them. A repeat-aware reading of long VNTRs is the instrument's
  one known failure.
- **The one-target-per-element limit is the compact table's, not the sweep's
  (2026-09-22, lane-cache, `df6c5f4` registration, `bda8a83` fix).**
  `crispri.py`'s `deletion_drop` scored a **structural zero** whenever the
  measured gene was not the element's top predicted target — and a zero there
  means "the table had nothing to say", while every consumer read it as "the
  model predicts no effect". It was **8,991 of 9,237** covered K562 training
  pairs (97.3%) and **1,704 of 1,744** held out (97.7%); 78 of the 114 held-out
  positives were structural zeros. The same sweep's per-element response cache,
  written at `threshold=0.0` with every gene in the 1 Mb window signed and per
  cell line, answers **62.1%** and **65.3%** of them, so the fix cost **0
  AlphaGenome requests** and 107 seconds.
  **The registration expected the gain to shrink or hold, and committed that
  before the re-run. It nearly doubled on the arm that has the pairs and fell by
  two thirds on the arm that has fourteen**, and both are reported in one table
  at the same size. Held-out K562 AUPRC **0.633 → 0.691**, gain +0.083 (+0.031
  to +0.166) → **+0.141 (+0.082 to +0.231)**; elements not in training 0.587 →
  0.646; leave-chromosome-out 0.674 → 0.740. **GM12878 falls, 0.962 → 0.900,
  +0.097 → +0.035 with the interval crossing zero** — on 14 regulated pairs, 6
  of them structural zeros, an interval that already touched zero in September.
  The verdict stands on the same pre-registration, the same split and the same
  200 resamples; only where one feature reads its number changed.
  **Half the registered reasoning was right and the half that was wrong is the
  interesting half**: `top_target` is untouched and its weight barely moves
  (1.284 → 1.329), so the old gain was indeed carried by the indicator — but
  calling the newly visible magnitudes noise was wrong, and the weight on
  `deletion_drop` rises **15.94 → 25.30**.
  **Wherever the compact table and the cache both speak they agree exactly: 0
  disagreements over every covered pair on all 24 chromosomes.** The cache is
  the same run with the censoring removed, not a re-score. The 2026-09-16
  section of ATTRIBUTION.md is annotated rather than rewritten, including its
  "Limits" sentence, which said a full per-gene table "would cost the sweep
  again" — it would not have, and it was already on disk.
- **The CRISPRi result made showable, one cell type (2026-09-27, lane-publish, `fe61c36` registered,
  `42d7b1b` scored, 0 requests).** Held against Gschwind et al. 2026 (Nature, doi:10.1038/s41586-026-10781-4)
  Supplementary Table 3 on the same pairs with the benchmark's own estimator, which reproduces the
  published distance figures exactly (0.4359, 0.3631) before any model was scored. **Negative first: the
  held-out comparison flatters this project** — every one of the 190 held-out positives sits in an
  H3K27ac element, 1,438 of the 4,188 negatives do not, and the activity term here reads H3K27ac while
  the published held-out model reads DNase only; the registered "above ENCODE-rE2G" band is withdrawn.
  Training, leave-one-chromosome-out: baseline **0.507** (below ABC, 0.565) → frozen model with the
  deletion **0.724** (ENCODE-rE2G 0.662, Extended 0.737). Held-out, pooled and weighted as registered:
  0.567 → 0.677, gain +0.110 (+0.061 to +0.174). Post hoc, DNase only: 0.476 → 0.639 against
  ENCODE-rE2G 0.556 (0.468–0.631), **"in the range of ENCODE-rE2G"**, unpaired. The gain is invariant
  to coverage: all three registered arms pass (coverage-matched median +0.145, 1,000 of 1,000 draws
  above zero; a covered-or-not flag alone +0.001). Second cell type: only HCT116 has the n (34
  regulated, AlphaGenome carries its own tracks), costed at **705 requests** and registered with its
  pass rule, not run. Prior art: the AlphaGenome preprint already added an input-gradient score to
  ENCODE-rE2G on this dataset; new here is the deletion form and the frozen held-out test.
  Outside-reader summary: docs/CRISPRI-RESULT.md.
- **The last one-target consumers moved, and the unknown's control (2026-09-27, lane-onetarget,
  `b7352be` registered, `d717b28` result, 0 requests).** motif_transfer, syntax_tiling, candidates and
  loci `read_gene_input` read `ElementResponses` beside the compact head; 0 disagreements between the
  cache head and the table over 849,469 elements. **Negative first: the loci window rule is
  falsified** — crediting every coding gene at the bar drops strict hits 10 → 9 of 17 (HOXD lost to
  EVX2, a narrow table win before), so `gene_input` keeps the one-target sum and carries
  `window_reading` beside it. motif_transfer reproduced to the digit (AUROC 0.5422 stands);
  candidates held its falsifier (top node gene changed at 3 of 58). `constrained_unknown_targets`
  now carries its own length-matched random-window control: **real unknown 62.3% of 531 blocks name a
  target against 86.0% of 44,100 random windows, −23.7 points, below chance**; naming a coding gene
  29.8% against 67.0%. The neutral tier reads −19.3, so most of the gap is organiser blocks sitting in
  gene-poorer sequence than their control, not something peculiar to the unknown; inside length
  deciles the unknown is +4.6 points above neutral. The borrowed "87% of random windows" is withdrawn.
  This confirms, on a second script, the reading that holds milestone 1.3's second clause as not met.
- **How independent that measured confirmation is, component by component (2026-09-28, lane-nodeindep,
  `d3f51f0` registered, `366db9a` repair, `f752e5e` result, 0 requests).** Eleven components audited:
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
- **The Progress tab computed from files (2026-09-28, lane-progress, `084e092`).** `/api/state` answers
  four questions from files at request time, so the tab cannot lag the plan: **milestones** from section 6
  through `roadmap.parse_milestones` (6 of 7 reached; 1.3 held with the reason its own row states; every
  reached row carrying the caveats it emphasises, including 1.0's unpriced BioForge answer and 1.1's
  74%/65% being about 70% DNase depth); the **external review R1–R9** from section 5 item 11 with each
  acceptance test and the follow-up rows that closed it, an item counting as done only because a
  follow-up row says so; **README's Status as five separate claims**, with assay coverage read from
  `unknown_coverage.json` and the one independent prediction from `crispri_published.json`, each with its
  file's date and its file's own qualification, and the node independence audit's verdict carried beside
  the claim it qualifies; and **what waits on the owner**, from section 4's owner column, the blocked
  steps that name him and the board. **Agreement between README and each result file is computed and
  shown, so a README figure that stops matching its result reads as a disagreement rather than as the
  truth.** Four tests hold it, including that neither the server nor the page contains any of these
  figures as text, so a hand-written number fails the suite, and that the view renders with an empty
  board and no results with every figure absent rather than invented. Nothing is stored. *(The lane also
  found that section 6's table carried an unescaped `|` inside 1.1's proof cell, which truncated that
  cell for the parser and for any renderer; escaped here.)*
- **The reader under the one-gene tables, and the second consumer moved
  (2026-09-22, lane-reader2, `45e4f92` and `350c4c3`, 0 requests).**
  `attribution/targets.py` gained `ElementResponses`, the general form of
  `crispri.ElementCache`: ask what the sweep predicted for a gene at an element
  on a cell's track and get a signed value **or a named silence** — not in the
  window, not on that track, not cached — **never a zero standing in for a
  missing answer**. It reproduces the compact tables exactly where they speak,
  6,239 of 6,239 elements including the log2 fold change, and costs one
  chromosome at a time (chr21 0.6 s, chr1 5.7 s and a peak near 2.9 GB; the tree
  is 775 MB and is never read whole).
  **`eqtl.py` is the second consumer, and the registration was wrong in a way
  worth more than being right.** It predicted the eGene rate would FALL when the
  1,353 elements the threshold had been skipping were asked too. It moved by
  **one thousandth**: uniform 0.529 (2,372 elements) → **0.528 (3,613)**,
  constrained 0.448 → 0.449, VISTA 0.412 → 0.423. The decomposition shows it was
  wrong twice by the same amount — dropping the 75 elements with no answer
  available is +0.017, adding the 1,316 unnamed at 0.496 is −0.018.
  **So the falsifier fired and the threshold got priced instead: `MIN_EFFECT` =
  0.1 excludes a third of the elements to buy five points** (54.6 against 49.6,
  chance 12.7; 7.5 and 7.9 points on the other two sets). And the window the
  compact table could not express says something new: **the measured eGene's
  median rank is 1 of 34.1 genes, top-3 0.776, top-5 0.851** — when the model is
  wrong about which gene, it is usually wrong by one or two places.
  **A distinct defect class was separated on the way, and it is not the CRISPRi
  one:** a compact `predicted = null` is a **threshold**, not a silence. The
  sweep scored everything and nothing cleared 0.1 — all 3,045 null elements have
  a head below `MIN_EFFECT`. CRISPRi's zero meant "the table had no row"; this
  means "the table had a rule", and **any consumer conditioned on `predicted` is
  quoting a rate conditional on that gate**.
- **Next, open: the rest of the one-target consumers, all at 0 requests.** About
  twenty modules read the compact table through `attribution/targets.py:25`
  (`run_elements`). Highest value first: `eqtl.py:246`, which asks whether the
  model's target is an eGene of one gene per element and whose answer the web UI
  displays at `index.html:1661`; `target_calibration.py:173`, where the whole
  calibration is conditioned on `top_target`; `motif_transfer.py:1025` and
  `syntax_tiling.py:252`, where an element that moves a non-top gene reads as
  "did not move"; then a cache-backed reader beside `run_elements` itself, which
  is the broadest structural fix. `scripts/crispri_contact.py` needs the same
  re-run but costs 4DN Hi-C range requests rather than model requests, and its
  stored result now says which reader it used.
- **The calibration's `top_target` gate, counted and removed — and the number
  that matters is not the reliability but what the curve quotes off its own
  population (2026-09-22, lane-calib, `839d3c5` registration, `5881d45` result,
  0 requests).** The gate admitted **245 of 8,796** training pairs at a 76.7%
  base rate and excluded **8,551 at 2.91%**; two thirds of the excluded carry a
  number the sweep did predict, one third are genuinely outside the scorer's
  window and now say so by name. The gate turned out to act in **three** places,
  not one: the magnitude (`score()` called `crispri.annotate` with no cache, so
  today's fix never reached this module — `deletion_drop` was a structural zero
  on 97.6% of training pairs), the element (the fallback ranked candidates on
  the magnitude predicted for *whatever other gene the table named*), and the
  population (245 pairs band all **612,323** sweep targets).
  **The registered direction was wrong and is reported first: reliability
  improved and crossed its threshold** — 6/10 bins to **7/10** on the same
  held-out 1,715 pairs, AUPRC **0.559 → 0.677**, Brier 0.0418 → 0.0345. The
  registration's own improvement protocol is what stops that becoming an
  upgrade: a rise counts only beside a non-narrowing predicted range, a
  non-falling AUPRC **and a prevalence gap that closed**. The range widened and
  AUPRC rose, so it is a real ranking gain and not a flattening — but the gap
  did not close: the fixing shift goes **+0.6793 → +0.7025, larger**, and all
  three remaining bin failures are in the same direction as before. **The
  2026-09-17 verdict of failed stands.**
  **The clause registered as mattering most returned the sharpest result in the
  area.** The 245-pair curve that prices every predicted target, applied to the
  held-out pairs its own gate excludes, quotes **0.4126 against a measured
  0.0603 — ×6.85** — while landing at ×0.93 on the 40 it admits. Re-fitting
  through the window barely moves it (0.3855, ×6.40), so **the problem is the
  population, not the features**. That is the registered "the failure is
  happening" signature, unambiguously: a curve fitted on the calls a gate admits
  is honest on that population and nowhere else.
  **And it is wrong twice, because the excluded set is not flat.** On genes that
  are not the element's top target, a predicted rise or zero was called
  regulated **33 of 2,079 (1.6%)** and a predicted drop above 0.1 **58 of 67
  (86.6%)** — a fifty-five-fold separation on a population that used to score a
  structural zero. 67 pairs, 34 of them held out, and the row says so.
  **Coverage goes from 2.7% to 63.1% of training pairs and 2.3% to 66.1% held
  out**; on chr21 the answerable (element, gene) questions go from 7,846 to
  347,591. The published `score()` keeps its no-cache default as the control and
  was not re-run.
- **The genome re-banded, and 94.2% of its bands lose their number (2026-09-22,
  lane-band, `9a44b65` registration, `d519bdb` result, 0 requests). Two things
  this session wrote above are corrected here by the lane sent to act on them,
  and both stay visible.** First, **593,765** targets are banded, not 612,323:
  that larger figure is what the sweep NAMES, and 18,558 of them carry no band
  because their gene has no GENCODE v50 TSS. Second, **no swept target is off
  the gate** — `sweep_chromosome` bands the `predicted` and `predicted_coding`
  keys and `crispri.deletion_values` sets `top_target` on a match to either — so
  the ×6.85 is not an error in the shipped table but what that curve quotes for
  a population that became askable only when `ElementResponses` existed. The
  coordinator's guessed split was argued rather than assumed, and it does not
  exist in the sweep. **The real mismatch is inside the gate**: 31.8% of the
  fitted 245 pairs are K562-named against 5.3% of the sweep, and 34.3% carry a
  drop above 0.2 against 5.8%.
  **Re-banded on the registered axis, 559,607 of 593,765 targets (94.2%) lose
  their band** — 34,158 keep 0.9–1, 62 move up and **none moves down** — and all
  four registered predictions were exact, to the count. **Of the 55,957 targets
  published at 0.9–1, 21,861 (39.1%) are not calibrated at all.** A stratum
  keeps a number only with at least 30 pooled pairs and a Wilson interval inside
  one published band; otherwise it reads **"not calibrated here"** with its count
  and its interval, which was written into the registration before anyone knew
  how many strata would need it.
  **The registration's own falsifier fired, and it is recorded as written: the
  axis was wrong.** The track axis bands both its strata where the drop axis
  bands 5.8% of the genome — and on it **44,372 of the 55,957 published 0.9–1
  targets (79.3%) fall two bands**. Either way the conclusion is the same and it
  is the point of the lane: **the claim that 55,957 pairs are at least 90%
  likely does not survive**, a third of that band being uncalibrated or four
  fifths of it two bands lower.
  **On the gate every published band is honest** — 0.25–0.5 measures 0.3514,
  0.5–0.75 measures 0.6629, 0.9–1 measures 0.9754 — and after re-banding the top
  band means **105 of 105**. The off-gate population gets its own table for the
  first time: at drop 0, 4,716 pairs measure **0.0091**; at 0 < drop ≤ 0.1,
  3,783 pairs measure **0.0428** — against the 0.4126 the shipped curve would
  have quoted them. **No web UI hunk is needed**: nothing in `server.py`,
  `index.html` or `cli.py` reads this result or any band label, which was
  verified directly rather than assumed.
  **One thing was measured and deliberately not adopted**, and the reason is the
  lane's own discipline: the two axes are near-disjoint (Jaccard 14.5%) and
  their union bands **128 of 128** pairs at 1.0000. Adopting it now would be
  exactly the post-hoc choice this lane exists to correct, so it goes to the
  next registration instead.
- **Next, open, and the first is now the area's highest-value item: re-band the
  genome from an off-gate curve.** An off-gate confidence is quoted from an
  off-gate curve or it is not quoted. Then the model's **second** targets as a
  shortlist, held against an independent perturbation, since the 86.6% against
  1.6% separation is something no compact one-gene table could ever offer. Then
  the prevalence term, unchanged since 2026-09-17 and now confirmed on 23.7× the
  pairs (×1.23 to ×1.44 on every slicing): re-fit the intercept per screen with
  that screen's own base rate as an offset before any band is quoted.
  *(Written before the lane ran, and its first clause is answered and corrected
  above: the re-banding is done, and there was no off-gate population in the
  sweep to re-band. The list below supersedes this one.)*
- **The union axis, registered and then fitted (2026-09-22, lane-union,
  `ccb06e3`, 0 requests) — the separation confirms, the LEVEL does not, and the
  coordinator's decision below keeps it beside the drop axis rather than in
  front of it.** The registration did the work the situation demanded, because
  the axis had been chosen after seeing one result: it split the 128-of-128 into
  three questions and registered **the level as a description rather than a
  test**, since a hypothesis and its confirmation cannot be the same 128 pairs.
  **The separation is confirmable and was confirmed on populations that took no
  part in the choice**: off the gate, on 1,419 pairs, track-only reads **0.0719
  [0.0596, 0.0865]** against neither at **0.0222 [0.0193, 0.0255]**, ×3.2 with
  the intervals clearing wide; and on GTEx cis-eQTLs, an independent assay,
  0.5939 against 0.5062, ×1.17 and marginal. **The lane's own registered
  prediction there was a null, and it is recorded as wrong.** The search that
  found the axis was priced by permutation at p = 0.0001, which bounds the
  search without choosing the axis.
  **The same test kills the level: the union stratum reads 1.0000 on the gate
  and 0.0895 off it, ×11.** So 0.9–1 is a property of *union ∧ on-gate ∧ K562 ∧
  screen-tested*, not of the axis — the third time on a third axis, after ×6.85
  and ×7. Nothing on this disk can confirm the level, and the one independent
  population readable on the same axis says it is not there.
  **The sharpest objection is the lane's own**: the `track only` cell rests on
  **23 pairs**, and 23 pairs buy no band under the project's own thin-stratum
  rule — the union clears that rule only by pooling them with the 105 the drop
  axis had already banded, so **the 23,468 targets the union adds rest on
  evidence the rule refuses when read alone.**
  Banded: **57,626 of 593,765 (9.70%) at 0.9–1, 536,139 (90.30%) "not calibrated
  here"**, 20,293 up and none down; of the 55,957 published at 0.9–1, 37,333
  keep it and **18,624 (33.3%) are declared uncalibrated**. Coverage sits in the
  drop axis's refusal mode rather than the track axis's aggregation mode, which
  is the right failure to have and a small return for an axis that looked like
  it would band everything.
  **Decision, 2026-09-22: the banding stands as the registered rule computed it,
  and the union does NOT become the shipped default.** Overturning a rule after
  seeing its output is the move every lane today has been asked not to make, so
  the number stands; but two things bar promotion, and both are measured rather
  than felt. The level does not travel (×11), and the genome-wide targets it
  would be quoted for are on-gate yet mostly **neither K562 nor screen-tested** —
  the same inside-the-gate mismatch counted earlier at 31.8% against 5.3%. And
  the increment is unpriced. **So the union is reported beside the drop axis,
  every band carries the population it was fitted on, and promotion waits on
  pricing the 23-pair cell on pairs that can carry it.**
  **Two corrections to the re-banding of a few hours earlier, reported rather
  than quietly fixed**: its windowed off-gate table was training-only, because
  the driver called `add_features` on one arm, so off-gate is 10,226 here
  against 8,549 there; and its track axis returned `another track` for every
  off-gate pair, collapsing to a constant exactly where it needed reading.
  **And a defect the lane names in its own adoption rule**: clause (c) cannot
  distinguish "these sub-populations differ" from "too small to say", and here
  it is mostly the second — reported as a defect rather than claimed as a
  heterogeneity finding.
- **The shared logistic fit was undamped and seven call sites inherited it, and
  no published number moves (2026-09-22, lane-solver, `9abed00` registration,
  `b9a5489` repair, 0 requests).** `crispri.logistic_fit` took the whole Newton
  step, 25 times. Where the classes nearly separate one full step saturates
  every linear predictor, a saturated row on the wrong side then offers a
  curvature of **9.4e-14** against a gradient of order one, and the next solve
  sends the weights to **1e12**. On a random family of ill-conditioned designs
  the old solver returns a point **worse than the zero start it began from in 64
  of 3,738 fits (1.7%)**, worst at |w| = 5.0e12; the repaired one loses none.
  **What makes it a hazard rather than a bug is what divergence preserves: the
  ordering.** On the eight-row design now pinned in the tests the diverged fit
  hands seven of eight rows a probability of exactly 0 — three of them positives
  — and still scores **AUROC 1.0**, the same as the correct fit. The test that
  guarded the function asserted the ordering and a weight's sign, exactly the
  two properties a divergence leaves alone, so it passed on a broken solver. The
  new tests assert what divergence cannot conserve, and two of them fail on the
  old implementation.
  **Four call sites read only AUPRC and AUROC, where this would have been
  invisible; `target_calibration.fit` and its fold are read as probabilities and
  feed four published result files, where a diverged fit reads as "confident and
  wrong" rather than as a broken solver.**
  **Repaired, and nothing moved — checked against the committed copies rather
  than against figures in a document**, because those figures moved today for
  other reasons. Four files byte-identical, two identical but for a date stamp,
  one +550/−0 from another lane's same-day field. The 55,957 top band, the
  CRISPRi verdict and the held-out AUPRC of 0.6909 all reproduce exactly. The
  reason is measured rather than assumed: across all 72 real fits the two
  solvers disagree by at most **3.4e-7** against a largest weight of 27.5. The
  defect was live and inherited; it had simply never fired on this project's
  data. **Left knowingly undone**: the second copy of the line search in
  `target_calibration` is pinned by an equivalence test so it cannot drift, and
  merging it waits for the staging guard's two-day window.
- **The second copy of the median node-open threshold, and why it was not a
  defect (2026-09-22, lane-nodes, `1e57e3e` registration, `c6fc229` result).**
  `candidates.py:398` re-derived, line for line and floor for floor, the split
  that retired `nodes_open` the same day — but **both fields it feeds are
  write-only**: `node_open_element_closed` appears exactly once in the whole
  repository, on the line that writes it, and every score reads the absolute
  `open_on_element` instead. **So no published number rests on a
  between-biosample reading of it**, and the 69 candidates' classes,
  confidences and controls are computed without it.
  **The judgement was registered in advance and it is the useful part**: the
  *count* of a median split is half the nodes whatever the depth, but the
  *membership* is depth-invariant — scale every density in a cell and its median
  scales with it — so comparing memberships across cells is rank-normalised and
  legitimate. It simply can never say a node is shut, and the module's docstring
  had been calling it "openness" beside `open_on_element`, which is absolute.
  That juxtaposition was the whole of the problem.
  **All three registered falsifiers fired negative: 0 disagreements over 260,026
  node calls, 69 of 69 published blocks identical on both fields, and the
  `max(1.0, median)` floor — the one part of the rule that is an absolute
  threshold on a density rather than a rank — binds in 9 of 312 files, all
  chrY, and on none of the 19 chromosomes the blocks sit on.** Reported as the
  negative it is: **a maintenance hazard, not a defect.** The rule now has one
  definition, `node_open_threshold` in `genome/reader.py`, called from both
  sites; the basis string travels with the record; and a test pins that neither
  module can grow the literal back.
- **The prevalence term, registered then fitted — it does not transport, and
  the genome-wide band table is withdrawn as a quotable confidence
  (2026-09-22, lane-prevalence, `1a1a852` registration, `e4914ef` fit, 0
  requests).** The census that had to be written first settles the question on
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
- **The next item, ahead of the list below, which item 1 of that list is now
  answered by the bullet above: price the increment on the pairs that justify
  it.** The union's whole gain over the drop axis is the 23,468-target
  `track only` cell, resting on 23 pairs. Until that is priced the union stays
  beside the drop axis. Then a homogeneity clause that can separate difference
  from under-power, which the one used here cannot and says so.
- **Next, open.** 1. **The union axis, registered properly**: K562 track and
  predicted drop are near-disjoint, either alone gives 100% regulated, and their
  union would band the genome from two single-band populations rather than from
  a drop axis that bands 5.8% or a track axis that aggregates sub-populations
  reading 0.58, 0.60 and 0.87. 2. **The 52 off-gate pairs above a drop of 0.1**,
  carrying 84.6% and 83.3% and refused a band because 26 and 24 pairs cannot
  hold one — the model's second and third targets, and the cheapest way to widen
  that shortlist is more off-gate pairs rather than more features. 3. **No test
  guards the genome-wide totals**, which is how "612,323" travelled into four
  documents. Then the model's **second** targets as a
  shortlist, held against an independent perturbation, since the 86.6% against
  1.6% separation is something no compact one-gene table could ever offer. Then
  the prevalence term, unchanged since 2026-09-17 and now confirmed on 23.7× the
  pairs (×1.23 to ×1.44 on every slicing): re-fit the intercept per screen with
  that screen's own base rate as an offset before any band is quoted.
  *(Two of this list closed within the hour and the bullet above has them: the
  reader landed as `ElementResponses`, `eqtl.py` was re-run through it, and the
  Evidence view now carries both eGene columns. What remains is
  `target_calibration.py`, `motif_transfer.py` and `syntax_tiling.py`, each
  needing its own registration.)*
- **The unknown scoring read off the sweep's whole window, and the naming rate turns out to sit
  below its own chance level (2026-09-27, lane-scoring, `77de3bb` registration, `7f7c8c1` result, 0
  requests).** Three sites in `attribution/unknown_scoring.py` moved onto `ElementResponses`, with every
  pre-existing chr21 figure reproduced to the digit as the control. **The compact table was not
  censoring this module's target**: 5,174 of 12,139 elements name a coding gene from the table and the
  same 5,174 from the window, 0 disagreements, so reading the window cannot name one extra target —
  checked before the run rather than discovered after it.
  **The decisive number is a control this module never had.** On 22,200 length-matched random windows
  inside the same scored span, **unknown blocks name a coding gene at 49.1% against 65.3% at random —
  16.2 points below chance.** The one count the window did raise, blocks carrying a target and a cell
  line (90 → 99), rose *faster* at the control (−16.69 → −16.99 points), and is reported as an inflation
  in those words. Of chr21's 20 constrained-unknown blocks, 1 carries a named target and 0 a cell line.
  **What the lane bought instead.** A magnitude comparison no longer conditioned on the `MIN_EFFECT`
  gate, which fires at 49.5% of unknown-block elements against 68.2% of the rest: −0.0294 at p 0.048
  becomes **−0.0625 at p 0.0005**. A VISTA contrast with **no invented zeros** — `against_vista` had
  been scoring an element covered only by gated elements as exactly 0.0 and averaging it in; removing
  them moves the contrast 0.3306 → 0.3126. And "nothing" on the 69 syntax candidates split into its two
  meanings: 35 blocks hold no registry element at all, 9 were asked and answered below the bar.
  **The registered falsifier fired, and it is worth more than the fix.** Over all 445 covered
  lentiMPRA elements rather than the gated 329, unsigned rho *rose* on HepG2, **0.1699 → 0.2081** (p
  0.002 → 0.0005): the gate discards elements whose sub-threshold effects still track measured
  activity. That is the eQTL lane's `MIN_EFFECT` finding of 2026-09-22 shown against a measurement
  rather than a database. **Next here**: lift `matched_random_windows` into
  `scripts/constrained_unknown_targets.py`, which still carries the borrowed 87% sentence.
- **Owner.** genomeos-i1 since 2026-09-13 (this lane's files; earlier genomeos-f7 and genomeos-c6); the organism-level closure runs
  on genomeos-73's Body runtime and the block evidence on genomeos-fe's
  decoding results.

---

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
- **Code.** `attribution/variation.py` and `genomeos variation --chrom C`
  (this lane). `attribution/syntax.py` and `genomeos syntax --gene G --chrom C`
  (area I's lane, landed 2026-09-12 as dd2b620: per base of one gene,
  constrained across mammals against variable among the imported people,
  rs12913832 named as HERC2's second most constrained variable base). It builds
  on `attribution/bigwig.py` and
  `attribution/constraint.py` (the species axis, done genome-wide),
  `genome/mouse.py` (MGI orthology, node synteny), `attribution/compile.py`
  (the BioLang output), `molecules/graph.py` (the graph the new edges join).
- **Design.** GRAMMAR-BY-COMPARISON.md; GENOME-AS-CODE.md §7 (the mapping
  table with the construct per row); SEQUENCE-GRAMMAR.md §2 and §5 (the
  motif scan, `requires:`).
- **Data.** `variation_chr21`, `variation_chr15`, `variation_vista_genome_wide`
  (Gnocchi per block and element, the two-axis table, controls). To come:
  the other chromosomes, `origin_genome_wide` (stratum per gene, age per library),
  `duplication_chr*`, `motifs_chr*`, the decompiled programs.
- **Requirements.** Every comparison enters as evidence with a confidence
  (D4, D42, D43); within-human constraint is never read as proof of
  neutrality; curated duplication before heuristics (D20's order); the
  species axis and the human axis are reported side by side in the same
  table, never merged into one score; each step is scored on chr21 first and
  run genome-wide second, as the budget was.
- **Verified 2026-09-12.** The within-human track (gnomAD Gnocchi, 1 kb Z
  score) reads through the in-house bigWig reader: three chr21 windows, nine
  requests, 40 KB, 4.5 s, one of the three constrained in humans (max Z 4.3).
  Ensembl's homology endpoint returns 92 mammalian orthologues of OCA2 in one
  call; JASPAR 2026 serves 1,019 CORE vertebrate profiles; UCSC hosts the
  470-mammal phyloP, `genomicSuperDups` and the HPRC 90-way human alignment.
- **Step 1 done (2026-09-12).** `genomeos variation --chrom C`
  (`attribution/variation.py`, `variation_chr21`, `variation_chr15`,
  `variation_vista_genome_wide`): gnomAD Gnocchi (Z per kilobase, 76,156
  genomes; ≥ 2.18 the top decile) read over every UNKNOWN block and every
  scored element beside Zoonomia, each filed in one of four cases (syntax,
  relaxed, recent, tolerant) at confidence 0.3 to 0.6, with the canonical
  coding segments and VISTA as controls. Coding segments read 40% and 31%
  human-constrained (chr21, chr15) against 2% and 17% for the neutral tier,
  so the axis is real and noisy at block resolution; the constrained_unknown
  tier reads 1.7% and 3.2%, no more than neutral, 18 of chr15's 31 measured
  blocks `relaxed`. On the elements the axes add: those constrained on both
  name a gene 75% and 72% of the time, above every other case on both
  chromosomes. VISTA positives 20.2% against negatives 16.9% over 2,223
  elements. rs12913832 is constrained across mammals at the base and its
  kilobase is unscored by Gnocchi. The genome (24 chromosomes, 137 MB of
  ranges): coding 38.1%, introns 18.9%, regulatory tier 11.9%, neutral 2.7%,
  constrained_unknown 2.4%; elements by case name a gene 73.9% (syntax, 825),
  56.0% (relaxed), 50.0% (recent), 45.3% (tolerant), the pair beating either
  axis alone on every count; the constrained_unknown tier is 492 relaxed, 346
  tolerant, 77 syntax and 30 recent of 948 measured; chrY unmeasured, not
  zero. GRAMMAR-BY-COMPARISON.md §6.
- **Step 3 done (2026-09-12).** `genomeos duplications --chrom C`
  (`genome/duplications.py`, `duplication_chr*` for 24 chromosomes,
  `duplication_genome_wide`): UCSC's curated segmental duplications over
  every UNKNOWN block, tier, shared-k-mer pair and scored element. 64,216
  pairs, 5.4% of the genome; the constrained_unknown tier is 4.7% duplicated
  with 216 of 1,098 blocks (19.7%) mostly copies, concentrated on chrY (72%),
  chr21 (61%) and chr22 (42%); neutral 10.3%, regulatory 2.8%; 2.0% of scored
  elements sit inside a duplication; the `similar_to` heuristic fired once
  and was right. A copy reads mammal-constrained (paralogous alignment) and
  unscored by gnomAD (mapping), which is the cause of the two axes
  disagreeing. GRAMMAR-BY-COMPARISON.md §7. Paralogues per gene come from
  step 2's Compara stream (in progress: `knowledge/homology.py`,
  `scripts/origin_genome_wide.py`, `genomeos origin`).
- **Step 4 built (2026-09-12).** `genomeos motifs --chrom C [--gene G] |
  --library L` (`genome/motifs.py`, `motifs_chr*`, `motifs_genome_wide`):
  JASPAR 2026's 1,019 CORE vertebrate profiles scanned over every canonical
  promoter (TSS ± 1 kb) and scored element through a k-mer index of feasible
  motif cores (chr21 in 16 s), held against a dinucleotide-preserving shuffle;
  a gene's `requires:` is its enriched factors with a hit (eight per promoter);
  per library the enriched factors and the factor pairs recurring across
  members beyond the two shares (the operators). Lesson: a single-base
  shuffle made long zinc-finger matrices the only enriched class; the
  dinucleotide control restored MEF2, FOX, PBX, PKNOX1 and MLXIP on real
  promoters. The genome (20,067 promoters): at the factor level the
  libraries recover REST for the nervous system (7.7x), IRF2/3/7/9 for
  immunity and ZNF143 and THAP11 for the housekeeping cores without being
  told; at the pair level the recurring pairs are near-identical matrices
  (MEF2A with MEF2D) and GC-rich profiles co-hitting, the promoter's GC
  content read twice, so the operator table stays inferred at low
  confidence until profiles are clustered into families and the expectation
  is GC-matched. GRAMMAR-BY-COMPARISON.md §8.
- **Step 2 built (2026-09-12).** `genomeos origin [--gene G | --library L]`
  (`knowledge/homology.py`, `scripts/origin_genome_wide.py`,
  `origin_genome_wide`, `origin_presence_genome_wide`): 508 species placed
  once on a 25-clade ladder shared with human (Life to Homo; proxy names for
  the nodes Ensembl's classification omits; strain and assembly suffixes
  stripped for the lookup), the vertebrate and pan-taxonomic Compara dumps
  streamed once (3.9 and 2.6 million rows, 15 s together, species the dumps
  name placed on the way), the deepest clade with an orthologue per gene and
  its paralogues. 19,392 coding genes: 3,774 Life core, 8,250 eukaryote core,
  3,935 animal core, 1,815 vertebrate core, 905 mammal core, 281 amniote,
  200 tetrapod, 121 chordate, 92 primate, 18 ape, 1 human-specific in the
  set; 64,050 `paralogue_of` edges in the knowledge graph (394,068 edges).
  The 42 libraries aged: the oldest are core.splicing, replication,
  translation and timer.telomere (96 to 98% animal-wide or older), the
  youngest systems.ligand_receptor_protocol (55%) and organ_skin (15%
  amniote or younger). Caveats: a vertebrate-dump origin is a lower bound, a
  pan-dump origin can overshoot (family-level calls at the deep strata:
  HOXA1 reads eukaryote-wide); Mus musculus is absent from the human dump
  (its relatives place the grade); a dump with fewer than fifty species
  fails loudly. GRAMMAR-BY-COMPARISON.md §10.
- **Step 5 built (2026-09-12).** `genomeos profiling [--library L]`
  (`knowledge/profiling.py`, `profiling_genome_wide`): each library's
  presence profile over the 355 species, its gain curve along the ladder,
  and pairwise residual correlations as gained-and-lost-together candidates.
  The developmental blueprint libraries share one history (0.97 to 0.99),
  the cores against the vertebrate-born ligand-receptor protocol are the
  least alike (−0.87); shared members inflate pairs and 102 of the species
  are bacteria, so the correlation is `inferred` and a dependency claim
  needs a gene-level profile and a shuffled null. GRAMMAR-BY-COMPARISON.md §11.
- **Step 7 built (2026-09-12).** `genomeos across --locus ZRS | HERC2_OCA2`
  (`knowledge/across.py`, `across_ZRS`, `across_HERC2_OCA2`): the human
  element aligned to mouse, opossum, chicken, frog and zebrafish through
  Compara's pairwise LASTZ alignments, both constraint axes, and the JASPAR
  factors whose motif holds the *same aligned site* in every species (a first
  version counting presence anywhere called 170 factors conserved). The ZRS
  aligns fully to mouse, opossum and chicken (87 to 89%) and 10% to
  zebrafish; 23 of 606 factors hold their site in all four, the HOX, PBX and
  MEIS class with an ETS site, which is the enhancer's known logic, and HAND2
  holds in the three amniotes. The OCA2 enhancer aligns only to mammals,
  constrained among people (Z 3.46) far more than across mammals (7.7%);
  the MiT E-box class and SOX10 hold their sites. Matrix families still
  count many times over. GRAMMAR-BY-COMPARISON.md §12.
- **Refinement: families, composition, sites (2026-09-12).** Matrices are
  counted as curated TFClass families (C2H2 zinc fingers individually) and
  operator expectations are taken within GC-decile by repeat-share strata of
  all 20,067 promoters (`promoter_composition_genome_wide`). Result: no
  family pair recurs beyond its matched expectation in any of the 42
  libraries (84 naive pairs explained by GC or repeats; the last survivor,
  ZNF135 with ZNF460, sat on Alu: 19.4% of those promoters' bases against
  1.1%), while REST stays enriched in systems.nervous at 5.7x. Across
  species the site call was corrected the same day (a factor held at
  different places in different species had been counted as one site): at
  88% identity nearly every window of conserved sequence keeps some matrix,
  so a single locus's site list is uninformative and held-site density does
  not separate four limb enhancers from four matched negatives; the grammar
  question became a positives-against-negatives panel: 39 limb and 33
  neural VISTA enhancers against 76 constraint-matched negatives hold the
  same density of sites (10.2 and 10.7 against 10.6 per kb) and no family
  more often at a false discovery rate of 5% (275 families tested). The
  readout, not enhancer grammar, is what fails; the next readout is in-silico
  mutagenesis against the same panel. GRAMMAR-BY-COMPARISON.md §15.
  GRAMMAR-BY-COMPARISON.md §13.
- **Refinement: library histories against an age-matched null
  (2026-09-12).** Exclusive members, eukaryotic species alone, 200 random
  gene sets per pair matching each side's origin make-up: the developmental
  libraries' shared history (§11) and the cores against the ligand-receptor
  protocol are what their ages predict (z −0.9 to 2.2); 35 of 860 pairs are
  more alike than age predicts, led by immunity and signalling (immune with
  the protocol z 9.3), and seven of the top eight survive removing
  duplicated genes. GRAMMAR-BY-COMPARISON.md §14.
- **Step 6 built (2026-09-12).** `genomeos decompile GENE --chrom C`
  (`genomeos/decompile.py`): every layer read for one gene assembled into a
  BioLang-flavoured view with an evidence note per line and an `unknown { }`
  block naming the layers not yet read; assembly only, no inference.
  GRAMMAR-BY-COMPARISON.md §9.
- **Measured readout (2026-09-13).** `genomeos/knowledge/satmut.py`,
  `scripts/satmut_grammar.py`, result `satmut_grammar`: Kircher et al. 2019's
  saturation-mutagenesis MPRA (44,658 GRCh38 measurements over 21 elements)
  as ground truth for which bases matter. It explains §15's negative: at the
  0.85 match threshold JASPAR sites cover 99.9% of an element's bases, so a
  "held site" was a statement about conserved sequence; at 0.95 they cover
  65% and enrich functional bases 1.11x. Per base the motif score predicts
  function better than conservation (median AUC 0.58 against 0.54) and both
  are weak; together they are better than either alone (conserved bases in a
  strict site are functional 36.5% against 21.1% outside, p 1e-15). No factor
  family survives correction, and the ZRS has 6 strongly functional bases of
  485, so §12's reading is unsupported by measurement. Uses no AlphaGenome
  quota, so it runs beside the deletion chain. GRAMMAR-BY-COMPARISON.md §16.
- **Grammar as a test that could fail, and it failed (2026-09-16).** `attribution/motif_grammar.py`,
  results `motif_grammar` and `motif_grammar_preregistration`, **pre-registered in code and committed
  before the held-out chromosomes were read**. On 51,376 ENCODE4 lentiMPRA elements, three nested
  ridge feature sets — composition; plus strict family-collapsed JASPAR counts at 0.95; plus 330
  arrangement columns (gap bins, strand, helical phase) — with chr8, 9, 21 and 22 held out (6,146
  elements) and a shuffled-grammar falsifier. **The registered claim fails in all three cell lines**:
  arrangement adds +0.0021 [−0.0049, +0.0085] in K562, −0.0056 [−0.0120, +0.0014] in HepG2, +0.0022
  [−0.0050, +0.0097] in WTC11, and the shuffle does as well or better in two of three. **What worked
  is the counts**: +0.234 Spearman in K562 and +0.235 in HepG2 over composition. Four limits recorded
  rather than tuned away, the largest being that an episomal reporter has no nucleosome, which is the
  mechanism that would make helical phase matter. GRAMMAR-BY-COMPARISON.md §17.
  *(2026-09-28, R6: the ENCODE4 lentiMPRA reporter is integrated by lentivirus, not episomal (Agarwal et al.
  2025), so it is chromatinised though outside its native locus; the "no nucleosome" limit is weaker than
  stated. The null stands.)*
- **The count positive tested for transfer, and it survives — the one surviving sequence-to-function
  claim in area J (2026-09-17).** `attribution/motif_transfer.py`, pre-registration committed first.
  On 2,223 VISTA transgenic-mouse verdicts with nine chromosomes held out (643 elements, 351 positive),
  **conservation + counts 0.655 against conservation alone 0.596, +0.060 [+0.022, +0.099]**, on
  elements chosen for conservation; counts beat a dinucleotide shuffle by +0.067; 0 of 1,000 label
  permutations reached the observed value. The control that fails is kept: which *tissue* a positive
  drives is predicted by length and composition alone, and counts add nothing. On the 227 lentiMPRA
  elements inside the constrained-unknown blocks counts add +0.238 in K562 — but swept blind over
  152,556 windows of those blocks the same model **cannot tell real sequence from its own shuffle**,
  because a lentiMPRA element is a selected cCRE and a uniform tiling is not. GRAMMAR-BY-COMPARISON.md §18.
- **Missing.** *(Corrected 2026-09-27 by lane-records: the line that stood here was stale in all
  four of its clauses, and is kept below rather than overwritten.)* It read: "The motif scan and operator
  patterns; phylogenetic profiling between libraries; the decompiled locus view; one developmental
  locus run end to end." All four were built on 2026-09-12 by steps recorded in this area, and three
  then produced negatives: **the motif scan** (`genome/motifs.py`, 1,019 JASPAR CORE profiles over
  all 20,067 canonical promoters, §8); **operator patterns** — searched, and there are none: no
  TFClass family pair recurs beyond its GC-decile × repeat-share expectation in any of the 42
  libraries, the last survivor sitting on Alu at 19.4% of bases against 1.1% (§13), and arrangement
  adds +0.0021, −0.0056, +0.0022 Spearman on 51,376 measured elements, every interval containing zero
  (§17); **phylogenetic profiling** (§11), since tested twice, with 14 of 42 libraries at or below a
  prevalence-matched real-gene null and the fourteen being the ancient core (§19); **the decompiled
  locus view** (`genomeos decompile`, the Blocks-tab card, `/api/decompile`, §9); and **one
  developmental locus end to end**, the ZRS across five species (§12), its site call corrected (§13)
  and carried into a 148-element VISTA panel where held-site density does not separate (§15).
  **What is actually missing**: enhancer pair logic read outside promoters (GATA plus T-box in heart
  has never been read in enhancers); a species-tree null, the only thing that would let a pair be
  called a dependency; and a model that reads the sequence itself, scored on §17's own held-out
  elements — the one readout §17 names as able to overturn its negative.
- **Next.** 1. Done: the human axis over the genome and the case in the
  compiled `region` note, and **the Blocks-tab lane is done 2026-09-17**: every UNKNOWN block
  carries its case, both axes and the confidence, and the lane draws a bar beneath it — amber
  `syntax`, purple `relaxed`, blue `recent`, grey `tolerant`, and nothing where the block was not
  measured on both axes. The case had been computed and committed for weeks while appearing
  nowhere. 2. Done: origin per gene and age per library, and **`orthologue_of` as species counts per clade is done 2026-09-17**: every protein node carries `orthologues_by_clade` from the presence bitmasks, one number per clade rather than one over 355 species — TP53 and APP sit in all 65 Euteleostomi, while OR5H1 is Boreoeutheria-restricted, a distinction a single total cannot make. A clade a gene has no reading in is absent rather than 0. `paralogue_of` edges were already there. Formerly: `orthologue_of` still to add
  as species counts per clade on the graph nodes. 3. `genomicSuperDups` done; paralogues from the Compara stream; `paralogue_of`
  and `orthologue_of` edges in the knowledge graph (genomeos-fe's `molecules/graph.py`, hunk announced first). 4. Done, with families and a GC-by-repeat null: promoters carry family
  enrichments per library and no pair logic; the pair search moves to
  enhancers (GATA plus T-box in heart). 5. Done at library level, and **the gene-level profile with a
  shuffled-library null is done 2026-09-17** — it is a negative about the method. Against a
  prevalence-preserving shuffle all 42 libraries look coherent (median excess +0.307 mean pair
  Jaccard, not one draw above observed). Against a shuffled-library null built from **real genes
  matched on prevalence**, 21 are clear, 7 marginal and **14 at or below it** — and the fourteen are
  the ancient core (translation, energy, transcription, splicing, DNA repair, chromatin), whose genes
  sit in nearly every species and so co-occur with any other such gene. The survivors are patchy and
  lineage-restricted: developmental timing, ECM adhesion, segmentation. The prevalence shuffle
  measures phylogenetic non-independence and calls it co-evolution. Nothing is called a dependency; a
  species-tree null would settle a pair and this project holds none. GRAMMAR-BY-COMPARISON.md §19.
  6. Done: `genomeos decompile GENE --chrom C`, and **the Blocks-tab card is done 2026-09-17**: a
  Decompile button on any gene block reads `/api/decompile` and shows the program beside the map, with
  the layers that carry nothing named above it — a decompiler that hid its empty layers would read as
  completeness, which is the opposite of what the view is for. 7. Done, and tested on a VISTA panel with matched negatives: motif sites held
  across species do not separate enhancers from inactive conserved sequence.
  **The second readout is now done too and agrees (2026-09-17).** AlphaGenome
  in-silico mutagenesis, pre-registered before the instrument existed
  (`attribution/satmut_vista.py`, §20) and run to its full budget: 3,200
  windows of 25 bp over the panel's own arms, each window's predicted effect
  against its Zoonomia phyloP, read as the difference in mean within-element
  Spearman. **+0.0453 at one-sided p 0.1741 — below the registered weak bar of
  +0.05, which is failure and was the registered expectation.** By the median
  the negatives are the closer-tracking arm (+0.0333 against +0.0199). The run
  also censored itself and said so: `MIN_EFFECT` 0.1 is a whole-element floor,
  58.2% of 25 bp windows fell under it and scored exactly 0.0, and 47 elements
  went flat and left the reading (113 of 160 read). The floor cut negatives
  harder (61.4% against 55.1%), which widens the difference rather than
  narrowing it, so the failure is robust; a post-hoc no-floor recomputation off
  the request cache, costing no quota, reads +0.0480 at p 0.1045 on all 160.
  Two readouts, two routes, one answer: this panel does not separate.
  GRAMMAR-BY-COMPARISON.md §21.
- **Owner.** genomeos-bb (this lane's files: `attribution/variation.py`,
  `knowledge/homology.py`, GRAMMAR-BY-COMPARISON.md); the compiled output
  stays with genomeos-f7's `attribution/compile.py` and takes isolated hunks.

---

## 4. Data jobs to be done

Per-chromosome coverage on 2026-09-11 (25 human chromosomes including X, Y
and M). Each row is a job that can run unattended through the Progress tab
or the CLI; results land in `data/results` and are committed.

| Layer | Done | Job | Notes |
|---|---|---|---|
| Sequence, GENCODE rows, ENCODE elements, RepeatMasker, curated UNKNOWN pass, domains (`fetch_<chrom>`) | 25 of 25 (2026-09-11) | `genomeos data fetch --chrom C --analyse` | sequences 1.4 GB in `data/reference`, local only |
| Anatomy inventory | 25 of 25 | done | `anatomy_hg38_by_chromosome` |
| UNKNOWN classification | 25 of 25 with curated repeats; genome-wide 98.5% classified | done | `unknown_genome_wide` summary committed |
| CTCF domains | 24 of 25 (chrM has none) | done | |
| Reader (open nodes per cell type) | 25 of 25, eleven cell types | done | 42.0% (keratinocyte) to 77.1% (hepatocyte) of coding genes read per cell; resumable per cell type, any ENCODE biosample |
| Proteome compiled | 25 of 25 (2026-09-11) | done | 19,478 coding genes: sequence 99.1%, domains 99.0%, function 86.3%, interactions 81.5%, pathways 58.2%, experimental structure 45.2%, predicted structure 98.2%, disease 25.5% |
| Translation verified against UniProt | 25 of 25 (2026-09-11) | done | 19,249 genes: 90.9% canonical identical, 97.8% exact for some isoform; the remaining disagreements are triaged by mechanism, including hg38 frameshift and nonsense alleles detected automatically |
| Knowledge graph | genome-wide (2026-09-11) | done | 41,982 nodes, 330,018 edges after the 2,985 `modifies` edges, 19,283 compiled proteins, largest component 15,165, 3,924 components |
| Modifiable sites (post-translational state) | genome-wide (2026-09-11) | done | 96,362 sites on 13,083 proteins, 352 named writers; `ptm_genome_wide` |
| Enhancer targets, predicted (AlphaGenome) | 24 of 24 chromosomes (4,800 elements) | done | 2,994 elements move a gene, 2,291 name a coding gene; 90.2% of those sit inside the element's CTCF node, against 79.1% for random boundaries on the full archive (qualified 2026-09-14) (79.8% to 91.8% per chromosome) and 71.2% are exactly the nearest TSS; 1,157 behave as silencers. Needs a key; the daily quota ran out mid-run and the job waited and resumed rather than failing |
| All-elements deletion scoring (every ENCODE enhancer inside a node, AlphaGenome) | 24 of 24 (2026-09-16) | done, `scripts/enhancer_targets_all_chain.py` | 961,227 elements, 778,780 requests; 63.7% name a gene (the "87%" once quoted beside it is 52 of 60 random windows in the twelve-locus panel, scored by a different layer, Wilson 75.8–93.1%, and is not a control for this rate); the node beats as many random boundaries by +2.9 points (95% CI +2.03 to +3.81) over 440,377 coding-target elements against a count-matched control, +1.2 to +6.6 across four baselines, and +5.89 on 661 measured CRISPRi pairs (re-audited 2026-09-27). Summaries committed, element tables local and packed. See area I |
| Human panel (HPRC release 1, 90 assemblies, storage units and executor-ready values) | 23 of 24 (2026-09-16); chr2 unread | `scripts/human_panel.py`, one chromosome at a time | 1,786,068 storage units, 174,210 executor-ready (GTEx fine-mapped or MPRA allele pair); seven chromosomes fail a control, six of them on the neutral tier. chr2 (22.6 GB of alignment) needs about 20 GB free to start. See area I |
| HG002 twin | 22 of 22 autosomes | done | 4.05 M PASS variants applied, zero reference mismatches; chrX and chrY are not phased in the GIAB benchmark, so they are out of scope rather than pending |
| Curated repeats distilled | 25 of 25 | done | the 25 RepeatMasker BEDs (60 MB) stay local; summaries committed |
| Composition budget: constraint per UNKNOWN block, a tier per block | 24 of 24 (2026-09-12) | done | Zoonomia phyloP read per base over HTTP ranges, never stored, 2,252 MB for the genome in eight hours; 1,009 Mb tiered, constrained 1.93%, constrained-unknown 32 Mb; see area I |
| Segment parser scored against GENCODE | chr21 and chr22, twelve runs, series closed 2026-09-12 | `genomeos segments --chrom C [--predicted-sites] [--promoter-anchored] [--coding markov] [--rna-filter]` | every result kept side by side; RNA over the exons reaches 92.7% gene precision at 57.9% sensitivity; see area B |

One-off jobs, each needing a decision or a resource:

| Job | Needs | Owner |
|---|---|---|
| Telomere length from HG002 reads (stream a BAM range, no download) | **done 2026-09-12**: 218 MB in 48 range requests over the 601 GB BAM, about 2,676 bp per end at confidence 0.3 | genomeos-8e |
| AlphaGenome live predictions for UNKNOWN regulatory blocks | `ALPHAGENOME_API_KEY` | Albert |
| Measured Hi-C boundaries against the CTCF-only node edges (`genomeos domains --hic GM12878`) | **done 2026-09-12** with the 4DN key: five biosources on chr21, GM12878 supports 58% of the edges against 42% at random | genomeos-8e |
| Cohort expression for the other TCGA studies (`genomeos cancer expression --study`) | 30 s each, pick the studies | genomeos-f7 |
| Retrospective therapeutic benchmark (public tumours with approved targets) | choose the cases | genomeos-f7 |
| Packer 2019 disagreements resolved to one labelling depth | **done 2026-09-11**: 84% on the 164 single-tissue ids; the rest is depth | genomeos-73 |
| Mouse chromosome as the second mammal for node comparison | **done 2026-09-11**: mm10 chr19 and chr11 with MGI's curated orthology, 379 tested nodes at 93 to 94% in the same human neighbourhood | genomeos-fe |
| Sender & Milo per-tissue turnover as maintenance timers in `bio.std` | **done 2026-09-11**: `bio.std.human_turnover` | genomeos-73 |

---

## 5. Next steps, consolidated

**Where things stand at the end of 2026-09-17.** Both genome-wide sweeps are finished: every
chromosome scored by deletion (961,227 elements) and every chromosome read by the human panel
(24 of 24, 187,966 executor-ready units). The executor test and its hold-out are done, `share:`,
`commitment`, `competence` and the integrated reads are in the engine, and the worm's fate rules
read what a cell could know.

**What 2026-09-17 was, in one line: the day the project checked its own numbers and found four of
them wrong.** Not one was found by a failing test — every one preserved the quantity its test
watched.

| what was wrong | how it survived | what it cost to find | where the correction landed |
|---|---|---|---|
| the 882-block tier reading | compared against a tier average, unmatched | withdrawn after two corrections | **withdrawn in place**: ATTRIBUTION.md carries the section headed "the section below overstated in both directions", and +0.284 at p 0.00014 is still legible above +0.100 at p 0.18 |
| the bigWig summariser | preserved every ratio; only absolutes were wrong (×4.62 at element level) | a peer reading the loop by eye | **re-run**, 23 chromosomes against a stop rule fixed first: measured_bp 552.0 → 535.0 Mb, constrained_bp 38.46 → 35.46 Mb, both numbers in the record |
| the body program's germ-layer splits | preserved every layer total; three populations held zero cells, one of them 24.5% of the ectoderm | using `share:` for the first time | **restated in place** with `share:`, the last split absorbing the rounding and saying so in its own evidence line |
| the Blocks lane and `decompile` | poised genes arrive inside `silent_genes` by design | reading what the lane actually drew | **fixed**: poised is a third state the reader emits first, and `decompile` prints the cell instead of dropping it |

**Two of these had reached a conclusion before they were caught, and both conclusions were withdrawn
in the open rather than quietly edited.** That is the claim worth making about this record, and it is
checkable in the files: the wrong number is still printed beside the right one. "Nothing propagated"
would be a claim about luck, and it fails on the first file a reader opens.

Two more that were never wrong, only unsaid: **0.45% of the unknown space has ever been measured by
any assay held here**, and the compiled genome states **940,803 facts of which none is experimental**.
The experiment that would change the first is designed and unrun (284,001 oligos with mappability at
four read lengths). Area J's gene-level profiling closed as a negative about its own method: 14 of 42
libraries do not beat a null made of real genes.

The transferable lesson is in LESSONS.md and it is the day's only general one: **a defect that
preserves the quantity it is checked against is invisible to that check.** Both fixes were the same
shape — give the two quantities two names.

Two lanes run in worktrees under genomeos-79 (items 6 and 7); area G is closed.

**The result of the day is a negative, and it points area I somewhere else.** With the sweep
complete, the 882 blocks of the real unknown were read against every scored element
(`1d0a137`), and the reading was then **twice corrected and finally withdrawn** (`4a2c199`,
`66a45ac`): the per-element rates were taken against the neutral tier unmatched, and the tiers
differ in distance to a promoter by an order of magnitude (81 kb syntax, 240 kb relaxed, 418 kb
neutral, against 42 kb for the genome's elements). Standardised on length, GC and that distance,
**no tier differs from the neutral tier at all**. What survives is cruder and holds: **every
element inside an UNKNOWN block acts less than the genome's scored elements**, from -0.046
(recent) and -0.100 (syntax) to -0.433 (structural). The unknown space is quieter under deletion;
its tiers cannot be ordered against each other by this instrument. The same day, chr22 replicated chr21's ablation as
a negative on both chromosomes, and 44 of the 69 syntax candidates gained nothing from the
entire sweep. What is next, in order of what it decides:

1. ~~The executor claim's external endpoint is inconclusive, and extending it costs a
   pre-registration.~~ **Done 2026-09-17 (`2a07c80`): the extension ran on its 2,260 unspent pairs
   and returned +0.075, one-sided p 0.00037 — the WEAK band by its own bar, not the declared success,
   carried by GM12878 (+0.089 on 1,895 pairs) with Jurkat flat, and with its declared causality
   secondary passing for the first time in the series (+0.198 where DAP-G fine-maps against +0.057
   elsewhere). The ENCODE objection is narrowed, not closed.** E1 was widened genome-wide and run on 2026-09-14 (`cbff5f6`): 1,000
   of 3,260 allele pairs, strongest measured effect first, units agreeing 0.515 against
   matched nulls at 0.501, **+0.014, one-sided p 0.38, upper 95% bound 0.073**. Its own
   pre-registration named that band as deciding nothing, and its declared secondary failed
   (agreement did not rise where DAP-G fine-maps the variant: -0.152 on 57 pairs). So the
   ENCODE objection to E2 and E3 is still open. The 2,260 remaining pairs would tighten the
   bound, but spending them *because* the first 1,000 came out at p 0.38 is optional
   stopping run backwards: it needs a new pre-registration, written first, naming the
   budget, the looks and what each outcome would mean. **Both are now done**: the
   pre-registration is in ATTRIBUTION.md (`1da09a8`), the runner takes the declared schedule
   (`e290744`, with the 1,000-pair skip verified against the rows the first run scored), and the
   run of the remaining 2,260 pairs is under way tonight as the registry job
   `executor_e1_extension`. It is read when it lands, against that document and nothing else.
   Area I.
2. ~~The panel's chr2, then the 24-chromosome rollup.~~ **Done 2026-09-16** (`17a76d8`): chr2
   read with every control passing and moved nothing; eight claims hold on all 23 pooled
   chromosomes, five are named chromosome-specific. The constrained-unknown depletion ends at
   0.829 of the neutral tier on 19 of 23 (p 0.0013).
3. ~~Whether the neutral tier is neutral.~~ **Answered 2026-09-16** (`e9245d9`), and it was
   bigger than the question: *every* non-coding tier reads above the GC- and timing-matched
   background (regulatory 23 of 23, fossil 22 of 23, neutral 21 of 23), so the background is not
   a neutral baseline and carries a systematic +6 to +11 points that belongs to the control. Tier
   statements are read against the neutral tier from now on. What is left is measuring what the
   background actually holds — coding bases, conserved elements, segmental duplication — which is
   a day's work and changes no conclusion already drawn.
3a. ~~The 69 syntax blocks, pre-registered.~~ **Registered (`9b2a4ea`), then cancelled unrun and
   the observation withdrawn (`4a2c199`, `66a45ac`).** Assembling the arms showed they do not
   overlap in the covariate that decides the outcome, and the free version of the question answered
   it: the 40 elements read **-0.100 against the genome's own elements and +0.100 against neutral,
   neither distinguishable from zero**. The registration stands in ATTRIBUTION.md marked unrun, with
   the reason; 2,400 requests were not spent. Nothing is claimed for those blocks.
3b. **Measurement over the real unknown is the binding constraint, and it is now measured and
   designed for.** The coverage map is done (`3124bb4`): counting a block as measured if lentiMPRA,
   VISTA or the CRISPRi benchmark overlaps it by a single base, **the whole unknown space is 1,008.8 Mb
   of which 4.554 Mb has been measured — 0.45%**; the real unknown is 30.6 Mb with 0.52% measured, in
   202 of its 882 blocks; the best-covered tier, regulatory, is at 1.07% and structural at 0.005%.
   Every model built over that space has been extrapolating from half a percent of it, which is what
   the four corrections of 2026-09-16 and 2026-09-17 look like from inside.
   The library that would fix it is written out rather than proposed (`30b2e77`, `a4df8fe`,
   `e08edf9`): **284,001 oligos** after the unorderable ones are dropped — 91,919 tiling the real
   unknown at 300 bp with the untouched blocks first (41,037 of them) and the seven CRISPRi-bridge
   blocks next, 91,918 genomic negatives matched one to one on GC and promoter distance, 8,245
   lentiMPRA actives as in-batch positives, and 91,919 dinucleotide-preserving shuffles. **This row
   said 284,598 until 2026-09-21, which was `e08edf9`'s own output before `97af38c` charged both
   low-complexity rules as a union and `b7c9ff4` added the mappability columns; the rebuild dropped
   597 oligos and the figure was never carried forward here. `data/results/unknown_library.json`
   is the authority and reads 284,001.** The synthesis
   filter was recomputed independently of genomeos-79's lane and the two agree to **0.4%** (41,113
   orderable against 41,297). Two judgements are recorded as design, not measurement: the repeat flag
   is **not** an exclusion (31.3% of the lentiMPRA windows a reporter already read are majority
   interspersed repeat), and the 53 thin blocks are **kept** (median TSS distance 23.6 kb against
   100.9 kb, so dropping them would remove the blocks nearest genes — the covariate that confounded
   everything corrected this week). Manifest under `data/knowledge/library/`, local and git-ignored;
   a design, not an order.
   **Asked on 2026-09-21, because an order is a purchase and 284,001 is one: can it be distilled and
   still meet? It can, and the cut that pays is not the proportional one.** Two arms are sized 1:1 to
   the test arm for a pairing the analysis never consumes — `compare.py` is stratified direct
   standardisation and takes controls per stratum, not one per target — and 202 of the 878 blocks
   tiled are blocks that already hold a measured element, carrying 55% of the test arm. So the
   recommended distillation, **Design A at about 100k**, keeps every base of the 12.31 Mb nothing has
   ever measured tiled end to end (the 676 untouched blocks, 41,037 oligos, plus the bridge blocks)
   and cuts the controls to stratified rather than matched one to one: 30,000 genomic negatives,
   25,000 shuffles, 2,000 positives across activity deciles. It is 35% of the designed library and
   forfeits the already-measured blocks and a per-oligo pairing that only the extreme tail could have
   used. Detectable effect size goes from d 0.009 to **d 0.021 against negatives and 0.022 against
   shuffle**, still a factor of two inside the smallest effect this project's own evidence predicts
   (d 0.05 real-against-shuffle, predicted; d 0.13 real-against-neutral, measured at n=227).
   **The floor is about 12,000 oligos and it is a real floor**: below it the detectable d exceeds
   0.07, so a null would be produced by the sizing rather than by the biology — an experiment that
   cannot fail, which is the one outcome this project has spent the week refusing to buy. A 10k design
   sits on that line, not above it.
   **The arm to keep is the one it would be most tempting to cut.** The shuffle arm is expected to
   come out tied, because §18 reading 3 found the count model cannot separate real sequence from its
   own dinucleotide shuffle over the whole space. But that was a *model score* against a shuffle; the
   library would be the first *measurement* of real unknown sequence against its own composition, and
   a measured tie over 12 Mb is the strongest negative available here — the first evidence that the
   constrained-unexplained sequence holds no reporter activity above its own base composition. Drop
   the arm and that result cannot be bought at any price. The arm whose size is hardest to defend is
   the genomic negatives: the comparison it makes was already run and withdrawn once, because the
   neutral tier is not a control (median TSS distance 418 kb against 42 kb), and the matching leaves
   a 17% residual gap.
   **One decision is still open and it changes which oligos go on the plate rather than how many:
   mappability is carried beside the synthesis filter and folded into nothing.** Over the untouched
   blocks only 45,550 of 91,919 test oligos are unique at k24 (88,282 at k100), and 8,131 oligos the
   repeat proxy passes are unmappable at k24. If the follow-up has to find the locus again, roughly
   half the test arm cannot be attributed. The result file says so itself in `mappability_not_folded`.
   **The comparison itself is now a library function**, `genomeos/compare.py` (`e08edf9`): direct
   standardisation, the control arm carrying the n of the matched targets, the dropped targets counted,
   and both groups' covariate medians in the same result, with nine tests each of a mistake that was
   published this week. Any lane doing a two-group comparison should import it rather than write one.
   **What remains here is not ours to compute:** the library is an experiment, and until one is run
   the honest description of the 98% is that it is unmeasured rather than unknown.
   **Both of the pieces that did remain have landed (2026-09-17), and this row was still calling them
   open until it was audited against the results:**
   - **Which oligos cannot be measured even in principle** — `measurability_real_unknown`. Family A
     (assembly gap, GC extreme, homopolymer, tandem low complexity, non-unique in block) and family B
     (segmental duplication, interspersed repeat) with a stated precedence, and an independent
     reproduction of the peer's reading that matched it exactly: 680 untouched blocks, 13.77 Mb,
     45,900 oligos.
   - **Mappability measured rather than proxied** — `mappability_real_unknown`. Umap multi-read
     mappability (Karimzadeh et al. 2018) in place of the repeat proxy, the two arms reported
     separately, and the four k columns kept **beside** the synthesis filter and never added into it.
   - **The CRISPRi shortlist against the real unknown** — `shortlist_in_real_unknown`, and it is thin,
     which is the point: of the **882** real-unknown blocks, 163 hold any scored element and only
     **7 hold a shortlist element** (coding drop above 0.2, the band where the screens called 105 of
     105). A positive set whose expectation comes from a screen rather than from the model exists, and
     it is seven blocks.
   Area I.
4. **The known-locus benchmark with more loci**, not more readings of the one it has.
   **Stated-interval scoring landed 2026-09-17 and cost one request rather than a
   handful; more loci remains.** The scorer's input is 1 Mb centred on the element, so its reach is 524 kb each
   way, and `read_reach` now checks every published target against it for free: **SOX9 is
   1,450 kb from its element, SHH 979 kb, IRX5 1,164 kb — all three outside the model's input.**
   That retires a standing result: the ZRS's own deletion names LMBR1 and was graded a miss
   against SHH, which was never a candidate. Those loci are now **unaskable** rather than
   pending, and the pending line no longer quotes a price. The headline `target_derived` stays
   **15/17**; a new line beside it holding the two unaskable loci out reads **13/15 (0.867)** —
   slightly *worse*, because both were hits through other layers, which is the direction an
   honest correction goes. H19 was the one locus left askable and was scored: deleting the ICR
   names MIR675 (the miRNA inside H19), H19 itself at rank 3 and IGF2-AS strongest overall, all
   inside the imprinted cluster, and **never MRPL23, the nearest-TSS trap** — yet it scores a
   miss, because the published targets are IGF2 and H19 and what came back was their
   read-through and their antisense. The direction is *represses* where the published mechanism
   raises IGF2; the model carries no allele and no methylation state, which the registration
   said in advance. **Two things are reported rather than changed, both the owner's call:**
   which denominator is the real one, and that the deletion layer takes `predicted_coding`
   first and so can hide a non-coding published target behind a coding read-through (H19 is the
   first locus where that bites; it did not change this verdict). LOCI-BENCHMARK.md §18.
   **A third set of nine landed 2026-09-21 (`ca50b59`, `loci_third`, §21) for three requests, and
   the finding is about the BASELINE rather than the model.** Derived target reads 15/17, 8/9 and
   8/9 across three independently curated sets — a spread of **0.007**. The nearest-gene baseline
   over the same three reads 0.647, 0.333 and 0.778 — a spread of **0.445**. So a benchmark quoting
   a derived rate alone looks reproducible while the thing it is measured against depends entirely
   on who picked the loci, and the baseline's misses here are **silences**: at both long-range loci
   the node named no gene at all, so a baseline that abstains at distance and answers nearby is not
   the baseline its rate describes.
   Reach three frames deep: 3/17 unaskable, then 2/11, now **0/9**, registered as the expectation
   before it ran. The sharpest instance is that **SHH is unaskable through the ZRS at 979 kb and
   askable through SBE2 at 456 kb, where it scores a hit** — same gene, different published element,
   and the difference is the model's input rather than the biology.
   **The direction axis was asked both ways for the first time, and registered undecidable at n=2
   before it ran, so it is two readings and not a rate.** MYC: deleting the PVT1 promoter *raises*
   MYC (+0.139) as published, the opposite sign to the panel's own MYC element on the other side of
   the gene — the layer disagreed with itself across two elements at one target, correctly. HBG1:
   deleting the BCL11A −115 motif *lowers* HBG1 (−1.02), where destroying that motif is one of the
   best-replicated ways to **raise** HbF. Element in the reader's own cell type, target named first
   and correctly, 115 bp away — every input favoured the model and the sign came back backwards.
   The coverage artefact reproduces on this third set, and `input_presence` now says why it can: the
   deletion input is **bought** (9/9 targets, 17/31 controls) and the constraint track is **free**
   (31/31 both arms), so the direction claim's collapse under a coverage stratum is a test while the
   value-on-syntax null under the same stratum is arithmetic.
   **One design consequence, taken 2026-09-21:** `Context.score_region` substituted an overlapping
   annotated element at all three intervals this set paid for, so every reading it labelled "stated"
   was an annotated one. It now takes `substitute=False` and the stated-interval driver passes it, so
   a published coordinate is scored rather than pointed at; the overlap is reported on both branches,
   because an annotation over a published interval makes §10's premise false for that locus.
   **A fourth frame answered this item on 2026-09-21, and it is the hardest negative the benchmark
   has produced: the 0.88 the three curated frames agree on was measured almost entirely where
   proximity and the publication agree, and on twenty loci where they disagree the model reads 1 of
   20.** The frame was drawn by a RULE committed before the source file was read (`2103c5a`), not
   curated: every element in the held-out arm of the ENCODE CRISPR benchmark whose published target
   is not its own nearest coding TSS. The rule reads no model output, no effect size and no
   chromatin score, and it selects exactly the geometry on which the cheap answer is wrong by
   construction — it is adversarial to the benchmark's own headline, which is why it was worth
   drawing. Of 175 held-out elements with a regulated gene, **85 have their own nearest coding TSS
   as the published target**: at half of the field's measured pairs, naming the nearest gene is
   simply correct, and any frame drawn without a geometry rule inherits that. 60 passed the rule,
   the registered cap drew 24, and 4 died at reach (registered 1 to 4, and registered higher than
   the third set's 0 of 9 because a rule selecting elements that skip a nearer gene selects for
   distance). **One request** of a registered budget of 30: the finished sweep had already deleted
   an element inside 19 of the 20 graded intervals.
   **Derived target 6/24 (0.250), 6/20 where the model could answer (0.300)**, against 15/17, 8/9,
   8/9. Above its own chance floor of 0.137, so the layers are not guessing, and nowhere near 0.88.
   **The split under it is the finding**, because the derived rate is a union and the model's own
   layer is not what carries it: deletion **1/20**, GTEx eQTL 3/20, summed window 2/20, the node
   heuristic 0/20, lookup 0/20. What the deletion layer named instead: the published target 1,
   **the nearest coding TSS 14**, another gene 5. The node rule named the nearest coding TSS at 18
   of 20. Those two rows sit beside each other and a long way from the publication: **when
   proximity and the published answer disagree, the model goes with proximity.**
   **The positive control says the reading is not broken.** Drawn by the same rule from the same
   file — of the elements this frame rejected for having the shortcut right, the one with the
   smallest element-to-TSS distance, 1,111 bp — deleting it named HMGA1 first at -0.2271, by the
   deletion layer. The same layer that is right at 1.1 kb names the wrong gene at the twenty loci
   where the right answer is not the nearest one.
   **What it does not say**, and section 22 says so itself: not that the deletion layer is the
   nearest-gene rule in general, since at the curated frames it also gets directions and cell types
   right, which proximity cannot do; and a CRISPRi positive can be indirect, a column this frame
   deliberately did not use. What it does say is that the agreement of the first three frames could
   never have shown this, because none of them contained the geometry that tests it.
   **Two things fall out for free and are the next of this item**: 36 loci pass the same rule
   undrawn, so raising the cap in a commit that says so first takes n from 20 to about 56 at
   near-zero request cost; and the 14 elements dropped for a non-coding target are the first
   published, perturbation-backed non-coding frame this project has had — what section 21 said
   curation could not produce, and what the `predicted_coding`-first defect exposed at H19 still
   needs. LOCI-BENCHMARK.md §22, `188f4ce`.
   **The first of those two was done the same night and the frame is now exhausted: n = 50 graded
   of 60 drawn, and the negative holds in both halves** (`f407915` the raise, `aa1f4be` the run).
   The cap went from 24 to 60, which is every element the rule passes, so it no longer binds and
   nothing is left for a later hand to choose; `assess` and `select` were untouched, and the raise
   was registered before the run with its expected cost, with the interval the remaining 36 would
   have to read for the frame to be called stable, and with what a divergence would mean.
   **3 requests, exactly what `plan` costed in advance, and 4 for the whole frame against a budget
   of 30. Derived target 10/60 (0.167) over every locus the rule returned, 10/50 (0.200) where the
   model could answer, on a chance floor of 0.110** — against 15/17, 8/9, 8/9. Node heuristic 3/50,
   lookup 0/50.
   **The halves diverge, and the registration named the direction before the numbers existed.**
   The cap ran in genome order, so the first 24 are chr1–chr10 and the next 36 chr11–chrX, and the
   later chromosomes are the gene-dense ones, where the floor per locus is *lower* — so a weaker
   second half is the expected direction of a density effect and not a finding. First half 6/20
   (0.300) on a floor of 0.137; second half **4/30 (0.133), below the registered 0.15 to 0.45**, on
   a floor of 0.092. The floor did fall as predicted, **but the margin over it fell too, 2.2× to
   1.4×, so density does not account for the whole drop: the half scored first was the better of
   the two, and that is reported rather than smoothed.**
   **The per-layer split is unchanged in kind and worse in degree**: deletion **3/50**, GTEx eQTL
   3/50, summed window 4/50, node 3/50, lookup 0/50. What the deletion layer named instead over
   the 50: the published target 3, **the nearest coding TSS 26**, another gene 19, nothing 2 —
   against a pure proximity rule's 41 of 50. In the second half it named the nearer gene 12 times
   and a third gene 14, so it is not only proximity; what is stable across both halves is that
   **the published answer is the rarest of the three things the model says.** The control passed
   again (HMGA1 at 1,111 bp, -0.2271), and the coverage artefact reproduces a fourth time at the
   larger n: the target claim reads p 0.006 raw and p 0.188 with coverage held fixed, on 97 of 239
   control windows bought against 50 of 50 loci.
   **Two defects were found in the machinery by running it twice, and both would have cost
   something real.** `score_intervals` writes through `save_result`, which overwrites, and the
   second run asks only for what the first did not cover — so writing its rows alone would have
   **deleted the PCBP1 deletion the first run bought** and regraded that locus as unasked; a
   `merge()` now carries old rows forward by locus, and the result names the one it carried. And a
   reporting bug that read `chrom` off the scored row rather than off `expected` fired after the
   build and before the save, costing 25 minutes of reads; `run` now saves the built rows before
   computing any reading, so a defect in a derived block costs a free `--reaggregate` instead.
   `loci.STATED_INTERVAL_RESULTS` lists `loci_fourth_intervals` in the file now, with a test that
   the runtime patch and the file agree.
   **So this item is closed as a frame**: every element the rule passes is drawn, and taking n past
   60 means registering a different frame rather than raising a cap. What stays free is the 14
   non-coding targets. `f407915`, `aa1f4be`.
   **The non-coding frame was then run and it is two loci rather than fourteen, and it turned the
   H19 caveat into a verdict** (`70c1c10` the registration, `1220095` the result). The 14 elements
   the fourth frame dropped at its first branch are not a non-coding set: that branch compares gene
   SYMBOLS, and the ENCODE benchmark's symbol column is as old as its screens. Joining the
   `measuredGeneEnsemblId` the same file carries to GENCODE shows **twelve are stale symbols**
   (SSFA2 is ITPRID2 and accounts for ten of them, SARS is SARS1, WDR61 is SKIC8), leaving **n = 2**
   genuine lncRNA targets, LINC00885 and CCDC26. That was settled and committed before scoring, so
   it is arithmetic rather than a discovery.
   **Nearest-coding-TSS is not the baseline for a non-coding target**, and the registration said why
   before the run: two of the four derived layers cannot name a non-coding gene at all — `gene_input`
   reads `predicted_coding` only and `node` returns a coding TSS — so they are zero by construction
   and quoting them as beaten baselines would be arithmetic dressed as a result. The registered
   baseline is instead **the model's own any-gene prediction**, the `predicted` field the sweep
   already stored.
   **Scored for 0 requests. Derived 0/2, heuristic 0/2, reach 0/2 as registered — and the registered
   baseline reads 1/2.** At CCDC26 the stored sweep row names the published lncRNA at **-1.691**,
   twelve times the strongest coding effect in the window, and `loci.read_deletion` reports **GSDMC
   at -0.1433** instead — a coding gene 204 kb away that the silencing barely moves — because the
   reader walks `predicted_coding` first and breaks at the first key holding a gene. **So the
   `predicted_coding`-first defect that section 18 found at H19 is now confirmed with a positive:
   the model had the answer and the reader discarded it.** At H19 it did not change the outcome;
   here it is the whole outcome. The one-line fix is named and deliberately not made, because
   `read_deletion` is the shared reader every frame is scored through, and a test pins the current
   behaviour so whoever fixes it is told section 23 must be re-read.
   **Two consequences for other owners, stated and applied nowhere.** Under the corrected symbols
   **8 of the 12 pass every step of the fourth frame's own rule and would have been drawn into it**,
   so that frame's draw is smaller than its rule specifies and "exhausted" holds for the rule as
   implemented rather than as written; whether it is re-drawn, and what that does to 0.200 at n = 50,
   is a decision to take with the reader fix rather than separately. And the general form is worth
   more than either: **any rule matching a published gene symbol against GENCODE symbols can
   silently reject or mis-route a renamed gene**, and the source files already carry the Ensembl ids
   that fix it. LOCI-BENCHMARK.md §23.
   **Both defects were then fixed, registered before anything ran, for 0 requests** (`22cbc74` the
   registration alone, `284b946` the reader, `3374359` the join, `1553abd` the corrected draw,
   §24). The repaired reader ranks `predicted` and `predicted_coding` together and answers with the
   largest effect whatever the biotype, and the argument for it is that **the change is net
   unfavourable by construction**: |predicted| ≥ |predicted_coding| for one element, so a published
   *coding* target can be displaced and never promoted, and four of the five frames are coding-only
   — the fourth frame's rule requires a coding target — so in four of five it can only lower the
   rate. The favourable readings existed and were refused.
   **It cost the benchmark two of its three curated frames.** §9's seventeen hold at 15/17; §19 and
   §21 fall 8/9 → **7/9** each; §22 falls 10/50 → **9/50**; §23 rises 0/2 → **1/2**, which was
   disclosed as known in advance rather than predicted. Every number landed inside its registered
   bracket. The deletion layer alone, which the union hides: 9/17→6/17, 5/9→4/9, 7/9→4/9, 3/50→2/50.
   **0 of 400 negative-window claims moved**, so the matched-window tables of §8, §19, §21 and §22
   are untouched and no p-value changes.
   **The Ensembl join then corrected the draw and the cap binds again: the rule as written passes
   68, not 60**, and 68 is a superset of the 60, so the frame was scored uncapped and nothing §22
   scored was dropped. The corrected frame reads **10/58 (0.172)** graded and 10/68 (0.147) drawn.
   **It did not rise**, so the registered rise-handling never had to be used and §22's 0.200 is
   rewritten nowhere.
   **The sharpest result of the correction goes against this benchmark.** The nearest-coding-TSS
   baseline jumps 3/50 (0.060) to **10/58 (0.172)** and now ties the derived rate — and **7 of its
   10 hits are the seven ITPRID2 elements**, seven perturbations of one gene in one cell line
   recovered together because one symbol went stale. The corrected frame's largest single target is
   7 of 58, so it is less independent than its n suggests, and that was registered before the draw
   was scored.
   **The coordinator's decision of 2026-09-22, on the one thing the lane left open: §22 stays
   interpretable, and the control is restated rather than waived.** Under the repaired reader the
   frame's positive control names `ENSG00000288879` at -1.1735 where it named HMGA1 at -0.2271, and
   §22 registered that a control failure makes nothing in that run interpretable. But
   `ENSG00000288879`'s body **overlaps HMGA1's** in GENCODE, so what failed is the hit rule's
   symbol equality and not the reading the control exists to check: the layer still names the
   published target's locus at 1,111 bp, which is the whole claim the control was carrying. Three
   conditions come with that, because a criterion restated after it fails is exactly the move this
   project has caught itself making before. **One**, it is disclosed as post-hoc and judged against
   a fact that existed before the failure — the GENCODE overlap — and not against whether it
   rescues the run. **Two**, it does not travel: the hit rule stays symbol equality everywhere, so
   **no rate moves**, and the control is reported as failing the letter while passing the claim.
   **Three**, the general question it exposes gets its own registration rather than a ruling here.
   **That question is now the benchmark's most interesting open one**, and the lane measured it
   before anyone asked: of the 21 loci where the repaired reader changed the answer, 19 are
   locatable and **5 promote a gene whose body overlaps the published target** (IGF2-AS over IGF2,
   ENSG00000259006 over MC1R, ENSG00000288879 over HMGA1, ENSG00000240739 over SLC2A3, CCDC26 over
   CCDC26), with 7 more inside 100 kb (HAGLROS 3.2 kb from HOXD, SLC25A3P2 19 kb from ERP29,
   LINC02594 47 kb from SOST, MIR1204 inside PVT1). All are still scored as misses. **So the
   benchmark's "miss" conflates naming a gene somewhere else with naming the right place under
   another name**, and which of the two it is has never been measured. Also left undone, and named:
   `read_gene_input` still reads `predicted_coding` alone, deliberately excluded so two movements
   would not sit behind one number.
   **That question is now measured, classes and thresholds committed before a single miss was
   looked at** (`bf33233` the registration, `d8777fe` the measurement and §25, 0 requests). Of
   **98 strict misses** across the six frames scored through the repaired reader, **15 overlap the
   published target's body, 11 more are within 10 kb and 50 within 100 kb**, leaving 22 further
   away. **So the benchmark's "miss" is overwhelmingly a NEIGHBOUR rather than an alias** — which
   is neither of the two things the question posed. Controls passed before any number was read:
   42 of 42 strict hits classify as exact symbol matches, and the strict rate recomputed here
   reproduces all six published figures.
   **The decisive number was the baseline under the same tolerance.** An overlap-tolerant rate,
   reported beside the strict one as a description, reads 15/50 and 16/58 where the strict reads
   9/50 and 10/58 — and the nearest-coding-TSS rule under the same tolerance rises 3/50 → 7/50 and
   10/58 → **15/58**. **The gap between the model and the rule it is measured against does not
   move.** That column was added after the distribution was seen and is disclosed as late,
   admissible only because it can only make the tolerant reading look worse. §22's proximity
   finding survives and not narrowly: the deletion layer reads 4/50 and 4/58 tolerant against about
   twenty loci where it names the nearer gene.
   **The coordinator's decision, 2026-09-22: the registered case was met and a second rate is NOT
   adopted.** The thresholds cleared — overlap share 0.153 against 0.15, gap 0.120 against 0.05 —
   and the lane recorded that as computed rather than adjusting it, which is the only reason this
   decision can be trusted. Three measured facts decide it against adoption. **It cleared by
   0.003**, and the pool is not six independent frames: the two fourth-frame entries are two
   readings of one draw, and dropping either reads 0.143 or 0.163, so the threshold flips on which
   reading of one frame is counted. **The tolerance buys the model nothing against its baseline.**
   And **only 3 of the 8 distinct overlapping loci are the case the question was about** — a
   non-coding gene lying inside the target (IGF2-AS in IGF2, a pseudogene in SLC2A3, the control's
   lncRNA in HMGA1) — while three are tail overlaps between distinct protein-coding genes, LRRC23
   over ENO2 by 499 bp, which are not aliases in any sense. The class was **not** re-cut after the
   wider test cleared, because narrowing a class once it has passed is the move the registration
   exists to prevent; the composition is reported beside the count and the count stands.
   **Declining a threshold that was met is the mirror of adopting one that was not, so this
   decision is pinned to something checkable rather than to taste: the benchmark reports a
   comparison, and a tolerance that lifts both arms equally changes no comparison.** What would
   reverse it is registered here in advance — **a frame drawn FOR this question**, rather than
   overlaps found in frames drawn for something else, in which overlap tolerance moves the model's
   margin over the nearest-gene rule by at least 0.05. That frame is the next step of this item.
   LOCI-BENCHMARK.md §25.
5. **Albert's language decisions** (BIOLANG-v0.4-ECONOMY.md §10). **No decision blocks code any
   more, and this row said otherwise until 2026-09-21.** Decision 2, amounts or concentrations, was
   the one holding stage 2: it was **priced 2026-09-17** (the registered test could not be run,
   because `volume` is a fraction and a concentration needs an absolute one, so the blocker was a
   language field costing two citations at that scale), then **resolved by Albert on 2026-09-19** and
   **implemented the same day** — a compartment now states an absolute volume beside its fraction,
   and §9's table reads stage 2 as "not started; no longer blocked". Decision 1 is priced, decision 8
   was settled by measurement. Of what remains, decision 9 (refusal or rate) waits on Fukushige &
   Krause's per-stage conversion tables, 4, 5 and 6 are each named measurable with the gate that would
   settle them, and 3 and 7 are preference. **So the live item here is no longer a decision but a
   build: stage 2 — pool, cost, allocation — gated on burden (Ceroni 2015, Frei 2020) and absolute
   abundance against PaxDb, with the burden gate's own stated trap that a shared pool with
   proportional allocation rescales every gene equally and so cannot move a rank correlation.**
   *The sentence above is left as it was written and is wrong in one particular: that trap belongs
   to the ABUNDANCE gate, where §5.3 states it, not to the burden gate, which had no such problem
   and had already passed when this row was written.*
   **Stage 2 is built and gated as of 2026-09-21, and one of its two gates is a negative that was
   derived before the data were fetched.** Four commits after the language surface (`d334b48`) and
   the arithmetic (`eaa603d`): gate (a) `6929a40` and `26b76bc`, gate (b) `e6d4e4f`, `6ebc0b2`,
   `b57de43` and `a9f566f`.
   **Gate (a), burden: passed on three clauses of four.** The unrelated gene's factor falls 1.0,
   1.0, 1.0, 1.0, 0.50, 0.25, 0.05 across seven registered demand levels and returns to 1.0 at
   every one of them when the pools are multiplied by 1e6 — the clause a bug would have passed.
   The order-of-magnitude clause against Ceroni and Frei is **registered not askable**, because
   this project holds neither measured curve, and a test stops it contributing to a pass.
   **Gate (b), absolute abundance: FAILED, on 15,159 genes joining PaxDb's integrated human set to
   GTEx v8 fibroblast median TPM. log10 MAE 0.952005 for the one-to-one baseline and 0.952005 for
   the allocation arm — improvement 0.000000000, bootstrap [0.0, 0.0], against a 0.05 bar.** The
   failure is **algebraic, and it was committed before the transcript arm was fetched** (`b57de43`):
   `Economy._share` returns `capacity / wanted` to every demander of a pool, with no gene index, so
   the allocation arm differs from the baseline by one global constant and a fitted baseline
   absorbs it exactly. The run was done with the pools driven short on purpose, factor 1.754e-4,
   moving every prediction 3.8 decades; the error did not change in the sixth decimal.
   **§5.3's consequence was agreed in advance and stands: the pool layer is optional.** It is
   expressible, it charges a cost, it refuses six modelling errors, it reproduces a burden of the
   documented shape — and it improves no abundance prediction.
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
6. **Four lanes opened by genomeos-79 on 2026-09-16, running in worktrees and merged onto dev one
   piece at a time** (proposed by that session, folded here as the roadmap's rows). They follow from
   the CRISPRi benchmark, `138824f`: on held-out K562, adding the AlphaGenome deletion to activity
   over distance lifts AUPRC **0.550 to 0.633 (+0.083, 95% +0.031 to +0.166)**, as pre-registered in
   the code before the held-out pairs were scored; GM12878 agrees on 14 regulated pairs with an
   interval touching zero. Its denominator was then checked on both arms (`8359872`, after the
   coverage artefact this benchmark and the locus benchmark have each produced): coverage is
   all-or-nothing per element, 3,472 of 3,941 elements fully covered and none partly, and the arms do
   differ — **95.8% of regulated pairs covered against 88.9% of the rest** — but the covered subset is
   not the easier one: distance reads 0.438 AUPRC on all 10,356 valid pairs against 0.441 on the 9,237
   covered, activity over distance 0.517 against 0.519. **So the lift stands with its scope named: a
   claim about cCREs inside a node near a tested gene, not about enhancers at large.**
   - **Area I, measured contact.** 4DN or ENCODE Hi-C contact in place of 1/distance in the CRISPRi
     scorer, K562 and GM12878, rerunning the same pre-registration: does measured contact move the
     activity model the way the deletion did? `genome/hic_contact.py`. Note for it: the boundary
     loader already discards the strength column (`load_scored_boundaries` was added for that), and
     measured boundaries did not separate the one rearrangement the benchmark has.
   - **Area A, a methylation state machine in BioVM.** One CpG with a de novo writer (DNMT3A/B),
     maintenance at division (UHRF1, DNMT1, fidelity below 1) and an eraser (TET1-3), falsified
     against the solo-WCGW loss measured in K562, HepG2 and GM12878. Links to the twin's epigenetic
     age.
   - **Area F, synthetic lethality and selectivity.** A tumour's own losses against DepMap
     dependencies with BRCA1/2 to PARP1 as the control that must be recovered, and marker selectivity
     (A and B and not C) over the healthy-tissue atlas. Target discovery only, inside the non-goals.
   - **Area J, motif spacing and orientation** against lentiMPRA activity, beyond motif counts, so a
     grammar rule is a test that can fail. **Landed `fb2d83d`, and it FAILED as registered**: on
     45,228 elements with chr8, chr9, chr21 and chr22 held out, the grammar's gain over strict motif
     counts is K562 +0.002, HepG2 -0.006, WTC11 +0.002 Spearman, every interval containing zero and
     the shuffled-grammar control doing as well. The features were verified real before the null was
     accepted (83 per element at the median, no empty elements). Its by-product is the largest number
     in the table: composition to strict family counts is **+0.234** in K562 and +0.235 in HepG2.
     **Which factors have sites matters; where they sit does not, at 200 bp on episomal DNA.**
     *(2026-09-28, R6: integrated lentiviral reporters, not episomal DNA — Agarwal et al. 2025.)*

   **All four landed on 2026-09-16, and three of them are negatives.** Rows as their lane reported
   them:
   - **Area A methylation, `d3c5c52`.** The first chromatin mechanism the VM runs: a CpG state machine
     over divisions, with `param x = unknown` now in BioLang so the engine refuses to run until a
     program binds the six unmeasured rates. The solo-WCGW ordering holds on all eight methylomes in
     points (+0.160) and relative terms (+0.454) but thinly with region held fixed (margin 0.0075),
     and the coordinating session added the caveat the lane had not: the ordering was already known
     here, so this is consistency with a fact in hand rather than a prediction. Its own negative is in
     the section: the dense class erodes too, which context dependence alone cannot produce.
   - **Area I measured contact, `5a31c39`. Read the `reading` field, not the `verdict`.** Balanced 4DN
     Hi-C at 5 kb (K562, GM12878) in place of 1/distance **costs 0.081 AUPRC on 9,165 training pairs**
     (interval -0.118 to -0.011, clear of zero); the held-out gain is +0.012 with the interval
     straddling zero, so the pre-registered rule reads "passed" while the measurement settles nothing,
     and the result file now says so in a computed `reading` field beside the verdict. The coverage
     warning was acted on: the matrix reaches 451 of 451 regulated training pairs and 99.2% of the
     rest, both baselines unmoved, missing contacts counted by named reason and never as zeros. Why it
     fails: at CRISPRi distances raw contact is mostly the distance decay (0.379 against 0.442) and
     observed-over-expected alone collapses to 0.083.
   - **Area F synthetic lethality and selectivity, `2305bee`.** DepMap 24Q4, 19,749 pairs, **44 hits at
     q <= 0.10** with a permuted-label null finding nothing in 20 of 20 permutations, and five of seven
     declared controls recovered (ARID1A-ARID1B, SMARCA4-SMARCA2, MTAP-PRMT5, MTAP-MAT2A, WRN in
     MSI-high). The two PARP1 rows are significant, negligible (0.047 gene-effect units) and fall to
     nominal under the sweep's own correction, as the pre-registration said to expect of a knockout
     screen. **Marker selectivity FAILS its claim**: approved targets rank below random gates (median
     percentile 25.1 against 51.5), and because the placebo null is not flat (delta -0.496) the
     benchmark read against its null sits at the 87th percentile, short of its declared bar. The flaw
     is located: a 7-TPM threshold crossing mistaken for surface density.
7. **Two lanes opened 2026-09-17 by genomeos-79, both following from the day's findings rather than
   opening new ground.** Each reports a pre-registration before scoring.
   - **Area I, calibration. Landed `bdc2364`, and its pre-registered reliability claim FAILED — on the
     population, not the curve.** Out of sample on training chromosomes the curve is reliable (9 of 10
     bins, ECE 0.006); on held-out pairs 6 of 10 bins fall inside their Wilson intervals against a
     declared 7, and **all four failures are under-predictions in the same direction**. The cause is
     named: the held-out screens call **6.53% of pairs regulated against 4.97%** in the screens it was
     fitted on, 1.31 times the prevalence, and a single log-odds shift of +0.679 restores 9 of 10 bins
     — fitted on the held-out labels and reported as a description of the failure, not a test. What
     does hold is the table with no fitted model between it and the measurement, on 286 pooled K562
     pairs: at a predicted drop **above 0.2 the screens called 105 of 105 regulated**. Genome-wide only
     **24,114 coding targets (5.5%)** reach that drop, 3.2% (chr18) to 10.1% (chr19) per chromosome —
     a shortlist to hold against an independent perturbation, not a claim. Scope: K562, conditional on
     a screen having tested the pair, with four named falsifiers. Turn the CRISPRi benchmark into a
     calibration curve, so every predicted target in the finished sweep carries a confidence backed by
     measured precision rather than by the model's own effect size. `attribution/calibration.py`, `scripts/calibrate_targets.py`. It is the
     direct consequence of the structure-against-perturbation lesson: if perturbation is the
     discriminating variable, its output is what deserves calibrating. **The denominator to watch is
     the one the locus benchmark and the CRISPRi lane have each already been caught by** — a
     calibration curve conditioned on "was this element scored" inherits the coverage of both arms.
   - **Area J, transfer. Landed `b7e1405`: the motif-count positive transfers off the reporter, and
     the unknown space's binding constraint turns out to be coverage.** On VISTA's in-vivo assay with
     nine chromosomes held out, strict family counts add **+0.060 AUROC over conservation** (interval
     clear of zero, dinucleotide-shuffled null at 0.562, 0 of 1,000 label permutations reaching the
     observed value): conservation plus counts 0.655 against conservation alone 0.596 — weak, on an
     unmatched comparison, and the lane says so. Two negatives beside it: the tissue-group gains (heart
     0.846, neural 0.763) **vanish against a length-plus-composition control** (0.839, 0.744; counts
     add +0.006 and +0.019 with intervals containing zero), and swept blind over 30.6 Mb of whole
     blocks the model **cannot separate real sequence from its own dinucleotide shuffle** — a model
     fitted on selected candidate cCREs is not a genome-wide activity track. **The number that sets
     area I's next question: 161 of the 882 blocks hold a measured element and 721 hold none.** What
     the real unknown needs is an assay over it, not another model.
     The lane's original brief: whether the +0.234 Spearman on lentiMPRA transfers to VISTA and to the
     882 blocks.
     `attribution/motif_transfer.py`. **Two traps are already known here**: VISTA as an endpoint is
     weak by construction in this project (its negatives are not matched to its positives), and the
     882 blocks were read with the model's own deletion predictions as the outcome, so a motif reading
     that uses the same elements is not independent evidence about the same blocks. If it disagrees
     with the 0.248-against-0.293 reading, that disagreement is the result.
8. **Disk. Re-measured 2026-09-17, and both halves of this row were stale.** The Chrome
   code-sign clones are **gone**: 63 of those directories remain in the system temp tree and
   every one is 0 bytes, with the whole temp tree at 1.0 GB. Nothing here is Albert's call any
   more. The project is **12 GB, not 7.3** — `data/knowledge` 6.2 GB, `data/reference` 1.4 GB,
   `.git` 129 MB — so it nearly doubled while this row said otherwise. The disk is at **92%,
   36 GB free of 460 GB**, and the weight is outside the repository: `~/Library` 165 GB,
   `~/Documents` 50 GB, `~/SITES DEV` 32 GB. Inside the repository the recoverable item is
   `.claude/worktrees` at 2.8 GB across 8 agent worktrees, 4 of them locked, clearable when
   their lanes are idle. What remains a real question is whether distilling `data/knowledge`
   is worth it, and that is inside the project rather than Albert's.
9. **A coordinating session, from 2026-09-21, and the drift it exists to stop.** Albert's
   direction: the work of the last days was landing but not against a published plan, so what
   was going on could not be read off the Progress tab. One session (named in §8) now holds
   this document, assigns the lanes below, and keeps the board honest; every other session
   builds and sends finished row text rather than editing the roadmap itself.
   The drift was measurable rather than felt. On the day this row was written the work board
   carried **31 entries of which 2 were live**: eleven said "working" with their last update
   between five and eight days old, and the sessions that wrote them are gone. A reader opening
   the Progress tab saw eleven lanes in flight and two were. The board's own staleness label
   was doing its job and the tab was still misleading, because a stale entry from a dead
   session and a stale entry from a session about to come back look the same.
   **The open lanes, one owner each, and nothing owned by a session that no longer exists:**
   - **Area F, the last clause of milestone 1.2** — **done the same day (`e21d34a`)**, and the
     lane changed hands once on the way: it was assigned to genomeos-e6, whose session ended
     before it started, and ran as lane-cancerF. The pin is 0, the benchmark is 7/7 on all three
     readings, its seventh case is reached by an expression call rather than a DNA event, and it
     turned up a third defect of the same class in EML4-ALK. Area F's row has the detail.
     **The first test of this row's own rule, and it held**: the board said the lane was gone,
     so the work was reassigned in minutes rather than sitting under a name nobody was behind.
   - **Item 4, more loci** — **done the same day as lane-loci4 (`2103c5a`, `188f4ce`)**, with
     Albert's approval to spend AlphaGenome requests on it, and it cost 1 of a registered 30. It
     is a negative and it is the item's answer rather than another reading of it: a frame drawn
     by a rule committed before the file was read, on the geometry the first three frames did not
     contain, takes the derived rate from 0.88 to 0.300 and the model's own layer to 1 of 20. The
     item's row has it. **Two follow-ons need no new registration**: the 36 undrawn loci that pass
     the same rule, and the 14 non-coding targets it set aside.
   - **Item 5, stage 2 of the economy** — genomeos-0e, last seen seven hours before this row.
     Pool, cost and allocation are in the language and the arithmetic; the burden and abundance
     gates are what remain.
   - **Item 3b, the assay over the real unknown** — designed, unordered, and not a lane anyone
     can open without Albert: 284,598 oligos is a purchase, not a commit.
   The rule this row establishes, so the drift does not recur: **a lane exists on the board or
   it does not exist.** An entry whose session is gone is retired rather than left to read as
   work in flight, and the Progress tab separates what is running from what was run.

   **The lane table above is the state at the row's first hour and is left as written. This is
   where it stood at the end of 2026-09-22**, thirteen lanes later, every one of them registered
   before it measured and every one reporting its own negatives:
   - **Closed with a result**: area F's two mechanism defects and the CNA/SV/expression cases
     (milestone 1.2 at **9/9**); stage 2 of the economy with gate (b) a **derived** negative; the
     fourth locus frame at n=50, its reader defect, the non-coding pair, and the re-read of all
     five frames; area B's stricter site call, its twelfth and thirteenth biosamples, and the
     depth banding of the whole reader family; the calibration's `top_target` gate, the genome
     re-banding, the union axis and the prevalence term.
   - **Running**: `lane-solver` (the shared logistic fit is undamped and every caller inherits
     it), `lane-worm2` (area E's exposure read against a lineage mean, and the 332-to-414 move).
   - **Decisions taken here rather than deferred**: §22 stays interpretable on a restated control;
     no overlap-tolerant rate, because the tolerance lifts the baseline equally; the union axis
     banded as registered but **not** promoted; and the **593,765-target band table retired**
     rather than corrected.
   - **Still Albert's, and none of it blocks a lane**: mappability, whether the oligo library is
     ordered and at what size, whether a surface score should read the alteration, a permission
     rule for `scripts/commit_own.sh` so a registration can be committed on its own, and disk at
     **95%**.
   **What the day actually established, across four areas and one shared theme:** a number is
   trusted for the population it was measured on and no further. It showed up as a curve quoted
   ×6.85 off its own gate, a band table that is one screen's to within 0.37%, a reading predicted
   by assay depth to z = −0.11, a field that returns half the nodes by construction, and a
   direction that passes at 41 of 44 on downward effects and fails at chance upward. **Three
   guards were added so the recurrences are caught by a machine rather than by a reader**:
   `genomeos work retire` for lanes nobody is behind, a message check in `commit_own.sh` for the
   shared scratchpad, and `tests/test_band_totals.py` for the genome-wide totals.

10. **The execution plan of 2026-09-27: one picture, waves ordered so nothing is built twice, and a
   capacity model measured rather than assumed.** Written after five independent assessors read every
   area against its own tests and result files. This item is the plan the lanes execute; the rows
   above it are the record.

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

   **There is no breakthrough, and the assessment says so.** There is one result worth showing an
   outside reader — the AlphaGenome deletion feature lifting held-out CRISPRi AUPRC 0.550 → 0.633
   (+0.083, 95% +0.031 to +0.166) on 1,744 K562 pairs, frozen before the held-out pairs were scored —
   and it needs a second cell type with real n and one published baseline figure beside it. The
   unusual asset is the discipline, not the biology: registration that binds a single author by
   commit order. Roughly two of sixteen documented withdrawals bound someone else's claim; the rest
   bound this project's own. Moving that ratio is the strategic goal of the next stretch.
   *(2026-09-27, later: the CRISPRi result now stands beside the published figures on the same
   pairs — held-out 0.691 after the cache repair, "in the range of ENCODE-rE2G" on a DNase-only
   refit, still one cell type; HCT116 is the second, at 705 requests. Area I has the detail.)*

   **Two decisions were taken on 2026-09-27 rather than deferred.** *The oligo library is not
   ordered*: permanent non-goal D8 forbids wet-lab synthesis, and 1.3's own exit criterion asks for "a
   designed experiment for the rest, which is built" — it is built, so that clause is met, and the
   plan loses its only months-long dependency. What blocks 1.3 is clause 2, which was never split the
   way clause 3 was on 2026-09-17; it is split here the same way — attributed where measurement
   exists, a labelled *lead* elsewhere. *BioLang stays a generated package*: the export already runs
   36 files at 8/8 checks with the application unimportable, nothing outside this repository consumes
   it yet, and 2.0's clause "GenomeOS depends on it" is a means rather than an end; it is restated as
   "the engine is separately installable and independently tested", and a real split waits for a
   real consumer.

   **Capacity, measured on 2026-09-27.** Apple M3 Pro, 11 cores (5 performance), 18 GiB RAM with about
   6.5 GiB reclaimable under load, data volume at 95% with 22 GiB free. Twenty-one lanes ran on
   2026-09-22 at a mean of about **35 minutes each** and **2–3 concurrent** without a file collision.
   The binding limits are, in order: RAM (the per-element cache reader peaks near 2.9 GB on chr1),
   the shared checkout's index and shared documents, the single AlphaGenome quota, disk for anything
   that fetches, and the coordinator's review of each report against its result file. So:
   **at most four lanes at once, of which at most two read the per-element cache, at most one spends
   AlphaGenome quota, and at most one fetches more than 1 GB.** Only the coordinator runs the full
   `scripts/check.sh` (six minutes, CPU-bound); lanes run their own targeted tests. **The machine must
   be on mains power** for any run longer than an hour: it sleeps on battery or with the lid closed.

   **The waves, ordered so each one's result is an input and never a re-do.**

   | wave | lane | why here | lane-sessions |
   |---|---|---|---|
   | 0 *(running)* | the node-containment claim: audit its baseline and interval, then hold it against 661 measured CRISPRi pairs | quoted seven times with no interval and an undefined baseline; everything node-based waits on it | 2 |
   | 0 *(running)* | `unknown_scoring` off the one-gene table | it is behind 1.3's clause 2 and was on no list | 2 |
   | 0 *(running)* | records truth pass | stops lanes reading rows that are wrong | 1 |
   | 0 *(running)* | repository size against GitHub limits | the daily merge to `main` needs the headroom known | 0.5 |
   | 1 | make the deletion result publishable: a second held-out cell type, and the published ENCODE-rE2G / ABC figure on the same pair set | the one external contribution in reach, no wet lab | 3 |
   | 1 | the reader family in **one** lane: `nodes_open` given a between-biosample definition, `genes_poised` against H3K27me3, `genes_read` normalised | one recomputation of 264 cached blocks instead of three; closes 1.1 | 2.5 |
   | 1 | the therapeutic evidence tier (observed alteration > this patient's RNA > a STRING neighbour) | `BURIED_SURFACE_TARGETS` 2 → 0; must be a tier, because CD19 has no origin record and a no-origin penalty would demote the expression route | 2 |
   | 1 | 2.0 part one: the 27 engine-only test files (162 tests) travelling with `biolang`, and a version string of its own | the only missing clause of 2.0's proof | 2.5 |
   | 2 | 2.0 part two: grammar docs copied in, `bio repl` in the isolated check, a default for `import protein:` | release mechanics after the tests exist | 2.5 |
   | 2 | the remaining one-target consumers under **one** registration with a falsifier per module | five modules, one shared reader; one registration instead of five | 3 |
   | 2 | D: verify NIST/GIAB Q100 T2T-HG002 (trio-binned, so its haplotype is the parent), then phase the 1,430 candidates | a phased resource this project had never named; probe before spending | 0.5 + 2 |
   | 2 | C: reachability of the Ochoa 2020 phosphosite table and CPTAC through the cBioPortal client already shipped | turns "can be modified" into "observed modified in N experiments", never "occupancy" | 0.5 + 2 |
   | 3 | E2/E3 widened past chr21 and chr22 — 184,604 fine-mapped units on disk, 161 read | the largest unread population in area I; spends quota, so it runs alone on the key | 3 |
   | 3 | A: the three facts with evidence kind `none`; v0.4 stage 4, partitioning division | small correctness, then the next stage | 0.5 + 2 |
   | 3 | ↳ **A done 2026-09-27 (lane-parser, `61c5838` census, `b41ddd6` parser).** A key a block does not read is now an error naming the key, the block, the closest match and the accepted keys; the census before it found 4 hits in 428 program texts, all `rate:` on rules in the engine's own smoke and REPL programs, so the inhibition there was never applied. Two bugs behind it: a nested transcript never inherited its gene's evidence, and NKX2-5-201's loci matched no GENCODE record (now GENCODE v50 ENST00000329198.5, curated). `none` facts 3 → 2, each remaining one with its reason stated. Open: `genomeos evidence` drops facts silently when a program fails to parse (26,845 → 26,824, no error); `engine_package.json` regenerated 12/12 but uncommitted, its removal check needs Albert. v0.4 stage 4 waits behind item 11 R1–R5 | done | — |
   | 3 | ↳ **Parse failures reported 2026-09-30 (lane-evparse, night shift; `91b085c`)** | `evidence.py`, `web/server.py` (pass-through), CLI | No program fails to parse today (26,851 facts in 42 programs; 986,947 with the compiled chromosomes); `c22c8c6` (2026-09-28) had already named failures on stderr and exited 1. Added: `parse_errors` (path and message), `parse_error_count` and `complete: false` whenever any program fails, under every filter and on the csv branch; the CLI prints "incomplete: N program(s) failed to parse". Checked on a scratch copy with one broken program: the count falls by exactly 21 and is reported incomplete; every other program's rows are byte-identical; 8 JSON payloads and 4 CLI runs byte-identical on the real repository. Deviations: a clean run carries no `complete: true` (the accepted evidence-view test pins the whole `/api/evidence` payload's digest); the Evidence tab does not yet show the report | done | — |
   | 3 | ↳ **Completion 2026-09-30 (lane-evparse; `854c2a4`, additions only)** | `evidence.py`, `index.html`, tests | After the review: a clean read now states `complete: true`, `parse_errors: []` and `parse_error_count: 0` under every filter and on the csv branch; the Evidence tab shows "Incomplete: N program(s) failed to parse" with each file and error above the rows whenever a read is incomplete, and the CSV status line says so. The evidence-view contract test's digests were redefined on the line after the held ones, with a dated reason; each new payload with the three report keys removed hashes to its old digest. 73 evidence and web tests pass. **Remaining:** the superseded "clean read adds nothing" comment and two tests kept as strict expected failures stay until their lines clear the guard (from `91b085c`, 2026-10-01 22:40), then removed | done (clean-up listed) | — |
   | 3 | ↳ **v0.4 stage 4 done 2026-09-28 (lane-v04, `4ff3cd4` registered, `73af2d5` built)** | `lang/grammar.py`, `ir/model.py`, `runtime/division.py` | One new key, `partition: duplicate | contents | binomial` on an `event`, with the checkpoint expressed as the event's own `when` through `ir.model.matches` rather than a private mechanism. `runtime/division.py` conserves every species exactly, splits binomially below the regime threshold under `units: copies` and by exact halves above, and refuses half a molecule. The default stays `duplicate`, so no committed program changes. **The registered negative is the finding: the spec's own gate for this stage is refuted, not passed, and the refutation was derived from the arithmetic and committed before the instrument existed.** §8 asked that dilution separate a stable protein from a short-lived one; partitioning is gene-blind, so it contributes **0.0** of the 9.392136 decades between a 240 h and a 6 h protein (identical whether contents are partitioned or duplicated), and in a steady state dividing **compresses** the separation 8.80-fold, 40.0000 → 4.5455. What passes is the stage's own claim: a protein given no half-life falls by exactly 2^−8 partitioned and not at all duplicated; 2,000 of 2,000 binomial splits of 20 copies conserve exactly (mean 9.8895 against 10, variance 4.8993 against 5); variance is exactly 0 above the threshold; an unmet checkpoint fires 0 divisions. **8 of 8 registered clauses met**, 24 tests added, check.sh 221/221, engine package 12/12 (36 → 38 files, 177 → 192 packaged test functions). **Reported and excluded from the verdict** because the registration did not name it: the demo declares a division cost of 1e10 ATP against a 3e9 pool, so the resource checkpoint refuses every division — a pool is a standing stock and a cost is a draw over a cycle, and nothing converts between them until stage 3 supplies ATP over time. It cannot claim any measured half-life (none is held), nor that real partitioning is binomial (no partitioning-error data held) | the engine result file is the coordinator's to rebuild | done |
   | 3 | H: budget constraints in `design`; network knockouts in `experiment` | the two concrete items left in the weakest area | 2 |
   | 3 | ↳ **H done 2026-09-28 (lane-h, `faf75c4` registered, `6c3c66d` budget, `8e705f7` network knockouts, `66675d8` pin fix).** A `design` run under a `Budget` never answers over the limit (each knockout, addition or knob moved costs 1 unless priced; `at_most` never counted knob moves); over-limit candidates run once so the result names what the limit excluded; no budget is byte for byte the old behaviour. two_intestinal_founders answers -POP-1 at limit 1 and the wild type at limit 0, with -POP-1 named as excluded at loss 0. A network knockout holds nodes at zero or blinds one edge's target, and ranks the mean change over untouched nodes against matched random knockouts: CycD off turns the Fauré 2006 cycle into the G1 fixed point (Rb, p27, Cdh1 on), consequence 0.556, p = 0.1 over all nine matched knockouts. 20 tests; no falsifier fired. **Negatives:** the rank is coarse on small networks (0.1 is the floor at nine draws; one ring edge kills oscillation at p = 1.0 because every ring edge does); a hard limit hides the fix in the excluded list rather than pricing it; and **`NetworkRuntime.run(clamp=)` does not hold a species at zero inside the Runge–Kutta stages** (a clamped `a.mRNA` left protein A at 1.87 against 19.2, a knockdown not a knockout), so node knockouts zero what makes a species instead; the clamp is its own lane. Open: no CLI or grammar surface for network knockouts | done | — |
   | 3 | ↳ **Clamp diagnosed 2026-09-30 (lane-clamp, night shift; `fbb5840` test only)** | `runtime` (no change) | The open item above was already closed: `eb65f5d` (census) and `51c10b8` (fix) on 2026-09-28 hold a clamped species at its clamp inside every Runge–Kutta stage. Diagnosis, measured: the old 1.87 was an integrator leak, not slow decay (steady from 30 h to 60 h where pure decay gives 1e-129; scaling with the step, 0.48, 0.96, 1.87, 2.97 at dt 0.005 to 0.04; matching k_tl·r·dt/(2k) on a one-gene network). What `clamp=` promises, from its docstring and callers: it fixes one state variable (a signal such as Nodal, or an mRNA at zero); it does not remove protein, and no knockout operation was added. `fbb5840` tests the protein of an mRNA clamped at zero against the independently calculated P(t) = P0·e^(−kt), k = ln 2 / the declared 0.139 h half-life, at dt 0.005, 0.01 and 0.02, on the one-gene network and the ring: agreement within RK4's own error (4.5e-6 relative at 1 h); the pre-fix code fails it at the first step. No result or pin moves | done | — |
   | 3 | ↳ **Clamp fixed 2026-09-28 (lane-clamp, `eb65f5d` census, `51c10b8` fix).** `NetworkRuntime.run(clamp=)` now holds a clamped species at its value inside every Runge–Kutta stage (derivative zero, clamp written into each stage state); API unchanged, unclamped runs identical. The ring's protein A with its mRNA clamped at zero reads 1.1e-66 (was 1.87, against 19.2), so the row above's "knockdown" is historical. Census: two callers. **Gastrulation moved, away from its expectations**: default mesoderm 0.108 → 0.100 and endoderm 0.417 → 0.425 (stated 0.35 and 0.20); layer order holds; **its test's mesoderm check passes only by float rounding** (0.35 − 0.10 = 0.2499999… against a bar of 0.25), before and after the fix, so that pass is not evidence. The cell view cannot move. Seen for R3 in the same function: a missing regulator reads as 0 and a gene with no stated maximum gets max rate 0 | done | — |
   | 3 | ↳ **Gastrulation measured honestly 2026-09-28 (lane-gastrula, `a398e9c` test, `19423bd` labels).** The layer order is real and tested. The proportions miss: mesoderm 0.10 against 0.35, endoderm 0.433 against 0.20, ectoderm 0.467 against 0.45; the test now resolves a float tie against the model and marks mesoderm `xfail(strict=True)`. **The expected proportions cite no source** ("inferred, order of magnitude" since the first commit) and are labelled unsourced; the thresholds and three repression rules cite nothing either. Cause: each fate owns a fixed band of NODAL (SOX17 above about 1.8, TBXT only 1.3–1.8), set by uncited repression strengths, not by the gradient. Two free gradient numbers can fit any three proportions (one computed point gave 0.467 / 0.317 / 0.217), so tuning them was refused and the defaults stand. Open: a sourced germ-layer census for the targets and a sourced Tbxt–Sox17 repression before any retune | done | — |
   | 3 | ↳ **Gastrulation against a measured census 2026-09-28 (lane-census, `3c4e301` registered, `d55cb19` comparison).** The one measured human gastrula, Tyser et al. 2021 CS7 (E-MTAB-9388, 1,195 Smart-seq2 cells, one embryo, author labels), with the mouse atlas (Pijuan-Sala et al. 2019, 116,312 cells) as a no-verdict secondary. Registered before the run: two label mappings × two site filters, interval = their spread, tolerance 0.05; the registration says the lane had seen both the model's shares and the counts while choosing mappings, and predicted a mesoderm miss. **One default run: falsified on all three layers** — model 0.475 / 0.100 / 0.425 against CS7 ectoderm 0.046–0.172, mesoderm 0.694–0.773, endoderm 0.116–0.213. Bounded: a sampled census of one embryo is not a 1-D axis, and the model names no stage. **A rule has the wrong sign:** Lolas et al. 2014 show Brachyury *activates* Sox17; the module's inferred `Tbxt inhibits SOX17` contradicts it (Sox17 inhibiting TBXT is supported). No measured NODAL threshold maps to model units; zebrafish work (Dubrulle et al. 2015) argues induction kinetics, not a clamped level, set fates | done; the sign fix is its own lane | — |
   | 3 | ↳ **Tbxt–Sox17 sign corrected 2026-09-28 (lane-sign, `b1f3405` registered, `5cbce26` result).** `Tbxt activates SOX17` per Lolas et al. 2014 (Fig. 1A ChIP-seq, Fig. 3A knockdown); `Sox17 inhibits TBXT` now cites Fig. 3C; strengths and thresholds unchanged and labelled unsourced. **The registered expectation was half wrong:** default shares 0.475 / 0.100 / 0.425 → 0.633 / **0.000** / 0.367 — the mesoderm band vanishes and endoderm narrows rather than widens, because `runtime/grn.py` combines a gene's activators as a **mean**, so a second activator that is low where Nodal is high halves SOX17's drive. The CS7 comparison is still falsified on all three layers (mesoderm outside by 0.694). The three-layer order test is now a strict xfail with its reason; the corrected census run sits in its own file beside the 3c4e301 record. Open, each its own registered lane: whether the missing band comes from the unsourced numbers, the mean-of-activators rule (now named as a model assumption for R3), or a missing interaction | done | — |
   | 3 | ↳ *Correction to the row above (2026-09-28, lane-combine, `93caf61`): mesoderm is absent under all four activator rules tested (mean, capped sum, max, OR), so the mean did not remove it; it only sets how much endoderm the second activator costs. The row's causal sentence is wrong and kept as written.* | — | — |

   *(Correction to the wave-1 row above, from the lane that built it: the engine-only tests are **24 files and 157 tests, not 27 and 162** — two files import only the engine but load gate scripts from `scripts/`, so they stay with the application.)*

   **The stop list, which is as much the plan as the waves.** Nothing more is built on the node model
   until wave 0's measured test lands. No more biosamples: 1.1 needs its readings to mean what their
   names say. Nothing is built in area J: all seven steps exist, and only bookkeeping remains. No
   AlphaGenome requests on the 882 blocks: that returns predictions, and the matched table already
   showed this instrument cannot order the tiers. Area E's fate ceiling is not chased: 0 fates decided
   independently of the lineage, generation 7 at chance, and no public dataset resolves terminal
   division; measured contact atlases would strengthen the founder claims and not this one. BioForge's
   confidence is not calibrated: it is uncheckable by construction.

   **The rules every lane follows**, collected in one place: claim the board with `--who` before the
   first edit; write the commit message to a lane-unique file; register before producing any number,
   with its falsifier; commit with `scripts/commit_own.sh`, own files only; send roadmap rows to the
   coordinator as finished text; report a negative first. A lane whose files overlap a running lane
   waits rather than negotiates. The coordinator verifies every report against its result file before
   placing its row, runs the full check before the daily merge, and merges `dev` into `main` at the
   end of each working day and never more than two days apart.

   **Estimate.** About 36 lane-sessions across waves 0–3; at four concurrent lanes, bounded by review
   rather than by the machine, about **1.5 working days**, and the whole unblocked backlog about
   **3 days**. After that the project's limits are measurement it does not hold, not work it has not
   done.
   *(2026-09-29, beside the original: this says more than any result has shown. The coherence pilot
   established only that it did not improve prediction over simpler approaches, not that further work
   cannot improve the biology, nor that new measurement necessarily would.)*

11. **The external review of 2026-09-28, adopted as the next order of work.** A read-only review of the
   code, documents and stored results (tests not run) found the project's strongest asset is traceable
   evidence and its main gap predictive reliability, and set a priority: **context preservation,
   uncertainty reporting and one independent biological benchmark before any expansion of simulation
   scope.** Its nine items, each confirmed by it against a line of code, are adopted with its acceptance
   tests as written; its delivery order is binding: **R1–R5 before any increase in model complexity**,
   then R6–R7, then the joint engine R8, with R9 alongside. Consequences for the rest of this section:
   v0.4 stage 4 (partitioning division) and the E2/E3 widening wait behind R1–R5; lane-joint's pretest
   becomes the first step of R8, not a wave-3 lane of its own.

   | # | Item (review's priority) | Code it names | Acceptance (the review's) | Lane, order |
   | --- | --- | --- | --- | --- |
   | R1 | Context in executable rules (P1) | `attribution/compile.py:173`, `measured.py:970` | a K562-specific rule is inactive in HepG2; conflicting results from two cell types survive compilation and round-trip serialisation | after lane-split (same files) |
   | ↳ R1 | **Done 2026-09-28 (lane-context, `8dcdade` registered, `6a4a8a9` build, `25f81c6` chr21)** | `compile.py`, `measured.py`, `ir/model.py`, `lang/grammar.py` | Every compiled rule carries its cell as an executable `when: cell_type = <cell>`, through the existing `when` clause the runtime already enforces; one experimental rule per element, gene and cell, never the strongest across cells. `unknown` in a `when` matches no context, so a missing context is never universal (grammar regenerated, engine 12/12). **Acceptance met** (`tests/test_rule_context.py`): a K562 rule is inactive in HepG2; a two-cell conflict survives compile → parse → BioIR JSON → parse. Over 24 programs: rules 440,589 and experimental 212 unchanged (no measured link is regulated in two cells); rules without a `when` 440,589 → 0; a chr21 HepG2 run now uses 155 of 5,176 rules; two cell-level nulls the collapse had hidden are named (IL6ST, PPIF). **Negatives:** a run with no `cell_type` now uses no compiled rule; `Rule.applies` and `Event.applies` ignore `a|b`, `absent` and comparisons; the removal guard holds two chr21 lines (one held-out rule still ungated, the 155th in that HepG2 run) and two superseded split tests (strict xfail) until its two-day window passes | follow-up after the window | done |
   | ↳ R1 | **Rule conditions fixed 2026-09-28 (lane-when, `49aca5d` census, `1d8713a` fix)** | `ir/model.py` | Rules and events now read `when` through the one matcher decisions and timers use (`ir.model.matches`), so `a|b`, `absent` and `>=n <=n >n <n` hold on a rule; `unknown` still matches no context. Census first: **0 of 456,157 rule and event clauses used those constructs**, so nothing had failed silently yet and no program changed (135,698,709 old-against-new checks, 0 differ; the chr21 HepG2 run integrates 154 of 5,176 rules). 89 test cases, 31 failing before. Open: `k != v` is parsed but no matcher implements it (0 uses) — the coordinator's decision is that the parser refuses it until a program needs it, since implementing it means deciding what a missing key means; `Event.applies` has no caller | `!=` refusal queued | done |
   | ↳ R1 | **`!=` refused 2026-09-28 (lane-ne, `ff0d767`)** | `lang/parser.py` | `k != v` used to parse into the literal value `!=v`, which no matcher implements. It is now a parse error naming the clause and the alternatives (`a|b`, `absent`, comparisons) on every block that reads a `when`; **0 of 450,205 clauses in the repository's programs are refused**. The contradicting grammar bullets are resolved additively (the old one marked superseded). `!=` stays unsupported until a program needs a rule for what a missing key means | — | done |
   | R2 | CRISPRi negatives read as outcomes (P1) | `measured.py:362` | significant decrease, significant increase, non-significant, missing and invalid kept apart with units, power, study and context; a significant increase never becomes "no effect"; a low-powered null cannot strongly reject; an increase does not establish a silencer | taken into lane-split (it holds `measured.py` and already carries power) |
   | R3 | A validated bridge from annotation to dynamics (P1) | `compile.py:295`, `runtime/grn.py:109` | an end-to-end fixture shows a supported expression response; missing regulator state gives an explicit unresolved-model diagnostic; another citation for the same mechanism does not change its simulated strength; predicted and measured forms of one relation are not counted twice | after R1 |
   | ↳ R3 | **Done 2026-09-28 (lane-bridge, `06b7557` registered, `b986238` amended for A8, `eda2728` build)** | `runtime/grn.py`, `attribution/bridge.py`, `compile.py` | The runtime names its combination rules (activators by **mean**, inhibitors by product) and its declared zero (`@zero`), and reports every regulator without state and every regulated gene without a rate (a warning, `Trajectory.unresolved`, an error when strict); integrated values unchanged, and a full-suite run with the warning as an error passed, so no caller relied on a silent zero. The bridge maps one observation to one strength per (element, gene, cell) through the observable RHO = expression removed / intact (2^x predicted, 1 + f for a one-pair CRISPRi link), with the Hill threshold and coefficient declared assumptions; it does not claim a log2 fold change is a rate constant; measured takes precedence; a conflict or a summarised maximum is reported, never averaged. **Acceptance met** (17 tests): a fixture reproduces −1, −0.4, +0.7 log2 and a CRISPRi −0.3 within 1%; a missing regulator is diagnosed; a second citation leaves the strength unchanged. **Audit: 0 of 440,589 compiled rules can be simulated as compiled**; 76,469 of 237,613 gene–cell pairs (32%) have several mechanisms and cannot be parameterised from single-removal observations under the mean rule; 161,144 lack only declared rates; 39 K562 relations were carried twice in two units (13 now superseded). All 212 CRISPRi links rest on one pair, so the A8 maximum never summarised more than one. Next, each its own registration: the activator-combination rule; declared rates for the single-mechanism pairs | — | done |
   | ↳ R3 | **Activator rule decided on evidence 2026-09-28 (lane-combine, `a6fc5c4` census, `7f3e597` registered, `93caf61` result): mean retained, undetermined** | `runtime/grn.py` (`combine_activators`) | Census: of 41 hand-written programs only SOX17 in the gastrulation demo has two activators; 48,405 compiled gene–cell pairs do, none simulable; 4 tests reach such a gene; `bio test` passes under every rule. **Mesoderm is absent under all four rules**, and the CS7 comparison is falsified under each. Choice test (a case counts against a rule only if no choice of its numbers can produce it; nothing re-tuned; the lane read Bothma before writing the rule and said so): Bothma et al. 2015 enhancer pairs contradict mean and max (additive knirps-late and hunchback-central pairs) and every monotone rule (a snail pair below its stronger single); 81 K562 gene–cells with two or more significant CRISPRi decreases rule out max. **No rule survives**: measured pairs change with enhancer strength, which no fixed function reproduces. The bridge's non-identifiability (76,469 of 237,613 pairs) stands; nothing moved. Open, its own registration: a strength-dependent (promoter-competition) form, which needs numbers the project lacks | — | done |
   | ↳ R3 | **Rates declared 2026-09-28 (lane-rates, `ff2062a` registered, `e921fae` result)** | `attribution/bridge.py`, `runtime/grn.py` | Registered before any rate was read: a measured transcription rate is the gene's **total** rate in the one unperturbed cell it was measured in, never the `max_rate` ceiling no steady-state measurement observes; a half-life is a degradation constant; a copy number becomes nothing. What the split identifies is the element's own contribution, not a saturation constant, and `DECLARED_STRENGTH = 1` is a convention licensing nothing. **Simulable mechanisms 0 → 22,576; compiled rules 0 → 22,580 of 440,589** — inside the registered 15,000–60,000. Still unresolved: **138,543 genes nobody measured**, marked and never filled; **76,469 not identifiable** under the mean rule, byte-identical as the registered invariant required; 25 inhibitions past the split's bound. A borrowed genome median would make 160,392 runnable in relative units only and is never counted as measured. **Every rate is mouse NIH 3T3 fibroblast**, so each of the 22,576 human pairs borrows across species *and* cell type; the matched human-against-mouse half-lives are 1.56 h and 2.73 h, a factor the Spearman of 0.617 does not see. Sources cached, unlicensed, not redistributed | a human synthesis-rate source for the 138,543 | done |
   | ↳ R3 | **No human rate to borrow 2026-09-28 (lane-rates2, `4c3a225` registered, `8a0ef9d` result)** | `attribution/bridge.py`, `build_gene_rates.py` | Registered before any source was read for its numbers: a source must give a per-gene rate in molecules per cell per unit time **with its per-cell calibration named**; a half-life cannot supply a synthesis rate; and `T = N·ln2/t_half` was named in advance as the construction *not* made. Ten candidates read in the primary source. **The only absolute per-cell rate published is mouse** (Shao 2022, mouse ESC, "cell⁻¹ min⁻¹", not distributed per gene). **The only human genome-wide transcription rates are a mouse constant rescaled**: Hausser 2019's HeLa rates are an abundance times one global decay rate times ~225,000 mRNA per cell — a number the authors say they could not measure and obtained by scaling **the 180,000 mRNAs of a mouse 3T3 cell** by a volume ratio. Adopting it would have laundered the exact constant this work set out to escape. The construction was refused: no human per-gene copy number exists genome-wide, and where it can be checked it disagrees with the one measurement by a median 0.236; a test now fails if any rate in the table is that construction. **Every count is unchanged, exactly as registered.** The whole of the change is two new numbers: **0 pairs stand on a human rate, 22,576 on a borrowed-species one**, now counted per pair rather than stated once in a header, so no later result can quote the 22,576 without the species travelling with it. Next, each its own registration: Ietswaart 2024 measures human K562 **and** mouse NIH 3T3 by one method, which would replace the present species test (mouse MEF against human K562, confounding species with cell type) with a matched one | a matched species test | done |
   | ↳ R3 | *Correction, 2026-09-29, lane-s5: the "1.56 h against 2.73 h" quoted in the two rows above is Schofield's own pair of cells (human K562, mouse MEF), not the half-lives the rate table carries, which are Schwanhäusser's NIH 3T3 values (median 9.96 h). Per gene, the K562 half-life is a median 0.23 times the table's value, and **most of that gap is measurement method, not species**. The factor the rows attribute to species is therefore mostly a method difference.* | — | — |
   | R4 | Confidence, effect size and evidence quality separated (P1) | `compile.py:271`, `organism/forge.py:170`, BIOFORGE-CONFIDENCE.md | effect magnitude alone never raises certainty; every probability names its outcome, calibration population and method, and an unavailable probability stays unavailable; the calibration text corrected (a constant can be calibrated to a base rate and still not discriminate); tests that only preserve the old constant updated | after lane-h releases `forge.py`; compiler half after R1 |
   | ↳ R4a | **BioForge half done 2026-09-28 (lane-conf, `38087a1` tests, `16a0434` build, `3b29509` docs)** | `genomeos/certainty.py`, `organism/forge.py` | Evidence category, effect with unit, measurement uncertainty, model score and probability kept apart; a `Probability` cannot exist without a `Calibration` naming outcome, population and method, and a missing one stays `None` with its reason. The fixed `confidence=0.3` is retired: a design answer carries effect = wild-type loss − answer loss, score = loss, probability None ("no calibration record"); the emitted `experiment` states no confidence. Effect magnitude alone never raises certainty (four designs and an effect sweep). BIOFORGE-CONFIDENCE.md corrected beside the original: a constant can be calibrated to a base rate (error 0) and still not discriminate (AUC 0.5). Line-number pins that only preserved 0.3 replaced by property tests. **Negative: its first commit imported a module it had not committed and blocked every push for about an hour** (guard added, `26dc738`). Open for R4b: `compile.py:271`, `attribution/budget.py:60`, and `attribution/human_panel.py:1321`, where confidence rises with the size of a gap, the pattern the review names | R4b: the two attribution files now, compile.py after R1 | R4a done |
   | ↳ R4b | **Attribution half outside compile.py done 2026-09-28 (lane-conf2, `ef31028` census and strict-xfail tests, `3c5e122` build)** | `human_panel.py`, `variation.py`, `budget.py` | Census: **five confidences rose with effect size** (`block_class` lineage_restricted, polymorphic `0.3+min(0.3,gap*3)`, core, variable; `variation.case_of` `0.3+0.3*min(dm,dh)`). They now carry a Certainty record: effect in its unit, Poisson SD for counts, probability None with the reason; the CLI's "conf" column is gone. `budget.guess` (with `:60`) is a hand-set evidence-quality score that does not move with the constrained fraction within a label; it stays numeric because `compile.py:107` emits it, labelled, with a record beside it. The `test_variation` pins of 0.6 and 0.3 are replaced by a property test. No stored result rewritten; on their next run 49,078 panel rows and 29,481 variation cases (all 0.3–0.6, or 0.0 unplaced) lose `confidence` and gain the record. Open, for the compiler half after R1: `compile.py:276` (effect magnitude read as confidence, emitted at :294 and :304), `:260` default 0.4, `:269` constant 0.9, `:107`; outside attribution, the per-rule constants in `genome/unknown.py` | compile.py after R1 | done outside compile.py |
   | ↳ R4c | **Compiler half done 2026-09-28 (lane-bridge, `96bc5e5`)** | `attribution/compile.py` | Predicted elements and rules state no confidence; each carries its effect in log2 units and "probability unavailable", with the full Certainty record in every program header. The gene 0.9 and domain 0.4 constants are gone; the region budget score stays, labelled hand-set, not a probability. Property tests: effects from 0.01 to 5 give the same stated certainty. Over 24 programs: block confidences 486,255 → 45,878; rules stating a confidence 440,589 → 212 (the measured ones); 18,973 gene and 14,279 domain constants removed. **Negatives:** the evidence explorer's "weak" count rose 812,921 → 928,094 because it reads an unstated confidence as 0.0 (`evidence.py` must tell unstated from low); `confidence_calibration.py` still recomputes the old formula; the committed chr21 rule lines are held by the removal guard until 2026-09-30 | the two follow-ups queued | **R4 done** |
   | ↳ R4d | **Follow-ups done 2026-09-28 (lane-evid, `8d57243` calibration, `262ea4e` explorer)** | `evidence.py`, `ir/model.py`, `lang/parser.py` (two lines), `attribution/confidence_calibration.py` | A missing `confidence:` now parses to `UNSTATED` (equal to 0.0 in any calculation, told apart by `confidence_stated()`; a stated 0.0 exists in eight places, so value alone could not tell them). **The difference lives only in parsed programs: BioIR JSON writes it back as a stated 0.0.** The Evidence explorer (CLI, CSV, by-program, web) counts three categories: **986,941 facts = 19,293 stated weak + 48,114 stated strong + 919,534 unstated** (the old rule called 938,827 "weak"); means are over stated values only. The calibration calibrates nothing: no predicted fact states a confidence, and measured blocks are hand-set ranks. `model_score_curve_genome` describes the unclipped |log2 fc| against outcomes with probability None: lentiMPRA ordered over 6 bands (0.148 → 0.366), CRISPRi all-matched ordered (0.052 → 0.127, a bound); CRISPRi tested pairs and VISTA still refused. Found: `predict/enhancer_target.py` stores a confidence of |lfc| capped at 0.7 (differs from |lfc| on 39,540 of 43,681 chr1 links), which the compiled header calls equal to |lfc|; two scripts still count unstated facts as weak — both now in lane-evid2 and lane-ontology | lane-evid2 running | done |
   | ↳ R4e | **Leftovers done 2026-09-28 (lane-evid2, `a378182` census, `2085279` links, `37172b6` weak counts)** | `predict/enhancer_target.py`, two census scripts | New predicted enhancer–gene links carry the signed log2 fold change, its unit and a Certainty record (model score = |lfc| unclipped, probability None); no confidence is written. The 1,052,700 stored links keep their capped confidence (it differed from |lfc| on 955,106) and readers accept both forms; `decompile`, the CLI and the node audit read the new form. The census scripts count weak only on stated confidences: **14,088 of 960,096 facts are weak, not 928,094; 914,006 state none** (`evidence_compiled_genome_r4`, `measured_layer_genome_r4`). Open: compile.py `element_certainty` should read `link_score()`, or new links get no model score; magnitude-derived confidences remain in `alphagenome_adapter` (`VariantEffect`) and `genome/missense.py` | R4f queued | done |
   | ↳ R4f | **Last leftovers done 2026-09-28 (lane-conf3, `1219917` census, `7485df0` compiler, `6e833c3` adapter and missense)** | `compile.py` (`element_certainty` only), `predict/alphagenome_adapter.py`, `genome/missense.py` | `element_certainty` reads `link_score()`, so links written since R4e keep a model score (chr21 recompiles byte-identical). The variant adapter's `min(0.7, 0.2 + 0.25|lfc|)` and missense's `round(min(0.7, score), 3)` are replaced by Certainty records with probability None; adapter rules state none; the CLI and API write the record. Property tests over 307 effect sizes and 200 scores change nothing but effect and score. **No stored figure moved**: 40 confidences in `alphagenome_rs12740374` keep their values; the missense confidence was never stored. The web tooltip and ALPHAGENOME.md no longer say "capped at 0.7" (coordinator). Open: the compiled header still names the score as the capped stored field, rewritten at the next recompile | — | **R4 done** |
   | R5 | Evaluation independence through provenance (P1) | `measured.py:199` | immutable study, assay, source-file and partition identifiers on every record; validation outcomes kept out of features, candidate selection, starting labels and search objectives; automated overlap checks including related intervals and variants; reused benchmarks reported as reused | lane-split, running (its audit is the first half) |
   | ↳ R5, R2 | **done 2026-09-28 (lane-split: `0af7e1b` audit, `36fa298` split and provenance, `1cff24a` outcomes)** | `attribution/measured.py`, `compile.py` | **Audit: no reported held-out figure was exposed.** Every held-out claim reads the two files apart and fits on training only; the merged loader fed evaluation and display only. The defect was latent: 43 of 212 compiled measured rules were held-out-only links, now kept out of `targets:` and marked `HELDOUT_MARK`. One caveat recorded: `target_rebanding` and `union_axis` pool training and held-out K562 in their band rates by registration, so no later held-out K562 figure may be read against the re-banded table as out of sample. Every record now carries split, source file, study, assay and five power columns; `development_only()` refuses held-out pairs; the overlap check finds **0 identical pairs but 249 of 4,378 held-out pairs sharing bases with a training interval** (36 on the same gene); the held-out file is flagged a reused benchmark (≥10 scored results read it). R2: five outcomes, a null informative at PowerAtEffectSize20 ≥ 0.8 (6,169 of 9,810 negatives); an increase raises no rule or silencer. CRISPRi disagreements 31 → 23; tested-gene agreement 0.7578 (n 128) → 0.8083 (n 120); calibration ECE 0.4108 → 0.4507 and the **matched difference 0.137 (p 0.070) → 0.090 (p 0.156)**, weaker; the 10-of-12 band headline unchanged. Open: the committed chr21 program's two R2 lines wait for R1's recompile; `crispri.Pair` (the scorers' own record) still lacks power; the 2026-09-17 ATTRIBUTION.md sections still quote 128, 0.7578 and 4,613, superseded by the dated section | R1 now running (lane-context) | done |
   | R6 | Assay observations kept before aggregation (P2) | `measured.py:390` | one strong reporter tile cannot set a region's label; maxima only as labelled descriptive statistics; conflicting observations inspectable; lentiMPRA described as integrated, not episomal | after R2 |
   | ↳ R6 | **Done 2026-09-28 (lane-assay, `e0c911a` census, `bea3478` registered, `180ad41` + `f063eff` build)** | `measured.py`, `mpra.py` | Census: the lentiMPRA per-cell label was **the maximum over matched tiles, the only summary that decided a label**; of 17,869 elements 1,006 match several tiles, and in 237 cell readings the tiles disagree, all of which the maximum read as active. Now every tile is kept with interval, strand, overlap and value; a cell is labelled by the share of tiles at or above 1.0 (a tie is `reporter_tiles_conflict`, counted apart); `activity` is the median and the maximum survives only as a labelled descriptive field no verdict reads; every assay block carries `outcome_kind`; saturation-mutagenesis repeats are read per base. **Acceptance met:** one strong tile among inactive ones no longer makes an element active. Moved: lentiMPRA agrees 3,361 → 3,231, disagrees 14,508 → 14,512, 126 conflicts; calibration narrow arm ECE 0.0905 → 0.0963. The ENCODE files carry no replicate values, so no per-tile uncertainty exists to keep. "Episomal" corrected in code and annotated beside the original in the docs. Rebuilt results sit beside the old ones (`*_r6.json`), the old names held by the removal guard. Open (A8): a compiled experimental rule's strength is the largest training |effect| per element, gene and cell — passed to R3 | A8 in R3 | done |
   | R7 | A compositional ontology (P2) | `compile.py:283` | origin, role, activity, target and evidence status separated; a repeat-derived regulatory element keeps both; a predicted increase does not force a silencer label; alternatives survive serialisation; constraint read as selection, not mechanism | after R3 |
   | ↳ R7 | **Done 2026-09-28 (lane-ontology, `2ba726b` census, `8812e1a` registered, `329c2f9` build, `d2ab0ae` header fix)** | `compile.py`, `lang/grammar.py`, `lang/parser.py` | Elements and regions state `origin`, `molecular_role`, `activity`, `target_relation` and `evidence_status` from a closed vocabulary (`grammar.AXES`); `,` joins values that all hold and `|` unresolved alternatives; the parser refuses other values and BioIR JSON round-trips them. `class:` is derived from the ENCODE registry role; region `role:` kept verbatim (Module.unknowns reads it). Constraint is evidence of selection only, and no value means "no function". **Acceptance met** (`tests/test_ontology.py`): a repeat-derived dELS keeps both properties; a predicted increase is `represses_target` with `silencer|insulator_like|competing_promoter|unknown`, never a silencer label; alternatives survive serialisation. Census: all 459,449 elements had said `class: enhancer`, including 156,925 whose rule represses. Over 24 programs every non-axis line is unchanged; new: 101,011 predicted elements repeat-derived and 81,058 partly; 193,027 insulator-like beside enhancer-like; 21,283 regions now say no sign of selection instead of neutral. **The classifier calls a block regulatory before it reads RepeatMasker**: of 15,536 regions with a registry or CpG role, 8,709 are mostly repeat, which it never read. **Negatives:** the derived class moved nothing (every predicted element is pELS or dELS, so the old constant happened to match); the budget's "best guess neutral" stays on 2,546 regions; no predicted element has a constraint reading | budget relabel and classifier order queued | **R7 done** |
   | ↳ R7 | **Follow-up done 2026-09-28 (lane-budget, `b6e6d17` census, `bbce3d8` registered, `308484a` build)** | `genome/unknown.py`, `attribution/budget.py` | The classifier now reads RepeatMasker before any class rule (it had never read it for 15,320 of 26,806 blocks), so a regulatory block can also be repeat-derived: 8,709 of 15,536 are, 6,534 partly. New budget writes (`budget_axes_<chrom>`, `budget_axes_genome_wide`) state `evidence_status` (constraint as selection only), `origin` and `legacy_tier`; fossil → repeat_unconstrained, neutral → unconstrained_unknown, and eight labels renamed one to one — "best guess neutral" on 2,546 blocks now reads "no sign of selection, which is not a sign of no function". Every tier count, base pair, score and stored file is unchanged, and no consumer figure moves, since all still read `budget_<chrom>`. **Negatives:** origin differs from the compiled programs on 11 blocks (coverage just under one half, rounded by the compiler); 16 blocks without a constraint reading keep `constrained_unknown`; the legacy keys stay until about 20 consumers move; the census commit's subject says 24,337 blocks where the count is 23,882 (corrected in the next commit's message) | move consumers to `budget_axes`; regenerate `unknown_<chrom>` | done |
   | ↳ R7 | **Consumers moved 2026-09-28 (lane-consumers, `775d4be` registered, `690a71c` build)** | `budget.read_axes` and 18 direct readers | Every budget reader goes through one reader (18 direct; 7 more through `organise.blocks` or lexicon's context map); consumers join on `legacy_tier`, so stored keys stay. The CLI, the web budget card and 23,708 compiled role lines use `repeat_unconstrained` / `unconstrained_unknown` and the R7 labels; no program's role line says neutral, fossil, dead or best guess. **Registered invariant held: 0 numeric differences** across 24 chromosomes and 19 consumer outputs (9.9 million numbers in the programs, 4.4 million leaves in results). Behaviour change: `genomeos budget` without `--chrom` prints the stored R7 sum instead of overwriting `budget_genome_wide`. Open after 2026-09-30: regenerate `unknown_<chrom>` so the compiler reads the classifier's origin (changes origin text on 11 blocks); two guard-held "best guess" strings in variation.py; renaming legacy keys in stored results is a separate schema change | after 2026-09-30 | done |
   | R8 | The joint inference engine, bounded (P2) | `attribution/closure.py`; genomeos-8a's brief | pretest first (coupling against independent scoring on validation folds; a negative closes it); then fixed-boundary search solving a synthetic two-change case, and gains on independent evidence over unchanged annotation and independent-block baselines; convergence never reported as validation | after R1–R5 |
   | ↳ indep | **Independent benchmark probed 2026-09-28 (lane-indep, `042a4aa` registered, `0924a5c` result): NOT READABLE** | `attribution/indep.py` | Census: every CRISPRi, reporter and eQTL set the attribution layer has read is listed; nothing from the IGVF portal had been. Chosen: IGVF K562 MHC CRISPRi Perturb-seq (Gersbach, IGVFFI4093WUVB, CC BY 4.0, 34,279 pairs), registered before scoring (element-bootstrap AUPRC against distance, the compiled target, the node rule and scE2G v1.2; 0 requests; uncovered pairs missing, never imputed). Coverage was fine (68% within reach), but **only 13 covered distal decreases against a floor of 20**: 74 of the screen's 93 decreases are promoters, and its distal hits are 181 increases against 19 decreases, which anyone reusing it needs the authors' reading of. **The project still has no independent target benchmark; every target claim rests on the reused ENCODE file.** Next, by cost: HCT116 (705 requests, the owner's decision), the IGVF K562 many-loci Perturb-seq (0 requests, needs single-cell processing), five ENCODE FlowFISH loci (0 requests, too few alone) | the owner's HCT116 decision | open |
   | ↳ indep | **Second probe 2026-09-28 (lane-indep2, `12ad319`): nothing scorable at 0 requests** | `scripts/indep_igvf_probe.py` | IGVF has released no processed table for the 16 many-loci K562 Perturb-seq sets (raw reads only, 691 GB against 19 GiB free; their analysis sets IGVFDS3624DGBH and IGVFDS3481OUHP exist but are unreleased). Their guide library, read without labels, is covered 0.58 by the sweep (161 of 279 distal elements, 5,676 pairs), so it clears the coverage floor the day a table is released. The one other released K562 table (DC-TAP, IGVFFI0957PYTA) is Ray 2025, already in the ENCODE held-out file, and a remainder of a study already read was not registered as independent. Every other released element-level table is in WTC11, WTC11-derived cells or HCT116, where the sweep has no track. **So the cheapest independent benchmark is HCT116 at 705 requests (the owner's decision); the free one waits on IGVF's release** | pre-register on IGVFDS3624DGBH when released | open |
   | ↳ indep | **HCT116 arm run 2026-09-28 (lane-hct116, `6257b95` protocol, `a39073d` result; the 705 requests approved by the owner): passes, NOT replicated** | `scripts/crispri_hct116.py` | The fe61c36 registration applied as written with the frozen weights; 705 requests, 0 refused, 0 retried, each asking for ALL_FOLDS. On 363 covered HCT116 pairs (34 regulated) the deletion gain is **+0.022, 95% −0.058 to +0.145**; weighted over all 396 pairs +0.014 (−0.042 to +0.126). It passes the registered rule, but the interval includes zero, so **the CRISPRi result stays replicated in one cell type**. Why it could not decide: the 34 positives are five genes at five loci (SSFA2, FAM3C, MYC, KITLG, CCND1), so the bootstrap has five units, and 228 of the 363 covered pairs lie beyond the model's 1 Mb window (37% answered against 66% in held-out K562). Model versions labelled "mixed", but on these 705 elements the ALL_FOLDS answers equal the stored sweep on 54,760 of 54,760 values, so the mix explains nothing. A decisive second cell type needs more genes within 1 Mb, not more requests | the IGVF release | done |
   | ↳ 1.3 | **Clause 2 re-tested against a matched control 2026-09-28 (lane-clause2, `e1dbcf3` registered, `0c8b82d` result): still not met** | `scripts/clause2_matched_control.py` | Registered before any matched figure existed, then run once on 24 chromosomes, 0 requests; the gate passed (with matching off it reproduces `d717b28`'s pooled counts to the digit). Matched on block length **and local coding-gene density**, the real unknown names a coding gene **−27.25 points** below its own windows (95% over blocks −30.91 to −23.58; over chromosomes −30.35 to −24.40; n = 531 blocks, none undrawable). **The registered falsifier fired**: −27.25 is below −18.62, so on the same estimator density closes about 11 points of the unmatched −38.62 gap. The neutral tier sits exactly as far down (−27.15), and real unknown minus neutral is −0.10 (−4.55 to +4.30). Sensitivities agree (−28.48 matching on distance to the nearest TSS; −20.06 rejecting only on the target set). **A residual confounder is named, not hidden**: the TSS covariate is matched (4.11 against 4.05) but a carrying window holds 13.3 scored elements against 6.2 per carrying block, and more elements mean more chances to name a gene. The measured-arm route stays registered and unrun | a measured arm per block | clause 2 not met |
   | ↳ R8 | **Pretest done 2026-09-28 (lane-joint, `dcc277b` registered, `4098753` result): NEGATIVE, R8 closes** | `attribution/joint_pretest.py` | CRISPRi training split only (held-out never opened, a test spies on the file opener); leave-one-chromosome-out folds; R2 outcomes; 3,713 scored K562 pairs (423 decreases, 3,290 well-powered nulls); 0 requests. Independent predicted drop AUPRC 0.766; **element competition 0.754 (gain −0.013, 97.5% component interval −0.024 to +0.001); gene budget 0.747 (−0.020, −0.027 to −0.011)**; the self-share control 0.758. Neither interval lies above 0 and neither beats the control; on pairs with a scored partner, where coupling could matter, both lose more. 4,529 pairs missing from the cache, not imputed. **No joint search engine is built**; the synthetic two-change case and the search are cancelled. Scope: coupling of the deletion score as held, K562 only, on pairs the sweep answers (0.766 is not comparable with all-pairs figures) | — | closed |
   | ↳ R8 | *Correction, 2026-09-28, second review (item 12, S3): the row above overstates its negative. The pretest rejected **two particular score transformations** (element competition and gene budget) on the available K562 pairs; it did not show that joint inference cannot help. "No joint search engine is built" and "cancelled" are withdrawn as conclusions; the two transformations are retired, the negative stands for them, and a bounded pilot reopens the objective.* | — | reopened as S3 |
   | R9 | Result manifests and public claims (P2) | `genomeos/results.py:23`, README | every result carries source accessions and versions, checksums, assembly, coordinate convention, code revision, parameters, exclusions and partitions; a second environment rebuilds one representative result or names what is missing; public summaries separate sequence, annotation and assay coverage, software correctness and independent prediction | now; README corrected 2026-09-28 by the coordinator |
   | ↳ R9 | **manifests done 2026-09-28 (lane-manifest, `fe0880a`, `4b8e8b9`; evidence fix `c22c8c6`)** | `genomeos/manifest.py`, under `result_manifest` | Census before: **0 of 955 results carried a manifest**, 718 said nothing on any field, none recorded a code revision or coordinate convention. The writer stamps git sha, dirty flag and argv; a new result name must be complete (written, then `ManifestError`); rewriting a historical result warns; nothing historical rewritten. `crispri_direction` rebuilt **byte-identical** from its manifest in a second clean worktree with a fresh offline venv, all 3 input sha256 matching. On another machine the missing dependency is the unpublished AlphaGenome all-element sweep (model version unpinned); the CRISPRi tables are pinned only by sha256 (upstream `main`). `genomeos evidence` now fails loudly when a program does not parse. Side finding fixed: uv.lock named genomeos 0.9.0 against 1.0.0 | open: the other writers still need manifests | done, one writer |
   | ↳ R9 | **Headlines rebuilt 2026-09-28 (lane-manifest2, `42e786f`, `190a144`, `65d2d23`, `f48b909`)** | the five README headline writers | **Every headline reproduced: 0 values differ** over two reruns each in fresh worktrees with offline environments and 0 requests — node containment +2.90 (3,004 values) and +5.89 (2,607), real unknown 62.3% vs 86.0% (3,059), benchmark 9/9 (324), assay coverage 160,447 of 30,602,182 bp (144). `manifest_rebuild.py` matched 96 of 96 inputs by sha256 on node_containment_audit. The one mismatch was not a number: a note added to the benchmark result by hand in `2805552`, which the script now writes. `constrained_unknown_targets` used to check its 2026-09-16 reproduction against whatever file was on disk, which after one rewrite compared it with itself; it now reads that run from git. Open: `crispri_published` (now unblocked); `code.dirty` counts result files written by an earlier writer; the rebuild should ignore "seconds" timings | crispri_published next | done, five of six |
   | ↳ R9 | **Sixth headline rebuilt 2026-09-28 (lane-manifest3, `0126346`, `86e3c52`)** | `crispri_published`, `manifest.py`, `manifest_rebuild.py` | **All six README headline results now reproduce from a clean checkout.** The CRISPRi result equals the committed file in all 256 values over two reruns in fresh offline environments (0 requests); a third checkout by `manifest_rebuild` matched all 4 inputs (1.43 GB) by sha256. Contract gaps closed: `code.dirty` ignores data/results/ and lists those paths apart; the rebuild ignores wall-clock timing keys at any depth and lists them. `crispri.Pair` carries the five power columns, with no figure moved. Found in passing: the post hoc DNase-only 0.639 sits 0.008 above ENCODE-rE2G's published interval (0.631) while the summary says "in the range of"; the coordinator added the three reasons beside the number in CRISPRI-RESULT.md rather than changing either | — | done, six of six |
   | ↳ R9 | **Model dependency recorded 2026-09-28 (lane-pin, `837d5da`)** | `manifest.py`, five headline writers | The AlphaGenome all-element sweep ran on client alphagenome 0.9.0 (uv.lock, unchanged since df1a184) against the v1main DnaModelService (RNA_SEQ GeneMaskLFC scorer, 1 Mb window), 2026-09-12 to 2026-09-16. **The model version is unpinned: none was requested, and no response or file records one** (ALL_FOLDS is only the client's documented default). The track table was never stored; a fingerprint of the cache stands in (371 tracks modal across 966,615 answers, 316 tissue names, sha256 bdf63a52). Recorded as `model_dependency("alphagenome")`, carried by the five headline results that read the sweep (additive, 0 lines removed) and their writers; every new cached answer records client version, requested model version, API, date and a track checksum. 0 requests. **Coordinator's decision:** future runs request the model version explicitly, so what a new answer came from is known; the stored sweep's version stays unknown and is labelled so | explicit model version, queued | done |
   | ↳ R9 | **Model version pinned 2026-09-28 (lane-pin2, `981b5b0`)** | `predict/alphagenome_adapter.py`, `enhancer_target.py`, `manifest.py` | Every AlphaGenome client is made by `create_client`, which requests `ALL_FOLDS` through the installed 0.9.0 client's `model_version` parameter (read from its source; its enum names ALL_FOLDS as the default when none is set); a test fails if a client is created anywhere else. Cached answers carry a `model_version` label, "unrequested" for the stored sweep or "ALL_FOLDS", and are never asked again; a result reading both states "mixed" in its manifest, because the two are not guaranteed to be the same model. 0 requests | — | done |

   **What the project already agreed with.** R5's merge defect and R2's missing power column were
   reported independently the same night by genomeos-8a and confirmed by the coordinator; lane-split
   was running on them when the review arrived. BIOFORGE-CONFIDENCE.md already said the fixed 0.3 is
   not a probability; the review corrects its reason.

12. **The second external review, adopted 2026-09-28: the original objective is unfinished, and the
   status measures used so far do not measure it.** The review credits explicit uncertainty,
   context-specific rules, assay distinctions and reproducible headline results, and then states what
   this plan had lost sight of: **the objective was to revise genomic labels jointly and show that the
   revisions improve biological accuracy, and neither milestone percentages nor test counts measure
   that.** Its eight items are adopted with its delivery order. Four of its code claims were checked
   before adoption and hold: a confidence the program did not state is written to BioIR JSON as a stated
   0.0 (`ir/model.py`); a new result that fails its manifest check is still written, and a retry then
   finds the file, treats it as historical and only warns (`results.py`); the revision stamp ignores
   untracked source (`manifest.py`, `--untracked-files=no`); and CI never runs on a pull request's
   later pushes (`ci.yml` omits `synchronize`). A fifth, found in checking the fourth: **GitHub runs the
   daily scheduled check on the default branch, which is `main`, so the daily second opinion has been
   testing a branch 170 commits behind the one being worked on.**
   *(Correction, 2026-09-28, lane-s7, `b50c873`: **the fifth claim above is wrong and was the
   coordinator's.** The repository's default branch is `dev` (GitHub API), and all 13 scheduled runs since
   2026-09-16 tested the tip of dev; the claim came from this checkout's `origin/HEAD`, which a clone sets
   once and which was never updated. **What holds instead is worse: the daily run has failed every day
   since 2026-09-22 and nothing waited on it**, on one test (`test_cancer_alterations.py::
   test_a_five_prime_partner_keeps_its_own_ectodomain`, 1 failed of 1,992) that reads git-ignored caches
   without the `@needs_caches` skip its eight neighbours carry. The local pre-push check links this
   machine's caches into its worktree, so it passed the same commits — **"every push passed the full
   suite" was true of the local check and false of CI for seven days.** And main moved twice (`36f36d6`,
   `45a5aa6`) without a green run on either sha, because main's rule requires the check but does not bind
   admins.)*

   **Three status measures, kept apart from now on:** software delivered; biology independently
   validated; release readiness. Discontinued investigations are recorded separately from implemented
   capabilities.

   *(Built 2026-09-28, lane-status, `9f18960`: `/api/state` now returns the three measures and a separate
   discontinued list, computed from files and never summed. **Software delivered** (software, not knowledge):
   6 of 7 milestones; item 11 **8 of 9 — R8 reopened, not done**, since an item's state is now read from its
   follow-up rows and a row recording a negative never closes its item; item 12 3 of 8; engine checks 12 of 12,
   built from a dirty tree. **Biology independently validated**, seven categories: **one counted**, the CRISPRi
   result in README's wording with its qualifications (HCT116 passes, not replicated; not a paired test; a
   reused benchmark); node containment shown, not counted (not independent); clause 2 not met; the therapeutic
   benchmark built knowing its answers; **none** in function, activity, context, origin, therapeutic target and
   label revision. A claim whose required qualification cannot be read is withheld, and the claim panel's bare
   band field ("above the published interval") is removed. **Release readiness** lists what is missing and
   what is met, reading CI from a cached file that shows its age; the promotion gate is only ever dry-run.
   **Discontinued**, counted by no measure: nine investigations, each quoted from its own row. The CI cache is
   refreshed with `scripts/ci_status_cache.py`.)*

   | # | Item | Code it names | Acceptance (the review's) | Order |
   | --- | --- | --- | --- | --- |
   | S1 | Uncertainty survives serialisation | `ir/model.py:43` | unavailable confidence is explicit and preserved through parsing, JSON, exports and summaries; the test that accepts the loss is replaced by a preservation test | wave 1 |
   | ↳ S1 | **Done 2026-09-28 (lane-s1, `1307fee` census, `92654f5` build)** | `ir/model.py`, `lang/parser.py`, `lang/tools.py`, `bio.py`, `runtime/{uncertainty,located,debugger,methylation}.py`, `cli.py`, `web/server.py`, `index.html`, `organism/{experiment,forge}.py` | BioIR JSON writes an unstated confidence as `null` and marks the file `records_unstated_confidence: true`; a stated 0.0 stays 0.0; every IR default is `UNSTATED`; means are over stated values with the unstated counted beside, nothing unstated is weak, and weakest-link scores name unstated links. **919,533 of 919,533 unstated confidences now survive BioIR JSON (0 before).** Check weak-rule lists 445,072 → 4,695 (440,377 unstated counted apart); check means moved in 26 programs (live reports only). Evidence explorer counts, its CSV and the uncertainty reports of 50 programs are unchanged; no stored result moved. The loss-accepting test is replaced by a preservation test; the 18 census negatives pass. A file written before the mark cannot be told apart and is read as the stated 0.0 it carries; none exists. Found in the census: the parser gave a transcript that stated 0.0 its gene's value (`cconf or conf`); no figure depended on it. Engine package 12/12 | — | done |
   | S2 | Calibrate the power model before any sample size is used | `scripts/clause2_design_power.py:460` | the simulation reproduces its own anchor's quantity (the review reports that it takes a 13.45% window-level anchor and produces 16.80% element detection and 56.96% positive windows; reproduce first); sensitivity and false positives used consistently; observed and latent correlation distinguished; the planned analysis simulated under the null and alternatives, with shared controls and clustering; designs compared at equal total assay cost; sizes above the eligible block population flagged. **Until then 20–5,000 is exploratory, not an experimental specification** | wave 1 |
   | ↳ S2 | **Done 2026-09-28 (lane-s2, `b903e1e` reproduction, `019a9ab` registered, `cab70d9` result)** | `scripts/clause2_design_power.py`, `clause2_design_power_calibrated.json` | The review's 16.80% and 56.96% reproduce exactly from the committed code, each traced to its cause. Recalibrated at the window level with observed-scale correlation, per-element sensitivity and false positives on both sides, shared controls and chromosome clustering, designs costed in tested elements and capped at their eligible population; **the anchor now reproduces** (13.459% against 13.45%). **At 80% probability of detection, a three-quarters rate is infeasible in all 24 designs** (best 0.466 at every eligible block). Shared controls are anti-conservative in all 105 cells. lane-design's reference 750 / 100 / 30 / 20 / 20 becomes infeasible / 187 / 75 / 50 / 30, and 50 of its entries exceeded their own eligible population. **The measured arm's own reading rule reads "model_failed" in 53–82% of experiments at a three-quarters rate** — the rule, not only the sample size, needs revisiting. Open: the result was produced through a scratch driver outside the repository, so its manifest's argv cannot be rebuilt; a rerun through `--calibrated` is queued | rerun via the committed entry point | done |
   | ↳ S2 | **Follow-ups done 2026-09-29 (lane-rule, `1cb7559` registered, `1ef8990` build, `d85fc5d` rebuild, `4f094bf` result)** | `scripts/clause2_design_power.py`, `scripts/clause2_measured_arm.py`, `clause2_reading_rule.json` | **Rebuilt through `--calibrated`: 25,003 of 25,003 values, `manifest_rebuild` 201 of 201 inputs, 0 differences** — the row above's open item is closed. **Reading rule revisited: no rule on the estimator decides at the registered error rates within the eligible blocks** (best 0.376 and 0.328 against a registered 0.8 at both poles); the committed rule, which never says "cannot_decide" above the floor, is kept as the record with the revised one beside it | — | done |
   | S6 | The evidence contract enforced | `results.py:52`, `manifest.py:83` | a failed new result goes to a separate location and only validated results enter the registry; an explicit legacy allowlist; untracked executable source counted in the revision stamp | wave 1 |
   | ↳ S6 | **Done 2026-09-28 (lane-s6, `0e50036` census, `e1464b2` registered, `f0dbfd0` build)** | `results.py:save_result`, `manifest.py:code_revision`, `data/results_legacy.txt` | Census before: 221 `save_result` sites in 148 files; 1,015 results, 955 historical at `fe0880a` (656 tracked, 299 git-ignored); all 60 new names complete in all 69 committed versions, so no committed result shows the retry path — the defect was real and unexercised. Now a failed new result goes to git-ignored `data/quarantine/results/` with its reason, never to `data/results/`, on every attempt; historical means a committed list of 955 names generated once from git history (`--check`: 0 differences) and never extended by code; a missing list counts as empty, so enforcement fails closed; an allowlisted result still regenerates with a warning; untracked code under genomeos/, scripts/, tests/ or any `*.bio` marks the stamp dirty and is named, untracked data does not. Census after unchanged (1,015 / 71 / 67; 0 quarantined), as registered. 19 new tests; the 2 tests that asserted the defect now assert the fix. **Open: 13 writers write or edit results in `data/results/` without `save_result`** (11 write, 2 edit in place), outside the contract | the 13 bypassing writers | done |
   | ↳ S6 | **Every writer under the contract, 2026-09-29 (lane-contract, `d6a77dd` census, `07caa6c` clause 2 corrections, `8934a72` build)** | `scripts/results_writer_census.py`, `tests/test_results_writers_guard.py`, `organism/provenance.py` | A static taint census found **36 sites in 28 files** writing into data/results outside `save_result`; lane-s6's 13 were a lower bound by 7. All 20 JSON writers now go through `save_result` with complete manifests; the two in-place editors became results of their own (`crispri_hct116`, `clause2_statistical_corrections`), committed blocks untouched. **No figure moved: 22 outputs run before and after in clean worktrees, 0 differences.** **The HCT116 arm's manifest had lacked sources and exclusions all along** — nothing checked it while it sat inside another result — and the contract quarantined it until they were stated. A guard test fails on any direct write under genomeos/ or scripts/ (23 at base, 0 now, a planted write caught); 2 exceptions and 11 fetched-data sites listed out of scope. **Side findings, not fixed: `celegans_fate_reads` no longer reproduces its committed file (82 values differ); `celegans_commitment` depends on Python's string-hash seed (12 values move in the 4th decimal).** Also: two census scripts headed Apache-2.0 imported AGPL application code, contrary to LICENSING.md; the coordinator set both headers to AGPL-3.0-or-later, as the table assigns `scripts/`. Open: fetched data kept under data/results; 3 registry files with no code writer | the two C. elegans results | done |
   | ↳ S6 | **Two worm results reproduce again, 2026-09-29 (lane-repro, `b587bb1` code, `236fa2b` results)** | `organism/commitment.py`, `scripts/celegans_fate_reads.py`, `tests/test_celegans_determinism.py`, `tests/test_fate_rules.py` | **`celegans_fate_reads` moved at `9d42485` (§7.5 built), an intended runtime change whose result was never rerun — stale, not wrong.** `git bisect`, 9 steps: the whole writer at the parent `4a1dc0d` reproduces the committed file (0 values differ), `9d42485` gives all 82, and `9d42485` equals HEAD. Three changes inside it: crossings recheck a cell's decision, the mean over a zero-length window is the value in force, and a same-type decision is not a revision. **exposure(cell) >= 1: 5 cells decided at birth, 170 fired** — 165 re-assertions, 0 change a type, 68 refused by the lock; with rechecks off, 5 and identical types. Regenerated through `save_result` with both counts side by side; 81 values changed; the shipped read's fates unchanged (522 / 483 / 928). **`celegans_commitment` followed PYTHONHASHSEED** through an unsorted set in `competence_test`; sorted, it is identical under seeds 0, 1, 2 and 12345 and reproduces neither committed value (9 competence values move in the fourth decimal, 0.0781 → 0.0774 bits; conclusions unchanged). A determinism test runs both writers under two hash seeds. Stale passages annotated beside the originals by the coordinator: this plan's area E, ORGANISM-FROM-ONE-CELL.md (the two `cell` rows and "it never does") and the embryo_factors.bio comment | — | done |
   | ↳ S6 | **Backup gap found 2026-10-01 (coordinator, read-only)** | git-ignored data | **No recoverable backup exists**: Time Machine has no destination configured (its last destination fails to mount) and the work OneDrive does not cover this folder. About 9.2 GB exists only on the owner's Mac: `data/knowledge` 6.5 GB (the AlphaGenome per-element answers about 1.4 GB, about 778,780 requests to recreate; the epigenome, panels, compiled programs, CRISPRi benchmark, contact and therapeutics stores), `data/reference` 1.4 GB, `data/cache` 945 MB, `data/individuals` 329 MB. Integrity today: of the 241 local inputs that committed results recorded, 239 match their sha256 and 0 have changed (the other two are descriptive labels, not files). **When the owner chooses a destination:** copy the model cache and the essential inputs first; `data/individuals` (a person's genome) only to encrypted local media, never a cloud service; verify every copied file's sha256 against the source, and restore a small sample to a scratch directory and re-hash it, since a copy alone does not establish recoverability | open (owner's choice) |
   | ↳ S6 | **Backup made and verified 2026-10-01 12:16 (coordinator, on the owner's authorisation)** | `/Volumes/Albert Canfield/GenomeOs/GenomeOS-snapshot-2026-10-01T1216` | The owner authorised an unencrypted local snapshot, `data/individuals` included, to an external USB APFS volume (238 GiB free; no earlier backups there to preserve). The whole project folder was copied from a clean tree at `8fcd37e`, git history included: **49,937 files, 10,376,834,634 bytes, 0 errors**, originals only read. Each file was hashed while copied and re-read; because that re-read could come from memory, **every file was then re-verified against the snapshot's `MANIFEST.sha256` with the OS cache disabled: 49,937 of 49,937 match, none missing**. **Restore test:** a fixed-seed sample of 105 files (1.50 GB: the ten largest, five or more from every group, `data/individuals`, the AlphaGenome answers, the other knowledge stores, reference, caches, `.git`, `.claude`, ignored results and tracked files, and fifty at random) read back from the drive uncached into a scratch folder: all match the manifest and the originals; the scratch copies, which included personal genomic data, were deleted. **Not backed up, by choice:** `.env` (the AlphaGenome and 4DN API keys), kept off an unencrypted drive pending the owner's decision; and regenerable files (`.venv`, test and lint caches, `build/`, bytecode, Finder metadata). The snapshot holds `MANIFEST.sha256` and `backup_report.json` at its root | done |
   | S7 | Release checks follow the revision being promoted | `.github/workflows/ci.yml:8` | CI runs on later pushes to a pull request; the scheduled check exercises the branch being worked on; checks are required on the revision promoted to `main`; local hooks stay, but not as the only enforcement | wave 1 |
   | ↳ S7 | **Done 2026-09-28 (lane-s7, `b50c873`)** | `.github/workflows/ci.yml`, `scripts/promote_main.sh` | `synchronize` added; a push to main runs CI (one run per promotion, a record, not a gate); the scheduled run asserts it is on `refs/heads/dev` (the checkout keeps `github.sha` so the recorded sha is the tested one). **Gate: `scripts/promote_main.sh`**, dry run by default, refuses a sha that is not on dev, not a fast-forward of main, or whose latest `test` run did not conclude `success` (it rejects `skipped` and `neutral`, which the branch rule would accept); 18 tests against a temporary remote. Dry run on the real repo: **no sha between main and dev has a green run**, so nothing can be promoted until the daily run is green. CONTRIBUTING's claims that every push to dev runs CI and that nothing reaches main unless CI is green were both false and are corrected. Recommended to the owner, not changed: turn on `enforce_admins` for main. Reported, not fixed (not a lane's file): the checkout guard refuses a push to `main` or `:main` but not `<sha>:refs/heads/main`; the old PR-path `scripts/promote.sh` still exists | the red test fixed by the coordinator the same hour | done |
   | ↳ S7 | **Gate found insufficient 2026-09-29 (coordinator)** | `scripts/promote_main.sh`, `CONTRIBUTING.md`, the checkout guard | Both recommendations closed first, as the owner ruled: **`enforce_admins` is on** (the owner's setting, verified by reading main's protection), and the guard now refuses every main-targeting push form and every `gh api` write to main's ref or protection (tested on 50 commands: the old hook 28 correct, the new one 50). Then `promote_main.sh --push 9a59faf` passed every check of the gate and **GitHub refused the push** (GH006, "Required status check \"test\" is expected"), although GitHub Actions had recorded a successful `test` check run on that exact sha (run 36550391987, `workflow_dispatch`, check suite on `dev`). **A successful manually triggered check is insufficient to establish GitHub eligibility**: the gate's claim that it lets through what the rule would is false, and is corrected beside it in the script and in CONTRIBUTING. Every earlier promotion went through as an admin, whom the rule did not bind. Why GitHub did not count the run is not established. Not bypassed. **Route now:** a branch frozen at the chosen sha (`promote-9a59faf`), a pull request (#11) whose checks run on it, merged by the owner; the coordinator verifies `main`'s file tree and ancestry against the chosen sha. Open: make the gate check what GitHub accepts, or prepare the pull-request branch itself; retire the old `scripts/promote.sh` | open |
   | ↳ S7 | **Cause narrowed 2026-09-29 (the external reviewer's correction, checked against the source by the coordinator)** | — | Beside "why GitHub did not count the run is not established" above, and the coordinator's guess that the run was recorded against `dev`: **the event type is the documented explanation.** GitHub's troubleshooting page (docs.github.com, "Troubleshooting required status checks", section "Checks from some workflow jobs are not evaluated") says checks from a workflow run are evaluated for a pull request, or for a branch ruleset's required checks, only when the run was triggered by `push`, `pull_request`, `pull_request_review`, `pull_request_target`, `deployment` or `deployment_status`, and names `workflow_dispatch` as not evaluated; `schedule` is not on the list. Run 36550391987 was `workflow_dispatch` (measured: event, sha and app read from the API). **Measured:** PR #11's `pull_request` run passed and the merge was accepted. **Inference:** the same event rule decided the refused direct push; the page covers pull requests and rulesets, not a classic rule's direct push. **Consequence for the gate (a proposal, not built):** a green run from `workflow_dispatch` or the daily `schedule` is our own pre-check only; the check that makes a sha eligible is the one on its pull request, so the gate should read that run's event, or stop at preparing the branch | open |
   | ↳ S7 | **Gate made honest 2026-09-30 (lane-gate, night shift; `3a893fc`)** | `scripts/promote_main.sh`, `scripts/promote.sh`, CONTRIBUTING.md | The gate prints the event of the `test` run it relies on: `pull_request` or `push` reads "eligible for protected promotion"; any other event, or a run it cannot find, reads "CI pre-check passed … a pre-check only", with the inference limitation kept wherever the event rule is cited. `--push` is refused before arguments are read; `--prepare SHA` makes the same checks and pushes only `promote-<sha7>` (refusing a branch at another sha) and prints the compare URL; `--verify SHA` checks main's tree and ancestry with git only; `scripts/promote.sh` exits 1 with a pointer. 18 existing and 22 new tests; they check the script's decisions against a local remote and do not simulate GitHub's enforcement. All additive: no `--force` (lane-s7's `b50c873` lines are held until 2026-09-30 23:03), so the old `--push` code stays, unreachable. **Clean-up after 23:03:** promote_main.sh lines 6–7, 54, 86, 243–246 and the stand-in refusal at 238–242; outdated text at 4, 13–14, 22–23, 216; CONTRIBUTING 87–88, 93, 99–100; the two "would run" test assertions. Not checked against the live API (no downloads) | open (clean-up) |
   | ↳ S7 | **Promoted 2026-10-02: `main` = `a1ce607` (PR #13, opened by the supervisor under Albert's delegation of 2026-10-01, merged by Albert; verified independently by the coordinator)** | — | The third promotion, and the first the coordinator did not prepare: Albert delegated opening PRs to main to the supervisor session in his own words on 2026-10-01, which reverses CONTRIBUTING's "sessions open no PRs" for that session only and is recorded in the commit-cycle memory. Candidate **`747b7ffb53d05a7d407fea7d6537535f3ccf2ff0`** was frozen on `promote-747b7ff`; the branch was updated with main (head `054eb46`, parents `747b7ff` and `f7016ff`, file tree identical to `747b7ff`); `test` and `live-readers` passed on the pull request. **Checked by the coordinator against the remote rather than taken from the report:** `origin/main` is `a1ce607` as claimed; `747b7ff` **is** an ancestor of it; the two file trees are **byte-identical** (`b0eaf2abcba23cedc744423efc5518b8eac87d8b` both sides); and `main`'s protection reads `enforce_admins: true` with `test` required (GitHub API, read after the merge). **Excluded and held for the next promotion:** everything after `747b7ff` — lane-cell2's eligibility work, the `bc28f79` estimator fix, lane-re2g, lane-s8 (`bee8f20`) and lane-rebuild. **Cadence set by the supervisor: at most one pull request to main a day, the next no earlier than 2026-10-03** | — | done |
   | ↳ S7 | **Qualifying event, not eligibility 2026-09-30 (lane-gate; `03497e7`, additions only) — incomplete** | `scripts/promote_main.sh`, CONTRIBUTING.md, tests | After the review: a `pull_request` or `push` run is a qualifying CI event only; the gate now prints, directly after its held line, that protected promotion eligibility is not established (the exact revision, required checks and protection rules decide; GitHub's decision is authoritative), and says the same on the pre-check branch, in a dated header comment and in CONTRIBUTING. 43 gate tests pass. **Incomplete:** the held line still prints "eligible for protected promotion" until its lines clear the guard (from `3a893fc`, 2026-10-01 22:39): then reword line 215, line 220, the comments at 35 and 40–41, the test at 305/309 and the check at 274, and CONTRIBUTING 116–120; the `b50c873` `--push` lines listed above clear on 2026-09-30 23:03. Direct pushes stay refused; promotion stays the owner's | open (blocked by the guard) |
   | ↳ S7 | **Gate finding, 2026-09-30 (narrowed by the external reviewer)** | `scripts/promote_main.sh` | The dry run refused `bdee236`: "does not contain origin/main (7ff4e37); moving main there is not a fast-forward". In the current workflow the gate rejects candidates that lack main's ancestry, which a pull-request merge commit on main produces; it does not necessarily reject every future candidate. No gate change requested. The promotion goes through a branch frozen at the candidate (`promote-bdee236`, tree `5d5973d6251dd78bbfc0eb7900ecaa48f946669b`), the owner's pull request, and a report of the updated head's and the synthetic merge commit's full tree hashes and required checks before merge | open (finding) |
   | ↳ S7 | **Promoted 2026-09-30: main = `f7016ff` (PR #12, merged 19:33 UTC)** | — | The candidate `bdee236` was frozen on `promote-bdee236`; the pull request was brought up to date with main (head `9c64371`, parents `bdee236` and `7ff4e37`), `test` and `live-readers` passed on it, and it was merged with every protection on. Verified afterwards: main `f7016ff42c42a3ca2237a9b1c6f163edbba3e0c4` has parents `7ff4e37` and `9c64371`; its file tree is `5d5973d6251dd78bbfc0eb7900ecaa48f946669b`, identical to `bdee236`'s (no file differs); `bdee236` is an ancestor of main; protection unchanged (admins bound, `test` required, up to date required, no force push, no deletion). Main now holds both judge versions (v2 direction, v3 target), the 2026-09-30 rebuild with the tenth therapeutic case, the evidence completeness report, the gate's qualifying-event wording, the clamp test, the CI fix, and the response map as API and CLI capability. **The response map's Evidence-tab page is not in it** and still awaits its separate permission. **Biological validation is unchanged**: nothing in this promotion validates a biological claim. Commits on dev after `bdee236` (from `155be1d`) wait for a later promotion | done |
   | ↳ S7 | **Promoted 2026-09-29: main = `7ff4e37` (PR #11, merged by the owner)** | — | The owner opened PR #11 from `promote-9a59faf` (frozen at `9a59faf`) and merged it at 12:48 UTC once `test` and `live-readers` passed on it, with every protection on. Verified by the coordinator afterwards: main's file tree is `f2b2af35`, identical to `9a59faf`'s (no file differs); `9a59faf` is a parent of main, beside the old main `45a5aa6`; protection unchanged (admins bound, `test` required, up to date required, no force push, no deletion). 232 commits reached main, the first promotion since `45a5aa6`. **It does not include the versioned judge** (`5931609` onward), promoted separately once its exact revision passes CI. The next promotion branch will not contain the merge commit, so under "up to date" the owner updates it first; the tree check stays the same | done |
   | S3 | Reopen the coherence objective with one bounded pilot | `attribution/joint_pretest.py:70` | R8's negative is kept and its two transformations retired, but its code and this plan stop claiming more than it showed; a small model with explicit roles, target links and cell context solves a synthetic case that needs two simultaneous corrections, then is scored on untouched evidence against unchanged labels and independent-block scoring; it returns ranked alternatives, supporting and conflicting evidence, unresolved alternatives and the next measurement that would separate them; **success is better prediction on untouched evidence at useful coverage, never a higher internal coherence score**; expansion stops if it adds nothing | wave 2, after wave 1 |
   | ↳ S3 | **Discontinued 2026-09-29 (lane-pilot; `a9660ff` registered, `1a492c8` gate 1, `a89bf5c` gate 2 registered, `197c560` result, `6bbf1d6` + `40b66e9` merge recount)** | `attribution/pilot.py`, `attribution/pilot_bio.py`, `attribution/joint_pretest.py` (R8's reading corrected beside the original; element competition and gene budget RETIRED) | **Gate 1 (synthetic) passed:** 303 of 303 planted pairs of simultaneous corrections found and committed; 101 of 101 indistinguishable pairs kept as one family with a separating measurement; the greedy control found 0 of 303. But the independent-block control also solved split-and-retarget and swap (101 of 101 each), and **only context hand-off (0 of 101) needed coupling**. **Gate 2 failed, `beats_unchanged_only`:** on chr1/6/8/9/10/11/19 the pilot beat the unchanged labels (Gasperini2019 +0.165 [0.074, 0.266], Morris +0.431, Xie +0.277, also at their coverage) but **passed 0 of 4 CRISPRi decrease endpoints**; minus independent-block, Gasperini2019 **−0.065 [−0.125, −0.012]**, the rest span 0; minus distance, above 0 on Xie only (+0.105 [0.015, 0.261]). Prior-only (distance, compiled target, H3K27ac) beat distance on Gasperini2019 (+0.193) and Xie (+0.086): descriptive, needs its own registration (H3K27ac selected the held-out file's positives). Metric: 6 validated, 6 errors, 249 untested of 261 committed, 118.34 of each per CPU-hour (0.0507 h), none tested by a held-out CRISPRi screen; 1,834 abstained changes. S4 at tiny coverage: Gasperini2019 target 36 of 44 (coverage 44 of 2,558). Defect found after the run and disclosed: evidence-free merges. An internal development result on withheld sources | — | **discontinued investigation (stop rule fired)** |
   | S4 | What a correct attribution means | — | origin, function, activity, target and context kept separate; not naming a coding gene is not absence of function; reporter activity, perturbation, contact and conservation never substitute for one another; "the observation model is inadequate" kept as an explanation; target accuracy reported apart from role accuracy and from coverage | wave 3 |
   | ↳ S4 | **Done 2026-09-29 (lane-s4, `b87f59e` registered, `ab62d99` result)** | `attribution/correctness.py`, `scripts/s4_correctness_run.py` | A correct attribution is judged claim by claim on five axes kept apart (origin, molecular role, activity, target, cell context). A table of 25 observation kinds × 5 axes says what each establishes, refutes, cannot read ("the observation model is inadequate"), only suggests, or cannot establish: a reporter establishes reporter activity and never a target or a direction in place; CRISPRi establishes target, direction and context and never a role or an origin; a well-powered null refutes one gene in one cell and never establishes absence of function; conservation, contact, association and the model output establish nothing. Target accuracy, role accuracy and coverage are three shares that refuse addition; every judging source passes `holdout.check_provenance`. **Applied once to the unchanged labels: origin 0 of 470,077 and molecular role 0 of 659,403 claims can be judged; target 98 of 440,377, all 98 established, but only 39 were refutable** (all K562, all held), because none of the 23 well-powered nulls lies in its rule's stated cell — so 98 of 98 is not evidence the targets are right; activity 93 of 98, **all 4 judged repressions wrong in direction**; context 39 of 39. 0 unresolved, 0 model-inadequate (the labels make no reporter claim). The ablation's 128 measured links split exactly into 98 / 23 / 7 — consistency with it, not new evidence. 46 tests. An internal development benchmark reading, not a validation | the pilot reports through `judge(claims, labels, sources=[S])` | done |
   | ↳ S4 | **Traced 2026-09-29 (lane-repress; `5d14d37` classes registered, `54c47f1` code, `0762073` tests, `3cc59bf` doc, `641909e` result)** | `scripts/repression_trace.py` | The four repressions S4 judged wrong, and the CCND1 comparison, followed from the model's cached answer to the deciding benchmark row. **No sign or convention error on either side.** Each compiled direction is the single most extreme of 371 tracks, and the model moves every one of these genes both ways (the winner by 0.0095 to 0.0711). **Cell mismatch is primary for CD83, HEMGN, BEX4 and CCND1**: the claim names the winning track's cell, the screen measured another, and for the first three the model's own value on the measured cell has the screen's sign. **ID1 is the one genuine contradiction, and a separate finding**: the model's cached K562 value (+0.1433) disagrees with the K562 screen (−0.158); that does not by itself refute the compiled whole-blood claim. The judge follows its registration: on activity `refutable` is never computed, and the registered rule judges direction from another cell when the stated one has no response, which is how all 5 wrong directions (and 54 correct ones) were decided. Nothing fixed, 0 labels or verdicts moved; `by_cell` keeping one track per cell name recorded, not fixed. 0 model requests. An internal development trace of five claims chosen for being wrong: not a rate, not a validation. **The owner's decision the same day:** the v1 result stands as history; a versioned rule judges direction only in the stated cell and returns "not assessed in this context" otherwise, symmetrically for all 59 cross-cell verdicts, with cross-cell disagreement kept as its own finding (lane-judge2, running) | done |
   | ↳ S4 | **Judge v2 done 2026-09-29 (lane-judge2; `5931609` registered, `172572c` code, `4567219` result, `8f3689e` doc)** | `attribution/correctness.py` (`RULE_V1`, `RULE_V2`), `scripts/s4_correctness_run.py --rule v2` | The owner's decision applied: the direction rule is versioned. **v1 is S4's rule unchanged** and reproduces `attribution_correctness.json` key by key, which stays as the historical result; every script that reproduces a v1 result pins it. **v2, now the default, judges a direction only in the stated cell**; otherwise `not_assessed_in_this_context`, with every other cell's response kept beside the verdict as a cross-cell finding. Re-judged once on S4's claims and sources: **activity 98 judged (93 / 5) → 39 judged (39 / 0), 59 not assessed in this context** (54 agreeing and 5 disagreeing cross-cell findings, the 5 being the four repressions and CCND1); target, context, origin and role unchanged; no unexplained count change. ID1 stays a separate finding. **A finding for the owner, not acted on:** 59 of the 98 established targets were decided only from another cell. Not regenerated, derived: 19 activity verdicts in the prior-only test's S4-beside block and 16 + 3 in the pilot gate's would move under v2. 0 model requests | done |
   | ↳ S4 | **The owner's rulings, 2026-09-29** | — | **The direction correction is resolved.** **Neither 39-of-39 figure is general accuracy**: each describes the small subset of claims assessable in the stated cell (direction under v2, target under v3). **Target, v3:** judge the compiled target claim only in its stated cell, for positive and negative verdicts alike; keep the 59 other-cell relationships as supporting evidence elsewhere; report 39 supported in context, 59 unassessed in context, 98 supported somewhere; version it and keep the earlier results (lane-judge3, running). **Main:** `9a59faf` is promoted as an interim verified snapshot once `enforce_admins` is verified on, through the protected path with no bypass on rejection; it does not include the versioned judge, which is promoted separately once its exact revision passes CI | running |
   | ↳ S4 | **Judge v3 done 2026-09-29 (lane-judge3; `2c95f66` registered, `b7cb814` code, `d1e09ae` result, `d9c9ea9` doc)** | `attribution/correctness.py` (`RULE_V3`, now the default), `scripts/s4_correctness_run.py --rule v3` | The owner's target ruling applied: a target claim is judged only in its stated cell, for supported and refuted verdicts alike, and every other-cell relationship is kept beside it as evidence elsewhere. Re-judged once on S4's claims and sources, from a clean worktree; v1 and v2 reproduce their committed files, which stand unchanged beside it. **Target: 39 supported in context, 0 refuted in context, 59 unassessed in context, 98 supported somewhere; direction as under v2 (39 / 0 / 59 unassessed / 93 supported somewhere); context 39 / 0.** **Every 39 of 39 describes only the small subset of claims assessable in the stated context, mostly K562, and none is general accuracy.** Exactly the 59 target claims moved from v2, each established only from another cell before and now carrying its agreeing other-cell responses as supporting evidence elsewhere; no count change unexplained; every count matched the registration. v1's target precondition (a response elsewhere turned a stated-cell null into a context error) is dropped, as registered. **Findings for the owner, not acted on:** one stated-cell null can now count against both target and context (0 cases in S4's claims), and the context axis keeps its own condition that the gene responds in some cell. Not regenerated, derived: the prior-only and pilot-gate S4-beside blocks would move as the doc lists (e.g. prior unchanged labels 43 of 43 → 24 supported in context, 19 unassessed). 0 model requests | done |
   | S5 | Runnable simulation separated from validated human kinetics | `attribution/bridge.py:174` | 22,580 of 440,589 rules (5.1%) simulable on mouse NIH 3T3 rates supports assumption-dependent simulation only; transferred rates labelled; parameter ranges propagated and the predictions stable across them identified; one human cell context with matched measurements preferred over more rules run on defaults | wave 3 |
   | ↳ S5 | **Done 2026-09-29 (lane-s5, `0ea0a90` registered, `83013e6` label and pre-run amendment, `a54ea12` result, `74e6a2f` audit)** | `attribution/bridge.py`, `runtime/grn.py`, `scripts/s5_rate_ranges.py`, `scripts/bridge_audit.py` | Every output built on a transferred rate says **assumption-dependent simulation, not measured human dynamics** and names, per gene and per trajectory, its source, the species and cell it was measured in, its tier and the context it runs in; a run left on the default half-life names it (mixed units). **22,576 of 22,576 pairs labelled**; labelled = simulable in every tier, no count moved. Ranges read from measured spreads: replicate disagreement 2^±1.12; K562 over NIH 3T3 half-life 0.077–0.639 (2,262 genes); the rate bracketed by rate-carries-over against copy-number-carries-over; the combination rule applies to no simulable pair. Runtime equals the closed form at 13,000 of 13,000 sampled points. **Stable over the range: direction and fold for all (they are the input); level, effect size and response time for none; on/off for 16,101 (71.3%, all on); switch for 16,295 (all "not a switch"); both for 14,683 of 22,580 rules; rank order for 0.85% of level comparisons.** Stable is not validated. One registered expectation failed (10 rank comparisons under R3). The one human context is **K562** (4 of 7 measurement kinds, 347 simulable pairs, 199 with a K562 half-life); a validated kinetics test there lacks a readout that is not an input, a time course after an acute perturbation, and an absolute rate or copies per cell. Outside the lane: `runtime/abundance_gate.py` uses a borrowed mouse protein half-life the label does not cover | K562 time course after acute perturbation; absolute K562 calibration | done |
   | S8 | The next experiment's information value | — | a representative arm (how common a function is) and a disagreement-selected arm (which model is right) kept apart, selection probabilities recorded; independent loci, assay feasibility and distinguishable hypotheses before volume; HCT116 stays inconclusive, and re-reading held-out data is not fresh validation | wave 3 |
   | ↳ S8 | **Built 2026-10-01, pushed `bee8f20` (lane-s8): the gate that asks first, checked by the coordinator against the files** | `genomeos/attribution/experiment_value.py`, `tests/test_experiment_value.py` | Four gates in a fixed order; a design is refused at the first it fails, always with the reason and **never with a sample size**. **(1) A blinded eligibility count:** the unit of analysis named, then units, label prevalence and independent loci. Blinded by construction — the records counted carry no field for a prediction, a score or a hypothesis, and the gate refuses outright if handed records that do. This is the check N1, the wiring diagnostic, lane-contrast, lane-entex and lane-indep each lacked. **(2) Assay feasibility:** an unavailable assay drops a candidate, while an unstated detection limit or an **undeclared readout comparability** refuses the whole design — the paired-enhancer survey's failure was comparability never declared. **(3) Distinguishability:** each hypothesis's predicted outcome in the assay's own units beside the detection limit; two hypotheses predicting the same outcome are refused as indistinguishable whatever the assay can see, and a separation below the detection limit is refused because the assay could not report it. **(4) Freshness:** a candidate already read is development evidence for good, a fresh-validation claim is refused when any surviving candidate has been read, and only a context whose registered reading is `replicated` counts as a replication, so an inconclusive second cell type stays inconclusive and its registered reading travels word for word. The **representative** arm (how common a function is) and the **disagreement-selected** arm (which model is right) are reported apart, each selection carrying the probability that selected it; pooling them raises. **Volume comes last** and says of itself that it is the registered target allocated over the eligible independent loci, **not a power calculation** — there is no variance model and no minimum detectable effect. The **independent-locus convention** is in code with its reason and its limits: same target gene, or within 1 Mb on one chromosome, chained; an operational grouping for counting, **not established biological independence**, and a count under it an upper bound on the independent evidence a design holds. **Verified by the coordinator, not taken from the lane's report:** 46 tests pass; **0** dataset or cell-line names in the module; **0** default values in `Floors`, so every floor is the design's own and registered with it and none can be moved after an outcome. The lane also reports nine deliberate mutations each caught by a test. No experiment selected, no dataset read, 0 model requests, no money | lane-s8 | done |
   | ↳ S8 | **Four things S8 deliberately left to a decision, 2026-10-01 (lane-s8, recorded by the coordinator; none chosen by the lane)** | (1) **The floor values.** `Floors` ships no defaults: a real design must fix `min_units`, `min_positive_units`, `min_label_prevalence`, `min_independent_loci`, `min_distinguishing_candidates`, `target_loci` and `units_per_locus` before any count. The lane invented no number, which is right; **the numbers are the owner's**. (2) **One shared locus constant.** `experiment_value.py` and `attribution/cell2.py` registered the same 1 Mb chained convention on the same day (`7e45429`), independently and with a comment in each saying the two must not diverge. **Coordinator's decision: they become one constant, after lane-cell2 finishes**, because that module is live and a shared file must not be edited under a working lane; until then the comment is the guard. (3) **Where an incomplete design fails.** An undeclared readout comparability or an unstated detection limit refuses the whole design, while an unavailable assay drops only the candidate; the lane's reason is in the module and the alternative is a two-line change plus a test. (4) **Volume is an allocation, not power.** A detectable-effect calculation would need a variance assumption registered per assay, with data and money implications, so the lane did not start it. Also left undone: no design is instantiated — no script, CLI, API surface or result file — and `independent_loci` treats positions as points, not intervals, so wide elements need checking before they are passed to it | — | owner decides (1); coordinator holds (2) |

    **The 2026-09-30 rebuild: input list and acceptance checks** (prepared 2026-09-29 by the
    coordinator from a read-only inventory, spot-checked; nothing run). The removal guard
    (`scripts/check_staged.py`, `--since 2.days` on committer dates: exactly 48 hours, per path) held
    each item below; each clears at the time given, BST on 2026-09-30. One lane after 04:22, one commit
    per item. **Expected values are regression expectations, not targets: any figure outside an item's
    acceptance stops the lane for explanation.**

    1. **Two results under their own names** (clears 00:24). `scripts/measured_layer.py --write-programs`,
       then `scripts/confidence_calibration.py --retired-formula` (since `8d57243` a bare run writes
       `model_score_curve_genome`). Accept: `measured_layer_genome` equals `measured_layer_genome_r6` except
       `evidence` (weak 928,094 → 14,088, unstated 914,006), date, seconds and manifest; lentiMPRA 3,231
       agree, 14,512 disagree, 126 conflicts; CRISPRi 97 and 23; `programs_written` 25 entries; calibration
       n 17,743, observed 0.1821, ECE 0.0963. Not yet confirmed that `--retired-formula` reproduces `_r6`
       exactly; the `_r6` files stay until it does. No test pins either file.
    2. **The chr21 program** `data/organisms/human/noncoding_chr21.bio` (clears 03:04; path clear 04:22),
       written by the same run. Accept: 5,175 rule lines change (5,174 `confidence:` lines from `25f81c6`
       and the held ICOSLG rule at L79999 from `36fa298`), 1 basis line (SUMO3 at L80225, `36fa298`), and
       the header (7 lines out, 1 in: L32 and L39–43 from `96bc5e5` go with L33, `d2ab0ae`'s correction);
       0 role or axis lines change; `# test: rules == 5176` holds. Never copy
       `data/knowledge/compiled/noncoding_chr21.bio` over it (older than `690a71c`, pre-R7 role text). Tests:
       `test_measured_layer.py::test_the_committed_chromosome_program_states_experimental_facts`, `bio test`.
    3. **`attribution/compile.py`'s overridden lines** (clears 03:27): L60–67 (the first `MODEL_SCORE_NAME`
       and its correction) and L377–381 (`score = pc.get("confidence")` and the lines overriding it) become
       one definition each. The header at L609 still says "capped at 0.7" and a recompile alone does not
       change it, so it is rewritten here. Accept: compiled outputs identical except that header;
       `test_compile_certainty.py`, `test_r4f_certainty.py` green.
    4. **`attribution/variation.py`'s two "best guess" strings** (L73, L104–105; clears 00:58). Wording
       only. Accept: `variation.build` shows 0 numeric differences.
    5. **The grammar bullet pair** (`lang/grammar.py` L483–490, regenerated into
       `docs/BIOLANG-GRAMMAR.md` L17–24; clears 02:06): keep the "Current behaviour" bullet and "matches no".
       Tests: `test_grammar.py::test_grammar_document_is_current`,
       `test_rule_context.py::test_the_grammar_says_what_unknown_means_in_a_when_clause`,
       `test_engine_package.py`. Held by the same guard; not confirmed it was on the original list.
    6. **Two strict xfails** in `tests/test_crispri_split.py` (L122–126, L137–141; clears 01:01), kept
       "until the removal window lets it be rewritten". Accept: rewritten to assert the corrected behaviour
       and passing; no assertion loosened.
    7. **The tenth therapeutic case into the headline result** (clears 00:59; `manifest_headlines.json`
       03:54). No usable copy remains, so it is a fresh run at a clean revision
       (`scripts/manifest_rebuild.py data/results/therapeutic_benchmark.json --where DIR`; public network
       sources, no paid requests). In the same commit: `tests/test_manifest_headlines.py:43` (9, 9, 9) →
       (10, 10, 10), `manifest_headlines.json` `rebuilt_at`, README.md:122 "9/9", docs/THERAPEUTICS.md:1057.
       Accept: rows 1–9 identical except `seconds`; 10, 10, 10; surface 5 of 5; copy number 2; ERBB2 rank 1
       at 0.494; MYC 0.246 at 8 copies; `gate_moved_the_target` and `magnitude_tiebreak_changed_order` empty;
       the allele-fraction list `["ERBB2"]`; the four tests at `test_therapeutic_benchmark.py` L432–518,
       skipped until the case is present, run and pass.
    8. **Not in this rebuild, each its own lane:** regenerating `unknown_<chrom>` (24 runs) with
       `compile.py` reading the classifier's origin (origin text on exactly 11 blocks; 0 numeric
       differences; clears 02:59); the R7 legacy-key rename (no guard hold, a schema change); and checking
       whether `measured_layer_genome_r4` and `evidence_compiled_genome_r4` (`37172b6`) were held too.

    Around the whole rebuild: `scripts/check.sh` green before and after; CI on the exact revision; a
    promotion only on the owner's specific authorisation, through a pull request.

    **Result, 2026-09-30 05:12 BST (lane-rebuild; `1a11a7b`, `d13e50b`, `e801b6b`, `e9f8526`, `1163e4d`,
    `5ab1b95`, `aa7f92c`; checked by the coordinator against the result files).** Items 1–7 done, no guard
    refusal, no override; `check.sh` green before (2,704 passed) and after (2,710; bio 221/221).
    Rebuilt and verified: (1) `measured_layer_genome` equals `_r6` except evidence (weak 14,088 stated at or
    below the line, 914,006 unstated counted apart, the old 928,094 kept beside), date, seconds, name and
    code block; lentiMPRA 3,231 / 14,512 / 126, CRISPRi 97 / 23, 25 programs; `--retired-formula`
    reproduces `_r6`: n 17,743, observed 0.1821, ECE 0.0963; clean stamps; the `_r6` files kept. (2) The
    chr21 program's 5,825 `confidence:` matches reconcile as 5,174 predicted rules + 1 measured + ICOSLG +
    446 region fields + 200 measured fields + 3 header comments; 5,175 rule lines, 1 basis line and the
    header change, no role or axis line; `rules == 5176`. (3) one model-score name and one score line;
    every recompiled program changes by exactly one header line. (4) the variation strings: 0 numeric
    differences over 14,288 leaves. (5) the grammar pair, not refused. (6) the two strict xfails rewritten,
    no assertion loosened. (7) the tenth therapeutic case: 10 / 10 / 10, rows 1–9 identical in all 288 shared
    fields except `seconds`, the four tests pass; README and the pins moved in the same commit.
    Legitimate changes, explained: items committed in the order 3, 5, 2, 1, 4, 6, 7 so the chr21 header was
    written once in item 3's wording; the chr21 generation date changes; item 7 was run by its writer at a
    clean revision because `manifest_rebuild.py` checks out `65d2d23`, before the tenth case, so rows 1–9
    gain the three fields `64f0afe`'s code writes, and the local therapeutics store differs from the bytes
    the old manifest recorded; `manifest_headlines.json` also moved `value` and `quoted_as` to ten. Not
    done: CI on `aa7f92c`; the `_r6` files; item 8.

    **Reconciliation, 2026-09-30 08:30 (the reviewer's questions, answered from existing records).** Item 7
    is an **updated result, not a reproduction**, with complete provenance. Its manifest is complete (clean
    stamp at `5ab1b95`) and names 24 inputs against the old 22: two added, the tenth case's hand-written
    `data/demo/benchmark/erbb2_myc_gastroesophageal.vcf` and `.cnv`; one changed, the local store
    `data/knowledge/therapeutics` (sha256 `12bf0ffff535…`, 48,840,306 bytes, 977 files → `fbda4e49b635…`,
    48,853,641 bytes, 984 files). The change is additions only, proved by hashing today's store without the
    seven files written on 2026-09-28 09:26 (Open Targets and HPA records for STARD3, GRB7, MIEN1, MYC, the
    tenth case's genes): it gives exactly the old digest. Both input sets are available on this machine (the
    store is git-ignored, so on no other). Cases 1–9: no scientific value changed; three fields were added
    (`magnitude_tiebreaks`, `quantities_in_this_tumour`, `rank_without_gate`, written by `64f0afe`'s code)
    and `seconds` changed. The 10 / 10 / 10 counts ten benchmark cases (target gene recovered, modality
    verdict correct, top mechanism defensible), not ten validated therapies. **CI:** the same `pytest -q` runs
    in both places. CI collected 2,706 tests, as did the rebuild's clean-worktree pre-push run (2,695 passed,
    6 skipped, 5 xfailed); the local 2,721 included 15 tests of the then-uncommitted
    `tests/test_response_map.py`. CI skips 91 more than a clean local run: 90 need git-ignored local data or a
    key (reference sequences 25 + 2, HG002 VCFs 21, VEP caches 12, benchmark and compiled files, JASPAR, the
    worm's inputs) and 1 is the history test, skipped only in a shallow checkout; nothing is deselected. What
    that skip gives up in CI is the check that the rule version the discovery view prints, `b7e4bf0`,
    resolves to a commit here; the check still runs in every full checkout and in the pre-push hook, whose
    temporary worktree shares the full history. A fixture repository could test the lookup but not that
    fact about this repository's history, so none was added.
    *Accepted by the external reviewer, with a portability limitation: rerunning the therapeutic result
    currently depends on inputs that exist only on the owner's machine (the git-ignored
    `data/knowledge/therapeutics` store). CI passes the tests it can run; local and pre-push checks cover the
    additional data-dependent behaviour; the one shallow-history skip has the narrow consequence stated above
    and keeps its coverage in full checkouts.*

   Delivery: wave 1 first; then one bounded pilot on cached evidence, ordinary CPUs and local regulatory
   neighbourhoods, chromosomes as validation partitions, supported connections kept across neighbourhood
   boundaries; genome-wide deployment only after the pilot shows improvement.

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
   | C4 | **Hide evidence, predict it** | Hold out one assay or study, infer labels from the rest, predict what was held out; rotate. Split by locus and study, track shared provenance (R5), model each assay's endpoint separately (R2, R6). Plus the diagnostic run first: **remove AlphaGenome predictions entirely — which conclusions survive on experimental evidence alone?** An internal development benchmark only: evidence already used repeatedly cannot become a fresh external validation set, and it is labelled so | cached data; R5 provenance | **A, now: lane-c4** — it is the measuring stick every other approach is scored by |
   | ↳ C4 | **Done 2026-09-28 (lane-c4, `cc824f4` registered, `b7e4bf0` result)** | `attribution/holdout.py`, `attribution/ablation.py` | **Ablation:** without AlphaGenome, **93 of 440,377 predicted links survive** — a CRISPRi pair on the same element and gene establishes the stated direction; 28 are contradicted (23 by a well-powered null, 5 by the opposite sign), 7 inconclusive, and 440,249 never measured by the same endpoint (438,872 elements never screened). 39 predicted rules survive in their own cell. The 212 experimental rules and all sequence and registry labels were never model-dependent. Of 6 headline results, node containment survives only on its measured arm (+5.893 [3.183, 8.458]), 3 were never model-dependent, constrained_unknown_targets is inconclusive and crispri_published is about the model. The model's element selection kept 449 of 661 measured CRISPRi decrease links out of the programs. **Harness:** `holdout.score(labels, source)` and `compare`, over 16 sources and 15 scored endpoints, locus bootstrap, masked by study and locus, the CRISPRi held-out file never evidence (`LeakError` otherwise); about 52 s for all 16 sources on one CPU. **First scores, negatives first:** the unfitted `rest` rule is below distance to TSS on every scorable CRISPRi decrease endpoint and on GTEx; **the unchanged labels beat distance on none of the four CRISPRi decrease endpoints** (Morris −0.366 [−0.522, −0.187]); lentiMPRA, VISTA and saturation mutagenesis are predicted near chance by every labelling; **14 of 20 CRISPRi endpoints are unscorable**. An internal development benchmark reading only, never a validation | the pilot scores against unchanged labels, independent-block scoring **and distance to TSS** | done |
   | ↳ C4 | **Correction 2026-09-29, beside the row above (the external reviewer's reading, checked against the result)** | `attribution/ablation.py` `placement()` | "The model's element selection kept 449 of 661 measured CRISPRi decrease links out of the programs" overstates. `c4_alphagenome_ablation.json` `placement` says: **212 placed** in a compiled element (169 training, 43 held-out); **62 not placed that an ENCODE registry element would carry** (43, 19): only these are the model's selection; **387 that no registry element meets under the overlap rule** (259, 128), which no selection among registry elements could have placed. The 449 are placement failures, not 449 proven absent enhancers, and why the 387 fail (build or coordinate convention, interval boundaries, the overlap rule, candidate exclusion) is not established. Audit A (lane-place, running) measures it on all screened pairs, nulls included, with ranking frozen | open |
   | C5 | **Natural variation as experiments** | Different haplotypes are alternative versions of regulatory sequence; where genotype and a molecular measurement are paired (allele-specific expression, splicing, chromatin), ask whether an explanation predicts the allelic difference. Mapping bias, imprinting and linked variants handled explicitly. VCFs alone supply no outcome | a paired dataset, if one is reachable | **A, now: lane-c5 probe** (what paired data exists, its licence, size, and bias controls); build in phase C |
   | ↳ C5 | **Natural variation probed 2026-09-28 (lane-c5, `3156e5b`): open paired allelic data reaches about a fifth of the attribution layer** | `scripts/c5_paired_variation_probe.py`, DATA.md "Natural variation as experiments: the C5 probe" | EN-TEx (4 donors, reads mapped to each donor's own phased genome, allele counts provided, open) reaches 207,479 of 961,227 elements and 85,571 of 440,377 pairs in the same donor and tissue. Geuvadis chr21 (449 samples) has 3,608 of 5,174 pairs contrastable, 395 in the LCL context. GTEx phASER is open but cannot be paired: its genotypes are controlled-access (not requested). **Linked variants decide what a result can say**: a median 1,580 to 2,012 heterozygous SNVs lie within ±1 Mb of a target's TSS, so a within-person imbalance cannot be credited to one element — EN-TEx supports a coordination test, not an attribution, and the between-individual Geuvadis contrast is the design that can address linkage. "Reached" means reach, not power: read depth sits in columns not read. Imprinted and immunoglobulin/HLA targets counted for exclusion; one binomial test needs 47 reads to tell 0.70 from 0.50 and 783 for 0.55. Counts only; no allelic outcome read; 0 requests | phase C only if the pilot helps: EN-TEx within-donor coordination first, Geuvadis between-individual second | probe done |
   | ↳ A–C | **Three bounded audits started 2026-09-29 (proposed by the external reviewer, run as coordinator lanes; not a continuation of the discontinued pilot)** | **A, candidate coverage** (lane-place): why 449 of 661 measured CRISPRi decrease links are not placed, by cause (build and coordinates, boundaries, overlap rule, candidate exclusion), on all screened pairs with nulls, ranking frozen, at most one outcome-blind policy registered before it is scored; stops if a change only moves the denominator or inflates overlap. **B, cached context contrast** (lane-contrast): feasibility of one feature, same-cell signed deletion minus the median other-cell value, given that `by_cell` keeps the last same-named track; counts and a duplication check first, a registration with a frozen threshold (proposed +0.02 AP, paired locus interval above zero, no coverage loss) only on a clear go, no scoring in this lane. **C, EN-TEx allele coordination** (lane-entex): go or no-go on same-donor ATAC and RNA allelic counts, phasing, depth, exclusions and exposure, with power from total depth only and no signed outcome opened; a coordination test, not an attribution; independence from GTEx unresolved. All three: no model requests, no paid services, previously read outcomes are development evidence, a negative finding completes the checkpoint | C4, S4 harness; C5 probe | running |
   | ↳ A | **Checkpoint 2026-09-29 (lane-place; `dd8c49c` policy registered, `0183ef5` code, `6c0d39f` result)** | **No coordinate bug; the placement rule, not the model, loses most measured links.** All inputs are GRCh38, 0-based half-open (sha256 in the result); all 440,377 compiled elements sit at their cCRE's coordinates; tested intervals touch a cCRE 91.2% of the time against about 42% when moved 20 kb. C4's 212 / 62 / 387 reproduce exactly. Of C4's 387 "no registry element" links, **365 do overlap a cCRE** and fail the reciprocal-0.5 overlap rule on width (290 tested intervals are wider than 700 bp, where no 150–350 bp cCRE can meet it; 72 more fail at a reachable width), 3 overlap only partly, **21 are absent from both element sets**, and 1 is a zero-width source row. **The current rule places positives and nulls at the same rate (32.1% vs 32.8%).** One outcome-blind policy, registered before it was scored (overlap ≥ half the smaller interval): placement 75.2% of positives vs 54.8% of nulls; rescues 285 of 449 positives vs 3,093 of 9,460 nulls, excess +0.308 [0.258, 0.356], +0.26 [0.214, 0.305] over 20 kb-shifted intervals, +0.140 [0.068, 0.215] within 75–700 bp. No registered stop fired, **but 235 of the 285 rescues come from intervals over 700 bp**, where half the placements are carried by more than one element and about 28% of positive placements would also occur 20 kb away. Distance-to-TSS AUROC 0.873–0.880 on every subset. **Not applied anywhere.** Internal development evidence: these outcomes were read before (seven commits listed in the result). 9.5 s, 1.14 GB, 0 downloads, 0 requests. **Open, not measured:** C4's survival census (`CrispriIndex.of`) uses the same reciprocal rule, so its "93 of 440,377 supported" likely has the same width blind spot | C4 correction above | checkpoint |
   | ↳ A | **Correction 2026-09-29, beside the row above (the external reviewer's reading, checked against `scripts/placement_audit.py`)** | Three things in the coordinator's row say more than the audit does. (1) **The breakdown of the 387:** 362 fail the reciprocal rule on width (290 wider than 700 bp, 72 at a reachable width), 3 overlap only partly, 21 are absent from both element sets, 1 is a zero-width row; 365 (362 + 3) overlap a cCRE, and 365 is not the number of width failures. (2) **"Nulls" are a non-decrease comparator**, not evidence of no regulation: the audit counts every pair that is not a significant decrease, so its 14,073 include 4,113 underpowered observations and 159 increases beside 9,801 well-powered nulls. (3) **Equal placement rates (32.1% vs 32.8%) do not establish that placement carries no signal** about real links. The registered result stands as written; a sensitivity restricted to well-powered nulls would be additional development analysis. Next, bounded and read-only: a census of C4's evidence attachment under the current and the registered interval policy, recorded apart, replacing no headline and no production matching (lane-census) | checkpoint 2026-09-29 | open |
   | ↳ A | **Census 2026-09-29 (lane-census; `1253d69` registered, `32af2cb` result)** | Read-only, recorded apart, stopped after one pass: C4's CRISPRi evidence attachment counted **by observation** under the current reciprocal-0.5 rule and under audit A's registered policy, each response of a wide perturbation held as one observation with a candidate set, never multiplied into several validations. Denominators unchanged and checked (440,377 predicted links, 14,734 screened observations); C4's committed census and audit A's placements reproduce. **Current rule: 93 supported, 23 refuted by a well-powered null, 5 by the opposite sign, 7 unassessed (underpowered); all 128 attached observations are unique-element assignments and none of the 3,551 intervals over 700 bp attaches.** **Policy: 208 supported, 50 and 16 refuted, 18 split in direction, 16 unassessed; of 308 attached observations 163 are unique and 145 ambiguous; of the 208 supported, 120 unique and 88 ambiguous.** The policy adds 180 observations, 132 of them ambiguous and 152 from intervals over 700 bp, where only 26 of 152 attach uniquely; Gasperini2019 gives 115 of the 180; Klann and Morris attach only ambiguously. Per link, 93 / 28 / 7 would become 261 / 100 / 19, but the 261 rest on 226 observations and 141 of them only on ambiguous ones. Nothing adopted, no matcher changed, C4 not regenerated; adopting a matching change waits on review of what a set-valued observation can establish. Internal development evidence. 4.9 s, 959 MB, 0 requests | audit A | checkpoint |
   | ↳ A | **Scope of the census counts, beside the row above (the external reviewer's reading, checked)** | **The census reproduces C4's any-cell attachment**: a predicted link stores no claim cell and the observation is matched on gene, so 93 and 208 are historical attachment and direction-agreement counts, **not support or refutation in the stated context under the v3 judge**, which is the evaluator for any context-specific claim. **"Unique element" means one candidate among the cCRE registry elements the overlap rule admits**, not an experimentally resolved single cause. Transitions, recomputed in memory with nothing written: 116 observations newly supported under the policy (37 unique, 79 ambiguous), 92 stay supported (83 unique, 9 ambiguous), 1 supported becomes split in direction; so +115 and the +27 on unique support are net changes | census | note |
   | ↳ A | **Discovery-only review started 2026-09-29 (lane-discover), under amended rules** | All 14,734 observations under both rules with their IDs; the changed set reported continuously (tested-interval coverage from the union of candidate overlaps, per-candidate coverage, candidate multiplicity; no 700 bp boundary); discovery screened by 20 kb-shifted enrichment as a diagnostic, not a false-attribution estimate; **verdict eligibility stays under the current rule**, broader eligibility undecided without a defensible observation model; a previously unique assignment that becomes ambiguous keeps its historical verdict under the old rule and is unresolved under the new one; the 24 / 13 stated-cell split recorded reproducibly (provisional until then); wording "all admitted registry candidates"; no production matcher change; stops after its registered report | census | running |
   | ↳ A | **Discovery-only review done 2026-09-29 (lane-discover; `a9ab0db` registered and code, `66051ee` result)** | Reading `merits_a_separately_reviewed_proposal_only_in_the_passing_strata`; nothing implemented, no matcher changed, historical verdicts kept, verdict eligibility under the current rule. **Discoverable:** 184 observations (180 newly attached, 4 with added links), 48 unique among all admitted registry candidates and 136 ambiguous (median 2 candidates, up to 12); union coverage of the tested interval median 0.42 for the unique (58% of the median interval in no admitted candidate), 0.67 for the ambiguous; 5,995 more gain candidates but stay unattached; 0 lost. **Screen (diagnostic):** pooled 184 against 76 at ±20 kb, 2.42×; 8 strata pass, 7 of them on 16 or fewer observations; WTC11_DC_TAP/WTC11 fails at 0.86; Nasser2021/GM12878 undefined. **Ambiguity introduced:** 136 of the 184, and 13 previously unique assignments become ambiguous, each keeping its historical verdict (10 supported, 2 refuted by a well-powered null, 1 underpowered) and marked unresolved under the new rule. **Same cell vs other cell:** the 24 / 13 split of the 37 newly supported unique observations, 39 of 79 ambiguous with a stated-cell match and the 39 / 54 check now recorded reproducibly; of the 184, 74 in the stated cell and 110 in another cell only. Internal development evidence; the screen is a point ratio with no interval; the counts are C4's any-cell attachment, not stated-context support. 10.5 s, 1.06 GB, 0 requests; the result file is 7.2 MB (a compaction needing `--force` was refused by the permission check and not made) | census | done |
   | ↳ A | **Evidence view, read-only, approved on conditions 2026-09-29 (lane-evview, running)** | The external reviewer accepted the discovery review as a completed audit and approved one bounded addition to the existing evidence view, `/api/evidence/discovery`, reading the committed `discovery_review.json` by pinned sha256 and recomputing nothing. One record per observation, referenced by id from each section: tested interval, measured gene, assay cell, the admitted registry candidates with union coverage, the attached claims each with its context (match, other cell, no stated cell), and a mutually exclusive row category (all attached claims match, some match, none match, no assessable claim context, no attached claims). Sections: the discovery table for Gasperini2019/K562 only (116 of 5,299 screened observations, both read from the file; 41 all-match, 9 some, 66 none); all 13 unique-to-ambiguous transitions across studies (3 all-match, 10 none), each with its historical verdict and rule version; the other strata as labelled counts only. Fixed wording: "The stated-cell subset is selected using model predictions and is not representative of all screened pairs"; every verdict marked as the historical C4 any-cell reading, not a v3 stated-context verdict; "the discovery view adds access to evidence, not new validation". **Why the 37 discovered nulls are all other-cell-only, provisional (in memory, per link, not yet recorded):** the stated cell is the model's strongest track; for agreeing links the model's K562 effect is near its strongest (median share 0.876, 70 of 143 links state K562), for null links it is small (0.112, 0 of 49), so the stated-cell subset is selected by the prediction; source selection and gene-only matching are not shown to contribute nothing. No matcher, scoring, label or result change | discovery review | running |
   | ↳ A | **Evidence view accepted 2026-09-29 at `7ae0681` (the external reviewer)** | `/api/evidence/discovery` and one section of the Evidence tab serve 128 observation records (116 in the Gasperini2019/K562 discovery table, 13 resolution changes, one shared), reconciled with the committed `discovery_review.json` by its sha256; 20 targeted tests and the 60 web and evidence tests pass; the v1, v2 and v3 result files and the existing `/api/evidence` output are unchanged. The lane that wrote it was stopped before committing; the coordinator reviewed it against the specification, took two screenshots of the real page (the mixed-context and none-match rows; all 13 resolution changes with their historical verdicts), and committed it as written. **Visual acceptance covers the supplied cropped views, not the complete page layout.** It improves access to evidence and its uncertainty and adds no biological validation. Closed; no further changes or reruns requested | discovery review | done |
   | ↳ MAP | **Cellular control and response map, increment 1 started 2026-09-30 (the reviewer's scoped task; lane-map)** | Architecture mapped first (read-only): intervals (`coords.Locus`, no assembly field; GRCh38 and half-open recorded in result manifests), elements and their five axes (`ir/model.py`, `grammar.AXES`), genes and proteins (`ir/model.py`, committed `lib/data/proteome.json.gz` with UniProt function and Reactome), cell context (`when: cell_type`), perturbations (`measured.CrispriPair`, five outcomes), relations (six rule actions; contact in `genome/hic_contact.py`, local cache only, not a loaded observation), evidence and verdicts (`evidence.py`, `correctness.py` v1–v3, `Certainty`, unresolved alternatives), runtime (`grn.py`); `decompile` reads no CRISPRi, contact or measured layer; no human time course, dose response or persistence data exists locally. **Example: the β-globin locus in K562**, chosen for traceability: committed CRISPRi decreases (Reilly HBG1 −0.7356 and −0.6908 in `crispri_direction.json`, where the model named HBE1 instead; Schraivogel2020 and Nasser2021 HBE1 and HBG2 decreases in `discovery_review.json`, one interval lowering two genes; many underpowered observations), a context mismatch (v3 holds `EH38E2941908` → HBE1 stated in HT1080 as not assessed in this context, the K562 decrease kept as other-cell evidence), and protein function from the committed proteome; these outcomes were examined before. **Increment:** one read-only view that assembles typed assertions by reference (entity identity, functional participation and evidence status kept apart; relations encodes, binds, physically contacts, changes accessibility, changes measured RNA or protein, participates in a process, associated with a state), shows where the explanation stops, and marks dynamics questions not assessable; no schema, compiler, verdict, label or result change | architecture review | running |
   | ↳ MAP | **Correction 2026-09-30, beside the row above (found by lane-map against the committed sources)** | Two statements in the row above were wrong. For the proximal Reilly HBG1 interval (chr11:5253147-5253547) the model named **HBG2**, not HBE1; it named HBE1 only for the distal interval (5275847-5276247). And the committed proteome gives HBG1, HBG2 and HBE1 only Reactome R-HSA-983231 (megakaryocyte and platelet); **no committed annotation says oxygen transport** for them, which the map therefore shows as a stop point | increment 1 | note |
   | ↳ MAP | **Increment 1 done 2026-09-30 (lane-map; `ed734fe`), the Evidence-tab section blocked by the permission system** | `genomeos/response_map.py` (read-only builder, a rule check that refuses a payload breaking the safeguards, cycles without order, a CLI `python -m genomeos.response_map globin_k562 [--summary]`) and `/api/evidence/response-map?example=globin_k562`; reuses `evidence_discovery`, `correctness` (cell comparison, `NOT_ASSESSED`), `measured` outcomes and `hic_contact`; reads v3 verdicts from the committed file and evaluates nothing. 107 assertions, each with its source, record key and file sha256: 22 observed (CRISPRi effects in K562, the HBB_LCR DNase and lentiMPRA readings, 10 K562 Hi-C contacts from the local cache), 14 predicted (8 compiled claims, not assessed in their stated cells, the K562 decreases kept as other-cell evidence), 26 inferred (proteome annotations, not measured in K562; the erythrocyte program's adult HbA rules, labelled not K562), 45 unknown (underpowered observations, not negative). Measured effects keep the relation unresolved between RNA and protein because no committed record says which was read out; 12 stop points; every dynamics question not assessable with current data. 16 tests; the evidence and web suites pass with `/api/evidence` and `/api/evidence/discovery` unchanged. **Establishes software capability only**: the evidence and its limits are inspectable; no biological interpretation, no validation. **Open:** the Evidence-tab section (the owner's decision after the permission refusal), hence no ordinary-page visual check; 13 of 23 significant local benchmark rows not assembled; the source-specific readers (about 990 lines) to be separated from the reusable part before any second example | increment 1 | done (UI blocked) |
   | ↳ MAP | **Wording corrected 2026-09-30 (the external reviewer)** | "Establishes software capability only … no biological interpretation" is wrong: the 26 inferred assertions are interpretations, inherited from existing sources. Correct: **the map organises observations, predictions, inferences and missing evidence, and introduces no new biological validation or evaluator verdicts.** "Observed" describes only what an experiment measured, never an indirect effect upgraded into a mechanism; unknown entries are questions or missing relationships, not facts. Accepted provisionally as a backend; **delivered as API and CLI capability**; the Evidence-tab section and its ordinary-page visual check wait for the existing permission gate; no further expansion | increment 1 | note |
   | ↳ MAP | **Page done 2026-10-01 (`d2043c3`), accepted at the report level (the external reviewer)** | After the owner approved the narrow Evidence-tab edit, made through the normal permission process: a read-only section reading `/api/evidence/response-map`, listing assertions by status (observed: what an experiment measured, never a mechanism; predicted; inferred: interpretations inherited from existing sources; unknown: questions or missing relationships), chains as ordered references, the stop points, the dynamics questions and the sources; 66 lines added, none removed; 4 page tests; 72 web and evidence tests pass. Checked on the ordinary rendered Evidence tab (a full-page capture, nothing hidden, cropped to the card); screenshots retained locally, git-ignored, in `data/cache/review_screenshots/2026-10-01_response_map/` (`rmap_fullpage.png` sha256 `e9930441…`, `rmap_page_part1.png` `89f118a8…`, `rmap_page_part2.png` `0e5aaba9…`), beside the discovery view's in `2026-09-29_discovery_view/` (`583ccc56…`, `e30dd58e…`). Increment 1 complete; no expansion | increment 1 | done |
   | ↳ N1 | **Next independent test: revised feasibility proposal, 2026-10-01 (not started; nothing downloaded, fitted or scored)** | **Data:** Replogle et al. 2022, genome-scale Perturb-seq in K562 (CRISPRi knockdown of each gene; transcriptome-wide response). **No prior project exposure found**: no file, script, result or document here reads it; whether the upstream model saw related data cannot be established, so independence is not complete. **Claim, narrowed:** the frozen motif-to-gene predictions anticipate downstream transcriptional responses to a factor's knockdown better than a proximity baseline. It does not validate individual enhancer–gene connections: a factor's knockdown can act indirectly. **Step 1, metadata only:** file listings, sizes, licences, perturbation identities and the measurement fields offered; the potential factor count from the perturbation list intersected with the frozen JASPAR mapping and the existing prediction coverage. No differential-expression result is opened and no factor is selected by its responders. Prefer the authors' processed pseudobulk and summary statistics; the single-cell collection is not downloaded by default. **Step 2, registration before any measurement:** freeze the motif mapping, the prediction score, the TSS window of the proximity baseline, the gene universe, the control-based expression matching, the definition of "responds", and any knockdown-quality threshold, which is fixed before the knockdown measurements are read and whose later assessment is recorded as outcome exposure even when it only decides eligibility. Identical eligible genes for both methods; the perturbed gene excluded from its own target scoring; factors with no measurable responders or an undefined AUROC reported in coverage, never dropped. **Analysis:** the paired gain (attribution minus proximity) per factor, its uncertainty and the evaluable coverage; a factor-level bootstrap is limited because factors share motifs and downstream programmes, and that is stated. 30 eligible factors is a feasibility floor, not shown power; +0.02 AUROC stays the proposed useful gain, without a claim that the study can detect it. **Cost:** metadata requests and one processed download of a size to be confirmed; minutes of compute on cached programs; 0 model or paid requests | — | proposal |
   | ↳ N1 | **Step 1 done 2026-10-01: metadata only (coordinator; no expression value, response statistic or knockdown measurement read)** | **Source:** Figshare+ article 20029387, DOI 10.25452/figshare.plus.20029387.v1, published 2022-06-09, **CC BY 4.0**; 12 AnnData files (K562 genome-wide day 8, K562 essential day 6, RPE1 day 7): single-cell files 8.70–65.83 GB, pseudobulk files 79,766,954 (K562 essential), 95,350,546 (RPE1) and 374,587,922 bytes (K562 genome-wide, normalized md5 a3dfaa94…). The companion article 20022944 holds only sequencing-run manifests. **Read by HTTP range requests, about 6 MB in all:** the `.obs` index (perturbation identities), the names of the `.obs` and `.var` fields, and the measured-gene identities. Never `.X` or any field's values. The `.obs` fields include knockdown and response statistics (`control_expr`, `fold_expr`, `pct_expr`, `energy_test_p_value`, leverage scores): their names only were read, and any later reading of them is outcome exposure. **Perturbations:** genome-wide 9,867 genes (11,258 pseudobulk rows, 585 non-targeting) × 8,248 measured genes; essential 2,058 genes (2,285 rows) × 8,563. The genome-wide library includes factors not expressed in K562 (PAX6, NEUROD1, MYOD1), so a perturbed factor is not yet an eligible one. **Candidate coverage against the frozen JASPAR 2026 mapping (1,019 profiles, 887 factor names) and the committed predictions (`motifs_chr*.json`):** genome-wide, 870 JASPAR factors perturbed; 141 have promoter predictions and 137 element predictions; with the perturbed gene excluded and targets restricted to measured genes, 66 factors have at least 10 element-predicted targets (48 at 20) and 103 at least 10 promoter-predicted genes. Essential screen: 45 factors, 5 with at least 10 element-predicted measured targets, below the floor. **Limits found:** the element arm rests on the committed sample scan of 300 elements per chromosome (7,063 elements, 2,582 target genes, 777 of them measured), not the 440,377 compiled elements, so any full scan would be a new frozen prediction set; the 66 is an upper bound before the knockdown-quality threshold, which reads measurements and belongs to step 2. **Smallest suitable subset:** `K562_gwps_normalized_bulk_01.h5ad` (374.6 MB), not the single-cell collection. Stopped before step 2 | — | step 1 done |
   | ↳ N1 | **Step 2 stopped before the freeze, 2026-10-01 (lane-n1; nothing committed, no value read, nothing downloaded)** | Point 3 first, as required: **the pseudobulk file alone cannot support a per-gene "responds" endpoint.** `.X` (dense, 11,258 × 8,248) is the per-perturbation mean of per-cell gemgroup z-scores (STAR Methods, "Filtering and internal normalization", doi 10.1016/j.cell.2022.05.013; the Figshare description; the 2019 producer code `CellPopulation.average()`, the 2022 code not public). `control_expr`, `fold_expr`, `pct_expr` are undocumented by name; the documented knockdown (mean unnormalized target expression in perturbed cells over non-targeting cells) fits `fold_expr` with `control_expr` as denominator, an inference. `energy_test_p_value` (a permutation test on 20 principal components), the Anderson–Darling and Mann–Whitney counts and the leverage scores are perturbation-level and cannot label a gene; the per-gene test results are not distributed; the file has no layers, `uns`, `obsm` or `varm`. **Smallest additional requirement:** a derivation, no new file: SE = 1/√`num_cells_filtered`, T = X·√n, assuming X averages exactly those cells (unconfirmed; the field is float64), unit variance under no response, independent cells and a normal approximation that fails for sparse genes (calls too liberal on low-expression genes, which may correlate with motif scores). Checking it needs measurement access but no outcome row: T across the 585 non-targeting rows (514 core), selectable from the identity index, about 19 MB by range reads, with its pass rule frozen before any factor row is read. Fallbacks: the authors' per-gene results (on request only) or the 66 GB single-cell file | — | waiting (owner) |
   | ↳ N1 | **Frozen package committed 2026-10-01 (lane-n1; `b359867` code, `f9a9130` registration, `46c10d7` doc), reviewed once by the coordinator; no measurement read** | `genomeos/attribution/n1_perturb_response.py`, `scripts/n1_register.py`, `scripts/n1_run.py` (refuses without a named `--authorisation` and refuses a second run), `data/results/n1_registration.json` (clean stamp at `b359867`, manifest complete, 28 frozen inputs hashed including the 24 motif files, JASPAR, TRANSFAC and the GENCODE set; Figshare md5 of both pseudobulk files). 26 synthetic tests; 43 with the writer guard pass. **Endpoint:** T = X·√`num_cells_filtered`, assumptions stated; responds at |T| ≥ 3 (nominal two-sided rate 0.0027). **Gate first,** on all 585 non-targeting rows: pass only if the |T| ≥ 3 share is at most 1.5× nominal overall and in each of 10 control-expression deciles (each expecting at least 20 exceedances); a failure stops the run, with no fallback inside it. **Eligibility:** at least 25 cells, at least 10 expected target UMIs, `fold_expr` ≤ 0.40 (inferred meaning, quoted), `pct_expr` = `fold_expr` − 1 or the run stops; multiple rows pooled by cells; heterodimers excluded; families as bootstrap clusters. **Analysis:** per-factor AUROC over the 860-gene shared universe (the perturbed gene and genes within 10 kb of its TSS excluded), undefined AUROCs counted by reason, equal-weight paired mean, 10,000 family-clustered bootstrap draws; +0.02 and a lower bound above 0 reported apart; 30 factors a floor. **Coverage:** 860 of 8,248 measured genes (322 without a scanned promoter, 7,066 with a promoter but no sampled element); 56 factors (59 rows, 39 clusters) before eligibility. **Observation from the review:** the arms flag very different shares of the universe (for CGGBP1, 40 genes by element attribution against 598 by promoter proximity), a frozen property of the committed predictions that bears on how a gain is read. **Measurement access the run needs:** both pseudobulk files (374.6 MB each, md5-checked), then the gate's rows, then the knockdown fields at the 59 candidate rows, then normalized X at eligible rows only | — | frozen (owner decides on the run) |
   | ↳ N1 | **Amendment 1, 2026-10-01 (lane-n1 on the external reviewer's code review; `5d1710a` code, `8ebaa9a` registration, `7a875cb` doc), checked by the coordinator; no measurement read** | Versioned and additive (497 lines added, none removed); `f9a9130` and `b359867` preserved; the old runner refuses and points to `scripts/n1_run_v2.py`. (1) **AUROC:** a constant score with both classes present is 0.5, not undefined (`auroc_v2`, ties count half; undefined only when a class is empty). (2) **Knockdown units unsupported:** the STAR Methods define knockdown as the ratio of mean unnormalized target expression in perturbed against non-targeting cells, but the fields are undocumented, the methods also produce depth-adjusted counts, and the 2022 producer code is not public; cells × `control_expr` is not established as expected UMIs and the `pct_expr` check shows algebra only. **v2 stops before the eligibility rule** (`ELIGIBILITY_SUPPORTED = False`), with no replacement threshold. (3) **Framing:** N1 is an exploratory comparison against an operational response label; the gate does not validate per-gene responses, cell independence, reference-mean uncertainty or calibration across cell counts and pooled rows; the 585 controls include the 514 core ones used in the normalization. (4) **Expression matching:** the departure is recorded (`f9a9130` banded expression only in the gate); the primary is now an AUROC stratified by control-expression decile, pair-weighted, the unstratified figure reported without a criterion; frozen and tested, not run while (2) stands. (5) **Freeze check:** v2 refuses before opening any file unless the original registration (`8478dd76…`), the module (`0ea47ad9…`), the runner (`b33ff4f3…`) and the constants match, which they do at `5d1710a`; its reader decodes only the 860 selected columns, and the original reader's whole-row decoding is recorded. 50 tests pass, 1 skipped with its reason. **An authorised v2 run would download both files, run the gate on the 585 controls, and stop at `eligibility_unsupported`, reading no candidate row** | — | frozen; blocked on (2) |
   | ↳ N1 | **Blocked for full execution, 2026-10-01 (the external reviewer accepted amendment 1)** | No further code until the measurement definitions are resolved; eligibility is not relaxed, no measurement is opened and no new threshold is set while waiting. **Inquiry drafted for the owner to authorise:** to the dataset's listed contact (Figshare, `joseph@wi.mit.edu`) and the corresponding author who co-led the Perturb-seq data analysis (`normantm@mskcc.org`), asking how `fold_expr`, `pct_expr` and `control_expr` in `K562_gwps_normalized_bulk_01.h5ad` are calculated (normalization, control selection), whether `num_cells_filtered` counts the cells averaged into each normalized pseudobulk row, and whether a documented knockdown-quality flag or recommended eligibility rule exists, with a request for producer code or a data dictionary. **Separately proposable:** the gate alone as a control-tail diagnostic; passing it would show only that the selected controls satisfy the registered exceedance rule, not that the statistic is generally calibrated or that the experiment is ready, and it cannot resolve the knockdown definitions | — | blocked (owner) |
   | ↳ N1 | **Amendment 2, 2026-10-01 (lane-n1; `29e7a88` code, `fd3d928` registration, `8764bc6` doc), checked by the coordinator; no p-value, expression or knockdown value read** | On the owner's direction (no author contact; the authors' own per-gene results): **Figshare+ 21632564** (Replogle and Weissman, CC0, doi 10.25452/figshare.plus.21632564.v1), `anderson-darling p-values, BH-corrected.csv.gz`, 488,720,141 bytes, md5 `abb0310e…`; per STAR Methods, a per-gene Anderson–Darling test of normalized expression in perturbed against non-targeting cells, **adjusted by Benjamini–Hochberg** (the adjustment family undocumented), the documented level p < 0.05. **The question changes:** assigned perturbations, no knockdown filter; ineffective perturbations may weaken the signal, an absent response does not show a factor has no regulatory role, and no factor is selected by its responses. **Responds:** published adjusted p < 0.05; empty, non-numeric, non-finite or out-of-range entries, and genes absent or duplicated, are excluded and counted, never non-responders. **Coverage from header and identities only:** the file's 11,258 columns are the pseudobulk's row labels; all 56 candidate factors present (CGGBP1, FOXD3, LHX3 on their P1 rows, the paper reporting P2 perturbations generally without effect); gene coverage of the 860-gene universe is counted in a labels-only pass at run time, before any value is parsed. **Strata:** deciles of the raw pseudobulk mean over the 585 non-targeting rows, rule frozen, values read only in the authorised run. **Kept:** the 860-gene universe, the 10 kb exclusion, amendment 1's stratified AUROC primary, the 30-factor floor, the bootstrap and the two separate criteria. **Dropped:** the derived |T| label, its gate and the knockdown rule. **Freeze:** `scripts/n1_run_v3.py` refuses unless both earlier registrations, the module (`9be800c9…`), itself (`406b60da…`), the constants and the published header digest match, which they do; the earlier runners refuse. 58 tests pass, 1 skipped. Additions only (461 lines); `f9a9130` and `8ebaa9a` unchanged. **A run needs:** the published file (488.7 MB) and the raw pseudobulk (374.6 MB) with md5 checked and sha256 recorded; raw X at the 585 control rows for the strata; a labels-only pass over the published file; then the 56 principal columns at covered universe genes | — | frozen; owner decides on the run |
   | ↳ N1 | **Run once under amendment 2, 2026-10-01: insufficient coverage, no estimate (coordinator, on Albert's authorisation of 15:19 BST; `data/results/n1_result_amendment_2.json` from `1cb790a`, clean)** | Freeze check passed (registrations `8478dd76…`, `ecfb55fc…`; module `9be800c9…`; runner `406b60da…`; constants; header digest); both downloads match Figshare's md5 (raw `4570b53c…`, published `abb0310e…`); identities match. Registered reading, word for word: "12 defined paired differences, below the floor of 30; no estimate is reported and no rule is relaxed." **No paired gain, no interval; neither criterion assessed** (+0.02 point gain; lower bound above 0). **Coverage:** genes 496 of 860 in the published file (364 absent, none duplicated, no missing entry at the 496); factors 56 available, 56 analysed, **12 evaluable**, 44 undefined in both arms for `no_responders` (no gene below adjusted p 0.05 among the 496); families 39 analysed, **12 represented** among the evaluable, one each; 56 responders in all (1 each for 5 factors, 2–3 for 5, MEF2A 11, NRF1 27), compared within 1–9 control-expression deciles of 49–50 genes. Per-factor stratified AUROCs for both arms and the unstratified mean over the 12 (−0.061, no criterion, below the floor, not an estimate of the gain) are in docs/ATTRIBUTION.md. **Supports:** the frozen pipeline (candidates, universe and endpoint together) cannot answer N1's question with this endpoint; a limitation of the pipeline, not a finding about the predictions. **Cannot establish:** better or worse prediction; that the 44 factors lack a regulatory role in K562 (assigned perturbations, no knockdown check, expression in K562 not examined); adjusted p < 0.05 marks a detected distribution difference with no magnitude or direction; the BH family is undocumented, the published level used as frozen, not re-adjusted, no 5% FDR claimed for the subset; the 30 floor is not power. **Reporting gap found in the pre-run check:** genes inside the used deciles are not a result field (derivable from decile sizes). **Cost:** no money, no model requests; 863,308,063 bytes of CC0 downloads in 25 s (kept in git-ignored `data/cache/n1/`); run 14.3 s wall, 8.9 s CPU, 85 MB peak. Note to the amendment 2 row: its "58 tests" counts the N1 file with the results-writer guard; the N1 file alone is 41 passed, 1 skipped | — | done: stopped below the floor; no further N1 amendment proposed (a new test needs data not yet read; owner decides) |
   | ↳ N1 | **Closed 2026-10-01 under its registered stop rule: completed, inconclusive (coordinator, on the external reviewer's advice relayed by Albert)** | Results preserved; no threshold relaxed, no factor selected, universe not enlarged. **Wording correction beside the run row and the record:** "no primary aggregate estimate or decision", not "no estimate" or "no evidence either way"; the per-factor scores and the secondary unstratified mean (−0.061) are descriptive results below the reporting floor, not a reliable conclusion about comparative performance. **Coverage diagnosis (identities and documentation only, 15:59–16:03 BST, no response or expression value read):** identifier mismatch ruled out (the published file's 5,530 gene rows are unversioned Ensembl IDs, unrepeated, all among the pseudobulk's 8,248); the pseudobulk's gene set is documented (">0.01 UMI per cell", Figshare+ 20029387), the published test's is not; the reason for the 364 absent universe genes is unresolved, narrowed to an undocumented restriction by the producers of the genes they tested. **The larger limit is the candidates:** 39 of 56 factors' own genes are absent from the pseudobulk's gene list (by ID and by name), 5 present but untested, 12 tested; the 44 factors without a responder split 33/4/7, the 12 evaluable 6/1/5; NRF1 is absent from the list yet has 27 responders, so the literal reading is recorded, not relied on. **For another design:** draw the universe from the 5,530 tested genes and check each factor's own gene against the expressed-gene list, both from identities before any value is read; a new design on this deposit is no longer blind for these 56 factors | — | closed: completed, inconclusive |
   | ↳ After N1 | **Next, 2026-10-01: the wiring diagnostic, with paired enhancers held at metadata-only feasibility** | No earlier wiring diagnostic is recorded in ROADMAP or ATTRIBUTION.md, so it starts now. **Question:** does the assigned element-to-gene wiring carry predictive information beyond what it shares with simpler alternatives? **Candidate benchmark:** the already-examined K562 CRISPRi pairs (10,356; `scripts/crispri_published.py`, `data/results/crispri_published.json`), with scores and evaluation rules fixed: the real per-pair deletion feature against rewired assignments matched on distance, coverage and connectivity; registered before any null is computed; no model requests. **Stop:** if distance, coverage and connectivity cannot be controlled adequately, report that and stop. Development evidence, not independent validation. Paired enhancers stay at metadata-only feasibility until a concrete, adequately covered test exists; no new framework; the coherence pilot is not reopened | lane-wiring | feasibility first |
   | ↳ Wiring | **Closed 2026-10-01 as infeasible on the K562 CRISPRi benchmark (lane-wiring; `31e624e` thresholds and code before any diagnostic, `f58a2a3` descriptive counts after the decision, `9b88526` result and record); no registration written, no rewired assignment scored** | The primary rewiring keeps each element and attaches its cached deletion value to another gene the element was measured against at a matched distance; pairs, labels, features, split and estimator unchanged; a link that cannot move keeps its real value. Four schemes (distance bins 0.2 or 0.3 log10, with or without K562 expression quartiles) were fixed in order, with thresholds set before diagnostics; **none met them**. **Three quantities kept apart (training, S1 / S3):** *links that can move:* 60 / 67 of 423 regulated, 37.0 / 47.0 moved per rewiring against a floor of 100, and only 38 / 42 can change class (the rest swap only with another regulated link of the same element); *distinct alternative assignments:* 88 / 107 over 62 / 74 genes, 52 / 58 links with one or two alternatives, the same alternative in 91% of the draws in which a link moves (S1); *independent genes or loci:* 42 / 49 genes, 28 / 32 loci (1 Mb chaining, rough) against 290 genes and 216 loci for all regulated links. Expression matching leaves 6 to 9 moving. Held-out K562 and the secondary rewirings also fail. **The 100-link floor was a feasibility requirement, not demonstrated power.** Why: 188 of the 423 regulated links are their element's only link with a cached value, and 360 are its nearest. The expression covariate reuses N1's cached control expression (disclosed in the record). Nothing loosened, no predictions bought, no extra rewirings. **Says nothing about the wiring either way**; a property of this benchmark and the cache's coverage, not of the model. Reopening needs a benchmark with several measured genes per element at comparable distances, under a new registration. 0 model requests | — | closed: infeasible |
   | ↳ Wiring | **Clarification closed 2026-10-01 with three wording corrections (coordinator, on the owner's questions; no scoring, no randomisation search, no new registration)** | `scripts/wiring_counts.py` lists every count read-only, fitting no model and computing no metric: `uv run --frozen python scripts/wiring_counts.py`. **Scheme:** primary rewiring `element_kept`, strata (element, distance bin), bins 0.20 log10 at S1 and 0.30 at S3 with a random offset per rewiring, no expression classes, a derangement on the gene so in-degree and link counts are exact, 20 seeded draws. **S1 alone:** 60 of 423 regulated links have an alternative; moved per draw mean 37.0, smallest observed 30, largest observed 42; 59 unique links moved across the draws; 2.0 per draw receive an identical feature block; distinct genes 42, elements 46, chromosomes 13; locus clusters 28. **S3 alone:** 67 with an alternative; moved per draw mean 47.0, smallest observed 39, largest observed 51; 63 unique links moved; 2.8 identical per draw; distinct genes 49, elements 50, chromosomes 14; locus clusters 32. **Which scheme gave the earlier figures:** 60 is S1's links-with-an-alternative and 47 is S3's mean moved per draw, so the two were quoted together across schemes and measure different things. **Correction 1:** "independent" is withdrawn from the rows above; these are distinct genes, elements and chromosomes, and 1 Mb chaining gives operational locus clusters, establishing no biological independence. **Correction 2:** "at most 47 can change per rewiring" is withdrawn; 47.0 is S3's mean and the largest observed at S3 is 51, with no maximum over all rewirings proved. The closure rests on the tested schemes failing the frozen feasibility rules. **Correction 3:** "so the score would change" is withdrawn; what was measured is feature identity (the attached `top_target` and `deletion_drop` match the link's own in 1 of 23 assignments at S1, 0 of 28 at S3), and a differing block does not establish a different score because no model was fitted. The owner's reading holds: those swaps never test the assignment against a non-target and leave the class unchanged, so admissibility did not make an alternative informative | — | closed: infeasible; clarification closed |
   | ↳ Paired enhancers | **Metadata feasibility done 2026-10-01; one local proposal set out and recommended against, for the owner's review (coordinator; nothing implemented, no outcome scored, no money, no model requests)** | Generalisation stays closed: no public dataset has designed enhancer pairs across many loci. Surveyed: Lin 2022 (MYC), Hsiung 2024 (same locus), Xie 2017 (Mosaic-seq), Gasperini 2019 (pairs incidental to high multiplicity, not designed), Pacalin 2024 (CRISPRa+CRISPRi, T cells), mouse deletion series (no overlap). **Only candidate: Lin et al. 2022, GSE160768** — K562 dCas9-KRAB, 87,025 sgRNA pairs over all single and pairwise combinations of the 7 MYC enhancers e1–e7 across 1.8 Mb, 30 doublings (quoted from the series design); 21 enhancer pairs, 1 locus; controls 40 non-targeting, 34 at 3 negative regions, 31 at the MYC promoter; 2 biological replicates; counts in Zenodo 10.5281/zenodo.6823833 (MIT, ~330 kB, **not opened**, they hold outcomes); **no prior exposure anywhere in ATTRIBUTION or ROADMAP**, so it could be registered unexposed. **Two verified blockers:** (1) **no expression readout exists in the series** — its 58 samples are 16 ATAC-seq, 22 ChIP-seq, 20 other, no RNA-seq — while the project's quantity is a log2 expression change, so a fitness-to-expression mapping would be an assumption under test; (2) **of the project's elements overlapping e1–e7 (26 by this check against `all_elements/chr8.json`, the lane's 23; overlap convention, unresolved), none predicts MYC** — the targets are MIR1206, MIR1205, MIR1207/1208, MIR1207, CCDC26 and ENSG00000285108, MYC appearing only at the 3 promoter-overlapping elements beside CASC11. **No paired capability either:** the project predicts single-element expression change, and R8's coupling pretest already closed with no coupling term beating independent deletion scoring, so a paired prediction needs a new model. **Decision value low:** a fitness change cannot separate a correct assignment to a microRNA from an unassigned MYC effect without an expression measurement. A singles-additive baseline is constructible in the same library (a single is a targeting plus a non-targeting sgRNA), the one requirement well met. Chromatin arm covers WT, 4 singles and **2 pairs** only. **Recommended against on local terms, not for failing to generalise:** no quantity on both sides to compare. **Would change it:** an expression readout on single and paired perturbation of these enhancers in the same context, or a dataset whose elements carry a project prediction for the gene it measures; Hsiung 2024 (GSE260832) is the natural second check but is also a fitness screen | — | proposal recorded; owner decides |
   | ↳ Paired enhancers | **No-go accepted 2026-10-01 by the owner for the proposed comparison with current predictions, with three corrections to the row above (coordinator; search not reopened)** | Accepted on the two grounds that hold: the comparison mixes a growth readout with expression predictions, and the project has no frozen paired-response model. **(1)** "No public dataset has designed enhancer pairs across many loci" is too broad; what is supported is **no suitable dataset was identified in this survey** — Lin 2022, Hsiung 2024, Xie 2017, Gasperini 2019, Pacalin 2024 and three mouse series, by literature and repository metadata in about 90 minutes. Not a proof of absence. **(2)** Listing "none predicts MYC" as a blocker is **withdrawn**: with a suitable expression measurement that disagreement would be a prediction failure the benchmark exposes, and **selecting only experiments whose implicated gene the model already predicts would bias validation**. The verified count stands; it does not disqualify the dataset. The disqualifying reasons are the readout mismatch and the absent paired model. **(3)** "No RNA measurement at all" is narrowed: of the deposit's 58 samples none is an RNA-seq library (16 ATAC-seq, 22 ChIP-seq, 20 other), but a published expression measurement such as supplementary qPCR could exist outside them and the supplementary material was not retrievable. **Also:** Pacalin 2024 combines activation at one element with repression at another, a different intervention from paired silencing with different comparability requirements, not assessed against them. **Carried forward:** a research proposal states what the project predicts, what the experiment measures and how the two compare, then verifies coverage and controls, before any analysis code; eligibility never depends on whether the model already predicts the implicated gene | — | closed: no-go accepted |
   | ↳ E2/E3 | **Widening further is a no-go at 0 requests, 2026-10-02 (lane-e23; nothing built, nothing committed by the lane, no request sent)** | **The widening already happened.** E2R and E3R ran on 2026-09-15 over the 16 replication chromosomes (chr6–chr20 and chrX) for 3,000 requests, and their registered readings stand word for word: **E2R +0.091** on 1,166 pairs, one-sided p **0.0043**, upper 95% bound 0.146, verdict `negative` under the rule written before the run — *"a difference between 0 and 0.10 resolves nothing"* — the discovery run's **+0.268 falling far outside** that interval; **E3R −0.030** on 775 pairs, upper bound 0.042, verdict `futility`. **What remains on disk under the same bars, none moved:** 1,505 unspent E2R pairs of 2,671 assembled, of which **654 lie at 184 independent loci where nothing has been read — 1,308 requests, 7.1 per locus**; and 269,783 unspent E3R pairs of 270,558, of which 135,738 lie at **144** never-read loci for 271,476 requests. **E3R's 141,673 unspent units chain into 206 loci in all** under the imported 1 Mb convention, so its entire remainder buys at most 144 places in the genome; that is what the rule says about this population, not an argument to relax it. **Nothing is free:** the prediction cache covers **1 of 3,010 E2R sides and 0 of 539,566 E3R sides**, and no unspent pair's outcome has been seen. **The only genuinely new population cannot be counted at all:** chr1–chr5 hold 682,295 catalogued storage units, every catalogue written *after* the 2026-09-15 assembly, but the genome-wide GTEx distillation `data/knowledge/human_panel/executor/gtex_wide/` is **absent from disk** — only the chr21/chr22 distillation survives — so a 1.4 GB re-stream of `GTEx_Analysis_v8_eQTL.tar` precedes any chr1–chr5 pair count, and the lane refused to quote its order-of-magnitude extrapolation as one. **No-go on three named requirements:** that absent distillation; **a new pre-registration**, because spending the remainder *after* seeing +0.091 fall short of the pre-registered 0.10 is optional stopping run backwards — the same objection §5 item 1 raised against E1's extension and resolved only by registering first; and the limit the registration states itself, *"no outcome here makes E2 and E3 independent of the model"*, AlphaGenome's training being kin to GTEx, so the whole arm is **instrument evidence** and cannot answer the ENCODE objection, which stands **narrowed, not closed** after E1's extension returned +0.075 (p 0.00037) in its own WEAK band. **Prior exposure:** every one of the 16 chromosomes has been read, so nothing here could be fresh validation; 275 of E2R's 459 loci and 62 of E3R's 206 already contain a read pair. If the owner funds it anyway the smallest honest version is **1,308 requests** for 654 pairs at 184 never-read loci, gated behind the new pre-registration and labelled instrument evidence everywhere it appears. **Recommendation: spend nothing.** 0 requests sent, no money, no byte downloaded, no cached value read | — | no-go |
   | ↳ Defect | **Fixed 2026-10-01: the superseded node-containment control no longer appears in the genome-wide script's report (coordinator; reporting only, no biological rerun, no rescoring, 0 requests)** | The defect this plan records at the node-containment claim: `scripts/enhancer_targets_all_genome_wide.py` printed "(random boundaries 0.791, +2.6 points)" on every run, two figures from the 113,399-element deletion archive shown beside a 440,377-element measurement as if they were its control, superseded by the 2026-09-27 audit. **Now** `node_audit_lines()` reads `data/results/node_containment_audit.json` at run time and reports its figures as that dated audit's, printing its date, its path and "not recomputed here", so no later fold can drift from it or present them as newly computed; an unreadable audit quotes no figure. It states the headline **+2.90 percentage points** on the share of elements whose coding target lies inside the element's own CTCF node **against its named comparator, uniform random boundary placement** (the audit's published control: as many uniform positions as the caller has edges, no 50 kb merge) over 440,377 modelled pairs; the **95% CI +2.03 to +3.81** as the interval of that excess against that control, from 2,000 bootstrap resamples **over the 24 chromosomes, not the 440,377 pairs**; ahead on **19 of 24** chromosomes; **separately** the four defensible baselines spanning **+1.21 to +6.58**, the headline one of them and not their summary; and the qualifier that this is **internal benchmark evidence on the model's own reading, not observed enhancer-gene pairs and not independent validation**. The superseded entry stays in `CONTROLS` because it records what was superseded, and a test fails if it ever sits there without its `superseded_by` note. **Tests:** 23 (15 before, 8 new) covering the removal and each of those meanings; the removal test was confirmed to fail with the old line restored and to pass once reverted | — | done |
   | ↳ Defect | **Fixed 2026-10-01: per-cell deletion values are no longer collapsed to one track, with the legacy data left legacy (coordinator; `genomeos/predict/enhancer_target.py`, `tests/test_enhancer_target.py`; 0 model requests, no committed result touched)** | The limitation S4 recorded on 2026-09-29 as found and not fixed: `aggregate()` wrote `by_cell[name] = value` per emitted track, so where several tracks carried one cell's name the last was kept and the rest discarded at write time — and the versioned judge rules on direction only in the stated cell. **Now** each cell name also carries `values` (the multiset received, sorted so the summary is independent of arrival order), `emitted`, `signs_disagree` (computed apart from any average, since opposite values average to zero), `mean_of_emitted_log2fc` **labelled as a descriptive statistic of emitted values only — the tracks are not known to be biological replicates and nothing reads it as a cell-level effect**, and a `completeness` note. **`emitted` is a lower bound, not a track count:** the adapter passes a value only when |log2 fold change| > its threshold, so values at the threshold and exact zeros never arrive; the cell's total track count is not recorded and not inferred, and no track identity is invented. **Prediction behaviour deliberately unchanged:** `by_cell` still holds the last emitted value, `predict_target` still selects by largest drop or rise and reads neither field, so this does **not** fix claim-context selection; about 20 modules and 13 scripts were traced and **none switched** to the new field, a mean or a different direction rule — that needs its own registration. **Versioned:** new rows carry `by_cell_schema` 2 and a summary; `per_cell()` returns an older row as `legacy_last_track` with no count and no aggregate status. **Unrecoverable:** the raw per-track values are discarded in `score_element()` and nothing retains them; all 24 chromosome archives were checked and carry no summary, so every stored per-cell value is legacy and only new requests could recover the rest — none were made. **How often a cell name carried several tracks cannot be answered from disk** (the track table was never stored); known: 371 emitted values per element at the mode against 316 distinct tissue names, so names do repeat, and the 2026-09-29 trace showed individual cases. `emitted` records it from now on. **Tests:** 21 (10 new, behavioural); reducing the summary to one value per cell fails four | — | done |
   | ↳ Defect | **Correction 2026-10-01: the judge never read the collapsed per-cell value, so the fix above implies no re-judging (coordinator, on the owner's reading of the code)** | The row above justified the fix partly by a dependency that does not exist, and that justification is **withdrawn**. Checked: `genomeos/attribution/correctness.py` references `by_cell` **nowhere**; `compiled_claims()` takes the gene, action and cell from the compiled rule's own fields; `_direction_v2()` (used under v3) compares that claim against **experimental observations** in the stated cell, never against the model's per-cell numbers; and `predict_target()` names the winning track's tissue, which **becomes** the claim's stated context, so the two are the same cell in this path and cannot disagree. The 2026-09-29 trace's "cell mismatch" was between a claim's cell and the cell the screen measured; treating it as a mismatch between the collapsed value and the stated cell conflated the two. **No judge registration, no re-judging; current judges and every committed result unchanged.** **What the fix is worth, precisely:** future stored answers preserve more evidence and an explicit reader identifies legacy representations; **existing consumers gained no protection automatically**, since `by_cell` keeps its meaning and none was switched. **Conditions recorded for any future abstention policy on track disagreement, before such work starts:** it belongs before a claim is issued, not in the evaluator, because an evaluator that skipped already-issued claims on model uncertainty would drop the difficult ones and make reported performance misleading; it must report coverage beside performance; legacy disagreement is **unknown**, not absent; and a numerical sign difference must be distinguished from a materially conflicting effect, since with track identities unavailable disagreement alone does not establish biological inconsistency | — | closed |
   | ↳ Defect | **Done 2026-10-01: a reusable discovery-overlap primitive, with the evidence interface left alone (coordinator; `genomeos/attribution/measured.py`, `genomeos/attribution/ablation.py`, `tests/test_measured_rules.py`; nothing adopted, re-scored or re-judged, 0 model requests)** | The limitation: the production rule needs the overlap to be at least half of **both** widths, so it pairs intervals only when their widths are within a factor of two, and audit A's alternative (dd8c49c, `min_side_half`: half of the **smaller** width) lived only in the census script. **Delivered:** `ATTACHMENT_RULES` names both predicates, `attaches(rule, ...)` applies one by name beside the new `min_side_overlap()`. **An overlap under the broader rule is a candidate, not a finding:** it does not establish which element produced a measured response. **The evidence interface is untouched:** `CrispriIndex.of()` takes **no** rule argument and applies the production rule only, because returning one observation independently for several elements would make it decisive evidence for each; a discovery pass needs observation ids, whole candidate sets and rule provenance, and **that interface is not built here**. **Bug found and fixed before commit:** the first draft validated the rule name only on the path that computed an overlap, so an unknown name could be answered with a plain "no overlap" when the chromosome was absent or the slice empty — the "always refused" claim was false as written; the name is now looked up before any coordinate is touched, with the empty cases tested. **The implication test establishes overlap inclusion only** (20,000 random pairs): **not** preservation of unique attribution and **not** verdict eligibility, since a broader rule can retain every observation while resolving each to several elements, lowering resolution rather than adding evidence. **The counts are the historical census's** (`placement_census.json`, 2026-09-29: 180 newly overlapping, 128 → 308 of 14,734) and **their reproduction was not performed**; the acceptance scope is narrowed to the predicates' definitions, the inclusion property, the width asymmetry, the threshold, the interface's silence on rules, and name refusal whatever the intervals. **The census's resolution figures do not identify elements:** of the 180, 48 unique, 77 "one link" — **which is not a uniquely identified element**, since several candidates can share one target link — and 55 several, so most stay unresolved as to which element produced the response. Adoption for attribution remains a separate decision with its own registration. **Tests:** 10; 63 pass across this file, the CRISPRi suite and the correctness suite | — | done |
   | ↳ Defect | **Fixed 2026-10-01: a deletion gain is refused where the feature was never available, and every interval states how it was made (coordinator, on a peer session's audit, checked against the files first; `genomeos/attribution/crispri.py`, `tests/test_crispri_gain_reporting.py`; the committed result is not rewritten, 0 model requests)** | **Defect 1:** `score_published` computed `deletion_gain` for every held-out cell type and recorded `deletion_available` separately, so `crispri_published.json` (2026-09-27) carries a gain and an interval for **three strata with no deletion value**: HCT116 +0.0004 [-0.0002, +0.004], Jurkat +0.0038 [0.0, +0.0059], WTC11 -0.0049 [-0.0174, 0.0]. Outside `MODEL_CELLS` no deletion value was ever read, so the difference is an artefact of the model form; **all three must be read as "no gain measured", not as gains**, and two have intervals quotable as a small positive. The peer named WTC11; the file has three. **Now** `gain_where_available()` does not compute the gain at all where the feature is missing, returning `gain: None`, `ci95: None` and a reason that names the artefact, so a number cannot leak by being blanked afterwards. **Defect 2:** every interval rested on 200 chromosome resamples with nothing said about it — at 200 draws each tail is a handful of values. Both estimators now report `clusters`, `draws_requested`, `draws_dropped` and `met_minimum` against a stated `MIN_RESAMPLES = 1000`, and `BOOTSTRAPS` rises 200 → 2000 for new runs. **No existing file is recomputed**, so no committed interval changes. **Deliberately not done:** the committed result is **not rewritten**; rebuilding it under a new name would change every stratum's interval and the quotable headline, which is the owner's decision, not a reporting fix. **Tests:** 8, including a spy proving the gain is not computed in an unavailable stratum; restoring the old behaviour fails two; 38 pass across this file and the two CRISPRi suites | — | done |
   | ↳ rE2G | **Paired against ENCODE-rE2G on identical held-out pairs, 2026-10-02 (lane-re2g; `f22136b` `7cc4b3b` `2dd6cb5` `d2ca94d` registration before any comparator score, `bdba855` implementation and the second registration, `95da5e4` `20ddd52` `948b1e6` `36dd6de`)** | **THE LIKE-FOR-LIKE COMPARISON IS THE BINDING ONE AND IT IS NEGATIVE.** Registered in advance with its falsifier: `dnase + distance + deletion` against reconstructed ENCODE-rE2G on the K562 primary (1,918 pairs, 118 positives) is **+0.0593, 95% [−0.001, 0.1485] — "no difference detected"**. So **README may not say the deletion model ranks better than ENCODE-rE2G**, and the H3K27ac comparison may not be cited alone. **The registered H3K27ac comparison**, kept with its qualifiers and never alone: +0.0689 **[+0.0022, 0.1605]**, 22 clusters, 2,000 of 2,000 draws — "ranks better than ENCODE-rE2G on these pairs"; pooled (4,378 pairs) +0.1376 [0.0734, 0.2067]. **The margin is not robust:** both comparisons are borderline, the point estimates differ by 0.0096 and the lower bounds by 0.0032, and with matched inputs the interval includes zero. **No causal claim is made about H3K27ac** — an adjustment that moves a borderline bound across zero does not identify its cause; the selection is stated as the reason the matched comparison is the cited one, which is that **all 190 of 190 held-out positives carry H3K27ac** (52 plus 138, none in No-H3K27ac, CTCF or H3K27me3) against 1,438 non-regulated pairs lacking it, where our activity term is `sqrt(DNase × H3K27ac)` and rE2G's held-out model reads DNase only. **Gate passed on its registered terms:** rE2G reconstructed from the five ENCODE portal files Supplementary Table 12 names (1,569,167,142 bytes, md5-verified) reproduces its published pooled weighted AUPRC **0.5393 against 0.556151, miss 0.0168**, inside the registered 0.02; the join is the benchmark's own, pinned at `v1.0.0` (`50587422`, byte-identical to main), 3,956 of 4,378 pairs joined and 422 scored 0 by the benchmark's own fill rule. **GM12878 (7 clusters) reports a gain and no interval, so no reading**, and its earlier "+0.2713 [0.0418, 0.5293] → ranks better" is **withdrawn** rather than replaced by a hand-written interval. HCT116, Jurkat and WTC11 report **feature unavailable**, never a number, and carry **no** deletion value at all — 0 of 396, 0 of 75, 0 of 1,921. **Identity:** 1,918 K562 pairs and all 4,378 compared, **0 differing** deletion values and 0 differing scores, so "ours" is the frozen headline model; K562 AUPRC 0.7272 and pooled 0.677 equal the headline's. **A free exact cross-check:** against a common comparator the two deltas differ by the deletion gain, giving **+0.136** against the headline's independently derived **0.1361**. **Post hoc, descriptive, one of several unregistered arms:** `dnase + distance` alone sits below reconstructed rE2G on these pairs (−0.1227, interval [−0.2379, −0.0204] excludes zero); `activity + distance` alone is −0.0671 [−0.1768, 0.0291] on K562 and +0.0281 [−0.0506, 0.1004] pooled, both including zero, so no directional word is used for them. The earlier coordinator claim that the pooled delta is largely not a test of the deletion feature is **refuted by measurement**: a count of zero-valued columns does not measure a feature's contribution. **Corrected in passing:** "3,849 pairs with a deletion value" was element coverage mislabelled, overstating reach threefold; the coordinator introduced that figure and relayed it. **Limits:** a reused benchmark, two frozen models, not fresh validation; rE2G was trained on K562 pairs from the same compendium and held these out itself; HCT116's biosample is **chosen, not established** (1 of 16 candidates) and the residual 0.0168 may live there. 0 model requests, no money | lane-re2g | done |
   | ↳ rE2G | **Two corrections to the row above, 2026-10-02 (the reviewer withdrawing its own ruling, and a breach neither it nor the coordinator caught at first)** | **(1) The `dnase + distance`-alone arm is NOT post hoc, and the row above mislabels it.** It is in the registration code at `bdba855` lines 119–123 as `dnase_plus_distance_alone`, with its own reading, **committed at 00:41 before any figure of it existed at 01:47**. So it is a **pre-specified secondary arm within the registered K562 population**, and its reading may stand with that label: −0.1227, interval [−0.2379, −0.0204], excludes zero. The reviewer called it post hoc at 01:36 and has withdrawn that; the coordinator relayed the error into the row above and withdraws it here. The `activity + distance` arm is unchanged: its interval includes zero, so no directional word is used. **(2) A registered reading word was attached to a descriptive population.** `re2g_likeforlike.json` carried `paired_delta.pooled_descriptive.reading = "ranks better than ENCODE-rE2G on these pairs"`, and the top-level `reading` dict listed it beside the primary — while the lane's own registration (`scripts/re2g_likeforlike.py:16`, `bdba855`) says "Registered population: the 1,918 held-out K562 pairs. The pooled population is descriptive context." Lifted on its own that sentence reads as "with matched inputs it ranks better", **the exact overstatement this lane exists to have withdrawn**, three lines below a primary reading of "no difference detected". **Neither the coordinator's two reports nor the reviewer's first review caught it.** Fixed by the lane: descriptive populations carry `reading: null` with "descriptive: no registered reading", the top-level `reading` holds the registered population only, and a test fails if any population not named as registered carries a reading word — written over the result's populations rather than naming the one that broke it. Only output labels changed; **no number moved** | lane-re2g | corrected |
   | ↳ Rebuild | **The benchmark result rebuilt under its own name, 2026-10-02 (lane-rebuild; `9dedf52`, with `156de87` `cd263bc` `c67bbfe` before it); the 2026-09-27 file is kept untouched** | Albert approved the rebuild; the coordinator and the supervisor both recommended it. **What moved:** every interval is now 2,000 draws instead of 200 and states its clusters, draws requested and draws dropped. **K562 +0.1361, 95% [0.0734, 0.2270], 22 clusters, 0 dropped — the headline gain survives the stricter interval.** **What is now withheld rather than reported:** GM12878's gain (+0.0146) and the carried HCT116 arm's (+0.0222) carry **no `ci95` at all**, on 7 and 5 clusters against the ten-cluster minimum committed at `57681b2`; neither was replaced by a hand-written interval. **The HCT116 arm is carried, not dropped:** it was **absent** from the first rebuild, which would have read as never measured, and it is now carried from `crispri_published.json` at its sha with its **705 delivered requests** recorded, so no later reader re-buys them. Strata with no deletion value report null with the reason, never a number. Old and new figures sit side by side with each population named. 0 model requests — every deletion value is the sweep's own cache | lane-rebuild | done |
   | ↳ Decision | **2026-10-02, Albert: the CRISPRi published-benchmark result is rebuilt under a new name (lane-rebuild, started 00:15)** | Approved after the coordinator and the supervisor both recommended it, because the committed `crispri_published.json` of 2026-09-27 carries two defects fixed in code at `bc28f79`: every interval rests on 200 chromosome resamples with nothing recorded about how it was made, and three strata report a deletion gain and an interval for cell types with no deletion value (HCT116 +0.0004, Jurkat +0.0038, WTC11 -0.0049). **Conditions, all of them Albert's through the supervisor:** a new result name with the old file untouched; 2,000 draws; a stratum whose feature is unavailable reports null with the reason; every interval states its clusters, draws requested and draws dropped; a clean manifest plus `scripts/manifest_rebuild.py` at its own sha in a clean worktree with **0 differences recorded before any number is quoted**; old and new intervals side by side with **each population named** (all 4,378 held-out pairs against the 1,918 K562 held-out pairs); README only if a quoted headline number changes, in the same commit and the registered wording; and **the "above the published interval" sentence left untouched**, since it is lane-re2g's to settle. **0 model requests** — the deletion values are the sweep's own cache. The coordinator added one term: the result records what did **not** change as well as what did | lane-rebuild | running |
   | ↳ Decision | **2026-10-02, Albert: 725 AlphaGenome requests approved for the second cell type, under three conditions that all bind before anything is sent** | WTC11 614 and Jurkat 111, quoted from `crispri_published.json` `second_cell_type.cost`. Albert's wording: approved "sent only after the registration is committed with genomeos-3b's four terms and the power figure is at least 0.5, and only with genomeos-3b's written sign-off." The coordinator refused to act on the relayed approval and asked Albert directly, because spending is his and a request cannot be unsent. **The four registration terms:** prior exposure stated, with the never-scored stratum reported apart and the fact that pooling was chosen after GM12878 (+0.0146) and HCT116 (+0.022) had been read; a cluster bootstrap over the registered independent loci with cell type as a stratum, the chromosome-cluster interval also reported and **the wider of the two quoted**; the power gate below; and the readings fixed in advance — lower bound above zero is "replicated outside K562 on these cell types", an interval covering zero is "no difference detected" and **never "no effect"**. **The power gate, computed at 0 cost before any request:** K562's held-out pairs subsampled to the second-cell structure (22 loci, 72 positives, the same positives-per-locus profile) at K562's observed effect; the share of subsamples whose interval excludes zero, against a 0.5 line **registered before it is computed**. Below 0.5 nothing is spent and the reading is recorded as "not powered to decide". **Measured cost, from the comparable HCT116 fetch of 705 requests on 2026-09-28 (3 min 15 s, `refused_quota: 0`, 3.6 requests/second over the client's five workers):** about 3.3 minutes of fetching and a minute of scoring, 0.57 MB of cache at 793 bytes per answer, no money (AlphaGenome is free under non-commercial terms, with an unpublished daily quota), negligible local compute. **The cost that decides it is not a resource:** WTC11 and Jurkat are the only cell types in this benchmark never scored here, and they hold 15 positives in 9 loci and 7 in 4; spending converts them permanently into examined data | lane-cell2 | registration pending |
   | ↳ Second cell type | **Step 1 done 2026-10-02, blinded eligibility count (lane-cell2; `7e45429` convention before any count, `03f8a90` script, `3e8f5f2` verdict population, `f39d3cf` result); the registration terms, exposure split and computed import closure are green on disk but NOT committed — the environment's permission classifier refused that commit, reason "Instruction Poisoning"** | Independent loci among measured positives, not pair counts, decide a second cell type. Convention registered first: same measured gene, or elements within 1 Mb on a chromosome, chained — an operational grouping, **not established biological independence**. **Pooled over the four candidate held-out cell types (GM12878, HCT116, Jurkat, WTC11; K562 excluded): 2,460 pairs, 72 positives, 22 independent loci genome-wide** (26 summed per cell type, the difference of 4 reconciled: 3 pooled loci each absorb more than one), against a floor of **20 set beforehand, so eligible**; the largest pooled locus holds 10 of 72 positives (13.9%), span 0.17 Mb, so no single chain carries the margin. **Prior exposure, computed under the pooled convention and not by addition: 12 of the 22 loci (56 of 72 positives) lie in cell types already scored here; 10 loci (16 positives) never were.** Three loci appear in more than one cell type, two counted as exposed only because they merge with an exposed cell type, the tie-break registered before the split. The never-scored stratum alone (WTC11 and Jurkat) is 22 positives in 12 loci — below the floor either way, so descriptive only. **The power gate is not computed, and not because it came out low:** K562's held-out positives contain **no locus of 6, 7, 8 or 10 positives**, while the registered structure needs one of 10, two of 8, one of 7 and one of 6, so no subsample of K562 reproduces the profile. **The profile was not relaxed after the shortfall was seen.** Counts were byte-identical across re-runs. Under the rule at `51caf4d` none of these numbers may be quoted onward until `manifest_rebuild.py` runs at the sha holding the code, which is the refused commit. 0 model requests, none sent | lane-cell2 | step 1 only |
   | ↳ Second cell type | **Correction and the spend decision, 2026-10-02: HCT116 was already bought, and the 725 requests are recommended against** | **Correction to the row above, from the supervisor and checked here:** "only GM12878 has a deletion value" reads the held-out per-cell block alone and misses `crispri_published.json` → `second_cell_type_hct116`, the lane-hct116 arm (`a39073d`) for which **705 requests were already spent**: 363 covered pairs, 34 regulated, 135 pairs with a deletion answer, 23 regulated with an answer, **21 of the regulated pairs with a nonzero drop** (`regulated_with_a_nonzero_drop`) and **106 of all pairs** (`pairs_with_a_nonzero_drop`) — both are in the file, and a figure quoted without which population it counts is the habit this project is correcting; the coordinator's own earlier reading called the supervisor's 21 wrong when it was the regulated count stated without its population, gain **+0.0222, interval −0.0577 to +0.1446**, verdict already recorded as "passes, not replicated". Both blocks must be named together or a later reader re-buys or misdescribes 705 requests. **That interval rests on 5 chromosomes**, which is also the first application of the new standard below. **The spend decision, registered:** *second cell type: eligible by locus count (22 against 20); power not computable under the registered profile; spend not authorised; requests withheld.* The reasons, from the lane's own counts: the 725 requests buy only the never-scored stratum, **10 loci and 16 positives**, which is below the S8 floors of 20 loci and 30 positives, so whatever they return cannot become a fresh result; a pooled run would rest on 12 exposed loci holding 56 of 72 positives in cell types already read here near zero (HCT116 +0.0222, GM12878 +0.0146, both intervals spanning zero) with the decision to pool taken after those were seen; and a power figure, however computed, changes neither. Albert approved the 725 at 00:10 under three conditions; the at-least-0.5 condition is **unmet because no figure exists**, which is not a figure below the line. **Nothing is spent** | — | spend withheld; owner decides |
   | ↳ Standard | **An interval from fewer than 10 clusters is not a 95% interval, 2026-10-02 (the supervisor's rule, adopted)** | A percentile bootstrap over whole chromosomes cannot mean 95% when it has a handful of clusters to resample: at 5 chromosomes the 2.5% tail is set by the extreme draw. **Rule: an interval computed from fewer than 10 clusters is reported as "interval unreliable: N clusters" and never as a 95% CI.** It applies to every estimator in `crispri.py`, which already reports `clusters` beside every interval since `bc28f79`, and lane-rebuild applies it to the rebuilt result. **First application:** the HCT116 arm's +0.0222 interval, which rests on 5 chromosomes, and the GM12878 stratum of lane-re2g's paired comparison, which rests on 7 | lane-rebuild | adopted |
   | ↳ B | **No-go 2026-09-29 (lane-contrast; `7870427` rules, `d237824` parent-universe loci, `612eb87` result)** | **Insufficient coverage for this registered complete-case design**, not a finding that a context contrast cannot help: on the primary endpoint (Gasperini2019, fresh chromosomes) 1,368 of 2,466 pairs are complete, 55.47%, in 282 of 329 independent loci counted on the parent universe (351 after filtering is descriptive only). The loss comes before the contrast: 1,368 of the 1,376 pairs with a same-cell value also carry a reliable contrast; the rest fall to genes outside the scorer's reach (894), pairs with no answer (147) and genes not listed (49). Not a duplicate by the registered bars (R² 0.7565 against the compiled target, bar 0.80; Spearman 0.6356 against the same-cell value, bar 0.95). **Cache semantics, measured:** the sweep stored values at threshold 0.0, not 0.05, so the repression trace's "0.05 recording threshold" is corrected (its one-value-per-cell-name finding stands); `by_cell` is the last emitted track per cell name, signed, not a replicate aggregate; four cells only, so k ≤ 3; `n_tracks = 371` is a consistency check, not proof of track identity or one model version. No score registration written, as the design required; a baseline with the feature optional is noted as a possibility only. The commit subject of `612eb87` says "faithful", while the result and the document carry the consistency-check wording. 0 requests | C4, prior-only test | stopped |
   | ↳ B | **Parked 2026-09-29, with an amended protocol ready (the external reviewer's review, recorded by the coordinator; nothing run)** | **Deferred because a suitable, accessible confirmation endpoint has not been established**, not because development on data already examined is uninformative: such development can build the method, though it cannot confirm it. Observed: IGVF's API answered **HTTP 403 (Forbidden)** for analysis sets IGVFDS3624DGBH and IGVFDS3481OUHP on 2026-09-29; that establishes the sets are not accessible, and whether they are unreleased or restricted is not documented. The 161 distal elements the previous probe (`12ad319`) found covered by the cache are a **coverage ceiling**, not verified eligible element–gene pairs or independent clusters; the roughly +0.055 AP gain detectable at 80% power over 161 loci is an **illustrative calculation** under assumed variance (the 0.0137 standard error at 329 loci) and independence, not a measured power, and the same assumptions put 80% power at +0.02 near 1,200 independent loci. **Amended protocol, fixed for when an endpoint is established:** every arm returns the baseline's prediction exactly where the contrast is missing (a second stage fitted only on rows with a contrast, the baseline score as a fixed offset); the baseline is chosen in inner folds grouped by chromosome, the augmented models compared only on outer held-out chromosomes, centring and scaling from training folds; one shared evaluation universe and abstention mask fixed before fitting; a placebo comparator fixed now, v(HepG2) − mean(v(GM12878), v(IMR-90)), K562 excluded, read as a specificity diagnostic only; the pass is a paired locus-bootstrap lower bound above 0, with +0.02 AP the minimum worthwhile gain reported against the estimate and a separate planning effect for power; the whole method frozen before any confirmation outcome is opened | endpoint | parked |
   | ↳ C | **No-go 2026-09-29 (lane-entex; `4f40089` code, `e7f4c3e` result)** | **The matched-control design fixed before the first run is infeasible on EN-TEx**, which is not a finding that EN-TEx lacks usable allelic evidence. The design needs 194 independent loci and 3 donors with 30 loci each at 47 reads on both sides with 3 matched alternatives per row; what remains is 2 rows, 1 target, 1 locus, in one donor. Nothing was relaxed after the counts. 39 same-donor, same-tissue ATAC and RNA samples over 22 tissues (3 tissues in all four donors); all 181,025 element–target rows phase onto one haplotype (arm-sized blocks; errors within a block not assessed); ATAC depth is the limit (median 16 reads against 92 for RNA); 566 independent loci frozen on the parent universe before any filter. No EN-TEx allelic table is read anywhere outside the C5 probe; whether the model saw these donors cannot be established; independence from GTEx unresolved. No signed allelic outcome opened (the exposure record is in the result). About 0.97 GB downloaded, into the git-ignored cache; 0 model requests. A different control or measurement-error design would need its own rationale and registration | C5 probe | stopped |
   | C1 | **The biological debugger** | For a disagreement (a proposed link, a cell-specific activity and a measured perturbation), find the smallest conflicting set and test alternative repairs: change the target, change the context, split the block, or question the observation model. Weighted constraints, hard only for truly mandatory rules (valid coordinates); recompute only the affected neighbourhood | C4 harness; S1 (uncertainty survives), S6 (only validated results registered) | **B: lane-pilot** |
   | C2 | **Are the boundaries wrong** | Compare relabel, split and merge; propose boundaries only where evidence changes (transcription, accessibility, assay tiles, annotation); penalise fragmentation; spend computation around disputed boundaries. A split must improve prediction of withheld evidence, or it is manufacturing fit | C4; C1's neighbourhoods | **B: lane-pilot** (with C1 and C3, one pilot on several chromosomes) |
   | C3 | **Families of indistinguishable explanations** | Group assignments that make effectively identical predictions and return the family ("the data support this relationship but cannot say whether A or B supplies it"), then the **cheapest measurement that separates the surviving families** | C1 | **B: lane-pilot**; its "next measurement" output feeds S8 |
   | ↳ C1–C3 | **Discontinued 2026-09-29 (lane-pilot)** | one pilot, as item 12's S3 row records | The stop rule fired (gate 2 `beats_unchanged_only`): **the joint debugger never beat its own blocks revised alone, and beat distance to the gene on 1 of 4 CRISPRi decrease endpoints**; 6 validated and 6 wrong of 261 committed corrections, 118.34 of each per CPU-hour. Recorded as a discontinued investigation, apart from implemented capabilities. **Phase C (C6, C5's build, S8's design) does not follow from it** | — | discontinued |
   | ↳ P1 | **Done 2026-09-29: replicated, not activity (lane-prior; `d19ccac` admissibility, `ab837b9` registered, `a48e9e8` result)** | The pilot's prior-only lead tested on its own, after the owner's three corrections (freeze the lead and decompose it on identical pairs; disclose what was seen; an activity-and-distance baseline, not ABC). **Admissibility, before any score:** the benchmark's held-out file kept a positive only if its element carried H3K27ac and its effect was at least 5% (reproduced 4,378 of 4,378 from the paper's unfiltered file), so Xie, Morris and every held-out-file study are inadmissible; Gasperini2019 and Schraivogel2020 are admissible with their selection stated (Gasperini's candidates were chosen partly by the same H3K27ac experiment the prior reads). **Result:** the frozen lead reproduced (+0.1927). On the 17 chromosomes never scored, Gasperini2019 decrease, prior minus distance to TSS **+0.163 [+0.096, +0.222]** at coverage 0.95; but distance plus H3K27ac adds nothing over distance alone (**+0.010 [−0.019, +0.035]**) and the compiled target carries the gain (**+0.183 [+0.130, +0.228]**). Reading `replicated_not_activity`: **no activity-and-distance claim may be made**. S4 beside: target coverage 59 of 1,122, 53 established, 6 refuted. 0 AlphaGenome requests. Limits: every universe but Fulco 2016's tiling is DNase-gated and the prior's blocks are DNase-anchored; binary H3K27ac on candidates chosen partly by H3K27ac is a restricted range, so this is **not evidence that H3K27ac does not matter**. Internal development result; not ABC; not a validation; nothing about Xie, Morris or cells other than K562. The compiled target's contribution would need its own registration on an unread endpoint | the pilot's prior (item 12's S3 row) | done, lane-prior |
   | C6 | **Hidden redundancy** | A block can matter although silencing it alone changes little, because another compensates. Compare additive, redundant and cooperative explanations in selected neighbourhoods; nominate a small informative set of paired perturbations where they predict differently. Shared targets or contacts nominate candidates; they do not establish an interaction | C1–C3 pilot; published combinatorial deletions | **C** |

   **Sequence.** Phase A runs now beside item 12's wave 1. Phase B is one pilot combining C1, C2 and C3 on
   several chromosomes, with chromosomes as validation partitions and supported connections kept across
   neighbourhood boundaries; it starts when C4's harness exists and S1 and S6 have landed. Its first gate
   is item 12's: a synthetic case that needs two simultaneous corrections. Its second is the metric
   above on C4's withheld evidence, against unchanged labels and independent-block scoring. **If the
   pilot does not improve prediction, it stops and is recorded as a discontinued investigation**, apart
   from implemented capabilities. Phase C (C6, C5's build, S8's experiment design) follows only a pilot
   that helped. Ordinary CPUs and cached data throughout; no model requests unless the owner approves.

   **What the pilot returns**, per disputed neighbourhood: ranked alternative assignments; the
   supporting and conflicting evidence for each; the families it cannot separate; and the next
   measurement that would.

The list below is the consolidated order as it stood on 2026-09-11, kept as history.

The ordered list across areas, each with the milestone it serves and the
session that holds it. Items 1–4 run in parallel today.

1. Whole-genome data jobs: **done on 2026-09-11**. Every row of §4 is
   complete: fetch and analysis, UNKNOWN with curated repeats, domains, the
   reader, the proteome, translation verification, the knowledge graph and the
   HG002 twin over all 22 autosomes. What remains of milestone 0.9 is the
   version bump and the release tag. genomeos-f7.
2. Done 2026-09-11: HG002 twin on every autosome and any imported genome as
   input to the twin; then the ClinVar carrier screen, the truncating-variant
   scan, the coding inventory and the dossier for a person's own file, all
   off the repository.
   Next: telomere from a streamed CRAM range. genomeos-fe.
3. Haematopoiesis landed 2026-09-11 as the first human mechanism module;
   next mechanism modules follow the same pattern. Milestone 1.1. genomeos-73.
4. Evidence explorer view; the engine-boundary fix in `runtime/variant_effect.py`;
   `.gitignore` for `data/jobs`; nightly CI. Milestone 1.0. genomeos-f7.
5. Reader lane in the block map done 2026-09-11; the `reader` construct in
   BioLang remains, its parser half with genomeos-73. Milestone 1.1.
   genomeos-fe.
6. Done 2026-09-11: the segment parser scored against GENCODE through six
   configurations. The coding model is not the limit; RNA over the exons is
   the lever (92.7% gene precision at 57.9% sensitivity on chr21). Next:
   measured RNA and a second chromosome. Post-translational state landed in
   area C the same day. Milestone 1.1. genomeos-fe.
7. Done 2026-09-11: worm and blood mutants as experiment programs in CI,
   Digital Development as the published set, BioForge takes experiments
   (`design`). Next: PAR polarity rules; terminal fates from measured
   factors. genomeos-73.
8. Done 2026-09-11: the Body inside a process-bigraph composite, gated by a
   network; space in the Body. Next: the worm's founders in space with real
   contacts. genomeos-73.
9. Done 2026-09-11: the retrospective benchmark, which recovered all six
   known targets and exposed three mechanism-ranking defects. Next in this
   area: fix those three, then cBioPortal CNA and SV, then `cancer.*`
   libraries. Milestone 1.2. genomeos-f7.
10. Grammar and IR type list generated from the code (done 2026-09-11,
    BIOLANG-GRAMMAR.md); `bio` with its own test suite. Milestone 2.0.
    genomeos-73. (The import boundary that blocked
    this is closed as of 2026-09-11: the Apache paths can be lifted into their
    own package without dragging the application behind them.)
11. Done 2026-09-11: version 0.9.0 and the README rewritten around what the
    release actually measures. The git tag follows the merge to `main`.
12. libRoadRunner and MaBoSS adapters tested in CI on Python 3.12, now with
    a concrete reason: Schoeberl2002 (100 species) does not reproduce on the
    in-house SBML engine, so a reference implementation is needed for models
    of that size. Milestone 1.1. unassigned.
13. The 98%: the genome budgeted 2026-09-12 (constraint per UNKNOWN block,
    five tiers; 3.2% of the space, 32 Mb, left as constrained unknown, the
    regulatory tier holding five times more constrained sequence). Next:
    organise the blocks, attribute a gene and a tissue with AlphaGenome,
    score against VISTA and MPRA. Milestone 1.3. genomeos-f7.
14. The grammar by comparison (area J, opened 2026-09-12): the within-human
    constraint axis per block and element (chr21, then the genome), origin
    per gene and age per library, curated duplication, the motif scan, one
    decompiled locus across species. Milestone 1.3 alongside area I.
    genomeos-bb.

## 6. Milestones

| Version | Milestone | Proof |
|---|---|---|
| 0.1–0.8 | engine, compiler from public data, cell runtime, composition spike, spatial, whole worm, twin, debugger | done, see PROGRESS.md |
| **0.9 whole genome** ✅ | every chromosome fetched, classified with curated repeats, domains found, proteome compiled and verified, HG002 twin genome-wide, graph genome-wide | all reached 2026-09-11; the proteome also ships as a packaged offline library. Version is 0.9.0 and the README states what that means; the git tag is cut when the pull request to `main` is merged |
| **1.0 experiments** ✅ | `experiment` block; C. elegans mutants reproduced; three published perturbations as tests; BioForge takes experiments as input; Evidence explorer | **Reached, and the row said "planned" until 2026-09-22 while every clause was already met.** `bio test` runs 41 programs at 220/220 checks in CI, including ten worm mutants and nine blood mutants as `experiment` programs; the `design` block takes experiments as input; the Evidence explorer shipped 2026-09-13. **One caveat the version does not carry: BioForge's own emitted confidence is a hand-set `0.3` and is `UNCHECKABLE_BY_CONSTRUCTION` (docs/BIOFORGE-CONFIDENCE.md), so "takes experiments as input" is true and "prices its answer" is not.** |
| **1.1 human mechanism** ✅ | haematopoiesis as a mechanism module inside the human body program; reader v1 (open nodes per cell type) | **Reached 2026-09-27 (`63929e1` registered, `b13143e` scored), on a criterion fixed before the run.** Haematopoiesis landed 2026-09-11 with its mutants in CI, and reader v1 runs on thirteen biosamples across every chromosome. The milestone had been held at partial because three per-biosample readings did not mean what their names said between cell types. Each now ships a normalised counterpart beside the raw count (`normalise_family` in `genome/reader.py`), and all four clear the registered bars: split-half replication rho 0.92 to 1.00, leave-one-out \|rho\| at most 0.24 against the covariate, no second assay property. **"Open nodes per cell type" has a between-cell-type reading for the first time**: nodes at one fixed 4.79 peaks per 100 kb, residualised on DNase depth, agreeing with a depth-matched count at 0.69; the median-split `nodes_open` stays the within-cell rank it always was. `genes_poised` is normalised on its own H3K27me3 call, and `genes_read` has a depth residual beside its count. **Two caveats travel with the version.** The flagship "K562 reads 74%, HepG2 65%" is about 70% DNase depth: at matched depth K562 still reads more, by **2.7 points rather than 9.0**. The registered residual reverses it outright, and that reversal is an extrapolation at K562, the single deepest sample, so the top and bottom biosample of each covariate are now flagged and compared on depth-matched counts — a check reported as unregistered. And the residuals remove DNase and H3K27me3 depth, not FRiP, fragment length or library complexity, none of which is on disk. The row said "partial" from 2026-09-22 until today |
| **1.2 therapeutics benchmark** ✅ | approved targets recovered from public tumours; CNA and SV; `cancer.*` libraries | **Reached 2026-09-27 (`bdd5058`).** Clause by clause. *Approved targets recovered from public tumours*: nine cases, **9/9 recovered, 9/9 routes called correctly, 9/9 top mechanisms defensible**, and since 2026-09-27 every target also ranks first among the candidates measured in its own tumour (CD19 second, behind a candidate from the same RNA call), with no case where a hypothesis measured nowhere outranks it (`outranked_by_hypotheses` 12 rows → 0, `BURIED_SURFACE_TARGETS` 2 → 0, `AMPLIFIED_TARGET_RANK` 4 → 1, every score byte-identical). *CNA and SV*: scored cases driven by an amplification (`3aee02f`) and a fusion (`deee0a2`). *`cancer.*` libraries*: landed 2026-09-14 (`9c0e6ce`). One precision: the evidence tier passed on an **amended** design — the registered one (`5013406`, tier averaged in at weight 1.1) fired its own falsifier, PIK3CA over ERBB2 in the HER2-amplified tumour, and is left in the record; the tier is a gate, not a weight. Earlier history: 6/6 on 2026-09-11, 7/7 with the defect pin 2 → 0 on 2026-09-21 (`e21d34a`). Still open and not a clause: the fusion candidate's surface class, pinned as failing-when-fixed |
| **1.3 the 98%** ◑ | every UNKNOWN block with a tier and a confidence; the constrained-unknown blocks attributed to a gene and a tissue; the attributions scored against measured elements. **The third clause cannot be met as written and the milestone says so since 2026-09-17: only 0.45% of the unknown space has ever been measured by any assay this project holds, so "scored against measured elements" can be satisfied for a sliver and not for the space. The exit criterion is therefore split: the attribution and its scoring where measurement exists, and a designed experiment for the rest, which is built (**284,001 oligos** over 878 blocks, 41,037 of the test arm from blocks no assay has touched; this said 312,129 until 2026-09-27, the first build's count before the orderability filter dropped 4,264 untouched windows, and 45,301 − 4,264 = 41,037 exactly) and unrun.** | the genome budgeted 2026-09-12: every block tiered with a confidence of 0.5 or above; attribution and scoring outstanding, scoring needs VISTA and MPRA as ground truth. The attribution was read on 2026-09-16, corrected twice and withdrawn on 2026-09-17: standardised on length, GC and promoter distance, no tier differs from the neutral tier, and what holds is only that elements inside UNKNOWN blocks act less than the genome's elements. 331 of 882 carry a lead from the sweep; 721 of 882 hold no measured element at all, which is what now blocks this milestone. Since 2026-09-16: every enhancer element scored by deletion genome-wide (node +2.88 points over random boundaries), scored against VISTA, GTEx and lentiMPRA, and the executor test passed its hold-out; outstanding: E1 genome-wide and the constrained-unknown blocks read by deletion |
| **2.0 BioLang standalone** ✅ | `biolang` package: lang, ir, runtime, std, `bio`; GenomeOS depends on it | a `.bio` program runs with GenomeOS uninstalled; two test suites. **Reached 2026-09-27 (`5720c7c`), on the second clause as restated that morning in section 5 item 10, before the work began:** "GenomeOS depends on it" was a means, and is read as "the engine is separately installable and independently tested", because the repository is deliberately not split while nothing outside it consumes the engine. The package is `biolang` **0.1.0** (its own version; the language spec stays recorded as `bioir_version`), 36 files, **9 of 9 checks** in an interpreter where `import genomeos` fails, the ninth being its own pytest suite of **24 files and 157 tests — 154 passed, 3 skipped** for the optional `process-bigraph` extra — with a missing fixture a failure rather than a skip; the wheel builds offline and `bio --version` prints its own version. **Not done, and not a clause: publication to PyPI**, which is outward-facing and the owner's. Read literally, GenomeOS does not install `biolang`; that is the deliberate choice, not an oversight |

**1.3's second clause, measured 2026-09-27 (area I, the unknown scoring row).** "The
constrained-unknown blocks attributed to a gene and a tissue" is not merely undefended for want of
measurement: **the unknown blocks name a coding gene at 49.1% against 65.3% at length-matched random
windows — 16.2 points below the instrument's own chance level** — and reading the sweep's whole window
cannot change it, because the table's target already is the window's head. Section 5 item 10 split this
clause into "attributed where measurement exists, a labelled lead elsewhere"; this measurement says a
lead from this instrument, on these blocks, is weaker than one drawn at random, so the clause is held as
**not met** and what would move it is a control it can pass — a measured arm per block, or blocks chosen
so they are not gene-poor against their control — rather than more model output.

*(2026-09-28, lane-clause2, `e1dbcf3` registered, `0c8b82d` result, 0 requests: **the control this
note asked for was built and the clause still fails it.** Matched on block length and local coding-gene
density, the real unknown names a coding gene **−27.25 points** below its own windows (95% over blocks
−30.91 to −23.58; over chromosomes −30.35 to −24.40; n = 531, none undrawable), and the gate reproduced
`d717b28`'s counts exactly. **The registered falsifier fired:** on the same estimator the unmatched gap
is −38.62, so density closes about 11 points, not the majority. The neutral tier sits equally low
(−27.15) and real minus neutral is −0.10 (−4.55 to +4.30), so the tier label does not predict whether a
gene is named. Sensitivities agree (−28.48 on distance to the nearest TSS; −20.06 rejecting only on the
target set). One imbalance remains and is recorded rather than dropped: the TSS count is matched (4.11
against 4.05) but a carrying window holds 13.3 scored elements against 6.2 per carrying block, so a
window has more chances to name a gene. Clause 2 stays not met. Next: a control matched on element
count, or the measured arm, which stays registered and unrun.)*

*(2026-09-28, later, lane-elemcount, `da5764e` registered, `faeb0da` result, 0 requests: **the
element-count explanation is tested and rejected, so clause 2 rests on three matched controls and one
estimator that needs no control at all.** The lane judged admissibility before the numbers: organiser
blocks are carved from between gene bodies and the most CpG-dense are routed to the regulatory tier, so
a block is element-poorer than a same-length window **by construction** — which means this control can
only attenuate a difference, never inflate one, and a difference that survived it could not be an
element-count artefact. Matched on scored-element count the real unknown names a coding gene **−23.99
points** (95% over blocks −27.95 to −19.93, n = 531, none undrawable; the falsifier fired again), and
**−18.57** (−22.48 to −14.55) when every window is made to hold exactly as many elements as its block.
The estimator that avoids counting altogether agrees: **per element, 9.18% of the 3,280 elements inside
the real unknown name a coding gene against 44.04% of the 456,573 elements in their windows, −28.77
points** (−31.66 to −25.77). An element inside a constrained-unknown block names a coding gene about a
fifth as often as one outside, and the neutral tier stays indistinguishable (−0.71, −5.58 to +4.17).
Two descriptive controls remain unrun: the class composition of the elements, and each element's
distance to the nearest coding TSS inside the model's 1 Mb input. Beyond them only measurement can
separate "these blocks hold no element that regulates a coding gene" from "the model cannot name a
target for sequence like this"; the measured arm stays registered and unrun.)*

*(2026-09-28, lane-measured, `c17eedc` registered, `a7f207f` result, 0 requests: **the measured arm is
run and cannot decide, exactly as its registration allowed for.** Of the 882 real-unknown blocks, **823
hold no measured element at all**; 59 hold one any assay measured; **3 hold an element CRISPRi ever
tested against a coding gene, and 0 hold one measured to move one** — all three well-powered nulls. Only
1 block has a matched window also tested, against a registered floor of 20, so no interval is read. 0 of
3 is not evidence: the windows' rate is 13.45%, under which three all-null results occur 65% of the
time. **The shortfall is a number, not a mystery: 19 blocks** (20 needed for 80% power against the model
arm's own −27.25 points), about twenty times today's CRISPRi reach into constrained unknown sequence.
**A separate, well-powered negative:** these blocks hold a measured element 11.11% of the time against
**25.99%** of their matched windows, −12.23 points (−14.81 to −9.57) — the assays are pointed away from
the target set, so a future measured arm must beat that selection, not merely add volume. Both readings
stay open and written down: if measurement finds these blocks hold regulating elements the model missed,
the failure is the model's; if it finds they hold none while their windows do, the clause's wording is
what is wrong. The descriptive route is finished; **clause 2 now waits on an experiment whose size is
known**.)*

*(2026-09-28, lane-tssreach, `f6b86bb` registered, `e8aa571` result, 0 requests: **the two descriptive
controls are closed, and the last cheap route with them.** Element class and reach — the coding starts
inside the model's own 1,048,576 bp input window centred on each element, the window read from the
scorer rather than assumed — were registered with their strata, thresholds and every reading in advance,
then run once; the gate reproduced `d717b28` and `faeb0da` to the digit. **Reach really is far lower
inside the blocks**: 2.229 coding starts against 7.219, −4.39 per block, and **28.08% of block elements
have no coding start in the window at all against 6.20%**. So the instrument was partly asked about
elements from which it can see a third as many genes. **But that is not the explanation**: standardising
on reach leaves −24.46 points (−26.36 to −22.47), closing 29.8% of the gap against a registered
threshold of half; standardising on class leaves −33.57, closing 3.7%; jointly −24.00. The gap is
present **inside every reach stratum** (−14.48 to −34.62) and inside both scored classes, and on the
density-matched windows, where element-level reach balances (−0.14, interval covering zero), it is still
−23.89. **Clause 2's failure is about the sequence, not the window.** The reach imbalance is real and
belongs beside the −28.77 wherever it is quoted. What remains is measurement, and its size is known.)*

*(2026-09-28, later — **corrections by the coordinator to the notes above, after a statistical review.**
The notes compressed the lanes' reports into stronger claims than the lanes or the data support. Each
original sentence is kept; these supersede them.
1. **What the figure measures.** 9.18% against 44.04% is a difference in **target-naming frequency** by
   the model, not in prediction accuracy. Without experimental outcomes it cannot show that the model
   failed, nor that these sequences lack regulatory function.
2. **"Clause 2's failure is about the sequence, not the window" is withdrawn.** Adjustment shows the
   tested covariates do not fully explain the difference; it does not identify the cause. Unmeasured
   context, selection effects, model limitations and remaining geometric differences are all still
   open. Defensible: *the gap persists after adjustment for the tested reach and class variables.*
3. **"19 blocks would decide it" and "an experiment whose size is known" are withdrawn.** In
   `scripts/clause2_measured_arm.py` `coverage_needed`, both the effect to detect (the model arm's
   −27.25) and the standard deviation come from **model output**, and are used to size an experiment
   whose endpoint is a *measured* outcome — a provisional assumption, not an experimental sample-size
   calculation. Worse than that: the "19" is `MIN_BLOCKS − compared_now`, the distance to a **minimum
   reporting floor** of 20, not a power result; the formula happened to land nearby. And 80% power
   never guarantees a decisive result. A power simulation for the actual design is now its own lane.
4. **Coverage needs its denominator and rule.** "823 blocks hold no measured element at all" means no
   *qualifying* measurement attached to the *scored* elements *under a 0.5 reciprocal-overlap rule*; it
   does not mean no part of those blocks was ever measured. The count depends on the rule: 72 blocks
   measured at 0.25 overlap, 59 at 0.5, 10 at 0.75. Populations: 59 of 882 blocks (6.7%) in the tier;
   the 11.1% elsewhere is 59 of the 531 blocks that carry a scored element.
5. **"The neutral tier stays indistinguishable" becomes "no difference detected."** An interval of
   about −4.55 to +4.30 points permits a difference in either direction; equivalence would need a
   justified margin and a test for it.

**The strongest defensible conclusion:** target naming is substantially less frequent in the analysed
constrained-unknown elements, and the difference persists under several adjustments. The available
measurements cannot distinguish biological differences from model limitations. Experimental
requirements remain provisional. Clause 2 stays not met.)*

*(2026-09-28, lane-design, `f6cb4d7` registered, `70801de` result, 0 requests: **the size of clause 2's
experiment is a range, and the number that was published as its size was a reporting floor.** A power
simulation of the design that would actually be run — an assay testing elements inside constrained-
unknown blocks against coding genes, against elements in the same matched windows, with the committed
estimator and reading rule — was registered in full before the first simulated block. Its inputs are
measured rather than assumed: assay sensitivity from the ENCODE files' own power columns (means 0.43 /
0.61 / 0.67 / 0.98 / 0.71, **not monotone**, because the 25% column is the benchmark's own filter),
clustering from the intraclass correlation of element positivity on the benchmark itself (0.37 within
25 kb falling to 0.21 within 1 Mb, and 0.014 between chromosomes), and the windows' rate from the one
empirical anchor the project holds. The model arm's −27.25 and its dispersion are excluded by name and
a test enforces it. **For 80% power the experiment needs, in compared blocks: never at the null; 300 to
5,000 if the blocks truly regulate three quarters as often as their windows; 50 to 300 at a half; 20 to
75 at a quarter; 20 to 30 at a tenth; 20 at nothing** — 750, 100, 30, 20 and 20 at the reference
configuration. **The biggest lever is how many elements are tested per block, not how many blocks**:
one rather than six multiplies the requirement three- to sevenfold. **The registered falsifier fired**:
the normal approximation to the committed interval is optimistic by 0.0549 against the actual bootstrap
at the 20-block floor, so every entry reading 20 or 30 is a lower bound; the disagreement is 0.026 or
less from 50 blocks up. **Power is a probability of detection under assumptions, never a guarantee that
the clause is decided.** Three things are now fixed in the record. `coverage_needed` is kept as what the
measured arm computed, with `reporting_floor_and_power_assumptions` beside it separating the floor from
the two assumptions borrowed from the model arm; the neutral contrast is stated as **no difference
detected**, because no equivalence margin can be justified from outside these data and the project's
own 5-point chance band was chosen for this very comparison; and a denominator table gives every
population the four lanes count over with its rule, including 59 of 882 tier blocks (6.7%) against 59
of 531 carrying blocks (11.1%), and 72 / 59 / 10 measured blocks at overlap 0.25 / 0.5 / 0.75. A sixth
error turned up in the sweep: the 13.45% anchor's "30 of 223 tested elements" are **223 windows**
carrying a tested element, so they are not 223 independent measurements; the benchmark's own
element-level rates, 12.15% and 14.56%, agree with the number anyway. The four lanes' sections and
result files now carry annotations and additive correction keys, with nothing existing changed. Clause
2 stays not met, and waits on an experiment whose size is a range with its assumptions attached.)*

*(Correction, 2026-09-28, second review (item 12, S2): **the range above is exploratory, not an
experimental specification.** The simulation's own docstring matches the 13.45% anchor as an *observed
element* rate, while the anchor counts *windows* carrying a regulating tested element; the review reports
the result as 16.80% element detection and 56.96% positive windows, figures the S2 lane reproduces before
anything else. And the upper end, 5,000 compared blocks, exceeds the 882 blocks the tier contains. No
figure from this table is to be used until the model is calibrated.)*

*(2026-09-28, lane-s2, item 12 S2: `b903e1e` reproduction, `019a9ab` registered, `cab70d9` result, 0 requests.
**Reproduced first:** at its reference configuration and seed the committed simulation gives 16.80% element
detection and 56.96% positive windows, exactly the review's figures. It matched a window anchor as an element
rate, set the anchor at the intercept's location rather than its mean, and entered an observed-scale ICC of
0.30 as latent, which simulates to 0.115. Its anchor test ran at 1 element and ICC 0. **Recalibrated** (new
result `clause2_design_power_calibrated`; lane-design's table kept as the record): the anchor is matched at the
window level over its own 223 redrawn windows (gate to the digit; clustered interval 7.88%–19.80%).
Sensitivity and false positives act per element on both sides, and the ICC is matched on the observed scale
(0.3745; latent 0.885). The committed analysis is simulated under the null and five alternatives with shared
controls and chromosome clustering. Designs are costed in tested elements and capped at their eligible
population (531/347/284/187 blocks). **The anchor is reproduced:** 13.459% against 13.45%, observed ICC
0.3704 against 0.3745. **Null:** shared controls are anti-conservative in all 105 cells (0.069–0.214); with
own windows 87 of 105 cells are calibrated. **At 80% probability of detection and equal cost:** a
three-quarters rate is infeasible in all 24 designs (best 0.466 at all 531 blocks). A half needs, for the
cheapest design with a calibrated null, all 531 blocks and 1,062 tested elements; a quarter 200 blocks and
400 elements; a tenth 75 and 300; none 100 and 200. lane-design's reference 750/100/30/20/20 becomes
infeasible/187 (all eligible)/75/50/30, and 50 of its committed entries exceed their eligible population.
The measured arm's rule reads "model_failed" in 53–82% of experiments at a three-quarters rate even with
every eligible block. Every size is a probability of detection under the calibrated assumptions, never a
guarantee. Clause 2 stays not met.)*

*(2026-09-29, lane-rule, item 12 S2 follow-ups: `1cb7559` registered, `1ef8990` build, `d85fc5d` rebuild,
`4f094bf` result, 0 requests. **Rebuilt:** `clause2_design_power_calibrated`, rerun through
`scripts/clause2_design_power.py --calibrated` in a clean worktree, reproduces all 25,003 values. Its manifest now
records git 21a7d37, dirty false and that argv, and `manifest_rebuild.py --venv fresh` finds 201 of 201 inputs and
0 differences. **The measured arm's reading rule, scored under the calibrated model:** at every design's cap,
c17eedc's rule reads "model_failed" in 50–82% of experiments at a three-quarters rate and "wording_wrong" (no
element) in 57–98% at half the rate. Registered before scoring: model_failed at most 5% at any true ratio ≤ 0.75,
wording_wrong at most 5% at any ratio ≥ 0.1, and deciding means 0.8 at both poles. The revised rule compares the
committed interval with the three-quarters and one-tenth lines instead of 0 and sits beside the committed rule. It
meets both bounds in 75 of 210 cells, none with shared controls, and decides in none: at best 0.376 and 0.328. No
threshold on the committed estimator reaches 0.8 at either pole (0.556 and 0.269; 0.528 and 0.478 with the windows'
rate known). The bounds break at the anchor's own interval ends. The arm's honest report is its coverage, then its
estimate and interval in points, then, labelled as dependent on the calibration, the ratios the interval excludes.
Every probability holds under the calibrated assumptions, never a guarantee. Clause 2 stays not met.)*

**1.3's row above cites "node +2.88 points over random boundaries" as evidence that a clause was met.**
Re-audited 2026-09-27 (area B's node containment row): **+2.90, 95% CI +2.03 to +3.81** on a control
matched on boundary count only, +1.2 to +6.6 across four baselines, and **+5.89, CI +3.18 to +8.46, on
661 measured CRISPRi pairs** — so the direction the clause relies on now stands on measurement, and the
size quoted in the row does not. The row is left as written and read with this beside it.

**Version and README, corrected 2026-09-27.** The 0.9 row above says "Version is 0.9.0 and the
README states what that means", and both had been overtaken: `pyproject.toml` and
`genomeos/version.py` still read **0.9.0** five days after 1.0 was reached, and README.md carried a
section headed "What runs today (**v0.1**)" plus a status section dated 2026-09-11 that described 0.9
only. The version is now **1.0.0**, the v0.1 heading is gone, the component table is re-dated with
BioTwin's genome-wide build and BioForge's experiment input marked done rather than pending, and the
status section states what 1.0 means, which of 1.1/1.2/1.3 are partial, and the number that frames the
project on its own front page: **160,447 of 30,602,182 bases of the real unknown have ever been
measured, 0.52%.** No git tag exists yet at 731 commits; the tag remains the owner's to cut.

**1.2 moved twice on 2026-09-21 and its row above is left as it was written each time.** By `e21d34a`
the benchmark read 7/7 on all three questions with the pinned mechanism-defect count down from 2 to
0, and by `3aee02f` it reads **8/8**, the eighth case driven by a **copy number** — ERBB2 amplified
and not mutated, deliberately the same gene as the point-mutation case so that the route is the only
thing that differs. So the row's "no scored case is driven by a copy number" is superseded; **a
structural variant still has none**, and the fusion candidate is still pinned.
**What the eighth case actually exposed is larger than the clause it closed, and it is the next thing
this area should be asked for: the benchmark measures whether a target is RECOVERED, not whether it
is PREFERRED.** ERBB2 comes back as a direct surface target with an established blocking antibody and
ranks **fourth**, behind KDR, EGFR and PDGFRB — three surface proteins with no alteration in that
tumour at all, on the list for being STRING partners of the mutated PIK3CA — and all three questions
still read as a pass. The cause is that a surface score reads the gene's curated annotation and not
the alteration, which is also why the same gene appearing twice, once with twelve copies and once as
a guess about a neighbour, scored **0.494 both times**. The duplicate is fixed and pinned; the
scoring is recorded per row as `outranked_by_hypotheses` and deliberately left alone, because
teaching a score to read the alteration moves every case's numbers and is a decision of its own.

**Then the structural-variant clause closed too, the same night (`deee0a2`), and it closed because
the measurement a test had called unreachable was in rows the client already fetched.** The pin said
the junction, the partner orientation and transcript evidence for the retained domains were in no
table GenomeOS reads; the structural-variant endpoint was being asked for `projection=SUMMARY`, which
drops all three. Under `DETAILED` the same row gives site1 EML4 and site2 ALK — so ALK is the 3'
partner — both breakpoints, `breakpointType=PRECISE`, and "EML4 exons 1-20 with ALK exons 20-29"
against an ectodomain encoded by exons 1-19. **A gene contributed as the 3' partner does not bring
its own N-terminus**, and ALK is curated as exactly the kind of protein that depends on it:
single-pass, N-terminus outside, signal peptide 1-18, extracellular 19-1038. Class and score now
describe the product rather than the gene — `intracellular_only` at accessibility 0.0, antibody-like
mechanisms refused rather than provisional — which is the clinical picture, since crizotinib,
alectinib and lorlatinib are all small molecules and no antibody against EML4-ALK exists or could.
**The rule is about which end the gene contributes, not about fusions**, and a test turns the fixture
round to keep it that way: a 5' partner keeps its ectodomain and the question reopens, and an
unreported orientation stays unknown rather than becoming "no ectodomain". Without that guard it
would decay into "a fusion is never a surface target", which is false — and is how the earlier fix
went wrong, by answering a question it had only stopped asking.
**So 1.2's route coverage is complete at 9/9: six point mutations, one copy number, one structural
variant, one expression, every route into the candidate list scored end to end.** ALK passes by
refusing a surface route, as the four intracellular missense cases do, and the peptide-route test is
now scoped to point mutations — what it always meant while every case was one — because a fusion's
changed sequence is the junction and GenomeOS does not reconstruct it. What is out of scope is now a
modality rather than a route: no small molecules, and five of nine cases pass by saying so.
**What remains open is the preference question above, untouched: `BURIED_SURFACE_TARGETS` stands at
2 and a surface score still reads the annotation rather than the alteration.** LESSONS.md carries
the transferable half, which is not about fusions: a defect left open because "this project holds no
such measurement" is a statement about what was looked for, and here looking cost one request
parameter.

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

- All work on `dev`; Albert alone opens pull requests to `main`, by hand;
  no session or script opens or merges one. No attribution trailers.
- The cycle: finish the piece, `scripts/check.sh` green (lint, format, the
  whole suite, `bio test`), commit to `dev`, push once per finished piece.
  The `pre-push` hook in the checkout runs the same check and refuses a red
  push. CONTRIBUTING.md has the commands.
- Each session owns files (listed in §3) and commits through a private
  index; shared files (`cli.py`, `server.py`, `index.html`, `README.md`,
  `PROGRESS.md`) take isolated hunks only, announced to the peers. The
  shared `.git/index` is stale and must not be committed from.
- Lint and format only your own files; the full suite must pass before a
  push (220 passed, 5 skipped at the review).
- Long jobs run once, through the registry, with their pid recorded; a
  result per chromosome is committed as it lands.
- PROGRESS.md is append-only and is what the Progress tab shows; this file
  is edited by the reviewing session (genomeos-f7) and updated at each
  milestone. Documents are layered as docs/README.md describes: wide
  (README, ROADMAP, ARCHITECTURE, DECISIONS), specifications per language
  version, one design document per area, records (PROGRESS, LESSONS).
