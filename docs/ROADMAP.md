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
boundaries by +2.88 points, no unknown tier can be ordered against another once the
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
  `docs/ARCHITECTURE.md` §10.1, milestone 2.0): 31 files, Apache-2.0, no
  dependencies, run in a Python with no site-packages where `import genomeos`
  fails outright — all four verbs, every module, the standard library and an
  organism, with the application absent.
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
- **Next (epigenome and nodes).** 1. Done as an alternative; the next test is a
  stricter site call (best-hit strand, stronger motifs) judged on all four
  measurements at once, and whether Hi-C questions and enhancer-to-gene
  questions need two node sets. 2. The poised call in the Blocks lane:
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
- **Data.** `proteome_chr*` (4 of 25: chrM, 21, 22, Y),
  `translation_vs_uniprot_chr*` (3), `graph_chr21`, `graph_genome`,
  `ptm_genome_wide`, compiled definitions in
  `data/knowledge/proteins` (18 MB).
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
- **Owner.** genomeos-d3 (from 2026-09-15; genomeos-d2, genomeos-d1 and
  genomeos-73 before).

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
  program, whose mean confidence is 0.27. Measured and predicted are almost
  equal in number, which is what a project that compiles public data and then
  predicts on it should look like; the weak half is concentrated in exactly
  one place, area I's attributions, rather than spread through the
  hand-written biology.
- **Missing.** One view promised in the UI track is still unbuilt:
  **Cell** (a cell type's active rule set and graph, run and watch). No git tag and no PyPI package; no nightly
  real-data CI; no "known phenotypes it must reproduce" list per library from
  a biologist. Done 2026-09-11: version 0.9.0 in `pyproject.toml` and the
  README rewritten around what the release measures. Done 2026-09-11: the jobs
  registry de-duplicates by pid liveness, so server restarts no longer spawn
  copies of a job. Done 2026-09-11: `data/jobs` metadata and logs are
  git-ignored, the registry kept in code.
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
  decision. Still one scheduled run a day. 5. Cell view: **done 2026-09-17**. The reader has run on eleven cell types over
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
- **Code.** `forge.py`, `design.py`, `runtime/debugger.py`.
- **Design.** DESIGN-MINIMAL-CELL.md, ACTION-PLAN.md Phase 5.
- **Data.** `design_neuron`.
- **Missing.** No calibration of the predicted confidence against outcomes;
  knockouts of network species (the `experiment` block perturbs organisms;
  network perturbations are covered by tests: the repressilator's `# test:`
  claims and the Fauré cell-cycle knockout in `test_boolean.py`). Done
  2026-09-11: the `design` block (BioForge takes experiments as input;
  knockouts, additions and knobs searched toward targets under constraints,
  the answer emitted as a predicted experiment; POP-1, PIE-1, GFI1 and IKZF1
  found), eight worm mutants and nine blood mutants as experiment programs
  in CI, Digital Development as the published set.
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
  chromosomes (sign test one-sided p 0.011)**. Real and small; the headline 90.2% of
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
- **Missing.** The motif scan and
  operator patterns; phylogenetic profiling between libraries; the
  decompiled locus view; one developmental locus run end to end.
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
| All-elements deletion scoring (every ENCODE enhancer inside a node, AlphaGenome) | 24 of 24 (2026-09-16) | done, `scripts/enhancer_targets_all_chain.py` | 961,227 elements, 778,780 requests; 63.7% name a gene (matched random windows 87%); the node beats as many random boundaries by +2.88 points over 440,377 coding-target elements, ahead on 18 of 24 chromosomes. Summaries committed, element tables local and packed. See area I |
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
| **1.0 experiments** | `experiment` block; C. elegans mutants reproduced; three published perturbations as tests; BioForge takes experiments as input; Evidence explorer | `bio test` passes the mutant programs; benchmark tests in CI |
| **1.1 human mechanism** | haematopoiesis as a mechanism module inside the human body program; reader v1 (open nodes per cell type) | lineage choices and counts reproduced with confidence above "low" |
| **1.2 therapeutics benchmark** ◑ | approved targets recovered from public tumours; CNA and SV; `cancer.*` libraries | benchmark built and in CI 2026-09-11: 6/6 targets recovered, 6/6 routes correct, 4/6 top mechanisms defensible after three fixes it prompted. **2026-09-21 (`e21d34a`): 7/7 recovered, 7/7 routes, 7/7 defensible, the pinned mechanism-defect count down from 2 to 0, and the seventh case reached by an expression call rather than a DNA event.** The `cancer.*` libraries landed 2026-09-14 (`9c0e6ce`). Outstanding: **no scored case is driven by a copy number or a structural variant** — both routes have unit controls only — and the fusion candidate keeps a surface class its product may not have, blocked on a junction or orientation measurement this project does not hold and pinned as a failing-when-fixed test |
| **1.3 the 98%** ◑ | every UNKNOWN block with a tier and a confidence; the constrained-unknown blocks attributed to a gene and a tissue; the attributions scored against measured elements. **The third clause cannot be met as written and the milestone says so since 2026-09-17: only 0.45% of the unknown space has ever been measured by any assay this project holds, so "scored against measured elements" can be satisfied for a sliver and not for the space. The exit criterion is therefore split: the attribution and its scoring where measurement exists, and a designed experiment for the rest, which is built (312,129 oligos) and unrun.** | the genome budgeted 2026-09-12: every block tiered with a confidence of 0.5 or above; attribution and scoring outstanding, scoring needs VISTA and MPRA as ground truth. The attribution was read on 2026-09-16, corrected twice and withdrawn on 2026-09-17: standardised on length, GC and promoter distance, no tier differs from the neutral tier, and what holds is only that elements inside UNKNOWN blocks act less than the genome's elements. 331 of 882 carry a lead from the sweep; 721 of 882 hold no measured element at all, which is what now blocks this milestone. Since 2026-09-16: every enhancer element scored by deletion genome-wide (node +2.88 points over random boundaries), scored against VISTA, GTEx and lentiMPRA, and the executor test passed its hold-out; outstanding: E1 genome-wide and the constrained-unknown blocks read by deletion |
| **2.0 BioLang standalone** | `biolang` package: lang, ir, runtime, std, `bio`; GenomeOS depends on it | a `.bio` program runs with GenomeOS uninstalled; two test suites |

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

---

## 7. Improvement and creativity pool

Ideas judged worth keeping, not yet scheduled. Each would be judged by the
ambition in §1 before it is started.

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
