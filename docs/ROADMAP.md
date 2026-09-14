# GenomeOS: scope, architecture, areas, plan

Reviewed and reorganised on 2026-09-11 against the code, the tests, the
result files and the running jobs. This is the one document that says what
GenomeOS is, how it is built, what each area is for, what is done, what is
missing and what comes next. Evidence per finished task stays in
[PROGRESS.md](PROGRESS.md); the original phased plan is kept as history in
[ACTION-PLAN.md](ACTION-PLAN.md); lessons the code carries are in
[LESSONS.md](LESSONS.md).

Measured state on 2026-09-11 (evening): 85 commits in two days, 44 `genomeos`
commands plus the `bio` toolchain, 16 web views, 255 tests passing and 4
skipped, lint clean, three `bio.std` modules, seven organism programs, every
human chromosome fetched and analysed, the whole human proteome compiled and
verified against UniProt, a genome-wide knowledge graph, the whole worm grown
from one cell, a human body grown as populations, haematopoiesis as mechanism,
a therapeutic design dataset with optional peptide/HLA binding, and the engine
separable from the application under its own licence.

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
| **The 98%** (GenomeOS) | genome budgeted | every UNKNOWN block of the 24 chromosomes tiered with Zoonomia constraint read per base (1,009 Mb): structural 23%, fossil 33%, regulatory 34%, neutral 7%, constrained-unknown 3.2% (1,098 blocks, 32 Mb, 1% of the genome); `genomeos budget`, the Progress tab card; the attributions compiled to a BioLang program per chromosome (`--bio`); scored against VISTA, eQTL, lentiMPRA, GWAS and ClinVar; the organiser reads copies first and leaves a real unknown of 882 blocks, 69 of them constrained on both axes | AlphaGenome chromatin tracks on the 15,536 regulatory blocks; closure at gene, cell and organism level |
| **Cancer and therapeutics** (GenomeOS) | 1.0 of the design dataset | tumour-only pipeline, cohort expression, altered-protein reconstruction, 15 mechanisms, design dataset with negative set | CNA/SV profiles; benchmark against approved targets; peptide/HLA predictor |
| **Web UI and CLI** (GenomeOS) | 17 views, 45 commands | every layer visible with evidence pills; jobs with progress; gene dossier; the Progress tab says what is going on and what is planned; the Evidence explorer reads the project's own confidence (22,419 facts, mean 0.57) | Cell view; release tags; nightly CI |

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
- **Next.** 1. `bio` packaged as an extra entry point with its own test set.
  2. A `population` type in BioIR (today a counted `Cell`). 3. PAR polarity
  rules for the worm's first divisions, so par-2 and par-3 knockouts are
  predicted rather than stated.
- **Owner.** genomeos-c1 (the engine: language, IR, runtime, toolchain, v0.4).

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
  nearest TSS. The node model has now been tested on the whole genome rather
  than argued for: the boundary predicts where an element acts nine times out
  of ten. Two findings worth as much as the headline: the nearest-gene
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
- **Fold-change profiles.** chr2, chr7, chr12, chr17, chr21 and chr22 carry
  them for every cell type and mark; every other chromosome carries peaks only,
  and the record says so.
- **Next (epigenome and nodes).** 1. An orientation-aware boundary caller
  (reverse-to-forward flips among motif-bearing CTCF sites) as an alternative
  node set, compared on the same measured boundaries and on the enhancer
  deletions before any switch. 2. The poised and read-by-marks calls in the
  Blocks lane, whose file belongs to area B: today it shows poised genes as
  silent through `silent_genes` but has no poised colour. 3. A neural, gonadal
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
  decide the embryonic fates they can (below); nothing in the language
  models commitment, competence or hysteresis; the PAR polarity rules
  are not written (par-2 and par-3 are the misses against Digital
  Development); the external engines (libRoadRunner, MaBoSS, CompuCell3D)
  are deferred for lack of Python 3.14 wheels, although CI runs 3.12 and
  could test them as optional extras. Done 2026-09-11: haematopoiesis,
  the Body as a process-bigraph process gated by a network, space (sites,
  fields, division in place, contact inhibition, migration; the French flag
  grown from one cell), the measured reader, Packer 2019 at 84% on
  single-tissue ids (the rest is labelling depth), Digital Development as
  the published knockout set (7 of 11 modelled transformations).
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
  decisions fire.
- **Next.** 1. When genomeos-c1 lands `commitment` and the exposure and
  mean reads (v0.4 §7.2a), rewrite `fates.bio` and `exposure.bio` with them
  in place of the generated `_integrated` factors, and score again. 2. PAR polarity rules for the first divisions.
  3. The worm's founders in space with the real contacts replacing named
  senders (needs the contact runtime requested from genomeos-c1), with the
  AC/VU pair scored as an equivalence group. 4. Glia from factors: sheath
  and socket against their sister neurons. 5. libRoadRunner and MaBoSS
  adapters tested in CI on 3.12.
- **Owner.** genomeos-d1 (from 2026-09-14; genomeos-73 before).

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
  distinctive claim. A third question is scored apart and the pipeline fails
  it, though less than it did: 4 of 6 now, after two fixes the benchmark
  prompted. What remains is that a cytoplasmic protein (BRAF, PIK3CA) still
  heads its list with a mechanism needing an extracellular epitope, now capped
  at the provisional ceiling of 0.25 and labelled with its unmet requirement,
  so it misleads nobody but should not be first. Pinned at 2, may fall, must
  not rise.
- **Three fixes the benchmark prompted.** (1) A small-molecule precedent
  counted as evidence that a mechanism needing an extracellular epitope could
  reach the target, so a kinase inhibitor "supported" a radioligand against
  cytoplasmic BRAF; such precedent now requires positive evidence of
  reachability. (2) An agonist antibody was ranked first against an activating
  driver, which would push the pathway the tumour already over-drives;
  agonism now needs positive evidence that triggering the target is the
  intent, and is refused on ignorance rather than offered on it. (3) A
  mechanism whose hard requirement merely went unanswered could outrank one
  whose requirements were established; such a mechanism is now capped at 0.25
  and carries the reason. All three are the same mistake in different
  clothes: treating an unanswered question as permission.
- **Missing.** Copy number and structural variants from cBioPortal profiles
  (only patient-supplied CNV today); a `cancer.*` library layer; the three
  mechanism-ranking defects above; expression-driven targets such as CD19 and
  BCMA, which need patient RNA rather than a variant; the demo runs at data
  level 1 only.
- **Next.** 1. Retrospective benchmark on public cases with known targets;
  a test that fails if a known target drops out of the top ranks.
  2. cBioPortal CNA and SV profiles. 3. `cancer.*` libraries with driver
  evidence. 4. A `NeoantigenProvider` behind a licence-checked optional
  extra.
- **Owner.** unassigned since genomeos-b9 left; genomeos-f7 takes the
  benchmark.

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
  data. 5. Cell view.
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
  (genomeos-h1): the regional stage is fixed and chr22 replicates (above); next either the executor test on the 1,053
  early-replicating storage units with an eQTL or an MPRA allele pair (allele-dependent
  read-out, which needs GTEx by allele or in-silico swaps and so waits on quota), or the
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
- **Missing.** The human axis as a Blocks-tab lane; the motif scan and
  operator patterns; phylogenetic profiling between libraries; the
  decompiled locus view; one developmental locus run end to end.
- **Next.** 1. Done: the human axis over the genome and the case in the
  compiled `region` note; next a Blocks-tab lane. 2. Done: origin per gene and age per library, `orthologue_of` still to add
  as species counts per clade on the graph nodes. 3. `genomicSuperDups` done; paralogues from the Compara stream; `paralogue_of`
  and `orthologue_of` edges in the knowledge graph (genomeos-fe's `molecules/graph.py`, hunk announced first). 4. Done, with families and a GC-by-repeat null: promoters carry family
  enrichments per library and no pair logic; the pair search moves to
  enhancers (GATA plus T-box in heart). 5. Done at library level; next a gene-level profile with a shuffled-library
  null before any pair is called a dependency.
  6. Done: `genomeos decompile GENE --chrom C`; next a Blocks-tab card. 7. Done, and tested on a VISTA panel with matched negatives: motif sites held
  across species do not separate enhancers from inactive conserved sequence;
  next a different readout, AlphaGenome in-silico mutagenesis per element on
  the same panel (are the bases that matter the ones that are kept).
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
| Enhancer targets, predicted (AlphaGenome) | 24 of 24 chromosomes (4,800 elements) | done | 2,994 elements move a gene, 2,291 name a coding gene; **90.2% of those sit inside the element's CTCF node** (79.8% to 91.8% per chromosome) and 71.2% are exactly the nearest TSS; 1,157 behave as silencers. Needs a key; the daily quota ran out mid-run and the job waited and resumed rather than failing |
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
| **1.2 therapeutics benchmark** ◐ | approved targets recovered from public tumours; CNA and SV; `cancer.*` libraries | benchmark built and in CI 2026-09-11: 6/6 targets recovered, 6/6 routes correct, 4/6 top mechanisms defensible after three fixes it prompted. Outstanding: the two remaining ranking defects, CNA and SV, the library layer |
| **1.3 the 98%** ◑ | every UNKNOWN block with a tier and a confidence; the constrained-unknown blocks attributed to a gene and a tissue; the attributions scored against measured elements | the genome budgeted 2026-09-12: every block tiered with a confidence of 0.5 or above; attribution and scoring outstanding, scoring needs VISTA and MPRA as ground truth |
| **2.0 BioLang standalone** | `biolang` package: lang, ir, runtime, std, `bio`; GenomeOS depends on it | a `.bio` program runs with GenomeOS uninstalled; two test suites |

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
