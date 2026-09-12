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
| **BioTwin** (individual) | 0.9 | any GRCh38 genome imported locally and used by carrier lookup, the gene report, the twin and the gene-by-gene walk; HG002 across all 22 autosomes with measured telomere and epigenetic age; fork, run, diff; predicted effect of a person's own regulatory variants; ClinVar carrier screen, genome-wide truncating-variant scan, coding inventory by consequence and a one-page dossier, all off the repository; the GIAB trio read for inheritance, 96.4% of HG002's calls inherited | telomere from real reads; true de novo variants out of the trio; polygenic scores |
| **BioForge** (design) | search only | random-restart search under constraints with predicted labels; minimal-cell design estimate | experiments as input; published perturbations reproduced; constraint language |
| **Genome decoding** (GenomeOS) | 0.9 reached | every chromosome fetched and analysed (2026-09-11): inventory, UNKNOWN classified genome-wide with curated repeats (98.5%), domains on 24, curated repeats and ENCODE elements on 25, reader on 25 with eleven cell types; the node model tested genome-wide (90.2% of enhancers act inside their node) and on two mouse chromosomes; the segment parser on chr21 at 92.7% gene precision once RNA over the exons counts as evidence | the parser on every chromosome with measured RNA; enhancer targets and boundaries from Hi-C |
| **Molecules** (GenomeOS) | 0.9 reached | the whole human proteome compiled from seven public databases (19,478 coding genes, 25 chromosomes) and packaged as a 2.3 MB offline library, translation verified on 19,249 of them, genome-wide knowledge graph, RNA layer with GTEx, pathways as reachability and as kinetics where a curated ODE model exists; 96,362 modifiable sites with their 352 writers as `modifies` edges in the graph; the dominant isoform per tissue from GTEx (APP695 in the cerebellum) | measured modification state; isoform expression per cell type or stage; an engine that handles 100-species models |
| **The 98%** (GenomeOS) | genome budgeted | every UNKNOWN block of the 24 chromosomes tiered with Zoonomia constraint read per base (1,009 Mb): structural 23%, fossil 33%, regulatory 34%, neutral 7%, constrained-unknown 3.2% (1,098 blocks, 32 Mb, 1% of the genome); `genomeos budget`, the Progress tab card; the attributions compiled to a BioLang program per chromosome (`--bio`), chr21 with 747 entities and 155 rules passing its own checks; scored against VISTA, eQTL, GWAS and ClinVar | AlphaGenome chromatin tracks on the 15,536 regulatory blocks; closure at gene, cell and organism level |
| **Cancer and therapeutics** (GenomeOS) | 1.0 of the design dataset | tumour-only pipeline, cohort expression, altered-protein reconstruction, 15 mechanisms, design dataset with negative set | CNA/SV profiles; benchmark against approved targets; peptide/HLA predictor |
| **Web UI and CLI** (GenomeOS) | 16 views, 43 commands | every layer visible with evidence pills; jobs with progress; gene dossier | Evidence explorer and Cell views; release tags; nightly CI |

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
- **Next.** 1. `bio` packaged as an extra entry point with its own test set.
  2. A `population` type in BioIR (today a counted `Cell`). 3. PAR polarity
  rules for the worm's first divisions, so par-2 and par-3 knockouts are
  predicted rather than stated.
- **Owner.** genomeos-73 (parser, IR, body runtime, std); genomeos-fe
  (`bio.py`, central dogma).

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
  | + measured RNA, six ENCODE cell lines pooled | 143 | 65.9% | 62.9% | 62.9% |
  | predicted panel ∪ measured lines | 154 | 63.7% | 63.0% | 66.1% |
  | predicted panel ∩ measured lines | 71 | 71.6% | **97.2%** | 54.3% |

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
  to 62.9% of genes, and pays for it in precision, 89.1% to 62.9%, because
  measured transcription includes the non-coding transcription the model's
  gene-level tracks exclude. A model that knows genes plus a measurement that
  knows the cell is the honest filter. Results are
  `segments_chr21_predicted_sites_measured_<cell|panelN>`. A measurement
  lesson found by the closure test (area I) and fixed here on 2026-09-12:
  ENCODE labels a strand track by the read, and IMR-90's total RNA-seq reads
  antisense, so `MeasuredRna` now probes sixty first exons per strand in both
  orientations and swaps when the swapped one carries twice the signal (K562
  26.4 against 1.8, HepG2 19.7 against 1.4, GM12878 36.3 against 6.4, A549 29.8
  against 5.6, MCF-7 36.2 against 2.2, IMR-90 0.85 against 22.4, swapped). The
  six-line panel numbers stand as five lines' worth until the rerun; SK-N-SH
  and HeLa-S3 have no total RNA-seq bigWig on ENCODE.
- **The two filters combined (2026-09-12).** `--rna-combine union|intersection`
  applies when both the predicted panel and the measured lines are given.
  The union is the most sensitive call of the series, 66.1% of genes; the
  intersection is the strictest, 69 of 71 candidates are genes (97.2%
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
  | + measured RNA, six lines | 289 | 60.1% | 58.8% | 88.6% | 65.1% |
  | predicted ∩ measured | 142 | 63.7% | **94.4%** | 79.9% | 62.0% |

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
- **Next.** 1. The `reader` construct in
  BioLang, so context gating comes from
  chromatin (the parser half belongs to the language owner). 2. Orthology
  from Ensembl Compara. 3. Measured Hi-C boundaries for the node edges: built
  2026-09-12 (`genomeos domains --chrom C --hic BIOSOURCE`, `genome/hic.py`,
  the biosource's deepest 4DN boundary file streamed once into
  `data/knowledge/hic` and compared with the CTCF-only boundaries under the
  same random control as the predicted map; tested offline) and blocked on a
  free 4DN account key that only Albert can create, since every 4DN download
  returns 403 without one, ENCODE has no Hi-C domain BEDs and the 3D Genome
  Browser's hg38 archive is gone; one command once `FOURDN_KEY` and
  `FOURDN_SECRET` are in `.env` (DATA.md "Optional: 4D Nucleome boundary
  calls"). The L1 ORF2
  and Alu sequence signatures are superseded by the RepeatMasker pass and
  stay as the fallback for chromosomes not yet distilled.
- **Owner.** genomeos-fe.

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
- **Owner.** genomeos-fe.

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
- **Missing.** Telomere length from a real 30x BAM (100 GB) not done; no
  polygenic scores; a missense variant is ranked by annotation, not predicted,
  and the trace still uses the canonical transcript rather than the one the
  tissue makes; the trio's 12,296 candidates are not yet reduced to the 60 to
  100 that a child really carries.
- **Next.** 1. Telomere from a streamed CRAM range instead of a BAM download.
  2. The missense effect predicted rather than ranked, on the transcript the
  tissue makes. 3. De novo candidates normalised across the three call sets
  (representation, then phasing by parent) until the count lands where the
  literature says.
- **Owner.** genomeos-fe.

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
  haematopoiesis; the worm's terminal fates below the founders are the
  observed lineage program, although every cell now carries its measured
  transcription factors (Ma 2021 atlas, `express`); the PAR polarity rules
  are not written (par-2 and par-3 are the misses against Digital
  Development); the external engines (libRoadRunner, MaBoSS, CompuCell3D)
  are deferred for lack of Python 3.14 wheels, although CI runs 3.12 and
  could test them as optional extras. Done 2026-09-11: haematopoiesis,
  the Body as a process-bigraph process gated by a network, space (sites,
  fields, division in place, contact inhibition, migration; the French flag
  grown from one cell), the measured reader, Packer 2019 at 84% on
  single-tissue ids (the rest is labelling depth), Digital Development as
  the published knockout set (7 of 11 modelled transformations).
- **Next.** 1. Terminal fates from measured factors: rules on the atlas
  (ELT-2, HLH-1, PHA-4 and the rest) with precedence over the lineage
  lookup, scored by the diff. 2. PAR polarity rules for the first divisions.
  3. The worm's founders in space with the real contacts replacing named
  senders. 4. libRoadRunner and MaBoSS adapters tested in CI on 3.12.
- **Owner.** genomeos-73.

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
- **Missing.** Two views promised in the UI track were never built:
  **Evidence explorer** (filter every rule by evidence kind and confidence,
  export the weak ones) and **Cell** (a cell type's active rule set and
  graph, run and watch). `pyproject` still says 0.1.0 and there are no
  release tags or PyPI package; the README "Project family" table still
  says BioLang and BioLib are "designed, not built"; no nightly real-data
  CI; no "known phenotypes it must reproduce" list per library from a
  biologist; `data/jobs` metadata and logs are tracked in git and change
  on every run; the jobs registry de-duplicated only in memory, so server
  restarts spawned five copies of the proteome job (fixed 2026-09-11 with
  a pid liveness check).
- **Next.** 1. Evidence explorer view over every loaded module.
  2. README refresh and a version bump to 0.9.0 with a git tag when the
  whole-genome milestone lands. 3. `.gitignore` for `data/jobs`, keeping
  the registry in code. 4. Nightly CI job that runs the real-data tests
  against cached reference data. 5. Cell view.
- **Owner.** genomeos-f7 (views, packaging, CI); genomeos-fe (jobs).

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
  in which cell types, silent in which). chr21: 747 entities (446 regions,
  155 elements, 83 genes, 63 domains), 155 rules, 20 unknowns, three `# test:`
  lines that `bio test` passes in 0.13 s; predicted evidence capped at 0.7,
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
  and variant mapping fails; the two axes disagreeing had a cause. The budget's
  "real unknown" therefore splits: copies first, attribution second, and the
  organiser reads `duplicated_fraction` at 0.5 or more as the copy flag. The
  classifier's own `similar_to` heuristic fires almost never (1 pair on chr21,
  0 on chr15); the curated pairs replace it.
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
  so sub-kb elements read their window.
- **Next.** 1. Every element in a gene's node scored per cell, not a sample:
  the closure needs a gene's whole regulatory input before it can judge it,
  which is the self-hosted model's first job (about 6,600 elements on chr21).
  2. The rest of the organising evidence
  on the compiled blocks: the copy flag and the human axis on regions (both
  measured now), repeat family, the reader's openness on the element itself
  rather than its node; the constrained_unknown tier re-read with copies set
  apart, since on chr21 three quarters of its blocks are duplicated. 3. The
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
  compact vertebrate, mouse for transfer.
- **Owner.** genomeos-f7 (this lane's files); the organism-level closure runs
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
  kilobase is unscored by Gnocchi. GRAMMAR-BY-COMPARISON.md §6.
- **Missing.** The human axis genome-wide and as a Blocks-tab lane; the
  origin per gene and library; curated duplication; the motif scan and
  operator patterns; phylogenetic profiling between libraries; the
  decompiled locus view; one developmental locus run end to end.
- **Next.** 1. The human axis over the other 22 chromosomes as a resumable
  job (a few MB each) and the case per block in the compiled BioLang
  `region` as a second evidence line. 2. Origin per gene through
  Ensembl homology at eight taxa; age distribution per library; the graph
  gains `orthologue_of`. 3. Paralogues and `genomicSuperDups`; `paralogue_of`
  edges; `similar_to` demoted to fallback. 4. JASPAR scan of promoters and
  constrained elements; `requires:` with evidence; recurring motif
  combinations per library. 5. Phylogenetic profiling between libraries.
  6. `genomeos decompile LOCUS`. 7. OCA2/HERC2 and the SHH ZRS across human,
  mouse, chicken and zebrafish.
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
| Telomere length from HG002 30x reads (stream a CRAM range, no download) | design of the range streaming | genomeos-fe |
| AlphaGenome live predictions for UNKNOWN regulatory blocks | `ALPHAGENOME_API_KEY` | Albert |
| Measured Hi-C boundaries against the CTCF-only node edges (`genomeos domains --hic GM12878`) | a free 4DN account: `FOURDN_KEY`, `FOURDN_SECRET` in `.env`; the code is built and tested | Albert |
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
  not a VCF.
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
