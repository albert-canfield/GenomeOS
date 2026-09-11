# GenomeOS: scope, architecture, areas, plan

Reviewed and reorganised on 2026-09-11 against the code, the tests, the
result files and the running jobs. This is the one document that says what
GenomeOS is, how it is built, what each area is for, what is done, what is
missing and what comes next. Evidence per finished task stays in
[PROGRESS.md](PROGRESS.md); the original phased plan is kept as history in
[ACTION-PLAN.md](ACTION-PLAN.md); lessons the code carries are in
[LESSONS.md](LESSONS.md).

Measured state at the review: 68 commits in two days, 42 `genomeos`
commands plus the `bio` toolchain, 16 web views, 220 tests passing and 5
skipped, lint clean, three `bio.std` modules, six organism programs, 25
chromosomes inventoried, four proteomes compiled, the whole worm grown from
one cell, a human body grown as populations, a therapeutic design dataset.

---

## 1. Main scope

GenomeOS is an executable model of biology for geneticists: a compiler that
turns a real genome plus public scientific knowledge into an executable
model, and a virtual machine that runs that model through time, from one
cell to an organism, with evidence and uncertainty on every fact.

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
| Engines import nothing from genome, knowledge or web | **broken twice**: `runtime/central_dogma.py` and `runtime/variant_effect.py` import the genome layer. Fix: move the codon tables and translation into the engine and have the genome layer call them, not the reverse |
| `lang` and `ir` depend only on themselves | `Locus` lives in `genomeos/coords.py`; move it under `ir` when the packages split |

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
- **Missing.** The `experiment` block (perturb, run, compare, report; in
  progress); a written grammar (BNF) and a versioned language spec; the
  two boundary leaks above; a `biolang` package with its own tests;
  imports from a registry; `bio.std` has no signalling, human development
  or timer modules beyond the worm.
- **Next.** 1. `experiment` block with knockouts checked against known
  C. elegans mutants (genomeos-73, in progress). 2. Fix the two engine
  imports. 3. Grammar spec in one file, generated from the parser's block
  table. 4. `bio` packaged as an extra entry point with its own test set.
  5. `bio.std.signalling` from CellPhoneDB and `bio.std.timers.human`
  from the Carnegie stages.
- **Owner.** genomeos-73 (parser, IR, body runtime, std); genomeos-fe
  (`bio.py`, central dogma).

### B. Genome decoding (reverse-engineering the sequence)

- **Goal.** Every base of a real genome accounted for: gene, transcript,
  exon, regulatory element, repeat, domain, or an UNKNOWN with a class and
  a confidence; the delimiters the cell reads learned from the sequence.
- **Code.** `genome/annotation.py`, `blocks.py`, `unknown.py`, `repeats.py`,
  `regulation.py`, `domains.py`, `signals.py`, `anatomy.py`, `fetch.py`,
  `lookup.py`.
- **Design.** GENOME-AS-CODE.md, GENOME-ANATOMY.md, SEQUENCE-GRAMMAR.md,
  UNKNOWN.md, NODES-READER-WRITER.md.
- **Data.** `anatomy_hg38_by_chromosome`, `unknown_chr*` (25),
  `rmsk_chr*` (18), `ccres_chr*` (24), `domains_chr*` (2),
  `signals_chr21`, `gencode_v50_chr*` rows (18).
- **Requirements.** Any chromosome on demand (`genomeos data fetch --chrom`);
  classification order fixed (curated repeats before ORFs); every class
  names its evidence layer; the block map shows it.
- **Missing.** Domains on 2 of 25 chromosomes; curated-repeat UNKNOWN pass
  on 14 of 24 (job running); the L1 ORF2 signature (long_orf is 23.6% of
  UNKNOWN and mostly LINE-1); silencers have no source; enhancer → gene is
  inferred from the CTCF domain, never measured (Hi-C or a predictive model);
  the sequence grammar has signals but no segment parser (Viterbi over the
  grammar) and no JASPAR promoter scan; the **reader** of
  NODES-READER-WRITER.md (which nodes are open in which cell type, from
  DNase/ATAC/methylation) has no code; no second mammal for node comparison.
- **Next.** 1. The fetch job runs the curated UNKNOWN pass and the domains
  for the chromosome it fetched, so a new chromosome arrives fully analysed;
  the genome-wide job finishes the rest. 2. Reader v1: ENCODE DNase per
  biosample distilled to "open nodes per cell type", shown as a lane in the
  block map (after the background jobs, same network budget). 3. Segment
  parser: Viterbi over the grammar with the learned PWMs, evaluated against
  GENCODE per chromosome. 4. Mouse chr19 as the comparison genome. The L1
  ORF2 and Alu sequence signatures are superseded by the RepeatMasker pass
  and stay as the fallback for chromosomes not yet distilled.
- **Owner.** genomeos-fe.

### C. Molecules (RNA, proteins, pathways, the knowledge graph)

- **Goal.** From gene to transcript to protein isoform to modified state,
  compiled from public sources into one definition per protein, verified
  against the curators, and executable as pathways.
- **Code.** `molecules/compiler.py`, `proteome.py`, `verify.py`, `rna.py`,
  `graph.py`, `reactome.py`, `uniprot.py`, `alphafold.py`; `flow/trace.py`.
- **Design.** PROTEIN.md, FLOW.md.
- **Data.** `proteome_chr*` (4 of 25: chrM, 21, 22, Y),
  `translation_vs_uniprot_chr*` (3), `graph_chr21`, compiled definitions in
  `data/knowledge/proteins` (18 MB).
- **Requirements.** UniProt accession is the id; Gene → Transcript →
  Protein isoform, never Gene → Protein; predicted never equals
  experimental; STRING association is not interaction; never call an API
  per simulation tick.
- **Missing.** 21 chromosomes of proteome (about 20 minutes per 400 genes,
  so roughly 17 hours of source time, running); translation verification on
  22 chromosomes; the knowledge graph exists for chr21 only; pathways run as
  reachability, not kinetics; no post-translational state model beyond the
  UniProt feature list; no isoform-level expression.
- **Next.** 1. Let the proteome job run to the end, one process (fixed
  today). 2. `genomeos verify` after each chromosome compiles, as part of
  the job. 3. Genome-wide graph once the proteome finishes. 4. Kinetic
  pathways where BioModels has an SBML version of the Reactome pathway.
- **Owner.** genomeos-fe.

### D. The individual (BioTwin)

- **Goal.** A person's genome plus measured state plus environment, forked,
  run and compared.
- **Code.** `twin/twin.py`, `twin/clocks.py`, `genome/telomere.py`,
  `genome/variants.py`, `runtime/variant_effect.py`, `calibrate`.
- **Design.** ACTION-PLAN.md Phase 3, STORAGE.md, LESSONS.md (timers).
- **Data.** `hg002_chr21`, `hg002_telomere_stream`,
  `hg002_methylation_stream`, `clock_GSE41169`, `clinvar_chr21_agreement`,
  `calibration_hematopoietic_stem_attrition`.
- **Requirements.** HG002 is the test human; measured state is labelled
  with its caveats (cell-line signature); every inferred parameter carries
  its confidence.
- **Missing.** The twin holds chr21 only (genome-wide streaming of the HG002
  VCF in progress); telomere length from a real 30x BAM (100 GB) not done;
  AlphaGenome live predictions need a key; no path for a user's own VCF
  beyond HG002; no pedigree or trio; no polygenic scores.
- **Next.** 1. Twin on every fetched chromosome (genomeos-fe, in progress).
  2. `genomeos twin build --vcf ANY` for a user's own genome with the same
  reports. 3. Telomere from a streamed CRAM range instead of a BAM download.
  4. AlphaGenome adapter live once the key exists (Albert).
- **Owner.** genomeos-fe (build, lookup); genomeos-f7 (user-genome path).

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
- **Missing.** The human program is counts, not mechanism (low confidence by
  construction); no runnable haematopoietic differentiation module (the
  0.3 proof in the old roadmap was never built; only the cell type and the
  library exist); Packer 2019 disagrees with the lineage at 26% of internal
  nodes (labelling depth); the Body runtime does not yet run inside a
  process-bigraph composite with a GRN or SBML model; the spatial engine has
  no division or movement; the external engines (libRoadRunner, MaBoSS,
  CompuCell3D) are deferred for lack of Python 3.14 wheels, although CI
  runs 3.12 and could test them as optional extras.
- **Next.** 1. `experiment` block and mutant phenotypes (genomeos-73).
  2. First human mechanism module: haematopoiesis from HSC to the eight
  lineages, decisions from cited transcription-factor rules, checked against
  the Sender & Milo counts and the known lineage choices. 3. Body runtime
  as a bigraph process beside a GRN. 4. Cells that divide and move in the
  spatial engine (the worm's 28 founders in space). 5. libRoadRunner and
  MaBoSS adapters tested in CI on 3.12.
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
- **Missing.** Copy number and structural variants from cBioPortal profiles
  (only patient-supplied CNV today); a `cancer.*` library layer; a
  peptide/HLA binding predictor (licence decision pending); no
  retrospective benchmark showing the pipeline recovers approved targets
  (ERBB2, CD19, BCMA, EGFR) from their tumours; the demo runs at data
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
- **Missing.** Search only; no experiment language (arrives with the
  `experiment` block); no published perturbation reproduced; no
  calibration of the predicted confidence against outcomes.
- **Next.** 1. Reuse the `experiment` block as the BioForge input. 2. Pick
  three published perturbations with quantitative outcomes (repressilator
  variants, a Boolean cell-cycle mutant, a C. elegans founder mutant) and
  make them tests. 3. Constraint language for design (must express, must
  not express, budget in genes and bases).
- **Owner.** genomeos-73 after the experiment block.

---

## 4. Data jobs to be done

Per-chromosome coverage on 2026-09-11 (25 human chromosomes including X, Y
and M). Each row is a job that can run unattended through the Progress tab
or the CLI; results land in `data/results` and are committed.

| Layer | Done | Job | Notes |
|---|---|---|---|
| Sequence, GENCODE rows, ENCODE elements, RepeatMasker (`fetch_<chrom>`) | 18 of 25 | `genomeos data fetch --chrom` | missing chr15–20, chrX; about 4 min each |
| Anatomy inventory | 25 of 25 | done | `anatomy_hg38_by_chromosome` |
| UNKNOWN classification | 25 of 25, curated repeats on 14 of 24 | `unknown_genome_wide` (running, one process) | redo of chr15–22, X, Y after their fetch |
| CTCF domains | 2 of 25 | `genomeos domains --chrom` | run after each fetch |
| Proteome compiled | 4 of 25 | `proteome_genome_wide` (running) | roughly 17 h of source time remaining |
| Translation verified against UniProt | 3 of 25 | `genomeos verify --chrom` | fold into the proteome job |
| Knowledge graph | 1 of 25 | `genomeos graph` | after the proteome |
| HG002 twin | 1 of 25 | `genomeos twin build` per chromosome | in progress (genomeos-fe) |
| Curated repeats distilled | 18 of 25 | part of fetch | |

One-off jobs, each needing a decision or a resource:

| Job | Needs | Owner |
|---|---|---|
| Telomere length from HG002 30x reads (stream a CRAM range, no download) | design of the range streaming | genomeos-fe |
| AlphaGenome live predictions for UNKNOWN regulatory blocks | `ALPHAGENOME_API_KEY` | Albert |
| Cohort expression for the other TCGA studies (`genomeos cancer expression --study`) | 30 s each, pick the studies | genomeos-f7 |
| Retrospective therapeutic benchmark (public tumours with approved targets) | choose the cases | genomeos-f7 |
| Packer 2019 disagreements resolved to one labelling depth | design | genomeos-73 |
| Mouse chromosome as the second mammal for node comparison | Ensembl GRCm39 fetch through the same code | genomeos-fe |
| Sender & Milo per-tissue turnover as maintenance timers in `bio.std` | none | genomeos-73 |

---

## 5. Milestones

| Version | Milestone | Proof |
|---|---|---|
| 0.1–0.8 | engine, compiler from public data, cell runtime, composition spike, spatial, whole worm, twin, debugger | done, see PROGRESS.md |
| **0.9 whole genome** | every chromosome fetched, classified with curated repeats, domains found, proteome compiled and verified, HG002 twin genome-wide, graph genome-wide | every row of §4 at 25 of 25; numbers in PROGRESS.md; README and version bumped |
| **1.0 experiments** | `experiment` block; C. elegans mutants reproduced; three published perturbations as tests; BioForge takes experiments as input; Evidence explorer | `bio test` passes the mutant programs; benchmark tests in CI |
| **1.1 human mechanism** | haematopoiesis as a mechanism module inside the human body program; reader v1 (open nodes per cell type) | lineage choices and counts reproduced with confidence above "low" |
| **1.2 therapeutics benchmark** | approved targets recovered from public tumours; CNA and SV; `cancer.*` libraries | benchmark test in CI |
| **2.0 BioLang standalone** | `biolang` package: lang, ir, runtime, std, `bio`; GenomeOS depends on it | a `.bio` program runs with GenomeOS uninstalled; two test suites |

---

## 6. Improvement and creativity pool

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

---

## 7. How the work is organised (several sessions, one checkout)

- All work on `dev`; Albert opens pull requests to `main`. No attribution
  trailers.
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
  milestone.
