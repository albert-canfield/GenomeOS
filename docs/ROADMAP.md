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
| **BioTwin** (individual) | 0.9 | any GRCh38 genome imported locally and used by carrier lookup, the gene report, the twin and the gene-by-gene walk; HG002 across all 22 autosomes with measured telomere and epigenetic age; fork, run, diff; predicted effect of a person's own regulatory variants | telomere from real reads; pedigree and trio; polygenic scores |
| **BioForge** (design) | search only | random-restart search under constraints with predicted labels; minimal-cell design estimate | experiments as input; published perturbations reproduced; constraint language |
| **Genome decoding** (GenomeOS) | 0.9 nearly there | every chromosome fetched and analysed (2026-09-11): inventory, UNKNOWN classified genome-wide with curated repeats (98.5%), domains on 24, curated repeats and ENCODE elements on 25, reader v1 (open nodes per cell type) on 2; signals learned | the segment parser (in progress); reader on every chromosome; a second mammal |
| **Molecules** (GenomeOS) | 0.9 reached | the whole human proteome compiled from seven public databases (19,478 coding genes, 25 chromosomes), translation verified against UniProt on 19,249 of them, genome-wide knowledge graph, RNA layer with GTEx, pathways as reachability | chrY verification; kinetics where BioModels has an SBML pathway; post-translational state beyond the UniProt feature list |
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
- **Missing.** A written grammar (BNF) and a versioned language spec; the
  two boundary leaks above; a `biolang` package with its own tests;
  imports from a registry; `bio.std` has no signalling, human development
  or timer modules beyond the worm.
- **Next.** 1. Fix the two engine imports (done for `experiment`: knockouts
  against the wild type reproduce the classic founder mutants, 2026-09-11).
  2. Grammar spec in one file, generated from the parser's block
  table. 3. `bio` packaged as an extra entry point with its own test set.
  4. `bio.std.signalling` from CellPhoneDB and `bio.std.timers.human`
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
- **Done (2026-09-11).** Enhancer to gene (feature b): deleting an element in
  its 1 Mb window and reading which gene moves turns "nearest coding gene in
  the CTCF domain, inferred 0.4" into a named gene with a tissue and a
  magnitude. On 200 sampled chr21 distal enhancers, 63.5% move some gene by at
  least 0.1 log2, and where a coding gene is named it is the nearest TSS in
  the node 67.8% of the time and inside the node 87.4%. Two findings worth
  more than the headline: a third of the time the nearest-gene heuristic names
  the wrong gene, and 43 of 127 elements the registry calls enhancer-like
  behave as silencers.
- **Next.** 1. Reader lane in the block map and a `reader` construct in
  BioLang so context gating comes from chromatin. 2. Segment parser
  (`genomeos segments --chrom C`): Viterbi over the grammar with the learned
  PWMs, evaluated against GENCODE per chromosome. 3. Fetch chr15–20 and chrX.
  4. Mouse chr19 as the comparison genome. The L1
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
- **Missing.** Pathways run as reachability, not kinetics; no
  post-translational state model beyond the UniProt feature list; no
  isoform-level expression. chrY verification and the triage of the
  translation disagreements by mechanism both landed on 2026-09-11.
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
- **Missing.** Telomere length from a real 30x BAM (100 GB) not done; no
  pedigree or trio; no polygenic scores.
- **Next.** 1. Telomere from a streamed CRAM range instead of a BAM download.
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
  haematopoiesis (first mechanism module, from the stem cell to the blood
  lineages gated by the established factors, landed 2026-09-11); Packer 2019 disagrees with the lineage at 26% of internal
  nodes (labelling depth); the Body runtime does not yet run inside a
  process-bigraph composite with a GRN or SBML model; the spatial engine has
  no division or movement; the external engines (libRoadRunner, MaBoSS,
  CompuCell3D) are deferred for lack of Python 3.14 wheels, although CI
  runs 3.12 and could test them as optional extras.
- **Next.** 1. First human mechanism module: haematopoiesis from HSC to the eight
  lineages, decisions from cited transcription-factor rules, checked against
  the Sender & Milo counts and the known lineage choices. 2. Body runtime
  as a bigraph process beside a GRN. 3. Cells that divide and move in the
  spatial engine (the worm's 28 founders in space). 4. libRoadRunner and
  MaBoSS adapters tested in CI on 3.12. (The `experiment` block and the six
  founder mutants landed on 2026-09-11.)
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
| Sequence, GENCODE rows, ENCODE elements, RepeatMasker, curated UNKNOWN pass, domains (`fetch_<chrom>`) | 25 of 25 (2026-09-11) | `genomeos data fetch --chrom C --analyse` | sequences 1.4 GB in `data/reference`, local only |
| Anatomy inventory | 25 of 25 | done | `anatomy_hg38_by_chromosome` |
| UNKNOWN classification | 25 of 25 with curated repeats; genome-wide 98.5% classified | done | `unknown_genome_wide` summary committed |
| CTCF domains | 24 of 25 (chrM has none) | done | |
| Reader (open nodes per cell type) | 25 of 25 | done | K562 and HepG2; `reader_genome_wide` |
| Proteome compiled | 25 of 25 (2026-09-11) | done | 19,478 coding genes: sequence 99.1%, domains 99.0%, function 86.3%, interactions 81.5%, pathways 58.2%, experimental structure 45.2%, predicted structure 98.2%, disease 25.5% |
| Translation verified against UniProt | 25 of 25 (2026-09-11) | done | 19,249 genes: 90.9% canonical identical, 97.8% exact for some isoform; the remaining disagreements are triaged by mechanism, including hg38 frameshift and nonsense alleles detected automatically |
| Knowledge graph | genome-wide (2026-09-11) | done | 41,982 nodes, 327,024 edges, 19,283 compiled proteins, largest component 15,165, 3,924 components |
| Enhancer targets, predicted (AlphaGenome) | chr21 sampled (200 of 6,618) | `enhancer_targets_<chrom>` | needs a key; about eight seconds per element, so it is a sampling job, not a sweep |
| HG002 twin | 22 of 22 autosomes | done | 4.05 M PASS variants applied, zero reference mismatches; chrX and chrY are not phased in the GIAB benchmark, so they are out of scope rather than pending |
| Curated repeats distilled | 25 of 25 | done | the 25 RepeatMasker BEDs (60 MB) stay local; summaries committed |

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

## 5. Next steps, consolidated

The ordered list across areas, each with the milestone it serves and the
session that holds it. Items 1–4 run in parallel today.

1. Whole-genome data jobs: **done on 2026-09-11**. Every row of §4 is
   complete: fetch and analysis, UNKNOWN with curated repeats, domains, the
   reader, the proteome, translation verification, the knowledge graph and the
   HG002 twin over all 22 autosomes. What remains of milestone 0.9 is the
   version bump and the release tag. genomeos-f7.
2. HG002 twin on every chromosome, then `twin build --vcf` for any genome.
   Milestone 0.9. genomeos-fe.
3. Haematopoiesis landed 2026-09-11 as the first human mechanism module;
   next mechanism modules follow the same pattern. Milestone 1.1. genomeos-73.
4. Evidence explorer view; the engine-boundary fix in `runtime/variant_effect.py`;
   `.gitignore` for `data/jobs`; nightly CI. Milestone 1.0. genomeos-f7.
5. Reader lane in the block map and a `reader` construct in BioLang.
   Milestone 1.1. genomeos-fe.
6. Segment parser scored against GENCODE. Milestone 1.1. genomeos-fe.
7. Three published perturbations as `experiment` programs in CI; BioForge
   takes experiments as input. Milestone 1.0. genomeos-73.
8. Body runtime inside a process-bigraph composite with a network model.
   Milestone 1.1. genomeos-73.
9. Retrospective therapeutic benchmark; cBioPortal CNA and SV; `cancer.*`
   libraries. Milestone 1.2. genomeos-f7.
10. BIOIR-v0.3 spec and the BioLang grammar in one file; `bio` with its own
    test suite. Milestone 2.0. genomeos-73. (The import boundary that blocked
    this is closed as of 2026-09-11: the Apache paths can be lifted into their
    own package without dragging the application behind them.)
11. README refresh, version 0.9.0 and the first git tag when item 1 lands.
    Milestone 0.9. genomeos-f7.
12. libRoadRunner and MaBoSS adapters tested in CI on Python 3.12.
    Milestone 1.1. unassigned.

## 6. Milestones

| Version | Milestone | Proof |
|---|---|---|
| 0.1–0.8 | engine, compiler from public data, cell runtime, composition spike, spatial, whole worm, twin, debugger | done, see PROGRESS.md |
| **0.9 whole genome** ✅ data | every chromosome fetched, classified with curated repeats, domains found, proteome compiled and verified, HG002 twin genome-wide, graph genome-wide | all reached on 2026-09-11; the proteome also ships as a packaged offline library. Outstanding for the release itself: version bump to 0.9.0, git tag, README refresh |
| **1.0 experiments** | `experiment` block; C. elegans mutants reproduced; three published perturbations as tests; BioForge takes experiments as input; Evidence explorer | `bio test` passes the mutant programs; benchmark tests in CI |
| **1.1 human mechanism** | haematopoiesis as a mechanism module inside the human body program; reader v1 (open nodes per cell type) | lineage choices and counts reproduced with confidence above "low" |
| **1.2 therapeutics benchmark** | approved targets recovered from public tumours; CNA and SV; `cancer.*` libraries | benchmark test in CI |
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
- **A self-hosted AlphaGenome.** The weights are published; self-deployment
  removes the per-request quota and turns enhancer-to-gene, the predicted
  reader and in-silico mutagenesis from per-chromosome jobs into genome-wide
  ones (docs/ALPHAGENOME.md).

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
