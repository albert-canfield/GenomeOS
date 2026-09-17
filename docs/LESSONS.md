# Lessons the core has learned from data

Every real-data run leaves a summary (docs/STORAGE.md). This page is the
distilled knowledge itself: what the engine now carries, where it came from,
and how sure we are. When a lesson changes a parameter or a rule in the code,
the code cites it. Numbers are from the 2026-09-10 runs; see data/results/.

## Genome and annotation

| Lesson | Evidence | Embedded as |
|---|---|---|
| Human mitochondrial genes start at AUU/AUA as well as AUG; MT-ND2 begins with AUU read as Met | 13 of 13 mtDNA proteins at published lengths only after adding initiator codons | `VERTEBRATE_MITOCHONDRIAL_START_CODONS`, `translate(initiator=True)` |
| ~0.1% of GENCODE "complete" coding models do not begin with ATG (NCAM2-201, DSCAM-203) | 2 of 2,705 chr21 transcripts | tolerance in the annotation test, not a rule |
| Ensembl and GENCODE GFF3 differ in feature types and id prefixes but carry the same information | C. elegans III: >2,000 genes, >97% of coding transcripts translate cleanly | one parser handling both dialects |
| The GIAB HG002 benchmark is unphased; 55,210 chr21 variants, 0 reference mismatches, net length change -1.8 kb / +0.6 kb per haplotype | `hg002_chr21.json` | unphased-het-to-hap2 policy, documented |
| ClinVar labels truncating frameshifts "nonsense" | agreement >90% for all four consequence classes once this is allowed | `compatible` mapping in the ClinVar test and distiller |

| The genome's block delimiters are statistical signals, not tokens: donor GT 99.1%, acceptor AG 99.8%, start ATG 220/221 with Kozak A at −3 in 55%, stops TGA>TAA>TAG, polyA signal in the last 40 nt 66%, CpG-island promoters 58% | `signals_chr21.json` | `genomeos signals`; the sequence-grammar design in docs/SEQUENCE-GRAMMAR.md |
| A learned donor matrix alone gives ~90% recall at ~7 false hits per kb: local signals must be combined by a grammar with context, the way gene finders and the spliceosome do | `signals_chr21.json` separability | scan reports relative scores; segments are never asserted |

| Genomes come in three organisations, measured: compact (mtDNA: 68.5% coding, no introns, 785 genes/Mb, one strand), dense (worm III: 25.8% coding, 95 bp introns, 192 genes/Mb) and sparse (human chr21: 0.8% coding, 50% intron, 2 kb introns, 5.5 genes/Mb, ~4 transcripts per gene). Exons per transcript (5–6) are conserved between worm and human; intron length is what differs | `anatomy_comparison.json` | `genomeos anatomy`, docs/GENOME-ANATOMY.md design budgets |

| Whole human genome: 20,107 coding genes, 1.21% of 3.088 Gb coding, 61.3% intronic, 27.9% intergenic; median 8 exons per transcript; gene-dense chromosomes (chr19 25/Mb) have the shortest introns (842 bp) and sparse ones the longest (2.7 kb): bases per gene is one budget split between introns and spacing | `anatomy_hg38_by_chromosome.json` | docs/GENOME-ANATOMY.md; Progress tab genome-wide table |

| The UNKNOWN 45% of chr21 is, by sequence alone: 30% assembly gaps, 18% unique intergenic, 15% long-ORF blocks (transposon ORFs and pseudogenes), 11% centromere, 6% satellite arrays, 5% mixed, 4% interspersed repeat, 4% promoter-like, 6% unclassified; 94% classified | `unknown_chr21.json` | `genomeos unknown`, Blocks tab labels |
| Sequence alone cannot see regulation: ENCODE's chromatin registry (13,356 elements on chr21, distilled to 150 KB) reclassifies 17.5% of the UNKNOWN space as regulatory with curated evidence; the CTCF motif marks node boundaries | `encode_ccres_chr21.json`, `unknown_chr21.json` | `regulatory` class, ENCODE blocks in the map, docs/NODES-READER-WRITER.md |
| Exact k-mer counting only sees young repeats; diverged Alu/L1 copies need mismatch-tolerant motif search (Alu core, ≤7 of 36 mismatches, 0 random hits per 200 kb) | first vs second run of `genomeos unknown` on chr21 | `alu_density`, window composition |

## Libraries (the genome's reusable modules)

| Lesson | Evidence | Embedded as |
|---|---|---|
| 45 libraries have data-computed membership: core.replication 235 genes, signalling toolkit >1,000, immune 1,531, ligand-receptor 2,417 | GO annotations + Reactome all-levels hierarchy | `genomeos/lib/data/members.json.gz` (114 KB), used when raw files are absent |
| Human GO annotation is thin for developmental master regulators and stem-cell markers (they are annotated as "regulation of transcription") | 20 of ~440 catalogue genes miss their library term; overall agreement 95.5% | catalogue keeps literature-curated genes; verification threshold 94% |
| Two GO ids I wrote from memory were obsolete (GO:0006306, GO:0014065) | caught by the ontology check | every `go_terms` entry is validated against the ontology in tests |
| Reactome's Ensembl mapping lists only leaf pathways | "Signaling by WNT" had no members until membership was propagated up the hierarchy | `Reactome.from_files` with the relations file |

## Timers and ageing

| Lesson | Evidence | Embedded as |
|---|---|---|
| Horvath 2013 predicts age in whole blood with r > 0.8 and MAE < 10 years from 353 CpGs; Hannum similar | GEO GSE41169, `clock_GSE41169.json` (per-sample table) | coefficient tables shipped in `genomeos/twin/data/`; intercept 0.6955 |
| Horvath's age transform is log below 20 years, linear above | round-trip test | `horvath_transform` / `_inverse_horvath` |
| Net telomere attrition of ~25-35 bp/year in blood and skin requires telomerase compensation in stem cells (colon crypt ~0.98, HSC ~0.6) or the model senesces by 50 | `genomeos age` calibration against cohort attrition | `telomerase_compensation` parameters with evidence |
| Senescence must be a hazard on the shortest telomere, not a wall on the mean | with a hard wall almost no cell senesces before 90 | exponential hazard with `telomere_hazard_scale_bp` |
| TelSeq-scale telomere length for HG002 (a man in his 40s) is 2,948 bp from 2 M streamed reads; TelSeq runs ~half of Southern-blot lengths | `hg002_telomere_stream.json` | `southern_equivalent_bp_inferred` = 2× (inferred, confidence 0.5) |

| Fitting the stem-cell division rate to HG002's streamed telomere at an assumed 45 years gives 3.8 divisions/yr (net 91 bp/yr on the Southern scale): higher than cohort attrition, as expected for a lymphoblastoid cell line sample | `calibration_hematopoietic_stem_attrition.json` | `genomeos calibrate`; twin runs use it with confidence 0.4 |
| GIAB HG002 DNA is a lymphoblastoid cell line: Horvath reads 97.5 y and Hannum 62.7 y from 409 clock CpGs at 36x, far above the donor's age, the known culture signature | `hg002_methylation_stream.json` | recorded in the result note; the twin stores it as measured state with that caveat |

## Dynamics and development

| Lesson | Evidence | Embedded as |
|---|---|---|
| The repressilator oscillates with Elowitz's dimensionless parameters (beta = 5, alpha ≈ 216) in both BioLang and the BioModels SBML | 5-8 peaks in 60 h; SBML assignment rules evaluate to the published values | `data/demo/repressilator.bio`, SBML engine test |
| The Fauré cell cycle has a G1-arrest fixed point without CycD and a 7-state cycle with it | Boolean engine | `data/models/mammalian_cell_cycle.bnet` |
| Segment count = frozen cells / (wavefront speed × period); cells freeze ~18 cells behind the FGF front | segmentation runtime | `SegmentationResult.expected` |
| A tristable switch read against a clamped morphogen gives endoderm → mesoderm → ectoderm; half its rules are inferred and the report says "low" | gastrulation runtime | `clamp` argument on the network runtime; uncertainty report |
| Sulston's early lineage follows from per-founder cycle times: 4 cells at 22 min, ~24 at 100 min, E slowest | lineage engine | `data/demo/celegans_lineage.bio` |

## The UNKNOWN space, genome-wide

- Streaming one chromosome at a time kept peak disk under 300 MB for a
  1 Gb job; the summaries total a few hundred KB. Stream, distil, discard
  scales.
- The largest class of intergenic DNA by chromatin evidence is regulatory
  (34%), not repeat: ENCODE's elements are dense enough that a 10 kb window
  with two of them is the common case, not the exception.
- An ORF finder run on repeat-rich DNA finds LINE-1 ORF2 everywhere; a class
  named "long ORF" is only honest once the repeat signatures are subtracted
  first. Order of classification matters more than the classifiers.

- gnomAD's Gnocchi is one Z score per kilobase, and the range reader counts
  every kilobase an interval touches: a 170 bp exon measures 1,000 track
  bases. Report shares of kilobases, never of bases, for anything shorter
  than the step (`aggregate()` in `attribution/variation.py` says
  `track_bases` for that reason).
- Within-human constraint at one kilobase agrees with the mammalian axis
  where both are strong and adds little on its own: elements constrained on
  both axes name a gene 75% and 72% of the time (chr21, chr15), above every
  other case; the constrained_unknown tier reads 1.7% and 3.2%
  human-constrained, no more than the neutral tier, and chr15's neutral tier
  reads 17%. Coding segments read 40% and 31%, introns 26% and 18%, so the
  axis is real and its block-level form is noisy; the case is `inferred` at
  0.3 to 0.6, never a verdict (`variation_chr21`, `variation_chr15`). Over
  the genome the ordering holds on every count: elements constrained on both
  axes name a gene 73.9% (825), mammals only 56.0%, people only 50.0%,
  neither 45.3% (`variation_genome_wide`). A track with no coverage (Gnocchi
  on chrY) is recorded as unmeasured on every block, never as zero: absence
  of a score is not evidence of tolerance.

- The most constrained intergenic blocks are often copies: UCSC's curated
  segmental duplications cover 4.7% of the constrained_unknown tier genome-wide
  but 61% of it on chr21 (15 of 20 blocks), 72% on chrY and 42% on chr22, and
  one block in five of the tier is mostly a copy. A duplicated block reads as
  constrained across mammals because the alignment lands on its paralogue and
  is unscored by gnomAD because variant mapping fails there, so the two axes
  disagreeing is a duplication flag first. The classifier's shared-20-mer
  `similar_to` fired once genome-wide (and was right): curated pairs replace
  it, the heuristic stays as fallback (`duplication_genome_wide`).

- A motif control must keep the dinucleotides: shuffling single bases
  destroys the CpG and GC runs promoters have, and long GC-rich zinc-finger
  matrices then read as the only enriched class; a dinucleotide-preserving
  shuffle (Altschul and Erickson) restores MEF2, FOX, PBX and PKNOX1 on real
  promoters. And a factor-pair statistic that assumes independence finds the
  GC content twice (MEF2A with MEF2D, CGGBP1 with GC-rich zinc fingers);
  operators need profiles clustered into families and a GC-matched
  expectation first (`motifs_genome_wide`).
- Library-level motif enrichment recovers textbook associations without
  being told them: REST 7.7x in `systems.nervous`, IRF2/3/7/9 in
  `systems.immune`, ZNF143 and THAP11 in `core.translation` and
  `core.replication`; the heart's GATA plus T-box pair is not in promoters,
  where the scan looked, but in enhancers.

- Ensembl's "vertebrate species" list omits the outgroups its gene trees
  use: yeast, fly and worm are in the vertebrate Compara dump (50,000 rows)
  but not in `info/species?division=EnsemblVertebrates`, so a ladder built
  from the species list stopped at Chordata for every gene. Place the species
  the dump names, not the species the list names (`_stream_placing`), and
  read the pan-taxonomic dump for the strata below the animals. Ensembl's
  taxonomy classification also omits Amniota, Tetrapoda, Theria and
  Boreoeutheria, so a duck reads as a fish unless Aves stands in for Amniota
  (`PROXIES`). At the deepest strata a Compara "orthologue" is a
  family-level call (HOXA1 eukaryote-wide through plant homeobox proteins):
  read the grade, not the exact stratum (`origin_genome_wide`).

- A motif-pair statistic has three confounds, and each one produced
  "operators" on its own: near-identical matrices (MEF2A with MEF2D), GC
  (CGGBP1 with ZNF93 in 27 libraries) and transposons (ZNF135 with ZNF460:
  promoters with both are 19.4% Alu against 1.1%). Count TFClass families
  (C2H2 zinc fingers individually), take expectations within GC by
  repeat-share strata, and the promoter operator table goes from 38
  libraries to none while single-family enrichments such as REST in the
  nervous system survive (`motifs_genome_wide`,
  `promoter_composition_genome_wide`).
- Across species, a list of "held" motif sites at one locus mostly restates
  that the sequence is conserved: at 88% identity nearly every window keeps
  some of 1,019 matrices, and held-site density does not separate VISTA
  enhancers from matched negatives. A site must be one genome position held
  by the same factor in every well-aligned species (the first call mixed a
  factor's positions across species), and grammar claims need negatives.
  Record positions in genome coordinates from the start; a position in a
  concatenated alignment depends on which species' blocks were joined.
- Correlated presence profiles are mostly shared age: random gene sets of
  the same origin make-up correlate at 0.55 to 0.73, which is where the
  "one developmental history" sat. Test library pairs against age-matched
  sets on exclusive members, and remove duplicated genes before calling a
  shared history (`profiling_genome_wide` null).

- Motif sites held across species do not separate enhancers from inactive
  conserved sequence: 39 limb and 33 neural VISTA enhancers against 76
  constraint-matched negatives hold the same density of factor-strict sites
  and no family more often after correcting for 275 families
  (`across_panel_vista`). Do not read a single conserved locus's motif list
  as its grammar (the ZRS looked like textbook HOX, PBX and MEIS logic), and
  change the readout rather than sweeping thresholds until a test passes.

- A motif-site call needs calibrating against measurement before it is used
  as evidence: with 1,019 JASPAR profiles at 85% of the matrix range, some
  site covers 99.9% of a regulatory element's bases, so "in a site" and
  "held across species" say nothing; at 95% sites cover 65% and enrich
  measured functional bases 1.11x. Use 0.95 for coverage questions and 0.85
  only for per-gene best hits (`satmut_grammar`, Kircher 2019).
- Measured against 9,834 saturation-mutagenesis bases, the motif score
  predicts a functional base slightly better than mammalian conservation
  (median AUC 0.58 against 0.54), both weakly, and the two together beat
  either alone: conserved bases inside a strict site are functional 36.5% of
  the time against 21.1% outside. Significance alone is depth-dependent (1%
  of FOXE1's bases, 77% of IRF4's), so report an effect-size definition
  beside it.
- Motif *arrangement* is not detectable over motif counts in a reporter assay
  (2026-09-16, `motif_grammar`, pre-registered, held out on chr8, chr9, chr21
  and chr22 of ENCODE4's lentiMPRA). Strict family-collapsed site counts lift
  Spearman with measured activity from 0.17-0.35 to 0.40-0.43; adding pairwise
  spacing bins, relative orientation and helical phase for 55 family pairs moves
  it by +0.002, -0.006 and +0.002, every interval across zero, and the control
  that permutes site labels and strands within the element does as well. Which
  factors have sites matters; where they sit relative to each other does not, at
  200 bp and off the chromosome. Fix the claim in code and commit it before the
  held-out fold is read: the honest answer to "is there a grammar" is worth as
  much as a positive one.
- Held across mammals and invariant among people does not sort genome units
  into syntax and value slots, even at base resolution (2026-09-14,
  `lexicon_axes_chr21`, `lexicon_axes_chr22`).
  - Per-base phyloP sees JASPAR 0.95 sites beyond their letters in every
    context: sites beat column-permuted decoys of the same letters in 549 of
    743 and 571 of 749 factors.
  - Per-base diversity in 89 HPRC assemblies sees at most a 2 to 3%
    depletion.
  - Among sites held across mammals, the share also depleted among people is
    the same for motifs and decoys (9.0% against 9.2%, 4.5% against 7.9%).
  - The panel resolves coding constraint (CDS at 0.46 to 0.49 of expected)
    and not motif constraint.
  - Every motif claim needs a decoy with the same letters, because AT-rich
    homeobox sites looked like syntax until their decoys did too.
  - A value slot cannot be told from a mutation-rate hotspot with diversity
    alone: the dELS class replicates as "more conserved than its flanks, more
    variable among people", and open chromatin's mutation rate reads the
    same.
  - The fossil tier failed a fourth time. Matched on replication timing, its
    copies are as variable as the same subfamily's intronic copies (0.994 and
    1.086), and the unmatched excess was timing.
- Compression prices knowledge in bits, and three rules make the price honest (2026-09-14,
  `compress_chr21`, `compress_chr22`, `compress_chr18`).
  - A layer's model must be inactive where its annotation says it does not apply. Left active,
    a JASPAR-site model saved 13 kb inside its 22 kb of sites and lost 87 kb elsewhere, because
    a model that says "uniform" still takes weight in a mixture.
  - Price an annotation against a baseline that holds the same sequence. Generic models primed
    with the duplication partners recover two thirds of the duplication layer's gain; the
    alignment alone is worth 0.22 to 0.33 bits per claimed base.
  - Base-level compression sees copies, not function. Unique sequence costs 1.89 to 1.94 bits
    per base in every tier, coding exons included, and constrained_unknown differs from neutral
    by 0.01 bits per base. A curated repeat label (31 bits) is worth about what the copies the
    models have already seen are worth: it loses to blind copy finding when subfamilies are
    fitted on one chromosome, and just pays (+3.8 kb per Mb) when they are fitted on two.

- **Three measurements of chromatin structure, three nulls, and one perturbation that works
  (2026-09-17).** Asked which element acts on which gene, every structural reading this project has
  made has failed, and they were made by different lanes on different data:
  1. **CTCF orientation** is real at measured boundaries (two thirds of sites point away in H1, K562
     and HepG2) and predicts nothing about our edges (convergent edges at their strand-shuffle
     median).
  2. **Measured boundary calls and their strength** do not distinguish a published rearrangement from
     a random cut of the same kind, in any of five cell types, called or scored.
  3. **Balanced 4DN Hi-C contact at 5 kb**, put in place of 1/distance in the CRISPRi scorer, *costs*
     0.081 AUPRC on 9,165 training pairs (interval clear of zero); its held-out +0.012 straddles zero.
     At CRISPRi distances raw contact is mostly the distance decay (0.379 against 0.442 for distance
     alone) and observed-over-expected alone collapses to 0.083.
  What does add signal at the same task is the **perturbation**: the AlphaGenome deletion lifts
  held-out K562 AUPRC 0.550 to 0.633 over activity-over-distance, and the node's containment excess is
  +2.88 points genome-wide, small but measured against random boundaries. The reading is not that
  structure does not exist — it is measured, and it is real where it was measured — but that **at these
  distances structure is not the discriminating variable and perturbation is**. A fourth structural
  reading needs a reason why it would differ from these three, stated before it is run.

- **A matched comparison that pools repeated controls invents its own sample size (2026-09-17).**
  Written the obvious way — for each target, append every control in its stratum to a list, then
  compare the two lists — a stratum with 50,000 controls is counted once per target that lands in it.
  Large strata dominate the estimate, and the standard error is taken on a list far longer than the
  data, so the p-value is manufactured: this produced "+0.284 at p 0.00014" where the standardised
  answer is +0.100 at p 0.18, and "+0.025 at p 0.014" where it is -0.031 at p 0.97. Both were
  committed before the defect was found. The correct form is direct standardisation: compute each
  stratum's rate once over the pool, average it over the targets that fall in those strata, and take
  the error on the number of matched targets. **The tell was speed, not statistics** — the pooling
  version is quadratic and ran for thirteen minutes on 961,227 elements before it was stopped, twice;
  the standardised version answers in 3.7 seconds. A matched estimate that gets slower as the arms
  grow is doing something other than standardising.

- **Standardising on covariates is not a substitute for conditioning on coverage — and a coverage
  stratum only bites a claim that had to be BOUGHT (2026-09-17).** Two halves, found a day apart on
  two panels. First: an apparent separation survived standardising on GC, distance and constraint
  (+0.41, p 0.026) and died only when "has anything been measured inside this window" entered as a
  stratum (−0.11, the controls scoring *more* often than the targets). None of the usual covariates
  says whether anybody ever spent a measurement on a window, so no amount of matching on them
  substitutes for asking. Second, the corollary: the same manoeuvre on a different claim moved it
  0.3529 → 0.3578, which is nothing. The mechanism is what generalises — **a direction had to be
  bought** (it needs a deletion spent on that window, so windows differ in whether anyone paid),
  **a value on constrained sequence is free** (read from a trio's variants and a public phyloP track,
  which cover every window whether or not anyone chose it).

  **You can tell which kind you have before running anything: count how many windows carry the
  input at all.** In the same result the bought input is present at 60 of 85 control windows and the
  two free inputs at 85 of 85. A claim whose input is missing nowhere cannot have a coverage
  artefact, and the stratum is a null manoeuvre on it; a claim whose input is missing anywhere is
  where the artefact lives. So report the difference both ways as a matter of course — that is what
  tells a reader which kind of claim they are being shown — and never quote the phrase "conditioning
  on coverage changed nothing" without saying that the input was universal, because on a free claim
  it was never going to.

- **The number was never wrong; the population it describes was never stated (2026-09-17).** Three
  instruments found the same defect on one day, and it is not an arithmetic error in any of them:

  - a project-wide `mean_confidence` of **0.617**, pooling predicted facts at 0.242 with experimental
    ones at 0.872 — a figure lying strictly between two populations half a point apart and describing
    neither;
  - one stated confidence tested against four measured populations, reading **+0.472, −0.291, −0.137
    and +0.518** against them. No constant fixes four disagreeing signs, because a level is a
    property of a population;
  - every per-block `core` call in the human panel, which is a depletion measured against a baseline
    that turned out to be **71% gene body** — no value changed when that was discovered, only the
    reading.

  Each number is computed correctly and each is unreadable on its own. That is what makes the class
  hard: there is nothing to catch, because nothing is broken. The defect is a missing sentence, not a
  wrong value, and a test asserting the value would pass.

  **So before a number is quoted, name the population it is over in the same breath, and if that
  cannot be done in one clause the number is not ready to be quoted.** This is why the answer all
  week was "report both and name each" rather than "pick the right one": where two populations exist,
  a single figure is not a summary of them but a fact about neither.

## The reference is one haplotype

- hg38 carries loss-of-function alleles at dozens of loci (olfactory
  receptors, CASP12, FCGR2C, IFNL4, CYP2D7, KIR2DS4, SIGLEC16): the
  reference reading frame is shifted or stopped where UniProt describes the
  working protein. A translation engine that disagrees with UniProt there is
  right, and the disagreement is a fact about the genome to report, not to
  fix.
- Make the signature specific before trusting it: "CDS length not divisible
  by three" holds for 5% of coding transcripts (5'-incomplete models); the
  same test restricted to canonical transcripts without an incomplete tag
  flags 29 genes, and the six mitochondrial genes among them that translate
  perfectly prove it measures the annotation, not the engine.
- Every batch of disagreements so far hid one real bug (a synonym match,
  selenocysteine, a stop-word symbol). Read them; do not round them away.

## Several sessions in one checkout

- The git index is shared by everyone working in the same checkout: a
  plain `git commit` ships whatever anyone has staged, and `git add -A`
  sweeps half-finished work into someone else's commit. Commit through a
  private index (`GIT_INDEX_FILE`, `read-tree HEAD`, add explicit paths,
  apply your own hunks of shared files, `write-tree`, `commit-tree`,
  `update-ref` with the old value) and the shared index stops mattering.
- Format and lint only the files you changed; a package-wide `ruff format`
  rewrites the other person's uncommitted lines.
- Say which files are yours before you start; the overlap (cli.py,
  server.py, index.html) is small when each feature is one subcommand, one
  route and one card.
- Committing through a private index solves the clobbering problem and
  creates a second one: `update-ref` does not touch `.git/index`, so the
  shared index falls behind by every commit made that way. After 125 commits
  it reported 200 phantom changes to an editor and carried 64 staged
  deletions of files that were on disk and committed. `git read-tree HEAD`
  after each commit, and the drift never starts.
- The shared index goes stale on its own: twice in one day it held staged
  deletions of every committed result file (248,000 lines), one plain
  commit away from leaving the repository. `git status` shows it as `D `
  rows; a mixed `git reset -q` fixes the index without touching the tree.
- Filtering shared-file hunks by keyword drops the hunk that does not
  contain the keyword: an `individual predict` subparser went out in a later
  fix-up because its hunk said "predict" and "regulatory", never
  "individual". After every push of a shared file, run the new command once
  from the pushed tree (the pre-push hook checks tests, not that the CLI
  surface is complete), or filter by line range instead of by word.
- **A liveness check that can match itself is not a liveness check
  (2026-09-15).** Three times in two days a process has looked for another
  process by name and found its own command line instead. A
  `ps | grep | kill` pipeline matched its own shell and killed itself
  mid-command, twice, the second time from a heredoc written to avoid the
  first; and a watcher polling
  `pgrep -f "enhancer_targets_all_chain.py"` matched a sibling watcher whose
  invocation contained that same string, so it reported the chain alive for
  hours and would have reported it alive forever. The false positive is the
  dangerous direction: the watch goes quiet and quiet looks like healthy.
  The pattern you search for must be one your own invocation cannot contain,
  and in this repository it need not be a pattern at all — every registry job
  writes `data/jobs/<name>.heartbeat`, and `jobs.last_activity(name)` against
  `jobs.STALL_AFTER` answers the question without naming a process. Ask the
  job whether it is alive; do not ask the process table whether something
  spelled like it exists.
- A job that stops is not always a job that died: the chain was paused
  deliberately at 03:09 so the AlphaGenome key could go to another lane, and
  from outside that is indistinguishable from a crash. `jobs.key_holder()`
  says who holds the key and what for, and it is the first thing to read
  before reporting a stall.
- **A stale base is a silent revert, and care is not the fix
  (2026-09-15).** Three times in one day a lane staged a whole shared file
  from a copy it had read before HEAD moved, publishing the older text over
  what landed in between. Twice it was a paragraph; the third took out the
  entire `share:` construct — parser clause, IR field, runtime, grammar row,
  documentation and 155 lines of tests — and nobody noticed until the lane
  that wrote it came back to its own feature. The last of the three was
  committed by a session that knew the rule, had a tool for it, and still
  lost the race: it read the base blob in one shell call and moved the ref in
  a later one. So the interval between reading and writing must not exist —
  read the base and write the index in one call, with the parent pinned at
  `read-tree` time — and the check must be mechanical, not remembered:
  `scripts/check_staged.py` asks the staged tree whether it deletes lines a
  recent commit added, names the peer commit when it does, and refuses. It
  costs a second and it replays the real accident correctly. A rule nobody
  can follow perfectly under concurrency belongs in a program.

- **A warning you read this morning does not survive being in a hurry; `set -o pipefail` does
  (2026-09-17).** `check_staged.py | tail -2 && git commit-tree ...` reports the refusal and exits 0,
  because a pipeline's status is its last command's, so the `&&` meant to stop the commit runs it.
  This is written in the tool's own docstring in bold, by the session that wrote the tool, after it
  had happened twice — and it happened a third time the same day, to that session, in a compound
  command that staged, checked and committed in one line to save a round trip. The commit was
  harmless and that was luck, not care. **If the guard must be piped, run the shell with
  `set -o pipefail`**, which turns the refusal back into exit 2; verified both ways. The general
  form is the one this file keeps arriving at from different directions: a rule that depends on
  remembering it under pressure is not a rule, it is a hope, and the fix is a mechanism — here one
  shell option, one word long.

- **A claim relayed between sessions loses its evidence and keeps its confidence, so the check
  belongs at the reader (2026-09-17).** One session wrote a summary of a finding; a second session
  half-drafted an edit to a shared document from that summary. The claim was wrong — it inverted a
  table three hundred lines below it in the same file — and it would have entered the record through
  the second session's hands with the first session's name nowhere near it. **A correct claim and a
  confidently phrased one are indistinguishable at the receiving end**, because everything that would
  distinguish them, the table and the file and the line, is exactly what the summary replaced.
  Writers checking their own work is necessary and cannot be sufficient: the failure mode is that the
  writer already believes it. So **before acting on a peer's claim about the record, open the file it
  describes** — and when sending one, relay the commit, path and line rather than the conclusion, so
  the reader can do that cheaply. Of the six instances of this class in one day, four were a claim
  checked against a summary of the record instead of the record; the only thing that ever caught one
  was opening the file.

  **The same mechanism hides a false PREMISE better than a false finding, because a premise is
  inherited rather than asserted.** One session proposed a length-matched control set, a second wrote
  it into a brief, a third registered a test around it — and the controls had been matched on element
  length exactly, zero tolerance, all 85 of them, since before any of it. `candidate_windows` says so
  in its first docstring line and the committed result agrees for every window. Nobody checked,
  because each of the three could reasonably assume an earlier one had; a premise arrives already
  believed, with no claimant to interrogate. Worse, one of the three had printed those very lengths
  on screen hours before while checking something else and had not read them: **evidence that answers
  a question nobody has asked yet passes unnoticed.** So when a piece of work exists to settle a
  question, read the construction of the thing it is about before running it — the cost here was two
  commands against a file committed for days, and it retired the whole exercise.

- **The staging guard earns its keep on your OWN lines, which is the case that feels like friction
  (2026-09-17).** A session correcting its own published sentence had `check_staged.py` refuse the
  commit, because the correction removed two lines the same session had added an hour earlier. That
  is the moment the guard looks like bureaucracy — the author knows why the line is going. It is also
  precisely the moment worth the friction: **a self-correction and a self-serving revert are
  indistinguishable from the inside**, and the flag costs one look and one `--force` that prints what
  it takes out. A guard that only stopped other people's mistakes would be a guard nobody needed.

  **And later the same day it caught what no test could have.** A lane re-ran a layer without
  `--write-programs`, so `programs_written` came back empty; the guard refused the commit because it
  removed 25 lines. Nothing about the result's findings was wrong, and **no test covered that field** —
  no test ever would have, because a test asserts what somebody thought to assert, and nobody had
  thought about a record of which programs a run wrote. **A diff guard asserts on removal itself**,
  which is the one property that does not require knowing in advance what matters. That is why it
  belongs beside the test suite rather than inside it: tests cover the fields with owners, the guard
  covers the fields without. Of three refusals that day, two were legitimate rewrites and one was a
  silent degradation — a ratio that argues for the flag costing a look rather than a confirmation.

- **A blind must be a boundary, not a list of places (2026-09-17).** A lane was set up to test a
  claim without seeing the previous set's answer, and was forbidden the new section and the new
  result files. The answer was also in the *old* section, put there deliberately an hour earlier so
  that a correction would sit beside the claim it corrected — `1 of 9` appears six times in that
  file, twice outside the forbidden section. A lane told to test a claim reads the claim, and the
  claim lives in the old section, so the blind forbade where the answer was written for the new
  section and permitted where it was written for honesty. **A list of places is a promise that a
  section will not be scrolled; a boundary is a file nobody may open, and anyone can check it with
  one grep.** The fix is never to unwrite the true sentence to protect a test — that trades a real
  result for a clean one — but to move the boundary and hand the lane a written brief carrying only
  what it is allowed to know.

  The general form is worth more than the instance: **two good record-keeping rules can pull opposite
  ways on the same paragraph, and the one already in the file wins silently.** "Put the correction
  next to the claim it corrects" was written an hour before "blind the lane that will test the claim"
  existed, so nothing announced the conflict. When a new rule arrives, the question to ask is not
  whether it is right but which existing rule it now contradicts.

## Engineering

- **Before grading an instrument's answer, check that it could have given the right one
  (2026-09-17).** The known-locus benchmark asked whether deleting a published enhancer names its
  published target, and recorded the ZRS as a miss: deleting it names LMBR1, the gene it sits
  inside, not SHH. That stood as the flagship long-range failure for weeks. The scorer resizes its
  input to 1 Mb around the element, so its reach is 524 kb each way, and **SHH is 979 kb away**. It
  was never in the model's input. Two more of the panel's targets are in the same position, SOX9 at
  1,450 kb and IRX5 at 1,164 kb. The check costs nothing — GENCODE gene bodies against arithmetic,
  no request — and it had never been written, because a reading that comes back looks like an answer
  whether or not the question reached the instrument. **A wrong answer and an unasked question
  arrive in the same shape**, and only the second one is free to detect. So for any instrument with
  a window, a detection limit, a vocabulary or a dynamic range, assert that the expected answer is
  *inside* it before scoring the result, and when it is not, say unaskable rather than pending, so
  nobody prices the requests that cannot buy it. **State that range in the units the instrument
  actually uses**, which is the part that is easy to get wrong: here it is gene-body overlap with
  the window, not distance to a promoter. On chr21, 79 of 5,174 predicted targets have a TSS beyond
  524 kb and every one is a long gene — mostly RUNX1, 1.22 Mb — whose body reaches in. A reach test
  written in the convenient units would have been wrong in the safe-looking direction, calling
  answerable questions unanswerable. The correction made the model look slightly worse
  here (0.882 to 0.867, because the out-of-reach loci were hits through other layers), which is the
  direction that should raise least suspicion.

- **A defect that preserves the quantity it is checked against is invisible to that check
  (2026-09-17).** Found twice in one day, in unrelated code, by two sessions.
  1. The bigWig summariser credited whole bins to an interval, so every absolute base count from a
     1 kb track was overstated — 4.62x at element level, where a cCRE is shorter than a bin. It
     survived because it **preserved every ratio**: numerator and denominator scaled together, so the
     means and fractions everything was concluded from were exact.
  2. The human body program stated germ-layer splits as fractions of the remainder, printed at four
     decimals, so the last split of each layer took 1.0000 of nothing and three populations held zero
     cells — one of them Glia, 24.5% of the ectoderm. It survived because it **preserved the layer
     totals**: the twenty-year count moved by 1.8e-6 when the three came back, and every assert the
     program carries passed before and after.
  In both cases the wrong number was the one nothing depended on, which is exactly why nothing caught
  it, and in both cases the fix was the same: **give the two quantities two names**. `bases` and
  `overlap_bases`; `fraction` and `share`. A field that answers two questions will be right for one of
  them and summed for the other. When a check and a defect are both conservative in the same
  quantity, the check cannot see the defect — so ask what a wrong version would *not* conserve, and
  assert that instead.

- **The fifth instance of that class, and it is a search: a grep that returns nothing looks
  identical to a thing that is not there (2026-09-17).** Asked what was left to build, I ran
  `ls genomeos/biovm/*methyl*`, got no matches, concluded "no module at all", told the user that the
  methylation state machine was the one unbuilt lane, and began building it. It has been running
  since 2026-09-16 as `genomeos/runtime/methylation.py` — U/H/M dyad states, replication, maintenance
  at fidelity below 1, de novo, TET erasure, every rate cited, with the solo-WCGW falsifier already
  run on eight methylomes. The directory I searched has never existed. **The empty result preserved
  every quantity anyone would have checked**: the command succeeded, the exit status was ordinary,
  the output was consistent with the conclusion, and nothing in it mentioned the path I had guessed.
  That is the same defect shape as the other four, wearing a search instead of an arithmetic.

  The general form: **a negative search result is a statement about the search, not about the world**,
  and the two are reported in the same shape — silence. So a search that informs a decision has to
  carry its own denominator. Say what space was searched and how big it was, not only what came back;
  when the answer is "none", show that the space was non-empty and correctly addressed. In a result
  file the rule is that **"measured and came back empty" and "never looked for" must never be the
  same output** — report the number assessed beside the number that had the property, and assert the
  assessed count is non-zero.

  **The same holds for a count, which is the commoner case:** a peer reported seven occurrences of a
  string from `grep -c 'A\|B'`, and the true count of A was six. `grep -c` counts *lines that match*,
  not occurrences, and the alternation folded a `B` line into the total; `grep -o A | wc -l` gives
  six. Nothing looked wrong — the command succeeded and returned a plausible integer. **A count is a
  statement about the query, not about the file**, so when a number will be quoted, check that the
  command counts the thing being claimed: occurrences or lines, one pattern or several, and whether
  two on a line collapse into one.

  **A figure typed beside a computed one will be wrong eventually, and the same change gave us the
  control arm (2026-09-17).** A lane reported a range as "0.014 to 0.871". The document in that very
  commit reads **0.0000 to 0.8714**, and is right, because that paragraph interpolates its bounds
  from the result file. The wrong number existed only in the commit message and a peer message —
  prose that was typed. *Same change, same minute, same session: the computed path was correct and
  the remembered path was wrong, and they diverged exactly where one of them stopped touching the
  data.*

  The typed figure was wrong twice over, which is the instructive part. The true minimum is **0.0** —
  the inert element, the row the whole change existed to surface — so quoting a non-zero floor
  silently asserts that zero is not a value, when here zero was the finding. And 0.014 was not the
  non-zero minimum either; that is **0.0064**. It was the smallest number visible in a twelve-row
  listing printed earlier: not a minimum over the data, nor over the non-zero data, but **over what
  happened to be on screen**. An extremum quoted from a screenful is a statement about the screen,
  which is the search lesson above wearing its third costume. So interpolate figures from the
  artefact wherever they are repeated — including in the commit message, which is three inches to the
  left of the paragraph that already does it.

  Memory is not a defence against this, and the timing is the proof: this arrived about twenty
  minutes after I had explained the same principle to another session and been thanked for it. Four
  of the five corrections this week needed an outside reading to catch; so did the fifth. The only
  thing that has worked all week is checking a statement against the artefact rather than against a
  recollection of it, **because the artefact does not round in our favour.**



- **A backtick in a double-quoted commit message is a command, and the word disappears
  (2026-09-17).** `git commit-tree -m "the method took ` + "`" + `compiled` + "`" + ` and the dispatcher never sent it"` runs
  `compiled` as a command; the shell prints "command not found" among the git output, where it reads as
  noise, and the message lands with a hole where the word was. One of 56 messages that day lost a word
  this way, and it was only caught because the stray error line was noticed. The repository's record is
  the thing being damaged, and it is damaged silently. Write commit bodies through a quoted heredoc
  (`<<'EOF'`), or a file, or without backticks; the same applies to `$(` and `!` in double quotes. The
  general form is the one this project keeps rediscovering: a shell will quietly reinterpret text you
  believed you were merely passing along. The same tool can be defeated two ways by how it is
  invoked: piping it keeps the reason and loses the exit code, redirecting it to /dev/null keeps the
  exit code and loses the reason, and this session did both to `scripts/check_staged.py` on
  consecutive days.



- SBML species given as `initialAmount` inside a compartment of size ≠ 1 are concentrations in the maths, and a kinetic law's value is substance per time: integrate rate/volume. Skipping that ran a 100-species model into zeros within a second and looked like instability; it was units.
- process-bigraph registers processes with `core.register_link` and applies float updates additively; the composite ran our two engines on one clock.
- libRoadRunner, MaBoSS, CompuCell3D and biolearn (torch) have no Python 3.14 wheels yet; the in-house engines keep their interfaces so adapters can replace them.
- Streaming beats downloading: the HG002 telomere estimate cost 12 seconds and 0 bytes, the methylation clocks 54 seconds and 0 bytes for 1.2 GB of bedMethyl; the same information from BAM/CRAM would have cost 100-200 GB of disk.
- A strand-specific RNA-seq track is labelled by the read, not the transcript: ENCODE's IMR-90 total RNA-seq (ENCSR424FAZ) reads antisense, so APP sits on the "plus strand" file, and a reader that takes an exon's signal from the file named after the gene's strand sees the cell as silent everywhere. The closure test caught it (every chr21 gene at zero in one cell of four); `MeasuredRna` and `attribution/closure.py` now probe a few dozen exons in both orientations and swap when the swapped one carries twice the signal. Never trust a strand label without a housekeeping gene to check it against.
- **Ensembl's human homology dump has no mouse (2026-09-12).** `homologies/homo_sapiens/Compara.116.protein_default.homologies.tsv.gz` lists 199 species, among them *Mus caroli*, *Mus spretus* and the rat, but not *Mus musculus*; the mouse–human pairs are in the mouse dump (`homologies/mus_musculus/...`, `homology_species == homo_sapiens`). A reader that filters on the species it expects and finds zero rows should say so loudly; the first Compara pass wrote an empty orthology file and the comparison silently fell back to MGI.
- **ENCODE's per-strand WGBS "methylation state at CpG" bigWigs carry every cytosine (2026-09-14).** On chr21 the K562 plus-strand file (ENCFF459XNY) has 346 positions in 2 kb, at CCCTT and AACTT as well as CG, most at 0%; a mean over them reads a CpG fraction diluted by unmethylated CHH. The per-CpG bedMethyl bigBed of the same experiment carries coverage and percent per strand call and is the source to read (`BigBed` in `genome/epigenome.py`). Check a file's positions against the reference sequence before averaging it.
- **One biosample term is several cells (2026-09-14).** The portal's "cardiac muscle cell" is H7-derived for two marks and RUES2-derived for three, and "CD14-positive monocyte" is a female ENCODE donor for all five marks and a Roadmap male for three: a rule that picks the deepest file per mark assembles a cell that does not exist. Pick per cell type the biosample that covers the most marks, then the file (`pick_marks`), and record the biosample on every field.
- **A model trained on a track is not an independent test of that track (2026-09-14).** Measured marks predict AlphaGenome's deletion direction better than the registry class does (AUC 0.62 against 0.53 on chr21 and chr22), but AlphaGenome learned from ENCODE's histone ChIP-seq in the same lines, so part of the agreement is the model reading back its inputs. Controls that keep the marginals (another line's marks, marks shuffled within class and DNase call) say the cell's own marks matter; only measured perturbations can say the element does.
- **A label that is right on average can carry no information about your calls (2026-09-14).** CTCF orientation is real at measured boundaries (two thirds of sites point away from them in H1, K562 and HepG2), yet splitting our CTCF-only edges by orientation predicts nothing (convergent edges at their strand-shuffle median). The edges were placed without the strand, so the strand cannot rescue them afterwards; used to place boundaries it works (reverse-to-forward flips beat random 1.4 to 1.6 times). Test a mechanism where it acts, not only as a filter on calls made without it. And check a boundary set's density first: GM12878's calls put a random position within 20 kb of one 54% of the time.
- **Consumers read absence from a list as the other answer (2026-09-14).** The reader saved `silent_genes[:200]`, and the Blocks lane, the gene report and the decompiler all call a gene "read" when it is not in that list, so every silent gene past the 200th on a large chromosome was shown as read. A list that stands for a set must be complete, or its consumers must ask the positive list.
- **A containment metric needs a placement control, not only a size one (2026-09-14).** "An enhancer acts inside its node" rises when nodes are fewer or bigger, so compare it with as many boundaries placed at random, and compare callers at the same boundary count. Shuffling node lengths along the chromosome is the wrong null: real boundaries sit in gene-dense regions, so shuffled ones cut fewer pairs and the null reads above the observed. On the full deletion archive the CTCF-only nodes beat random placement by 2.6 points, and the orientation-aware nodes, better on Hi-C, fall 14.5 points below it.
- **A calibration curve's shape can transfer while its level does not (2026-09-17).** The CRISPRi
  screens were turned into a measured confidence for the sweep's predicted targets. The shape held:
  leave-chromosome-out the curve is reliable in 9 of 10 bins at an expected calibration error of
  0.006. The pre-registered claim still failed on held-out pairs, in 4 of 10 bins, and all four
  failures were under-predictions in the same direction, because the held-out screens call 6.53% of
  their pairs regulated against 4.97% in the training screens. A single log-odds shift of +0.679,
  fitted on the held-out labels and therefore a description of the failure rather than a test,
  restores 9 of 10 bins. The level is a property of how a screen chose which pairs to test, not of
  the element: **a confidence band quoted without its population's base rate measures nothing**, and a
  calibration carried to a new assay needs that assay's own intercept before its numbers mean
  anything. The same reasoning names what the band cannot say genome-wide: every predicted target
  lands in the top bands because the calibration's population is regulated 77% of the time, so the
  discriminating axis is the predicted drop, and only 5.5% of coding targets clear the 0.2 where the
  screens called 105 of 105.
