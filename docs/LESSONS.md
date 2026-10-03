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
| Re-measured 2026-09-28: the order is right, the proportions are not. Mesoderm 0.10 against a stated 0.35, endoderm 0.433 against 0.20; the old test passed mesoderm only because 0.35 - 0.10 = 0.24999999999999997 < 0.25. A tolerance check must resolve float ties against the model, and an expectation with no source is not a bar | gastrulation runtime, `run_gastrulation(cells=60, hours=30)` | `tests/test_gastrulation.py` `_within`; mesoderm `xfail(strict=True)` |
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
  +2.88 points genome-wide, small but measured against random boundaries *(re-audited 2026-09-27: +2.90, CI +2.03 to +3.81, on a control matched on boundary count only; +1.2 to +6.6 across four baselines; and **+5.89, CI +3.18 to +8.46, on 661 measured CRISPRi pairs** — the direction survives measurement, the size was never identified)*. The reading is not that
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

- **Base overlap does not predict effect overlap, so a union of covariates is not bracketed by its
  own singles (2026-09-17).** Four covariates were removed from a background one at a time, and I
  twice described the combined answer as lying "somewhere between the largest single and the sum".
  Measured, it does not: on **three of five tiers the union falls outside that interval**, and in
  both directions — the sum overstates on two tiers and understates on two others.

  The tempting shortcut is the one the base counts invite. The four overlap by only **11.5%**
  (168.4 Mb double-counted, sum over union 1.1148), so "they barely overlap, therefore the sum is
  about right" looks safe and is wrong on four of the five tiers. **Removing a set of covariates does
  not subtract their effects; it changes the population the comparison is against**, so the union is
  a different comparison rather than a combination of the singles. Overlap in bases and overlap in
  effect are separate quantities and the first says nothing about the second.

  The extreme case is worth its own flag, and the lane gave it one (`union_overshoots_the_offset`):
  one tier "explained" **141%** of its offset, which is not more than fully explained — the ratio
  crosses 1, so with all four covariates gone the tier reads *below* the background the exclusion
  leaves. A share above 100% there is a statement that the comparison has moved to a different
  population, wearing the costume of progress.

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

  **A sixth instance, 2026-09-22, and this one is about the guard rather than the defect.**
  `crispri.logistic_fit` took the whole Newton step, and where classes nearly separate one full step
  saturates every linear predictor; a saturated row then offers a curvature of 9.4e-14 against a
  gradient of order one and the next solve sends the weights to **1e12**. Seven call sites inherited
  it. **What divergence preserves is the ordering**: on the eight-row design now pinned in the
  tests, the diverged fit gives seven of eight rows a probability of exactly 0 — three of them
  positives — and still scores **AUROC 1.0**, the same as the correct fit. The test guarding the
  function asserted the ordering and the sign of a weight, which are exactly the two properties a
  divergence leaves alone, so it passed on a broken solver. The repaired tests assert what a
  divergence cannot conserve: that the weights stay finite, and that the fit never returns a point
  whose objective is worse than the zero start it began from (the old solver: −4.4e12 against
  −5.55). **A test written from the quantity a result is read through will pass on the failure that
  leaves that quantity alone.**

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

- **A test that pins a defect names what would close it, and that name is a claim to be checked
  like any other (2026-09-21).** Half of the EML4-ALK defect was pinned as a test asserting
  today's wrong answer on purpose — ALK classed `direct_surface` at accessibility 1.0, both read
  off the full-length gene — with a docstring naming the three measurements that would close it:
  the fusion junction, the partner orientation, or transcript evidence for the retained domains,
  "none of the three in any table GenomeOS reads". **All three were in rows the client already
  fetched.** cBioPortal's structural-variant endpoint was asked for `projection=SUMMARY`, which
  drops them; the same row under `DETAILED` gives `site1HugoSymbol=EML4` and `site2HugoSymbol=ALK`
  — site1 is the 5' partner, so ALK is the 3' — both breakpoints, `breakpointType=PRECISE`, and an
  annotation reading "EML4 exons 1-20 with ALK exons 20-29" against an ectodomain encoded by exons
  1-19. The pin was right, the honesty was right, and the sentence that made it durable was wrong:
  a blocker recorded as "this project holds no such measurement" is a statement about what was
  looked for, and looking cost one request parameter. **Whenever a defect is left open on the
  grounds that the data does not exist, the cheapest next act is to check the source already in
  the code for the field, before the note hardens into a reason.** The lesson is not about fusions
  and not about cBioPortal: a `SUMMARY` projection, a default column set and a summary endpoint
  all drop fields silently, and none of them report what they withheld.

  **A second instance the next day, from a different direction, which is why this is a rule and
  not an anecdote (2026-09-22).** The direction question could only be asked of 44 element-gene
  pairs because "the sweep stores one target per element", and closing the other arm was costed at
  30 requests and authorised. It cost none. That sentence was true of the compact
  `all_elements/<chrom>.json` — a **derived table** — and false of the measurement: the same run
  wrote a per-element response cache at `threshold=0.0` holding every gene in the scorer's window
  with a signed per-cell value. The answerable set went from 44 pairs to 152, and the extra 108 are
  where the claim fell over. **A limit inherited from a summary of the data is not a limit of the
  data**, and the cheapest thing to do with any sentence of the form "we only have one X per Y" is
  to look for what the run wrote before it was summarised.

- **When nothing appears on both sides of a split, a correction fitted on one side cannot be
  checked for the other (2026-09-22).** The calibration's level error was blamed on prevalence, and
  the repair was to carry each screen's own base rate as an offset. The census written before the
  fit ended the question: **the screens it was fitted on and the screens it is read on are disjoint
  sets.** So the correction exists only after a screen has been run, and the targets that need it —
  everything in no screen at all — are exactly the ones it can never be quoted for. The same
  arithmetic shows how little the pooled number meant: one shift of +0.6793 restored 9 of 10 bins
  and is **a mean over per-screen shifts of −1.303 to +4.228**, with the wrong sign for the screen
  holding 63% of the pairs. **Before fitting a per-group correction, count the groups on each side
  of the split; if the intersection is empty, the thing being built is a description of the groups
  you have.** The honest output then is the part that transports — here the ordering, which moved
  by at most 0.022 AUPRC — and "no number" for the rest.

- **A reading that a new sample cannot surprise is a reading of the instrument (2026-09-22).** The
  epigenome layer reports `enhancers_active` per biosample, and across eleven biosamples it spans
  4.1× while `genes_read` spans 1.3×. The wide one correlates with DNase peak count at rho 0.8818,
  which on its own says only that they move together. The test that settles it costs one fetch:
  fit the reading on the biosamples you have, then bring in a lineage and a material the layer has
  never held — here testis and ovary, the first bulk tissues in a set of eleven cultured
  populations — and see whether the fit predicts them. **Testis landed at z = −0.11 and ovary at
  −1.29.** A number that a genuinely new sample cannot move away from its assay-depth prediction is
  measuring the assay. `genes_read` is the control that makes the point rather than a caveat
  against it: testis leaves the eleven's entire range at 15,078 of 20,094, on mid-range depth, in
  the direction its biology predicts. **Before a per-sample reading is used as evidence, fit it
  against the cheapest measure of how hard the sample was sequenced and report both.**

  *Two corrections from the lane sent to act on this entry. `enhancers_active` is produced by the
  **DNase reader**, `genomeos/genome/reader.py`, not by the epigenome layer as written above. And
  the family check it prompted found something worse than a confound in the same file:
  **`nodes_open` is a median split of each biosample against itself**, so it returns about half the
  nodes for every one of them — 9,988 to 10,016 of 20,002 across thirteen, a span of 1.00 where
  depth spans 6.81. A reading that cannot distinguish two samples even in principle is a stronger
  version of the same lesson, and it had been printed in the CLI and in two views as a comparison.*

- **A range built from two non-overlapping populations cannot be failed (2026-09-22).** The same
  lane registered that germline methylation would fall outside the somatic range, and it did not.
  The more useful finding is that it could not have: the "somatic range" is transformed lines at
  0.074–0.604 and non-transformed at 0.805–0.852 with **nothing in between**, so the interval it
  spans contains no value any sample could take. A test against an aggregate range is only a test
  if the aggregate is one population, and the check is to plot the members before quoting the span.

- **A curve is honest on the population it was fitted on and nowhere else, and the gate that chose
  that population is usually invisible (2026-09-22).** The confidence attached to every predicted
  enhancer target was fitted on the **245** training pairs where the measured gene happened to be
  the element's top predicted target — a population whose base rate is **76.7%** — and then used to
  band all **612,323** targets the sweep names. Applied to the held-out pairs its own gate excludes,
  it quotes **0.4126 against a measured 0.0603, ×6.85**, while landing at ×0.93 on the 40 it admits.
  Re-fitting it through the uncensored window barely moved that (×6.40): **the error is the
  population, not the features.** The reliability diagram never showed it, because a reliability
  diagram is drawn on the fitted population by construction. What shows it is one question, and it
  costs nothing: **quote the curve for the cases the gate threw away and compare with what was
  measured there.** Any confidence, band or calibrated probability in this project should carry the
  population it was fitted on beside it, and a consumer applying it off that population is quoting a
  number the curve was never asked about.

  **Two corrections to this entry as it was first written, both found by the lane sent to act on
  it, and both left visible because the entry is about not trusting a number past its population.**
  The figure **612,323** is the count of targets the sweep NAMES; **593,765** are banded, and the
  other 18,558 carry no band because their gene has no GENCODE v50 TSS. And the ×6.85 was not an
  error in the shipped table: `sweep_chromosome` bands the `predicted` and `predicted_coding` keys,
  and `crispri.deletion_values` sets `top_target` on a match to either, so **no swept target is off
  the gate**. The ×6.85 is what that curve quotes for a population that could not be asked about
  until `ElementResponses` existed. The lesson stands and is sharper stated correctly: **the curve
  was honest on every population it could reach, and the exposure appeared the moment a new
  population became reachable.** The mismatch that is real sits INSIDE the gate — 31.8% of the
  fitted 245 pairs are K562-named against 5.3% of the sweep, and 34.3% carry a drop above 0.2
  against 5.8%.

## Summaries are where claims grow

On 2026-09-28 four lanes tested milestone 1.3's clause 2 and reported carefully: the measured arm
said *cannot decide*; the reach control said its registered threshold was not met. The coordinator
compressed those reports into the plan and the progress log, and in compressing them made five claims
none of the lanes had made — that the failure was "about the sequence, not the window", that "19 more
blocks would decide it", that 823 blocks held "no measured element at all", that the neutral tier
"behaves identically", and, by the framing, that a difference in how often the model names a gene was a
verdict on the model's accuracy. A statistical review caught all five the same day.

Each is a different slip, and each is common:

- **An adjustment that fails to explain a difference does not explain it either.** "The tested
  covariates do not account for the gap" is what the data say; naming the cause is not.
- **A sample size planned from the model's own output, for an experiment measuring a different
  endpoint, is circular.** Both the effect and the variance came from the thing under test. And the
  "19" was not even that calculation's answer: it was the distance to a minimum-reporting floor.
- **A count needs its denominator and its rule.** 59 of 882 and 59 of 531 are different statements;
  59 at one overlap rule is 72 or 10 at another.
- **An interval that covers zero is absence of evidence of a difference, not evidence of none.**
  Equivalence needs a margin justified from outside the data.
- **Frequency is not accuracy.** How often an instrument says something is a property of the
  instrument until an outcome is measured.

**The rule adopted:** when a lane's result is carried into the plan or the log, its registered reading
travels with it word for word (*cannot decide*, *persists after adjustment*, *no difference detected*),
and any stronger sentence the summariser wants to write is a new claim that needs its own evidence. The
lanes' own discipline did not fail here; the step after it did.

## A verification tool can fail in the direction that flatters

On 2026-10-02 `scripts/manifest_rebuild.py`, the tool the night's compliance claims rested on, resolved
every declared input by `entry["path"]`. For an input recorded by `manifest.files_entry` — a group of
files read together, a person's per-chromosome calls, a census's tracked programs — that field holds a
**label**, not a path. So every grouped input was reported absent while its files sat on disk: the
rebuild opened and hashed not one of them, the entry never reached the `inputs` list the report shows,
and the run stopped with the dependency named unavailable, which review item R9 accepts as an outcome in
its own right. Five committed results were excused that way, among them `context_evidence.json`, where 12
of 13 declared inputs, standing for 873 files, were never opened. Two things made it quiet rather than
merely wrong: the `inputs` list carried no denominator, so the 12 unchecked entries left it and the one
that remained read as though all had matched; and a manifest declaring no inputs at all got a full
verdict having opened nothing.

The error ran one way. A result whose inputs were grouped could not fail the check, because the check
never reached its numbers — it was let off as resting on something unavailable. Nothing in the report
said what had not been looked at. It was found by a lane that hit it, worked around it locally by
recording each file under its own path, and said so rather than keeping the workaround to itself.

**What now makes it impossible.** `files_entry` names every file in a group, each with its own sha256 and
byte count, and marks the entry a group, which also tells a group apart from `input_entry` on a directory.
`manifest.group_digest` is the one implementation of the group digest, so the writer and the checker
cannot answer differently about which bytes it covers. The rebuild opens and hashes every member on its
own, and reports `inputs_declared`, `inputs_checked`, `inputs_unchecked` and `files_opened_and_hashed`,
with one entry per declared input, so an input nobody looked at cannot leave the list. An input the tool
cannot resolve — a group written before its members were named, an absolute path outside the linked
stores, a manifest with no inputs — stops the rebuild by name, and is never a pass. Ten tests pin it;
nine of the ten fail against the code as it stood, and the tenth asserts the group digest did not move,
so a digest recorded earlier still verifies.

**The shape to watch for**, which this project has met before: an instrument that preserves the quantity
its test watches. The question to ask of a verification tool is not whether it reports zero differences,
but whether it opened the files it says it checked, and whether it says how many it did not.

## A commit's message can describe one lane's work and carry another's (2026-10-02)

Four commits tonight carry one lane's files under a different lane's message, and a fifth
did the same on 2026-09-22 (`45e4f92`). History is **not** rewritten — in a checkout four lanes are committing to, a rebase to correct a
message risks far more than the wrong message costs. The mapping is recorded instead, so
a later reader is not misled by `git log`:

| commit | its message describes | its content actually is |
| --- | --- | --- |
| `ec8536d` | lane-re2g restamping its registration | **only** lane-rebuild's `docs/ATTRIBUTION.md` section |
| `dd5a223` | lane-re2g's identity-check loop fix | **only** lane-rebuild's `docs/ATTRIBUTION.md` section |
| `cc677e8` → relabelled `a3f1c04` | lane-mrfix's grouped-input fix | **only** lane-map2's five files (`response_map2.py`, `web/server.py`, `web/static/index.html`, `scripts/response_map_increment2.py`, `tests/test_response_map2.py`). **Corrected on the machine before it was pushed**: identical tree, same parent, `update-ref` with the old value pinned, so nothing was rewritten that anyone else had |
| `d0bc88a` | lane-rebuild's 2,000-draw benchmark rebuild — the message of `7fac349`, five hours older, word for word | **only** lane-198's `data/results/placement_cause_198.json` and its 155-line `docs/ATTRIBUTION.md` section. Already pushed, so not repairable the way `a3f1c04` was |

**How it happens, and the first explanation here was wrong.** The cause is the **shared
scratchpad**, not a shared file. Lanes write their commit message to a file and pass it to
`commit_own.sh -F`. When the guard refuses a call that both writes the heredoc *and*
commits, the heredoc never runs — so the message file is never written. If a **stale file
of the same name** is already there from another lane, the next call finds it and commits
that lane's text. `msg1.txt` is the name that collided. The content committed was correct
and complete each time, and each lane verified its own text byte for byte afterwards; only
the message was wrong, and no work was lost.

The same refusal-then-rewrite sequence hit this coordinator repeatedly tonight; it escaped
the collision only because its message files carried distinctive names.

**Read the direction of the error before naming the cause.** The lane that wrote `d0bc88a`'s
content reported it as a shared-index race between its run and a peer's. It is not: the
commit holds **exactly** that lane's own two paths and nothing else, which is what its own
`git add` of its own paths into its own private index produces. What came from elsewhere was
the **message**, and a message reaches `git commit-tree` from one place only, the file given
to `-F`. An index race would have produced the opposite signature — foreign *content* under
the lane's own message. Both faults lose a message, so both feel identical from inside the
lane; only the content tells them apart, and they have different fixes.

The lane then settled it in one command, and its own account of why it had not is the part
worth keeping: it reasoned from "my message is gone" to the failure mode it had read about,
instead of running `head -1` on the file it had just handed to `-F`. **Before naming a cause,
read the artefact you gave the tool.** What that one command showed: the file was dated
00:32, five hours and eleven minutes before the lane's own changes, and the directory each
session is told is "session-specific" in fact held some seventy `msg*.txt` files from every
lane of the session at once. A generic name in a shared directory was a collision waiting.

**What was built, 2026-10-02 06:00, and how much of the hole each part closes.**
`commit_own.sh` now refuses the message **file** on three grounds, each with its own exit
code and each overridable with `--force` only after the file has been opened:

| | refuses | exit | what it actually catches |
| --- | --- | --- | --- |
| 1 | a basename that is not `msg-<lane>-<purpose>-<epoch>.txt`, or that omits the lane given by `-L` / `GENOMEOS_LANE` | 67 | **the whole class.** The epoch is the working part: write the file in the same call that commits, and a lost heredoc makes the retry name a file that does not exist, so the script stops at "no such message file" instead of reading what was there |
| 2 | a message file more than 30 minutes older than the oldest uncommitted change among the paths being committed | 68 | only the hours-old case. `d0bc88a`'s message predated its files by five hours; but one of the five reused a file **15 minutes** old, well inside the slack, so this check would have passed it |
| 3 | a first line equal to the first line of any of the last 200 commits | 66 | the subset whose reused text was already committed — which is most of them, since the file being reused is usually a message that already landed |

The check that existed before this caught **none of the four after the first**: it compared
the whole message to `HEAD`'s only, and the reused text was older than `HEAD` every time. A
guard aimed one commit deep is a guard aimed at the case that does not happen.

Two things were deliberately **not** done. The slack in check 2 was left at 30 minutes with
its known hole written beside the number, because tightening it below about a quarter of an
hour starts refusing the legitimate order of work — write the message, then rerun the
generator that rewrites a result. And the stronger fix proposed earlier, refusing a staged
diff of a shared document with hunks outside the committing lane's declared section, is
still unbuilt; lanes stage by section (`scripts/stage_section.py`) by convention instead.

The three checks landed **while five lanes were mid-commit**, which the same proposal was
held back for a night earlier. The difference is that these three refuse on the message file
alone: a lane that has already staged nothing loses nothing by being told to rename a file,
whereas a refusal that reads the staged diff can only fire after the work is staged.

**The shape to watch.** A commit message is evidence about a commit, and like any other
record it can be wrong while every number inside it is right. `git log --stat` tells you
what a commit did; its subject tells you only what someone meant it to do.
## A nine-minute check is not a verdict on any one tree (2026-10-02)

`scripts/check.sh` on this project takes about nine minutes. In a checkout five
lanes commit to, that is long enough for the tree to change underneath it, and the
failures that result exist in **neither** the tree it started on nor the tree it
ended on.

The case. A coordinator run started 06:05 and reported **4 failed, 3,367 passed** at
06:14:33, every failure in `tests/test_promote_main.py`, each expecting the retired
`would run: git push ...:refs/heads/main` or the old "eligible for protected
promotion" wording. A lane committed `8df8799` at **06:14:39** — `promote_main.sh`
and `tests/test_promote_main.py` changed together, 94 insertions. pytest had
collected the **old** test code at 06:05 and executed it against the **new** script
on disk nine minutes later. Rerunning the committed pair: **45 passed**.

So the run tested a combination that was never committed: one lane's new script
against another snapshot's old assertions. Neither lane did anything wrong, and the
lane that landed it was right to finish rather than stop a step before its commit.

**What told the difference, and it was not the rerun.** `git diff` on the failing
file was **empty** — so the script was committed, not someone's working copy — and
`git log -- <file>` named the commit and its timestamp. A rerun alone would have
been "it passed the second time", which is not a diagnosis and is exactly how a real
intermittent failure gets waved through.

**The rule.** Before believing a failure in a shared checkout, read `git log` and
`git diff` for the file that failed. A failure whose file was committed by a peer
during your run is an artefact of your run.

**The mechanism** (queued 2026-10-02, supervisor's ruling): a check whose verdict
will be cited runs on the exact tree being committed, not on the live checkout.
Build the tree from the private index with `git write-tree`, materialise it in a
temporary worktree with the data stores linked as `scripts/pre-push.sh` already
does, run the check there, then commit that same tree. The verdict then names a tree
hash, and a peer's mid-run commit cannot reach it.

## A function that prints its own file by line number (2026-10-02)

`scripts/promote_main.sh` had `usage() { sed -n '8,10p' "$0" | sed 's/^# //' >&2; }`
— it printed its own usage by **line number**, so any edit above line 10 silently
changed what `--help` said. It had already rotted once: a second definition was added
below the first rather than the first being fixed, and bash keeping the last
definition is the only reason `--help` was right at all.

Removing six dead lines above it moved the usage text to 6-8. The repair is not the
new range; it is that a test now asserts the exact three lines `--help` prints, so the
range cannot rot silently again. **A test on a line range is not a test of what the
function says.**

The same mistake in a different place: the two clean-up lists in docs/ROADMAP.md gave
line numbers for every item, and after one intervening commit every number but one was
wrong — the refusal had moved from 251 to 238, a case from 93 to 54. Every item still
existed and still read as described, so the lists were sound and only their addresses
had rotted. **Line numbers in a written plan are labels, not addresses: match by
content, and quote the line's text beside its number.**

## A guard that has never been shown to fire is not a guard (2026-10-02)

Adopted after a lane planted its own counterexample unprompted. Before trusting a new
check, plant something that **must** trip it and something that **must not**, confirm
both, delete the probe, and say in the commit message that this was done.

The near-miss is the half that is usually skipped, and it is the half that catches a
check written too broadly. The lane testing "no shell script pushes to main" planted
the push line twice in one file, once as a comment and once executable, and showed the
check reported only the executable line. Without the commented copy it would not have
known whether it had written a test or a grep for a string.

For the three message-file refusals above, done in the live tree before the commit that
added them: a generic name refused 67, a first line copied from HEAD refused 66, a file
backdated to the previous day refused 68, and a correctly-formed name carrying another
lane's name refused 67 — with `git log` unchanged after all four, so nothing was
committed by the demonstration. The near-misses that must not trip: a message 20
minutes older than the work (inside the slack) and a correctly named fresh file, both
of which went through.

## A count is not a measurement of the thing you want to count (2026-10-02)

Twice in one night, in different hands.

**Once as a count of zeros read as a contribution.** 3,528 rows of the CRISPRi
benchmark carry a zero in the deletion column, and the coordinator concluded the pooled
gain was "largely not a test of the deletion feature". Measurement refuted it:
`activity+distance` alone scores **-0.0671** on K562 and **+0.0281** pooled,
`dnase+distance` alone **-0.1227**. Neither clears zero, so the deletion features carry
the advantage after all. A count of absent values says nothing about what the present
ones contribute.

**Once as a count of silence read as a count of error.** 198 of the 473
measured-perturbation loci attach to no compiled element, and that figure was written
into docs/ROADMAP.md and passed on as "what the model misses despite measurement". The
diagnosis of those loci reports what the measurements there actually said: **98
significant decreases, 8 increases and 593 well-powered nulls**, with only **66 of the
198** holding a decrease. So 132 of them are places the screen tested and found
nothing, where the model's silence agrees with the measurement. The located omission is
66 loci, not 198 — a third of what was claimed.

**The rule.** Before a count is quoted as evidence that something is wrong, get the
outcome breakdown at its denominator. "The model says nothing here" and "the model is
wrong here" are the same count only where the measurement said something at every
place. Both errors above survived review for hours because the arithmetic was correct;
what was wrong was the noun.


## Ten copies of one function, and nothing comparing them (2026-10-02)

`result_manifest.code_cleanliness` is what every published result says about which code
was uncommitted when its numbers were written: which files the writing lane had not
committed, which files other lanes had not committed, and whether any of those lay on
the counting path — the import closure the result's figures could have been read
through. Thirty blocks in `data/results` carry it.

It was copy-pasted into **ten** files, and the import closure behind it into **eight**.
Nothing compared the copies. Run against the same entry script on the same tree, they
gave **three different answers**:

| copies | paths | what they followed |
|---|---|---|
| six | 41 | only the `genomeos` package |
| `crispri_published_v2`, `crispri_benchmark_v2` | 68 | also `scripts.*` and each parent `__init__.py` |
| `cell2_eligibility` | 44 | the package, with a plain `Path.is_file()` |

Two of the three were wrong in ways the copies hid from each other.

**`data/results/cell2_eligibility.json` publishes three paths that name no file.**
`genomeos/genome/Genome.py`, `genomeos/genome/Annotation.py` and
`genomeos/knowledge/Reactome.py` are the *classes* those modules export; the files are
`genome.py`, `annotation.py` and `reactome.py`. `scripts/cell2_eligibility.py` was the
one copy that never took the cd263bc fix, so it asked `Path.is_file()` instead of
matching each path component against what its directory actually lists — and on this
Mac's case-insensitive APFS the answer was yes. That result's counting path is 45
entries here and 42 where the filesystem is case-sensitive, and a counting-path
difference is a **real** difference, not an environment field
(`tests/test_manifest_rebuild_environment.py`). Its recorded "0 differences" rebuild was
run on the case-insensitive machine.

**A file whose code does enter the numbers was reported as off the path.**
`scripts/placement_cause_198.py` reaches `scripts/response_map_coverage.py` through
`sys.path.insert(0, str(ROOT / "scripts"))` and reuses its closure;
`response_map_coverage.attaches_to` is what defines the 198 the result is about. Its
published counting path of 45 entries did not name that file. So the block's own note —
a foreign file outside the counting path cannot have entered the count — did not cover
it. **No copy caught this, including the broad one.** Resolving bare imports against the
directories the source literally puts on `sys.path`, as well as the repository root, is
what does.

**What was not wrong, and this bound is the point.** The 24 blocks written by the narrow
copies **under-state** their counting path, but their reported assurance still holds on
their recorded facts. Measured by an over-approximating test across all 30 published
blocks: **zero** has a `foreign_uncommitted_code` entry that falls inside the broad
closure but outside the narrow one. So no published
`foreign_uncommitted_code_on_the_counting_path` value is wrong, and nothing published
has to be withdrawn. The narrow copies also still reproduce their published lengths
exactly, so this was divergence between algorithms and never drift over time.

**The fix, and the rule.** One `counting_path` and one `code_cleanliness` in
`genomeos/manifest.py`, beside `code_revision` and `stamp`, both taking what a caller
may legitimately vary as arguments — the entry whose closure is the counting path, and
the paths the writing lane is answerable for — and guessing neither. The unified closure
is the broadest of the three, because a counting path is a claim about what *could* have
entered a number and the superset is the honest one. `save_result` refuses a registry
write whose manifest lacks the block, on the quarantine path that was already there; the
legacy allowlist is respected, so a name written before the contract still warns rather
than fails. No committed result was retro-edited: the two already written with the broad
closure reproduce their published counting paths exactly (67 and 68), and the narrower
entries pick the correct value up when their owning lanes next regenerate them.

**The rule.** A figure that several files each compute their own way is not one figure,
and the agreement nobody tests is the agreement that is not there. Where a number is
published as a property of the repository, one function computes it and a test fails if
a second definition appears — `tests/test_code_cleanliness_shared.py` fails if any file
under `scripts/` or `genomeos/` defines its own `counting_path` or `code_cleanliness`.
A re-export keeps a name that other writers already reach through; a `def` is a second
implementation, and a second implementation is drift waiting to happen.

## `pkill -f` and `pgrep -f` match every session on the machine (2026-10-02)

Twice, five hours apart, an unanchored pattern reached work that was not its own.

**First, a waiter that matched itself.** A loop waiting for a check to finish used
`pgrep -f check.sh`, which matched the waiting shell's own command line. Nine waiters
each saw "a check is running" — themselves — and waited on one another. A lane was
blocked about twenty minutes. Fixed by anchoring: `pgrep -f '^bash scripts/check.sh'`.

**Then a kill that reached a peer's push.** A lane ran, meaning to stop one stale
background check of its own:

```
pkill -f "scripts/check.sh" ; pkill -f "pytest -q"
```

Three `check.sh` runs were live in this checkout. One was inside another session's
pre-push worktree. Killing it made `pre-push.sh`'s subshell return non-zero, pre-push
removed the worktree immediately, and `bio test` — a child of the killed `check.sh`,
still running — lost its working directory. That push reported **red** with pytest
already **green at 3,485 passed**: the failure was an absent path under TMPDIR and then
a *relative* path to a tracked, present file, which is what a deleted working directory
looks like from inside.

**The rule.** Kill a background job by its own **PID**. Where a pattern is unavoidable,
anchor it to the start of the command line and make it specific enough that no other
session's process can match. `-f` matches the whole command line of every process on the
machine, including the one doing the matching.

**Two things this also teaches about reading a failure.** The push's own log was enough
to tell it was not a test failure: pytest had completed, and the error was in a later
leg, on a path that cannot go missing. And the lane that caused it volunteered the fact
unprompted, which is the only reason it was diagnosed at all rather than filed as an
intermittent `bio test` fault — the worst outcome, because an unexplained intermittent
failure teaches everyone to re-run and ignore.

**A mechanism is owed and not yet built.** `pre-push.sh` removes its worktree the
instant the check returns non-zero, without waiting for the check's children. That is
what turned one stray kill into a confusing failure rather than a clean one, and the
cleanup should wait. A guard in the shared-checkout hook refusing an unanchored
`pkill -f` is requested but not written: hooks are the owner's.


## The verdict moved twice without the mechanism moving once (2026-10-02)

`test_the_cleanliness_block_sets_the_uncommitted_files_against_the_counting_path`, in
`tests/test_crispri_benchmark_v2.py` and `tests/test_crispri_published_v2.py`, asserted
that `own_code_is_committed` was `True` and that
`foreign_uncommitted_code_on_the_counting_path` was empty. Those assert **the state of
the live checkout, not the behaviour of any code**. They were readiness checks for
writing a result at that moment, filed as unit tests.

In a checkout several sessions commit to, that makes the verdict a function of who else
is editing. Measured, from a run neither the fixing lane nor the coordinator produced: a
peer's status artefact shows 4,186 passed and exactly **two** failures, these two tests,
tripped by that lane's own uncommitted files. While the fixing lane worked, the field read
`['genomeos/attribution/compile.py', 'genomeos/attribution/direction_v2.py']` and
`own_code_is_committed` was `False`; its own 219-test run passed in that state; then the
peer committed and **the old assertion would have gone green with no change to any code it
tests**. The verdict moved twice and the mechanism moved once in neither direction.

**The cost was not the red itself.** The project's rule is "check.sh green on your own
files before any commit". With four lanes, at least one always has uncommitted code on the
counting path, so the rule is unachievable as literally stated, and every lane must decide
for itself whether a red suite is its own fault — every time. That is a standing invitation
to force past a red that *is* its fault. It also blocked a money-gating result for half an
hour while its lane correctly waited on a peer.

**The replacement, and why it is a tightening rather than a loss.** The two assertions
became invariants that hold in any tree state: every file in the foreign-on-path list is in
both the counting path and the foreign list; `own_uncommitted_code` is a subset of
`OWN_CODE`; `own_code_is_committed == (own_uncommitted_code == [])`. Then a **hermetic**
test in a temporary git repository covers the three states the live assertions only
pretended to check, plus a fourth planting an uncommitted foreign file *off* the path,
which must not be named — holding the field to being the filter it claims rather than a
copy of the dirty list. The signal given up was a duplicate: `scripts/check_staged.py`
refuses the commit of any staged result whose block has either value wrong, at publication,
which is the point that matters.

**Three things this teaches beyond the one test.**

A test that reads the live working tree is testing the tree, not the code. The giveaway is
that its verdict can change while nothing it covers changes — so when a failure appears and
disappears without an edit to the thing under test, suspect the assertion before the
mechanism.

**Check how many copies exist, and do not change the ones that are not the defect.** Four
copies of this test were found, not the two assumed. Two had the defect; two were already
invariant-style. Three *further* sites assert the same two field values and were
deliberately left alone, because each reads a **committed** `data/results/*.json` and
asserts what that published file recorded at write time — a fixed historical fact in git,
and the opposite of the defect. An over-eager reading of the ruling would have changed all
seven and silently weakened three.

**The replacement nearly passed for the wrong reason.** `manifest.is_code` returns true
only for paths under `CODE_ROOTS = ("genomeos/", "scripts/", "tests/")`, so a planted
`entry.py` at a temporary repository's root is not code to the revision stamp, and all
three hermetic cases would have passed while testing nothing. Both planted files had to sit
under `scripts/`. The counterfactual is what caught it: `inspect.getsource` of the real
function with its one path-filter line replaced, asserting exactly one occurrence so a
rewrite fails rather than passing vacuously, and shown **non-vacuous** by running the same
claim against an unstripped copy, which passes. A new test built to satisfy the
counterfactual rule was itself nearly worthless without one.

**And a premise the coordinator got wrong, corrected by the lane.** The brief said
`save_result`'s quarantine enforces these two values. It does not: `save_result` routes
through `manifest.cleanliness_problems`, whose docstring says it is "Checked by key set,
because that is what can be checked", and `genomeos/results.py` mentions neither value. The
value-level refusal is `scripts/check_staged.py` alone, deliberately — "Writing such a
result locally stays possible, deliberately: trial runs need it… The refusal belongs at the
commit that publishes the file, not at the write." The conclusion survived; the mechanism
named in support of it did not.

## When a classification is uncertain, err toward comparing (2026-10-02)

Two exemptions in `scripts/manifest_rebuild.py` were narrowed on the same day, and the
fix that was right for the first was **wrong** for the second. The difference is worth
keeping, because the obvious generalisation is the mistake.

**The first: one leaf, so an exact path.** The tracer's `opens` count could not be
reconciled across a rebuild. A suffix or name match would have been the natural rule, and
it is wrong: `data/results/response_map_increment3.json` asserts
`/per_element_response_cache/opens: 0` as a **computed** value, so a name rule would have
swallowed a real committed assertion and left nothing to fail on. Exempting the one exact
path leaves that one standing.

**The second: 536 paths, so a structural narrowing.** `is_timing` matched a key *name* at
any depth. On the committed `gene_row_locus.json` it set aside 4 leaves, of which only
`/seconds = 336.6` was a timing — the other three were the key-frequency counts 963,406 /
3,209 / 705 under `schema_keys["seconds"]`, spelled `seconds` because the key vocabulary
tallies key *names* found in the data. **This was a live silent pass, not a latent risk**:
altering 963,406 to 963,405 on the real 17,000-leaf result returned `[]`. A quantity the
run computed could be changed and the comparison said nothing.

An explicit list cannot fix that, and the number is the reason: **a timing name appears at
536 distinct paths across the ~1,100 results here, and every new result invents more.** So
the name test was kept and narrowed *structurally* — it stops at any key **the data
supplied**, detected by a `keys` or `vocabulary` token anywhere in the ancestry and
**sticky downward**, because a tally's values nest and testing only the immediate parent
would still swallow a count one level deeper.

**The rule, now standing for every verification tool: when a classification is uncertain,
err toward COMPARING.** A false-loud difference is read by someone; a false-silent pass is
read by no one. So the container test is deliberately generous: a container wrongly called
a tally gets a genuine timing compared, which is noise somebody resolves; a container
wrongly called ordinary sets a computed count aside in silence, which is the failure the
tool must not have. That asymmetry has its own test.

The same reasoning decides what **not** to narrow. `MUST_HOLD` matches a key at any depth
and was deliberately left alone: a tallied `own_code_is_committed` causes a spurious
**failure**, never a pass, and a rule whose error is loud is not the one to fear.
`is_resource` appears at only 3 distinct paths and *could* have been exact-path; it was
narrowed instead, because **the defect is in the matching rule, not in either list**.

**Two further things this taught.**

**The audit, not the fix, is the deliverable.** Two name-matched exemptions had been found
by accident in one day, which is weak evidence that a file was written with name matching
as a habit. Asking "how many others" returned a definite answer — two rules could end in a
silent pass, and there is no third — now pinned by
`test_every_exemption_is_either_an_exact_path_or_a_narrowed_name_test`. Two accidents are
a reason to count, not to fix twice.

**A report must print the set-aside leaf's VALUE, not only its path.** A path alone cannot
show that a leaf set aside as a reading is not one:
`/key_vocabulary/…/schema_keys/seconds` reads as a timing until its 963,406 sits beside it.

## Two instruments on one run, and the flattering one was the clean bill of health (2026-10-02)

`data/results/constrained_unknown_targets.json` is a result README relies on. Rebuilt
with the committed tool at `53f3b33`, it came back: **193 of 193 inputs declared and
checked, 0 unchecked, 193 files opened and hashed, 3,858 of 3,896 leaves compared, no
`must_hold` failures**, and the reading *"0 differences AFTER SETTING ASIDE 1 leaves (1
timing fields)"* — the one set-aside leaf a genuine `/seconds` = 236.7 wall-clock
reading. By every figure the rebuild produces, that result is reproducible.

An open-tracer watching the same run recorded **195 reads under `data/` against 193
declared**. The undeclared one is `data/results/unknown_chr21.json`, reached at
`scripts/constrained_unknown_targets.py:344-346` where a **self-check** calls
`unknown_scoring.unknown_blocks("chr21")`, which loads it at
`genomeos/attribution/unknown_scoring.py:89`.

**The two instruments disagree about the same run, and the one that reads as a clean bill
of health is the one that is wrong.** The file is git-**tracked** and 355,767 bytes, so a
clean worktree has it and the rebuild runs green — while its bytes are pinned by **no
declared sha256**. So a change to that file would change the result, and the rebuild would
still print "0 differences".

**The rule this establishes, measured on a live gate result rather than inferred from
code: a passing `manifest_rebuild` is not sufficient evidence for a published result,
because it is blind by construction to any file the manifest does not name.** "Rebuilt, 0
differences" means *"the declared inputs reproduce the file"*, not *"the result is
reproducible"*. The subset check has to run beside it, and a result's own report should not
be able to say the stronger thing.

**Three shapes of the hazard were found, and the one that was being described for months
was not among them.** The standing account was that `targets.run_elements` reads
`elements_where` out of a result file, so a pointer file names a table no grep can see. At
`190a1442`, `organise.inputs` **does** follow `elements_where` and declares the table it
points at. The real shapes:

1. **A declaration that silently shrinks.** `inputs()` appends the pointed-at table only
   `if w.exists()`, so a missing input is dropped from the declaration without a word —
   worse than never following the pointer, because the declaration still claims
   completeness.
2. **A code path no declaration function inspects** — here a self-check helper. No amount
   of care inside `inputs()` can see it, which is why the answer is a tracer at write time
   rather than a better declaration function.
3. **An undeclared input that is machine-local.** `therapeutic_benchmark` reads
   `data/knowledge/vep/vep_cache.jsonl`, git-ignored, declared nowhere in its 24-input
   manifest, under a README figure. Unlike (1) and (2) this one cannot even be rebuilt
   elsewhere: another machine would miss the cache and diverge.

**And the tracer that found all this had two defects of its own, both of which had already
produced a pass it had not earned** — the only direction of error that matters in an
instrument:

- it patched `builtins.open` alone, so `io.open` and therefore all of `pathlib` and pandas
  were invisible: **25 reads recorded where there were 3,259**, and one declared input
  reported "never read" when it had been read 3,258 times. `sys.addaudithook` on the
  `"open"` event sees them all, which is why it is the instrument and a list of entry
  points to patch is not;
- it dumped its record once per run at `atexit`, so under `--workers 2` **a pool child
  that read nothing overwrote the record**, and a result came back with **0 paths opened**
  — which satisfies "every read is declared" **vacuously**. Fixed by one append-only
  record written with `os.write`, and by withholding `verdict` **by name** unless reads
  were actually seen.

In the lane's own words: *"Had I not re-run after the first fix, I would have reported two
false passes."* An instrument's silence is not evidence; a verdict it cannot justify should
be withheld rather than defaulted to a pass. State the limit too — an audit hook cannot see
a C library opening a file directly, so htslib and pysam reads are invisible, and the claim
is about Python-level reads only.

## A stub may substitute a dependency's behaviour, never its existence (2026-10-02)

The AstroREG-2 sender passed **122 of its own tests**, a line-by-line code review by the
supervisor, and a dry run that exercised every refusal. It was launched against a live
paid API and failed **before sending a single request**:

    AttributeError: module 'genomeos.predict.enhancer_target' has no attribute
    'live_scorer_and_fetch'      scripts/astroreg2_send.py:234

`enhancer_target.Scorer` is a **type alias** — `Callable[[str, int, str, str], …]` — not a
class or a factory. The live client is `alphagenome_adapter.create_client` and
`AlphaGenomeAdapter`. The sender called a convenience helper that **had never been
written**.

**It failed safe, and that part was built right.** The exception landed one line after the
budget object was created and one line before any charge: `run_already_completed` false,
`charged_elements` 0, no ledger file, no completion event — so the authorisation was still
the same single run rather than a consumed one. Every guard on the way in had fired
correctly: all five of the approver's clauses checked, the plan digest matched the reviewed
list, cap 1,232 with 0 charged.

**A SECOND defect surfaced in the same twenty minutes, and it had been making the whole
tree red.** The push after the failed launch was refused by the pre-push check — 1 failed,
4,334 passed — on
`test_results_writers_guard.py::test_nothing_under_genomeos_or_scripts_writes_into_data_results_but_save_result`:
the sender's ledger was at `data/results/astroreg2_ledger.jsonl`, written directly rather
than through `save_result`. **The guard was right and the fix was not an exemption**: a
ledger is not a result — it is an append-only record of charges, written *before* each
charge precisely so it cannot be reconstructed from outcomes — so it moved to
`data/ledgers/`. Adding it to the allowlist would have papered over a file being in the
wrong place, and that guard's claim is worth keeping absolute.

**Why each defect was invisible, which is the lesson rather than either fix.**

- **Every test used a stub client and no network** — correctly, because testing it live
  spends money. So **every test substituted the entry point, and not one of them could
  prove it exists.**
- **The dry run passed because the scorer is constructed only after the `--send` check**,
  so the one line that only a real send could reach was the one line never reached.
- **The writers guard lives in another file**, so nobody running the sender's own tests
  would ever see it fail.

Two distinct blind spots: **a stubbed dependency, and a test in a file nobody in the loop
was running.**

**The rules adopted, free and permanent.** A stub may substitute a dependency's
**behaviour**, never its **existence**: every live entry point a paid path names gets an
existence-and-signature test — real import, `inspect.signature`, no network, no key — plus
a sweep of the path for attribute access on imported modules, asserting each resolves. And
a sign-off on a paid path requires **the full-suite verdict from the status file of the
committed tree**, not the module's own tests.

**The sign-off itself became blob-bound.** A review is of code, so
`SUPERVISOR_SIGNOFF` now records the **git blob sha** of the sender and of the module
holding the refusals, and `may_send` refuses if either differs — planted with a
one-character edit. The consequence is deliberate and will recur: the fix commit voids the
sign-off and requires a re-sign. **No approval can carry over to code the reviewer never
read.**

**The reviewer recorded its own miss as its own, which is why this entry exists.** In its
words: it **read** the line calling `live_scorer_and_fetch()` and **did not check that it
exists**, and it ran only the two module test files rather than the full suite, so the
writers-guard failure was invisible to it. "Dry run reviewed" was given on a review with
exactly the two blind spots above. A review that reads a call without resolving it is
reading prose, not code.

**And the chain held where it mattered.** No red tree reached the remote: the pre-push
check refused, did not retry, and said *"a red check is a red check: fix it and push
again."* Nothing was spent, nothing was published, and the authorisation survived intact.

## A declaration decoupled from the read cannot even be wrong in a way the filesystem would show (2026-10-02)

Three defects in input declarations were found in one day, each worse than the last, and
the third names the general case.

**First: a declaration that silently shrinks.** `organise.inputs()` appended a pointed-at
table only `if w.exists()`, so a missing input left the declaration shorter **with no
record**. The decisive argument came from the code rather than from principle:
`genomeos/attribution/targets.py:50` **already raises** on exactly that condition — *"A
summary whose table is missing raises, rather than reading as a run with no elements"* — so
**the declaration was laxer than the read it describes**, and a manifest was buildable for
a join that could not be read.

**Second: another file answering in the absent one's place.** `therapeutic_benchmark`
declared `data/knowledge/Ensembl2Reactome.txt`. That file is absent here, so
`genomeos/lib/membership.py:98` falls back to `Membership.distilled()`, which reads
`data/results/library_members.json` — 643,666 bytes, distilled **from** the very file that
is gone. **The manifest named the input that was discarded and omitted the input that was
read**, and the string `library_members` did not occur in it. A rebuild re-reading the
right file cannot tell. Found only because `save_result` **refused** the regenerated result
and named it.

**Third, and the general case: a declaration built without looking at anything.** The
`inputs` block written into every `organised_<chrom>` result was **four f-strings**. Being
decoupled from the read, it could not be falsified by the filesystem at all. Measured on
chr21:

- it named `budget_chr21`, while `read_axes` restated the tiers from **`budget_axes_chr21`**
  — and did so on **all 24 chromosomes**, because both files exist for all 24. The field
  recording which branch ran was discarded;
- it omitted **`enhancer_targets_all`** — **12,139 of the 12,439 elements** joined on
  chr21, and **the only run whose summary is a pointer at a ~1.4 GB local table.** The
  largest input, and the only unrecoverable one, was the one left out;
- it named `variation` and `duplication` **unconditionally** — the mirror defect, a
  declaration that will not shrink when it should.

**The rule that covers all three: a declaration must be produced BY the read, not written
alongside it.** A declaration that shrinks in silence is wrong and findable. A declaration
that names the wrong file is wrong and hard to find. **A declaration that never consulted
the filesystem is not even the kind of thing the filesystem can contradict.**

**The remedy, three-way and not a blanket raise.** An earlier ruling here — *a missing
declared input must raise with the path* — **was wrong applied verbatim**, and a peer proved
it on the very file it was handed: `Ensembl2Reactome.txt` is **genuinely optional in the
code with a documented fallback**, so raising would have made a README-backed benchmark
unrunnable. What holds:

1. **required** input missing → **raise, naming the path**;
2. **conditional** input missing → **declare WHICH BRANCH RAN**: the fallback *and* the file
   the fallback reads are declared, the stand-in is **hashed in the absent input's place**,
   and the swap is recorded (`inputs_optional_absent`);
3. **neither** present → **raise**, naming the packaged copy of last resort that would
   otherwise have answered **unpinned**.

**It is the house style, not three accidents.** Counted rather than estimated: **147
`exists()` calls inside input-declaring functions across 71 files; 91 of them actually
decide whether an entry enters an input list; 73 of those are over a fixed named path,
which is where required-versus-conditional lives; 16 enumerate a candidate set and are
legitimately tolerant; 7 already raise.** Of the 73, **exactly one has been classified
correctly.** And **16 must be left alone** — discovery over a candidate set, where a raise
would make a partial genome unwritable.

**A larger gap sits beside all of this:** `organise.run_and_save()` calls `save_result`
**with no manifest at all**, so `organised_<chrom>` carries *"manifest incomplete: missing
sources; missing inputs"* — and `inputs()` is never called by its own writer, only by three
downstream ones. **A declaration that does not exist is a separate and larger gap than one
that shrinks.**

## A red verdict that was about nothing in any tree (2026-10-02)

Several red verdicts that day were statements about no tree at all, and each read exactly
like a real failure until somebody opened the log.

Two were `git check-ignore` exiting 128 past a symlinked data store, so **20 tests ERRORED
at setup** — nothing was asserted wrongly; the tests could not start. **One of those
refused a push and held `origin/dev` thirty commits back for an hour.** One more was
`ruff` handed `data/results/manifest_headlines.json` as a lint argument: **96 "errors" in
a result file, the run dead at the `ruff-check` leg with `counts: null`, no test run at
all.** A fourth was seen and is **no longer in the record**, which is its own lesson: the
status store is **keyed by tree hash alone**, so a peer checking the same tree overwrites a
verdict, and **it cannot be used to count anything over time.** The coordinator reported
four and the record holds three; a sweep found exactly one of the ruff kind, and the second
was not written up as fact.

**The cause of the symlink pair was not what it looked like.** `.git/hooks/pre-push` was a
**COPY** taken at 04:34 that still symlinked the stores, while `scripts/pre-push.sh` had
moved at 17:52 to read-only clones — and `scripts/install-hooks.sh` **copies**, so the hook
was not the script and the correct fix never reached what the push actually ran.

**Two confident diagnoses were measured and refuted, and both were planted as negative
tests** — because a fix aimed at either would have passed its own test while every push
stayed red:

- **`$TMPDIR` behind `/var -> private/var` is innocent.** A worktree there with real store
  directories answers `check-ignore` 0, measured twice. Git resolves a linked worktree's
  root physically and a relative pathspec is taken from the process's physical cwd, so the
  leading symlink never enters the question.
- **`realpath` before `check-ignore` fails in the DANGEROUS direction.** An ignore rule is
  about a NAME, so resolving the path asks a different question: it can exit 1, **"not
  ignored"**, reporting a correct marker as **wrong** rather than unanswerable. `--no-index`
  changes neither reading; it was checked rather than assumed.

**The fix is to ask git about the name the marker concerns:** truncate to the first
component that IS a symlink and ask about that, since an ignore rule on a directory covers
its subtree. Trusted only when the symlink's own name is ignored; otherwise it refuses and
says it stops at the symlink rather than guessing past it.

**The lesson is that a verdict must say what KIND of thing made it red, or a tooling fault
and a broken codebase are the same artefact.** A test red and a tooling red had the same
`verdict` ("red"), the same `exit_code` (1) and a **byte-identical** `reason` ("the check
exited 1"). Every status file now carries `error_class`: `test` (an assertion failed),
`test_setup` (pytest red with errors and no failures — tests could not start, where
environment and tooling faults land), `tooling` (a non-test leg), `tree_moved`, `unknown`.
The symlink reds classify as `test_setup`; the ruff-on-JSON red as `tooling`. **It changes
nothing about acceptance — a merely-tooling red is still refused for a push, asserted by
test — and everything about how long it takes to know whose problem it is.**

**And the hook itself became a mechanism rather than a note**, since the stale copy is the
same defect class as a dirty verification tool: the installed hook is now a fixed thin
wrapper that runs `git show HEAD:scripts/pre-push.sh` and **refuses by name when the
working-tree copy differs from HEAD's**, because an uncommitted verification script gives
no verdict. **What it does not close, recorded in the wrapper and asserted by a test: a
commit can still weaken its own pre-push check**, because the hook runs the `pre-push.sh`
of the commit being pushed. That is the commit review's case — a visible diff in a tracked
file. **The wrapper closes the SILENT case.**

## The handler could not be reached, and its plant proved it against a patched body (2026-10-02)

A guard was fixed, planted, reviewed and committed, and **the branch it added could not fire.** It
had been reasoned about rather than measured, by the person who asked for it.

`genomeos/results.py: _is_registry` decided whether a write was going into the result registry, and
the flag it returned gated **four** refusals. It read `results_dir.resolve() == RESULTS_DIR.resolve()`
inside a `try`, with `except OSError: return False` — failing OPEN, so an exception turned all four
refusals off at once. That was found, and the fix was `except OSError: return True`: a question that
cannot be answered is treated as the registry, so the write refuses. It was planted both ways and it
went in at `8560b82`.

**Measured afterwards, on CPython 3.12 and APFS:**

| condition | `Path.resolve()` non-strict | `os.path.samefile` |
| --- | --- | --- |
| symlink loop | **`RuntimeError`** | `OSError` errno 62 |
| unsearchable parent | **returns, no raise** | `OSError` errno 13 |
| 300-character name | **returns, no raise** | `OSError` errno 63 |

`issubclass(RuntimeError, OSError)` is **False**. So on a symlink loop the writer *raised* instead of
refusing, and on the other two `resolve` handed the path back unchanged and the spelling comparison
answered `False` — the fail-open it was supposed to have closed. **The `except OSError` branch was
unreachable through the filesystem.** It fired in exactly one place: a test that patched the module's
source to raise.

**That is the part worth keeping.** The counterfactual was real and it passed, so nothing looked
wrong. But it was planted against a **patched body**, and what it proved was that the handler's logic
was correct *if reached* — never that anything could reach it. A plant that constructs its own trigger
tests the branch; it does not test the path to the branch. The question "can this condition actually
arise from the thing the code touches?" was never asked, and it is a different question from "does
the handler do the right thing".

The fix is identity by file — `os.path.samefile`, which asks the filesystem which file a path names —
and it closes the hole and makes the handler reachable **for the first time**, because `samefile`
raises `OSError` for all three conditions where `resolve` raised the wrong class or nothing at all.
A reader seeing `except OSError` merely moved would otherwise assume it had been widened.

**Whose error it was.** The coordinator's, not the lane's: it pressed for `OSError -> True` after
reasoning about the handler, and reviewed and approved a plant that patched the source rather than
asking whether the filesystem could produce the exception. The lane that implemented the real fix
measured all three conditions, found the gap, and reported it as a third hole nobody had briefed.

**The rule: a plant that supplies its own trigger proves the branch, not the reachability.** When a
handler exists for a condition the environment is supposed to produce, produce it — a real symlink
loop, a real unreadable directory, a real over-long name — and if that cannot be done, say in the
test that the branch is unreachable as written rather than letting a patched body stand in for the
world.

## The answer was on the next line, and a committed registration had already forbidden the method (2026-10-03)

Two lanes were briefed that night on questions **already measured**, by the two sessions whose job
was to stop exactly that. Both briefs cited a note and neither read the paragraph after it.

**The mechanical cause is worth keeping, because it is not carelessness and a better pattern does not
fix it.** `docs/ROADMAP.md` records the open question and its answer in consecutive paragraphs, and
the text wraps:

```
  which of the two it is has never been measured.
  That question is now measured, classes and thresholds committed before a single
```

Any single-line `grep` for the question's own phrasing lands on **"never been measured"** and the
answer sits on the following line, outside the match. Both of us searched for the phrasing we
remembered and both of us got the stale sentence. The lane that was sent to do the work found the
answer on its second tool call, by grepping **all of `docs/`** rather than one file and by using
`alias|synonym` rather than the question's wording — and then by **reading the section**.

**The guard is not a cleverer regex. It is `grep -A5`, or reading the section.** A pattern tuned to
the sentence you remember will keep finding the sentence you remember.

**The second half is worse than the duplication.** The coordinator's brief defined a tolerance tier —
"within 100 kb counts as the right place under another name" — and asked for it as a reported rate.
The governing registration, committed before a single miss was classified, forbids exactly that, with
its reason written down:

> The distance classes are DESCRIPTION ONLY and feed no rate, ever. The reason is not taste: a rule
> that admitted `within_10_kb` would score the nearest-gene baseline as a hit, and that baseline is
> the control the whole fourth frame was drawn to beat. A benchmark whose hit rule tolerates distance
> cannot then report that the model follows proximity.

Had the lane complied, it would have broken a committed prohibition and produced a rate that flattered
the model by scoring its own control as a hit. **It is the same act as using a track share as a
control after the registration forbade it in writing** — eighth instance in one session of reasoning
from a quantity without reading what governed it.

**And the measurement, which the briefs were asking to redo:** of 98 strict misses, **15 overlap the
published target's body, 11 more within 10 kb, 50 within 100 kb, 16 elsewhere** — so 76 of 98 name a
gene inside 100 kb. The headline is neither of the two things the question posed: **the miss is
overwhelmingly a NEIGHBOUR, not an alias.** And the control decides what to do about it: the
tolerance moves the model 10 → 16 and moves the nearest-TSS-in-node baseline 10 → 15, so **the
comparison the benchmark reports does not move.** Tolerating aliases lifts the proximity rule as much
as it lifts the model, which is why no second rate was adopted.

**One smaller thing the lane did right and should be the norm.** It ran no acceptance suite and said
so rather than reporting a green: it had changed no file, so there was nothing to accept, and
"running the suite to decorate a close-out would be a green with no change behind it." A green
attached to no change is not evidence; it is furniture.

## A run hashed its own previous output as a declared input, then replaced it (2026-10-03)

A census script's first write globbed its own output directory, hashed **the previous run's copy of
its own result** as a declared input, and then overwrote that file. The manifest it produced named a
path with a sha256 that the bytes at that path no longer matched — **an input pin that nobody,
including the writer, could ever verify.**

It is worth recording because nothing about it looks wrong while it is happening. The declaration is
honest: those bytes really were read. The hash is correct for what was read. The only defect is the
ORDER — the file is read, declared, and then replaced by the same process, so the declaration
describes a state that the run itself destroyed. A reader checking the pin afterwards finds a
mismatch and cannot tell it from a tampered input.

**The fix is an exclusion with its reason on the artefact's own face:** the script excludes its own
output name while writing it, states the exclusion in the result's `exclusions` and in an `itself`
field, and a test asserts both halves. The bad file was deleted and regenerated from the committed
script rather than patched.

**The general form: a writer may not declare an input it is about to replace.** Any process whose
output lives in the directory it reads has this hazard, and the project has several — censuses,
organisers, anything that globs `data/results`. The cheap check is to ask, of every declared input,
whether this run writes to that path.

The lane found this in its own first write and fixed it before the result was committed, which is
the only reason there is nothing to correct.

## A tool that works on one interpreter or OS is not a working tool (2026-10-03)

Two faults in one day, unrelated in mechanism and identical in shape: each passed on every developer
Mac and made a script **unrunnable in CI**, so the thing that broke was not the logic but the
assumption that one machine's behaviour is the behaviour.

**`mktemp -t <prefix>` is not portable.** GNU mktemp, which is what CI runs, REFUSES a `-t` template
with fewer than three trailing X's: the run of 2026-10-02 12:08 UTC printed `mktemp: too few X's in
template 'genomeos-index'` and stopped. BSD mktemp accepts a bare prefix and invents the suffix
itself, so the same line had worked locally for as long as it had existed. `scripts/push_own.sh` could
not run on Linux at all, which is the only platform that matters for a push gate. The fix is a full
path with six X's and a comment saying which mktemp refuses what.

**argparse's positional matching is version-dependent across CPython MICRO versions.** Two optional
positionals with an option between them are not parsed the same way by 3.12.x as by 3.12.y, and that
produced two full-suite reds that existed **only in CI** — nothing to reproduce locally, no file to
bisect, and the temptation to blame the runner. The repair is `class StableParser(argparse.ArgumentParser)`
in `genomeos/cli.py`, which carries CPython's own guard so the parse does not depend on which micro
version happens to be installed.

**The lesson is the class, not the two instances.** A shell builtin's flags, a standard library's
argument matching, a filesystem's case sensitivity, a locale's sort order: each is a place where "it
works" is a statement about one interpreter and one OS. Where a script is run by CI, the question is
not whether it runs but **where it has been shown to run**, and a green that has only ever been
produced on a Mac is evidence about a Mac.

## A red run whose first job dies early reports almost nothing (2026-10-03)

A failure summary was read as the whole suite's verdict when it came from **one job of several**. The
line quoted was "5 errors"; the real line was `25 failed, 3844 passed, 129 skipped, 5 xfailed, 5
errors`, and it came from `test-bare` ALONE — the `test` job had died at shellcheck and never reached
pytest at all. So the figures that were discussed described a fraction of the suite, and the jobs
nobody had looked at had not run.

What makes it hard to see is that the summary is perfectly accurate about what it covers. Nothing in
it is wrong. It simply does not say which jobs produced it, and a job that dies in an early step
produces no test output to be missing from the total — the absence has no line of its own.

**Read which jobs RAN before reading what they said.** A red run is first a question about the matrix:
which jobs started, which reached their test step, which exited before it. Only then does a count mean
anything, and a count quoted without its job name is a count without a denominator.

## A resource guard inside a test reads the test process, not the work (2026-10-03)

Memory ceilings were written against `resource.getrusage(resource.RUSAGE_SELF).ru_maxrss` from inside
a test. The reasoning was sound and the instrument could not carry it: `ru_maxrss` is a **process
high-water mark**, so it is monotonic and never comes down. A ceiling on it is therefore a ceiling on
**every test that ran before it in the same process**, and it fails in a full suite for memory that
some other test used — green alone, red together, with nothing about the code under test having
changed.

The repair is injection: the figure under test is passed in and asserted on, so the assertion is about
the work rather than about the process that happens to be hosting it. In one case the honest quantity
was different again — bytes read from `data/` — and measuring that instead removed the dependence on
the process entirely. `tests/guard_injection.py` and `tests/test_loci_gene_input.py` carry the
reasoning beside the code.

A latent instance of the same ceiling survived in `tests/test_argmaxcell.py` **only because of
alphabetical collection order**: it happened to run early enough that the high-water was still low. It
had never passed for a reason.

This is the third plant for this class, and the plant is the part that keeps being got wrong. Showing
that the guard fires on a figure above the ceiling proves the comparison, not the measurement.
**Verify against a guard that IGNORES the injected figure** — feed it a number that must trip it and
confirm it does not, so a guard still reading the process cannot pass itself off as one reading the
work.
