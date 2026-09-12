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
- Across species, families still overcount homeodomains (HOX, NK,
  paired-related and HD-LIM share the TAAT core): count sites, hits mapped
  to genome coordinates and merged where they overlap. Record positions in
  genome coordinates from the start; a position in a concatenated alignment
  depends on which species' blocks were concatenated.

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

## Engineering

- SBML species given as `initialAmount` inside a compartment of size ≠ 1 are concentrations in the maths, and a kinetic law's value is substance per time: integrate rate/volume. Skipping that ran a 100-species model into zeros within a second and looked like instability; it was units.
- process-bigraph registers processes with `core.register_link` and applies float updates additively; the composite ran our two engines on one clock.
- libRoadRunner, MaBoSS, CompuCell3D and biolearn (torch) have no Python 3.14 wheels yet; the in-house engines keep their interfaces so adapters can replace them.
- Streaming beats downloading: the HG002 telomere estimate cost 12 seconds and 0 bytes, the methylation clocks 54 seconds and 0 bytes for 1.2 GB of bedMethyl; the same information from BAM/CRAM would have cost 100-200 GB of disk.
- A strand-specific RNA-seq track is labelled by the read, not the transcript: ENCODE's IMR-90 total RNA-seq (ENCSR424FAZ) reads antisense, so APP sits on the "plus strand" file, and a reader that takes an exon's signal from the file named after the gene's strand sees the cell as silent everywhere. The closure test caught it (every chr21 gene at zero in one cell of four); `MeasuredRna` and `attribution/closure.py` now probe a few dozen exons in both orientations and swap when the swapped one carries twice the signal. Never trust a strand label without a housekeeping gene to check it against.
- **Ensembl's human homology dump has no mouse (2026-09-12).** `homologies/homo_sapiens/Compara.116.protein_default.homologies.tsv.gz` lists 199 species, among them *Mus caroli*, *Mus spretus* and the rat, but not *Mus musculus*; the mouse–human pairs are in the mouse dump (`homologies/mus_musculus/...`, `homology_species == homo_sapiens`). A reader that filters on the species it expects and finds zero rows should say so loudly; the first Compara pass wrote an empty orthology file and the comparison silently fell back to MGI.
