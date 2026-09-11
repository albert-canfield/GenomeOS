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
