# Investigating the UNKNOWN blocks

**Status on 2026-10-10.**
- The classification below holds against `unknown_genome_wide.json` (2026-09-10): 26,806 UNKNOWN blocks, 1,008,756,271 bases (33% of hg38), 98.53% carrying a class. A class here comes from sequence, RepeatMasker and ENCODE chromatin evidence; it names what kind of sequence a block is, not what it does.
- The next layer is in docs/ATTRIBUTION.md: every block carries a constraint tier (`budget_genome_wide.json`; of the 1,008.8 Mb, structural 22.8%, repeat_unconstrained 32.6%, regulatory 34.2%, unconstrained_unknown 7.3%, constrained_unknown 3.2%; since review R7 of 2026-09-28 the stored keys `fossil` and `neutral` read under those names). The "real unknown" is the 882 constrained-unknown blocks that are not mostly copies: 30,602,182 bases.
- Measured: counting a block as measured if lentiMPRA, VISTA or the CRISPRi benchmark overlaps it by a single base, 4,554,153 of the 1,008,756,271 unknown bases (0.45%), and 160,447 of the real unknown's 30,602,182 bases (0.52%, in 202 of 882 blocks) (`unknown_coverage.json`, 2026-09-17).
- Nothing in this document is withdrawn. The reading of 2026-09-16 that the tiers differ in function was withdrawn on 2026-09-17 (ATTRIBUTION.md, "The matched table, computed correctly").
- Milestone 1.3 (the 98%), clause 2, is not met: on chr21, 140 of 285 unknown blocks (49.1%) name a coding gene against 10,908 of 16,706 length-matched random windows (65.3%) (`unknown_scoring_chr21.json`); this is how often a target is named, and no measured outcome enters it. Clause 3's designed experiment, 284,001 oligos over 878 blocks, is built and unrun (`unknown_library.json`).
- The last section, "What it cannot do yet", carries dated notes for what has been built since.

Forty-five percent of chromosome 21 is UNKNOWN in the block map: the space
between annotated genes. It is not empty; it is where our knowledge is
thinnest. `genomeos unknown` classifies it from the sequence itself, largest
block first, with the evidence that produced each label.

## How it classifies

| class | signal | evidence |
|---|---|---|
| gap | run of N | curated (assembly) |
| telomere | TTAGGG arrays covering the block | curated pattern |
| centromere | CENP-B boxes at > 1 per kb, or ~171 bp periodicity of the sequence | curated pattern / inferred |
| tandem_repeat | short tandem repeats over half the block | predicted |
| low_complexity | low 3-mer entropy over half the block | predicted |
| interspersed_repeat | more than 60% of positions carry a 16-mer that occurs ≥ 40 times in the chromosome | predicted (learned from the chromosome's own k-mers) |
| interspersed_repeat_SINE | as above, with the Alu 5' core present | curated pattern |
| long_orf | an open reading frame of 250–2,500 aa in non-repetitive sequence: an unannotated gene, a pseudogene, or a transposon's own ORF (L1 ORF2 is ~1,275 aa) | predicted |
| satellite_array | a stop-free frame longer than 2,500 aa: an array of a repeat unit with no stop in one frame | inferred |
| promoter_like | a CpG island in a block ≤ 20 kb, or ≥ 2 islands per 100 kb | predicted |
| mixed_intergenic | blocks ≥ 20 kb are read in 5 kb windows; a mix of interspersed repeats and unique sequence | predicted |
| unique_intergenic | ≥ 70% of windows unique and non-repetitive: the space where regulation is expected | inferred |
| regulatory | ENCODE candidate regulatory elements (promoter-like, enhancer-like) at ≥ 2 per 10 kb: experimental chromatin evidence, checked before every sequence-only rule | curated |
| gene_desert | ≥ 500 kb with nothing above | inferred |
| unclassified | none of the above dominates | none |

Two design choices matter. Repeats are **learned, not assumed**: the 16-mer
table is built from the chromosome being investigated, so any element with
many copies shows up without a repeat library. Names are **only ever set by
patterns**: the file `data/patterns/known_motifs.tsv` (name, class, regex,
evidence, confidence, note) is the place to declare what a motif is and how
sure you are; the tool reports pattern hits per block and uses them to name
classes (Alu core, CENP-B box, telomeric repeat). Add your own rows to teach
it more.

Among the 60 largest blocks, pairs that share many 20-mers are reported as
duplication candidates (`similar_to`).

## Curated repeats: RepeatMasker through the UCSC API

Sequence patterns only find young repeats. `genomeos repeats --chrom chr21`
fetches RepeatMasker's annotation of the chromosome from the public UCSC
track API (20 MB of JSON once), keeps a 0.7 MB BED with class, family, name
and divergence, and a class summary (`rmsk_chr21.json`: 64,812 copies,
20.5 Mb; LINE/L1 6.3 Mb, SINE/Alu 3.5 Mb, alpha satellite 2.5 Mb). The
UNKNOWN classifier then grades each block by curated coverage before any
sequence rule: ≥ 50% interspersed repeat → `interspersed_repeat_<class>`
(curated, 0.9); ≥ 50% satellite → centromere or satellite array; ≥ 50%
simple or low-complexity → tandem repeat; and a block that is ≥ 30%
interspersed repeat can no longer be called `long_orf`.

On chromosome 21 this moved 2.1 of the 2.3 Mb of `long_orf` into
`interspersed_repeat_LINE` (3.2 Mb in total, the LINE-1 ORF2 frames), added
0.45 Mb of LTR, and raised the classified share from 96.3% to 97.6%. What
remains `long_orf` (0.13 Mb) is non-repetitive and worth a look.

## Chromosome 21 (85 s, with ENCODE evidence)

| class | blocks | bp | share of UNKNOWN |
|---|---|---|---|
| gap (assembly N runs) | 15 | 6.16 Mb | 29.6% |
| regulatory (ENCODE cCRE-dense) | 177 | 3.65 Mb | 17.5% |
| unique intergenic | 52 | 2.95 Mb | 14.2% |
| centromere (CENP-B boxes, 81% high-copy k-mers) | 1 | 2.35 Mb | 11.3% |
| long ORF (transposon ORFs, pseudogenes, unannotated) | 38 | 2.27 Mb | 10.9% |
| satellite array | 6 | 1.22 Mb | 5.9% |
| mixed intergenic | 9 | 0.76 Mb | 3.6% |
| interspersed repeat | 7 | 0.41 Mb | 2.0% |
| promoter-like (CpG island) | 17 | 0.21 Mb | 1.0% |
| unclassified | 113 | 0.76 Mb | 3.7% |

96.3% of the UNKNOWN bases now carry a class. Chromosome 21 has 13,356 ENCODE
elements (9,981 distal enhancer-like, 2,158 proximal, 482 promoter-like, 385
CTCF-only, 350 open-chromatin), kept as a 150 KB file distilled from the
64 MB registry.

Two lessons learned on the way. Exact 16-mer counting finds only *young*
repeat copies; most Alu and L1 copies are 10–20% diverged, so the first
version saw the intergenic space as 99% unique, and mismatch-tolerant
matching of the Alu 5' core (up to 7 substitutions in 36 bases, zero hits in
200 kb of random sequence) recovered them. And sequence alone cannot see
regulation: adding one experimental layer (chromatin) moved 17.5% of the
UNKNOWN space from "unique" to "regulatory" with curated evidence.

## The whole genome (84 min of compute, one chromosome on disk at a time)

`unknown_genome_wide` streamed every chromosome from UCSC, fetched its
RepeatMasker annotation and ENCODE elements, classified its UNKNOWN blocks
with the curated evidence first, kept the summary and deleted the rest.
chrM has no UNKNOWN space.

| total | value |
|---|---|
| UNKNOWN blocks | 26,806 |
| UNKNOWN bases | 1.009 Gb (33% of hg38) |
| classified | 98.5% (was 96.6% before the curated repeat pass) |

| class | bases | share of UNKNOWN |
|---|---|---|
| regulatory | 338.8 Mb | 33.6% |
| interspersed repeat LINE | 287.1 Mb | 28.5% |
| gap | 159.7 Mb | 15.8% |
| unique intergenic | 74.8 Mb | 7.4% |
| centromere | 48.2 Mb | 4.8% |
| interspersed repeat LTR | 22.5 Mb | 2.2% |
| satellite array | 21.5 Mb | 2.1% |
| interspersed repeat SINE | 17.8 Mb | 1.8% |
| unclassified | 14.8 Mb | 1.5% |
| mixed intergenic | 10.0 Mb | 1.0% |
| promoter like | 6.4 Mb | 0.6% |
| interspersed repeat | 5.0 Mb | 0.5% |
| long orf | 1.2 Mb | 0.1% |
| interspersed repeat DNA | 0.6 Mb | 0.1% |

*(Header row added 2026-10-10; the rows are unchanged and agree with `unknown_genome_wide.json` → `chromosomes`.)*

| chromosome | UNKNOWN blocks | UNKNOWN bases | classified | largest class |
|---|---|---|---|---|
| chr1 | 2,425 | 77.2 Mb | 99% | regulatory |
| chr2 | 1,943 | 64.1 Mb | 98% | regulatory |
| chr3 | 1,298 | 50.4 Mb | 99% | regulatory |
| chr4 | 1,398 | 62.2 Mb | 99% | interspersed repeat LINE |
| chr5 | 1,363 | 54.7 Mb | 99% | regulatory |
| chr6 | 1,396 | 48.0 Mb | 99% | regulatory |
| chr7 | 1,400 | 45.2 Mb | 98% | regulatory |
| chr8 | 1,180 | 40.3 Mb | 98% | regulatory |
| chr9 | 1,080 | 49.1 Mb | 99% | gap |
| chrX | 1,338 | 71.7 Mb | 99% | interspersed repeat LINE |
| chrY | 375 | 46.7 Mb | 99% | gap |
| chr10 | 1,190 | 38.7 Mb | 98% | regulatory |
| chr11 | 1,453 | 39.0 Mb | 97% | regulatory |
| chr12 | 1,165 | 35.3 Mb | 99% | regulatory |
| chr13 | 771 | 52.2 Mb | 99% | gap |
| chr14 | 942 | 40.7 Mb | 99% | gap |
| chr15 | 684 | 33.8 Mb | 98% | gap |
| chr16 | 862 | 29.4 Mb | 98% | gap |
| chr17 | 1,022 | 20.1 Mb | 98% | regulatory |
| chr18 | 612 | 29.6 Mb | 99% | regulatory |
| chr19 | 1,066 | 14.9 Mb | 98% | regulatory |
| chr20 | 762 | 21.6 Mb | 99% | regulatory |
| chr21 | 446 | 20.8 Mb | 98% | gap |
| chr22 | 635 | 22.9 Mb | 99% | gap |

What the table says. A third of the space between genes is regulatory by
ENCODE's chromatin evidence; another third is LINE-1 and its relatives by
RepeatMasker's curated library, which is where the earlier "long ORF"
class went (it read the ORF2 frames of LINE-1 copies); gaps are the
unassembled 16%; one in thirteen bases is unique intergenic sequence with
no evidence of anything yet, which is the honest residue this tool exists
to shrink. Only 1.5% stays unclassified.

## Where it shows

- `genomeos unknown --chrom chr21` runs it and saves `data/results/unknown_chr21.json`
  (per block: class, evidence, confidence, features, pattern hits, similar blocks).
- The **Blocks** tab labels UNKNOWN blocks with their class and colours them
  by it; the "Investigate UNKNOWN" button starts the job, which appears in
  **Progress**.

## What it cannot do yet

Enhancers, silencers and insulators leave no signal a k-mer or a regex can
read; ENCODE's chromatin evidence now places them (`regulatory` class,
`genomeos regulation` for their targets), but which gene each one reaches
is inferred from the CTCF domain, not measured. Where the AlphaGenome job
has run (`enhancer_targets_chr21`: 200 chr21 enhancers deleted one by one,
docs/ALPHAGENOME.md feature b) the element carries a `predicted` target with
a tissue and a magnitude next to the inferred one; 87% of the predicted
coding targets fall inside the inferred node. Hi-C would make it a
measurement.
*(2026-10-10: measured since, in part. 4DN boundary calls are read by `genomeos domains --hic`
(docs/NODES-READER-WRITER.md); on 661 CRISPRi Regulated=TRUE element-gene pairs the measured pair
sits inside one CTCF node 5.89 points above the published uniform control, 95% +3.18 to +8.46 over
chromosomes (`node_containment_measured.json`). The AlphaGenome sweep now covers all 961,227
registry elements of 24 chromosomes, not 200 (`enhancer_targets_all_genome_wide.json`). Outside the
elements a CRISPRi screen tested, an element's target is still predicted or inferred, not measured.)*
Repeat families are
named by RepeatMasker's curated library where a chromosome has been
distilled with `genomeos repeats`; the sequence patterns remain the fallback
for a chromosome that has not. Silencers have no curated source here yet.
*(2026-10-10: one exists since 2026-10-03, ReSE's assay-validated silencer records, which validate
in K562 and HepG2 only. Against them the predicted `silencer_like` label, a direction of effect, is
not enriched: 420 of 1,181 joined elements, 35.563% against a stratum-matched 37.707%, difference
−2.144 points, 95% [−5.897, +1.158], registered band NO_SIGNAL (`silenceragree.json`).)*
