# Investigating the UNKNOWN blocks

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

## Chromosome 21, first run (55 s)

| class | blocks | bp | share of UNKNOWN |
|---|---|---|---|
| gap (assembly N runs) | 15 | 6.16 Mb | 29.6% |
| unique intergenic | 73 | 3.82 Mb | 18.4% |
| long ORF (250–2,500 aa; transposon ORFs, pseudogenes, unannotated) | 75 | 3.07 Mb | 14.7% |
| centromere (CENP-B boxes, 81% high-copy k-mers) | 1 | 2.35 Mb | 11.3% |
| satellite array (stop-free frame > 2,500 aa) | 6 | 1.22 Mb | 5.9% |
| mixed intergenic | 14 | 1.06 Mb | 5.1% |
| interspersed repeat | 19 | 0.84 Mb | 4.0% |
| promoter-like | 46 | 0.83 Mb | 4.0% |
| unclassified | 175 | 1.29 Mb | 6.2% |

Inside blocks of 20 kb or more, 28% of 5 kb windows are interspersed repeat
and 72% unique; Alu density has a median of 0.18 per kb and peaks at 2 per kb.
The centromere of chromosome 21 (chr21:10.65–13.0 Mb) is found from its CENP-B
boxes alone.

A lesson learned on the way: exact 16-mer counting finds only *young* repeat
copies. Most Alu and L1 copies are 10–20% diverged, so their 16-mers are no
longer shared and the first version saw the intergenic space as 99% unique.
Mismatch-tolerant matching of the Alu 5' core (up to 7 substitutions in 36
bases, zero hits in 200 kb of random sequence) is what recovers the old copies.

## Where it shows

- `genomeos unknown --chrom chr21` runs it and saves `data/results/unknown_chr21.json`
  (per block: class, evidence, confidence, features, pattern hits, similar blocks).
- The **Blocks** tab labels UNKNOWN blocks with their class and colours them
  by it; the "Investigate UNKNOWN" button starts the job, which appears in
  **Progress**.

## What it cannot do yet

Enhancers, silencers and insulators leave no signal a k-mer or a regex can
read; they need cross-species conservation, chromatin data or a predictive
model (AlphaGenome, `predicted` evidence). Repeat *families* are found but not
named beyond Alu; a learned clustering of high-copy k-mers into families, and
Dfam consensus sequences as patterns, are the next steps.
