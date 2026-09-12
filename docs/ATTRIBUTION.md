# The 98%: attributing function to the non-coding genome

Area I of ROADMAP.md. Code: `genomeos/attribution/` (`bigwig.py`, `constraint.py`,
`budget.py`), `scripts/budget_genome_wide.py`, `genomeos budget`. Results:
`data/results/budget_<chrom>.json`, `budget_genome_wide.json`.

## The question, and the house

The coding 2% is decoded: every gene translated, verified against UniProt,
compiled into a protein definition, packaged. The other 98% is classified by
sequence class (`genomeos unknown`: repeat families, regulatory clusters,
centromere, satellite, unique intergenic) but nobody has said what each block
*does for the organism*. Albert's framing on 2026-09-11: the DNA is the
blueprint, the living organism is the built house, and a plan that calls half
the house pipes is wrong before a brick is laid because the quantities do not
fit. Attribute a function to each block, run the organism forward, and let the
quantities and the outcome falsify the attributions.

This document is the plan for that, and the record of what has been built.

## What the quantities already say

Two published constraints act exactly like "50% can never be pipes".

- **Constraint across mammals.** Zoonomia aligned 241 placental mammals
  (Christmas et al. 2023, Science) and scored every human base with phyloP.
  A base at phyloP ≥ 2.27 is constrained at 5% FDR; about a tenth of the
  genome passes. Everything else has been free to mutate for tens of millions
  of years, and a base that can change without consequence is not carrying a
  function the organism depends on.
- **Mutational load.** If much more than a tenth of the genome mattered, each
  child's roughly 70 new mutations would include too many harmful ones for the
  population to persist (Graur 2017). This is a closure constraint at the level
  of the species, the same argument as the house budget.

So the honest split is not 2% known and 98% unknown. Roughly 85% of the genome
is repeat-derived or unconstrained, and there the best guess "fossil, fill,
structural" is a good guess with evidence attached. Roughly 10% is constrained
non-coding, and that is where the attribution problem lives.

## The tiers

Every UNKNOWN block gets a tier, a label and a confidence from its sequence
class plus its constraint. The tiers are the vocabulary of the next steps, not a
verdict.

| tier | what it means | rule (constrained fraction of the block's bases) |
|---|---|---|
| structural | assembly gap, centromere, satellite or tandem array: a mechanical role, not read as genes | by class |
| fossil | transposable-element remains without constraint | repeat class, < 3% (0.8) or 3–5% (0.6) |
| regulatory | ENCODE-backed element cluster; the target gene is the open question | regulatory or promoter-like class; ≥ 5% raises the confidence to 0.7 |
| constrained_unknown | under selection, not coding, not regulatory by the registry: the real unknown | ≥ 5%, or ≥ 3% with three or more conserved elements; also repeat-derived blocks ≥ 5% (possibly exapted) |
| neutral | unique sequence with no constraint: best guess, nothing | < 3% (0.6); 3–5% weakly constrained (0.4) |

Two thresholds carry the rules. Below 3% of bases constrained a block reads as
unconstrained: neutral sequence reaches 1 to 2% by chance at a 5% FDR, and the
whole of chromosome 21's UNKNOWN space sits at 1.3%. From 5% up the block
carries constraint worth a name. Both are constants in `budget.py`, stated in
every result, and open to revision once more chromosomes are in.

## Sources, and how they are read

- **Zoonomia phyloP 241-way** is a 9.6 GB bigWig at UCSC. It is never
  downloaded. `bigwig.py` is a standard-library reader that fetches the
  chromosome tree and the R-tree index with HTTP range requests, coalesces the
  data sections overlapping the blocks into a few ranged fetches, inflates
  them, and keeps per block: bases seen, mean, maximum, and how many reach the
  threshold. Chromosome 21's 431 sequence-bearing blocks cost 294 requests
  and 33 MB in under four minutes. Assembly gaps are not asked for.
- **100-vertebrate conserved elements** (`phastConsElements100way`) come
  through the UCSC REST API in 5 Mb windows, one request per second, and are
  cached once per chromosome under `data/knowledge/constraint` (git-ignored,
  rebuilt on demand). Per block: elements overlapping, bases covered, highest
  lod.

Both are evidence of purifying selection, which is the nearest thing to a
measurement of "this matters" that exists for every base. Neither says what
the base does.

## Chromosome 21, the first budget (2026-09-11)

| tier | blocks | Mb | of UNKNOWN | of chromosome | constrained kb |
|---|---|---|---|---|---|
| structural | 27 | 9.75 | 46.9% | 20.9% | 0.8 |
| fossil | 104 | 3.91 | 18.8% | 8.4% | 39.9 |
| regulatory | 194 | 3.93 | 18.9% | 8.4% | 50.6 |
| constrained_unknown | 20 | 0.29 | 1.4% | 0.6% | 18.0 |
| neutral | 101 | 2.93 | 14.1% | 6.3% | 39.8 |

Of 20.8 Mb of UNKNOWN blocks on chromosome 21, 1.3% of the measured bases
are constrained, against 10.7% for the genome as a whole: constraint lives in
and around genes, which is where the UNKNOWN pass does not look. Twenty blocks
totalling 287 kb are the real unknown of this chromosome. The most constrained
of them, 47 kb at 16.82 Mb, is unique intergenic sequence with 179 conserved
elements and a lod of 2,163, which is what an unannotated regulatory region or
an unannotated gene looks like. The registry's regulatory class shows only
1.2% constrained bases, a reminder that a cCRE is a class of chromatin state
and mostly lineage-specific, not a conserved element.

## The genome (2026-09-12)

The job ran the 24 chromosomes in eight hours, smallest first, 2,252 MB of
range requests over a track that was never downloaded.

| tier | blocks | Mb | of UNKNOWN | of genome | constrained Mb |
|---|---|---|---|---|---|
| structural | 238 | 229.7 | 22.8% | 7.4% | 0.05 |
| fossil | 7,302 | 328.5 | 32.6% | 10.6% | 4.13 |
| regulatory | 15,536 | 345.2 | 34.2% | 11.2% | 8.10 |
| constrained_unknown | 1,098 | 32.0 | 3.2% | 1.0% | 1.66 |
| neutral | 2,632 | 73.3 | 7.3% | 2.4% | 1.07 |

1,009 Mb of UNKNOWN blocks, 780 Mb of them holding sequence, 15.0 Mb of
those constrained: 1.93%, against 10.7% for the genome as a whole. No class
of the UNKNOWN space comes near the genome's average (centromere 0.6%,
unique intergenic 2.5%), which is the quantitative form of "constraint lives
in and around genes". The real unknown is 1,098 blocks and 32 Mb, one percent
of the genome; the regulatory tier holds five times more constrained sequence
and its target genes are the next attribution. That attribution has started:
ranking chromosome 21's distal enhancers by constraint and deleting the 100
most constrained in AlphaGenome names a gene for 73% of them, against 63.5%
of a uniform sample, and the constrained ones act beyond their CTCF node twice
as often (ALPHAGENOME.md, feature b, "The constrained ones first";
`data/results/constrained_targets_<chrom>.json`).

The chromosomes sort themselves without being told. chrY and chrX are the
controls: 0.66% and 0.88% constrained, fossils 27% and 73% of their space.
The gene-dense chr17, chr19 and chr20 put half their space in the regulatory
tier. The gene-poor chr4, chr13 and chr18 carry the most neutral and
constrained-unknown sequence. Per-chromosome budgets are in
`data/results/budget_<chrom>.json`, the sum in `budget_genome_wide.json`, and
the Progress tab shows both.

## What comes next, in order

1. **Organise the blocks** by the evidence that already exists: CTCF node
   membership, the reader's open fraction per cell type, the AlphaGenome
   enhancer target, repeat family, segmental duplication, and now constraint.
2. **Attribution as BioLang `element` blocks**: a function, a target gene, a
   tissue and a confidence labelled `predicted`, from priors per tier.
   AlphaGenome's predicted DNase, histone and CAGE tracks say whether a block is
   active and where; deleting it in silico says which gene moves; in-silico
   mutagenesis says which bases inside it matter. The hosted service handles
   thousands of predictions, not a million blocks, so the self-hosted model in
   the roadmap pool gates the move from chromosome 21 to the genome.
3. **Score against measured ground truth**, the way the segment parser was
   scored against GENCODE: VISTA enhancers (about 3,300 tested in mouse
   embryos), ENCODE4 lentiMPRA (hundreds of thousands of elements in K562,
   HepG2, WTC11), GTEx eQTL, the GWAS catalog and ClinVar. Precision and
   sensitivity reported plainly per tier.
4. **Closure tests at three levels.** Gene level: open elements with their
   predicted effects must reproduce GTEx expression per tissue (19,000 genes by
   54 tissues). Cell level: an element deleted in silico moves its gene and
   changes the cell type's program; compare with IMPC knockouts and DepMap.
   Organism level: the Body runtime with its turnover timers must still
   produce a human's cell counts and proportions over decades. The last one
   validates the mechanism modules, not single blocks, and stays as the sanity
   check Albert described.
5. **Simpler genomes as the comparison.** C. elegans first (GenomeOS already
   grows it; the lineage is the finished house), Drosophila for enhancer
   grammar (genome-wide STARR-seq), fugu as the compact vertebrate with the
   same gene count and an eighth of the sequence, mouse (already in) for
   transfer. The minimal synthetic cell JCVI-syn3.0 stays as the humbling
   benchmark: 473 genes, 149 of them still without a known function.

## A note on the alphabet

Four letters, not ten, and that helps at every layer. Sixty-four codons for
twenty amino acids and a stop; a 3-periodic fifth-order Markov model has 4⁶
contexts per frame position, which is small enough to learn from one
chromosome; the parser's signals are matrices over four symbols. For the
UNKNOWN blocks the same fact is already at work in the classifier: k-mer
composition (the high-copy 16-mers, GC and CpG) is what separates a LINE
fossil from unique sequence. The alphabet does not tell a block's function,
but it makes the space of hypotheses small enough to enumerate, which is what
the tiers do.
