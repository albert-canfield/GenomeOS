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
of a uniform sample; genome-wide, 2,305 elements later, 75.7% against 62.4%,
with the strongest effects among the constrained (31.9% at or above 0.3 log2
against 24.5%) and the node model holding as well on them as on the rest
(89.8% inside their node). Constraint predicts function on 22 of 24
chromosomes (ALPHAGENOME.md, feature b, "The constrained ones first";
`data/results/constrained_targets_<chrom>.json` and
`constrained_targets_genome_wide.json`).

Scoring against measured ground truth started with VISTA's 2,442
transgenic-mouse enhancer elements, read blind through the registry, the node,
constraint and the AlphaGenome deletion. Over 23 chromosomes and 2,223
elements every layer leans the right way: an ENCODE element sits on 84.5% of
positives against 66.1% of negatives, the deletion moves a gene for 72.7%
against 64.9% with the difference concentrated in the strong effects (52%
against 39%), and the predicted tissue falls in the measured tissue group for
68% of the judged cases against 35.5% by chance. Constraint separates least,
because VISTA chose its elements for conservation to begin with. The deletion
effect is a target finder first and an activity call second (ALPHAGENOME.md
"Against measured enhancers"; `vista_<chrom>.json`, `vista_genome_wide.json`).

GTEx eQTLs test the target itself. With the match defined by the eQTL's gene
and the nearest coding TSS as the control, the deletion model's coding target
is an eGene in 71.7% of the uniform elements that carry an eQTL against 66.2%
for the nearest TSS, and where the two disagree the model is right 364 times
to the node's 309. On the most constrained elements the two are level, and
fewer of them carry an eQTL at all, since constrained sequence carries fewer
common variants (ALPHAGENOME.md "Against measured targets";
`eqtl_targets.json`).

Activity was scored last and closes the series: ENCODE4's lentiMPRA library,
51,376 elements assayed in K562, HepG2 and WTC11, against the registry,
constraint, the measured reader and AlphaGenome's predicted DNase on the cell's
own track. The registry orders activity (promoter-like 27.5% active against 6.9%
with no class in K562) but a distal enhancer-like element is barely more active
than none, since the reporter measures the sequence out of its node; constraint
adds a point or two; the predicted reader on the cell's own track correlates
with activity at 0.44 in K562 and 0.32 in HepG2 against 0.03 and 0.12 on the
other cell's track, and picks the active line for 84.6% of elements active in
exactly one. The model knows the cell, not the amount (ALPHAGENOME.md "Against
measured activity"; `mpra_genome_wide.json`).

Consequence is the weakest test by construction. GWAS lead variants sit in
regulatory elements barely above a shifted control (1.1 to 1.2 times), a named
target and a strong effect move it a few points the right way, and constrained
elements hold fewer because selection removes the common variants a GWAS
needs; ClinVar's non-coding pathogenic variants number a few dozen per set.
The measured layers, VISTA for tissue, eQTL for target, lentiMPRA for
activity, the reader for openness, are where the attribution is tested
(ALPHAGENOME.md "Against consequence"; `consequence_targets.json`).

The chromosomes sort themselves without being told. chrY and chrX are the
controls: 0.66% and 0.88% constrained, fossils 27% and 73% of their space.
The gene-dense chr17, chr19 and chr20 put half their space in the regulatory
tier. The gene-poor chr4, chr13 and chr18 carry the most neutral and
constrained-unknown sequence. Per-chromosome budgets are in
`data/results/budget_<chrom>.json`, the sum in `budget_genome_wide.json`, and
the Progress tab shows both.

## The attributions as a program (2026-09-12)

`genomeos budget --chrom chr21 --bio` compiles the chromosome's attributions into
BioLang (`genomeos/attribution/compile.py`), and the result is committed as
`data/organisms/human/noncoding_chr21.bio`, where `bio test` checks it with the
other organism programs.

- Every UNKNOWN block is a `region` whose role is its tier and label, with the
  constraint numbers and the sequence class as evidence. The constrained_unknown
  tier keeps `role: unknown`, so the program's own count of unknowns is the
  chromosome's real unknown, and a `# test:` line pins it.
- Every regulatory element with a predicted coding target, from the
  constrained-target run and the uniform enhancer-deletion run, is an
  `element` with its class, locus, domain, target and basis (the deletion's
  gene, magnitude and tissue, the constrained fraction, the node verdict), and
  one `rule` that activates or inhibits the target with the predicted magnitude
  as strength. Predicted evidence is capped at 0.7. Target genes are declared
  once as stubs.
- The CTCF nodes those elements sit in are `domain` blocks, each preceded by a
  comment with the reader's view: in which cell types the node is open and in
  which it is silent.

Chromosome 21 compiles to 747 entities (446 regions, 155 elements, 83 genes,
63 domains) and 155 rules with 20 unknowns, and passes its three checks in
0.13 s. Mean confidence is 0.62 for regions and 0.25 for elements, which is
what the evidence supports and is stated rather than rounded up. A program that
imports it can knock an element out, run the network, and see which gene moves,
which is what the closure tests in the next section will do.

## Closure at the gene level (2026-09-12)

`genomeos closure --chrom chr21` is the first of the three closure tests, the one
with enough bandwidth to falsify. For each coding gene and each of four ENCODE
cell lines that have both a DNase reader and a total RNA-seq track (K562, HepG2,
GM12878, IMR-90) it reads three things independently:

- the promoter's openness, a DNase peak within 1 kb of the TSS (the reader);
- the regulatory input from the elements attributed to the gene: an element is
  active where a DNase peak overlaps it, and the input is the signed sum of the
  predicted magnitudes over the active elements;
- the measured expression, the fraction of the canonical exons' bases carrying
  RNA-seq signal on the gene's strand, read through the bigWig range reader
  (2.9 MB for the chromosome).

| cell | expressed, promoter open | expressed, promoter closed | open with activating input | open, no element |
|---|---|---|---|---|
| K562 | 64% (135) | 7% (86) | 57% (7) | 65% (122) |
| HepG2 | 73% (118) | 8% (103) | 62% (8) | 76% (107) |
| GM12878 | 84% (89) | 22% (132) | 100% (3) | 85% (85) |
| IMR-90 | 76% (130) | 2% (91) | 60% (15) | 77% (108) |

The reader's claim holds in every cell. The elements add nothing at this
resolution: within a cell they do not raise the expressed fraction, and across
cells, for the 39 genes whose input and expression both vary, the cell where the
elements are most active is the most expressed cell 15% of the time, against 25%
by chance and 33% under shuffled cells. Twelve attributions are rejected by name
(MAP3K7CL in HepG2 with input 0.77 and no expression; SIM2 and EVA1C in K562;
OLIG1 in IMR-90), and 43 expressed genes have a closed promoter and no active
element, most of them in GM12878, whose reader has the fewest peaks.

The reading is not that the targets are wrong. VISTA and the eQTLs say the
target and the tissue are usually right. It is that a target plus "active where
a DNase peak overlaps" plus a magnitude scored in whichever tissue moved most
does not reproduce a cell. The deletion scored on each cell's own track was
the obvious repair, and it ran the same day: across cells, over 71 genes, the
most active cell is the most expressed 24% of the time, the shuffled rate
exactly. Within a cell the sign carries a little, since open-promoter genes with
a repressing input are expressed less often than those with an activating one
in all four lines (by 5 to 21 points): the model's direction is right where its
magnitude is too small to rank cells. The negative stands and is sharper for
having the right magnitudes. Only 83 of the 221 genes have any attributed
element and most have one, with effects of 0.1 to 0.2 log2; a gene's expression
in a cell is set by its promoter and by the elements the two sampled runs never
scored. An attribution of "this element, this gene" is not yet an attribution of
"this cell". Scoring every element in a gene's node per cell, about 6,600 on
chromosome 21 alone, is the self-hosted model's first job.

One measurement lesson on the way: IMR-90's ENCODE total RNA-seq carries its
strand labels inverted (APP reads on the "plus" file), so the module probes each
cell's orientation over forty exons per strand and swaps the tracks when the
swapped orientation carries more than twice the signal, recording the choice.

## Syntax against values (2026-09-12)

Albert's framing: a language has syntax, operators and values, and a variable
trait should show the same syntax in every person with only the value changing.
`genomeos syntax --gene HERC2 --chrom chr15` reads one gene three ways at once:
Zoonomia constraint per base (what selection held still across 241 mammals, the
syntax), the variants every imported person carries there (where people differ,
the values), and the gene's features (canonical exons, cCREs). Each variable
position is then a value or a value in syntax, the rare change on protected
sequence.

| gene | span | syntax bases | variable positions (SNV) | values in syntax | if variation ignored syntax |
|---|---|---|---|---|---|
| HERC2 | 211 kb | 11,817 (5.6%) | 160 (0.76 per kb) | 3 (1.9%) | 5.6% |
| OCA2 | 345 kb | 2,995 (0.9%) | 544 (1.58 per kb) | 3 (0.6%) | 0.9% |

Two things the table says. Variation avoids syntax: in both genes the persons'
variants land on constrained bases less often than they would by chance, which is
the language claim in one number. And the values in syntax are the ones worth
asking about: in HERC2 the second most constrained variable position of the whole
gene, phyloP 3.41 in an intron, is rs12913832, the base inside the HERC2 enhancer
that sets OCA2 expression and with it eye colour. HG002 and HG003 carry the
alternative allele twice, HG004 once. The tool did not know the answer; the
literature's value fell out of the reading as the position a language-minded
reader would ask about first. The other two HERC2 values in syntax are coding.

The command works for any gene and any set of imported genomes (`--name`), reads
about 1 to 2 MB of ranges, and saves `data/results/syntax_<GENE>.json`. With more
genomes imported the value slots fill in; with the reader and the attributions the
next question, what a slot controls, is asked of the right positions.

## What comes next, in order

1. **Every element in a gene's node scored per cell**, not a sample, so the
   closure judges a gene's whole regulatory input; the self-hosted model's
   first job.
2. **The rest of the organising evidence** on the compiled blocks: repeat
   family and segmental duplication on regions, the reader's openness on the
   element itself rather than its node, and the human axis (area J's Gnocchi
   score) beside the mammalian one. The two axes disagreeing is itself a
   flag: chr15:84,395,903-84,398,315 is 23% constrained across mammals and
   carries far more human variation than expected (Gnocchi Z of -8.5 over its
   kilobases), which reads as a mutation hotspot or a mapping artefact rather
   than a value slot, and the constrained_unknown tier holds 18 such blocks on
   chr15 to look at first.
3. **The remaining attribution**: AlphaGenome's predicted DNase, histone and
   CAGE tracks say whether a block is active and where, and in-silico
   mutagenesis says which bases inside it matter; both on the 15,536
   regulatory blocks. The hosted service handles
   thousands of predictions, not a million blocks, so the self-hosted model in
   the roadmap pool gates the move from chromosome 21 to the genome.
4. **Score against measured ground truth**: done for tissue (VISTA), target
   (eQTL), activity (lentiMPRA), openness (the readers) and consequence
   (GWAS, ClinVar), each recorded above with its precision and its limit.

5. **Closure tests at three levels.** Gene level: open elements with their
   predicted effects must reproduce GTEx expression per tissue (19,000 genes by
   54 tissues). Cell level: an element deleted in silico moves its gene and
   changes the cell type's program; compare with IMPC knockouts and DepMap.
   Organism level: the Body runtime with its turnover timers must still
   produce a human's cell counts and proportions over decades. The last one
   validates the mechanism modules, not single blocks, and stays as the sanity
   check Albert described.
6. **Simpler genomes as the comparison.** C. elegans first (GenomeOS already
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
