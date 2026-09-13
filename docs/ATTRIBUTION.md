# The 98%: attributing function to the non-coding genome

Area I of ROADMAP.md. Code: `genomeos/attribution/` (`bigwig.py`, `constraint.py`,
`budget.py`, `organise.py`, `candidates.py`, `unknown_scoring.py`, `lexicon.py`),
`scripts/budget_genome_wide.py`, `scripts/syntax_candidates.py`, `scripts/lexicon.py`,
`genomeos budget`. Results: `data/results/budget_<chrom>.json`, `budget_genome_wide.json`,
`syntax_candidates_genome_wide.json`, `unknown_scoring_chr21.json`, `lexicon_<chrom>.json`.

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

Chromosome 21 first compiled to 747 entities (446 regions, 155 elements, 83
genes, 63 domains) and 155 rules with 20 unknowns, from the two samples; with
every element scored it compiles to 5,972 entities (446 regions, 5,174 elements,
193 genes, 159 domains) and 5,174 rules, still 20 unknowns, 2.7 MB, and passes
its three checks in 0.26 s. Mean confidence is 0.62 for regions and 0.25 for elements, which is
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

### The closure passes on the whole input (2026-09-13)

Every one of chromosome 21's 12,139 enhancer elements was then scored in
AlphaGenome on the four cell lines' own tracks, and the closure reran on the same
table with a fair judge: ties among cells broken at random rather than by cell
order, and a null that permutes each gene's inputs across cells 1,000 times.

| element input across the four cells | genes | most active cell is the most expressed | null | p |
|---|---|---|---|---|
| tissue-agnostic magnitude, DNase-gated | 138 | 30% | 25% | 0.07 |
| the deletion on the cell's own track | 146 | 36% | 25% | 0.002 |
| the cell's own score, DNase-gated | 141 | 41% | 25% | 0.001 |

A gene's whole regulatory input, scored in the cell's own track and read where
the cell's chromatin is open, picks the cell that makes the gene well above
chance. Within a cell it carries too: open-promoter genes with activating input
are expressed more often than those with no element in three of the four lines,
and the unexplained genes fall from 43 to 9. The tissue-agnostic magnitude still
fails, so the earlier negative was read correctly: one sampled element per gene
cannot carry a cell, the whole input can. The effect is modest, a mean rank
correlation of 0.19, and 72 attributions are rejected by name, S100B in K562
first. Those are the list to read next.

On chromosome 22, 19,708 elements scored the same way, the across-cell pass
replicates at a smaller effect: 32% of 375 genes for the cell's own DNase-gated
score against a 25% null (p 0.003), 30% (p 0.013) ungated, and 31% (p 0.003) for
the tissue-agnostic input, which passes here though not on chromosome 21. Within
a cell, repressing input lowers the expressed fraction of open-promoter genes in
three of the four lines, while in K562 and HepG2 activating input does not raise
it above genes with no element: the model's repressors carry a cell better than
its activators. The 235 rejected attributions of the tissue-agnostic reading open
with lineage genes the model activates in the wrong line, VPREB1 and VPREB3,
B-cell genes, in erythroid K562.

Read from the scorer's side, those three are an artefact of the tissue-agnostic
input: their magnitude comes from Jurkat and skeletal muscle, and their K562
values sit near zero. The closure's headline is therefore the cell's own score
counted where a DNase peak sits on the element. In that reading the rejected set
splits in two. Between 42% and 59% of it are genes no line of the panel makes,
neural (OLIG1, OLIG2, GRIK1, S100B), keratin, lens (CRYBB1, CRYBB2) and liver
(CYP2D6) genes whose tissues the four lines do not cover: not wrong targets. The
rest, 22 genes on chromosome 21 and 96 on chromosome 22, are made by another line
and are the genuine wrong-cell cases. Against the attributions the cell honours
they carry about twice as many elements, a higher summed input, and promoters
open in fewer of the four lines. The closure is summing many weak elements into a
large input for genes whose promoters are only marginally available.

That reading was then tested, and it does not hold as a fix. The strongest open
element, the mean over open elements, and the sum weighted by the promoter's DNase
percentile all score at or below the plain sum on both chromosomes, and the
rejected count barely moves. The control settles what the elements are worth: the
promoter's own DNase percentile, with no elements, picks the expressing cell 27% and
28% of the time, no better than chance, while the summed element input picks it 41%
and 32%. The elements carry which-cell information that the promoter's openness
does not.

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

## The organiser (2026-09-12)

`genomeos organise --chrom C` joins, per UNKNOWN block, the three readings that
exist: the budget's tier and mammalian constraint, area J's human axis and the
case the two axes make, and the copy flag from curated segmental duplications.
It recomputes nothing and it reads copies first. A block half or more of which is
a duplication is a copy before anything else is said about it, because its
mammalian alignment is paralogous and its human variants do not map, so both
axes lie about it. A constrained_unknown block that is not a copy is then read by
its case.

| tier | blocks | copies | after copies | syntax | relaxed | recent | tolerant | unmeasured |
|---|---|---|---|---|---|---|---|---|
| constrained_unknown | 1,098 | 216 | 882 | 69 | 437 | 29 | 329 | 18 |
| regulatory | 15,536 | 689 | 14,847 | 612 | 1,011 | 3,831 | 8,932 | 461 |
| fossil | 7,302 | 1,240 | 6,062 | 0 | 0 | 732 | 5,130 | 200 |
| neutral | 2,632 | 693 | 1,939 | 0 | 0 | 306 | 1,537 | 96 |
| structural | 238 | 56 | 182 | 0 | 0 | 1 | 37 | 144 |

The real unknown after copies is 882 blocks and 30.6 Mb, and it is three things
rather than one. Sixty-nine blocks, 387 kb in all, are constrained on both axes
with no annotation: the sharpest candidates for an unannotated element or gene,
and the first to read one by one. Four hundred and thirty-seven, 13.4 Mb, are held
across mammals yet variable among people: a frame whose value varies, or a
function this lineage has lost. Three hundred and twenty-nine, 16.4 Mb, are free
on both scales at kilobase resolution, and want their alignment checked before
anything is attributed. Chromosome 21's twenty unknown blocks become five after
copies, none of them syntax. The compiled program carries the copy flag as well:
a copy's region says so in its role and its note, with the program's checks and
unknown count unchanged. Results are `organised_<chrom>.json` (tallies, copies,
candidates; the full join is recomputed on demand) and `organised_genome_wide.json`.

## The syntax candidates read one by one (2026-09-13)

`genomeos/attribution/candidates.py` (`scripts/syntax_candidates.py`,
`syntax_candidates_genome_wide.json`, six minutes) reads each of the organiser's 69
blocks, 387 kb constrained on both axes and not copies, through every layer the
repository already holds, with no model call and no remote track. Per block: the
GENCODE v50 gene on each side and which of its ends faces the block, the CTCF node
and its inferred target (nearest coding TSS inside it) and the targets the deletion
runs named for elements of the same node; how many of the Gnocchi kilobases it
touches hold a neighbour's exon bases; the 100-vertebrate conserved elements merged
into segments and the repeats they sit in; the registry and the reader (ENCODE DNase
in eleven cell types) on the conserved bases themselves; every conserved segment of
45 bp or more scored for an open frame with the chromosome's codon usage, held
against 300 canonical internal exons of the same chromosome, the segment parser over
the block with 10 kb either side, and open frames translated against the local
UniProt proteome; and VISTA, lentiMPRA, ClinVar, GWAS and GTEx where they overlap.
Every layer is read again on up to forty same-length windows in nearby intergenic
space, so each quantity has its chance level beside it.

The controls changed the reading before any block was classed. A registry element
lies in 54% of same-length windows in these neighbourhoods, so "a cCRE in the block"
says little: the registry counts only where it sits on the conserved bases, and
alone it adds 0.05 of confidence. The reader counts only when the element is open in
at least two cell types and twice the control mean. The parser recovers 104 of 760
known internal exons (13.7%) given the same flank, so its silence is not evidence;
and one coding-like segment among many is what the control rate (5.6% of control
segments) gives, so a coding reading needs a binomial tail of 0.05 or less, with
segments inside transposon remnants not counted.

| reading | blocks | kb | mean confidence | what carries it |
|---|---|---|---|---|
| regulatory element | 23 | 193.1 | 0.34 | a dELS, pELS or CTCF-only element on the conserved bases, or the reader alone; the reader above controls for 8, 14 of the 23 registry-only at 0.3 |
| promoter-like | 8 | 20.6 | 0.39 | PLS or DNase-H3K4me3 on the conserved bases, the reader above controls for 3; beside FAM237A, CDKN1A, DHH, CRCT1, FAM167B, HOXC13-AS |
| coding exon candidate | 5 | 21.5 | 0.36 | open frames above the median exon beyond chance, or the parser joining the block to a neighbour's exons |
| extended 3' end, unmeasured | 9 | 23.5 | 0.19 | conserved sequence within 3 kb of a coding gene's 3' end, no chromatin mark (PGRMC1, NR2E1, INSM2, CACNA1F) |
| unexplained | 24 | 127.8 | 0 | no layer supports a class; 15 of the 69 are closed in all eleven cells and carry no registry element |
| copy of a known protein, structured-RNA copy, primate-repeat artefact | 0 | 0 | | the checks ran and found none: 330 of 47,523 conserved bases lie in primate-specific repeats |

So of the 69: 31 plausible regulatory elements (0.2 to 0.55), 5 that look like
unannotated coding sequence and 9 that may be unannotated RNA (3' extensions), and
24 unexplained; mean confidence 0.21 over the set. The strongest single readings are
chr2:206,641,282 (promoter-like at FAM237A, open in four cells, and the one lentiMPRA
element that is active in K562, log2 1.04), chr12:49,094,807 (PLS next to DHH, open
in five cells), chr7:25,257,129 and chr22:38,571,213 (enhancer-like, open in 11 and
10 of 11 cells against control means of 1.2 and 4.3), and, on the coding side,
chr7:55,588,336 (the block the classifier had called long_orf: 9 of 11 conserved
segments open-framed with exon-like codon usage, binomial p 3e-10, and a seven-exon
parser structure of 909 bp joining exons of the lncRNAs VOPP1-DT and ENSG00000233977,
with no known human protein sharing its peptide) and chr8:43,271,960 (5 of 8
segments, p 3e-5, a 903 bp parser structure beside the POTEA node). Ground truth
reaches few blocks: lentiMPRA tested 7 (2 active, the other at TCF7L2's node in
HepG2), VISTA one (hs1103, negative, inside the 84 kb PBX3 block), ClinVar none, and
the GWAS and GTEx rows touch 2 and 1 blocks, but those two tables were distilled only
near scored elements, so their absence says nothing.

The set as a whole, against its control windows (the r-th window of every block is
replicate r, 40 replicates):

| layer | 69 candidates | control windows | p |
|---|---|---|---|
| a registry element anywhere in the block | 69.6% | 53.9% | 0.049 |
| cell types open on the block (mean of 11) | 2.78 | 2.33 | 0.15 |
| a parser exon on conserved sequence | 17.4% | 4.3% | 0.024 |
| conserved segments coding-like | 13.2% of 227 | 5.6% of 2,711 | |
| conserved bases inside a registry element | 7.6% | 19.8% of control conserved bases | |
| conserved bases open, per cell type | 2.2% | 3.7% of control conserved bases | |

This is the finding the table of readings should be read through. The candidates'
constrained bases are less often a registry element and less often open than
conserved bases a few tens of kilobases away, which is partly how they were chosen
(a block dense in registry elements is tiered regulatory, not constrained_unknown)
and partly what they are: sequence held on both axes that the chromatin assays of
these cells do not mark. The coding signal is the one that stands above its control,
more than twice the rate of coding-like segments and four times the parser exons, so
unannotated coding or transcribed sequence is a better-supported hypothesis for part
of this set than the regulatory count suggests, and the regulatory count is mostly
registry-only readings at 0.3. Nine blocks have every touched Gnocchi kilobase
holding a neighbour's exon bases (the human axis may be borrowed; their confidence is
lowered by 0.1). Of the 69 nodes, 49 have a coding TSS to infer a target from and 26
hold elements the deletion runs named a target for, but none of the scored elements
lies on a candidate block, so no candidate has a deletion reading of its own.

What is weak: codon log-odds partly measure GC content, and the coding-like peptides
of chr8:43.27 Mb are proline and arginine rich, the mark of GC-rich sequence; the
extended-3'-end reading is a rule of adjacency with no RNA behind it; the cell panel
has no neural, gonadal or embryonic tissue beyond H1 and SK-N-SH, where much
constrained non-coding sequence is active. The layers that would settle most of the
69, measured RNA over the blocks (ENCODE total RNA-seq, which the closure already
reads) and AlphaGenome's predicted tracks and deletion on each block, were not used
here, the first because this run kept to local data and the second because the
all-elements chain holds the request quota.

## What the whole-chromosome scoring bought the 98%, decided on chr21 (2026-09-13)

Albert's question before more quota goes to the sweep: chromosome 21's 12,139 elements
have all been deleted in AlphaGenome with the effect per gene and per cell line (5,409
requests), so what did that buy the UNKNOWN space, does it sharpen the syntax reading,
and does measurement back it? `attribution/unknown_scoring.py`
(`unknown_scoring_chr21.json`, 33 s, no model call) answers from the tables on disk.

2,334 elements (19.2%) overlap an UNKNOWN block, edges counted, and 1,155 of them move a
gene: 1,854 over 180 of the 194 regulatory blocks, 225 over 54 of 104 fossil, 203 over 44
of 101 neutral, 35 over 5 of 27 structural, and 17 over 5 of the 20 constrained_unknown
blocks. Of the 5,174 elements in the closure's input (those naming a coding gene), 665
lie over an unknown block.

**The ablation.** The closure was rebuilt from the committed gene by cell table and
reproduces it exactly (40.7% against 24.8% shuffled, p 0.001, rho 0.191), then rerun
without and with only those 665, each against fifty random removals of the same number
of elements per gene, so a fall cannot be the mere amount of input taken away.

| element input (the cell's own score, DNase-gated) | genes | most active cell is the most expressed | p | rho | rejected |
|---|---|---|---|---|---|
| all elements (the committed closure) | 141 | 40.7% | 0.001 | 0.191 | 73 |
| without the unknown-block elements | 140 | **40.5%** | 0.001 | 0.213 | 67 |
| random removals of the same size, 50 draws | 139.8 | 39.5% (36.2% to 41.4%) | | | |
| only the unknown-block elements | 61 | **28.7%** (shuffled 25.2%) | 0.26 | 0.072 | 35 |
| random subsets of the same size, 50 draws | 53.3 | 35.8% (23.8% to 49.1%) | | | |

Dropping the elements over unknown space does not move the one test that works: the
closure keeps 40.5%, above 86% of random removals of the same size, the correlation
rises and six rejections disappear. On their own those elements do not pick the cell (p
0.26) and do worse than 88% of random subsets of the same size, though that last gap is
inside the draws' spread. The model's reading over unknown space adds nothing to the
closure on this chromosome.

**Are they different?** Against the other 9,805 elements, raw and inside 480 strata of
length, GC and distance to the nearest coding TSS (40 distance bins, because eight left
the unknown elements at 260 kb against 173 kb; with 40 the matched means are 254 and 246
kb), with the flag permuted inside strata 2,000 times:

| per element | over unknown blocks | the rest, matched | p |
|---|---|---|---|
| moves a gene | 49.5% | 68.0% | 0.0005 |
| names a coding gene | 28.6% | 43.1% | 0.0005 |
| abs log2 of the best target, movers | 0.309 | 0.338 | 0.048 |
| silencer-like (expression up on deletion), movers | 37.7% | 38.8% | 0.49 |
| cell lines at 0.1 log2 or more, of four, movers | 0.45 | 0.80 | 0.0005 |
| distance to the coding target's TSS, kb | 45.0 | 49.3 | 0.11 |

They move genes less often, name a coding gene less often and act in about half as many
cell lines, with a marginally smaller effect; the silencer share and the reach to the
target do not differ. One confound no stratum removes: 99.9% of the other elements sit
inside a gene body (introns) against 3.9% of these, so the comparison is intergenic
against intragenic, which is what the unknown space is by construction.

**The syntax reading.** Chromosome 21 has no block constrained on both axes at any
tier (the organiser found none), so the test the question asks for has no members here.
The nearest one, mammal-constrained blocks (relaxed, 5 blocks, 23 elements) against the
rest: moving 47.8% against 47.5% matched (p 0.96), abs log2 0.315 against 0.281 (p
0.74), cell lines 0.91 against 0.35 (p 0.02 per element, but 0.28 with blocks as the
unit, which is the honest null since elements share their block's label). Nothing
survives. Of the 69 genome-wide syntax candidates none lies on chr21; two lie on
chromosomes whose scoring is complete. chr19:53,943,950 holds one scored element that
moves nothing. chr22:38,571,213 holds one that names DMC1 at 0.108 log2 (psoas muscle),
with no cell line reaching the bar. So the sweep gave the candidates one weak target
and no cell.

**Against measurement.** lentiMPRA covers 329 scored elements in K562 and HepG2 (the
library was drawn from those lines' DNase peaks): the signed deletion effect correlates
with measured activity at rho 0.22 in K562 and 0.24 in HepG2 (both p 0.0005 by
permutation), the unsigned magnitude at 0.10 (p 0.057) and 0.17 (p 0.002). Only 59 of
the 329 lie over unknown blocks (rho 0.01 in K562, 0.25 in HepG2), too few to say the two
groups differ. VISTA has 19 elements on chr21: 10 of 13 positives against 2 of 6
negatives hold an element that moves a gene, reported and not read. ClinVar and the GWAS
table are not used: the GWAS rows were distilled only near the sampled elements.

**Coverage and cost.** Before the sweep the two samples named a coding target in 11 of
the 446 unknown blocks; after it, 140, and 90 of those have a cell line where the
deletion reaches 0.1 log2. By tier: regulatory 9 to 117 of 194, neutral 2 to 15 of 101,
fossil 0 to 7 of 104, constrained_unknown 0 to 1 of 20 (no cell), structural 0 of 27.
The blocks so touched span 3.9 Mb of the 20.8 Mb of unknown space, but the elements
naming a target cover 185 kb of it, 0.9%, against 3.1 kb before. The runs cost 0.45
requests per element on chr21, 0.74 on chr19 (21,590 requests, 603 quota waits) and 1.0
on chr22, 0.76 over the complete runs. With 83,128 elements scored so far and about
884,000 in a node left in the registry's 1.06 million, finishing the sweep is about
674,000 requests, some 125 chromosome 21s. Scoring the constrained unknown directly is
of another order: the 882 blocks that are not copies and the 69 candidates are about a
thousand deletions (a few thousand if the long blocks are cut into their conserved
segments), under a fifth of chr19.

**Verdict.** For the 98% the sweep does not pay. Its elements over unknown space are
removable from the closure without loss, carry no cell on their own, act in fewer cell
lines than intragenic elements, barely reach the constrained unknown (17 elements, one
block named, no cell), and gave the syntax candidates one weak target. What it did buy
is bookkeeping: 129 more unknown blocks, mostly regulatory tier, with a named target that
the one falsifying test does not use. The model itself is not the problem, since its
signed effect tracks measured activity (rho 0.22 to 0.24): the registry's elements are
simply not where the unknown space's constraint is. Recommendation: stop the sweep as an
instrument for area I after chr20, and spend the next quota scoring the 882 non-copy
constrained-unknown blocks and the 69 candidates directly, their conserved segments
deleted with matched control windows, which is what would test the coding and regulatory
readings of the candidate section. The closure's replication on a third chromosome needs
no quota, since chr19 is complete. What is weak: one chromosome, and a small one with no
syntax block; the cell panel of four lines; and an element counted as over unknown space
when any base of it touches a block.

## The lexicon: the genome's units at every scale, with their nulls (2026-09-13)

Albert's brief: treat the genome as a maths problem. Take its units from small to big, put
them in one index with statistics, and make every count carry its expectation. Then ask two
questions of it. Is the fossil tier a library for its node? Does the index separate syntax
from value slots? Code: `attribution/lexicon.py`, run by `scripts/lexicon.py chr21 chr22`,
with no model call and nothing streamed. The index itself is
`data/knowledge/lexicon/lexicon_<chrom>.json.gz` (git-ignored). The committed summaries are
`data/results/lexicon_chr21.json` and `lexicon_chr22.json`.

**Schema.** One index per chromosome, so the genome is added one chromosome at a time.

- **Levels.**
  1. `kmers`: every 6-mer counted in each of 13 contexts, against an order-2 Markov chain
     fitted inside the same context.
  2. `seeds`: the data-driven vocabulary. These are the 2,000 most frequent exact 16-mers
     occurring at least 50 times, extended by consensus and joined into families when their
     occurrences co-occur at a fixed offset.
  3. Curated units: RepeatMasker subfamilies and families (simple repeats and low complexity
     left out), cCRE classes, GENCODE CDS and exons, segmental duplications, and the JASPAR
     sites the motif runs recorded.
  4. UNKNOWN blocks with class, budget tier and human-axis case.
  5. CTCF nodes.
  6. Libraries with at least five coding genes on the chromosome.
- **The 13 contexts.** Each base gets one context, painted in rising priority: tier of the
  UNKNOWN block, gene body, cCRE class, exon, CDS. The contexts are cds, utr_or_exon,
  exon_noncoding, promoter_like, enhancer_like, ctcf_only, intron, one per budget tier
  (unknown_structural, unknown_fossil, unknown_regulatory, unknown_constrained,
  unknown_neutral) and other.
- **What every unit records.**
  - Occurrences (sampled to 400 for a seed, by one pseudo-random key per 256-base bin shared
    by all seeds, so the seeds of one segment keep the same loci) and the context of each.
  - Composition: the repeat share.
  - Conservation: the share of its bases inside a 100-vertebrate conserved element.
  - The human axis: the share of its measured bases in blocks whose Gnocchi kilobases pass
    variation.py's bar, and the share measured at all.
  - The tiers and cases of the blocks it falls in, and its overlap with the 69 candidates.
  - A context-enrichment table.
- **The background model.** Windows of 10 kb carry a GC bin (5) and a repeat-share bin (4),
  giving 20 strata.
  - A unit's expected count in a context is its own density inside each stratum spread over
    that context's bases in the stratum.
  - Its expected conservation is the conserved share of the same context inside the same
    strata.
  - Its expected human share is the share over measured bases of the same strata.
  - Tails are Poisson and are taken on occurrence-equivalents (bases divided by the mean
    occurrence length), because the bases of one occurrence are not independent draws.
  - Every family of tests is corrected by Benjamini-Hochberg at 5% and its size is reported.
  - A unit that paints a context (a cCRE class, the CDS) is not tested in it.
  - A JASPAR site is also read against the rest of its own element, because the motif runs
    only record sites inside the constrained and enhancer targets.

**Two errors caught before any number was read.**
1. The first Poisson tail summed the lower tail from exp(-expected). That underflows to zero
   above an expectation of about 745, so every large count came back with p = 1.0.
2. The first seed sample was independent per seed. That broke the co-occurrence the families
   are found by: families went from 29 to 288 with it and back to 32 with the shared key.

Both are now tests.

**The numbers, chr21 then chr22.**

- **Size.** 3,143 and 3,098 units: 1,022 and 962 subfamilies, 46 families each, 67 and 82
  JASPAR factors, 2,000 seeds each. 446 and 635 blocks, 228 and 247 nodes, 16 and 33
  libraries.
- **Context enrichment.** 22,076 and 22,795 tests; 2,172 and 2,279 pass at a threshold of
  0.0049. The strongest with at least ten occurrences are known biology or known structure:
  - U6 snRNA copies in non-coding exons (18.6x and 17.6x);
  - the (CATTC)n satellite in the structural tier (17x);
  - (GA)n in promoter-like elements (23.8x on chr21);
  - L1MB2 in chr21's constrained unknown (29 against 1.24);
  - JASPAR sites in enhancer-like elements (8x to 10x). Those sites were recorded only
    inside elements, so that enrichment is their selection, not a finding.
- **k-mers.** 52,411 and 52,568 tests over the contexts; 10,377 and 10,096 pass BH, and
  1,168 and 1,317 are also twofold over the Markov expectation. The top five per context
  contain CG in 60 of 65 words on chr21 and 49 of 65 on chr22, in every context. An order-2
  chain fitted over a whole context cannot absorb CpG islands, which cluster, so at this
  depth the k-mer level records that heterogeneity and no vocabulary. Poisson on counts of
  thousands is overdispersed besides. **Negative as a lexicon; kept as a level.**
- **Seeds.** The vocabulary found from sequence is the repeat library.
  - 1,965 of chr21's 2,000 seeds and 1,980 of chr22's lie at least half in RepeatMasker, and
    every other one is a tandem word of period six or less. On neither chromosome is there a
    recurring 16-mer that is neither.
  - The largest family, 1,360 and 1,004 seeds, is the alpha satellite (span capped at 135 by
    the extension reach). The next ones are Alu and L1 pieces.
  - **Negative:** at 50 occurrences and the top 2,000 words, no unannotated recurring
    segment exists on these two chromosomes.

**Question 1, the fossil tier as a library for its node.**

- **Within-node similarity.** 104 fossil-tier blocks with at least 500 repeat bases on
  chr21 (19 nodes holding two or more) and 106 on chr22 (20 nodes).
  - Within a node the mean cosine of subfamily composition is 0.089 against 0.071 between
    nodes (chr21), and 0.178 against 0.100 (chr22).
  - Against a label permutation inside repeat-base strata: null 0.069 and 0.092, p 0.044
    and 0.001, Fisher-combined p 0.0005.
  - Against the same partition shifted along the chromosome, which keeps node sizes and
    contiguity: null 0.088 and 0.142, p 0.47 and 0.11, combined p 0.21.
  - Fossils of one node are alike because neighbours are alike. The CTCF partition adds
    nothing to proximity. **Negative.**
- **Link to the libraries of the node's genes.**
  - 772 and 995 library by subfamily tests; 4 and 5 pass BH.
  - Random gene sets of the same sizes pass 7 to 16 (chr21) and 8 to 23 (chr22) in ten
    draws each.
  - The real libraries pass fewer pairs than random gene sets. **Negative.**

**Question 2, syntax against value slots.** A unit is syntax when it is held across species
(conserved share at least 0.10) and constrained among people (human share at least 0.25).
It is a value slot when held but free among people, or free on both axes, while present in
at least four contexts with 100 occurrences and conserved above its neighbourhood.

- **Classes.** chr21 then chr22: 2,445 and 2,395 units tested; syntax 21 and 22, relaxed 16
  and 8, recent 271 and 443, tolerant 764 and 1,277, unmeasured 1,373 and 645 (the
  satellite seeds of the structural tier, which Gnocchi does not cover).
- **Syntax against the nulls.**
  - 12 and 12 survive the context-matched conservation null (26 and 23 pass it at all, of
    2,445 and 2,395 tests).
  - 10 and 10 of those are JASPAR factors that then fail against their own elements: 0 of
    22 factors pass on either chromosome.
  - Pooled, the sites are 0.335 conserved against 0.308 for the rest of their elements on
    chr21 (p 0.050) and 0.336 against 0.310 on chr22 (p 0.045). Overlapping sites make both
    p optimistic.
  - What survives every null is two units per chromosome, none of them new: the
    promoter-like cCREs (0.16 against 0.09 and 0.17 against 0.11), an LSU rRNA copy on
    chr21's acrocentric arm (0.58 against 0.02), and chr22's (GATG)n tandem in
    regulatory-tier blocks (0.11 against 0.025).
- **Value slots.** 1 on chr21 (the hAT family, 0.126 conserved against 0.033 expected) and
  0 on chr22.
- **By tier.** Almost all syntax units sit in the regulatory tier (20 of 21, 22 of 22).
  Their human case follows the case of their blocks, because the human axis is block-level.
- **Candidates.** chr22's one candidate holds Alu, L1 and L2 pieces, all tolerant or recent.
- **Verdict.** The index does not separate syntax from value slots at unit scale. It finds
  the syntax already known (promoters), shows that a motif site's apparent conservation is
  its element's, and has no resolution on the human axis. **Negative.** It is the resolution
  of the axes that fails, not the idea. The human axis is a kilobase score attached to
  blocks. Conservation is element membership, not Zoonomia's per-base phyloP; the brief
  asked for phyloP and the index does not yet read it.

**Cost of the genome-wide index.**
- **Time and memory.** The chr21 and chr22 run took 135 s on one core with a peak resident
  set of 1.03 GB: 1.38 to 1.39 s per Mb including both questions. That puts the 3.1 Gb
  genome at about 71 minutes one chromosome at a time. Peak memory grows with the largest
  chromosome, about 5 GB for chr1 by linear extrapolation.
- **Disk.** The index is 8.9 to 9.4 kB per Mb gzipped, about 29 MB genome-wide under
  data/knowledge. The committed summaries are 34 to 37 kB per chromosome, under 1 MB for
  the 24.
- **What would give the axes their resolution, not spent here:**
  - per-base Zoonomia phyloP streamed from the bigWig by range requests (a local one-bit
    constrained flag would cost about 390 MB uncompressed genome-wide);
  - per-kilobase Gnocchi in place of the block summary;
  - a genome-wide JASPAR scan at 0.95 (motifs.py's calibrated threshold) so sites stop
    being a sample of selected elements.
- **Per chromosome, not pooled.** Each chromosome's seed vocabulary is its own top 2,000,
  so a genome-wide vocabulary needs the seed counts merged across chromosomes. That is
  another pass of the same cost.

**What is weak.** Two small chromosomes. Conservation by element membership rather than
per-base phyloP. A human axis that is a property of blocks rather than units. Seeds capped
at the top 2,000 per chromosome. JASPAR sites that are a sample of selected elements.
Poisson tails that remain optimistic wherever occurrences cluster, which is why the library
link carries its random-gene-set calibration.

## What comes next, in order

1. Done 2026-09-13: the whole-input closure passing on chromosomes 21 and 22,
   the rejected set read, and four corrections of the input tested and
   rejected against a promoter-only control (the closure section above).
2. Done 2026-09-13: **the 69 syntax candidates read one by one** (the section
   above): 31 plausible regulatory, 5 coding candidates, 9 possible 3'
   extensions, 24 unexplained. What remains of the step: measured RNA over the
   69 blocks and their controls, to test the coding and 3'-extension readings,
   and AlphaGenome's predicted tracks and deletion per block once the chain frees
   the quota. The two axes disagreeing is itself a
   flag: chr15:84,395,903-84,398,315 is 23% constrained across mammals and
   carries far more human variation than expected (Gnocchi Z of -8.5 over its
   kilobases), which reads as a mutation hotspot or a mapping artefact rather
   than a value slot, and the constrained_unknown tier holds 18 such blocks on
   chr15 to look at first. The cause is now measured: curated segmental
   duplications cover 60.6% of chromosome 21's constrained_unknown tier (15 of
   20 blocks mostly duplicated; the most constrained block of the chromosome is
   100% duplicated with 8 partners) against 7.8% of the regulatory tier, so on
   that chromosome most of "deep constraint with no human score" is a copy
   whose alignment is paralogous and whose variants do not map. Genome-wide the
   tier is 4.7% duplicated by bases, the chr21 figure being its acrocentric
   short arm, but 216 of its 1,098 blocks, one in five, are mostly copies. The
   organiser reads a block's duplicated fraction first (`duplication_<chrom>.json`,
   area J) and attributes second; the real unknown is what remains after the
   copies.
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
