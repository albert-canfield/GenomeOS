# The 98%: attributing function to the non-coding genome

Area I of ROADMAP.md. Code: `genomeos/attribution/` (`bigwig.py`, `constraint.py`,
`budget.py`, `organise.py`, `candidates.py`), `scripts/budget_genome_wide.py`,
`scripts/syntax_candidates.py`, `genomeos budget`. Results:
`data/results/budget_<chrom>.json`, `budget_genome_wide.json`,
`syntax_candidates_genome_wide.json`.

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
