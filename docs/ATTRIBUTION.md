# The 98%: attributing function to the non-coding genome

Area I of ROADMAP.md. Code: `genomeos/attribution/` (`bigwig.py`, `constraint.py`,
`budget.py`, `organise.py`, `candidates.py`, `unknown_scoring.py`, `lexicon.py`, `human_panel.py`),
`scripts/budget_genome_wide.py`, `scripts/syntax_candidates.py`, `scripts/lexicon.py`, `scripts/human_panel.py`,
`genomeos budget`. Results: `data/results/budget_<chrom>.json`, `budget_genome_wide.json`,
`syntax_candidates_genome_wide.json`, `unknown_scoring_chr21.json`, `lexicon_<chrom>.json`, `human_panel_chr21.json`,
`element_types.json` (`element_types.py`, `scripts/element_types.py`).

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

> **Annotated 2026-09-27**, figures unchanged. Every number in this section is read from the compact
> one-gene-per-element table. Two of them were computed over structural zeros and one set of six was
> conditioned on the sweep's `MIN_EFFECT` gate without saying so; the section below,
> "the unknown scoring off the one-gene table", says which and reports the corrected figures beside
> them. The naming rates here are not among them: the compact table's coding target was already
> chosen from the whole window, and the re-run reproduces 5,174 of 12,139 exactly. What this section
> lacked, and could not have had, is a matched-random control on its own instrument: the unknown
> blocks name a coding target at 49.1% against 65.3% at length-matched random windows on chr21.

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

## The lexicon at base resolution: Question 2 re-asked, chr21 and chr22 (2026-09-14)

The first lexicon named why its negatives might be soft. The sharper instrument is
`attribution/lexicon_axes.py`, run by `scripts/lexicon_axes.py chr21 chr22` with no model
call. It writes `lexicon_axes_chr21` and `lexicon_axes_chr22`. The byte tracks and the full
unit rows go under `data/knowledge/lexicon`. Each soft axis was replaced.

- **Mammals.** Zoonomia phyloP over 241 mammals, every base.
  - The track is streamed once through `attribution/bigwig.py`, the reader the budget uses:
    97 and 95 MB in coalesced ranges, 4 minutes per chromosome. It is kept as one signed byte
    per base, with 2.27 exactly on a bin edge (23 MB gzipped).
  - 80% and 72% of bases are measured, and 2.4% and 3.7% of the measured bases are
    constrained.
- **People.** Which measure to trust was part of the question.
  - genomeos-h1 showed that Gnocchi and the 89-assembly panel disagree (-0.08 per kilobase)
    and that per-kilobase panel counts are 15.6 times overdispersed.
  - A unit a few bases long can only be read by a per-base measure. The axis is therefore
    the panel's recurring substitutions: a column where at least two assemblies carry
    something other than the commonest state and one of them is a base. Only bases inside
    blocks with at least 80 of 89 assemblies present count (34 million informative bases
    per chromosome, 0.57% and 0.60% variable). These are read from h1's stores through
    h1's `Panel`.
  - The expectation is the base's strand-collapsed trinucleotide by context by GC bin, so
    CpG mutability is inside it, and again with the base's phyloP bin added.
  - The overdispersion is met where it arises. Every test's variance is built from the
    unit's own occurrences (a sandwich estimate, floored at Poisson). Every family is also
    read on the same units shifted as one piece around the chromosome.
  - Gnocchi is joined per base for comparison only.
- **Motifs.** Every JASPAR 2026 CORE vertebrate profile was scanned over the whole
  chromosome at 0.95, both strands.
  - The scan is indexed: one sort of the 8-mers, then feasible cores as ranges. It checks
    identical to `motifs.all_hits` and takes 20 s.
  - It gives 13.5 and 11.6 million merged sites in 895 and 894 distinct site sets, 21
    names being aliases of another's sites.
  - Every profile's columns were also permuted and scanned the same way. A decoy keeps its
    motif's letters and information and loses its order, so an effect a decoy shares
    belongs to the letters.

**The tests.** A row is a unit in one context with at least 20 occurrences. Units are 1,125
curated units (the first lexicon's, minus its 0.85 sites), about 890 motif site sets, about
880 decoy sets and 2,000 seeds. That gives 10,831 and 10,587 real rows.

1. Mammals against the same context in the same GC by repeat stratum: 9,268 and 9,427
   tests, 1,152 and 1,267 passing.
2. People against trinucleotide by context by GC: 5,445 and 6,597 tests, 343 and 267
   passing.
3. People against that expectation with phyloP added: 5,445 and 6,597 tests, 421 and 529
   passing.
4. Each axis against the occurrences' own flanks (as long as the occurrence, 10 bases
   away):
   - mammals: 7,490 and 7,771 tests, 649 and 555 passing;
   - people: 5,567 and 6,370 tests, 263 and 302 passing.

The shifted copies read with the real thresholds pass mammals 94 and 160 times and people
52 and 34 times, but the flank tests only 67 and 109 times (mammals) and 35 and 37 (people).
The flank tests are the ones to read.

**Motif against its own decoy.** Each factor is compared with its decoy in the same context,
both as site over flank.

- **Mammals.** Motifs win in every context on both chromosomes: 549 of 743 and 571 of 749
  factors, median difference +0.145 and +0.153, largest in enhancer-like elements (+0.275
  and +0.252).
- **Against the rest of their own cCRE.** Motif sites pass for 334 of 690 and 405 of 723
  factors; decoys for 37 of 679 and 69 of 711. At 0.95 and at base resolution, mammals do
  see motif sites, beyond their letters, their flanks and their element.
- **People.** The same comparison gives a median difference of -0.018 (315 of 577, sign
  p 0.03) and -0.024 (318 of 572, p 0.008), with no context consistent across the two
  chromosomes.
  - Enhancer-like: -0.041 then -0.021 (n.s.).
  - Fossil: -0.050 then +0.038.
  - Unknown-regulatory: -0.013 (n.s.) then -0.041.
  - Against their own cCRE, motif sites carry 0.986 and 0.975 of the variation the rest of
    the element predicts; decoys 0.999 and 1.006. Factors passing: 0 of 555 on chr21 and
    17 of 599 on chr22, against 5 and 1 decoys.
  - Factors with similar matrices share sites, so every sign-test p here is optimistic.

**Question 2, re-asked.** Syntax needs a row held across mammals and among people in the same
context.

- Rows held against both flank nulls: 57 and 24 motif rows against 23 and 11 decoy rows, a
  surplus that comes from the mammal side.
- Among rows held across mammals against flanks, the share also held among people is 9.0%
  for motifs against 9.2% for decoys on chr21, and 4.5% against 7.9% on chr22.
- Against the strata, decoys pass the people test 192 times where motifs pass 233 (chr21),
  and 189 where motifs pass 193 (chr22).
- The rows that look like syntax are AT-rich homeobox, NFAT and STAT sites in introns. Their
  people depletion is shared by decoys of the same letters.
- Value slots (held across mammals, more variable among people than their flanks): 3 motif
  rows against 7 decoys on chr21, 0 against 1 on chr22.
- The only row that replicates is the dELS class. It is more conserved than its flanks
  (0.043 against 0.029, z 12.5; 0.052 against 0.040, z 12.7) and more variable among people
  than its flanks (O/E 1.01 against 0.93, z 6.7; 1.00 against 0.91, z 7.9). That is what a
  value slot would look like. It is also what open chromatin's mutation rate would look
  like, and diversity alone cannot tell the two apart.
- **Verdict. Negative, and stronger than the first.** At base resolution the mammal axis
  separates motif sites from their letters in every context. The people axis at 89
  assemblies sees at most a 2 to 3% depletion and adds nothing the mammal axis does not
  already carry, so units do not separate into syntax and value slots. The panel is not
  blind: canonical CDS reads 0.49 and 0.46 of expected variation. It is that motif-site
  constraint among people, if present, is an order of magnitude weaker than coding
  constraint, below what 89 haplotypes resolve per unit.

**Which human measure.** The panel, for the reason above. Gnocchi does not agree with it at
unit scale either: rank correlation with panel depletion is -0.067 on chr21 and 0.013 on
chr22, and with phyloP -0.128 and 0.100. Nothing in the re-ask depends on Gnocchi.

**The fossil tier on the sharper axes.**

- **Tier-wide.** The fossil tier is 0.97% and 1.03% constrained, against 1.5% and 1.7% for
  intron. The tier was selected for low block phyloP, so that comparison is circular.
- **People, tier-wide.** Observed over the trinucleotide by GC expectation is 1.145 and
  1.160, against intron 1.019 and 0.945 and CDS 0.491 and 0.458.
- **People, matched.** Matched subfamily by subfamily and within replication-timing
  terciles (h1's Repli-seq per kilobase), the tier's copies against the same subfamily's
  copies elsewhere:
  - against intronic copies: 0.994 (38 of 62 subfamilies more variable, z -0.22) and 1.086
    (27 of 34, z 3.0);
  - against unknown-regulatory copies: 0.946 and 1.001;
  - against neutral-tier copies: 1.111 on chr21 (chr22 has only 2 subfamilies to pair).
- The unmatched comparison (fossil copies more variable in 54 of 67 subfamilies) was
  replication timing, and matching removes it.
- **Motif sites inside fossils.** The chr21 hint (people -0.05, p 0.0008) reverses on chr22
  (+0.038).
- **Negative.** The fossil tier is no less variable among people than its own subfamilies
  elsewhere. This is the fourth independent failure on that tier: the node libraries, the
  methylation (genomeos-e1), h1's block variability, and now per-base diversity with timing
  matched. The tier is what its name says.

**The syntax candidate on these chromosomes.** chr22:38,571,213-38,574,930 is 5.5%
constrained per base and carries 4 recurring substitutions against 19.0 expected (23.2 with
phyloP). That is one block, against h1's 59 of 69 candidates variable.

**Cost, and why the other chromosomes are not built.**
- **Run time.** Both chromosomes took 400 s in one process, and memory peaked at 6.9 GB.
  Most of the memory is per-base arrays, which scale with length, so chr1 would need about
  18 GB without chunking.
- **Streaming.** phyloP costs 97 MB and 4 minutes per chromosome to stream, and 23 MB to
  keep: about 1.5 GB kept for the genome. The panel stores are h1's to build (5.3 hours of
  streaming for the genome).
- **Decision.** Nothing separates, so none of that is spent.

**What is weak.**
- **Human depth.** 89 haplotypes, recurring substitutions only (indels and structural
  alleles left out).
- **Motif calls.** Sites are 0.95 predictions, not bound sites.
- **Flanks.** The flank null assumes the flanks share the site's mutation environment.
- **Sign tests.** Motif-against-decoy sign tests count factors that share sites.
- **Replication timing.** It is matched for the fossil comparison only.
- **Next.** Measurement is what could still sort a site into syntax or slot: saturation
  mutagenesis and MPRA allele effects at motif sites, or per-base gnomAD frequencies for
  population depth.

## Many human genomes: fixed, storage and cannot place, chr21 (2026-09-14)

Albert's question: compare as many human genomes as can be read, and catalogue what is common
(the syntax), what varies among a few recurring values (the storage, "like a multiple-select
db table column") and what cannot be placed, the unknown space above all and the known parts as
the control. Until now the human axis was Gnocchi, a Z per kilobase from 76,156 genomes, SNVs
only; the only genomes read were the GIAB trio. `attribution/human_panel.py`
(`scripts/human_panel.py`, `human_panel_chr21.json`, no model call) reads the genomes.

**What is read, and how.**
- **The panel.** HPRC release 1 aligned to hg38 by Cactus (UCSC `hprc90way`): hg38, T2T-CHM13
  and 88 haplotype assemblies of 44 people, so 89 genomes beside the reference. The alignment
  is uncompressed MAF, 3.28 GB for chr21 and 266.6 GB for the genome. UCSC's REST API returns
  the byte offset of every block; the blocks are streamed in 32 MB ranges and distilled on the
  spot into three things kept under `data/knowledge/human_panel/chr21` (git-ignored, 8.4 MB):
  each block's span with the status of every assembly (aligned, deleted with contiguous
  flanks, replaced by other sequence, missing), every column where some assembly differs from
  hg38 with each assembly's allele, and the bases each assembly inserts. chr21 cost 99
  requests and 236 s: 21,495 blocks, 1,054,284 variable columns, 424,828 variants once
  deletion runs and insertions are collapsed into events.
- **Structure.** The HPRC arrangement tracks (`hprcArrV1`: duplications, inversions,
  deletions, insertions, double rearrangements, with the number of genomes carrying each) and
  the release-2 v2.1 structural variants from 233 assemblies (`hprc2v21Sv`, AC over AN), read
  as bigBed by range. The classes use the alignment's own events; these tracks are annotation
  (`structure_by_tier`).
- **Frequencies beyond 89 haplotypes.** gnomAD v4.1.1 genomes for substitutions and short
  indels, gnomAD v4.1 structural variants, and the TRExplorer catalogue's tandem-repeat
  allele-size histograms (TenK10K, typically 3,850 alleles), all bigBed by range. Each carries
  its own sample size and none is merged with a panel count: every event reports the panel's
  minor count out of 89 and, beside it, gnomAD's AF with its AN.
- **Replication timing** (genomeos-9c's control: late-replicating DNA mutates more, so
  "variable" may be timing rather than selection). ENCODE UW Repli-seq wavelet-smoothed
  signal for 11 lines (BG02ES, BJ, GM12878, HeLa-S3, HepG2, HUVEC, IMR-90, K562, MCF-7, NHEK,
  SK-N-SH) exists only on hg19, so each hg38 kilobase's centre is lifted through UCSC's
  `hg38ToHg19` chain and the 11 signals are averaged: 35,916 kilobases, 545 requests, 6.5 MB.
  Tertiles of the chromosome's values are late, middle and early.
- **Read-out evidence** for the later executor test: GTEx v8 fine-mapped eQTLs (DAP-G) and
  MPRAVarDB allele pairs, both bigBed.

**The two readings.** A *recurring* variant is one where two or more of the 89 assemblies
depart from the commonest state; singletons (assembly errors among them) are left out.
- **Blocks.** Presence, identity, π, the share of bases touched by a recurring deletion,
  replacement or insertion, and the count of recurring variants against an expectation from
  the chromosome's own rate outside coding exons, matched per kilobase on GC (seven strata)
  and timing (tertiles). Classes: *unplaced* (under half aligned, a quarter of assembly-bases
  missing, or fewer than 5 expected variants), *lineage-restricted* (presence under 0.5),
  *polymorphic* (presence under 0.95, or a tenth of the bases under a recurring structural
  event), *core* (half the expected rate or less, Poisson lower tail 0.01 or less) and
  *variable*. Every block is also compared with two control windows of its own length
  elsewhere on the chromosome, GC within 0.02 and the same timing tertile.
- **Units.** Every block, every gene's canonical CDS and the whole chromosome (the background)
  tiled into 200-bp units. Each assembly's value is its combination of states at the unit's
  recurring variants (base, gap, inserted length, absence). *Fixed*: one value in 95% of
  assemblies. *Storage*: recurring values cover 80% of the assemblies and there are at most
  8 of them. *Hypervariable*: otherwise, a value per haplotype rather than a column.
  *Unplaced*: under half aligned or under 80% of the assemblies informative. The same rule
  applied to the length allele alone (net bases gained or lost) reads copy numbers.

**The controls, first.** Coding exons against the tiers, pooled; ratio is recurring variants
over the matched expectation:

| set | presence | identity | π | recurring per kb | ratio, GC matched | ratio, GC and timing |
|---|---|---|---|---|---|---|
| canonical CDS (221 genes, 318 kb) | 0.9991 | 0.99930 | 0.00073 | 3.77 | 0.351 | **0.411** |
| neutral tier | 0.9978 | 0.99874 | 0.00137 | 7.20 | 0.958 | **0.907** |
| fossil tier | 0.9976 | 0.99842 | 0.00155 | 7.77 | 1.059 | 0.988 |
| regulatory tier | 0.9946 | 0.99855 | 0.00161 | 8.94 | 1.030 | 1.087 |
| constrained_unknown | 0.9972 | 0.99846 | 0.00242 | 7.99 | 1.029 | 0.649 |
| structural (8.5% aligned) | 0.6248 | 0.98792 | 0.00259 | 12.91 | 1.578 | 1.064 |

| blocks | core | variable | polymorphic | lineage-restricted | unplaced | core of callable (GC only) | matched windows |
|---|---|---|---|---|---|---|---|
| CDS, per gene | 89 | 64 | 5 | 0 | 63 | **58.2%** (66.0%) | 14.3% |
| neutral | 6 | 70 | 3 | 0 | 22 | 7.9% (9.2%) | 3.7% |
| fossil | 1 | 81 | 2 | 0 | 20 | 1.2% (2.4%) | 4.2% |
| regulatory | 6 | 166 | 11 | 0 | 11 | 3.5% (9.3%) | 6.9% |
| constrained_unknown | 4 | 6 | 0 | 0 | 10 | 40% (30%) | 10% |
| structural | 0 | 2 | 2 | 1 | 22 | | |

| 200-bp units, placed | fixed | storage | hypervariable | fixed expected, GC | fixed expected, GC and timing |
|---|---|---|---|---|---|
| CDS | **70.3%** | 29.5% | 0.1% | 43.5% | 44.1% |
| neutral | 48.2% | 50.4% | 1.4% | 49.5% | 47.6% |
| fossil | 44.4% | 54.4% | 1.2% | 48.7% | 47.4% |
| regulatory | 41.7% | 56.5% | 1.9% | 46.6% | 46.8% |
| constrained_unknown | 56.8% | 42.5% | 0.8% | 48.3% | 52.9% |
| background (176,885 units) | 48.0% | 50.5% | 1.5% | | |

The instrument passes the ordering the brief made its condition: coding exons carry 0.41 of the
matched rate of recurring variation and the neutral tier 0.91, coding units are fixed 26 points
above their matched expectation and neutral units sit on theirs, and callable genes are core
58% of the time against 8% of neutral blocks and 14% of control windows of the same length, GC
and timing. It does not pass "overwhelmingly core": 64 of 153 callable genes are variable, and
30% of 200-bp coding units hold a common alternative, because synonymous and common coding
variation are real (the checks, in `control_verdict`, were written after the first read of
these numbers and describe the ordering seen). Only the coding exons stand clearly off the
background. No tier of the unknown space does: regulatory and fossil units are slightly less
fixed than matched windows, the neutral tier matches them, and the constrained_unknown tier's
four core blocks include the two copies at 13.76 Mb (100% segmental duplication, where the
alignment is paralogous) among ten callable blocks.

**Replication timing, the control.** Within a GC stratum the background rate is up to 1.5
times higher in late than in early DNA (GC 35-40%: 8.4 against 5.9 recurring variants per
kb; GC 45-50%: 11.1 against 7.2; GC 50-55%: 11.7 against 8.2; near 1 in the most GC-poor and
GC-rich strata). How much of each difference
survives matching on timing as well as GC:
- CDS against neutral, pooled ratio: a gap of 0.61 before (0.35 against 0.96), 0.50 after
  (0.41 against 0.91). **82% survives.**
- CDS units' excess of fixed over matched: 26.8 points before, 26.2 after. **98% survives.**
- Genes called core: 66% before, 58% after. The regulatory tier's core blocks fall from 9.3%
  to 3.5% of callable, so most of its apparent depletion was early-replicating sequence;
  constrained_unknown rises from 30% to 40% because its blocks are late.
- Across the chromosome's 176,885 background units the classes are spread evenly over
  timing: fixed 33% late, 34% middle, 33% early; storage 34%, 34%, 33%; hypervariable 32%,
  32%, 35%. Variable units are not enriched in late DNA at this scale. What is late is what
  cannot be placed: 74% of unplaced units.

So fixed against variable does not separate only before the control; the coding signal is
selection, not timing. The same control shows that most of the unknown space's small
departures from the background were timing or GC.

**The three catalogues** (200-bp units; the unknown space is 20.8 Mb):

| catalogue | unknown space: units | Mb | coding exons: units | kb |
|---|---|---|---|---|
| fixed (one value in 95%) | 24,272 | 4.85 | 1,143 | 227 |
| storage (a small value domain) | 26,969 | 5.40 | 480 | 97 |
| cannot place: hypervariable | 774 | 0.15 | 2 | 0.5 |
| cannot place: not aligned or missing | 21,281 | 4.26 | 3 | 0.6 |

- **Value domains.** Of the 26,969 storage units, 11,673 hold two recurring values, 9,024
  three, 3,642 four and 2,629 five to eight. The median effective number of values is 1.68,
  the median commonest value holds 73% of the assemblies. 20,236 are substitutions only; the
  rest carry deletion, insertion or absence values (1,995 deletions with substitutions, 1,706
  insertions with substitutions, 1,105 all three, 1,600 indels alone, 137 a block absent in
  some genomes); 5,955 overlap a TRExplorer tandem-repeat locus.
- **Frequencies.** 49,764 of the 60,750 substitution and indel events (81.9%) have a gnomAD
  site. Panel minor share against gnomAD minor allele frequency: Spearman 0.69. The two are
  kept in separate fields.
- **An example.** chr21:33,160,383, in a constrained_unknown block, is an (AC)n repeat of six
  copies in hg38. In the panel, 55 of 89 assemblies carry 6 more bases (nine copies), 33 carry
  hg38's six, one carries 8 more. TRExplorer's TenK10K histogram has nine copies in 1,646 of
  3,850 alleles and six in 190.
- **Hypervariable units.** 774, and 758 of them sit on a tandem repeat. Read by length allele,
  359 become storage columns of copy number and 406 stay hypervariable.
- **The executor shortlist** (Albert's third class, programs that read a stored value, is not
  tested here). 2,185 storage units hold a GTEx fine-mapped eQTL variant and 107 an MPRAVarDB
  allele pair: 2,216 in all. 1,794 are regulatory tier, 225 fossil, 194 neutral, 1
  constrained_unknown. 1,509 are early-replicating, 768 two-value, and 1,053 are early with
  at most three values: the strongest candidates, since early DNA that is still diverse is not
  explained by mutational input. The full catalogue is local
  (`data/knowledge/human_panel/chr21/storage_catalogue.json.gz`); the result keeps counts and
  a showcase.

Storage is not rare. Half of all placed 200-bp units anywhere on chr21 hold a recurring
alternative, and no unknown tier holds more than its matched background. At this resolution
a column with a few recurring values is the default state of human sequence under neutral
drift; what marks function is depletion against matched windows (the coding exons) and,
possibly, a read-out that depends on the value (the shortlist).

**Calibration: slots whose values are known.** Each locus is a regional store of the same
alignment (±20 kb, the flanks as the local rate). Positions were looked up and verified
against the sources (rsIDs in gnomAD, CDS from GENCODE, repeat loci from TRExplorer by motif
size).

| locus | unit | class | recurring values | commonest | effective | by length | timing | rate vs flanks |
|---|---|---|---|---|---|---|---|---|
| rs12913832, HERC2 enhancer of OCA2 | 200 bp | **storage** | 2 (G in 13 of 89; gnomAD G 0.487, AN 152,216) | 0.85 | 1.33 | fixed | middle | 0.92 |
| rs4988235, MCM6 enhancer of LCT | 200 bp | **storage** | 2 (A in 8; gnomAD A 0.397, AN 152,102) | 0.91 | 1.20 | fixed | early | 0.86 |
| ABO | canonical CDS | hypervariable (storage at 16) | 10 | 0.21 | 6.93 | **storage**: 64 at hg38's length, 25 a base longer; rs8176719's frame-restoring C, absent from hg38's O allele, is in at least 31 of 89 (gnomAD 0.354, AN 151,710) | early | 2.03 |
| HLA-A | canonical CDS | hypervariable (storage only at 32) | 19 | 0.10 | 19.6 | fixed | early | 1.80 |
| DRD4 exon 3 48-bp VNTR | locus ±50 bp | **storage** | 8 | 0.43 | 4.09 | **storage**: 4 copies in 58, +144 bp (7 copies) in 22 | early | 5.04 |
| SLC6A3 3' UTR 40-bp VNTR | locus ±50 bp | **storage** | 7 | 0.61 | 2.41 | **storage**: 10 copies in 62, -39 bp (9) in 24, -77 bp (8) in 3 | early | 2.98 |
| INS promoter 14-bp VNTR | locus ±50 bp | hypervariable | 12 of 66 | 0.07 | 44.8 | hypervariable | early | 11.3 |

Five of seven come out as Albert means. The two eye-colour and lactase switches are two-value
columns, rs12913832 in a kilobase Gnocchi does not score. ABO is a two-value column by length
and a ten-value one on substitutions. The DRD4 and SLC6A3 repeats are
copy-number columns with the textbook alleles. HLA-A is a twenty-value slot and is filed as
one, not as a small column. The failure is the INS VNTR: the alignment cannot represent the
2-kb class III expansions, breaks them into replaced blocks and reads a value per haplotype.
Long tandem repeats need the assemblies read directly or a repeat-aware caller. All seven loci
are early or middle-replicating and more variable than their flanks, which is the pattern of a
held value rather than of mutational input.

**Against Gnocchi.** Per kilobase outside coding exons, 90% or more aligned (34,392 kb):
"constrained" is Gnocchi's mean Z at 2.18 or more; "depleted" is the panel's Poisson lower
tail at 0.05 or less against the GC and timing matched rate.

| kilobases | n | panel ratio | Gnocchi Z | late | early | GC | duplicated | missing assemblies | in unknown space |
|---|---|---|---|---|---|---|---|---|---|
| both | 755 | 0.25 | 3.21 | 3% | 68% | 0.485 | 0.8% | 0.0 | 19% |
| Gnocchi only | 2,235 | 0.91 | 2.98 | 2% | 65% | 0.454 | 1.5% | 0.0 | 18% |
| panel only | 3,765 | 0.23 | -0.53 | 42% | 22% | 0.390 | 3.1% | 0.0 | 31% |
| neither | 16,796 | 0.93 | -0.75 | 38% | 22% | 0.384 | 3.3% | 0.0 | 34% |
| Gnocchi unscored, panel depleted | 1,854 | 0.20 | | 51% | 40% | 0.465 | 58% | 45% | 44% |
| Gnocchi unscored, other | 8,987 | 1.71 | | 31% | 41% | 0.426 | 23% | 6% | 33% |

- **Agreement is weak.** Spearman between Z and the panel ratio is -0.08. A Gnocchi-constrained
  kilobase is panel-depleted 25% of the time, against 18% for the rest (1.4 times).
- **The Poisson tail cannot be trusted per kilobase.** 6,374 kilobases are called depleted
  where 1,091 are expected, and the counts are 15.6 times overdispersed. Recurring variants
  share genealogies, so a kilobase with a shallow local tree carries few of them. The panel's
  per-kilobase depletion is mostly that, which is why blocks are judged against matched
  windows and not by their tail.
- **Where the two disagree.** Gnocchi's constrained kilobases are early-replicating and
  GC-rich, next to genes. The panel-only ones are late and GC-poor, where 89 haplotypes' tree
  variance dominates. Gnocchi's rare-variant model sees purifying selection the panel's common
  variants cannot: 2,235 Gnocchi-constrained kilobases carry a normal load of common
  variation. Where Gnocchi is silent and the panel reads depleted, 58% of the kilobases are
  segmental duplications with 45% of assembly-bases missing: copies, where neither axis
  speaks (the organiser's copies-first rule again).
- **By block.** None of chr21's blocks is Gnocchi "syntax". Of 13 core blocks with a case, 3
  are human-constrained ("recent") and 10 tolerant; of 286 variable blocks, 85 are
  human-constrained. Human-lineage invariance at block scale and Gnocchi's aggregate
  constraint do not agree on this chromosome.

**The 69 syntax candidates against their own flanks.** 6 core, 59 variable, 4 unplaced. Their
recurring variation is 1.21 times their ±20 kb flanks pooled (median 1.03). Blocks constrained
across mammals and in Gnocchi carry no less common human variation than the sequence around
them. That is the same negative genomeos-i1's lexicon reached with coarser axes, now at base
resolution: the candidates are not invariant among people. By reading: regulatory 3 core of
22 callable, promoter-like 1 of 7, unexplained 1 of 24, 3' extension 1 of 7, coding exon 0 of
5.

**Sensitivity.**
- **Unit width dominates** the fixed and storage split. Fixed share of placed CDS against
  neutral units: 82.8% against 68.3% at 100 bp; 70.3% against 48.2% at 200; 45.9% against
  19.2% at 500; 29.6% against 5.8% at 1 kb. The gap holds at every width; the absolute
  catalogue sizes do not.
- **The fixed bar** moves both sets together: at 90%, CDS 76.1% and neutral 56.3%; at 99%,
  58.7% and 36.1%.
- **The value-count bar** of 4, 8 or 16 moves only the hypervariable share (neutral 6.0%,
  1.4%, 0.3%).
- **Core ratio** of 0.3 instead of 0.5: genes 51% core of callable, neutral 5.3%, regulatory
  0%, constrained_unknown 10%.

**Cost of the genome.**
- **Stream.** 266.6 GB of MAF: 5.3 hours at the measured 13.9 MB/s, about 1 GB of stores.
- **Memory.** 5 to 6 GB for chr1 and chr2 with their events held in memory.
- **gnomAD.** The storage units' frequencies cost 200 MB on chr21, about 16 GB genome-wide.
- **Timing.** Replication timing costs 6.5 MB per chromosome.
- **Wall time.** The chr21 run took 8.6 hours end to end. Under ten minutes was the
  chromosome itself; the rest was the regional calibration and candidate stage, unprofiled,
  which is the part to fix before a genome run. The genome is affordable in streaming and
  disk; it waits on that fix and on a decision about whether a panel of 89 haplotypes is worth
  reading everywhere after the negatives above.

**What is weak.**
- **The sample.** 44 people. The panel places common alleles but does not estimate their
  frequency: rs12913832's G is 15% in the panel and 49% in gnomAD.
- **Tail probabilities.** The Poisson tails are optimistic by the overdispersion.
- **Missing data.** The expectation scales with aligned bases, not with informative
  assemblies, so partly missing regions read as depleted.
- **One alignment.** Presence in a duplicated block is paralogous, and long VNTRs break.
- **The structural tracks.** gnomAD SVs include rare megabase calls (any gnomAD SV item covers
  93% of the regulatory tier, common ones 3.7%); HPRC insertion items are points.
- **One chromosome**, small, with no Gnocchi syntax block.

## Element types from deletion behaviour, against the registry's classes (2026-09-14)

Code: `genomeos/attribution/element_types.py`, `scripts/element_types.py`. Result:
`element_types.json`. No model call.

**The question.** The registry's classes are a weak guide to what an element does. The class
predicts whether a deletion lowers or raises the gene at AUC 0.532 (the epigenome layer's table
in NODES-READER-WRITER.md). Do categories read off measured behaviour do better than the labels
we inherited?

**The vector.** Each element gets eleven numbers, taken only from the deletion archive: every gene
in the 1 Mb window, on the RNA-seq tracks of K562, HepG2, GM12878 and IMR-90. They record:

- how strong the largest effect is, and its sign;
- the net direction over all genes and lines;
- how many genes move by 0.1 log2 or more, and what share of the window that is;
- how many lines act, and whether acting lines disagree in sign;
- the spread between the lines' strongest effects;
- whether the most-moved gene is the nearest TSS, and how far its start is;
- the mean effect over the window.

No class, mark, DNase, sequence or constraint value enters. The data are the 20,954
enhancer-like elements of chr21 and chr22 (the fit set). The transfer set is the 37,790 of chr19
and chr20, and chr17, chr18 and chrY are assigned for shares only: 101,278 elements in all.

**The choice of k, fixed before any result.** k-means with k from 2 to 10. The rule: take the
highest mean silhouette on a fixed 5,000-element subsample, among the k whose smallest cluster
holds at least 2% of the elements. It picks k = 4, at silhouette 0.277. The rule is weak: k = 2
scores 0.273, and on chr19 and chr20 the same rule picks 2.

**The four types (chr21 and chr22).**

| type | share | what the deletion does | registry mix (dELS / pELS) | CTCF-bound | open (DNase, per line) | with conserved bases |
|---|---|---|---|---|---|---|
| 0 inert | 51.5% | nothing reaches 0.1 in any line (strongest effect about 0.05) | 80 / 20 | 41% | 15% | 35% |
| 1 silencer-like | 16.5% | expression rises; about 40% of lines act; all agree | 67 / 33 | 48% | 15% | 37% |
| 2 activator-like | 26.2% | expression drops; about 47% of lines act; all agree | 65 / 35 | 49% | 23% | 47% |
| 3 line-discordant | 5.7% | rises in some lines and drops in others; 68% of lines act, most genes moving | 57 / 43 | 58% | 28% | 55% |

The types are behavioural. They are not the registry renamed: Cramér's V between type and class
is 0.17, against 0.39 between class and measured chromatin state. They are not geometry either.
Length, GC and distance to the nearest gene predict membership at AUC 0.55 to 0.63. A k-means
on those three features agrees with the types at adjusted Rand 0.017. Of the eleven features,
only the distance to the most-moved gene is much explained by geometry (out-of-fold r² 0.41);
the rest reach at most 0.09. The budget tier says little. 83 to 89% of every type sits in
annotated sequence outside the UNKNOWN blocks, and type against tier gives V = 0.04.

Type 3 is the interesting one. The registry cannot see it, and it stands out on every axis the
vector was blind to. It is the most CTCF-bound, most open and most conserved type, with the
nearest genes (median 3.5 kb). Beyond what its class mix predicts, its chromatin state carries
7.7 points less "no mark" and 3.6 points more active promoter. It is also what one would expect of
promoter-proximal elements with several genes within reach, which is what makes the model's
direction differ from line to line: a cluster to read as an observation about the model's window,
not yet as a kind of element.

**Controls.**

1. *Features permuted one by one*, keeping every marginal and destroying the dependence between
   features: silhouette falls from 0.277 to 0.090 at k = 4, inertia rises from 122,118 to 188,755,
   and the criterion picks k = 8 on the permuted cloud. The joint structure is real.
2. *Whole vectors shuffled between elements of the same class, CTCF flag and number of open lines*
   (20 strata: two classes by the CTCF flag by nought to four open lines). This cannot change the cloud, so the silhouette is untouched by construction;
   what it destroys is which element carries which behaviour. Type against measured state falls
   from V = 0.098 to 0.063, type against conserved bases from 0.126 to 0.029, type against tier
   from 0.040 to 0.016, and type against class stays at 0.174 (it is held fixed by the strata).
   So the associations with chromatin and with constraint are carried by the element, not by its
   class and openness; they are also small in absolute terms.
3. *Transfer.* Scale and centroids fitted on chr21 and chr22, applied unchanged to chr19 and
   chr20: adjusted Rand 0.996 against a refit on those chromosomes, 99.9% of elements in the
   matched cluster, silhouette 0.275. Shares are stable over all seven chromosomes: inert 46 to
   55%, silencer-like 15 to 17%, activator-like 24 to 30%, discordant 4.9 to 9.4%. The partition
   transfers; the *number* of types does not, since the criterion picks 2 on the transfer set.
4. *Geometry alone*: above.

**Does a discovered type beat a curated one?** The same three questions, the same units (element
by cell line), the same metric as the epigenome layer. A type is assigned to an element for a
given line from the other three lines' behaviour only, so the outcome never enters its own
feature. Held out two ways: five folds grouped by element within chr21 and chr22 (which
reproduces the layer's published numbers exactly), and fitted on chr21 and chr22 then scored on
chr19 and chr20.

| features | acts (AUC) | rise among acting (AUC) | magnitude (Spearman) |
|---|---|---|---|
| registry class | 0.626 / 0.639 | 0.532 / 0.523 | 0.183 / 0.220 |
| + DNase in the line | 0.645 / 0.658 | 0.594 / 0.584 | 0.217 / 0.251 |
| + the line's marks and methylation | **0.676** / — | **0.616** / — | **0.333** / — |
| + the line's mark peaks and methylation | 0.668 / 0.678 | 0.610 / 0.604 | 0.307 / 0.326 |
| discovered type (from the other three lines) | 0.644 / 0.659 | 0.669 / **0.686** | 0.195 / 0.206 |
| discovered type on top of the curated features | 0.702 / 0.717 | 0.698 / 0.705 | 0.342 / 0.368 |
| control: type of another element of the same class and openness | 0.551 / 0.534 | 0.521 / 0.514 | 0.105 / 0.094 |
| the same behaviour vector uncompressed | 0.772 / 0.809 | 0.842 / 0.870 | 0.433 / 0.506 |
| two of its columns: the direction in the other three lines | 0.575 / 0.578 | 0.838 / 0.865 | 0.122 / 0.130 |

Each cell is the within-chromosome cross-validation and the transfer to chr19 and chr20. The
marks-with-fold-change row has no transfer figure: no signal profile is cached beyond chr21 and
chr22, so that model meets constant columns there and its transfer score (0.08 for magnitude) is
an artefact. Mark peaks plus methylation is the curated comparison on held-out chromosomes.

95% intervals over 200 element resamples of the transfer set, with the share of resamples showing
no gain:

| comparison | acts | rise among acting | magnitude |
|---|---|---|---|
| type - curated (peaks) | -0.024 to -0.015 (1.00) | **+0.072 to +0.092 (0.00)** | -0.127 to -0.113 (1.00) |
| type + curated - curated | +0.037 to +0.041 (0.00) | +0.095 to +0.107 (0.00) | +0.040 to +0.044 (0.00) |
| type - shuffled type | +0.121 to +0.129 (0.00) | +0.165 to +0.180 (0.00) | +0.107 to +0.117 (0.00) |
| type - the two direction columns | +0.076 to +0.086 (0.00) | **-0.184 to -0.173 (1.00)** | +0.070 to +0.082 (0.00) |

**The answer, and why it is not the win it looks like.** On direction, the one question the
registry answers at chance, discovered types beat the curated features by 0.07 to 0.09 AUC
(0.686 against 0.604 on held-out chromosomes, 0.669 against 0.616 in the layer's own
cross-validation). On the other two questions they lose: 0.02 worse on whether an element acts
and 0.12 worse on how large the effect is. The gain survives the shuffle control, so it belongs
to the element and not to its stratum. But it is not a discovered *category* doing the work. Two
columns of the same vector — the sign and the net direction of the deletion in the other three
lines — score 0.865, 0.18 above the four types, and the whole vector 0.870. Compressing
behaviour into four names throws away most of what makes it predictive, and the k the criterion
picks is not even the best k for the task: types alone on the transfer set score 0.515 at k = 2,
0.686 at k = 4 and 0.831 at k = 6, without a rule that would have chosen 6 in advance.

So the honest reading is that measured behaviour predicts measured behaviour: what "beats" the
registry class is the model's own direction in three sibling cell lines, and the types are a
lossy summary of it. As labels for the attribution the four types are worth having (they are
stable, they transfer, and they name the discordant 5.7% the registry cannot see); as a
classifier they are neither the cheapest nor the best use of the same data.

**The caveat that limits all of it.** Every number in the vector is AlphaGenome's prediction, and
AlphaGenome was trained on ENCODE histone ChIP-seq, DNase and RNA-seq of these same four lines.
Where a type agrees with ENCODE-derived marks, part of that agreement is the model reading its
own inputs back. And the direction result is a statement about the model's consistency across its
own tracks, not about an element measured in a cell: "rise on deletion" is the model's silencer
call, and nothing here is a deletion experiment.

**Cost.** 765 s wall clock, 3.1 GB peak, 353 MB read over 49 files, 101,278 elements and 234,970
element-cell units. Building the vectors and the clusters — the classification itself — is 33 s
of that, 0.32 s per thousand elements; the rest is the epigenome reading (9 s) and the scoring
with its bootstrap (709 s). Genome-wide, scaling linearly by sequence: about 709,000 elements,
2.5 GB read, under 4 minutes of clustering, and peak memory set by the largest chromosome's
archive rather than by the genome, since chromosomes are read one at a time. That is the cheap
half. The deletion scoring these vectors are made of cost 395,510 s of AlphaGenome time for these
101,278 elements (3.9 s each), which extrapolates to about 770 hours for the genome — one
request per element, already spent on these ten chromosomes by the chain. A curated class needs
no model call at all, which is the real cost difference between the two ways of naming an
element.

## Many human genomes, replicated on chr22, and the slow stage fixed (2026-09-14)

Two things had to hold before the human panel could be trusted at scale: the analysis had to
stop costing hours, and chr21's results had to survive a second chromosome. chr22 is fully
worked by every other layer and about the same size (3.39 GB of alignment, streamed in 484 s on
a loaded link: 25,531 blocks, 1,107,993 variable columns).

**The slow stage.** The chr21 run took 8.6 hours, 8.5 of them in calibration and candidates.
A profile showed the analysis was never the cost: 0.1 s per locus. The time went to range
requests, each of which opened a new TLS connection. On UCSC's download host a handshake cost
0.5 s when quiet and 4 to 17 s when loaded; a chr7 subset of seven candidates spent 541 s on 99
requests.
- **The fix.** One kept-alive connection per host and thread, with patient retries and a 60 s
  timeout. A failed request alternates with UCSC's second download host, which serves the same
  files (sizes checked for every file used), and a size mismatch between hosts raises. Each
  bigBed reader is opened once per track; chains and GENCODE are parsed once per chromosome.
- **Replication timing** now comes from the 11 Repli-seq bigWigs fetched once whole (108 MB,
  git-ignored), cached per chromosome. chr22's timing took 6 s.
- **Stalls.** The script times every stage and records a stalled stream as pending, and
  `--fill-pending` reruns only those stages. chr22 needed it once: its gnomAD and Gnocchi stages
  timed out and were filled in 319 s.

| chr21, same store | before | after |
|---|---|---|
| blocks, controls, units | about 90 s | 52 s |
| storage catalogue (200 MB of gnomAD ranges) | 366 s | 1,147 s on a slower link |
| calibration loci | about 8.5 hours with the candidates | 91 s |
| the 69 candidates | | 2 s |
| whole run | 31,099 s | 1,392 s |

Every committed number reproduced exactly: controls, catalogues, tiers, units, matched windows,
calibration and candidates. What is left is transfer: gnomAD's 200 MB for the storage units
dominates and scales with the link. A genome run is now about the 5.3 hours of streaming plus
roughly 16 GB of gnomAD ranges, not days of handshakes.

**chr22 against chr21, claim by claim.**

| claim | chr21 | chr22 | holds? |
|---|---|---|---|
| CDS pooled ratio to the GC and timing matched rate | 0.411 | 0.423 | yes |
| neutral tier pooled ratio | 0.907 | 1.045 | yes (both near 1) |
| CDS 200-bp units fixed, against matched expectation | 70.3% against 44.1% | 71.4% against 46.8% | yes |
| callable genes core | 58.2% (153 genes) | 54.9% (359 genes) | yes |
| share of the CDS-neutral gap surviving the timing match | 82% (0.61 to 0.50) | 71% (0.87 to 0.62) | yes, smaller |
| late against early rate within a GC stratum | 1.1 to 1.5 times | 1.1 to 1.6 times | yes |
| storage the default: placed background units holding a recurring alternative | 52.0% | 53.1% | yes |
| no unknown tier more fixed than its matched expectation (units) | none | none: fossil 42.0/45.3, regulatory 40.0/45.8, constrained_unknown 41.4/44.7, neutral 41.1/42.6 | yes |
| unknown blocks called core, against matched windows | fossil 1.2/4.2, regulatory 3.5/6.9, CU 40/10, neutral 7.9/3.7 | fossil 5.5/15.2, regulatory 7.5/10.1, CU 16.7/28, neutral 10.0/20.7 | yes, and chr21's CU excess (two copies) does not recur |
| unit classes even over timing | flat (storage 34% late, fixed 33%) | storage 37% late against fixed 31% | **no**: on chr22 variable units lean late |
| Gnocchi against panel per kilobase, Spearman | -0.08 | -0.02 | yes, weaker |
| panel-depleted given Gnocchi-constrained, against the rest | 25% against 18% | 28% against 26% | yes, no agreement on chr22 |
| Gnocchi's constrained kilobases early and GC-rich | 65-68% early, GC 0.45-0.49 | 32-44% early, GC 0.51 | **partly**: GC-rich yes, early much weaker |
| the panel's depleted-only kilobases late and GC-poor | 42% late, GC 0.39 | 53% late, GC 0.44 | yes |
| overdispersion of recurring counts per kilobase | 15.6 | 26.2 | yes, larger |
| Gnocchi-unscored, panel-depleted kilobases are copies | 58% duplicated, 45% missing | 47% duplicated, 29% missing | yes |
| core blocks Gnocchi-constrained, against variable blocks | 3 of 13 against 85 of 286 (23% against 30%) | 18 of 40 against 130 of 461 (45% against 28%) | **no**: on chr22 they agree at block scale |

**What replicates.** Coding exons sit at 0.41 to 0.42 of the matched rate on both chromosomes,
and most of that gap is not replication timing. Storage is the background state on both: about
half of all units carry a common alternative. No tier of the unknown space holds more fixed
sequence than matched background on either chromosome, and on chr22 its blocks are less often
core than matched windows in every tier. The two human axes barely agree per kilobase. Counts
are overdispersed well beyond a Poisson, and Gnocchi's silent, panel-depleted kilobases are
copies.

**What does not replicate.**
- **Timing and unit class.** On chr22 variable units lean late (37% against 31% for fixed), so
  timing is part of what separates fixed from storage units there, even though the coding
  signal survives it.
- **Gnocchi and replication timing.** The early-replication skew of Gnocchi's constrained
  kilobases is a chr21 feature.
- **Gnocchi at block scale.** Here the two measures agree modestly on chr22 and not on chr21.
  Block-level agreement with Gnocchi is not a conclusion either way from two chromosomes.

**chr22's catalogues and value domains** (200-bp units, unknown space):
- **Fixed** 3.29 Mb, **storage** 4.79 Mb, **cannot place** 3.75 Mb (148 kb hypervariable).
- **Domains.** 8,702 two-value and 7,890 three-value storage units; median effective values 1.71.
- **Frequencies.** 69.8% of events matched in gnomAD (chr21 81.9%); panel minor share against
  gnomAD MAF, Spearman 0.70 (chr21 0.69).
- **Hypervariable units.** 742, and 391 become copy-number columns when read by length.
- **The executor shortlist.** 4,034 units: 3,968 with a GTEx fine-mapped eQTL, 245 with an
  MPRAVarDB allele pair.
- **Presence** is lower in chr22's unknown tiers (fossil 0.968, constrained_unknown 0.927,
  against 0.998 and 0.997 on chr21). Its low-copy repeats likely make more of the unknown space
  structurally polymorphic, and its tiers read above the background rate (fossil
  1.22, constrained_unknown 1.74).

**chr22's own candidates.** The organiser's real unknown on chr22, 10 constrained-unknown blocks
that are not copies, each read against its own ±20 kb flanks: all 10 variable, pooled ratio 0.92,
median 0.77. The one genome-wide syntax candidate on the chromosome, chr22:38,571,213, is the
most depleted: 5 recurring variants where its flanks predict 13.1, ratio 0.38, a Poisson lower
tail of 0.0101, just above the 0.01 bar. It is the single block of the set worth a second look,
and with a 15- to 26-fold overdispersion that tail is not evidence on its own. chr21's three listed
real-unknown blocks read the same way, all variable against their flanks (pooled 1.14).

**Which human measure to trust.** For "is this sequence held among people", trust Gnocchi over
the panel. Gnocchi counts depletion of mostly rare variants in 76,156 genomes against a mutation
model. The panel counts common variants in 89 haplotypes, whose number per kilobase follows the
local genealogy (15- to 26-fold overdispersed) more than selection. The panel's recurring-variant
depletion is interpretable only pooled over many kilobases, as for coding exons and tiers, or
against matched windows, never per kilobase or per short block. For single bases use phyloP.

The panel is the measure to trust for what Gnocchi cannot see:
- presence and structure;
- indels and copy number;
- the alleles themselves and how they combine into values;
- sequence Gnocchi leaves unscored, where it mostly reveals copies.

**What is weak.**
- **Two chromosomes, both small** and both acrocentric, with short arms the panel cannot align.
- **The second host.** The mirror is trusted on file size, not on content.
- **The before timing** is from a run the infrastructure stall may have lengthened. The
  profile's per-request handshake cost is the measured cause either way.

## Many human genomes, twenty chromosomes: what holds and what was the chromosome (2026-09-15)

Two chromosomes could not carry the conclusions, and three of chr21's readings had already failed
to replicate on chr22. The panel has now run on chr5 to chr22, chrX and chrY: 20 of 24
chromosomes, 186.3 GB of alignment streamed in 13,418 s at 13.9 MB/s and 9.0 hours of run time
over two lanes. chr1 to chr4 are **not** read: they are 79.5 GB of the alignment, 30% of it, and
the sweep was stopped four times by free disk, never by anything the panel costs. Every claim
below is stated with its
spread over the chromosomes and the chromosomes it fails on, never as one pooled number; the
reader is `genome_wide()` in `attribution/human_panel.py`
(`scripts/human_panel.py --genome-wide`, `human_panel_genome_wide.json`).

**chrY is read and pooled with nothing.** Only 19 of the 89 HPRC haplotype assemblies carry it, so
a recurring value there is two of 19 rather than two of 89, a ten times higher allele bar. 94% of
its placed units read fixed against 48% elsewhere, its coding exons sit at 1.66 of the matched
rate instead of below it, no callable gene is core and no storage unit carries an eQTL or an MPRA
pair. That is a property of the panel's depth on a male-only chromosome, not of the sequence.
Every pooled figure below is over the 19 chromosomes that are not chrY.

**What holds on every chromosome read.**

| claim | min | median | max |
|---|---|---|---|
| coding exons, ratio to the GC- and timing-matched recurring rate | 0.292 | **0.379** | 0.444 |
| coding 200-bp units fixed above their matched expectation | +22.1 pts | **+26.2** | +29.7 pts |
| callable genes called core | 42.0% | **59.0%** | 70.7% |
| genes core above their own matched windows | +34.3 pts | **+49.0** | +61.4 pts |
| share of the coding-neutral gap surviving the timing match | 71% | **90.9%** | 98% |
| storage share of placed background units | 35.9% | **48.5%** | 53.4% |
| overdispersion of recurring counts per kilobase | 7.4 | **9.8** | 26.2 |

The instrument's ordering is not a chromosome. Coding exons carry between 0.29 and 0.44 of the
matched rate of recurring variation on all nineteen, coding units are fixed 22 to 30 points above
matched expectation on all nineteen, and most of that gap survives the replication-timing match on
all nineteen (the 71% floor is chr22, which was the worst case at two chromosomes too). Storage is
the default state of human sequence everywhere: on the eighteen autosomes read it is between 47.2%
(chr5) and 53.4% (chr19) of placed 200-bp units, and chrX alone sits lower.

**What was the chromosome.**

| claim | chr21, chr22 said | over 19 | verdict |
|---|---|---|---|
| Gnocchi's silent, panel-depleted kilobases are copies | 58%, 47% duplicated | 8% to 58%, median 25%, fails on 12 | **chromosome-specific** |
| core blocks agree with Gnocchi more than variable blocks | chr21 -6.6 pts, chr22 +17 pts | +5 to +29 pts, fails only on chr21 | **chr21 was the outlier: they agree** |
| unit classes spread evenly over replication timing | flat on chr21, late-leaning on chr22 | 0.4 to 6.6 pts, fails on 5 | **no, and the direction is not fixed** |
| no unknown tier more fixed than its matched expectation | none on either | fails on 5 | **no: constrained_unknown is** |
| the neutral tier within a quarter of its matched rate | 0.91, 1.05 | 0.85 to 1.35, fails on 4 | **mostly, and always by reading high** |
| Gnocchi's constrained kilobases are early-replicating | 65-68% on chr21, 32-44% on chr22 | 36% to 70%, median 51%, fails on 9 | **weak everywhere, strong nowhere but chr21** |

- **The copies reading does not generalise.** chr21 and chr22 both said that where Gnocchi is
  silent and the panel sees depletion, the sequence is duplicated: 58% and 47%. Over 18
  chromosomes the median is 25% and twelve fall below 30% (chr18 0.08, chr12 0.086, chr19 0.086,
  chr8 0.10, chrX 0.105). It is high where segmental duplication is high (chr9 0.57, chr16 0.45,
  chr22 0.47, chr21 0.58) and low elsewhere. It was a fact about two acrocentrics.
- **Gnocchi and the panel do agree at block scale.** chr21 said they did not (-6.6 points) and
  chr22 said they did (+17). On eighteen of nineteen chromosomes core blocks are Gnocchi-
  constrained more often than variable blocks, by 5 to 31 points, median 20. chr21 is the only
  chromosome where the sign goes the other way. Per kilobase they still barely agree: Spearman
  between -0.02 and -0.21, median -0.08, on every chromosome.
- **Gnocchi's constrained kilobases and replication timing.** chr21 said its constrained
  kilobases were 65 to 68% early-replicating and chr22 said 32 to 44%. Over eighteen the early
  share runs from 0.362 (chr17) to 0.698 (chr8) with a median of 0.509. A tertile would give 0.33,
  so the skew towards early DNA is there on every chromosome; what is not there is the strength
  chr21 showed, which is the top of the range. Half the chromosomes do not reach 50%. The honest
  statement is a weak early skew everywhere, not a property that separates constrained kilobases.
- **Replication timing and unit class is unsettled.** The gap between the late share of storage
  units and of fixed units runs from 0.4 to 6.6 points and exceeds 3 points on five chromosomes.
  Two chromosomes could not settle it and nineteen do not either; the coding signal survives the
  timing match regardless, which is the part that matters.

**A new reading the two chromosomes could not give: the constrained-unknown tier.** Pooled ratio
to the matched rate, per chromosome, median over the nineteen:

| tier | median ratio | below 1 on | sign test |
|---|---|---|---|
| canonical CDS | 0.379 | 19 of 19 | — |
| constrained_unknown | **0.916** | 15 of 19 | p = 0.0096 |
| fossil | 1.054 | 1 of 19 | — |
| regulatory | 1.070 | 0 of 19 | — |
| neutral | 1.108 | 2 of 19 | — |

This is the first thing in the unknown space that reads as held, and it is small. Read against 1
it is an 8% depletion; read against the neutral tier, which is the honest comparison because the
matched background under-predicts every unknown tier by about a tenth, it is 0.83 where coding is
0.34. chr21's 0.649 overstated it and chr22's 1.742 inverted it; neither chromosome could see the
effect for what it is. The tier is also more fixed at the unit level than matched expectation on
twelve of nineteen chromosomes (chr19 +6.5 points, chr11 +6.4, chr8 +6.1, chrX +5.0). The fossil
and regulatory tiers read *above* their matched background on 18 and 19 of 19: whatever the
panel's matched windows are, the unknown space is not less variable than them.

**chrX is a third instrument.** Storage is 35.9% of its placed units against 47 to 53% everywhere
else, and fixed 61.4% against about 50%; only 3.1% is unplaced, so this is not a coverage
artefact. Its coding exons still sit at 0.342 and its genes are core 69.3% of the time, so the
controls hold; it is the diversity that is lower, as a chromosome with three quarters of the
autosomal effective population size should be. It is reported with the autosomes and flagged here.

**The genome-wide storage catalogue** (19 chromosomes, 200-bp units of the unknown space):

| catalogue | Mb |
|---|---|
| fixed (one value in 95% of assemblies) | 268.2 |
| storage (a small recurring value domain) | 274.3 |
| cannot place (hypervariable, or not aligned) | 61.0 |

- **1,370,932 storage units.** Value-set sizes: 607,208 hold two recurring values, 452,660 three,
  181,827 four, 69,584 five, 30,888 six, 17,355 seven, 11,266 eight. The median effective number
  of values is **1.65** and the median commonest value holds **75.0%** of the assemblies, the same
  to two places as chr21's 1.68 and 73% and chr22's 1.71. A storage column is a two- or
  three-value column almost always: 77% of them hold two or three values.
- **By tier:** regulatory 616,947, fossil 542,233, neutral 122,060, constrained_unknown 54,135,
  structural 35,557.
- **By event:** 1,016,839 substitutions only; 92,097 deletion with substitutions, 84,274 insertion
  with substitutions, 56,435 all three, 34,266 deletions alone, 33,722 insertions alone, 26,217
  indels together, 13,974 with a structural value. 308,831 sit on a TRExplorer tandem repeat.
- **Frequencies.** 82.0% of the substitution and indel events have a gnomAD site, close to chr21's
  81.9% and above chr22's 69.8%.
- **Hypervariable units:** 34,022, of which 15,280 become copy-number columns when read by length
  allele alone and 18,197 stay hypervariable.

**The executor shortlist is no longer chromosome-shaped.** **134,297** storage units carry a
measured value event: 131,891 with a GTEx fine-mapped eQTL and 6,436 with an MPRAVarDB allele
pair. chr21 and chr22 together gave 6,250, so this is **128,047 new executor-ready units**, a
twentyfold shortlist. Per chromosome it runs from 2,083 (chrX) to 12,875 (chr7); chr19 gives 9,156
from 5.4 GB, the best return per byte streamed, and chr11 alone gives 11,101. The catalogues stay
local and gzipped under `data/knowledge/human_panel/<chrom>/storage_catalogue.json.gz`, which is
where `attribution/executor.py` already reads them from, so the executor lane needs no new format:
it needs a longer `chroms` list.

**Cost, measured.**

| | |
|---|---|
| alignment streamed | 186.3 GB in 13,418 s, 13.9 MB/s on a shared link |
| run time | 9.0 h over two lanes |
| per chromosome | 574 s (chr22) to 3,608 s (chr7); the median is 1,352 s |
| the slow stage | the storage catalogue's gnomAD ranges, 300 to 2,300 s and growing with the chromosome |
| local store | about 40 MB per chromosome, git-ignored |

**What is missing and why.** chr1 to chr4 are unread. Free space on the machine fell to 13 or 14 GB
four separate times, none of it the panel's doing — the whole checkout including every cache is
about 12 GB of a 460 GB disk — and the sweep was held at the 15 GB floor every lane was given. The
panel's own footprint per chromosome is about 40 MB, so this is a disk the panel shares, not a disk
the panel spends. chr4 and chr5 were each killed once in mid-run for it; chr5 was finished on the
fifth attempt and chr4 was not. The four missing chromosomes are the four largest and 30% of the
alignment; the claims above are over 70% of it and every chromosome from chr5 down. Nothing in the
readings varies with chromosome size — chr5, the largest read, sits at the median of every claim —
so the four are expected to extend rather than change them, and that is a prediction, not a
result.

**What is still weak.**
- **The tier ratios are not centred.** Every unknown tier but constrained_unknown reads above its
  GC- and timing-matched background, and the neutral tier, the control, reads 1.12. Either the
  matched windows are drawn from sequence that is quieter than the tiers or the tiers are genuinely
  freer. Until that is understood, the tier ratios should be read against the neutral tier and not
  against 1.
- **The constrained-unknown depletion is a sign test on 19 chromosomes**, not an effect size with
  an interval. It is 0.83 of neutral where coding is 0.34 of neutral.
- **The calibration loci were read once, on chr21**, and are referred to rather than re-read; they
  do not depend on which chromosome is being swept, but they have not been re-measured since.
- **chrY's 19 haplotypes** mean the panel has no reading at all for 57 Mb of the genome.

## The library after the exclusions: 284,598 oligos, and two judgements that are design and not measurement (2026-09-17, later)

genomeos-79's lane measured which oligos cannot be read even in principle, and the manifest is rebuilt
on it. The filter is **recomputed here from sequence alone rather than imported**, so the two counts
check each other:

| unorderable, over the untouched blocks | this build | genomeos-79 |
|---|---|---|
| homopolymer ≥ 10 | 3,421 | 3,394 |
| GC outside 25–75% | 749 | 751 |
| low complexity, by sequence (a 6-mer covering half the window) | 220 | — |
| low complexity, by annotation (half the window simple/low-complexity/satellite) | 163 | 105 |
| not unique within its own block | not implemented | 23 |
| **total** | **4,264** | **4,273** |
| **orderable** | **41,037** | **41,297 (90.6%)** |

Two independent implementations, **nine windows apart** on the number a synthesis quote would use.

**Both low-complexity rules are charged, as a union, and neither lane's is the one to keep.** The
annotation rule inherits whatever RepeatMasker did or did not annotate, so an unannotated simple repeat
passes it — and on unannotated sequence that is exactly the failure mode to assume. The sequence rule
fires on the oligo itself, which is closer to what a synthesiser fails on, and misses what a curator
saw in a diverged repeat. The union costs about 0.2% of the library and removes the need for either
side to defend an arbitrary threshold. The within-block uniqueness test is a **lower bound** on
mappability and is not presented as a mappability arm by either build.

**The library is now 284,001 oligos**: 91,919 test (41,037 of them from untouched blocks, 2,277 in the
seven bridge blocks), 91,918 matched genomic negatives, 8,245 positives, 91,919 scrambles.

**Judgement one: the repeat flag is not an exclusion.** Interspersed repeat (17,112 windows) and
segmental duplication (210) are carried **flagged and kept**. The reason is measured, not argued:
**31.3% of the lentiMPRA windows a reporter already measured successfully are themselves majority
interspersed repeat.** Dropping them would drop a class the assay demonstrably reads. It is an
attribution caveat — a hit inside a recent duplication cannot be assigned to one locus — and the
manifest carries the flag so the analysis can condition on it rather than the design pre-empting it.

**Judgement two: the thin blocks stay, and this one is the more dangerous.** 53 of the 680 untouched
blocks fall below three attributable oligos, and they are not a random 53: **median length 1,402 bp
against 7,996, and median distance to a coding TSS 23.6 kb against 100.9 kb.** Dropping them would
systematically remove the blocks nearest genes — which is the covariate that confounded every
comparison this project corrected this week. A library that quietly drops its promoter-proximal
sequence would reproduce the artefact in its own design, and no amount of standardising afterwards
recovers sequence that was never ordered.

**The bridge, and it is thin.** The CRISPRi-calibrated shortlist reaches the real unknown in **seven
blocks and eight elements** — three relaxed, four tolerant, none of the 69 syntax and none of the 29
recent. The strongest moves NXPH1 by 1.15. The library tiles those seven deliberately (2,283 oligos)
and labels them, so the screens' drop band can be tested off their own population; for the other 875
blocks the honest expectation is that no CRISPRi evidence speaks to them at all. The near-absence is
structural rather than biological: the sweep scored registry elements inside a node, and a
real-unknown block is unannotated by construction, so most of its sequence was never eligible.

**Still unmeasured and worth stating before anyone quotes the attributable figure:** genome-wide
mappability. The uniqueness test compares a window only against its own block, so it is a lower bound;
a real read needs a Umap or Bismap track or a whole-genome k-mer index, and neither is held locally.
Repeat and segdup content stand in as proxies, and they are proxies.

## The library, written out: 312,129 oligos that would measure the real unknown (2026-09-17)

The coverage reading below says the sequence this project most wants measured is small enough to
measure. `scripts/unknown_library.py` emits that library rather than proposing it — every oligo with
its coordinates, sequence, block and arm, into a manifest under `data/knowledge/library/` (101 MB,
local, git-ignored). It takes 3.5 minutes and sends nothing anywhere.

| arm | oligos | median GC | median TSS distance | what it is for |
|---|---|---|---|---|
| test | **101,295** | 0.373 | 307 kb | all 30.6 Mb of the real unknown, tiled end to end at 300 bp, **untouched blocks first (45,301)** |
| genomic negative | **101,294** | 0.373 | 355 kb | neutral-tier oligos matched to each test oligo on GC and order of magnitude of TSS distance |
| positive | 8,245 | 0.520 | 44 kb | lentiMPRA elements already measured active, so the batch calibrates against a known answer |
| scrambled | 101,295 | 0.373 | — | each test oligo's own dinucleotide-preserving shuffle: composition held, arrangement destroyed |

**312,129 oligos over 878 blocks**, by case: tolerant 54,590, relaxed 44,588, syntax 1,250, recent
797, unmeasured 70.

**Three things about the design that are the point of writing it out.**

1. **The negatives are matched one to one**, 101,294 against 101,295, with the same median GC to four
   decimal places. The first build sampled the neutral pool four times too coarsely and matched only
   54% of the test arm — a library with a hole where its control should be, and invisible until the
   counts were printed side by side. Every comparison this project got wrong this week was a control
   that was not matched on GC and promoter distance; a designed library is the one place those can be
   fixed in advance rather than standardised afterwards.
2. **The scrambled arm doubles the cost and earns it.** Composition is the one thing the motif lane
   found that transfers (+0.234 to strict family counts, +0.060 AUROC on VISTA), so a test oligo that
   reads active must be read against its own composition, not against the library average. A cheaper
   design scrambles a fifth of the test arm and comes to about 232,000; tiling only the untouched
   blocks with all four arms comes to about 110,000.
3. **The untouched-first ordering is the budget.** If only part of the library can be built, the first
   45,301 test oligos are the sequence that no assay this project holds has ever touched, and the
   design spends down from there.

**What it is not.** A design, not an order: no synthesis, no vendor, and no claim that any of this
sequence does anything. Its value is that it converts "the 98% is unmeasured" from a complaint into a
file with 312,129 rows in it, and that the arithmetic behind every arm is in the script rather than in
a sentence.

## Known defect: every element-level Gnocchi base count in the committed results is about 4.4x too high (2026-09-17)

Recorded before the fix rather than after it, because the wrong numbers are committed and readable now.

genomeos-79 found that `attribution/bigwig.py` credits a whole bin whenever a track's step is greater
than 1. For a step-1 track — Zoonomia phyloP, the Umap mappability tracks — it is exact and nothing
here is affected. For gnomAD Gnocchi, which is one value per kilobase, an interval is credited every
bin it touches in full.

**Measured from the committed results, capping each interval's claimed bases at its own length, which
is an upper bound on the truth — so these are lower bounds on the over-credit:**

| | intervals | claimed | capped at length | at least | claiming more than their own length |
|---|---|---|---|---|---|
| chr21 blocks | 309 | 7,420,000 bp | 7,383,606 | 0.5% over | 50 |
| **chr21 elements** | 231 | **282,000 bp** | **63,912** | **77.3% over** | **231 of 231** |
| chr22 elements | 186 | 219,000 | 49,707 | 77.3% | 186 of 186 |
| chr19 elements | 191 | 218,000 | 49,828 | 77.1% | 191 of 191 |

A 350 bp element in `variation_chr21` reports `"bases": 1000`. A cCRE is shorter than a Gnocchi bin,
so at element level the field is not slightly wrong: it counts bins and calls them bases.

**What is affected, exactly.** In every `variation_chr*.json`: `gnocchi.bases` per block and per
element, and the two fields `attribution/variation.py` derives from it — `by_tier.measured_bp`, and
`by_tier.human_constrained_bp`, which is `fraction_above × bases` and therefore inherits the inflation
although `fraction_above` itself is sound. **Ratios survive**: every mean Z, every `fraction_above`,
every case assignment and every mammalian-constraint figure is untouched, because the numerator and
denominator scale together.

**What is not affected.** No document quotes `measured_bp` or `human_constrained_bp`. What the
documents quote from Gnocchi are means, fractions, and kilobase counts from the human panel, which is
a different reader. The conclusions of area J and of the tier readings rest on ratios and do not move.

**Fixed and measured the same day (`aee3569`), and it is half a bug rather than a bug.** The
touched-bin quantity is legitimate and documented — for a 1 kb track, "bases in every bin this
interval touches" is a real thing to count. The defect was that **one field served two questions** and
the absolute one was summed. So the fix is a name, not a comment: `IntervalStats` now carries both
readings from one pass, `bases`/`above` unchanged, and `overlap_bases`/`overlap_above` counting only
the bases inside the interval, which can never exceed its own length. The phyloP-derived numbers were
exact throughout and "the reader was wrong" would have made them look suspect.

**The after column, chr21, one read of the track for both (7 requests, 0.22 MB):**

| | own length | touched bins | own bases | overstated |
|---|---|---|---|---|
| UNKNOWN blocks (310) | 9,514,145 | 7,444,000 | 7,205,452 | **×1.033** |
| registry elements (9,758) | 2,678,600 | 11,645,000 | 2,520,003 | **×4.62** |

The 4.4× lower bound measured here was ×4.62 in fact. **Every registry element is shorter than a
Gnocchi bin**, so the 350 bp element reporting 1,000 bases was the general case rather than an
example. Constrained bases at element level fall from 3,202,000 to 683,878; the *fraction* moves only
0.275 → 0.2714, which is the ratio surviving as predicted.

**Deltas in `variation_chr21`** (`measured_bp` and `human_constrained_bp` now count the interval's own
bases; the old quantity survives as `touched_bin_bp` under its own name): regulatory 2,654,000 →
2,531,725 with the constrained fraction 0.1692 → 0.1619; fossil 2,661,000 → 2,610,966, 0.0165 →
0.0152; neutral 1,980,000 → 1,919,278, 0.0207 → 0.0189; constrained unknown 121,000 → 116,760, 0.0165
→ 0.0171. No case, no mammal fraction, no mean and no `fraction_above` moved beyond edge-bin
weighting.

**The uncomfortable half of the good news.** The corrected fields are exactly the set nobody had drawn
a conclusion from — which is why no document had to be struck, and also why the number could stay
wrong for weeks. *A figure nothing depends on is a figure nothing checks.*

**Re-run complete, all 24 files agree (2026-09-17).** 9,105 range requests, 134.11 MB, 27.8 minutes,
and `touched_bin_bp` is populated in every file so an old reading and a new one can be compared with
both quantities present. **Genome-wide: `measured_bp` 552.0 Mb → 535.0 Mb (3.2% overstated),
`constrained_bp` 38.46 Mb → 35.46 Mb (8.4%).**

The stop rule — halt the chain if any tier's constrained *fraction* moves more than edge-bin weighting
can explain — was not tripped, and the distribution is the prediction rather than a pass by luck: the
largest move in any tier above 100 kb is **0.0153** (chr19 regulatory), then 0.0109, 0.0106, 0.0102,
0.0095, with chr22 and chrY moving exactly 0. The observed maximum is three quarters of the bar the
rule allowed, which is the argument against loosening it if this runs again.

## The compiled genome states 940,803 facts and not one of them is experimental (2026-09-17)

Area I compiles a chromosome's non-coding space into a BioLang program where every region carries a
role, an evidence kind and a confidence. Only chr21 had ever been compiled, because the inputs for the
rest were unfinished; both sweeps closed this week, so `scripts/compile_genome_programs.py` compiles
all 24 — **227.7 MB of generated program, in under a minute** — and reads the Evidence explorer over
them.

| evidence kind | facts | share |
|---|---|---|
| experimental | **0** | **0.0%** |
| curated | 42,040 | 4.5% |
| predicted | **880,754** | **93.6%** |
| inferred | 18,009 | 1.9% |
| **weak (confidence ≤ 0.5)** | **812,921** | **86.4%** |

**940,803 facts, and the strongest evidence kind is absent entirely.** Per chromosome the weak share
runs from 0.777 (chr19, the most gene-dense) to 0.898 (chr18), with chrY at 0.383 because most of it
is structural and curated rather than attributed.

This is the Evidence explorer's own promise turned on the project's largest artefact, and it agrees
with the coverage reading from the opposite direction: that one measured how little of the unknown
space any assay has touched (0.45%), this one measures how little of what the project *states* about
it rests on measurement (none of it). The two are the same fact seen from the data and from the model.

**A performance decision with a number behind it.** The compiled programs hold 940,803 of the
project's 967,426 stated facts and cost 21 seconds to parse, against 1.8 seconds for everything
hand-written. `evidence.collect` therefore excludes them unless asked (`compiled=1` on the API), and
the exclusion happens before parsing rather than after, because filtering afterwards would have cost
the same twenty seconds. A default nobody waits for is not a default.

**What it does not say.** Nothing here is a claim that the attributions are wrong. `predicted` is an
honest label on a model's output and the project has spent two days measuring how far those outputs
carry — the node's +2.88 points *(Re-audited 2026-09-27, NODES-READER-WRITER.md: +2.90, 95% CI +2.03 to +3.81, on a control matched on boundary count only; +1.2 to +6.6 across four defensible baselines; and **+5.89, CI +3.18 to +8.46, on 661 measured CRISPRi pairs** — the direction stands on measurement, the size quoted here was never identified.)*, the calibration's failure on prevalence, E1's +0.075. It says that
the compiled genome is a hypothesis in the shape of a program, and that it says so itself, fact by
fact, which is what the evidence field was built for.

## The experimental layer: 5.98% of the compiled elements were ever eligible, and 19,282 facts were raised (2026-09-17, later)

The section above ends on a count of zero: the compiled genome stated 940,803 facts and not one of them
was `experimental`. That zero was true of the compiler, not of the project, which already holds real
measurements over some of the same sequence. `genomeos/attribution/measured.py` attaches them, and
`scripts/measured_layer.py` counts what that buys. **It buys 2.01% of the program** — but that number
on its own invites the wrong conclusion, so it is stated third, behind two denominators.

**Three denominators, in order, because a small census is a statement about the search and not about
the genome:**

| | | |
|---|---|---|
| **1. read** | compiled facts parsed at all | **960,094** in 24 programs |
| **2. eligible** | compiled elements inside *any* assay's footprint, at one shared base | **26,316 of 440,377 — 5.98%** |
| **3. raised** | elements a measurement is *of*, under the overlap rule | **19,070 — 4.33% of all, 72.5% of the eligible** |

**414,061 compiled elements have never been covered by a CRISPRi screen, a lentiMPRA library or a
VISTA test at all.** They are not elements the assays contradicted; they are elements nothing asked
about, and they sit in their own bucket (`elements_no_assay_ever_covered`) rather than with the ones
that were measured and disagreed. Of what *was* eligible, **72.5% was raised** — so the overlap rule
is not the bottleneck and neither is the matching: the footprint is. The result of this lane is a
statement about coverage, and only within that 5.98% is it a statement about the predictions.

| evidence kind | before | after | share after |
|---|---|---|---|
| **experimental** | **0** | **19,282** | **2.01%** |
| curated | 42,040 | 42,049 | 4.38% |
| predicted | 880,754 | 880,754 | 91.74% |
| inferred | 18,009 | 18,009 | 1.88% |
| **facts** | **940,803** | **960,094** | |

Nothing moved between kinds. **The predicted fact is never overwritten, because two names answer the
two questions:** `<id>` is what AlphaGenome's deletion predicts and `<id>_measured` is what an assay
measured over the same DNA, each with its own evidence, confidence and source, side by side in the same
program. The 9 new `curated` facts are gene stubs a screen named and the model had not.

**The overlap rule, and it is a constant.** A measurement counts as a measurement *of* a compiled
element only at `RECIPROCAL_OVERLAP = 0.5` — the overlap covers at least half the element **and** at
least half the tested interval. Proximity, similarity and containment never raise a fact. Reported at
three values, genome-wide, with the containments the rule refuses:

| reciprocal overlap | elements matched | CRISPRi | lentiMPRA | VISTA |
|---|---|---|---|---|
| ≥ 0.25 | 22,170 | 2,754 | 19,876 | 252 |
| **≥ 0.50 (the rule)** | **19,070** | **1,505** | **17,869** | **98** |
| ≥ 0.75 | 2,964 | 61 | 2,863 | 40 |
| *contained only, not upgraded* | | *1,495* | *0* | *1,941* |

The last row is the rule doing its job. A VISTA sequence has a median length of 1.4–2.3 kb and a
compiled element 293–299 bp, so **1,941 elements sit inside a VISTA test and are not measured by it**:
a positive over 2 kb says the 2 kb drives expression, not that this 300 bp does. The same refusal costs
CRISPRi 1,495 elements. Between 0.5 and 0.75 the count falls sevenfold, which says the matches are real
overlaps of similar intervals rather than near-identities; below 0.5 it rises by a sixth, mostly VISTA.

**Coverage, as a named number.** `coverage_elements_measured` = **19,070 of 440,377 compiled elements,
4.33%**. Per chromosome it runs from 3.52% (chr13) to 5.46% (chr19), and it is flat: no chromosome is
measured and none is blind. That flatness matters, because it means the floor is a property of the
assays and not of one chromosome's luck. The eligible share is flat for the same reason — 4.29%
(chr13) to 7.69% (chr19) — and the share of the eligible that was raised is flatter still, 63.8%
(chr8) to 82.0% (chr13).

| chromosome | compiled | eligible | share eligible | raised | of all | of eligible |
|---|---|---|---|---|---|---|
| chr1 | 43,681 | 2,923 | 6.69% | 2,123 | 4.86% | 72.6% |
| chr2 | 32,917 | 1,922 | 5.84% | 1,397 | 4.24% | 72.7% |
| chr3 | 27,682 | 1,502 | 5.43% | 1,119 | 4.04% | 74.5% |
| chr4 | 20,646 | 1,007 | 4.88% | 774 | 3.75% | 76.9% |
| chr5 | 22,772 | 1,151 | 5.05% | 845 | 3.71% | 73.4% |
| chr6 | 23,769 | 1,522 | 6.40% | 1,124 | 4.73% | 73.9% |
| chr7 | 22,596 | 1,391 | 6.16% | 982 | 4.35% | 70.6% |
| chr8 | 18,552 | 1,220 | 6.58% | 778 | 4.19% | 63.8% |
| chr9 | 19,272 | 1,046 | 5.43% | 778 | 4.04% | 74.4% |
| chr10 | 21,023 | 1,269 | 6.04% | 923 | 4.39% | 72.7% |
| chr11 | 25,731 | 1,740 | 6.76% | 1,149 | 4.46% | 66.0% |
| chr12 | 22,814 | 1,178 | 5.16% | 888 | 3.89% | 75.4% |
| chr13 | 9,711 | 417 | 4.29% | 342 | 3.52% | 82.0% |
| chr14 | 13,385 | 689 | 5.15% | 503 | 3.76% | 73.0% |
| chr15 | 14,474 | 762 | 5.27% | 581 | 4.01% | 76.2% |
| chr16 | 14,799 | 911 | 6.16% | 669 | 4.52% | 73.4% |
| chr17 | 21,608 | 1,407 | 6.51% | 1,019 | 4.72% | 72.4% |
| chr18 | 9,250 | 550 | 5.95% | 358 | 3.87% | 65.1% |
| chr19 | 17,265 | 1,328 | **7.69%** | 942 | 5.46% | 70.9% |
| chr20 | 12,702 | 819 | 6.45% | 607 | 4.78% | 74.1% |
| chr21 | 5,174 | 259 | 5.01% | 200 | 3.86% | 77.2% |
| chr22 | 9,612 | 571 | 5.94% | 421 | 4.38% | 73.7% |
| chrX | 10,817 | 723 | 6.68% | 542 | 5.01% | 75.0% |
| chrY | 125 | 9 | 7.20% | 6 | 4.80% | 66.7% |
| **pooled** | **440,377** | **26,316** | **5.98%** | **19,070** | **4.33%** | **72.5%** |

Eligible by assay, pooled: lentiMPRA 21,133 elements, CRISPRi 3,908, VISTA 2,375. VISTA is the assay
the overlap rule costs most — 2,375 elements eligible, 98 raised — because its sequences are long.

| chromosome | compiled elements | measured | coverage | CRISPRi | lentiMPRA | VISTA | agree | disagree |
|---|---|---|---|---|---|---|---|---|
| chr1 | 43,681 | 2,123 | 4.86% | 134 | 2,022 | 12 | 408 | 1,634 |
| chr2 | 32,917 | 1,397 | 4.24% | 73 | 1,340 | 6 | 228 | 1,127 |
| chr3 | 27,682 | 1,119 | 4.04% | 79 | 1,056 | 2 | 210 | 854 |
| chr4 | 20,646 | 774 | 3.75% | 15 | 763 | 5 | 126 | 647 |
| chr5 | 22,772 | 845 | 3.71% | 36 | 818 | 5 | 149 | 678 |
| chr6 | 23,769 | 1,124 | 4.73% | 73 | 1,062 | 6 | 178 | 900 |
| chr7 | 22,596 | 982 | 4.35% | 49 | 941 | 7 | 168 | 786 |
| chr8 | 18,552 | 778 | 4.19% | 103 | 689 | 7 | 133 | 566 |
| chr9 | 19,272 | 778 | 4.04% | 15 | 765 | 5 | 151 | 623 |
| chr10 | 21,023 | 923 | 4.39% | 46 | 883 | 4 | 158 | 737 |
| chr11 | 25,731 | 1,149 | 4.46% | 224 | 980 | 3 | 187 | 817 |
| chr12 | 22,814 | 888 | 3.89% | 72 | 826 | 6 | 139 | 699 |
| chr13 | 9,711 | 342 | 3.52% | 0 | 341 | 1 | 63 | 279 |
| chr14 | 13,385 | 503 | 3.76% | 22 | 487 | 3 | 105 | 387 |
| chr15 | 14,474 | 581 | 4.01% | 16 | 566 | 5 | 101 | 473 |
| chr16 | 14,799 | 669 | 4.52% | 55 | 625 | 3 | 139 | 496 |
| chr17 | 21,608 | 1,019 | 4.72% | 88 | 947 | 5 | 186 | 770 |
| chr18 | 9,250 | 358 | 3.87% | 28 | 334 | 5 | 68 | 272 |
| chr19 | 17,265 | 942 | **5.46%** | 211 | 777 | 1 | 244 | 542 |
| chr20 | 12,702 | 607 | 4.78% | 50 | 569 | 4 | 103 | 478 |
| chr21 | 5,174 | 200 | 3.86% | 18 | 188 | 0 | 44 | 146 |
| chr22 | 9,612 | 421 | 4.38% | 14 | 409 | 3 | 103 | 310 |
| chrX | 10,817 | 542 | 5.01% | 84 | 475 | 0 | 136 | 341 |
| chrY | 125 | 6 | 4.80% | 0 | 6 | 0 | 0 | 6 |
| **pooled** | **440,377** | **19,070** | **4.33%** | **1,505** | **17,869** | **98** | **3,527** | **14,568** |

**The disagreement is the number worth reading, and it is reported per assay.**

| assay | elements | agrees | disagrees | not asked |
|---|---|---|---|---|
| CRISPRi | 1,505 | 97 | 31 | 1,377 |
| lentiMPRA | 17,869 | 3,361 | 14,508 | — |
| VISTA | 98 | 69 | 29 | — |

**Two agreement rates, with two names, because one of them is a trap.** Over the 128 elements where a
CRISPRi screen tested the very gene the model named, **75.8% agree**. Over all 1,505 elements a screen
touched, **6.45% agree**. The gap is not noise: in 1,377 of them the screen measured nearby genes and
never measured the one the deletion predicted. The first number is the one a reader quotes and it is
computed over a subset the prediction itself selected; the result stores both under
`agreement_rate_where_the_predicted_gene_was_tested` and
`agreement_rate_over_all_matched_elements`, with the denominator beside each, so the narrow one cannot
be mistaken for the wide one.

**4,613 element-gene pairs were measured as not regulated, and every one is kept.** "No effect in this
screen" is experimental evidence about the absence of an effect, so it is written into the measured
block by name (`measured no effect on SOD1`) rather than dropped or turned back into UNKNOWN. It gets
no `rule`, because a rule would state a relation the assay says is not there — the absence lives on the
element, the relation lives on the rule, and neither field answers the other's question. The 212
regulated pairs do become experimental rules, with the direction read from the sign of the effect: the
screen silences the element, so a gene that falls was being activated by it.
*(2026-09-28: superseded by R2, `1cff24a`. With significant increases and underpowered nulls read as their own
outcomes, 4,613 → 4,565 pairs measured as not regulated, n 128 → 120, observed 0.7578 → 0.8083, ECE 0.4108 →
0.4507. See "CRISPRi negatives read as outcomes" below; the figures here are kept as they were.)*

**What lentiMPRA's 14,508 disagreements do and do not mean.** An episomal reporter measures whether a
200 bp sequence drives transcription out of its chromosome; the compiled claim is that deleting the
element in its chromosome moves a gene. Only 18.8% of the compiled elements the library tested are
active in any of K562, HepG2 or WTC11 at log2(RNA/DNA) ≥ 1.0. That is a weaker contradiction than a
CRISPRi negative on the named gene, and it is reported under its own assay name for exactly that
reason. It is also consistent with what the lentiMPRA lane already found: the reporter and the locus
are different questions.

*R6 note (2026-09-28): "episomal" above is wrong for this assay. The ENCODE4 lentiMPRA integrates
its reporter by lentivirus (Agarwal et al. 2025, Nature, https://www.nature.com/articles/s41586-024-08430-9). The reading still holds, because an integrated reporter also sits
outside the sequence's native locus. After R6 the lentiMPRA disagreements number 14,512, and 126
elements whose tiles split evenly are counted apart (see the R6 result below).*

**Is the measured 4.33% a fair sample?** Both comparisons go through `genomeos/compare.py`, standardised
on length, GC and distance to the nearest coding TSS, with the medians of both arms printed:

| comparison (chr21) | targets | controls | matched difference | p | length (t/c) | GC (t/c) | TSS (t/c) |
|---|---|---|---|---|---|---|---|
| measured vs unmeasured | 200 | 4,974 | +0.040 | 0.21 | 339 / 296 | 0.506 / 0.488 | 30.1 kb / 29.0 kb |
| disagrees vs agrees | 146 (15 unmatched) | 44 | **−0.255** | 1.00 | 337 / 342 | 0.499 / 0.545 | 31.5 kb / 32.2 kb |

| comparison (chr22) | targets | controls | matched difference | p | length (t/c) | GC (t/c) | TSS (t/c) |
|---|---|---|---|---|---|---|---|
| measured vs unmeasured | 421 | 9,191 | +0.082 | 0.009 | 338 / 290 | 0.556 / 0.549 | 17.3 kb / 13.5 kb |
| disagrees vs agrees | 310 (27 unmatched) | 103 | −0.047 | 0.87 | 333 / 341 | 0.549 / 0.568 | 18.4 kb / 13.8 kb |

The hit both arms are scored on is "the deletion moves the target by at least 0.2 log2", and `p` is
`compare.py`'s one-sided p for *targets above controls*, so a value near 1 is a difference in the other
direction rather than a null. The measured
elements are **slightly** more often ones the model predicts strongly (+0.040 on chr21, not significant;
+0.082 on chr22, p 0.009), and the imbalance columns say why to be careful: the measured arm is longer
(338 against 290 bp) and sits further from a coding TSS on chr22 (17.3 kb against 13.5 kb), which is the
assays' own design showing through. On chr21 the elements an assay **disagrees** with are 25 points
*less* often ones the model predicts strongly, and their GC is lower (0.499 against 0.545); on chr22 the
same difference is −0.047 and does not hold. One chromosome is not a finding, which is why both are
printed.

**A defect this section was built to avoid.** `evidence.collect` skips the compiled tree *before*
parsing it, for speed, so a census that forgets `compiled=True` reports zero experimental facts and the
zero reads as a finding rather than as a bug. The result therefore records `compiled_facts_read`
(960,094, from 24 programs) beside every count, `genomeos evidence --measured` prints it first of the
three denominators, and a test pins the distinction: a census that read nothing and a census that
measured nothing must never look the same. The compiled programs also carry the empty case in words — a
chromosome with no measurement says "this program states no experimental fact" rather than falling
silent, and every program's header states its eligible count before its raised count. This is the same
move the locus benchmark made in separating "outside the model's input" from "the model got it wrong",
and the general form is the lesson written the same day: *a negative search result is a statement about
the search, not about the world*.

**Cost.** The whole genome — 24 chromosomes read, matched, compared and recompiled, plus the 960,094-fact
evidence read — is **106 seconds** and no requests. `uv run python scripts/measured_layer.py
--write-programs`; `uv run genomeos evidence --measured` prints the census, and `--compiled` gives the
Evidence explorer the rows behind it.

**What this is not.** It is not a claim that 95.7% of the compiled elements are wrong; it is the
statement that **94.02% of them were never eligible to be asked**. The layer is thin because the
assays' footprint is thin, not because the matching is strict: 72.5% of everything that was eligible
was raised. The number to argue about is 5.98%, and 4.33% and 2.01% follow from it.

## The fourth assay measures bases, not elements: saturation mutagenesis covers 15 compiled elements and measures 4 (2026-09-17, later)

The section above wired three assays and said what it had left undone: saturation mutagenesis
(`genomeos/knowledge/satmut.py`) was not wired in and was the obvious fourth. It is now, and the
reason it was worth adding is the reason its census is so small. The other three ask whether an
ELEMENT does something. Kircher et al. 2019 (GSE126550) asked which BASES inside one matter —
44,658 single-base substitutions over **21 regulatory elements**, which is the entire assay. Twenty-one.

**Three denominators, in the order the layer always states them:**

| | | |
|---|---|---|
| **1. read** | compiled facts parsed at all (`compiled=True`) | **960,096** in 24 programs |
| **2. eligible** | compiled elements inside satmut's footprint, at one shared base | **15 of 440,377** |
| **3. raised** | elements a satmut experiment is a measurement *of*, under the overlap rule | **4** |

**The footprint, beside the three assays already wired.** This is the result of the lane, and the
lane stops here rather than being padded out:

| assay | measurements available | elements eligible | elements raised | share of the compiled genome raised |
|---|---|---|---|---|
| lentiMPRA | 53,986 elements | 21,133 | 17,869 | 4.058% |
| CRISPRi | 14,734 pairs | 3,908 | 1,505 | 0.342% |
| VISTA | 2,442 elements | 2,375 | 98 | 0.022% |
| **saturation mutagenesis** | **21 elements** | **15** | **4** | **0.0009%** |

satmut's footprint is **158 times smaller than VISTA's** and 1,409 times smaller than lentiMPRA's,
measured in eligible elements. Adding it moved the genome-wide eligible count from 26,316 to 26,324
(+8; seven of its fifteen were already eligible through another assay) and the raised count from
19,070 to 19,072 (+2; two of its four were already raised by lentiMPRA or VISTA). The experimental
layer grew from 19,282 facts to 19,284, and its share of the compiled genome from 2.0083% to 2.0085%.
**Two facts is the whole of it, and the honest reading is the eligible count, not the raised one.**

Its 21 experiments sit on 13 chromosomes (chr1 and chr10 three each; chr2, chr6, chr8 and chr11 two
each; one each on chr5, chr7, chr9, chr19, chr20, chr22 and chrX). **chr21 holds no satmut experiment at all**, and
chr22 holds one (GP1BA) whose interval contains no compiled element — so on the two chromosomes this
project usually reads first, the fourth assay's honest output is zero, of the never-looked kind.

**What "agrees" and "disagrees" mean for a base-level assay against an element-level prediction.**
The compiled claim is that deleting the element moves a gene. Saturation mutagenesis never deletes
the element and never measures that gene, so its verdicts are asymmetric by construction:

- **agrees** — at least one measured base inside the element is functional (some substitution at it is
  significant at the data portal's own defaults: p < 1e-5 with ≥ 10 barcodes). The element contains
  bases whose identity changes activity, which is measured support of the same weak kind lentiMPRA's
  "active" is.
- **it cannot disagree.** `disagrees["satmut"]` is **zero by construction, not by result**, and the
  census carries that sentence as a string (`SATMUT_CANNOT_DISAGREE`) so the zero can never be read as
  "never contradicted". An element every one of whose measured bases is inert is counted under its own
  name, `bases_measured_none_functional`. **An element whose bases mostly do not matter is not an
  element that does nothing.** Single-base substitution cannot see a function carried redundantly
  across a site; a null is depth-dependent besides (the functional share of these same 21 elements
  runs from 1% in FOXE1 to 77% in IRF4 with barcode depth); and none of it touches the deletion or
  the gene. A base-level null is evidence about those bases and about nothing larger. It keeps
  `experimental` as its evidence kind and is written into the measured block by name, exactly as a
  CRISPRi negative is.
- **`bases_not_measured`** is the fourth outcome and the never-looked one: the experiment met the
  overlap rule but covers none of the element's own bases. Measured-and-inert and never-looked are
  two zeros with two names, and `elements_in_the_footprint_never_raised` is a third.

It raises **no rule**, and not for CRISPRi's reason. A CRISPRi negative raises none because a rule
would state a relation the assay says is absent; satmut raises none because it names no relation at
all. Its experiments carry gene names given by their authors (SORT1, IRF4), and that name is curation,
not a measured target — three of the four matched elements happen to predict the gene their experiment
is named for, and this lane does not count that as agreement, because the name is not a measurement.

**The four elements, each one named, because a median over four elements is a worse summary than the
four rows.** Every one of them is covered end to end (100% of the element's bases measured):

| element | chrom | length | GC | TSS | experiment | overlap | bases measured | functional | strong | predicted target | other assays |
|---|---|---|---|---|---|---|---|---|---|---|---|
| EH38E1374646 | chr1 | 349 | 0.547 | 9.1 kb | SORT1.2 | 0.582 | 349 (100%) | 266 (76%) | 232 | PSRC1 −0.41 | lentiMPRA |
| EH38E2863563 | chr1 | 350 | 0.577 | 9.9 kb | IRF6 | 0.583 | 350 (100%) | 127 (36%) | 93 | IRF6 −0.54 | — |
| EH38E2046168 | chr2 | 350 | 0.634 | 69.7 kb | UC88 | 0.593 | 350 (100%) | 48 (14%) | 27 | TBR1 −0.12 | lentiMPRA, VISTA |
| EH38E2438393 | chr6 | 350 | 0.426 | 4.6 kb | IRF4 | 0.776 | 350 (100%) | 305 (87%) | 235 | IRF4 −0.22 | — |

1,399 bases measured, 746 functional (53.3%), 587 strongly so (42.0%). All four agree, none is inert,
and none is a disagreement — nor could any of them have been.

**The one element three assays and a base-level assay all reached, which is the whole distinction in a
single row.** EH38E2046168 is VISTA-positive in forebrain in two transgenic experiments, active in
lentiMPRA in WTC11 (log2 1.49, silent in K562 and HepG2), and **86% of its bases individually do not
matter** — 48 functional of 350, 27 strongly. An element can be a demonstrated in-vivo enhancer while
almost none of its single bases is load-bearing. Had the two questions been collapsed into one verdict,
this row would have read "14% functional, therefore mostly dead", which the embryo contradicts.

**The rule at three values, for satmut alone**, with the containments it refuses:

| reciprocal overlap | elements matched |
|---|---|
| ≥ 0.25 | 11 |
| **≥ 0.50 (the rule)** | **4** |
| ≥ 0.75 | 1 |
| *contained only, not upgraded* | *3* |

The fall from 11 to 4 is the same shape the other assays show: a satmut experiment is 187–601 bp and a
compiled element 293–350 bp, so the two are commensurable, and the eleven at 0.25 are mostly elements
sitting half in and half out of the tested interval. One element (chr9, FOXE1) misses the rule at
0.488 and is not raised, which is the rule doing its job rather than an unlucky rounding.

**Was the comparison run? No, and the refusal is in the result.** Both groups are described through
`genomeos/compare.py` — `input_presence` first, which says the claim is **bought** (the input, "bases
of this element measured one at a time", is present on 4 rows and missing on 440,373), and then
`imbalance`, which prints length, GC and median distance to a coding TSS for both arms per chromosome
along with the coverage of each. On chr1 the two targets are longer (349.5 against 297 bp), more GC
(0.562 against 0.489) and half as far from a coding TSS (9.5 kb against 18.0 kb) as the chromosome's
other 43,679 elements. The standardised difference is **not** computed and
`comparison_refused` says why in the result: four targets is below the lane's bar of twenty, and a
stratified difference over four elements is noise carrying a p-value. The footprint is the finding.

**What could not be assessed.** Whether these elements' bases behave the same way in their own
chromosome: the assay is episomal, so it lands at `REPORTER_CONFIDENCE` 0.75 with lentiMPRA and not
with CRISPRi, however fine its resolution. Whether an inert element exists: none of the four is inert,
so the `bases_measured_none_functional` branch is exercised by tests and not yet by data — the
element that would most likely have shown it (chr1, 68 bases measured, none functional) sits at
overlap 0.113 and is not raised. And nothing at all about the 99.9991% of compiled elements this
assay has never been pointed at.

*R6 note (2026-09-28): "the assay is episomal" was not checked against Kircher et al. 2019 for each
experiment. The conclusion depends only on the reporter being outside the locus, which holds whether
the construct is a plasmid or integrated.*

**Cost.** 24 chromosomes read, matched, compared and recompiled, plus the 960,096-fact evidence read:
**159 seconds** and no requests, against 106 before the fourth assay. The satmut table is 2.6 MB
gzipped, parsed once per process and cached. `uv run python scripts/measured_layer.py
--write-programs`; `uv run genomeos evidence --measured` prints the census with the fourth assay's
counters and the sentence that explains its zero.

## Are the project's own confidences calibrated? The stated level sits outside its own interval in 10 of 12 bands, and the ordering holds in 1 of the 4 reliability tables (2026-09-17)

Every compiled fact carries a confidence and the Evidence explorer counts them. Nobody had asked whether a fact stated at 0.4 is right about 40% of the time. The measured layer made the question askable, and the answer is reported per assay because *agreement* is a different event in each of the four.

**What the number is.** `round(min(PREDICTED_CAP, max(0.05, |log2 fold change|)), 2), as attribution/compile.py writes it on the element block and on the rule beside it; PREDICTED_CAP = 0.7, so no compiled element can state more than that however large the predicted effect`

**What a failure would look like, written down before the rates were computed.** The four patterns are named constants in `genomeos/attribution/confidence_calibration.py` and the bands are imported unchanged from `target_calibration`, where another lane fixed them on the same day for a different quantity. A band holds a verdict only at 30 measured elements or more, an assay gets a verdict only at 3 judged bands or more, and the level clears the bar at 70% of judged bands consistent.

- **`uninformative`** — the observed agreement rate does not rise with the stated confidence: the highest judged band's rate is at or below the lowest judged band's, or their 95% Wilson intervals overlap. If this is the pattern, the number on a compiled fact orders nothing, and a reader who preferred a 0.6 fact to a 0.2 one gained nothing by it.
- **`ordered_but_miscalibrated`** — the rate rises - the highest judged band is above the lowest and their intervals are disjoint - but the band's mean stated confidence lies outside the 95% Wilson interval of its observed rate in more than 30% of the judged bands. The number ranks and its level is wrong. This is the ordinary outcome for a scoring system read as a probability, and it is a useful result: a rank is worth having. The offset is reported as the median signed gap, and a level is a property of the population and not of the scorer - target_calibration's curve transferred while its level did not, because the base rate of the new population differed.
- **`calibrated`** — the rate rises AND the band's mean stated confidence lies inside the 95% Wilson interval of its observed rate in at least 70% of the judged bands. Then the number may be read as a probability, within the scope below and nowhere else.
- **`level_withheld_ordering_only`** — the assay's agreement is not the compiled claim's own event, so the ordering may be read and the level may not: a rate of 0.19 against a stated 0.25 is not a miscalibration when the rate being measured is 'a 200 bp sequence drives a reporter' and the claim is 'deleting the element in its own chromosome moves this gene'. The level verdict is withheld rather than computed and ignored.

**The answer, per assay and never pooled.** The rate rises with the stated confidence in `lentimpra:where_the_predicted_gene_was_tested`; it does not in `crispri:over_all_matched_elements`; the pre-registered rule refuses a verdict in `crispri:where_the_predicted_gene_was_tested`, `vista:where_the_predicted_gene_was_tested`, `satmut`.

| table | n | mean stated | observed | median offset | ECE | pattern |
|---|---|---|---|---|---|---|
| `crispri:where_the_predicted_gene_was_tested` | 128 | 0.347 | 0.7578 | **+0.4724** | 0.4108 | `refused` |
| `crispri:over_all_matched_elements` | 1,505 | 0.3268 | 0.0645 | **-0.2913** | 0.2623 | `uninformative` |
| `lentimpra:where_the_predicted_gene_was_tested` | 17,869 | 0.2786 | 0.1881 | **-0.1371** | 0.0905 | `level_withheld_ordering_only` |
| `vista:where_the_predicted_gene_was_tested` | 98 | 0.3262 | 0.7041 | **+0.5179** | 0.3779 | `refused` |

**The offsets have opposite signs, and that is the result of record.** `crispri:where_the_predicted_gene_was_tested` +0.4724, `crispri:over_all_matched_elements` -0.2913, `lentimpra:where_the_predicted_gene_was_tested` -0.1371, `vista:where_the_predicted_gene_was_tested` +0.5179. The same stated number, measured against two populations, is off in OPPOSITE DIRECTIONS. That is not an inconsistency in the tables; it is what a level is. target_calibration found the same thing for a different quantity - the curve transferred and the level did not, because the base rate of the new population differed - and it is the reason no single offset can be added to the compiler's formula to fix it. An offset quoted without the population it was measured on is not a number.

**Three denominators, in order.** Of **440,377** compiled elements that state a confidence, **19,072** are measured by any assay (**4.33%**) and **421,305** are not. 0 have no covariates and enter no comparison.

**The measured slice against the rest, on the covariates and on the confidence itself:**

| | measured | never measured | ratio |
|---|---|---|---|
| length (bp) | 340.0 | 293.0 | 1.16 |
| GC | 0.5014 | 0.4793 | 1.046 |
| distance to a coding TSS (bp) | 27826.0 | 22717.0 | 1.225 |
| **stated confidence** | 0.21 | 0.19 | 1.105 |

**Coverage by band — the selection this lane cannot undo, as a number:**

| band | measured | never measured | coverage |
|---|---|---|---|
| 0.1-0.25 | 11,227 | 281,434 | 3.84% |
| 0.25-0.5 | 4,879 | 93,399 | 4.96% |
| 0.5-0.75 | 2,966 | 46,472 | 6.00% |

### crispri

*silencing this element in its own chromosome, does the predicted gene fall?* — the screen perturbs the element where it lives and measures the very gene the deletion named, which is the compiled claim asked directly. A negative here is a strong contradiction and the rate may be read as a probability.

3,908 compiled elements in this assay's footprint and 1,505 measured; 128 of them are ones where the screen tested the very gene the deletion named.

**where_the_predicted_gene_was_tested** — n = 128, observed 0.7578, mean stated 0.347, expected calibration error 0.4108.
*(2026-09-28: superseded by R2, `1cff24a`. With significant increases and underpowered nulls read as their own
outcomes, 4,613 → 4,565 pairs measured as not regulated, n 128 → 120, observed 0.7578 → 0.8083, ECE 0.4108 →
0.4507. See "CRISPRi negatives read as outcomes" below; the figures here are kept as they were.)*

| band | compiled | in footprint | measured | coverage | mean stated | agrees | observed | 95% CI | inside | length | GC | TSS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1-0.25 | 292,661 | 2,142 | 62 | 0.02% | 0.1566 | 39 | **0.629** | 0.505–0.738 | no | 347.0 | 0.5215 | 22951.5 |
| 0.25-0.5 | 98,278 | 992 | 29 | 0.03% | 0.3486 | 25 | **0.8621** | 0.694–0.945 | no *(not judged)* | 348.0 | 0.4942 | 21183.0 |
| 0.5-0.75 | 49,438 | 774 | 37 | 0.07% | 0.6649 | 33 | **0.8919** | 0.753–0.957 | no | 349.0 | 0.5 | 12723.0 |

**Verdict: refused.** 2 bands hold at least 30 measured elements, below the 3 this lane fixed before computing anything. The table is printed and no verdict is drawn from it: a curve through fewer than three points is a line by construction. The bands that are populated and below the bar are named with the size of the shortfall, because a bar missed by one is still a bar and moving it after seeing the counts is how a pre-registration stops being one. Short of the bar: band 0.25-0.5 holds 29, 1 below 30. The gaps are still described: 0.1-0.25 +0.4724, 0.25-0.5 +0.5134, 0.5-0.75 +0.227, median +0.4724.

Is the band reading the confidence or the covariate? Bought or free first: A stated confidence on the compiled fact is **free** (128/128 of the measured arm, 438,872/438,872 of the unmeasured); A verdict from crispri over this element is **bought** (128/128 of the measured arm, 0/438,872 of the unmeasured).

Standardised on length, GC and distance to a coding TSS: 37 high-confidence elements against 91 low, 34 matched and 3 dropped for want of a control; raw +0.1886, matched **+0.1374** at one-sided p 0.06975.

**over_all_matched_elements** — n = 1,505, observed 0.0645, mean stated 0.3268, expected calibration error 0.2623.

| band | compiled | in footprint | measured | coverage | mean stated | agrees | observed | 95% CI | inside | length | GC | TSS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1-0.25 | 292,661 | 2,142 | 739 | 0.25% | 0.1587 | 39 | **0.0528** | 0.039–0.071 | no | 346.0 | 0.5059 | 21583.0 |
| 0.25-0.5 | 98,278 | 992 | 417 | 0.42% | 0.3512 | 25 | **0.06** | 0.041–0.087 | no | 346.0 | 0.5115 | 12146.0 |
| 0.5-0.75 | 49,438 | 774 | 349 | 0.71% | 0.6533 | 33 | **0.0946** | 0.068–0.130 | no | 347.0 | 0.5447 | 6227.0 |

**Verdict: `uninformative`.** Ordering: the point estimates rise but the intervals overlap — 0.1-0.25 at 0.0528 to 0.5-0.75 at 0.0946, intervals overlapping. Level: the stated confidence lies inside its band's interval in 0 of 3 judged bands (0%, bar 70%); median signed offset -0.2913, range -0.5587 to -0.106.

Is the band reading the confidence or the covariate? Bought or free first: A stated confidence on the compiled fact is **free** (1,505/1,505 of the measured arm, 438,872/438,872 of the unmeasured); A verdict from crispri over this element is **bought** (1,505/1,505 of the measured arm, 0/438,872 of the unmeasured).

Standardised on length, GC and distance to a coding TSS: 349 high-confidence elements against 1,156 low, 349 matched and 0 dropped for want of a control; raw +0.0392, matched **+0.0364** at one-sided p 0.034782.

### lentimpra

*out of its chromosome, does a 200 bp copy of this sequence drive a reporter?* — episomal: it measures the sequence and not the locus, and it never sees the predicted gene. A silence is a weak contradiction, and the base rate of 'active in a reporter' has no reason to equal the base rate of 'deleting this moves that gene'. Ordering only.

*R6 note (2026-09-28): this assay is an integrated lentiviral reporter (Agarwal et al. 2025, Nature, https://www.nature.com/articles/s41586-024-08430-9), not episomal; the rebuilt result states it that way. The reporter is still outside its native locus.*

21,133 compiled elements in this assay's footprint and 17,869 measured; 17,869 of them are tested by it, which is every one, so the two denominators coincide - this assay never asks about the predicted gene at all, and the wide rate is kept as a named key so the narrow one is never the only rate on the page.

**where_the_predicted_gene_was_tested** — n = 17,869, observed 0.1881, mean stated 0.2786, expected calibration error 0.0905.

| band | compiled | in footprint | measured | coverage | mean stated | agrees | observed | 95% CI | inside | length | GC | TSS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1-0.25 | 292,661 | 12,621 | 10,627 | 3.63% | 0.1584 | 1,618 | **0.1523** | 0.145–0.159 | yes | 334.0 | 0.4957 | 37713.0 |
| 0.25-0.5 | 98,278 | 5,377 | 4,570 | 4.65% | 0.3421 | 937 | **0.205** | 0.194–0.217 | no | 342.0 | 0.5029 | 23825.5 |
| 0.5-0.75 | 49,438 | 3,135 | 2,672 | 5.41% | 0.6483 | 806 | **0.3016** | 0.284–0.319 | no | 344.0 | 0.5171 | 14798.0 |

**Verdict: `level_withheld_ordering_only`.** Ordering: rises with the stated confidence — 0.1-0.25 at 0.1523 to 0.5-0.75 at 0.3016, intervals disjoint. Level: the stated confidence lies inside its band's interval in 1 of 3 judged bands (33%, bar 70%); median signed offset -0.1371, range -0.3467 to -0.0061.

Is the band reading the confidence or the covariate? Bought or free first: A stated confidence on the compiled fact is **free** (17,869/17,869 of the measured arm, 422,508/422,508 of the unmeasured); A verdict from lentimpra over this element is **bought** (17,869/17,869 of the measured arm, 0/422,508 of the unmeasured).

Standardised on length, GC and distance to a coding TSS: 2,672 high-confidence elements against 15,197 low, 2,672 matched and 0 dropped for want of a control; raw +0.1335, matched **+0.1293** at one-sided p 0.0.

### vista

*in a transgenic mouse embryo at e11.5, is this sequence an enhancer?* — a different organism, one developmental stage and a reporter construct. A negative is a negative in a mouse embryo, which is not a measurement of a human cell's transcription, and the assay's own positives are enriched by how its sequences were chosen. Ordering only.

2,375 compiled elements in this assay's footprint and 98 measured; 98 of them are tested by it, which is every one, so the two denominators coincide - this assay never asks about the predicted gene at all, and the wide rate is kept as a named key so the narrow one is never the only rate on the page.

**where_the_predicted_gene_was_tested** — n = 98, observed 0.7041, mean stated 0.3262, expected calibration error 0.3779.

| band | compiled | in footprint | measured | coverage | mean stated | agrees | observed | 95% CI | inside | length | GC | TSS |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.1-0.25 | 292,661 | 1,644 | 55 | 0.02% | 0.1676 | 39 | **0.7091** | 0.579–0.812 | no | 346.0 | 0.4711 | 63594.0 |
| 0.25-0.5 | 98,278 | 491 | 14 | 0.01% | 0.3393 | 12 | **0.8571** | 0.601–0.960 | no *(not judged)* | 347.5 | 0.5602 | 39688.5 |
| 0.5-0.75 | 49,438 | 240 | 29 | 0.06% | 0.6207 | 18 | **0.6207** | 0.440–0.773 | yes *(not judged)* | 347.0 | 0.5536 | 21559.0 |

**Verdict: refused.** 1 band holds at least 30 measured elements, below the 3 this lane fixed before computing anything. The table is printed and no verdict is drawn from it: a curve through fewer than three points is a line by construction. The bands that are populated and below the bar are named with the size of the shortfall, because a bar missed by one is still a bar and moving it after seeing the counts is how a pre-registration stops being one. Short of the bar: band 0.25-0.5 holds 14, 16 below 30; band 0.5-0.75 holds 29, 1 below 30. The gaps are still described: 0.1-0.25 +0.5415, 0.25-0.5 +0.5179, 0.5-0.75 +0.0, median +0.5179.

Is the band reading the confidence or the covariate? Bought or free first: A stated confidence on the compiled fact is **free** (98/98 of the measured arm, 440,279/440,279 of the unmeasured); A verdict from vista over this element is **bought** (98/98 of the measured arm, 0/440,279 of the unmeasured).

Standardised on length, GC and distance to a coding TSS: 29 high-confidence elements against 69 low, 27 matched and 2 dropped for want of a control; raw -0.1184, matched **-0.1395** at one-sided p 0.871033.

### satmut

*which single bases inside this element change a reporter's activity?* — it cannot disagree at all: substituting one base at a time never deletes the element and never measures the predicted gene, so its agreement rate is 1.0 by construction wherever any base is functional. A reliability table over it would measure nothing twice over, and at four matched elements it would not measure it. See measured.SATMUT_CANNOT_DISAGREE.

**No table.** 15 elements in the footprint, 4 measured, 4 agreeing. 4 matched elements, and the assay cannot disagree: zero by construction, not by result: saturation mutagenesis substitutes one base at a time in a reporter and never deletes the element or measures the predicted gene, so it can support the compiled claim and cannot contradict it. An element whose measured bases are all inert is counted under bases_measured_none_functional, which is evidence about those bases only. A reliability table needs an event that can come out false, and this one cannot, so no table is computed. The group is described in `the_measured_slice` with the rest.

**Why there is no pooled number.** One agreement rate over all four assays would be a rate whose event changes with the denominator. A CRISPRi negative on the named gene contradicts the compiled claim; a lentiMPRA silence says a 200 bp copy does not drive a reporter out of its chromosome; a VISTA negative says a mouse embryo did not stain at e11.5; saturation mutagenesis cannot disagree at all and so contributes agreements and never disagreements. Pooling them would also pool their sizes - lentiMPRA's 17,869 matched elements against CRISPRi's 128 tested pairs is 140 to 1 - so the pooled number would be the reporter's number wearing the perturbation's name. Every table here is per assay and the result holds no pooled rate.

**What this cannot say.** This is a reliability table over MEASURED compiled elements, and that is not a random sample of the compiled genome. 4.33% of compiled elements are measured at all (19,072 of 440,377) and only 5.98% were ever eligible, because an assay's footprint had to reach them. The measured slice is longer, more GC-rich and closer to a coding TSS than the rest - the assays chose candidate regulatory sequence, which is the same property the deletion scores highly - so the bands' rates are conditional on an element having been chosen for an assay. The table says how often the project's confidence is borne out ON MEASURED ELEMENTS. It says nothing directly about the 414,053 elements no assay ever covered, and the imbalance between the two arms is printed in the result rather than described here, so the size of the extrapolation is a number.

**What would falsify the transfer.** The transfer from measured elements to the compiled genome fails if any of these is observed. (a) A screen that tiles a region without choosing candidate sequence - an unbiased tiling rather than a cCRE list - finds the top band's rate far below what it is here, which is what selection on testability looks like from outside. (b) The agreement rate moves with the measured slice's own covariates once they are held fixed: if the standardised difference between the high and low confidence arms collapses toward zero when length, GC and distance to a coding TSS are matched, then the band was reading the covariate and not the confidence, and on unmeasured elements the covariate distribution is different. (c) A new assay extends the footprint into elements that look unlike the measured ones - far from a promoter, AT-rich, short - and the bands do not hold there. (d) The prevalence shift target_calibration already found: a population whose base rate of agreement differs from this one's needs its own intercept, so a band quoted without its population's rate beside it has stopped being a measurement of anything.

**Cost.** 24 chromosomes read, binned and compared in **70.6 seconds**, no requests and no model. `uv run python scripts/confidence_calibration.py`; `--markdown` prints this section back out of `data/results/confidence_calibration_genome.json`, which is how every figure above got here rather than being typed beside a computed one.

## 99.5% of the unknown space has never been measured, and the part that matters would fit in one library (2026-09-17)

Three lanes hit the same wall today from different directions: the deletion sweep's per-tier readings
did not survive standardisation, a calibration of those predictions failed on the prevalence of the
screens that measured it, and a motif model that transfers to VISTA cannot beat its own shuffle when
swept blind. The transfer lane put a number on the wall — 161 of the 882 blocks hold a measured
element. `scripts/unknown_coverage.py` puts it on every tier, for no requests, counting a block as
measured if **any** of lentiMPRA, VISTA or the CRISPRi benchmark overlaps it by a single base.

| tier | Mb | measured Mb | share of bases | blocks measured | lentiMPRA / VISTA / CRISPRi |
|---|---|---|---|---|---|
| regulatory | 345.2 | 3.698 | **1.07%** | 6,422 of 15,536 | 5,857 / 343 / 1,241 |
| fossil | 328.5 | 0.561 | 0.17% | 1,541 of 7,302 | 1,455 / 71 / 93 |
| structural | 229.7 | 0.012 | **0.005%** | 39 of 238 | 39 / 0 / 0 |
| neutral | 73.3 | 0.122 | 0.17% | 332 of 2,632 | 314 / 17 / 13 |
| real unknown, tolerant | 16.4 | 0.053 | 0.32% | 92 of 329 | 82 / 16 / 1 |
| real unknown, relaxed | 13.4 | 0.105 | 0.78% | 99 of 437 | 75 / 38 / 2 |
| real unknown, syntax | 0.39 | 0.002 | 0.54% | 8 of 69 | 7 / 1 / 0 |
| real unknown, recent | 0.24 | 0.001 | 0.33% | 3 of 29 | 3 / 0 / 0 |
| constrained unknown, copies | 1.42 | 0.001 | 0.09% | 6 of 216 | 5 / 0 / 1 |

**The whole unknown space is 1,008.8 Mb and 4.554 Mb of it has been measured: 0.45%.** The real
unknown is 30.6 Mb with 0.16 Mb measured, **0.52%**, in 202 of its 882 blocks. Every model this
project has built over that space — the deletion sweep, the calibration, the motif counts, the tier
readings — has been extrapolating from half a percent of it, and the corrections of the last two days
are what that extrapolation looks like when it is finally checked.

**One number is actionable, and it is the reason to state the rest.** Tiled at 300 bp:

- the **untouched part of the real unknown, 13.77 Mb, is 45,900 oligos** — one MPRA library, at a size
  that is routine;
- **all 30.6 Mb of the real unknown is 102,000 oligos** — one large library;
- the whole unknown space is **3.4 million oligos**, which is not a library but a programme.

So the asymmetry that matters is not between tiers, which this project can no longer rank, but between
the real unknown and everything else: **the sequence we most want measured is small enough to measure
in one experiment.** The median untouched block is 2.7 kb in the syntax case and 6.3 kb in the relaxed
case, so an oligo library tiles whole blocks rather than sampling them.

**What this section is not.** It is not a claim that anything in those 13.77 Mb is functional; it is
the statement that we do not know, that no assay this project holds has asked, and that asking is
affordable. Coverage is the precondition for a claim, not a claim.

## The matched table, computed correctly: no tier differs from the neutral tier, and the section below overstated in both directions (2026-09-17, later)

The section below reported two matched numbers — syntax against neutral at **+0.284, p 0.00014**, and
relaxed against neutral at **+0.025, p 0.014** — and both were wrong, from the same defect in how the
matching was implemented. The control arm was built by copying every control of a stratum into a list
once per target, so a stratum with many controls counted many times over, and the standard error was
computed from that inflated list. Large strata dominated the estimate and the p-values were
manufactured by an n that did not exist.

Written the standard way — each stratum's rate computed once, averaged over the targets that fall in
it, with the error taken on the number of matched targets — the whole table changes, and the
run that produced it takes 3.7 seconds instead of the thirteen minutes that first exposed the defect.

**Against the neutral tier, matched on length, GC and distance to the nearest coding TSS:**

| label | elements | raw | matched | p |
|---|---|---|---|---|
| real unknown, syntax | 40 | +0.283 | **+0.100** | 0.18 |
| real unknown, recent | 22 | +0.344 | +0.136 | 0.18 |
| constrained unknown, copies | 38 | +0.260 | +0.105 | 0.18 |
| fossil | 23,248 | +0.011 | -0.023 | 1.0 |
| real unknown, relaxed | 1,430 | -0.038 | -0.031 | 0.97 |
| real unknown, tolerant | 1,784 | -0.064 | -0.040 | 1.0 |
| regulatory | 144,977 | +0.169 | -0.050 | 1.0 |

**Not one of them differs from the neutral tier once standardised.** Every raw gap, in either
direction, is within noise of zero after length, GC and promoter distance are held fixed. The
regulatory tier's +0.169 — the largest and most obvious of the raw gaps — becomes -0.050.

**What does survive, and it is the one statement this table supports:** every element inside an
UNKNOWN block acts less than the genome's scored elements, whatever tier it sits in. Matched, the
deficits run from **-0.046 (recent) and -0.100 (syntax)** through **-0.225 (regulatory)**, **-0.302
(neutral)**, **-0.304 (fossil)**, **-0.337 (relaxed)**, **-0.357 (tolerant)** to **-0.433
(structural)**. The unknown space is quieter than the rest of the genome under deletion; how its
tiers rank inside it is not measurable this way.

**So both readings of 2026-09-16 and 2026-09-17 are withdrawn.** "The real unknown is the least active
tier" was measured against an unmatched control. "The relaxed case reads above neutral once matched"
was measured with a broken matcher. The third statement is the one to carry: **the tiers cannot be
ordered against each other by this instrument at all.** The 69 syntax blocks read +0.100 above neutral
and -0.100 below the genome's elements, on 40 elements, neither of them distinguishable from zero.

**The engineering lesson is worth more than the reading.** A matched comparison written by pooling
repeated controls inflates its own sample size and manufactures significance; the same code was also
quadratic, which is how it was caught — it ran for thirteen minutes and was stopped, twice, before the
arithmetic was written properly. A matched estimate that gets slower as the arms grow is doing
something other than standardising.

## The syntax observation dissolves under matching, the tiling run is cancelled, and yesterday's headline needs a correction (2026-09-17)

Assembling the pre-registered tiling run found its own reason not to spend the budget, before a single
request: **the arms do not overlap in the covariate that dominates the outcome.** On chr22 the syntax
windows sit a median **2.6 kb** from a coding transcription start; the relaxed arm sits at 51 kb and
the neutral arm at 62 kb. Genome-wide, over elements the sweep already scored, the medians are
**81 kb (syntax), 240 kb (relaxed), 418 kb (neutral)** — against **42 kb for the genome's scored
elements as a whole.** Matching arms that differ by an order of magnitude means matching on their thin
tails, so the cheaper question was asked first, over elements already scored and with no request
spent (`scripts/syntax_blocks_matched.py`).

**Against the genome's own elements, the 40 syntax elements are ordinary.**

| comparison | raw | matched on length, GC and TSS distance |
|---|---|---|
| syntax (40) against every other scored element | 0.575 vs 0.640, **-0.065** | **-0.033, p 0.67** |
| syntax (40) against the neutral tier | +0.283 | **+0.284, p 0.00014** |
| relaxed (1,430) against the neutral tier | -0.038 | **+0.025, p 0.014** |

So the observation that prompted the registration is real only against one yardstick. The 0.575 is
what an ordinary element does — the genome's scored elements move a gene 0.640 of the time — and it
looked remarkable only because the neutral tier was the comparison, and the neutral tier's elements
sit ten times further from a promoter than the genome's.

**The registered run is therefore not run, and that decision is part of the record.** Its neutral arm
is the yardstick that has just been shown to be the wrong one; it would most likely have returned
"success" against that arm and measured promoter distance more precisely. The pre-registration stays
in this document and in `attribution/syntax_tiling.PRE_REGISTRATION` as written, marked unrun, with
this as the reason. Cancelling a registered run because a cheaper analysis answered its question is
allowed; changing its bars after seeing data is not, and neither happened here.

**The correction to yesterday's reading.** The 882-block section says the real unknown is the least
active tier in the genome, from 0.248 per element against the neutral tier's 0.293. That comparison is
unmatched, and matching reverses part of it: **the relaxed case, which is 1,430 of the tier's 3,280
elements, reads +0.025 ABOVE the neutral tier once length, GC and TSS distance are held fixed**, not
below it. What survives matching is the larger and cruder statement: every unknown tier sits far below
the genome's scored elements — neutral -0.247 and relaxed -0.284, matched — so the unknown space is
much less active than the genome's elements, and the ordering *within* it is not established. The
tolerant and recent cases hold the remaining 1,850 elements and have not been matched; until they are,
"the least active tier" should be read as "far below the genome's elements, with the ordering against
neutral unsettled".

**The lesson, which is the third of its kind this week.** A tier average is not a control. The neutral
tier was built to be the background for *constraint*, and it is a background for nothing else: its
elements differ from the genome's in distance to a promoter by a factor of ten, and that is the
variable the outcome tracks. Every tier-against-tier number in this document is now suspect in the
same way, and the ones that matter should be re-read with the strata applied.

## Pre-registration: the 69 syntax blocks, tiled and deleted against two matched arms (2026-09-17)

Written before the instrument existed and before any window was scored; the machine-readable copy is
`PRE_REGISTRATION` in `attribution/syntax_tiling.py`, committed in the same state.

**What prompts it.** Reading the finished sweep over the 882 blocks of the real unknown made the tier
the least active in the genome per element (0.248 against 0.293 at the neutral tier, P 7e-9 below).
One slice read the other way: the 69 blocks of the `syntax` case — held across mammals *and*
constrained among people — at **0.575 of 40 elements**. Forty elements, chosen by where ENCODE
happened to call a cCRE inside those blocks, compared with a tier average rather than matched
sequence, and sliced after the reading. That is an observation, not a result, and this is the
instrument that can fail it.

**The design, fixed now.** Each block of each arm is tiled in non-overlapping **300 bp** windows from
its start; a block yielding more than **12** contributes 12 evenly spaced across it, so a 40 kb block
cannot outvote a 2 kb one. Each window is one deletion request, as an element is in the sweep. Three
arms:

- **syntax**, the 69 blocks;
- **relaxed**, the 437 blocks held across mammals but variable among people — this arm differs from
  the first in the human axis alone, which is the variable the case claims;
- **neutral**, the tier the whole reading was taken against.

Each syntax window is matched on its own chromosome to one relaxed and one neutral window with **GC
within 0.04** and **distance to the nearest coding TSS within 35% relative**; length is equal by
construction. **The arms are deliberately not matched on constraint**, because mammalian constraint
defines the tier and human constraint is what the relaxed arm isolates: matching on either would
remove the variable under test. Unmatched syntax windows are dropped from that comparison and counted.

**Outcome and bars.** A window moves a gene at |log2| ≥ 0.1 with a named gene, the sweep's threshold.
Primary: the difference in that share, syntax minus neutral and syntax minus relaxed, each with a
one-sided p and a 95% upper bound. **Success is at least +0.10 against both arms at p 0.01 or better.
Weak is 0.05 to 0.10 against both, or +0.10 against one arm only. Failure is below 0.05 against
neutral — and failure means the syntax case stops being carried in the roadmap as the sharpest
candidate.** Budget 2,400 requests; one futility look at half, no interim success look.

**Declared secondary.** Inside the syntax arm, windows overlapping an element the sweep already scored
should act about as often as those 40 did; windows overlapping none are the new information. If only
the overlapping windows act, the tiling has added nothing and the arm's rate is the old observation in
a new coat.

**What it cannot do, stated first.** This is the same model that produced the observation, so it
cannot be independent evidence about it. What it fixes is the sampling: 40 ENCODE-chosen elements
against a tier average become a few hundred windows fixed by the blocks themselves against matched
sequence. A positive would say the model reacts more to deleting these sequences than to deleting
matched sequence — not that the sequence is functional. Measured evidence on these blocks remains the
thing that would, and lentiMPRA's coverage of them is the next question after this one.

## E1's extension ran: a weak positive by its own bar, carried by one cell line, with its declared secondary passing (2026-09-17)

The 2,260 pairs were scored as the pre-registration below declares them — a fresh sample, the first
1,000 skipped and verified to be the same 1,000, one futility look at 1,130 and no interim success
look. 4,518 requests, 71 minutes, no look fired.

**The result, against the three outcomes named before it ran:**

| | units | controls | difference | one-sided p | upper 95% |
|---|---|---|---|---|---|
| **extension, 2,260 pairs** | 637/1,096 = 0.581 | 490/968 = 0.506 | **+0.075** | **0.00037** | 0.111 |
| first run, 1,000 pairs (`cbff5f6`) | 205/398 = 0.515 | 191/381 = 0.501 | +0.014 | 0.38 | 0.073 |

**+0.075 is the weak band: 0.05 to 0.10, "settles nothing on its own".** It is not the declared
success, which required 0.10 *in size* at p 0.01 — and the p of 0.0004 does not convert it into one,
because the bar was set on size for exactly this reason. Nor is it the declared failure. The
registration's own words stand: a weak positive.

**The declared secondary passed, and this is the first prediction in this series to be met rather than
missed.** Agreement should be higher where DAP-G fine-maps the variant than where it does not. In the
first run it failed, -0.152 on 57 fine-mapped pairs at power too low to mean much. Here, at three
times that power: **fine-mapped +0.198 (132 units, 90 controls, p 0.0018) against +0.057 for the rest
(p 0.008)**. The split E2 and E3 found in the model, and the two-instruments reading found in two
measurements, appears again in the endpoint built to be independent of the first.

**And one cell line carries it, exactly as the two-instruments reading was carried by one cell line.**
GM12878 reads **+0.089 on 1,895 pairs (p 0.0001)**; Jurkat reads **-0.001 on 365 pairs**. So the
ENCODE objection to E2 and E3 is **narrowed, not closed**: in lymphoblastoid cells the model tracks a
reporter assay's measured direction weakly and detectably, and outside them this endpoint has 365
pairs in one line and nothing to say.

**Two things stated because they are uncomfortable.**
- **The extension's pairs are weaker by construction** — the order was by measured effect size, so
  these are the 2,260 after the strongest 1,000 — and yet they give the larger difference (+0.075
  against +0.014). That is the opposite of what the ordering predicts. The two samples' intervals are
  not inconsistent (the first run's upper bound was 0.073), so the honest reading is that the first
  1,000 were an unlucky draw rather than that weak effects are easier to detect; it is recorded as a
  surprise, not explained away.
- **The pooled figure over all 3,260 pairs is +0.059** (842/1,494 against 681/1,349). The
  pre-registration allows that to be described and forbids it as a test, because the decision to run
  the second sample was taken after seeing the first. It is written here once, as a description.

**What this leaves.** The executor claim's external endpoint is no longer silent: it is weakly
positive, significantly so, in one cell type, with its causality secondary met. E2's +0.268 remains far
outside this endpoint's interval, so whatever the model does at a fine-mapped eQTL it does not do at
that size against a reporter assay. The next thing that would move this is not more pairs of the same
kind — 3,260 is the database — but a second cell type with enough fine-mapped variants to test, and
MPRAVarDB does not currently hold one.

## Pre-registration: E1's remaining 2,260 pairs, read as a replication and not as a continuation (2026-09-16)

Written before any request of the extension and before any of its pairs were scored. The run it
extends is `cbff5f6`: 1,000 of E1W's 3,260 pairs, agreement 0.515 against matched nulls at 0.501,
**+0.014, one-sided p 0.38, upper 95% bound 0.073**, which its own pre-registration called a result
that decides nothing.

**Why this needs a document at all.** The obvious move -- spend the remaining pairs and read all 3,260
together -- is optional stopping run backwards: the decision to continue was taken *because* the first
look came out at p 0.38, so a pooled p-value from that sequence does not mean what a p-value means.
There are two honest routes and this takes the second: abandon the endpoint, or treat the unspent
pairs as a fresh sample with its own pre-registered reading, whose result stands whether or not it
agrees with the first.

**The sample.** All 2,260 pairs of E1W not scored in `cbff5f6`, in the order the original
pre-registration fixed (strongest measured MPRA effect first), with no further selection. They are
weaker by construction than the first 1,000, since the order was by effect size: this reading is
therefore of *weaker measured effects*, and that is a difference from the first sample, not a flaw in
it, but it must be stated wherever the two are put side by side.

**The primary reading**, identical in form to the first: the share of units whose model read-out moves
in the direction the reporter measured, minus the same share for the matched nulls drawn from the same
library, with its one-sided p and its 95% upper bound. Success is a difference of **0.10 or more at
one-sided p 0.01 or better**; **0.05 to 0.10** is a weak positive that settles nothing on its own;
**below 0.05, with an upper bound under 0.10**, is the endpoint failing at a size that would matter,
and is to be reported as the ENCODE objection standing unanswered rather than as an absence of
evidence.

**The secondary**, also carried over: agreement should be higher where DAP-G fine-maps the variant
than where it does not. It failed in the first run (-0.152 on 57 fine-mapped pairs) at power too low
to mean much; it is declared again here, and a second failure is a failure of the prediction, not of
the pairs.

**One look, no more.** At 1,130 pairs, for futility only: if the difference is at or below 0 with an
upper bound under 0.05, the run stops and is reported as stopped for futility. There is no interim
success look, so nothing can be gained by watching it.

**What no outcome here can do.** It cannot make E2 and E3 independent of the model -- only a
measurement that does not pass through AlphaGenome can, and the one this project has (two instruments,
2026-09-15) is carried by a single cell line. A success would say the model's agreement with a
reporter assay survives at the size the executor claim needs; a failure says the model does not track
a reporter assay's direction at that size, which is compatible both with the executor claim being
wrong and with a plasmid fragment being a poor proxy for a gene in its chromosome. The pooled figure
over all 3,260 pairs may be reported as a description, never as a test.

## The 882 blocks of the real unknown, read against the completed sweep (2026-09-16)

Milestone 1.3 asks for the constrained-unknown blocks attributed to a gene and a tissue. The sweep
finished today, so every ENCODE element inside a node has a deletion answer on every chromosome, and
the question is arithmetic: `scripts/constrained_unknown_targets.py`, no model request, 7 minutes over
tables on disk. The target is the organiser's real unknown — the constrained-unknown tier with the
copies removed, **882 blocks**, the number area I has been carrying as its remaining work.

**Per element, the real unknown is the least active tier in the genome, not the most.** The rate that
does not depend on how long a block is:

| tier | elements inside its blocks | move a gene | rate |
|---|---|---|---|
| regulatory | 145,002 | 66,915 | **0.462** |
| fossil | 23,487 | 7,122 | 0.303 |
| neutral | 5,706 | 1,670 | 0.293 |
| **real unknown (882 blocks)** | **3,280** | **812** | **0.248** |

Against the neutral tier that is z = -5.7 (P 7e-9); against fossil, -6.9; against regulatory, -24.6.
An element sitting in the sequence this project calls its sharpest attribution target moves a gene
**less often** than one sitting in sequence the budget calls neutral. The blocks are denser in elements
than neutral blocks (3.7 per block against 2.2), which is why the per-block rate goes the other way
(+4.6 points inside length deciles) — and that per-block excess is not consistent, running from -4.3 to
+6.9 points across the ten deciles. **Density of cCREs is not evidence of function; the per-element
rate is the one to quote, and it is negative.**

**331 of the 882 carry at least one element that moves a gene**, and they are leads rather than
findings: the same model names a target at 87% of matched random windows, and these rates sit below
the tiers they are read against.

**One slice reads the other way, and it is 40 elements.** Split by the case the two axes make, the
**syntax** case — constrained across mammals *and* among people, the 69 blocks named as the sharpest
candidates — has 40 elements inside 27 of its blocks, of which 23 move a gene: **0.575, 95% CI 0.42 to
0.73**, above the neutral tier at P 0.0002 and indistinguishable from the regulatory tier (P 0.10).
The `recent` case reads 0.636 on 22 elements. Both are small, both are post hoc in the sense that the
cases were defined before this reading but this slicing was not pre-registered, and the honest summary
is: **the tier as a whole fails, and the 69 blocks it was narrowed to do not — on 40 elements.** That
is the number to grow, and the pre-registration to write, before anything is claimed for them.

## What scoring every element bought the 98%, replicated on chr22 and read genome-wide (2026-09-16)

The reading of 2026-09-13 was decided on chromosome 21, the only chromosome then fully scored. With
the sweep complete, chr22 has its own closure table and can be asked the same questions, and the 69
syntax candidates can be read against every chromosome at once. No model request: 95 s of arithmetic
over tables already on disk.

**The falsifiable test replicates, and it is a negative on both chromosomes.** The closure recomputed
with the elements over UNKNOWN blocks dropped from every gene's input barely moves — chr21 0.407 to
0.405, chr22 0.321 to 0.302 — and both sit inside random removals of the same size per gene (50 draws:
0.86 of chr21's draws at or below the observed, 0.14 of chr22's). Recomputed with *only* those
elements it is **0.287 at p 0.26 on chr21 and 0.257 at p 0.42 on chr22**, against shuffled inputs at
0.25: no better than chance on either chromosome, and no better than random subsets of the same size
(0.12 and 0.06 of draws at or below). Both runs reproduce their committed closure first, so this is a
test of the elements and not of the arithmetic.

**Elements over unknown blocks act less, not more, and that replicates too.** Within strata of length,
GC and distance to the nearest coding TSS, an element over an UNKNOWN block moves a gene **0.185 less
often on chr21 and 0.116 less often on chr22** (both p 0.0005 by permutation within strata). The
unknown space is gene-poor sequence by construction, and the strata are there to take that out; what
is left is still negative.

**The 69 syntax candidates, at last read against a complete genome.** The chr21 run could reach 2 of
them; now all 69 sit on fully scored chromosomes:

| what the scoring gives a candidate | blocks |
|---|---|
| nothing — no element over it moves a gene | **44** |
| a target and a magnitude, no cell line at the bar | 17 |
| a target, a cell and a magnitude | **8** |

**25 of 69 gain a named target**, 11 of them a coding gene — GJD3 at chr17:40,365,372-40,367,951 with
log2 -1.39 in memory T cells, ZNF503 at chr10:75,412,413-75,428,875 in monocytes, GJD2 on chr15 in
endodermal cells. Nine of the 25 were read as "unexplained" before the deletion was scored, which is
the one category where the sweep changed a block's reading rather than confirming it.

**What this does not say.** A named target is a prediction from one model, and the same model names a
target at 87% of matched random windows; the number that matters is not 25 of 69 but whether these
differ from matched controls, which the sweep's own fold answers for elements at large (+2.88 points) *(Re-audited 2026-09-27, NODES-READER-WRITER.md: +2.90, 95% CI +2.03 to +3.81, on a control matched on boundary count only; +1.2 to +6.6 across four defensible baselines; and **+5.89, CI +3.18 to +8.46, on 661 measured CRISPRi pairs** — the direction stands on measurement, the size quoted here was never identified.)*
and which nobody has answered for these 69. The two-thirds that gain nothing are the honest headline:
**scoring every element of the genome leaves 44 of 69 candidate blocks exactly where they were.**

## Many human genomes, every chromosome: the panel closes at 24 of 24 (2026-09-16)

chr2 was the last one, read in 3,308 s with every control passing, and it moved nothing: coding exons
at 0.344 of the matched rate, genes core 0.603, the constrained-unknown tier at 0.880. The panel has
now read **24 of 24 chromosomes**, about 267 GB of Cactus alignment of 90 human assemblies streamed
and never stored, no model call anywhere in it.

**1,925,587 storage units** over 385.3 Mb, 3,442,635 of their events matched in gnomAD, **187,966
executor-ready** (184,604 with a fine-mapped GTEx eQTL, 8,611 with an MPRA allele pair), 47,157
hypervariable of which 45,878 sit on a tandem repeat. By tier: regulatory 866,918, fossil 769,793,
neutral 177,771, constrained unknown 74,109, structural 36,996.

**Eight claims are genome-wide, holding on all 23 pooled chromosomes** (chrY is read and pooled with
nothing): coding exons below 0.6 of the matched rate (median 0.379), coding units fixed above matched
(+26.2 points), callable genes core (59.0%), genes core above their own matched windows (+49.0),
the coding-neutral gap surviving the timing match (90.4%), the storage share of placed background
units (48.4%), and overdispersion per kilobase (9.65).

**Five are chromosome-specific, and the rollup names the chromosomes** rather than pooling them away:
the neutral tier within a quarter of its matched rate fails on chr8, chr17, chr20 and chrX; no unknown
tier above matched fails on five; unit classes even over replication timing fails on five; Gnocchi's
constrained kilobases being early-replicating holds on 14 of 23 and is weak wherever it holds; and
Gnocchi's silent panel-depleted kilobases being copies now fails on **14 of 23**, which is the claim
that looked strongest on the two acrocentrics and is the one this sweep destroyed.

**What closing the panel settles, and what it does not.** It settles the size of the storage
catalogue and the executor-ready set, which is what the rest of area I reads from. It does not settle
the constrained-unknown depletion by itself — that rests on the tier baseline below, and the honest
form of it is 0.829 of the neutral tier, on 19 of 23 chromosomes. And it leaves the panel's own
weakest point untouched: a repeat-aware reading of long VNTRs, which the instrument still does not do.

## The matched background is not a neutral baseline, and that is why the controls fail (2026-09-16)

Seven of the 23 chromosomes the panel has read fail one of its six controls, and on all six autosomes
and X the failing check is the same one: the neutral tier's fixed rate missing the GC- and
replication-timing-matched background it is measured against. Six special chromosomes would be one
explanation. `scripts/panel_tier_baseline.py` asks the committed results whether it is instead one
thing happening everywhere, and it is.

| tier against the matched background | above 1 on | one-sided p | median ratio |
|---|---|---|---|
| regulatory | **23 of 23** | 1e-7 | 1.077 |
| fossil | **22 of 23** | 3e-6 | 1.062 |
| neutral | **21 of 23** | 3e-5 | 1.105 |
| constrained unknown | 6 of 23 | 0.99 | 0.916 |
| coding exons | 0 of 23 | 1.0 | 0.379 |

**Every non-coding tier fixes less than its own control.** The background windows are matched on GC
and on replication timing, and on those two axes they match; what they are not is neutral sequence.
They carry coding exons, conserved elements and everything else that falls in a window of the right
composition, so they fix more than the unknown sequence they are the control for, by about 6 to 11%.
A tier's distance from 1 therefore contains that offset, and the six control failures are the tail of
it rather than six chromosomes with something wrong: the deviation shrinks as the tier gets bigger
(median 0.163 on the eleven chromosomes with the least neutral sequence, 0.083 on the eleven with the
most, Spearman -0.33), so the tier crosses the 0.05 band wherever it is both offset and small.

**What survives, and it is the reading the panel lane had already arrived at.** Against the neutral
tier of the same chromosome, rather than against the background:

- **constrained unknown 0.829, below the neutral tier on 19 of 23** (p 0.0013) — the only unknown tier
  that is depleted, and it stays depleted under either baseline (below the background on 17 of 23,
  p 0.017; on the 17 chromosomes that pass every control, 14 of 17 either way, p 0.0064).
- **fossil 0.977 (17 of 23) and regulatory 0.994 (12 of 23)** — indistinguishable from neutral, which
  against the background both appeared to exceed.
- **coding exons 0.351, below it on 23 of 23** — the positive control, unchanged by the choice.

So the offset changes no conclusion about the constrained-unknown tier and removes two that were
never claimed but were visible in the table: fossil and regulatory sequence is not *more* variable
than neutral, it only looked so beside a background that is not neutral. The honest form of every
tier statement is against the neutral tier, with the background offset stated; quoting "x% above the
matched rate" for a non-coding tier carries a systematic +6 to +11 points that belongs to the control,
not to the tier.

**What it does not say.** It does not say the matched background is badly built for what it was built
for -- matching GC and timing is what makes the coding ratio (0.387) meaningful, and coding is the
one tier the offset does not touch, because a real depletion of that size swamps it. It does not
measure what the background carries; that is the next step, and the candidates are the obvious ones
(coding bases, conserved elements, segmental duplication). And it rests on 22 chromosomes of one
panel. chr2 landed after this was first written and moved nothing: every figure above is the
23-chromosome reading with it in.

> **The next step was taken the following day (2026-09-17, the section below), and two of those three
> candidates were wrong — one of them backwards.** Coding bases account for **none** of the offset in
> any tier (0.6%, 1.3%, 0.3%), and I could have ruled them out here without measuring anything: the
> panel cuts canonical CDS out of the background by construction, so what is left is 1.0% coding. I
> named a candidate the construction had already excluded, which is this week's recurring mistake —
> read how a thing is built before theorising about what it contains. Segmental duplication is the
> worse one: I listed it as an ingredient of the offset and it **masks about a third of it**, because
> the background's duplicated kilobases are noisy (9.94 recurrent events per kb against 6.44
> elsewhere), so removing them *widens* the gap by 29 to 34%. Only conserved elements survived, and
> they are the largest single contributor at 31%, 25% and 18%. Everything measured together accounts
> for 41%, 56% and 25%, leaving **44% to 75% unexplained**. The finding above stands untouched; the
> mechanism I guessed for it did not.

## What the matched background is made of, and why composition is at most half of its offset (2026-09-17)

The section above established that the panel's GC- and replication-timing-matched background is not a
neutral baseline and named three candidates for what it carries: coding bases, conserved elements,
segmental duplication. It did not measure them. `genomeos/attribution/panel_background.py` and
`scripts/panel_background_composition.py` do, on 24 chromosomes, from local data only.

The background is rebuilt exactly as `human_panel.build_background` builds it -- the kilobases of the
hg38 grid inside the panel's alignment blocks, with canonical CDS bases and the events inside them
cut out, binned by GC and by replication timing. The rebuild is checked against the committed result
before anything is read from it: on all 24 chromosomes every rebuilt `ratio_gc_rt` matches the
committed one to within **0.0005**. Every kilobase is then annotated from GENCODE 50, phastCons
100-way, RepeatMasker, genomicSuperDups, ENCODE cCREs and the panel's own replication timing, and the
offset is decomposed by rebuilding the background with each candidate taken out of it and recomputing
the tier ratios through the panel's own code.

**What the background holds.** Medians across the 23 chromosomes the panel's medians rest on; the
last five columns are shares of the group's bases. "Coverage" is the share of the panel's assemblies
informing the group, which is effort and not composition, and is reported beside the covariates
because it has to be held separately from them.

| group | bases | unit | GC | median distance to a coding TSS | coverage | recurring/kb | coding (any transcript) | exons (any gene) | conserved elements | segmental duplication | repeat |
|---|---|---|---|---|---|---|---|---|---|---|---|
| **the background** | 2,808 Mb | kilobase | 0.399 | 87 kb | 1.000 | 6.63 | 1.0% | 7.0% | 5.3% | 4.7% | 51.1% |
| the background outside every block | 2,006 Mb | kilobase | 0.402 | 76 kb | 1.000 | 6.39 | 1.3% | 8.6% | 5.9% | 3.9% | 49.0% |
| coding exons (the positive control) | 34 Mb | 124 bp | 0.470 | 16 kb | 0.995 | 2.93 | 100% | 100% | 75.5% | 7.1% | 3.4% |
| structural | 198 Mb | 72.7 kb | 0.418 | 213 kb | 0.740 | 10.41 | 0.0% | 0.0% | 0.2% | 4.0% | 82.4% |
| fossil | 316 Mb | 19.4 kb | 0.390 | 56 kb | 0.996 | 7.07 | 0.0% | 0.2% | 2.8% | 7.0% | 62.4% |
| regulatory | 345 Mb | 14.3 kb | 0.440 | 24 kb | 0.998 | 7.20 | 0.0% | 0.4% | 3.9% | 2.8% | 55.9% |
| constrained unknown | 32 Mb | 7.7 kb | 0.390 | 80 kb | 0.999 | 6.26 | 0.0% | 0.4% | 8.4% | 4.0% | 44.5% |
| neutral | 72 Mb | 7.9 kb | 0.405 | 58 kb | 0.987 | 7.40 | 0.0% | 0.4% | 3.1% | 10.2% | 45.4% |

Every cell above was assessed on every chromosome: all 24 carry GENCODE, phastCons, RepeatMasker,
genomicSuperDups and cCREs, so no share in this table is a zero that means "never looked". Sampled
phyloP (300 kilobases per chromosome, 7,200 in all, read by range and never downloaded) puts the
background at **3.5% constrained bases** (standard error 0.4 points over kilobases) against the
neutral tier's 1.52%, fossil's 1.40%, regulatory's 2.41% and constrained unknown's 5.23%, the last
four taken from the budget's block-by-block phyloP with the blocks it never measured counted, not
zeroed (fossil 6,756 of 7,302; regulatory 15,155 of 15,536; neutral 2,431 of 2,632).

So the background is a little coding, a little conserved, half repeat, and 71% of it lies outside
every unknown block. Within the background, the kilobases that hold a conserved element are the
quietest thing in it -- 5.33 recurring events per kilobase against 6.73 where there is none -- and
cCREs (6.18 against 6.73), exons (6.37 against 6.66) and promoters (6.37 against 6.62) are quieter
too. Segmental duplication and satellite go the other way, hard: 9.94 against 6.44 and 10.70 against
6.57.

**What each covariate accounts for.** Each row is the tier's median ratio when that thing is taken
out of the background, and what that move is as a share of the tier's distance from 1.

| taken out of the background | fossil (1.062) | regulatory (1.077) | neutral (1.105) |
|---|---|---|---|
| CDS of any transcript | 1.062, accounts for none | 1.076, accounts for none | 1.104, accounts for none |
| exons of any gene | 1.054, **13%** | 1.052, **33%** | 1.093, **11%** |
| conserved elements | 1.043, **31%** | 1.058, **25%** | 1.086, **18%** |
| within 2 kb of a coding TSS | 1.062, accounts for none | 1.073, 6% | 1.102, 3% |
| all three of those together | 1.037, **41%** | 1.034, **56%** | 1.079, **25%** |
| segmental duplication | 1.080, *widens by 29%* | 1.104, *widens by 34%* | 1.139, *widens by 33%* |
| low-coverage kilobases (not a covariate) | 1.075, *widens by 15%* | 1.083, *widens by 12%* | 1.123, *widens by 31%* |

Read down the table: the coding bases the earlier section led with account for **none** of the offset
through other transcripts' CDS and 11% to 33% through exons at large; conserved elements account for
18% to 31%; promoter proximity for 0% to 6%. Taking all three out together closes **25% (neutral),
41% (fossil) and 56% (regulatory)** of the offset and leaves **44% to 75% of it standing**. Segmental
duplication -- the third candidate the earlier section named -- moves it the *other way*: the
background's duplicated kilobases are noisy, so they were masking about a third of the offset, and a
duplication-free control reads further from the tiers, not nearer. The three exclusions are assessed
on 23 of 23 chromosomes; the coverage row on 22, because on chrX the cut leaves under half the
background's kilobases and the arm reports that rather than a ratio built from the scraps.

**Coverage, held separately.** A covariate describes the sequence; coverage describes what was done
to it, and standardising on the first while the arms differ in the second manufactures a difference
out of effort. This lane has one effort axis and it is named as a constant: `COVERAGE_MEASURE`, the
share of the panel's 89 assemblies informing a window, since a recurring event needs two assemblies
carrying the minority state and a window half the panel reached cannot show what a whole-panel window
can. It is conditioned two ways. First by rebuilding the background without the kilobases below 0.95
(the row above, 22 chromosomes). Second inside the two-group comparison: 100-base windows, the hit
being "does this window carry a recurring event", each tier against the panel's own control -- every
non-coding window -- through `compare.standardised`, stratified three ways.

| tier | windows | GC + timing | + distance to a coding TSS | + coverage | median TSS, target vs control |
|---|---|---|---|---|---|
| fossil | 265,454 | +0.0231 | +0.0223 | **+0.0198** | 248 kb vs 89 kb |
| regulatory | 317,835 | +0.0184 | +0.0187 | **+0.0206** | 55 kb vs 89 kb |
| neutral | 66,203 | +0.0150 | +0.0152 | **+0.0195** | 291 kb vs 89 kb |
| constrained unknown | 29,696 | -0.0182 | -0.0142 | **-0.0114** | 262 kb vs 89 kb |
| structural | 20,077 | +0.0154 | +0.0149 | **+0.0338** | 545 kb vs 89 kb |

No target window was dropped for want of a control in any of the three, on any chromosome. The
difference does not die on coverage; on three of the five tiers it grows. That is what licenses the
stronger of the two sentences the result file carries as data rather than as prose, and the file
carries the one that is available and not the one that sounds better:
**"nothing I measured explains it, including coverage"**. Had the coverage stratum not been run, the
only sentence available would have been "the covariates I measured do not explain it", which is a
different and weaker claim; the two are separate string constants so that neither can be upgraded to
the other by paraphrase.

**What this changes.** Nothing about any tier, and one thing about the explanation offered for the
offset. The section above wrote that the background windows "carry coding exons, conserved elements
and everything else that falls in a window of the right composition, so they fix more than the
unknown sequence they are the control for". That is the right direction and the wrong magnitude, and
it leads with the wrong term: canonical coding bases are already cut out of the background by
construction, what remains is 1.0% coding and 7.0% exonic, and the three named candidates account for
between a quarter and a half of the offset with duplication pulling the other way. The honest form is
that composition explains part of the offset and not most of it. Every tier statement still reads
against the neutral tier, the +6 to +11 points still belong to the control rather than to the tier,
and the constrained-unknown depletion is untouched: it is the one tier below the background in the
ratio (0.916, 6 of 23 above 1) and the one tier below the control in the window comparison (-0.0114
with coverage held), as it was.

**Coverage of the measurement itself.** 2,870,132 kilobases lie inside the panel's alignment blocks
over the 24 chromosomes; 2,856,365 were measured. The 13,767 that were not are two named categories:
10,503 held under half a kilobase once the canonical-CDS cut was applied, which is
`build_background`'s own rule, and 3,264 have no GC value in the store. A further 35,876 were
measured but carry no replication timing and are matched on GC alone, exactly as the panel matches
them. No track was absent on any chromosome and GENCODE was read for all 24, so every "0%" in the
composition table is a measurement. The window comparison subsamples long chromosomes by a fixed step
between 2 and 20, which is recorded per chromosome. phyloP was sampled, not read whole: 300
kilobases of each chromosome's background, with the unsampled kilobases counted.

**What is left undone.** The largest single part of the background by bases is not any of the
covariates above: 2,006 Mb of 2,808 Mb lies outside every unknown block, and it reads 6.39 recurring
events per kilobase against the tiers' 7.07 to 7.40. What that sequence *is* -- gene bodies and
introns against intergenic sequence the budget never called a block -- was not measured here, and it
is the next covariate rather than a zero in the table above. It is also the last one worth trying
before concluding that the offset is a property of how the budget cuts blocks rather than of what the
background contains. Written to `data/results/panel_background_composition.json`.

## The 71% of the background nobody had looked at: it is gene bodies, and the cut is the class (2026-09-17, later)

The section above measured what the panel's matched background holds, found that composition explains
at most half of the offset every non-coding tier reads above it, and named its own leftover: **2,006
Mb of the 2,808 Mb background -- 71% -- lies outside every unknown block**, reading 6.39 recurring
events per kilobase against the tiers' 7.07 to 7.40, and nobody had measured what that sequence is.
`genomeos/attribution/panel_leftover.py` and `scripts/panel_leftover_classes.py` measure it on 24
chromosomes from local data only. It was the last covariate the lane had left before concluding that
the offset is a property of how the budget cuts blocks rather than of what the background contains,
and the answer is that the question, put that way, cannot be answered -- for a reason worth more than
the answer would have been.

**What was eligible.** Over the 23 chromosomes the panel's medians rest on: 3,031 Mb of sequence,
2,846,945 kilobases inside the panel's alignment blocks, 28,468,458 100-base cells inside those, five
tiers and eight declared classes. Counting chrY as well -- measured and reported, excluded from every
median as everywhere in this panel -- 2,870,132 kilobases lie inside the alignment blocks and
2,856,365 were measured, which are the section above's two numbers to the kilobase, because the
background is rebuilt through the same `build_background`. The leftover comes back at **2,005,628,435
bases**, the figure that section committed, to the base. The 13,767 kilobases not measured are its two
named categories: 10,503 hold under half a kilobase after the canonical-CDS cut and 3,264 have no GC
in the store. GENCODE was read on all 24. Six of the eight classes have bases somewhere; the other two
are measured zeroes and say so rather than being absent.

**How a block is cut, read rather than inferred.** `genome/unknown.unknown_blocks` has one rule:

    cursor = 0
    for every gene in the annotation, by start:
        if gene.start - cursor >= min_size:  cut a block from cursor to gene.start
        cursor = max(cursor, gene.end)
    if length - cursor >= min_size:          cut a last block

So a base is never cut into a block for exactly one of three reasons, and they are exhaustive: it is
inside the span of an annotated gene; it is in an inter-gene gap shorter than `min_size` (1,000
bases); or it is in a tail shorter than that after the last gene. The first reason is wider than it
looks, because `Annotation` loads `gene`, `ncRNA_gene` and `pseudogene` alike, so lncRNA bodies and
pseudogene bodies are skipped along with the coding exons the rule is aimed at. `replay_the_cut` walks
that loop keeping the bookkeeping it throws away, and before anything is read from it the blocks it
recovers are checked against the committed budget: **26,806 blocks for 26,806, identical on 24 of 24
chromosomes**, none in the replay only and none in the budget only.

**What the 2,006 Mb is.** Over the 23 chromosomes. The classes are made disjoint by subtraction in
the order shown, and the bases a nesting sent to a higher-priority class are counted rather than left
to the reader: 175.8 Mb of non-coding gene body and 29.6 Mb of pseudogene body lie inside a coding
gene's span and are counted there. No base is in a block and in a class.

| what it is | bases | share of the never-cut | why the budget did not cut it |
|---|---|---|---|
| **coding-gene introns** | **1,187.3 Mb** | **58.5%** | inside an annotated gene span |
| **non-coding gene bodies** | **713.6 Mb** | **35.1%** | inside an annotated gene span: `ncRNA_gene` is a gene |
| UTR and the non-coding exons of coding genes | 92.6 Mb | 4.6% | inside an annotated gene span |
| CDS of a non-canonical transcript | 29.7 Mb | 1.5% | inside an annotated gene span |
| pseudogene bodies | 6.3 Mb | 0.31% | inside an annotated gene span: `pseudogene` is a gene |
| inter-gene gaps under 1,000 bases | 1.33 Mb | **0.066%** | the length rule, and nothing else |
| a tail under 1,000 bases | 0, measured | 0% | the length rule, and nothing else |
| in none of the eight classes | 0, measured | 0% | the cut and the annotation would have to disagree |

**99.93% of the never-cut sequence is declined for being inside an annotated gene span.** That line
is the finding, and it decides what can be asked. "The budget did not cut it" and "it is a gene body"
are not two hypotheses about those bases; they are one predicate written twice. No sample size
separates them, because it is a property of the cutting rule and not of the sample. Only 1.33 Mb --
sequence declined for its length alone -- is ordinary intergenic sequence where the two come apart,
and that is 0.066% of the mass.

Every group prints the covariates it has to be compared on. Coverage is the share of the panel's 89
assemblies informing the group, which is effort and not composition; the unit is the kilobase the
panel bins by, and a kilobase belongs to the class holding at least half of it.

| group | kilobases | bases | GC | median distance to a coding TSS | coverage | recurring/kb |
|---|---|---|---|---|---|---|
| the background | 2,833,407 | 2,807.6 Mb | 0.399 | 87 kb | 1.000 | 6.63 |
| the never-cut background | 2,030,746 | 2,005.6 Mb | 0.402 | 76 kb | 1.000 | 6.39 |
| coding-gene introns | 1,210,315 | 1,192.4 Mb | 0.405 | 57 kb | 1.000 | **6.16** |
| non-coding gene bodies | 713,606 | 713.5 Mb | 0.392 | 143 kb | 1.000 | **6.90** |
| UTR and non-coding exons | 82,698 | 79.2 Mb | 0.421 | 10 kb | 1.000 | **5.43** |
| mixed: no class holds half | 15,865 | 12.7 Mb | 0.494 | 8 kb | 1.000 | 5.97 |
| pseudogene bodies | 6,264 | 6.3 Mb | 0.421 | 68 kb | 1.000 | 8.04 |
| CDS of a non-canonical transcript | 1,161 | 0.8 Mb | 0.581 | 5 kb | 1.000 | 6.20 |
| inter-gene gaps under 1,000 bases | 837 | 0.8 Mb | 0.487 | 4 kb | 1.000 | 6.85 |

**What each class accounts for.** The background is rebuilt with each class taken out of it through
the panel's own `build_background` and the tier ratios recomputed, exactly as the section above did
for exons, conserved elements and promoter proximity. Every row is assessed on 23 of 23 chromosomes.

| taken out of the background | fossil (1.062) | regulatory (1.077) | neutral (1.105) | constrained unknown (0.916) | structural (1.163) |
|---|---|---|---|---|---|
| **coding-gene introns** | 1.027, **56%** | 1.027, **66%** | 1.073, **30%** | 0.889, *widens by 32%* | 1.170, *widens by 4%* |
| UTR and non-coding exons | 1.057, **8%** | 1.060, **23%** | 1.100, 4% | 0.914, *widens by 3%* | 1.150, **8%** |
| non-coding gene bodies | 1.062, accounts for none | 1.091, *widens by 18%* | 1.106, accounts for none | 0.918, accounts for none | 1.181, *widens by 11%* |
| CDS of a non-canonical transcript | 1.062, none | 1.076, none | 1.104, none | 0.916, none | 1.162, none |
| pseudogene bodies | 1.061, none | 1.079, none | 1.104, none | 0.916, none | 1.160, none |
| inter-gene gaps under 1,000 bases | 1.062, none | 1.077, none | 1.105, none | 0.916, none | 1.163, none |

**Coding-gene introns are the largest single covariate anyone has found for this offset**: 56%
(fossil), 66% (regulatory) and 30% (neutral), against 41%, 56% and 25% for the three of the section
above *together*. They are also 58.5% of the leftover by bases and the quietest large class in the
background after the UTRs. Non-coding gene bodies run the other way for the regulatory tier -- they
are noisier than the leftover around them (6.90 against 6.39) and removing them *widens* the offset by
18% -- which is the same shape as segmental duplication in the section above and the same warning: a
class large enough to move an average can move it either way, and "it is in the background" is not an
argument about direction.

**Two arms of that decomposition are degenerate, and the result says so rather than this paragraph.**
Take *every* gene body out of the background and what is left is 99.8% block sequence (median over 23
chromosomes, never below 99.2%); take every never-cut base out and it is 100%. That is a tier average
standing in for a control, which this project does not accept, so both arms keep their numbers, carry
`"degenerate": true` beside the share of the remaining background that is blocks, and are excluded
from every explained share on all 23 chromosomes. They are reported because the fact that *you cannot
take the gene bodies out of this background and still have a background* is the same finding as the
99.93% above, reached from the other end.

**The comparison that carries the weight.** 100-base windows, the hit being "does this window carry a
recurring event", through `compare.standardised`, stratified three ways, with canonical-CDS windows
out of both arms. Neither arm of any row is a tier average: the first row sets block sequence against
the sequence the cut declined, and every other sets one class against the *rest* of the never-cut, so
no tier appears in either arm anywhere in this table. Medians over 23 chromosomes.

| target | targets | controls | GC + timing | + distance to a coding TSS | + coverage | above 0 on | p | dropped |
|---|---|---|---|---|---|---|---|---|
| **every block, against the never-cut** | 699,309 | 1,844,353 | +0.0272 | +0.0269 | **+0.0255** | **22 of 23** | 3e-6 | 3,012 |
| CDS of a non-canonical transcript | 3,561 | 1,840,792 | -0.0434 | -0.0445 | **-0.0421** | 4 of 23 | 1.0 | 3 |
| UTR and non-coding exons | 92,164 | 1,752,189 | -0.0406 | -0.0373 | **-0.0362** | 0 of 23 | 1.0 | 187 |
| coding-gene introns | 1,075,024 | 769,329 | -0.0155 | -0.0182 | **-0.0167** | 1 of 23 | 1.0 | 1,293 |
| non-coding gene bodies | 649,872 | 1,194,481 | +0.0229 | +0.0243 | **+0.0230** | 21 of 23 | 3.3e-5 | 8,676 |
| pseudogene bodies | 9,298 | 1,835,055 | +0.0396 | +0.0386 | **+0.0365** | 17 of 23 | 0.017 | 48 |
| inter-gene gaps under 1,000 bases | 14,434 | 1,829,919 | +0.0111 | +0.0188 | **+0.0140** | 18 of 23 | 0.0053 | 79 |

`input_presence` was run on every one of those rows *before* the coverage stratum, and on all 23
chromosomes both inputs come back **free**: no window lacked a coverage value and none lacked a local
annotation, so no coverage artefact can live in the assignment to an arm. Free presence is not a
constant level, which is what a stratum bites on, and the level does vary -- 82,673 of the 699,309
block windows and 123,381 of the 1,844,353 never-cut windows are below full panel coverage -- so the
stratum had something to bite. The difference survives it on all three rungs, which is what licenses
the stronger of the two sentences the result carries as data rather than as prose: **"nothing I
measured explains it, including coverage"**.

**The verdict, which is about mass.** The result carries one of four sentences and it carries this
one: *the cutting rule and the biological class are the same predicate on all but a fraction of the
never-cut bases, so these data cannot assign the offset to one rather than the other; what they do say
is that the panel's background is mostly gene-body sequence and that gene-body sequence carries fewer
recurring events than block sequence at matched GC, timing, distance to a coding TSS and coverage.*
The brief that opened this lane asked whether the offset tracks "sequence the budget declined to cut
into a block" rather than any biological class. It tracks both, because on 99.93% of the bases those
are one predicate. That is not a dodge; it is the strongest true thing here, and the two shorter
readings -- "the tier comparison is an artefact of block-cutting", "the offset is introns" -- are both
quotable from it and both wrong as stated.

**The mechanism, a different question on a different sample, and it is on the line.** On the 0.066%
where the two readings do come apart -- ordinary intergenic sequence declined for its length alone --
the sequence is *not* quiet: it reads +0.0140 above the rest of the never-cut where the blocks read
+0.0255, which puts it 55% of the way from the never-cut sequence to the blocks. The reading is that
the quietness travels with the gene body rather than with the decision not to cut. It is decided by a
margin of **0.00125 points** against a threshold chosen by hand, so the result carries
`mechanism_on_the_line: true` and the sentence that goes with it: this is the direction the numbers
point, not a result that survives its own threshold. Those windows are also gene-proximal by
construction (median 4.8 kb to a coding TSS against 77 kb) and GC-rich (0.51 against 0.40), and while
the comparison holds both, a 14,434-window arm carrying a 0.066% class is not where this gets settled.

**What this changes.** No tier conclusion, and the reading of five claims:

1. **Every "x% above the matched background" for a non-coding tier.** The +6 to +11 points were
   assigned to the control on 2026-09-16; this says what the control is. That comparison is, on 71%
   of its bases, block sequence against gene bodies, and coding-gene introns alone carry 56% to 66%
   of the gap for fossil and regulatory. The honest form of a tier statement is still against the
   neutral tier, as it has been since that section.
2. **Every per-block `core` call in `human_panel`.** `block_class` calls a block core at half its
   expected rate or less and the expectation comes from this background through `bg_rt.expected`, so
   a core block is depleted relative to a baseline that is 71% gene body, not relative to intergenic
   sequence. That is the "core of callable" column in the chr21 and chr22 tables.
3. **The 200-bp unit expectations** (`expected_gc_rt_matched`). The reference units are every piece
   of the grid that is not canonical CDS, which is the same 71%.
4. **The matched control windows** drawn per block. `forbidden` is canonical CDS plus the structural
   tier, so a control window may land inside a gene body, and most of the chromosome it can land in
   is one.
5. **The section above's "44% to 75% unexplained".** Most of it now has a name. Exons, conserved
   elements and promoter proximity account for 41%/56%/25%; coding-gene introns alone account for
   56%/66%/30%. The two sets overlap and are not additive, and their union was not measured here,
   which is the first thing left undone below.

The constrained-unknown tier is untouched, as through every correction this week: still the one tier
below the background (0.916, above 1 on 6 of 23), and taking the introns out pushes it *further*
below (0.889), which is the opposite of an explanation.

**Coverage of the measurement.** Eligibility is above; this is what was done with it. 2,833,407
kilobases measured over the 23 chromosomes, 2,030,746 of them outside every block, and 2,596,722
windows after a fixed subsampling step between 4 and 20 per chromosome. 13,035 windows were taken
outside that step because they carry one of the two length-rule classes, which are the discriminator
and would otherwise have had almost no power; they are 0.5% of the windows, and being the noisy ones
they make the block-against-never-cut difference *smaller* rather than larger. 52,677 windows holding
canonical CDS were excluded from both arms. 738 windows were measured and placed in neither arm,
being split between a block and a class with neither holding half, and are counted rather than pushed
into one. No window lacked a coverage value on any chromosome. Every class that exists exists on all
23, so no share in the tables above is a zero that means "never looked". Nothing was fetched: the
stores, the committed budgets and GENCODE were all on disk, and the top-level distillation can be
redone from the committed result with `--from-result` without reading a store again.

**What is left undone.**

- **The union of the covariates, which is the number a reader will want and this lane does not have.**
  Coding-gene introns account for 56%/66%/30% and the section above's three for 41%/56%/25%, and
  nobody has taken all four out of the background at once. They overlap heavily -- conserved elements
  and exons sit inside gene bodies -- so the union is somewhere between the larger of the two and
  their sum, and quoting either end as "the offset is explained" would be wrong. It is one more
  exclusion arm and it is the first thing to do here.
- **Whether an intron is quiet for being an intron.** Nothing here separates transcription-coupled
  repair, alignment quality inside genes, and the correlation between gene bodies and early
  replication that survives a three-tertile stratum. This lane measured *where* the quiet sequence
  is, not why it is quiet.
- **The mechanism arm.** Decided by 0.00125 points. A purpose-built comparison -- intergenic sequence
  matched to gene-body sequence on distance to a TSS rather than stratified on it, which is what the
  GC 0.51-against-0.40 imbalance is really saying -- would settle it or show that it cannot be.
- **The `mixed` kilobases.** 15,865 of them, 12.7 Mb, the quietest group in the covariate table at
  5.97 per kb, and excluded from every class comparison by the half-a-kilobase rule. Too small to
  move a share and too large for anyone to call them nothing.
- **chrY.** Measured and excluded from the medians by the panel's own convention; its numbers sit in
  `per_chromosome` and nobody has read them.

Written to `data/results/panel_leftover_classes.json`.

## All four covariates out at once: 99% of the fossil offset, 141% of the regulatory, which is an overshoot rather than an explanation (2026-09-17, later still)

The two sections above each took a piece out of the offset by which every non-coding tier reads
above the panel's matched background, and neither took the other's. Exons of any gene, conserved
elements and promoter proximity account for 41% / 56% / 25% (fossil / regulatory / neutral)
together; coding-gene introns alone account for 56% / 66% / 30%. Nobody had taken all four out of
the background at once, and the last section named that as the first thing left undone.
`genomeos/attribution/panel_union.py` and `scripts/panel_union_arms.py` do it on 24 chromosomes from
local data only, by the same method, through the same `build_background`, with the same degeneracy
flag and the same two permitted sentences.

**What was eligible.** Over the 23 chromosomes the panel's medians rest on: 3,031 Mb of sequence,
2,846,945 kilobases inside the panel's alignment blocks, 28,468,458 100-base cells inside those, 5
tiers and 4 declared covariates, of which 4 have bases somewhere. Nothing below would mean anything
if the rebuilt background were not the panel's own, so before anything is read from it the rebuilt
ratios are checked against the committed panel's: 144 tier ratios on 24 chromosomes, largest
absolute difference 0.0005.

**The cheap half first: how much the four share.** The union has to sit between the largest single
arm and the sum of the singles, and where in that interval is decided by how much sequence the
covariates hold in common. Over the 23 chromosomes the four cover 1,466.9 Mb between them and
1,635.3 Mb when each is counted on its own: 168.4 Mb counted more than once, a sum 1.1148 times the
union.

| pair | bases in common | share of the first | share of the second |
|---|---|---|---|
| coding-gene introns and exons of any gene | 15.9 Mb | 1.3% | 7.6% |
| coding-gene introns and conserved elements | 55.4 Mb | 4.7% | 35.5% |
| coding-gene introns and within 2 kb of a coding TSS | 23.5 Mb | 2.0% | 30.3% |
| exons of any gene and conserved elements | 46.2 Mb | 21.9% | 29.6% |
| exons of any gene and within 2 kb of a coding TSS | 29.8 Mb | 14.1% | 38.5% |
| conserved elements and within 2 kb of a coding TSS | 8.4 Mb | 5.4% | 10.9% |

**Every group, with the covariates it has to be compared on.** Coverage is assemblies informing the
window, as a fraction of the panel's assemblies, which is effort and not composition. The unit of a
background row is the kilobase the panel bins by; the unit of a tier row is its block. The tier rows
are here to be described, not to be used: no tier average is a control anywhere in this section.

| group | units | bases | median length | GC | median distance to a coding TSS | coverage | recurring/kb |
|---|---|---|---|---|---|---|---|
| the background | 2,833,407 | 2,807.6 Mb | 1,000 | 0.399 | 87 kb | 1.0 | 6.626 |
| **the background free of all four** | 1,307,462 | 1,307.4 Mb | 1,000 | 0.391 | 145 kb | 1.0 | 7.117 |
| kilobases holding any of the four | 1,525,945 | 1,500.2 Mb | 1,000 | 0.407 | 55 kb | 1.0 | 6.23 |
| kilobases holding coding-gene introns | 1,244,198 | 1,222.6 Mb | 1,000 | 0.406 | 55 kb | 1.0 | 6.154 |
| kilobases holding exons of any gene | 297,246 | 280.1 Mb | 1,000 | 0.436 | 24 kb | 1.0 | 6.372 |
| kilobases holding conserved elements | 226,332 | 215.5 Mb | 1,000 | 0.39 | 85 kb | 1.0 | 5.34 |
| kilobases holding within 2 kb of a coding TSS | 86,116 | 82.6 Mb | 1,000 | 0.456 | 1 kb | 1.0 | 6.398 |
| tier: fossil | 7,085 | 316.0 Mb | 19,422 | 0.3897 | 56 kb | 0.9963 | 7.072 |
| tier: regulatory | 15,516 | 344.5 Mb | 14,279 | 0.44 | 24 kb | 0.9984 | 7.195 |
| tier: neutral | 2,543 | 71.7 Mb | 7,874 | 0.4048 | 58 kb | 0.9872 | 7.4 |
| tier: constrained unknown | 1,068 | 31.7 Mb | 7,676 | 0.3902 | 80 kb | 0.9987 | 6.261 |
| tier: structural | 219 | 198.1 Mb | 72,680 | 0.4181 | 213 kb | 0.7395 | 10.41 |

**The arms.** The background is rebuilt with each covariate taken out of it and then with all four
taken out, through the panel's own `build_background`, and the tier ratios recomputed. The ratio the
panel itself builds is in each column heading; the cell is the ratio after the exclusion, and what
share of its distance from 1 that closed. Every row is assessed on 23 of 23 chromosomes.

| taken out of the background | fossil (1.062) | regulatory (1.0774) | neutral (1.1045) | constrained unknown (0.9163) | structural (1.1627) |
|---|---|---|---|---|---|
| coding-gene introns | 1.027, **56%** | 1.0266, **66%** | 1.0731, accounts for 30% of the offset | 0.8891, widens the offset by 32% | 1.1695, widens the offset by 4% |
| exons of any gene | 1.0538, accounts for 13% of the offset | 1.0519, accounts for 33% of the offset | 1.0931, accounts for 11% of the offset | 0.912, widens the offset by 5% | 1.1498, accounts for 8% of the offset |
| conserved elements | 1.0425, accounts for 31% of the offset | 1.0582, accounts for 25% of the offset | 1.0862, accounts for 18% of the offset | 0.9008, widens the offset by 19% | 1.1457, accounts for 10% of the offset |
| within 2 kb of a coding TSS | 1.0616, accounts for none | 1.0726, accounts for 6% of the offset | 1.1015, accounts for 3% of the offset | 0.9158, accounts for none | 1.187, widens the offset by 15% |
| **all four at once** | 1.0008, **99%** | 0.9683, **141%** | 1.0465, **56%** | 0.866, widens the offset by 60% | 1.1405, accounts for 14% of the offset |

**Where the union sits, which is the number this lane exists to produce.**

| tier | largest single arm | sum of the singles | **the union** | where in that interval | what that reads as |
|---|---|---|---|---|---|
| fossil | 56.5% | 101.8% | **98.7%** | 0.9323 | the union is near the sum of the singles |
| regulatory | 65.6% | 129.6% | **141.0%** | 1.1778 | the union explains MORE than the sum of its singles |
| neutral | 30.0% | 61.3% | **55.5%** | 0.8134 | the union sits between the largest single arm and the sum, as overlapping covariates do |
| constrained unknown | -0.6% | -56.8% | **-60.1%** | - | the union explains LESS than its largest single arm |
| structural | 10.4% | -0.7% | **13.6%** | - | the union explains MORE than the sum of its singles |

**Neither end of the interval was the answer, and on 3 of the 5 tiers the union fell outside the
interval altogether (regulatory, constrained unknown, structural).** The fossil union is 98.7%
against a largest single of 56.5% and a sum of 101.8%: quoting the largest single would have
understated it by 42.3%, and the sum happens to land within 3.1% of it. The neutral union is 55.5%
against 30.0% and 61.3%, sitting 0.8134 of the way from the first to the second, which is the shape
an ordinary set of overlapping covariates makes.

**The regulatory tier is the exception, and it is the one figure in this section that must not be
quoted as good news.** Its union reads 141.0%, and a share above 100% is not an offset more than
explained. It is the ratio carried past 1 and out the other side: the tier reads 1.0774 against the
background as the panel builds it and 0.9683 against the background with all four gone, so the
exclusion did not land on 1, it went through it. The result carries that as a flag and a sentence
rather than leaving it to a reader: *the union does not close this tier's offset, it crosses it:
with all four covariates gone the tier reads BELOW the background the exclusion leaves, so the share
above 100% is over-correction and not explanation, and the covariates removed were quieter than the
tier rather than merely as quiet*. It is true of regulatory and of no other tier. Over-correction on
1 of 5 tiers is also why the arithmetic below stops at what the arms did and does not go on to call
the offset accounted for.

**The degeneracy check, which is what both earlier lanes hit and what could have made this
unanswerable.** Take every gene body out of this background and the remainder is 99.8% block
sequence, which is a tier average standing in for a control; two arms of the last section are
reported without being counted for exactly that reason. The question here was whether the union arm
goes the same way. It does not. Each arm carries the share of the background it leaves that is block
sequence, against the same line the last section drew (0.95), and the union arm leaves a background
that is 0.5261 block sequence at the median and 0.7093 on its worst chromosome. It is degenerate on
0 of 23, so the question can be answered this way and the table above is evidence rather than a
reported artefact.

| arm | bases taken out of the background | block share of what it leaves, median | worst chromosome | degenerate on |
|---|---|---|---|---|
| coding-gene introns | 1,224.9 Mb | 0.457 | 0.6606 | 0 of 23 |
| exons of any gene | 210.8 Mb | 0.2836 | 0.4659 | 0 of 23 |
| conserved elements | 166.1 Mb | 0.2711 | 0.4553 | 0 of 23 |
| within 2 kb of a coding TSS | 106.2 Mb | 0.2677 | 0.4476 | 0 of 23 |
| **all four at once** | 1,466.9 Mb | 0.5261 | 0.7093 | 0 of 23 |

**The comparison that carries the weight, run both ways round one control.** 100-base windows, the
hit being "does this window carry a recurring event", through `compare.standardised`, stratified
three ways, with canonical-CDS windows out of both arms. The control is the same in every row and it
is sequence, never a tier average: the windows that carry not one base of any of the four and lie in
no block. The first row asks whether the covariates are quieter than what removing them leaves,
which is what an exclusion arm needs to be true; the rest ask what each tier reads against a
background with the four already gone. Medians over 23 chromosomes.

| target | targets | controls | GC + timing | + distance to a coding TSS | + coverage | above 0 on | p | dropped |
|---|---|---|---|---|---|---|---|---|
| **the four covariates, against the rest** | 1,331,006 | 500,729 | -0.0411 | -0.0397 | **-0.0386** | 0 of 23 | 1 | 30,738 |
| tier: fossil | 265,472 | 500,729 | +0.0092 | +0.0091 | **+0.0097** | 16 of 23 | 0.0466 | 2,975 |
| tier: regulatory | 317,851 | 500,729 | -0.0045 | -0.0036 | **-0.0014** | 11 of 23 | 0.661 | 7,272 |
| tier: neutral | 66,202 | 500,729 | +0.0125 | +0.0131 | **+0.0132** | 16 of 23 | 0.0466 | 995 |
| tier: constrained unknown | 29,705 | 500,729 | -0.0227 | -0.0239 | **-0.0216** | 4 of 23 | 1 | 247 |
| tier: structural | 20,079 | 500,729 | +0.0049 | +0.0029 | **+0.0216** | 14 of 23 | 0.202 | 2,301 |

`input_presence` was run on every one of those rows *before* the coverage stratum, which is the
order that stops "conditioning changed nothing" from reading as a test passed when no test was
administered. On the first row: the panel's alignment over the window is **free** on 23 of 23; the
local annotation of the window is **free** on 23 of 23. Free presence is not a constant level, and
the level does vary -- 77,351 of the 1,331,006 covariate windows and 44,594 of the 500,729
covariate-free windows are below full panel coverage -- so the stratum had something to bite, and
the difference survives it on all three rungs.

**The sentence this result carries, as data rather than as prose:** *"nothing I measured explains
it, including coverage"*. It is one of two constants and there is no third: without a coverage
stratum the strongest available sentence is "the covariates I measured do not explain it", and
"nothing explains it" is not available at any sample size, because the next covariate nobody has
measured is always still there.

**What this changes.** One number, and one reading:

1. **The last section's "44% to 75% unexplained" and its successor.** With all four out at once, per
   tier, in the result's own words: fossil accounts for 99% of the offset; regulatory accounts for
   141% of the offset; neutral accounts for 56% of the offset. That is the number to quote, and it
   is neither the largest single arm nor the sum.
2. **The sum was never the number.** For the regulatory tier the sum of the singles is 129.6% --
   more than the whole offset, which is what adding overlapping covariates buys, and 168.4 Mb of the
   sequence behind it is counted twice. Adding the two sections' numbers would have been wrong on
   this tier by 11.4% too low, and the direction it is wrong in is not fixed: the sum overstates the
   union on fossil, neutral, constrained unknown and understates it on regulatory, structural.
3. **No tier conclusion moves.** The constrained-unknown tier is still the one tier below the
   background (0.9163), and with all four out it reads 0.866: the arm widens the offset by 60%,
   which is the same shape the last two sections found and the opposite of an explanation.

**Coverage of the measurement.** Eligibility is above; this is what was done with it. 2,833,407
kilobases measured over the 23 chromosomes, 1,307,462 of them free of all four covariates, and
2,583,687 windows after a fixed subsampling step between 4 and 20 per chromosome. The kilobases not
measured are the two named categories: 10,472 below half a kilobase after the coding cut, 3,066 no
gc in the store. Every covariate that exists exists on all 23, so no cell in the tables above is a
zero that means "never looked"; where one would have been, the cell says `not assessed` and carries
the reason. chrY is measured and reported and kept out of every median, as everywhere else in this
panel, and its numbers sit in `per_chromosome`. Nothing was fetched: the stores, the committed panel
results and GENCODE were all on disk, and the top-level distillation can be redone with
`--from-result` without reading a store again.

**What is left undone.**

- **Segmental duplication, the fifth covariate, which runs the other way.** The first section
  measured it and it *widens* the offset when removed -- fossil 29%, regulatory 34%, neutral 33%. It
  is not in this union, so the union above is a union of the four that close the offset and says
  nothing about what a fifth that opens it would do to the total. A five-way arm is the obvious next
  exclusion and it is one `COVARIATES` entry away.
- **The union arm is not degenerate, but it is not neutral either.** What it leaves is 0.5261 block
  sequence at the median against 0.7093 at the worst. That is far below the line, and it is also far
  above where a background that shared nothing with the tiers would sit. Nobody has asked what a
  background enriched for blocks does to a ratio short of degeneracy.
- **Why an intron is quiet.** Carried forward unchanged from the last section: nothing here
  separates transcription-coupled repair, alignment quality inside genes, and the correlation
  between gene bodies and early replication that survives a three-tertile stratum. This lane
  measured how much the four covariates overlap, not why any of them is quiet.
- **The remainder has no name, and on one tier it is not a remainder.** What all four together
  leave: fossil 1.3% still standing, regulatory over-corrected by 41.0%, neutral 44.5% still
  standing. This lane adds no candidate for what is left, and none for why the regulatory arm goes
  past the mark. The honest form of both is the sentence the result carries, not a fifth covariate
  chosen because it was to hand.
- **chrY.** Measured and excluded from the medians by the panel's own convention; its numbers sit in
  `per_chromosome` and nobody has read them.

Written to `data/results/panel_union_arms.json`.

## Compression: what each piece of knowledge is worth in bits, chr21, chr22 and chr18 (2026-09-14)

Albert's probe: try different ways of classifying the genome and find which buys the most
understanding per unit of effort. This one is information-theoretic and needs no labels to score
against. The best model of a sequence is the one that describes it in the fewest bits, so each
piece of knowledge the project holds is turned into a model that predicts the next base. The
chromosome is coded under that model, and the gain is read in bits against generic baselines. A
layer that does not shorten the code is not knowledge about the sequence, however it reads in a
table.

Code: `attribution/compress.py`, run by `scripts/compress.py chr21:chr22 chr22:chr21 chr18:chr22`.
Results: `compress_chr21.json`, `compress_chr22.json` and `compress_chr18.json`, plus the
sensitivity run `compress_chr21_trained_chr20_chr22.json`. No model API is
called. The only requests are UCSC's CpG islands, one per chromosome, cached under
`data/knowledge/compress`.

**How it is kept honest.**
- **Code lengths.** They are ideal code lengths, -log2 of the probability given to the base that
  came. An arithmetic coder reaches them within two bits per chromosome.
- **Held-out fitting.** Every static model is fitted on both strands of another chromosome (chr22
  for chr21 and chr18, chr21 for chr22). It is applied to the chromosome being coded, and its
  parameters are not charged.
- **Adaptive models.** These are fitted on the chromosome itself and learn as they code: counts
  are updated after each base, and a reverse-strand context is added once its bases are known. The
  decoder can rebuild them, so they pay for their parameters.
- **Annotations are paid for.** A layer's annotation is side information the decoder needs: where
  each repeat is and its subfamily and divergence, exon coordinates and phase, the islands, the
  cCREs, the sites, the blocks and their tiers. Its description length is subtracted from the
  layer's gain. Positions and lengths use an adaptive bit-length code, labels an adaptive
  categorical code. Label dictionaries belong to the fitted models and are reported, not charged
  (266 kb for RepeatMasker's names).
- **The mixture.** Models are combined by a causal windowed Bayesian mixture: a model's weight at
  a base is its likelihood over the previous W bases. W and the temperature come from a grid of 12
  on the generic stack, and the log2(12) bits of that choice are charged; W = 8 was chosen on chr21
  and chr22, W = 16 on chr18.
- **A layer's models are inactive outside the region it claims.** The decoder knows that region,
  because the annotation is transmitted. The first run left them active everywhere, and each layer
  then cost more outside its claim than it saved inside. The motif sites saved 13 kb in their 22 kb
  and lost 87 kb elsewhere on chr21, because a model that says "two bits" everywhere still takes
  weight in a mixture. That was the mixing scheme failing, not the knowledge.
- **CpG islands, decided in advance.** genomeos-i1's lexicon was swamped by CG words because a
  chain fitted over a whole context cannot absorb island clustering. Here GC and CpG structure is
  its own layer: a causal composition state over the previous kilobase, active everywhere, plus the
  island annotation. The cCRE, motif and tier layers are also read left out of the full stack, so
  what they add is measured on top of that structure.
- **Repeats dominate, as expected.** Every later layer is therefore also reported over the
  repeat-aware stack.
- **A control for the duplications.** The duplication layer reads partners on earlier chromosomes,
  which the chromosome-local generic models cannot see. The control primes the adaptive models with
  the same partner sequence but without the alignment.

**Baselines, bits per base** (A/C/G/T only; the assembly gaps cost 1.5 to 1.8 kb to describe):

| model | chr21 | chr22 | chr18 |
|---|---|---|---|
| order 0, adaptive | 1.976 | 1.997 | 1.970 |
| order 2, adaptive | 1.924 | 1.932 | 1.918 |
| best single order, adaptive (k) | 1.743 (15) | 1.685 (14) | 1.745 (15) |
| best short order, held out (k) | 1.783 (10) | 1.726 (10) | 1.822 (9) |
| bzip2 -9, 2-bit packed | 1.824 | 1.786 | 1.804 |
| xz -9e, 2-bit packed | 1.706 | 1.673 | 1.700 |
| xz -9, letters | 1.700 | 1.667 | 1.746 |
| naive: static orders 1 to 12 mixed, fitted on the other chromosome | 1.692 | 1.633 | 1.745 |
| generic: naive plus adaptive orders 12 to 24 (blind copy finding) | 1.590 | 1.545 | 1.650 |
| generic primed with the duplication partners (control) | 1.528 | 1.475 | 1.636 |
| every layer, models only | 1.470 | 1.418 | 1.607 |
| every layer, annotation paid | 1.532 | 1.499 | 1.665 |
| best net: primed generic plus duplications | **1.500** | **1.446** | **1.630** |

**Which order stops helping.**
- **Short-range statistics saturate at order 9 to 10.** Held out, order 10 is the best short
  context (1.783 on chr21). Adaptive single orders get no better from 9 to 12, as contexts outrun
  the data.
- **From order 13 the adaptive model becomes a copy finder.** It is best at 14 or 15 and flat to
  16.
- **In-sample fitting cheats.** The in-sample entropy falls to 0.10 bits per base at order 16: the
  chromosome memorised, not modelled.
- **Against the general compressors.** The generic mixture beats xz by 0.12, 0.13 and 0.05 bits per base on
  chr21, chr22 and chr18.

**The layers.** Figures are net of annotation, per claimed base and per megabase of the whole
chromosome, given as chr21 / chr22 / chr18.

| layer | claims | annotation, bits per claimed base | saved over generic, gross | **net per claimed base** | **net kb per Mb** | net kb per Mb, left out of the full stack | pays? |
|---|---|---|---|---|---|---|---|
| duplications | 12.6 / 13.3 / 2.0% | 0.016 / 0.027 / 0.030 | 0.69 / 0.71 / 0.94 | **+0.68 / +0.68 / +0.91** | **+85 / +91 / +18** | +73 / +76 / +17 | yes, everywhere |
| GC and CpG | all bases | 0.0003 / 0.0005 / 0.0002 | 0.0075 / 0.0071 / 0.0057 | +0.007 / +0.007 / +0.006 | +7.2 / +6.6 / +5.6 | +3.5 / +3.2 / +3.9 | yes, everywhere |
| tiers and blocks | 36 / 29 / 37% | 0.001 / 0.002 / 0.001 | 0.020 / 0.020 / 0.001 | +0.019 / +0.018 / -0.000 | +6.7 / +5.3 / -0.1 | -1.1 / -1.8 / -1.2 | only over generic, on two chromosomes |
| repeats | 52 / 53 / 51% | 0.097 / 0.112 / 0.096 | 0.072 / 0.074 / 0.043 | **-0.025 / -0.038 / -0.054** | -12.6 / -20.4 / -27.4 | -25.3 / -35.8 / -29.9 | no; +3.8 kb per Mb on chr21 when fitted on two chromosomes |
| cCREs | 9 / 15 / 8% | 0.081 / 0.080 / 0.081 | 0.006 / 0.011 / 0.004 | -0.075 / -0.068 / -0.078 | -6.8 / -10.1 / -6.5 | -7.1 / -11.2 / -6.7 | no |
| coding (codon phase) | 0.8 / 1.8 / 0.7% | 0.144 / 0.137 / 0.136 | 0.034 / 0.049 / 0.028 | -0.11 / -0.09 / -0.11 | -0.9 / -1.6 / -0.7 | -0.9 / -1.9 / -0.8 | no |
| motif sites | 0.06 / 0.09 / 0.04% | 1.67 / 1.68 / 1.61 | 0.53 / 0.53 / 0.47 | -1.14 / -1.15 / -1.14 | -0.6 / -1.1 / -0.4 | -0.7 / -1.2 / -0.4 | no |

- **Duplications pay, and most of that is access to the partner's sequence.** Generic models
  primed with the partners recover two thirds of the gain. The alignment itself still pays: +0.22,
  +0.22 and +0.33 bits per claimed base over the primed stack (+28, +29 and +7 kb per Mb). The
  copy model agrees with its partner at 89.5%, 89.0% and 84.0% of the bases it predicts. A
  genome-wide generic coder holds more than the partners, so over it the alignment's value is at
  most that.
- **GC and CpG pay a little, everywhere.** On island bases the layer takes 0.06 to 0.09 bits
  (1.577 to 1.518 on chr21, 1.860 to 1.773 on chr22, 1.890 to 1.817 on chr18). It is the only
  layer that needs no annotation to do most of its work.
- **Repeats do not pay over blind copy finding when fitted on one chromosome. Negative, and set by
  the training.**
  - 65,000 to 128,000 labelled copies cost 31 bits each and save 0.04 to 0.07 bits a base over
    models that find copies with no annotation. Over the naive stack they are about level: +7.3
    and +4.0 kb per Mb on chr21 and chr22, -0.9 on chr18.
  - **Fitted on two chromosomes** (chr21 fitted on chr20 and chr22,
    `compress_chr21_trained_chr20_chr22`), the gross gain over generic rises from 0.072 to 0.104
    bits per claimed base. The net turns from -0.025 to +0.007 (+3.8 kb per Mb) over generic, and
    from +7.3 to +24.5 kb per Mb over naive.
  - Left out of the full stack it still loses 12.5 kb per Mb, where the duplications and tiers
    already hold much of what it knows.
  - Every other verdict is unchanged by the doubled training: duplications +86.4, GC and CpG
    +7.9, cCREs -6.7, coding -0.8 and motifs -0.7 kb per Mb.
  - A repeat library fitted on the whole genome would likely pay. RepeatMasker's labels are worth
    about what the copies the models have already seen are worth.
- **Coding, cCREs and motif sites do not pay. Negative.**
  - Codon phase is worth 0.03 to 0.05 bits per coding base and costs 0.14 to point to.
  - A cCRE class is worth about 0.01 bits per base of its element or less.
  - A recorded JASPAR site does predict its bases (0.5 bits each, since sites were selected on
    their score), but pointing to it costs 1.6 bits per base. Motif knowledge in this form is
    worth less than its address.
- **Tiers and blocks are repeat content by another name. Negative.** They save 0.02 bits per
  claimed base over generic on chr21 and chr22 and nothing on chr18. Over the repeat-aware
  stack they save +0.005, +0.002 and -0.002, and left out of the full stack they lose bits.

**The tiers, and the prediction.** Bits per base under the naive, generic, repeat-aware and full
stacks. "Unique" is outside RepeatMasker and segmental duplications. chr18 first: its 2.87 Mb of
constrained_unknown in 43 blocks is the most the genome offers for its size, and only 1.5% of it is
duplicated.

| chr18 | Mb | interspersed | duplicated | naive | generic | +repeats | full | interspersed bases, full | unique, full | unique, full, GC-standardised |
|---|---|---|---|---|---|---|---|---|---|---|
| coding exons (calibration) | 0.52 | 0.3% | 5.1% | 1.959 | 1.945 | 1.945 | 1.861 | | 1.923 | 1.938 |
| structural | 5.52 | 1.0% (tandem 98.6%) | 1.6% | 1.220 | 0.179 | 0.186 | 0.192 | | | |
| fossil | 6.41 | 57.9% | 4.4% | 1.772 | 1.735 | 1.705 | 1.674 | 1.548 | 1.912 | 1.909 |
| regulatory | 9.65 | 51.0% | 1.5% | 1.755 | 1.725 | 1.698 | 1.685 | 1.501 | 1.915 | 1.911 |
| constrained_unknown | 2.87 | 42.7% | 1.5% | 1.828 | 1.811 | 1.795 | 1.768 | 1.626 | **1.921** | **1.918** |
| neutral | 4.88 | 43.3% | 1.6% | 1.821 | 1.795 | 1.778 | 1.760 | 1.609 | **1.911** | **1.911** |
| chromosome | 80.09 | 42.7% | 2.4% | 1.745 | 1.650 | 1.628 | 1.607 | 1.518 | 1.916 | 1.913 |

| full stack | fossil | regulatory | constrained_unknown | neutral | coding exons | unique, constrained_unknown | unique, neutral |
|---|---|---|---|---|---|---|---|
| chr21 (duplicated share) | 1.473 (23%) | 1.626 (8%) | 1.144 (61%) | 1.475 (22%) | 1.788 | 1.934 | 1.927 |
| chr22 (duplicated share) | 1.238 (37%) | 1.569 (7%) | 1.258 (42%) | 0.935 (64%) | 1.837 | 1.918 | 1.888 |

- **Constrained_unknown compresses like neutral sequence. Negative, stated loudly.**
  - **All bases, chr18:** 1.768 against 1.760 under the full stack. **Unique bases:** 1.921
    against 1.911.
  - **Block bootstrap, unique bases, 38 blocks against 82:** the difference is +0.0006 (-0.0010 to
    0.0022) under the naive stack, +0.0075 (0.0048 to 0.0112) under the generic stack and +0.0105
    (0.0073 to 0.0146) under the full stack.
  - The sign is the predicted one: constrained sequence carries fewer unannotated near-copies than
    neutral sequence. The size is half a percent of a two-bit alphabet. At base level the 32 Mb
    thought most interesting holds no more information than unconstrained unique sequence, and no
    less.
  - On chr21 and chr22 the tier compresses easier than neutral, because 61% and 42% of it are
    copies. The organiser's copies-first rule, measured in bits. With 5 and 11 blocks the unique
    comparison has no power there: +0.007 (-0.003 to 0.015) and +0.030 (-0.017 to 0.069).
- **Fossils do not "compress hard". Negative for the prediction.**
  - Their repeat bases cost 1.55 bits per base under the full stack on chr18 (1.68 under the
    naive one), 1.34 on chr21 and 1.10 on chr22. Old copies with a mean divergence around 20% leave a quarter of the letters to
    chance.
  - The fossil tier as a whole sits 0.09 bits below neutral on chr18 and level with it on chr21.
  - Its repeat bases are 0.06 bits cheaper than the neutral tier's (1.61) and dearer than the
    regulatory tier's (1.50).
- **No model moves unique sequence.** Every stack codes it at 1.89 to 1.94 bits per base on every
  tier and chromosome, coding exons included: 1.94 on chr18 against 1.91 for neutral once
  GC-standardised. The compressible part of a chromosome is its copies, tandem arrays (structural
  0.18 bits per base on chr18) and old repeats. Function without repetition is invisible to this
  instrument, which is the limit of the method, not only of the tiers.

**Efficiency.**

| | chr21 (40.1 Mb coded) | chr22 (39.2 Mb) | chr18 (80.1 Mb) |
|---|---|---|---|
| wall clock | 1,132 s | 1,125 s | 2,356 s |
| seconds per Mb | 28.2 | 28.7 | 29.4 |
| peak memory | 7.5 GB | 8.4 GB (same process as chr21) | 10.0 GB (model codes spilled to disk, 5.8 GB) |
| read from local disk | 2.94 GB | 2.99 GB | 2.77 GB |
| requests, bytes streamed | 2 requests, 224 KB, for chr21 and chr22 on the first run | cached | 1 request, 97 KB |
| model calls | 0 | 0 | 0 |

- **Where the time goes (chr21).**
  - 45 mixtures, 524 s: 11.7 s each for 12 to 32 models.
  - The mixing grid, 95 s.
  - The order 0 to 16 table, 173 s.
  - The primed control, 115 s; the adaptive models, 79 s.
  - The layers: repeats 49 s, GC and CpG 20 s, tiers 18 s, duplications 7 s, coding 5 s, cCREs 2
    s, motifs 1 s.
  - Almost all of the disk read is the earlier chromosomes' FASTA for the duplication partners.
- **Result per unit of cost.**
  - Duplications buy +85 kb per Mb (+28 over the primed control) for 7 s and a read of the earlier
    chromosomes.
  - GC and CpG buy +7 kb per Mb for 20 s and one 100 KB request.
  - Repeats cost 49 s and lose 13 to 27 kb per Mb.
  - The rest cost under 20 s each and lose bits.
- **A genome-wide run.**
  - **As run here** (47 mixtures, the table and the control): about 29 s per Mb, so 2,900 Mb of
    bases is about 23 hours on one core.
  - **Memory, as run here.** It grows with the chromosome (7.5 GB at 40 Mb, 10 GB at 80 Mb with the
    models on disk), so chr1 and chr2 extrapolate past this 19 GB machine. That is what the reduced
    pass below is for.
  - **Disk and network.** Model spill is 2 bytes a base a model (17 GB for chr1, deleted after).
    Reading the earlier chromosomes for partners is about 3 GB per chromosome, 70 GB of local reads
    unless cached. Network is 24 requests, about 3 MB. Summaries are about 150 KB a chromosome.

**The reduced pass: the same verdicts in bounded memory** (`scripts/compress.py --pass`,
`compress_pass_<chrom>`; the chain is `scripts/compress_genome_wide.py`).

- **What it produces.** Every stack of the full run, every layer's gain net of its annotation over
  the naive, generic and repeat-aware stacks, cumulatively and left out of the full stack, the tier
  tables and the per-block sums the bootstrap needs. **What it drops:** the order 0 to 16 table, xz
  and bzip2, the mixing grid (the window and temperature chosen on chr21 and chr22 are fixed,
  W = 8, beta 1) and the partner-primed control.
- **One mixture, not forty-seven.** A model's weight in the windowed mixture depends only on its
  own code lengths, so a stack's mixture is the sum of its layers' weighted probabilities over the
  sum of their weights. All 33 stacks come out of one pass over the models, to within 1e-4 of
  mixing each separately (`test_one_pass_mixes_every_stack_as_separate_mixtures_would`).
- **Chunked context hashing.** Contexts are built by doubling rather than one base at a time, and a
  count is routed into hash groups of at most 16 M rows, each sorted on its own, so memory is
  bounded by a group whatever the chromosome; rows beyond 3 GB wait in temporary files. This is
  what chr1 and chr2 needed.
- **Segments.** The chromosome is mixed 32 Mb at a time. A segment reads the bases before it, so
  its models and its composition state are the whole chromosome's; only the mixing window (8 bases)
  restarts at a boundary, and each annotation item is charged to the segment holding its start, so
  its coordinate code restarts too. Both are conservative and both are small: on chr21, two
  segments, the annotation bill rises by 0.15% for repeats, 0.8% for coding and 1.2% for GC and
  CpG.
- **It reproduces the committed numbers.** `compress_pass_check_chr21` against the committed
  `compress_chr21`: all 33 shared stacks agree to four decimals (1.6918 naive, 1.5903 generic,
  1.4701 full), every layer's claimed bases and item counts are identical, and the tier bootstrap
  agrees (+0.0071, -0.0023 to 0.0144, against +0.0071, -0.0031 to 0.0149). The net per claimed base
  moves only where the segmented annotation bill moves it: repeats -0.0245 to -0.0247, coding
  -0.1106 to -0.1118, motifs -1.1401 to -1.1446. No verdict changes.
- **Cost, chr21 (40.1 Mb).** 641 s, 16.0 s per Mb, 4.6 GB peak memory, 0.32 GB of temporaries,
  2.94 GB read, 0 requests (cached), 0 model calls. The 23-hour full run becomes about 13 hours on
  an idle core. Stages: the four adaptive orders 111 s, the duplication layer 44 s, the per-segment
  models 401 s, the mixing 84 s.
- **The disk floor.** Temporaries are 34 bytes a base (the adaptive and copy code maps at 10, an
  adaptive model's count rows at 24), 8.5 GB at chr1. Free space is checked before every step and
  every 30 s, and the pass stops and deletes its temporaries rather than let free space fall below
  15 GB, which the other lanes need.

## Compression genome-wide: the tier gap holds over 901 blocks, and everything else gets worse (2026-09-14)

`scripts/compress_genome_wide.py` runs the reduced pass one chromosome to a process, then sums the
genome: `compress_pass_<chrom>` per chromosome, `compress_genome_wide` for the rollup. **Twenty-two
of twenty-four chromosomes, 2,467 Mb of the 2,875.** chr1 and chr2 are missing, as a negative: see
the disk paragraph below. Every stack, every layer's gain and the tier tables are summed over
bases; the bootstrap resamples blocks within their own chromosome.

**The question this was run to answer.** chr18's 43 constrained-unknown blocks carried the claim
that the tier holds no more information per base than neutral sequence. The genome holds 1,098 such
blocks over 32 Mb, and 901 of them (26.7 Mb, 13.9 Mb of it unique) are in these 22 chromosomes.

**The tier gap holds.** Constrained_unknown minus neutral, unique bases, bits per base, resampling
blocks within chromosome:

| stack | chr18 alone (43 blocks) | genome, 22 chromosomes (736 blocks with unique bases) |
|---|---|---|
| naive | +0.0009 (-0.0006 to 0.0024) | **+0.0030 (0.0017 to 0.0047)** |
| generic | +0.0084 (0.0055 to 0.0123) | **+0.0096 (0.0074 to 0.0122)** |
| repeat-aware | — | **+0.0096 (0.0074 to 0.0120)** |
| full | +0.0099 (0.0064 to 0.0136) | **+0.0117 (0.0095 to 0.0145)** |

- **It neither grows nor vanishes.** With seventeen times the blocks the estimate moves from
  +0.0099 to +0.0117 and the interval narrows from ±0.0037 to ±0.0025. chr18 was not a fluke and
  was not the whole story either.
- **It is still half a percent of the alphabet. Negative, and this is the finding.** 1.9299 against
  1.9182 bits per base on unique sequence; GC-standardised, 1.9265 against 1.9182. The 32 Mb the
  project calls its most interesting holds about one part in 170 more information per base than
  unconstrained unique sequence. Deep constraint across mammals is visible to this instrument only
  as a rounding error.
- **Per chromosome it is positive nearly everywhere and clear of zero in two thirds.** Twenty-one
  of the 22 chromosomes give a positive point estimate; chr15 alone is negative (-0.0002, -0.0069
  to 0.0053). Fourteen exclude zero. The two largest estimates, chrY (+0.178) and chr16 (+0.056),
  have the widest intervals and the fewest blocks, which is why the pooled figure is the one to
  quote.
- **Fossils are level with neutral. Negative for the prediction, now genome-wide.** +0.0013
  (-0.0013 to 0.0044) over 5,350 blocks. On chr18 it read as fossils sitting below neutral; over
  the genome there is no difference at all.

**Which layers pay, over 2,467 Mb.** Net of annotation, against the generic stack, and left out of
the full stack.

| layer | claims | annotation, bits per claimed base | net per claimed base | net kb per Mb | net kb per Mb, left out | pays? |
|---|---|---|---|---|---|---|
| duplications | 4.0% | 0.020 | **+0.541** | **+21.6** | +20.1 | yes |
| GC and CpG | all bases | 0.0003 | +0.0044 | +4.4 | +3.1 | yes |
| tiers and blocks | 29.9% | 0.0012 | **-0.0003** | **-0.1** | -1.1 | **no, and this is new** |
| repeats | 52.4% | 0.105 | -0.0625 | -32.7 | -34.8 | no |
| cCREs | 9.7% | 0.081 | -0.0753 | -7.3 | -7.6 | no |
| coding (codon phase) | 1.1% | 0.135 | -0.0986 | -1.1 | -1.2 | no |
| motif sites | 0.1% | 1.608 | -1.073 | -0.6 | -0.6 | no |

- **The tiers stop paying once the genome is looked at. Negative.** On chr21 and chr22 they were
  +0.019 and +0.018 bits per claimed base over generic, and the honest reading then was "repeat
  content by another name". Genome-wide they are -0.0003 over generic, -0.0024 over the
  repeat-aware stack and -0.0038 left out of the full stack. The chromosomes that made them look
  positive are the acrocentrics, whose tiers are mostly copies.
- **Duplications still pay and still pay most, but a third of what chr21 suggested.** +21.6 kb per
  Mb against +85 on chr21, because chr21's 12.6% duplicated share is an acrocentric short arm and
  the genome's is 5.7%. Per claimed base the layer is as strong as ever, +0.541.
- **GC and CpG still pay everywhere,** +4.4 kb per Mb for an annotation of three ten-thousandths of
  a bit per base.
- **Repeats lose more, not less, with scale.** -32.7 kb per Mb against -13 to -27 on the three
  chromosomes. 31 bits per labelled copy is not recovered when blind copy finding already has the
  sequence.
- **Nothing moves unique sequence, on any tier or chromosome.** 1.91 to 1.93 bits per base
  everywhere, coding exons included (1.9107). The whole compressible part of the genome is copies,
  tandem arrays (structural 0.41 bits per base) and old repeats.

**Stacks, bits per base, 2,467 Mb:** naive 1.7374, generic 1.6494, repeat-aware 1.6273, repeats
plus duplications 1.6061, full 1.6021.

**Efficiency, measured.**

| | genome pass, 22 chromosomes |
|---|---|
| wall clock | 16,863 s of coding, 4 h 50 m elapsed |
| seconds per Mb | 6.8 (chr18 7.0, chr21 14.0 where the duplication layer is heavy) |
| peak memory | 10.2 GB, on one chromosome at a time |
| temporary disk, peak | 1.58 GB at once for a finished chromosome |
| free disk, minimum reached | 18.8 GB against a 15 GB floor |
| read from local disk | 48.9 GB, almost all earlier chromosomes' FASTA for duplication partners |
| network | 18 requests, 4.0 MB (CpG islands) |
| model calls | 0 |

- **What each layer cost.** Unchanged in shape from chr21: duplications buy the most for the least,
  GC and CpG buy a little for 20 s a chromosome, and repeats are the most expensive layer and the
  most negative.

**chr1 and chr2 did not run. Negative, and the reason is the machine, not the method.** The pass was
interrupted at the disk floor four times, and the cost is now measured rather than estimated. chr2
is 242 Mb; counting **one** adaptive order over it consumed **17.5 GB of free space** (35.1 GB
free at 21:14, 17.6 GB when order 12 finished at 21:18), and the guard then refused order 16, which
reserves another 5.8 GB. Only 5.8 GB of that 17.5 is the count rows themselves — three rows a base,
the queries, the forward events and the reverse-strand events, at 8 bytes each; `du` attributes
under 1 GB to the spill directory and killing the process returns the rest, so the remainder is the
process's own paging onto the same volume. **chr2 would need about 38 GB free to start**, and this
460 GB machine at 94% full, with five other lanes writing, offers 30 to 35. `--rows-memory-gb`,
added for exactly this, does not help: it moves the 5.8 GB from disk into 19 GB of RAM that is
already 11 GB spoken for. Both chromosomes are skipped, the chain skips finished chromosomes so it
resumes on them when there is room, and the rollup names them as missing. Between them they hold
197 of the 1,098 blocks and 5.3 of the 32 Mb, which cannot overturn +0.0117 with an interval of
0.0095 to 0.0145.

**What is weak.**
- **The mixer.** The windowed Bayesian mixer is simple, and a stronger one (logistic mixing, as in
  GeCo3 or cmix) would lower every stack. The per-layer comparisons share the mixer, so they are
  fair to one another. The absolute numbers would move.
- **Training data.** Static models are fitted on a single other chromosome. Doubling it lowers
  every stack by 0.01 to 0.025 bits per base and turns the repeats verdict (above); nothing else
  moves.
- **The copy model** anchors on 12-mers nearest the linear expectation. It is right for most bases
  of a pair and wrong inside internal repeats.
- **One annotation code per layer.** The codes are universal rather than tuned. A better code for
  the coordinates could make repeats level over generic, but not the motif sites, whose address
  costs three times what they save.
- **Chromosomes.** chr21 and chr22 have acrocentric copies that dominate their tiers. chr18 alone
  carried the tier result until the genome-wide pass above; it now rests on 901 blocks over 22
  chromosomes, with chr1 and chr2 still missing.

## The executor test: does a block execute the value stored beside it? (2026-09-14)

Albert's third functional class is a program that reads a stored value: a unit whose read-out
changes with the value held beside it. The human panel catalogued the storage — 200-bp windows where
the 89 assemblies carry a few recurring values — and for some of those values a read-out has been
measured, as a GTEx eQTL or an MPRA allele pair. The test swaps the stored allele in silico and asks
whether the model's predicted read-out of the measured gene moves in the measured direction, more
often than the same swap at a matched unit whose value nobody has measured.

Code: `attribution/executor.py`, run by `scripts/executor_test.py`. Results:
`executor_shortlist.json` (no model call), `executor_test.json` (E1 and E2) and
`executor_holdout.json` (E3). The criterion, the alpha spent at each look, the order of the
requests, the stopping rule, the power and the ENCODE caveat are in `CRITERION`, written before the
first request; three amendments and the hold-out are in `AMENDMENTS` and `HOLD_OUT`, each dated and
each written before the numbers it governs.

**The design.** 15,927 storage units on chr21 and chr22 whose recurring value is a measured variant.
Each unit is tested at one variant, with two controls from the same chromosome within 250 kb (median
146 kb), the same GC stratum and, where the neighbourhood allows, the same replication-timing
tertile and number of recurring values, panel r-squared under 0.2 with the tested variant, and no
GTEx pair and no MPRA row anywhere in the control unit. A control borrows the unit's gene, tissue or
cell and its measured sign, so it asks how often the model's allele effect on that gene has that
sign when the variant is not a measured regulator. Both alleles go in one request over the 1 Mb
window, the read-out is the log2 fold change on the named gene, and a predicted effect under 0.001
is a no-call in either arm. Three endpoints, in the order they were run: E1 MPRA allele pairs, E2
eQTL values DAP-G fine-maps (PIP at or above 0.5), E3 the significant eQTL values DAP-G does not
fine-map.

**E1, the cleaner test, is uninformative.** 12 pairs, 2 units answered. chr21 and chr22 hold only 20
significant MPRA allele pairs in storage units, and the criterion said in advance that at that size
only a difference of 0.5 could be detected. Negative, and not evidence of absence.

**E2 is a success by the criterion.** 161 pairs: units agree 38 of 56 (0.679), matched controls 30
of 73 (0.411), a difference of **+0.268**, one-sided Fisher p 0.0021, the same direction on both
chromosomes (chr21 +0.149 at p 0.34, chr22 +0.301 at p 0.0024). It is driven by chr22.

**Then the hold-out, and it is what the test was for.** E2 alone cannot separate two readings: the
value is the cause of the read-out, or the model simply agrees with GTEx in eQTL-rich
neighbourhoods. E3's variants are significant eQTLs (median p 5e-34) that DAP-G does not fine-map:
572 of the 770 units carry a DAP-G row at the very variant with a median posterior of 0.044, so
these are values linked to a cause rather than the cause. Written before any E3 row was read: if a
unit executes its value, E3's difference should be clearly smaller than E2's; if E3 matches E2, the
difference is a property of the neighbourhood and not of causality, which would have been the more
important result. "Clearly smaller" was given a number before the rows were read: E2's difference
minus E3's with a one-sided 95% bound (`compare_endpoints`, `HOLD_OUT['read']`).

**E3 is diluted, as the executor claim predicts.** 1,415 pairs over 770 units, 2,000 requests: units
agree 239 of 450 (0.531), controls 351 of 749 (0.469), a difference of **+0.062**, p 0.021, upper
95% bound 0.111; chr21 +0.054 and chr22 +0.066, the same direction but neither alone significant. By
the criterion that is a negative endpoint — it needs 0.10 — and against E2 the gap is **+0.206 with
a one-sided 95% lower bound of +0.058**, so the gap clears zero. At the pre-registered budget of
1,000 pairs, read separately, the same: E3 +0.057 at p 0.066, gap +0.211, lower bound +0.059.

**The reading.** The hold-out went the way the executor claim predicts and not the way that would
have destroyed it. The agreement difference is carried by the values fine-mapping calls causal and
falls by three quarters at values that are only linked to a cause, in the same neighbourhoods, at
the same distances, with the same matching — which a property of eQTL-rich neighbourhoods could not
do. That is the strongest form this claim has had, and with the deletion target work it is one of
the two attribution claims still standing. What it is not: E3's difference is small and positive
rather than zero, so linked values agree a little too, as linkage alone would give; the whole test
rests on one model's predictions; and the model learned from ENCODE RNA-seq while GTEx expression is
the same kind of measurement, so agreement with an eQTL is partly the model reading back its inputs.
The clean way to break that objection was E1, and E1 has no power on two chromosomes.

**The widened E1, pre-registered and waiting on quota (2026-09-14).** E1 was the endpoint that
answers the objection to E2 and E3, and it had 12 pairs — because the shortlist was built from the
human panel's storage units and the panel has run on chr21 and chr22 only. The measurement is not so
limited. MPRAVarDB read genome-wide through the panel's own bigBed reader (`mpra_genome_wide`, 8
seconds, no model request) holds **231,336 rows over 129,442 variants; 18,038 rows are significant**
by the pre-registered rule, and **7,178 of those are in a cell line the model has RNA-seq tracks for**
(GM12878 3,169, HepG2 2,799, K562 739, Jurkat 386, HeLa 85), covering 6,196 distinct variants on 23
chromosomes. That is the widening: 20 rows became 7,178.

The plan is in `E1_WIDE` and `executor_mpra_wide.json`, written before any request.
- **Two strata, the headline named first.** **E1W_mpra_wide**, 3,260 pairs over 1,515 separate
  neighbourhoods in GM12878 and Jurkat, is the headline. **E1S_mpra_saturation**, 847 pairs, is the
  saturation-mutagenesis study alone: it mutates whole elements and its usable pairs sit at three
  loci in HepG2, so it is a within-element test, reported apart and never the headline.
- **The control is a measured null**, not an unmeasured unit: a variant the same study measured in
  the same cell line and found not significant (FDR above 0.2), no significant row anywhere, 200 bp
  to 250 kb away, nearest taken, borrowing the unit's gene, cell and measured sign. Median distance
  909 bp. Off chr21 and chr22 there is no panel to draw an unmeasured control from, and the
  same-library null is the better match anyway: same assay, same cell, same neighbourhood.
- **The read-out gene is chosen without the model**: DAP-G's gene at the variant where there is one
  (2,629 tests), else the nearest protein-coding TSS (3,567). The original E1 let the model pick the
  gene with the largest predicted effect when no eQTL named one, which favours the unit arm at the
  no-call threshold; the widened endpoint refuses the pair instead.
- **What each outcome means** is written down: a difference of 0.10 or more at p 0.01 with GM12878
  and Jurkat agreeing is the external answer to the ENCODE objection; a difference at or below zero,
  or a futility stop, says the model does not track an external measurement of direction and confines
  E2 and E3 to eQTL-shaped read-outs; anything between decides nothing.

**The direction convention, checked without the model.** The reporter's direction and the GTEx
eQTL's direction agree 0.518 of the time over the tested GM12878 variants and 0.670 where DAP-G
fine-maps them: the convention is the right way up, and the split between causal and linked values
appears in two measurements with no model between them. It is read on its own below, "Two
measurements agree without the model".

**The widened E1 ran, and by its own pre-registration it resolves nothing.** 1,000 pairs, **1,999
requests**, 33 minutes, the headline stratum in the pre-registered order (strongest measured effect
first, so every pair tested has |log2FC| at or above 1: the subset where the measurement is most
certain). Units agree **205 of 398 (0.515)**, matched nulls **191 of 381 (0.501)**, a difference of
**+0.014**, one-sided p 0.38, upper 95% bound **0.073**. No look fired: neither success at 150 or 400
pairs nor futility, and the run ended at the budget. The pre-registration named three outcomes and
this is the third: a difference between 0 and 0.10 decides nothing. It is not a confirmation, and it
is not the contradiction either — but the bound is worth stating, because it does exclude what the
endpoint was built to find: **a difference of 0.10 or more is outside the 95% bound, and E2's +0.268
is far outside it.** Whatever the model is doing when it agrees with a fine-mapped eQTL, it is not
doing it against a reporter assay's direction at this size.

- **The declared secondary failed.** If the model tracked the measurement rather than its own
  training data, its agreement should have split like the two measurements do, higher on fine-mapped
  variants. It does not: fine-mapped 17 of 33 against 12 of 18 controls (**-0.152**), the rest 188 of
  365 against 179 of 363 (**+0.022**). Only 57 of the 1,000 pairs are fine-mapped, so this is a
  failed prediction at low power rather than a measured absence, and it is recorded as a failure
  either way.
- **The arms are balanced, which the old E1's were not.** Choosing the read-out gene without the
  model left 398 units and 381 controls answered of 1,000 each, with the same median absolute
  predicted effect (0.00085 in both). The endpoint's null result is not an artefact of one arm being
  filtered harder than the other. 135 test sides returned no value at all: no track for the cell, or
  the named gene outside the window.
- **One cell line carried it.** The order by measured strength filled the budget with GM12878 (995 of
  1,000 pairs) and left Jurkat 5, so the pre-registered requirement of the same direction in both
  cells could not be tested. The saturation-mutagenesis stratum was not granted quota and has not run.
- **The comparison against E2 in the result file is descriptive only.** `compare_endpoints` was
  written for E2 against E3, two readings of the same kind of measurement, where a gap means dilution
  by linkage. E2 and this endpoint measure different outcomes, so the gap of +0.254 between them is
  not that dilution and is not read as one.
- **What it does not settle.** The reporter measures a plasmid fragment; the model predicts a
  chromosomal gene; and the two measurements agree with each other only 0.670 of the time at best
  (below). A null here is consistent with the model being wrong about direction, and equally with the
  two assays being too far apart for the question. The honest summary is that the external endpoint
  is now built, genome-wide, at 3,260 available pairs — and the first 1,000 of them did not answer.

**E2 and E3 asked again away from chr21 and chr22, pre-registered and waiting on quota
(2026-09-15).** The human panel reached 19 chromosomes overnight, so the two results this lane rests
on can be asked where the answer can come out differently. GTEx v8 was distilled over every
catalogued chromosome's storage units in 138 seconds (71.5 M pairs scanned, 8,995,460 kept inside
1,247,338 units, 1.4 GB streamed, nothing of the tar stored), and the pairs were assembled chromosome
by chromosome in 33 minutes with no model request: **E2R, 1,441 units and 2,671 matched pairs on 16
chromosomes, against the 81 units and 161 pairs the original had on two; E3R, 142,097 units and
270,558 pairs**. chr21 and chr22 are the discovery set and are excluded from the replication claim
outright; chrY has no significant GTEx pair in a storage unit and drops out.

Everything the first run fixed stays fixed — the read-out, the measured sign, the no-call at 0.001,
two matched controls per unit within 250 kb with the same GC stratum and, where the neighbourhood
allows, the same timing tertile and value count, panel r-squared under 0.2, each unit counted once.
Two things are new, both in `E2_REPLICATION` and both written before any request.

- **The order is round robin over chromosomes**, strongest evidence first inside each, so the first
  requests answer whether the effect exists away from the discovery chromosomes rather than filling
  up on whichever chromosome has the most units: the first 400 pairs carry 21 to 28 from every one of
  the 16.
- **Success needs a leave-one-out.** The pooled difference must reach 0.10 at the look's alpha *and*
  survive the removal of any single chromosome above 0.05. A pooled difference that only clears the
  bar because one chromosome carries it is reported as "not replicated as a genome-wide claim", never
  as a success. That is the instrument the four chromosome-specific readings this project found on
  the acrocentrics would have failed.

**What a replication would and would not buy.** It would say E2 is not chr21-and-chr22-shaped, and
E3R would re-ask the dilution contrast at power. It would not move the evidence outside the model:
this endpoint is one model agreeing with fine-mapping, the external endpoint (the widened E1) came
back at +0.014 with an upper bound of 0.073 and resolved nothing, and the benchmark lane has since
shown that a model producing a direction is evidence of having been scored rather than of an element
doing anything. The only part of this work whose evidence does not pass through AlphaGenome remains
the measured-against-measured reading below.

**The budget, fixed before the key:** 600 units of E2R (1,200 pairs, 1,800 requests) and 400 of E3R
(800 pairs, 1,200 requests), 3,000 requests in all, with the looks at 150 and 400 pairs and the
stopping rule unchanged.

**The replication ran: the direction holds away from chr21 and chr22, the size does not, and by the
pre-registered bar it resolves nothing (2026-09-15).** 3,000 requests exactly, under the shared-key
lock, in the round-robin order.

- **E2R, the replication.** 1,166 pairs over 16 chromosomes, 1,801 requests. Units agree **231 of 381
  (0.606)**, matched controls **255 of 495 (0.515)**: a difference of **+0.091**, one-sided p
  **0.0043**, upper 95% bound 0.146. The pre-registration asked for 0.10 and it did not get there, so
  the outcome is the third one, written before the run and not bent now: **a difference between 0 and
  0.10 resolves nothing**. What did hold is everything around the number. The difference is positive
  and significant at the final alpha, it survives **leave-one-out on every chromosome** (the lowest is
  +0.070, with chr6 removed, above the 0.05 floor the rule set), and 11 of the 16 chromosomes lean the
  same way (sign test p 0.105). Neither look fired: at 150 and at 400 pairs it was neither a success
  at its alpha nor futile.
- **The size did not replicate.** The discovery run read +0.268 on chr21 and chr22; **+0.268 is far
  outside this run's 95% interval**, whose upper bound is 0.146. The honest reading is that the first
  estimate was inflated about threefold by the two chromosomes it was found on, and that the effect
  which survives at genome scale is about a third of it and just under the bar the pre-registration
  set for calling it real. chr6 (+0.415) and chr7 (+0.383) are the two chromosomes that would have
  carried a bolder claim, and the leave-one-out is what keeps them from doing so.
- **E3R, the dilution, replicates cleanly.** 775 pairs, 1,199 requests: units 107 of 216 (0.495),
  controls 178 of 339 (0.525), a difference of **−0.030**, upper bound 0.042, which is a **futility**
  verdict by the rule — linked values show no agreement excess at all. Against E2R the gap is
  **+0.121 with a one-sided 95% lower bound of +0.031**, so the split between values fine-mapping
  calls causal and values merely linked to a cause is there again, on 16 chromosomes the discovery
  never saw. The dilution is the part of this result that replicated best.

**What it adds up to.** The executor claim's *shape* survived a genuine replication — causal values
agree, linked values do not, and no single chromosome carries it — while its *effect size* shrank by
two thirds and fell below the line drawn in advance. That is the fifth claim in this area to shrink
under a control tonight, and it is recorded as such rather than as a success. It also remains, as the
pre-registration says in the code itself, one model agreeing with fine-mapping: the endpoint whose
outcome is an external measurement resolved nothing at +0.014, and the only evidence here that does
not pass through AlphaGenome is the measured-against-measured reading below — 0.670 agreement where
DAP-G fine-maps against 0.518 overall, two instruments and no model between them, which is also the
one number the model's own 0.606 against 0.515 now sits closest to.

**Caveats kept with the result.**
- **Effect sizes are tiny.** The median absolute predicted log2 fold change is 0.0012 in units and
  0.0011 in controls, so the no-call threshold of 0.001 decides which pairs are scored at all. The
  grid is reported with every reading: at no-call 0, E2 is +0.131 and E3 +0.039; at 0.01 both
  collapse to a few dozen pairs (E3 +0.448 on 30 units and 7 controls).
- **A small contamination of E3.** 22 of its 770 units do carry a DAP-G posterior above 0.5
  somewhere, for a gene whose sign disagrees across tissues or which is not the variant's top eQTL.
  Dropping them: +0.052 at p 0.049, the same reading.
- **Two of the three endpoints are eQTL endpoints**, and both carry the ENCODE caveat.
- **The counting was corrected** after E1 and E2 were scored: a unit with two controls had been
  counted once per pair. E2's difference barely moved (+0.269 to +0.268) but its p-value fell from
  0.00036 to 0.0021, and the live E3 run's own file was rewritten from its cache under the corrected
  tally, which is what `executor_holdout.json` holds.
- **Requests.** 447 for E1 and E2 (a two-request smoke test, 308 before the control amendment, 118
  after it, 19 for the last 23 pairs) and 2,000 for E3, all through the existing adapter path with a
  1 Mb window and both alleles in one request. The rerun that produced the committed hold-out
  numbers spent none: every variant was already in the cache.

## Two measurements agree without the model, and they split the same way (2026-09-15)

Every claim in this area rests on one model's predictions, and the model learned from ENCODE tracks
that are kin to the outcomes it is asked about. This reading does not. It compares two measurements
of the same alleles, made by different laboratories with different instruments, with nothing between
them.

- **A reporter assay.** MPRAVarDB's log2 fold change for the alternative allele over the reference: a
  few hundred bases of human sequence on a plasmid, the two alleles side by side, transcripts counted.
- **An expression quantitative trait locus.** GTEx v8's slope for the same alternative allele: the
  gene's expression in hundreds of donors, against the allele they carry in their own chromosomes.
- **Fine-mapping to separate cause from company.** DAP-G's posterior at that very variant: at or
  above 0.5 the allele is probably the cause; below it, usually a passenger of one.

It was first read on the 179 fine-mapped variants the widened E1 happened to carry in GM12878, where
it gave 0.670 against 0.518. `TWO_INSTRUMENTS`, written before the genome-wide table was built, asked
it at scale: every variant with a significant MPRA row and a significant GTEx pair at the same
position with the same alleles, in every cell line, the eQTL taken from DAP-G where it fine-maps one,
reported per cell line and per study with the sample size beside every figure, against a control of
variants the same libraries measured and found not significant. GTEx was distilled at the MPRA-tested
positions rather than over storage units, so nothing here depends on the human panel either. 18
minutes, no model request.

**The split is real and the control is clean.** Of 5,252 rows both instruments measured (309 distinct
fine-mapped variants), the fine-mapped stratum agrees **221 of 350, 0.631**, against **2,525 of 4,902,
0.515** for the rest: a difference of **+0.116 at one-sided p 1.5e-5**. The control — 58,302 rows
whose MPRA measurement was not significant — agrees **0.502 in the fine-mapped stratum and 0.502 in
the rest, a difference of exactly 0.000**. Whatever produces the split, it is not the way the data
were assembled: when one instrument measured nothing, the two agree at chance in both strata.

**But one cell line carries it, and by the pre-registration that means it resolves nothing.** The
reading asked for four things and got three: the fine-mapped stratum above 0.60, a margin of at least
0.08 at p 0.01 or better, and a control at chance. The fourth was the same direction in the two
largest cell lines separately, and that fails. **GM12878 reads 0.668 against 0.508 (p 1.2e-7) on 289
fine-mapped rows; every other cell line together reads 0.459 against 0.524 on 61.** No second cell
line has more than 15 fine-mapped variants, so this is one cell line's result with no corroboration
rather than a contradicted one — the same shape as a claim carried by one chromosome, and it is
reported the way that would be.

**Why it is probably a cell type and not an accident.** Directions can only correspond if the
reporter's cell and the eQTL's tissue are the same tissue, and GM12878 against EBV-transformed
lymphocytes is the one such pair that exists in bulk. Where the cell and the tissue do match, the
fine-mapped stratum agrees **29 of 33, 0.879**; where they do not, **0.606**. By study, Tewhey's
lymphoblastoid library reads 0.754 of 138 and Abell's 0.592 of 147, while every non-lymphoblastoid
study has fewer than 16 fine-mapped rows and no consistent direction. The claim that survives is
therefore narrower than the one that was asked: **two instruments agree about which way a causal
allele pushes, in lymphoblastoid cells, where both are measuring the same tissue** — 0.631 pooled,
0.668 in GM12878, 0.879 where the tissue matches, against 0.515 for alleles that are only linked to a
cause.

**The 0.879 is the one figure here whose control does not behave, and it is stated as such
(2026-09-16).** The control is at chance overall — 0.502 fine-mapped against 0.502 for the rest, a
difference of 0.000 over 58,302 rows — but read inside the matched-tissue stratum alone it splits the
same way the result does: **65 of 104, 0.625, against 0.509 for the rest** (`control_not_significant`,
`matched_tissue_only`, in the result). Alleles whose reporter measurement was *not* significant should
carry no direction to agree about, so a split there says something in the matched-tissue overlap
favours agreement regardless of the measurement — 104 rows of it, against the 33 that carry the 0.879.
The pooled 0.631 and the GM12878 0.668 keep their clean control and are unaffected; the matched-tissue
figure is not evidence on its own until that control is understood.

**What it cannot mean**, kept from the pre-registration.
- **The overlap is not a sample of the genome.** MPRA libraries are built from eQTL and GWAS
  variants, so these are alleles already suspected of doing something; no statement about the genome
  follows.
- **A plasmid fragment is not a gene in its chromosome.** Agreement says the allele does something in
  both assays and in the same direction, not that the eQTL's gene is the reporter element's target.
- **Fine-mapping is itself an instrument.** The split is between what DAP-G calls causal and what it
  does not.
- **One cell type**, and 350 fine-mapped rows over 309 variants carry the whole reading.

It remains the only evidence in this area that does not pass through AlphaGenome, it is now measured
at twice the size with a null control that behaves, and it is smaller and narrower than the first
reading suggested.

## The all-elements sweep folded: what 605,137 deleted elements say, and what they do not (2026-09-15)

The chain (`scripts/enhancer_targets_all_chain.py`) is deleting every registry element of every
chromosome in AlphaGenome, one chromosome at a time, and each finished chromosome commits its own
summary. `scripts/enhancer_targets_all_genome_wide.py` folds the complete ones into one reading and is
re-run as each lands; every reading it has taken is kept in the result's `history`, so what follows can
be checked against what the same fold said one chromosome ago. Eighteen are complete (chr2, chr8 to
chr22, chrX and chrY), 605,137 elements, 425,526 requests.

| reading | genome-wide, 18 chromosomes | across chromosomes | the control it has to be read against |
|---|---|---|---|
| names a gene at all | 64.1% | 58.1% (chr13) to 70.8% (chr19), median 62.8% | 87% of matched random windows in the locus benchmark |
| the named coding gene is the nearest TSS | 62.0% | 48.8% (chrY) to 65.3% (chr22), median 61.9% | none measured; the nearest-TSS heuristic itself scores 8 of 12 published loci |
| the named coding gene is inside the element's own CTCF node | 77.2% | 68.4% (chr13) to 91.2% (chr19), median 75.1% | 74.3% for as many boundaries placed at random, measured on these same elements. *(2026-09-27: that random set is matched on boundary count only, not on the 50 kb floor every real node passes, so it is biased low; across four baselines the excess runs +1.2 to +6.6, and on 661 measured CRISPRi pairs it is +5.89. NODES-READER-WRITER.md, the random-boundary audit.)* |

**The node control, re-measured where the claim is now being made.** dbad593 qualified the old
"90.2% of enhancers act inside their own CTCF node" to 81.7% against 79.1% for randomly placed
boundaries — 2.6 points — over the archive as it stood. That control belonged to the archive's element
set, not to this one, so the fold now places its own: the same 176,650 elements naming a coding gene,
the same most-moved gene, the same `infer_domains` caller, against as many boundaries drawn uniformly
at random (200 draws, seed 7); at eighteen chromosomes that is 279,231 elements. It reproduces the
scorer's own containment figure to four decimals on every chromosome, which is the check that the two
readings are the same reading.

Pooled, the node keeps **77.6% against 74.9% at random: +2.73 points**, so dbad593's correction holds
at this size. But the per-chromosome control, once it was asked how big its own null's scatter is,
says something that supersedes most of what was written here at fourteen and fifteen chromosomes.

| | chr10 | chr18 | chr11 | chr17 | chr2 | chr12 | chr16 | chr14 | chr9 | chr19 | chr20 | chrX | chr22 | chr15 | chr13 | chr21 | chrY |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| excess, points | +7.2 | +6.9 | +6.2 | +4.1 | +3.8 | +3.5 | +1.8 | +1.7 | +1.1 | +0.8 | +0.7 | −0.1 | −0.7 | −1.0 | −1.4 | −2.2 | −22.1 |
| the null's own scatter | 1.6 | 2.4 | 1.4 | 1.2 | 1.2 | 1.3 | 1.4 | 1.8 | 1.5 | 0.8 | 1.6 | 1.4 | 1.8 | 1.9 | 2.2 | 2.9 | 6.2 |
| p | **.01** | **.01** | **.01** | **.01** | **.01** | **.01** | .23 | .35 | .51 | .36 | .70 | .94 | .67 | .57 | .46 | .45 | **.02** |

**Ten of the seventeen chromosomes have no measurable excess at all.** Against 200 random-boundary
draws rather than 20, six chromosomes beat their null (chr10, chr18, chr11, chr17, chr2, chr12, all
p 0.01) and one loses to it (chrY, p 0.02). Everything between +1.8 and −2.2 points sits inside the
scatter the null produces on its own, which on these chromosomes is 1.2 to 2.9 points wide.

This retracts the strongest form of the claim made earlier in this section. "The node is ahead on ten
and behind on five" was reporting the sign of numbers most of which are zero: **every chromosome the
node appeared to lose on except chrY is inside the null**, chr21 and chr22 included. The honest
statement is that the node beats a random partition of the same resolution on six chromosomes, loses on
one, and on ten the measurement cannot tell. That is a smaller claim than "the excess changes sign
between chromosomes", and it is the one the evidence supports.

**What survives is a clean line, and it is not about sign.** Sorting the eighteen by how finely the
caller cuts them separates the measurable from the unmeasurable almost perfectly:

| at or above 6.55 boundaries per Mb | chr8 7.65 | chr20 7.40 | chr11 7.00 | chr18 6.96 | chr12 6.91 | chr10 6.88 | chr2 6.83 | chr17 6.55 |
|---|---|---|---|---|---|---|---|---|
| | **+4.6** | +0.7 | **+6.2** | **+6.9** | **+3.5** | **+7.2** | **+3.8** | **+4.1** |
| p | **.01** | .70 | **.01** | **.01** | **.01** | **.01** | **.01** | **.01** |

Seven of the eight at or above 6.55 per Mb carry an excess distinguishable from their null; none of the
nine autosomes below it does. The two exceptions are chr20 at 7.40, above the line and unmeasurable,
and chrY at 0.84, far below it and the one clearly measurable *negative* in the set — 48 boundaries
over 57 Mb, nearly all in the short euchromatic arm, so uniformly placed boundaries leave the
element-dense region uncut and beat the real nodes by 22 points on 125 elements. Sixteen of eighteen
fall on the right side of the line.

**A second rule, fixed now, because the first one failed for a reason this one avoids.** The sign rule
called the *direction* of the excess and came unstuck on chrX, where there was no direction to call.
This one claims only *where an excess exists at all*, which the p-values say is the answerable
question: **an excess is measurable if and only if the chromosome is cut at 6.55 boundaries per Mb or
finer.** It is fitted on these eighteen, scores 16 of 18 with both misses named above, and says nothing
about sign. Every one of the six chromosomes still to land is cut at 6.58 or finer (chr3 7.26, chr7
7.10, chr4 6.94, chr6 6.90, chr5 6.85, chr1 6.58), so it predicts all six measurable and has six wrong
calls available to it. It is recorded in `MEASURABILITY_RULE` with its own scoreboard, and the sign
rule is **not** amended by it and keeps its failure.

**The pooled excess drifts in both directions, and that is still the finding.** It has read +1.85,
+2.47, +3.02, +2.85, +2.57 and +2.73 points at twelve through seventeen chromosomes — up as finely cut
chromosomes arrived, down as chr9 and chrX did — while touching none of the per-chromosome numbers
underneath. Now that most of those per-chromosome numbers are known to be zero, the reading is stronger
than before: **the genome-wide excess is an average over six chromosomes that carry a real effect and
ten that carry noise**, and its value depends on the mix. No single figure for it is a property of the
genome — dbad593's 2.6 included, and this paragraph's +2.73.

Nor does the excess follow the density relation that the containment level follows. chr19 is the
densest chromosome and highest on containment at 91.2%, and its excess over random is +0.7, among the
smallest; chr18 is sparse, low on containment at 69.1%, and carries the largest excess at +7.5. Being
easy to contain and being better-than-random at containing are different things.

**What the sign is a property of.** Eight candidate accounts were measured on the fourteen, each
computable from what the control already loads, and ranked against the excess:

| | boundaries per Mb | median node | node length CV | element-target span | coding genes per Mb | elements per Mb | boundaries follow the elements | the null's own level |
|---|---|---|---|---|---|---|---|---|
| Spearman with the excess | **+0.79** | −0.69 | −0.46 | −0.20 | +0.35 | +0.45 | −0.59 | −0.46 |
| with chrY dropped | **+0.75** | −0.62 | −0.41 | −0.03 | +0.21 | +0.33 | −0.50 | −0.34 |

Five chromosomes on, the boundary-density account has held its ground while the alternatives have
decayed: without chrY, gene density has stayed in the +0.15 to +0.32 band while boundary density has
held +0.74 to +0.77 throughout. But this whole table has to be read knowing that ten of the seventeen
excesses it correlates against are indistinguishable from zero, so most of what it ranks is noise. The
part that is not noise needs no correlation to state: **the six chromosomes whose excess beats its null
all carry between 6.55 and 7.00 boundaries per Mb**, and the one whose excess loses to its null carries
0.84. Fine cutting is necessary for a measurable excess in this sample and is not sufficient — chr20 at
7.40 per Mb, the finest cut of all, is not distinguishable from random.

One thing tracks the sign, and it is how finely the caller cuts. Every chromosome the node wins on
carries 5.77 boundaries per Mb or more; every one it loses on carries 5.78 or fewer. The median node
length is the same quantity inverted, and the rest fade when chrY is dropped. The mechanism that would
have been the interesting answer — that real CTCF boundaries crowd into the element-dense stretches and
so spend their cuts where the pairs are, while uniformly scattered boundaries waste most of themselves
on empty sequence — is there in the sign (−0.55, −0.41 with chrY out) but is weaker and partly the same
measurement. chrY is that mechanism at its limit: 48 boundaries over 57 Mb, essentially all of them in
the short arm, so random placement leaves the element-dense region uncut and beats the real nodes by 22
points.

This makes the node claim a statement about resolution rather than about CTCF. Where CTCF-only elements
are dense enough to cut the chromosome finely, the partition they produce holds element and target
together better than chance; where they are sparse, the nodes carry nothing that as many boundaries
thrown at the chromosome would not carry, and on chr13, chr15, chr21, chr22 and chrY they carry less.

**A rule fixed before the evidence, so it can fail.** The relation is fitted, not derived, so it is
recorded as a rule with a threshold and scored out of sample. Fixed on 2026-09-15 from these twelve:
*the node beats as many random boundaries above 5.8 boundaries per Mb and loses below it; within 0.3
per Mb of the line the call is refused.* In sample it separates 11 of 12, its one error being chr19 at
5.77 per Mb, predicted to lose and winning by 0.7 points. Boundary density needs no element table, so
the standing predictions for every chromosome the sweep has not reached could be computed and committed
before their elements existed:

| | chr8 | chr3 | chr7 | *chr11* | chr4 | chr6 | *chr10* | chr5 | *chr2* | chr1 | *chr9* | *chrX* |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| boundaries per Mb | 7.65 | 7.26 | 7.10 | *7.00* | 6.94 | 6.90 | *6.88* | 6.85 | *6.83* | 6.58 | *5.87* | *5.26* |
| predicted | + | + | + | *+* | + | + | *+* | + | *+* | + | *refused* | *−* |

*Italics mark the chromosomes that have since landed: chr11 +6.2 (right), chr10 +7.2 (right), chr2 +3.8
(right), chr9 +1.1 (refused), chrX +0.0007 (**wrong**). Seven remain, all predicted positive, and the
rule stands as written.*

So the rule stakes something: every remaining autosome should show a positive excess, chr9 is too close
to call and is refused rather than guessed, and **chrX alone should come back negative**. The fold
scores each chromosome against this table as it lands and records right, wrong and refused, so a
failure is visible in the committed result rather than absorbed.

**Held out so far: two right, none wrong, none refused.** chr11 finished at 00:24 and chr10 at 03:09
on 2026-09-15, both after the rule was fixed and committed, so neither was fitted on. The committed
result tags every chromosome with one of five words — `held out: right`, `held out: WRONG`,
`held out: refused`, `fitted on, not a test`, `still to land` — and carries a one-line scoreboard and a
verdict, so which chromosomes were fitted and which were tests can be read without reconstructing
either set.

| | predicted | boundaries per Mb | measured excess | |
|---|---|---|---|---|
| chr11 | positive | 7.00 | +6.2 points, p 0.01 | held |
| chr10 | positive | 6.88 | +7.2 points, p 0.01 | held |
| chr9 | *refused, too close to call* | 5.87 | +1.1 points, p 0.51 | refused, and rightly |
| **chrX** | **negative** | **5.26** | **+0.0007, p 0.94** | **FAILED** |
| chr2 | positive | 6.83 | +3.8 points, p 0.01 | held |
| chr8 | positive | 7.65 | +4.6 points, p 0.01 | held |

Held out: four right, one wrong, one refused. Every chromosome the rule has called positive and been
right about carries an excess distinguishable from its null, and the one it got wrong is one where
there is no effect to get right. That is the shape of a rule that has found something narrow and states
it too broadly, which is what the second rule below is for.

**The rule failed on chrX, and the failure stands.** chrX was the one chromosome the rule called
negative and the only standing prediction that could falsify it. It came back at +0.0007 — seven
ten-thousandths, positive — so the rule is wrong, and the committed result's `verdict` field reads
`FAILED on chrX`.

The failure has a tail that matters more than the failure. Asking the control for its null's scatter
was what produced the table above, and under the better 200-draw null chrX reads **−0.0007**: the same
chromosome, the opposite sign, and the rule would score right. The quantity is zero (p 0.94, scatter
1.4 points) and its sign belongs to the random seed. Converting a failure into a pass by improving the
instrument an hour after the failure is exactly the move this section promised not to make, so the
verdicts on held-out chromosomes are now frozen in `HELD_OUT_LOG` as they were measured when each
landed, the result records that chrX "recomputes differently now", and a test asserts the chrX verdict
stays `WRONG`. Rescuing it would have been indistinguishable from refitting.

What the failure teaches is not that the threshold is in the wrong place. It is that the rule made a
confident sign call on a chromosome where the sign is not measurable, and the ±0.3 per Mb refusal band
was too narrow to catch it — chrX sits 0.54 below the line and is still at zero. **The band is not
being widened.** Widening it now would be fitted to the chromosome that broke the rule, and the eight
chromosomes still to land have to test the rule as written.

chr10's is the largest excess of the fifteen and chr11's the third largest, from the third and fourth
finest boundary spacing. Both are the easy half of the rule's range, and two chromosomes drawn from the
half where it is confident is weak evidence for the rule and no evidence at all for the threshold.

**The refusal band is calibrated, as far as it goes.** It was drawn as a hedge — within 0.3 per Mb of
the line, decline to call — and the four chromosomes inside it, chr14, chr9, chr19 and chr13, have
excesses averaging 1.4 points in magnitude against 5.2 for the twelve outside, and not one of the four
is distinguishable from its null. The band does mark where the quantity vanishes. Its failure is that
the region where the quantity vanishes is far wider than the band: ten chromosomes are indistinguishable
from random and only four of them are inside it.

What the held-out chromosomes support meanwhile is the rule's account of the pooled number. The excess
rose +1.85 → +2.47 → +3.02 as two finely cut chromosomes arrived and fell to +2.85 when chr9's +1.1
joined, which is what a resolution account predicts and a property-of-CTCF account does not.

The naming rate is the same kind of number. 65.4% of elements name some gene, and "some derived layer
names a target" already fires at 52 of 60 matched random windows in the locus benchmark (303d0b0,
ef80075). The two instruments are not the same — the benchmark pools every layer over a whole window,
this is one element's own deletion — but the direction is unambiguous: naming a target is what this
kind of reading does almost everywhere, and a naming rate quoted without that baseline says nothing.

**What does vary, and it is not chromosome size.** Against hg38 chromosome length the three rates have
Spearman −0.05, −0.20 and −0.04 over the fourteen. Against registry elements per Mb — the element
density this run already holds, standing in for gene density — they are +0.67, +0.74 and +0.86. Length
correlates weakly and unstably (−0.07, +0.11, −0.05 at twelve, thirteen and fourteen chromosomes)
because the dense chromosomes happen to be the short ones, while the density relation has stayed put
for containment and nearest-TSS agreement. The reading is that these rates measure how much gene there
is to name near an element, not how well the genome is being understood, and a target within a node is
easier to find where targets are dense.

**The density relation, scored out of sample, fails for naming.** The correlation is a statement about
ranking fourteen chromosomes; what it claims in practice is that an arriving chromosome's rates can be
predicted from its element density. That is testable with no new data: fit each rate on density over
the twelve the relation was read from (chrY dropped as the outlier it is), then predict chr11 and
chr10, which arrived afterwards.

| | chr11 pred | chr11 actual | chr10 pred | chr10 actual | chr9 pred | chr9 actual |
|---|---|---|---|---|---|---|
| names a gene | 65.2% | **69.4%** | 65.1% | **60.9%** | 63.5% | 62.8% |
| nearest TSS | 61.6% | 62.5% | 61.5% | 62.9% | 60.3% | 61.9% |
| inside the node | 78.1% | 78.8% | 77.9% | 74.4% | 74.2% | 75.0% |

chr11 and chr10 have almost the same element density — 358 and 354 per Mb — so the relation predicts
almost the same naming rate for both, 65.2% and 65.1%. They came back **8.5 points apart**, at 69.4%
and 60.9%, residuals +4.2 and −4.2 points: as large as most of the spread the relation was supposed to
explain, and in opposite directions from two nearly identical inputs. chr9 then fitted it well, at
−0.7, +1.6 and +0.8 points on the three rates. Three predictions in, the naming errors are −0.7, +4.2
and −4.2, a mean absolute error of 3.0 points against a total spread of 12.7.

So the corrected reading: element density ranks chromosomes on all three rates and predicts an
individual chromosome's naming rate only sometimes, missing by a third of the whole spread when it
misses. Its correlation with naming fell from +0.75 to +0.67 when chr10 arrived and has recovered only
to +0.69. Something other than how much gene is nearby sets how often an element names one, and chr10
is where to look for it: the lowest naming rate of any autosome bar the sparsest, and simultaneously
the largest node excess of the fifteen. Whatever it is suppresses naming without touching placement.

**The acrocentric trap did not spring here.** chr21, chr22 and chrY pooled give 65.7% / 63.7% / 80.0%
against 65.4% / 61.6% / 78.4% for the other eleven, so unlike the compression lane's per-chromosome
layers, this reading was not bought by the small chromosomes. The standing prediction is the density
one: chrX, chr4 and chr5 are sparse, and if the density relation held at the pooled level the rates
should fall as they land. Three arrivals in, it has one point in its favour. chr11 and chr10, at 358
and 354 elements per Mb, are each about a tenth denser than the 323.6 per Mb of the twelve they joined
and should have nudged the pooled rates up; naming went 65.4 → 65.9 → 65.4%, up once and down once.
chr9 at 301 per Mb is the first arrival *sparser* than the pool it joined (330.6), and the pooled rates
duly fell, naming to 65.1% and containment 78.5 → 78.1%. chrX then landed at 147 elements per Mb, less
than half the pool's density and by far the sparsest arrival yet, and the pooled naming rate fell again
to 64.9%; chr2 at 329 followed and it fell again, to 64.4%. So the pooled form of the prediction has
three observations in its favour, all small, because one chromosome can only move a pooled rate by a
few tenths.

The per-chromosome test has now run five times, fitting on the original twelve and predicting each
arrival's own naming rate from its density: errors of +4.2 (chr11), −4.2 (chr10), −0.7 (chr9), +1.6
(chrX) and −3.2 (chr2), a mean absolute error of **2.8 points** against a spread of 12.7 across
chromosomes. Density carries something — it is not predicting at random — but it misses by a quarter of
the whole spread on average and by a third at its worst, and it cannot tell chr11 from chr10 at all.
chrX is its best case: sparsest arrival, second-lowest naming rate, right for the stated reason. chr4
and chr5 remain to come.

`data/results/enhancer_targets_all_genome_wide.json`; the per-element tables stay local under
`data/knowledge/alphagenome/all_elements/<chrom>.json`. Cells K562, HepG2, GM12878 and IMR-90.
Evidence: predicted (AlphaGenome deletion effect per element, per gene, per cell line) over inferred
(CTCF-only nodes). The node control re-runs with the fold, so its sign per chromosome is re-measured
every time a chromosome lands rather than being quoted from this table.

## Measured perturbations: the deletion held against CRISPRi screens (2026-09-16)

> **Annotation, 2026-09-22.** Every figure in this section is the 2026-09-16 record and is left
> exactly as it was scored. It was scored with `deletion_drop` read from the compact
> one-target-per-element table, which gave a structural zero to 97.7% of the held-out K562 pairs.
> The feature now reads the sweep's own per-element cache instead, and the held-out K562 gain is
> +0.141 (+0.082 to +0.231) rather than the +0.083 below, while GM12878's falls to +0.035. The
> re-scored figures and what moved are in *The deletion feature read from the sweep's own cache, not
> from the compact table (2026-09-22)*, near the end of this document; `data/results/crispri_benchmark.json`
> now holds that later run.

Every enhancer-to-gene target so far came from a predicted deletion or the nearest TSS in the
node; none had been held against an element silenced in a cell and its genes measured. The ENCODE
enhancer-gene benchmark (EngreitzLab/CRISPR_comparison, Gschwind et al. 2025) is that
measurement: 10,356 valid K562 pairs from Nasser 2021, Gasperini 2019 and Schraivogel 2020 for
training, 4,378 held-out pairs in five cell types. `scripts/crispri_score.py`
(`attribution/crispri.py`, `crispri_benchmark`, 14 s, no model request) joins each pair to the
all-enhancer deletion table by overlap: 9,237 training pairs (451 regulated) and 3,849 held-out
pairs sit on a deleted element.

Three predictors on the same pairs. Activity is sqrt(DNase × H3K27ac) measured in the screen's
own cell, from the benchmark's columns; contact is 1/distance, so no Hi-C enters. The deletion is
the predicted drop in the screen's cell line when the pair's gene is the element's top target,
zero otherwise, because the table keeps one gene per element.

| K562 training, same 9,237 pairs | AUPRC | AUROC |
|---|---|---|
| distance alone | 0.441 | 0.896 |
| activity over distance | 0.519 | 0.925 |
| deletion alone | 0.461 | 0.707 |
| logistic, activity + distance, leave-chromosome-out | 0.511 | 0.922 |
| logistic, activity + distance + deletion, leave-chromosome-out | **0.674** | 0.934 |

The deletion adds +0.163 AUPRC (95% interval +0.121 to +0.210, chromosomes resampled). The model
set was chosen after the single predictors had been read, so the claim was then fixed in the code
(`PREREGISTERED`) before the held-out pairs were scored: fitted on the K562 training pairs, the
model with the deletion has the higher AUPRC on held-out K562 and on held-out GM12878.

| held out | pairs (regulated) | activity + distance | + deletion | gain (95%) |
|---|---|---|---|---|
| K562 | 1,744 (114) | 0.550 | 0.633 | +0.083 (+0.031 to +0.166) |
| K562, elements not in training | 1,580 (91) | 0.497 | 0.587 | +0.090 (+0.037 to +0.173) |
| GM12878 | 62 (14) | 0.865 | 0.962 | +0.097 (0.000 to +0.227) |
| HCT116, Jurkat, WTC11 | | refused: no AlphaGenome line in the table | | |

**Passed**, with the caveat that GM12878 rests on 14 regulated pairs and its interval touches zero.

**Coverage on both arms, after genomeos-8a's denominator lesson (`coverage_by_arm`).** The sweep
deleted registry elements inside a CTCF node, so coverage selects on the element, and it is
all-or-nothing: of 3,941 training elements, 3,472 are covered and none partly. The arms are not
covered equally — regulated pairs 95.8% against 88.9% not regulated (held out: 95.8% against
87.6%) — so the covered subset holds a slightly higher share of positives (4.88% against 4.55%).
It is not, however, an easier subset: both baselines read the same on all valid pairs as on the
covered ones (distance 0.438 against 0.441, activity over distance 0.517 against 0.519; held out
0.320 against 0.320 and 0.519 against 0.524). The comparison between predictors is fair whatever
the coverage, because they are scored on the same pairs; what the coverage governs is what the
subset is a sample of, which is cCREs inside a node near a tested gene, not enhancers at large.

**What the deletion adds is which elements, not which gene.** Given one call per element, the
deletion calls only where its top target was tested and the predicted drop clears a threshold. At
0.1 it makes 124 calls on 3,472 elements and 119 are right (96%); at 0.2, 84 of 84. On those same
elements the closest tested gene is right 117 and 83 times. Choosing the gene is mostly distance on
this set (activity is one value per element, so within an element it picks the closest gene too). The deletion's contribution is high precision about which elements act at all, at low
recall: at 0.1 it reaches 30% of the 400 elements with a regulated gene. The nearest TSS inside the
CTCF node names a tested gene on 488 elements and is right on 223 (46%). The node boundary costs
recall here without buying precision.

**Limits.** One gene per element was kept by the sweep, so the deletion cannot rank a second
target; a full per-gene table would cost the sweep again. *(Annotation, 2026-09-22: this sentence
is wrong in its second half and is the reason for the section at the end of this document. One gene
per element was kept by the compact table, not by the sweep; the sweep's per-element cache holds
every gene in the 1 Mb window, so the full per-gene table cost nothing — it was already on disk.)*
The benchmark's activity columns are
measured in the screen's cell, while the table's per-element DNase is not used. The screens test
elements chosen near expressed genes, mostly in K562.

**Next, proposed as roadmap rows** (for the roadmap editor to fold; not written into ROADMAP.md):

1. Area I: measured contact. Replace 1/distance with 4DN or ENCODE Hi-C contact in K562 and
   GM12878 and rerun the same pre-registration; the question is whether measured contact moves
   the activity model the way the deletion did.
2. Area A: a methylation state machine in BioVM, one CpG with a de novo writer (DNMT3A/B),
   maintenance at division (UHRF1 and DNMT1 with a fidelity below 1) and an eraser (TET1-3). Its
   falsifiable prediction is the solo-WCGW loss already measured in K562, HepG2 and GM12878
   (NODES-READER-WRITER.md, the Alu lead), and its link is the twin's epigenetic age.
3. Area F: synthetic lethality from a tumour's own losses against DepMap dependencies, with
   BRCA1/2 to PARP1 as the control that must be recovered; and marker selectivity as A and B and
   not C over the healthy-tissue atlas. Target discovery and selectivity only, inside the
   non-goals: no payload, construct or claim that a therapy works.
4. Area J: motif spacing and orientation against lentiMPRA activity, beyond motif counts, so a
   grammar rule is a test that can fail.

## Measured 3D contact in place of 1/distance: the substitution fails (2026-09-16)

The contact term above was 1/distance, a power law standing in for a measurement that exists. 4DN
publishes the matrices, so the substitution can be made and the same pre-registration rerun.
`scripts/crispri_contact.py` (`attribution/crispri.py` with the new `genome/hic_contact.py`, result
`crispri_contact`, no model request) reads the contact between the element's bin and the TSS's bin
from a released in-situ Hi-C .hic of the screen's own cell line at 5 kb: K562 4DNFITUOMFUQ (6.7 GB,
merged replicates, MboI, Aiden lab) and GM12878 4DNFI1UEG1HD (22.6 GB, same protocol). Neither file
is downloaded. A .hic is an addressable container — header, master index of compressed blocks,
balancing vectors, expected-value vectors — so the reader asks for the ranges a question needs over
HTTP, decompresses each block, keeps the cells it came for and drops the rest: the 10,613 distinct
K562 bin pairs the benchmark asks about lie in the 534 blocks the file's index offers for them,
132 MB of its 6.7 GB, and what remains on disk is 1.6 MB of answers (data/knowledge/hic_contact,
local, never committed; the whole first pass over both matrices took 19 minutes).
The reader is standard library only (struct, zlib, urllib), version 8 and version 9. It was checked
against the file itself before anything was scored: a 10 kb raw count equals the four 5 kb counts
inside it, cell by cell, and the file's expected-value vector reproduces the distance decay of its
own blocks. Normalisation is the file's own balancing vector, the first of SCALE, KR, VC_SQRT and VC
that defines both bins — KR for 9,127 training pairs, VC_SQRT for 1,139 where KR left a bin
undefined, and 90 pairs no vector reaches, which are reported, never scored as a contact of zero.

**Coverage, both arms, before any lift is read** (`contact_coverage_by_arm`). The matrix reaches
451 of 451 regulated training pairs (100%) and 8,714 of 8,786 not regulated (99.2%); held-out K562
114 of 114 and 1,625 of 1,630. Both baselines are unmoved by the subset (distance 0.4413 on all
9,237 covered pairs against 0.4424 on the 9,165 that also have a contact; activity over distance
0.5190 against 0.5195), so this is the same problem as the section above, not an easier one. 1,158
of the 10,266 measured pairs have a contact of exactly zero at 5 kb — a sparse cell in a deep
matrix, kept as the zero it is — and 53 training pairs have element and TSS in one bin, where the
number is the bin's diagonal, its self-contact.

The pre-registration, fixed in `PREREGISTERED_CONTACT` before the held-out pairs were scored: fitted
on the K562 training pairs that have a measured contact, activity × measured contact has a higher
AUPRC than activity + distance on held-out K562.

| K562 training, same 9,165 pairs (the logistic rows leave-chromosome-out) | AUPRC | AUROC |
|---|---|---|
| distance alone | 0.442 | 0.897 |
| measured contact alone | 0.379 | 0.891 |
| measured contact over expected alone | 0.083 | 0.668 |
| activity × measured contact (unfitted product) | 0.441 | 0.914 |
| logistic, activity + distance | **0.512** | 0.922 |
| logistic, activity × contact | 0.431 | 0.910 |
| logistic, activity + distance + contact | 0.522 | 0.928 |
| logistic, activity + distance + deletion | **0.675** | 0.934 |
| logistic, activity + contact + deletion | 0.650 | 0.928 |

On training the substitution loses: −0.081 AUPRC (95% −0.118 to −0.011, chromosomes resampled).
Replacing distance with measured contact inside the deletion model loses too, −0.025 (−0.048 to
+0.008). Adding contact to distance rather than replacing it gains +0.010 (−0.010 to +0.027).

| held out | pairs (regulated) | activity + distance | activity × contact | gain (95%) |
|---|---|---|---|---|
| K562 | 1,739 (114) | 0.550 | 0.563 | +0.012 (−0.024 to +0.046) |
| GM12878 | 61 (13) | 0.897 | 0.922 | +0.025 (−0.029 to +0.067) |
| HCT116, Jurkat, WTC11 | | refused: no 4DN matrix chosen for these cell types | | |

**The pre-registered claim passes on the point estimate and means nothing on its own.** Held-out
K562 gives +0.012 AUPRC with an interval straddling zero, while the same comparison on five times
as many training pairs is clearly negative. Read together the answer is no: measured contact at
5 kb does not replace 1/distance, and nothing here resembles what the deletion did (+0.163 on
training, +0.083 held out, both intervals clear of zero). The closest thing to a positive is contact
as an extra feature beside distance, +0.023 on held-out K562 (0.000 to +0.063) and +0.010 on
training (−0.010 to +0.027), and with the deletion 0.633 → 0.644 (+0.010, −0.010 to +0.039). Every
one of those intervals touches or crosses zero. GM12878 rests on 13 regulated pairs and decides
nothing either way.

**Why, as far as the numbers say.** Raw contact at these distances is mostly the distance decay
itself: contact alone (0.379) reads close to distance alone (0.442), and once the decay is divided
out, observed over expected alone collapses to 0.083 AUPRC at a base rate of 4.9%. The benchmark's
pairs are close together — that is what CRISPRi can test — and there the power law is already a
good estimate of contact, while a 5 kb Hi-C cell is a noisy one. The model that fits both terms
reads it the same way: given distance, it puts weight on observed over expected (+2.46) and a
negative one on raw contact (−0.61), using the measurement only as a residual.

**Limits.** One matrix per cell line, MboI in situ Hi-C, at 5 kb: a Micro-C matrix, 1 kb bins, or a
loop-call feature instead of a single cell might each behave differently, and this result does not
speak for them. The element's midpoint and the benchmark's `startTSS` define the two bins, so a pair
whose enhancer or promoter straddles a bin edge is read one bin off. GM12878's held-out pairs are too
few to test anything. What is settled is narrower than the question: this contact, at this
resolution, does not improve this model.

## A measured confidence for a predicted target, and the prevalence that breaks it (2026-09-17)

> **Annotation, 2026-09-22.** Every number in this section stands as measured and none of it is
> rewritten; what follows says what was later found to be conditional on the compact table's one-gene
> projection. The `deletion_drop` feature below is a structural zero on 8,589 of the 8,796 training
> pairs (97.6%) and 1,680 of the 1,715 held-out ones (97.9%), because the table keeps one gene per
> element — so the magnitude is only ever spoken on the 245 pairs of the predicted-target row. Asked
> through `attribution/targets.ElementResponses`, two thirds of the excluded pairs carry a number the
> sweep did predict and one third are the named silence "the gene is not in the scorer's window".
> Re-fitting through that window on the same held-out pairs moves reliability 6 of 10 bins to 7 and
> AUPRC 0.559 to 0.677, but the prevalence gap this section names does not close (the +0.679 shift
> becomes +0.703), so **the verdict of failed below stands**. The predicted-target curve of the row
> at 245 pairs and 76.7%, whose weights band all 612,323 sweep targets, quotes 0.41 for the held-out
> pairs the gate excludes against a measured 0.0603 — ×6.85. See "The gate removed" (2026-09-22,
> later) and its pre-registration.

Three measurements of chromatin structure have now failed to say which element acts on which gene —
CTCF orientation, measured boundary strength, measured Hi-C contact at 5 kb — while the predicted
deletion works (LESSONS.md, "Three measurements of chromatin structure, three nulls"). The deletion's
output is therefore what downstream work will trust, and the number the sweep writes beside it is the
model's own effect size: `confidence = min(1, |log2 fold change|)`. That is a rescaled effect size,
not a probability, and nothing measured says what it means. `scripts/target_calibration.py`
(`attribution/target_calibration.py`, result `target_calibration`, 61 s, no model request, no
network) asks the CRISPRi screens to supply the missing scale, and then applies it to all 961,227
deleted elements.

Only features the sweep also has may enter, or the fit cannot be applied to it: the predicted drop in
K562 (the negated log2 fold change, floored at zero, non-zero only when the pair's gene is the
element's top predicted target), the top-target flag itself, the log distance from the element's
midpoint to the gene's **GENCODE v50** TSS — computed the same way on both sides, so the feature
transfers — whether the gene is the nearest TSS inside the element's CTCF node, and the registry
class. The benchmark's measured activity (DNase × H3K27ac in the screen's own cell) is deliberately
left out: the project has cached DNase for one chromosome, so a calibration using it could not be
applied to the sweep. What leaving it out costs is reported below.

**Coverage of both arms first, then any rate** (genomeos-8a's denominator lesson, and genomeos-79's
sharper form of it: a calibration bin is a rate conditioned on the pair having been scored at all).
Of 9,237 covered training pairs, 8,796 carry every feature: 437 of 451 regulated (96.9%) and 8,359 of
8,786 not regulated (95.1%). The 491 exclusions are one named category — 477 not-regulated and 14
regulated pairs whose gene has no GENCODE v50 gene entry — counted, never scored as a distance of
zero. The subset is not an easier problem: distance reads 0.4413 AUPRC on all covered pairs against
0.4436 on the feature-complete ones, activity over distance 0.5190 against 0.5265, the base rate
4.88% against 4.97%. On held-out K562 the arms are level — 112 of 114 regulated (98.3%) and 1,603 of
1,630 not (98.3%) — and the baselines again do not move (distance 0.3753 → 0.3733, activity over
distance 0.5502 → 0.5540, base rate 6.54% → 6.53%). The GENCODE gene-level TSS differs from the
benchmark's own `distanceToTSS` by a median 95 bp on those pairs, within 1 kb for 83% of them and
within 10 kb for 92% (on the training pairs a median 148 bp, 77% and 89%); the rest are genes whose
screen used another TSS.

Two populations, and the population matters more than the fit:

| population | training | rate | held out (K562) | rate |
|---|---|---|---|---|
| a pair a K562 screen tested, on a cCRE inside a node, every feature present | 8,796 | 4.97% | 1,715 | 6.53% |
| the same, and the tested gene **is** the element's top predicted target (the sweep's own pair) | 245 | 76.7% | 40 | 90.0% |

The pre-registration, fixed in `PREREGISTERED_CALIBRATION` before any held-out pair was scored:
fitted on the K562 training pairs with only the sweep's features, the calibrated probability is
reliable on held-out K562 and better calibrated than the effect size read as a confidence — (a) in at
least 7 of the 10 equal-count bins the bin's mean predicted probability lies inside the 95% Wilson
interval of the observed rate, and (b) the expected calibration error and the Brier score are both
lower than those of `min(1, |log2 fc|)` on the same pairs. The Hosmer–Lemeshow statistic is reported
everywhere but deliberately gates nothing: read on the training pairs before the held-out set was
touched, it rejects the fit's own in-sample curve (χ² 22.5, 8 df, p 0.004) at an expected calibration
error of 0.007, so at 8,796 pairs it is answering a question about counts of one and two positives per
bin (`CRITERION_NOTE`).

| held-out K562, 1,715 pairs (112 regulated) | ECE | MCE | Brier | bins inside | AUPRC |
|---|---|---|---|---|---|
| calibrated, the sweep's features | **0.0231** | 0.0843 | **0.04184** | 6/10 | 0.559 |
| the drop alone, fitted | 0.0256 | 0.0649 | 0.04639 | 7/10 | 0.428 |
| the effect size as a confidence (today's number) | 0.0575 | 0.1435 | 0.05522 | 0/10 | 0.428 |

| bin | pairs | coverage, regulated | coverage, not | predicted | observed | 95% interval |
|---|---|---|---|---|---|---|
| 1 | 171 | 0.974 | 0.983 | 0.0046 | 0.0292 | 0.0126–0.0666 **outside** |
| 2 | 172 | 0.974 | 0.983 | 0.0059 | 0.0058 | 0.0010–0.0322 |
| 3 | 171 | 0.974 | 0.983 | 0.0070 | 0.0175 | 0.0060–0.0503 |
| 4 | 172 | 0.974 | 0.983 | 0.0084 | 0.0058 | 0.0010–0.0322 |
| 5 | 171 | 0.974 | 0.983 | 0.0103 | 0.0175 | 0.0060–0.0503 |
| 6 | 172 | 0.974 | 0.983 | 0.0129 | 0.0058 | 0.0010–0.0322 |
| 7 | 171 | 0.974 | 0.983 | 0.0169 | 0.0234 | 0.0091–0.0586 |
| 8 | 172 | 0.974 | 0.983 | 0.0249 | 0.0640 | 0.0361–0.1109 **outside** |
| 9 | 171 | 0.974 | 0.983 | 0.0443 | 0.0936 | 0.0584–0.1466 **outside** |
| 10 | 172 | 0.983 | 0.983 | 0.3053 | 0.3895 | 0.3198–0.4641 **outside** |

> **Annotation, 2026-09-22 (sixth).** The `+0.679` below stands as computed and its meaning is
> narrower than the sentence around it. Split by screen, the held-out pairs need shifts of
> **-1.303** (K562_DC_TAP, 1,084 pairs), +0.945 (Xie, 416), +2.050 (Morris, 177) and +4.228
> (Klann and Reilly, 38): the one constant is their mean, and its sign is wrong for the screen
> that carries 63% of the population. "The held-out screens call 6.53% against 4.97%" is true
> and is not a fact about the two tables -- the spread inside each is larger than the gap
> between them, 3.5x across the fitted screens and 99x across the read ones. Not one screen
> here is in both tables. See "The prevalence term, registered before it was fitted"
> (2026-09-22, sixth).

The coverage columns come before the rate columns on purpose, and they are flat: no bin's rate is
inflated by having collected the better-covered elements. On training the leave-chromosome-out curve
is reliable by the same rule (9 of 10 bins, ECE 0.0061, Brier 0.02812, AUPRC 0.635).

**The pre-registered claim FAILS,** on clause (a) with 6 of 10 bins: the calibration is better than
today's number on every error it was compared on, and it is not reliable. The shape of the curve
transfers and its level does not, and the reason is visible in one line: the calibration's mean
predicted probability is 0.0441 while the held-out rate is 0.0653, a prevalence 1.31 times the
training screens'. All four failing bins fail in the same direction — the observed rate above the
predicted one — and none in the other. Add the single constant that matches the prevalence —
+0.679 in log odds, fitted on the held-out labels themselves, so a description of the failure and not
a test — and the same curve gives 9 of 10 bins inside their intervals, ECE 0.0100, Brier 0.04065,
Hosmer–Lemeshow p 0.08. A calibration is a rate conditioned on a population, and the base rate of a
CRISPRi screen is a property of how that screen chose its pairs. Nothing here was retuned after the
held-out set was read; the verdict stands as failed.

**The table to quote, with no model between it and the screens.** For the sweep's own population —
one element, its top predicted target — the drop bands need no feature but the drop, so no pair is
dropped for a missing value; the coverage columns read 1.00 on both arms everywhere but the zero band,
where one regulated pair of 26 lacks a GENCODE TSS (0.96). Each band's pairs are distinct elements.

| predicted K562 drop | pairs | regulated | rate | 95% interval |
|---|---|---|---|---|
| = 0 | 44 | 26 | 0.591 | 0.444–0.723 |
| 0 < drop ≤ 0.1 | 92 | 55 | 0.598 | 0.496–0.692 |
| 0.1 < drop ≤ 0.2 | 45 | 39 | 0.867 | 0.738–0.937 |
| 0.2 < drop ≤ 0.5 | 46 | 46 | 1.000 | 0.923–1.000 |
| 0.5 < drop | 59 | 59 | 1.000 | 0.939–1.000 |

Training and held-out pairs pooled, which is legitimate only because the held-out test above was
scored first. Read it as the conditional it is: *if* a K562 CRISPRi screen tests an element's top
predicted target, a drop above 0.2 was called regulated 105 times out of 105, and a drop at or below
0.1 about three times in five. The fitted version of the same population (`deletion_drop`,
`log_tss_distance`, `node_target`, 245 training pairs) is reliable in 3 of 3 bins in sample (ECE
0.014) and in 3 of 3 on the 40 held-out pairs (ECE 0.066), which decides nothing at that size and is
excluded from the pre-registration for exactly that reason.

**The sweep, read through the calibration.** 961,227 deleted elements, 612,323 with a predicted
target; 18,558 of those name a gene GENCODE v50 has no entry for and are refused rather than banded,
leaving 593,765 (440,377 for the coding target). The registry classes present are only dELS and pELS,
because the sweep scored those two classes and no others — three of the five fitted class levels are
empty in both the benchmark and the sweep, and the class term is one contrast (pELS −0.48 against
dELS).

| confidence band | any gene | coding gene | any gene, dELS | any gene, pELS |
|---|---|---|---|---|
| 0.05–0.1 | 45 | 0 | 34 | 11 |
| 0.1–0.25 | 3,837 | 1 | 2,591 | 1,246 |
| 0.25–0.5 | 178,387 | 184,256 | 166,181 | 12,206 |
| 0.5–0.75 | 293,568 | 175,655 | 217,661 | 75,907 |
| 0.75–0.9 | 61,971 | 42,552 | 39,538 | 22,433 |
| 0.9–1 | 55,957 | 37,913 | 26,869 | 29,088 |

Every predicted target lands at 0.25 or above, and that is the result, not a bug: the population the
calibration was measured on — an element's top predicted target, tested by a screen — is regulated
77% of the time, so a model fitted on it hands every sweep target a high probability. What separates
the sweep's targets from each other is the drop, and on that axis the genome is thin: 46.7% of coding
targets have a K562 drop of exactly zero, 39.3% sit at or below 0.1, and only 24,114 (5.5%) clear
0.2, the band where the screens called 46 of 46 and where no tested pair was ever called unregulated
(the not-regulated arm of the top two bands is empty, 0 pairs, which is the finding and not a gap in
coverage). The 0.9–1 band holds 37,913 coding targets. Per
chromosome the share clearing 0.2 runs from 3.2% (chr18) to 10.1% (chr19), — chr19 highest, chr18 lowest,
chrX at 7.9% — and chrY contributes 125 coding targets of which 12 clear 0.2. The full per-chromosome and per-class tables are in `data/results/target_calibration.json`;
the committed sweep summaries were not touched.

**What the missing activity would have bought.** With the benchmark's measured DNase × H3K27ac added,
the leave-chromosome-out ECE improves from 0.0061 to 0.0041 and AUPRC from 0.635 to 0.667; on held-out
K562, Brier 0.04184 → 0.03772 and AUPRC 0.559 → 0.640, with the same 6 of 10 bins. The activity is
worth having and cannot be applied: it exists for the benchmark's pairs and for one chromosome of the
genome.

**Scope, stated as narrowly as it holds.** The calibration is measured in one cell line on one kind of
element: K562 CRISPRi pairs whose element is an ENCODE cCRE inside a CTCF node and whose gene a screen
chose to test, which means a gene expressed in K562 within about a megabase. The probability is
conditional on that test happening. It says how often a K562 screen calls such a pair regulated; it
does not say how often an element regulates a gene, and for a sweep element whose strongest predicted
tissue is not K562 the drop entering the calibration is still the K562 drop. The ENCODE-training
objection is narrowed, not closed: E1's 2,260 fresh pairs (commit 2a07c80) put units at 637 of 1,096
agreeing (0.581) against matched nulls 490 of 968 (0.506), +0.075 with an upper 95% bound of 0.111 at
one-sided p 0.00037 — the weak band, since the pre-registration asked for 0.10 in size, and a small p
does not convert a weak effect into the declared one. One cell line carries it (GM12878 +0.089 on
1,895 pairs; Jurkat −0.001 on 365). So in lymphoblastoid cells the model weakly and detectably tracks
a reporter assay's direction, and outside them that endpoint says nothing. E1's declared secondary did
pass, and it is the first met prediction of that series: agreement is higher where DAP-G fine-maps the
variant, +0.198 on 132 units against 90 controls (p 0.0018) against +0.057 for the rest — a reason to
expect the transfer where the causal variant is known, not evidence that it happens.

**What would falsify the transfer** (`FALSIFIES_TRANSFER`). A CRISPRi screen in another cell line with
an AlphaGenome line in the deletion table (GM12878, HepG2, IMR-90 at this coverage) whose observed
rate per band falls outside the band in most populated bands, or whose ECE exceeds the effect size
read as a confidence. A screen that tiles elements without regard to K562 activity finding the top
band's rate far below the band, which is what the selection on testability would look like from
outside. The bands holding no better on a new screen's fine-mapped subset than off it, which would
remove the one reason E1 gives to expect any transfer. And the prevalence failure above, generalised:
a screen whose base rate differs from these screens' 5% needs its own intercept, so a band quoted
without its population's rate beside it measures nothing.

**Limits.** Held-out GM12878 has 62 covered pairs and 14 regulated and is refused as too small to
carry a reliability table; HCT116, Jurkat and WTC11 have no AlphaGenome line in the deletion table and
are refused with that reason. The sweep keeps one gene per element, so the calibration cannot rank a
second target. The GENCODE gene-level TSS stands in for the screen's own TSS. The 0.2-and-above band
rests on 105 pairs and its interval's lower bound is 0.92, not 1.0.

**Next, proposed as roadmap rows** (for the roadmap editor to fold; not written into ROADMAP.md):

1. Area I: the prevalence term. Re-fit the intercept per screen with the screen's own base rate as an
   offset, and pre-register that a band holds across the three K562 datasets (Gasperini, Schraivogel,
   Nasser) whose rates differ, before it is quoted for anything.
2. Area I: the drop above 0.2 as a shortlist. 24,114 coding targets clear it genome-wide; hold that
   list against an independent perturbation in a second cell line rather than against a model.
3. Area I: DNase per element genome-wide (one bigWig per cell line, the reader exists), which is the
   one cheap feature the calibration had to leave out and which improved every error it was allowed
   to touch.

## The library's bridge to the screens: seven blocks of 882 (2026-09-17)

**The library's positive set whose expectation comes from a measurement rather than from a model is
seven blocks.** Of the 882 blocks of the real unknown, 163 hold an element the sweep scored at all
and 7 hold one in the band where the CRISPRi screens called 105 of 105 — 2,283 oligos of 284,598.
Anyone who reads 882 should read 7 in the same breath, and the reason they are not in conflict is
that a thin positive set is what 0.45% coverage predicts: the screens tested annotated elements near
expressed genes, and a real-unknown block is unannotated by construction, so most of its sequence was
never eligible for the measurement that would have given it a prior. The library is not weakened by
this; it is the instrument that would change it.

The calibration leaves 24,114 coding targets in the band where the CRISPRi screens called 105 of
105. genomeos-9c's MPRA library tiles the 882 blocks of the real unknown. If the two overlap, those
oligos carry a prior from a measurement that is not the model's own effect size, and the library
tests directly whether the band holds off the population the screens chose.
`scripts/shortlist_in_real_unknown.py` (`shortlist_in_real_unknown`, committed results and the local
element tables only, no network) intersects them. It reproduces the 24,114 exactly, which is the
cross-check that it is the same shortlist.

**Seven blocks of 882, eight elements.** Three relaxed, four tolerant; none of the 69 syntax blocks
and none of the 29 recent ones. The strongest is chr7:8,752,961-8,938,786 (tolerant), whose element
moves NXPH1 by 1.15; then chr10:73,783,652-73,785,548 → ZSWIM8 (0.67), chr1:153,205,120-153,217,583
→ PRR9 (0.63), chr16:8,180,769-8,252,085 → RBFOX1 (0.53), chr4 → RASGEF1B (0.41), chr6 → TFAP2B
(0.40, two elements) and chr1:69,563,715-69,567,490 → LRRC7 (0.22).

**The context that makes seven readable.** 163 of the 882 blocks hold a scored element with a coding
target at all (relaxed 81, tolerant 62, syntax 11, recent 7, unmeasured 2); of those, 7 hold one
that clears the band. The best element in the next block down reads 0.17. So the shortlist and the
real unknown barely overlap, and the reason is structural rather than biological: the sweep scored
registry elements lying inside a CTCF node, and a real-unknown block is unannotated by construction,
so most of its sequence was never eligible for a deletion. 520,850 elements genome-wide have no
coding target at all and 416,263 sit below the band.

**What this means for the library.** It carries an internal bridge, but a thin one: eight oligo
positions with an independent prior, not a positive set. A library that wants the drop band tested
off the screens' population has to tile those seven blocks and say so, and the honest expectation
for the other 875 is that nothing in the CRISPRi evidence speaks to them. (This intersection counts
a block as scored only when its element has a coding target with a GENCODE TSS, which is stricter
than the 531 tested blocks and 331 leads counted for the tier reading; the two numbers answer
different questions and are not in conflict.)

## The proxy replaced by a measurement: mappability read at four read lengths, and the 8,131 oligos the repeat flag passes that a 24-mer cannot place (2026-09-17, later)

The section below charged 17,112 of the untouched blocks' windows to `interspersed_repeat` and 210 to
`segmental_duplication`, and called the remainder **attributable to one locus**. Both rules were a
**proxy**. They ask how much of a window is annotated repeat, not whether the window's sequence can
be told apart from the rest of the genome, and the difference showed in two places: the attributable
count swung **19,084 to 28,563 oligos** as `INTERSPERSED_MAX` moved from 0.25 to 0.75, a 50% swing on
a number nothing measures; and `non_unique_in_block`, the one rule that did ask about uniqueness,
compared each window only against its own block and so fired on **23 windows** — a lower bound with
no relation to the genome.

`genomeos/attribution/mappability.py` and `scripts/mappability.py` read the quantity itself.

### The track, and what reading it cost

**Umap multi-read mappability**, GRCh38 (UCSC hg38): for each base, the share of the length-*k* reads
overlapping it that align to exactly one place in the assembly. Karimzadeh, Ernst, Kundaje and
Hoffman 2018, *Nucleic Acids Research* 46:e120, doi:10.1093/nar/gky677; read from UCSC's copies at
`https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability/k{k}.Umap.MultiTrackMappability.bw`.

**k is part of every number below and the four columns are never pooled.** A 24-mer and a 100-mer
answer different questions about a 300 bp oligo, and the answers here differ by 31,170 oligos.

| track | *k* | file bytes | bytes read | range requests | wall time |
|---|---|---|---|---|---|
| Umap k24 | 24 | 2,028,863,031 | 21,038,867 | 1,903 | 331 s |
| Umap k36 | 36 | 1,647,208,796 | 18,134,671 | 1,747 | 306 s |
| Umap k50 | 50 | 1,331,314,604 | 16,234,544 | 1,687 | 295 s |
| Umap k100 | 100 | 864,604,710 | 15,971,162 | 1,567 | 268 s |
| **total** | — | **5,871,991,141** | **71,379,244** | **6,904** | **1,199 s** |

Nothing was downloaded: the four files weigh 5.87 GB and **1.22% of their bytes were fetched**,
through `genomeos/attribution/bigwig.py`, the same range reader Zoonomia and Gnocchi use — header,
chromosome tree and R-tree first, then only the data sections overlapping the oligos. `BYTE_CAP` is
96 MB per track and was not reached. The kept cache is the per-oligo summaries only,
**11.3 MB in `data/knowledge/mappability/`** against a `CACHE_CAP` of 100 MB; no track bytes are
stored. The whole set was read twice, once cold and once with `--refresh`, and every count in this
section is identical between the two runs.

### The rule, with its one justification

    a base is uniquely readable when its multi-read mappability is at least BASE_CUT = 1.0, and an
    oligo is mappable when at least MAPPABLE_FRACTION of its 300 bases are uniquely readable.

`BASE_CUT` is 1.0 rather than a fraction because anything under 1.0 already means some read covering
that base has a second home in the genome, which is exactly the failure being avoided; the tolerance
the design really has is over *how much of the oligo* may fail, and that is `MAPPABLE_FRACTION`,
reported at **0.50, 0.90 and 0.99** and not tuned. The universe is the library's **test arm as
written: 91,919 oligos of 300 bp over 878 blocks**, which is already the set that passed the
synthesis rules.

### The two arms, reported separately

The repeat proxy does not depend on *k*, so it is one number in every row. Nothing here adds the two
arms together.

| *k* | oligos | attributable by the **repeat proxy** | attributable by the **mappability track** | both arms | arms agree |
|---|---|---|---|---|---|
| 24 | 91,919 | 51,789 | 57,346 | 43,658 | 76.3% |
| 36 | 91,919 | 51,789 | 70,276 | 48,112 | 71.9% |
| 50 | 91,919 | 51,789 | 78,276 | 49,960 | 67.2% |
| 100 | 91,919 | 51,789 | 88,516 | 51,603 | 59.6% |

And over the **680 untouched blocks alone** — the published proxy's own universe, which contributes
41,037 of the test arm's oligos across 676 blocks:

| *k* | oligos | repeat proxy | mappability track | both arms | proxy passes, track says no |
|---|---|---|---|---|---|
| 24 | 41,037 | 23,863 | 26,081 | 20,028 | **3,835** |
| 36 | 41,037 | 23,863 | 31,806 | 22,120 | **1,743** |
| 50 | 41,037 | 23,863 | 35,145 | 22,980 | **883** |
| 100 | 41,037 | 23,863 | 39,541 | 23,776 | **87** |

### The cross-tabulation, which is the result

Over the whole test arm at `MAPPABLE_FRACTION` 0.90. The disagreement is not a problem to reconcile;
it is what the measurement had to say.

| *k* | both flag | proxy flags, **track says mappable** | **proxy passes, track says unmappable** | neither flags |
|---|---|---|---|---|
| 24 | 26,442 | 13,688 | **8,131** | 43,658 |
| 36 | 17,966 | 22,164 | **3,677** | 48,112 |
| 50 | 11,814 | 28,316 | **1,829** | 49,960 |
| 100 | 3,217 | 36,913 | **186** | 51,603 |

Both directions are large, and they move opposite ways with *k*.

- **The dangerous direction.** 8,131 oligos at k24, 3,677 at k36, 1,829 at k50 and 186 at k100 pass
  the repeat proxy and fail the measurement. As a share of the 51,789 the proxy calls attributable
  that is **15.7%, 7.1%, 3.5% and 0.36%**. These are oligos the design would order believing a
  positive could be placed. They are not an annotation artefact: the within-block rule found 23 of
  them and the genome-wide track at k24 calls 34,573 of the 91,919 test oligos unmappable, so the old
  lower bound was low by three orders of magnitude.
- **The other direction, which is larger.** 13,688 oligos at k24 rising to 36,913 at k100 are flagged
  by the proxy and fully mappable. At k100 that is 36,913 of the 40,130 proxy-flagged oligos: **92% of
  the repeat flag is not an attribution problem at all** once the read is 100 bp (70.6% at k50, 34.1%
  at k24). This is the same conclusion the lentiMPRA calibration reached from the other side, now
  measured: interspersed repeat is a caveat about annotation, not about whether a locus can be found.

### Does the flag carry information once the covariates are held fixed?

`genomeos/compare.py`, direct standardisation on block length, GC and distance to a coding TSS.
Targets are the 40,130 proxy-flagged oligos, controls the 51,789 it passes, outcome "the track calls
this unmappable". **40,128 targets matched, 2 dropped for want of a control, 0 excluded for a missing
covariate.** Imbalance before matching, identical for every *k*: block length 101,366 against 87,393
(ratio 1.16), GC 0.3867 against 0.3667 (1.055), TSS distance 296,976 against 313,568 (0.947).

| *k* | flagged unmappable | passed unmappable | raw difference | **matched difference** | 95% upper |
|---|---|---|---|---|---|
| 24 | 0.6589 | 0.1581 | +0.5019 | **+0.5009** | 0.5067 |
| 36 | 0.4477 | 0.0769 | +0.3767 | **+0.3708** | 0.3763 |
| 50 | 0.2944 | 0.0401 | +0.2591 | **+0.2543** | 0.2591 |
| 100 | 0.0802 | 0.0039 | +0.0766 | **+0.0763** | 0.0790 |

So the proxy is informative and never sufficient: at k24 a flagged oligo is 50 percentage points more
likely to be unmappable at equal block length, GC and TSS distance, and by k100 the gap is 7.6
points. The matched difference barely moves off the raw one, which is the honest reading of an
imbalance of 1.16 and 1.055 — small, and reported rather than assumed away.

### Sensitivity: the rule at three values

| *k* | mappable at 0.50 | at **0.90** | at 0.99 | dangerous cell at 0.50 / **0.90** / 0.99 |
|---|---|---|---|---|
| 24 | 75,782 | **57,346** | 45,969 | 439 / **8,131** / 15,352 |
| 36 | 82,496 | **70,276** | 68,559 | 142 / **3,677** / 4,523 |
| 50 | 86,163 | **78,276** | 77,090 | 71 / **1,829** / 2,315 |
| 100 | 89,884 | **88,516** | 88,302 | 25 / **186** / 251 |

k24 is the column that is sensitive to the cut-off (75,782 to 45,969, a factor of 1.65); at k50 and
above the rule barely matters (86,163 to 77,090; 89,884 to 88,302). The proxy's own swing is reported
in the same result and the track does not follow it: as `INTERSPERSED_MAX` moves 0.25 → 0.50 → 0.75
the proxy's count goes 41,013 → 51,789 → 62,177 while every track column stays exactly where it is,
and the dangerous cell grows 4,421 → 8,131 → 12,943 at k24.

### Every group, with its covariates

Oligo length is 300 bp by construction, so the length distribution reported is that of the blocks the
oligos were tiled from (quartiles, oligo-weighted). No tier is used as a control anywhere in this
section: the groups are cut by the two rules themselves.

| group | *k* | oligos | block length q25 / med / q75 | GC | median TSS distance | no TSS |
|---|---|---|---|---|---|---|
| both flag | 24 | 26,442 | 50,739 / 102,355 / 181,805 | 0.393 | 290.5 kb | 0 |
| proxy flags, track mappable | 24 | 13,688 | 48,767 / 97,001 / 178,533 | 0.377 | 309.6 kb | 0 |
| proxy passes, track unmappable | 24 | 8,131 | 40,858 / 83,723 / 158,892 | 0.370 | 296.0 kb | 0 |
| neither flags | 24 | 43,658 | 43,799 / 89,270 / 164,906 | 0.367 | 316.4 kb | 0 |
| both flag | 100 | 3,217 | 50,375 / 115,335 / 178,533 | 0.400 | 275.1 kb | 0 |
| proxy flags, track mappable | 100 | 36,913 | 49,944 / 101,277 / 181,805 | 0.383 | 300.5 kb | 0 |
| proxy passes, track unmappable | 100 | 186 | 31,699 / 83,650 / 164,641 | 0.393 | 261.9 kb | 0 |
| neither flags | 100 | 51,603 | 43,673 / 87,393 / 164,641 | 0.367 | 313.8 kb | 0 |

The dangerous cell comes from the **shortest** blocks at every *k* (median 83.7 kb against 89.3 kb at
k24, 83.7 against 87.4 at k100) and, at k24, is the group second-closest to a gene. Dropping it would
repeat the bias the thin-block decision was taken to avoid, so it is reported as a flag and not as a
deletion. All eight groups, and the four cut over the untouched blocks alone, are in the result with
these three covariates.

### Coverage: what could not be assessed, by name

- **0 oligos unassessed at every *k*.** All 91,919 have a mappability reading and a proxy reading; no
  chromosome of the arm is absent from a track and `BYTE_CAP` was never reached.
- **Oligos with no track value anywhere: 691 (k24), 405 (k36), 251 (k50), 83 (k100).** These are
  measured as unmappable, not as unassessed: Umap writes nothing where mappability is zero. That
  reading is licensed by a check in the output rather than by assumption — the test arm contains
  **not one N base** and no oligo is other than 300 bp, so a missing value cannot be an assembly gap.
- **4 of the 680 untouched blocks have no test oligo to assess**, all on chrY
  (`chrY:12920478-12930164`, `chrY:18650279-18655519`, `chrY:56954169-57015104`,
  `chrY:57062405-57067746`): the library's test arm stops at chrX. They are named, not zeroed, and
  they are why this lane's untouched-block universe is 676 blocks and 41,037 oligos where the section
  below reports 680 blocks and 45,570 windows — the rest of the difference being the 4,264 windows
  the synthesis rules had already removed before the manifest was written.
- **50,882 of the 91,919 test oligos come from blocks an assay has already touched** and are outside
  the published proxy's universe; they are counted separately, never pooled with the 41,037.
- **No proxy layer was missing**: RepeatMasker and genomicSuperDups are cached for all 23
  chromosomes of the arm, so no oligo is silently repeat-free.

### What changes for the library's quotable numbers, and what does not

- **The orderable count does not move. 41,297 of the untouched blocks' 45,570 windows can be
  synthesised and read back, and the test arm is 91,919 oligos.** Mappability is not a synthesis
  rule; it changes nothing a vendor quotes, and no number in the section below is withdrawn.
- **"Attributable to one locus" can no longer be quoted without a *k*.** Over the untouched blocks it
  is 26,081 oligos at k24, 31,806 at k36, 35,145 at k50 and 39,541 at k100, against the proxy's
  23,863. The proxy's 19,084–28,563 swing on an arbitrary repeat cut-off is gone; what replaces it is
  a spread that is a property of the follow-up assay's read length, which the design can choose and
  state.
- **Which *k* to quote is a question about the follow-up, not about the reporter.** An MPRA reads its
  members back by barcode, so mappability does not decide whether a number comes out; it decides
  whether a positive can be pinned to this locus afterwards. A CRISPRi follow-up places ~20 bp
  guides, so k24 is its column; a 100 bp resequencing or allele-specific read is k100's. The library
  should carry all four values per oligo and quote the one belonging to the experiment it promises.
- **The repeat flag stays a flag.** 92% of the proxy's flagged oligos are fully mappable at k100 and
  70.6% at k50, so dropping them for want of attribution would throw away sequence that can in fact
  be placed. The 8,131 oligos at k24 that the proxy passes and the track cannot place are the other
  half of the same correction, and they are the reason this arm exists.

Result: `data/results/mappability_real_unknown.json`. Tests: `tests/test_mappability.py`, on a bigWig
built in memory, no network.

## What the library cannot measure: 45,570 oligos really tile the 13.77 Mb, 41,297 can be ordered, and the repeat flag is not a reason to drop them (2026-09-17)

The section below established that 680 of the 882 real-unknown blocks — 13.77 Mb — have never been
touched by lentiMPRA, VISTA or the CRISPRi benchmark, and that at 300 bp this is **45,900 oligos**.
45,900 is arithmetic: 13,770,000 divided by 300. `scripts/measurability.py` asks the opposite of every
other question this project has asked about the unknown space. Not *what is in it* — no model, no
extrapolation — but **which parts of it a reporter assay could not measure even in principle**, so the
design excludes or flags them before oligos are paid for. It sharpens an artefact of the instrument.

`genomeos/attribution/measurability.py` fixes what "cannot be measured" means, in code, as named
constants with their reasons, in two families that are reported separately and never added together:

- **Family A, not synthesisable or not resolvable.** `assembly_gap` (more than `N_MAX` = 10% of the
  window is N: there is no sequence to order), `gc_extreme` (GC outside `GC_LOW` 0.25 – `GC_HIGH` 0.75:
  array synthesis and pool PCR drop out at the extremes), `homopolymer` (a run of one base at least
  `HOMOPOLYMER` = 10 long: coupling slips and the read-back miscounts the run), `tandem_low_complexity`
  (at least `TANDEM_MAX` = 50% RepeatMasker Simple_repeat, Low_complexity or Satellite),
  `non_unique_in_block` (at least `DUP_KMER_MAX` = 50% of the window's 40-mers recur inside its own
  block, so two windows of the tiling are the same molecule).
- **Family B, synthesisable but not attributable to one locus.** `segmental_duplication` (at least
  `SEGDUP_MAX` = 50% of the window inside a curated genomicSuperDups pair) and `interspersed_repeat`
  (at least `INTERSPERSED_MAX` = 50% LINE, SINE, LTR, DNA, Retroposon or RC). The assay returns a
  number; the number belongs to a sequence rather than to a place.

Two more categories are arithmetic, named so the base accounting closes to the last base:
`shorter_than_one_oligo` and `tiling_remainder`. **The set of 680 blocks and 13.77 Mb is reproduced
independently by this lane before anything is measured on it**, which is the check that the two
readings are about the same sequence.

### Excluded by reason, pooled over 680 blocks and 13,772,419 bp

Windows are 300 bp at step 300 — **no overlap**, the tiling whose arithmetic gave 45,900. Each window
is charged to exactly one reason in the precedence order above; the raw column is how many windows
fail that rule at all, so the two differ where a window fails several.

| reason | family | windows charged | raw | Mb | share of tiled |
|---|---|---|---|---|---|
| `assembly_gap` | A | **0** | 0 | 0.000 | 0.00% |
| `gc_extreme` | A | 751 | 751 | 0.225 | 1.65% |
| `homopolymer` | A | 3,394 | 3,438 | 1.018 | 7.45% |
| `tandem_low_complexity` | A | 105 | 162 | 0.032 | 0.23% |
| `non_unique_in_block` | A | 23 | 35 | 0.007 | 0.05% |
| `segmental_duplication` | B | 210 | 245 | 0.063 | 0.46% |
| `interspersed_repeat` | B | 17,112 | 19,111 | 5.134 | 37.55% |
| `shorter_than_one_oligo` | arithmetic | — | — | 0.000 | 0.00% |
| `tiling_remainder` | arithmetic | — | — | 0.101 | 0.74% |
| **usable, both families** | — | **23,975** | — | **7.192** | 52.61% |

**Base accounting: 13,772,419 of 13,772,419 bp, 0 unaccounted.** No silent zeros: `assembly_gap` and
`shorter_than_one_oligo` are genuinely zero and are printed as zero. There is **not one N base in the
whole 13.77 Mb**, because the budget's blocks were cut from called sequence; and no block is shorter
than one oligo, so nothing is lost to being too small to tile.

### The number the design needs, beside 45,900

- **45,570 oligos actually tile**, not 45,900. Exactly 13,772,419 / 300 = 45,908 windows would fit if
  blocks divided evenly, and 338 are lost to block tails shorter than 300 bp (101,419 bp of remainder).
  The arithmetic over-counted by 0.7%.
- **41,297 oligos can be ordered and read back** — everything Family A does not exclude. That is
  **90.6% of the tiling**, a 9.4% haircut, and it is the number a synthesis quote should be based on.
- **23,975 oligos are also attributable to one locus** — 52.6%. The gap between 41,297 and 23,975 is
  almost entirely `interspersed_repeat`.
- **53 of 680 blocks fall below three attributable oligos**, so they cannot carry a tiling as opposed
  to a sample. They are not a random 53: median length 1,402 bp against 7,996 for the rest, GC 0.440
  against 0.386, and **median distance to a coding TSS 23.6 kb against 100.9 kb**. Dropping them
  therefore drops the blocks closest to genes, which is a choice the design has to make knowingly.

### The repeat flag calibrated against sequence that was already measured

Family B is where the loss looks catastrophic, so it is the one that had to be checked against reality
rather than against another tier. The same rules were run over 18,991 300 bp windows centred on
**lentiMPRA elements — sequence an episomal reporter has already ordered and got numbers out of**
*(R6 note, 2026-09-28: the reporter is integrated by lentivirus, not episomal; Agarwal et al. 2025, Nature, https://www.nature.com/articles/s41586-024-08430-9.)*
(both columns are charged, mutually-exclusive counts, so Family B is counted after Family A):

| | untouched real unknown | lentiMPRA, already measured |
|---|---|---|
| windows | 45,570 | 18,991 |
| Family A excludes | 9.38% | **4.67%** (homopolymer 825, GC 30, tandem 31, N 1) |
| Family B flags | 38.01% | **31.31%** (interspersed 5,336, segdup 610) |
| median GC | 0.386 (blocks) | 0.497 |
| median distance to a coding TSS | 100.9 kb | 45.9 kb |

**A third of the elements a real MPRA library measured successfully are themselves more than half
interspersed repeat.** So `interspersed_repeat` is not a measurability exclusion at all: it is an
attribution caveat, and a library that dropped it would be dropping sequence of exactly the kind
lentiMPRA already reads. The recommendation is to **carry the 17,112 repeat-derived oligos flagged,
not excluded**, and to exclude only Family A. Family A's 4.67% false-exclusion rate on already-measured
sequence also says the Family A cut-offs are in roughly the right place rather than punitive.

This is not a matched comparison and is not offered as one: the lentiMPRA elements are ENCODE cCREs
and differ from the untouched blocks in GC (0.497 against 0.386) and in distance to a coding TSS
(45.9 kb against 100.9 kb), which is why both covariates are in the table. It is a plausibility check
on the cut-offs, not a test of a hypothesis, and no tier average is used as a control anywhere here.

### Per case, with the covariates that decide every comparison in this project

| case | blocks | Mb | oligos tiled | orderable (A) | attributable (A+B) | blocks < 3 | median length | GC | median TSS |
|---|---|---|---|---|---|---|---|---|---|
| tolerant | 237 | 7.60 | 25,211 | 22,850 | 13,396 | 4 | 16,104 | 0.374 | 176.1 kb |
| relaxed | 338 | 5.70 | 18,837 | 17,078 | 9,769 | 27 | 6,322 | 0.388 | 97.6 kb |
| syntax | 61 | 0.22 | 687 | 612 | 403 | 11 | 2,679 | 0.440 | 18.1 kb |
| recent | 26 | 0.15 | 496 | 447 | 271 | 3 | 2,890 | 0.448 | 14.9 kb |
| unmeasured | 18 | 0.10 | 339 | 310 | 136 | 8 | 1,701 | 0.407 | 41.8 kb |

The cases are listed, not compared: they differ by an order of magnitude in length and by more than
ten-fold in distance to a coding TSS, so no ratio between two of these rows means anything. What the
table is for is budgeting — **the 61 untouched blocks of the 69 the project cares most about are 687
oligos, about 1.5% of the library, and 403 of them are attributable** — and for seeing that 11 of those
61 and 8 of the 18 unmeasured blocks cannot carry a tiling at all.

Each excluded group also carries its own covariates, because a reason that selects on GC or on
distance to a gene biases what the library measures:

| charged reason | windows | blocks | median GC | median TSS |
|---|---|---|---|---|
| usable | 23,975 | 673 | 0.363 | 278.0 kb |
| `gc_extreme` | 751 | 212 | 0.237 | 511.9 kb |
| `homopolymer` | 3,394 | 516 | 0.377 | 248.0 kb |
| `tandem_low_complexity` | 105 | 66 | 0.377 | 310.3 kb |
| `non_unique_in_block` | 23 | 10 | 0.393 | 282.1 kb |
| `segmental_duplication` | 210 | 52 | 0.370 | 91.7 kb |
| `interspersed_repeat` | 17,112 | 601 | 0.387 | 273.5 kb |

(Window length is 300 bp by construction, so the length reported per group is the length of the blocks
the windows came from; those medians are window-weighted and therefore larger than the per-block ones.)
`gc_extreme` is the one reason that selects hard: the windows it removes sit at GC 0.237 and a median
half-megabase from any coding gene. The exclusion is AT-rich, gene-poor sequence, and saying so is the
point of printing the column.

### Composition of the 13.77 Mb, per block

Repeat-derived fraction per block: quartiles **0.265 / 0.392 / 0.473** — the untouched real unknown is
about 39% repeat by base, of which almost all is interspersed (median 0.373) and very little tandem
(median 0.012). Segmental duplication is negligible: **63 of 680 blocks touch a curated duplication at
all, 0.078 Mb in total, and only 6 blocks touch one at ≥ 0.98 identity** — unsurprising, since the
organiser already set aside any block half-duplicated as a copy before the real unknown was defined.

### Every arbitrary threshold at three values, so no number rests on a choice

| threshold moved | orderable (A) | attributable (A+B) | blocks < 3 usable |
|---|---|---|---|
| `n_max` 0.0 / 0.10 / 0.50 | 41,297 / 41,297 / 41,297 | 23,975 / 23,975 / 23,975 | 53 / 53 / 53 |
| GC 0.20–0.80 / 0.25–0.75 / 0.30–0.70 | 41,901 / 41,297 / 36,838 | 24,369 / 23,975 / 21,044 | 53 / 53 / 67 |
| `homopolymer` 8 / 10 / 12 | 37,807 / 41,297 / 42,622 | 21,909 / 23,975 / 24,756 | 63 / 53 / 48 |
| `tandem_max` 0.25 / 0.50 / 0.75 | 40,897 / 41,297 / 41,367 | 23,650 / 23,975 / 24,044 | 55 / 53 / 53 |
| `dup_kmer_max` 0.25 / 0.50 / 0.90 | 41,218 / 41,297 / 41,319 | 23,963 / 23,975 / 23,984 | 53 / 53 / 53 |
| `segdup_max` 0.01 / 0.50 / 1.00 | 41,297 / 41,297 / 41,297 | 23,961 / 23,975 / 23,999 | 54 / 53 / 51 |
| `interspersed_max` 0.25 / 0.50 / 0.75 | 41,297 / 41,297 / 41,297 | 19,084 / 23,975 / 28,563 | 94 / 53 / 32 |

`n_max` does nothing at any value, because there are no N bases. The orderable count moves by at most
13% across every Family A threshold, so **"about 41,000 orderable oligos" survives the choice of
cut-off**; the attributable count swings from 19,084 to 28,563 on the interspersed threshold alone,
which is one more reason to carry that flag rather than act on it.

### What this lane could not assess, and what it is not

- **Genome-wide mappability was not measured.** `non_unique_in_block` compares a window only against
  its own block, which is a lower bound on non-uniqueness and is why it fires on just 23 windows. A
  real mappability read needs a 36-mer or 50-mer uniqueness track (Umap/Bismap) or a whole-genome
  k-mer index, neither of which this project holds, and the disk could not take one. `segmental
  duplication` and `interspersed_repeat` are the proxies standing in for it, and they are proxies.
- **Coverage is complete otherwise: 680 of 680 blocks assessed, 0 unassessed, 0 without a distance to
  a coding TSS.** Every layer needed — the organiser's blocks, RepeatMasker, genomicSuperDups, the
  bgzip-indexed reference, GENCODE — was cached for all 24 chromosomes.
- **The thresholds were fixed in the module before any block was read and were not moved afterwards.**
  Where a threshold is arbitrary it is reported at three values rather than chosen.
- **This says nothing about whether any of the 13.77 Mb is functional**, and nothing about whether an
  MPRA is the right assay for it. It says how much of it an MPRA could physically ask about: about
  41,000 oligos of the 45,900 the arithmetic promised, of which about 24,000 would also point at one
  place in the genome.

`data/results/measurability_real_unknown.json`; `uv run python scripts/measurability.py`. No network.


### Beside the census, not inside it: the elements whose bases were measured (2026-09-17, later)

The cross-assay census asks reciprocal overlap — is the tested interval *this* element — and that is
the right question for an assay that perturbs an element and the wrong one for an assay that perturbs
bases. The layer had already accepted that satmut speaks about bases: a base-level null is evidence
about *those bases*, it keeps `experimental` as its kind, and it stays out of `assays_disagreeing`.
Then interval geometry decided which elements were allowed to hear it, which is where the seam showed
(genomeos-0e's ruling, 2026-09-17).

So the count now exists under its own name, `base_level_beside_the_census`, and is **never added to
the raised total**: 15 compiled elements have at least one of their
own bases substituted, of which 4 are the ones the element rule
raises. The other 11 are readings the census discards on interval shape rather than on what
was perturbed — including 1 whose measured bases are all inert, the first genuinely inert
elements the layer has, and the reading a base-level assay is best placed to make.

**Every row carries the share of the element actually substituted, and that is the whole guard.** The
measured fraction runs from 0.063 to
1.000, and the functional share of the measured bases from
0.0000 to
0.8714 across the same twenty-one experiments.
A count with no fraction beside it would let a thin measurement upgrade a fact, which is the trap the
reciprocal rule was built to stop: 15 bases of a 237 bp element and 350 of a 350 bp one are not the
same evidence and must not print as the same row. No threshold is imposed; the fraction is printed
instead, because a cut chosen after looking at twenty-one experiments is a cut chosen to fit them.

## The deletion feature read from the sweep's own cache, not from the compact table (2026-09-22)

The CRISPRi benchmark above scored `deletion_drop` as zero whenever the screen's measured gene was
not the element's single top predicted target. That zero said "the compact table had nothing to say
about this pair", and the logistic model consumed it as "the model predicts no effect". Those are
different statements and only one of them is a measurement.

The limit belongs to the derived table, not to the run. `data/knowledge/alphagenome/all_elements/<chrom>.json`
is a compact projection that keeps one gene per element (`compact()`, `scripts/enhancer_targets_all.py`),
but the same sweep also wrote a per-element response cache at
`data/knowledge/alphagenome/elements/<chrom>.json.gz` with `threshold=0.0` (`worker_scorer`), and that
cache carries **every gene in the scorer's 1 Mb window with a signed log2 fold change on each of
K562, HepG2, GM12878 and IMR-90's own track, uncensored**. `crispri_direction_both.py` already reads
it; this section moves `crispri.py` onto it. **Zero AlphaGenome requests**: the answer was on disk.

### The registration, written and committed before the benchmark was re-run

**How large the defect is, counted first.** This is the check that decides whether the change is
worth making, so it was run before anything was re-scored (`crispri.deletion_census`, folded into the
result file so it can be read back):

| arm | covered pairs | gene is the top target | structural zeros | of those, the cache scores | of those, a predicted fall |
|---|---|---|---|---|---|
| K562 training | 9,237 (451 reg.) | 246 | **8,991 (97.3%)** | 5,581 (62.1%) | 3,841 |
| K562 held out | 1,744 (114 reg.) | 40 | **1,704 (97.7%)** | 1,112 (65.3%) | 769 |
| GM12878 held out | 62 (14 reg.) | 8 | 54 (87.1%) | 53 (98.1%) | 33 |

So the answer to "does almost nothing change?" is **no, almost everything changes**. `deletion_drop`
was non-zero on 207 of 9,237 training pairs and on 35 of 1,744 held-out K562 pairs; after the fix it
carries a value the sweep actually predicted on 5,581 and 1,112 of them. Among the regulated pairs,
78 of the 114 held-out K562 positives were structural zeros and the cache answers 67 of them.

**Which direction the AUPRC is expected to move, and why.** Down, or to no change — not up. The
reason is in the same counts: `top_target` is a near-perfect separator on this set (36 of the 40
top-target pairs in held-out K562 are regulated, 90%, against 78 of 1,704, 4.6%), and the
2026-09-16 gain of +0.083 was carried by that indicator far more than by the magnitude beside it,
which was non-zero on 35 pairs. The fix leaves `top_target` exactly as it was and gives the other
1,704 pairs a magnitude that may be signal or may be noise. If the gain shrinks, the reading is that
the deletion's contribution was always "which elements the sweep selected" and never "how much the
gene moves"; if it grows, the cache's per-gene drops carry signal the compact table was censoring.
**Both outcomes are reported here with equal prominence**, in the same table, in the next section,
whichever way it goes. A correctness fix that lowers a headline is still a correctness fix.

**What is not touched.** The 2026-09-16 figures above are annotated, never rewritten: the
pre-registration `PREREGISTERED`, the feature set, the fitted model, the resampling and the held-out
split are all unchanged, and the only thing that changes is where `deletion_drop` reads its number.
`top_target` keeps its old meaning and stays a separate feature. `scripts/crispri_contact.py` and
`target_calibration.py` still read the compact table and their stored results are unchanged by this;
they are listed as limited consumers rather than silently re-run.

### The result: the registered direction was wrong on K562 and right on GM12878

`scripts/crispri_score.py` re-run the same way (107 s, **0 AlphaGenome requests**, same tables, same
`PREREGISTERED`, same split, same 200 resamples). The registration expected the gain to shrink or
hold. **On held-out K562 it nearly doubled**; on held-out GM12878 it fell by two thirds and its
interval now crosses zero. Both are below, in the same table, at the same size.

| held out | pairs (reg.) | activity + distance | + deletion, 2026-09-16 | + deletion, 2026-09-22 | gain then | gain now |
|---|---|---|---|---|---|---|
| K562 | 1,744 (114) | 0.550 | 0.633 | **0.691** | +0.083 (+0.031 to +0.166) | **+0.141 (+0.082 to +0.231)** |
| K562, elements not in training | 1,580 (91) | 0.497 | 0.587 | **0.646** | +0.090 (+0.037 to +0.173) | **+0.149 (+0.076 to +0.248)** |
| GM12878 | 62 (14) | 0.865 | 0.962 | **0.900** | +0.097 (0.000 to +0.227) | **+0.035 (−0.026 to +0.169)** |

K562 leave-chromosome-out on the training pairs moves the same way: 0.674 → **0.740** AUPRC with the
deletion, against an unchanged 0.511 without it, so the gain goes +0.163 (+0.121 to +0.210) →
**+0.228 (+0.185 to +0.277)**. The verdict is still **passed**, and it is the same pre-registration
that passed in September: nothing about the test was changed, only where one feature reads its number.

**The registered reasoning was half right.** The expectation rested on `top_target` carrying the old
gain, and that part holds — the indicator is untouched and its weight barely moves (1.284 → 1.329).
What the registration got wrong was calling the newly-visible magnitudes noise. They are not: the
fitted weight on `deletion_drop` rises 15.94 → **25.30**, which is the model leaning harder on a
column it can now read on 1,112 of the 1,704 held-out pairs it previously read as zero. The
predicted fall for a gene that is *not* the element's top target carries real signal about whether
silencing that element lowers that gene, and the compact table was throwing it away.

**GM12878 fell, and that is the honest caveat.** 0.962 → 0.900, gain +0.097 → +0.035 with an
interval from −0.026 to +0.169. It rests on 14 regulated pairs, 6 of them structural zeros of which
the cache answers 5, so a handful of newly-valued pairs move it either way; it was already the arm
whose interval touched zero in September. The arm that carries the weight of evidence is K562, with
114 regulated pairs of 1,744, and it rose. Neither arm is the whole answer and both are stated.

**What still cannot be read, and is not claimed.** Two diagnostics in the same result are gated on
`top_target` by construction and so barely move: the single predictor `top_target × (1 + drop)`
(AUPRC 0.4613 → 0.4647) and the one-call-per-element table (at threshold 0.1, 124 calls / 119 right
→ 128 / 123). Those numbers describe "which elements did the compact table select", which is a real
question with a correct answer, and they were not changed. The whole of the gain above lives in the
logistic model, where `deletion_drop` is a free column.

**A consistency check worth recording.** Over every covered pair on all 24 chromosomes, wherever the
compact table and the per-element cache both carry a value for the same (element, gene, cell), they
agree exactly: **0 disagreements**. The cache is not a different measurement or a later re-score; it
is the same numbers with the censoring removed. Where a top-target pair's drop changes at all, it is
because a *second* overlapping element also scored that gene, and the existing rule — the largest
fall across overlapping elements — can now see it.

**The other consumers of the one-target-per-element table**, found while doing this and left alone,
so the list exists: `crispri_direction.py:118` (`model_value`, superseded by
`crispri_direction_both.py`, which already reads the cache), `target_calibration.py:173`
(`element_for_pair`, the whole calibration is conditioned on `top_target`), `eqtl.py:246` (asks "is
the model's target an eGene" of one gene per element, and is directly improvable from the cache),
`motif_transfer.py:1025` and `syntax_tiling.py:252` (an element that moves a non-top gene reads as
"did not move"), `unknown_scoring.py:510`, `candidates.py:734`, `closure.py:69`, `compile.py:112`,
`organise.py:91`, `confidence_calibration.py:244`, `vista.py:278`, `decompile.py:113`,
`benchmark/loci.py:1137`, `cli.py:3832`, and the loader `attribution/targets.py:25` (`run_elements`)
that about twenty of them go through. `scripts/crispri_contact.py` is the nearest neighbour of this
fix and is deliberately not re-run here: its stored result still reads the compact table and says so.

## The eGene question asked of every element, not of the one gene a table kept (2026-09-22)

The fix above moved one consumer off the compact table. This section does two things: it puts the
reader under all of them — `ElementResponses` in `attribution/targets.py`, beside `run_elements` and
`attributed`, which about twenty modules load their elements through — and it moves the second
consumer, `eqtl.py`, with its own registration.

**The reader.** `ElementResponses.response(chrom, element_id, gene, cell)` returns a `Response` that
is either a signed log2 fold change or a `None` with the reason named: the gene is not in the
scorer's window, the gene was scored but not on that cell's own track, or the element is not in the
cache. It never returns a zero for a question it cannot answer, and `bool(response)` asks "did the
sweep answer", never "is the number non-zero", so a scored 0.0 is truthy and a silence is not.
`ranked(..., by="effect")` returns every gene in the window in the order
`predict.enhancer_target.predict_target` ranks them, so the head is the gene the compact table keeps
and the tail is what it dropped. `crispri.ElementCache` is the same reader written first, for the
CRISPRi benchmark; its `value()` and this one's agree by construction.

**What it costs**, measured: the tree is 775 MB gzipped over 24 chromosomes and must never be read
whole. One chromosome is held at a time and switching drops the previous one, so the bill is the
largest chromosome — chr21 (9.8 MB gzipped, 12,158 elements) 0.6 s and about 0.6 GB resident, chr1
(77 MB, 88,302 elements) 5.7 s and a peak near 2.9 GB. A caller therefore walks its elements
chromosome by chromosome; hopping costs a full decompression each time.

**That the cache is the same run, checked before any of it was used.** Over the 6,239 uniform,
constrained and VISTA elements that carry a prediction, the cache's effect-ordered head reproduces
the compact table's `predicted` gene **and** its log2 fold change exactly: 6,239 of 6,239, 0
disagreements. Where both carry per-cell values they agree exactly as well (800 comparisons, 0
disagreements). And every one of the 3,045 elements whose compact prediction is `null` has a cached
head below `MIN_EFFECT` = 0.1 — so `null` means "the best gene in the window moved by less than
0.1", not "nothing was scored". The cache is the censoring removed, not a re-score.

### The registration, written and committed before the eQTL score was re-run

`eqtl.py:246` asks, per element, "is the gene the model named among the genes GTEx ties to this
element?" — of one gene, chosen out of a window the sweep scored in full, and only where that one
gene cleared 0.1. Counted before anything is re-scored:

| set | elements with an eQTL | the question is asked of | never asked (no compact target) | genes in the window, median | eGene mentions in the window | elements with no eGene in the window |
|---|---|---|---|---|---|---|
| uniform | 3,725 | **2,372 (63.7%)**, rate 0.529 | **1,353 (36.3%)**, all cached | 30 (mean 34.1) | 15,292 of 17,834 (85.7%) | 112 (3.0%) |
| constrained | 1,394 | 1,054 (75.6%), rate 0.448 | 340, all cached | 29 (mean 33.1) | 4,993 of 5,881 (84.9%) | 65 (4.7%) |
| VISTA | 1,519 | 1,119 (73.7%), rate 0.412 | 400, of which 361 cached | 21 (mean 24.9) | 5,004 of 6,595 (75.9%) | 127 of the 1,480 cached |

Two things are wrong with the current number and they pull in opposite directions. **The denominator
is gated**: a third of the elements that carry an eQTL are dropped from it because the sweep's best
gene moved by less than 0.1, so the published 0.529 is the rate *given that the model was confident
enough to speak*, quoted as if it were the rate. **The negatives are conflated**: a `False` today
means either "the model named a different gene" or "the measured gene was never in the scorer's
window", and 14.3% of the uniform eGene mentions are of the second kind, with 112 elements that have
no answerable eGene at all and are counted as misses.

**The direction expected, and why.** Down. Adding the 1,353 ungated elements should lower the
headline, because on those the sweep's best gene moves by less than 0.1 and a top-1 pick out of a
30-gene window at that size is close to guessing; removing the 112 unanswerable elements should
raise it a little, and 1,353 outweighs 112. So `cache_target_is_an_egene` over every answerable
element is registered to come in **below 0.529 on uniform**, and the same way on the other two sets.
The floor to compare against is picking a gene at random from the window, about 4.1 eGenes in a
window of 34, near 0.12; the number to beat is the nearest-coding-TSS rule, 0.662, which already
beats the model on this test.

**The falsifier.** If the ungated elements are right about as often as the gated ones — if the
headline holds at 0.529 or rises — then `MIN_EFFECT` is not selecting the elements the model can
answer, and the gate is not earning the conditioning it imposes. That would be a finding about the
threshold and it is to be reported as one.

**Reported first if it falls.** A fall is the expected outcome and goes at the top of the result
section, in the same table as anything that rises.

### The result: the headline did not move, and that is the finding

`scripts/eqtl_targets.py` re-run over the same elements and the same distilled hits (125 s, **0
AlphaGenome requests**, nothing fetched). The registration expected the headline to fall. **On the
uniform set it moved by 0.001**, from 0.529 to 0.528, while the denominator grew by half, from
2,372 elements to 3,613. On VISTA it rose, 0.412 → 0.423. So the registered direction is wrong, and
the falsifier written beside it is the one that fired: **the elements the compact table refused to
name a gene for are answered about as well as the elements it named one for.**

| set | published, gated | asked of every answerable element | where the table named a gene | where it named none | chance in the same window |
|---|---|---|---|---|---|
| uniform | 0.529 (2,372) | **0.528 (3,613)** | 0.546 (2,297) | **0.496 (1,316)** | 0.127 |
| constrained | 0.448 (1,054) | **0.449 (1,329)** | 0.467 (1,010) | **0.392 (319)** | 0.116 |
| VISTA | 0.412 (1,119) | **0.423 (1,353)** | 0.441 (1,046) | **0.362 (307)** | 0.146 |

The movement decomposes into two corrections that nearly cancel, and both are real. Uniform:
0.529 over 2,372 → **0.546 over 2,297** when the 75 elements whose every eGene lies outside the
scorer's window are dropped, because the model could not have named one and they were being counted
as misses (+0.017) → **0.528 over 3,613** when the 1,316 elements the table left unnamed are added
at 0.496 (−0.018). Net −0.001. Quoting the old 0.529 as the model's accuracy was wrong twice over,
and the two errors happened to be the same size.

**What `MIN_EFFECT` = 0.1 is worth, now that it can be priced.** The gate excludes a third of the
elements that carry an eQTL. On the elements it keeps, the sweep's best gene is an eGene 54.6% of
the time; on the elements it throws away, 49.6% — five points, against a chance floor of 12.7% in
the same windows. So the threshold is a weak confidence signal, not the difference between an
answer and a non-answer, and a table built on it is discarding 1,316 usable answers to raise a rate
by five points. On the constrained and VISTA sets the same gap is 7.5 and 7.9 points.

**What the compact table could not express at all.** The measured eGene's rank in the sweep's own
ordering of the window, not merely whether it is the head of it: median rank **1** of 34.1 genes on
uniform, with the eGene in the model's top 3 for **77.6%** of elements and the top 5 for **85.1%**
(constrained 0.673 / 0.737 of 33.1 genes, VISTA 0.694 / 0.793 of 24.9). That is a much stronger
statement about what the sweep knows than the 0.528 top-1 rate, and it is not a comparison with the
nearest-coding-TSS rule, which gets one guess and takes 0.662; it is a statement that when the model
is wrong about which gene, it is usually wrong by one or two places.

**The negatives that were never answers.** 2,542 of the 17,834 uniform eGene mentions (14.3%) are
of genes outside the scorer's 1 Mb window, and on 112 elements *every* eGene is outside it. Those
elements now carry `cache_silence`, "no eGene of this element is in the scorer's window", instead of
a `False`. On VISTA 39 elements are not in the cache at all and say so rather than being scored.

**The control held exactly.** Every pre-existing field reproduces to the digit: uniform 0.529 on
2,372, coding 0.717 on 1,851, nearest TSS 0.662 on 3,725, and the disagreement table 516 / 842 of
1,379. The only reason any number in the result file differs from the September run is the block of
new fields beside it (3.9 MB → 6.2 MB).

**A correction to the record.** The commit that carries the reader and the registration above,
45e4f92, was committed with the wrong message text: it repeats the message of df6c5f4, from the
CRISPRi lane earlier the same day, because a stale message file of that name was picked up. Its
contents are `attribution/targets.py`, `tests/test_element_responses.py` and the registration
section here, and the history was left unrewritten rather than rewritten under a shared branch.

**What is not touched.** Every field the current result carries keeps its name and its meaning, so
the figures above stay comparable and the web view keeps working; the new fields are added beside
them. `predicted_target_is_an_egene` still means "of the elements where the compact table named a
gene", and the re-run must reproduce it exactly — that is the control. Zero AlphaGenome requests:
every number here is on disk. The rank of the measured eGene in the model's own ordering, which the
compact table could not express at all, is reported as new information rather than as a comparison.

## Pre-registration: the calibration is conditioned on `top_target`, counted before it is re-fitted (2026-09-22)

The published calibration ("A measured confidence for a predicted target, and the prevalence that
breaks it", 2026-09-17) attaches a probability to every predicted enhancer-to-gene target, and it is
the module with the most at stake in the compact table's one-gene projection: it is what downstream
work quotes when it says how sure the project is that an element acts on a gene. Its pre-registered
reliability claim already failed once, and it failed on the population rather than on the curve — the
held-out screens call 6.53% of pairs regulated against 4.97% in the screens the fit saw, and one
log-odds shift of +0.679 restores 9 of the 10 bins. This section counts the gate before anything is
re-fitted; `PREREGISTERED_GATE` in `attribution/target_calibration.py` carries the same text, fixed
in the module before the re-fit ran. `scripts/target_calibration_gate.py --census`, 137 s, 0
AlphaGenome requests, nothing fetched.

**Where the gate acts — three places, not one.** `top_target` is the flag that a pair's measured gene
is the single gene the compact `all_elements` table kept for that element.

1. *The magnitude.* `deletion_drop` is read from that table, so it is zero for every pair whose gene
   is not that one gene — a statement about the projection, consumed as a statement about the model.
   On the fitted population that is 8,589 of 8,796 K562 training pairs (97.6%) and 1,680 of 1,715
   held-out K562 pairs (97.9%).
2. *The element.* `matched_element` (`target_calibration.py:164`, the function the assignment names
   `element_for_pair`) picks the overlapping element whose top predicted target is the pair's gene,
   and failing that the element with the largest predicted magnitude **for whatever other gene the
   table named** — and the class and distance features are then taken from that element.
3. *The population.* The predicted-target calibration is fitted on the 245 training pairs the gate
   admits and read on 40 held-out ones, and it is that fit whose weights band all 612,323 sweep
   targets.

**What the gate excludes, and what the excluded pairs are when asked properly.** Through
`targets.ElementResponses`, which carries every gene in the scorer's 1 Mb window with a signed change
on the cell's own track:

| | training, every feature present | held-out K562 |
|---|---|---|
| pairs the calibration is fitted or judged on | 8,796 (437 regulated, 4.97%) | 1,715 (112, 6.53%) |
| admitted by the gate | 245 (188, **76.7%**) | 40 (36, **90.0%**) |
| excluded by the gate | 8,551 (249, **2.91%**) | 1,675 (76, **4.54%**) |
| of the excluded, the sweep did predict a change for this gene | 5,571 (65.1%) | 1,112 (66.4%) |
| of the excluded, the gene is not in the scorer's window at this element | 2,980 (34.9%) | 563 (33.6%) |
| of the excluded, not on this cell's track, or not cached | 0 | 0 |

So the gate is two thirds a censoring the sweep's own cache can undo and one third a real limit, and
the calibration can be fitted with a magnitude that means something on **5,816 training pairs instead
of 245**, a factor of 23.7. The base rate on the two sides differs by a factor of 26, which is the
whole difficulty: `top_target` separates this set almost perfectly on its own.

**How far down the window the screens test.** New, and not expressible in a one-gene table: the
measured gene's rank in the sweep's own ordering of the window. Median rank **24** of a median 50
genes in the window on the training pairs (53 held out); rank 1 on 128 of 5,814, in the top 3 on 359
and the top 5 on 566. A CRISPRi screen tests genes the model is not confident about, which is exactly
why a calibration fitted only on the ones it *is* confident about is the wrong curve for them.

**What is re-fitted and what is not.** The published curve is annotated, never rewritten. `score()`
keeps its default of no cache and reproduces the 2026-09-17 numbers exactly — that is the control —
and the re-fit is a second arm that passes the cache. `SWEEP_FEATURES` and `TARGET_FEATURES` do not
change, because a feature the sweep does not have cannot be fitted; only what `deletion_drop` *means*
changes, from "the compact table's entry for the one gene it kept" to "what the sweep predicted for
this pair's own gene at this element".

**The direction expected, and why.** Reliability is expected **not** to improve. The 2026-09-17 claim
failed on the population: 6 of 10 bins, all four failures in the same direction, mean predicted
0.0441 against an observed 0.0653. An uncensored magnitude is a better feature; it is not an
intercept, and it cannot move a prevalence that is a property of how a screen chose its pairs. The
registered expectation is 6 of 10 bins give or take one, the prevalence ratio still near 1.31, and a
rise in AUPRC as the only movement a better feature buys.

**What a confidence means for a pair the gate would have excluded.** This is the clause that decides
whether the re-fit is worth anything. The predicted-target curve is fitted on 245 pairs whose base
rate is 76.7% — the model's most confident calls, the one gene per element it was surest of — and the
sweep quotes it for 612,323 targets. A user who now asks about any *other* gene in the window, which
the reader makes askable and which is about fifty genes per element rather than one, would be quoted
a curve fitted on a population whose base rate is twenty-six times theirs. That is the classic form
of a calibration that looks reliable and is not. The registered measurement is direct: score the
shipped 245-pair curve on the held-out pairs the gate excluded but the sweep did answer, and report
its mean predicted probability against their observed rate.

- **What would show the failure is happening:** the shipped curve's mean predicted probability on
  those pairs sits far above their observed rate — a gap of the order of the 0.767-against-0.029
  base-rate gap, not of the 1.31 prevalence ratio already found — and its bins fall outside their
  intervals in one direction.
- **What would show it is not:** the shipped curve lands near their observed rate, which would mean
  the three features carry the population difference and the gate was only selecting on them.

Either way the number is reported, and a confidence for an off-gate pair is quoted from a curve
fitted on off-gate pairs or it is not quoted at all.

**If reliability improves.** That is the outcome that would tempt a lane to stop checking, so the
checks are fixed here and run whether it improves or not. A curve that predicts the base rate
everywhere is trivially inside every equal-count bin, so a rise in `bins_consistent` is reported only
beside (a) the width of the predicted range across the ten bins, which must not shrink, (b) AUPRC on
the same held-out pairs, which must not fall, and (c) the prevalence ratio and the log-odds shift,
which must have moved towards 1 and 0. If `bins_consistent` reaches 7 or more while the predicted
range narrows or AUPRC falls, the improvement is recorded as a **flattening** and the 2026-09-17
verdict of failed is not upgraded. The verdict is upgraded only if the prevalence gap itself closes,
and nothing in this change acts on the intercept.

## The gate removed: the registered direction was wrong, and the confidence quoted off the gate was seven times the measured rate (2026-09-22, later)

> **Annotation, 2026-09-22 (sixth).** The registered condition for upgrading this section's verdict
> was that the prevalence gap itself close. It has now been attacked directly with per-screen
> intercepts and does not close: the repair fails 3 of its own 4 registered clauses, and the gap
> was never the two-population fact it was stated as -- it is a per-screen level that moves in both
> directions inside each table. See the sixth section below.

> **Annotation, 2026-09-22 (third).** Every number in this section stands as measured. Two sentences
> of its prose are narrowed by the inventory in the section below, and neither changes a measurement.
> (1) "its weights band all 612,323 sweep targets" — 593,765 are banded; the other 18,558 carry
> `no_gencode_tss_for_the_predicted_gene` and get no band. (2) The ×6.85 is read on pairs the gate
> excludes, and **no swept target is off the gate**: `crispri.deletion_values` sets `top_target` for
> a match on either `predicted` or `predicted_coding`, and `sweep_chromosome` bands exactly those two
> keys. So ×6.85 is the error the module would make on the newly askable off-gate population, not an
> error in the published band table, and the correction to that table is a different one.

The registration above expected reliability not to improve. It improved, from 6 of 10 bins to 7,
which is exactly the pre-registered threshold, so `judge()` reads **passed** where it read failed.
That is the outcome the registration named as the one that would tempt a lane to stop checking, so
the checks fixed in advance were run, and what they show is that the improvement is real but is in
**ranking, not in level**: the prevalence failure the 2026-09-17 section diagnosed is untouched. The
larger result is elsewhere. `scripts/target_calibration_gate.py`, result `target_calibration_gate`,
518 s, **0 AlphaGenome requests**, nothing fetched.

**The control held exactly.** The shipped arm annotates from the compact table, as `score()` does,
and reproduces the published numbers to the digit: 8,796 training and 1,715 held-out K562 pairs, 6 of
10 bins, ECE 0.0231, Brier 0.04184, AUPRC 0.559, mean predicted 0.0441 against an observed 0.0653, a
log-odds shift of +0.6793 restoring 9 of 10. The held-out population is identical in the two arms
(1,715 pairs, 112 regulated), so the bin counts are comparable.

| held-out K562, the same 1,715 pairs | bins inside | ECE | Brier | AUPRC | `judge()` |
|---|---|---|---|---|---|
| shipped, the compact table (published 2026-09-17) | 6/10 | 0.0231 | 0.04184 | 0.559 | failed |
| re-fitted, the sweep's own window | **7/10** | 0.0214 | **0.03453** | **0.677** | passed |

**Why the pass is not an upgrade.** All three bins still outside fail in the same direction as before,
the observed rate above the predicted one (bins 1, 8 and 10; bin 10 predicts 0.3570 against an
observed 0.4593). The mean predicted probability moves 0.0441 → 0.0469 against an unchanged observed
0.0653, so the prevalence ratio only goes 1.48 → 1.39, and the single log-odds shift that fixes it
goes **+0.6793 → +0.7025 — larger, not smaller** — and still restores 9 of 10 bins. The registered
condition for upgrading the verdict was that the prevalence gap itself close; it did not, and nothing
in this change acts on the intercept. **The 2026-09-17 verdict of failed stands.** What moved is the
ranking: AUPRC +0.119 and Brier −17.5%, and the extra bin is bin 9 crossing its interval because the
better ranking reshuffled which pairs land in it. The registered flattening test clears the
improvement of the other explanation — the predicted range across the ten bins *widened*, 0.3007 →
0.3521, and AUPRC rose, so this is not a curve collapsing onto the base rate — but a real ranking
gain and an unmoved level are two different things and only the first happened.

The two deletion terms are not redundant once the magnitude is real: `deletion_drop` takes 20.72 →
**31.57** while `top_target` keeps its own weight (1.234 → 1.355) and the intercept falls 7.16 → 4.79.
The flag says the table named this gene; the magnitude says how far the sweep thinks it moves, and
before today the second was only ever spoken on the 2.8% of pairs where the first was true.

**The registered clause, measured: what a confidence means off the gate.** The predicted-target curve
is fitted on 245 pairs at a base rate of 76.7%, and its weights band all 612,323 sweep targets. Asked
about the 1,112 held-out pairs the gate excludes but the sweep does answer:

| curve | pairs | quoted | observed | 95% interval | ratio |
|---|---|---|---|---|---|
| the predicted-target curve, on the pairs the gate excludes | 1,112 | **0.4126** | **0.0603** | 0.0477–0.0758 | **×6.85** |
| the same curve, on the 40 pairs it admits | 40 | 0.8355 | 0.9000 | 0.7695–0.9604 | ×0.93 |
| the tested-pair curve, on the same excluded pairs | 1,112 | 0.0359 | 0.0603 | 0.0477–0.0758 | ×0.60 |

This is the signature the registration named for the failure happening, and it is not marginal. On
its own population the curve is honest to within 7%; one step off it, it quotes seven times the
measured rate. Re-fitting it through the window barely helps (0.4126 → 0.3855, ×6.40), because the
problem is the population it was fitted on and not the features it was fitted with. A confidence for
a pair the gate would have excluded, quoted from this curve, is not a probability of anything.

**And it is wrong twice, because the off-gate set is not flat.** With the sweep's own magnitude
restored, the pairs the gate excludes separate by more than fifty-fold. Training and held-out pooled,
legitimate only because the held-out arm above was scored first:

| predicted K562 drop, on genes that are *not* the element's top target | pairs | regulated | rate | 95% interval | elements |
|---|---|---|---|---|---|
| = 0 (the sweep predicts a rise or no change) | 2,079 | 33 | **0.0159** | 0.0113–0.0222 | 1,286 |
| 0 < drop ≤ 0.1 | 4,537 | 209 | 0.0461 | 0.0403–0.0526 | 2,363 |
| 0.1 < drop ≤ 0.2 | 34 | 30 | 0.8824 | 0.7338–0.9533 | 33 |
| 0.2 < drop ≤ 0.5 | 29 | 24 | 0.8276 | 0.6545–0.9240 | 26 |
| 0.5 < drop | 4 | 4 | 1.0000 | 0.5101–1.0000 | 3 |

A drop above 0.1 on a gene the compact table never named was called regulated **58 times out of 67
(86.6%)**, against 33 of 2,079 (1.6%) where the sweep predicted a rise or nothing — a fifty-five-fold
separation on a population that until today scored a structural zero and would have been quoted 0.41
by the sweep's own curve. The model names second and third targets, and the projection was throwing
them away. The caution is the size: 67 pairs above 0.1 off the gate, 34 of them held out, so the top
three rows are a direction with an interval and not a rate to quote to three figures.

**The base-rate problem, re-read on the wider set.** It does not go away and it is not an artefact of
the gate. The held-out screens call pairs regulated more often than the training screens on every
slicing: 6.53% against 4.97% on the published population (×1.31), 8.94% against 7.24% on the answered
population (×1.23), 6.03% against 4.18% on the answered pairs the gate excludes (×1.44). Widening the
population by a factor of 23.7 changed the prevalence ratio by less than a fifth of itself. The
2026-09-17 conclusion stands unchanged and is now measured on 5,816 training pairs rather than 245: a
band quoted without its population's base rate beside it measures nothing.

**Coverage: what fraction of the questions can now be answered at all.**

| | askable before | askable now | factor |
|---|---|---|---|
| training K562 pairs on a deleted element | 246 of 9,237 (2.7%) | 5,827 (63.1%) | ×23.7 |
| held-out K562 pairs on a deleted element | 40 of 1,744 (2.3%) | 1,152 (66.1%) | ×28.8 |
| chr21 (element, gene) questions | 7,846 | 347,591 | ×44.3 |
| chr22 (element, gene) questions | 13,108 | 873,596 | ×66.6 |

The benchmark's remaining third is a real limit and now says so: 2,980 training and 563 held-out
pairs carry the named silence "the gene is not in the scorer's window at this element", and not one
pair is a "not on this cell's own track" or a "not cached". The sweep rows are two chromosomes, not
genome-wide — the response tree is 775 MB and one archive is held at a time — and they count every
gene in the scorer's window with a K562 value, a median of 28 and 41 genes per element against the
one the compact table kept.

**What is not touched.** `score()` keeps its default of no cache and its published output; the
2026-09-17 result file and the genome-wide band table are unchanged and were not re-run. The feature
columns are the published ones. `matched_element` takes the reader as an optional argument and
reproduces its old choice exactly without one, so the control above is the shipped code path and not
a reconstruction of it.

**Next, proposed as roadmap rows** (for the roadmap editor to fold; not written into ROADMAP.md):

1. Area I: the sweep's bands are produced by the 245-pair curve, which overstates by ×6.85 one step
   off its population. Re-band the genome from the curve fitted on the answered pairs, and quote an
   off-gate confidence only from an off-gate curve.
2. Area I: 67 off-gate pairs above a drop of 0.1 carry 86.6% regulated. That is the shortlist to
   hold against an independent perturbation — the model's second targets, which no compact table has
   ever offered.
3. Area I: the prevalence term, unchanged from 2026-09-17 and now confirmed on 23.7 times the pairs.
   Re-fit the intercept per screen with the screen's own base rate as an offset before any band is
   quoted.

## Pre-registration: band the genome from a curve fitted on the population each band is quoted for (2026-09-22, third)

> *Annotation, 2026-09-22 (fifth), added later and changing nothing below.* This registration's
> falsifier fired on the axis, and the axis that replaced it was in turn superseded by their union,
> which is registered separately in `PREREGISTERED_UNION` and measured two sections down. Nothing
> in this registration is withdrawn: its direction, its count of 34,158 and its thin-stratum rule
> all survive and are reused unchanged by the union registration.

The section above found that the curve pricing every predicted enhancer target is fitted on 245
pairs and quotes ×6.85 one step off its own gate, and proposed re-banding the genome as the first
roadmap row that follows. This registers that re-banding before any target is re-banded, and it
begins by establishing that the repair everybody expected is not the repair the sweep needs.
Inventory and counts: 0 AlphaGenome requests, everything read from disk.

**Where the band table is, and what it is.** `sweep_chromosome()` and `sweep()` in
`genomeos/attribution/target_calibration.py` produce it; `scripts/target_calibration.py` drives them
and writes `data/results/target_calibration.json`. It is **not** a per-element row table — the
per-(element, gene) probability is computed in memory and never persisted. What is stored is an
aggregated counter, keyed on **(chromosome) × (`any gene` | `coding gene`) × band label → count**,
with a second cut by ENCODE registry class, and `genome_wide` is the same structure summed over the
24 chromosomes. `band_of()` buckets a probability into the eight half-open `CONFIDENCE_BANDS`.

**How many targets carry which band today**, from `genome_wide.predicted_target_bands`:

| band | `any gene` | `coding gene` |
|---|---|---|
| 0.05–0.1 | 45 | — |
| 0.1–0.25 | 3,837 | 1 |
| 0.25–0.5 | 178,387 | 184,256 |
| 0.5–0.75 | **293,568** | 175,655 |
| 0.75–0.9 | 61,971 | 42,552 |
| 0.9–1 | **55,957** | 37,913 |
| total banded | **593,765** | 440,377 |

**The 612,323 is not the banded count and never was.** 593,765 banded + 18,558
`no_gencode_tss_for_the_predicted_gene` = 612,323, and 612,323 + 348,904 `no_predicted_target` =
961,227 elements. So **18,558 of the 612,323 (3.0%) carry no band at all**, and the sentence "the
weights band all 612,323 sweep targets" — in this document twice above, in `LESSONS.md`, in
`ROADMAP.md` and in the module's own prose — overstates by that many. Nothing downstream consumes
the table: `genomeos/web/server.py` and `genomeos/web/static/index.html` contain no band label and
no reference to this result, `genomeos/cli.py` has no command that prints it, and the only reader is
`scripts/target_calibration.py`'s own stdout. **No web UI hunk is needed for the correction itself.**

**The expected repair does not exist in the sweep.** The obvious split — band on-gate and off-gate
targets from two curves — assumes some swept targets are off-gate. None are.
`crispri.deletion_values` sets `top_target` to 1.0 when the pair's gene matches **either**
`predicted` **or** `predicted_coding` at an overlapping element, and `sweep_chromosome` bands
exactly those two keys. Of the 440,377 coding-arm targets, 331,209 are the element's window head and
109,168 are the coding head only — a real distinction, and not this one, because both are compact
table entries and both carried `top_target` = 1 inside the fit. **The ×6.85 is therefore not an
error in the published table.** It is the error the module would make on the population
`targets.ElementResponses` has just made askable — the median 28 to 41 other genes per element — and
which the sweep does not yet band at all.

**What the population mismatch actually is.** The fitted 245 and the banded 593,765 are both inside
the gate and are still different populations, on two axes measured before this registration:

| | fitted 245 | swept targets | factor |
|---|---|---|---|
| predicted K562 drop above 0.2 | 34.3% | 5.8% | ×5.9 |
| predicted K562 drop exactly 0 | 15.5% | 46.5% | ×0.33 |
| target named on **K562**, the one cell this curve speaks for | 78 (31.8%) | 32,597 of 612,323 (5.3%) | **×6.0** |

The other 94.7% were named on placenta, CD14-positive monocyte, HepG2, testis, psoas muscle and
some three hundred other tracks, and enter a K562 curve through a K562 drop that is zero for most
of them.

**The populations, and how a target is assigned to one.** The axis is the predicted K562 deletion
drop, in the five `DROP_BANDS` strata the module already reports `measured_by_drop_band` on. It is
chosen over the top-target axis because it is the axis the two populations demonstrably differ on,
it is computable for every swept target and every benchmark pair with no extra data, and it is the
axis along which the fitted rate moves. Registered in advance, `any gene` arm:

| stratum | swept targets | training pairs | held-out | observed rate (pooled) |
|---|---|---|---|---|
| = 0 | 275,821 (46.5%) | 38 | 5 | 0.5814 |
| 0 < drop ≤ 0.1 | 239,487 (40.3%) | 83 | 9 | 0.5978 |
| 0.1 < drop ≤ 0.2 | 44,299 (7.5%) | 40 | 5 | 0.8667 |
| 0.2 < drop ≤ 0.5 | 24,281 (4.1%) | 35 | 11 | 1.0000 |
| 0.5 < drop | 9,877 (1.7%) | 49 | 10 | 1.0000 |

**86.8% of the genome's bands rest on 121 training and 14 held-out pairs.** The `coding gene` arm is
205,541 / 173,234 / 37,488 / 18,338 / 5,776 on the same strata.

**The rule.** A stratum keeps a numeric band only if it holds at least `MIN_POOLED_FOR_A_BAND` = 30
pooled pairs **and** the 95% Wilson interval of its pooled observed rate lies inside a single one of
the eight `CONFIDENCE_BANDS`. Otherwise its targets are banded `not calibrated here` and carry the
stratum's pair count, observed rate and interval instead of a number. Pooling training with held-out
pairs for the rate is declared here rather than discovered later: the 40 held-out pairs alone give
intervals of width 0.65 and would mark every stratum uncalibrated, which is true and says nothing,
so the interval is read on all 285 on-gate pairs, and what is being interval-bounded is a base rate
rather than a fitted curve's error. **If a stratum is too thin to fit at all it is banded `not
calibrated here` with its count and no curve is fitted for it** — the expected outcome at drop = 0,
whose 43 pooled pairs give 0.5814 [0.4335, 0.7160], spanning two published bands.

**The direction expected.** Most of the genome loses its number: the two bottom strata (515,308 of
593,765, 86.8%) and the 0.1–0.2 stratum (44,299, 7.5%) all fail the interval clause and become `not
calibrated here`, leaving only the two top strata (34,158 targets, 5.8%) with a numeric band, which
will be 0.9–1 for both. That is 80% to 95% of the table losing its band. For every target that keeps
one the direction is unchanged or up, never down — because **within the gate the curve was measured
honest in every stratum**: quoted 0.4932 against an observed 0.5789 at drop = 0, 0.6139 against
0.5663, 0.8746 against 0.8750, 0.9808 against 1.0000, 0.9999 against 1.0000. So the correction is
**not** that the published levels are wrong. It is that the table states eight-way bands to two
decimal places for 593,765 targets on the evidence of 285 pairs.

**The falsifier.** The split is wrong if the two bottom strata's pooled intervals do fit inside a
single published band, because the drop axis would then separate nothing the curve has not already
absorbed and the published table would stand as written. It is also wrong if the alternative axis
measured here — named on K562 against named on another track — yields single-band intervals where
the drop axis does not; both are computed, both are reported, and if the cell axis is cleaner it
replaces the drop axis and this registration is recorded as wrong on the axis while right on the
direction. Third: if re-banding moves fewer than half the targets, the registered magnitude was
wrong and is reported as wrong rather than rounded towards.

**What is not touched.** The published band table is annotated and kept, never overwritten. `sweep()`
keeps its output and its keys and the re-banding is a second block under a new key, so both can be
read side by side from the same result. The eight `CONFIDENCE_BANDS` and the fitted weights do not
change. Registered in full as `PREREGISTERED_BANDS` in `genomeos/attribution/target_calibration.py`.

## The genome re-banded: every registered count correct, and the registration's own falsifier fired on the axis (2026-09-22, fourth)

> **Annotation, 2026-09-22 (sixth).** The published band table this section reproduces to the digit
> -- 45 / 3,837 / 178,387 / 293,568 / 61,971 / 55,957 -- is, to within **0.37%**, the band table of
> one screen. Banded under Gasperini2019's own on-gate base rate it reads 32 / 3,843 / 176,025 /
> 295,305 / 62,394 / **56,166**, and that screen is 207 of the 285 on-gate pairs. Under the other
> seven screens of the same benchmark the top band runs from 27,514 to 198,475. So this section's
> 94.2% that lose their band were never carrying a property of the genome. See the sixth section
> below.

The registration above was met on every quantity it predicted and defeated on the one thing it was
least sure of. `scripts/target_rebanding.py`, result `target_rebanding`, **0 AlphaGenome requests**,
nothing fetched.

**The control held.** The re-banding walk recomputes the published band for every target as it goes,
and reproduces the published table to the digit: 593,765 `any gene` targets at 45 / 3,837 / 178,387 /
293,568 / 61,971 / 55,957 and 440,377 `coding gene` targets at 1 / 184,256 / 175,655 / 42,552 /
37,913. So the two bands below are read off the same walk over the same elements.

**The populations, measured.** Pooled over the 245 training and 40 held-out on-gate pairs:

| predicted K562 drop | pairs | observed | 95% interval | band it can carry |
|---|---|---|---|---|
| = 0 | 43 | 0.5814 | 0.4333–0.7162 | **not calibrated here** (spans 2 bands) |
| 0 < drop ≤ 0.1 | 92 | 0.5978 | 0.4957–0.6922 | **not calibrated here** (spans 2 bands) |
| 0.1 < drop ≤ 0.2 | 45 | 0.8667 | 0.7382–0.9374 | **not calibrated here** (spans 3 bands) |
| 0.2 < drop ≤ 0.5 | 46 | 1.0000 | 0.9229–1.0 | 0.9–1 |
| 0.5 < drop | 59 | 1.0000 | 0.9389–1.0 | 0.9–1 |

**The re-banding, on the registered axis.** `any gene`: **559,607 of 593,765 targets (94.2%) lose
their band**, 34,158 (5.8%) keep one and it is 0.9–1 — 34,096 unchanged, **62 up, none down**.
`coding gene`: 416,263 of 440,377 (94.5%) lose theirs, 24,114 keep 0.9–1, 24,100 unchanged, 14 up,
none down. The registration predicted 80–95% losing their band, the two top strata keeping 0.9–1,
and the count 34,158; all three are exact, and so is the direction clause.

**What the bands were actually called by the screens, per population.** This is the number the
correction is for. The published curve, scored on the pairs it bands and on the pairs it does not:

| published band | on the gate: pairs → observed | off the gate: pairs → observed | ratio |
|---|---|---|---|
| 0.25–0.5 | 37 → **0.3514** | 9,599 → **0.0188** | ×18.7 |
| 0.5–0.75 | 89 → **0.6629** | 627 → **0.2313** | ×2.9 |
| 0.75–0.9 | 37 → 0.8919 | never reached off the gate | — |
| 0.9–1 | 122 → **0.9754** | never reached off the gate | — |

**On the gate every published band is honest** — 0.3514 inside 0.25–0.5, 0.6629 inside 0.5–0.75,
0.8919 inside 0.75–0.9, 0.9754 inside 0.9–1. The ×6.85 of the section above is therefore not a level
error in the shipped table; it is what the same label means one population over, where 0.25–0.5 buys
1.9% and 0.5–0.75 buys 23%. **After re-banding, the top band means 105 of 105** (1.0000, 0.9647–1.0):
the screens called every pair in the two top drop strata regulated.

**The off-gate population now has a table of its own**, read through the sweep's own window so its
magnitude is real — the thing the shipped curve had no honest number for:

| off-gate, predicted K562 drop | pairs | observed | band |
|---|---|---|---|
| = 0 | 4,716 | 0.0091 | **0–0.02** |
| 0 < drop ≤ 0.1 | 3,783 | 0.0428 | **0.02–0.05** |
| 0.1 < drop ≤ 0.2 | 26 | 0.8462 | not calibrated here (26 pairs) |
| 0.2 < drop ≤ 0.5 | 24 | 0.8333 | not calibrated here (24 pairs) |
| 0.5 < drop | 2 | 1.0000 | not calibrated here (2 pairs) |

An off-gate target is banded 0–0.02 or 0.02–0.05 where the shipped curve quoted **0.4126**. And the
registered thin-stratum rule does its work on precisely the interesting rows: the 52 off-gate pairs
above a drop of 0.1 carry 84.6% and 83.3% and are still refused a number, because 26 and 24 pairs
cannot carry one. That is the model's second and third targets, and the honest answer there today is
a direction with an interval, not a band.

**The falsifier fired, on the axis.** The registration named one rival — the track the sweep named
the target on — and said that if it gave single-band intervals where the drop axis did not, it
replaces the drop axis and the registration is wrong on the axis while right on the direction. It
did: K562-named, 93 pairs, 1.0000 [0.9603, 1.0] → 0.9–1; named on another track, 192 pairs, 0.6823
[0.6134, 0.744] → 0.5–0.75. **Both strata carry a band, so the track axis bands 100% of the genome
where the drop axis bands 5.8%.** The registered verdict is recorded as written: **wrong on the
axis.**

**And on that axis the published top band is the band that falls.** Re-banded on the track,
`any gene` goes to 561,927 at 0.5–0.75 and 31,838 at 0.9–1: 194,125 targets up, 298,948 unchanged
and **100,692 down**. The fall is concentrated exactly where the table was most confident —
**44,372 of the 55,957 targets published at 0.9–1 (79.3%) drop to 0.5–0.75**, and 56,320 of the
61,971 at 0.75–0.9 (90.9%) drop with them, while 170,069 of the 178,387 at 0.25–0.5 rise to
0.5–0.75. `coding gene` is the same shape: 29,052 of 37,913 at 0.9–1 (76.6%) fall. The registered
axis reaches the same targets by a different route — of the 55,957 at 0.9–1 it leaves 34,096 there
and declares **21,861 (39.1%) not calibrated at all**. Whichever axis is read, **the claim that
55,957 element–gene pairs are at least 90% likely does not survive contact with the population it
was quoted for**; between a third and four fifths of that band is either uncalibrated or two bands
lower.

**Why neither axis is the answer, which is new and was not registered.** The two axes are nearly
independent, not two readings of one thing. Genome-wide, 32,597 targets are K562-named and 35,217
carry a drop above 0.2, and only **8,572 are both** — a Jaccard overlap of 14.5%. Crossed on the 285
on-gate pairs:

| named on K562 | drop > 0.2 | pairs | regulated | rate |
|---|---|---|---|---|
| no | no | 157 | 96 | 0.6115 |
| no | yes | 35 | 35 | **1.0000** |
| yes | no | 23 | 23 | **1.0000** |
| yes | yes | 70 | 70 | **1.0000** |

Either signal alone is sufficient and the two are close to disjoint, so the union is the better
split: **128 of 128 pairs with either signal were called regulated** (1.0000, 0.9709–1.0 → 0.9–1),
against 96 of 157 with neither (0.6115, 0.5334–0.6842 → 0.5–0.75). That bands the whole genome from
two populations that each carry a single-band interval on more than 128 pairs — of the 612,323 named
targets, 59,242 carry either signal and 553,081 carry neither, before the 18,558 without a GENCODE
TSS are removed from the banded set. **It is not adopted here**, because it was not registered and
picking it now would be the same post-hoc choice this lane exists to correct. It is the next
registration, and it is the row sent to the roadmap.

> *Annotation, 2026-09-22 (fifth), added later and changing nothing above.* That next registration
> was written and run. The union axis **separates** on two populations that took no part in
> choosing it, so the refusal here was right to hold it open rather than to kill it — but the same
> test shows the `1.0000` in the table above is a property of the population and not of the axis:
> off the gate the identical union stratum reads **0.0895**, a factor of 11 lower. The counts here
> are over all 612,323 named targets; over the 593,765 actually banded they are 31,838 K562-named,
> 34,158 above a drop of 0.2, 8,370 both and **57,626 in the union**. See the section below.

**What is not touched.** `sweep()` keeps its output and its keys and the published band table is
unchanged on disk; the re-banding is a separate result carrying the published band beside the new
one, on both axes, so a reader can take either and compare. The weights, `CONFIDENCE_BANDS` and
`TARGET_FEATURES` do not change. No web UI or CLI hunk is needed: `genomeos/web/server.py`,
`genomeos/web/static/index.html` and `genomeos/cli.py` contain no reference to this result and no
band label, and the only other consumer in the tree, `scripts/shortlist_in_real_unknown.py`, already
keys off the drop band at 0.2 and cites the 105 of 105 this section re-measures.

## The union axis, registered before it was fitted: the separation survives, the level does not (2026-09-22, fifth)

> **Annotation, 2026-09-22 (sixth).** Nothing here changes. The 9.70% this section bands and the
> refusal of the rest are now the honest genome-wide coverage in a stronger sense than when they
> were written: the prevalence term the roadmap held open as the repair for the other 90% has
> been registered and fitted, and it does not transport. See the sixth section below.

The section above measured that the union of the two axes bands 128 of 128 on-gate pairs and
**refused to adopt it**, because an axis chosen after seeing it win is priced by the search that
found it. This lane adopts it the right way round: registered first, in
`target_calibration.PREREGISTERED_UNION`, before a single rate was computed; then fitted.
`scripts/union_axis.py`, result `union_axis`, **0 AlphaGenome requests**, nothing fetched.

**What the registration had to do that the earlier ones did not.** The axis was chosen after its
result was seen, and no amount of care undoes that, so the registration separates three questions
the 128 of 128 runs together and gives each the evidence it can actually have.

| | what it asks | evidence available | status |
|---|---|---|---|
| **level** | is the union stratum inside 0.9–1 | only the 285 pairs that chose the axis | registered as a **description**, not a test |
| **search** | could two empty axes hand the search a stratum this clean | the same 285, permuted | a real test of a null already known false |
| **separation** | does either arm carry information where nobody looked | off-gate pairs; GTEx eQTLs | the only real test, and it is the one that ran |

The populations were counted before any rate, by functions that never touch `regulated`: on the
gate 70 / 35 / 23 / 157 over both signals, drop only, track only and neither; off the gate through
the sweep's own window 19 / 14 / 1,419 / 8,774; genome-wide, of 593,765 banded `any gene` targets,
8,370 / 25,788 / 23,468 / 536,139, so the union is **57,626 (9.70%)** and neither is 536,139
(90.30%). 72.6% of the on-gate evidence for this axis is one screen (Gasperini2019, 207 of 285).

**The search is not what produced it.** Permuting the labels with the crossing held at its sizes,
over all 14 groupings of the four cells the search really ran over: **1 of 20,000 permutations**
reached a Wilson lower bound of 0.9, p = 0.0001, and the best grouping in the real data is the
union itself. Registered in advance as nearly worthless, and it is: the null it rejects — that
neither axis carries information — was already known false on 2026-09-22.

**The separation is real, and the registration predicted the wrong way.** The registered prediction
was that the off-gate track arm would **not** separate, on the argument that "named on K562" is a
marker of what a K562 screen could test rather than a mechanism. It separated, on 1,419 pairs that
took no part in choosing the axis:

| off the gate, through the sweep's own window | pairs | regulated | observed | 95% interval |
|---|---|---|---|---|
| track only | 1,419 | 102 | **0.0719** | 0.0596–0.0865 |
| neither | 8,774 | 195 | 0.0222 | 0.0193–0.0255 |
| union (either signal) | 1,452 | 130 | **0.0895** | 0.0759–0.1053 |

The intervals clear each other by a wide margin — ×3.2 for the track arm alone, ×4.0 for the union.
**Recorded as registered: wrong on the prediction, and the arm is not a testability marker.** The
independent assay agrees, weakly: on the 3,331 banded targets whose elements carry a distilled GTEx
cis-eQTL, the element's predicted target is among its eGenes for **155 of 261** union targets
(0.5939, 0.5333–0.6517) against 1,554 of 3,070 neither (0.5062, 0.4885–0.5239). The intervals clear
by 0.0094 — a ×1.17 effect on a pre-selected index of elements, which is support for a direction and
nothing more. The third arm was registered in advance as unusable and is: the five non-K562 held-out
cell lines contribute **8 scored on-gate pairs** out of 2,063, which fail the gate rather than the
cell.

**And the same test destroys the level.** The union stratum reads **1.0000 on the gate and 0.0895
off it** — one axis, one definition, two populations, a factor of **11**. The 0.9–1 band is
therefore not a property of the union axis. It is a property of *union ∧ on the gate ∧ K562 ∧
tested by a screen*, which is the third time this project has found the same thing on a third axis,
after the ×6.85 and the ×7. **No population on this disk can confirm that the union stratum sits at
0.9–1 genome-wide**, and the one independent population that can be read on the same axis says it
does not sit there.

**The sharpest thing against the banding below, stated because it is ours.** The `track only` cell
carries 23 pairs, and under the registered rule 23 pairs **buy no band at all** — its own interval,
0.8569–1.0, spans two. The union stratum clears the rule only by pooling those 23 with the 105 the
drop axis had already banded. So the 23,468 genome-wide targets the union adds over the drop axis
rest on evidence the rule refuses when it is read alone. The increment is not priced by the pairs
that justify the increment, and the next registration should price it that way.

**The verdict, applied as registered rather than reinterpreted.** (a) the search prices below 0.01:
**yes**. (b) a population that took no part in the choice separates: **yes**, on both the off-gate
pairs and the eQTLs. (c) the `neither` stratum is one population: **no** — its drop sub-strata read
0.5714 on 42 pairs, 0.5889 on 90 and 0.7600 on 25, and none of them carries a band of its own. So
the registered fallback fires: **the union stratum is banded and `neither` is banded 'not calibrated
here'.** One honest qualification on clause (c): as written it cannot tell *these sub-populations
differ* from *these sub-populations are too small to say*, and here it is mostly the second. The
outcome is the conservative one either way, but the clause is a defect in this registration and is
reported as one rather than claimed as a heterogeneity finding.

**The banding.** `any gene`, of 593,765 targets: **57,626 (9.70%) banded 0.9–1**, 536,139 (90.30%)
`not calibrated here`. Moves: 536,139 lost their band, 37,333 unchanged, **20,293 up, none down**.
`coding gene`, of 440,377: 44,775 (10.17%) at 0.9–1, 395,602 (89.83%) refused; 26,375 unchanged,
18,400 up, none down. Against the two readings already published:

| reading | banded 0.9–1 | refused | of the 55,957 published at 0.9–1 |
|---|---|---|---|
| published 2026-09-17 | 55,957 | 0 | — |
| drop axis (2026-09-22, fourth) | 34,158 (5.8%) | 559,607 | 34,096 keep it, 21,861 uncalibrated |
| track axis (the falsifier that fired) | 31,838 | 0 | 44,372 fall two bands to 0.5–0.75 |
| **union axis (here)** | **57,626 (9.70%)** | **536,139** | **37,333 keep it, 18,624 (33.3%) uncalibrated** |

The union's gain over the drop axis is exactly the `track only` cell, 23,468 targets, and its whole
disagreement with the track axis is that it refuses to price the 536,139 the track axis banded at
0.5–0.75.

**Coverage, and which failure mode this is nearer.** The question the last two sections leave is
whether a method refuses too much or aggregates too much. The drop axis banded 5.8% and refused
94.2%; the track axis banded 100% by giving one number to populations reading 0.58, 0.60 and 0.87.
**The union axis is squarely in the drop axis's mode**: it bands 9.70% and refuses 90.30%. It buys
3.9 percentage points of the genome over the drop axis and declines to say anything about the rest,
which is the right failure to have but is a small return for an axis that looked, on the 285 pairs
that chose it, like the one that would band everything.

**What would settle the level**, none of it reachable at 0 requests: a CRISPRi screen outside the
EngreitzLab benchmark with on-gate union pairs; the benchmark's other five cell lines rescored so
their pairs reach the gate at all (8 of 2,063 do now); or the prospective form — publish the union
stratum's top targets as a list and have them tested.

**What is not touched.** The published 2026-09-17 table, the drop-axis re-banding and the track-axis
re-banding all keep their values and keys; the union is a third axis beside them under its own
result, so all four readings can be held against each other from disk. `CONFIDENCE_BANDS`, the
fitted weights and `TARGET_FEATURES` do not change. No web UI or CLI hunk: `server.py`,
`index.html` and `cli.py` hold no reference to this result and no band label. One consumer is worth
naming — `scripts/shortlist_in_real_unknown.py` keys off a drop of 0.2, and its 24,114 coding
targets are exactly this census's coding `drop only` + `both` (17,521 + 6,593), so its shortlist is
the union's magnitude arm and nothing else. Two corrections to the section above, both reported
rather than quietly fixed: its windowed off-gate table was **training-only**, because the driver
called `add_features` on one arm, which is why the off-gate population is 10,226 here and 8,549
there; and its track axis returned `another track` for every off-gate pair, so `pair_track` now
reads the track at the element and the axis means one thing on both populations. 21 tests.

## The prevalence term, registered before it was fitted: the one constant was a mean over shifts that run both ways, and the genome's band table is one screen's (2026-09-22, sixth)

**The correction cannot be applied where it is needed, and that is the result.** The repair the
roadmap has carried since 2026-09-17 — re-fit the intercept per screen with that screen's own base
rate as an offset — was registered in `PREREGISTERED_PREVALENCE` (commit 1a1a852) after the screen
census was counted and before a single offset was fitted. The census alone settles the scope: **not
one screen the calibration is read on is a screen it was fitted on.** The intersection of the two
screen sets is empty. A per-screen intercept is therefore a thing that exists after a screen has been
run and never before it, so it says nothing about a target in no screen at all, which is the whole
genome. What the fit then found is worse for the published table than the failure it was meant to
repair: the single log-odds shift of **+0.679** that the 2026-09-17 section reported as the
description of its failure is a mean over per-screen shifts running from **−1.303 to +4.228**, and
its sign is wrong for the largest held-out screen, 1,084 of the 1,715 pairs. And the genome-wide band
table that prices 593,765 targets is, to within **0.4%**, the band table of one screen.
`scripts/target_prevalence.py`, result `target_prevalence`, 55 s, **0 AlphaGenome requests**, nothing
fetched.

**The census, counted before anything was fitted.** On the scored population — covered, on a deleted
element, every feature present — the two tables are made of eight screens and share none of them.

| | screen | pairs | regulated | rate | on the gate | on-gate rate | own term |
|---|---|---|---|---|---|---|---|
| fitted | Gasperini2019 | 4,776 | 335 | 7.01% | **207** | 0.7874 | yes |
| fitted | Nasser2021 | 2,931 | 80 | 2.73% | 28 | 0.6786 | yes |
| fitted | Schraivogel2020 | 1,089 | 22 | 2.02% | 10 | 0.6000 | yes |
| read | K562_DC_TAP | 1,084 | 11 | 1.01% | 5 | 0.4000 | yes |
| read | Xie | 416 | 41 | 9.86% | 19 | 0.9474 | yes |
| read | Morris | 177 | 34 | 19.21% | 8 | 1.0000 | yes |
| read | Klann | 32 | 20 | 62.50% | 7 | 1.0000 | pooled |
| read | Reilly | 6 | 6 | 100% | 1 | 1.0000 | pooled |

The three fitted screens pool to 8,796 at 4.97% and the five read ones to 1,715 at 6.53%, which is
the ×1.31 the failure was named for. **The spread inside each table is larger than the gap between
them**: 3.5× across the fitted screens and 99× across the read ones. So "the held-out screens have a
higher prevalence" was never one fact; it is a pooled summary of five populations that have almost
nothing to do with each other. The registered size rule — at least 100 scored pairs and at least 10
regulated — admits six screens and pools Klann and Reilly into 38 pairs at 68.42%.

**On the gate, which is the population whose curve bands the genome, per-screen terms are
unfittable.** There are 285 on-gate pairs. Gasperini2019 is **207 of them (72.6%)** and no other
screen reaches 30. Whatever size rule is written, one screen carries the population, so the
correction the roadmap asked for degenerates there to the fit that is already shipped.

**The same held-out pairs, before and after.** 1,715 K562 pairs, 112 regulated, observed 6.53%. The
offset model is fitted on the training screens with their offsets and read on the held-out screens
with each held-out screen's own base rate, so like the published pooled shift it is **fitted on the
labels it is read against and is a description, not a test.**

| on held-out K562 | bins inside | ECE | Brier | AUPRC | mean predicted |
|---|---|---|---|---|---|
| the shipped curve | **6/10** | 0.0231 | 0.04184 | 0.559 | 0.0441 |
| one pooled shift, +0.6793 (the published description) | **9/10** | 0.0100 | 0.04065 | 0.559 | 0.0653 |
| per-screen offsets | **8/10** | 0.0158 | **0.03082** | **0.737** | 0.0551 |

The registration expected the per-screen offsets to reach 9 or 10 of 10 and called that a description
that could not fail. **It reached 8, below the one constant it was supposed to beat**, and that is the
first thing the split explains. Per screen:

| held-out screen | pairs | observed | shipped mean | shift it needs | offset mean | shift it still needs |
|---|---|---|---|---|---|---|
| K562_DC_TAP | 1,084 | 0.0101 | 0.0269 | **−1.303** | 0.0085 | +0.262 |
| Xie | 416 | 0.0986 | 0.0611 | +0.945 | 0.0877 | +0.216 |
| Morris | 177 | 0.1921 | 0.0701 | +2.050 | 0.1466 | +0.516 |
| Klann and Reilly, pooled | 38 | 0.6842 | 0.2271 | **+4.228** | 0.6038 | +0.524 |

**The published +0.679 is a mean over shifts of −1.303, +0.945, +2.050 and +4.228.** The screen that
carries 63% of the held-out pairs needs the correction in the *opposite* direction, so applying the
published constant to K562_DC_TAP makes it worse, and the pooled table only looks repaired because
the 38 pairs needing +4.2 are averaged against the 1,084 needing −1.3. The 2026-09-17 diagnosis that
"the held-out screens call 6.53% against 4.97%" is arithmetically true and mechanistically empty: the
level does not move between the two tables, it moves between every screen and every other screen,
inside each table as much as across them. Per-screen ECE falls on all four strata — 0.0183 → 0.0048,
0.0442 → 0.0342, 0.1220 → 0.0455, 0.4571 → 0.0971 — while the pooled bin count drops from the
constant's 9 to 8, because the constant is tuned on exactly the pooled table the bins are read on.

**The test that could fail: leave one screen out.** The shape is fitted on the other five screens with
their offsets in place; the held-back screen contributes nothing of its own but its base rate. The
comparator is the identical split with no offset anywhere, through the same solver.

| screen held back | own rate | rest's rate | pooled: bins, ECE, AUPRC | offset: bins, ECE, AUPRC | residual shift |
|---|---|---|---|---|---|
| Gasperini2019 | 7.01% | 3.73% | 9/10, 0.0096, 0.717 | 7/10, 0.0154, 0.715 | −0.442 |
| Nasser2021 | 2.73% | 6.19% | 5/10, 0.0190, 0.410 | **9/10**, 0.0094, 0.410 | +0.173 |
| Schraivogel2020 | 2.02% | 5.59% | 9/10, 0.0147, 0.397 | 9/10, 0.0141, 0.419 | +0.521 |
| K562_DC_TAP | 1.01% | 5.71% | 9/10, 0.0223, 0.268 | 7/10, **0.0049**, 0.264 | +0.205 |
| Xie | 9.86% | 5.03% | 7/10, 0.0422, 0.725 | **9/10**, 0.0339, 0.720 | +0.159 |
| Morris | 19.21% | 4.98% | 0/3, 0.1204, 0.727 | **2/3**, **0.0507**, 0.733 | +0.475 |

**Verdict: FAILED, on three of the four registered clauses.** Bins improve on **3 of 6** screens
against the registered 4. The pairs-weighted expected calibration error falls 0.01724 → 0.01383, a
fall of **0.00341** against the registered 0.005. AUPRC falls on three screens — Gasperini 0.7174 →
0.7152, K562_DC_TAP 0.2678 → 0.2640, Xie 0.7245 → 0.7204 — and the flattening clause is absolute, so
it fires. The one clause that passed is the third and it is the informative one: the median absolute
residual shift after the offset is **0.3235**, under the registered bar of 0.34 and less than half of
the +0.6793 the single pooled correction needed. **So a screen's own base rate accounts for rather
more than half of the level error and leaves the rest**, which is a real finding and not the repair
that was asked for.

**Why it fails where it fails, which was not registered.** ECE improves on five of the six screens,
sometimes by a factor of four; the bin count does not follow, and the screen that loses most is
Gasperini2019. Its own rate is 7.01% and the other five screens pool to 3.73%, so the offset pushes
its mean predicted from 0.0649 up to 0.0856 against an observed 0.0701 and it now **over**-predicts,
residual −0.442. The mechanism is that Gasperini is 54% of the pooled pairs, so "the rest" is a
different population from the one whose intercept the shape absorbed, and the offset and the fitted
intercept double-count. A correction that needs the held-back screen not to dominate the pool cannot
be applied to a population that is 72.6% one screen, which is the on-gate population.

**What it does to the genome-wide band, which is the first paragraph's claim.** The same 593,765
banded targets, walked once, banded under each screen's own on-gate base rate as an offset. The
`published` row is the same walk at a shift of zero, so the table is its own control and reproduces
the 2026-09-17 table to the digit.

| offset taken from | 0.02–0.05 | 0.05–0.1 | 0.1–0.25 | 0.25–0.5 | 0.5–0.75 | 0.75–0.9 | **0.9–1** |
|---|---|---|---|---|---|---|---|
| **published** | — | 45 | 3,837 | 178,387 | 293,568 | 61,971 | **55,957** |
| K562_DC_TAP | 3,590 | 23,723 | 339,732 | 147,518 | 36,427 | 15,261 | **27,514** |
| Schraivogel2020 | 712 | 2,990 | 118,803 | 321,280 | 88,191 | 25,973 | **35,816** |
| Nasser2021 | — | 2,764 | 30,714 | 343,756 | 139,609 | 34,999 | **41,923** |
| Reilly | — | 791 | 4,250 | 247,643 | 240,350 | 50,219 | **50,512** |
| Gasperini2019 | — | 32 | 3,843 | 176,025 | 295,305 | 62,394 | **56,166** |
| Xie | — | — | — | 3,758 | 136,452 | 311,554 | **142,001** |
| Klann | — | — | — | 3,415 | 75,423 | 340,928 | **173,999** |
| Morris | — | — | — | 3,138 | 45,702 | 346,450 | **198,475** |

**The top band runs from 27,514 to 198,475, a factor of 7.2**, on the identical 593,765 targets with
the identical weights; for the coding target it is 17,642 to 116,456, a factor of 6.6. Under
K562_DC_TAP's rate 367,045 targets (61.8%) fall below 0.25 and under Morris's rate none do. **And the
published table is Gasperini2019's table**: 55,957 against 56,166 on `any gene`, a difference of
0.37%, and 37,913 against 38,065 on the coding target, 0.40%. That is not a coincidence and it is the
finding — Gasperini2019 is 72.6% of the on-gate pairs the predicted-target curve was fitted on, so the
intercept the genome inherits is that screen's intercept. The band a swept element's target receives
is, to within a fifth of a percent, an answer to the question *what would Gasperini2019 have called
this*, and the other seven screens of the same benchmark in the same cell line would have answered
between half and three and a half times as often.

**What a band means for a target in no screen at all.** Exactly what the registration said before the
fit, and the fit did not move it: **it cannot be quoted.** The offset's input is a screen's own base
rate; a swept target has no screen, and the only substitute is a guess at some future screen's base
rate, which across six K562 screens of one benchmark runs from 1.01% to 100% and moves the top band
by a factor of 7. The term is not transportable even between *tested* pairs, because the two screen
sets do not intersect: Xie's intercept cannot be fitted on Gasperini, Nasser and Schraivogel, only
read off Xie's own labels. **So the prevalence term does not unlock the 593,765 targets the 2026-09-17
table priced. It gives the reason that table was never quotable**, and it leaves standing the two
things that need no intercept and that this project already has: the **shape**, which leave-one-screen-out
leaves almost untouched (per-screen AUPRC 0.264 to 0.733 under the offset, moving by at most
0.022 against the pooled model and by less than 0.005 on four of the six), reported as an ordering
and never as a probability; and the **measured band on a
stratum**, which is the re-banding's answer — a rate read directly on a stratum with its population
and its interval beside it, and a refusal where the stratum is thin. The re-banding's **5.8%** and the
union axis's **9.70%** of the genome remain the honest coverage, and this section removes the last
reason to think the other 90% was waiting on an intercept.

**What is not touched.** Additive throughout. `score()`, `sweep()` and `reband()` keep their behaviour
and their keys, `data/results/target_calibration.json` and `data/results/target_rebanding.json` are
unchanged on disk, and `SWEEP_FEATURES`, `TARGET_FEATURES` and `CONFIDENCE_BANDS` do not move. The
one change inside the shared fit is a backtracking line search in the new offset solver, which is
required because an undamped Newton step diverges on this near-separable problem; with every offset
set to zero it reproduces `crispri.logistic_fit`'s weights to 1e-6, which is a test, so the
comparator differs from the offset model in the offset and in nothing else. 13 tests.

**Limits.** Every base rate used as an offset is measured on the same pairs the reliability is read
on, so E1 is a description and only the leave-one-screen-out arm is a test. Klann and Reilly are 38
pairs pooled and carry one stratum between them. The on-gate rates that drive the genome table are
smoothed as (k + 0.5)/(n + 1) because four of the eight screens are at 100% on the gate, so the
factor of 7.2 is if anything an understatement of the spread. Six screens is a small number to read a
"4 of 6" clause on, and three of them are the screens the shape was fitted on.

**Next, proposed as roadmap rows** (for the roadmap editor to fold; not written into ROADMAP.md):

1. Area I: retire the 593,765-target band table rather than correct it. The prevalence term is now
   measured and it does not transport; the published confidence is one screen's intercept to within
   0.4%. Quote the re-banding's strata and the union axis, and say "no band" for the rest.
2. Area I: a screen's base rate as a declared input. If a band must be quoted for a planned
   experiment, quote it as a function of that experiment's expected base rate with the range shown,
   not as a single number — the machinery for it is `genome_under_each_screen`.
3. Area I: the shape, tested on its own. Leave-one-screen-out moves AUPRC by less than 0.005 in
   either direction, so the ranking is the part that survives; register a ranking claim (top-k
   precision on a screen that took no part in the fit) and stop registering level claims until a
   screen exists on both sides of the split.

## Pre-registration: the shared logistic fit is undamped, and every caller inherits it (2026-09-22, seventh)

Written before the repair is made and before any result is re-derived.

**The defect.** `genomeos/attribution/crispri.py:407`, `logistic_fit`, takes an undamped Newton
step: it forms the gradient and the ridge-penalised Hessian, calls `solve`, and adds the whole step
to the weights, twenty-five times or until the step is small. Newton's method for logistic
regression is not globally convergent, and on a design where the classes nearly separate the step
can carry the iterate so far that every linear predictor saturates. The sigmoid clamps its argument
to [-30, 30], so a saturated row contributes a curvature of about 9.4e-14 and, if its label is on
the wrong side, a gradient of order one: the next solve divides a gradient of order one by a
curvature of order 1e-13 and the weights leave for 1e12. The iterate never comes back, the step
never falls under 1e-6, and the function returns the twenty-fifth divergent point without a word.

**Why it is dangerous rather than merely wrong.** A diverged fit still *orders* the rows. On the
eight-row design in `tests/test_crispri.py` added with this repair the diverged weights reach
2.25e12 and give seven of eight rows a predicted probability of exactly 0 — including three of the
four positives — and its AUROC is still 1.0, identical to the converged fit's. Every downstream
figure in this project that is read off `logistic_fit` as a *ranking* (AUPRC, AUROC, a gain
interval, a top-k precision) is therefore blind to the failure, and every figure read off it as a
*level* (a reliability bin, a calibration curve, a confidence band) would be silently ruined. The
existing test, `test_logistic_fit_orders_a_separable_feature`, asserts the ordering and the sign of
one weight, which are exactly the two properties divergence preserves; it passes on a diverged fit.

**The inventory of callers**, each with what it fits and what a divergence there would look like:

| Caller | What it fits | Exposure |
|---|---|---|
| `crispri.py:611` | leave-one-chromosome-out folds of `FEATURES` on the K562 training pairs, 23 folds x 3 models | read as AUPRC/AUROC only: a divergence would be **invisible**, the ranking survives |
| `crispri.py:617` | the three `FEATURES` models on all 9,237 covered training pairs, applied to the held-out pairs | read as AUPRC/AUROC and a bootstrap gain interval: **invisible** |
| `crispri.py:843` | leave-one-chromosome-out folds of `CONTACT_FEATURES` (seven models incl. measured Hi-C contact) | read as AUPRC/AUROC: **invisible** |
| `crispri.py:849` | the seven `CONTACT_FEATURES` models on the training pairs with a measured contact | read as AUPRC/AUROC: **invisible** |
| `target_calibration.py:658` (`fit`) | every calibration weight vector in the project: `sweep`, `drop`, `activity`, `target`, and through `scripts/target_calibration.py`, `scripts/target_calibration_gate.py`, `scripts/target_rebanding.py`, `scripts/union_axis.py` | read as a **probability**: a divergence gives predicted 0 or 1 and a reliability table of empty bins, which would look like a legitimate "the model is confident and wrong" result |
| `target_calibration.py:1028` | leave-one-chromosome-out folds of `sweep` and `activity` inside `calibration_result` | read as a probability through `predict`: **plausible-looking**, as above |
| `tests/test_crispri.py:75` | a one-feature separable toy | asserts only the ordering, so it **cannot see** the failure |

`target_calibration.logistic_fit_offset:2089` is the same solver with a per-row offset and already
carries the backtracking line search; it is not exposed. `motif_grammar.fit` and `motif_transfer.fit`
are a different, linear, ridge solve and are not callers.

**The repair, registered.** Move the line search into the shared function: `crispri.logistic_fit`
gains `damped_step`, which halves the Newton step until it stops lowering the penalised
log-likelihood and returns a zero step if twenty halvings do not, which ends the iteration through
the caller's own convergence check. The edit is purely additive — the existing lines of
`logistic_fit` keep their text and one line is inserted between the solve and the update — because
`scripts/check_staged.py` refuses a commit that deletes lines a commit of the last two days added,
and `target_calibration.py`'s copy of the same rule was added today by e4914ef. That copy is
therefore **left in place, knowingly**, and pinned instead: a test asserts that
`crispri.damped_step` and `target_calibration.damped` return the same step on the same input, so
the duplication cannot drift into two rules. Removing it is a roadmap row, not this change.

**What must hold.** The damped fit must reproduce the undamped fit wherever the undamped fit
converges, to **1e-6 in every weight**, which is the tolerance the offset solver was already held
to. The reason it holds is that the line search accepts the full step whenever the full step does
not lower the penalised log-likelihood, which is every iteration of a converging run, so the
iterates are identical and not merely close.

**The results that are re-derived, and that must not move.** Each is compared against the **file
committed in this repository today**, not against a figure quoted in this document, because several
of these numbers moved this morning for unrelated reasons:

1. `data/results/crispri_benchmark.json` — the pre-registered CRISPRi verdict, the held-out AUPRC
   per cell type and the deletion gain intervals (`scripts/crispri_score.py`).
2. `data/results/crispri_contact.json` — the measured-contact comparison (`scripts/crispri_contact.py`,
   offline: every contact is already in `data/knowledge/hic_contact`).
3. `data/results/target_calibration.json` — the reliability table, the band assignment and the
   sweep summary (`scripts/target_calibration.py`).
4. `data/results/target_calibration_gate.json` — the gate arm.
5. `data/results/target_rebanding.json` — the re-banded curve and its strata.
6. `data/results/union_axis.json` — the union axis's 9.70%.
7. `data/results/target_prevalence.json` — the leave-one-screen-out arm, whose comparator is
   `fit` and whose offset arm is already damped.

**What happens if one of them moves.** It is reported first, as the finding, before anything else
in this section, and the old value is **withdrawn**, not reconciled. A number that changes when the
solver stops diverging was computed at a diverged point and was wrong when it was published; the
size of the change is not an argument about whether it mattered, because a diverged fit has no
error bar to be inside. The falsifier is stated here so it cannot be softened afterwards: if any
figure in the seven files above differs from its committed value by more than the rounding the file
itself applies, this section's title becomes that finding, the affected claim is marked withdrawn
where it was made, and the re-derived value is published with the date of the repair beside it.

### The outcome: nothing moved, and the reason is that the real fits were never near the edge

Written after the repair, against the registration above.

**No published number moved.** All seven files were re-derived on the repaired solver and compared
against the copy committed in this repository:

| File | After the repair |
|---|---|
| `crispri_benchmark.json` | byte-identical (verdict `passed`, held-out K562 AUPRC 0.6909 against 0.5501, deletion gain 0.1407 [0.082, 0.2313]) |
| `crispri_contact.json` | every number identical; the `date` stamp alone changed, and `reads: {}` — no matrix was opened, nothing was fetched |
| `target_calibration.json` | every number identical; the `date` stamp alone changed. The top band reads 55,957 of 593,765 as before |
| `target_calibration_gate.json` | byte-identical |
| `target_rebanding.json` | 550 insertions, **0 deletions**: the `by_stratum` counts that `target_calibration.py:1648` began emitting earlier today, on a file last saved before that. No existing value altered |
| `union_axis.json` | byte-identical (44,775 coding targets on the union, 395,602 losing their band) |
| `target_prevalence.json` | byte-identical (3 of 6 screens improving, weighted ECE 0.01724 to 0.01383, median residual 0.3235, verdict `failed`; the band table 27,514 to 198,475) |

**Why, measured rather than assumed.** Before the repair, every one of the 72 fits inside
`crispri.score` — the three `FEATURES` models on the full training set and on each of 23
leave-one-chromosome-out folds — was computed twice, once undamped and once with the line search.
The largest disagreement in any weight over all 72 was **3.4e-7**, and the largest weight reached
was 27.5. The real designs are not near the edge: 9,237 pairs with 451 positives, no separating
feature, and a Hessian that stays well conditioned. The defect was live and inherited by seven call
sites; it had simply never fired on the data this project has.

**How near the edge is near enough.** On a random family of ill-conditioned designs (8 to 80 rows,
1 to 5 columns, column scales drawn over five orders of magnitude, labels separable, near-separable
or random), the undamped solver returns a point **worse than the all-zero start it began from** in
**64 of 3,738 fits, 1.7%**, the worst of them at |w| = 5.0e12. The repaired solver loses on none.
So the exposure was real and roughly one fit in sixty of that family — it is the conditioning of
this project's particular feature matrices, not the solver, that kept the published numbers safe.

**What the repair leaves.** `crispri.logistic_fit` now damps, so all seven call sites inherit the
line search. `target_calibration.damped` and `penalised_log_likelihood` are still a second copy of
the same rule, kept because `check_staged.py` refuses to delete lines committed the same day, and
pinned by `test_the_two_copies_of_the_line_search_are_one_rule`, which asserts they return the same
step and the same objective on the same input. Merging them is a roadmap row for tomorrow, when the
window has passed; the pinning test is then the regression test for that merge.

**The lesson, which is about the test and not the solver.** `test_logistic_fit_orders_a_separable_feature`
had guarded this function since it was written, and it asserts the ordering of the scores and the
sign of one weight. Those are exactly the two properties a diverged fit keeps. On the design now
pinned beside it, the old solver returns weights of 2.25e12, gives seven of eight rows a predicted
probability of exactly 0 including three of the four positives, and scores **AUROC 1.0** — the same
AUROC the correct fit gives. A test that reads a model the way its headline figure reads it cannot
see a failure that the headline figure is blind to; the two tests added here assert instead that
the weights stay finite, that the probabilities separate, and that the solver never returns a point
whose penalised objective is below the start's, which is the invariant divergence actually breaks.

## Pre-registration: the unknown scoring off the one-gene table, and what a named target means against matched random windows (2026-09-27)

`attribution/unknown_scoring.py` is the module behind milestone 1.3's clause "the
constrained-unknown blocks attributed to a gene and a tissue". It reads the compact
one-gene-per-element table at three sites, and the two lanes of 2026-09-22 that moved
`crispri.py` and `eqtl.py` onto the sweep's own per-element response cache left it on the table.
This registers what the three sites do, what the cache can and cannot add, and the one comparison
that decides whether reading the whole window is an improvement at all. **0 AlphaGenome requests:
every number below and every number in the result section is on disk.**

### The three sites, and whether a silence is read as "no effect"

| site | what it reads | what a silence there becomes | is it a zero? |
|---|---|---|---|
| `unknown_scoring.py:65` (`scored_elements`, consumed by `annotate:111`) | `run_elements("enhancer_targets_all", chrom)` | `predicted` null → `abs_log2`, `silencer`, `cells_acting`, `target_tss_distance` all `None`; `moves_gene` and `names_coding` `False` | **no, a `None` that `stratified:242` then drops** — so the four magnitude metrics are silently *conditioned on the gate* |
| `unknown_scoring.py:497` (`candidates_on_scored_chromosomes`) | the same table, per candidate block | no element clears the gate → `gain: "nothing"`, `target: null`, `log2: null` | **no, but it reads as one**: "nothing" also means "no registry element lies here", and the two are not the same silence |
| `unknown_scoring.py:681` (`coverage`) | the same table plus the two samples' tables for `before_ids` | a gated element cannot make a block "named"; `cells_acting` `None` is falsy | no |
| `unknown_scoring.py:657` (`against_vista`, reached from site 1's rows) | `max((r["abs_log2"] or 0.0 for r in over), default=None)` | a VISTA element covered only by gated elements gets **`max_abs_log2` = 0.0** and is averaged into `mean_max_abs_log2_positive` / `_negative` | **yes. This is the zero.** |
| `unknown_scoring.py:588` (`against_mpra`) | `"abs_log2": r["abs_log2"] or 0.0` | the same zero is written into each pair | yes, but the field is never consumed: the correlation reads `predicted_by_cell`, which is `None` for a gated element, so those elements are dropped instead |

So the answer to the question that decides everything is: **one zero, at `unknown_scoring.py:657`,
in the module's only enhancer-level measured validation**; one written-and-unused zero at
`unknown_scoring.py:588`; and everywhere else a `None` that is dropped rather than zeroed — which
is not the CRISPRi defect but is the eQTL defect, a comparison conditioned on the gate and
published as a comparison.

### How much of the window the module reads, counted first

Over chr21's 12,139 scored elements (the one chromosome this result is computed on), against
`data/knowledge/alphagenome/elements/chr21.json.gz` read through `ElementResponses`
(0.7 s, one chromosome held at a time, 0 requests):

| | count | share |
|---|---|---|
| cached element-gene rows available | 347,614 (28.6 genes per window, median 28) | |
| gene readings the module takes (`predicted` and `predicted_coding`) | 13,020 | **3.75%** |
| per-cell readings available (4 lines × every gene) | 1,389,752 | |
| per-cell readings the module takes (4 lines × the top gene of a mover) | 31,384 | **2.26%** |
| elements the cache answers a magnitude for | 12,139 | **100%** |
| elements whose `abs_log2` is `None` today, dropped from every magnitude comparison | 4,293 | **35.4%** |
| of those, the cache answers | 4,293 | **100%** |
| elements not in the cache | 0 | |

The gate's rate is not the same on the two arms of the module's own comparison: **49.5% of the
2,334 elements over an unknown block clear it against 68.2% of the other 9,805**, an 18.5-point
difference. That is what makes the four gated metrics a conditioned comparison rather than a
comparison: the arms are selected differently by the very quantity the table gates on.

### What the cache cannot add, which is the more important half

**The compact table is not censoring the target.** `predict.enhancer_target.score_element`
computes `predicted_coding` by running `predict_target` over *every coding gene in the window* and
keeps its head, so the table's one coding gene was already chosen from the whole window. Counted:
**5,174 of 12,139 elements name a coding gene from the compact table, and 5,174 of 12,139 name one
from the cache gated at the same `MIN_EFFECT` = 0.1 — the same number, and the same elements.**
Reading the window cannot name one additional target at the threshold. Ungated it would name
11,372 (93.7%), which is not a gain but the removal of the only thing keeping the rate from being
"almost every element names something".

So this is registered as a check rather than a prediction: **any change at all in `names_coding`,
`moves_gene`, or the block counts in `coverage` would mean the compact table was censoring after
all**, and the re-run must reproduce them to the element.

What the cache can add is three things, each of which is either a de-conditioning or a silence
given its name: the head magnitude for all 12,139 elements rather than the 7,846 that clear the
gate; the per-cell question asked of every gene in the window rather than of the top gene only;
and the 4,293 gated elements' measured heads in place of the zeros at `against_vista`.

### The clause that decides it: what a named target means at matched random windows

`constrained_unknown_targets.json` carries the sentence that sets 1.3's standing — "a named target
here is a lead and not a finding: the same model names a target at 87% of matched random windows".
That 87% is borrowed from the locus benchmark (303d0b0, ef80075: some derived layer names a target
at 52 of 60 matched random windows), which is a different instrument — it pools every layer over a
whole window, this is one element's own deletion — and this module has never had a matched-random
control of its own. One is added here, on this instrument, with 0 requests: for each of chr21's 446
UNKNOWN blocks, `RANDOM_DRAWS` = 50 windows of the same length placed uniformly at random inside
the span of the chromosome's scored elements and rejected if they overlap any UNKNOWN block, with
the same block-level question asked of each. 444 of the 446 blocks admit such a window; the two
that do not are 5.01 Mb and 2.35 Mb long.

Measured during this census, on the compact table, so reported here as a census figure and not as a
prediction test:

| block-level question, of blocks or windows carrying at least one scored element | UNKNOWN blocks (285 of 446) | matched random windows (16,706 of 22,200) |
|---|---|---|
| carries an element naming a coding gene | 140, **49.1%** | 10,908, **65.3%** |
| carries an element naming a coding gene with a cell line at the bar | 89, **31.2%** | 8,005, **47.9%** |

**The unknown blocks name a target less often than length-matched random windows on the same
chromosome, by 16.2 points.** That is the honest reading of 1.3's clause on this instrument, and it
is worse than the borrowed 87% suggested, not better: it is not "naming is easy everywhere", it is
"naming is harder here than at random".

**This is why the direction of the change is registered as it is.** Because the target set is
identical before and after, no reading of the window can move the 49.1% or the 65.3%. The one
block-level figure the window can move is the cell line, and if it moves the matched-random rate by
as much or more than it moves the blocks, then reading the window has inflated a count without
improving what can be said — and that is to be reported as an inflation, in those words, not as a
gain.

### The directions registered, before the module is changed

Nothing below has been computed.

1. **`head_abs_log2` over all 12,139 elements** (new, beside the untouched `abs_log2`): the
   unknown-minus-rest difference is expected to be **larger in magnitude** than the gated
   −0.0294 (p 0.048), because the gate removes the sub-0.1 tail from 31.8% of the rest arm and
   50.5% of the unknown arm, so ungating pulls the unknown arm down further. Falsifier: a
   difference that narrows or flips sign, which would mean the gated magnitude gap was itself an
   artefact of the gate in the opposite direction.
2. **`cells_acting_window`** (any gene in the window at the bar on a cell's own track, against the
   top gene only): the raw share of elements acting in at least one line goes 3,260 of 7,846
   movers to 3,693 of all 12,139, and per element 15.8% → 17.9% over unknown blocks against
   29.5% → 33.4% elsewhere. The matched difference and its permutation p are not computed and are
   registered as unknown. The relative change is **larger on the rest arm than on the unknown
   arm**, which already says the window's extra cells are not a property of unknown sequence.
3. **VISTA**: `mean_max_abs_log2_positive` 0.4259 and `_negative` 0.0953 are computed over
   structural zeros — 4 of the 6 negatives and 3 of the 13 positives are covered only by gated
   elements. Replacing those zeros with the measured heads (median 0.064 on gated elements) must
   raise the negative arm more than the positive arm, so **the contrast is expected to shrink**.
   Falsifier: a contrast that widens.
4. **lentiMPRA**: the correlation is computed on the 329 elements whose top gene cleared the gate;
   the window's strongest per-cell effect exists for every covered element. The ungated rho is
   expected to **fall in magnitude or hold**, because the added elements are all small-effect.
   Falsifier: a rise, which would say the gate was discarding elements whose sub-threshold effects
   still track measured activity — the same finding the eQTL lane reached about `MIN_EFFECT` on
   2026-09-22, and it would be worth more than the fix.
5. **The 69 syntax candidates**: 35 have no scored element inside them at all, 9 have elements of
   which none clears the gate, and 25 have a mover. The current stored table reports only 2 rows
   because it was computed on 2026-09-13 when two chromosomes were complete; 19 are complete now
   and all 69 candidates lie on them. The re-run therefore changes that table for a reason that
   has nothing to do with the cache, so **both re-runs are reported**: the compact table over 69
   candidates, which is the control, and the cache over the same 69.
6. **Everything already published is annotated, never rewritten.** `abs_log2`, `silencer`,
   `cells_acting`, `target_tss_distance`, `moves_gene`, `names_coding` and every field of
   `coverage`, `ablation`, `mpra` and `vista` keep their names and their meanings, and the re-run
   must reproduce them to the digit; the new fields are added beside them. The 2026-09-13 section
   above is annotated in place.

## The result: the window named no target the table had not, and the one count it raised it raised faster at random (2026-09-27, later)

`unknown_scoring.run("chr21")` re-run with the response cache handed in (133.1 s against 32.7 s, the
extra time being 17 chromosome archives decompressed for the syntax candidates, **0 AlphaGenome
requests, nothing fetched**). Every pre-existing field computed from chr21's own tables reproduces
**to the digit** — the ablation, all six stratified metrics, `mpra`, `vista`, every count in
`coverage` — which is the control. Thirteen figures differ and all thirteen are in the two blocks
that read *other* chromosomes' results: `syntax_candidates` (2 complete chromosomes on 2026-09-13,
19 now) and `cost` (chr1 and chr20 finished in between). That confound was registered, and both
readings of the candidates are computed in the same run so the cache's effect can be told from it.

### The registered check: the compact table was not censoring, confirmed

| | compact table | the window, same threshold | disagreements |
|---|---|---|---|
| elements naming a coding gene, chr21 | **5,174** of 12,139 | **5,174** of 12,139 | **0** |

The same elements, not merely the same count. Ungated the window names one at 11,372 (93.7%), which
is the threshold doing the work and not the table. **So reading the whole window cannot name a
single target this module was not already naming**, and the naming rates of the 2026-09-13 section
stand exactly as published.

### The decider: what a named target means against matched random windows

New, and on this instrument rather than borrowed from the locus benchmark: 22,200 windows, 50 per
UNKNOWN block, of that block's length, placed uniformly at random inside the 41.6 Mb span of chr21's
scored elements and rejected if they overlap any UNKNOWN block; 444 of 446 blocks drawable, the two
that are not being 5.01 Mb and 2.35 Mb long; an element belongs to the window its midpoint falls in.
Rates over the 285 blocks and 16,706 windows that carry at least one scored element.

| block-level question | UNKNOWN blocks | matched random windows | difference |
|---|---|---|---|
| names a coding gene — **before and after, unchanged** | 140, **49.12%** | 10,908, **65.29%** | **−16.17 points** |
| a coding gene and a cell line, from the named gene (before) | 89, 31.23% | 8,005, 47.92% | −16.69 points |
| a coding gene and a cell line, from the whole window (after) | 97, **34.04%** | 8,524, **51.02%** | **−16.99 points** |

Two things follow and they are the answer to the task.

**The 1.3 clause is weaker than the borrowed 87% made it sound, not stronger.** On this instrument
the unknown blocks name a coding target at 49.1% where length-matched random windows on the same
chromosome name one at 65.3%. A named target at an unknown block is not merely unremarkable against
chance, it is **16.2 points rarer than at random sequence of the same size** — because the unknown
space is gene-poor by construction, which is the same confound the element-level strata correct for
and which the block-level claim never had a control for at all.

**The one count the window raised, it raised faster at random.** The cell line is the only
block-level figure the window can move, and it moves it: 89 → 97 blocks with a target and a cell by
the midpoint rule, 90 → 99 by `coverage`'s largest-overlap rule (regulatory 81 → 89, neutral 6 → 7,
fossil 3 → 3, constrained_unknown 0 → 0, structural 0 → 0). At the matched random windows the same
change is 8,005 → 8,524, and the gap against the control **widens** from −16.69 to −16.99 points.
That is an inflation and it is named as one: the extra cell lines are a property of reading a 1 Mb
window's worth of genes, not a property of unknown sequence, and any block count quoted off them
must be quoted against the control beside it. The registration said this outcome was to be reported
in those words, and it is.

### What the change did buy: three comparisons that were conditioned on the gate

**The gate selects the two arms unequally**, which is why this matters: 1,155 of the 2,334 elements
over an unknown block clear `MIN_EFFECT` (49.5%) against 6,691 of the other 9,805 (68.2%). Every
metric computed only on the movers was therefore conditioned on a variable that differs between the
arms by 18.5 points.

| per element, matched inside 443 shared strata | conditioned on the gate | over every element | |
|---|---|---|---|
| magnitude, unknown | 0.3091 (n 1,155) | **0.1815** (n 2,334) | |
| magnitude, rest | 0.3384 (n 6,691) | **0.2440** (n 9,805) | |
| difference | **−0.0294**, p 0.048 | **−0.0625**, p 0.0005 | registered to grow: **it did, 2.1×** |
| cell lines of four, unknown | 0.4523 (movers) | 0.2677 (all) | |
| cell lines of four, rest | 0.7962 (movers) | 0.6072 (all) | |
| difference | −0.3439, p 0.0005 | −0.3395, p 0.0005 | barely moves |

The magnitude result is the one that changes standing. Published, the unknown elements' effects were
0.309 against 0.338 matched, a 0.029 gap at p 0.048 — a difference at the edge of its null, quoted
over the 62% of elements the gate had kept. Over every element the sweep answered for, the gap is
0.0625 at p 0.0005. The direction was always right and the size was **less than half of what it
is**, because the gate was removing the sub-0.1 tail from half the unknown arm and a third of the
rest. Raw, the share of elements acting in at least one cell line goes 3,260 of 7,846 movers to
3,693 of all 12,139, and per element 15.8% → 17.9% over unknown blocks against 29.5% → 33.4%
elsewhere — larger on the rest arm, which is the same inflation the control catches at block level.

### The zero, and the measurement it was in

`against_vista` took `max((r["abs_log2"] or 0.0 ...))`, so a VISTA element covered only by elements
the gate had silenced was scored as a predicted change of **exactly zero** — a statement the sweep
never made — and averaged into the module's only enhancer-level measured contrast.

| VISTA on chr21 | published | over what the sweep predicted |
|---|---|---|
| mean max effect, 13 positives | 0.4259 | **0.4393** |
| mean max effect, 6 negatives | 0.0953 | **0.1267** |
| the contrast | 0.3306 | **0.3126** |
| rows that had to be given a zero | 3 positives, 2 negatives | none |

Registered to shrink, and it shrank, by 0.018. **One correction to the registration**: it said four
of the six negatives were covered only by silenced elements. Two were; the other two have no scored
element covering them at all, so `max_abs_log2` was `None` and `_mean` dropped them rather than
zeroing them. The new `positives_scored_as_a_structural_zero` and
`negatives_scored_as_a_structural_zero` fields are what told the two silences apart, which is the
distinction that was missing in the first place. With 19 elements the note above it still holds:
reported, not read.

### The falsifier fired on lentiMPRA, and it is worth more than the fix

The correlation against measured activity was computed on the 329 elements whose named gene cleared
the gate. The window's strongest per-cell effect exists for every covered element, so the same
correlation can be taken over all 445.

| lentiMPRA, chr21 | n | signed effect vs activity | unsigned magnitude vs activity |
|---|---|---|---|
| K562, gated | 329 | 0.2208 (p 0.0005) | 0.1043 (p 0.057) |
| K562, every covered element | **445** | **0.1885** (p 0.0005) | **0.0893** (p 0.060) |
| HepG2, gated | 329 | 0.2391 (p 0.0005) | 0.1699 (p 0.002) |
| HepG2, every covered element | **445** | **0.2095** (p 0.0005) | **0.2081** (p 0.0005) |

Registered: fall in magnitude or hold, with a rise named in advance as the outcome worth more than
the fix. The signed correlation fell on both lines as registered, 0.221 → 0.189 and 0.239 → 0.210,
which is what adding 116 small-effect elements to a signed measure does. **The unsigned magnitude
rose on HepG2, 0.1699 → 0.2081, and its permutation p went 0.002 → 0.0005** while falling slightly
on K562. So on one line of two the 116 elements the gate was discarding carry *more* information
about measured activity per element than the ones it kept. That is the same finding the eQTL lane
reached about `MIN_EFFECT` on 2026-09-22 from a completely different direction — the threshold is a
weak confidence signal and not the line between an answer and a non-answer — and it is the first
time it has been shown against a measurement rather than against a database.

### "Nothing" was two silences, and the syntax candidates say which

All 69 genome-wide syntax candidates now lie on completely scored chromosomes (19 of them; 17 caches
decompressed, chr12's and chr14's candidates needing none). Both readings, in the same run:

| the 69 candidates | compact table | the whole window |
|---|---|---|
| **"nothing"** | **44** | — |
| no registry element lies inside the block: never asked | — | **35** |
| asked, and answered below the bar | — | **9** |
| target and magnitude, no cell line at the bar | 17 | 15 |
| target, cell and magnitude | 8 | **10** |

The 44 "nothing" rows were two different statements. Thirty-five of those blocks hold no ENCODE
registry element at all, so the sweep was never asked about them and no amount of quota would have
changed the answer — a fact about the element registry's coverage of the unknown space, not about
the model. Nine were asked and answered, and the answer is a number: the strongest gene in the
window moves by 0.0204 to 0.0947 (chr1:29,326,813 LINC01756 0.0558; chr1:156,070,514 MEX3A 0.0759;
chr10:128,696,482 ENSG00000289400 0.0204; chr15:96,415,047 ENSG00000242295 0.0789; chr16:46,575,094
ENSG00000260251 0.0551; chr19:53,943,950 CACNG7 0.0628; chr3:171,938,715 TMEM212-AS1 0.0573;
chr6:108,188,805 NR2E1 0.0947; chrX:41,077,750 USP9X 0.0368). Two candidates move from "no cell line
at the bar" to a cell line, both through a gene that is not their element's top target:
chr17:40,365,372 (GJD3, −1.3884) and chrX:134,166,932 (CCDC160, −0.1194), each reaching the bar on
all four lines.

### The standing of milestone 1.3's second clause, stated plainly

Before this the clause's result said "a named target here is a lead and not a finding" and priced
that against 87% of matched random windows from another instrument. It can now be priced on its own:
**140 of the 285 unknown blocks carrying a scored element name a coding gene, 49.1%, against 65.3%
at length-matched random windows on the same chromosome.** Reading the sweep's whole window does not
move either number and cannot, and the one count it does move it moves faster at the control. So the
clause is not 0% defensible for want of a wider reading; it is 0% defensible because the quantity it
reports is below its own chance level, and what would change that is a control the attribution can
pass, not more of the model's output. The three things this lane did deliver are the control itself,
a magnitude comparison no longer conditioned on the gate (−0.0294 at p 0.048 → −0.0625 at p 0.0005),
and a measured contrast with no invented zeros in it.

**What still reads the compact one-gene table**, unchanged and listed rather than silently left:
`candidates.py:734`, `closure.py:69`, `compile.py:112`, `organise.py:91`,
`confidence_calibration.py:244`, `vista.py:278`, `decompile.py:113`, `motif_transfer.py:1025`,
`syntax_tiling.py:252`, `benchmark/loci.py:1137`, `cli.py:3832`, `target_calibration.py:173`,
`crispri_direction.py:118`, and `scripts/constrained_unknown_targets.py:42`, which is the script that
carries this milestone's own result and whose 87% sentence this section prices. On the evidence here,
moving any of them will not change a naming rate: the table's target is the window's head. What it
will change is any figure they compute over movers only, and any zero they write for an element the
gate silenced. `annotate` reads the cache only when it is handed a reader, so
`scripts/confidence_calibration.py` and `scripts/measured_layer.py`, which call it for covariates,
pay nothing and see `NOT_READ` in the window fields — a fourth named silence beside the reader's
three, because "this caller did not look" is not "the sweep did not answer".

## Pre-registration: the CRISPRi result against the published baseline on the same pairs, a second cell type costed, and the gain held against coverage (2026-09-27)

The one result in this project an outside reader could be shown is the CRISPRi one: adding the
AlphaGenome predicted-deletion feature to an activity-over-distance baseline raises held-out K562
AUPRC from 0.550 to 0.691 (gain +0.141, 95% +0.082 to +0.231, 1,744 pairs, 114 regulated). A
sceptical reading of 2026-09-27 named what it lacks: a second cell type with real n, a published
baseline on the same pairs, a coverage-matched arm for the *gain*, and the feature-selection
caveat. This section registers all four before anything new is scored. The registration is also in
code, `crispri.PREREGISTERED_PUBLISHED` and `crispri.PUBLISHED`, committed before the scoring code.
**AlphaGenome requests: 0 prediction requests; one metadata call (`output_metadata`) to read the
model's track list.**

### A second cell type: the model has the tracks, the cache does not, and it costs 705 requests

The held-out refusal said "no AlphaGenome line for this cell type". That is not true of the model.
Its own track list, read on 2026-09-27, carries RNA-seq tracks for **HCT116 (EFO:0002824, 3
tracks), Jurkat (CLO:0007045, 1) and WTC11 (EFO:0009747, 3)**, and DNase for all three. No proxy
is needed and none is proposed. The refusal is a fact about the sweep's per-element cache, which
kept per-gene values for four cell lines only (`enhancer_target.CELLS`).

Whether a cell type could carry a result is settled by its regulated pairs, not by the model:

| held-out cell | pairs | regulated | on a deleted element | regulated there | registry elements to score |
|---|---|---|---|---|---|
| K562 | 1,918 | 118 | 1,744 | 114 | (cached) |
| HCT116 | 396 | 34 | 363 | 34 | **705** |
| WTC11 | 1,921 | 15 | 1,616 | 14 | 614 |
| GM12878 | 68 | 16 | 62 | 14 | (cached) |
| Jurkat | 75 | 7 | 64 | 6 | 111 |

Only HCT116 has enough positives to say anything. The frozen feature is the largest predicted drop
over every registry element overlapping the pair, so scoring HCT116 costs one deletion request per
overlapping registry element: **705**. Deleting the 363 CRISPR elements themselves would cost 363
but changes the feature, so it is not the frozen model. The lane's budget is 20. **HCT116 is
registered here and not run**; the number goes to the coordinator. When it is bought: same code
path with HCT116 added to the cells kept per gene, weights unchanged; it passes if the gain on the
363 covered pairs is above zero, and counts as a replication only if the chromosome-bootstrap
interval excludes zero. A gain at or below zero makes this a K562 result, to be described as one.

### The published baseline: same pairs, same split, figures from the primary source

Gschwind et al., *An encyclopedia of human enhancer–gene regulatory interactions*, Nature 2026
(doi:10.1038/s41586-026-10781-4, PMC13471189), Supplementary Table 3, benchmarks on **the same two
files this module reads**. The training set is 10,356 K562 pairs with 471 positives, scored by
hold-one-chromosome-out for supervised models and unweighted AUPRC. The held-out set is 4,378 pairs
with 190 positives over five cell types, pooled, with AUPRC weighted by each pair's
`direct_vs_indirect_negative` (the 190 positives sum to 157.39, the paper's "157.39 weighted").

| set | predictor | AUPRC [95%] |
|---|---|---|
| training, unweighted | ENCODE-rE2G Extended | 0.737 [0.693, 0.775] |
| | ENCODE-rE2G | 0.662 [0.616, 0.706] |
| | ABC (DNase × H3K27ac, ENCODE Hi-C) | 0.613 [0.558, 0.657] |
| | ABC (DNase, average Hi-C) | 0.565 [0.511, 0.610] |
| | DHS & DNase RPM × distance (norm.) | 0.527 [0.479, 0.575] |
| | distance to TSS | 0.436 [0.387, 0.482] |
| held-out, weighted | ENCODE-rE2G | 0.556 [0.468, 0.631] |
| | ABC (DNase, average Hi-C) | 0.465 [0.378, 0.541] |
| | distance to TSS | 0.363 [0.281, 0.438] |

**The comparison this module has published so far is not comparable to any of these.** It scores
only the pairs on a deleted element (9,237 training, 1,744 held-out K562), unweighted, with
average precision rather than the benchmark's trapezoid, and held-out K562 alone rather than the
pooled five. So the numbers are re-scored on the published pair sets with the published estimator.
The estimator was checked on the one predictor that needs no model, before this registration: raw
distance to TSS gives **0.4359** unweighted on training (published 0.4359) and **0.3631** weighted
on held-out (published 0.3631).

Registered comparisons:

1. **Training, published split.** Hold-one-chromosome-out on all 10,356 pairs, fitted on the
   covered pairs of the other chromosomes (the frozen recipe); an uncovered pair carries deletion
   features of zero, the benchmark's own minimum-score rule. Unweighted, beside the training rows
   above. The pairs and the split are comparable. One thing is not: this module's deletion features
   were chosen after reading the single predictors on these same pairs. ENCODE-rE2G's features were
   also chosen by sequential selection on this set (its Fig. 4c), so both numbers carry a selection
   on the data, and that is stated with them.
2. **Held-out, published pairs and weighting.** The frozen weights applied to all 4,378 pairs,
   pooled, weighted. The handicap is stated before the number: 2,392 of the 4,378 pairs (HCT116,
   Jurkat, WTC11; 56 regulated) carry no deletion value. Bands against ENCODE-rE2G's interval
   [0.468, 0.631]: above it, "higher than the published ENCODE-rE2G figure on the same pairs";
   inside it, "in the range of ENCODE-rE2G"; below it, "below ENCODE-rE2G: the gain is over a
   baseline weaker than the published state of the art". `activity + distance` is placed against
   ABC's [0.378, 0.541] by the same bands, which says whether this project's own baseline is weak.
   The K562-only weighted figure is reported as this module's alone, because the source gives no
   per-cell-type held-out figure.

### The coverage-matched arm, for the gain

Held-out K562; the covered rate is 96.6% for regulated pairs against 90.6% for the rest.

- **Arm 1:** all 1,918 pairs, uncovered ones with deletion features of zero. The gain, with its
  chromosome-bootstrap interval.
- **Arm 2:** the covered pairs, with the regulated arm thinned at random to the non-regulated
  arm's coverage rate, over 1,000 draws. The median gain and the share of draws above zero.
- **Arm 3, the control:** `activity + distance` plus a covered-or-not indicator in place of the
  deletion features, fitted on all training pairs and scored on all 1,918. If coverage alone buys
  the gain, this arm shows it.

The gain is **invariant to coverage** only if arm 1's interval excludes zero, arm 2's share is at
least 0.95, and arm 3's gain is below arm 1's. Otherwise the headline must say "on the pairs the
sweep covered".

### What each outcome allows the result to be called

- Held-out pooled above ENCODE-rE2G's interval: "matches or exceeds the published model on its own
  benchmark". This would be a surprise, since 55% of the pairs carry no deletion.
- Inside the interval: "a sequence-model feature brings a simple baseline into the published
  model's range".
- Below the interval: "improves this project's own baseline; not competitive with the published
  state of the art".
- The coverage bar failing: the claim is restricted to covered pairs whatever the band.
- HCT116 unrun: the result is **one cell type** and is described as such. GM12878 (14 regulated) is
  reported but carries no weight.

## The result: the CRISPRi gain survives coverage and the published pair sets, the held-out comparison is flattered by a filter, and the result is still one cell type (2026-09-27, later)

The registration above was committed first (fe61c36). The code is `crispri.score_published`, the
script `scripts/crispri_published.py` (90 s), and the result `data/results/crispri_published.json`.
**AlphaGenome requests: 0.** The deletion values are the sweep's own cache.

### What goes first: three limits the numbers do not remove

1. **The held-out comparison with ENCODE-rE2G flatters this module, for a reason the registration
   did not foresee.** The held-out file's positives were selected on chromatin at the tested
   element (Gschwind et al. 2026, Methods: "filters based on effect size and chromatin state").
   All 190 positives sit in an `H3K27ac` or `High H3K27ac` element. Of the 4,188 negatives, 1,438
   do not (552 CTCF, 297 H3K27me3, 589 no H3K27ac). This module's activity term is
   sqrt(DNase × H3K27ac), read from the same file, so on the held-out set it partly reads the rule
   that made a pair positive. The published held-out figure is ENCODE-rE2G's DNase-only model, which
   does not read H3K27ac. The registered band (below) came out "above the published interval", and
   **that wording may not be used.** A post hoc diagnostic, not registered and labelled so, refits
   both models with DNase alone.
2. **It is one cell type.** HCT116 is the only other held-out cell type with enough positives, and
   AlphaGenome has its own HCT116 tracks, but the frozen feature costs 705 deletion requests against
   a budget of 20. It was not run. GM12878's 14 regulated covered pairs give a weighted gain of
   +0.015 (95% −0.094 to +0.205), which says nothing either way.
3. **The idea is not new.** The AlphaGenome preprint (Avsec et al. 2025, bioRxiv
   10.1101/2025.06.25.661532, Fig. 4j) already added an AlphaGenome input-gradient score to
   ENCODE-rE2G-extended on this CRISPRi dataset. What is this project's own is the deletion form of
   the feature and a held-out test frozen before scoring. No figure from that preprint is quoted here,
   because none was read from the source.

### 1. Training, the published pair set and split (the cleanest comparison)

All 10,356 pairs, hold-one-chromosome-out, unweighted, scored by the benchmark's estimator:

| model | AUPRC | published, same pairs |
|---|---|---|
| distance (this module) | 0.423 | distance to TSS 0.436 [0.387, 0.482] |
| activity + distance | **0.507** | ABC (DNase, avg Hi-C) 0.565 [0.511, 0.610] → **below** |
| activity + distance + deletion | **0.724** | ENCODE-rE2G 0.662 [0.616, 0.706] → above; ENCODE-rE2G Extended 0.737 [0.693, 0.775] → **inside** |

Deletion gain +0.217 (95% +0.177 to +0.261). **This project's own baseline is weaker than
published ABC.** It sits just under ABC's interval, level with the paper's simple DHS-and-DNase
baselines (0.506 to 0.527). With the deletion, the frozen model reaches the range of the paper's
best model. The caveat stays attached: the deletion features were chosen after reading the single
predictors on these same pairs. ENCODE-rE2G's features were also selected on this set (its
Fig. 4c), so neither figure is free of selection.

### 2. Held-out, the published pairs and weighting

All 4,378 pairs, five cell types pooled, weighted by direct-effect probability, frozen weights:

| model | weighted AUPRC | registered band |
|---|---|---|
| distance | 0.363 (published 0.363) | — |
| activity + distance | 0.567 | above ABC's [0.378, 0.541] |
| activity + distance + deletion | 0.677 | above ENCODE-rE2G's [0.468, 0.631] |

Deletion gain +0.110 (95% +0.061 to +0.174). This is lower than on K562 alone, as registered,
because 2,392 pairs carry no deletion value. Per cell type (weighted; the deletion exists only for
K562 and GM12878):

| cell | pairs | positives (weighted) | activity + distance | + deletion | gain [95%] |
|---|---|---|---|---|---|
| K562 | 1,918 | 118 (100.8) | 0.591 | 0.727 | +0.136 [+0.077, +0.228] |
| GM12878 | 68 | 16 (14.3) | 0.785 | 0.799 | +0.015 [−0.094, +0.205] |
| HCT116 | 396 | 34 (22.6) | 0.495 | 0.495 | no deletion value |
| WTC11 | 1,921 | 15 (13.4) | 0.520 | 0.515 | no deletion value |
| Jurkat | 75 | 7 (6.3) | 0.575 | 0.579 | no deletion value |

The baseline landing *above* ABC here but *below* it on training is the filter at work.

**Post hoc, DNase in place of sqrt(DNase × H3K27ac)** (not registered; found after the numbers
above were read):

| | DNase + distance | + deletion | gain [95%] |
|---|---|---|---|
| training, all pairs, LOCO | 0.496 | 0.720 | +0.224 [+0.182, +0.268] |
| held-out pooled, weighted | 0.476 (ABC 0.465: inside) | 0.639 (ENCODE-rE2G upper bound 0.631) | +0.164 [+0.102, +0.237] |
| held-out K562 covered, headline estimator | 0.502 | 0.681 | +0.179 [+0.102, +0.287] |

With H3K27ac taken out, the baseline falls onto ABC, the deletion's gain gets *larger*, not
smaller, and the combined model lands just above ENCODE-rE2G's published interval. That is an
unpaired comparison against a published interval, not a test, so the defensible phrase is **"in
the range of the published ENCODE-rE2G model on the same pairs"**. It is not "exceeds".

### 3. The coverage-matched arms, held-out K562

| arm | result | bar | met |
|---|---|---|---|
| 1. all 1,918 pairs, 174 uncovered at zero | gain +0.136 [+0.081, +0.225] | interval above zero | yes |
| 2. regulated arm thinned from 96.6% to 90.6% covered (107 of 114 kept), 1,000 draws | median +0.145, range +0.124 to +0.160, 100% above zero | ≥ 95% above zero | yes |
| 3. control: a covered-or-not indicator in place of the deletion | gain +0.001 [−0.001, +0.003] | below arm 1 | yes |

**The gain is invariant to coverage.** Being on a deleted element buys nothing (arm 3). The gain
holds when the arms are matched (arm 2) and when uncovered pairs are counted (arm 1).

### What the result may now be called

- Allowed: "On the ENCODE CRISPRi benchmark, adding an AlphaGenome predicted-deletion feature to an
  activity-and-distance model, frozen before the held-out pairs were scored, raises held-out K562
  AUPRC from 0.55 to 0.69 (+0.14, 95% +0.08 to +0.23; 1,744 pairs, 114 regulated). The gain does
  not depend on which pairs the sweep covered, and the combined model sits in the range of the
  published ENCODE-rE2G model on the same pairs. Tested in one cell type."
- Not allowed: "exceeds ENCODE-rE2G" (held-out band flattered by the H3K27ac filter; unpaired);
  "generalises across cell types" (HCT116 unrun, GM12878 uninformative); "a strong baseline" (this
  project's own activity baseline is below published ABC on training).
- Open, and handed to the coordinator: HCT116 at 705 requests, the only purchase that could make it
  two cell types.

The outside-reader summary is `docs/CRISPRI-RESULT.md`.

## Pre-registration: the last one-target consumers under one registration, and the constrained-unknown headline against its own matched random windows (2026-09-27, third)

Four modules still read one target gene per element from the compact
`all_elements` tables: `motif_transfer.tier_blocks` (the sweep lead per block),
`syntax_tiling.moved`, `candidates.node_context` and the locus benchmark's
`read_gene_input`. They are moved onto `attribution.targets.ElementResponses`
together, under this one registration, whose code form is
`PRE_REGISTRATION` in `scripts/onetarget_consumers.py`. The old reading is
kept beside the new one in every case, so each run carries its own control.
0 AlphaGenome requests; one chromosome's cache archive held at a time.

**The shared premise, registered as expected false.** The roadmap line said an
element that moves a non-top gene reads as "did not move". The compact head is
the maximum over the window by `predict_target`'s own rule, so "some gene is at
the bar" and "the head is at the bar" are one statement. The survey measured it
on chr21 before this registration: 7,846 elements yes/yes, 4,293 no/no, 0
disagreements. What a one-target reader loses is how many other genes move and
which, never the yes or no. Shared control: the window's head at the bar names
the table's `predicted` gene and its coding head names `predicted_coding`, with 0
disagreements required; an uncached element keeps its compact reading and is
counted by name.

| Consumer | Must reproduce to the digit | May change, and which way | Falsifier |
| --- | --- | --- | --- |
| motif_transfer | real unknown 531 tested / 331 with a lead; neutral 1,181 / 768; leads by case 11, 167, 19, 131, 3; hence the agreement AUROC 0.5422 by construction (count model not re-run) | a new per-block count of genes at the bar, expected at or above the distinct head genes, about 2.3 per moving element | any block's lead or tested count differs |
| syntax_tiling | nothing: the registered run was cancelled unrun | its rows gain the genes at the bar from the cache when it runs; `moved` unchanged | a cached element where `moved` and "some gene at the bar" disagree (unit test) |
| candidates | the chr22 candidate's node fields (the only candidate chromosome complete on 2026-09-13) | the other 68 differ from the stored file by input, counted and attributed to it; per gene, window count >= table count, strictly; the top node gene two-sided, counted | any gene whose window count is below its table count |
| loci gene_input | all 17 loci's stored gene_input fields; strict hits 10 of 17, among 15 of 17 | strict hits hold or rise; the misses' median rank holds or improves (MYC_8q24 and PMP22_CMT1A, rank 2, likeliest to flip) | strict hits below 10, or any of SHH_ZRS, HERC2_OCA2, FTO_IRX3, SOX9_PierreRobin (where gene_input alone carries the headline) losing its hit |

The benchmark's `target_derived` 15/17 is not recomputed in place: swapping a
layer of the shared benchmark is its owner's call, so the counterfactual rate
with the window layer is reported beside it. The candidates' labels and set
tests stand as stored.

**The constrained-unknown headline gets its own control**
(`CONTROL_REGISTRATION` in `scripts/constrained_unknown_targets.py`).
`matched_random_windows` is lifted into the script and first checked against
chr21: handed the same 446 blocks and the same predicate it must reproduce
7f7c8c1 to the digit (22,200 drawn, 2 undrawable, 16,706 carrying an element,
10,908 naming; 140 of 285 blocks). Any difference and no genome-wide figure is
reported. Then every field the script wrote on 2026-09-16 must reproduce (331 of
531 blocks with an element carry a mover, 0.6234; 0.2476 per element; tiers,
cases, deciles, the 331 named blocks), and the primary is new: the real
unknown's 0.6234 against 50 windows per block of its length, drawn in the
chromosome's scored span and rejected on overlap with any organiser block,
pooled over 24 chromosomes. Readings, with a band of 5 points: **below chance**
is a live possibility and the expected one (chr21, coding targets: 49.1% against
65.3%) and would say a lead at a real-unknown block is rarer than at random
sequence outside unknown space, which is what the 331 leads are worth, not a
failure; **at chance** says a lead carries no block-level information; **above
chance** would be the first above-chance reading of the leads. Secondary: the
same control on "names a coding gene", the neutral tier as a second target set,
and genes at the bar per block. Caveat registered now: the random windows avoid
all unknown space and so sample gene-richer sequence; the tier comparison in
length deciles remains the within-unknown control. The borrowed 87% sentence is
replaced by whatever this measures.

## The result: the last one-target readers moved onto the cache, the locus layer lost HOXD, and the real unknown names a target 24 points below chance (2026-09-27, fourth)

The run of the registration above. `data/results/onetarget_consumers.json` and
`data/results/constrained_unknown_targets.json`; 0 AlphaGenome requests; one
chromosome's cache archive held at a time. Negatives first.

**The locus benchmark's window rule is falsified.** Crediting every coding gene
at the bar, not only each element's coding head, takes the gene-input layer's
strict hits from 10 of 17 to 9: HOXD is lost, gained by none. HOXD was a narrow
table win (HOXD1 13.859 against EVX2 13.199, summed |log2|), and the window
reading hands EVX2 the lead as rank 2. No sole-carrier locus lost its hit
(SHH_ZRS, HERC2_OCA2, FTO_IRX3 and SOX9_PierreRobin all hold), and the misses'
median rank of the first published target improved from 6 to 5, but the
registered falsifier is "strict hits below 10", and it fired. The one-target sum
stays the better reader of this layer; `read_gene_input` keeps it and carries
the window reading beside it as `window_reading`, never in place of it. The
counterfactual `target_derived` rate with the window layer is 15 of 17, the same
as stored, because HOXD is also hit by another derived layer.

**The constrained-unknown headline reads below chance, as registered.** The
lifted `matched_random_windows` reproduced 7f7c8c1 on chr21 to the digit
(22,200 drawn, 2 undrawable, 16,706 carrying an element, 10,908 naming; 285 and
140 blocks), so the genome-wide figure is reported. Of the 531 real-unknown
blocks carrying a scored element, 331 carry one that moves a gene, 62.3%,
against 86.0% of 44,100 length-matched random windows outside the organiser's
blocks (27,254 of 31,676 carrying an element): -23.7 points, **below chance**.
On "names a coding gene" the gap is -37.2 points (29.8% against 67.0%). The
neutral tier reads the same way, -19.3 and -36.0 points (65.0% against 84.3%;
28.8% against 64.8%), so the gap is the organiser's blocks against gene-richer
sequence outside them, the caveat registered, rather than something special to
the real unknown; inside length deciles the real unknown stays 4.6 points above
neutral, the within-unknown control, unchanged. A lead at a real-unknown block is rarer than at random sequence, and
the borrowed 87% sentence is replaced in the result by this measurement.

**What reproduced to the digit.** Every field `constrained_unknown_targets`
wrote on 2026-09-16 (331 of 531, 0.6234; 0.2476 per element; tiers, cases,
deciles, the 331 named blocks). motif_transfer: real unknown 531 tested / 331
with a lead, neutral 1,181 / 768, leads by case 11, 167, 19, 131, 3, so the
agreement AUROC 0.5422 stands by construction. candidates: the chr22
candidate's node fields reproduced; the other 68 differ from the stored file by
input (93 elements in their nodes then, 7,285 now), as registered. loci: all 17
stored `gene_input` fields reproduced on the input the stored benchmark had;
on today's input MYC_8q24 differs because stated intervals were added after
7a1c751.

**The shared premise was false, as registered.** Over every element read from
both sources (849,469 for the candidates' chromosomes, 4,794 in the locus
windows) the window's head at the bar names the table's `predicted` gene, and
its coding head names `predicted_coding`, with 0 disagreements; no element in a
motif-transfer block was uncached, and no block's mover count moved. What the
table lost is the other genes: over the 331 led real-unknown blocks, 849 genes
at the bar against 416 head genes (179 blocks name more); neutral 2,139 against
932. The registered "about 2.3 per moving element" was the survey's per-element
figure; counted distinct per block it is 849 over 812 moving elements, 1.05,
because neighbouring elements share genes. candidates: no gene's window count
fell below its table count (the falsifier), node genes 223 to 267, and the top
node gene changed at 3 of 58 candidates with one (SYP to FOXP3 on chrX, MAP3K12
to SP7 and HOXC6 to HOXC11 on chr12). syntax_tiling: not re-run, as registered;
its rows now carry `genes_at_bar` from the cache the scorer writes, and the unit
test that `moved` equals "some gene at the bar" passes on 2,000 random records.

## The CRISPRi split audited: no held-out figure read the held-out pairs through the merged loader or the compiled program, and 43 held-out links sat in the compiled programs unmarked (2026-09-28)

Two defects were reported in `attribution/measured.py`. First, `CrispriPair` dropped the
benchmark's power columns, so a well-powered negative and a blind one counted the same. Second,
`load_crispri` merged the training file with the held-out file, and `compile.py _measured_blocks`
wrote the regulated links from both into the compiled BioLang programs as experimental facts. The
audit came first and was committed before any code changed. It traced every consumer of
`load_crispri`, of the measured layer's CRISPRi data and of the `<id>_measured` blocks through
`genomeos/` and `scripts/`. It then checked every other CRISPRi reader behind a held-out claim, to
confirm that none of them reaches the two paths indirectly. The table is in
`data/results/crispri_split_audit.json`.

**Exposed: none.** No scorer, fit or feature behind a held-out or validation figure read held-out
pairs through either path. Every held-out claim comes from a module that reads the two files
separately with `crispri.load` and fits on training rows only. The defect is latent, not realised:
the compiled programs carried **43 held-out-only links out of 212 measured rules** across the 24
chromosomes. The committed `data/organisms/human/noncoding_chr21.bio` carried 1 of its 2 (ICOSLG,
GM12878). The first feature to read a compiled program would have seen the test set.

| consumer | path | figures | status |
|---|---|---|---|
| `crispri.score` (`scripts/crispri_score.py`) | `crispri.load` per file | held-out K562 AUPRC 0.550 → 0.691, +0.141 [+0.082, +0.231] (`crispri_benchmark.json`, CRISPRI-RESULT.md) | CLEAN: fit on training, held-out scored frozen |
| `crispri.score_published` (`scripts/crispri_published.py`, 42d7b1b) | `crispri.load` per file | training LOCO 0.507 / 0.724; held-out pooled 0.567 / 0.677 registered, 0.476 / 0.639 DNase-only; coverage arms +0.136 | CLEAN; the post hoc DNase choice is already disclosed and is not a leak |
| `crispri.score_contact` | `crispri.load` per file | `crispri_contact.json` held-out arm | CLEAN |
| `target_calibration` and its gate and prevalence scripts | `crispri.load` per file | `heldout_k562`, `heldout_other_cells`, `held_out_before_and_after` | CLEAN: every `fit` and `fit_with_offset` runs on training |
| `target_rebanding`, `union_axis` band levels | `crispri.load` per file | band rates | CLEAN for these paths; the rates pool training and held-out K562 by registration, so no later held-out K562 figure may be read against the re-banded table as out of sample |
| `crispri_direction`, `crispri_direction_both` | `crispri.load` per file | per-cell sign agreement | EVALUATION-ONLY: nothing fitted |
| `loci_fourth`, `loci_noncoding` | `crispri.load(HELDOUT)` only | locus rates | EVALUATION-ONLY: registered 2026-09-21 before the file was read |
| `node_containment_audit.py` stage 2 | `crispri.load` of both files | 661 regulated pairs, measured share 0.758 | EVALUATION-ONLY: the pairs are ground truth for a caller built from CTCF cCREs, nothing is fitted, and the figure is not quoted as held-out |
| `measured.rows`/`census` (`scripts/measured_layer.py`) | `load_crispri`, merged | 19,070 of 440,377 elements measured; agreement over 1,505 and over 128 | EVALUATION-ONLY: fixed rules, pooled by design, not a held-out claim |
| `confidence_calibration` | `load_crispri`, merged | 128 elements, observed 0.758 against stated 0.347; 10 of 12 bands outside | EVALUATION-ONLY: bands imported unchanged, nothing fitted |
| `compile._measured_blocks` → `noncoding_chr*.bio` | the compiled program | 43 of 212 rules held-out-only | REPORT-ONLY today: its readers are the evidence census, the CLI, the explorer index, `bio test` and three syntax censuses |
| `measurability.py`, `unknown_coverage.py` | own file lists, intervals only | coverage footprints | REPORT-ONLY: no label read |

**What the fix is therefore allowed to change: no figure.** `load_crispri` keeps its default of
both files, because no exposed caller exists to move. The fix, committed after this audit, gives it a `split` argument
(`"training"`, `"heldout"`, `"all"`). Each `CrispriPair` gains its `split`, taken from the
file it came from, and the five power columns exactly as the headers name them:
`PowerAtEffectSize10`, `15`, `20`, `25` and `50`. In both files every valid row has a value. The
training file holds 9,885 negatives, 6,193 of them with power at least 0.8 at a 20% effect (the
brief's "6,169 of 9,810" was an approximation of this). The held-out file holds 4,188 negatives,
3,675 of them well powered.

The compiled program **keeps the held-out links, marked**, and this is the least invasive choice. Removing them would change the rule
counts every program tests itself on, and it would make the program omit a real measurement. The
chosen mark is carried by the structured fields a feature would read, not only by a comment:

- the `targets:` of a `_measured` element lists only genes regulated in a training pair;
- a link found only in held-out pairs is written with `measured.HELDOUT_MARK` in its evidence
  source, and `measured.is_heldout(source)` is the test any future reader must apply;
- a link found in both splits takes its action and strength from the training pairs.

A test is to assert that no held-out-only gene reaches `targets:` and that every held-out-only rule
carries the mark.

That test now exists: `tests/test_crispri_split.py`.

**The fix, as built.** In the committed chr21 program, `targets: ICOSLG` left the held-out element
and its rule gained the mark. The rule count stayed 2. Genome-wide, 1 of the 169 training-supported
rules had taken its strength from a stronger held-out pair; it now takes it from the training pair.
No reported figure moves. The `rules`, `entities` and `unknowns` counts each program tests itself
on do not move either, because every measured link is still stated. Only the chr21 lines that the
fix touches were rewritten in the committed copy. A full recompile would also have pulled in
unrelated drift in the reader comments since 2026-09-17. The git-ignored programs under
`data/knowledge/compiled/` pick up the mark on the next
`scripts/measured_layer.py --write-programs`.

**Evaluation independence (review item R5).**
- *Provenance.* Every `CrispriPair` now carries its partition (`split`), its `source_file`, its
  `assay` and its `study` (the benchmark's `Dataset` column). All four are frozen at parse time.
- *Guard.* `measured.development_only` raises on a held-out pair rather than filtering it out
  silently.
- *Overlap check.* `measured.split_overlap` compares the two partitions by interval as well as by
  identical pair, and a test recomputes it from the files against
  `data/results/crispri_split_audit.json`.

It finds **no identical pairs, but 249 of the 4,378 held-out pairs share at least one base with a
training interval**:

- 36 are on the same gene: 4 near-identical (reciprocal overlap ≥ 0.9) and 32 overlapping. By
  cell, 23 are in HCT116 and 13 in WTC11.
- 213 are on another gene. 66 of those sit on a near-identical element, and 164 of the 213 are in
  K562.

No held-out K562 pair repeats a training gene at an overlapping interval. `crispri.score` already
reports an arm restricted to held-out elements absent from training; every other held-out figure
includes the related pairs. The related-variant half of the check has nothing to compare here,
because the benchmark tests intervals, not alleles.

**Reused benchmark.** At least ten scored results have read the held-out file: crispri_benchmark,
crispri_contact, crispri_published, crispri_direction, crispri_direction_both, target_calibration,
target_calibration_gate, target_prevalence, target_rebanding and union_axis. loci_fourth and
loci_noncoding read its held-out arm too. Each registration froze its model before its own read.
The choices between registrations, however, were made by people who had seen the earlier reads. In
the audit JSON, every row on this file is therefore flagged `reused_benchmark`, and its figures are
to be quoted as a reused benchmark, not as a first touch of a test set.

**Correction to the count above.** The brief's "6,169 of 9,810" is exact, not an approximation.
9,810 is the number of *non-significant* training pairs, and 6,169 of them have
`PowerAtEffectSize20` ≥ 0.8. The 9,885 quoted above also counted the 75 significant increases,
which the benchmark marks `Regulated` FALSE.

**CRISPRi negatives read as outcomes (review item R2).** In the benchmark file, `Regulated` means a
significant *decrease*. A significant increase is therefore marked FALSE as well, and until today
`for_element` listed such a gene as "measured no effect". Each pair now has one of five outcomes,
kept apart before anything is pooled. The rule was fixed before any verdict was recomputed:

- `significant_decrease`;
- `significant_increase`;
- `not_significant_well_powered`, meaning at least 0.8 power at a 20% effect;
- `not_significant_underpowered`;
- `missing`.

Invalid comparisons (`ValidConnection` FALSE) are still counted and never read as outcomes. Both
files hold 0 missing and 0 invalid rows. The training file has 471 decreases, 75 increases and
9,810 non-significant pairs, 6,169 of them well powered. The held-out file has 190, 84 and 4,104,
with 3,632 well powered. The 25% column cannot separate the nulls: every non-significant pair in
both files is at least 0.8 there, which is the benchmark's own filter.

The agreement rules change as follows:

- A gene takes the strongest outcome among its pairs.
- Only a well-powered null disagrees with the compiled claim.
- An underpowered null is `predicted_gene_null_underpowered`.
- An increase is `predicted_gene_increased`. Neither of these two counts as agreement or as
  disagreement.
- An increase raises no rule and no silencer label.

The figures that move are these, old → new:

| figure | old | new |
|---|---|---|
| CRISPRi disagreements (`measured_layer_genome`) | 31 | 23 |
| elements where the predicted gene was tested | 128 | 120 |
| agreement rate over them | 0.7578 | 0.8083 |
| pairs measured as not regulated | 4,613 | 4,565 |
| predicted gene significantly increased (new) | — | 1 |
| predicted gene null but underpowered (new) | — | 7 |
| pairs by outcome (new) | — | decrease 212, increase 48, well-powered null 3,130, underpowered null 1,435 |
| confidence calibration, CRISPRi where tested: observed / mean stated | 0.7578 / 0.347 | 0.8083 / 0.3576 |
| expected calibration error | 0.4108 | 0.4507 |
| median gap over populated bands | +0.4724 | +0.5432 |
| does the stated confidence carry information? matched difference (p) | 0.1374 (0.070) | 0.0898 (0.156) |
| the same, raw (p) | 0.1886 (0.0035) | 0.1208 (0.040) |
| committed chr21 program, pairs "measured as not regulated" in its header | 32 | 31 |

Some figures do not move: the band headline (the stated level lies outside its own interval in 10
of 12 bands), every pattern verdict, the all-matched rate of 0.0645 over 1,505, the coverage of
19,070 of 440,377, and the `rules`, `entities` and `unknowns` counts in every program. The
figures of 2026-09-17 in the sections above (128, 0.7578, 4,613 and the calibration table) are
superseded by this table and are left as they were written. The 24 git-ignored programs were
recompiled (`scripts/measured_layer.py --write-programs`), so the held-out mark and the five
outcomes are now on disk. In the committed chr21 copy, only the header line and the one basis line
that changed were rewritten.

**Not yet in the committed chr21 copy.** Two lines are held back from this commit: the header
count (32 → 31) and the SUMO3 basis line. The basis line replaces one this lane wrote in `36fa298`
the same night, and the checkout's removal guard refuses that without `--force`, which this lane
does not use. The programs in `data/knowledge/compiled/` already carry both changes. The next
recompile of the committed copy will pick them up.

## Pre-registration: context in executable rules, one rule per element, gene and cell (R1, 2026-09-28)

The external review of 2026-09-28 (ROADMAP section 5 item 11, R1) found two things. The compiler
wrote the cell a CRISPRi screen silenced an element in only into evidence text ("silenced in
K562"), so a rule measured in K562 ran in every simulated cell. And `measured.rule_links` kept one
rule per element and gene, the strongest result across every cell that measured it, so two cells
that disagreed became one number.

**Route.** BioLang already has the condition form: a rule's `when: k = v` clause, which the
network runtime honours through `Module.active_rules` (a rule whose `when` does not match the
context is not integrated). Nothing new is added to the grammar's keys. The compiler writes
`when: cell_type = <cell>`, the context key `measured.CONTEXT_KEY`. One addition is needed in
the runtime: `cell_type = unknown` (`measured.CONTEXT_UNKNOWN`) must match no cell. Today a
`when` value is compared by plain equality, so `unknown` would match a context that happened to
be called "unknown"; the addition makes it match nothing, and the grammar's text says so.

**Fixed before the build.**

1. The review's acceptance tests, as it words them (`measured.R1_ACCEPTANCE`): a K562-specific
   rule is inactive in HepG2; conflicting results from two cell types survive compilation and
   round-trip serialisation (compiled text, parse, `Module.to_dict`, JSON, `from_dict`).
2. A rule with no recorded context is written `when: cell_type = unknown`, never without a
   `when`. That holds for measured rules (a pair with an empty cell) and for the predicted
   AlphaGenome rules, which name the track's biosample in their evidence ("AlphaGenome deletion,
   K562") and are gated on it the same way. A predicted rule whose track was not recorded
   ("strongest track") gets `unknown`.
3. The unit of a compiled experimental rule is one (element, gene, cell) observation
   (`measured.RULE_UNIT`). Within one cell, the split rule of the audit above is kept: a link
   with training support in that cell takes its action and strength from the training pairs
   only, and a link held out in that cell is marked `HELDOUT_MARK`. An element-gene link
   regulated in K562 (training) and in WTC11 (held-out) was one training rule and becomes two
   rules, the WTC11 one marked.
4. What still makes no rule: a measured null and a significant increase (the benchmark's
   `Regulated` label excludes positive effects, and the audit above pinned that an increase
   raises no rule). Where one cell's null or increase sits beside another cell's regulated link
   for the same gene, the measured element's basis text names it cell by cell, because
   `for_element` today keeps only the strongest outcome per gene and drops the other cell's.

**Expected counts, with the direction.** Predicted rules: unchanged in number (one per element),
all gated. Experimental rules: the same or more per program, never fewer, because one rule per
(element, gene) becomes one per (element, gene, cell). Every program's `# test: rules ==` line
moves by the same amount. A simulation run with no `cell_type` in its context now activates no
compiled rule; before, every one ran. Old and new counts are reported with the result.

**Falsifiers.** Any compiled rule without a `when`; an experimental rule count below the old one;
a training rule whose strength changes; a HepG2 run that integrates a K562 rule.

## The result: every compiled rule now runs only in its own cell, and no measured link in the genome was regulated in two cells, so no rule count moved (R1, 2026-09-28, later)

Built as registered. `measured.rule_links` yields one link per (gene, cell), with the split rule
applied within each cell; the compiler writes `when: cell_type = <cell>` on every experimental and
every predicted rule; `Rule.applies` and `matches` treat `unknown` as matching no context, and the
grammar's `when` text says so (BIOLANG-GRAMMAR.md regenerated; `package_engine.py` 12 of 12).
`for_element` keeps each cell's own outcome per gene (`outcomes_by_cell`), and the measured
element's basis names every cell whose outcome the pooled per-gene lists hide. Tests:
`tests/test_rule_context.py`, 9, written first and failing 8 of 9 before the build.

**Counts, all 24 programs recompiled (old → new).** Rules 440,589 → 440,589; experimental rules
212 → 212; held-out-marked experimental rules 43 → 43; predicted rules 440,377 → 440,377. Rules
without a `when`: 440,589 → 0. Rules gated on `unknown`: 0 (every predicted rule named its track's
biosample, and every CRISPRi pair its cell). Every program's `# test: rules ==` line is unchanged
and still holds. Measured element blocks 19,072 → 19,072. Basis lines that name a cell-level
difference: 0 → 2, both held-out, both a link in one cell beside an underpowered null in the
other (IL6ST on chr5, regulated in Jurkat, underpowered in GM12878; PPIF on chr10, regulated in
GM12878, underpowered in Jurkat). Before, the null in the second cell was dropped from the element
without a word. Across the genome the rules fall into 317 contexts. K562 has the most, 27,445.
chr21 in a HepG2 run integrates 154 of its 5,176 rules; before, it integrated all 5,176.

**What the zero means.** The experimental count could only move if one element-gene link was
regulated in two cells, and in this benchmark, intersected with the compiled elements, none is.
The collapse the review found was real in the code and had no case to act on in the data. The
synthetic tests are the only place a two-cell conflict exists so far, and they pin that it
survives compilation and the BioIR round trip (text, parse, `to_dict`, JSON, `from_dict`).

**Negatives.** A compiled program run with no `cell_type` in its context now integrates none of
its rules. That is intended, but any caller that ran one without a context got every rule before
and now gets none. `Rule.applies` and `Event.applies` still compare a `when` value by plain
equality, so the grammar's alternatives (`a|b`), `absent` and comparisons work for decisions and
timers, which call `matches`, but not for rules or events. That was already the case before R1 and
is not changed here. Two tests in
`tests/test_crispri_split.py` pinned one rule per element and gene. They are now strict xfail, not
deleted, because the checkout's removal guard refuses to take out lines that lane-split wrote the
same night. Their invariant is re-pinned per cell in the new file.

## Pre-registration: one independent benchmark, the IGVF MHC CRISPRi screen in K562 (review item 11, 2026-09-28)

Review item 11 asks for one independent biological benchmark before any widening of scope, and for
generalisation shown on evidence that is truly independent before any claim of discovery. lane-split's
audit (`0af7e1b`) found the ENCODE CRISPRi held-out file reused by at least ten scored results, so it
cannot serve as that benchmark any more. This section fixes a test on a screen that no script here has
read. The test is registered before scoring and costs 0 AlphaGenome requests.

**Census of what the attribution layer has read** (grep of `data/knowledge/`, `scripts/`, `genomeos/`,
`docs/DATA.md`). CRISPRi: both ENCODE EPCrisprBenchmark files, i.e. Nasser 2021 (Fulco 2016/2019,
Ulirsch, Klann 2017, Xie 2017 and the rest), Gasperini 2019, Schraivogel 2020 in training, and Xie 2019,
Morris 2023 (STING-seq), Klann 2021, Reilly 2021, Ray 2025 DC-TAP (K562, WTC11), Guckelberger 2024
(HCT116), Nasser GM12878/Jurkat held out. Activity: ENCODE4 lentiMPRA ENCSR106SZM (K562, HepG2, WTC11)
and Kircher 2019 saturation mutagenesis (GSE126550). Other sources: VISTA, GTEx v8 eQTLs with DAP-G,
the GWAS Catalog, ClinVar non-coding, and the reporter-allele sets behind E1 (Tewhey, Abell). Features:
the ENCODE cCRE registry V3, the ENCODE epigenome, and 4DN Hi-C (GM12878, K562). Nothing from the IGVF
portal has been read.

**Candidates checked on the primary portals** (ENCODE and IGVF APIs, 2026-09-28):

| Candidate | Measures | n | Cell | Licence | Coverage with 0 requests |
| --- | --- | --- | --- | --- | --- |
| **IGVF K562 MHC CRISPRi Perturb-seq** (Gersbach; IGVFDS7132YVKO, file IGVFFI4093WUVB, 0.96 MB, released 2026-06-08) | element → gene, signed effect and FDR (FRACTEL) | 34,279 pairs, 581 DHSs, 59 genes, chr6 | K562 | CC BY 4.0 | 4,472 of 6,544 eligible pairs within 500 kb (68%); **chosen** |
| ENCODE Flow-FISH K562 genes outside the benchmark (TRIR, OTUD5 Engreitz; FADS2, PVT1, HBS1L Sabeti) | element → one gene | 5 genes | K562 | ENCODE, open | yes, but five loci are too few for a readable test |
| IGVF K562 CRISPRi Perturb-seq "many loci" (Gersbach, 16 sets, 2026-07-27) | element → gene | not processed | K562 | CC BY 4.0 | no released element-level table; raw reads only |
| IGVF HCT116 FlowFISH (Engreitz: CCND1, MYC, KITLG, FAM3C, ITPRID2) | element → one gene | 156 sets | HCT116 | CC BY 4.0 | the sweep has no HCT116 track; this is the 705-request item |
| IGVF HepG2 electroporated MPRA (Mohlke) and CYP3A4 MPRA (Ahituv) | variant allelic effect | about 10⁴ variants | HepG2 | CC BY 4.0 | the sweep deletes whole elements and never scores alleles; not covered |
| IGVF lentiMPRA IGVFDS9175DMRR | element activity | genome-wide | HEK293T | CC BY 4.0 | the sweep has no HEK293T track |
| Morris 2023 STING-seq | element → gene | — | K562 | — | already read (in the held-out file) |

**The published model.** scE2G v1.2 in K562 (IGVF prediction set IGVFDS5428HHMB, unfiltered file
IGVFFI1706PNVV, 326 MB, CC BY 4.0). Its chr6 27–36 Mb rows are distilled to
`data/cache/indep/scE2G_K562_chr6_MHC.tsv.gz`. It was trained on the ENCODE CRISPRi benchmark, which
does not contain this screen.

**What was seen before this registration.** The column names, the table's shape (34,279 rows, 581
elements, 59 genes, chr6 only), the count of rows at FDR < 0.05 in either direction (385), and a
label-blind coverage count. Of 20,228 eligible pairs, the sweep covers 4,651. Within 500 kb of the TSS
it covers 4,472 of 6,544. Beyond 500 kb it covers almost nothing, because the scorer read only the genes
inside its 1 Mb window. No label was joined to any prediction.

**The registered test** (`genomeos/attribution/indep.py`, `PREREGISTERED`; `scripts/indep_benchmark.py`
is fixed in the same commit):

- **Pairs.** One DHS and one gene. Excluded and counted: elements within 1 kb of any GENCODE v50 TSS
  (12,921 pairs), genes absent from GENCODE v50 chr6 (1,086), and pairs whose gene and an overlapping
  element also appear in either ENCODE benchmark file (44). The primary population is the pairs within
  500 kb of the TSS, the sweep's reach; the 13,684 pairs beyond it are counted and not scored.
- **Labels (R2).** A decrease is FDR < 0.05 with effect < 0, and it is the positive. An increase is kept
  apart and reported alone. "Not detected" (FDR ≥ 0.05) is the negative. The screen gives no power, so
  "not detected" does not mean "no effect".
- **Coverage.** A pair is covered when a sweep element overlaps the DHS and the per-element cache has
  a K562 value for the gene. Pairs that are not covered are missing: never imputed and never scored as
  zero. The test is not read if fewer than 20 covered positives remain, or if the covered share falls
  below 0.5. As a selection check, the distance AP on uncovered pairs is reported beside the AP on
  covered pairs.
- **Predictor and baselines.** The predictor is the deletion's largest predicted fall on K562's own
  track. The baselines are:
  - distance: minus the log distance from the element to the TSS;
  - the current annotation: 1 when the gene is the element's compiled top target;
  - the node rule: 1 when the gene is the nearest TSS in the element's CTCF node;
  - scE2G's `Score`, compared on the pairs that both it and the sweep cover.
- **Metric.** AUPRC (average precision, with tied scores counted as one threshold). The 95% interval
  comes from 2,000 bootstrap resamples of elements (seed 20260928), with the paired difference taken
  on the same resamples. A resample of the 59 genes is reported as a sensitivity check. The base rate
  is reported beside every AP.
- **Pass.** The 95% interval of AP(deletion) − AP(distance) lies wholly above 0.
- **Falsifier.** That interval lies wholly below 0, or the upper bound of AP(deletion) does not exceed
  the base rate. In either case the result is FAIL.
- **Anything else is NOT ESTABLISHED.**
- **Readings.**
  - PASS: the ranking transfers to one new K562 screen at one locus, and to nothing wider.
  - FAIL: the ranking that the ENCODE file supported does not transfer. The target bands and the
    calibration must then be restated as unvalidated outside that file.
  - NOT ESTABLISHED: this screen cannot separate the deletion from distance, and that is not support.
  - NOT READABLE: coverage is the whole result.

  One locus bounds every reading. The element bootstrap does not sample variation from locus to locus,
  and the MHC is unusual. AlphaGenome was trained on ENCODE K562 tracks, so the objection that it has
  seen K562 chromatin remains open.

## The result: the independent MHC screen cannot be read, because only 19 distal pairs in it are measured decreases and 13 of those fall within the sweep's reach (review item 11, 2026-09-28, later)

The test registered in `042a4aa` was scored once (`scripts/indep_benchmark.py`, result
`data/results/indep_mhc_crispri.json`, partition `independent`, source IGVFFI4093WUVB, 0 AlphaGenome
requests). **Verdict: NOT READABLE.** No AUPRC was computed, and none is quoted here, because the
registration forbids reading a test below its floors.

**Why it cannot be read.** The registration requires at least 20 covered positives; the screen supplies
13. Of the screen's 93 significant decreases, 74 lie on promoters: the element sits within 1 kb of a
GENCODE v50 TSS, so these are the TSS controls and promoter hits the registration excludes. That leaves
**19 distal decreases in 21,358 distal pairs**. Of those 19, 14 lie within 500 kb of the gene's TSS, and
13 of the 14 are covered. Coverage is not the limit: 4,472 of 6,544 pairs within reach are covered
(68%, above the 0.5 floor). The shortfall is in the measurement: in this locus the screen finds almost
no distal element that lowers a gene.

**The measurement leans the other way.** Among distal pairs the screen calls **181 significant
increases against 19 decreases**; within reach the count is 72 against 14. Under R2 these are kept as
their own outcome. A CRISPRi screen whose distal hits are mostly increases is measuring something other
than enhancers being silenced, such as indirect effects, silencers, or the FRACTEL test's pooling of the
two directions. Before this screen is used as a target benchmark in any form, its authors' reading of
those increases is needed. The sweep predicts no fall for 16 of the 54 covered increases; that is
descriptive, not a test.

**The selection check**, reported as registered and read as descriptive only: the distance baseline's
AP is 0.119 on covered pairs (base rate 0.0029) and 0.25 on uncovered pairs within reach (base rate
0.0005, one positive). With one positive, that second number means nothing.

**What this means for item 11.** The project still has no independent target benchmark, and that has to
be said rather than papered over. Every target claim continues to rest on the reused ENCODE file, so it
remains in-sample by R5's definition. The census found no other study that the stored sweep can score
with 0 requests and at sufficient size. What remains:

1. the IGVF K562 CRISPRi Perturb-seq over many loci (Gersbach, 16 measurement sets, 2026-07-27). It
   needs 0 model requests, but no element-level table has been released: using it means processing raw
   single-cell reads, or waiting for the IGVF analysis set;
2. HCT116 (Guckelberger 2024, IGVF FlowFISH). The sweep has no HCT116 track, so this costs the 705
   requests already waiting on the owner's approval;
3. the five ENCODE K562 FlowFISH genes outside the benchmark (TRIR, OTUD5, FADS2, PVT1, HBS1L). These
   cost 0 requests but give five loci, and each would need its own registration with a floor it can
   meet.

## Probe: no second independent benchmark exists yet, because the 16 many-loci K562 sets are raw reads and their analysis set is unreleased (review item 11, 2026-09-28, lane-indep2)

**Negative first.** The IGVF portal holds no processed element-gene table for the Gersbach K562 CRISPRi
Perturb-seq over many loci, and none for any other K562, HepG2, GM12878 or IMR-90 CRISPRi screen that
the project has not already read. No benchmark was registered and nothing was scored. The project still
has no independent target benchmark. `scripts/indep_igvf_probe.py` records the census and the two
label-blind counts below (`data/results/indep_igvf_probe.json`, 0 requests, portal read through
`https://api.data.igvf.org` on 2026-09-28, IGVF data CC BY 4.0).

**The 16 many-loci K562 sets** (IGVFDS5637QJRP, 8708IDCC, 5440YEBR, 6162UVFJ, 0567GXFW, 5791LZOT,
5632GXBM, 2463TEJP, 1341HYXT, 8094NBFV, 0146NONB, 6790CUPS, 6013FBZQ, 6260UZSB, 1474DQOO, 8491ILFO;
released 2026-07-27) each carry only fastq and seqspec files: 597 GB of reads, plus 94 GB in their
16 gRNA-capture auxiliary sets. Each set is `input_for` an analysis set that exists but is not released
(HTTP 403): IGVFDS3624DGBH for 15 of them and IGVFDS3481OUHP for one. The processed table is therefore
coming, and has not arrived. Three further K562 sets from 2026-02-06 (IGVFDS1915EGTT, 2136PGKM,
5321WVPQ) feed a third unreleased analysis set, IGVFDS1588RHDM, and use the MHC library, so they are the
same locus again.

**What the library says, label-free.** The guide file (IGVFFI9377CEMM, 39 kB, md5 2c1975d5…, released
2026-07-27) describes "350 DHS genome-wide with transcriptome-wide Perturb-seq and growth readout":
1,587 targeting guides on 342 elements across 23 chromosomes. Of those, 63 lie within 1 kb of a
GENCODE v50 TSS, 9 overlap an ENCODE benchmark element, and 161 of the 279 distal elements (58%) are
covered by the stored sweep, giving 5,676 element-gene pairs with a K562 value. This clears the 0.5
coverage floor before any label is seen. Whether it clears the floor of 20 covered distal positives is
unknown until the table is released.

**Processing it here is not feasible.** It would mean 691 GB of downloads against 19 GiB free disk and a
2 GB cap, so it cannot be done even one set at a time. The steps would be the IGVF Perturb-seq pipeline
(Cell Ranger or kallisto|bustools for gene and guide counts, then guide assignment, then SCEPTRE or
PerTurbo per element and gene), with several CPU-hours and 32–64 GB RAM per set. This lane stops here,
as its brief says.

**The one other released K562 table has already been read.** The Engreitz "K562 Random DC-TAP-seq"
(IGVFDS7288SJVF; file IGVFFI0957PYTA, 0.57 MB, md5 758e1162…, released 2026-08-08, replacing the revoked
IGVFFI2436YUJO and the archived IGVFFI9586LSUT) is Ray et al. 2025, and that study sits in the ENCODE
held-out file. A label-blind overlap check read only the coordinate, gene, type, power and promoter
columns. It found 1,289 of the 7,281 targeting pairs (920 elements) matching a K562_DC_TAP benchmark pair for the same
gene, and 88 more matching pairs from other ENCODE studies. Another 1,222 distal pairs with power ≥ 0.8 at a 15% effect are absent from the ENCODE file, but
they come from the same screen, cells and pipeline whose positives already informed the project's
thresholds. A remainder of a study already read is not independent evidence, so it was not registered.

**Everything else on the portal is out of reach at 0 requests.** The sweep scores only K562, HepG2,
GM12878 and IMR-90. The remaining released element-level CRISPRi tables are in WTC11 (Gersbach
IGVFDS0523KGZP; Engreitz IGVFDS4617SJVM and IGVFDS3911MOCN; Hon IGVFDS7340YDHF and four cardiomyocyte
sets), in WTC11-derived neurons and endothelial cells, or in HCT116 (the 705-request item). HepG2 has
only a phenotype FACS screen with guide counts (Sherwood IGVFDS6504OLWV). GM12878 and IMR-90 have none.

**Next.** Watch IGVFDS3624DGBH. When it is released, fetch its element-level file, count label-blind
coverage, and pre-register on it exactly as `genomeos/attribution/indep.py` did for the MHC. It needs 0
requests, and the covered share is already known to be 0.58.

## Census: every summary in the measured layer, what reads it, and how often it bites (R6, 2026-09-28)

Review item R6 says the measured layer keeps the highest reporter activity over matching intervals and
the strongest CRISPRi effect, and that such summaries favour extremes and hide disagreement. Before any
change, `scripts/measured_aggregation_census.py` (result `data/results/measured_aggregation_census.json`,
0 requests) lists every aggregation in `attribution/measured.py` at `0924a5c`, with the consumers
`git grep` finds, and counts genome-wide how often each one is taken over more than one observation.

**One summary decides a label: the lentiMPRA maximum (A1).** `Layer.for_element` keeps, per cell, the
largest log2(RNA/DNA) of every tile that meets the overlap rule, and `cells_active` is read from that
maximum. It feeds the lentiMPRA verdict, the pooled verdict and `assays_agreeing`, the census counts in
`measured_layer_genome.json`, the calibration table in `confidence_calibration_genome.json`, the
agree/disagree header of every compiled programme, and the `basis:` line of every `<id>_measured`
block. Of the 17,869 elements lentiMPRA measures, 16,863 match one tile, 975 match two and 31 three or
more. That gives 3,018 cell readings over more than one tile, and in **237 of them the tiles disagree
about the threshold**: 224 are ties (one of two tiles active), 7 have a minority of tiles active and 6 a
majority. Under the maximum, all 237 read active.

**The rest, and why they are out of scope or harmless.**
- A2, `max_activity`, and A3, the `overlap` of every assay: maxima that no code reads. They are kept but
  should be named as descriptive.
- A4, `lentimpra.elements`: tile names without values, strand or overlap. The disagreement above cannot
  be seen in a row today.
- A5, `mpra.parse`: averages the rows that share an interval in one cell's file. The average never
  fires, because each of the three ENCODE files (ENCFF802FUV, ENCFF475FKV, ENCFF769REH) holds one row
  per interval. The strand column is dropped (28 rows per file are on `-`). The p and q columns are -1
  throughout, and no file has a replicate-level value, so **no per-tile uncertainty exists to keep**.
- A6 and A7, the pooled CRISPRi gene lists and `outcomes_by_cell`: strongest outcome across cells and
  across datasets. R1 (`6a4a8a9`) left every pair and each cell's own outcome in the row and names
  differing cells in `basis:`, so these are left alone.
- A8, `rule_links`: a compiled rule's `strength:` is the largest |effect| among the training pairs of
  one element, gene and cell. The compiler reads it. It is reported to the coordinator, not changed here.
- A9, VISTA: an element agrees when any matched VISTA element is positive. 5 of the 98 VISTA elements
  match more than one, and **none matches a positive and a negative**.
- A10, saturation mutagenesis: only a locus's primary experiment is read, and its repeats are counted but
  never looked at. 2 repeats sit under the 4 matched elements, both in SORT1's group (SORT1.2 primary,
  plus SORT1 and SORT1-flip, which is the same element in the other orientation). This is the one
  orientation measurement the layer holds, and it is discarded.
- A11, `confidence_of`: the strongest assay present sets the row's confidence. The docstring names this,
  and R4's lanes own confidences.
- A12, the pooled verdict: already keeps conflict as `mixed`.

**The word "episomal".** It describes the ENCODE4 lentiMPRA in `measured.py`,
`confidence_calibration.py`, `scripts/unknown_coverage.py`, four places in this file,
GRAMMAR-BY-COMPARISON and ROADMAP. The cited assay integrates its reporters by lentivirus (Agarwal et al.
2025, Nature, https://www.nature.com/articles/s41586-024-08430-9). An integrated reporter still sits
outside its native locus, so the conclusion drawn from the word ("it measures the sequence and not the
locus") stands, but the word itself is wrong.

## Pre-registration: reporter tiles kept, the label from a declared share rule, the maximum descriptive only (R6, 2026-09-28)

Fixed in `attribution/measured.py` (`REPORTER_UNIT`, `REPORTER_LABEL_RULE`, `REPORTER_SUMMARY`,
`REPORTER_MAX_FIELD`, `REPORTER_UNCERTAINTY`, `TILES_CONFLICT`, `OUTCOME_KIND`, `R6_ACCEPTANCE`) before
any row is rebuilt, on the census above.

**The observation.** One lentiMPRA tile in one cell: its interval, strand, reciprocal overlap with the
compiled element and log2(RNA/DNA). Every matched tile is listed in the row under `tiles`. Uncertainty
is recorded as absent, with the reason, because the files hold none.

**The aggregation model.** Per cell, a label comes from the share of tiles on each side of
`MPRA_ACTIVE` (1.0): **active** when more than half are at or above it, **silent** when more than half
are below, **conflicting** when exactly half are. The per-cell `activity` reported is the median of the
tiles, and it is descriptive: the label is not read from it. The maximum stays only as
`activity_max_descriptive` (per cell) and `max_activity_descriptive` (over cells), and no verdict reads
either. The lentiMPRA verdict: **agrees** when some cell is active, **disagrees** when every cell is
silent, and **`reporter_tiles_conflict`** otherwise. The last counts as neither agreement nor
disagreement, and the census reports it under its own key. Across cells the rule is unchanged: activity
in any cell is support, and `cells_active`, `cells_silent` and the new `cells_conflicting` list each cell.

**Outcomes kept apart.** Each assay block carries `outcome_kind`: endogenous gene regulation (CRISPRi),
reporter activity from an integrated lentiviral construct outside the native locus (lentiMPRA),
developmental activity in the mouse embryo (VISTA) and base-level sensitivity in a reporter (saturation
mutagenesis). The pooled verdict still counts assays, and never averages their values.

**Saturation mutagenesis.** The repeat experiments of a locus are read over the element's bases and
listed per experiment (bases measured, bases functional), with the count of bases whose functional call
differs between experiments. The verdict stays on the primary experiment, and this is stated. For SORT1
this puts the reverse-orientation experiment, SORT1-flip, in the row.

**Acceptance (tests/test_measured_aggregation.py).** (1) An element matched by five tiles, one at +3.0
and four at +0.2 in K562, is silent in K562, its lentiMPRA verdict is `disagrees`, all five tiles are
listed and `activity_max_descriptive` is 3.0. (2) An element matched by two tiles at +1.5 and +0.5 is
conflicting in K562, its verdict is `reporter_tiles_conflict`, neither counter moves, and `basis:` names
both values. (3) A one-tile element reads exactly as before.

**What may move, stated before the build.** Only the 237 cell readings where tiles disagree can change
label: at most 224 to conflicting and 7 to silent, while 6 stay active. So the lentiMPRA `agrees` of
3,361 can only fall, `disagrees` of 14,508 can only rise, and no element can move from disagreeing to
agreeing. CRISPRi and VISTA counts, the 19,072 experimental element blocks and 212 rule blocks, and every
compiled rule stay the same. Programme text changes only on `basis:` lines and the header's agree and
disagree counts. The calibration table's lentiMPRA arm loses the elements that become conflicting. Any
other change is a fault of the build.

## The result: no reporter label now rests on one tile, and 130 elements the maximum had called active are not (R6, 2026-09-28, later)

Built as registered (`attribution/measured.py` `reporter_block`, `reporter_label`, `reporter_counts`;
tests `tests/test_measured_aggregation.py`). The layer and the calibration were rebuilt in a clean
worktree at `49aca5d`, with only this lane's files changed, through `save_result` with complete
manifests. 0 requests.
The results are saved beside the old ones rather than over them, as
`data/results/measured_layer_genome_r6.json` and `data/results/confidence_calibration_genome_r6.json`,
because the removal guard protects lines of the old files committed in the last two days. The old
names still hold the pre-R6 figures, and the next rebuild after 2026-09-30 can overwrite them with
`--result-name` left at its default.

**What moved, old -> new, and nothing else did** (element by element against a rebuild without the
change at the same commit):

| count | before | after |
| --- | ---: | ---: |
| lentiMPRA agrees | 3,361 | 3,231 |
| lentiMPRA disagrees | 14,508 | 14,512 |
| lentiMPRA `reporter_tiles_conflict` (new, counted apart) | — | 126 |
| elements matched by more than one tile | (not recorded) | 1,006 |
| cell readings where tiles disagree / of them conflicting | (not recorded) | 237 / 224 |
| pooled verdict agrees -> `reporter_tiles_conflict` | — | 125 |
| pooled verdict agrees -> disagrees | — | 4 |
| satmut bases read by more than one experiment / functional calls differ | (not read) | 349 / 84 |
| CRISPRi and VISTA agrees/disagrees, element blocks 19,072, rule blocks 212, every rule and every confidence | — | unchanged |

So 130 elements that the maximum had called reporter-active now read otherwise: 126 because their tiles
split evenly, and 4 because a minority of their tiles were active. That is inside the registered bound
(at most 224 conflicting and 7 silent cell readings; no element moved towards agreement). `basis:`
changed on the 1,006 multi-tile elements, and on chr21 on 12 lines plus the header count (44 -> 43
agreeing). The committed `noncoding_chr21.bio` takes those 13 lines only; it stays held as its own
header says.

**The calibration table.** The lentiMPRA arm's narrow denominator loses the 126 conflicting elements
(17,869 -> 17,743), the observed rate falls from 0.1881 to 0.1821, and ECE rises from 0.0905 to 0.0963.
Its wide denominator is now a different set of rows and is reported (17,869, observed 0.1808): a
conflict counts there as not agreeing. The level stays withheld, as before, because this assay's
agreement is not the compiled claim. Across all tables, 15 bands are populated and 1 has the stated
level inside its interval (it was 12 and 2).

**Not from this change.** The rebuilt `measured_layer_genome.json` also reports 928,094 weak compiled
facts where the stored copy said 812,921. That comes from `96bc5e5` (compiled predictions no longer
state a confidence), since the stored file had not been rebuilt after it. R6 moved no confidence.
(The rebuilt file is `measured_layer_genome_r6.json`; see above.)

**What R6 did not do.** No per-tile uncertainty exists in the ENCODE element files, so none is carried.
The row names this rather than inventing one. A8, the strength of a compiled experimental rule, is still
the largest training |effect| within one element, gene and cell. The pairs sit in the row, but the
compiler reads the maximum. That belongs to the compiler's owner. The saturation-mutagenesis verdict is
still read from the primary experiment; the repeats, SORT1-flip included, are listed beside it with
their disagreement counted.

## Pre-registration: a bridge from compiled annotation to dynamics, and what it will not claim (R3, 2026-09-28)

Registered before any of it was built; the constants are in `genomeos/attribution/bridge.py`.

**The problem.** A compiled program is executable annotation. Its rule `strength` is an
observation's magnitude in that observation's own unit (|log2 fold change| clipped at 1 for an
AlphaGenome deletion, |fractional change| for a CRISPRi screen), its target genes are stubs with no
transcription parameters, and the element it names as a source is not a species the runtime
holds. `runtime/grn.py` read a missing source as 0.0 and a missing `max_rate` as 0.0 without a
word, so a compiled rule ran and moved nothing.

**The mapping.** One observation: removing element E changes gene G's steady-state expression in
cell C by the fold RHO = expression(removed) / expression(intact); RHO = 2^x for a predicted log2
fold change x read from the rule's evidence note (a strength of exactly 1 with no note is
censored), RHO = 1 + f for a CRISPRi EffectSize f. One model: E's state is its presence (1 intact,
0 removed, held by clamp, never defaulted); G must declare `basal_rate` b > 0 and `max_rate` V > 0
(a.u./h); Hill threshold 1 and coefficient 2 are declared assumptions, so H(intact) = h = 0.5. For
one mechanism alone on its gene in its cell, activates: s = b(1/RHO - 1)/(V h); inhibits:
s = (b + V)(1 - 1/RHO)/(V h). A fitted s outside (0, 1] means the declared parameters cannot
produce the observed response, and is reported, not clipped.

**What it will not claim.** A deletion log2 fold change is an observed effect, not a rate constant,
not an affinity, not a dose response, and says nothing about any other perturbation. The bridge
fits one dimensionless strength so that the model reproduces that one response under parameters
the caller declares; anything else the simulation says is the model's.

**The combination rules, named.** The runtime combines a gene's activators as a mean of
s_i H(x_i) and its inhibitors as a product. Under the mean, adding an activator that is low where
another is high halves the drive there, so adding an activator can lower expression (lane-sign,
`b1f3405`, `5cbce26`: correcting Tbxt to activate SOX17 removed mesoderm from the gastrulation
model). The bridge states this as an assumption and does not change it; a change would be its own
registration, with the gastrulation and design tests as dependents. Because a gene with two or
more regulating mechanisms in one context cannot be parameterised from single-deletion
observations under these rules, the bridge reports it as not identifiable.

**The diagnostic.** Every reason a mechanism does not become a parameter is reported by name:
regulator state missing, gene parameter missing, observation missing, observation censored, not
identifiable, conflicting observations, response out of range. The runtime reports a regulator
with no state and a regulated gene with no `max_rate` on every run (`UnresolvedModel` in strict
mode, a warning and a list on the trajectory otherwise); the values it integrates do not change.
A source renamed `<id>@zero`, network_experiment's edge knockout, is a declared zero and is not
unresolved.

**One mechanism, one parameter.** A mechanism is (element without `_measured`, gene, cell). When
a predicted and a measured rule describe it, the measured one is simulated and the predicted one
is recorded as superseded; another citation of the same kind with the same observation is the same
parameter; with a different observation it is a conflict, reported and not averaged.

**The double-counting audit's question, asked before it is answered.** In the 24 compiled
programs, how many (element, gene, cell) mechanisms carry both a predicted and a measured rule that
one run in that cell would integrate side by side, and how many (gene, cell) pairs have more than
one active regulatory rule, so that no single-deletion observation identifies a strength under the
mean rule?

**Acceptance, fixed now.** (1) End to end: a program compiled from one element whose deletion is
observed at log2 fold change -1 in K562, parameterised with declared b and V, reproduces
expression(removed)/expression(intact) = 0.5 within 1% with the element clamped at 1 and at 0.
(2) A run whose regulator has no state raises `UnresolvedModel` naming the missing state in strict
mode and records it otherwise; an `@zero` source raises nothing. (3) A second citation of the same
mechanism, a duplicate rule citing another source or a predicted rule beside a measured one,
leaves the simulated strength unchanged. Negatives are written first.

**Amendment, same day, before the build (lane-assay's census item A8).** A compiled experimental
rule's strength is the largest |EffectSize| among the regulated pairs of one element, gene and cell
(`measured.rule_links`), a maximum used as a parameter whenever there is more than one pair. The
compiler will write each link's pair count on the measured element's evidence note (`<gene> in
<cell> from <n> pairs`). The bridge maps a link that rests on one pair as RHO = 1 + f and reports a
link that rests on several as `observation_summarised`; it never fits to the maximum. No summary
rule (median, mean, a meta-analysis) is chosen here; choosing one is its own registration. The
audit will also count the experimental links that rest on more than one pair.

## Pre-registration: R8's coupling pretest, and a negative closes R8 (R8, 2026-09-28)

Review item R8 allows a joint inference engine only after a pretest shows that a coupling term carries
information that independent per-element scoring does not (genomeos-8a's brief,
`.claude/briefs/joint-annotation-search.md`). This section fixes that pretest before any number is
produced. The code is `genomeos/attribution/joint_pretest.py`, the run `scripts/joint_pretest.py`, the
result `joint_pretest`. No search is built in this lane, whatever the outcome.

**Data and independence (R5).** The CRISPRi benchmark's training file only, through
`measured.load_crispri(split="training")` and `development_only()`; the held-out file is never opened,
not for fitting and not for looking (a test spies on the file opener). lane-split's `split_overlap()`
found 249 held-out pairs sharing bases with training intervals; that is irrelevant here, because the
held-out file is not read. Folds are leave-one-chromosome-out within the training split. The training
file is K562 only (10,356 valid pairs), so per-cell reporting (R1) is K562; any other cell would be
counted and computed only at 20 positives or more.

**Outcomes (R2).** A positive is a significant decrease; a negative is a non-significant pair with
PowerAtEffectSize20 at or above 0.8. Significant increases and underpowered or power-unknown nulls are
excluded from the metric and counted apart.

**Scores.** From the stored per-element response cache (`targets.ElementResponses`), 0 new requests:

- *Independent* (the baseline): a pair's predicted drop `d`, the largest predicted fall of the measured
  gene on K562's own track over the cached elements overlapping the tested interval, floored at 0
  (`crispri.deletion_drop`). A pair no cached element answers is missing and never imputed.
- *Element competition*: `d / (sum of d over every gene the screen tested against the same interval in
  the same cell + 0.1)`.
- *Gene budget*: `d / (sum of d over every interval the screen tested against that gene in that cell + 0.1)`.

The floor 0.1 is the scorer's own `MIN_EFFECT` bar, not fitted. Partners are read by their predictions
whatever their measured outcome; their labels are never read. Each coupling variant enters a logistic
model on `[d, share]` fitted on the other chromosomes; the independent score is `d` itself. The control
is the self-only share `d / (d + 0.1)` in the same model: the share's nonlinearity with no partner.

**Metric.** AUPRC on the out-of-fold scores pooled over the folds. The gain of each variant over the
independent score carries a 97.5% two-sided interval (Bonferroni over two variants) from 2,000
resamples of connected components of the perturbation-gene graph (seed 0); the same interval over
chromosomes is reported beside it. Secondary, never deciding: the same comparison restricted to pairs
with at least one scored partner.

**Pass rule.** A variant passes if the lower end of its component interval lies above 0 and its point
gain exceeds the self-share control's. R8 proceeds if at least one variant passes.

**Falsifier.** A pass is void if the self-share control gains at least as much: the gain would then be a
reshaping of `d`, not coupling. A pass whose chromosome interval includes 0 is reported as weak and
names the chromosome that carries it.

**Readings.** Neither passes: a negative. Coupling carries no information about CRISPRi targets beyond
independent per-element scoring on the data held, so a joint search over this score cannot beat
per-block scoring; **R8 closes and no search engine is built**. One passes: that term is kept, frozen
as registered, and R8 proceeds to its build (the synthetic two-change case first). Both pass: both are
kept. Void: reported as a negative with the control's gain beside it; R8 closes.

## The result: neither coupling term beats independent deletion scoring, so R8 closes and no search is built (R8, 2026-09-28, later)

Run once as registered above (`scripts/joint_pretest.py`, result `joint_pretest`, 0 AlphaGenome
requests, 84 s, peak resident memory 5.4 GB with one cache reader and one chromosome at a time).
Training split only; the held-out file was never opened.

**Negatives first.** On 3,713 scored and labelled K562 pairs (423 significant decreases, 3,290
well-powered nulls; 534 components, 23 chromosomes), out-of-fold AUPRC:

| score | AUPRC | gain over independent | 97.5% interval, components | 97.5% interval, chromosomes |
| --- | --- | --- | --- | --- |
| independent (predicted drop) | 0.766 | | | |
| element competition | 0.754 | −0.013 | −0.024 to +0.001 | −0.022 to −0.000 |
| gene budget | 0.747 | −0.020 | −0.027 to −0.011 | −0.029 to −0.008 |
| self-share control | 0.758 | −0.008 | −0.015 to −0.003 | −0.013 to −0.003 |

Neither variant's interval lies above 0, and neither beats the control's point gain: the registered
reading is **neither passes**. Gene budget is worse than the independent score with an interval that
excludes 0; element competition is worse at the point and its component interval just reaches 0.
Restricted to pairs that have at least one scored partner, where the coupling is not trivial, both
lose more: element competition 0.678 against 0.701 (2,444 pairs, 235 positives, −0.024, −0.045 to
−0.002), gene budget 0.676 against 0.707 (2,894 pairs, 270 positives, −0.031, −0.040 to −0.021). The
share read alone, with no fit, is no better: 0.758 for element competition and 0.705 for gene budget.

**Counted apart.** Of 10,356 training pairs (all K562; no other cell is in the training file, so the
per-cell result is the pooled one), 1,119 have no cached element overlapping the tested interval and
3,410 have an overlapping element that did not score the gene in K562: 4,529 missing, never imputed.
Of the 5,827 scored, 46 significant increases and 2,068 underpowered nulls are kept out of the metric.
The 0.766 is therefore an AUPRC on the pairs the sweep answers, which excludes genes outside its
window; it is not comparable with the benchmark's all-pairs figures.

**Reading.** On the data held, a gene's rank among the genes tested against the same perturbation, and
an element's share of its gene's predicted regulators, carry no information about which pairs a
CRISPRi screen calls regulated beyond each pair's own predicted drop; both lower precision. A joint
search over this score cannot beat per-block scoring, so **R8 closes as a negative and no search
engine is built**. What this does not say: it tests coupling of the deletion score as held, in K562
alone, on the pairs the sweep answers; it does not test coupling of a different score or in a cell
with its own training data. lane-split's overlap check (249 held-out pairs share bases with training
intervals) did not bear on this run, because the held-out file was not read.

## The result: no compiled rule can be simulated as compiled, one gene-cell pair in three has several mechanisms and cannot be parameterised from single removals, and 39 relations were carried twice (R3, 2026-09-28, later)

Built as registered in `06b7557` and amended in `b986238`. `runtime/grn.py` now names its
combination rules (`ACTIVATOR_COMBINATION = "mean"`, `INHIBITOR_COMBINATION = "product"`) and its
declared zero (`EXPLICIT_ZERO_SUFFIX = "@zero"`), and `NetworkRuntime.diagnose` reports every
regulator with no state and every regulated gene with no `max_rate`: on the trajectory's
`unresolved` list and as an `UnresolvedModelWarning`, or as `UnresolvedModel` when the runtime is
strict. The integrated values are unchanged, and no test in the suite raised the warning, so no
existing caller depended on an undeclared zero; network_experiment's `@zero` edge knockout is
recognised as declared. `attribution/bridge.py` reads one observation per mechanism (2^x from a
predicted rule's effect note, 1 + f from a CRISPRi rule resting on one pair), picks one
representation by precedence, and fits one strength per mechanism under the registered mapping.
The compiler writes each measured link's pair count on its measured element and states in every
program's header that it is executable annotation and not a simulation. Tests:
`tests/test_grn_bridge.py`, 17, written first; all three acceptance tests pass, including the
end-to-end fixture (log2 fold changes -1, -0.4 and +0.7, and a CRISPRi f = -0.3, each reproduced
within 1%) and the fitted strength differing from the compiled |log2 fold change| in every case.

**The audit's answer (`scripts/bridge_audit.py`, `data/results/bridge_audit.json`, 24 programs,
every cell context a rule names).** Of 440,589 compiled regulatory rules, **0 can be simulated as
compiled**: every target gene is a stub. They form 440,550 mechanisms over 237,613 (gene, cell)
pairs. **76,469 of those pairs (32%) have more than one mechanism**, 279,406 mechanisms between them,
and cannot be parameterised from single-removal observations under the mean rule; the other
161,144 pairs have one mechanism whose observation reads cleanly and lack only declared rates.
**39 relations carry both a predicted and a measured rule in the same cell**, all in K562, which a
plain run would integrate side by side, averaging a |log2 fold change| with a |fractional change|
(SMIM1: 0.848 beside 0.472); the bridge supersedes the predicted one where it can act (13) and the
other 26 sit on genes that are not identifiable anyway. All 212 experimental links rest on one
CRISPRi pair, so the maximum lane-assay's census found (A8) never summarised more than one pair in
these data; no predicted rule is censored, because every one now carries its effect note.

**Negatives.** The bridge fits one strength to one observation and nothing else: it cannot say what
a partial deletion, a dose or a second perturbation would do, and a gene with two mechanisms is
refused rather than solved. The mean rule stays, so adding an activator can still lower
expression; changing it is its own registration. The committed chr21 copy holds its 5,176 rule lines
at their earlier text under the removal guard until 2026-09-30; the 24 git-ignored programs carry
everything.

## Pre-registration: five axes instead of one label, class kept as a derived summary (R7, 2026-09-28)

Review item R7. The census (`2ba726b`, `data/results/ontology_census.json`) found every one of the
459,449 compiled elements saying `class: enhancer`, the constant written by `compile.py`, including
the 156,925 whose rule inhibits its gene; 101,011 predicted elements at least half interspersed repeat
with nothing saying so; and region roles that join a tier, a budget label and a copy flag in one
string, 9,499 of which read an absence of constraint as a fossil, a dead frame or neutral sequence.
The 15,252 regulatory blocks never had their repeat coverage read, because the block classifier names
a block `regulatory` before it reads RepeatMasker.

**The axes.** An `element` and a `region` may state five keys, each with a closed vocabulary held in
`genomeos.lang.grammar.AXES` (every axis has `unknown`; a value may take a qualifier after `/`):

- `origin`: `unique`, `repeat_derived/CLASS` (>= 0.5 of the interval interspersed repeat, CLASS the
  largest), `partly_repeat_derived/CLASS` (> 0, < 0.5), `satellite`, `tandem_repeat`,
  `segmental_duplication` (>= 0.5 duplicated), `assembly_gap`.
- `molecular_role`: `promoter_like`, `enhancer_like`, `insulator_like`, `open_chromatin`, `silencer`,
  `competing_promoter`, `structural`, `coding_candidate`.
- `activity`: `activates_target`, `represses_target`, `no_effect_measured`, `active_in_reporter`,
  `inactive_in_reporter`.
- `target_relation`: `predicted_deletion_target`, `nearest_tss_in_domain`,
  `measured_perturbation_target`, `tested_no_effect`, `unassigned`.
- `evidence_status`: `curated_annotation`, `registry_biochemical`, `predicted_model`, `measured`,
  `measured_negative`, `conflicting`, `under_selection` (>= 5% of bases at phyloP >= 2.27),
  `selection_weak` (3-5%), `selection_not_detected` (< 3%), `selection_not_measured`.

`,` joins values that all hold (overlapping roles: a CTCF-bound dELS is `enhancer_like,
insulator_like`); `|` joins unresolved alternatives, one of which holds and none of which is chosen.
The parser refuses a value outside the vocabulary and keeps the axes in the entity's `attrs` as
`{"ontology": {axis: [[alternative, ...], ...]}}`, which `Module.to_dict` and `from_dict` carry.

**Constraint.** Constraint is evidence of selection and is stated only on `evidence_status`. No axis
has a value meaning "no function"; an unconstrained fossil or unique block says `molecular_role:
unknown` and `evidence_status: ..., selection_not_detected`. The predicted elements carry no
constraint field in their runs, so they say `selection_not_measured` (nothing is fetched to fill it).

**What each block says.** A predicted element: origin from RepeatMasker over its interval; role from
the registry entry of its own id (`CCRE_ROLE`, plus `insulator_like` when CTCF-bound); activity from
the rule's direction; target relation `predicted_deletion_target`, plus `nearest_tss_in_domain` when
the run's verdict says they agree; status `registry_biochemical, predicted_model,
selection_not_measured`. When the prediction is an increase on deletion, the role gains the
unresolved group `silencer|insulator_like|competing_promoter|unknown` and the activity says
`represses_target`: a direction of effect never becomes a silencer label. A `_measured` block keeps
the sequence facts (origin, role) and states its own activity, relation and status from the assays:
CRISPRi decrease `activates_target`, increase `represses_target` (with the same group), well-powered
nulls only `no_effect_measured` and `tested_no_effect`; lentiMPRA or VISTA active in any cell
`active_in_reporter`, inactive everywhere `inactive_in_reporter`, a lentiMPRA tie the two as
alternatives; `measured`, `measured_negative` when a gene was measured unchanged, `conflicting` when
the agreement verdict is disagrees. A region: origin from its sequence class, or for a block whose
class is not a repeat class from RepeatMasker over the block, plus `segmental_duplication` from the
copy flag; role from the registry classes it contains for a `regulatory` block, `promoter_like|unknown`
for a CpG-island block, `structural` for a centromere, `structural|unknown` for a satellite or tandem
array, `coding_candidate|unknown` for a long ORF, `unknown` otherwise; activity `unknown`; relation
`unassigned`; status the annotation it rests on plus the selection reading.

**Backward compatibility.** Two choices, both the least disruptive the census allows. `class:` stays,
now derived from the registry role (`ROLE_CLASS`, the first role) and labelled as that summary in the
grammar; since the census shows every predicted element is pELS or dELS, the derived class is
`enhancer` for all of them and no class count moves, while an element outside the registry would now
say `unknown` rather than inherit `enhancer`. No runtime, evidence-explorer or web path reads `class:`
from a compiled program (census, consumers), so nothing downstream changes. A region's `role:` stays
verbatim: `Module.unknowns` counts `role: unknown`, and every program's `# test: unknowns ==` line and
the organiser read it, so it is kept as the budget's tier summary, labelled so in the grammar, with
the five axes beside it as the authoritative reading. Rewriting the budget's labels is a separate
step (24 stored budget results).

**Acceptance (the review's, as tests).** (1) A repeat-derived regulatory element keeps both: an
element over a LINE with a dELS registry entry compiles to `origin: repeat_derived/LINE` and
`molecular_role: enhancer_like`, and a regulatory region over repeats keeps `repeat_derived` beside
its registry roles. (2) A predicted expression increase on deletion does not force a silencer label:
the element's `class:` is the registry summary, its activity `represses_target`, and `silencer`
appears only inside the unresolved group. (3) Unresolved alternatives survive serialisation: compile
-> parse -> `Module.to_dict` -> JSON -> `from_dict` returns the same groups, and an out-of-vocabulary
value is a parse error.

**Expected counts over the 24 programs.** Unchanged: 440,377 + 19,072 `class: enhancer`; 440,589
rules; 1,098 unknowns; every region `role:`. New: 459,449 elements and 26,806 regions each state five
axes. Predicted elements: `repeat_derived` 101,011, `partly_repeat_derived` 81,058, `unique` 258,308;
`insulator_like` beside `enhancer_like` on 193,027; `represses_target` with the unresolved group on
156,925 and `silencer` alone on 0. Regions with `molecular_role` meaning no function: 0. The count of
repeat-derived `regulatory` blocks is not known in advance (never read) and is reported as found.
Known costs: programs grow by five lines per element and region; the committed chr21 copy's changed
lines may be held by the removal guard until 2026-09-30.

## The result: every compiled element and region states five axes, no count the programs test moved, and a repression is never a silencer (R7, 2026-09-28, later)

Built as registered (`genomeos.lang.grammar.AXES`, the parser's `_ontology`, `compile.py`'s
`element_axes`, `measured_axes` and `region_axes`; `tests/test_ontology.py`, 7 tests). The 24
programs in `data/knowledge/compiled` were recompiled and every line that is not an axis line or the
new header paragraph is byte-identical to the programs before (checked per file), so nothing
registered as unchanged moved: 440,377 + 19,072 `class: enhancer`, 440,589 rules, 1,098 unknowns,
every region `role:`. Counts after (`data/results/ontology_census_after.json`, beside the census):

- **Predicted elements (440,377, all five axes):** origin `repeat_derived` 101,011, `partly_repeat_derived`
  81,058, `unique` 258,308, exactly as registered; role `enhancer_like` 440,377 with `insulator_like`
  beside it on 193,027; `represses_target` 156,925, every one with the role left as
  `silencer|insulator_like|competing_promoter|unknown`, and `silencer` alone on 0; target relation
  `predicted_deletion_target` 440,377, `nearest_tss_in_domain` beside it on 268,940; status
  `selection_not_measured` on all, because the element runs carry no constraint.
- **Measured blocks (19,072):** activity `inactive_in_reporter` 14,539, `active_in_reporter` 3,296, the
  two as alternatives (a lentiMPRA tie) 126, `no_effect_measured` 976, `activates_target` 181,
  `represses_target` 43 (the role left open on all 43), `unknown` 243; relation
  `measured_perturbation_target` 145, `tested_no_effect` 1,284, `unassigned` 17,643; status
  `conflicting` 14,537, `measured_negative` 1,406.
- **Regions (26,806):** origin `repeat_derived` 16,192, `partly_repeat_derived` 9,958, `unique` 418,
  `segmental_duplication` beside another origin on 2,894, `assembly_gap` 118, `satellite` 113,
  `tandem_repeat` 7. Of the 15,536 regions with a registry or CpG role, **8,709 are at least half
  interspersed repeat and 6,534 partly**: the repeat-derived regulatory blocks the classifier had
  hidden, the count the registration could not predict. Role `unknown` 11,061; no region carries a
  value meaning no function; `selection_not_detected` 21,283, `selection_weak` 2,962,
  `under_selection` 2,410, `selection_not_measured` 33.

**Negatives and limits.** The derived `class:` moved nothing because every predicted element is pELS
or dELS in the registry: the old constant was right for the registry and wrong in kind, and the
change is that it is now read, not assumed. A region's `role:` keeps the budget's strings, including
"no evidence of function, best guess neutral" on 2,546 regions; the axes beside it say `unknown` and
`selection_not_detected`, and rewriting the budget's labels (24 stored results) is its own step. No
element has a per-element constraint reading, so the selection axis is unmeasured for all predicted
elements; filling it needs phyloP over 440,377 intervals. Segmental duplication is read for regions
only. The committed chr21 copy received the axis lines and the header paragraph as pure additions
(29,108 lines, none removed), leaving its guard-held rule lines as they were until 2026-09-30.

## Pre-registration: every budget reader through one reader, names only (R7 follow-up, consumers, 2026-09-28)

The budget's R7 restatement (`budget_axes_<chrom>`, commits b6e6d17, bbce3d8, 308484a) renames two
tiers (fossil to `repeat_unconstrained`, neutral to `unconstrained_unknown`) and eight labels, and
adds `evidence_status`, `origin` and `legacy_tier`. Every consumer still read `budget_<chrom>`. This
step moves them all onto one reader, `budget.read_axes(chrom)`, registered in
`budget.CONSUMERS_REGISTERED` before any consumer moves.

- **What moves.** The reader returns the axes record, with the measurements the restatement dropped
  (class evidence, phyloP mean and max, conserved elements, cost) filled back from the stored budget
  after a block-by-block agreement check. A consumer whose own result joins on the old key reads
  `legacy_tier`, so organise, compress, human_panel, lexicon, variation, unknown_scoring,
  duplications, element_types and the epigenome keep their stored keys. The CLI, the web budget panel
  and a compiled region's `role:` line show the new tier names and labels.
- **The invariant.** No numeric figure moves: tier counts and base pairs are identical by
  construction, and every number a consumer writes must be identical before and after. Only names
  move, and only where a person reads them. Run times and a byte counter of files read are exempt.
- **The proof.** Each consumer that writes a result is run before and after into scratch, never over
  a committed result, and every numeric leaf is compared, tier keys mapped through `LEGACY_TIER`.
  Any other difference is a failure of the move and is reported with its cause.
- **Negatives first.** Seven modules the census counted as budget readers read a tier through
  `organise.blocks` or lexicon's context map, not the budget file; they change through organise and
  keep the legacy key. The censuses (`budget_axes_census.py`, `ontology_census.py`) and the
  restatement itself keep reading `budget_<chrom>`, because they describe it.
- **Not in this step.** Regenerating `unknown_<chrom>` so the classifier's origin is stored, and
  compile.py reading that origin instead of re-deriving it, go ahead only if cheap and only if the
  comparison shows no numeric change; otherwise after 2026-09-30.

## The result: every budget reader now goes through one reader, and across 24 chromosomes and 19 consumer outputs no number moved (R7 follow-up, consumers, 2026-09-28, later)

Built as registered above. `budget.read_axes(chrom)` returns the `budget_axes_<chrom>` record and
fills back the stored measurements after its agreement check, which holds on all 24 chromosomes
(26,806 blocks). `read_axes_genome_wide()` returns the R7 genome-wide sum and `tier_name()` gives a
person the R7 name of any tier key.

**Negatives first.** Nothing had to be explained away: 0 numeric differences in every comparison.
The only non-exempt differences are names and one input list. Seven modules the census counted read a
tier through `organise.blocks` or lexicon's context map, so they moved with organise and needed no
edit. The consumers' stored results keep the legacy keys `fossil` and `neutral`, as do result names
such as `epigenome_fossil` and lexicon's `unknown_fossil` context; renaming those keys changes the
stored schemas of the consumers and every reader of them, not the budget, and is not done here. variation.py's
case reading still says "a best guess" (a reading of the two axes, not of the budget tier): its two
lines were added on 2026-09-27 and are held by the removal guard until 2026-09-30.

| Consumer (run before and after, into scratch) | Numeric leaves compared | Differences |
| --- | ---: | ---: |
| organise.organise, 24 chromosomes | 12,862 | 0 (23 run times exempt) |
| organise.blocks, 24 chromosomes (also what syntax_tiling, measurability, motif_transfer, unknown_coverage, syntax_blocks_matched, constrained_unknown_targets read) | 265,301 | 0 |
| compress.tier_blocks, 24 | 53,612 | 0 |
| unknown_scoring.unknown_blocks, 24 | 131,389 | 0 |
| element_types block tiers, 24 | 53,612 | 0 |
| duplications.run_and_save (not saved), 24 | 108,137 | 0 |
| variation.build, chr21, chr22, chrY, human axis stubbed (no request) | 22,223 | 0 |
| human_panel.build, chr21, human axis stubbed | 3,621,823 | 0 (1 run time) |
| lexicon context map, chr21 | 1,797 | 0 |
| panel_background.rebuild, chr21 | 2 | 0 |
| epigenome summary, fossil and alu actions, chr21 | 3,185 | 0 |
| budget distil, legacy and R7 | 610 | 0 |
| web budget card | 347 | 0 |
| CLI budget (chr21 and genome), duplications, variation, organise text | 491 numbers | 0 lines whose numbers differ |
| compiled programs, 24 chromosomes | 9,896,402 numbers | 0 lines whose numbers differ |

`organise.inputs` now names `budget_axes_<chrom>` beside `budget_<chrom>`, so a manifest built from it
lists both files the reader opens. The compiled programs changed on role lines only: 397 of 446 regions
on chr21 and 23,708 of 26,806 genome-wide, each now `role: repeat_unconstrained, repeat-derived sequence, no
sign of selection, which is not a sign of no function` and the like; no role line in any of the 24
programs says neutral, fossil, dead or best guess. The committed chr21 copy received exactly those 397
role lines (all from 2026-09-12, none guard-held); its held rule lines are untouched. The CLI and the web
budget card print the R7 tier names, and their notes say constraint is evidence of selection only.

**Left for after 2026-09-30.** Regenerating `unknown_<chrom>` so the classifier's origin is stored,
and compile.py reading that origin instead of re-deriving it: the comparison shows no numeric change,
but the step is not cheap (24 classifier runs and 24 rewritten results) and it changes origin text on
11 blocks whose interspersed coverage sits just under one half (commit 308484a), so it is not names only.

## Pre-registration: milestone 1.3's second clause against a control matched on length and gene density (2026-09-28)

**What is compared.** Clause 2 of milestone 1.3 (the constrained-unknown blocks attributed to a gene
and a tissue) is held as not met because the real unknown names a coding gene far below length-matched
random windows: in d717b28, 29.8% of 531 blocks against 67.0% of windows, -37.2 points. The audit then
found the neutral tier -36.0 points below the same control, and the real unknown +4.6 points above
neutral inside length deciles, which names local gene density as the likely confounder. This test asks
the same question once more against a control that accounts for it. Script:
`scripts/clause2_matched_control.py`, whose `PRE_REGISTRATION` holds the constants below; tests in
`tests/test_clause2_matched_control.py`. No matched figure had been computed on any chromosome when
this was committed.

**Admissibility.** The control was chosen after the audit named the confounder. That is admissible
because the audit only named the variable; the matching variable, its bins, the estimator, the
interval, the pass rule and the falsifier are fixed here before the comparison runs, and the same
windows are asked of the neutral tier so the control is not tuned to the target.

**Targets.** The real unknown: constrained-unknown blocks with copies out (882 blocks, 531 carrying a
scored element), all 24 chromosomes. The neutral tier (2,632 blocks, 1,181 carrying) is the secondary
target set.

**Control.** For each block, windows of the block's exact length drawn as d717b28 drew them (uniform
inside the span of scored-element midpoints, rejected on overlap with any organiser block, same seed
per chromosome), and also rejected unless the count of GENCODE protein-coding TSSs within 524,288 bp
of the window's midpoint falls in the block's own decile. Deciles are nearest-rank cuts of each target
set's own block counts pooled genome-wide. 50 accepted windows per block, at most 20,000 tries; a block
with no accepted window leaves the comparison as undrawable.

**Question and primary.** Does a block or window carrying a scored element carry one whose deletion
names a protein-coding gene. Primary statistic: for the real unknown, the mean over compared blocks of
(block yes, 0 or 1) minus (share of that block's carrying windows saying yes), in points.

**Interval.** 95% percentile bootstrap over blocks, 10,000 resamples, seed 20260913. A bootstrap over
the 24 chromosomes is reported beside it as a sensitivity for correlation along a chromosome, with the
pooled rates as d717b28 printed them.

**Secondary and sensitivities.** Secondary: the neutral tier under the same matching; real unknown
minus neutral on the matched difference; the question "moves any gene" on both sets. Sensitivities:
matching on distance to the nearest coding TSS instead of the count; rejecting windows only on overlap
with the target set itself. A balance check (covariate and elements per block against per window) is
reported, not gated.

**Gate.** Before any matched figure is read, the same function with the density test off and
d717b28's 4,000-try cap must reproduce d717b28's pooled counts to the digit: real unknown 44,100
windows drawn, 31,676 carrying, 21,217 naming a coding gene, 531 and 158 blocks; neutral 131,600,
88,317, 57,212, 1,181 and 340. If it does not, no matched figure is reported.

**Pass rule.** Clause 2, as a labelled lead with measurement still outstanding, passes on this
instrument only if the primary's 95% interval over blocks lies wholly above 0. An interval covering 0
reads as "at chance against a matched control" and the clause stays not met; an interval wholly below
0 keeps it not met, and "below chance" if the point is also below -5 points. The clause's third part
(scored against measurement) is untouched by any outcome.

**Falsifier.** The audit's account (most of the gap is gene density) is falsified if the matched
difference keeps more than half of the unmatched gap, that is, falls below -18.62 points (half of
-37.23).

**Registered, not run.** A measured arm: the designed experiment (284,001 oligos over 878 blocks)
tiles each block and a density-matched window set, so an MPRA readout would give measured activity per
block against its own control. Nothing measured exists for it (0.45% of unknown space has any assay),
and it is not run.

**Cost.** 0 AlphaGenome requests; the response cache is never opened.


> **Annotation (2026-09-28, lane-design, after a statistical review; the text above is unchanged).**
> Two readings of this registration need narrowing. (a) The question it fixes, "does a block or window
> carrying a scored element carry one whose deletion names a protein-coding gene", is a **target-naming
> frequency of the deletion model**. A difference in it is not a difference in prediction accuracy, and
> on its own it shows neither that the model failed nor that the blocks' sequence lacks regulatory
> function; no experimental outcome enters it. (b) "0.45% of unknown space has any assay", in
> *Registered, not run*, is a coverage figure and needs its rule and denominator like every other: the
> corresponding block-level counts are 59 of 882 tier blocks (6.7%) or 59 of the 531 carrying a scored
> element (11.1%) at reciprocal overlap 0.5, and 72 or 10 blocks if that rule is set to 0.25 or 0.75.
> The full denominator table is in the lane-design result section below.

## The result: in HCT116 the frozen deletion model gains +0.022 with an interval across zero, so the CRISPRi result passes its second cell type but is not replicated there (2026-09-28, lane-hct116)

The registration is `crispri.PREREGISTERED_PUBLISHED['second_cell_type']` (fe61c36, 2026-09-27,
above): the same code path with HCT116 added to the cells kept per gene, weights unchanged; it
passes if the gain on the 363 covered HCT116 pairs is above zero, and counts as a replication only
if the chromosome-bootstrap interval excludes zero. The owner approved up to 760 requests on
2026-09-28. The fetch and scoring code (`scripts/crispri_hct116.py`) was committed in 6257b95
before any HCT116 value was read. The result is `data/results/crispri_published.json`, key
`second_cell_type_hct116`, inserted beside `second_cell_type` without changing any existing line.

**AlphaGenome requests: 705,** the costed figure exactly. All 705 were answered; there were no quota
refusals and no retries. Each answer is stored in the sweep's cache format under its own root
(`data/knowledge/alphagenome/elements_hct116`, local), with a ledger entry written before every
request was sent.

### The registered verdict first

| covered HCT116 pairs (363; 34 regulated; 5 chromosomes) | AUPRC (average precision) |
|---|---|
| distance | 0.334 |
| activity + distance | 0.517 |
| activity + distance + deletion | 0.539 |
| **gain** | **+0.022, 95% −0.058 to +0.145** (200 chromosome resamples) |

**Passes, not replicated.** The gain is above zero, so the rule that would have made this "a K562
result" does not fire. The interval includes zero, so the registration forbids calling it a
replication. score_published's per-cell block (all 396 pairs, weighted, the benchmark's estimator)
agrees: 0.495 → 0.509, gain +0.014 (95% −0.042 to +0.126). The frozen weights reproduce the K562
headline to the last digit on the same run: +0.1407 (95% +0.082 to +0.231), 1,744 pairs.

The registration left one choice open, and it was fixed before scoring. It names the pairs, the
sign rule and the bootstrap, but not the AUPRC estimator. The primary figure uses the estimator of
the K562 headline this arm replicates.

### Why this arm could say so little: five genes, and most pairs out of the model's reach

These limits were not stated in the registration. They are reported here and change no verdict.

- **The 34 positives are five genes at five loci:** SSFA2 (10), FAM3C (8), MYC (7), KITLG (5) and
  CCND1 (4). The pairs fall on chr8 (257), chr11 (74), chr2 (34), chr7 (26) and chr12 (5). A
  chromosome bootstrap over five units is a bootstrap over five loci, which is why the interval is
  wide. The registration's "34 regulated" read as more evidence than it is.
- **228 of the 363 covered pairs get no deletion value in any cell,** 11 of the 34 regulated among
  them. Their gene lies outside the scorer's 1 Mb window: the median distance is 2.16 Mb for these
  pairs against 245 kb for the answered ones, and they are mostly the long-range MYC and CCND1
  screens. Only 135 pairs (23 regulated) carry a deletion answer, and 106 (21 regulated) carry a
  predicted fall. That is 37% of covered pairs answered, against 66% in held-out K562 (1,152 of
  1,744). "Covered" in the registration meant "on a deleted element", not "within the model's
  reach". The two had not been told apart.

### The model version: mixed by rule, identical by measurement

Every new request asked for ALL_FOLDS through `alphagenome_adapter.create_client`. The frozen weights
and the `top_target` flag read 6,097 answers from the 2026-09 sweep, which asked for no version.
So the manifest says **mixed** (`result_manifest.model_dependencies[0].mixed`), as the rule requires.

The purchase also measured what the rule could only assume. Each of the 705 elements was answered
twice, once by the sweep with no version requested and once now with ALL_FOLDS requested. Across
54,760 gene-by-cell values on K562, HepG2, GM12878 and IMR-90, **all 54,760 agree to four decimal
places**, and the top predicted target is the same for all 705 elements
(`second_cell_type_hct116.versus_the_sweep`). For these elements, the unrequested sweep answered as
ALL_FOLDS answers now. The flag stays because the rule counts labels, and 705 elements are not the
whole sweep. For this claim, though, the mixed versions do not explain the weak HCT116 gain: the
weights were fitted on the same model's answers.

### What the CRISPRi result may now be called

- Allowed: the K562 sentence of 2026-09-27, unchanged, then: "In a second cell type, HCT116
  (363 pairs, 34 regulated, five genes), the frozen model's gain is +0.02 (95% −0.06 to +0.14):
  positive, not significant, and limited by the 1 Mb window, which reaches only 135 of the pairs."
- Not allowed: "replicated in a second cell type", "generalises across cell types".
- The registration's failure clause did not fire, so the claim need not be withdrawn to K562. It
  remains **replicated in one cell type**. An informative second cell type needs more genes, not
  more requests for these five. It also needs pairs inside 1 Mb, or a longer-range feature
  registered for them.

## Pre-registration: milestone 1.3's second clause against a control matched on length and scored-element count (2026-09-28, lane-elemcount)

**What is compared.** Clause 2 of milestone 1.3 is still held as not met. d717b28 put the real unknown
-37.23 points below length-matched random windows on "names a coding gene"; e1dbcf3/0c8b82d matched
those windows on local coding-gene density as well and the gap stayed -27.25 points (95% over blocks
-30.91 to -23.58, n = 531 blocks), with the neutral tier equally low and every sensitivity agreeing.
That run matched its covariate (4.11 coding TSSs per compared block against 4.05 per window) and
reported exactly one residual imbalance, naming it the last live explanation for the gap besides "the
instrument does not work": a carrying window holds 13.27 scored elements against 6.18 per carrying
block, so a window has more chances to name a gene. This test closes that imbalance. Script:
`scripts/clause2_element_count_control.py`, whose `PRE_REGISTRATION` holds the constants below; tests in
`tests/test_clause2_element_count_control.py`. `scripts/clause2_matched_control.py` is not rewritten;
its draw, gate, estimator and reading bands are imported and reused. No element-count-matched figure had
been computed on any chromosome when this was committed.

**Admissibility, decided before the numbers.** The imbalance is expected by construction, and matching
on it is therefore only half admissible. An organiser block is an intergenic gap: `unknown_blocks`
carves blocks out of the sequence left between annotated gene bodies, so no block contains a gene,
while a window placed uniformly at the same length usually straddles one. The budget's classifier then
sends CpG-island blocks to the `regulatory` tier, so the most promoter-like -- hence most cCRE-dense --
intergenic blocks are removed from `constrained_unknown` before this comparison starts. Scored elements
are ENCODE SCREEN cCREs, which cluster at promoters and inside gene bodies. A constrained-unknown block
is element-poorer than a random window of its length by definition, not by accident, and matching on
element count conditions on a variable the target set's own definition sets. Two consequences are
registered here rather than chosen afterwards.

1. The matched window set is pulled toward the target set itself. Windows as element-poor as a block,
   at the block's length, are disproportionately the same intergenic deserts the blocks are. This
   control can therefore only attenuate a difference, never inflate one, so a difference that stays far
   below 0 is informative and cannot be an element-count artefact.
2. A difference that moves to 0 is weak. It is equally consistent with the element-count artefact and
   with the control having become a copy of the target. On its own it does not restore clause 2: it
   replaces "below chance" with "unresolved by this instrument" and hands the question on.

Because of (2), a second estimator is registered that removes the number-of-trials imbalance
arithmetically instead of by matching, and it is read first if the matched difference moves: the
per-element rate, which does not depend on how many elements anything holds, computed on d717b28's own
unmatched windows.

**Targets.** Unchanged from e1dbcf3. The real unknown: constrained-unknown blocks with copies out (882
blocks, 531 carrying a scored element), all 24 chromosomes. The neutral tier (2,632 blocks, 1,181
carrying) is the secondary target set.

**Control.** For each block, in start order, windows of the block's exact length drawn as d717b28 drew
them (uniform inside the span of scored-element midpoints, rejected on overlap with any organiser
block, one random number per try, same seed 20260913 per chromosome), and additionally rejected unless
the number of scored elements whose midpoint falls inside the window lies in the same decile of
scored-element count as the block's own number. 50 accepted windows per block, at most 20,000 tries,
e1dbcf3's cap unchanged. A block for which no window is accepted inside the cap is counted as
undrawable and leaves the comparison; a block that carries an element but gets no accepted window
carrying one also leaves the comparison. Both counts are reported. Both rules are e1dbcf3's.

**Banding.** Deciles of scored-element count: the nearest-rank 10th to 90th percentiles of the target
set's own per-block element counts, pooled genome-wide over every block of the set including the blocks
that hold none, ties collapsed so there may be fewer than ten bins. This is the same `decile_edges`
e1dbcf3 used for coding-TSS count, applied to the new variable. Each target set is banded on its own
distribution. A block holding no element falls in the lowest bin and can only draw windows holding
none, so it leaves the comparison as it already did.

**Question, statistic, interval.** Unchanged. The question is whether a block or window carrying at
least one scored element carries one whose deletion names a protein-coding gene. The primary is the
per-block matched difference: over blocks that carry an element and have at least one accepted window
carrying one, the mean of (block yes, 0 or 1) minus (share of that block's carrying windows saying
yes), in points. The interval is a 95% percentile bootstrap over blocks, 10,000 resamples,
`default_rng(20260913)`, with a bootstrap over the 24 chromosomes reported beside it and the pooled
rates as d717b28 printed them.

**Secondaries.** The per-element rate described above, on the unmatched windows and on the
element-count-matched ones; the neutral tier under the same matching with its own deciles; real unknown
minus neutral; and the question "moves a gene" on both sets. Sensitivities: matching jointly on the
scored-element decile and e1dbcf3's coding-TSS decile, so both imbalances are closed at once; matching
on the exact element count rather than its decile; and rejecting windows only on overlap with the
target set itself.

**Gate.** Before any element-count-matched figure is read, the same draw function with the element test
switched off and d717b28's 4,000-try cap must reproduce d717b28's pooled counts for both sets to the
digit: real unknown 44,100 windows drawn, 31,676 carrying, 21,217 naming a coding gene, 531 and 158
blocks; neutral 131,600 / 88,317 / 57,212, 1,181 and 340. If it does not, no matched figure is reported.

**Balance check, reported and not gated.** Scored elements per carrying block against per carrying
window -- the quantity being matched, which must come close to equal for the matching to have worked --
and the mean coding-TSS count of the compared blocks against their windows', which is free again here
and may drift back apart.

**Pass rule.** Clause 2 passes on this instrument only if the primary's 95% interval lies wholly above
0: a real-unknown block then names a coding gene more often than sequence of the same length holding as
many scored elements. The clause's third part (scored against measurement) is untouched by any outcome.

**Readings, fixed here.**

- Interval wholly above 0: clause 2 passes as a labelled lead against an element-count-matched control,
  with measurement still outstanding.
- Interval covering 0, that is the difference moving close to 0: the element-count imbalance was the
  gap, the below-chance readings of d717b28 and 0c8b82d were an element-count artefact, and clause 2's
  status must be re-examined -- it stops being "below chance" and becomes "unresolved by this
  instrument", not "met". By the admissibility judgement above this outcome is weak on its own, because
  the matched control is pulled toward the target set; the per-element secondary is read next, and the
  clause stays not met until a measured arm answers it.
- Interval wholly below 0: element count is not the cause either. With length, coding-gene density and
  the number of scored elements all matched and the difference still far below 0, no covariate named so
  far explains the reading, and clause 2 stays not met and below chance on its existing ground.

**Falsifier.** The element-count account of the gap is falsified if the element-count-matched difference
keeps more than half of d717b28's unmatched gap, that is if it is below -18.61 points (half of -37.23).
This is e1dbcf3's threshold, unchanged so that the two controls are read on one scale. The share of the
unmatched gap and the share of the density-matched -27.25 points that element count closes are both
reported.

**The measured arm is not run.** It stays registered for its own lane, and no figure from it is read
here. Nothing measured exists for these blocks: 0.45% of unknown space has any assay.

**Cost.** 0 AlphaGenome requests. The per-element response cache is never opened; the stored all-element
archive is read one chromosome at a time.

## Pre-registration: how independent the node containment result is, and what a fresh CRISPRi set could see (area B, 2026-09-28, lane-nodeindep)

**What is audited.** The node containment claim (area B): on 661 CRISPRi `Regulated=TRUE` pairs the
element and its measured gene share a node 5.89 points more often than under the published
random-boundary control, 95% +3.18 to +8.46 (`data/results/node_containment_measured.json`,
`scripts/node_containment_audit.py` stage 2). lane-split (`0af7e1b`) classified that arm
EVALUATION-ONLY and flagged the ENCODE held-out file as a reused benchmark. EVALUATION-ONLY says a
scorer fitted nothing to the pairs; it does not say the scorer was built without them. This audit
asks the second question of every component of the node caller and reports it under lane-split's
legend: CLEAN, EXPOSED, EVALUATION-ONLY. Script: `scripts/node_independence_audit.py`; tests in
`tests/test_node_independence_audit.py`. No figure of the existing claim is changed by any outcome;
a qualification, if one is owed, is written beside it.

**Checks, fixed here.** (1) The default caller is imported in a fresh interpreter and every
`genomeos` module it pulls in is listed, so "the caller never reads CRISPRi" is a fact about the
import closure. (2) `git log -S` dates each constant of the caller — the 50 kb node floor, the 5 kb
boundary merge, the CTCF-only class filter, the node confidence, the control's seed and draw count —
against the first commit naming the benchmark. (3) The stage 2 containment counts are recomputed on
every covered chromosome from the CRISPRi tables and the caller alone, with the AlphaGenome archive
never opened, and held against the committed result: identical counts show the measured arm does not
depend on the model output stage 1 scores. (4) The 661 pairs are decomposed by benchmark file and
cell. (5) The six rejected node callers are scored on the same pairs. (6) Every CRISPRi set this
project has not scored for containment is counted and its power stated.

**The sensitivity, and how it will be read, before it is run.** The default caller `ctcf_only` was
kept over six orientation callers on 2026-09-14 (`26da99a`) and again on 2026-09-21 (`56e2c50`), on
four measurements: 4DN Hi-C boundary support, node content over the AlphaGenome archive, mouse
synteny and the HOXD interval. Node content is the containment statistic in its modelled form, so the
default is the maximum of seven callers on a statistic correlated with the one it is now scored on.
All seven are therefore scored on the 661 measured pairs, each against its own controls (`uniform`,
`uniform_merged`, `circular`; `count_matched` is left to the committed result, being the slow draw
and unnecessary for a within-caller comparison). Reading rule, fixed here: if the rejected callers
also clear zero, the selection cost is small and the claim's excess is close to what any CTCF-based
caller would give; if only the default clears zero, the measured excess carries a selection premium,
and the spread across the seven callers is the size of it. This is a sensitivity on pairs already
read, not a new test, and it is reported as one.

**Power, and the rule for whether an independent arm runs.** An exact one-sided binomial at 0.05, null
share = the pooled `uniform` control share, alternative = the measured share. An arm runs only if its
regulated pairs reach 80% power at the claim's own effect. Candidates: the IGVF MHC Perturb-seq screen
(`data/results/indep_mhc_crispri.json`), the IGVF DC-TAP K562 remainder (`IGVFFI0957PYTA`, the pairs
not in either ENCODE file, scored by no result here), and the IGVF many-loci K562 library, unreleased.

**Seen before this was committed.** For the DC-TAP file, label-blind: 7,475 rows, 6,512 distal
targeting pairs, 948 elements, 263 genes, K562 only, 15 chromosomes, 29 loci at 2 Mb clustering, 5,069
distal pairs not in either ENCODE file. Then, to do the power arithmetic and for no other purpose, the
count of its positives: 92 distal significant pairs, 51 of them decreases, 36 of those not in either
ENCODE file, from 35 elements. Their containment — which pair falls inside which node — has not been
computed and is not computed by this audit. From the MHC result, already committed: 19 eligible distal
decreases. No containment figure for either set existed when this was committed.

**Why the arms are counted and not run.** Stated in advance: 19 pairs and 36 pairs are the whole of
what is readable at 0 requests, and the arithmetic above is what they can see. Reading their
containment now converts the only unread pairs this project holds into read ones for a test that
cannot separate the claim from chance; they are left unread and logged, so a pooled arm stays possible
when IGVF releases the many-loci table.

**Cost.** 0 AlphaGenome requests; no cache reader, nothing fetched.

## The result: the node caller was built before the pairs existed and reads nothing from them, and the one thing that is not independent is the size of the number (area B, 2026-09-28, lane-nodeindep, later)

The registration above was committed as `da5764e`, the script's split decomposition was repaired in
`366db9a` (the benchmark pair carries no `source_file` field, so the registered decomposition would
have raised on its first regulated pair), and the audit then ran once: 138 caller-chromosome cells,
0 AlphaGenome requests, nothing fetched, no per-element cache opened.
`data/results/node_independence_audit.json`.

**The table, under lane-split's legend.**

| Component | Status | What decides it |
| --- | --- | --- |
| The 661 pairs | EVALUATION-ONLY | Both benchmark files pooled. Nothing is fitted to them, so neither is a training set and neither is a test set; but the held-out file is a reused benchmark (lane-split `0af7e1b`: at least ten scored results read it), so this is not a first touch |
| Which split | **471 training / 190 held-out** | 71% of the claim rests on the *training* file, which is K562 only; the held-out 190 are 118 K562, 34 HCT116, 16 GM12878, 15 WTC11, 7 Jurkat |
| CTCF cCRE input | CLEAN | ENCODE SCREEN V3, external, names no gene and no perturbation |
| Boundary rule (CTCF-only midpoints merged within 5 kb) | CLEAN | `MERGE_BOUNDARIES_WITHIN = 5_000` written 2026-09-10, never changed; the class filter `c.cls == "CTCF-only"` written 2026-09-10 and last touched 2026-09-14 |
| The 50 kb floor | CLEAN | `MIN_DOMAIN = 50_000` written 2026-09-10 in `0127b3f`, never changed since. Inherited as a TAD size floor. It is not neutral for the number — the published control does not pass through it, which is why `uniform_merged` is reported beside `uniform` — but it is not tuned on the pairs |
| Thresholds inside the caller | CLEAN | The default caller has none beyond the class filter and those two lengths. `STRICT_RELATIVE = 0.95` belongs to the rejected orientation callers |
| The annotation | CLEAN | GENCODE 50; the gene is the benchmark's own `measuredGeneSymbol` |
| The random-boundary control | CLEAN | Seed 7 and 20 draws inherited unchanged from `scripts/oriented_domains.py` (2026-09-14); the audit's own `SEED = 7` (2026-09-27) exists to reproduce it, not to choose it |
| The statistic | EVALUATION-ONLY | Nothing fitted; the caller's output is compared with a measurement it never saw |
| Dependence on model output | CLEAN | Stage 2 recomputed on all 23 covered chromosomes with the AlphaGenome archive never opened: counts **identical**, 0 differences |
| **Choice of default caller out of seven** | **EXPOSED** | Not to the pairs — to the same statistic in its modelled form |

**The three facts that make "CLEAN" a measurement and not an assertion.** The default caller
imported in a fresh interpreter pulls in **17 `genomeos` modules** and not one of them names a
CRISPRi source; the only file under `genomeos/genome/` that mentions CRISPRi at all is
`hic_contact.py`, which the closure does not contain. Every constant of the caller was written on
**2026-09-10** (the 50 kb floor, the 5 kb merge, the CTCF-only class, the 0.4 node confidence) or
**2026-09-14** (the control's seed and draw count); the first commit naming the benchmark is
`138824f`, **2026-09-16**, and the first CRISPRi containment figure is `2b77663`, **2026-09-27**. So
every choice inside the caller predates the arrival of the pairs by two to six days and the arrival
of any containment number on them by thirteen to seventeen. And stage 2 reproduces to the count
from the cCREs, GENCODE, the chromosome length and the pairs alone, so the measured arm borrows no
credit from, and no contamination through, the model archive that stage 1 scores.

**What is exposed, and how much it costs.** The default was kept over six orientation callers on
2026-09-14 (`26da99a`) and again on 2026-09-21 (`56e2c50`), on four measurements, one of which —
node content over the deletion archive — is this same containment statistic computed on model
output. No CRISPRi figure existed on either date, and `scripts/oriented_domains.py` reads no CRISPRi
file, so the exposure is not to the pairs. It is to the statistic. The registration fixed the
reading rule before the run: if the rejected callers also clear zero the selection cost is small; if
only the default clears zero the excess carries a selection premium and the spread is its size.

The result is the second case, and not narrowly. Excess over each caller's own `uniform` control on
the same 661 pairs: **ctcf_only +5.89**, oriented_ctcf_only +1.86 (95% -0.54 to +3.93),
oriented_strong -3.97, oriented_strict -3.97, oriented_best_hit -13.50, oriented -14.98. **Zero of
the six rejected callers clear zero; the spread across the seven is 20.87 points.** The default is
the maximum of the seven on the measured statistic, exactly as it was the maximum of the seven on
the modelled one. The direction survives this — no alternative caller reverses it, they simply fail
— but the *size* +5.89 is a selected maximum and is not an unbiased effect estimate.

Two things keep this from being worse than it is. The default was not chosen out of the seven: it
was written on 2026-09-10 and the six challengers were built four days later to displace it and did
not. And the ranking was fixed before any CRISPRi number existed. A selection premium is still a
premium.

**The interval, re-drawn.** This audit's own chromosome bootstrap of the default gives +5.89, 95%
**+3.32 to +8.64**, against the committed **+3.18 to +8.46**: same point, same sign, a different
resample stream (one generator seeded 11 walks all seven callers here; the committed result draws
its own). The claim's interval reproduces to about a sixth of a point.

**Negatives first, on the independent arms.** No arm was run, and the registration decided that in
advance with an 80%-power gate. At the claim's own shares — null 0.699, alternative 0.7579 — an
exact one-sided binomial needs **358 pairs** for 80% power. The IGVF MHC Perturb-seq screen brings
**19** (power 0.127, and all 19 sit in one 4 Mb locus on chr6, so its independent-unit count is 1,
not 19). The IGVF DC-TAP K562 remainder brings **36** (power 0.101). Together, 55 against 358. The
IGVF many-loci K562 library is still unreleased (HTTP 403). So the answer to "is there an
independent arm" is **no, not at 0 requests**, and their containment was left uncomputed on purpose:
reading it now would convert the only unread pairs this project holds into read ones for a test that
cannot separate the claim from chance. They stay logged for a pooled arm when IGVF releases the
many-loci table.

**One correction to the registration's own inventory.** The registration recorded 6,512 distal
targeting DC-TAP pairs; recounted by the script it is **6,446** (7,475 rows, 7,281 targeting, 194
positive control, 6,446 of the targeting rows carrying no promoter gene). Every figure derived from
it is unchanged: 5,069 distal pairs not in either ENCODE file, 51 distal significant decreases, 36
of those fresh, from 35 elements.

**What the claim may be called.** That the node containment direction is **confirmed on a
measurement the node caller never read, and could not have read**: no fitting, no threshold from the
pairs, every constant older than the benchmark, and a stage-2 reproduction that never opens the
model archive. What it may **not** be called: an *independent* confirmation, an *out-of-sample* test
or a *held-out* result — 71% of the pairs are the training file, the held-out file is a reused
benchmark, and the caller that scores them was selected as the best of seven on the same statistic
in modelled form and is the only one of the seven that clears zero. And **+5.89 may not be quoted as
an unbiased effect size**; it is the maximum of a selected family, with the four-baseline range
(+4.5 to +10.2 across `uniform`, `uniform_merged`, `circular` and `count_matched`) and the 89.1%
K562 composition still beside it. The existing claim is not changed by any of this; this section is
the qualification, written beside it.

## The result: matched on the number of scored elements the real unknown still names a coding gene 24 points below its windows, and the per-element rate, which needs no matching, is 29 points down (2026-09-28, lane-elemcount, later)

The registration above was committed as da5764e before the run. The run then happened once, over all 24
chromosomes, with 0 AlphaGenome requests and the per-element cache never opened.

**Gate.** Passed. With the element test switched off and d717b28's 4,000-try cap the same draw function
gave back d717b28's pooled counts to the digit for both target sets: real unknown 44,100 windows drawn,
31,676 carrying, 21,217 naming a coding gene, 531 and 158 blocks; neutral 131,600 / 88,317 / 57,212,
1,181 and 340.

**Bands.** The real unknown's scored-element deciles collapse to the cuts 0, 1, 3, 6, 11 (five bins,
because 351 of the 882 blocks hold no element at all); the neutral tier's to 0, 1, 2, 7.

**Primary.** The element-count-matched difference is **-23.99 points**, 95% over blocks **-27.95 to
-19.93**, over chromosomes -27.71 to -20.56, **n = 531 blocks**. No block was undrawable and no carrying
block failed to draw a carrying window, so the comparison is the same 531 blocks as before. 29.76% of
the blocks name a coding gene against 53.75% of their matched carrying windows. The interval is wholly
below 0, the point is far below the five-point chance band, and the falsifier fired: element count
closes 35.6% of d717b28's unmatched gap and 12.0% of the density-matched -27.25 of 0c8b82d. By the
registered reading, **element count is not the cause of the gap either and clause 2 stays not met**.

**How much of the imbalance the deciles actually closed, and the sensitivity that closes all of it.**
The registered banding is coarse where it matters: with five bins and a top bin open above 11 elements,
a carrying window still holds 11.94 scored elements against 6.18 per carrying block, against 17.90 in
this run's own unmatched draw and 13.27 in 0c8b82d's density-matched one. So the primary closed a little
over half the imbalance, not all of it. The registered sensitivity that closes it exactly -- matching on
the block's own element count rather than its decile, so every accepted window holds the same number of
elements as its block, 6.177 against 6.177 -- gives **-18.57 points, 95% -22.48 to -14.55, n = 531**,
with no undrawable block. With the number of chances made identical, the real unknown is still about
eighteen and a half points below. Matching jointly on the element decile and e1dbcf3's coding-TSS decile
gives -18.07 (-21.73 to -14.36) with both covariates balanced (4.11 TSSs against 4.02); rejecting only
on overlap with the target set gives -18.12 (-22.14 to -13.99). All three sensitivities agree with the
primary in sign and magnitude, and the exact-count point sits a hundredth of a point on the shallow side
of the registered falsifier line, which is fired on the primary as registered.

**The per-element rate, which removes the imbalance arithmetically.** On d717b28's own unmatched windows,
with nothing matched away at all, the share of a block's scored elements that name a coding gene minus
the share of its windows' elements that do is **-28.77 points, 95% over blocks -31.66 to -25.77**, over
chromosomes -31.82 to -25.91, n = 531. Pooled, 301 of the 3,280 elements inside the real unknown name a
coding gene (9.18%) against 201,072 of 456,573 window elements (44.04%). The element-count-matched
windows give -18.81 (-21.70 to -15.82). This is the registered answer to the imbalance that needs no
control at all, and it says the same thing: a scored element inside a constrained-unknown block names a
coding gene about a fifth as often as a scored element in comparable sequence outside one. The gap is not
a number-of-trials effect.

**Secondaries.** The neutral tier under its own element-count matching sits at -23.28 (-26.03 to -20.47,
n = 1,181), and real unknown minus neutral is -0.71 (-5.58 to +4.17): the two tiers remain
indistinguishable, as in 0c8b82d. Per element and unmatched, the neutral tier is -26.43 (-28.52 to
-24.23). On the broader question "moves a gene" the real unknown is -15.89 (-20.18 to -11.77).

**What this settles and what it does not.** The imbalance that 0c8b82d named as the last live
explanation is now tested and rejected. Three quantities have been matched -- block length, local
coding-gene density and the number of scored elements -- and a fourth estimator avoids the third by
arithmetic, and the real unknown stays between 18 and 29 points below its control on every one of them.
No covariate named so far explains the reading. What the run does not settle is why: it cannot separate
"these blocks hold no element that regulates a coding gene" from "the deletion model cannot name a
target for sequence like this". Two further descriptive controls remain testable from the archive alone,
and neither has been run: the class composition of the elements (a block's cCREs may be dELS and
CTCF-only where a window's include PLS and pELS), and the distance from each element to the nearest
coding TSS inside the model's 1 Mb input, which is the sharpest remaining candidate because the model
can only name a gene it can see. Beyond those, only measurement separates the two readings, and the
measured arm stays registered and unrun in its own lane.

**Cost.** 0 AlphaGenome requests. 4,182 s over 24 chromosomes; one chromosome's archive in memory at a
time. Result: `data/results/clause2_element_count_control.json`.


> **Annotation (2026-09-28, lane-design, after a statistical review; the text above is unchanged).**
> Three narrowings. (a) **9.18% against 44.04%** is the share of scored elements for which the deletion
> model names a protein-coding gene, in each arm. It is a **difference in target-naming frequency, not
> in accuracy**: no measured outcome enters it, so "a scored element inside a constrained-unknown block
> names a coding gene about a fifth as often" is a statement about the model's output on two sets of
> sequence and not about either the model's correctness or the sequence's function. (b) "No covariate
> named so far explains the reading" is the right form, and the stronger form some later notes took
> from it is not: adjustment shows that **the tested covariates do not explain the difference**, never
> what does. (c) "the two tiers remain indistinguishable" should read **no difference detected**: the
> contrast is -0.71 points with a 95% interval of -5.58 to +4.17, which permits a real difference of
> over five points in either direction. No equivalence margin can be justified from outside these data,
> so no equivalence test is run; the reasoning is in the lane-design sections below.

## Pre-registration: the two descriptive controls clause 2 has left -- element class and reach inside the model's own window (2026-09-28, lane-tssreach)

Milestone 1.3's clause 2 is settled as not met across three matched controls and one control-free
estimator. The sharpest number needs no matching at all: per element, **9.18% of the 3,280 elements
inside the real unknown name a coding gene against 44.04% of the 456,573 elements in their windows**, a
pooled **-34.86 points**, **-28.77** as a mean over the 531 compared blocks (`faeb0da`). lane-elemcount
named two descriptive controls it did not run and called the second the sharper. This registers both and
runs them once. **Nothing here restates, weakens or re-opens any existing claim**, and neither control can
make the clause pass: both are descriptions of the same elements, and the best either can do is move the
-34.86 toward 0, which replaces "below chance" with "explained by where the blocks sit" -- not with
"attributed". The clause's third part, scored against measurement, is untouched by every outcome below.

**(a) Class composition.** A block's cCREs may be dELS or CTCF-only where a window's include PLS and
pELS. Per element, its ENCODE SCREEN class from `data/results/ccres_<chrom>.bed.gz` field 5, joined by
element id; an id absent from the registry is counted `unclassified` and reported, never dropped.

**(b) Reach.** *The model can only name a gene it can see.* Per element, on its midpoint
`m = (start + end) // 2`:

1. `reach_count` -- GENCODE protein-coding TSSs in `[m - 524288, m + 524288)`, half-open. The window is
   the scorer's own and is read from the scorer, not assumed: `AlphaGenomeAdapter._live_scorer` resizes
   the variant's reference interval to `dna_client.SEQUENCE_LENGTH_1MB` = **1,048,576 bp** before every
   deletion score, and every element in the archive was read in that window.
2. `reach_distance` -- `|m - nearest coding TSS|` in bp, uncapped.
3. `genes_in_window` -- protein-coding genes whose annotated span intersects the same window, as the
   sensitivity closer to what the RNA-seq gene scorer enumerated: a gene whose body reaches the window
   but whose TSS does not is still scored.

The earlier TSS matching (`0c8b82d`) was at block-midpoint level, which is not this quantity.

**Comparison and draw.** Elements inside real-unknown blocks against the elements of that block's own
accepted windows, over the same 531 blocks, on `d717b28`'s **unmatched** draw unchanged -- 50 windows of
the block's exact length per block, uniform inside the span of the chromosome's scored-element midpoints,
rejected on overlap with any organiser block, 4,000 tries, seed 20260913 per chromosome. No matching is
added, because the quantity to be described is the same in every window the -34.86 was computed on. The
coding-TSS-density-matched draw of `0c8b82d` and the neutral tier are registered sensitivities.

**Estimators.** Primary for (b): the per-block difference in mean `reach_count` (block elements minus its
carrying windows' elements), in coding TSSs, with the pooled means, the zero-reach shares and the same
statistic for `log10(1 + reach_distance)` and `genes_in_window` beside it. The statistic that **decides
between the readings** is the per-element rate of `names_a_coding_gene` **directly standardised** on
reach: both arms' elements fall in the fixed strata **0, 1, 2, 3-4, 5-6, 7-9, 10-14, 15-24, 25+** -- cut
points fixed here, not read from the data -- a stratum enters when it holds at least 30 elements in each
arm, block elements in every excluded stratum are reported with their share and their own rate, and the
standardised difference is the block rate minus the block-weighted sum of the window rates. It answers:
*if the block's elements had sat in windows of the same reach, how much of the gap would be left?* The
same standardisation on the six classes answers (a), with the per-class rates in both arms so it can be
seen whether the gap lives between classes or inside them; class x stratum jointly is a sensitivity.

**Interval.** 95% percentile bootstrap over the compared blocks, 10,000 resamples,
`default_rng(20260913)`, with the bootstrap over the 24 chromosomes beside it -- the two bootstraps of the
matched control, unchanged. For a standardised difference the resample is over blocks, recomputing the
block weights and rates and holding the window per-stratum rates at their pooled values, because the
window arm is the reference population and has 100 times the elements.

**Gate.** Before any reach or class figure is read, this script's draw must reproduce `d717b28`'s pooled
counts for both target sets to the digit (real unknown 44,100 / 31,676 / 21,217, 531 and 158 blocks;
neutral 131,600 / 88,317 / 57,212, 1,181 and 340) **and** `faeb0da`'s per-element counts (3,280 block
elements, 301 naming; 456,573 window elements, 201,072 naming). If either fails, nothing is reported.

**Thresholds, fixed before the run.** The half-the-gap convention of `e1dbcf3` and `da5764e`, on the
pooled per-element scale: a description **explains** the gap when standardising on it leaves the pooled
difference above **-17.43 points**, and does not when it stays at or below. "Far fewer" for reach is
fixed as either the pooled mean `reach_count` inside the blocks being below half the window mean, or the
zero-reach share inside the blocks exceeding the window share by more than 10 points. "Comparable" is
fixed as the per-block reach difference's 95% interval lying wholly inside +/-1.0 coding TSS of 0 and the
two zero-reach shares differing by less than 2 points.

**The readings, in advance.**

- **Reach far lower and standardising closes more than half.** The instrument was asked to name a gene it
  could not see, and the -28.77 is partly an artefact of reach. What clause 2's failure *means* changes:
  it stops being a statement about the blocks' sequence and becomes a statement about where the blocks
  sit relative to coding genes -- which is what the organiser selected them for. The clause still does not
  pass; the standardised difference becomes the number to quote, and it is still negative.
- **Reach far lower but the gap survives standardisation.** The reach imbalance is real and must be
  stated, but it is not the explanation; the reach-standardised difference becomes the number to quote
  and the failure stays about the sequence.
- **Reach comparable.** The model could see as many coding genes from inside the blocks as from their
  windows, the instrument was not asked to do the impossible, and the failure is about the sequence, not
  the window. Nothing about the clause's status changes and the -28.77 stands as written.
- **Reach comparable yet the gap closes.** Incoherent on its face: to be reported as a likely
  implementation error in this script, not as a finding about clause 2, until it is traced.
- **(a) explains** when the class-standardised difference closes more than half: the blocks hold a
  different mix of cCRE classes and the per-element gap is largely composition, and the per-class rates
  become the number to quote. **(a) does not explain** when it keeps at least half: inside every class
  with enough elements the block elements still name a coding gene far less often, so class composition
  is a real difference between the two element sets and not the cause of the gap.

**Cost.** 0 AlphaGenome requests; the per-element response cache is never opened; one chromosome's
archive in memory at a time. No matching is added, so no window's contents are evaluated on a rejected
try: lane-elemcount's 3,811 s chr2 cost came from exactly that and does not arise here. Script:
`scripts/clause2_reach_control.py`; result: `data/results/clause2_reach_control.json`.

## The measured arm of clause 2: 3 of the 882 real-unknown blocks hold an element a screen ever tested against a coding gene, and the arm cannot decide (area I, 2026-09-28, lane-measured)

Registration: `clause2_measured_arm.PRE_REGISTRATION` (`c17eedc`, committed before any assay was
counted), on top of `clause2_matched_control.PRE_REGISTRATION` (`e1dbcf3`, above), whose target sets,
window draw, decile matching, seed, estimator, both bootstraps and reproduction gate are imported
rather than restated. The model question rides through the same draw, so the gate is free: this run
gave back `0c8b82d`'s committed matched primary to the digit -- -27.25 points, 95% -30.91 to -23.58,
531 compared blocks, 44,081 windows, 158 blocks yes, 12,876 window-yes of 22,228 -- and only then was
a measured figure read. Script: `scripts/clause2_measured_arm.py`, tests in
`tests/test_clause2_measured_arm.py`; result `data/results/clause2_measured_arm.json`. 0 AlphaGenome
requests; the per-element response cache is never opened.

**What was left registered and unrun, and what this is.** Every lane on clause 2 ended with the same
sentence: beyond the descriptive controls, only measurement separates "these blocks hold no element
that regulates a coding gene" from "the deletion model cannot name a target for sequence like this".
The arm they registered (`second_route_not_run`) is a **designed** experiment -- 284,001 unsynthesised
oligos over 878 blocks -- and it correctly said nothing measured exists for it. That is still true and
it is still unrun. This is the other arm: the measurements the project already holds (CRISPRi,
lentiMPRA, VISTA, saturation mutagenesis, through `genomeos/attribution/measured.py`), asked of the
same elements, in the same blocks, against the same windows, with the same statistic.

**Coverage, which is the first result and was registered as such.** Of the 882 real-unknown blocks
(3,280 scored elements in the 531 that carry one):

| | blocks | elements |
| --- | ---: | ---: |
| eligible: some assay's interval shares a base with an element | 89 | 151 |
| measured: some assay measures an element at reciprocal overlap 0.5 | 59 | 78 |
| — lentiMPRA | 57 | 75 |
| — CRISPRi | 3 | 3 |
| — VISTA | 1 | 1 |
| — saturation mutagenesis | 0 | 0 |
| CRISPRi tested it against a protein-coding gene | **3** | **3** |
| measured to move a coding gene (decrease or increase) | **0** | **0** |
| measured a well-powered null on every coding gene tested | 3 | 3 |
| measured active with no gene named (lentiMPRA or VISTA) | 8 | 9 |

**823 of the 882 blocks hold no measured element at all**, and 3 blocks -- 0.34% -- hold the only kind
of element that can answer a clause asking for a gene. Of the 3, one has a matched window that is also
CRISPRi-tested, so the primary compares one block. The registered floor is 20
(`measured.MIN_FOR_A_COMPARISON`, adopted unchanged), so **no interval is read and the registered
reading is `cannot_decide`**: the measured arm cannot separate the two accounts of clause 2, and that
was registered in advance as a finished lane and a real result rather than a failure.

**Why 0 of 3 is not an answer.** All three tested elements came back well-powered nulls, which points
the same way the model does. It is not evidence: 13.45% of the matched windows' CRISPRi-tested
elements move a coding gene (30 of 223), and under that rate three tested elements come back all-null
65% of the time. The neutral tier is thinner still -- 1 block tested, 1 regulated, 0 comparable.

**The coverage that would be needed**, by the formula fixed before the count (d = 0.2725, the model
arm's own matched difference; s = 0.4324, the standard deviation of this run's per-block model-arm
differences): **20 compared blocks for 80% power**, which is also the floor, against the 1 there is --
19 short. In assay terms that is a CRISPRi screen reaching roughly twenty times as much constrained
unknown sequence as every published screen this project holds reaches today: about 2.3% of the real
unknown's blocks, against 0.34%.

**A well-powered negative that is not about clause 2 at all.** On the same windows, with the same
estimator and 531 compared blocks, a real-unknown block holds a measured element 11.11% of the time
against 25.99% for its length- and gene-density-matched windows: **-12.23 points, 95% -14.81 to
-9.57**. Measurement is not missing from these blocks at random with respect to the control. The
assays were pointed at candidate regulatory sequence near expressed genes in a handful of cell lines,
which is the same property the deletion model scores highly, so the sequence clause 2 is about is
exactly the sequence the screens avoided. Any future measured arm has to beat that selection, not
just add volume. (The same table at reciprocal overlap 0.25 gives 72 measured blocks and at 0.75 gives
10, with 3 and 0 CRISPRi-coding blocks: the rule is visible, and no setting of it produces a usable
arm.)

**What this does and does not settle.** It does not decide clause 2, and it was registered not to
pretend otherwise. Both readings stay open and both are written down: if measurement ever says these
blocks hold elements that regulate coding genes, the failure is the model's and the clause's stated
reason changes from "these blocks hold nothing" to "the deletion model cannot name a target for
sequence like this"; if measurement says they hold none while their matched windows do, the model was
right and clause 2's **wording** is what is wrong, because a milestone cannot ask for a gene and a
tissue from sequence that regulates no gene. What the run does settle is that the project cannot reach
either reading from the measurements it holds, that the shortfall is 19 blocks rather than a few
thousand, and that the measured layer is biased against the target set by 12 points before any of that
begins. The descriptive route is finished; clause 2 now waits on an experiment, and the size of the
experiment is now a number.

**Cost.** 0 AlphaGenome requests. 527 s over 24 chromosomes; one chromosome's element archive and one
chromosome's measured layer in memory at a time.


> **Annotation (2026-09-28, lane-design, after a statistical review; the text above is unchanged).**
> Four narrowings, three of them about the same paragraph. (a) **"The coverage that would be needed ...
> 20 compared blocks for 80% power ... 19 short" is withdrawn as a power result.** Both inputs of that
> formula are model output -- the effect is the model arm's own matched difference in target-naming
> frequency and the dispersion is the standard deviation of the model arm's per-block differences -- and
> they were used to size an experiment whose endpoint is a **measurement**. The "19" quoted beside it is
> `MIN_BLOCKS - compared_now` = 20 - 1, the distance to `measured.MIN_FOR_A_COMPARISON`, a
> **minimum-reporting rule**; that the formula also returned 20 is a coincidence, not a confirmation.
> (b) Consequently **"the shortfall is 19 blocks" and "the size of the experiment is now a number" are
> withdrawn.** A power simulation of the actual design gives a **range**: 20 to 75 compared blocks if
> the blocks truly regulate a quarter as often as their windows, 50 to 300 at a half, and 300 to 5,000
> at three quarters, with no configuration reaching 80% power at the null. (c) **"13.45% of the matched
> windows' CRISPRi-tested elements move a coding gene (30 of 223)"**: the 223 and the 30 are **windows**
> carrying such an element, not elements, so they are not 223 independent measurements. The
> element-level rates over the benchmark's own distinct tested elements are 12.15% (479 of 3,941,
> training) and 14.56% (247 of 1,697, held-out), which agree with it. (d) **"823 of the 882 blocks hold
> no measured element at all"** means no qualifying measurement attached to the block's *scored
> elements* under the 0.5 reciprocal-overlap rule; it does not mean no part of those blocks was ever
> measured, and at 0.25 or 0.75 the measured count is 72 or 10 blocks against 59. The committed result
> file carries the same corrections as an additive key, and `coverage_needed` is kept in the code as the
> record of what was computed, with `reporting_floor_and_power_assumptions` beside it.

## The result: the model really can see three times fewer coding genes from inside a constrained-unknown block, and that accounts for under a third of the gap (2026-09-28, lane-tssreach, later)

Both controls of the registration above ran once, genome-wide, 280 s, **0 AlphaGenome requests**, the
per-element response cache never opened. The gate passed: the draw gave back `d717b28`'s pooled counts
(44,100 windows drawn, 31,676 carrying, 21,217 naming, 531 and 158 blocks; neutral 131,600 / 88,317 /
57,212, 1,181 and 340) and `faeb0da`'s per-element counts (3,280 block elements, 301 naming; 456,573
window elements, 201,072 naming) to the digit.

**The negative first: neither control explains the gap.** The registered threshold was that standardising
on a description must leave the pooled per-element difference above **-17.43 points**. Standardising on
reach leaves **-24.46 points** (95% over blocks -26.36 to -22.47; over chromosomes -26.65 to -22.33;
n = 531 blocks, 3,252 of the 3,280 block elements). Standardising on class leaves **-33.57**
(-35.40 to -31.49). Reach closes **29.8%** of the -34.86, class **3.7%**, the two together **31.2%**
(-24.00, -25.86 to -22.06). The registered outcome is `reach_short_but_gap_survives`.

**Reach is nonetheless far lower inside the blocks, and this is the first time it has been measured where
the model reads it.** Per element, over the scorer's own 1,048,576 bp window centred on that element:

| | inside the blocks | in their windows |
| --- | --- | --- |
| coding TSSs in reach (mean) | **2.229** | **7.219** |
| elements with **no** coding TSS in reach | **28.08%** | **6.20%** |
| coding genes whose span reaches the window (mean) | 2.55 | 7.78 |
| mean log10(1 + bp to the nearest coding TSS) | 5.381 | 4.855 |

Per block, the difference in mean reach is **-4.39 coding TSSs** (95% over blocks -4.95 to -3.81, over
chromosomes -5.37 to -3.65), in genes whose span reaches the window -4.54 (-5.11 to -3.96), and in
log10 distance +0.349 (+0.281 to +0.415) -- a block's elements sit about 2.2x further from the nearest
coding gene start than the elements of sequence of the same length elsewhere. The earlier midpoint-level
matching could not have shown this: it matched the count around a block's midpoint, not what each element
could see.

**The gap lives inside every stratum of reach, not between them.** Block rate against window rate, in
points, by coding TSSs in reach: **0: 0.98% against 19.37% (-18.39)**; 1: 4.82 / 30.70 (-25.87); 2: 18.06
/ 39.19 (-21.13); 3-4: 9.70 / 44.32 (-34.62); 5-6: 21.58 / 46.98 (-25.41); 7-9: 17.57 / 48.23 (-30.66);
10-14: 38.00 / 52.48 (-14.48). The two richest strata hold 19 and 9 block elements and fall below the
registered minimum cell of 30; they are excluded and reported here rather than dropped. Note that the
reach = 0 stratum is **not** a hard ceiling: 19.37% of window elements with no coding TSS in reach still
name a coding gene, because the scorer reads genes whose body enters the window without their start.

**(a) Class composition is a real difference and a small one.** The all-element sweep scored only two
registry classes in this comparison, so the mix is nearly fixed by the sweep's own selection: **dELS
96.59% / pELS 3.41%** inside the blocks against **87.20% / 12.80%** in their windows, a per-block
promoter-like share difference of **-0.0298** (-0.0545 to -0.0042). Inside dELS the block elements name a
coding gene 8.46% against 42.27% (-33.81 points) and inside pELS 29.46% against 56.07% (-26.61): the gap
is inside both classes, not between them.

**The sharpest sensitivity.** On `0c8b82d`'s coding-TSS-density-matched windows, element-level reach comes
almost into balance -- per-block difference **-0.14 coding TSSs, 95% -0.31 to +0.05**, means 2.229 against
3.279 -- and the per-element gap is still **-27.72 crude and -23.89 reach-standardised** (-25.81 to
-21.88). Where the model can see about as many coding genes from both arms, it still names one three
times less often from inside the blocks. The neutral tier behaves the same way (reach 2.42 against 6.64,
crude -34.01, reach-standardised -20.60, -21.93 to -19.17), and on the looser question `moves_a_gene`
standardising changes almost nothing (reach -34.31, class -36.87).

**What this settles and what it does not.** Reach is a real imbalance in the instrument and must be
stated wherever the -28.77 is: the model is being asked about elements from which it can see about a
third as many coding genes, and 28% of them have no coding start in the window at all. But it is not the
explanation -- standardising on it keeps seven tenths of the gap, and the gap is present inside every
stratum including the richest ones. Clause 2's failure is about the sequence, not the window, and nothing
above changes the clause's status: it stays not met, on the ground `d717b28`, `0c8b82d` and `faeb0da` put
it on. One thing is named and not run: `genes_in_window` was registered as a descriptive statistic only,
and standardising on it rather than on TSS reach was not registered, so it is left for whoever takes it
up rather than computed after the fact. The clause's third part, scored against measurement, is
untouched. Result: `data/results/clause2_reach_control.json`.


> **Annotation (2026-09-28, lane-design, after a statistical review; the text above is unchanged).**
> **"Clause 2's failure is about the sequence, not the window" is withdrawn**, and so is the same
> sentence in the registered reading `reach_short_but_gap_survives` and in the class reading. What the
> run establishes is narrower and still worth having: **the gap persists after adjustment for the
> tested reach and class variables**. Standardising on a covariate and finding the difference survives
> says that covariate does not explain it; it does not identify what does. Unmeasured context,
> selection effects, limitations of the model itself and remaining geometric differences between the
> two element sets are all still open. The measured quantities in this section -- 2.229 coding starts
> against 7.219, 28.08% against 6.20% with none in the window, the per-stratum rates -- are unaffected
> and stand as written. The committed result file carries the same correction as an additive key.

## Pre-registration: a power simulation for the experiment clause 2 actually needs, and an equivalence decision taken before the contrasts are read (2026-09-28, lane-design)

A statistical review of the four clause 2 lanes found five places where the notes claim more than the
data carry, and the coordinator's corrections are in docs/ROADMAP.md beneath milestone 1.3. This lane
takes the two that are matters of design rather than of wording -- the sample size and the equivalence
claim -- and settles them with a simulation registered in advance. Script:
`scripts/clause2_design_power.py`, whose `PRE_REGISTRATION` holds every grid and justification below;
tests in `tests/test_clause2_design_power.py`. No simulated block had been drawn at the registered
settings when this was committed. **Nothing here re-opens, weakens or restates a committed figure**,
no new descriptive statistic about the blocks is computed, and the run makes **0 AlphaGenome
requests**.

**What is wrong with the number this replaces.** `clause2_measured_arm.coverage_needed` published
`n_for_80_percent_power = 20` blocks and "19 short". Two separate faults, and they compound.

1. Both inputs of that formula are **model output** used to size an experiment with a **measured**
   endpoint. The effect is the model arm's own matched difference in how often a deletion model names
   a target, and the dispersion is the standard deviation of the model arm's per-block differences.
   Neither is an estimate of anything about measured regulation. They are named here so that they can
   be excluded by name, and the script carries a test asserting that neither appears as a number in
   it.
2. The "19" that was reported is `MIN_BLOCKS - compared_now` = 20 - 1: the distance to
   `measured.MIN_FOR_A_COMPARISON`, a **minimum-reporting floor**, which is a rule about when a number
   may be printed. The power formula returned 20 as well, and the coincidence made a reporting rule
   read as a power result.

**The design being powered.** A CRISPRi screen (or an assay of equivalent power that names a gene)
tests elements inside constrained-unknown blocks against protein-coding genes, and elements inside the
same matched windows the four lanes drew. Per element the endpoint is the committed one: a significant
change in a coding gene's expression on silencing, either sign. Per block it is 1 if at least one
tested element is found to regulate a coding gene. The estimator and the reading rule are the
committed ones, unchanged: the mean over compared blocks of (block endpoint) minus (share of that
block's tested windows whose endpoint is 1), read as a detection when the 95% interval over blocks
excludes 0.

**The grids, each with what fixes it.**

- **Window regulation rate.** The one empirical anchor the project holds: 30 of 223, 13.45%, with the
  exact binomial interval 9.26% to 18.64% as the grid. Two independent anchors are reported beside it
  and not pooled in: 12.15% of the 3,941 distinct tested elements of the benchmark's training file and
  14.56% of the 1,697 of the held-out file have at least one significant pair. **A correction is
  recorded here rather than left to be found**: the 223 and the 30 are *windows carrying a tested
  element*, not tested elements, because they come from a count of windows; several windows of one
  block can carry the same element, so those 223 are not 223 independent measurements and the exact
  interval understates the uncertainty. The element-level anchors are the element-level ones.
- **Block regulation rate.** Registered as a *ratio* to the window rate: 1.0 (the null the design must
  be able to sit at), 0.75, 0.5, 0.25, 0.1, 0.0 (the clause's strongest form). **The model arm's
  -27.25 points is not used**, here or anywhere: it is a difference in target-naming frequency, not an
  estimate of a difference in measured regulation, and no published estimate of the latter exists,
  which is why the answer is a range over this grid rather than a number.
- **Assay sensitivity.** The means of the five power columns of the ENCODE benchmark files
  (`PowerAtEffectSize10..50`), read from the tables rather than assumed, and reported by split because
  the held-out file is better powered. The columns are not monotone in effect size -- the 25% column
  is the benchmark's own inclusion filter and sits near 1 for nearly every pair -- so they are used as
  five separate sensitivity levels and the non-monotonicity is reported rather than smoothed.
- **False positive rate.** 0 (the benchmark controls its own) and 0.01 per tested element.
- **Tested elements per block.** 1 (the minimum that makes a block testable, and what the three tested
  blocks that exist today have), 2, and 6 (6.177, the scored elements per carrying block of the
  matched control).
- **Tested windows per block.** 1 (today's 223 windows over 134 blocks is about 1.7), 3, 10.
- **Clustering.** Measured, not guessed. The one-way ANOVA estimator on the benchmark's 5,638 distinct
  tested elements puts the intraclass correlation of "this element regulates a coding gene" at 0.375
  within 25 kb, 0.350 within 50 kb, 0.300 within 100 kb, 0.273 within 250 kb and 0.209 within 1 Mb;
  the same estimator with chromosomes as clusters gives 0.0142. The grids are 0, 0.21, 0.30, 0.35
  within a block and 0, 0.014, 0.03 between chromosomes, and the script recomputes all of them at run
  time.

**How the two levels of clustering are handled, said plainly.** Elements within a block are a
**hierarchical simulation**: each block, and each window separately, carries a logit-normal random
intercept whose latent-scale variance is fixed from the registered ICC, and its elements are
conditionally independent given it. Blocks within a chromosome are a **design effect** applied to the
variance of the mean, `1 + (N / 24 - 1) * rho`, not simulated. That asymmetry is deliberate and is the
conservative thing to be explicit about: the committed interval is a bootstrap over blocks, which
assumes blocks independent, so wherever the chromosome ICC is above 0 the committed test is
anti-conservative and the design effect is what says by how much.

**Reference configuration and sweep.** The reference is the middle or the measured value of each grid:
window rate 13.45%, sensitivity `PowerAtEffectSize20` pooled, no false positives, 6 elements per
block, 3 tested windows per block, element ICC 0.30, chromosome ICC 0.014. The primary sweep is the
block-rate ratio at that reference; the sensitivity analysis is one factor at a time. The grids are
not fully crossed, because a range whose ends are combinations nobody would defend is not a range.

**How the numbers are obtained.** For each configuration, 200,000 simulated blocks with their windows
give the mean and standard deviation of the per-block difference; power at every sample size then
follows from the normal approximation to the committed interval, with the chromosome design effect.
That the approximation is safe is **checked and not assumed**: registered cells (ratio 0.5 and 0.25,
at 20, 50 and 200 blocks, chromosome ICC set to 0 so the comparison is of the approximation alone) are
also run with the actual percentile bootstrap, 2,000 resamples over 1,000 simulated experiments each.
**Falsifier:** if the largest disagreement in power exceeds 0.05, the normal-approximation table is
reported as unreliable and the bootstrap numbers are the ones read.

**What is reported.** For every configuration the smallest sample size reaching 50%, 80% and 90%
power, and across configurations the **range** at each ratio with the assumption at each end.
Sample sizes are in *compared blocks*: blocks with a tested element in both arms. The
minimum-reporting floor of 20 is applied as a floor on the answer and is **never** reported as a power
result. At ratio 1.0 the two arms have the same true rate, so what the table calls power is the test's
false-positive rate, and it is reported under its own name as the null calibration. Every number is a
probability of detection under its configuration's assumptions, and the result says so in its own
text: 80% power means failing to detect a real effect of the assumed size one time in five, and no
sample size guarantees that an experiment decides clause 2.

**Equivalence, decided before the contrasts are read.** The committed real-minus-neutral contrast is
-0.10 points, 95% -4.55 to +4.30 (`names_a_coding_gene`, density-matched) and -0.71, -5.58 to +4.17
element-count-matched. "The neutral tier behaves identically" is not what those intervals say. An
equivalence claim needs a margin justified from **outside** these data, and the admissibility rule is
fixed here before looking: a margin may come from a published threshold for a meaningful difference in
regulatory annotation rate, from a stated decision cost, or from a requirement fixed elsewhere in this
project before this comparison. The project's own five-point chance band is **inadmissible**, because
it was chosen in d717b28 for this very comparison. If no admissible margin exists, **no equivalence
test is run** and the finding is stated as "no difference detected". The script then reports, labelled
as a description and never as a margin, the smallest margin at which a two-one-sided test on the
committed interval would pass -- a statement about what the data would need, not about what anyone has
justified.

**The denominator table.** Every population the four lanes count over, each count with its denominator
and the rule that produced it, read from the committed result files rather than retyped, including the
overlap-rule sensitivity at 0.25, 0.5 and 0.75 and the two shares that share a numerator (59 of 882
tier blocks, 6.7%; 59 of the 531 carrying blocks, 11.1%).

**Registered honestly about its own history.** The script was executed once at non-registered settings
(2,000 simulated blocks, 20 check experiments) as a code check before this commit, and one reported
quantity was added afterwards: the null calibration described above. No grid, threshold, reference
value, reading or falsifier was changed after any output was seen.

**Cost.** 0 AlphaGenome requests. The model is not called; the per-element response cache is never
opened; the committed result files are read and never rewritten.

## The result: the measured experiment behind clause 2 needs between 20 and 5,000 compared blocks depending on an effect nobody has estimated, and "19 short" was a reporting floor (2026-09-28, lane-design, later)

The registration above was committed as `f6cb4d7` before the registered run. The run then happened
once, in 11 s, with **0 AlphaGenome requests**. Result: `data/results/clause2_design_power.json`.

**The negatives first, and there are four.**

1. **The registered falsifier fired.** The check on the normal approximation to the committed interval
   allowed a disagreement of 0.05 in power. The largest is **0.0549**, at ratio 0.25 with **20 blocks**
   -- the reporting floor -- where the bootstrap gives 0.724 and the approximation 0.779. So the
   approximation is **optimistic at the floor** and, by the rule fixed in advance, the bootstrap
   numbers are the ones to read there. A post-hoc extension of the same check, named as post-hoc in
   the result, locates the boundary: above the floor the largest disagreement is 0.0399 (ratio 0.25 at
   30 blocks) and it falls to 0.026 or less from 50 blocks upward. **Every entry in the tables below
   that reads 20 or 30 blocks is therefore a lower bound**, and anyone designing at that size should
   simulate rather than read it off.
2. **No sample size can be given without an assumed effect, and there is no estimate of that effect.**
   The ratio of the blocks' true regulation rate to the windows' is the one quantity that decides the
   answer, and nothing measures it. That is not a gap in this lane; it is the reason clause 2 is
   undecided.
3. **The measured arm's own null is not powered at all by anything the project can afford.** At ratio
   0.75 -- the blocks regulating three quarters as often as their windows -- **no configuration
   reaches 80% power below 300 compared blocks**, and the reference configuration needs 750. Against
   the 3 blocks that carry a tested element today, and the 1 the primary could compare, a difference
   of that size is out of reach by three orders of magnitude.
4. **No equivalence margin can be justified**, so "the neutral tier behaves identically" cannot be
   rescued into a claim. It becomes "no difference detected".

**The sample-size range, with its assumptions.** Compared blocks -- blocks carrying a tested element in
**both** arms -- for 80% correct-direction power, across the 17 registered configurations. The
reference configuration is in brackets.

| true block rate, as a multiple of the windows' | 80% power at | reference | what widens it | what narrows it |
| --- | --- | --- | --- | --- |
| 1.0 (the null) | never, by construction | — | — | — |
| 0.75 | **300 to 5,000** | 750 | 1 tested element per block (5,000) | no clustering inside a block (300) |
| 0.5 | **50 to 300** | 100 | 1 tested element per block (300) | no clustering inside a block (50) |
| 0.25 | **20 to 75** | 30 | 1 tested element per block (75) | no clustering inside a block (20) |
| 0.1 | **20 to 30** | 20 | 1 tested element per block (30) | — |
| 0.0 (the clause's strongest form) | **20** | 20 | — | — |

The assumptions behind every row: the windows' measured rate is 13.45% on the observed scale, assay
sensitivity is the pooled mean of `PowerAtEffectSize20` (0.6674), there are no false positives, six
elements are tested per block and three tested windows accompany it, elements within a block correlate
at 0.30 and blocks within a chromosome at 0.014. Varying each of those one at a time is what produces
the ranges. **The single largest lever is not the effect but the number of elements tested per block**:
testing one element rather than six multiplies the requirement by three to seven, because a block's
endpoint is "any tested element regulates" and one element is a much blunter instrument than six.
Clustering is the second lever, in the other direction: pretending elements inside a block are
independent makes the experiment look between a third and a half cheaper than it is.

**Power is a probability.** Every number above is the probability that an experiment of that size
produces an interval excluding 0, *if* its configuration's assumptions hold. At 80% power a real
effect of the assumed size is missed one time in five. No sample size makes an experiment guaranteed
to decide clause 2, and none of these numbers says that an experiment of that size would.

**The anchors, all measured rather than chosen.** Assay sensitivity comes from the benchmark's own
power columns: means of 0.4324, 0.6075, 0.6674, 0.9780 and 0.7052 for effects of 10, 15, 20, 25 and
50%, which is **not monotone** -- the 25% column is the benchmark's own inclusion filter and sits near
1 for nearly every pair -- so the five are used as five sensitivity levels and the non-monotonicity is
reported rather than smoothed. The held-out file is better powered than the training file throughout
(0.8497 against 0.5903 at 20%). Clustering comes from the one-way ANOVA estimator on the benchmark's
5,638 distinct tested elements: the intraclass correlation of "this element regulates a coding gene"
is **0.3745 within 25 kb, 0.3495 within 50 kb, 0.2998 within 100 kb, 0.2727 within 250 kb, 0.2090
within 1 Mb**, and **0.0142 between chromosomes**. The window rate anchor is 30 of 223, 13.45%, exact
interval 9.26% to 18.64%, and the benchmark's own element-level rates agree with it without being
pooled into it: **12.15%** of 3,941 tested elements in the training file and **14.56%** of 1,697 in
the held-out file, 12.88% over all 5,638.

**A correction to the anchor's own wording.** The 223 and the 30 are **windows carrying a
CRISPRi-tested element**, not tested elements: they come from a count over windows. Several of a
block's 50 windows can carry the same tested element, so the 223 are not 223 independent measurements
and the exact interval above understates the uncertainty. The number survives its mislabel -- the
element-level rates land in the same place -- but the unit is now written down.

**Null calibration.** At ratio 1.0 the two arms have the same true rate, so what the table would call
power is the test's false-positive rate. It is **0.0500 at the reporting floor and 0.0500 to 0.0508 at
500 blocks** in all 17 configurations, so the table's null is calibrated and the rest of it can be
read. The Monte Carlo error of the simulated null mean is under 0.13 points.

**The reporting floor, separated from every power statement.** `measured.MIN_FOR_A_COMPARISON` is 20
compared blocks and is a rule about when a number may be printed. "19 short" is 20 minus the 1 block
compared today and is that rule's distance, nothing else. The formula published beside it in
`clause2_measured_arm.coverage_needed` took **both** its inputs from the model arm -- the effect to
detect and the dispersion -- to size an experiment with a measured endpoint, and returned 20, so the
two appeared to confirm each other. They are unrelated quantities that happened to agree.
`coverage_needed` is **kept exactly as it was**, because it is the record of what that run computed
and what `a7f207f` reported; `reporting_floor_and_power_assumptions` is added beside it and separates
the floor, the two model-derived assumptions (each labelled as an assumption, with its source) and the
provisional calculation, which states in its own field that it is not a power result for the measured
experiment. A test asserts that neither 0.2725 nor 0.4324 appears as a number anywhere in the
simulation script.

**Equivalence: no margin, so no test.** The admissibility rule was fixed before the contrasts were
read. Nothing supplies a margin from outside these data: no published threshold states what difference
in regulatory annotation rate is meaningful, the project states no decision cost for treating the
neutral tier as the real unknown, and the only in-project candidate -- the five-point chance band --
was chosen in `d717b28` for this very comparison and is inadmissible by that rule. **So no equivalence
test is run and the finding is "no difference detected."** The committed contrast of -0.10 points,
95% -4.55 to +4.30, permits a real difference of up to about four and a half points in either
direction. For the record, and labelled in the result as a description and never as a margin, the
smallest margin at which a two-one-sided test on the committed interval would pass is 4.55 points for
`names_a_coding_gene` density-matched, 5.58 element-count-matched, 7.87 and 10.53 on `moves_a_gene` --
what the data would need, not what anyone has justified.

**The denominators, every population in one table.** Read from the committed result files, not
retyped. Each count carries its denominator and the rule that produced it.

| population | count | denominator | share | rule |
| --- | ---: | --- | ---: | --- |
| constrained-unknown blocks, copies out | 882 | the tier itself | — | `organise.blocks`, 24 chromosomes, duplicates removed |
| ... carrying a scored element | 531 | 882 blocks | 60.2% | a cCRE midpoint inside the block |
| scored elements inside them | 3,280 | elements, not blocks | — | every scored element of a carrying block |
| blocks eligible for an assay | 89 | 882 blocks | 10.1% | one shared base with a tested interval |
| blocks eligible for an assay | 89 | 531 carrying blocks | 16.8% | same numerator, other denominator |
| blocks with a measured element | 59 | 882 blocks | **6.7%** | reciprocal overlap 0.5 both ways |
| blocks with a measured element | 59 | 531 carrying blocks | **11.1%** | same numerator, other denominator |
| blocks CRISPRi-tested against a coding gene | 3 | 882 blocks | 0.34% | overlap 0.5, gene protein-coding |
| blocks measured to move a coding gene | 0 | 882 blocks | 0% | significant either sign |
| blocks with no measured element at all | 823 | 882 blocks | 93.3% | **no assay measures any of its scored elements at overlap 0.5** -- not "no part of the block was ever measured" |
| blocks the measured primary could compare | 1 | 882 blocks | 0.11% | tested element in the block AND in one of its windows |
| blocks with a measured element, overlap 0.25 | 72 | 882 blocks | 8.2% | the rule set to 0.25 |
| blocks with a measured element, overlap 0.75 | 10 | 882 blocks | 1.1% | the rule set to 0.75 |
| blocks CRISPRi-tested, overlap 0.25 / 0.75 | 3 / 0 | 882 blocks | 0.34% / 0% | no setting of the rule produces a usable arm |
| neutral-tier blocks | 2,632 | the tier itself | — | the secondary target set |
| ... carrying a scored element | 1,181 | 2,632 blocks | 44.9% | as above |
| matched windows drawn | 44,081 | windows | — | 50 per block, exact length, coding-TSS decile |
| ... carrying a scored element | 31,554 | 44,081 windows | 71.6% | a cCRE midpoint inside the window |
| ... carrying a CRISPRi-tested element | 223 | 44,081 windows | 0.51% | **windows, not elements** |
| ... holding one measured to move a coding gene | 30 | 223 windows | **13.45%** | the anchor, in windows |
| blocks the model arm compares | 531 | 882 blocks | 60.2% | carries an element and has a carrying window |
| window-elements in the per-element estimator | 456,573 | elements in the unmatched windows of the 531 | — | 9.18% against 44.04% is over these |
| block elements entering the reach standardisation | 3,252 | 3,280 block elements | 99.2% | strata with 30 elements in each arm |

**The sweep of the four lanes' own write-ups.** Their dated sections above now carry an annotation
paragraph each, added beneath the original text and never in place of it, and the four committed
result files each carry an additive
`corrections_after_the_statistical_review_2026_09_28` key naming the field or the string it corrects.
No existing key or value in any of them was changed, and the annotator refuses to write if one would
be. What the sweep found, beyond the five points the review named: the 30-of-223 unit error above,
which no lane would have caught from its own numbers.

**What this settles and what it does not.** It settles that the size of clause 2's experiment was
never known: what was published as its size was a reporting floor with a model-derived calculation
beside it. It replaces that with a range that is honest about depending on an unestimated effect, and
it says which design choice moves the range most (how many elements per block are tested, not how many
blocks). It does not make clause 2 decidable, it does not estimate the effect, and it does not say
that an experiment of any size here would decide the clause -- power is the probability of detection
under assumptions. Clause 2 stays not met, and now waits on an experiment whose size is a range with
its assumptions attached rather than a number.

**Cost.** 0 AlphaGenome requests. 11 s; no genome archive is opened, and the committed clause 2
results are read and never rewritten.

## Reproduced before anything changed: the committed power simulation, set from a 13.45% window anchor, produces 16.80% element detection and 56.96% positive windows (item 12 S2, 2026-09-28, lane-s2)

The second external review (docs/ROADMAP.md section 5, item 12, row S2) reports that the simulation
behind the 20-to-5,000 range above does not reproduce the quantity it was anchored on. **Both figures
it reports come back exactly from the committed code**, run unchanged at its registered reference
configuration (window anchor 0.1345, sensitivity 0.6674 from `PowerAtEffectSize20`, no false
positives, 6 tested elements per block and per window, 3 windows per block, ICC 0.30) on the committed
seed 20260928, whose first draw in `sweep` is exactly this call: **element detection 0.16802 and
positive windows 0.5696** (positive blocks 0.5695, as the null requires), identical to the committed
result's reference row at ratio 1.0. `scripts/clause2_design_power.py --reproduce-review` prints it;
`reproduce_the_review` and four new tests in `tests/test_clause2_design_power.py` hold it.

**Where each figure comes from**, computed in closed form (Gauss-Hermite quadrature, no simulation) by
`quadrature_of_the_committed_reference` and agreeing with the simulation to within its Monte Carlo
error (0.16824 and 0.5693):

1. **The level.** The anchor is 30 of 223 *windows* carrying a CRISPRi-tested element, and `moments`
   reads a window as positive if *any* of its six tested elements is. Six elements each detected at
   13.45% give 1 - 0.8655^6 = 57.97% positive windows on their own; the anchor was matched as an
   element rate and the simulation's windows are more than four times as often positive as the windows it was
   anchored on.
2. **The element rate is not the anchor either.** `moments` puts the logit-normal intercept's
   *location* at logit(0.1345 / 0.6674) = logit(0.2015). With a spread of 1.187 on the logit (the
   latent identity applied to 0.30), the *mean* true rate is 0.2521, not 0.2015, so the observed
   element rate is 0.6674 x 0.2521 = 16.8%.
3. **Observed and latent correlation were not kept apart.** The 0.30 is a one-way ANOVA ICC of 0/1
   calls on the benchmark's elements: an observed-scale figure. It was entered into the latent-scale
   identity, and the observed-scale correlation the simulation then produces is **0.115**, well
   under half the correlation that was measured. The registration said the realised binary-scale ICC
   would be reported beside the latent one; the committed result carries only the latent one.
4. **The test that watched the anchor looked where the definitions coincide.**
   `test_the_window_arm_reproduces_the_observed_anchor` runs `moments` at 1 element, 1 window and ICC
   0, the one configuration in which window, element and location rates are the same number. It
   passed, and the reference configuration was never checked against its anchor.

**The upper end is above the population.** `N_GRID` searches up to 5,000 compared blocks. The tier
holds 882 blocks and 531 of them carry a scored element; 1,000, 1,500, 2,000, 3,000 and 5,000 exceed
the tier, and 750 exceeds the carrying blocks, so the committed reference entry of 750 at ratio 0.75
and the 5,000 at the top of the range are sizes no experiment on this tier can have.

**What this does and does not change.** Nothing committed is altered: the 20-to-5,000 table stays as
the record of what that run computed, read as exploratory, as the correction under milestone 1.3 says.
This section establishes only that the review's two figures are the committed code's own output, and
why. The recalibrated model is registered next, before it is run. 0 AlphaGenome requests; 0.5 s.

## Pre-registration: clause 2's power model recalibrated to its own anchor, with the planned analysis simulated end to end and designs compared at equal assay cost (item 12 S2, 2026-09-28, lane-s2)

Registered before the anchor's windows are redrawn and before any configuration of this model is
simulated. The full text is `S2_REGISTRATION` in `scripts/clause2_design_power.py`, with its constants
beside it; the tests are in `tests/test_clause2_design_power.py`. The model is added beside
lane-design's, whose table (70801de) stays as the record and is read as exploratory; the new result is
`clause2_design_power_calibrated`. **0 AlphaGenome requests.** The model arm's -27.25 points and its
dispersion stay excluded, and the test that enforces it now covers this code too.

**The review's acceptance, and how each part is met.**

1. **The endpoint, as the committed estimator reads it.** Per tested element the assay calls a
   significant change in a protein-coding gene, either sign. A block is positive if any of its tested
   elements is called, a window likewise; per compared block d = Y - (positive tested windows / tested
   windows). A rate of this endpoint belongs to a unit and to how many elements it had tested, so every
   rate in the result names both.
2. **The anchor at the level it was measured.** 30/223 is a rate of positive *windows*. The windows are
   redrawn with the committed draw (`clause2_matched_control.matched_windows`, committed seed and decile
   edges, the measured arm's own verdicts) and each window's tested elements are recorded; the gate is
   the committed counts to the digit (44,081 drawn, 31,554 carrying, 223 tested, 30 regulating, 134 and
   27 blocks, 882 and 531). The calibration then requires the observed positive-window rate *over those
   windows, each with its own number of tested elements* to come back at 30/223. The anchor's
   uncertainty is a cluster bootstrap over windows linked by a shared tested element and over blocks;
   the wider interval's ends are the sensitivity grid, and the exact binomial interval is printed beside
   them as the as-if-independent one.
3. **Sensitivity and false positives, consistently.** Both act on each element, and the same s and f
   that deconvolve the anchor simulate the planned screen. s: the mean of a benchmark power column
   (reference `PowerAtEffectSize20`, 0.6674; the other four as variants). f: reference 0.01, the rate a
   false discovery rate of 5-10% implies at a call rate near 13% (an assumption, not read from the
   tables); 0 and 0.03 as variants.
4. **Observed and latent correlation, kept apart.** The benchmark ICCs are ANOVA ICCs of 0/1 calls,
   observed-scale figures. The latent spread is solved so that the calls' observed-scale ICC equals the
   measured one (the bin nearest the median carrying-block length; the others as variants), and the
   chromosome spread so that the calls' chromosome ICC equals 0.0142. Latent ICCs are reported beside
   the observed ones and never used in their place.
5. **The planned analysis, simulated end to end.** 2,000 experiments per cell. Blocks are drawn
   without replacement from the eligible population, so their chromosomes are the real ones; one
   chromosome intercept is shared by the blocks and windows on it; with g = 4 the blocks of a
   chromosome share one set of tested windows (**shared controls**). The committed estimator, the
   interval over blocks (computed exactly, the n-fold convolution that the committed 10,000-resample
   percentile interval estimates; a registered check on 200 experiments at 8 cells reports how often the
   two readings differ) and the interval over chromosomes (1,000 resamples) are computed on each. At
   the null the share excluding 0 is the **false-positive rate**; a cell is calibrated at 0.05 plus two
   Monte Carlo standard errors or less. Written in advance so it can be scored: **shared controls are
   expected to be anti-conservative**, because the committed interval treats blocks that share windows
   as independent. The probability of each registered reading is reported too: with 20 or more
   compared blocks the measured arm reads "model_failed" whenever the interval reaches 0.
6. **Designs at equal assay cost.** k tested elements per block and per window (1, 2, 3, 6: the same in
   both arms, so the null stays a null), m windows per block (1, 3, 10), g blocks per window set (1,
   4): 24 designs. Cost is counted in tested elements, N k for the blocks plus (window sets) m k for the
   windows. Designs are compared by the realised cost of the smallest searched N reaching 80%
   probability of detection, the cheapest feasible design named per ratio, and by probability of
   detection at budgets of 250 to 8,000 tested elements.
7. **The eligible population as a hard cap.** A k-element design can use only the real-unknown blocks
   holding at least k scored elements. N is never simulated above that count; a target not reached
   within it, and a budget that buys more blocks than it, are printed **"infeasible"**, never as a
   number. Every entry of lane-design's committed table above its own design's eligible population is
   listed beside it as infeasible.

**The check S2 exists for.** After calibration the anchor's own windows are *simulated* (4,000
replicates) and so are 20,000 clusters of four calls. The anchor is reproduced if the simulated
positive-window rate is within 0.5 points of 30/223 and the simulated observed-scale ICC within 0.02
of its target. If the reference calibration fails either, no size from the run is read as calibrated
and the result says so first. The implied element-level rate is set beside the benchmark's 12.15% and
14.56% as an out-of-sample comparison, not a gate.

**The block arm.** Its mean true per-element rate is ratio x the window arm's (1.0, 0.75, 0.5, 0.25,
0.1, 0.0), at the same spreads; lane-design applied the ratio to the location, which is not the mean
once there is spread. The estimand is computed in closed form and each interval's coverage of it is
reported.

**Sensitivity analysis.** One factor at a time, recalibrated and re-checked each time, at lane-design's
reference design (6, 3, 1) and at (1, 3, 1), which every carrying block can enter.

**Assumptions stated, not tested.** The anchor's tested elements were chosen by the benchmark's
designers and the planned experiment's would be chosen for the design; the window arm's per-element
rate is carried across that difference. One sensitivity per element, from pair-level power. Matched
windows holding k scored elements exist for every block (their supply is reported). Blocks and windows
share the latent spread and differ only in the mean.

**Words.** Every size is a probability-of-detection statement under the calibrated assumptions, never a
guarantee that an experiment decides clause 2.

## The result: recalibrated, the power model reproduces its 13.45% anchor, no design reaches 80% at a three-quarters rate within the eligible blocks, and shared controls are anti-conservative in every cell (item 12 S2, 2026-09-28, lane-s2, later)

The model registered in `019a9ab` was run once (1,562 s, **0 AlphaGenome requests**) and written to
`data/results/clause2_design_power_calibrated.json` with a complete manifest; the stamp is dirty only
with peers' files. The run went through a scratch driver that calls `collect_calibrated` and
`save_result` exactly as `--calibrated` does and also kept a copy of the output in case the save
failed; the manifest's argv names it. lane-design's table stays as the record.

**The gate.** The redraw reproduced the committed window draw to the digit: 44,081 windows, 31,554
carrying, 223 with a tested element, 30 with a regulating one, 134 and 27 blocks, 882 and 531.

**The anchor at its level.** 159 of the 223 windows carry one tested element and 64 carry 2 to 11.
Beneath them are 241 distinct tested elements, 23 of which regulate a coding gene (9.54%, an
element-level rate below the 13.45% window rate, as it has to be). 127 windows share a tested element
with another. Clustered, the anchor's interval is **7.88% to 19.80%** over shared-element components
(design effect 1.74, about 128 effective windows) and 8.82% to 18.81% over blocks, against the
as-if-independent 9.26% to 18.64%; the sensitivity grid used 7.88% and 19.80%.

**The calibration, latent and observed kept apart.** The median carrying block is 31.6 kb, so the
registered rule took the 25-kb benchmark ICC, **0.3745 on the observed scale**. Solved at sensitivity
0.6674 and false-positive rate 0.01: total logit spread 5.04, chromosome part 1.28. **The latent ICC
this implies is 0.885, against 0.3745 observed**; for chromosomes, 0.057 latent against 0.0142
observed. The mean true per-element rate of the window arm is 14.94%, while the rate at the
intercept's location is 0.39%: most units almost never regulate and a few often do, which is what a
0.37 observed correlation at 67% sensitivity requires. The implied observed element rate is 10.82%,
set beside 9.54% beneath the anchor's windows and the benchmark's 12.15% and 14.56%.

**The check S2 exists for: the anchor is reproduced.** Simulated, the anchor's own 223 windows come
back **13.459% positive against 13.45%** (+0.009 points; tolerance 0.5) and clusters of calls give an
**observed ICC of 0.3704 against 0.3745** (tolerance 0.02). The committed reference under the same check
gives 56.96% and 0.115. Every feasible sensitivity variant reproduces its own anchor too. One variant
has no solution: at `PowerAtEffectSize10`'s sensitivity no spread produces an observed ICC above 0.280,
so the measured 0.3745 is out of reach of an assay that weak.

**The null, simulated with its shared controls and clustering.** 210 null cells, 2,000 experiments
each, nominal 0.05. **With shared controls (g = 4) all 105 cells are anti-conservative, false-positive
rate 0.069 to 0.214**, and under the alternatives the interval covers the true difference in 65% to
93% of experiments: the registered prediction is confirmed. Without them 87 of 105 cells are calibrated; the failures sit at 20
to 30 blocks with 3 or 10 windows (up to 0.1205), and from 50 blocks up 76 of 81 are calibrated
(largest 0.066). One window per block is calibrated in all 35 of its cells (at most 0.054). The interval
over chromosomes runs 0.045 to 0.139. The exact interval over blocks and the committed 10,000-resample
one read 0 to 3 of 200 simulated experiments differently per checked cell.

**What the endpoint can show.** Even if the blocks regulated nothing (ratio 0), the expected difference
is -9.82 points with one tested element per unit and -22.02 with six; at a three-quarters rate it is
-2.46 and -4.83.

**The sample-size table at equal cost** (80% probability of detection under the calibrated
assumptions; cost in tested elements; "infeasible" = not reached within the design's eligible
population: 531 blocks for one element, 347 for two, 284 for three, 187 for six):

| Blocks regulate, relative to windows | Designs reaching 80% (of 24) | Registered cheapest (its null rate) | Cheapest whose null is calibrated |
| --- | --- | --- | --- |
| 0.75 | **0: infeasible in every design** | infeasible | infeasible; at all 531 blocks the best is 0.466 |
| 0.5 | 16 | 1 element, 3 shared windows: 400 blocks, 727 elements (0.142) | 1 element, 1 window: **all 531 blocks**, 1,062 elements, 0.809 |
| 0.25 | 24 | 1, 3 shared: 150 blocks, 290 elements (0.125) | 1, 1: 200 blocks, 400 elements, 0.824 |
| 0.1 | 24 | 1, 3 shared: 75 blocks, 158 elements (0.122) | 1, 3: 75 blocks, 300 elements, 0.852 |
| 0 | 24 | 1, 3 shared: 50 blocks, 114 elements (0.107) | 1, 1: 100 blocks, 200 elements, 0.893 |

The registered cheapest design at every ratio shares its windows, and its null is anti-conservative,
so by the registered rule its sizes read as optimistic. The last column is derived after the run from
that rule and nothing else (`cheapest_with_a_calibrated_null`, `--read-calibrated`), and is labelled so
in the code. At lane-design's reference design (6 elements, 3 windows) the committed 750 / 100 / 30 /
20 / 20 become **infeasible (0.263 at all 187 eligible blocks) / 187, all of them / 75 / 50 / 30**.

**At fixed budgets.** With 250, 500, 1,000 and 2,000 tested elements the best design with its own
windows detects a half-rate difference with probability 0.283, 0.482, 0.788 and 0.942. At 4,000 and
8,000 elements 19 and 22 of the 24 designs are infeasible, because the budget buys more blocks than the
tier holds. No design at any budget or size detects a three-quarters rate with probability above 0.47.

**The registered reading, simulated.** With 20 or more compared blocks the measured arm reads
"model_failed" whenever the interval over blocks reaches 0. At a three-quarters rate, with every
eligible block tested, that reading comes out in 53% to 82% of experiments depending on the design: the
rule reads an interval that could not decide as the blocks regulating no less often than their windows.

**Sensitivity, one factor at a time.** At a three-quarters rate both reference designs are infeasible
in every variant. At a half, one element with three windows needs 200 blocks (anchor at 19.80%) to 531
(false-positive rate 0.03), and is infeasible at the anchor's low end; six elements with three windows
need 150 to 187 or are infeasible (anchor's low end; false-positive rate 0.03).

**lane-design's table against the eligible population.** 50 of its committed entries exceed their own
design's eligible population, among them the reference 750 at 0.75 against 187 and the 5,000 at the top
of that ratio's range; they are listed in the result as infeasible and the committed file is unchanged.

**What this settles and what it does not.** It settles that the power model now reproduces the quantity
it is anchored on, at the level the anchor was measured, with its correlation on the right scale, and
that the planned analysis has been simulated under the null with the structure a real experiment would
have. It shows that the committed interval cannot be used with shared controls, and that, under these
assumptions, no searched design on this tier detects blocks regulating three quarters as often as their
windows with probability above 0.47. It does not estimate the
effect, which remains the reason every size is conditional, and the window arm's per-element rate is
carried from elements a benchmark chose to test. Every size here is a probability of detection under
these assumptions, never a guarantee that an experiment decides clause 2. Clause 2 stays not met.

## Item 13 C4 registered: hide evidence and predict it, with the AlphaGenome ablation run first (2026-09-28, lane-c4)

This section registers the measuring stick for the coherence programme (ROADMAP section 5, item 13):
the leave-one-source-out harness that C1, C2 and C3 will all be scored by, and the diagnostic that
runs before it, **the AlphaGenome ablation**. Both are fixed here before either is run. The constants
are in `genomeos/attribution/ablation.py` and `genomeos/attribution/holdout.py`; the runners are
`scripts/c4_ablation.py` and `scripts/c4_holdout_run.py`. Nothing calls the model (0 requests). The
per-element response cache is never opened: the model's labels are read from the compiled programs,
which are what the project states.

**What was read before this registration, so the reader can judge it.** While designing the census the
lane read the published values of the six headline results in `manifest_headlines.json`, including the
measured arm of the node-containment result (+5.893 points, chromosome interval [3.183, 8.458]) and the
`cannot_decide` reading of clause 2's measured arm. It also counted the units each source holds per
endpoint (below), which decides which endpoints clear the floors, and it parsed the compiled programs
once to confirm that they hold 440,377 predicted links, the published count. No metric of any labelling
on any source had been computed when this section was written.

**What both may and may not be called.** The harness is an **internal development benchmark only**.
Every source in it has been read by this project before, several of them repeatedly. The ENCODE CRISPRi
held-out file alone is read by at least ten scored results (lane-split), so it is a reused benchmark
and cannot become fresh validation. The ablation reads the same files. A score may be called: "on
evidence this project has already read, labelling L predicts held-out source S's own endpoint with
metric M [95% locus-bootstrap interval] at coverage C, beside distance to TSS and the unchanged labels."
It may not be called any of the following:

- an external, independent or fresh validation;
- evidence of biological accuracy beyond the sources held here;
- a test of regulation, when the endpoint is reporter activity, in-vivo activity, base sensitivity or
  association (activity is not regulation, and association is not perturbation);
- a comparison between sources or endpoints;
- out of sample for anything tuned on the CRISPRi benchmark;
- a verdict on what the evidence can predict when `rest` scores low. `rest` is an unfitted rule, so a low
  score says only that this rule transfers poorly.

### Part 1, run first: the AlphaGenome ablation

**Removed.** The all-element deletion sweep (`enhancer_targets_all_<chrom>`) and its per-element response
cache; the sample target runs of the same model (`constrained_targets_*`, `enhancer_targets_*`); the
compiled predicted layer (every `element` block without `_measured`, every `rule` with `evidence:
predicted`); and any value computed from these. An input path containing `alphagenome`,
`enhancer_targets`, `constrained_targets` or `knowledge/compiled/` counts as model output.

**Censused.** There are three families.

1. The six headline results rebuilt in `manifest_headlines.json`.
2. Every element block of the 24 compiled programs, predicted and measured, axis value by axis value
   (`class`, `targets`, `origin`, `molecular_role`, `activity`, `target_relation`, `evidence_status`),
   with the region blocks counted.
3. Every compiled rule, split by evidence kind after R1.

**The classification rule.** Each conclusion is put in exactly one of three classes.

- **Was never model-dependent.** No ablated input enters the stated value or selects the units it is
  computed on. A model file read only for a quantity the conclusion does not state does not count, and it
  is named. Examples are chromosome names and lengths, or a column that the stated figure does not use.
  For a per-element label, the label's *value* is classified. The fact that an element carries a block at
  all, because the model named a target, is counted separately as placement.
- **Survives.** An ablated input enters the conclusion, and with every ablated input removed the same
  statement is still established on experimental evidence alone.
  - A predicted link (element, gene, direction) survives when a CRISPRi pair on the same element meets
    `measured.RECIPROCAL_OVERLAP` (0.5), names the same gene, and is significant in the stated
    direction: activates means `significant_decrease`, inhibits means `significant_increase`. Any cell
    and either benchmark file counts.
  - A compiled rule survives on the same test, measured in the rule's own `when` cell (R1).
  - A result survives when the record holds an experimental arm of the same statistic, computed without
    ablated input, with the model figure's sign and a 95% interval that excludes zero. The interval is
    the result's own headline interval. The surviving figure is the arm's, never the model's.
- **Does not survive.** Everything else, with one reason:
  - `about_the_model`: the statement's subject is a model output (its value, its gain or its agreement);
  - `contradicted`: the same endpoint measured the same element and gene, and found either a
    well-powered null or a significant change in the opposite direction;
  - `inconclusive`: the same endpoint measured it without deciding. This covers an underpowered null, a
    missing effect, and an experimental arm that is below its floor or has an interval across zero;
  - `never_measured`: nothing of the same endpoint was measured. The element and gene were never
    perturbed together (for a rule: not in its cell), or the result has no experimental arm.

Reporter activity, VISTA, saturation mutagenesis and GTEx association are other endpoints. They never
make a link survive (review S4). Where they touch a link that does not survive, they are counted beside
it.

On a predicted block, the axis values that carry the link take the link's class. These are `targets`,
`activates_target` or `represses_target`, `predicted_deletion_target`, `nearest_tss_in_domain`, and the
repression-alternatives role group. `predicted_model` (on predicted blocks) and `conflicting` (on
measured twins) are `about_the_model`. Every other value is sequence, registry or assay content and was
never model-dependent. The same holds for a measured twin's activity and relation. Placement is counted
from the measured side: of the CRISPRi significant decreases (element, gene, cell), how many have a
model-compiled element that meets the overlap rule, and, for the rest, how many have an ENCODE cCRE that
would carry them with no model involved.

**The registered reading of each headline result.** Each reading comes from its writer, read before the
run, and is checked at run time against the result's manifest inputs.

| result | model input | reading |
| --- | --- | --- |
| node_containment_audit (+2.90) | value: the archive's most-moved coding gene | experimental arm `node_containment_measured`, uniform control, chromosome interval |
| node_containment_measured (+5.89) | bookkeeping only: stage 1's rows supply chromosome name and length | never model-dependent |
| constrained_unknown_targets (62.3% vs 86.0%) | value: the model's mover | experimental arm `clause2_measured_arm`, its registered `primary_reading` |
| therapeutic_benchmark (9/9) | none (checked against the manifest) | never model-dependent, unless the manifest names a model input |
| unknown_coverage (0.52%) | a column the figure does not state (organise's `attributed_elements`) | never model-dependent |
| crispri_published (deletion feature over ENCODE-rE2G) | value: the deletion feature | `about_the_model` |

### Part 2: the leave-one-source-out harness

**Sources and endpoints.** Each assay is a separately holdable unit, and each endpoint is modelled apart.

| source | unit | endpoint(s) | primary metric |
| --- | --- | --- | --- |
| `crispri:<study>` (the benchmark's `Dataset`, R5; 11 studies) | element-gene-cell pair | `decrease`: significant decrease (1) against well-powered null (0); `increase`: significant increase (1) against well-powered null (0). Underpowered nulls and missing effects are counted, never scored (R2) | average precision |
| `lentimpra:<cell>` (K562, HepG2, WTC11) | one tile in one cell (R6, never aggregated) | `activity`: log2(RNA/DNA); `active`: at or above 1.0 | Spearman; average precision |
| `vista` | element | positive in any tissue (1) or negative (0) | average precision |
| `satmut` | measured base of each locus's primary experiment | functional (1) or inert (0) | average precision |
| `gtex` | element x protein-coding gene with a TSS within 1 Mb, at the elements GTEx ties to a gene | associated in any of 49 tissues (1) or not (0; untested and null are not told apart) | average precision |

AUROC is the secondary metric everywhere (for `activity`, the AUROC of `active`). The prevalence is
reported beside average precision. The GTEx elements are the ones the eQTL distillation indexed, which
are elements that an earlier target run predicted for, so that unit set was chosen by those runs.

**Unit counts, read before this registration.**

- **CRISPRi decrease, positives/negatives:**
  - Gasperini2019: 360/4,893
  - Schraivogel2020: 23/1,276
  - Xie: 42/393
  - Morris: 35/132
  - K562_DC_TAP: 12/1,221
  - WTC11_DC_TAP: 15/1,886
  - HCT116: 34/0
  - Nasser2021: 111/0
  - Klann: 21/0
  - Reilly: 8/0
- **CRISPRi increase:**
  - Gasperini2019: 42/4,893
  - WTC11_DC_TAP: 20/1,886
  - every other study is under 20 positives or has no informative negatives
- **lentiMPRA:** about 54,000 tiles per cell, and 2,864 loci
- **VISTA:** 1,267/1,175, and 948 loci
- **saturation mutagenesis:** 2,694/7,140 bases, and 19 loci
- **GTEx:** 21,166/124,183, and 421 loci. Of 30,210 element-gene associations, 8,333 are with a gene that
  is not protein-coding and 711 with a gene absent from GENCODE v50. Both groups are counted and never
  scored.

**The split.** The split is by study and by locus. Holding out S does three things.

1. It removes S and every provenance sibling of S. CRISPRi studies that share a value of the
   benchmark's `Reference` column are siblings, which makes K562_DC_TAP and WTC11_DC_TAP siblings
   (Ray et al. 2025). The three lentiMPRA cells are one library.
2. It masks every record of S's endpoint family whose interval shares a base with an S unit. These are
   the related intervals of lane-split's `split_overlap`, so an element measured by two studies cannot
   leak across the split. The relation counts per source reuse `measured.split_overlap` itself.
3. It never admits a pair from the CRISPRi benchmark's held-out file as evidence (R5), whatever is held
   out.

Records of another family are evidence in full: a reporter predicting a screen is the question, and a
second screen of the same pair is a replicate. A labelling that declares it read S, a sibling, the
held-out file, or a same-family source outside `evidence(without=S)` is refused (`LeakError`).

**Metric and interval.**

- **Interval.** A 95% percentile bootstrap over loci, never over pairs. The run uses 1,000 resamples,
  seed 20260928, and the same resamples for every labelling of a source.
- **Loci.** A locus is a connected component of a source's units, joined by a shared 1 Mb bin of the
  midpoint or, for pair endpoints, a shared gene. Each resample draws as many loci as the source has.
- **Ranks under resampling.** Spearman under resampling keeps the full-sample ranks, weighted.
- **Abstention.** An abstention ranks below every stated score, and coverage is reported beside every
  value.
- **Floors.** An endpoint is described and not scored below 20 positives or 20 negatives
  (`MIN_FOR_A_COMPARISON`), below 20 units for `activity`, or below 10 loci. By the counts above, 15
  endpoints are scored and 16 CRISPRi endpoints are described:
  - HCT116, Nasser2021, Klann and Reilly have no well-powered null at all;
  - K562_DC_TAP is under 20 positives on both endpoints;
  - WTC11_DC_TAP is under 20 positives on `decrease`;
  - Morris, Schraivogel2020 and Xie are under 20 positives on `increase`.

**The labellings of the first run.** Three labellings are scored for every source.

- `distance`: minus the distance from the unit's midpoint to its gene's TSS. CRISPRi uses the
  benchmark's own TSS columns, and GTEx uses GENCODE v50. For element and base units it is the nearest
  GENCODE v50 protein-coding TSS.
- `unchanged`: the compiled predicted layer as it stands, matched by one shared base. For pair endpoints
  it is the largest strength among matched links that name the unit's gene in the endpoint's direction:
  activates for `decrease`, inhibits for `increase`, and either for `associated`. It is 0 when the
  matched links name another gene or direction. For element and base endpoints it is the largest strength
  of any matched link. It abstains where no predicted element shares a base, and it never reads the
  measured twins.
- `rest`: the prediction from the remaining sources, an unfitted rule. The score is 10 x G + E +
  closeness.
  - G counts the view's sources that support this gene at this interval: a CRISPRi training pair on the
    same gene, significant in the endpoint's direction, or a GTEx association with the same gene.
  - E counts the sources that say the interval does something: lentiMPRA active by R6's rule in some cell,
    a VISTA positive, or a functional saturation-mutagenesis base. For non-regulation endpoints it also
    counts a significant CRISPRi training pair and a GTEx association on any gene.
  - Closeness is 1 / (1 + distance / 100 kb).

The three paired differences are also scored on the same resamples: rest minus distance, rest minus
unchanged, and unchanged minus distance.

**The call the pilot makes.** It is `holdout.score(labels, held_out_source)`, and `compare(a, b,
source)` gives the paired difference. It is deterministic, with a fixed seed and resample count. It is
cached in memory and under `data/cache/holdout`, keyed by the predictions themselves. It runs on
ordinary CPUs, uses numpy over at most about 145,000 units, and makes no request.

**Expectation, written before the run.** Distance should be the strongest simple predictor of CRISPRi
decreases, as the benchmark literature finds. `rest` should add little there, because after the locus
mask no screen of the same pair remains. The unchanged labels should have low coverage of the screens'
pairs. None of this changes a rule above.

*Correction, 2026-09-28, lane-c4, after the run.* The prose above miscounts the sources, and the rules
are unaffected. The benchmark holds **10** CRISPRi studies, not 11: Nasser2021 spans both files and is one
study. So there are **16** holdable sources, not the 17 of the registration commit's subject, and **14**
CRISPRi endpoints are described without a score, not 16. The 15 scored endpoints were counted correctly.
The error was in the counting prose, not in any constant.

## The result: without the model, 93 of 440,377 predicted links survive; in the first held-out scores, neither the unfitted rest rule nor the unchanged labels beats distance to TSS on any CRISPRi decrease endpoint (2026-09-28, lane-c4, later)

Both runs were made once, as registered in the section above, in a clean worktree at `cc824f4` (the
revision stamp is not dirty). Results: `data/results/c4_alphagenome_ablation.json` and
`data/results/c4_holdout_scores.json`. Nothing called the model (0 requests), and the per-element
response cache was never opened. The ablation took 9 s and the harness 52 s, both on one laptop CPU.
**Every score below is an internal development benchmark reading.** Each is on evidence this project
has already read, and none is a fresh or external validation. The negatives come first.

### Negatives

1. **The prediction from the remaining sources (`rest`, an unfitted rule) is below distance to TSS on
   all four scored CRISPRi decrease endpoints and on GTEx.** Every one of these intervals excludes zero.
   The paired differences in average precision are:

   | held-out source | rest minus distance |
   | --- | --- |
   | Gasperini2019 | -0.278 [-0.331, -0.221] |
   | Morris | -0.305 [-0.429, -0.091] |
   | Schraivogel2020 | -0.309 [-0.608, -0.140] |
   | Xie | -0.344 [-0.440, -0.210] |
   | GTEx | -0.155 [-0.171, -0.141] |

   As registered, this says the rule transfers poorly. It says nothing about what the evidence can
   predict. The rule puts element-level activity (an active reporter tile, a VISTA positive) ahead of
   distance, so every gene paired with an active element moves up together. Most of those genes are that
   element's negatives.
2. **The unchanged labels beat distance on no CRISPRi decrease endpoint.** The compiled predicted layer
   names one gene per element. It is below distance on three endpoints, each with an interval excluding
   zero:
   - Morris decrease: -0.366 [-0.522, -0.187];
   - GTEx: -0.232 [-0.243, -0.220];
   - WTC11_DC_TAP increase: -0.013 [-0.033, -0.004].

   On Gasperini2019 (-0.068 [-0.140, +0.003]), Schraivogel2020 (-0.229 [-0.584, +0.072]) and Xie (-0.145
   [-0.333, +0.071]) the interval spans zero. Its coverage of the screens' pairs is 45% to 77%. This is
   **not** the deletion feature that `crispri_published` scored. That feature read every gene in the
   scorer's window from the response cache; the compiled labels keep one target per element. The two
   results do not contradict each other, and this one does not re-test that one.
3. **Nothing scored here predicts reporter activity, in-vivo activity or base sensitivity well.**
   - lentiMPRA: the AUROC of `active` is 0.503 to 0.561 across the three cells and three labellings. The
     largest Spearman is 0.126 [0.117, 0.135] (unchanged, HepG2).
   - VISTA: AUROC 0.529 to 0.571.
   - Saturation mutagenesis: 19 loci. Every paired interval spans zero, and the best AUROC is distance's,
     0.646 [0.524, 0.738].
4. **14 of the 20 CRISPRi endpoints cannot be scored.**
   - HCT116, Nasser2021, Klann and Reilly hold no well-powered null at all.
   - K562_DC_TAP is under 20 positives on both endpoints.
   - WTC11_DC_TAP is under 20 positives on `decrease`.
   - Morris, Schraivogel2020 and Xie are under 20 positives on `increase`.

   The benchmark's held-out file therefore contributes only three scored endpoints: Morris and Xie
   `decrease`, and WTC11_DC_TAP `increase`.
5. **Ablation: the predicted target layer does not survive the removal of the model.**
   - Links (element, gene, direction): of 440,377 predicted links, **93 survive** (0.021%). A same-endpoint
     measurement establishes each of these 93. The other 440,284 do not survive:
     - 440,249 were never measured by the same endpoint. For 438,872 of these the element was never
       screened; for 1,377 it was screened, but the predicted gene was not tested.
     - 28 are contradicted: 23 by a well-powered null, and 5 by a significant change in the other
       direction.
     - 7 are inconclusive (underpowered nulls).
   - Where a screen tested the predicted gene on the element, 93 of 128 links survive.
   - Rules: of the 440,377 predicted rules, **39 survive in their own cell**. The other 440,338 were never
     measured in the rule's cell, and none is contradicted there.
   - Other endpoints: 6,434 of the never-measured links have evidence from another endpoint on the element.
     This is 3,216 active lentiMPRA readings, 3,202 GTEx associations with the same gene, 69 VISTA
     positives and 14 functional saturation-mutagenesis bases. None of it makes a link survive (S4).
   - Headlines: 2 of the 6 do not survive.
     - `constrained_unknown_targets` (62.3% against 86.0%) is inconclusive: its measured arm reads
       `cannot_decide`.
     - `crispri_published` is about the model.
6. **The model's element selection kept 449 of 661 measured regulatory links out of the executable
   programs.** These are the CRISPRi significant decreases, counted as element, gene and cell. Only 212
   of the 661 are carried, because a twin is written only beside a model-compiled element. The 212 are
   exactly the 212 experimental rules: 169 from training and 43 marked held-out. Of the 449 left out, 62
   have an ENCODE cCRE that meets the same overlap rule with no model involved, and 387 have none.

### What survives on experimental evidence alone

- **Headline results, 4 of 6.**
  - `node_containment_audit` survives, on its experimental arm only. The surviving figure is **+5.893
    points [3.183, 8.458]** on 661 measured pairs against the uniform control. It is not the model's
    +2.90.
  - `node_containment_measured` (+5.89), `therapeutic_benchmark` (9/9) and `unknown_coverage` (0.52%)
    were never model-dependent. The first reads model files only for chromosome names and lengths, and the
    third only for a column that its figure does not state.
- **Per-element labels.** Every sequence and registry label was never model-dependent. That is origin
  (440,377), class (440,377), the registry roles (633,404 values), and the registry and selection status
  values (880,754). All content of the 19,072 measured twins was never model-dependent either. The model's
  own `predicted_model` status (440,377) and the twins' `conflicting` status (14,537) are about the model.
- **Compiled rules.** The 212 experimental rules were never model-dependent in value. Their placement was
  model-selected, as item 6 above says.

### The first held-out scores beside the baselines

Average precision is given with its 95% locus-bootstrap interval, and the prevalence is in brackets
after the source. For lentiMPRA the value is the Spearman correlation of `activity`.

| held-out source (endpoint) | rest | distance | unchanged (coverage) |
| --- | --- | --- | --- |
| Gasperini2019 decrease (0.069) | 0.215 [0.178, 0.266] | 0.493 [0.439, 0.560] | 0.425 [0.367, 0.485] (0.55) |
| Morris decrease (0.210) | 0.447 [0.306, 0.669] | 0.752 [0.584, 0.883] | 0.386 [0.255, 0.521] (0.77) |
| Schraivogel2020 decrease (0.018; 12 loci) | 0.082 [0.029, 0.136] | 0.390 [0.203, 0.672] | 0.162 [0.018, 0.316] (0.45) |
| Xie decrease (0.097) | 0.245 [0.153, 0.384] | 0.589 [0.407, 0.747] | 0.444 [0.297, 0.587] (0.72) |
| Gasperini2019 increase (0.009) | 0.011 [0.007, 0.022] | 0.011 [0.007, 0.017] | 0.032 [0.006, 0.076] (0.54) |
| WTC11_DC_TAP increase (0.011; 22 loci) | 0.018 [0.009, 0.049] | 0.023 [0.011, 0.049] | 0.010 [0.005, 0.017] (0.65) |
| lentiMPRA K562 activity (Spearman) | 0.061 [0.051, 0.070] | 0.061 [0.052, 0.070] | 0.068 [0.058, 0.077] (0.40) |
| lentiMPRA HepG2 activity (Spearman) | 0.039 [0.029, 0.048] | 0.040 [0.031, 0.050] | 0.126 [0.117, 0.135] (0.40) |
| lentiMPRA WTC11 activity (Spearman) | 0.036 [0.026, 0.047] | 0.038 [0.027, 0.049] | 0.072 [0.063, 0.082] (0.40) |
| VISTA positive (0.519) | 0.550 [0.518, 0.587] | 0.543 [0.514, 0.575] | 0.574 [0.545, 0.605] (0.42) |
| saturation mutagenesis functional (0.274; 19 loci) | 0.370 [0.206, 0.531] | 0.380 [0.277, 0.497] | 0.274 [0.182, 0.448] (0.31) |
| GTEx associated (0.146) | 0.281 [0.263, 0.297] | 0.436 [0.417, 0.456] | 0.204 [0.189, 0.220] (0.59) |

**The positives, in their registered wording.** On evidence this project has already read, the unchanged
labels predict reporter activity slightly better than distance to TSS, and the intervals exclude zero:
- HepG2 activity: Spearman +0.086 [+0.076, +0.097];
- WTC11 activity: +0.034 [+0.024, +0.045];
- VISTA: average precision +0.031 [+0.010, +0.050].

K562 spans zero. These are activity endpoints: they are not a test of regulation, and the effects are
small next to the endpoints' own spread.

**The split did what it was registered to do.** Holding out Gasperini2019 masked 125 Nasser2021 and 76
Schraivogel2020 training pairs at related intervals. By `measured.split_overlap`, 357 of Gasperini2019's
5,299 pairs are related to another study's pair, and none is an identical pair. The three lentiMPRA
cells measure the same tiles (53,988 to 53,990 per cell), so each removes the other two entirely. The 4,378 pairs of the
held-out file were never evidence for any held-out source.

**What it may be called.** It is an internal development benchmark, and a measuring stick. The pilot
(C1 with C2 and C3) is scored against `unchanged` and `distance` with `holdout.compare`, on the same
resamples. **What it may not be called:** a validation of anything, or a finding that the sources cannot
predict one another. `rest` is one unfitted rule, and it is now the reference labelling to beat, not a
ceiling.

## Pre-registration: the measured arm's reading rule revisited, with the error rates a rule must meet fixed before any rule is scored (item 12 S2 follow-up, 2026-09-29, lane-rule)

c17eedc's rule reads "model_failed" whenever the 95% interval over blocks reaches 0 with at least 20
compared blocks, and its text says the blocks regulate coding genes "at a rate not below their matched
windows". cab70d9 simulated it: at a true three-quarters rate, with every eligible block tested, it reads
so in 53% to 82% of experiments, so a failure to detect is read as a finding. Its "wording_wrong" says the
blocks hold no element that regulates a coding gene, and it is read whenever the interval lies below 0,
which a half or a three-quarters rate also produces. The committed rule stays the record
(`measured_reading`, unchanged). The revised one is `revised_reading` beside it, registered in
`READING_RULE_2026_09_29` in `scripts/clause2_measured_arm.py` before any probability of it, or of any
rule other than c17eedc's, was computed.

**What each reading claims.** model_failed: a true ratio of 1 (the blocks' mean true per-element rate
over their windows'). wording_wrong: a true ratio of 0. cannot_decide: nothing.

**The error rates a rule must meet.** model_failed in at most 5% of experiments at every true ratio of
0.75 or below: the null of this reading is that the blocks regulate less than their windows, and 0.75 is
the smallest departure on the grid. wording_wrong in at most 5% at every true ratio of 0.1 or above,
since a tenth of the windows' rate is not "no element". A rate meets its bound within 0.05 + 2 Monte
Carlo standard errors (cab70d9's null rule), in every (design, size) cell where the rule may read; a cell
that fails either bound reads cannot_decide. cannot_decide is never an error, and its cost is reported as
the probability of deciding. The rule **decides** clause 2 at a design and size only if it also reads
model_failed with probability at least 0.8 at a true ratio of 1 and wording_wrong with probability at
least 0.8 at a true ratio of 0.

**The rule.** The committed estimator, interval and 20-block floor are unchanged. Only the lines move.
model_failed if the interval lies wholly above the expected difference at a three-quarters rate (-2.456,
-3.589, -4.168 and -4.832 points for 1, 2, 3 and 6 tested elements per unit). wording_wrong if it lies
wholly below the expected difference at a tenth (-8.84, -13.267, -15.757, -19.175). cannot_decide
otherwise. The lines come from cab70d9's reference calibration and are never re-estimated from an
experiment's own data. An arm whose units carry varying numbers of tested elements, as c17eedc's does,
has no line. No other line is tried. Reported beside the reading, and never as one: the registered
ratios the interval excludes.

**The best any rule on this estimator could do**, reported beside the rule: the Neyman-Pearson threshold
on the estimator's simulated distribution, and the same on the blocks' own rate as if the windows' rate
were known exactly. Descriptively, the same at the weaker margins (0.5 for model_failed; 0.25 and 0.5
for wording_wrong), to show which weaker claims the arm could support.

**The simulation.** cab70d9's, unchanged. The calibration is recomputed through the same redraw and gated
against the committed one. 24 designs, every size up to each design's cap, every ratio, 2,000
experiments per cell, seed 2026092913; the committed and the revised readings are computed on the same
experiments. Robustness, descriptive only: the lines stay at the reference while the truth comes from
each feasible one-factor variant, at the caps of designs (1, 10, 1), (1, 3, 1) and (6, 3, 1).

**Expected before running.** Both bounds are met in every own-window cell from 50 blocks up; some cells
at 20 to 30 blocks and some with shared controls fail. No design decides. The revised model_failed at
ratio 1 is about 0.5 at best, and wording_wrong at ratio 0 about 0.2 at best. The best rule on the
estimator also stays below 0.8 (about 0.6 and 0.3). The bounds do not survive every variant: an anchor at
its low end, or no false positives in truth, puts a reading's error above 5%. If no design decides, the
arm reports its coverage first, then its estimate and interval with the ratios the interval excludes,
then the revised reading where the interval clears a line, with its probabilities. It never reports
"model_failed" from an interval that only reaches 0. Every probability is a probability under cab70d9's
calibrated assumptions, never a guarantee. 0 AlphaGenome requests.

## S4 registered: what a correct attribution means, one table of what each observation can establish, refute or say nothing about, and a scorer that keeps target accuracy, role accuracy and coverage apart (item 12 S4, 2026-09-29, lane-s4)

The second external review (ROADMAP section 5, item 12, S4) asked the project to **define what a correct
attribution means** before it scores one. This section is that definition, fixed before the scorer was
applied to any real labelling. The constants are in `genomeos/attribution/correctness.py`; the runner is
`scripts/s4_correctness_run.py`; the tests are `tests/test_attribution_correctness.py`. Nothing calls the
model (0 requests).

**What was read before this registration.** The ablation's published counts (C4, above): 93 of the 128
measured links survive, 28 are contradicted (23 by a well-powered null, 5 by the opposite sign), 7 are
inconclusive, and 39 rules survive in their own cell with none contradicted there. R2's five CRISPRi
outcomes and R6's reporter rule (`measured.py`), R7's five axes (`lang/grammar.py`), the C4 harness
(`holdout.py`) and the C5 probe's linkage counts. The runner was exercised once on a synthetic program
against the real held-out units, to check that it runs; no compiled label was read.

### The definition

**Five claims, kept apart.** An attribution is correct or not claim by claim, and a claim is about one
of five things:

| Axis | What the claim says | Where the compiled labels state it |
| --- | --- | --- |
| origin | where the sequence came from | R7 `origin` |
| molecular role | what the element is biochemically | R7 `molecular_role` |
| activity | the direction of its effect on its target in place, or its activity in a reporter | R7 `activity`, and the rule's action |
| target | which gene it acts on | the gene a rule names (R7 `target_relation` says how it was named) |
| cell context | the cell in which the relation holds | the rule's `when: cell_type` (R1) |

`evidence_status` states what a claim rests on (a registry, a model, a measurement, constraint). It is
not a claim about the biology and is never judged.

**Not naming a coding gene is not absence of function.** No axis has a value meaning "no function".
`unknown`, `unassigned` and a group of unchosen alternatives (`silencer|insulator_like|...`) are not
claims: they are counted beside, and never judged. A well-powered CRISPRi null refutes one named gene in
one cell for an effect of 20% or more, and nothing else.

**Five verdicts.** A claim is *correct* when an observation the table allows establishes it and none
refutes it; *incorrect* in the reverse case; *unresolved* when both are present; *the observation model is
inadequate* when the observations exist but the assay's own reading rule cannot produce a label (a
lentiMPRA cell whose tiles split evenly under R6's share rule; one element, gene and cell significant in
both directions); and *not judged* otherwise, with one reason. **Competing explanations are kept, never
resolved by the scorer.** "The observation model is inadequate" stays live beside every verdict an
observation decides. An unresolved claim keeps four: the claim holds and the refuting observation is
wrong; the claim is wrong and the establishing one is wrong; the answer differs between studies,
conditions or contexts; the observation model is inadequate. A refuted target also keeps "a redundant
element compensates, which only a combinatorial perturbation separates".

### The table: what each observation can say about each axis

**E** establishes, **R** refutes, **M** the observation model is inadequate (the assay's rule cannot read a
label), **S** suggests only (counted beside, never decides a verdict), **—** cannot establish (says
nothing, and is never consulted to decide one). *Not loaded* means no holdout source supplies it, so it
judges nothing in the run registered below.

| Observation | Origin | Molecular role | Activity | Target | Cell context |
| --- | --- | --- | --- | --- | --- |
| CRISPRi, significant decrease (R2) | — | S (enhancer or promoter) | E activates; R represses, no effect | E, any cell | E, the cell screened |
| CRISPRi, significant increase (R2) | — | S (silencer, insulator, competing promoter; never chosen) | E represses; R activates, no effect | E, any cell | E, the cell screened |
| CRISPRi, well-powered null (R2) | — | — | E `no_effect_measured` only; never a direction | R, only in the cell screened | R, when the gene responds in another cell |
| CRISPRi, underpowered null (R2) | — | — | — | — | — |
| CRISPRi, missing (R2) | — | — | — | — | — |
| Two elements silenced together, a change (*not loaded*) | — | — | E a direction; R no effect | E | E |
| Two elements silenced together, a null (*not loaded*) | — | — | E no effect | R, removing redundancy as the explanation of a single null | — |
| lentiMPRA, one cell, active by R6 | — | S | E active in reporter; R inactive; never a direction | **—** | S |
| lentiMPRA, one cell, silent by R6 | — | — | E inactive in reporter; R active (that cell) | — | — |
| lentiMPRA, one cell, tiles split evenly | — | — | **M** | — | — |
| VISTA, positive (mouse e11.5) | — | S | E active in reporter; R inactive | — | — |
| VISTA, negative | — | — | E inactive in reporter; R active (e11.5 only) | — | — |
| Saturation mutagenesis, a functional base | — | S | S (the bases, not the element) | — | — |
| Saturation mutagenesis, every base inert | — | — | — | — | — |
| GTEx eQTL, associated | — | — | S (association, not perturbation) | S (linkage not excluded) | — (the units pool 49 tissues) |
| GTEx, not associated | — | — | — | — | — |
| Chromatin contact present (*not loaded*) | — | — | — | S | S |
| Chromatin contact absent (*not loaded*) | — | — | — | — | — |
| Conservation, constrained (*not loaded*; the labels' input) | — | — | — | — | — |
| Conservation, not detected (*not loaded*) | — | — | — | — | — |
| Allele-specific readout, imbalance (C5, *not loaded*) | — | — | S | S (a haplotype, not an element) | S |
| Allele-specific readout, balance (C5, *not loaded*) | — | — | — | — | — |
| AlphaGenome deletion prediction (a model output) | — | — | — | — | — |
| Sequence annotation (*not loaded*; the labels' input) | E / R | — | — | — | — |
| ENCODE cCRE registry (*not loaded*; the labels' input) | — | E / R (promoter-, enhancer-, insulator-like, open chromatin) | — | — | — |

**The cells that cannot establish, stated plainly.**

- A reporter (a lentiMPRA tile, a VISTA embryo, a saturation-mutagenesis base) can establish activity in
  a reporter. It can never establish a target, and never a direction of effect in place: activity is not
  regulation.
- CRISPRi can establish a target, a direction and a cell context. It can never establish a biochemical
  role (an increase is not a silencer) or an origin.
- A well-powered null refutes one gene in one cell. It never establishes absence of function, and it
  says nothing about direction.
- Conservation can suggest selection, which is recorded on `evidence_status`. It never establishes a
  role, an activity, a target, a context or an origin; an exapted repeat can be constrained.
- Contact and association (GTEx, an allelic imbalance) can suggest a target and never establish one. A
  within-person imbalance cannot be credited to one element when a median 1,580 to 2,012 heterozygous
  SNVs lie within 1 Mb of the TSS (C5).
- A model output is a claim to be judged, never an observation that judges one.
- An annotation or the registry establishes origin or signature only for a labelling that did not read
  it. The compiled labels read both, so neither is their check.
- Nothing in the table establishes a silencer, a competing promoter, a structural role or a coding
  candidate.
- **Combinatorial perturbation** is the route when single perturbations cannot separate additive,
  redundant and cooperative explanations. A joint null removes redundancy between the two as the
  explanation of a single null (C6, phase C). No paired perturbation is among the cached screens.

### How each axis is judged

- **Target.** Established by a significant change of the named gene on silencing the element, in any
  cell. Refuted by a well-powered null of that gene in the cell the claim states, and only when the gene
  responds in no cell. A response in another cell with a null in the stated cell leaves the target
  established and refutes the context. A response and a null both in the stated cell is unresolved. A
  null in another cell is not a refutation; it is counted (`well_powered_null_only_in_another_context`).
  A claim with no stated cell can be established and never refuted. `refutable` counts the claims whose
  stated cell was screened on the named gene with an outcome that can establish or refute.
- **Cell context.** Judged only where the gene responds to the element in some cell. Established by a
  response in the stated cell; refuted by a well-powered null there. A claim whose target never responded
  is not judged on context, so one failure is never counted on both axes.
- **Activity.** The direction (activates, represses) is judged only where the named gene responded: in
  the stated cell when it responded there, otherwise in every cell it responded in. Both directions in
  one cell is *the observation model is inadequate*; different directions in different cells are
  *unresolved*. A null is read on the target axis, never as a wrong direction. Reporter values are judged
  by reporter observations in the stated cell when one was measured there, otherwise as R7 defines them
  (active in at least one reporter context; inactive in every one).
- **Molecular role.** Judged only by an observation the table lets establish a biochemical signature.
  Activity in a reporter or in place suggests, and never decides.
- **Origin.** Judged only by a sequence annotation the labelling did not read.

**Matching.** An observation is of the element when its interval meets the measured layer's rule
(`measured.measures`, reciprocal overlap at least 0.5) for CRISPRi pairs, lentiMPRA tiles and VISTA
elements. A saturation-mutagenesis base must lie inside the element, and a GTEx unit's element must share
a base with it. lentiMPRA is read per cell from the matched tiles by R6's share rule. Genes are matched by
symbol, and cells after lowercasing and dropping every character that is not a letter or a digit (the
ablation's rule).

### The three quantities, never combined

`judge(claims, labels)` returns a report with three quantities, each per axis, each a `Share` with its own
numerator and denominator:

- **target accuracy**: target claims established over target claims established or refuted;
- **role accuracy**: the same, on molecular role and on activity, reported apart and never pooled across
  the two;
- **coverage**: claims that some observation the table allows establishes, refutes or cannot read, over
  claims stated with one definite value, on all five axes.

Unresolved and observation-model-inadequate counts are reported beside each accuracy, never folded into
it. `refutable` is reported beside target and context. The reasons a claim was not judged are reported
beside coverage, one per claim, the first that applies in this order: judging observations outside the
claim's scope (another gene, another cell), not judged until the target responds, only observations that
suggest, only observations that cannot judge this axis, no observation. Origin and context accuracy are
reported under their own names beside the three. **A `Share` refuses addition, and no field of the report
combines two quantities or two axes.**

**The split.** The scorer uses C4's discipline and nothing else. Every held-out source that judges a
labelling is first passed to `holdout.check_provenance`, and its `LeakError` is never caught. A labelling
built from `holdout.evidence(without=S)` is judged with `sources=[S]` alone. An observation of a kind
holdout does not load is refused with the same error when the labelling read its source. The pilot (S3,
lane-pilot) is scored by the same call.

### Applied once, to the unchanged compiled labels

**Scope.** The compiled programs as they stand, the same files C4 reads. Origin and molecular role come
from every predicted `element` block and every `region` block. Target, activity and context come from
every `rule` with `evidence: predicted`. The measured twins are the evidence, never the labelling, and
are not read. The labelling reads the model, the cCRE registry, RepeatMasker, the segmental duplication
track, phyloP constraint, GENCODE v50 and the unknown-block classes: none of them is a holdout source.
Every holdout source judges it: ten CRISPRi studies, three lentiMPRA cells, VISTA, saturation mutagenesis
and GTEx.

**Expected, written before the run.** Most claims on every axis are not judged. Origin and molecular role
are judged by no loaded source: 0 claims, by construction. Target is judged on the order of 100 of the
440,377 rules and refuted on about none, because a refutation needs a null in the stated cell; the claims
whose only null is in another cell (about 20) are counted beside. Activity is judged on about the same
claims as target, a few of them refuted by the opposite sign. Context is established on about 40 and
refuted on about none.

**What it may be called:** an internal development benchmark reading: of N claims a labelling states on
axis A, J can be judged by an observation the S4 table lets establish or refute A (coverage J / N); of
the J, C are established and I refuted (accuracy C / (C + I)), with U unresolved and M
observation-model-inadequate beside. **What it may not be called:**

- a validation, fresh or external: every source here has been read by this project before;
- the accuracy of the labelling on claims no observation judges: the judged claims are the ones screens
  chose to test, near the genes they chose, mostly in K562;
- a single score, or a comparison between axes: each quantity has its own denominator and evidence;
- evidence that an unjudged element has no function.

**Tests** (`tests/test_attribution_correctness.py`, synthetic units only):

- the three quantities are never summed: a `Share` raises on addition, `sum` fails, and no key in the
  serialised report names a combined score;
- for every observation kind, axis and value, observations the table marks *cannot* or *suggests* never
  decide a verdict, even several of them in the claim's scope;
- *the observation model is inadequate* is a verdict, and it survives JSON serialisation of a verdict and
  of a whole report;
- a labelling that read a source, or the CRISPRi held-out file, raises `LeakError`, and every source passes
  through `holdout.check_provenance`;
- a response elsewhere with a null in the stated cell is one context error, not also a target error.

## The result: the evidence held judges no origin and no molecular role of the unchanged labels, and 98 of 440,377 targets; all 98 are established, but only the 39 stated in K562 could have been refuted (item 12 S4, 2026-09-29, lane-s4, later)

Run once, as registered in the section above, in a clean worktree at `b87f59e` (the revision stamp is not
dirty; manifest complete). Result: `data/results/attribution_correctness.json`. It took 24 s on one
laptop CPU. Nothing called the model (0 requests), and the per-element response cache was never opened.
All 16 holdout sources judged the labels, and each passed `holdout.check_provenance`. **Every count below
is an internal development benchmark reading.** Each is on evidence this project has already read, and
none is a validation. The negatives come first.

### Negatives

1. **Almost nothing can be judged.** Coverage per axis, with the reason each unjudged claim was not
   judged (one reason per claim, in the registered order):

   | Axis | Judged / claims | Outside the claim's scope | Until the target responds | Suggest only | Cannot judge this axis | No observation |
   | --- | --- | --- | --- | --- | --- | --- |
   | origin | **0 / 470,077** | 0 | 0 | 0 | 24,495 | 445,582 |
   | molecular role | **0 / 659,403** | 0 | 0 | 3,519 | 36,112 | 619,772 |
   | activity | **98 / 440,377** | 1,377 | 23 | 3,194 | 18,747 | 416,938 |
   | target | **98 / 440,377** | 1,400 | 0 | 3,194 | 18,747 | 416,938 |
   | cell context | **39 / 440,377** | 1,436 | 23 | 3,115 | 18,826 | 416,938 |

   Origin and molecular role are judged on 0 claims, as registered. No loaded source can establish either
   of them, and the two sources that could (the sequence annotation, the cCRE registry) are the labels'
   own inputs. Beside these, the labels state 157,393 groups of unchosen role alternatives and 11,061
   `unknown` roles; the 26,806 region blocks state no activity and no target. None of these is a claim,
   and none is counted as absence of function.
2. **Target accuracy is 98 of 98, and this is not evidence that the targets are right.** The rule
   refutes a target only with a well-powered null in the cell the rule states. None of the 23 targets
   with a well-powered null has it in its rule's stated cell; those 23 are counted apart
   (`well_powered_null_only_in_another_context`). A rule that ignored the stated cell would have counted
   them as wrong. **Only 39 target claims were refutable**, meaning their stated cell was screened on the
   named gene. All 39 are established, and all 39 are stated in K562.
3. **Every judged repression claim is wrong in direction.** Activity (role accuracy on the direction
   axis) is **93 of 98**. All 4 judged `represses_target` claims are refuted, because the gene fell in
   each. One `activates_target` claim is refuted by a rise: CCND1, whose deciding pair is in the HCT116
   study of the benchmark's held-out file. The 23 claims whose gene responded nowhere are not judged on
   direction: a null has no direction.
4. **Role accuracy on molecular role has no denominator.** It is 0 established of 0 judged. 3,519 role
   claims have only suggestive evidence: an active reporter, a VISTA positive, a CRISPRi response or a
   functional saturation-mutagenesis base. By the table, none of these can decide a biochemical role.
5. **Cell context is established on 39 and refuted on none**, all in K562 (context accuracy 39 of 39,
   refutable 39). For 59 claims the gene responds only in cells other than the stated one, and the stated
   cell is untested.
6. **No claim was unresolved, and none reached "the observation model is inadequate".** On every axis
   both counts are 0. This is a statement about these labels, not about the evidence. 181 target claims
   sit on an element with a lentiMPRA cell whose tiles split evenly. But the unchanged labels make no
   reporter claim, and by the table a reporter cannot judge the direction claims they do make. Every decided verdict
   keeps "the observation model is inadequate" among its explanations.

### What else it shows

- **The scorer splits the ablation's measured links by axis and reproduces them.** The 128 predicted
  links the ablation found measured on their own gene are here 98 established targets (97 decreases, 1
  increase), 23 targets with a null only in another cell, and 7 touched only by underpowered nulls
  (cannot judge). The ablation's "93 survive, 5 contradicted by the opposite sign" becomes 98 established
  targets, of which 93 are right in direction and 5 are wrong. Its "23 contradicted by a well-powered
  null" becomes 23 targets not refuted, because the null is not in the stated cell. The registered
  expectations held on every axis. They were written from the ablation's counts, so they show
  consistency with the ablation, not new evidence.
- **Where the judging comes from.** 71 of the 98 target verdicts are decided by a Gasperini2019 pair
  and 87 in K562. 14 of the 98 target and activity verdicts, and 1 of the 39 context verdicts, rest on a
  pair of the benchmark's held-out file, which is used to judge and never as evidence.

**What it may be called:** of 440,377 target claims the unchanged labels state, 98 can be judged by an
observation the S4 table lets establish or refute a target (coverage 98 / 440,377). All 98 are
established (accuracy 98 / 98), with 0 unresolved and 0 observation-model-inadequate beside, and 39
refutable. Activity is 93 / 98 at coverage 98 / 440,377. Context is 39 / 39 at coverage 39 / 440,377.
Origin is judged on 0 / 470,077 and molecular role on 0 / 659,403. **What it may not be called:** a
validation; the accuracy of the labels on the 440,279 unjudged targets, which the screens did not choose
to test; a single score; or evidence that any unjudged element has no function.

## Item 13 pilot registered: one debugger over local neighbourhoods, and a synthetic gate that needs two simultaneous corrections before anything biological is scored (item 12 S3, item 13 C1-C3, 2026-09-29, lane-pilot)

This section registers the pilot's model and search, and its first gate, before the gate is run and
before any score of any pilot labelling on any held-out source is read. The biological gate (the
chromosomes, the pass rule against the three baselines, the metric and the compute budget) is
registered in its own section before it runs. The code is `genomeos/attribution/pilot.py`; the gate's
script is `scripts/pilot_synthetic_gate.py`. No model request; the per-element response cache is not
opened.

### R8's reading, corrected beside the original

R8's pretest registered, for a negative, the reading "coupling carries no information about CRISPRi
targets beyond independent per-element scoring ... so a joint search over this score cannot beat
per-block scoring. R8 closes; no search engine is built". That states more than the run showed. The
run rejected **two particular score transformations**, element competition and gene budget (each a
share d / (sum d + TAU) of the per-element deletion score), on the K562 training pairs the cache
covered. It did not show that joint inference cannot help. The original text stays in
`joint_pretest.READINGS`; the corrected reading is `joint_pretest.READING_CORRECTIONS`, and the two
transformations are `joint_pretest.RETIRED`. The pilot uses neither: no share of a score over genes or
over elements is computed anywhere in it.

### The model

- **Blocks** are ENCODE cCREs (v3). A block's label is one target gene or no target, and the set of
  contexts it is active in (K562, HepG2, WTC11, `invivo` for VISTA, `other` for saturation
  mutagenesis, and any other cell an observation names).
- **Observations** are weighted soft constraints over the blocks their interval overlaps, read as an
  OR. A CRISPRi significant decrease of gene g in cell c is explained when some overlapping block
  targets g and is active in c; a well-powered null is contradicted when that holds. A reporter tile,
  a VISTA element or a saturation-mutagenesis locus is explained when some overlapping block is active
  in its context. A GTEx association is explained when some overlapping block targets its gene.
  CRISPRi increases, underpowered nulls, GTEx non-associations and missing effects are not used.
- **Priors**, never from a held-out source: distance to the candidate gene's TSS as a power law,
  -ln(1 + d / 5 kb); the compiled predicted target (the unchanged labels) with a bonus of
  1 + min(strength, 1); and activity from ENCODE H3K27ac replicated peaks (a block's own cell for K562
  and HepG2; for other contexts the share of 13 biosamples with a peak). A part's prior is the
  share-weighted sum of each block's own prior over the interval it covers, so a split conserves a
  block's prior, a merge neither adds nor drops one, and no block lends its compiled target to its
  neighbour.
- **Gene-level coupling** (the pilot only): a gene that the view saw decrease under CRISPRi at
  another element gets +1.0 as a target; a gene with at least three well-powered nulls elsewhere and
  no decrease gets -0.5. The element asked about is excluded, so direct evidence is never counted twice.
- **Hard constraints** are only the mandatory ones: a part lies inside its block, parts do not
  overlap, and no part is shorter than 50 bp.

**Weights, in nats, registered and not fitted.**

| term | weight |
| --- | --- |
| CRISPRi significant decrease left unexplained | 3.0 |
| CRISPRi well-powered null contradicted | 1.5 |
| GTEx association left unexplained | 0.5 |
| reporter tile active / inactive | 1.0 / 0.3 |
| VISTA positive / negative | 1.0 / 0.5 |
| saturation mutagenesis functional / inert | 1.0 / 0.3 |
| activity prior (logit) | -1.0 |
| own-cell H3K27ac peak / no peak (logit) | +2.0 / -1.0 |
| breadth, other contexts (logit) | 2.0 x share of the 13 biosamples |
| the compiled link's model cell is this context (logit) | +0.5 |
| no target (log-weight) | -3.0, that of a gene at about 100 kb |
| split (fragmentation) / merge / an observation group marked inadequate | 2.0 / 1.0 / 2.0 |

### C1, C2 and C3 as built

- **The debugger (C1).** From the starting labels (the compiled target, and the activity the
  annotation implies), the violated observations name the conflicting blocks. Repairs are proposed
  only on those blocks' variables:
  - change the target, to a gene an observation there names, the compiled target, or no target;
  - change an activity context an observation there names;
  - split the block where evidence changes;
  - merge it with an adjacent block within 500 bp;
  - mark one observation kind at that block as having an inadequate observation model. This is
    allowed for reporter, VISTA, saturation mutagenesis and GTEx, and never for CRISPRi, the endpoint
    of regulation, whose lower null weight already carries incomplete knockdown.

  Every single repair is evaluated, and every pair whose second repair lies in the conflict the first
  leaves. Only the touched blocks' energy is recomputed. The best repair is applied while it lowers
  the energy by more than 0.05 nat, for up to six rounds. Each violated observation is reduced to a
  smallest conflicting set of observations and label assumptions (deletion-based, exact by
  enumeration over the involved parts). It is reported beside the conflict among the observations
  themselves, which is where a block boundary shows up as an assumption.
- **Boundaries (C2).** A split is proposed only at an observation's edge or an H3K27ac peak's edge
  inside the block. Splits and merges are counted as corrections and must improve prediction of
  withheld evidence like any other correction.
- **Families (C3).**
  - Alternatives within 1.0 nat of the best survive. Survivors that predict the same value for every
    observation of the neighbourhood form one family.
  - A change of the best alternative is committed only if every survivor shares it. Otherwise the
    pilot abstains.
  - The next measurement is the cheapest candidate that separates the best alternative from its
    strongest surviving rival. The candidates are a reporter tile (cost 1) or a CRISPRi pair (cost 4)
    over a block or a part of one. Ties are broken by expected information over all survivors. When
    none separates those two, the candidate with the most information per unit cost is chosen.
- **What a neighbourhood returns:** the ranked alternatives with their energies, the supporting and
  conflicting observations for each, the families, the committed and abstained changes, and the next
  measurement.

### Gate 1: the synthetic gate

**Four motifs.** Each is built on a synthetic chromosome under the weights above. Instance 0 of each is
hand-built, and 100 more are generated from seed 20260929.

1. **Split and retarget** (C2 with C1). One block's two halves are screened against two genes with
   opposite results. The truth is to split at the evidence edge and retarget one half.
2. **Swap** (C1, across blocks). Two blocks sit under two screens that each span both. The truth
   exchanges their targets together.
3. **Context and hand-off** (C1, context and target). The linked block is not active in K562 after all,
   and its neighbour carries the link instead.
4. **Indistinguishable** (C3). One screen spans two blocks, and nothing else tells which of them
   regulates the gene.

**Validity, by brute force** over every single repair and every pair of repairs in the full move
space: every gene and no target, every context, every breakpoint, merge and excuse. A two-correction
instance is kept only if all four conditions hold:

1. no single repair fits the evidence;
2. no single repair lowers the energy by more than 0.05 nat, so a one-repair-at-a-time search is stuck
   at the starting labels;
3. the planted pair fits;
4. the planted pair's energy is below every other alternative of at most two repairs by more than
   1.0 nat.

A family instance is kept only if both single alternatives fit, predict the same for every
observation, and lie within 1.0 nat of each other and of the best alternative of at most two repairs.

**Checks done before this registration.** These checks were construction only, and the pilot's search
was not run on any gate instance. All four hand-built instances pass the validity check. On 200 draws
from another seed, these are the acceptance rates:

| motif | draws accepted |
| --- | --- |
| split and retarget | 188 of 200 |
| swap | 81 of 200 |
| context and hand-off | 200 of 200 |
| indistinguishable | 160 of 200 |

**The pass rule.**

- **The three two-correction motifs.** This part must hold on each hand-built case and on at least 95%
  of the 100 generated instances of each motif. The pilot's returned labels must equal the planted
  truth, and both planted corrections must be committed.
- **The family motif.** This part must hold on the hand-built case and on at least 95% of the 100
  instances. The two alternatives must both survive in one family, and the pilot must abstain on which
  block supplies the link. Its next measurement's predicted outcome must also differ between them.

**If the gate fails, the pilot stops and nothing biological is scored.**

**Controls (descriptive, not in the pass rule):**
- the same search without pairs (greedy);
- the independent-block variant, in which each block is its own neighbourhood with its own copy of
  every observation over it and there is no coupling;
- a noisy copy of each instance with one extra random observation.

**Falsifier.** The gate is void if the greedy control recovers the planted truth on any kept instance:
the validity check would then not be doing its job.

**Compute.** At most 10 CPU minutes, measured with `time.process_time` in the run's own process.

**What each outcome means.**

- **Pass.** The search finds repairs that need two simultaneous changes, and it keeps indistinguishable
  alternatives apart as a family. That is all it means: it is a solver test on constructed cases, and
  says nothing about biology.
- **Fail.** The pilot is not fit to be scored. It is recorded as a discontinued investigation, and
  gate 2 is not run.
- **Void.** Reported as a failure of the gate's construction, with gate 2 not run.

## The result: the pilot finds every planted pair of simultaneous corrections and keeps every indistinguishable pair as one family, but the independent-block control also solves two of the three motifs (item 12 S3, 2026-09-29, lane-pilot, later)

The gate was run once, as registered in the section above, in a clean worktree at `a9660ff` (the
revision stamp is not dirty). The script is `scripts/pilot_synthetic_gate.py` and the result is
`data/results/pilot_synthetic_gate.json`. It took **6.84 CPU seconds** (8.08 s wall, 31 MB peak),
inside the registered budget of 600. No data was read and no model was called. **Verdict: pass.** In
its registered wording: the search finds repairs that need two simultaneous changes, and it keeps
indistinguishable alternatives apart as a family. That is all it means. It is a solver test on
constructed cases, and it says nothing about biology.

### What the pass does not show, first

1. **The independent-block control, which has no coupling, also recovers the planted truth on two of
   the three motifs.**
   - **Split and retarget: 101 of 101.** Both corrections lie inside one block, where the control may
     pair its own repairs.
   - **Swap: 101 of 101.** Here it is right for a reason the pilot forbids. Each block's private copy of
     the two spanning screens counts them twice, so each block alone sees a single repair that helps.
     The pilot reads a spanning screen once, as an OR over the blocks it covers, and for the pilot no
     single repair helps.
   - **Context and hand-off: 0 of 101.** Only this motif needs the coupling.

   So the gate shows that the search finds pairs. It shows the value of coupling on one motif of three.
   Whether coupling helps on real evidence is gate 2's question: the pilot against the independent-block
   variant.
2. **One extra random observation is not in the pass rule.** With one added (the noisy control), the
   pilot still returns the planted truth on 99, 101 and 84 of the 101 instances of the three motifs.
   Context and hand-off, the motif that needs the coupling, is also the least robust to noise.
3. **The greedy control found the planted truth on 0 of 303 instances.** This is the falsifier, and
   it did not fire: the validity check did its job.
4. **The generator rejected these draws before keeping 100 valid instances per motif:**
   - split and retarget: 6 of 106 draws, because a single repair lowered the energy;
   - swap: 130 of 230, for the same reason;
   - context and hand-off: 0 of 100;
   - indistinguishable: 12 of 112, because the two alternatives were not both within the margin of the
     best.

   The kept instances are a constructed subset. The per-motif instance digests are in the result.

### The gate, motif by motif

| motif | hand-built | generated passed | greedy control | independent-block control | noisy copy | median evaluations |
| --- | --- | --- | --- | --- | --- | --- |
| split and retarget | pass | 100 of 100 | 0 of 101 | 101 of 101 | 99 of 101 | 26 |
| swap | pass | 100 of 100 | 0 of 101 | 101 of 101 | 101 of 101 | 52 |
| context and hand-off | pass | 100 of 100 | 0 of 101 | 0 of 101 | 84 of 101 | 115 |
| indistinguishable (family) | pass | 100 of 100 | — | — | — | 34 |

- **The three two-correction motifs.** The pilot resolved all 303 instances, committing both planted
  corrections with no abstention.
- **The family motif.** The pilot abstained on all 101 instances, and both alternatives survived in one
  family.
  - On the hand-built case, its next measurement is a CRISPRi pair of the first block alone for the
    gene in K562. It costs 4 and carries 0.998 bits over the survivors.
  - That measurement separates the best alternative from its strongest rival.

### What it may be called, and what follows

It may be called a passed solver gate on constructed neighbourhoods. It may not be called evidence that
the pilot, or coupling, improves biological labels. Gate 2 is registered in its own section before it
runs. That gate scores the pilot's revised labels on withheld sources against three baselines: the
unchanged labels, the independent-block variant and distance to TSS.

## Item 13 pilot, gate 2 registered: the revised labels scored on withheld sources against the unchanged labels, the independent-block variant and distance to TSS, on seven validation chromosomes (item 12 S3, item 13, 2026-09-29, lane-pilot)

This section was written before any score of any pilot labelling, any validation of a correction and
any S4 judgement was computed, on any chromosome. The model and the search are exactly those
registered in `a9660ff` and passed by gate 1; no weight or search constant changed after gate 1. The
code is `genomeos/attribution/pilot_bio.py`, the script `scripts/pilot_biological_gate.py`, and the
result will be `data/results/pilot_biological_gate.json`. The result may be called **an internal
development result on withheld sources**. It is never a validation: every source here has been read
by this project before, and the CRISPRi held-out file is a reused benchmark that is never evidence.

### What runs

The harness's scorable sources are held out in turn:

- Gasperini2019, Morris, Schraivogel2020 and Xie (CRISPRi);
- the three lentiMPRA cells;
- VISTA;
- saturation mutagenesis;
- GTEx.

For each held-out source S, the pilot reads only `holdout.evidence(without=S)`. It builds neighbourhoods
of ENCODE cCREs on the validation chromosomes around the blocks S's units overlap. S's units give it
coordinates only, and their outcomes are never read by a labelling. The pilot then solves them, and the
harness scores the result on S with `holdout.evaluate(..., chroms=VALIDATION)`. That call uses a locus
bootstrap with 1,000 resamples, seed 20260928, and the same resamples for every labelling. WTC11_DC_TAP
is not held out: its only scorable endpoint (increase, 20 positives genome-wide) is below the floors on
any chromosome subset.

**Chromosomes are the validation partitions.** The validation chromosomes are chr1, chr6, chr8, chr9,
chr10, chr11 and chr19. They are the smallest set found by unit counts alone that keeps all four CRISPRi
decrease endpoints above the harness's floors (20 positives, 20 negatives, 10 loci):

| endpoint | positives / negatives on the validation chromosomes | loci |
| --- | --- | --- |
| Gasperini2019 decrease | 150 / 2,637 | 201 |
| Morris decrease | 21 / 71 | 33 |
| Schraivogel2020 decrease (chr8 and chr11 hold all its units) | 23 / 1,276 | 12 |
| Xie decrease | 26 / 148 | 22 |

Odd or even chromosomes lose Schraivogel2020 and one of Morris or Xie. Of the 15 endpoints C4 scored,
14 stay scorable. The other 17 chromosomes are the development partition. The pilot was timed and
debugged there, and nothing is scored there.

### The labellings

| labelling | role | what it is |
| --- | --- | --- |
| `pilot` | scored | The joint debugger: neighbourhoods are the blocks joined by an observation over two or more of them or by adjacency within 500 bp, cut into windows of at most 60 with the outside blocks held fixed, so an observation crossing a cut is kept on both sides. Gene-level coupling is on |
| `independent_block` | baseline | The same model and search with each block its own neighbourhood, holding its own copy of every observation over it, and no gene-level coupling: each block revised alone |
| `unchanged` | baseline | `holdout.unchanged_labels()`, the compiled predicted layer |
| `distance` | baseline | `holdout.distance_labels()` |
| `prior_only` | descriptive | The pilot's priors alone (distance, the compiled target, H3K27ac), reading no holdable source |
| `pilot_at_unchanged_coverage` | descriptive | The pilot, abstaining wherever the unchanged labels abstain |

**How a labelling scores a unit.**

- **What it reads of the unit:** only the query, which is the unit's interval, gene, TSS and cell.
- **A pair unit** (CRISPRi, GTEx) gets the block's link probability for that gene. For CRISPRi the gene
  must also be active in the unit's cell. The probability is contrasted with "no target" (and
  inactive), with the rest of the neighbourhood as in each surviving alternative. It is never a share
  over other genes.
- **An element unit** gets the block's activity probability in the unit's context: the lentiMPRA cell,
  `invivo` for VISTA, `other` for saturation mutagenesis.
- **Several alternatives:** the score is mixed over the surviving alternatives with weights exp(-ΔE).
- **Several blocks:** the largest score over the blocks the unit overlaps is taken.
- **No solved block:** the labelling abstains (None), which the harness ranks last and reports as coverage.
- **Direction:** the pilot does not model the sign of a regulation, so its `increase` score is its link
  score.
- **Reads:** the pilot and independent-block labellings declare every source in the view,
  `built_without=S`, and the non-holdable inputs (GENCODE v50, the cCRE registry, the compiled layer and
  so the model, H3K27ac).

### The pass rule, the falsifier, the readings and the stop rule

**Primary endpoints:** the four CRISPRi significant-decrease endpoints (Gasperini2019, Morris,
Schraivogel2020, Xie) on the validation chromosomes. **Secondary, reported beside and not in the rule:**
Gasperini2019 increase, lentiMPRA activity and active in each cell, VISTA, saturation mutagenesis and
GTEx.

**Pass.** The pilot passes if two conditions hold:
1. On at least 2 of the 4 primary endpoints, its paired difference in average precision has a 95%
   locus-bootstrap interval above zero against each of the three baselines, at a coverage of at least
   0.80 of the endpoint's units.
2. On no primary endpoint is its interval against any baseline entirely below zero.

The rule tests three baselines on four endpoints and asks for two endpoints against all three, so a
chance pass is unlikely without any correction for multiplicity.

**Falsifier.** A pass is void in either of two cases:

- The prior-only labelling, which reads no holdable source, passes against the unchanged labels and
  distance on the same endpoints, and the pilot's interval against it includes zero there. The gain is
  then the priors', not the debugger's.
- The pilot restricted to the unchanged labels' coverage does not beat them. The pass against the
  unchanged labels is then reported as coverage.

Any LeakError is fatal and never caught.

**Readings, fixed now.**

| outcome | reading |
| --- | --- |
| pass | An internal development result on withheld sources: joint revision predicts the held-out CRISPRi decreases better than the unchanged labels, the same revision one block at a time, and distance to TSS, on the endpoints named, at the coverage given. Phase C may follow; genome-wide deployment waits for a fresh set |
| beats the unchanged labels and distance, not the independent blocks | A failure of the coherence pilot. Revision helps but coupling adds nothing measurable. The joint debugger stops as a discontinued investigation; the per-block revision is a separate finding, to be registered on its own before any use |
| beats the unchanged labels only | A failure. The pilot stops as a discontinued investigation |
| beats none | A failure. The pilot stops and is recorded as a discontinued investigation, apart from implemented capabilities |
| void | A failure of the kind the falsifier names. The pilot stops |

**Stop rule (item 13, binding).** If the pilot does not improve prediction on withheld evidence at
useful coverage, it stops. It is then recorded as a discontinued investigation, separately from
implemented capabilities. A higher internal coherence score is never success.

### The metric: validated corrections per compute-hour

A **committed correction** is a target, activity, split or merge change of the best alternative that
every surviving alternative shares. Each is validated against the held-out source's own binary endpoint:

- CRISPRi: decrease;
- lentiMPRA: active;
- VISTA: positive;
- saturation mutagenesis: functional;
- GTEx: associated.

The test runs over the units on the validation chromosomes that the correction touches. For each
touched unit, δ = +1 if the correction newly predicts it (its new target's pairs, a context switched on)
and -1 if it stops predicting it (the old target's pairs, a context switched off). A split or merge is
judged by its parts' labels against the block's starting label. With π the source's prevalence on
those chromosomes, the correction is:

- **validated** when sum((y - π) δ) > 0;
- an **error** when the sum is below zero;
- **neutral** at zero;
- **untested** when it touches no unit.

**Compute-hours** are the pilot's own CPU time for building neighbourhoods and searching, over all
held-out sources, measured with `time.process_time`. The harness's scoring time and the loading time
are reported beside. Reported beside the metric:

- errors;
- neutral and untested corrections;
- abstentions: changes not committed, and neighbourhoods by status;
- observation groups marked inadequate.

Counts are given for the independent-block variant as well.

### S4 beside, not in the rule

Two sets of claims are judged by `correctness.judge(..., sources=[S])`:

- the pilot's committed labels on the blocks it solved, as target claims (one per active cell, or
  unstated) and context claims;
- the unchanged labels' claims on the same blocks.

Target accuracy, role accuracy and coverage are reported apart, with coverage beside every accuracy.
The pilot states no molecular role and no direction, so it has no role claims. A revision that improves
accuracy on the few judgeable claims while coverage stays near zero is an internal development result
at tiny coverage, and is reported as such.

### Compute budget

The pilot's own work may take at most **2 CPU-hours** over all held-out sources. Per neighbourhood the
search caps of `pilot.SEARCH` apply, and a capped neighbourhood is counted. The run uses one process and
cached data. It opens no per-element response cache and makes no model request, and peak memory is
measured.

### What was done before this registration, and what was expected

**On the development partition, and nowhere else.** Nothing done there computed a held-out score, a
correction's validation or an S4 judgement.

- **Timing.** The joint pilot was timed on chr2 and chr12 with Gasperini2019, Morris, GTEx, lentiMPRA
  K562 and VISTA held out: at most 28.4 CPU seconds a source, 2 neighbourhoods capped, peak 487 MB.
- **Disputes.** Most disputed neighbourhoods there are one reporter tile against the activity H3K27ac
  implies. Most abstain under the registered margin: a flip within 1 nat survives, and so does a no-op
  merge at exactly 1 nat.
- **Smoke run.** The whole script was run on chr21 and chr22 with every call that reads a held-out
  outcome stubbed. That run fixed one crash: a reused outcome's merged parts were not registered in the
  new neighbourhood.

**Expected, before the run.**

- Distance should stay hard to beat on the CRISPRi decrease endpoints, as C4 found.
- The pilot's priors contain distance and H3K27ac activity, which is an ABC-like combination, so any gain
  over distance should appear in `prior_only` too.
- The independent-block variant should be close to the pilot: committed corrections were few on the
  development chromosomes, and most disputes abstain.

None of this changes a rule above.

## The result: the pilot beats the unchanged labels but neither distance to TSS nor its own blocks revised alone on the CRISPRi decrease endpoints, and on Gasperini2019 the coupling makes it worse; the stop rule fires and the pilot is a discontinued investigation (item 12 S3, item 13, 2026-09-29, lane-pilot, later)

**How it was run.**
- **Once**, as registered in the section above.
- **Where:** a clean worktree at `a89bf5c`; the revision stamp is not dirty.
- **Files:** the script is `scripts/pilot_biological_gate.py` and the result is
  `data/results/pilot_biological_gate.json`.
- **Compute:**
  - the pilot's own work took **0.0507 CPU-hours** (182.6 CPU seconds, inside the registered budget of 2
    CPU-hours);
  - the harness's scoring took 45.6 s and loading 4.8 s, for 243 CPU seconds in all;
  - peak memory was 1.83 GB;
  - no model request was made, and the per-element response cache was never opened.

**What it may be called.** An internal development result on withheld sources, never a validation.

**Registered reading: `beats_unchanged_only`.** In the registered wording, "a failure: the revision
improves on the compiled labels but not on distance to TSS; the pilot stops as a discontinued
investigation". **The stop rule of item 13 fires.** The coherence pilot is recorded as a discontinued
investigation, apart from implemented capabilities.

### Negatives, first

1. **No primary endpoint passes.** The pass needed 2 of 4 endpoints above all three baselines.
2. **The coupling never helps, and on Gasperini2019 it hurts.** The pilot minus the independent-block
   variant, in average precision on the validation chromosomes:

   | endpoint | pilot minus independent-block |
   | --- | --- |
   | Gasperini2019 | **-0.065 [-0.125, -0.012]** |
   | Morris | +0.028 [-0.014, +0.085] |
   | Schraivogel2020 | +0.009 [-0.015, +0.109] |
   | Xie | +0.021 [-0.047, +0.130] |

   The Gasperini2019 interval lies entirely below zero. That alone breaks the rule's second condition,
   and means the same revision one block at a time predicts that screen better than the joint debugger.
3. **Against distance to TSS, the pilot is above zero on one endpoint of four.**

   | endpoint | pilot minus distance |
   | --- | --- |
   | Xie | +0.105 [+0.015, +0.261] |
   | Gasperini2019 | +0.087 [-0.002, +0.169] |
   | Morris | +0.003 [-0.174, +0.173] |
   | Schraivogel2020 | -0.016 [-0.440, +0.171] |

4. **What beats distance is not the debugger but its priors.** The prior-only labelling reads no
   holdable source: distance, the compiled target and H3K27ac only.
   - It is above distance on Gasperini2019 (+0.193 [+0.134, +0.249]) and Xie (+0.086 [+0.022, +0.181]).
   - On Gasperini2019 it is above the pilot itself, by 0.106 [+0.041, +0.168].
   - On the other three primary endpoints, the pilot minus prior-only interval includes zero.

   Two things limit this.
   - The priors include H3K27ac. The CRISPRi benchmark's held-out file, where Morris and Xie live,
     selected its positives on chromatin at the element (noted 2026-09-27 above). On Xie, therefore,
     the prior partly reads the rule that made a pair positive.
   - The prior-only gain is descriptive and not in the rule. It is not a finding of this pilot: it
     would need its own registration before any use.
5. **Validated corrections per compute-hour: 78.85, and errors per compute-hour: 78.85.**
   - Of 261 committed corrections, **4 are validated, 4 are errors, 0 are neutral, and 253 (97%) are
     untested**.
   - The 261 count each held-out run. The three lentiMPRA runs share one view, so they repeat the same
     58.
   - All 8 tested corrections are from the GTEx run:
     - 3 validated and 3 errors are target changes;
     - 1 validated and 1 error are merges.

     The validated target changes read a GTEx "not associated" as evidence, which the harness says is
     weaker than a well-powered null.
   - **No committed correction was validated or refuted by any held-out CRISPRi screen.** Morris and Xie
     had 0 committed corrections, Gasperini2019 22 and Schraivogel2020 2, all untested.
6. **Abstentions dominate.**
   - 1,834 changes were not committed, and 416 observation groups were marked inadequate.
   - Over all ten runs, neighbourhoods were: 46,737 undisputed, 2,117 abstained, 129 resolved,
     154 excused and 72 kept.
   - Most disputes are a single reporter tile against the activity H3K27ac implies. A flip within one
     nat survives, so these abstain.
7. **Merges without evidence, a defect of the model found after the run.**
   - 109 of the 261 committed corrections are merges (19 of Gasperini2019's 22).
   - Of the 20 resolved neighbourhoods returned as examples, 6 were resolved by a merge alone that
     explains no observation.
   - The mechanism: the prior prefers any gene within about 100 kb to "no target" (`no_target` = -3.0),
     so the starting labels are not prior-optimal. A merge is the one repair that can hand a block a
     gene that no observation names.
8. **A counting defect in the metric, found after the run.**
   - `validate` compared a merge's label with its first block's starting label only, so the second
     block's change was never counted. The registered definition (each block against its own starting
     label) was not implemented for merges.
   - The registered run's figures above are kept. A recount with the definition as registered follows in
     its own commit.
   - A second defect cannot matter here: the falsifier's "prior-only above the unchanged labels" arm reads
     a comparison the run did not make, so it is always false. It matters only for a pass, and nothing
     passed.
9. **S4 (beside, not in the rule): coverage is tiny.**
   - **Gasperini2019.** The pilot's target claims on its solved blocks: 36 of 44 judged are established
     (coverage 44 of 2,558). The unchanged labels' claims on the same blocks: 28 of 28 (coverage 28 of
     1,652).

     The pilot states each target in every cell it thinks the block is active in, so its claims become
     refutable in K562, and 8 are refuted there. The unchanged labels state the model's cell.

     Context claims: the pilot 28 of 28 (coverage 28 of 2,149), the unchanged labels 12 of 12 (12 of
     1,652).
   - **Schraivogel2020.** Target claims: the pilot 3 of 4 (coverage 4 of 1,051), the unchanged labels 3
     of 3 (3 of 862).
   - **Xie.** Target claims: 1 of 1 each.
   - **Every other held-out source** (Morris, lentiMPRA, VISTA, saturation mutagenesis, GTEx) judges
     none.
   - **Role accuracy.** The pilot states no role or direction. The unchanged labels' direction on
     Gasperini2019's blocks is 27 of 28.

   This is an internal development result at tiny coverage.
10. **Secondary endpoints.**
    - The pilot, the independent-block variant and the prior-only labelling agree within 0.002 on every
      lentiMPRA endpoint.
    - Against the unchanged labels, the pilot is below on WTC11 activity (-0.027 [-0.045, -0.008]). At
      the unchanged labels' coverage, it is also below on VISTA (-0.041 [-0.075, -0.008]).
    - On GTEx, the coupling lowers the score: the pilot minus the independent-block variant is -0.005
      [-0.008, -0.002], and the pilot minus prior-only is -0.0065 [-0.0097, -0.0030].

### The primary endpoints

Average precision on the validation chromosomes is given with its 95% locus-bootstrap interval; the
coverage is in brackets.

| held-out source (positives / units, loci) | pilot | independent-block | prior-only | unchanged | distance |
| --- | --- | --- | --- | --- | --- |
| Gasperini2019 decrease (150 / 2,787, 201) | 0.554 [0.463, 0.677] (0.97) | 0.619 [0.527, 0.727] (0.97) | 0.660 [0.574, 0.757] (0.97) | 0.390 [0.304, 0.494] (0.57) | 0.467 [0.387, 0.581] (1.0) |
| Morris decrease (21 / 92, 33) | 0.866 [0.695, 0.974] (0.99) | 0.838 [0.648, 0.962] | 0.828 [0.625, 0.960] | 0.435 [0.255, 0.610] (0.87) | 0.863 [0.718, 0.965] |
| Schraivogel2020 decrease (23 / 1,299, 12) | 0.375 [0.098, 0.526] (0.85) | 0.366 [0.078, 0.530] | 0.407 [0.078, 0.574] | 0.162 [0.018, 0.316] (0.45) | 0.390 [0.203, 0.672] |
| Xie decrease (26 / 174, 22) | 0.745 [0.566, 0.924] (0.99) | 0.725 [0.541, 0.873] | 0.726 [0.539, 0.874] | 0.469 [0.312, 0.669] (0.83) | 0.641 [0.412, 0.809] |

**What holds, in the registered wording** (an internal development result on withheld sources):

- The pilot's revised labels predict these screens better than the unchanged compiled labels:
  - Gasperini2019: +0.165 [+0.074, +0.266];
  - Morris: +0.431 [+0.269, +0.587];
  - Xie: +0.277 [+0.087, +0.443];
  - Schraivogel2020: +0.213 [-0.052, +0.394], an interval across zero.
- The gain is not coverage alone. At the unchanged labels' own coverage the pilot is still above them:
  - Gasperini2019: +0.111 [+0.026, +0.209];
  - Morris: +0.272 [+0.119, +0.412];
  - Xie: +0.269 [+0.086, +0.434].
- As item 13's instruction to this lane said, beating the unchanged labels alone is not success.

### What the pilot returned, three neighbourhoods of the run

1. **Schraivogel2020 held out, chr8, EH38E3844209 and EH38E3844210** (compiled target CALB1).
   - **The evidence:** a Gasperini2019 CRISPRi decrease of LINC00534 over the region.
   - **The smallest conflicting set:** that decrease and "EH38E3844210 targets CALB1".
   - **The one family** holds two survivors:
     - merge both blocks under LINC00534 (energy 1.029);
     - retarget EH38E3844210 from CALB1 to LINC00534 (1.158).

     The data support the LINC00534 relationship but cannot say whether one block or both supply it.
   - **The pilot abstains.** Its next measurement is a CRISPRi pair of EH38E3844209 alone against CALB1
     in K562 (cost 4, 0.997 bits), which separates the two survivors.
2. **Gasperini2019 held out, chr1, EH38E2780707.**
   - **The dispute:** one K562 reporter tile is inactive, against the K562 activity its H3K27ac peak
     implies. The smallest conflicting set is that tile and "EH38E2780707 is active in K562".
   - **The survivors:**
     - keep the labels (10.30);
     - inactive in K562 (11.00);
     - also active in WTC11 (11.14).
   - **The pilot abstains.** The next measurement is a K562 reporter tile over the block (cost 1, 0.69
     bits). This is the commonest dispute of the run.
3. **Gasperini2019 held out, chr1, EH38E2840923,** "resolved" by merging it under MEF2D (energy -2.94 to
   -4.29). The merge explains no observation. It is the defect in negative 7, returned as a committed
   correction.

### What follows

The coherence pilot stops, as item 13's stop rule requires. It is recorded as a discontinued
investigation, apart from implemented capabilities. Phase C (C6, C5's build, S8's experiment design)
does not follow from it.

Three things are recorded for whoever reopens it. None of them is a result:

1. The prior-only labelling's lead over distance on Gasperini2019 and Xie is a separate observation. It
   needs its own registration, with the H3K27ac selection of the held-out file's positives in mind.
2. The merge defect and the metric's merge count are defects of this pilot's code.
3. On these sources the debugger commits almost nothing that a withheld screen can test (0 of 24 CRISPRi
   corrections tested).

### Correction, the merge count recounted as registered (2026-09-29, lane-pilot, later)

Negative 8's counting defect is fixed in `6bbf1d6`. `scripts/pilot_metric_recount.py` re-solved every
held-out source exactly as the registered run did, in a clean worktree at that commit (stamp not
dirty). With the run's counter it reproduced the run's counts exactly for all ten sources, and only then
counted each merged block against its own starting label, as registered. The result is
`data/results/pilot_metric_recount.json`. Nothing was scored again, so no score, pass or reading changed.

| count | as run (kept above) | as registered |
| --- | --- | --- |
| committed corrections | 261 | 261 |
| validated | 4 | **6** |
| errors | 4 | **6** |
| neutral | 0 | 0 |
| untested | 253 | **249** |
| validated per CPU-hour (0.0507 h) | 78.85 | **118.34** |
| errors per CPU-hour | 78.85 | **118.34** |

The two added validated corrections and the two added errors are GTEx merges: among the GTEx run's
merges, 3 are now validated, 3 errors and 3 untested. Every other run's counts are unchanged. **Still no
committed correction is validated or refuted by any held-out CRISPRi screen.** The stop rule's firing
does not depend on this count: it fired on the scores.

## Rebuilt through its committed entry point: all 25,003 values of clause2_design_power_calibrated reproduce, and a second checkout rebuilds it from its manifest with 0 differences (item 12 S2 follow-up, 2026-09-29, lane-rule)

cab70d9 wrote `data/results/clause2_design_power_calibrated.json` through a scratch driver outside the
repository. The driver called `collect_calibrated` and `save_result` as `--calibrated` does, and the
manifest's argv named the driver: a scratchpad path, at a1900285 with a dirty tree. So
`scripts/manifest_rebuild.py` could not run it.

**The rerun.** `python scripts/clause2_design_power.py --calibrated` was run in a clean worktree at
21a7d37 with the stores linked (1,958 s, **0 AlphaGenome requests**). Since cab70d9 the script had
changed only by the `--read-calibrated` view, which the calibrated path does not call. Every value was
compared with the committed file except the date, the one wall-clock key (`seconds`) and the manifest:
**25,003 values, 0 differences**. The two manifests agree on every field except `code`.

**What was written.** lane-manifest2's pattern. The committed file keeps every value, its date
(2026-09-28) and its timing. Only `result_manifest.code` is replaced, by the clean rerun's: git 21a7d37,
dirty false, argv `scripts/clause2_design_power.py --calibrated`.

**Rebuilt from the manifest.** `scripts/manifest_rebuild.py --venv fresh` then made a second checkout at
21a7d37, with a fresh environment from the committed uv.lock, offline. It found all 201 inputs with
their recorded sha256, ran the recorded argv (exit 0) and compared field by field:
**0 differences**. The bytes differ only in the date, the timing and the code stamp, which the rebuild
ignores. The finding of cab70d9 is unchanged; what changed is that its record can now be rerun.

## The result: the committed reading rule calls a three-quarters rate the model's failure in 50% to 82% of experiments and half the rate "no element" in 57% to 98%; the revised rule holds its error rates in 75 of 210 cells, decides in none, and no rule on this estimator reaches 0.8 at either pole (item 12 S2 follow-up, 2026-09-29, lane-rule, later)

The rule registered in `1cb7559` was scored once, through `python scripts/clause2_design_power.py
--reading-rule` in a clean worktree at `1ef8990` (2,092 s, **0 AlphaGenome requests**), and written to
`data/results/clause2_reading_rule.json` with a complete manifest (dirty false). The gate passed first:
recomputed through the same redraw, the reference calibration, the eligible populations (531, 347, 284,
187) and the registered lines all equal cab70d9's. The model is cab70d9's: 24 designs, every size up to
each cap, six true ratios, 2,000 experiments per cell, seed 2026092913. Both rules read the same
simulated experiments. Every probability below holds under cab70d9's calibrated assumptions and is
never a guarantee.

### Negatives first

- **No design decides clause 2.** In 0 of 210 (design, size) cells does the revised rule reach 0.8 at
  both poles. At its best it reads model_failed at a true ratio of 1 with probability **0.376** (1
  element, 3 own windows, all 531 blocks). It reads wording_wrong at a true ratio of 0 with
  **0.328** (6 elements, 10 own windows, all 187 blocks).
- **A threshold on the committed estimator does not reach 0.8 either.** The Neyman-Pearson threshold,
  read from the simulated experiments themselves and so slightly optimistic, reaches at most
  **0.556** for model_failed against a three-quarters rate and **0.269** for wording_wrong against a
  tenth. A rule that knew the windows' rate exactly and thresholded the blocks' own rate reaches at most
  0.528 and 0.478. Within the eligible blocks, at the registered error rates, no rule of either kind
  decides the clause.
- **The first registered expectation failed.** It predicted both bounds in every own-window cell from
  50 blocks up; they hold in 63 of 81. With 10 windows per block the committed interval lies wholly
  below the one-tenth line in 6% to 13% of experiments at a true tenth, against a nominal 2.5%. When
  positives are rare, the percentile interval over blocks is anti-conservative in its lower tail. With
  shared controls no cell meets both bounds: model_failed errors reach 0.286 and wording_wrong errors
  0.166.
- **The bounds do not survive the anchor's own uncertainty.** The lines are in points, and points move
  with the anchor. With the lines held at the reference and the truth drawn at the clustered interval's
  ends:
  - at 7.88% the rule reads model_failed at a three-quarters rate in 12% to 21% of experiments;
  - at 19.80% it reads wording_wrong at a tenth in 80% to 99.6%, and at a quarter in 21% to 65%;
  - with no false positives in truth, wording_wrong at a tenth rises to 23% to 48%.

  Of the two checked designs that meet both bounds under the reference, 1 element with 3 windows keeps
  them in 7 of 12 feasible variants and 6 elements with 3 windows in 4. Only chromosome ICC 0 keeps
  both.

### The committed rule (c17eedc), simulated at every design's cap

| True ratio | model_failed, range over 24 designs | wording_wrong | cannot_decide |
| --- | --- | --- | --- |
| 1.0 | 0.907 - 0.977 | 0.024 - 0.093 | 0 |
| 0.75 | **0.503 - 0.822** | 0.178 - 0.496 | 0 |
| 0.5 | 0.021 - 0.431 | **0.569 - 0.980** | 0 |
| 0.25 | 0 - 0.113 | 0.887 - 1.000 | 0 |
| 0.1 | 0 - 0.030 | 0.970 - 1.000 | 0 |
| 0 | 0 - 0.006 | 0.994 - 1.000 | 0 |

Above the floor it never says cannot_decide. At a three-quarters rate it reads the model's failure in at
least half the experiments; this seed reproduces cab70d9's 53% to 82% as 50% to 82%. At half the
windows' rate it reads that the blocks hold **no** regulating element in 57% to 98%. At a true ratio of 1
it reads "no element" in up to 9.3%.

### The revised rule (1cb7559), on the same experiments

model_failed / wording_wrong / cannot_decide, at two designs where the rule may read, every eligible block
tested; the committed rule's model_failed / wording_wrong at the first design beside them:

| True ratio | 1 element, 3 own windows, 531 blocks | 6 elements, 10 own windows, 187 blocks | committed, 1 element, 3 windows |
| --- | --- | --- | --- |
| 1.0 | 0.376 / 0 / 0.624 | 0.296 / 0 / 0.704 | 0.975 / 0.025 |
| 0.75 | 0.0245 / 0 / 0.9755 | 0.026 / 0 / 0.974 | 0.583 / 0.417 |
| 0.5 | 0 / 0 / 1 | 0 / 0 / 1 | 0.042 / 0.958 |
| 0.25 | 0 / 0.004 / 0.996 | 0 / 0.001 / 0.999 | 0 / 1 |
| 0.1 | 0 / 0.0555 / 0.9445 | 0 / 0.052 / 0.948 | 0 / 1 |
| 0 | 0 / 0.276 / 0.724 | 0 / 0.3275 / 0.6725 | 0 / 1 |

**Where it may read:** 75 of 210 cells, every one with own windows (75 of 105), none with shared
controls. Of the 30 own-window cells that fail, 26 have 10 windows per block and fail the wording_wrong
bound. The other four have 20 or 30 blocks: 1 element with 1 window at 20 and 30 (model_failed at a
three-quarters rate, 13% and 7%), and 2 and 6 elements with 3 windows at 20 (wording_wrong at a tenth,
6.5% and 6.3%). Those cells read cannot_decide whatever their interval. In the cells where the rule may
read, the largest error is 1.7% to 6.0% for wording_wrong (always at a tenth) and 1.4% to 5.9% for
model_failed (always at three quarters). That is inside the tolerance of 0.05 plus 2 Monte Carlo
standard errors, not at the nominal 2.5%.

### What the arm can honestly report instead

No rule decides at the feasible sizes. What the arm can honestly report, in this order:

1. its coverage, as c17eedc registered;
2. its estimate and interval in points: the difference between the blocks' and the windows' positive
   shares, which needs no calibration;
3. only under cab70d9's assumptions, and labelled so: the registered ratios the interval excludes, and
   the revised reading where the interval clears a line, with the probabilities above beside it.

It does not report model_failed from an interval that only reaches 0, and it does not report "no
element" from an interval that is only below 0.

Descriptively, and not as registered readings, weaker claims the estimator could carry at 0.8 are:

- "no less than half the windows' rate": 0.986 at best (2 elements, 10 windows, 347 blocks);
- "less than half": 0.999.

"Less than a quarter" stays below 0.8 at 0.768. These are Neyman-Pearson bounds: optimistic, and as
dependent on the calibration as the lines are. Applied to the committed measured arm (3 blocks), both
rules read cannot_decide.

### What this settles and what it does not

Under the calibrated assumptions, the revised rule removes both of the committed rule's misreadings. A
near-null difference is no longer read as the model's failure, and a half rate is no longer read as
"no element". The price is that the rule rarely decides. It also shows that within the eligible blocks
no rule on this estimator decides clause 2 at the registered error rates. The rule makes the arm decide
nothing, and its error control holds only if the anchor is right: the anchor's own clustered interval
breaks it. A rule whose lines scale with the experiment's own window rate would not move with the
anchor. It is not registered here. The committed rule stays the record, the revised one sits beside it
dated 2026-09-29, and clause 2 stays not met.

## The pilot's prior-only lead, step 1: the benchmark's held-out file kept a positive only if its element carried H3K27ac, so only Gasperini2019 and Schraivogel2020 may test a prior that reads H3K27ac (item 13 follow-up, 2026-09-29, lane-prior)

**Why this step exists.** Gate 2 of the coherence pilot left one descriptive lead. Its prior-only
labelling reads no holdable source: the pilot's distance power law, the compiled target and H3K27ac.
On the seven validation chromosomes it scored above distance to TSS:

- Gasperini2019: +0.193 [+0.134, +0.249];
- Xie: +0.086 [+0.022, +0.181].

The pilot recorded a caveat itself: the held-out file, which holds Xie, chose its positives on
chromatin. Before any score is read, this step asks of every CRISPRi study in C4's harness whether its
tested elements or its positives were chosen using chromatin at the element.

**Files.**

- The record is `STUDIES` and `BENCHMARK_EVIDENCE` in `scripts/prior_only_test.py`. Every quote there
  is verbatim, with its URL.
- The local counts come from `uv run python scripts/prior_only_test.py --admissibility`. That command
  counts the benchmark's own annotation and scores no labelling.
- No score of any labelling was computed for this step.

### The rule

**What the prior reads at an element:**

- H3K27ac, in its activity term (ENCODE replicated peaks; for K562, ENCFF532MMV of ENCSR000AKP);
- DNase, but only through the ENCODE cCRE registry, which defines its blocks and so where it
  abstains.

**The three classes:**

- **Inadmissible.** The endpoint's positives were selected or filtered using H3K27ac at the element.
  The prior would then read the rule that made a pair positive, and a gain over distance could not be
  told apart from the selection.
- **Admissible, with the selection stated.** The tested elements, positives and negatives alike, were
  chosen on chromatin before the outcome was measured.
  - A selection that does not depend on the outcome cannot make H3K27ac predict the outcome by
    construction.
  - It does restrict the population, so a result holds for chromatin-selected candidates only.
- **Admissible.** Neither of the above holds.

Only the decrease endpoint is considered. The prior states a link and no direction.

### What the benchmark did (Gschwind et al., Nature 2026)

**Held-out file: its positives were filtered on H3K27ac. Its negatives were not.**

- Supplementary Methods 22.1: "To obtain a clean set of positives likely to represent real enhancer
  interactions, we applied filters based on effect size and chromatin state at the tested element."
- Supplementary Methods 22.1.4: "we removed positives with elements in the CTCF element, H3K27me3
  element, or no H3K27ac categories."
- Methods: "Positives lacking H3K27ac signal (CTCF, H3K27me3 or no H3K27ac categories) were removed to
  focus on enhancer-mediated regulation."
- Negatives were filtered by power only.

**The rule was reproduced here exactly.**

- **The unfiltered file:** the paper's repository also holds the held-out file before this filter
  (`EPCrisprBenchmark.combined_heldout.annotated.tsv.gz` in EngreitzLab/ENCODE-rE2G-Paper, 4,437 rows,
  sha256 `2390273590c5…ddd`). It was read in a scratch copy and not added to the project.
- **One rule turns it into the project's held-out file, 4,378 rows of 4,378, with no exception:** keep
  every negative, and keep a positive only if its `elementChromatinCategory` is `H3K27ac` or
  `High H3K27ac` and its `EffectSize` is at most -0.05.
- **All 59 removed rows are positives.** 25 lie outside the H3K27ac categories, and 34 lie inside them
  with a smaller effect. None was relabelled a negative.
- **Morris lost 1 positive and Xie 2.** All three lie on validation chromosomes.

**Training file: no chromatin filter.**

- Methods: "Datasets were filtered to retain element–gene pairs with distances of 1 kb–1 Mb from TSS,
  excluding elements overlapping GENCODE v.29-annotated promoters or target gene bodies. Negative pairs
  required ≥80% power to detect a 15% decrease in expression."
- A copy "with same filters applied as for the filtered held-out dataset" exists in the paper's
  repository and was used for one comparison figure. The project holds the unfiltered file.
- The local check agrees. The project's training file keeps positives outside the H3K27ac categories
  (Gasperini2019: 3 CTCF, 2 no H3K27ac; Nasser2021: 2 CTCF), which the held-out rule would have removed.

**Both files: the elements are DNase peaks.** In the benchmark's re-analysis they are "the top 150,000
K562 DNase-seq peaks" for Gasperini2019, and the ABC pipeline's K562 DNase candidates for Xie, Klann and
Morris. This touches positives and negatives alike.

### Each study

| source | file | tested elements chosen by | positives | class | decrease endpoint in the harness |
| --- | --- | --- | --- | --- | --- |
| Gasperini2019 | training | DNase peaks intersected with **H3K27ac**, p300, GATA1 and Pol II ChIP (pilot library); DNase peaks ranked by a classifier over 170 ChIP tracks (at scale); 948 DNase peaks chosen for epigenomic diversity | readout only | **admissible, selection stated** | 360 / 5,253 genome-wide; 210 / 2,466 on the development chromosomes |
| Schraivogel2020 | training | DNase hotspots overlapping GenoSTAN active-enhancer states (histone marks, p300, DNase; H3K27ac among them not confirmed by a quote) | readout only | **admissible, selection stated** | 23 / 1,299, all on chr8 and chr11 (validation) |
| Nasser2021 | both (one harness source) | K562: every DHS within 450 kb of 30 genes (Fulco 2019), tiling (Fulco 2016), DHS (Klann 2017); GM12878 and Jurkat: accessible regions, ABC-selected elements | held-out part filtered on H3K27ac | inadmissible | unscorable: 111 positives, 0 well-powered negatives |
| HCT116 | held-out | DNase peaks around five genes | **H3K27ac above the median** in the study's own enhancer call, and the benchmark's filter | inadmissible | unscorable: 34 positives, 0 negatives |
| K562_DC_TAP | held-out | DNase peaks chosen at random | benchmark's filter | inadmissible | 12 positives, below the floor |
| WTC11_DC_TAP | held-out | DNase peaks chosen at random | benchmark's filter | inadmissible | 15 decreases, below the floor |
| Klann | held-out | every K562 DNase peak | benchmark's filter | inadmissible | unscorable: 21 positives, 0 negatives |
| Morris | held-out | fine-mapped variants in cCREs from DNase, **H3K27ac** and ATAC peaks | benchmark's filter | inadmissible | 35 / 167 genome-wide; 14 on the development chromosomes |
| Reilly | held-out | mostly tiling; every DHS at GATA1 | benchmark's filter | inadmissible | 8 positives |
| Xie | held-out | DNase peaks with H3K4me1, half chosen for p300 | benchmark's filter | inadmissible | 42 / 435 genome-wide; 16 on the development chromosomes |

The URLs and quotes for each row are in `STUDIES`. The sources are:

- Gasperini: PMC6690346;
- Schraivogel: PMC7610614;
- Fulco 2019: PMC6886585;
- Nasser: PMC9153265;
- Guckelberger: bioRxiv 10.1101/2024.07.12.603288;
- Ray: bioRxiv 10.1101/2025.09.16.676677;
- Klann 2021: bioRxiv 10.1101/2021.03.08.434470;
- Morris: PMC10518238;
- Reilly: PMC8925018;
- Xie: PMC6904118.

### What follows from it

- **Inadmissible: every study of the held-out file.** Its positives are the benchmark's held-out
  positives, kept only with H3K27ac at the element. **This covers Xie and Morris, and so half of the
  pilot's lead.** Nasser2021 is inadmissible as the harness pools it, and unscorable anyway.
  - Restoring the 59 removed positives from the unfiltered file would not rescue Morris or Xie. All
    three of their removed positives lie on validation chromosomes, where the prior was already
    scored.
  - On the development chromosomes they have 14 and 16 positives, below the harness's floor of 20.
- **Admissible, with the selection stated: Gasperini2019 and Schraivogel2020.** Their positives
  are the screens' readout alone. Their tested candidates were chosen on chromatin.
  - Gasperini2019's pilot library intersected DNase peaks with ENCSR000AKP. That is the same K562
    H3K27ac experiment whose replicated peaks the prior reads.
  - A result on Gasperini2019 therefore holds for chromatin-selected candidates only. It says nothing
    about unselected elements.
- **Where each admissible endpoint can be scored.**
  - Schraivogel2020 has every unit on chr8 and chr11. Both are validation chromosomes, where the pilot
    already scored the prior-only labelling.
  - Gasperini2019 is the only admissible endpoint with units where the prior was never scored: 210
    positives, 2,256 negatives and 329 loci on the 17 development chromosomes.

### What this lane read before registering, and what it could not read

**Read, and so already seen when the test is registered:**

- The benchmark's chromatin category of every positive and every pair, per study, on all chromosomes.
  - This shows H3K27ac's association with the outcome on every study. For Gasperini2019, 355 of 360
    positives lie in an H3K27ac category, against 4,911 of 5,299 pairs.
  - Schraivogel2020's 23 positives all lie in an H3K27ac category, against 912 of 1,306 pairs.
- Unit counts per partition.
- No score of any labelling.

**Not read:**

- the GenoSTAN mark list (S1 Appendix of Zacher et al. 2017);
- the ten small sources behind Fulco 2019's 429 borrowed pairs;
- the script that applied the held-out filter. It was reproduced here, not found.

## The pilot's prior-only lead, step 2 registered: tested once on Gasperini2019 on the 17 chromosomes where no pilot labelling was ever scored, and decomposed on identical pairs into distance, distance plus H3K27ac, and the full prior (item 13 follow-up, 2026-09-29, lane-prior)

Registered before any score of any labelling of this test was computed. The constants are in
`scripts/prior_only_test.py` (`registration()`), and the text here is the same.

**Status.** An internal development result on withheld sources, never a validation.

**Compute.** No model request, and the per-element response cache is not opened.

### The lead, frozen

The registered lead is the pilot's own figure:

- prior-only minus distance to TSS on Gasperini2019;
- on the seven validation chromosomes;
- **+0.1927 [+0.1335, +0.2491]**, 0.6601 against 0.4674 (197c560).

The run first recomputes that figure on those chromosomes. **If it does not come back to the fourth
decimal, the run is void and only the discrepancy is reported.** A mismatch would mean the prior, its
inputs or the harness changed after the lead was read.

Xie's half of the lead (+0.0858) is not tested: step 1 found Xie inadmissible.

### Where it is tested

**Fresh: the pilot's 17 development chromosomes.** No labelling of the pilot, the prior-only one
included, was ever scored there. The pass rule is read here and nowhere else.

**Primary endpoint: Gasperini2019 decrease.** It is the one admissible endpoint that clears the
harness's floors there: 210 positives, 2,256 negatives, 329 loci.

**Seen: the seven validation chromosomes.** They are reported beside and carry no weight. They hold:

- the reproduction check;
- the decomposition of the frozen lead on Gasperini2019;
- Schraivogel2020, whose every unit lies on chr8 and chr11.

### The labellings

**Three prior labellings on identical pairs.** Each is built by `pilot_bio.solve(mode="prior")` over the
ENCODE cCREs the held-out units overlap, and read by `pilot_bio.labels_from`.

- They share their blocks, their units and their abstentions (None outside every cCRE).
- They differ by one term at a time. What differs is which inputs are present; no weight of `pilot.W`
  is changed.
- A pair's score is sigmoid(t + 3) × sigmoid(a), the largest over the blocks the unit overlaps:
  - d runs from the block's midpoint to the benchmark's TSS of the pair's gene;
  - t = −ln(1 + d / 5,000) is the distance term;
  - a is the activity logit in the pair's cell.

The three are:

- **(a) Distance alone, `prior_distance`.** There are no compiled links and no H3K27ac peaks, so a = −1
  for every block and the score is monotone in d.
- **(b) The activity-and-distance baseline, `prior_distance_activity`.** Distance plus H3K27ac:
  a = −1, plus 2 with a K562 ENCODE H3K27ac replicated peak over the block, or −1 without. It has no
  compiled target.
  - It reads H3K27ac only. The pilot's prior reads no DNase except through the cCRE registry that
    defines its blocks.
  - **It is not ABC.** ABC multiplies activity by contact and normalises over every candidate element
    of the gene, with activity the geometric mean of DNase and H3K27ac. The project holds no
    experimental DNase-seq signal, only AlphaGenome's predicted DNase. So a true ABC-form score cannot
    be computed faithfully here, and it is left out.
- **(c) The prior-only labelling exactly as the pilot defined it, `prior_only`.** It is (b) plus the
  compiled target:
  - t gains 1 + min(strength, 1) where the block's compiled link names the pair's gene;
  - a gains 0.5 where that link's model cell is the pair's cell.

**Two baselines.** Both are C4's own:

- **distance to TSS:** minus the distance from the unit's midpoint to the benchmark's TSS;
- **the unchanged compiled labels.**

### The metric

**C4's harness, unchanged.**

- **Endpoint and score:** average precision of the CRISPRi significant decrease against the
  well-powered nulls.
- **Abstentions:** ranked below every stated score, with coverage beside.
- **Paired differences:** each comes from `holdout.compare`, which scores both labellings on the same
  units and the same 1,000 locus resamples (seed 20260928), with its 95% interval.

**The six comparisons:**

| comparison | what it asks |
| --- | --- |
| (c) − distance | the frozen lead: **the pass rule** |
| (b) − (a) | activity, on identical pairs |
| (c) − (b) | the compiled target, on identical pairs |
| (a) − distance | the pilot's machinery: its distance form over cCRE blocks and its abstentions |
| (b) − distance | the activity-and-distance baseline against distance |
| (c) − unchanged | the prior against the compiled labels |

### The pass rule, the falsifier and the readings

**Pass rule.** The lead replicates if, on the fresh partition, (c) − distance has a 95% locus-bootstrap
interval above zero on Gasperini2019 decrease, at a coverage of at least 0.80 of its units.

**Decomposition.** Each part is reported whatever the pass:

- activity carries the gain when (b) − (a) lies above zero;
- the compiled target adds to it when (c) − (b) lies above zero;
- the machinery carries it when (a) − distance lies above zero.

**Falsifier.** A replication is not an activity-and-distance reading:

- (i) if (b) − (a) does not lie above zero;
- (ii) if (c) − (b) lies above zero while (b) − distance does not. The gain is then the compiled
  target's, an AlphaGenome-derived term.
- (iii) if (a) − distance lies above zero while (b) − (a) does not. The gain is then cCRE coverage and
  the distance form.

The run is void if the reproduction check fails. Any LeakError is fatal and never caught.

**Readings.**

- **`replicated_activity`.** An internal development result on withheld sources:
  - on chromosomes where it was never scored, the frozen prior predicts Gasperini2019's CRISPRi
    decreases better than distance to TSS;
  - on identical pairs, H3K27ac activity adds to distance.

  The project's prior reproduces the known activity-and-distance relationship on CRISPRi decreases whose
  positives were not chosen by chromatin, among candidates that were. Not a new finding; not a
  validation of the attribution layer.
- **`replicated_not_activity`.** The lead replicates, but a clause of the falsifier fires. No
  activity-and-distance claim may be made, and what carries the gain is named from the decomposition.
- **`not_replicated`.** The lead stays a descriptive observation of the validation chromosomes and may
  not be used.
- **`reversed`.** The lead was specific to the chromosomes it was read on.
- **`no_admissible_endpoint`.** The lead cannot be tested with what the project holds.
- **`void_reproduction`.** Only the discrepancy is reported.

### What a pass would mean, and what it would not

**What a pass would mean.** The project's prior reproduces the known relationship: an element carrying
H3K27ac near a gene's TSS is a CRISPRi hit for that gene more often than its distance alone predicts.

- Its weights were registered before any score.
- The data are K562 CRISPRi decreases whose positives were not chosen by chromatin.
- The chromosomes are ones where this prior was never scored.
- The tested candidates were chosen on chromatin: DNase peaks, and for Gasperini2019 also H3K27ac, p300,
  GATA1 and Pol II ChIP. **The statement holds for chromatin-selected candidates only.**

**What a pass would not mean:**

- **A new finding.** The relationship is published on these screens (Fulco 2019, Nasser 2021, Gschwind
  2026), and this project read it on 2026-09-16 (next section).
- **A validation** of the attribution layer, of the compiled labels, or of the coherence pilot. The
  pilot stays discontinued.
- **A reproduction of ABC.**
- **An external, independent or fresh validation.**
- **Anything about** unselected elements, about cells other than K562, or about the held-out file's
  studies.

### What had already been seen when this was registered

Correction 2 of the coordinator's review: excluding a selected endpoint handles selection, not
independence. Everything below is development evidence.

1. **`crispri_benchmark` and `crispri_contact`** (138824f and 5a31c39, 2026-09-16) read "activity over
   distance" against distance on the pooled training file (Gasperini2019, Nasser2021 and Schraivogel2020,
   K562). **Average precision was 0.519 against 0.441** (9,237 pairs, 451 positives).
2. **The held-out file** has been read by at least ten scored results (`crispri_split_audit`, 0af7e1b).
3. **C4** (lane-c4, b7e4bf0, 2026-09-28) scored distance, the unchanged labels and rest on every CRISPRi
   endpoint on all chromosomes.
   - Gasperini2019 decrease: distance 0.493 [0.439, 0.560], and the unchanged labels 0.425 at coverage
     0.55.
   - Schraivogel2020: distance 0.390, the unchanged labels 0.162.
4. **The pilot's gate 2** (lane-pilot, 197c560) scored the prior-only labelling, the pilot, the
   independent-block variant, the unchanged labels and distance on the validation chromosomes.
   Prior-only minus distance was:
   - +0.1927 on Gasperini2019;
   - +0.0858 on Xie;
   - −0.0345 on Morris;
   - +0.0165 [−0.4385, +0.1934] on Schraivogel2020.

   Nothing was scored on the development chromosomes.
5. **This lane, before registering** (step 1, d19ccac), read the benchmark's chromatin category of every
   positive and pair per study. Gasperini2019 has 355 of 360 positives in an H3K27ac category, against
   4,911 of 5,299 pairs. It also read unit counts per partition. It computed no score.
6. **The literature.** ABC was built on Fulco 2019's screens, and ENCODE-rE2G was trained on the
   training file.

### S4 beside, and the budget

**S4 is descriptive and in no rule.**

- **What is judged:**
  - (c)'s labels on the blocks it solved on the fresh partition, as target and context claims
    (`pilot_bio.claims_of`);
  - the unchanged labels' claims on the same blocks.
- **How:** each is judged by `correctness.judge` with `sources=[crispri:Gasperini2019]`.
- **What is reported:** target accuracy, context accuracy and coverage, apart.
- (a) and (b) name no target, so they state no claim.

**Budget.** One process and at most 0.5 CPU-hour.

## The result: the prior's lead over distance replicates where it was never scored (+0.163 [+0.096, +0.222] on Gasperini2019), but H3K27ac does not carry it; the compiled target does, so no activity-and-distance claim may be made (item 13 follow-up, 2026-09-29, lane-prior, later)

**How it was run.**

- **Once**, as registered in the section above (ab837b9).
- **Where:** a clean worktree at that commit. The revision stamp is not dirty.
- **Files:** the script is `scripts/prior_only_test.py` and the result is
  `data/results/prior_only_test.json`.
- **Compute:** 8.9 CPU seconds and a peak of 1,248 MB (`ru_maxrss`). **0 AlphaGenome requests**, and the
  per-element response cache was never opened.

**What it may be called.** An internal development result on withheld sources, never a validation.

**Registered reading: `replicated_not_activity`.** In the registered wording: "the lead replicates where
it was never scored, but a clause of the falsifier fires: activity is not what carries it, and no
activity-and-distance claim may be made. What carries it is named from the decomposition."

**The frozen lead reproduced first, to the fourth decimal.** On the validation chromosomes, prior-only
minus distance is **+0.1927 [+0.1335, +0.2491]**, 0.6601 against 0.4674. That is the pilot's 197c560
figure, so the run is not void.

### Negatives, first

Every figure below is Gasperini2019 decrease on the 17 fresh chromosomes unless it says otherwise:

- 210 positives, 2,256 negatives, 329 loci.
- All three prior labellings state a score for 2,347 of 2,466 units (coverage 0.9517) and abstain on
  the same 119.

1. **Activity carries nothing measurable on identical pairs.** Distance plus H3K27ac minus distance
   alone, (b) − (a), is **+0.0098 [−0.0188, +0.0349]**. On the seen chromosomes it is +0.0159
   [−0.0069, +0.0427].
2. **The activity-and-distance baseline does not beat distance to TSS.** (b) − distance is −0.0194
   [−0.0604, +0.0185]. On the seen chromosomes it is +0.0084 [−0.0175, +0.0420].
3. **The pilot's machinery alone is slightly worse than distance to TSS.** (a) − distance is
   **−0.0292 [−0.0580, −0.0062]**, entirely below zero. Two things cost it:
   - it measures from the block's midpoint rather than the unit's;
   - it abstains on 119 units that lie outside every cCRE.

   On the seen chromosomes it is −0.0075 [−0.0206, +0.0070].
4. **Falsifier clauses (i) and (ii) fire, and (iii) does not.**
   - (i): activity does not carry the gain.
   - (ii): the compiled target does, and the activity-and-distance baseline is not above distance.

   So the registered activity-and-distance reading, "reproduces the known activity-and-distance
   relationship", **may not be made.**
5. **Schraivogel2020 says nothing.** It is admissible but lies on the seen chromosomes only: 23
   positives, 1,276 negatives, 12 loci. Every interval includes zero:
   - (c) − distance: +0.0165 [−0.4385, +0.1934];
   - (b) − (a): −0.0309 [−0.4027, +0.1096];
   - (c) − (b): +0.0387 [−0.0784, +0.0796].
6. **S4, beside and in no rule: coverage is about 5%.**

   | claims judged (`correctness.judge`, `sources=[crispri:Gasperini2019]`) | target: established of judged | target coverage | context: established of judged | context coverage |
   | --- | --- | --- | --- | --- |
   | the prior-only labelling (c) on its fresh blocks | 53 of 59 (6 refuted) | 59 of 1,122 | 40 of 40 | 40 of 991 |
   | the unchanged labels on the same blocks | 43 of 43 | 43 of 823 | 24 of 24 | 24 of 823 |

   - The prior restates each compiled target in every cell it calls active. That makes some claims
     refutable in K562, and 6 are refuted there. Gate 2 saw the same.
   - (a) and (b) name no target and state no claim.
7. **Half of the lead was never tested.** Xie is inadmissible (step 1).

### What holds, in the registered wording

The lead replicates on the chromosomes where the prior was never scored, at coverage 0.95 (≥ 0.80).

| Gasperini2019 decrease, fresh chromosomes | average precision [95% locus bootstrap] | coverage |
| --- | --- | --- |
| (c) prior-only, as the pilot defined it | **0.6815 [0.6141, 0.7487]** | 0.9517 |
| (b) distance plus H3K27ac | 0.4990 [0.4221, 0.5913] | 0.9517 |
| (a) distance alone, in the prior's form | 0.4892 [0.4126, 0.5823] | 0.9517 |
| distance to TSS (C4) | 0.5184 [0.4383, 0.6066] | 1.0 |
| unchanged compiled labels (C4) | 0.4562 [0.3911, 0.5226] | 0.5365 |

The paired differences (`holdout.compare`):

| comparison | fresh (primary) | seen (the frozen lead, beside) |
| --- | --- | --- |
| (c) − distance: **the pass rule** | **+0.1631 [+0.0956, +0.2223]** | +0.1927 [+0.1335, +0.2491] |
| (b) − (a): activity | +0.0098 [−0.0188, +0.0349] | +0.0159 [−0.0069, +0.0427] |
| (c) − (b): the compiled target | **+0.1825 [+0.1298, +0.2280]** | +0.1843 [+0.1225, +0.2361] |
| (a) − distance: machinery | −0.0292 [−0.0580, −0.0062] | −0.0075 [−0.0206, +0.0070] |
| (b) − distance | −0.0194 [−0.0604, +0.0185] | +0.0084 [−0.0175, +0.0420] |
| (c) − unchanged | +0.2253 [+0.1585, +0.2909] | +0.2705 [+0.1871, +0.3505] |

**The decomposition names what carries the lead: the compiled-link terms.**

- **The terms:** the compiled target's log-weight bonus, and the +0.5 activity where the link's model
  cell is the pair's cell. Both are AlphaGenome-derived.
- **The data:** the same on both partitions.
- **Why the compiled labels alone score low:** alone they score below distance, at 0.54 coverage.
  The gain comes from putting them together with distance inside the prior.

### What it may and may not be called

**It may be called** an internal development result on withheld sources:

- on chromosomes where it was never scored, the pilot's frozen prior predicts Gasperini2019's CRISPRi
  decreases better than distance to TSS;
- but activity is not what carries it, and no activity-and-distance claim may be made;
- the decomposition names the compiled target, an AlphaGenome-derived term.

**It may not be called any of these:**

- **A reproduction of the activity-and-distance relationship.** The registration's hoped-for reading
  failed.
- **ABC**, or a test of ABC.
- **Evidence that H3K27ac does not matter.** The prior reads binary peak overlap on candidates that
  were chosen partly by H3K27ac (Gasperini's pilot library used ENCSR000AKP, the prior's own K562
  experiment). In this sample 4,911 of 5,299 pairs lie in an H3K27ac category. One possible reason
  activity added nothing is that range restriction, and it was not tested. The project's earlier
  "activity over distance" (0.519 against 0.441, 2026-09-16) used a quantitative signal on the pooled
  training file.
- **A finding about the compiled target or AlphaGenome.** It was not the registered hypothesis, only
  what the decomposition names. The project has scored AlphaGenome deletion against these screens
  before (`crispri_benchmark`, `crispri_published`), and so has the AlphaGenome preprint.
- **A validation of any kind**, or anything about Xie, Morris or cells other than K562.

### Added after the run: the coordinator's audit of how the screens chose their elements

The coordinator's audit of six K562 screens arrived after the registration (ab837b9) and after the
run. It is recorded here, dated and labelled as added later. Nothing in the registration or the result
was changed for it.

**Where it agrees with step 1's record (d19ccac), which read the same sources:**

- **Gasperini 2019 used H3K27ac to select elements in both libraries.**
  - Pilot library: DHS intersected with H3K27ac, p300, GATA1 and Pol II.
  - At-scale library: a classifier over 170 K562 ChIP tracks.
- **Schraivogel 2020 selected through the GenoSTAN active-enhancer state.** That state is presumed to
  include H3K27ac, but this is not verified from S1 Appendix.
- **Fulco 2019 (most of Nasser2021) selected on DNase only.**
- **Fulco 2016 is tiling. Klann 2017 selected on DNase only.**
- **Nasser2021 has a minority slice of 429 pairs** from ten smaller sources whose selection is
  unverified.

**Two additions, stated as limits:**

1. **Xie did not use H3K27ac, but half its set was chosen for strong p300.** p300 binding correlates
   with H3K27ac, so Xie's tested universe is H3K27ac-enriched through p300.
   - Xie stays inadmissible for the reason step 1 gives.
2. **Every tested universe except Fulco 2016's tiling is gated on DNase.** An activity term that reads
   DNase faces that selection on every CRISPRi endpoint, and excluding datasets cannot remove it.
   - The prior's activity term reads H3K27ac, not DNase. But its blocks are ENCODE cCREs, which are
     DNase-anchored.
   - Its coverage and its abstentions therefore read the same gate. On Gasperini2019 it abstains on
     only 119 of 2,466 units, and its machinery term, (a) − distance, is slightly below distance.
   - **This is a limit of the test, not something it handled.**

**Where the audit's summary does not cover this test's endpoints.** "No screen defined its positives
with a chromatin mark" holds for the six K562 screens audited. It does not hold for two things:

- the benchmark's held-out file, whose positives the benchmark itself filtered on H3K27ac (step 1,
  reproduced row for row);
- Guckelberger 2024 (HCT116), whose own enhancer call required "an H3K27ac level above the median of
  all tested elements".

For the held-out file, then, the issue is a label chosen by chromatin, not only a selected universe.
That is why step 1 excluded it.

### What follows

- **The pilot's lead is accounted for.** On Gasperini2019 it is the compiled target combined with
  distance, not H3K27ac, and it holds on chromosomes the pilot never scored.
- **Its Xie half stays untested,** because Xie is inadmissible.
- **Any use of the compiled target's contribution needs its own registration.** That registration would
  make the compiled target the hypothesis, on an endpoint neither the pilot nor this lane has read.
- **The coherence pilot stays a discontinued investigation.**

## The four repressions S4 judged wrong, traced: no sign error; three are judged in a cell where the model itself predicts the fall the screen measured, and ID1 is contradicted in the measured cell (item 12 S4 follow-up, 2026-09-29, lane-repress)

This is an internal development trace of five existing claims. It is not a validation, and it changes no
label and no verdict. The cause classes were registered in `scripts/repression_trace.py` before any trace
code was written. The lane had already read the rows by then, and the registration says so. Class (f)
and the rule that picks the primary class were written after that reading. Result:
`data/results/repression_trace.json`. No model request was made; the per-element response cache was
read, never asked.

S4's registered reading stands as written: **every judged repression claim is wrong in direction.** The
trace adds where each direction came from and why it disagrees.

### Where each direction came from

Each compiled direction comes from one element's cached AlphaGenome deletion answer. The steps are:

- The answer holds the gene's largest drop and largest rise over 371 tracks.
- `enhancer_target.predict_target` keeps the larger of the two.
- That track's name becomes the rule's `when: cell_type`.
- `compile.py` writes `inhibits` for a rise, and S4 reads `inhibits` as `represses_target`.

For all five claims the model moves the named gene both ways on different tracks. Benchmark line numbers
are in the decompressed table, with the header as line 1. Compiled line numbers are in the local,
untracked programs, which the manifest pins by sha256.

| Gene | Compiled rule | Largest rise | Largest drop | Wins by | Deciding row | Model on the measured cell | Classes | Primary |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CD83 | chr6:53657, inhibits, K562 | +0.3049 K562 | -0.2954 OCI-LY7 | 0.0095 | held-out 67, Nasser2021, GM12878, -0.249, adj. p 3.4e-40 | GM12878 -0.0968, same sign | b, f | **b** |
| HEMGN | chr9:144732, inhibits, Brain_Cerebellar_Hemisphere | +0.1380 cerebellar hemisphere | -0.0962 K562 | 0.0418 | training 9323, Gasperini2019, K562, -0.223, 5.9e-11 | K562 -0.0948 and -0.0962, same sign | b, d, f | **b** |
| ID1 | chr20:73387, inhibits, Whole_Blood | +0.8119 whole blood | -0.7408 left lobe of liver | 0.0711 | training 6322, Gasperini2019, K562, -0.158, 1.8e-14 | K562 +0.1433, **opposite sign** | b, c, d, f | **c** |
| BEX4 | chrX:105892, inhibits, hair_follicular_keratinocyte | +0.1397 hair follicular keratinocyte | -0.1229 K562 | 0.0168 | training 10306, Gasperini2019, K562, -0.294, 6e-82 | K562 -0.1208 and -0.1229, same sign | b, d, f | **b** |
| CCND1 (comparison) | chr11:215344, activates, CD8_positive__alpha_beta_memory_T_cell | +0.1086 K562 | -0.1634 CD8 memory T cell | 0.0548 (the drop) | held-out 104, HCT116, +0.129, 1.0e-4 | no HCT116 value cached | b, f; c unknown | **b** |

### The classes

- **(a) Sign or convention error: absent in all five.** Every step carries the sign the step before it
  implies. The model's value is the alternate allele over the reference, and the alternate is the
  deletion. A rise means `represses`, which is written `inhibits` and read `represses_target`. On the
  screen side, a significant EffectSize below 0 is a significant decrease, the benchmark's `Regulated`
  agrees in all four, and `crispri_decrease` refutes `represses_target`. A convention flipped on one
  side would have turned the three same-sign cells in the table into disagreements.
- **(b) Cell mismatch: present in all five.** It is primary for CD83, HEMGN, BEX4 and CCND1. In each, the
  claim names the model's most extreme track, and the screen measured another cell. For CD83, HEMGN and
  BEX4, the model's own value for the measured cell has the screen's sign.
- **(c) Genuine contradiction: present only for ID1, and primary there.** In K562 the cached model value
  (+0.1433) and the screen (-0.158) disagree, so the disagreement survives aligning the cells. The cache
  keeps one K562 value for ID1, and any other K562 track is not seen.
- **(d) Indirect or ambiguous: present for HEMGN, BEX4 and ID1, and never primary.** The named gene itself
  fell in every case, so another gene cannot be what disagrees.
  - HEMGN: ANP32B also fell in the same screen (-0.214).
  - BEX4: TCEAL8 also fell in the same screen (-0.074).
  - ID1: the model's any-gene prediction is MIR3193 (-0.859 in HepG2), and the run's verdict is "another
    gene in the same domain". The element lies 1.3 kb from MIR3193 and 1.8 kb past ID1's annotated end.
- **(e) Judge inconsistency: absent as a defect.** `correctness.verdict_of` reproduces each verdict and its
  deciding observation. The next subsection explains what `refutable` means.
- **(f) Direction split across tracks: present in all five.** The rise beat the drop by 0.0095 to 0.0711
  in the four repressions. In CCND1 the drop beat the rise by 0.0548.

### What `refutable` means, and whether `incorrect` should be reachable when it is false

**Where the field is set:**

- `Verdict.refutable` is commented "target and context: the stated cell was screened on the gene".
- `_target` sets it when the gene has an establishing or refuting observation in the stated cell.
- `_context` sets it on every decided verdict.
- `_direction` never sets it, so every activity verdict carries the default `false`.

On the activity axis, `false` means "not computed". It does not mean "could not have been refuted". The
registered quantities report the field only beside target and context. The per-verdict records and the
per-axis tallies still carry the default on every axis.

**Where `incorrect` can occur with the field false:**

- **Target and context: never.** Neither can give `incorrect` with `refutable` false.
- **Activity: yes, by the registered rule.** The direction is judged "in the stated cell when it
  responded there, otherwise in every cell it responded in".

**What the committed S4 activity verdicts show:**

- All 5 wrong directions were decided only in another cell.
- All 39 decided in the stated cell are correct.
- The other 54 correct ones were also decided in another cell.

**What another rule would give (descriptive only):**

- If direction were refuted only in the stated cell, as a target is, the 5 would move to not judged and
  activity would read 93 of 93.
- If direction were judged only in the stated cell, activity would read 39 of 39, with 59 moved to not
  judged.

This is a tension between the S4 table, whose CRISPRi cells speak of "the cell screened", and the
activity rule. The judge follows its registration exactly. Changing the rule is a registration question
for the coordinator, not a bug fix, and no verdict was changed here.

### The comparison case

**Where it is.** The CCND1 claim that HCT116 refutes is in `attribution_correctness.json`, not in
`pilot_biological_gate.json`. The pilot's result names neither CCND1 nor HCT116, and its sources do not
include `crispri:HCT116`.

**What it shares with the four.** It shares (b) and (f), in mirror image: the drop won on a T-cell track,
and the screen measured a rise in HCT116. Whether the model agrees in HCT116 is unknown, because the cache
keeps no HCT116 value.

### Fixed, and found but not fixed

**Nothing is fixed.** No bug was found on the path from the model's answer to the label, or from the raw
row to the verdict.

**Found but not fixed: `by_cell` keeps one value per cell name.** `enhancer_target.aggregate` keeps one
`by_cell` value per cell name: the last track above the scorer's 0.05 recording threshold. In three cases
the largest drop or rise sits on another track with the same cell name:

- HEMGN: K562 -0.0948 is kept, against -0.0962 on another K562 track.
- BEX4: K562 -0.1208 is kept, against -0.1229.
- CCND1: K562 +0.1036 is kept, against +0.1086.

This changes no class here, and no S4 verdict reads `by_cell`. `crispri_direction` and `closure` do read
it as a cell's value. Which track a cell should keep is a decision, and the cached answers cannot change
without new model requests.

### What remains unknown

- The model's value for CCND1 on an HCT116 track.
- Whether another K562 track for ID1, or another GM12878 track for CD83, carries the other sign.
- How many of the 156,925 predicted repressions (of 440,377 predicted rules) share class (f). The run
  tables keep only the winning direction.
- Why the model's ID1 prediction swings from -0.74 in liver to +0.81 in whole blood at this element.

**What it may be called:** an internal development trace of five existing claims, showing where each
claimed direction came from and which registered class explains its disagreement with the observation
that decided it.

**What it may not be called:**

- a validation of the labels, the model or the judge;
- a rate, because the five claims were chosen for being judged wrong;
- evidence about claims S4 did not judge.

## The direction rule, versioned: S4's v1 result stays as history, and v2 judges a direction only in the cell the claim states (item 12 S4 follow-up, 2026-09-29, lane-judge2)

This section sits beside S4 and the repression trace and edits neither. S4's committed result
(`data/results/attribution_correctness.json`) was judged under the rule S4 registered. That rule is kept
unchanged as **v1**. A second rule, **v2**, is used for every judging call from this registration on.

### Registered before any claim was re-judged (2026-09-29)

The decision is Albert's, 2026-09-29, quoted in full (`correctness.RULE_DECISION`):

> "Judge rule: preserve the historical result, correct future judging. The current rule falls back to
> another cell when the claimed cell lacks measurements. That cannot establish that a cell-specific
> direction is wrong. Registration makes the analysis traceable; it does not make that interpretation
> valid. Introduce a versioned rule: judge direction in the stated cell; otherwise return 'not assessed
> in this context.' Keep cross-cell disagreement as a separate finding. Apply this symmetrically: the
> trace says 54 correct and five incorrect verdicts came from other cells. All 59 should become
> unassessed under the revised rule—not just the five errors. Preserve the original results beside the
> correction. ID1 remains a separate, useful finding: its cached K562 prediction disagrees with the K562
> experiment. That does not directly refute the compiled whole-blood claim."

**The v2 rule, word for word** (`correctness.DIRECTION_RULE["v2"]`):

> v2 (from 2026-09-29): the direction (activates_target, represses_target) is judged only in the cell
> the claim states, and only by the named gene's significant changes there: a change the claimed way
> establishes it, a change the other way refutes it, and both in the stated cell is
> `observation_model_inadequate`. When the claim states no cell, or the named gene responded in some cell
> but has no significant change in the stated cell, the claim is not judged, with the reason
> `not_assessed_in_this_context`, whatever the gene did in any other cell: agreement in another cell no
> longer establishes a direction, just as disagreement in another cell no longer refutes one. Each
> significant change of the named gene in another cell is kept on the verdict as a separate cross-cell
> finding (`cross_cell`), marked `agrees` or `disagrees` with the claimed direction; it never decides the
> verdict. A claim whose gene responded in no cell keeps v1's reason, since no cell decides it under
> either rule. Reporter values and no_effect_measured, and the target, context, origin and
> molecular-role axes, are judged exactly as under v1

**`refutable` under v2** (`correctness.REFUTABLE_RULE["v2"]`):

> v2: on a direction claim, true when the stated cell was screened on the named gene with an outcome
> that can establish or refute the direction (a significant change there; the table lets no null decide
> a direction), false otherwise; target and context keep v1's definition; where the field is not
> computed (origin, molecular role, and the activity values that are not a direction) v2 writes null,
> never false

**Which rule judged a result** (`correctness.RULE_RECORD`):

> a result judged under v2 names its rule: `judge_rule` in its manifest parameters and `rule` in its
> body. A result that names no rule was judged under v1, the only rule before 2026-09-29; the scripts
> that reproduce such results pin RULE_V1 and add nothing to their output, so their committed files
> still reproduce (scripts/s4_correctness_run.py without `--rule v2`, scripts/prior_only_test.py,
> scripts/pilot_biological_gate.py, scripts/repression_trace.py)

**What was seen before registering** (`correctness.SEEN_BEFORE_V2`). The 39 / 59 split was already known:

> seen before this registration, and used to write EXPECTED_V2: the repression trace
> (data/results/repression_trace.json, commits 54c47f1..641909e) had already counted the committed S4
> activity verdicts by where they were decided: 39 in the stated cell, all established; 54 established
> and 5 refuted only in another cell. It had also computed that judging direction only in the stated
> cell would give 39 judged, 39 established, with 59 moved to not judged. The 39 / 59 split and the S4
> counts expected below were therefore known before this rule was registered. Also read before
> registration: the committed S4-beside blocks of prior_only_test.json and pilot_biological_gate.json
> (their activity, target-refutable and context counts), from which the expected moves there are
> derived. The registration makes the re-judge traceable; it does not make the interpretation valid

**The counts expected under v2** (`correctness.EXPECTED_V2`):

> on the committed S4 claims and sources: activity is judged on 39 of 440,377 claims, all 39
> established (0 refuted, 0 unresolved, 0 observation-model-inadequate); the 59 verdicts v1 decided only
> in another cell (54 established, 5 refuted) are not judged, reason not_assessed_in_this_context, each
> with its other-cell observations kept as a cross-cell finding, 54 agreeing and 5 disagreeing with the
> claimed direction; every other activity reason keeps its v1 count (1,377 outside the claim's scope, 23
> until the target responds, 3,194 only suggest, 18,747 only cannot judge, 416,938 no observation);
> activity refutable is 39; no count moves on target (98 of 98), context (39 of 39), origin or molecular
> role (0 judged). ID1's represses_target claim in Whole_Blood (EH38E3426791) is not assessed in this
> context, with Gasperini2019's K562 decrease of ID1 kept as a disagreeing cross-cell finding. How many
> of the 39 judged verdicts also carry a cross-cell finding was not counted before registration. Beside
> S4, derived and not rerun: in the prior-only test's unchanged_on_the_same_blocks (Gasperini2019)
> activity would go from 41 of 43 to 24 of 24, 19 moving (17 established, 2 refuted); in the pilot
> gate's unchanged_on_same_blocks from 27 of 28 to 12 of 12 on Gasperini2019 (16 moving: 15
> established, 1 refuted), from 3 of 3 to 0 of 0 on Schraivogel2020 (3 moving, all established) and 1 of
> 1 unchanged on Xie; the pilot's and the prior's own labellings state no activity claim, so none of
> their counts move

The re-judge is run once, with the same claims and sources as S4, and written as a new result,
`data/results/attribution_correctness_v2.json`. `attribution_correctness.json` stays as it is, beside
it. If any count other than the 59 moves, the lane stops and reports it.

### The re-judge under v2, beside v1 (2026-09-29, later)

**S4's original reading stands as the historical v1 result.** `attribution_correctness.json` is
unchanged. In the same run, v1 was judged again and matched that file key by key, in order, leaving
out only what records the run (`result`, `date`, `seconds`, the manifest). The v2 result is
`attribution_correctness_v2.json`, judged once from commit 172572c in a clean worktree, with a clean
stamp. Its manifest records `judge_rule: v2`. No model request was made.

**Why v2 corrects future judging.** In Albert's words, the v1 rule "falls back to another cell when
the claimed cell lacks measurements. That cannot establish that a cell-specific direction is wrong."
The same holds the other way: agreement in another cell cannot establish that the direction is right.
So v2 treats both alike. All 59 verdicts decided only in another cell become "not assessed in this
context", the 54 correct ones as well as the 5 wrong ones. What the other cells show is kept beside
each claim as a separate finding.

**The counts, side by side** (`side_by_side` in the result):

| Axis | Claims | v1 judged | v1 established / refuted | v2 judged | v2 established / refuted | Not assessed in this context (v2) |
| --- | --- | --- | --- | --- | --- | --- |
| Origin | 470,077 | 0 | 0 / 0 | 0 | 0 / 0 | 0 |
| Molecular role | 659,403 | 0 | 0 / 0 | 0 | 0 / 0 | 0 |
| Activity | 440,377 | 98 | 93 / 5 | 39 | 39 / 0 | 59 |
| Target | 440,377 | 98 | 98 / 0 | 98 | 98 / 0 | 0 |
| Context | 440,377 | 39 | 39 / 0 | 39 | 39 / 0 | 0 |

Unresolved and observation-model-inadequate are 0 on every axis under both rules. Details on activity:

- **Other not-judged reasons:** every one keeps its count under both rules. That is 1,377 outside the
  claim's scope, 23 until the target responds, 3,194 only suggest, 18,747 only cannot judge, and
  416,938 no observation.
- **`refutable`:** v1 never computed it on activity. Under v2 it is 39, the directions decided in the
  stated cell.
- **Counts that move only with the 59:**
  - `decided_by_kind` goes from 97 decreases and 1 increase to 39 decreases.
  - Verdicts decided with a held-out-file pair go from 14 to 1.
- **Check:** every count change is accounted for by the 59 moved claims (`unexplained_count_changes` is
  empty).

**The 59 are the ones registered.**

- Under v1, all 59 were decided only in another cell: 54 established, 5 refuted.
- Under v2, each keeps its other-cell observations as a cross-cell finding: 54 agree with the claimed
  direction and 5 disagree.
- The 5 that disagree are the four repression claims the trace followed, plus the CCND1 comparison.
- 48 of the 59 other-cell observations are in K562.
- None of the 39 directions judged in the stated cell carries a cross-cell finding.
- Under v2, no repression claim is judged at all.

These are the counts registered before the run. The 39 / 59 split had already been seen in the
repression trace, so agreement with the registration is expected, not evidence for the rule. Albert's
reasoning justifies the rule; the registration only makes the re-judge traceable.

**ID1 is a separate finding.** Under v2, ID1's `represses_target` claim in Whole_Blood (EH38E3426791) is
not assessed in this context. Gasperini2019's K562 decrease of ID1 is kept beside it as a cross-cell
finding that disagrees. Apart from the verdict, the repression trace found that the model's cached K562
prediction disagrees with the K562 experiment ("model on K562: by_cell[K562] +0.1433; measured
EffectSize -0.1582"). That does not directly refute the compiled whole-blood claim.

**Other judging paths, audited and not changed.** Albert's decision covers direction only.

- **`_target`: a finding only, for the owner to decide separately.** It establishes a target from a
  response in any cell. So 59 of S4's 98 established targets were decided only from another cell. Only
  the other 39 had a response in the stated cell (`refutable`). A target is never refuted from another
  cell. If target were judged only in the stated cell, those 59 would move. Nothing on the target axis
  was changed here.
- **`_context`:** decides all 39 of its verdicts in the stated cell.
- **`_reporter`:** falls back to other cells when the stated cell has no reporter.
- **`_no_effect`:** pools every cell when the claim states none.
- S4 states no reporter or `no_effect_measured` claim, so 0 verdicts would move on either path.
- **Origin and molecular role:** no cell is read.

**Beside S4, derived from the committed counts and not rerun.** The prior-only test's and the pilot
gate's results were not regenerated. Their S4-beside blocks would move only on the unchanged labels'
activity counts:

- **Prior-only test, Gasperini2019:** from 41 of 43 to 24 of 24 (19 move).
- **Pilot gate, Gasperini2019:** from 27 of 28 to 12 of 12 (16 move).
- **Pilot gate, Schraivogel2020:** from 3 of 3 to none judged (3 move).
- **Pilot gate, Xie:** 1 of 1, unchanged.

The pilot's and the prior's own labellings state no activity claim. Their scripts pin v1, so their
committed files still reproduce.

**What it may be called:** an internal development benchmark reading of the same S4 claims under the
v2 direction rule, beside the v1 reading.

**What it may not be called:**

- a validation, fresh or external;
- evidence that the 59 directions are right or wrong: they are not assessed in the cell they name;
- a correction of S4's committed result, which stands as the v1 result.

## The target rule, versioned: v3 judges a target claim only in the cell it states, for supported and refuted verdicts alike, and keeps every other cell beside it (item 12 S4 follow-up, 2026-09-29, lane-judge3)

This section sits beside S4, the repression trace and the versioned direction rule, and edits none of
them. S4's result (`attribution_correctness.json`, v1) and the v2 re-judge
(`attribution_correctness_v2.json`) stay as they are. A third rule, **v3**, is v2 with the target claim
judged only in its stated cell. It becomes the rule a call uses when it names none once the code that
implements it is committed; v1 and v2 stay selectable.

### Registered before any claim was re-judged under v3 (2026-09-29)

The decision is Albert's, 2026-09-29, quoted in full (`correctness.RULE_DECISION_V3`):

> "Target axis: judge the compiled target claim only in its stated cell. Preserve the 59 other-cell
> relationships as supporting evidence elsewhere. Report 39 supported in context, 59 unassessed in
> context, 98 supported somewhere. Version the change and preserve earlier results. Apply the context
> restriction to both positive and negative verdicts." And: "Avoid presenting either 39/39 figure as
> general accuracy: it describes the small subset assessable in the stated context."

**The v3 target rule, word for word** (`correctness.TARGET_RULE["v3"]`):

> v3 (from 2026-09-29): a target claim is judged only in the cell the claim states, and only by the named
> gene's observations there that the table lets establish or refute a target: a significant change of the
> gene on silencing the element in the stated cell supports the claim in context (`correct`), a
> well-powered null of the gene in the stated cell refutes it in context (`incorrect`), and both in the
> stated cell is `unresolved`. What the gene did in any other cell neither establishes the target, nor
> refutes it, nor is a condition for either: a response in another cell no longer establishes a target,
> and no longer keeps a null in the stated cell from refuting one. When the claim states no cell, or the
> stated cell holds no such observation of the gene but another cell holds a response, the claim is not
> judged, with the reason `not_assessed_in_this_context`. On every target verdict the stated cell decides,
> and on every one not assessed in this context, each response and each well-powered null of the gene in
> another cell is kept beside the verdict as a cross-cell finding (`cross_cell`): a response `agrees`
> (supporting evidence elsewhere) and a null `disagrees` (contrary evidence elsewhere); it never decides
> the verdict. A claim the stated cell does not decide and whose gene responded in no cell keeps v1's
> reason and detail, with no cross-cell finding (a null only in another cell stays
> `well_powered_null_only_in_another_context`), since no cell decides it under any rule. `refutable` keeps
> v1's definition on target. Direction claims are judged by DIRECTION_RULE[v2] unchanged, and context,
> reporter values, no_effect_measured, origin and molecular role exactly as under v1 and v2

**The v1 condition on a refutation, and what v3 does with it** (`correctness.TARGET_V1_CONDITION`):

> The v1 target rule reads another cell on the negative side too: a well-powered null in the stated cell
> refutes the target only when the gene responds in no cell, and when the gene responds in another cell
> the target is established and the null is read on context instead. S4's registration gives the reason:
> 'A response in another cell with a null in the stated cell leaves the target established and refutes the
> context', with, on context, 'so one failure is never counted on both axes' (docs/ATTRIBUTION.md, S4
> registered, How each axis is judged; S4's test
> test_response_elsewhere_and_null_in_the_stated_cell_is_one_context_error_not_two). The reason is
> bookkeeping between the target and context axes. It is not a condition on what a null can show: S4's
> only such condition is the null's power (PowerAtEffectSize20 >= 0.8, R2), a property of the stated-cell
> observation itself, which v3 keeps; S4 names no requirement that the perturbation be shown to work, in
> the stated cell or elsewhere. The bookkeeping reason cannot be met within the stated cell alone: judged
> only there, the target claim asks what the context claim asks, of the same observation, and telling a
> target failure from a context failure needs another cell. v3 therefore does not keep it on the target
> axis, as the owner's decision requires for negative verdicts: a well-powered null in the stated cell
> refutes the target in context whatever other cells show. Where the gene responded in another cell and
> has a null and no response in the stated cell, that null now refutes the target in context and, under
> the unchanged context rule, the context claim; the two stay on their own axes and are never combined
> (NO_SUM). On S4's claims the case occurs 0 times (target refutable 39, all 39 established in the stated
> cell; context refuted 0). The context axis's own condition, judged only where the gene responds in some
> cell, is not changed by v3 and is left to the owner as a finding

The brief for this lane described the v1 condition the other way round: a refutation allowed only when
the gene responded somewhere, with `not_judged_until_the_target_responds` otherwise. That is the context
axis's condition, and the direction rule's, not the target's. `_target` never returns that reason.

**The owner's terms** (`correctness.IN_CONTEXT_SCOPE`, `correctness.IN_CONTEXT_TERMS`). They are computed
on:

> target claims, and direction claims (activates_target, represses_target) on the activity axis: the two
> claim kinds v3 judges only in the stated cell. Counts of claims, per axis, never combined across axes
> and never a rate

| Term | What it counts |
| --- | --- |
| `supported_in_context` | claims the stated cell supports: verdict correct, decided only by observations in the stated cell |
| `refuted_in_context` | claims the stated cell refutes: verdict incorrect, decided only by observations in the stated cell |
| `unresolved_in_context` | claims the stated cell both supports and refutes |
| `observation_model_inadequate_in_context` | direction claims whose gene moved both ways in the stated cell |
| `unassessed_in_context` | claims not judged, with the reason not_assessed_in_this_context: the stated cell holds nothing that decides them and another cell holds a response of the gene |
| `supported_elsewhere_only` | claims not supported in context that carry at least one cross-cell finding that agrees: a response of the gene in another cell (on direction, the claimed way) |
| `supported_somewhere` | supported in context plus supported elsewhere only: the claims with supporting evidence in some cell |

**Beside every such count, and beside every 39 of 39** (`correctness.IN_CONTEXT_CAUTION`):

> these counts describe only the small subset of claims assessable in the stated context, where a screen
> tested the named gene in the cell the claim names (mostly K562); an N of N here, such as 39 of 39, is
> not general accuracy, of the labelling or of any axis

**Which rule judged a result** (`correctness.RULE_RECORD_V3`):

> a result judged under v3 names its rule (`judge_rule: v3` in its manifest parameters, `rule: v3` in its
> body) and carries `in_context`, the target and direction counts in the owner's terms, with
> IN_CONTEXT_CAUTION beside them and beside each established-over-decided share S4's quantities print on
> target, activity and context. From the code that implements v3 on, v3 is the rule a call uses when it
> names none; v1 and v2 stay selectable and unchanged, and v2's rule record keeps the rules and the
> default it registered. scripts/s4_correctness_run.py without --rule pins v1 and with --rule v2 pins v2;
> --rule v3 judges the same claims and sources under all three rules, refuses to write unless v1
> reproduces attribution_correctness.json key by key and v2 reproduces the judged body and side-by-side
> counts of attribution_correctness_v2.json, and writes attribution_correctness_v3.json beside both, which
> it reads and never writes

**What was seen before registering** (`correctness.SEEN_BEFORE_V3`). The 39 / 59 split was already known:

> seen before this registration, and used to write EXPECTED_V3: the v2 re-judge's audit of the other
> judging paths (docs/ATTRIBUTION.md, 'The re-judge under v2, beside v1', commit 8f3689e) had counted S4's
> 98 established targets by where they were decided: 39 in the stated cell, the 39 refutable, and 59 only
> from another cell; the owner's decision quotes the 39 / 59 / 98 it expects. The split and the counts
> expected below were therefore known before this rule was registered. Also read before registration: the
> target, context and activity tallies of the committed attribution_correctness.json and
> attribution_correctness_v2.json, including v2's cross-cell findings on activity (54 agreeing, 5
> disagreeing), and the committed S4-beside blocks of prior_only_test.json and pilot_biological_gate.json
> (their target, target-refutable, context and activity counts, and the judged verdicts they list), from
> which the moves expected there are derived. No claim had been judged under v3. The registration makes
> the re-judge traceable; it does not make the interpretation valid

**The counts expected under v3** (`correctness.EXPECTED_V3`):

> on the committed S4 claims and sources, target: 39 supported in context, 0 refuted in context, 0
> unresolved, 0 observation-model-inadequate; 59 unassessed in context (reason
> not_assessed_in_this_context, detail 'the gene responded only in other cells'), each keeping its
> other-cell responses as cross-cell findings that agree; 98 supported somewhere (39 in context, 59
> elsewhere only). Every other target count keeps its v1 and v2 value: 1,400 outside the claim's scope
> (1,377 element screened but gene not tested, 23 well-powered null only in another context), 0 until the
> target responds, 3,194 only suggest, 18,747 only cannot judge, 416,938 no observation; refutable 39.
> Moving with the 59: decided_by_kind from 97 decreases and 1 increase to 39 decreases, and verdicts
> decided with a held-out-file pair from 14 to 1. How many of the 39 or the 59 also carry a disagreeing
> cross-cell finding (a well-powered null in another cell) was not counted before registration. Activity
> (direction) is exactly as under v2: 39 supported in context, 0 refuted in context, 59 unassessed in
> context (54 with an agreeing and 5 with a disagreeing cross-cell finding), 93 supported somewhere (39 in
> context, 54 elsewhere only). Context stays 39 established, 0 refuted; origin and molecular role 0
> judged. Between v2 and v3 exactly the 59 target verdicts move; between v1 and v3, those 59 and v2's 59
> directions. The 39 of 39 on target and on direction describe the small subset assessable in the stated
> context and are not general accuracy. Beside S4, derived from the committed counts and not rerun
> (neither result is regenerated, and neither gate's registered verdict reads an S4 block): in the
> prior-only test on Gasperini2019, unchanged_on_the_same_blocks target from 43 of 43 to 24 supported and
> 0 refuted in context, 19 unassessed in context, 43 supported somewhere, and the prior's own labelling
> from 53 of 59 to 40 supported and 6 refuted in context, 13 unassessed in context, 53 supported
> somewhere; in the pilot gate on Gasperini2019, unchanged_on_same_blocks from 28 of 28 to 12 supported
> and 0 refuted in context, 16 unassessed in context, 28 supported somewhere, and the pilot's own
> labelling from 36 of 44 to 28 supported and 8 refuted in context, 8 unassessed in context, 36 supported
> somewhere; on Schraivogel2020, unchanged_on_same_blocks from 3 of 3 to none judged, 3 unassessed in
> context, 3 supported somewhere, and the pilot's 3 of 4 unchanged; on Xie, 1 of 1 unchanged for both; no
> target is judged on any other source. Their activity counts move as EXPECTED_V2 derived; the pilot's and
> the prior's own labellings state no activity claim

The re-judge is run once, with the same claims and sources as S4, and written as a new result,
`data/results/attribution_correctness_v3.json`. The v1 and v2 files stay byte-identical beside it. If
any count other than the 59 target verdicts moves between v2 and v3, the lane stops and reports it
before explaining it.

### The re-judge under v3, beside v1 and v2 (2026-09-29, later)

**S4's reading stands as the v1 result, and the v2 re-judge as the v2 result.**
`attribution_correctness.json` and `attribution_correctness_v2.json` are unchanged, byte for byte. In the
same run, v1 was judged again and matched its file key by key, leaving out only what records the run.
v2 was judged again and matched its file's judged body and side-by-side counts. The v3 result is
`attribution_correctness_v3.json`, judged once from commit b7cb814 in a clean worktree, with a clean
stamp. Its manifest records `judge_rule: v3`. No model request was made.

**The owner's decision** (quoted in full in the registration above): judge the compiled target claim
only in its stated cell, for positive and negative verdicts alike, and keep the 59 other-cell
relationships as supporting evidence elsewhere.

**The counts in the owner's terms, side by side** (`in_context_side_by_side` in the result). Under v1
and v2 a verdict decided with another cell is counted apart.

| Claims | Rule | Supported in context | Refuted in context | Decided with another cell | Unassessed in context | Supported somewhere |
| --- | --- | --- | --- | --- | --- | --- |
| Target | v1 | 39 | 0 | 59, all established | 0 | 98 |
| Target | v2 | 39 | 0 | 59, all established | 0 | 98 |
| Target | v3 | 39 | 0 | 0 | 59 | 98 |
| Direction | v1 | 39 | 0 | 59: 54 established, 5 refuted | 0 | 93 |
| Direction | v2 | 39 | 0 | 0 | 59 | 93 |
| Direction | v3 | 39 | 0 | 0 | 59 | 93 |
| Context | v1, v2, v3 | 39 | 0 | 0 | not a term of this axis | not a term of this axis |

**Beside every figure in this table, and every 39 of 39 in this section** (`correctness.IN_CONTEXT_CAUTION`):

> these counts describe only the small subset of claims assessable in the stated context, where a screen
> tested the named gene in the cell the claim names (mostly K562); an N of N here, such as 39 of 39, is
> not general accuracy, of the labelling or of any axis

Unresolved and observation-model-inadequate are 0 on every axis under all three rules. Origin and
molecular role judge nothing under any rule.

**Target, in detail.**

- **Other not-judged reasons:** every one keeps its v1 and v2 count. That is 1,400 outside the claim's
  scope (1,377 element screened but gene not tested, 23 well-powered null only in another context), 0
  until the target responds, 3,194 only suggest, 18,747 only cannot judge, and 416,938 no observation.
- **`refutable`:** 39 under all three rules.
- **Counts that move only with the 59:**
  - `decided_by_kind` goes from 97 decreases and 1 increase to 39 decreases.
  - Verdicts decided with a held-out-file pair go from 14 to 1.
- **Check:** the moved claims account for every count change, between v2 and v3 and between v1 and v3
  (`unexplained_count_changes` is empty in both).

**The 59 are the ones registered.**

- Between v2 and v3, exactly 59 claims move. All are targets that v2 established only from another cell.
  All become not assessed in this context, with the detail "the gene responded only in other cells".
- Each keeps its other-cell responses beside it as cross-cell findings that agree: supporting evidence
  elsewhere.
- The registration left one count open. It is 0: none of the 59, and none of the 39, carries a
  disagreeing cross-cell finding (a well-powered null in another cell). The 39 carry no cross-cell
  finding at all.
- Between v1 and v3, 118 claims move: these 59 targets and v2's 59 directions.

Every count matches the registration. The 39 / 59 split had already been seen, so agreement is expected,
not evidence for the rule. Albert's decision justifies the rule; the registration only makes the
re-judge traceable.

**The direction correction is resolved.** v3 carries v2's direction rule unchanged. Every direction
count under v3 equals v2's.

**Left to the owner as findings, and not changed here.**

- **One null on both axes.** A gene can respond in another cell and have a well-powered null, but no
  response, in the stated cell. That null now refutes the target in context and, under the unchanged
  context rule, the context claim. The two stay on their own axes and are never combined. This happens 0
  times on S4's claims.
- **Context's own cross-cell condition.** The context axis is judged only where the gene responds in some
  cell. So another cell's response is still the precondition that lets a stated-cell null refute a
  context claim.

**Beside S4, derived from the committed counts and not rerun.** The prior-only test and the pilot gate
were not regenerated, and their scripts pin v1. Neither gate's registered verdict reads an S4 block. On
the target axis their S4-beside blocks would move as follows:

- **Prior-only test, Gasperini2019:**
  - unchanged labels: from 43 of 43 to 24 supported and 0 refuted in context, 19 unassessed;
  - the prior's own labelling: from 53 of 59 to 40 supported and 6 refuted in context, 13 unassessed.
- **Pilot gate, Gasperini2019:**
  - unchanged labels: from 28 of 28 to 12 supported and 0 refuted in context, 16 unassessed;
  - the pilot's own labelling: from 36 of 44 to 28 supported and 8 refuted in context, 8 unassessed.
- **Pilot gate, Schraivogel2020:** unchanged labels from 3 of 3 to none judged, 3 unassessed; the pilot's
  3 of 4 does not move.
- **Pilot gate, Xie:** 1 of 1 for both, unchanged. No target is judged on any other source.

Their activity counts move as the v2 section derived. The caution above applies to every one of these
figures: each N of N describes the small subset assessable in the stated context.

**What it may be called:** an internal development benchmark reading of the same S4 claims under v3, in
the owner's terms, beside the v1 and v2 readings.

**What it may not be called:**

- general accuracy: every N of N here, the 39 of 39 on target and on direction among them, describes the
  small subset assessable in the stated context;
- a validation, fresh or external;
- evidence that the 59 targets are right or wrong in the cell they name: they are not assessed there;
- a correction of S4's result or of the v2 re-judge, which stand as the v1 and v2 results.

## Audit A, candidate coverage: 365 of the 387 decrease links C4 counted without a registry element have one sharing their bases, 362 of them failing only the reciprocal rule's width requirement, and 21 are absent from both element sets (item 13 C4 follow-up, 2026-09-29, lane-place)

This section sits beside C4 and edits nothing in it. C4 (`b7e4bf0`) sorted the 661 measured CRISPRi
decrease links with `ablation.placement`: 212 meet the overlap rule in a compiled predicted element,
62 more would meet it in an ENCODE cCRE, and 387 meet it in neither. The 449 not placed are placement
failures, not 449 proven absent enhancers. This audit asks why each link fails, one cause at a time,
over every screened pair. Decreases and nulls are reported apart at every step, and the ranking and
the intervals stay frozen. `scripts/placement_audit.py`, `data/results/placement_audit.json`,
`tests/test_placement_audit.py`.

**Status.** This is internal development evidence and a bounded audit. It is not a biological result.
Every outcome it reads was read before: `138824f`, `ecc3452`, `f48b909`, `a39073d`, `b7e4bf0`,
`cc824f4` and `b7cb814` (`prior_exposure` in the result). No resplit restores independence.

### Registered before any per-outcome count (2026-09-29, `dd8c49c`)

- **The cascade.** Every link, keyed as C4 keys it (chrom, start, end, gene, cell, split), gets the
  first cause that holds:
  1. placed in a compiled element;
  2. a registry element meets the rule (candidate exclusion, split by whether the whole-chromosome
     sweep scored that element and named a target);
  3. a zero-width source row;
  4. a registry element holds at least half of the smaller interval, but the tested interval is
     outside 75 to 700 bp, so no 150 to 350 bp registry element can ever meet the reciprocal rule;
  5. the same, with a reachable width;
  6. a partial boundary overlap;
  7. no registry element shares a base.
- **The one policy**, fixed from the conventions alone: an element carries a tested interval when
  their overlap is at least half of the smaller of the two widths. The threshold is the measured
  layer's own 0.5. Nothing is padded, moved or merged.
- **Three stops.** The policy is reported as denominator or overlap inflation if any of these holds:
  - (a) its rescue of positives minus its rescue of nulls has a gene-clustered bootstrap 95% interval
    whose lower bound is at or below zero;
  - (b) that excess does not beat the same excess on the same links with every tested interval moved
    20 kb either way;
  - (c) the excess vanishes within the 75 to 700 bp intervals.
- The only outcome figures seen before the registration were C4's 212, 62 and 387.

### Each input's build and convention

The audit reproduces C4's placement exactly before it reads anything else.

| Input | sha256 | Declared | Found |
|---|---|---|---|
| CRISPRi training file | `9eddfe18…d79405` | GRCh38, BED 0-based half-open | see below |
| CRISPRi held-out file | `7783534a…5ffd38` | the same | see below |
| cCRE registry, 24 files | `95e79f14…20c4a6` (combined) | ENCODE cCREs v3 GRCh38, BED 0-based half-open | 1,063,878 elements, 150 to 350 bp |
| Compiled programs, 24 files | `20d2806e…ead0facf` (combined) | the registry's coordinates by id | all 440,377 predicted elements are registry ids at the registry's exact coordinates (0 moved, 0 missing) |

What the CRISPRi files show, read without outcomes:

- **Build.** The GRCh38 columns touch a registry element at 91.2% of the 5,638 distinct tested
  intervals. The same intervals moved 20 kb touch one at 41.6% and 41.9%. The name field's source
  coordinates, where they differ (3,749 lifted-over intervals), touch one at 37.2%, the chance level of
  another build.
- **TSS.** The benchmark's TSS sits within 1 kb of GENCODE v50's for 2,147 of 2,334 genes, with a median
  difference of 29.5 bp.
- **Liftover** kept the width of all but 66 rows.
- **One base.** Moving every tested start or end by one base changes the cause of 1 to 3 positives and
  23 to 25 nulls. The DC_TAP rows are 499 bp wide, a one-base end convention of the source.
- **One source defect.** The held-out Nasser2021 PTGER4 row is chr5:40490507-40490507, with zero width
  in its source name too.

### The attrition (current rule, training / held-out)

Positives are the 661 decrease links. Nulls are the 14,073 other links: 9,801 well-powered nulls,
4,113 underpowered nulls and 159 significant increases. The result keeps each null kind apart.

| Cause | Positives | Nulls |
|---|---|---|
| placed in a compiled element | 212 (169 / 43) | 4,613 (3,002 / 1,611) |
| registry would: candidate exclusion | 62 (43 / 19) | 3,182 (2,236 / 946) |
| of which: registry element not scored by the sweep | 11 | 467 |
| of which: scored, no target named | 22 | 1,565 |
| of which: scored, only a non-coding target named | 29 | 1,150 |
| zero-width source row | 1 (0 / 1) | 0 |
| rule and width: no registry element can reach the reciprocal rule at this width | 290 (183 / 107) | 3,174 (2,513 / 661) |
| rule and width: reachable width, these widths fail it | 72 (60 / 12) | 1,892 (1,240 / 652) |
| partial boundary overlap | 3 (1 / 2) | 111 (90 / 21) |
| no registry element shares a base | 21 (15 / 6) | 1,101 (804 / 297) |

**Reading.**

- **No convention bug.** There is no build, liftover or 0/1-based error. All three sets are GRCh38 BED
  half-open, and the compiled set is the registry's.
- **C4's label is inaccurate.** `not_placed_no_registry_element` is wrong for 365 of its 387 links: a
  registry element shares bases with each of them. 362 fail only because the reciprocal rule asks two
  intervals to be within a factor of two in width. Under that rule, a tested interval wider than 700 bp
  cannot be placed anywhere, whatever DNA it holds, and 290 of the 661 positives are that wide.
- **Placement is width, not DNA.** The current rule places positives and nulls at the same rate:
  32.1% against 32.8%, a difference of −0.7 points [−4.8, +3.4].
- **Genuine absence is small.** 21 positives have no registry element sharing a base, so none in
  either set. One more is the zero-width source row.
- **The 62 are the model's element selection or the sweep's scope.** No compiled element was dropped
  or moved.

### The registered policy, scored (`0183ef5`, clean stamp)

No registered stop fires.

| | Current rule | Policy |
|---|---|---|
| positives placed in a compiled element | 212 / 661 (32.1%) | 497 / 661 (75.2%) |
| nulls placed in a compiled element | 4,613 / 14,073 (32.8%) | 7,706 / 14,073 (54.8%) |
| placed or registry would, positives / nulls | 41.5% / 55.4% | 96.2% / 91.4% |
| placed and carried by more than one compiled element, positives / nulls | 0 / 0 | 131 (26.4%) / 1,593 (20.7%) |
| chance placement: intervals moved +20 kb and −20 kb, positives | 3.2%, 4.4% | 21.3%, 21.5% |
| chance placement: the same, nulls | 5.7%, 5.7% | 18.8%, 18.3% |

- **Rescue.** Among the links the current rule leaves unplaced, the policy rescues 285 of 449
  positives (63.5%) and 3,093 of 9,460 nulls (32.7%).
- **Stop (a).** The excess is +0.308 [+0.258, +0.356].
- **Stop (b).** The excess on moved intervals is +0.049 and +0.046, and the excess minus the moved
  excess is +0.26 [+0.214, +0.305].
- **Stop (c).** Within 75 to 700 bp the policy rescues 46 of 150 positives against 1,031 of 6,172 nulls,
  an excess of +0.140 [+0.068, +0.215].

**Errors and where the gain sits.** The following table was added after the first run and is outside
the registration.

- 235 of the 285 rescued positives, and 2,061 of the 3,093 rescued nulls, are tested intervals wider
  than 700 bp.
- In that class about half of the policy's placements are ambiguous: 117 of 235 positives and 1,081 of
  2,061 nulls. At reachable widths it is 14 of 258 and 512 of 5,644.
- About 28% of the policy's positive placements (21.3 / 75.2) would also occur 20 kb away, against
  about 12% under the current rule.

**Matched baseline: distance to TSS on the same links.**

| Links | AUROC | AP | Prevalence |
|---|---|---|---|
| all 14,734 | 0.877 | 0.403 | 4.5% |
| placed under the current rule (4,825) | 0.880 | 0.455 | 4.4% |
| placed under the policy (8,203) | 0.876 | 0.471 | 6.1% |
| rescued by the policy (3,378) | 0.873 | 0.508 | 8.4% |

The policy's placed set is no easier or harder for distance by AUROC. Its AP rises with prevalence, so
any score on it must be read against 0.471 on the same links, not against the current rule's figure.

**Independent loci.**

- Positives: 661 links on 415 genes and 594 tested intervals.
- Nulls: 14,073 links on 2,105 genes and 5,318 intervals.
- The intervals are gene-clustered bootstraps, 2,000 draws.
- Fold hashes: training `a679e35a…`, held-out `0cc6255d…` (the result's `folds`).

**Cost, all in.**

- The run: 9.5 s wall, 9.4 s CPU and 1.14 GB peak memory. Most of the memory is the whole-chromosome
  element tables, read for membership only.
- A first run, discarded because another session's untracked script made its stamp dirty, cost the
  same.
- Preparation was a handful of outcome-blind reads of widths and name fields, each a few seconds, plus
  three test runs of about 10 s.
- 0 downloads, 0 model requests and 0 paid services.

**What it may be called.** An internal development audit of C4's placement. It finds that the
reciprocal rule's width requirement, not absence, is why 362 of the 387 links fail, and that one
outcome-blind policy passes its three registered stops on development evidence.

**What it may not be called:**

- 449 absent enhancers, or 387 links with no registry element;
- a validation, or independent evidence: every outcome here was read before;
- evidence that a policy placement names the right element: a quarter of them are ambiguous, and about
  half are in the wide class;
- a change to C4. The policy is applied nowhere. C4's survival census uses the same reciprocal rule
  (`CrispriIndex.of`), so it carries the same width blind spot; that is not measured here. Adopting the
  policy anywhere would be its own registered change.

*2026-09-29, added afterwards beside this section (the external reviewer's reading, checked by the
coordinator): "null" here means every screened pair that is not a significant decrease, so the 14,073
include 4,113 underpowered observations and 159 increases beside the 9,801 well-powered nulls. Read it
as a non-decrease comparator, not as evidence of no regulation, and do not read the equal placement rates
under the current rule as showing that placement carries no signal about real links. The registered
result above is unchanged; a sensitivity restricted to well-powered nulls would be additional
development analysis.*

## Audit B, cached context contrast, checkpoint 1: no-go, insufficient coverage for this registered complete-case design (the external reviewer's audit B, 2026-09-29, lane-contrast)

**The question.** Can the per-element cached AlphaGenome deletion answer support one added feature: the
same-cell signed deletion minus the median other-cell signed deletion, per element, with a frozen
aggregation? This checkpoint asks only whether it can. It scores no labelling, makes no model request,
reads 4DN contacts from their local cache only, and downloads nothing. The rules, including the go rule,
were committed in `7870427` before any count (`scripts/context_contrast_feasibility.py`). The result is
`data/results/context_contrast_feasibility.json`, stamped `d237824`, clean.

**What the cache holds.** The sweep's scorer was built with `threshold=0.0`
(`scripts/enhancer_targets_all.py`, `worker_scorer`), not the adapter's default of 0.05. So every
gene-track value except an exact 0.0 is emitted, and values below 0.05 are stored. Of the 4,840 answers
overlapping an admissible pair, 4,810 store a `by_cell` value with |v| <= 0.05, as do 6,895 of the 6,928
pair-to-answer matches (0.9952). The statement beside the repression trace that `by_cell` keeps "the last
track above the scorer's 0.05 recording threshold" is corrected here for the sweep's answers: the
threshold was 0.0. Its finding stands: `by_cell[cell]` is the value of the **last** emitted track carrying
that cell name, signed and rounded to four places. It is not a replicate aggregate and not the largest,
and the track's identity is not stored. The track table, every other track of a cell name, and the
per-track values are not on disk. A track that returned exactly 0.0 is dropped, so its gene's `n_tracks`
falls below 371. Two genes of one name share a row, so `n_tracks` rises above 371. Only four cells are
kept: K562, the same cell for both admissible studies, and HepG2, GM12878 and IMR-90, so k is at most 3.
A missing cell is never read as zero, and the other cells are a context, not negative controls.
`n_tracks = 371` is a consistency check: every value was emitted, so no exact-zero track can have moved a
cell's value to an earlier track of that name. It does not prove independently which track `by_cell`
holds, or that one model version answered every element. All 4,840 answers are unrequested versions.

**Completeness on lane-prior's admissible endpoints (decrease endpoint).** Independent loci are counted
on the parent universe (each set's own units), frozen before any filtering. The recount after filtering
splits connected components (329 became 351 on the primary endpoint), so it is descriptive only. At every
k from 1 to 3 the counts are identical, because a reliable row holds all four cells.

| Endpoint set | Units (parent loci) | An answer overlaps | Gene listed, same-cell value | Reliable contrast, k = 1, 2, 3 | Complete with activity and distance (coverage) | Independent loci of the complete pairs (filtered recount) | With K562 contact |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Gasperini2019, fresh (primary) | 2,466 (329) | 2,319 | 1,376 | 1,368 | 1,368 (55.47%) | 282 (351) | 1,366 |
| Gasperini2019, seen | 2,787 (201) | 2,660 | 1,556 | 1,547 | 1,547 (55.51%) | 173 (217) | 1,547 |
| Schraivogel2020, seen | 1,299 (12) | 1,088 | 706 | 699 | 699 (53.81%) | 11 (18) | 699 |

All three sets together hold 6,552 pairs in 529 independent loci. The 3,614 complete pairs fall in 454 of
them, and 8 loci hold both studies. The two studies share no Reference and so are not provenance
siblings. Every unit is K562 and comes from the benchmark's training file. The activity columns (DHS.RPM,
H3K27ac.RPM) and the TSS are present for every pair. Cached K562 contact covers 99.64% of primary pairs,
with 0 cache misses and no matrix opened.

**Where the coverage is lost.** 1,368 of the 1,376 primary pairs that hold a same-cell value also hold a
reliable contrast, so most of the loss predates the contrast. It lies in the deletion answer itself:
- 894 pairs whose gene lies beyond the scorer's reach from every overlapping answer;
- 147 pairs that no cached answer overlaps;
- 49 pairs whose gene is in reach but not listed.
Only 8 are lost to the reliability conditions: 6 rows with exact-zero tracks and 2 merged rows.

**Duplication.** The registered bars were not reached, so the feature is not a duplicate by the rule.
- **D1.** The compiled target's description explains 0.7565 of the contrast's variance, against the bar
  of 0.80. That variance sits in the tail the compiled target names. It names 215 of the 3,614 pairs, and
  81 of the 82 K562-activating links have a negative contrast.
- **D2.** The Spearman correlation with the same-cell value is 0.6356, against 0.95.
- **Resolution.** 8.99% of contrasts lie inside the rounding bound of 1e-4. The median |contrast| is
  0.0011.

**Reading: `no_go_completeness`.** G1 (semantics) and G3 (not a duplicate) pass. G2 fails: the complete
pairs clear the floors (190 positives, 1,178 negatives) but cover 55.47% of the primary endpoint's units,
against the registered 80%. This is insufficient coverage for this registered complete-case design. It
does not say that a context contrast cannot help. As registered, no score registration is written, and
the lane stops here.

**A possibility only, not registered and not launched.** A baseline with the feature optional, using the
baseline's own score wherever the contrast is missing, could keep the baseline's coverage. Nothing here
tests it.

**Prior exposure. Every dataset read is development evidence.**
- `crispri_benchmark`, `crispri_contact` (`138824f`, `bda8a83`, `5a31c39`)
- `crispri_published` (`42d7b1b`, `a39073d`)
- `target_calibration` (`bdc2364`)
- `crispri_direction`, `crispri_direction_both` (`5c842ca`, `d492834`)
- `c4_holdout_scores`, `c4_alphagenome_ablation` (`b7e4bf0`)
- `pilot_biological_gate` (`197c560`)
- `prior_only_test` (`a48e9e8`)
- `repression_trace` (`641909e`)
- S4 v1, v2 and v3 (`ab62d99`, `4567219`, `d1e09ae`)

Before writing its rules, this lane probed chr21's archive and the 3,209 loose answer files for their
threshold and track counts, with no pair joined.

**Hashes.**
- **Folds** (sha256 of the sorted pair identities, no outcome):
  - primary `48f6895a…`, its complete pairs `5051f01d…`;
  - Gasperini2019 seen `d9515f04…`;
  - Schraivogel2020 seen `d642ec7d…`.
- **Sources:** training file `9eddfe18…`, the 3,209 loose answers `d07d978b…`, K562 contact cache
  `8b48075b…`. The 23 archives and 24 compiled programs are pinned one by one in the manifest.

**Cost.** Two runs of one process each:
- the provisional run on `7870427`, stamped `0183ef5`;
- the final run on `d237824`, which added only the parent-universe loci, the loss steps and these
  notes. Its gate, semantics and duplication numbers are identical to the provisional run's.

Each run took about 50 s wall and 50 s CPU. The streaming reader peaked at 1.30 GB resident at most,
against about 5.4 GB for a whole-archive reader; 23 archives were streamed and 3,209 loose files read. There were 0
downloads, 0 model requests and 0 network requests.

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
