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

## The matched background is not a neutral baseline, and that is why the controls fail (2026-09-16)

Seven of the 23 chromosomes the panel has read fail one of its six controls, and on all six autosomes
and X the failing check is the same one: the neutral tier's fixed rate missing the GC- and
replication-timing-matched background it is measured against. Six special chromosomes would be one
explanation. `scripts/panel_tier_baseline.py` asks the committed results whether it is instead one
thing happening everywhere, and it is.

| tier against the matched background | above 1 on | one-sided p | median ratio |
|---|---|---|---|
| regulatory | **22 of 22** | 2e-7 | 1.078 |
| fossil | **21 of 22** | 5e-6 | 1.058 |
| neutral | **20 of 22** | 6e-5 | 1.107 |
| constrained unknown | 6 of 22 | 0.99 | 0.917 |
| coding exons | 0 of 22 | 1.0 | 0.387 |

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

- **constrained unknown 0.827, below the neutral tier on 18 of 22** — the only unknown tier that is
  depleted, and it stays depleted under either baseline (below the background on 16 of 22, p 0.026).
- **fossil 0.974 (17 of 22) and regulatory 0.983 (12 of 22)** — indistinguishable from neutral, which
  against the background both appeared to exceed.
- **coding exons 0.356, below it on 22 of 22** — the positive control, unchanged by the choice.

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
panel, with chr2 still unread.

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
| the named coding gene is inside the element's own CTCF node | 77.2% | 68.4% (chr13) to 91.2% (chr19), median 75.1% | 74.3% for as many boundaries placed at random, measured on these same elements |

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

**What the deletion adds is which elements, not which gene.** Given one call per element, the
deletion calls only where its top target was tested and the predicted drop clears a threshold. At
0.1 it makes 124 calls on 3,472 elements and 119 are right (96%); at 0.2, 84 of 84. On those same
elements the closest tested gene is right 117 and 83 times. Choosing the gene is mostly distance on
this set (activity is one value per element, so within an element it picks the closest gene too). The deletion's contribution is high precision about which elements act at all, at low
recall: at 0.1 it reaches 30% of the 400 elements with a regulated gene. The nearest TSS inside the
CTCF node names a tested gene on 488 elements and is right on 223 (46%). The node boundary costs
recall here without buying precision.

**Limits.** One gene per element was kept by the sweep, so the deletion cannot rank a second
target; a full per-gene table would cost the sweep again. The benchmark's activity columns are
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
