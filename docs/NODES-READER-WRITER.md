# Nodes, reader, writer, executor: your model against the biology

You proposed that DNA, though linear, works in nodes; that a cell reads only
the nodes it needs (a nail cell and a neuron carry the same text and run
different parts); and that a software design for a human would have
libraries for the cell, a progression from one cell to billions, the organ
systems, an isolated reader, and an isolated writer. Here is how each idea
maps onto what is known, what GenomeOS already has for it, and what it
changes in the design.

## 1. The genome works in nodes: true, and measurable

The linear text folds. Chromosomes are organised into **topologically
associating domains** (TADs, typically 0.1–1 Mb), inside which a gene and
its enhancers touch each other while being insulated from the neighbours;
the boundaries are marked by **CTCF** sites held by cohesin loops, and
domains group into active (A) and inactive (B) compartments. A node, in
your words, is a domain: the genes plus the regulatory elements that reach
them, delimited by CTCF boundaries. Its content is read as a unit and can be
switched off as a unit (heterochromatin, B compartment).

What GenomeOS has: the block map at the gene level; the pattern file now
carries the CTCF core motif so boundary candidates can be found in
sequence; and, streamed from ENCODE, the experimental registry of candidate
regulatory elements (promoters, enhancers, CTCF sites) for chromosome 21.
What follows: a `Domain` block above `Gene` in BioIR, built from CTCF sites
and cCREs, so the map shows nodes, not only genes.

## 2. Same text, different cell: the reader is the epigenome

Neuron and nail cell differ not in DNA but in what is readable: DNA
methylation, histone marks and chromatin openness decide which nodes the
transcription machinery can enter, and a small set of master transcription
factors decides which programme runs inside them. That is the reader: the
machinery is shared, the *configuration* is per cell type, and it is
inherited through cell division. In GenomeOS the reader configuration is
`cell_type ... expresses:` (which silences the rest), the epigenetic clocks
(a reader setting that drifts with age), and the promoter `requires:` lists
derived from motifs (which factors must be present to read a node).

"Part of the DNA is ignored for a specific function" is exactly right: in
a typical cell type roughly half the genes are off, and whole domains sit in
the inactive compartment.

## 3. The writer: three different writers, one of them rare

- **Replication** copies the text; it is a copier, not an author.
- **Epigenetic writers** (DNMTs, histone modifiers) write the reader
  configuration during development and keep it through divisions: this is
  the writer that turns one cell into a neuron lineage or a skin lineage,
  and it writes on chromatin, not on the sequence.
- **Sequence writers** are rare and bounded: the germline (recombination,
  mutation across generations), V(D)J recombination in lymphocytes (the one
  programmed somatic rewrite), transposons, and damage. Mature cells do not
  rewrite their own instructions on purpose.

So the design has a copier, an epigenetic writer that sets reader state, and
an evolutionary writer that edits the text between generations. GenomeOS
models the copier (replication and its timers), the epigenetic writer as
events and cell-type context, and the evolutionary writer only as the
variants a genome carries.

## 4. The executor

Transcription, translation, folding, metabolism, signalling: the machinery
that turns a read node into proteins and behaviour. This is the runtime
(BioVM) and the core libraries; it is the part the genome must always keep
on, and the part the minimal-cell design counts as essential.

## 5. From one cell to billions: nodes plus timers

Development is the epigenetic writer running under timers: lineage trees
(the C. elegans engine), the segmentation clock, master switches turning
nodes on in sequence, and maintenance clocks afterwards. Your "instruction
when to multiply and maintain" is the timer library.

## 6. Your library list, mapped

| your library | biology | GenomeOS |
|---|---|---|
| cell: structure, components, rules, read/write/update/remove/kill, mitochondria | core libraries + organelle genomes | `core.*`, mtDNA, essential set |
| nodes: progression from one cell to billions | lineage, cell cycle, differentiation, timers | lineage engine, ageing runtime, `timer.*` |
| specific systems | blueprint and systems libraries | `blueprint.organ_*`, `systems.*` |
| isolated reader | transcription machinery + epigenetic configuration per cell type | `cell_type` context, `requires`, clocks |
| isolated writer | copier / epigenetic writer / germline | replication timers, events, variants |

## 7. How to figure out what each block is

Layered evidence, each layer with its own kind and confidence:

1. **Sequence patterns** (done): repeats, ORFs, islands, motifs, learned signals.
2. **Experimental chromatin data** (this step): ENCODE cCREs say where the
   reader actually binds; CTCF sites give the node boundaries.
3. **Conservation**: a block kept across species is under selection; our
   second genome (C. elegans) is too far for this, mouse is the next.
4. **Predictive models** (AlphaGenome, `predicted`): what a block does to
   expression in each tissue when it or its neighbours change.
5. **Perturbation data** (Virtual Cell models, CRISPR screens): what happens
   when the node is switched.

Each layer moves a block from UNKNOWN toward a named function with a
confidence, and the interface should show which layer said so.

## Do enhancers stay inside their node? A predicted test (2026-09-11)

The node model makes one testable claim: an enhancer reaches genes inside
its CTCF domain and not across the boundary. GenomeOS could not measure this
(no Hi-C), so it asked a predictive model. For 200 distal enhancers of
chromosome 21, spread along its length, AlphaGenome deleted the element in
its 1 Mb window and reported the predicted expression change of every gene
across 371 RNA-seq tracks (`scripts/enhancer_targets.py`, feature b in
docs/ALPHAGENOME.md, `data/results/enhancer_targets_chr21.json`):

- 90 elements move a coding gene by ≥ 0.1 log2; for 78 of them (86.7%) the
  gene that moves most is inside the node the element sits in, and for 61
  (67.8%) it is exactly the nearest TSS in the node, the target `genomeos
  regulation` had inferred at 0.4.
- 12 elements name a coding gene beyond the boundary. Some are the model's
  known long-range weakness, some may be boundaries the CTCF-only inference
  draws in the wrong place; each one is a candidate for a better node, which
  is what the `regulation` rows now show as `predicted` next to `inferred`.
- 110 elements move no coding gene by that threshold, and the share that
  names a gene falls from 52% under 20 kb to 32% beyond 100 kb: the model
  reads the near ones and is agnostic about the far ones, as it says of
  itself.
- 43 of the 127 named effects are rises, not drops: the registry's
  "enhancer-like" class, defined by chromatin marks, includes elements that
  behave as silencers for the gene they move.

This is a prediction agreeing with an inference, and both are labelled as
such; it is not a measurement. But it is the first evidence in the project,
beyond the placement of CTCF sites, that the node is the right unit for
regulation, and it gives every scored enhancer a named gene, a tissue and a
magnitude instead of a distance.

The same job then ran on every chromosome (4,800 elements, 200 per
chromosome, `data/results/enhancer_targets_genome_wide.json`): 2,291 move a
coding gene, and for 90.2% of them the gene is inside the element's node
(79.8% to 91.8% per chromosome, median 86.4%); 71.2% are exactly the
nearest TSS the inference had named; 1,157 of the 2,994 named effects are
rises, silencer-like. Twenty-four chromosomes giving the same shape is what
turns one chromosome's observation into a property of the genome as this
model reads it: the node bounds regulation, and nearest-gene is right two
times in three.

## The second mammal: mouse chr19 through the same code (2026-09-11)

If nodes are a property of mammalian genomes rather than of one annotation,
a mouse chromosome fetched and analysed exactly like a human one should give
nodes whose genes sit together in human too. `genomeos mouse --chrom chr19`
does that with no mouse-specific code beyond the three URLs: UCSC mm10
sequence (19 MB), GENCODE vM25 gene models (rows of chr19 kept), ENCODE
SCREEN's mouse cCRE registry (11,597 elements on chr19), the same
`infer_domains` (`genomeos/genome/mouse.py`,
`data/results/mouse_mm10_chr19.json`). Orthology is gene-symbol identity
(App ↔ APP), a cheap lower bound.

| mouse chr19 (61.4 Mb) | MGI orthology (curated) | symbol identity |
|---|---|---|
| genes | 1,394 (718 coding) | |
| nodes | 371 (256 with coding genes), median 109 kb; human chromosomes give 102–115 kb | |
| mouse nodes with ≥ 2 genes matched to a human node | 109 | 91 |
| all matched genes in one human node | 58 (53%) | 49 (54%) |
| in adjacent human nodes (same neighbourhood, boundaries drawn differently) | 45 (41%) | 39 (43%) |
| scattered over distant human nodes | 6 (6%) | 3 (3%) |
| same human neighbourhood | 94% | 97% |

The orthology is MGI's curated mouse–human homology report, streamed once
and distilled to 24,584 symbol pairs (`data/results/mgi_mouse_human_orthology.tsv.gz`);
symbol identity (App ↔ APP) stays as the fallback and is reported beside it.

A second, larger chromosome gives the same answer. Mouse chr11 (122 Mb,
1,622 coding genes, 27,432 elements, 836 nodes, median 102 kb; human 17 and
parts of 5, 7 and 22): 270 mouse nodes with ≥ 2 matched genes, 159 in one
human node (59%), 91 in adjacent nodes, 20 scattered, 93% in the same human
neighbourhood (`data/results/mouse_mm10_chr11.json`). Two chromosomes, 379
tested nodes, 93–94% same neighbourhood: the number is a property of the
comparison, not of the chromosome picked first.

Read two ways. Synteny is not the surprise: mouse chr19 is human 11q13 and
10q23–26 and everyone knows it. The number that speaks to the model is the
split: where a mouse node's genes fall into several human nodes, they fall
into *adjacent* ones 45 times out of 51. The neighbourhood is conserved; the
exact boundary is where the two CTCF registries (different depths, different
cell types assayed) disagree. Which is to say: the node is a real unit of
organisation in both genomes, and an inferred boundary from CTCF-only
elements is a resolution limit, not the biology. Curated orthology raised
the tested count from 91 to 109 and left the shape where it was, which is
what the symbol-identity run had predicted. The 262 mouse nodes with fewer
than two matches are mostly nodes with no or one coding gene (115 of 371
have none).

## The edges of the nodes, held against a predicted contact map (2026-09-12)

The CTCF-only boundaries were never more than a proxy. AlphaGenome's
predicted contact map gives an independent opinion: insulation minima at
2 kb over 28 cell types. On chr21, 37% of the 227 inferred boundaries sit
within 20 kb of a predicted minimum, against 28% for boundaries placed at
random, and the median distance is 34 kb (docs/ALPHAGENOME.md, use 4). That
is agreement a little above chance, not confirmation. The results that
carry the node model are therefore about node *content*: the enhancer
deletions landing inside the node on 24 chromosomes (90%) and the mouse
nodes landing in one human neighbourhood (93 to 94%). Where the edges sit
is a resolution question the project cannot settle without Hi-C or a
better predictor of insulation, and the node confidence stays at 0.4.
The measurement is one key away: 4DN publishes boundary calls for GM12878,
H1, HFFc6, K562, HCT116, HepG2 and IMR-90 as BED files, `genomeos domains
--hic GM12878` reads and compares them, and only the account key is
missing (docs/DATA.md).

## Does CTCF orientation say which node edges are real? (2026-09-14)

The node edges are CTCF-only registry elements, placed without reading which way
the CTCF motif points. The mechanism says the strand matters. Cohesin extrudes a
loop until it meets CTCF in the right orientation, so loops close between
convergent sites: a forward site upstream and a reverse site downstream (Rao et
al. 2014).

`scripts/ctcf_orientation.py` tests this. It scans every CTCF-only element with
JASPAR MA0139 at the project's 0.85 threshold. MA0139's forward strand reads
GCCACCAGGGGGCGC, the core Rao et al. call forward. That is 35,783 elements
genome-wide, 14,221 of them with a site. Each domain edge takes the strands of
the elements merged into it (`orient_boundaries` in `genome/domains.py`). A
domain is convergent when its upstream edge carries a + site and its downstream
edge a − site. It is divergent when the pair points away from it. The edges are
held against 4DN's measured boundary calls within 20 kb, with two controls:

- as many positions placed uniformly at random, the comparison's own control;
- strands shuffled 500 times among the elements that have a site, which keeps
  positions and which edges have a site and changes only orientation.

(`domains_ctcf_orientation`)

**The positive control: orientation is real at measured boundaries.**
Single-strand sites within 20 kb upstream of a measured boundary are mostly −,
and those downstream mostly +:

| biosample | − share upstream | + share downstream |
|---|---|---|
| H1 | 66% | 66% |
| K562 | 65% | 65% |
| HepG2 | 67% | 66% |
| IMR-90 | 57% | 56% |

A boundary sits between sites pointing away from it, as extrusion predicts.
GM12878's calls show no lean (50% and 51%). They are dense enough that a random
position lies within 20 kb of one 54% of the time genome-wide, so every
GM12878 comparison below is near saturation. The 58% on chr21 was mostly that
density (random 42% there).

**Our edges split by orientation: no information.**

| class, genome-wide (19,930 edges) | edges | H1 | K562 | HepG2 | IMR-90 | GM12878 |
|---|---|---|---|---|---|---|
| all | 19,930 | 24.7% | 15.4% | 15.8% | 16.3% | 59.3% |
| convergent pair | 2,360 | 31.0% | 21.7% | 20.0% | 18.7% | 58.6% |
| divergent pair (wrongly oriented) | 1,056 | 33.6% | 19.8% | 19.9% | 19.1% | 58.8% |
| a site, neither | 4,989 | 34.0% | 20.5% | 21.2% | 18.8% | 59.7% |
| no site | 11,525 | 18.6% | 11.4% | 12.1% | 14.5% | 59.4% |
| random positions | | 17.9% | 11.3% | 12.1% | 14.1% | 54.2% |

Convergent edges sit at their strand-shuffle medians in every biosample:
31.2%, 21.1%, 20.3%, 18.2% and 59.0%, with p from 0.17 to 0.69. Divergent
edges are no worse than their shuffles. On chr21, as asked, the 33 convergent
edges score 55% and the 14 divergent 86% against GM12878. Both counts are too
small to say anything, and they point the wrong way.

What does carry information is whether an edge has a site at all: edges with a
site are supported one and a half to two times as often as edges without one
in H1, K562 and HepG2. 58% of our edges (11,525) have no MA0139 site at the
project's threshold.

**Orientation used to place boundaries: it works.** Place a boundary between
consecutive single-strand sites (at least 5 kb apart) where the strand flips
from − to +, and take the + to − flips as the control:

| biosample | − to + flips supported | + to − flips supported | random |
|---|---|---|---|
| H1 | 25.2% | 16.2% | 18.0% |
| K562 | 17.5% | 10.2% | 11.1% |
| HepG2 | 19.5% | 9.9% | 12.3% |
| IMR-90 | 18.2% | 12.4% | 13.6% |
| GM12878 | 58.5% | 57.5% | 54.1% |

There are 3,093 and 3,109 of these flips genome-wide. The − to + flips beat
random by 1.4 to 1.6 times in four biosamples. The + to − flips, sites pointing
at each other inside a loop, fall below random, which is what the mechanism
predicts. They also beat our current edges' enrichment over random in K562,
HepG2 and IMR-90 (1.57, 1.59 and 1.33 against 1.36, 1.30 and 1.16), and match
it in H1.

**Reading.** Orientation carries information, but not as a label on the edges
we already have. Those edges are limited by where the caller puts them. It takes
whichever CTCF-only element comes first beyond a 50 kb minimum domain. More than
half of those elements have no CTCF motif. And a single site cannot say which
neighbour it bounds. An orientation-aware caller built on − to + transitions is
the better edge, and it is the next step for the node model. Replacing the
CTCF-only edges in `infer_domains` would change every node on every chromosome
and every enhancer-target result built on them, so it waits for a decision
rather than being slipped in.

## The orientation-aware caller, beside the CTCF-only one (2026-09-14)

The orientation measurement above showed that − to + flips of CTCF motifs mark
boundaries better than the CTCF-only edges do. The next step was a caller built
on that rule. It is `infer_domains(..., orientation=strands)` in
`genome/domains.py`, a named option. Without it the output is byte-identical to
before (checked on chr21, chr2 and chr19), and every committed node still comes
from the CTCF-only caller.

**The rule, fixed before any of the four tests was run.**

- **Sites.** Registry elements with CTCF ChIP support: the CTCF-only class plus
  any class flagged CTCF-bound. CTCF-only alone leaves out every CTCF site that
  overlaps a promoter or an enhancer, including all of those inside HOXD.
- **Strand.** A site needs MA0139 hits, at 0.85, on one strand only. An element
  with hits on both strands is skipped, as in the orientation measurement.
- **Boundaries.** One boundary midway between each reverse site and the forward
  site after it, merged within 5 kb, with the old caller's 50 kb minimum node.
- **Variant.** `oriented_ctcf_only` applies the same rule to CTCF-only sites,
  to show what the site set does.

`scripts/oriented_domains.py` runs all four measurements on the same
chromosomes (`domains_oriented_comparison`).

| measurement | CTCF-only (current) | oriented (new) | oriented, CTCF-only sites |
|---|---|---|---|
| nodes, genome | 20,002 | 17,995 | 3,273 |
| **1. Hi-C: edges within 20 kb of a measured boundary (enrichment over random)** | | | |
| H1 | 24.7% (1.38) | **31.0% (1.74)** | 26.1% (1.46) |
| K562 | 15.4% (1.36) | **22.5% (1.96)** | 18.1% (1.63) |
| HepG2 | 15.8% (1.30) | **24.3% (1.99)** | 20.0% (1.65) |
| IMR-90 | 16.3% (1.16) | **20.2% (1.45)** | 18.3% (1.35) |
| GM12878 (saturated) | 59.3% (1.10) | 59.8% (1.11) | 58.8% (1.09) |
| measured boundaries reached, H1 / K562 / HepG2 / IMR-90 / GM12878 | 36.6 / 36.0 / 34.5 / 31.1 / 29.0% | **41.2 / 47.5 / 47.9 / 34.8** / 26.3% | 6 to 7% |
| **2. Node content: scored elements whose most-moved coding gene is in their node** | | | |
| share of 113,399 | **81.7%** | 63.2% | 93.8% |
| as many boundaries placed at random | 79.1% | 77.7% | 94.5% |
| excess over random | **+2.6 points** | −14.5 points | −0.7 points |
| both callers cut to the same 5,542 boundaries | **82.2%** | 67.4% | not applicable |
| **3. Mouse synteny: mouse nodes in one human neighbourhood (MGI)** | | | |
| chr19 | **95.9%** of 122 (55.7% in one node) | 91.3% of 138 (31.2%) | 96.5% of 85 (60.0%) |
| chr11 | **92.4%** of 331 (58.3%) | 88.0% of 343 (27.1%) | 86.7% of 210 (59.0%) |
| **4. HOXD: the published HOXD11 to HOXD13 boundary** | | | |
| edge in 176,096,240 to 176,109,754 | no (nearest 106 kb away) | no (nearest 16 kb away) | no |
| nodes over the nine genes | 1 | 2, split HOXD8 / HOXD4 | 1 |

**Reading, measurement by measurement.**

1. **Hi-C.** The new caller wins in every biosample with resolving power, by a
   wide margin: enrichment 1.74 to 1.99 against 1.30 to 1.38 in H1, K562 and
   HepG2. It also reaches more measured boundaries (47.5% against 36.0% in
   K562). GM12878 is saturated, as before.
2. **Node content.** The new caller loses outright. Its nodes separate an
   element from the gene its deletion moves more often than randomly placed
   boundaries would: 63.2% inside against 77.7% at random. At the same
   resolution it keeps 67.4% against the old caller's 82.2%. So it is not
   finer nodes but where the boundaries go. Flips among CTCF sites at promoters
   and enhancers fall between an element and its gene.
   - The old caller's own excess over random is only 2.6 points, on all
     scored elements rather than the 4,800 sampled distal enhancers behind
     the 90.2%. The node-content claim is weaker on the full archive than it
     read on the sample.
   - The archive covers chr15 to chr22 and chrY completely, and chr1 and chr14
     in part.
3. **Mouse.** The new caller loses. The same neighbourhood holds at 88 to 91%
   against 92 to 96%, and "one node" halves (27 to 31% against 56 to 58%).
   These numbers use each caller's human nodes with full gene lists, so the old
   caller's chr19 and chr11 figures differ slightly from the committed 94.5% and
   92.6%, which came from a 12-gene-per-node index.
4. **HOXD.** Neither caller recovers the boundary. The new one does split the
   cluster into two landscapes, but between HOXD8 and HOXD4 rather than between
   HOXD11 and HOXD13. The rule itself places a flip at 176,103,509, inside the
   published interval between HOXD12 and HOXD11. The 50 kb minimum node
   inherited from the old caller then removes it, because the flip at
   176,080,249 came first. Lowering that floor after seeing HOXD would be
   fitting the test, and more boundaries would cost node content further.

**Verdict: keep both, and the default stays CTCF-only.** The new caller wins on
Hi-C, loses on node content and mouse synteny, and does not recover HOXD. The
disagreement is informative in itself. The boundaries insulation maps see
(convergent-loop anchors at CTCF sites in gene-dense regions) are not the
boundaries that bound the deletion model's enhancer reach, which stays within
larger, gene-sparse-edged nodes. A node set for Hi-C-like questions and a node
set for enhancer-to-gene questions may need to be different objects.

Two open points:

- Node content is the model's reading, not a measured enhancer-gene set, so
  the second test inherits AlphaGenome's window and training.
- Whether a stricter site call would move the verdict is not tested here. The
  scan counts an element as both-strand when a weak off-by-two hit sits on the
  other strand. Taking the best hit's strand instead, or requiring stronger
  motifs, is the refinement to try, judged on all four measurements at once
  rather than tuned on one.

## Reader v1 (built 2026-09-11)

`genomeos reader --cell-type K562 --versus HepG2 --chrom chr21` is the first
reader: ENCODE's DNase-seq peaks for a cell type (one 1-2 MB narrowPeak file
from the portal, rows of the chromosome kept) laid over the nodes, the
promoters and the enhancers. Per node: peaks and open fraction; per coding
gene: read (promoter open, TSS ± 1 kb) or silent; per enhancer: active or
not. Two cell types compared give the genes one reads and the other does
not.

| chr21 | K562 (blood) | HepG2 (liver) |
|---|---|---|
| peaks kept | 5,159 | 1,849 |
| coding genes read | 135 of 221 (61%) | 118 of 221 (53%) |
| enhancers active | 2,005 of 12,139 (17%) | 1,401 (12%) |
| silent nodes | 23 of 228 | 34 of 228 |
| read only here | RUNX1, ITGB2, S100B, GRIK1, KCNJ6 … | TFF1, TFF3, ABCG1, FTCD, MX2 … |

RUNX1 and ITGB2 open in the blood line and TFF1/TFF3 and FTCD in the liver
line is what the biology says; 106 genes are read in both. Evidence:
experimental for the peaks, inferred for "read" (an open promoter is
necessary for transcription, not proof of it). Results:
`reader_<cell>_<chrom>.json`, `reader_K562_vs_HepG2_chr21.json`.

### Eleven cell types, every chromosome (2026-09-11)

`scripts/reader_genome_wide.py` now takes any list of ENCODE biosamples (the
reader resolves the released GRCh38 DNase peak file by name) and resumes per
cell type; the default list adds nine to K562 and HepG2. Genome-wide, coding
genes whose promoter is open (read) out of 20,094:

| cell type | read | share | active enhancers | silent nodes |
|---|---|---|---|---|
| hepatocyte | 15,489 | 77.1% | 197,561 | 382 |
| cardiac muscle cell | 15,268 | 76.0% | 233,853 | 775 |
| H1 (embryonic stem) | 15,129 | 75.3% | 133,114 | 1,351 |
| SK-N-SH (neuroblastoma) | 14,883 | 74.1% | 143,197 | 850 |
| K562 (erythroleukaemia) | 14,828 | 73.8% | 193,255 | 2,085 |
| IMR-90 (lung fibroblast) | 14,057 | 70.0% | 213,753 | 1,411 |
| astrocyte | 13,885 | 69.1% | 203,756 | 743 |
| CD14-positive monocyte | 13,647 | 67.9% | 141,431 | 2,989 |
| HepG2 (liver cancer) | 13,015 | 64.8% | 124,478 | 1,465 |
| GM12878 (B lymphoblastoid) | 9,802 | 48.8% | 57,321 | 4,933 |
| keratinocyte | 8,432 | 42.0% | 63,062 | 4,259 |

Read with the caveat that peak counts differ by experiment depth as much as
by biology: GM12878 and keratinocyte have the fewest peaks in their files,
so their low shares are partly the assay. What holds across all eleven is
the shape: two thirds to three quarters of coding genes have an open
promoter in any one cell, the rest is the cell's identity, and the silent
nodes (whole CTCF domains without a peak) are where the reader is not
looking at all. Every cell type is a reader lane on the Blocks tab and a
column on the Progress tab; chr21's per-cell results are committed, the
other chromosomes stay local.

### The reader lane in the block map (2026-09-11)

The Blocks tab has a `reader` selector listing every cell type read on the
loaded chromosome. With one selected, each node in the domain lane is filled
in proportion to its open fraction in that cell (full at 5% of bases under a
DNase peak) and a silent node gets a red edge; every coding gene carries a
dot, filled when the reader calls it read (promoter open) and hollow when
silent. Switching K562 to HepG2 on chr21 is the picture of the whole idea:
the text does not change, the reading does. The node and gene attributes
(`K562_open_fraction`, `K562_node`, `HepG2_read`, …) travel with the blocks
(`/api/blocks`), so any view can ask which cell reads a gene. Inferred, as
the reader itself: an open promoter is necessary for transcription, not
proof of it.

### Every chromosome (`scripts/reader_genome_wide.py`, 110 s after the peaks are fetched)

| genome-wide | K562 (blood, female) | HepG2 (liver, male) |
|---|---|---|
| coding genes read | 14,828 of 20,094 (74%) | 13,015 of 20,094 (65%) |
| enhancers active | 193,255 | 124,478 |
| silent nodes | 2,085 | 1,465 |
| read in both | 11,987 | |

Two checks fall out of the numbers. K562 reads no gene on chrY (0 of 61)
because the line is female; HepG2, male, reads 9. And the two cell types
read three quarters and two thirds of the coding genes respectively, with
11,987 in common: the housekeeping core plus what each lineage adds.
Per-chromosome results stay local (`reader_<cell>_<chrom>.json`); the
summary is `reader_genome_wide.json`.

### Openness is not reading: poised and marked promoters (2026-09-14)

The epigenome layer caught the reader in two wrong answers. In H1 every HOXA
promoter is open, and ten of eleven carry H3K4me3 with H3K27me3. Those genes
are poised, not read. In keratinocytes and GM12878 the DNase files are shallow,
so promoters that carry H3K4me3 and H3K27ac, the marks of an active promoter,
come out closed. `read_chromosome` now uses the marks where the layer has
them, which is every chromosome for all eleven cell types:

- An open promoter with an H3K27me3 peak and no H3K27ac peak is **poised** and
  is not read.
- A closed promoter with H3K4me3 and H3K27ac peaks is **read by its marks**.
- `silent_genes` lists every gene not read, closed or poised, in full. The
  list used to be cut at 200 while the Blocks lane, the gene report and the
  decompiler take "not in the list" as read, so every silent gene past the
  200th on a large chromosome showed as read.

| biosample | read by openness alone | poised | read by marks | read now |
|---|---|---|---|---|
| K562 | 14,828 | 2,135 | 82 | 12,775 |
| HepG2 | 13,015 | 1,553 | 179 | 11,641 |
| GM12878 | 9,802 | 200 | 1,096 | 10,698 |
| H1 | 15,129 | 3,151 | 64 | 12,042 |
| IMR-90 | 14,057 | 1,436 | 116 | 12,737 |
| SK-N-SH | 14,883 | 1,467 | 140 | 13,556 |
| cardiac muscle cell | 15,268 | 2,422 | 105 | 12,951 |
| keratinocyte | 8,432 | 304 | 2,836 | 10,964 |
| hepatocyte | 15,489 | 1,414 | 60 | 14,135 |
| astrocyte | 13,885 | 1,348 | 102 | 12,639 |
| CD14-positive monocyte | 13,647 | 2,434 | 138 | 11,351 |

H1's 3,151 poised promoters are the size of the bivalent set reported for
embryonic stem cells. Keratinocyte and GM12878, the two lowest read shares
blamed on assay depth above, gain 2,836 and 1,096 genes from their marks.

**Checked against measured RNA** (`reader_poised_check`, scripts/epigenome.py
reader-check). A gene counts as expressed when ENCODE total RNA-seq covers 30%
of its exons on its strand, the segment filter's threshold. This covers every
chromosome in the four biosamples with such tracks:

| | K562 | HepG2 | GM12878 | IMR-90 |
|---|---|---|---|---|
| read by openness, expressed | 64.8% | 61.8% | 87.4% | 51.6% |
| read by marks, expressed | 58.5% (82) | 57.0% (179) | 91.2% (1,095) | 40.5% (116) |
| poised, expressed | 6.5% | 7.9% | 16.5% | 6.8% |
| closed, expressed | 8.0% | 7.8% | 19.5% | 7.5% |
| precision of "read": openness only → now | 56.4 → 64.8% | 55.4 → 61.8% | 86.0 → 87.8% | 47.0 → 51.5% |
| expressed genes called read: openness only → now | 94.8 → 93.8% | 91.9 → 91.7% | 75.2 → 83.8% | 93.1 → 92.4% |

Poised genes are expressed as rarely as closed ones, and genes read by their
marks about as often as genes read by openness. Removing the poised genes buys
6 to 8 points of precision for under one point of recall, and the marks buy
GM12878 8.6 points of recall. "Read" still means open or marked, which is
necessary for transcription, not proof of it. Half of IMR-90's read genes show
no RNA over their exons at this threshold.

## The epigenome layer: openness is not function (2026-09-14)

Reader v1 says whether a locus is open. It cannot say whether an open element
carries an active-enhancer mark in one cell and a Polycomb mark in another,
and the registry classes that stood in for that were built by someone else
from marks we never read. `genome/epigenome.py` reads them.

**What is measured.** For each of the reader's eleven biosamples, one released
GRCh38 ENCODE experiment per mark, chosen by a fixed rule and pinned in
`data/results/epigenome_manifest.json`: H3K4me3 (promoter), H3K27ac (active
enhancer or promoter), H3K4me1 (enhancer, primed or active), H3K27me3
(Polycomb repression), H3K9me3 (constitutive heterochromatin), and whole-genome
bisulfite sequencing. The rule keeps one biosample per cell type where the
portal offers several. Cardiac muscle cell is all RUES2-derived, where the
first pass had mixed in H7-derived H3K4me3 and H3K27me3. CD14-positive
monocyte is all one ENCODE donor, where the first pass mixed two Roadmap
donors. After that: complete files first (replicated peaks and pooled fold
change), then the largest peak file. Peaks are streamed once and kept per
chromosome. Signal is read by range into 200 bp bin means. Methylation is read
by range from the per-CpG bedMethyl bigBed into 200 bp bins (calls, calls with
5+ reads, their fractions and reads). All three caches live under
`data/knowledge/epigenome`; nothing raw is kept. The per-strand methylation
bigWigs ENCODE also releases carry every cytosine rather than CpGs, so a mean
over them would dilute the fraction. The bigBed is the clean source, read with
a bigBed reader built on the bigWig reader's R-tree.

**What the profiles cost, and why the broad job ran this time.** One mark's
fold change, for one cell type, on one chromosome, is 25 to 27 range requests
over **two** connections: chr21 is 233,550 bins, 18.5 MB and 22.0 s; chr20 is
322,221 bins, 21.8 MB and 34.3 s. The broad job was deferred once as too slow,
and the handshake was the reason: before the shared range reader kept its
connections open (46cb5d2) each of those 25 requests opened its own TLS
connection, at 0.5 s on a quiet server and 4 to 17 s under load, which is 12 to
425 s of handshake per profile on top of the transfer. Reusing the connection
leaves the bytes, and the bytes are what dominates now, at 0.6 to 0.8 MB/s a
stream and six to ten profiles a minute across eight workers. Nothing about the
job got cleverer; the reader stopped paying for the same connection 25 times.

Run to completion, that is all 1,320 profiles — eleven biosamples × five marks
× 24 chromosomes — for 2.2 GB of 200 bp bin means kept from roughly 55 GB
streamed. The last 593 took 59 minutes. **Every chromosome now carries fold
change as well as peaks, in every biosample and every mark**, so the layer no
longer answers UNKNOWN with "fold-change profile not read for this chromosome";
`genomeos epigenome coverage` says so. Methylation is complete on all 24
chromosomes for the eight biosamples that have GRCh38 WGBS; the other three
have none at all, which no amount of fetching fixes.

| biosample | DNase | H3K4me3 | H3K27ac | H3K4me1 | H3K27me3 | H3K9me3 | WGBS |
|---|---|---|---|---|---|---|---|
| K562 | yes | yes | yes | yes | yes | yes | yes |
| HepG2 | yes | yes | yes | yes | yes | yes | yes |
| GM12878 | yes | yes | yes | yes | yes | yes | yes |
| H1 | yes | yes | yes (Roadmap) | yes | yes (Roadmap) | yes | yes |
| IMR-90 | yes | Roadmap | Roadmap | Roadmap | Roadmap | Roadmap | Roadmap, low depth |
| SK-N-SH | yes | yes | yes | yes | yes | yes | yes |
| cardiac muscle cell (RUES2) | yes | yes | yes | yes | yes | yes | **none** (RRBS on hg19 only) |
| keratinocyte | yes | yes | yes | yes | yes | yes | **none** (no RRBS either) |
| hepatocyte (H9-derived) | yes | yes | yes | yes | yes | yes | yes |
| astrocyte | yes | yes | yes | yes | yes | yes | **none** (RRBS on hg19 only) |
| CD14-positive monocyte | yes | yes | yes | yes | yes | yes | yes, a different (Roadmap) donor |

Eleven of eleven carry all five marks and eight of eleven carry methylation.
The panel is thin where it matters most: there is no neural tissue beyond the
SK-N-SH neuroblastoma line, no gonad and no embryo beyond H1. Four of the
eleven are cancer or transformed lines (K562, HepG2, SK-N-SH, GM12878), and
all four are hypomethylated. Their fossil-tier means are 0.05, 0.18, 0.29 and
0.27, against 0.81 to 0.87 in H1, hepatocyte and monocytes. At 5+ reads,
IMR-90's WGBS measures only 8% of 1 kb windows and the monocyte WGBS 24%,
against 85 to 93% for the other six. The DNase file, the marks and the WGBS are not always from the same
donor or derivation (hepatocyte marks are H9-derived; the DNase biosample is
recorded in the manifest).

**The record.** `epigenome_at(chrom, start, end, cell_types)` returns one
record per cell type:

- `openness`: the DNase call.
- `marks`: per mark, the peak call, peak signal and fold change, with the
  experiment and file.
- `methylation`: the fraction over calls with 5+ reads, the number of calls
  and the mean coverage.
- `state`: inferred from the marks by fixed rules (bivalent, active promoter,
  active or poised or primed enhancer, Polycomb-repressed, heterochromatin, no
  mark).

Every field carries its evidence and confidence (0.9 for replicated peaks with
pooled signal, 0.9 or 0.6 for methylation by depth, 0.5 for the inferred
state). A field is UNKNOWN, with the reason, when the biosample has no
experiment or the chromosome was not read. `genomeos epigenome coverage` prints
the table above, `genomeos epigenome at --chrom chr7 --gene HOXA9` prints a
locus, and `genomeos epigenome summary --chrom C` writes `epigenome_<chrom>`:
marks on read and silent promoters, marks on each registry class, promoter
methylation, and methylation per budget tier. The layer is data only. No
BioIR construct reads it yet; the language change is the engine's decision.

Genome-wide (`epigenome_genome_wide`), 76 to 95% of the promoters the DNase
reader calls read carry H3K4me3 in every biosample, which is the check that
the two layers agree where they should. Of the silent ones, 3 to 4% carry
H3K27me3 in H1 and hepatocyte, 10 to 20% in SK-N-SH, cardiac muscle, astrocyte
and GM12878, and 25 to 34% in keratinocyte, IMR-90, K562, monocytes and HepG2.

HOXA9's promoter through it: bivalent in H1 (H3K4me3 peak and H3K27me3 at 34
times control), closed and Polycomb-repressed in K562, IMR-90 and GM12878,
bivalent and closed in keratinocytes.

### Does mark identity beat the registry class? (chr21 and chr22)

The deletion chain scored every enhancer-like element of chr21 and chr22 per
cell line. That is 83,814 element-cell pairs over 20,954 elements in K562,
HepG2, GM12878 and IMR-90, all four with every mark and WGBS, with no model
call needed. For each pair, a ridge model on standardised features was scored
out of fold (5 folds grouped by element). Three questions were asked: whether
the deletion moves a gene by 0.1 log2 or more, whether the effect is a rise
(silencer-like) among the 16,613 pairs that act, and how large the effect is.
Marks and methylation are read over the element midpoint ± 500 bp.
(`epigenome_direction_chr21_chr22`)

| features | acts (AUC) | rise among acting (AUC) | magnitude (Spearman) |
|---|---|---|---|
| registry class (pELS/dELS, CTCF-bound) | 0.626 | **0.532** | 0.183 |
| + DNase in that line (reader v1) | 0.645 | 0.594 | 0.217 |
| + the line's five marks and methylation | **0.676** | **0.616** | **0.333** |
| control: marks from another of the four lines | 0.654 | 0.605 | 0.253 |
| control: marks shuffled within chromosome × line × class × DNase | 0.645 | 0.593 | 0.216 |
| marks and methylation without class or DNase | 0.662 | 0.603 | 0.326 |

Gains of the line's own marks over registry plus DNase, as 95% intervals over
200 element resamples:

| question | gain in AUC or Spearman | gain over marks from another line |
|---|---|---|
| acts | +0.028 to +0.035 | +0.020 to +0.025 |
| direction | +0.013 to +0.030 | +0.005 to +0.019 |
| magnitude | +0.108 to +0.123 | not recorded |

No resample showed no gain in any of these comparisons.

The registry class says almost nothing about direction: 33% of acting pELS
and 38% of acting dELS rise. The measured state separates more:

| state in that line | rise among acting |
|---|---|
| active promoter | 24% |
| active enhancer | 25% |
| primed | 38% |
| no mark | 41% |
| poised | 44% |
| bivalent | 45% |
| Polycomb-repressed | 46% |

Elements with an H3K27me3 peak rise 45% of the time, against 35% without one.
So mark identity beats the class, and the gain is the cell's own marks:
another line's marks and shuffled marks both fall short. The gain in direction
is small (AUC 0.53 to 0.62), and most of it is openness plus the marks. The
larger gain is magnitude.

Two caveats limit what this means. The outcome is AlphaGenome's prediction,
and AlphaGenome was trained on ENCODE histone ChIP-seq, DNase and RNA-seq in
these same lines, so marks predicting its direction partly measure what the
model learned from those tracks. And "rise on deletion" is the model's
silencer call, not a measured one.

### Does the gain transfer, and does it survive the model's own sibling lines?

The comparison above was fitted and scored on the same two chromosomes. With
the fold-change profiles now genome-wide it can be fitted on chr21 and chr22
and scored on chromosomes it never saw, and the marks can be made to compete
with the thing that actually predicts the outcome: the same element's predicted
effect in the *other three* lines. A separate lane measured that at 0.865 AUC
for direction, against 0.616 for the marks and 0.686 for behaviour type. The
sibling features are the signed mean, the net sign and the mean absolute effect
over the other three lines; the line being predicted never enters its own
features (`scripts/epigenome.py direction-transfer`,
`epigenome_direction_transfer`).

Held out: chr14 to chr20, 564,824 element-line units over 141,209 elements, of
which 115,096 act. Three chromosomes were refused, each with its reason on the
record: chr1 (the sweep has scored 1.6% of it), chr13 (being written while this
ran) and chrY (two of the four lines are female, so their chrY marks and
predicted effects are not about a chromosome the cell has). A sweep in progress
is not a held-out chromosome; it is the first per cent of one.

| features, fitted on chr21 and chr22 | acts (AUC) | rise among acting (AUC) | magnitude (Spearman) |
|---|---|---|---|
| registry class + DNase in that line | 0.660 | 0.601 | 0.234 |
| + the line's five marks and methylation | 0.687 | 0.626 | 0.328 |
| control: marks from another of the four lines | 0.667 | 0.614 | 0.274 |
| control: marks shuffled within chromosome × line × class × DNase | 0.660 | 0.601 | 0.234 |
| **the model's own effect in the other three lines** | 0.787 | **0.907** | 0.475 |
| + registry and DNase | 0.791 | 0.912 | 0.470 |
| + the line's marks as well | 0.792 | **0.915** | 0.487 |
| control: + another line's marks instead | 0.788 | 0.913 | 0.471 |

95% intervals of each difference over 200 resamples of the held-out elements:

| difference | acts | direction | magnitude |
|---|---|---|---|
| marks over registry + DNase | +0.026 to +0.028 | +0.022 to +0.027 | +0.091 to +0.097 |
| marks over another line's marks | +0.019 to +0.021 | +0.010 to +0.014 | +0.052 to +0.056 |
| marks over shuffled marks | +0.025 to +0.028 | +0.022 to +0.029 | +0.091 to +0.096 |
| marks over the sibling model | +0.0003 to +0.0022 | +0.0027 to +0.0035 | +0.015 to +0.018 |
| the sibling model over the marks model | +0.102 to +0.106 | **+0.282 to +0.291** | +0.139 to +0.145 |

**The gain transfers, at the size it had in sample.** Direction rises from
0.601 to 0.626 on chromosomes the fit never saw, against 0.594 to 0.616 in
sample, and it is positive on all seven held-out chromosomes separately (+0.019
on chr20 to +0.038 on chr19). Shuffled marks score the fit without them, which
is the control working.

**Half of it is not the cell's own chromatin.** Another line's marks carry
0.614 of the 0.626, so the part that is the line's own marks rather than marks
at all is +0.010 to +0.014 AUC. The same pattern holds for acts and magnitude.

**And all of it is small beside the model reading its own lines.** Direction
from the other three lines alone scores 0.907; the entire marks model scores
0.626. Adding the line's marks on top of the sibling model moves direction by
+0.0027 to +0.0035 AUC and acts by +0.0003 to +0.0022, which at this sample
size is measurable and, as a claim about chromatin, is three thousandths.
Pooled AUCs hide two reversals worth stating: the marks made the sibling model
*worse* for acts on chr15 (0.7675 against 0.7680) and chr19 (0.8121 against
0.8165), and worse for magnitude on chr15 (0.4321 against 0.4322) and chr19
(0.5471 against 0.5494); against the sibling columns alone the chr19 losses are
larger (0.8121 against 0.8249, 0.5471 against 0.5667). The marks are not
uniformly additive once the model's own behaviour is in the comparison.

The honest reading is that mark identity does carry information about this
outcome, that it survives two controls and transfers off the training
chromosomes, and that it is nonetheless a rounding error against a feature that
is unambiguously the model reading back its own training tracks. What the first
comparison measured as "marks beat the registry class" is true and small; what
it could not see is that both are far behind the model's agreement with itself.
No measured outcome has been used anywhere in this comparison, so nothing here
is yet a claim about biology.

### Does H3K27me3 explain the loci that behave as repressed? The HOX clusters

The four HOX clusters were read in all eleven biosamples
(`epigenome_hox`): 39 HOX promoters among 4,415 coding promoters on chr2, chr7,
chr12 and chr17. The control is non-HOX promoters on the same chromosomes in
the same GC by CpG-density stratum, because Polycomb targets are CpG-island
promoters and an unmatched control would compare islands with everything.

| biosample | HOX promoters with H3K27me3 | matched control | silent HOX promoters with H3K27me3 | matched silent |
|---|---|---|---|---|
| K562 | 64% | 20% | 82% (17) | 45% |
| HepG2 | 90% | 23% | 96% (24) | 55% |
| GM12878 | 72% | 15% | 79% (33) | 35% |
| H1 | 100% | 26% | 100% (4) | 9% |
| IMR-90 | 67% | 17% | 79% (19) | 48% |
| SK-N-SH | 79% | 12% | 80% (10) | 8% |
| cardiac muscle cell | 95% | 21% | 95% (20) | 28% |
| keratinocyte | 56% | 20% | 65% (34) | 34% |
| hepatocyte | 100% | 13% | 100% (8) | 8% |
| astrocyte | 100% | 16% | 100% (25) | 43% |
| CD14-positive monocyte | 97% | 31% | 100% (21) | 66% |

Across cluster spans, H3K27me3 peaks cover most of the cluster in most cells,
at 1.5 to 34 times control against chromosome means of 0.45 to 1.04. The
exceptions are the textbook ones:

- K562 HOXB: 11% covered, seven of ten promoters active. K562 expresses HOXB.
- Keratinocytes HOXA and HOXC: 16% and 12% covered, seven active promoters
  each. Skin expresses both clusters.
- IMR-90 HOXA: six active promoters. Lung fibroblasts express anterior HOXA.
- H1: every HOXA and HOXD promoter is bivalent (H3K4me3 with H3K27me3 at 29
  to 34 times control). The DNase reader calls 11 of 11 HOXA promoters "read"
  there, because bivalent promoters are open. That is the openness-is-not-function
  error the layer exists to catch. Since the same day the reader calls them poised
  ("Openness is not reading", above).

So at HOX the answer is yes, with small numbers: 39 promoters, no formal test,
a matched control two to eight times lower in every biosample.

Across all silent promoters of those chromosomes, the answer is no.
H3K27me3 is the first mechanism in 4 to 39% of them. CpG methylation of 0.5 or
more covers 88% in H1, 90% in hepatocyte, 56% in SK-N-SH and 44% in monocytes.
"None of the three" is 1 to 5% where the methylome is intact, and 24 to 46% in
the hypomethylated lines (K562, HepG2, GM12878, SK-N-SH). A methylated closed
promoter with few CpGs is the default state rather than an explanation, so the
honest reading is narrower: H3K27me3 marks the developmental CpG-island loci
held silent, and it is not the general account of a closed promoter.

### Does methylation explain the fossil tier?

The genome's budget tiers its UNKNOWN space. The fossil tier is 328.5 Mb of
unconstrained transposable-element remains in 7,302 blocks. Every block of
every tier on all 24 chromosomes was tiled into fixed 1 kb windows: 817,079
windows, 320,656 of them fossil. Each window's methylation (calls with 5+
reads) was compared inside GC by CpG-density strata. The fossil-minus-tier
difference is weighted by fossil windows, and 99.8 to 99.96% of fossil
windows fall in shared strata. Blocks, not windows, were bootstrapped
(`epigenome_fossil`).

| biosample | fossil, raw | neutral, raw | fossil minus neutral, matched (95%) | fossil minus regulatory | fossil windows measured | fossil CpG calls covered |
|---|---|---|---|---|---|---|
| H1 | 0.869 | 0.874 | −0.004 (−0.007 to −0.001) | +0.010 | 91% | 62% |
| hepatocyte (H9) | 0.808 | 0.824 | −0.015 (−0.019 to −0.010) | −0.022 | 93% | 80% |
| CD14-positive monocyte | 0.865 | 0.858 | +0.008 (+0.004 to +0.012) | +0.003 | 24% | 9% |
| IMR-90 | 0.613 | 0.619 | +0.003 (−0.011 to +0.017) | −0.087 | 8% | 3% |
| K562 | 0.053 | 0.037 | +0.017 (+0.012 to +0.022) | −0.065 | 85% | 51% |
| HepG2 | 0.176 | 0.152 | +0.025 (+0.019 to +0.031) | −0.099 | 91% | 58% |
| GM12878 | 0.272 | 0.248 | +0.027 (+0.022 to +0.033) | −0.100 | 92% | 71% |
| SK-N-SH | 0.287 | 0.318 | −0.026 (−0.047 to −0.005) | −0.196 | 92% | 65% |

**Negative.** Matched for GC and CpG density, the fossil tier is no more
methylated than unique unconstrained sequence. In the three biosamples with an
intact, well-covered methylome, the difference is −1.5 to +0.8 points around a
background of 81 to 87%: fossils are methylated because the genome is. In the
hypomethylated lines the sign is split (+1.7 to +2.7 points in three, −2.6 in
SK-N-SH). The fossil tier is also no more covered by H3K9me3 peaks than
matched neutral sequence (−0.04 to +0.02 in every biosample). The coverage is
real: all 328.5 Mb of the tier was read, and 85 to 93% of its windows are
measured in six of the eight WGBS biosamples. It is not a coverage artefact,
but neither IMR-90 (8%) nor the monocyte WGBS (24%) can carry the conclusion
alone.

One family looked different. SINE remains (mostly old Alus, 14,000 to 15,000
windows) sit 15 to 18 points above matched neutral sequence in K562, HepG2,
GM12878 and SK-N-SH, and 0 to 3 points above it in H1, hepatocyte and
monocytes.

**The Alu lead, followed (`epigenome_fossil_alu`).** The comparison takes
125,219 RepeatMasker records of 200 bp or more lying wholly inside fossil-tier
blocks: 61,068 Alu and 64,151 L1. Each element's methylation is read from the
200 bp bins at least half inside it, so a flank can pull the value toward its
neighbourhood. The two families are not alike as sequence:

| family | median GC | median CpG per 100 bp | median divergence | median solo-WCGW share of CpGs |
|---|---|---|---|---|
| Alu | 0.52 | 2.24 | 0.11 | 0.13 |
| L1 | 0.35 | 0.53 | 0.18 | 0.43 |

A solo-WCGW CpG has no other CpG nearby and an A or T on each side. It is the
context hypomethylated genomes lose first (Zhou et al. 2018). Alu minus L1, in
points of methylation, weighted by Alu elements:

| design | K562 | HepG2 | GM12878 | SK-N-SH | H1 | hepatocyte | monocytes | IMR-90 |
|---|---|---|---|---|---|---|---|---|
| raw | +8.4 | +19.2 | +15.2 | +18.4 | +3.0 | +7.5 | +5.2 | +18.1 |
| GC × CpG matched | +7.1 | +5.0 | +4.0 | +10.8 | +2.5 | +5.6 | +4.7 | +9.8 |
| + divergence (age) | +8.8 | +6.6 | +4.6 | +12.6 | +1.5 | +4.1 | +1.7 | +8.0 |
| + solo-WCGW share | +8.3 | +6.7 | +4.6 | +12.7 | +1.6 | +4.3 | +1.8 | +9.9 |
| same 50 kb window, no sequence match | +2.1 | +11.0 | +8.5 | +9.4 | +3.4 | +6.1 | +4.2 | too few |

Matching on sequence removes most of the raw excess in HepG2 and GM12878
(19 to 7 points, 15 to 5). A residual of 4.6 to 12.7 points stays in the
hypomethylated lines, against 1.6 to 4.3 in the intact ones. But the matched L1 set
is a few thousand atypical, CpG-rich L1s (3,600 to 4,300 of 64,000), because
the families barely overlap in sequence space.

Inside each family, methylation falls steeply with the solo-WCGW share in the
hypomethylated lines: HepG2 Alu from 0.40 to 0.11, HepG2 L1 from 0.22 to 0.10.
It is nearly flat in H1 (0.90 to 0.86). At high solo-WCGW share Alu and L1 are
alike (0.107 against 0.102 in HepG2). The design that holds the region and the
sequence together leaves only 13 to 61 elements per line, too few to read.

So most of the Alu excess is sequence context, the known way hypomethylated
genomes lose CpG methylation, and some is region. A family effect beyond
context is not shown, and this data cannot separate it, because Alu and L1 do
not share enough sequence. **Not a family program. The lead is closed at this
resolution.** A per-CpG reading (the bigBed, not 200 bp bins) comparing
Alu-internal CpGs with L1-internal CpGs of the same flanking context would be
the next resolution if anyone needs it.

For the lexicon's negative (area I: the fossil tier is not a library for its
node), this closes the obvious mechanism too. Methylation keeps the fossil
tier quiet exactly as much as it keeps the rest of the unconstrained genome
quiet, and no more.

## The first chromatin mechanism the VM runs: CpG methylation across divisions (2026-09-16)

Everything above reads methylation. This runs it. `genomeos/runtime/methylation.py`
is an engine, not an annotation: per CpG dyad it holds three states
(unmethylated, hemimethylated, methylated) and steps them through one cell
division at a time — replication makes a methylated site hemimethylated,
UHRF1 plus DNMT1 restore it with a fidelity below 1 (Bostick 2007, Sharif 2007),
DNMT3A/B methylate de novo (Okano 1999), TET1-3 erase (Tahiliani 2009, Ito
2011) — with **both maintenance and de novo efficiency rising with the local CpG
density**, which is the published account of why an isolated CpG erodes and a
clustered one does not (collaborative models: Haerter 2014, Lövkvist 2016; the
genomic observation is Zhou 2018). The expectation and the steady state are
propagated exactly through the per-division transition matrix; individual sites
are sampled from a seeded generator when a lineage rather than a mean is wanted.

Six of the eight rates are not established per CpG per division, so **`param x =
unknown` is now BioLang**: `bio.std.methylation` states maintenance fidelity
(roughly 0.95–0.99: Laird 2004, Genereux 2005; above 0.99 clonally: Ushijima
2003) and the 35 bp solo definition, and says `unknown` for the two solo
fidelities, both de novo rates, TET turnover and the neighbour scale. The engine
**refuses to run** until a program binds them, and `bio check` lists what is
still UNKNOWN. `data/demo/methylation_erosion.bio` binds them as its own
labelled assumptions and tests what it predicts: from 0.95 at 100 divisions,
solo-WCGW falls to 0.17 while dense holds 0.77.

### The claim, registered before the data

Committed in 9a2860c, before any methylation value was read, with the window
classes fixed from **sequence alone**: solo-WCGW windows are 200 bp bins whose
every CpG is isolated (no other CpG within 35 bp) and A/T-flanked; dense windows
hold six or more CpGs, none isolated, and are not island-like (CpG obs/exp <
0.6, because islands are unmethylated in every cell type and would not be a
control). On chr19 to chr22 that is 2,976,042 CpGs, 231,295 of them solo-WCGW,
giving **11,257 solo-WCGW windows and 53,522 dense windows**. The prediction:
the dense-minus-solo gap is larger in the long-cultured lines (K562, HepG2,
GM12878) than in the three biosamples with an intact methylome (H1, hepatocyte,
CD14-positive monocyte). Falsified if any cultured line's gap is below any
intact one's. (`scripts/methylation_erosion.py`, `methylation_erosion`.)

| biosample | solo-WCGW | dense | gap | gap / dense | gap within 50 kb | solo windows measured |
|---|---|---|---|---|---|---|
| K562 | 0.172 | 0.408 | **+0.236** | 0.579 | +0.091 | 8,872 |
| HepG2 | 0.307 | 0.676 | **+0.368** | 0.545 | +0.231 | 9,970 |
| GM12878 | 0.277 | 0.631 | **+0.354** | 0.561 | +0.274 | 10,153 |
| H1 | 0.815 | 0.832 | +0.017 | 0.020 | +0.038 | 10,035 |
| hepatocyte (H9) | 0.760 | 0.836 | +0.076 | 0.090 | +0.084 | 10,365 |
| CD14-positive monocyte | 0.786 | 0.771 | −0.015 | −0.019 | +0.014 | 1,271 |
| SK-N-SH (outside the claim) | 0.381 | 0.723 | +0.342 | 0.473 | +0.216 | 9,944 |
| IMR-90 (outside the claim) | 0.555 | 0.739 | +0.184 | 0.249 | +0.111 | 371 |

**Not a blind test, and the brief is part of the record.** When this claim was
registered the ordering was already known here: the fossil-tier reading of
2026-09-14 had measured SINE and solo-WCGW loss in the same hypomethylated
lines, and the brief that opened this lane named the expected direction. The
registration still did its work — the window classes were fixed from sequence,
and the numbers were read once — but what the result establishes is that the
engine's context-dependent mechanism is *consistent with* a fact already in
hand, not that it predicted an unknown one. The region-matched arm is where a
test of the mechanism itself would live, and it is thin.

**The ordering holds, on all three metrics.** In points, the smallest cultured
gap (K562, +0.236) is 0.160 above the largest intact one (hepatocyte, +0.076).
Relative to the dense level the separation is wider (0.545–0.579 against
−0.019–0.090, margin 0.454). The region control — the gap recomputed inside 50 kb
windows that hold both classes, 878 to 3,012 such windows per biosample of the
claim (230 in IMR-90, whose WGBS is thin) — keeps
the ordering but **thinly**: K562's +0.091 is only 0.0075 above hepatocyte's
+0.084, so region (late-replicating, AT-rich PMD sequence) carries much of the
raw effect, and one line's margin is inside what the earlier fossil-tier
bootstraps would call noise. The two biosamples outside the pre-registered set
sit where the logic would put them: SK-N-SH, the fourth transformed line, with
the cultured group; IMR-90, a mortal strain grown to senescence, between them.
The monocyte's −0.015 is the only negative raw gap, and its WGBS measures 1,271
of the 11,257 solo windows (2,659 calls), the coverage problem already recorded
for that donor.

### What the fit implies — inferred, not measured

Starting each lineage at H1's own measured levels (0.832 dense, 0.815 solo),
with no TET turnover and de novo at 0.01 (dense) and 0.002 (solo) per strand per
division, the engine's inverses say:

- **At the literature fidelity of 0.97**, the observed dense level needs 123
  divisions in K562, 24 in GM12878, 17 in HepG2, 11 in SK-N-SH, 10 in IMR-90, 6
  in the monocyte and 0 in H1 and hepatocyte. These lines have doubled *thousands*
  of times, so 0.97 per division is far too lossy to describe them: real
  maintenance in these genomes must be nearer the top of the published range.
- Read the other way — fidelity solved at a fixed division count — that is what
  comes out. At 100 divisions: K562 0.968 dense against 0.960 solo (deficit
  0.008), HepG2 0.990/0.976 (0.014), GM12878 0.988/0.973 (0.014). At 500:
  0.972/0.980 (−0.008), 0.992/0.990 (0.002), 0.990/0.989 (0.001).
- Under the higher de novo scenario (0.05/0.01) the observed dense levels are
  **below anything a fidelity of 0.97 can reach**, and the engine returns nothing
  rather than a number — the refusal that `fidelity_for_level` exists to make.

So the ordering is confirmed and the rates are not. The data fixes where each
class sits, and any (divisions, fidelity, de novo) triple on a one-dimensional
family reproduces it; past a few hundred divisions the fidelity deficit the fit
needs falls below 0.008 and can invert, meaning the lower de novo rate of an
isolated CpG would suffice on its own. Everything in this section beyond the
table is inferred.

**And one thing the mechanism does not explain.** The dense class erodes too:
0.832 in H1 against 0.408 in K562 and 0.631 in GM12878. Context dependence
cannot produce that, because a CpG with six close neighbours is the case the
model says is well maintained. The cultured lines have lost methylation
genome-wide — a dosage or domain-level change (DNMT1/UHRF1, PMD-wide) that this
engine has no construct for. The engine's own demo run overshoots for the same
reason: a 0.60 gap at 100 divisions against a measured maximum of 0.37, because
it holds the dense class up while the real genomes let it fall. What was
measured here is that the *ordering* of the context effect tracks culture
history, which is the prediction; a mechanism for the global loss is the next
construct, not this one.

**Beside the twin's clocks.** `genomeos/twin/clocks.py` carries Horvath 2013 (353
CpGs) and Hannum 2013 (71 CpGs): elastic-net regressions on beta values that
predict chronological age and contain no mechanism at all. This engine is the
other kind of object — a mechanism with no fit to age — and the two meet at the
divisions. A clock reads a state; this says how the state changes when a cell
divides. Connecting them (a predicted drift at the clock's own CpGs, in their own
sequence context, against a twin's measured betas) is a real next step and is not
claimed here: none of the numbers above involve a clock.
