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
  error the layer exists to catch.

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

One family behaves differently. SINE remains (mostly old Alus, 14,000 to
15,000 windows) sit 15 to 18 points above matched neutral sequence in K562,
HepG2, GM12878 and SK-N-SH, and 0 to 3 points above it in H1, hepatocyte and
monocytes. Where a cell has lost methylation genome-wide, Alu remains keep
theirs. That is a property of a family under demethylation, not of the tier,
and it is the one lead here worth a matched follow-up (Alu against LINE of the
same age and CpG content).

For the lexicon's negative (area I: the fossil tier is not a library for its
node), this closes the obvious mechanism too. Methylation keeps the fossil
tier quiet exactly as much as it keeps the rest of the unconstrained genome
quiet, and no more.
