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
magnitude instead of a distance. The same job runs on any fetched
chromosome.

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
