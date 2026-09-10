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
