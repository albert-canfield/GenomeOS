# Can a small extract of the genome build a neuron that multiplies?

Short answer: the *program* that makes a neuron is small; the *machinery* a
neuron runs on is not, and the two properties you want, "neuron" and
"multiplies", are in tension in biology. Here is what the tools say.

## The experiments you are thinking of

- Cortical Labs' DishBrain (Kagan et al. 2022) and its CL1 system, and
  FinalSpark's Neuroplatform, grow human neurons from induced pluripotent stem
  cells on electrode arrays and use them for computation; Brainoware (Cai et
  al. 2023) does the same with brain organoids. The cells are not genome
  extracts: they are full human cells, differentiated, then kept alive.
- Direct conversion shows how small the identity switch is: three factors
  (Ascl1, Brn2, Myt1l) turn fibroblasts into neurons (Vierbuchen et al.
  2010); a single factor, NGN2, turns stem cells into functional neurons in
  days (Zhang et al. 2013).

## What the design tool counts

`genomeos design neuron` composes the libraries and the core-essential gene
set (Hart et al. 2017, 684 genes every human cell line needs):

| tier | genes | what it is |
|---|---|---|
| master switches | 5 | NEUROG2, ASCL1, NEUROD1, MYT1L, POU3F2: the identity program |
| identity libraries | ~2,000 | nervous-system development, synaptic machinery, signalling toolkit (data-driven membership) |
| essential | 684 | replication, transcription, translation, proteostasis: what keeps any cell alive |
| **minimal program** | **~2,700** | essential + identity + switches |
| broad program | ~13,700 | everything the 13 core libraries list, a superset |

and the base budget for the minimal program under the three measured
organisations (docs/GENOME-ANATOMY.md):

| layout | budget |
|---|---|
| compact, mtDNA-like (1.3 kb per gene) | 3.5 Mb |
| dense, worm-like (5 kb per gene) | 13 Mb |
| sparse, human-like (150 kb per gene) | 400 Mb |
| the actual human genome | 3,100 Mb |

So yes: the information needed to specify and run a neuron is on the order
of one percent of the genome, and in a compact layout a few megabases.

## Why it is still not "a small extract"

1. **Essential is a floor, not a design.** The 684 core-essential genes are
   what cell lines cannot lose one at a time in culture; a cell built from
   only them plus the neuron program has never been made. The smallest
   living genome ever built (JCVI-syn3A, 473 genes, 531 kb) is a bacterium
   with no nucleus, no organelles and no ability to become a neuron.
2. **A neuron does not multiply.** Mature neurons are post-mitotic: the same
   programme that makes them (NEUROG2 → NEUROD1) shuts the cell cycle down.
   What multiplies is the *progenitor* (SOX2, PAX6, NES; `genomeos design
   neural_progenitor`). The design is therefore two states and a switch:
   expand as progenitors, then differentiate. Every neural culture system,
   including the ones above, does exactly this.
3. **The genome is not the whole cell.** rRNA and tRNA arrays, centromeres,
   telomeres and replication origins are not in any gene list, and a genome
   must be booted inside an existing cell with ribosomes, membranes and
   mitochondria (docs/GENOME-AS-CODE.md, "the runtime that reads the code").
4. **Sparse layout is not waste we understand.** The human 150 kb per gene
   carries regulation we cannot yet read; a compact rewrite would discard it
   before we know what it does. This is the UNKNOWN space in the block map.

## What GenomeOS can do with this now

- Give the parts list and budget per cell type with evidence (`genomeos design`).
- Show the identity switch as BioLang: a `cell_type NeuralProgenitor` that
  expresses SOX2/PAX6, a `cell_type Neuron` that expresses the NEUROG2 targets,
  and an `event differentiate` guarded by NEUROG2 level, run and debugged like
  any other module.
- Point at what is missing: which essential genes lack a library, and which
  identity genes have only `inferred` evidence.

Building the cells is laboratory work outside this project. What the project
adds is the honest count of what a design has to contain and what it would
cost in bases, before anyone tries.

## The confidence on a BioForge design (2026-09-22)

Area H's third open item asked how often a design's answer matches the
published outcome, which presumes the number beside the answer is a
probability. It is not: `organism/forge.py:143` writes the literal `0.3` onto
every feasible answer, ignoring the loss, the feasibility, the evaluations and
the design's own stated confidence, and the same constant marks an unrelated
quantity in area I. A single-valued predictor has one reliability bin at every
sample size, so no calibration curve over it exists at any n. The repository
holds 2 designs with a published outcome cited in the file (4 counting two
recorded only in prose), all four solved, and the answers were written into the
module's rules from the same papers they are scored against. Verdict
`UNCHECKABLE_BY_CONSTRUCTION`, with the registration, the census, the interval
that was computed and withheld, and the five things that would change it, in
docs/BIOFORGE-CONFIDENCE.md.
