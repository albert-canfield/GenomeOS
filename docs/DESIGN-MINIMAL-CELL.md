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
probability. It is not: `organism/forge.py:170` writes the literal `0.3` onto
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

## Budgets in `design`, knockouts in network experiments (registered 2026-09-27)

Area H's two remaining concrete items, written down before the code. The
constants live in `genomeos/organism/forge.py` (`BUDGET_*`) and
`genomeos/forge/network_experiment.py` (`NETWORK_KNOCKOUT_*`).

**Budget in `design`.** The unit is the perturbation: a knockout, an addition
or a continuous knob moved off the program's own value each cost 1, unless the
caller prices a name (bases, edits, money). `at_most` already caps how many
factors are combined, but it never counted knob moves, so a design that varies
timers could move every knob at once. A budget is a hard limit on the answer,
not a loss term: the answer never spends more than it. Candidates over the
limit are still run once at the program's own knob values, so the result can
say what the spend bought and what the limit excluded, including whether an
excluded candidate reached the target more closely than the answer did.
*Falsifier:* on any shipped design and any limit from 0 to 3, an answer over
the limit, an excluded candidate missing from the list, or a cheaper-loss
excluded candidate not flagged. *Known case:* `two_intestinal_founders` in
`data/organisms/celegans/designs.bio` (Lin et al. 1995, pinned to POP-1 by
`tests/test_design.py`): at limit 1 the answer is -POP-1, solved; at limit 0
it is the wild type, unsolved, with -POP-1 listed as excluded at loss 0.

**Network knockouts in `experiment`.** An experiment on a network names nodes
held at zero for the whole run and edges `A>B` (B reads A as zero; every other
reader sees A as it is). The perturbed run is compared with the unperturbed
run from the same start: per node, the change in the fraction of the attractor
spent on (Boolean) or in the late mean level (continuous, as a log2 ratio), and
whether the dynamics changed kind (fixed point against cycle, oscillating
against not). The consequence is the mean change over nodes knocked out in
neither run, and it is ranked against knockouts of the same size, nodes and
edges matched separately, drawn from what was not chosen: 200 draws, or every
matched set when there are fewer; p = (1 + as-large-or-larger) / (1 + draws).
*Falsifier:* a knockout of a node no rule reads must score 0 with p = 1.0, and
no draw may contain a chosen node. *Known case:* the Fauré et al. 2006
mammalian cell cycle (`data/models/mammalian_cell_cycle.bnet`, pinned by
`tests/test_boolean.py`): from CycD on, the unperturbed network cycles, and
CycD held at zero rests in a fixed point with Rb, p27 and Cdh1 on and CycA,
CycB off.

## What the two items do now, and what they cannot tell you (2026-09-28)

Both items above are built and tested against the registration:
`tests/test_design_budget.py` (11 tests) and `tests/test_network_knockout.py`
(9 tests). Both known cases hold as registered — `two_intestinal_founders`
answers -POP-1 at limit 1 and the wild type at limit 0 with -POP-1 named as
excluded at loss 0, and CycD held at zero turns the Fauré cycle into the G1
fixed point with Rb, p27 and Cdh1 on, CycA, CycE and CycB off. The negatives
are the part worth reading.

**A hard limit hides the fix rather than pricing it.** At limit 0 the answer to
`two_intestinal_founders` is the wild type at loss 0.5, and the only trace of
the perturbation that would have solved it is the excluded list
(`limit_cost_the_target`). That is what "a limit, not a loss term" means, and a
caller who sets a limit below the cheapest fix gets an answer that looks like a
failed search unless they read the budget report.

**The chance rank is coarse, and on a small network it is nearly uninformative.**
With one node knocked out in a 10-node model there are nine matched knockouts,
so every p-value is a multiple of 0.1 and the strongest possible single-node
knockout still reads p = 0.1: the floor, not a small number. Cutting one edge
of the three-gene ring destroys the oscillation and still reads p = 1.0,
because every other edge of a ring does the same. The control answers "is this
change unusual among knockouts of this size", and on a small network the honest
answer is usually no.

**A clamp in the GRN runtime does not hold a species at zero.**
`NetworkRuntime.run(clamp=...)` re-applies the clamped value only between whole
steps, so inside the Runge-Kutta stages a clamped mRNA still takes its basal and
regulated rate and its protein is translated from those stages. Measured on the
three-gene ring: clamping `a.mRNA` to 0 left protein A at 1.87 instead of 0,
about a tenth of its unperturbed 19.2 — a knockdown, not a knockout. A network
node knockout therefore zeroes what makes the species (a gene's basal and max
rate, a protein's `produces` rules) instead of clamping it, which holds it at
exactly zero at every stage. Existing callers that clamp an external signal to
a non-zero level carry the smaller version of the same error; nothing in area H
depends on it, and it is not fixed here.
