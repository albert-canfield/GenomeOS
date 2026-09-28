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
*(2026-09-28: `16a0434` retired the literal. A design answer now carries a certainty record whose
probability is unavailable with its reason; see BIOFORGE-CONFIDENCE.md.)*

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

**Who reads `clamp=`, before it is fixed (census, 2026-09-28).** Two callers
pass it, and `data/results/grn_clamp_census.json` (written by
`scripts/grn_clamp_census.py before`) lists them with what depends on each.
`run_gastrulation` holds the NODAL signal at a non-zero level per cell; it
feeds `genomeos develop gastrulation`, the web develop view and the pin in
`tests/test_gastrulation.py` (layer order, and each proportion within 0.25 of
the module's stated expectation: mesoderm reads 0.100 against 0.35, inside the
bar by rounding). `Api.cell_run(silence=True)` holds silenced genes' mRNA at
zero, but the runtime already silences those genes and their derivative at
zero is zero, so its stages never left zero and its pin in
`tests/test_cell_view.py` cannot move. `Body._network_step` has its own
single-stage loop that skips clamped species and is not a caller;
`forge/network_experiment.py` avoids `clamp=` on purpose; no BioLang
`experiment` block and no script reaches it. Measured inside the stages of a
5 h run: a clamped-at-zero `a.mRNA` reaches 3.92, a protein clamped at 3
reaches 7.62, a signal clamped at 2 falls to 1.81. `tests/test_grn_clamp.py`
states the fix as strict expected failures: every stage and every recorded
point equals the clamp, and a protein whose mRNA is clamped at zero follows its
own decay and nothing else.

**The clamp now holds inside every stage (2026-09-28).** `NetworkRuntime.run`
masks a clamped species' derivative to zero and writes the clamp value into
every stage state, so it equals its clamp at every stage and every recorded
point; the public signature is unchanged and an unclamped run is bit-identical
(values in `data/results/grn_clamp_census_after.json`).
With `a.mRNA` clamped at zero, protein A on the ring now reads 1.1e-66 (its own
decay from 10) instead of 1.87. Among the census callers only gastrulation
moves, because NODAL no longer decays inside the stages it is read in: at the
test's 60 cells, 30 h, ectoderm 0.483 to 0.467 and endoderm 0.417 to 0.433,
mesoderm unchanged at 0.100, boundaries one cell further out; at the
command-line and web default of 120 cells, 40 h, mesoderm 0.108 to 0.100,
endoderm 0.417 to 0.425. Both move away from the module's stated expectation
(mesoderm 0.35, endoderm 0.20), and the test's mesoderm stays within the 0.25
bar only by floating-point rounding (0.35 - 0.1 = 0.24999999999999997). The
layer order is unchanged; the cell view does not move.

## Gastrulation against a measured census (registered 2026-09-28)

**What the module models.** `data/demo/gastrulation.bio` names human genes
and proteins (GENCODE; UniProt P48431, O15178, Q9H6I2, Q96S42) and cites
mouse and stem-cell papers for its interactions; it names no stage. It is a
generic amniote germ-layer switch read on a 1-D NODAL axis at steady state,
and its three `expected_*` proportions cite nothing.

**The census.** The only measured human gastrula with per-cell labels is
Tyser et al. 2021 (Nature 600:285, E-MTAB-9388): one Carnegie stage 7 embryo
(16–19 days), 1,195 Smart-seq2 cells after QC, FACS-sorted from three
dissected regions (665 caudal, 340 rostral, 190 yolk sac), eleven author
clusters. It counts sampled cells, not the embryo: the region shares are how
many cells were sorted, and n = 1. Its "endodermal cell" cluster holds
hypoblast and yolk-sac endoderm with definitive endoderm, its "ectodermal
cell" cluster holds amnion with non-neural ectoderm (subclusters are not in
the public SDRF), and its primordial germ cells sit inside "primitive
streak". The mouse atlas (Pijuan-Sala et al. 2019, Nature 566:490,
E-MTAB-6967; 116,312 10x cells, pooled whole embryos, E6.5–E8.5) is the
secondary, at E7.0, E7.25 and E7.5, the stages Tyser et al. match to human
CS7 epiblast and streak. Counts are distilled by
`scripts/gastrulation_census.py distill` into
`data/results/gastrulation_census_counts.json`.

**Is it comparable at all.** Only partly, and the registration says so
before the number. The model has three fates and a position; the census has
an undifferentiated epiblast, a primitive streak and extraembryonic tissue
that the model has no state for, and it changes with stage while the model
has none. The model's shares are fractions of axis length, set by two
unsourced gradient numbers (`decay_length`, `nodal_max`), so a pass would
show little; a large miss at default parameters is still informative,
because it says what the default gradient claims about an embryo is wrong.

**Registered comparison** (constants in `genomeos/runtime/gastrulation.py`):

- Model: one run of `run_gastrulation` with the command-line defaults, 120
  cells, 40 h, dt 0.05, `nodal_max` 6, `decay_length` 0.35
  (`CENSUS_MODEL_RUN`). No parameter is changed before or after.
- Human labels never mapped: yolk sac mesoderm, hemogenic endothelial
  progenitor, erythrocyte (extraembryonic or yolk-sac blood, off the
  embryonic axis).
- Two mappings. *differentiated*: ectoderm = ectodermal cell; mesoderm =
  nascent, emergent, advanced and axial mesoderm; endoderm = endodermal
  cell; epiblast and primitive streak left out. *marker*, the model's own
  reading: the SOX2 default holds the epiblast (ectoderm += epiblast cell)
  and TBXT holds the streak (mesoderm += primitive streak).
- Two site filters: all three regions, and the embryonic disc only
  (rostral and caudal), which removes the yolk-sac share of the endoderm
  cluster.
- Each variant gives three shares renormalised over its mapped cells. The
  census interval per layer is the lowest to highest share over the four
  variants; mapping ambiguity is carried by that spread, not by the
  tolerance.
- Tolerance 0.05 beyond the interval: the binomial standard error of a
  share at the smallest mapped total (about 600 cells) is at most 0.02, so
  0.05 is about 2.5 standard errors, and nothing larger is defensible for
  one embryo without inventing an embryo-to-embryo spread.
- Falsifier: any of the three model shares outside its interval by more
  than 0.05 on the human census. The verdict is then "falsified" for the
  default proportions; otherwise "not contradicted", never "validated".
- Mouse: the same two mappings (extraembryonic, blood, PGC, caudal
  epiblast and NMP labels never mapped; *marker* adds Epiblast to
  ectoderm, Primitive Streak to mesoderm, Anterior Primitive Streak to
  endoderm) at the three stages, reported beside it with no verdict.

**What was seen before this was written.** The model's shares (lane-gastrula,
`19423bd`: 0.475 / 0.10 / 0.425 at these defaults) and the distilled census
counts were both open while the mappings were chosen, so this registration
protects against choosing a mapping or tolerance after the verdict, not
against knowing the data. Expected outcome, stated now: falsified on
mesoderm, which dominates the sampled CS7 cells in every mapping.

**Sourced parameters for a later, separately registered lane** (none
applied here):

- `Tbxt inhibits SOX17` (inferred, strength 0.6) has the wrong sign for the
  published evidence: Lolas et al. 2014 (PNAS 111:4478, PMC3970479) find
  *Sox17* a direct Brachyury target that Brachyury activates, in mouse ES
  cell embryoid bodies (ChIP-seq), with Sox17 in turn repressing Brachyury
  (overexpression) and the two anti-correlated in E7.5 embryos. They report
  no dose, threshold or binding strength. Change: `Tbxt activates SOX17`,
  evidence experimental Lolas 2014, strength and threshold still unsourced;
  expected effect: the brake on SOX17 becomes an accelerator, endoderm
  widens and the mesoderm band narrows, away from the census.
- `Sox17 inhibits TBXT` (inferred) gains the same citation; no number.
- NODAL thresholds: the one measured dose ratio is Xenopus, not mammalian.
  Dyson & Gurdon 1998 (Cell 93:557) find about 100 and 300 bound activin
  molecules per animal-cap cell switch on *Xbra* and *Xgsc* (2% and 6% of
  receptors), and Shimizu & Gurdon 1999 (PNAS 96:6791) carry the same
  threefold step to nuclear SMAD2. The model's TBXT and SOX17 thresholds
  (1.0 and 3.0) happen to be threefold apart, but the numbers do not map:
  the units are occupied receptors per cell against model NODAL units with
  no stated relation, the measure is onset rather than a Hill half-maximum,
  and goosecoid is an organizer gene, not SOX17. Vincent et al. 2003 (Genes
  Dev 17:1646) give only an order in the mouse (lowering Nodal/Smad2 loses
  anterior definitive endoderm and prechordal plate first); D'Amour et al.
  2005 (Nat Biotechnol 23:1534) give an activin dose in a dish. Dubrulle et
  al. 2015 (eLife 4:e05042) find in zebrafish that thresholds on Nodal
  concentration do not predict target ranges, induction kinetics do: that
  bears on the model class (a clamped level read at 40 h), not on a number.
  No threshold here is citable as a value.

**The one comparison, run after the registration was committed (`3c4e301`).**
`scripts/gastrulation_census.py compare`, defaults, once; values in
`data/results/gastrulation_census_comparison.json`. Verdict: **falsified**,
on all three layers, not only mesoderm.

| layer | model (120 cells, 40 h) | human CS7 interval over four variants | outside by |
|---|---|---|---|
| ectoderm | 0.475 | 0.046 – 0.172 | 0.303 |
| mesoderm | 0.100 | 0.694 – 0.773 | 0.594 |
| endoderm | 0.425 | 0.116 – 0.213 | 0.212 |

Sampled CS7 cells are mostly mesoderm in every mapping; the model puts most
of its axis in ectoderm and endoderm. The mouse secondary (no verdict) says
the same for mesoderm and endoderm: across E7.0–E7.5 and both mappings,
mesoderm 0.46–0.94 and endoderm 0.03–0.12 against the model's 0.10 and
0.425; the model's ectoderm (0.475) falls inside the mouse range only
because the *marker* mapping counts the large E7.0 epiblast as ectoderm.
Descriptive, not registered: the module's own unsourced `expected_*`
values (0.45 / 0.35 / 0.20) would also miss the human interval on ectoderm
and mesoderm and sit inside it only on endoderm. What the miss means is
bounded by the comparability note above: a sampled cell census of one
embryo is not the shares of a 1-D axis, and the model has no stage. It does
not say the network is wrong; it says the default gradient's germ-layer
shares are not those of a measured gastrula, and that no sourced number
exists yet to set them.

## The Tbxt-SOX17 sign, corrected (registered 2026-09-28)

A correction to a cited fact, registered before the module is changed or
run. It is not a tuning: no number is chosen, and the change is expected to
move the model further from the CS7 census, not closer.

**The source, read in full** (Lolas M, Valenzuela PDT, Tjian R, Liu Z 2014,
PNAS 111:4478–4483, doi:10.1073/pnas.1402612111, PMC3970479; mouse ES-cell
embryoid bodies and E7.5 embryos):

- Brachyury binds and activates *Sox17*: "many key developmental genes
  (Fgf8, Foxa2, Foxj1, Sox17, and Dusp6) were directly targeted by
  Brachyury" (ChIP-seq, Fig. 1A, SI Fig. S8); "Brachyury depletion
  decreased Foxa2 and Sox17 expression at both protein and mRNA levels"
  (Fig. 3A, SI Fig. S9 A–B); the model they favour: "Brachyury first
  targets and promotes the expression of Foxa2 and Sox17".
- Sox17 represses *Brachyury*: "high levels of Sox17 would, in turn,
  directly or indirectly repress the expression of Brachyury (Fig. 3H) ...
  we overexpressed Sox17 after 2 days of differentiation ... we observed
  significant down-regulation of Brachyury at EB day 4 (Fig. 3C)"; in the
  embryo "Brachyury and Sox17 expression was negatively correlated in the
  definite endoderm region" (Fig. 3 D–F). Direct or indirect is left open.
- No dose, binding strength, threshold or Hill coefficient is reported for
  either interaction. Mouse, not human; embryoid bodies with activin, not a
  NODAL gradient.

**The change** (`SIGN_CORRECTION` in `genomeos/runtime/gastrulation.py`):

- `rule Tbxt inhibits SOX17 { strength: 0.6; threshold: 4.0; hill: 3;
  evidence: inferred "mesoderm dampens endoderm" }` becomes `rule Tbxt
  activates SOX17` with evidence experimental Lolas 2014; the old line stays
  in the module as a dated comment.
- `rule Sox17 inhibits TBXT` keeps its sign and numbers; its evidence
  becomes experimental Lolas 2014 (Fig. 3C).
- Strength, threshold and Hill stay exactly as they were (0.6 / 4.0 / 3
  and 1.0 / 2.0 / 3), and so do both confidences (0.4): the paper supports
  the sign and nothing else, and the module comment says so. No value is
  picked to improve any fit, now or after the run.

**Expected direction, stated before the run** (lane-census's prediction):
the brake on SOX17 in TBXT-high cells becomes an accelerator, so the
endoderm band widens and the mesoderm band narrows or vanishes; the model
moves further from the CS7 census, and the census verdict stays falsified.

**What will be measured, once, after this is committed:**

1. The default run's three shares, at `CENSUS_MODEL_RUN` (120 cells, 40 h;
   before: 0.475 / 0.100 / 0.425) and at the 60-cell, 30 h test fixture
   (before: 0.467 / 0.100 / 0.433), `SIGN_CORRECTION_BEFORE`.
2. `scripts/gastrulation_census.py compare`, with the mappings, site
   filters, tolerance (0.05) and model run registered in `3c4e301`,
   unchanged. The result file is rewritten for the corrected model and
   carries the pre-correction shares beside the new ones.
3. `tests/test_gastrulation.py` as it stands. A test is changed only where
   the corrected sign changes its outcome, with the reason and the old value
   in a comment; a claim the corrected model no longer meets becomes a
   strict xfail, never a wider tolerance. The mesoderm strict xfail from
   `a398e9c` stays: if the corrected model makes it pass, that is reported
   as a change, not hidden. Predicted: endoderm may leave the unsourced
   0.20 ± 0.25 band (above 0.45), and the layer-order test fails if
   mesoderm vanishes.

**The run, once, after the registration was committed (`b1f3405`).**
Negatives first. The registered direction was half wrong: endoderm did
**not** widen, it narrowed, 0.425 to 0.367; mesoderm did narrow, and all
the way to nothing, 0.100 to 0.000; the cells went to ectoderm, 0.475 to
0.633. The reason is the runtime, not the paper: `genomeos/runtime/grn.py`
takes the *mean* over a gene's activators, so a second activator that is
low where NODAL is high (Tbxt there is about 0.3) halves SOX17's drive in
the endoderm band (Sox17 at x = 0.09 falls from 35.8 to 18.0), and caps it
at 0.8 of what NODAL alone gave. In the former mesoderm band (x ≈ 0.5)
Tbxt no longer wins: Sox2 ends at 35.9 against Tbxt 0.24, where Tbxt had
been 21.7. The prediction assumed activators add; the runtime averages
them, and the registration did not look.

| layer | before (census run) | after | human CS7 interval | outside by, before → after |
|---|---|---|---|---|
| ectoderm | 0.475 | 0.633 | 0.046 – 0.172 | 0.303 → 0.461 |
| mesoderm | 0.100 | 0.000 | 0.694 – 0.773 | 0.594 → 0.694 |
| endoderm | 0.425 | 0.367 | 0.116 – 0.213 | 0.212 → 0.154 |

Census verdict, same mappings, site filters, tolerance and model run as
`3c4e301`, written to
`data/results/gastrulation_census_comparison_sign_corrected.json` (one
departure from the registration: the corrected run gets its own file, so
the committed record of the `3c4e301` run is kept, not rewritten; its test
is now a strict xfail on the fresh-run check): **falsified** on all three
layers, as before; endoderm moved
closer and still misses by 0.154, three times the tolerance; mesoderm and
ectoderm moved further. The 60-cell, 30 h test fixture gives the same
shares (0.633 / 0.000 / 0.367, before 0.467 / 0.100 / 0.433).

Tests: `test_three_layers_in_order` now fails (0 of 60 cells are mesoderm,
before 6) and is a strict xfail with that reason; its end-order and
uncertainty checks moved to a test of their own and pass. The molecular
uncertainty label rose from "low" to "medium" because two rules moved from
inferred to experimental evidence, although their numbers are still
unsourced; the test now asserts "medium" and says so. The mesoderm strict
xfail from `a398e9c` still xfails, now at 0.0 against 0.35. Ectoderm 0.633
and endoderm 0.367 are both still within 0.25 of the unsourced 0.45 and
0.20. No number was changed after the run.

What it says: with a sourced sign and unsourced numbers, this network in
this runtime has no mesoderm band at the default gradient. Whether that is
the numbers, the mean-of-activators rule, or a missing interaction is not
decided here; each would be its own registered lane.

## The activator combination rule: a choice test (registered 2026-09-28)

`genomeos/runtime/grn.py` combines a gene's activators by their mean
(`ACTIVATOR_COMBINATION = "mean"`), and `runtime/located.py` repeats the
formula inline. lane-sign found the mean made a correctly signed second
activator lower SOX17; lane-bridge found 76,469 of 237,613 compiled gene-cell
pairs cannot be parameterised from single removals under it. This section
asks whether measured evidence, not taste, picks a rule among four, each over
the terms t_i = s_i H(x_i) and all equal for one activator: **mean**
Σt/n; **sum_capped** min(1, Σt); **max**; **or** 1 − Π(1 − t_i).
`grn.combine_activators` computes each; the runtime still calls it with the
mean, so no value has moved.

**Census, computed before the registration** (`scripts/activator_combination.py`,
`data/results/activator_combination_census.json`; the per-rule runs use the
scratch instrument `scripts/activator_combination_patch.py`, which swaps the
rule in both runtimes and records every evaluation of a gene with two or more
activators):

- Every hand-written program (41 files under data/demo, data/organisms and
  genomeos/std, every context their rules name): **one gene, SOX17 in
  data/demo/gastrulation.bio (Nodal and Tbxt), has two activators.** The
  C. elegans, haematopoiesis, erythrocyte, body and std programs have none;
  data/organisms/human/noncoding_chr21.bio has 120 such genes, none with a
  `max_rate`, so none is simulated.
- The 24 compiled chromosomes: 48,405 gene-cell pairs with two or more
  activating rules (148,765 rules on them) of 183,304 with any; **0 simulable**
  (every compiled gene is a stub, as lane-bridge's audit found).
- The full test suite under the mean with the instrument: of about 2,050
  tests, **four** evaluate a gene with two activators, all SOX17 or the
  bridge's own pin: `test_gastrulation.py::test_census_comparison_as_run_once`,
  `::test_census_comparison_after_the_sign_correction`,
  `::test_three_layers_in_order` and
  `test_grn_bridge.py::test_under_the_mean_rule_adding_an_activator_can_lower_expression`.
  `bio test` over the three program trees: 220 of 220 checks pass under all
  four rules, and none evaluates a gene with two activators.
- The gastrulation model under each rule (census run 120 cells, 40 h; test
  fixture 60 cells, 30 h; no number changed):

| rule | census run ecto / meso / endo | fixture | CS7 verdict (3c4e301 mappings) |
|---|---|---|---|
| mean | 0.633 / 0.000 / 0.367 | 0.633 / 0.000 / 0.367 | falsified |
| sum_capped | 0.558 / 0.000 / 0.442 | 0.550 / 0.000 / 0.450 | falsified |
| max | 0.575 / 0.000 / 0.425 | 0.567 / 0.000 / 0.433 | falsified |
| or | 0.558 / 0.000 / 0.442 | 0.567 / 0.000 / 0.433 | falsified |

- The full suite under each alternative (the tests of this lane excluded):
  **sum_capped** fails three,
  `test_gastrulation.py::test_census_comparison_after_the_sign_correction`
  (the committed shares), `::test_ectoderm_and_endoderm_within_unsourced_expectation`
  (fixture endoderm 0.450, at the edge of 0.20 ± 0.25) and the bridge's mean
  pin; **max** and **or** fail two, the census pin and the bridge pin. No
  strict xfail passes under any rule (the mesoderm xfails stay at 0.0). These
  are every pin that would move.

  **Negative for the sign lane's reading:** the mesoderm band is absent under
  every rule, so the mean is not what removed it; it only sets how much
  endoderm the second activator costs (0.367 against 0.42–0.45).
