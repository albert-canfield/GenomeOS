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

**What evidence exists, and what it can discriminate.** Read before this
registration was written, so the outcome of the published cases below was
known to the lane when it wrote the pass rule:

1. *Published reporter pairs* (Bothma JP, Garcia HG, Ng S, Perry MW, Gregor T,
   Levine M 2015, eLife 4:e07956, PMC4532966): one transgene with the primary
   enhancer, the shadow enhancer, or both, live-imaged. With the construct's
   enhancers as the gene's whole activator set, each rule predicts the pair P
   from the singles S1, S2 without any number: mean P = (S1 + S2)/2, max
   P = max, sum_capped and or max ≤ P ≤ S1 + S2. The paper's own words:
   *knirps* "early ... super-additively" (P > S1 + S2) and later "a simple
   additive manner"; *hunchback* "sub-additively in anterior regions",
   "additively" in central ones ("the wild-type transgene produces
   significantly higher levels of expression than either of the transgenes
   driven by a single enhancer"); *snail* "the wild-type transgene displays
   significantly lower levels of expression than the mutant transgene
   containing only the shadow enhancer". **Can discriminate:** additive pairs
   contradict mean and max; the *snail* pair below its stronger single
   contradicts every monotone rule (sum_capped, max, or) and is the one
   behaviour only the mean produces; super-additivity contradicts all four.
2. *Published CRISPRi doubles* (Zhou J, Guruvayurappan K, Toneyan S et al.,
   bioRxiv 2023.04.26.538501, Cell Genomics 2024; Gasperini 2019 K562
   single-cell data): in 264 high-confidence enhancer pairs "the
   multiplicative model provided a better fit" than the additive one, and no
   interaction term was significant. **Cannot discriminate** these four: it
   compares two GLM links, and each rule's double-removal prediction depends
   on the gene's other activators and on saturation.
3. *The project's CRISPRi* (EPCrisprBenchmark training split, K562,
   `measured.load_crispri(split="training")` through `development_only`):
   single removals only, no doubles. **Can discriminate max only:** a
   gene-cell with two significant decreases is impossible under max
   (removing an activator below the maximum changes nothing). Under the mean
   the single-removal changes of a gene's activators sum to exactly zero
   (ΔA_i = (A − t_i)/(n − 1)), so decreases need balancing increases, but
   untested or undetected activators can carry them; the ratio is reported,
   not judged. Sum_capped and or allow any pattern of decreases.
4. *Textbook fixtures in tests*: the census above shows no C. elegans,
   haematopoiesis or cell-cycle program has a gene with two activators in
   this runtime (the Faure cell cycle runs in `runtime/boolean.py`, whose
   logic is written per gene), so **none can discriminate**. The one case,
   gastrulation, rests on unsourced strengths and thresholds: whether Tbxt
   removal lowers SOX17 (Lolas 2014 Fig. 3A) depends on them under the mean.

**The test** (`grn.COMBINATION_TEST`, run once by
`scripts/activator_combination.py choose`): a rule is contradicted by a
case only when no choice of its numbers can produce the observation. Cases:
the five Bothma behaviours above and the CRISPRi two-decrease count (max is
contradicted if at least 10 training gene-cells have two or more significant
decreases). A case that contradicts all four rules (knirps early) is
reported as outside the family and discriminates nothing. **Pass rule:**
switch from the mean to R only if the mean is contradicted by at least one
judged case and R by none; if no alternative survives, or several do, the
mean is retained and the choice recorded as undetermined. **No strength,
threshold, Hill coefficient or rate is re-tuned for any rule,** and the
census comparison is reported per rule but is not a criterion: picking the
rule that fits the census best would be a tuning with one discrete knob.

**Expected, stated before the run:** from the published cases alone the mean
is contradicted twice (knirps late, hunchback central) and every
alternative at least once (*snail*), so the registered outcome is "mean
retained, undetermined"; the CRISPRi count can only add a contradiction of
max, so it cannot change that. The run is therefore mainly the CRISPRi count
and its record; the lane says so rather than presenting a foregone verdict
as a discovery. The CRISPRi outcomes had not been counted when this was
written (only the table's schema was read); `tests/test_activator_combination.py`
pins the constants and the identities each rule's prediction rests on.

**The run, once, after the registration was committed (`7f3e597`).**
Written to `data/results/activator_combination_choice.json` (the computation
was executed again, unchanged, only to write that file beside the census
record rather than over it; it is deterministic and every figure was the
same). Negatives
first. **No rule survives, so the mean is retained and the choice
is undetermined**, as expected. Contradicted by: mean, knirps late and
hunchback central (additive pairs, where the mean gives their average);
sum_capped and or, snail (a pair below its stronger single, which no
monotone rule can give); max, all three of those and the CRISPRi count.
Knirps early (super-additive) was set aside as contradicting all four.

The CRISPRi training split (10,356 element-gene rows in K562, no conflicting
repeats, 2,113 gene-cells): 239 gene-cells have one significant decrease and
**81 have two or more** (47 two, 19 three, 7 four, 8 five or more; MYC nine),
so max is contradicted by 81 against the registered 10. Descriptive, not
judged: 14 of those 81 carry any significant increase, and the summed
increase magnitude is **0.075 of the summed decrease** on them, where the
mean predicts 1 if the screen had tested every activator of those genes. It
did not test every element, so this says only that the mean needs untested
or undetected activators whose removal raises expression by as much as the
measured decreases lower it.

What it means. The four candidates form one family (a fixed function of the
terms) and the measured pairs change behaviour with enhancer strength:
weak pairs add, strong pairs fall short of adding and, in snail, interfere,
which Bothma et al. attribute to competition for the promoter. No member of
the family reproduces both, so choosing one would be taste. The mean stays
the named assumption, and lane-bridge's finding stands unchanged: a gene
with two or more mechanisms in one context (76,469 of 237,613 compiled
gene-cell pairs) cannot be parameterised from single-removal observations,
because under the mean each removal's effect depends on every other
activator's term (ΔA_i = (A − t_i)/(n − 1)), and no double-removal data exist
in the project to separate them. A form with a strength-dependent
interaction (promoter competition) would be its own registration, and it
would need numbers the project does not have. Nothing moved: no rule
changed, no pin changed, the census verdict is the `5cbce26` one
(falsified). The per-rule gastrulation table above is the answer to "does
the rule explain the missing mesoderm": it does not.

## Declared rates for the compiled programs (registered 2026-09-28)

The R3 audit's starkest line is that **0 of 440,589 compiled rules can be
simulated as compiled**. 76,469 of 237,613 gene-cell pairs carry more than
one mechanism and are not identifiable from single removals under the mean
rule; that is settled and this section does not touch it. The other
**161,144 lack only declared rates**: a basal rate and a maximum
transcription rate per gene, which `runtime/grn.py` defaulted to zero and now
diagnoses. This section is the registration for supplying them. It was
written and committed before the rate table was built and before any count
was read.

**The sources.** `attribution/bridge.py`'s `MEASURED_RATE_SOURCES` holds them
with sizes and digests. In short: Schwanhausser et al. 2011 (Nature 473:337,
Supplementary Table 3, the file replaced in 2013 for the corrigendum, Nature
495:126) gives 4,338 transcription rates in molecules per cell per hour,
4,658 mRNA half-lives and 5,028 protein half-lives in hours, and 4,309 mRNA
copy numbers, all in **mouse NIH 3T3 fibroblasts**. Schofield et al. 2018
(Nat Methods 15:221, Supplementary Table 2) gives transcript half-lives in
**human K562** and **mouse embryonic fibroblasts** by one method, which is
the only matched cross-species pair here and so carries the species
falsifier. MGI's mouse-human homology report carries a mouse symbol to a
human one. Both supplementary files are Springer Nature material with no open
licence stated: they are cached in the git-ignored `data/cache/rates/`, read,
and never redistributed from this repository. The project's existing
`bio.std.human_turnover` (Sender & Milo 2021) is human cell **lifespans in
days** and is listed only so that nobody reaches for it: a cell's replacement
clock is not a molecular half-life and never becomes a decay constant.

**The unit mapping.** A measured half-life is a degradation constant, not a
maximum transcription rate: mRNA half-life becomes `delta_m = ln2 / t_half`
in 1/h and protein half-life becomes `delta_p`, and neither sets any level on
its own. A measured translation rate constant becomes `k_tl` and is not a
transcription rate. A measured mRNA copy number becomes **nothing**: it is
held for the acceptance check and never supplied as a parameter. And the
measured transcription rate becomes `T`, the gene's **total** transcription
in the one unperturbed cell it was measured in, `T = basal + max_rate * A *
R`. It does **not** become `max_rate`. `max_rate` is the ceiling of the
regulated term, reached only when every activator saturates, and no
steady-state measurement in any of these sources observes that state.

**The split, with no invented fraction.** For the case that matters, a gene
carrying one mechanism in this cell (a gene with more is already
`not_identifiable`), the measured `T` and the observed removal fold `RHO`
determine the basal rate exactly. Activator: intact is `b + V*s*h = T`,
removed is `b`, so `b = RHO * T` and `V*s*h = (1 - RHO) * T`. Inhibitor:
removed is `b + V = RHO * T`, so `V*s*h = (RHO - 1) * T` and `b = RHO*T - V`.
What is **identified** is `C = |1 - RHO| * T`, the element's own contribution
to transcription in molecules per cell per hour. `V` and `s` are **not**
separately identified; only their product is. `DECLARED_STRENGTH = 1` fixes
the split by convention and by nothing else, and it does not license the
claim that the element saturates its gene: every `s` in (0, 1] with
`V = C / (s*h)` reproduces this observation identically, and the simulation's
answer to any other perturbation depends on `V` and `s` separately.

**What a measured rate does and does not buy.** A fold is dimensionless and
the split is linear in `T`, so the fitted strength does not depend on `T` at
all. The measurement buys the absolute scale and the basal fraction, in
molecules per cell per hour; it buys nothing about the regulation. That is
why a gene with no measurement of its own is either an explicit
`gene_rate_unmeasured` marker with no number in its place, or, in a
separately counted lower tier, given the genome median as an openly borrowed
constant, which fixes no absolute number and makes the pair simulable in
relative units only. The borrowed tier is never added to the measured count.

**The species limit, stated not hidden.** Every transcription rate here is
mouse. A human gene simulated with one is borrowing a mouse fibroblast
constant, the parameterised genes carry `species: mouse`, and the audit
repeats it.

**Disclosed before the fact.** While surveying the source it was already
computed that Schwanhausser's `vsr` is not the identity `N * ln2 / t_half`
over its own columns: the median relative difference is 0.236 across the
4,309 genes carrying all three, because `vsr` comes from their ODE fit to
time courses and not from that ratio. That is a property of the source,
disclosed here, and is not presented as a passed test. The acceptance check
below therefore tests the runtime against `T / delta_m`, which the model does
determine.

**Acceptance** (`RATE_ACCEPTANCE`, fixed before the build, as lane-bridge's
fixture was). A one-gene program whose gene carries a measured `T` and a
measured mRNA half-life, and whose one measured regulator is observed at
`RHO`, parameterised through the bridge and integrated to steady state with
the element clamped at 1 and at 0, reproduces `RHO` within 1%; the same run's
intact steady-state mRNA equals `T / delta_m` within 1%, so the level reads
in molecules per cell and not in a.u.; a gene absent from the table is
reported `gene_rate_unmeasured` with no number put in its place, and a split
that leaves a non-positive basal rate is reported `rate_split_infeasible`.

One consequence to state plainly: because the declared split puts the
refitted strength exactly at 1, the bridge's acceptance window (0, 1] gains a
`STRENGTH_TOLERANCE` of 1e-9 at its upper end. That is a floating-point
allowance and nothing else; a strength above it is still
`response_out_of_range`.

**Falsifiers** (`RATE_FALSIFIER`, none of them computed when this was
written). *Species transfer*: over MGI one-to-one homology classes the
Spearman correlation between Schofield's mouse fibroblast and human K562
half-lives is computed once; if it is below 0.5 the mouse constants are
declared non-transferable and the registration states that the absolute human
numbers are supported by nothing measured in a human cell. *Inhibitor bound*:
under the declared split an inhibitory mechanism needs `RHO < 2` to leave a
positive basal rate; if more than half of the inhibitory mechanisms that
reach the split fail it, the binding limit is the declared Hill assumption
(`K = 1`, `n = 2`, so `H = 0.5` at an intact element) and not the data, and
that is the finding. *Invariant*: supplying rates must not move
`not_identifiable` (76,469) or any reason other than
`gene_parameter_missing`, because a rate is not an observation; a change
there is a defect.

**Expected count, with a direction** (`RATE_EXPECTED`). Of the 161,144
mechanisms now `gene_parameter_missing`, the number that becomes simulable on
measured rates will be strictly greater than 0 and strictly less than
161,144, and is predicted to fall between **15,000 and 60,000**: the source
carries 4,338 measured rates against roughly 20,000 protein-coding genes, and
those genes are the abundant ones, which may be over- or under-represented
among the targets of compiled non-coding elements. Simulable gene-cell pairs
rise from 0 by the same order. A count outside 15,000-60,000 is reported as a
miss.

## What the rates bought: 22,576 of 161,144, in molecules and not in a.u. (2026-09-28)

The registration above was committed in `ff2062a`; these are its counts, from
`scripts/build_gene_rates.py` and a re-run of `scripts/bridge_audit.py` over
all 24 compiled programs (`data/results/gene_rates.json`,
`data/results/bridge_audit.json`).

**The table.** Schwanhausser's 5,028 protein groups carry 4,338 transcription
rates. Requiring that every mouse symbol on a row point at the same human
gene through a 1:1 MGI homology class, and dropping the 19 human symbols that
would have received two different rows, leaves **3,603 human genes** with a
measured total transcription rate, a median of 1.76 molecules per cell per
hour and a range from 0.08 to 610.52. The strict alternative (rows naming
exactly one mouse symbol) would have left 1,945 rows and 1,424 genes; the cost
of that choice is recorded rather than hidden.

**The counts, before and after.** Nothing else changed: 440,589 active
regulatory rules, 440,550 mechanisms, 237,613 gene-cell pairs, 6,938 contexts.

| | before (no rates) | measured | borrowed median |
| --- | --- | --- | --- |
| simulable mechanisms | 0 | **22,576** | 160,392 |
| simulable gene-cell pairs (of 237,613) | 0 | **22,576** | 160,392 |
| compiled rules simulated (of 440,589) | 0 | **22,580** | 160,405 |
| `gene_parameter_missing` | 161,144 | 0 | 0 |
| `gene_rate_unmeasured` | - | 138,543 | 0 |
| `rate_split_infeasible` | - | 25 | 752 |
| `not_identifiable` | 76,469 | 76,469 | 76,469 |

**Against the registration.** The predicted interval was 15,000 to 60,000 and
the measured count is 22,576, inside it; the direction held, strictly above 0
and strictly below 161,144. The invariant held exactly: `not_identifiable`
stays at 76,469 in all three tiers, so no rate was mistaken for an
observation. The inhibitor bound was not the binding limit: 25 of the 5,245
inhibitory mechanisms that reached the split failed it, 0.48%, against the
registered "more than half" that would have convicted the declared Hill
assumptions. The species falsifier passed as registered, at a Spearman of
0.617 between Schofield's mouse fibroblast and human K562 half-lives over
2,133 one-to-one homologues, above the 0.5 threshold.

**What stays unresolved, and why.** 76,469 gene-cell pairs still carry more
than one mechanism and no single-removal observation identifies a strength
under the mean rule: unchanged, and no rate could have changed it. 138,543
mechanisms sit on a gene the source never measured, and they are marked, not
filled. 25 are inhibitions whose removal more than doubles the gene, which
the declared split cannot give a positive basal rate. That is the whole of
the 161,144.

**What the number does not mean.** Three things, plainly. First, the fitted
strength is unchanged by the measured rate, so the 22,576 is a gain in units
and in the basal fraction, not in knowledge of the regulation: what those
runs now report is molecules per cell per hour instead of arbitrary units.
Second, every one of those genes is a human gene in a human cell type running
on a **mouse fibroblast** constant, and in most of those contexts (K562,
placenta, glutamatergic neuron) it borrows across cell type as well as across
species. The matched Schofield measurement puts a number on the second
borrowing: the median transcript half-life is 2.73 h in mouse fibroblasts and
1.56 h in human K562, so the same gene's steady-state level moves by about a
factor of 1.75 between the two, a factor the rank correlation of 0.617 does
not see and the transferability verdict does not cover. Third, the borrowed
median tier makes 160,392 pairs runnable but fixes no absolute number at all:
those runs are in relative units and are never to be reported as measured.

**Where the acceptance was checked.** `tests/test_gene_rates.py`, eleven
tests: the split returns the measured total and the observed fold, the fold is
reproduced end to end within 1%, the intact steady state lands on `T / delta_m`
within 1% so the level reads in molecules per cell, doubling the measured
half-life doubles the level and moves no fold, an unmeasured gene gets a
marker and no number, an inhibition beyond the bound is reported and not
clipped, and a borrowed rate changes no strength.

## Registration: a second rate source, and a human one if there is one (2026-09-28)

Written before any source was applied, alongside the constants in
`genomeos/attribution/bridge.py`
(`HUMAN_RATE_SOURCE_SURVEY`, `HUMAN_RATE_CONSTRUCTION_NOT_TAKEN`,
`HUMAN_RATE_SURVEY_VERDICT`, `HUMAN_RATE_SURVEY_FALSIFIER`,
`RATE_ROW_PROVENANCE`, `HUMAN_RATE_SURVEY_EXPECTED`).

**The question.** The rates above leave two limits. 138,543 gene-cell pairs
sit on a gene nobody measured and carry a marker with no number. Every rate
that does exist is mouse NIH 3T3 fibroblast, so all 22,576 human pairs borrow
across species and cell type. Both have the same candidate fix: a second
source of **absolute** transcription rates, ideally human.

**What counts as a source.** A per-gene transcription or synthesis rate in
molecules per cell per unit time, or a quantity with a stated conversion to
it, with the per-cell calibration named: a cell count, an RNA mass per cell,
or a spike-in weight per cell. **A source that reports only half-lives cannot
supply a synthesis rate**, and recording that is the result, not a failure.
No rate is manufactured by pairing a half-life from one study with a copy
number from another unless this registration names that as an explicit,
labelled construction with its assumptions - and it does not.

**What a rate would NOT license,** carried over unchanged from the rates
registration: a measured total rate `T` becomes `basal_rate + max_rate * A * R`
at the observed state and **never becomes `max_rate`**; a half-life becomes
`delta_m` and sets no level; a copy number becomes nothing and is held only
for the acceptance check.

**How human and mouse would be kept apart.** Per row, not per note. Every row
of `data/results/gene_rates.json` carries `source` and `species_cell`, and the
audit counts `pairs_on_a_human_measured_rate` and
`pairs_on_a_borrowed_species_rate` for every tier. Two sources are never
averaged for one gene and never pooled into one median; a pair standing on a
mouse constant is never reported as a human result.

**The expected change, with a direction.** If a human absolute source is
found, simulable pairs on a measured rate rise above 22,576 and
`gene_rate_unmeasured` falls below 138,543. If none is found, **no count
changes at all**: measured stays 22,576, borrowed 160,392,
`gene_rate_unmeasured` 138,543, `not_identifiable` 76,469, and pairs on a
human rate stay 0. A survey is not a rate.

**The invariant.** The 76,469 not-identifiable pairs do not move, in any tier.
No existing measured figure moves: the 3,603 genes, the median of 1.76 and the
range 0.08 to 610.52 are unchanged, because nothing was added to the table.

**The falsifier.** A dataset distributing, per gene, for a human cell type, a
rate in molecules per cell per unit time with its per-cell calibration named.
`tests/test_human_rate_survey.py` holds it as a test that fails the day a
candidate is recorded as both human and absolute.

## The finding: there is no human absolute transcription rate to borrow (2026-09-28)

Ten candidates were read in the primary source. Three classes were searched:
metabolic labelling that reports synthesis and not only decay (TT-seq,
TimeLapse-seq, SLAM-seq, Bru-seq, single-cell new-RNA), absolute
transcriptome or proteome quantification paired with turnover, and any human
equivalent of Schwanhausser. Europe PMC full-text queries and the GEO
deposits were used; two articles (Schwalb 2016 in *Science*, Ietswaart 2024 in
*Molecular Cell*) were not readable by this lane and their records say exactly
what was read instead and claim nothing beyond it.

| source | species | cell | quantity | units | absolute | licence |
| --- | --- | --- | --- | --- | --- | --- |
| Schwalb 2016 (TT-seq) | human | K562 | synthesis rates, half-lives | not established here | not distributed per gene | subscription; GSE75792 public |
| Michel 2017 (TT-seq) | human | Jurkat | synthesis rate mu, decay lambda, 22,141 TUs | library units | no | CC BY 4.0 |
| Wachutka 2019 (TT-seq) | human | K562 | bond synthesis and cleavage | minutes | no | CC BY |
| Ietswaart 2024 (subcellular TimeLapse-seq) | human + mouse | K562 + NIH 3T3 | release, export, degradation constants | per hour | not established | subscription |
| Shao 2022 (TT-seq) | **mouse** | ESC | synthesis rate = labelled rate x copies per cell | **cell^-1 min^-1** | **yes** | CC BY 4.0 |
| Hausser 2019 | human | HeLa | transcription rate beta_m, thousands of genes | mRNA per cell per hour | **constructed** | CC BY 4.0 |
| Liu 2023 (SLAM-Drop-seq) | human | HEK293 | transcription, splicing, degradation, 399 genes | CPM per hour | no | CC BY 4.0 |
| Ramskold 2024 (NASC-seq2) | human + mouse | K562, fibroblasts | burst on/off rates, burst size | detected molecules | no | CC BY 4.0 |
| Schofield 2018 (cached) | human + mouse | K562 + MEF | transcript half-lives | hours | n/a | subscription, cached |
| Schwanhausser 2011 (in use) | **mouse** | NIH 3T3 | 4,338 transcription rates | molecules/(cell*h) | yes | subscription, cached |

Three of those rows decide it.

**The only absolute per-cell rate found is mouse.** Shao et al. 2022 state it
exactly: "RNA synthesis rate (cell^-1 min^-1, or copy/min per cell) was
calculated by multiplying labeled rate and transcript copy number", with the
copy number per cell coming from a spike-in weight model. So the method works
and the units are real - in mouse embryonic stem cells. Adopting it would swap
one borrowed species for the same species in another cell type, and its
per-gene rates are not distributed anyway: Table EV1 is elongation velocities,
and the rates would have to be re-derived from GSE168378.

**The one published set of human transcription rates is anchored on a mouse
constant.** Hausser et al. 2019 report `beta_m` for thousands of human HeLa
genes in mRNA per cell per hour, which is the right form. Reading the methods
settles it: the per-gene number is an abundance (from Eichhorn 2014 mRNA-seq
and ribosome profiling) rescaled by a total of about 225,000 mRNAs per HeLa
cell and by **one** global decay rate, 0.06 /h, from a median HEK293 half-life
of 11.4 h. And that 225,000 is not measured: they write "We could not find
direct measurements of the number of mRNAs per HeLa cell N_m", and obtain it by
taking the 180,000 mRNAs of a **mouse 3T3 cell** and scaling by a cell-volume
ratio 2500/2000. Adopting these rates would launder the very constant this
lane set out to escape into a human-labelled number, and per gene they are
proportional to expression and carry no kinetics of their own.

**The construction was available and was not made.** `T = N * ln2 / t_half`
would turn a human half-life into a human absolute rate. It needs a measured
per-gene human mRNA copy number; none was found genome-wide, and the strongest
evidence that none exists is the sentence above, from a group that needed it.
There is also a measured reason to distrust the identity: over the 4,309
Schwanhausser genes carrying a rate, a copy number and a half-life together,
the median relative disagreement between the measured `vsr` and
`N * ln2 / t_half` is **0.236**. So the construction has no human input and a
24% disagreement with the one measurement it can be checked against.
`tests/test_human_rate_survey.py` checks in the data that no rate in the table
is that construction.

**The counts, before and after: identical, as registered.** The re-run of
`scripts/bridge_audit.py` over all 24 compiled programs moved nothing.

| | before | after |
| --- | --- | --- |
| simulable gene-cell pairs, measured tier | 22,576 | 22,576 |
| simulable gene-cell pairs, borrowed tier | 160,392 | 160,392 |
| `gene_rate_unmeasured` | 138,543 | 138,543 |
| `not_identifiable` | 76,469 | 76,469 |
| pairs on a **human** measured rate | (not counted) | **0** |
| pairs on a borrowed-species rate | (not counted) | **22,576** |

The two new rows are the whole of what changed: the borrowing is now counted
per pair rather than stated once in a header, so no later result can quote the
22,576 without the species travelling with it.

**What stays unmeasured.** All of it. 138,543 pairs sit on a gene with no
measured transcription rate in any species, and this lane found no source that
reduces that number for a human cell. 76,469 remain not identifiable, which no
rate could ever change. The 22,576 that do run keep their mouse fibroblast
constant, and the honest description of the human rate table is that it is
empty.

**Where the next lane should look.** Ietswaart et al. 2024 measures human K562
and mouse NIH 3T3 - the very cell type the R3 rates come from - by one method.
It is named in the falsifier for a reason worth more than a rate: it would
replace the present species test (mouse MEF against human K562, which confounds
species with cell type) with a matched one. Shao et al. 2022's GSE168378 shows
the spike-in-per-cell calibration a human TT-seq dataset would need.

## Registration: runnable simulation kept apart from validated human kinetics (item 12 S5, 2026-09-28)

Written and committed before any label was built and before any pair was run
across a range, beside the constants in `genomeos/attribution/bridge.py`
(`SIMULATION_STATUS`, `TRANSFER_LABEL`, `S5_RANGES`, `S5_COMBINATION_RANGE`,
`S5_NOT_PROPAGATED`, `S5_PREDICTIONS`, `S5_SAMPLE`, `S5_EXPECTED`,
`S5_FALSIFIERS`, `S5_RUNTIME_GAP`, `S5_HUMAN_CONTEXT_KINDS`). The range widths
were read from the rate sources alone; no compiled program had been run and no
simulated level or classification existed.

**What the review found and what this does with it.** 22,580 of 440,589
compiled rules (5.1%) are simulable, every rate behind them is mouse NIH 3T3,
and that supports **assumption-dependent simulation, not measured human
dynamics**. That phrase is now a constant, and every output built on the bridge
carries it. A second category, validated human kinetics, is named and holds
nothing: a prediction checked against a human measurement that was not an input
to it, in the cell it is about. **Stable across the ranges below is not
validated**, and no count from this work is reported under that name.

**Negatives first.** Four things were already true before anything was run.
The runtime holds one mRNA half-life per module, and the bridge sets none, so a
compiled program run with measured rates reports a mouse transcription rate
divided by an a.u. default decay: mixed units, which the label will name. The
combination rule applies to no simulable pair, because the bridge fits only a
gene with one mechanism in its cell, and with one activator mean, capped sum,
max and OR are the same function; it is run anyway and every difference counted.
The Hill threshold and coefficient and `DECLARED_STRENGTH` are declared
assumptions, not transferred rates, and no measured spread for them can be
read, so they are not propagated: at the two states the bridge fits (element
present, element removed) no prediction depends on them, and the predictions
that do (a partial dose, whether an inhibitory pair is feasible at all) are
named undetermined. And the one matched human-against-mouse half-life figure
quoted so far, 1.56 h against 2.73 h, compares Schofield's human K562 with
Schofield's mouse MEF; the half-lives the rate table actually carries are
Schwanhausser's NIH 3T3 ones, median 9.96 h, so the distance between what the
model would use and the human measurement is larger than that figure says, and
most of it is method rather than species.

**Step 1, the label.** On the output, per gene and per run, never in a header.
`parameterize(rate_provenance=...)` writes the rate's source, the species and
cell it was measured in, and its tier on each gene whose rate it used, and
returns the same per gene; the runtime reads them back onto the trajectory with
the cell the run is in, and lists every parameter the run took from its a.u.
defaults. A rate whose species the caller did not state is labelled
`unstated`, never assumed mouse or human. The bridge audit counts labelled
pairs per tier and they must equal simulable pairs.

**Step 2, the ranges** (log2; percentiles by linear interpolation). Two
quantities are transferred into a simulable pair: `T`, the gene's total
transcription rate, and `t_half`, its mRNA half-life.

| range | read from | width |
| --- | --- | --- |
| R1 measurement | Schwanhausser's experiment and replicate columns: 95th percentile of \|log2(experiment/replicate)\| over every row with both (3,604 for `T`, 4,658 for `t_half`) | `T` x 2^±1.1186, `t_half` x 2^±1.1210 |
| R2 transfer to a human cell | per gene, log2 of Schofield's human K562 half-life over the Schwanhausser NIH 3T3 half-life the table carries, 2,262 genes | `t_half` x 2^[-3.6966, -0.6457] (0.077 to 0.639; median 0.229) |
| R2 for `T` | no matched human rate exists, so two readings of that one measurement bracket it: the mouse rate carries over, or the mouse copy number does (then `T` scales inversely with the half-life) | corners r in {0.077, 0.639} x reading in {rate, copy number} |
| R3 no transfer (sensitivity) | the across-gene central 90%: Schwanhausser `T` over 3,603 table genes, Schofield K562 `t_half` over 5,419 | `T` 0.39 to 10.02, `t_half` 0.40 to 7.03 h, the same box for every pair |

The primary range is the nominal mouse values together with every corner of R1
and R2. R3 asks what survives if the gene's own mouse numbers say nothing about
the human gene, and is reported beside the primary, never as it. K562 is the
only human cell whose half-lives are cached, so its spread stands in for every
human context and is a lower bound on the transfer to any other; cell-type
differences within one species at one method cannot be read from anything
cached.

**Step 3, the predictions and what stable means** (`S5_PREDICTIONS`). A
prediction is stable when it gives the same answer at every point of the
range. Direction and fold are the observation reproduced, and are listed so
that their stability is never mistaken for a result. Level and effect size in
molecules per cell, and response time (hours to cover half the way to the
removed level), are stable if their spread over the points is at most two-fold,
one log2 unit, the unit the observations are in. On/off (at least one molecule
per cell, the smallest level with a physical reading) and switch (removal takes
an expressed gene below one molecule) are stable if their answer never changes.
Rank is the order of two pairs in one cell context, each free anywhere in its
own range, reported as the share of within-context comparisons whose order
cannot change.

**The sample.** 1,000 measured-tier simulable pairs, a simple random sample
without replacement by `random.Random(20260928).sample` from the pairs sorted by
chromosome, cell, gene and element. Each is a one-gene program parameterised
through the bridge with its label and integrated by `NetworkRuntime` at every
point: to steady state with the element present, then from that state with it
removed. The same predictions for all 22,576 pairs come from the model's closed
form and are reported only if they agree with the runtime on every sampled
pair.

**Expected, with a direction.** The first four follow from the widths by
arithmetic and are consequences, not tests. Level and effect size are stable
for **no** pair: a multiplicative range of one width gives every pair the same
spread (R1 22.3-fold, R2 13.0-fold, primary 61.2-fold), so whether a level is
stable is decided by the width and not by the pair. Response time is stable for
no pair (4.7, 8.3 and 28.2-fold). Direction and fold are stable for all 22,576.
The combination rule changes no sampled pair. Then the tests: on/off is stable
for more than half of the pairs and fewer than all (stable on needs a nominal
level of at least 12.97 molecules per cell, stable off one below 0.212); switch
is stable for no more pairs than on/off; fewer than half of the within-context
rank comparisons are stable; under R3, on/off, switch and rank are stable for no
pair and only direction and fold survive.

**Falsifiers** (`S5_FALSIFIERS`). Runtime against formula on the sample: steady
states within 1%, response times within 2%, classifications identical except
within 1% of a threshold, or the census is withdrawn and only the sample with
its Wilson interval is reported. Label: every sampled trajectory names its
source, species, cell and context, and labelled pairs equal simulable pairs in
every tier of the audit, or step 1 is not met. Combination: any difference is a
bridge defect. Invariant: the audit re-run moves no existing count.

**Step 4, the one human context.** Chosen by a rule fixed now: among contexts
holding simulable pairs, the one where the project holds the most of seven
measurement kinds in that very cell (CRISPRi of elements, lentiMPRA, DNase,
human mRNA half-lives, an absolute human transcription rate, absolute mRNA
copies per cell, a time course after a perturbation), ties to the context with
more simulable pairs. It is named with what a validated human-kinetics test
there would need that is missing. Nothing is built.

**Amendment, before the registered run** (`S5_AMENDMENT`). Written after the
label was built and the propagation code was dry-run on chr21 (229 pairs, not
the registered run; none of its numbers is reported as a result). Two
corrections, neither of which moves a range, the sample, a threshold or a
criterion. First, the switch is not monotone in the level: removal switches an
activated gene off only when its intact level lies between one molecule and
1/RHO, a band, so an answer read only at the corners of a range can miss the
band's interior; it is evaluated over the whole level interval the corners span.
Every other prediction is monotone in the level or the half-life, so its corners
are its extremes. Second, the registered consequence "switch stable for no more
pairs than on/off" is false as written, and so is R3's "switch stable for no
pair": removal raises an inhibited gene, so an inhibitory pair can never be
switched off and its switch answer is stable whatever its on/off answer. Both
hold among activating pairs only. On/off and switch are therefore also reported
among activating pairs, and the two sentences as written are reported as wrong
by arithmetic, not by data.

## What survives the ranges: a sign and "expressed", never a level (item 12 S5, 2026-09-29)

The registered run of `scripts/s5_rate_ranges.py` under `0ea0a90` and its
pre-run amendment (`83013e6`), `data/results/s5_rate_ranges.json`, with
`data/results/bridge_audit.json` re-run for the label counts. The gate passed:
every range width re-read from the sources equals the registered number to four
decimals. 0 AlphaGenome requests.

**Everything below is assumption-dependent simulation, not measured human
dynamics. Stable across these ranges is not validated.** Nothing here was
compared with a human measurement it did not take as an input.

**The label reached every output.** All 22,576 measured-tier simulable pairs
(22,580 compiled rules) carry their label: source `schwanhausser2011`, measured
in *Mus musculus* NIH 3T3 fibroblasts, used in the pair's own human context;
0 unlabelled, 0 unstated, 0 on a human rate. In the audit, labelled pairs equal
simulable pairs in every tier (22,576 measured; 160,392 in the borrowed tier, of
which 137,816 are labelled `borrowed_median`), and no existing count moved. Every
one of the 26,000 trajectories the sample ran carried its transfer entry with
the right source, species, tier and context.

**The runtime and the formula agree.** On the 1,000 sampled pairs, at all 13
points each, the runtime's steady states were within 1% and its half-way times
within 2% of the closed form at 13,000 of 13,000 points, and all 1,000
classifications were identical, so the census over all 22,576 is reported. The
four activator rules gave identical runs on all 1,000 pairs: the combination rule
applies to no simulable pair, because the simulable tier is exactly the set of
genes with one mechanism.

**What is stable, over the primary range** (the nominal mouse values and every
corner of R1 and R2), census of 22,576 pairs:

| prediction | stable pairs | share | what the stable answer is |
| --- | --- | --- | --- |
| direction of the removal effect | 22,576 | 100% | the observation itself, reproduced |
| fold on removal | 22,576 | 100% | the observation itself, reproduced |
| level, molecules per cell | 0 | 0% | - |
| effect size, molecules per cell | 0 | 0% | - |
| response time after removal | 0 | 0% | - |
| on/off (at least one molecule per cell) | 16,101 | 71.3% | every one "on"; no pair is stably off |
| switch (removal silences it) | 16,295 | 72.2% | every one "not a switch"; no pair is a stable switch |
| on/off and switch together | 14,680 (14,683 rules) | 65.0% | |

The sample, through the runtime, agrees: on/off 702 of 1,000 (Wilson 95%
67.3%-73.0%), switch 720 (69.1%-74.7%), both 644 (61.4%-67.3%); each census
share lies inside its interval. Among the 17,356 activating pairs, on/off is
stable for 12,496 (72.0%) and switch for 11,075 (63.8%).

**Ranks.** Of 3,027,815 comparisons between two pairs in one cell context, the
order of their levels is fixed for 25,629 (0.85%), of their effect sizes for
59,628 (1.97%), and of their response times for none.

**Under R3** (the gene's own mouse numbers say nothing): on/off is stable for no
pair, switch only for the 5,220 inhibiting pairs that removal can never silence,
level and response-time ranks for no comparison, and effect-size ranks for 10 of
3,027,815. Only direction and fold survive, as expected.

**Against the registration.** Level, effect size and response time stable for
no pair: as derived, because one width gives every pair the same spread (61.2-
and 28.2-fold) and whether a level is stable is decided by the width, not by the
pair. Direction and fold stable for all: as derived, being inputs. Combination
rule: 0 differences, as expected. On/off stable for more than half and fewer
than all: **met**, 71.3%. Rank stable in fewer than half of comparisons:
**met**, 0.85%, 1.97% and 0%. The two sentences the amendment named are wrong
as written, by arithmetic: switch is stable for more pairs than on/off (16,295
against 16,101) and, under R3, for 5,220 pairs, all of them inhibitions; among
activating pairs both hold (11,075 against 12,496; 0 under R3). One registered
expectation **failed on the data**: under R3 rank is stable for "no pair", and
10 effect-size comparisons are stable, because effect size carries the observed
fold, which R3 does not vary.

**What the kinds of answer mean.** The answers that survive are the ones that
never depended on the transferred numbers (the direction and size of the fold,
which are the observation) and the threshold answers far from their threshold:
a gene whose nominal level is at least 13 molecules per cell stays expressed at
every point, and its removal does not silence it. Every number the simulation
adds to the observation - how many molecules, how large the change in
molecules, how fast the change - is unstable for every pair, and the order of
two genes by level survives for fewer than one comparison in a hundred. The
nominal level itself sits at the top of the transfer range: if the mouse rate
carried over, the human K562 half-lives would put the level 1.6 to 13 times
lower than the nominal for the central 90% of genes, and nothing measured says
which reading holds.

**The one human context: K562.** By the registered rule it is K562, the only
cell type holding four of the seven kinds in that very cell (a census of the
cached files): CRISPRi of elements (the ENCODE benchmark's 10,356 training and
1,918 held-out K562 pairs, and the IGVF MHC screen), lentiMPRA (51,375
elements), DNase, and human mRNA half-lives (Schofield, 5,419 transcripts). The
next, GM12878, WTC11 and HepG2, hold two or three; HepG2 has no CRISPRi set. 347
simulable pairs are in K562 (1.5% of the tier; 316 AlphaGenome, 31 CRISPRi), and
199 of them are on a gene with a human K562 half-life; 186 of those 199 fall
inside R2, which is in-sample by construction and not a test.

**What a validated human-kinetics test in K562 would need that is missing.**
No cell type holds the other three kinds, and each is what a test needs. First,
a readout that is not an input: the CRISPRi fold is what the bridge fits, so it
cannot validate anything the model then says. Second, **a time course after an
acute perturbation** in K562 (an inducible CRISPRi or a degron, sampled over
hours, for genes that carry a K562 half-life): with Schofield's human half-life
as the decay constant, the model's one dynamic prediction built wholly on human
K562 inputs is that the mRNA covers half the way to its new level in one K562
half-life, and nothing held measures it; the perturbation's own onset kinetics
would have to be measured too. Third, **an absolute transcription rate or
absolute copies per cell** in K562 (a spike-in-per-cell calibration of TT-seq or
RNA-seq), without which no level in molecules is human; none exists
(lane-rates2). Today 199 K562 pairs could run on a human degradation constant,
but their transcription rate would still be mouse, and a test of them would
still need the time course. Nothing was built here.

## The mesoderm the sign correction removed: a closed candidate list (registered 2026-10-02)

Correcting a rule to match its literature source made the model worse against the measurement. At
`5cbce26`, with `Tbxt activates SOX17` as Lolas et al. 2014 show (PNAS 111:4478, Fig. 1A ChIP-seq
and Fig. 3A knockdown), the gastrulation model lost its mesoderm band entirely, 0.10 to 0.00, while
the measured CS7 human range for mesoderm is 0.694 to 0.773. A program that became more faithful to
a published result lost the germ layer the real embryo is mostly made of, and that has a forced
consequence: something in the program was being compensated for by the wrong sign. Either another
rule is also wrong in a direction that cancelled it, or a required species or input is absent, or
the runtime semantics are wrong somewhere the existing tests did not reach.

**The combination rule is eliminated in advance and no run is spent on it.** `93caf61`, `267cc88`
and `a6fc5c4` record that mesoderm is absent under all four activator rules — mean 0.633/0.000/0.367,
sum_capped 0.558/0.000/0.442, max 0.575/0.000/0.425, or 0.558/0.000/0.442 — that only SOX17 has two
activators among the 41 hand-written programs, and that the averaging did not remove the band. The
earlier attribution to activator averaging, `3f7b0f5`, is withdrawn.

The candidate list is closed, fixed before any run, and may not grow between the single-component
pass and the pairwise pass. It is written as committed code, `REGISTRATION` in
`scripts/mesoderm_diagnosis.py`, and not as a data file, so it provably cannot grow without a commit
saying so. It has ten entries: W1 and W2, the two limbs of
the WNT3A-Brachyury loop (Yamaguchi et al. 1999, Genes Dev 13:3185; Martin & Kimelman 2008, Dev Cell
15:121); F1 and F2, the two limbs of the FGF-Brachyury loop (Isaacs et al. 1994, EMBO J 13:4469;
Ciruna & Rossant 2001, Dev Cell 1:37; Schulte-Merker & Smith 1995, Curr Biol 5:62); M1, a MIXL1
mesendoderm node (Hart et al. 2002, Development 129:3597); E1, an EOMES node (Arnold et al. 2008,
Development 135:501; Teo et al. 2011, Genes Dev 25:238; Costello et al. 2011, Nat Cell Biol 13:1084);
A1, a NODAL antagonist (Perea-Gomez et al. 2002, Dev Cell 3:745), registered and **not runnable**
because NODAL is not a modelled species here — its gene declares `max 0` and the runtime holds Nodal
at a clamped external level in every Runge-Kutta stage, so no rule the module can state will move it;
and R1, R2, R3, the three leave-one-out entries over rules already in the program (Lolas et al. 2014
Fig. 3C; Thomson et al. 2011, Cell 145:875, twice).

**Every entry cites a paper that says it acts in mesoderm specification, and the entries that could
not are excluded and named.** X1, leaving out `Nodal activates SOX17`, is excluded because
Kanai-Azuma et al. 2002 reports depletion of definitive gut endoderm and says nothing about
mesoderm; it is kept as a characterisation run outside the list and cannot clear the bound. X2, the
runtime's `a = 1.0` for a gene with no activators, is excluded because no paper about mesoderm
specification can source a runtime convention; it is reported as arithmetic. X3, the module's
unsourced thresholds, X4, the `expected_*` proportions, and X5, `nodal_max` and `decay_length`, are
excluded as interventions: adjusting any of them is tuning a free number, refused at `8bb9123` and
`19423bd`.

**The bound is deliberately weak and is not the CS7 range.** Mesoderm share at or above 0.01 at the
program's final stage, sustained over the last 10% of simulated time, under the default numbers —
120 cells, so two cells, not floating-point dust. Clearing it means mesoderm becomes non-zero at
all. Matching 0.694 to 0.773 would be tuning to the target and is not the goal.

**No number is chosen per candidate.** Every added gene takes the `max: 10; basal: 0.05` that TBXT
and SOX17 already carry, every added protein the `half_life: 2` all four already carry, every added
activating rule the `strength: 1.0; threshold: 1.0; hill: 2` of the module's own `Nodal activates
TBXT`, and every added inhibiting rule the `strength: 1.0; threshold: 2.0; hill: 3` of its
`Tbxt inhibits SOX2`. Nothing is swept.

**Method.** Leave-one-out and set-one-in over the closed list, every entry reported whether or not
it helps. Then, only if no single component clears the bound, a pairwise pass over the same list:
n = 10, n-choose-2 = 45, of which the nine pairs containing A1 are not runnable for A1's stated
reason and 36 are run. Every pair reported. A two-component compensation is exactly what
single-component runs miss, and it is this lane's leading hypothesis, because the program holds no
positive feedback on TBXT at all and a sourced autoregulatory loop has two limbs.

**What stays true whatever is found.** The sign correction at `5cbce26` stands and is not reverted:
`Tbxt activates SOX17` in every run, and a test asserts that no entry removes or flips it. The CS7
comparison stays falsified on all three layers, with `d55cb19`'s bound carried verbatim: "a sampled
census of one embryo is not a 1-D axis, and the model names no stage." If the diagnosis lands on a
threshold, `d55cb19`'s NODAL ceiling travels with it — no measured NODAL threshold maps to model
units, and Dubrulle et al. 2015 argue induction kinetics rather than a clamped level set fates, so a
threshold the model cannot source is a finding, not a fix. And whichever candidate restores the
band, the result says the program lacked it; it does not say it specifies mesoderm. This diagnoses a
program.

**Disclosed, because it shaped the hypothesis before the list was closed.** On the unmodified
committed module, with no candidate set in or left out, Tbxt's final level never exceeds its
competitors anywhere along the gradient: at clamped Nodal 6.0, 4.0, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5 the
final Tbxt reads 0.244, 0.274, 0.414, 0.775, 1.206, 0.256, 0.225, 0.215, and its highest value,
1.206 at Nodal 2.0, stands against Sox2 5.667 and Sox17 4.013 at the same point.

**Citations are unverified bibliography.** This lane is free: 0 model requests, no money, no
network. Every citation above is given from the lane's own knowledge and was not fetched or re-read;
volumes and pages are as recalled. The repository could not confirm them either — it holds WNT3A
only as CellPhoneDB ligand-receptor pairs in `genomeos/std/signalling.bio` and EOMES, MIXL1 and
WNT3A only as marker-gene names in `genomeos/lib/catalog.py`, and no regulatory citation for any of
them.

## The result: the program lacked any repressor of SOX17 but the rule whose sign was wrong (2026-10-02)

Twelve runs reported, from the closed list registered at `0790d51` before any of them: the ten
entries, the baseline, and the X1 characterisation that is outside the list and cannot clear the
bound. Result file `data/results/mesoderm_diagnosis_single.json`, written by committed code with
`own_code_is_committed` true and nothing of another lane's on the counting path. Ecto / meso / endo
at the final stage, with the sustained-tail minimum over the last 10% of simulated time:

| run | ecto | meso | endo | tail min | bound |
| --- | --- | --- | --- | --- | --- |
| baseline | 0.633 | 0.000 | 0.367 | 0.0000 | no |
| W1 Wnt3a activates TBXT | 0.625 | 0.000 | 0.375 | 0.0000 | no |
| W2 Tbxt activates WNT3A | 0.633 | 0.000 | 0.367 | 0.0000 | no |
| F1 Fgf4 activates TBXT | 0.625 | 0.000 | 0.375 | 0.0000 | no |
| F2 Tbxt activates FGF4 | 0.633 | 0.000 | 0.367 | 0.0000 | no |
| M1 MIXL1 node | 0.050 | 0.025 | 0.925 | 0.0167 | **clears** |
| E1 EOMES node | 0.000 | 0.000 | 1.000 | 0.0000 | no |
| A1 NODAL antagonist | — | — | — | — | no run |
| R1 leave out Sox17 inhibits TBXT | 0.367 | **0.633** | 0.000 | 0.6333 | **clears** |
| R2 leave out Sox2 inhibits TBXT | 0.283 | 0.000 | 0.717 | 0.0000 | no |
| R3 leave out Tbxt inhibits SOX2 | 0.642 | 0.000 | 0.358 | 0.0000 | no |
| X1 leave out Nodal activates SOX17 (outside the list) | 0.750 | 0.000 | 0.250 | 0.0000 | n/a |

**The pairwise pass did not run.** Registered condition 4 runs it only if no single component clears
the bound, and two did. n and n-choose-2 stay as registered, 10 and 45, and the list did not grow.

**The two null controls came back exactly on the baseline**, 0.633 / 0.000 / 0.367 to the digit, as
W2 and F2 were predicted to: nothing in the module reads Wnt3a or Fgf4, so the harness could not
have been moving the fates by adding a species.

**Baseline mesoderm is not a band that collapses; it is never there.** Over all 41 recorded samples
of the 40 h run the mesoderm share is 0.000 at every one, not merely at the end.

**Why the drive candidates changed nothing, as arithmetic.** At clamped Nodal 2.0, where Tbxt peaks,
the three genes decompose like this (activator drive `a`, repression factor `rep`, transcription):

| gene | a | activators | rep | repressed by | transcription |
| --- | --- | --- | --- | --- | --- |
| SOX2 | 1.0000 | 0 | 0.09037 | Tbxt 0.8202, Sox17 0.1102 | 1.4037 |
| TBXT | 0.8000 | 1 | 0.02866 | Sox2 0.2601, Sox17 0.1102 | 0.2793 |
| SOX17 | 0.0905 | 2 | **1.00000** | **nothing** | 0.9547 |

TBXT's drive is already 0.80 of its maximum and its repression factor is 0.029. The binding
constraint on TBXT is repression, not drive, so a sourced extra activator cannot restore the band:
W1 and F1 moved mesoderm by nothing at all, and the mesoderm share they produced is 0.000 under the
same bound that R1 clears at 0.633.

**The compensation, named.** `rep = 1.00000` for SOX17 is not a rounding: SOX17 has no inhibitor in
the program at all. The repression graph is Tbxt ⊣ SOX2, Sox2 ⊣ TBXT, Sox17 ⊣ TBXT, Sox17 ⊣ SOX2 —
SOX17 is a sink that represses two genes and is repressed by none. Before `5cbce26`, `Tbxt inhibits
SOX17` was SOX17's only brake, so **the program lacked any repressor of SOX17 other than the rule
whose sign was wrong.** That is what was being compensated for. Correcting the sign did not merely
remove the brake; it turned it into SOX17's second activator, while TBXT keeps two repressors, one of
them SOX17 itself, so TBXT now drives its own repressor with nothing on the other side. The module's
own header, "Mutual repression makes the choice sharp", is contradicted by its rule set: the
repression is not mutual, and for SOX17 it was mutual only through the rule that was wrong.

**R1 is the diagnosis and not a fix.** `Sox17 inhibits TBXT` is sourced by the same paper as the
correction — Lolas et al. 2014, Fig. 3C, Sox17 overexpression lowers Brachyury — so the one
single-component intervention that restores the band removes a cited rule. With its unsourced
strengths, thresholds and Hill coefficients the program cannot hold both of Lolas 2014's results at
once and still make mesoderm. That is a statement about this program, not about the embryo, which
holds both and makes mesoderm anyway.

**M1 clears the bound the other way round, and barely.** Its 0.025 comes with ectoderm collapsed to
0.050 and endoderm raised to 0.925, and its tail runs 0.0417, 0.0417, 0.025, 0.0167, 0.025 — three
cells, three, two and two of 120. A `Mixl1 ⊣ SOX2` term removes the ectoderm competitor where NODAL
is present; it does not rescue TBXT against SOX17, whose share grows. Clearing the weak bound is all
that is claimed for it.

**E1 and R2 both went to endoderm, and the corrected sign is why.** E1's second route into SOX17
takes endoderm to 1.000. R2, which lifts the Sox2 brake off TBXT, raises Tbxt — and a raised Tbxt
now activates SOX17, so endoderm rises from 0.367 to 0.717 and mesoderm stays at 0.000. Lifting a
brake off the mesoderm gene makes the model more endodermal; that is the corrected sign acting, and
it is the clearest single sign that the missing term is on SOX17 and not on TBXT.

**The pairs are bounded, by arithmetic and not by a run.** With W2 present, WNT3A's only activator is
Tbxt, which stays near zero, so the Wnt3a term on TBXT is no larger than it was in W1 alone, where
WNT3A ran unregulated at its full `max_rate`; W1 alone gave 0.000. The same holds for F1 and F2. The
two-limb loops are therefore bounded above by a zero. This is an inference from the numbers above,
not a measurement: the pairwise pass did not run, as registered.

**The ceilings travel with the answer.** The CS7 comparison stays falsified on all three layers,
carrying `d55cb19` verbatim: "a sampled census of one embryo is not a 1-D axis, and the model names
no stage." Nothing here was tuned: no strength, threshold, Hill coefficient, basal, max, `nodal_max`
or `decay_length` was varied in any run, the unsourced thresholds were excluded as X3 and are
reported as a finding rather than adjusted, and `d55cb19`'s NODAL ceiling stands — no measured NODAL
threshold maps to model units, and Dubrulle et al. 2015 argue induction kinetics rather than a
clamped level set fates. The sign correction at `5cbce26` is untouched in every run above. And the
finding is about a program: the program lacked a repressor of SOX17. It is not a claim that
repressing SOX17 specifies mesoderm.

**And the band that was lost rested on the module's weakest rule.** At `5cbce26^` SOX17 had exactly
one inhibitor, `Tbxt inhibits SOX17`, carrying `evidence: inferred "mesoderm dampens endoderm";
confidence: 0.4`, and its partner `Sox17 inhibits TBXT` carried `inferred "endoderm excludes mesoderm
program"; confidence: 0.4`. The mutual-repression toggle that held the 0.10 mesoderm band was built
out of two inferred rules at the module's lowest confidence, and the source contradicts one of them.
So the 0.10 was not a result the program had earned, and losing it is not a regression in anything
the program could show: `0.10 to 0.00` is a guess being withdrawn, not a measurement being lost.
