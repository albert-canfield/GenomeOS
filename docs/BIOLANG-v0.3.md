# BioLang v0.3: the organism layer

Adds to v0.2 (see [BIOLANG-v0.2.md](BIOLANG-v0.2.md)): a program can now
grow an organism from one cell. Six new block kinds, one runtime (the Body),
one command (`genomeos grow`). All v0.1 and v0.2 files compile unchanged.

## The idea

One genome, one bootstrap cell state, one environment; everything after that
is emergent. Each cell reads its context (name, lineage, generation, cell
type, stage, the factors it carries, the environment), the program's
**decisions** say what it does (divide, differentiate, migrate, quiesce, die),
**timers** say when, and **signals** from other cells change what it reads.
A cell for which no decision or no timer applies stops and is reported as
UNKNOWN. Nothing is invented.

## New constructs

```
import bio.std.development           # Zygote, Blastomere, Progenitor, GermCell, PostMitotic

organism Celegans {
  species: Caenorhabditis elegans
  genome: WBcel235
  root: P0                           # name of the first cell (default Zygote)
  resolution: cells                  # cells (named, default) or populations (counted)
  seed: 0                            # timers run with their measured spread from this seed; omit for means
  cell_type: Zygote                  # bootstrap type
  factors: SKN-1, PIE-1, PAL-1       # maternal factors present in the zygote
  environment: temperature = 20
  tempo: 1.0                         # multiplies every timer (species pace)
  observe: count, deaths, fates      # printed at checkpoints
  assert: count at 100 min in 22..34 # checked after the run
  reference: celegans                # ground truth to diff against
  evidence: experimental "Sulston et al. 1983"; confidence: 0.9
}

stage Gastrulation { from: 100 min; to: 350 min; evidence: ... }

timer ab_cycle { duration: 20 min; sd: 2; lengthening: 1.12; when: lineage = AB; evidence: ... }

signal Notch_P2_ABp {
  mode: contact                      # contact | gradient | systemic
  ligand: APX-1; receptor: GLP-1
  from: cell = P2                    # sender condition
  to: cell = ABp                     # receiver condition
  sets: Notch = received             # factor the receiver gains
  evidence: experimental "Mello et al. 1994"; confidence: 0.9
}

decision div_EMS {
  action: divide                     # divide | differentiate | migrate | quiesce | die
  when: cell = EMS, SKN-1 = present, Wnt = received
  daughters: MS, E                   # names; omitted = a/p then l/r suffixes
  lineages: MS = MS, E = E           # daughters that found a new lineage (generation 0)
  asymmetric: POP-1 -> MS            # which daughter keeps which factor
  timer: ems_cycle                   # omitted = first timer whose `when` matches
  evidence: experimental "Lin et al. 1998"; confidence: 0.9
}
decision fate_ABalaaaalal { action: differentiate; when: cell = ABalaaaalal; to: Neuron; name: AINL; ... }
decision die_ABalaaaalar  { action: die; when: cell = ABalaaaalar; after: 85 min; ... }
decision germline         { action: quiesce; when: cell = P4; ... }

domain APP_node { locus: chr21:25880000-26055000; genes: APP; boundaries: EH38E2135001, EH38E2135420 }

field Morphogen { diffusion: 0.4; decay: 0.02; source: 0,0 = 1.0 }

experiment pop1 {
  knockout: POP-1                    # factors never present; signals by id, ligand or receptor never sent
  until: 800 min                     # omitted = the organism's last stage
  expect: "MS adopts the E fate (Lin, Thorpe & Priess 1995)"
  assert: type EPrecursor at 100 min = 2
  evidence: experimental "Lin et al. 1995, Cell 83:599"; confidence: 0.9
}
```

`when` clauses also accept `absent` (`POP-1 = absent`: the factor is not
carried), which is how a mechanistic rule states what happens without a
factor, and what a knockout then triggers.

`when` clauses accept `any`, `absent`, alternatives `a|b`, and comparisons
`>=n`, `<=n`, `>n`, `<n` (`generation = >=3`). Repeated `assert` and `observe`
lines are allowed. `after` and `duration` take a unit (min, h, d, wk, yr).

## Semantics

- **Precedence.** For each action the first matching decision in module order
  wins; imported modules come first in import order. Write mechanism before
  lookup: `founders.bio` (maternal factors, Wnt, Notch) before the generated
  lineage program, and the mechanistic rule decides the cells it names.
- **Division.** The cell waits its timer (`duration × lengthening^(generation-1)
  × tempo`, plus `sd` noise when a seed is given), then two daughters are
  born with the parent's cell type and factors; `asymmetric` keeps a factor
  in one daughter and removes it from the other; `lineages` starts a new
  lineage at generation 0.
- **Signals** are applied when a receiver or a sender is born and both are
  alive; the receiver gains the factor and re-decides, which can replace a
  pending lookup division by a mechanistic one.
- **Populations.** `resolution: populations` in the organism block makes
  every node a counted population instead of a named cell: a `divide`
  without daughters grows the count in place by `fraction` (1.0 = doubling)
  and the population decides again; `differentiate` with `fraction` splits
  that share into a new node; `die` with `fraction` and `after` is a
  recurring loss (turnover) that runs while the decision applies; `quiesce`
  is re-read at every step, so `when: count = >=2.55e13` caps a tissue at
  its adult count and lets it regrow after losses. Populations wait from
  now, cells wait from birth. One runtime serves named cells (C. elegans)
  and counted tissues (human); `turnover_per_day` and `culled` are reported.
- **Spread.** A timer's `sd` is applied when the run has a seed: the organism's
  `seed:` by default, `--seed N` to override, `--means` to switch it off.
  The C. elegans program declares `seed: 0`; at their means the distilled
  timers put 194 cells at 350 min against 274 in the reference, with the
  spread restored 240 to 265 (seeds 0 to 3), so the spread is part of the
  measured biology, not noise added for effect.
- **Flows.** A `differentiate` decision with `fraction` and `after` on a
  population is a recurring flow: every `after`, that share of the pool moves
  into the target pool (created on first use, then merged into), for as long
  as the decision applies. This is how a stem-cell compartment feeds a
  lineage; `data/organisms/human/haematopoiesis.bio` is built from it.
- **Stages** are events: when a stage starts, populations and resting cells
  read the new stage and decide again.
- **Diamond imports** are merged once; identifier collisions between
  different files remain errors.
- **Uncertainty.** Decisions that fired contribute to the organism level,
  timers used to the cellular level, rules to the molecular level; every
  UNKNOWN stop counts against the organism level with confidence 0.

## Experiments

`genomeos grow FILE --experiments [names]` grows the wild type and each
perturbed organism with the same program, seed and horizon, and reports
every cell that decided differently: the founder-window cells first (which
decision was lost or gained, what type the cell took, how many descendants
that touches), then the counts and the experiment's own asserts. The
C. elegans mutants in `data/organisms/celegans/mutants.bio` remove POP-1,
SKN-1, PIE-1, APX-1, GLP-1 and PAL-1; the founder rules in `founders.bio`
are written from the genetics, so each knockout reproduces the published
phenotype at the founder level (MS becomes E-like without POP-1, EMS
daughters become C-like without SKN-1 through PAL-1, P2 becomes EMS-like
without PIE-1, ABp becomes ABa-like without the Notch signal, C and D lose
their fates without PAL-1). Terminal fates below the founders still come
from the observed lineage program, so the report is explicit about where
mechanism stops.

## Design (BioForge over organisms)

A `design` block asks which perturbation reaches a goal:

```
design lose_neutrophils_keep_the_rest {
  knockout_any_of: TAL1, GATA1, KLF1, SPI1, CEBPA, GFI1, IRF8, PAX5, NOTCH1, IKZF1
  at_most: 1                          # perturbations combined per candidate
  vary: timer cycle duration 10..120  # optional continuous knobs (timer duration, decision fraction or after)
  until: 3 yr
  target: type Neutrophil at 3 yr = 0 # loss: normalised distance from holding
  keep: type Erythrocyte at 3 yr >= 2e13   # constraints that must hold
  keep: type Monocyte at 3 yr >= 1e9
}
```

`genomeos grow FILE --design [names]` runs every candidate as an experiment
against the same program, seed and horizon, ranks feasible candidates by
loss and then by the number of perturbations, and emits the answer as an
`experiment` block with `predicted` evidence and low confidence: a design is
a hypothesis until the bench confirms it. Over the haematopoiesis module the
first question above answers GFI1 (C/EBPa would also take the monocytes,
since the GMP needs it); the second answers IKZF1 or TCF3 alone; over the C. elegans founders, "two intestinal
founders" answers POP-1 and "two EMS founders and no germline" answers PIE-1.

## Space

An organism may have a grid (`space: 30 x 1`, `origin: 0,0`, `sense: 30 min`)
and `field` blocks (diffusion, decay, point sources per hour). Cells then
hold a site; a division puts one daughter on the parent's site and the other
on the nearest free site (along `direction:` first), and a division that finds
no free site waits another cycle (contact inhibition, counted). A `signal`
with `mode: gradient` reads a field at the cell's site and sets its factor
while the value is at or above `threshold`; cells re-read gradients every
`sense` interval and decide again when a factor changes. `migrate` moves a
cell by `steps` sites along `direction` or up the gradient of `toward`, into
free sites only, and recurs with `after`. `data/demo/flag_organism.bio`
grows Wolpert's French flag from a single founder: the strip fills, the
gradient forms, and the three bands appear in order. A cell's `x` and `y`
are in its context, so `when: x = >=10` also works.

## Composition

With the `compose` extra, `genomeos.runtime.compose.build_body_composite`
puts an organism program into a process-bigraph Composite beside a
regulatory network. The network's species levels can gate the organism:
a gate names a factor, a species and a threshold, and while the species is
at or above it the factor is present in the organism's context, so any
`when: FACTOR = present` decision reads the network's state. That is the
bridge from molecular dynamics (rules, Hill functions) to the cell
decisions of the organism layer.

## Running

```
genomeos grow data/demo/celegans_lineage.bio --until 150 --depth 3
genomeos grow data/organisms/celegans/embryo.bio --until 800 --compare
genomeos grow data/organisms/celegans/embryo.bio --until 6000 --compare --json out.json
genomeos grow data/organisms/human/body.bio --until "20 yr" --depth 2
genomeos grow data/organisms/celegans/mutants.bio --experiments
genomeos grow data/organisms/celegans/mutants.bio --experiments pop1,skn1 --json mutants.json
```

The human program (`data/organisms/human/body.bio` with the generated
`tissues.bio`) is the same runtime at population resolution: Carnegie-stage
timers, germ layers by stated shares, sixteen tissue populations that grow
to the adult counts of Sender & Milo 2021 and then turn over at their
measured daily rates. Its uncertainty report says what it is: curated
counts and inferred shares, no mechanism.

`--compare` scores the body against `reference:` (cells born by name,
parent topology, division timing, terminal fates, deaths, count curve). The
reference for `celegans` is distilled by
`genomeos data distil --only celegans_lineage`, which also generates
`timers.bio`, `lineage_embryo.bio` and `lineage_larva.bio` under
`data/organisms/celegans/`.

## Standard prelude

| Module | Contents |
|---|---|
| `bio.std.development` | Zygote, Blastomere, Progenitor, GermCell, PostMitotic with CL ids |
