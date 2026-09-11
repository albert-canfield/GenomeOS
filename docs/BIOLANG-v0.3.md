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
```

`when` clauses accept `any`, alternatives `a|b`, and comparisons `>=n`,
`<=n`, `>n`, `<n` (`generation = >=3`). Repeated `assert` and `observe`
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
- **Populations.** A cell may carry `count` > 1; `fraction` on a
  differentiate decision splits that share into a new node. One runtime
  serves named cells (C. elegans) and counted tissues (human).
- **Diamond imports** are merged once; identifier collisions between
  different files remain errors.
- **Uncertainty.** Decisions that fired contribute to the organism level,
  timers used to the cellular level, rules to the molecular level; every
  UNKNOWN stop counts against the organism level with confidence 0.

## Running

```
genomeos grow data/demo/celegans_lineage.bio --until 150 --depth 3
genomeos grow data/organisms/celegans/embryo.bio --until 800 --compare
genomeos grow data/organisms/celegans/embryo.bio --until 6000 --compare --json out.json
```

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
