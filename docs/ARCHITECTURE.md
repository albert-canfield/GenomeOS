# GenomeOS Architecture

## 1. The problem GenomeOS solves

Existing tools approach biology from one of two ends. Genomics tools go
*DNA → molecular prediction* (annotation, variant effect). Physiology tools go
*equations → simulation* (SBML, CellML, agent-based tissue models). Nobody owns
the middle: a standard, executable abstraction that starts from a specific
genome and runs upward through cells, development and ageing, with the
uncertainty of current biological knowledge carried along explicitly.

GenomeOS is that middle layer.

## 2. The stack

```
┌──────────────────────────────────────────────────────────────────────┐
│ BioForge      design & experiment: fork a twin, perturb, compare      │
├──────────────────────────────────────────────────────────────────────┤
│ BioTwin       one genome + one bootstrap state + one environment      │
├──────────────────────────────────────────────────────────────────────┤
│ BioLang       the source language. Files: *.bio                       │
│ BioLib        reusable modules written in BioLang or imported         │
├──────────────────────────────────────────────────────────────────────┤
│ Compiler      FASTA/VCF/GFF3/SBML/CellML/GO/Reactome/AI → BioIR       │
├──────────────────────────────────────────────────────────────────────┤
│ BioIR         typed entities, contextual rules, evidence, UNKNOWN     │
├──────────────────────────────────────────────────────────────────────┤
│ BioVM         engines: central dogma · network (ODE/SDE) · cell       │
│               events · (planned) stochastic molecular · spatial       │
├──────────────────────────────────────────────────────────────────────┤
│ Genome engine lossless sequence store; linear now, graph-ready        │
└──────────────────────────────────────────────────────────────────────┘
```

Python package mapping (v0.1):

| Layer | Package | Status |
|---|---|---|
| Genome engine | `genomeos.genome` | working: FASTA (gz), Sequence, Locus, Genome.fetch |
| BioIR | `genomeos.ir` | working: Entity/Gene/Protein/Region/Rule/Parameter/Module, JSON round-trip |
| BioLang | `genomeos.lang` | working: v0.1 grammar, compile-time reference checks |
| BioVM | `genomeos.runtime` | working: central dogma, network runtime, cell ageing runtime |
| Compiler from external data | — | not started (Phase 1 of the action plan) |
| BioLib, BioTwin, BioForge | — | designed only |

## 3. The state model

The one equation everything follows:

```
S(t + Δt) = Runtime(G, S(t), E, Δt)
```

`G` is the genome. `S` is biological state. `E` is environment. The runtime is
the compiled semantics. Three consequences:

1. **A genome file is not a runnable program.** It needs `S0`, a bootstrap
   state (for a human: the zygote's cytoplasm, ribosomes, mitochondria,
   maternal RNAs). BioTwin always pairs `G` with a bootstrap.
2. **Age is state, not sequence.** Telomere length, methylation age, somatic
   mutations and senescence live in `CellState`, never in `Genome`.
3. **Two twins with the same `G` diverge** if `S0` or `E` differ. That is a
   feature: it is how BioForge compares interventions.

## 4. Execution regimes

Biology runs several kinds of computation at once. BioVM is therefore a set
of engines that a scheduler composes, not one `tick()` loop.

| Regime | Used for | v0.1 | Planned engine |
|---|---|---|---|
| Deterministic mapping | DNA → RNA → protein | `central_dogma` | keep |
| Continuous (ODE/SDE) | gene networks, signalling, metabolism | `grn.NetworkRuntime` (RK4, Euler-Maruyama) | wrap libRoadRunner for SBML-scale models |
| Discrete stochastic events | division, senescence, death, mutation | `cell.CellRuntime` (Poisson/hazard) | keep; add differentiation events |
| Stochastic molecular (SSA/rule-based) | promoter binding, complexes | — | PySB / BioNetGen adapter |
| Boolean / logical | large GRNs per cell | — | MaBoSS adapter |
| Spatial multicellular | morphogen gradients, development | — | CompuCell3D or Tissue Forge adapter |
| Composition / scheduling | all of the above in one organism | — | evaluate process-bigraph (Vivarium 2.0) |

**Composition decision (from the September 2026 research pass).** Vivarium
2.0's `process-bigraph` already provides typed hierarchical state, processes
with ports, multi-timescale scheduling and structural rewrites (division,
death), and the Covert lab's E. coli whole-cell model runs on it. GenomeOS
should not reinvent that scheduler. The plan is: BioIR stays GenomeOS's own
format (it carries evidence and UNKNOWN, which process-bigraph does not), and
BioVM engines are exposed as process-bigraph processes so that a GenomeOS
organism is a composite of them. Decision point is Phase 2 in the action plan.

## 5. Uncertainty as a type

Every `Rule`, `Parameter` and `Entity` has:

```
evidence:   experimental | curated | predicted | inferred | none
source:     citation / database id / model name
confidence: 0.0 .. 1.0
```

`Region.role` may be `UNKNOWN`. The compiler is forbidden from inventing a
role; a predicted role from an AI model is `predicted` with the model named.
`genomeos check` reports mean confidence per kind and lists UNKNOWN regions
and weak rules, so a simulation can say *where* it is guessing.

## 6. Context gating

Rules carry `when:` conditions (`cell_type = neuron`, `stage = gastrulation`).
The runtime activates only rules whose conditions match the current context.
This is how one module can hold knowledge for many tissues without every
rule firing everywhere, and how the same signalling library (WNT, BMP, Notch)
is reused across organs with different outcomes.

## 7. Genome representation

v0.1 stores one linear `Sequence` per chromosome and serves `fetch(locus)`.
The design target is a pangenome graph (HPRC Release 2, GBZ format) in which
an individual's chromosome is a path; `fetch(locus)` is the contract that
must survive that change. Diploid genomes are two `Genome` objects labelled
`maternal` / `paternal` until phased-graph support lands.

Large public files are read, not downloaded: two standard-library readers
fetch byte ranges over HTTP and decode only what a question needs,
`attribution/bigwig.py` for bigWig tracks (a 9.6 GB constraint track costs a
chromosome tens of megabytes) and `genome/bam_range.py` for indexed BAMs (a
600 GB alignment costs a telomere estimate 218 MB). Whatever is kept is the
distilled answer under `data/knowledge` or `data/results`, never the file.

## 8. Resolution levels

Simulating every molecule of every cell is impossible. BioVM will let each
tissue run at its own resolution:

```
L0 genome only        L3 cell population
L1 regulatory network L4 tissue (spatial)
L2 single cell        L5 organ    L6 organism
```

A heart can run at L4 while skin runs at L3, the way a physics engine uses
variable level of detail.

## 9. Debugger

The end-state developer experience: run, pause, step, breakpoints on
biological conditions (`break when telomere < 5000`,
`break when cell.type == cardiomyocyte`), and an explanation trace that
follows a differentiation event back through the rules and evidence that
produced it. The evidence fields in BioIR exist so that this trace is possible.

## 10. BioLang as a separable engine (goal, 2026-09-10)

GenomeOS builds the language, but the language should be able to leave home.
The picture to keep in mind:

```
today                                   later
─────                                   ─────
GenomeOS                                biolang        (engine + toolchain, like Node)
├── genomeos.lang    parser             ├── parser, IR, VM, std library, package manager
├── genomeos.ir      BioIR              ├── `bio run | check | compile | test | repl`
├── genomeos.runtime engines            └── embeddable API (Python first, others later)
├── genomeos.std     prelude
└── everything else  (data, twins,      GenomeOS       (an application on the engine:
    web UI, knowledge, results)         genome decoding, twins, anatomy, web UI)
```

Boundaries we keep now so that the split is possible without a rewrite:

1. **BioIR JSON is the contract.** Anything that talks to the engines does so
   through `Module.to_dict()` / `from_dict()`. Nothing in the language layer
   reads GenomeOS data paths (`data/…`) or results.
2. **Engines take a Module and return a Trajectory or a typed result.** No
   engine reaches into the CLI, the web server or the knowledge files.
3. **Standard library in BioLang, not Python.** What a program needs from the
   environment (cell types, ageing events, signals) is written as `bio.std.*`
   modules with evidence, so a future `biolang` package ships them unchanged.
4. **One thin entry point.** `genomeos bio FILE` (planned, v0.3) runs a file's
   `experiment` blocks; a `bio` alias is a one-line rename. The web UI and the
   other commands are GenomeOS features on top, never dependencies of the run.
5. **No dependencies in the language core.** Parser, IR and built-in engines
   stay pure Python; heavier engines (process-bigraph, roadrunner, MaBoSS) are
   adapters behind the same interfaces.

The step that turns the language into a runtime is the v0.3 program layer:
`experiment` and `organism` blocks (context, bootstrap state, clamps, run,
breakpoints, observe, assert, report, fork and compare), executed by a small
scheduler that maps clauses onto the engines. After that, the engine can be
packaged separately with its own tests, and GenomeOS becomes its first user.
