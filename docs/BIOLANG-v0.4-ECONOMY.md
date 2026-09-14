# BioLang v0.4: structure, economy and control (specification)

**Status: DRAFT for Albert's review, 2026-09-14.** Nothing past stage 1 is to
be implemented until he has steered the constructs below. Stage 1
(compartments and transport) is being built against its gate; its result is
recorded in §9 when it lands, pass or fail. Sections marked *open* are
decisions this document deliberately leaves to him.

Adds to v0.3 (see [BIOLANG-v0.3.md](BIOLANG-v0.3.md)). Every v0.1 to v0.3
program compiles and runs unchanged: the new layer is opt-in per program, and
a program that declares one `compartment` becomes a *located* program, to
which the stricter rules below apply in full.

## 1. Why

BioLang today can say that a gene produces a protein and that a rule fires
when a condition holds. It cannot say three facts every cell obeys:

1. **Everything happens somewhere.** A protein made in the cytosol is not in
   the mitochondrion until something carries it there.
2. **The machinery is finite.** Ribosomes, polymerases, tRNAs and import pores
   are counted, and what one gene takes another does not get.
3. **Everything costs something.** A peptide bond, a transcribed nucleotide
   and a pumped ion are paid for in ATP, and ATP is made by a metabolism that
   needs glucose and oxygen.

A model without them predicts a cell that expresses everything at once, at no
cost, everywhere. That is the main reason such models are decorative rather
than useful. A fourth fact comes with the other three: **a cell is a control
system**, keeping a small number of variables inside viable ranges while it
performs a role. The genome supplies capabilities and rules; what the cell
does next comes from its state.

## 2. The rule that governs this layer

An economy layer multiplies parameters, and unvalidated parameters produce
confident nonsense. So:

- Every quantity introduced here (a volume, a capacity, a cost, a set point,
  a gain) is a fact with `evidence` and `confidence`, exactly like a rule.
- Anything computed from a construct whose stage has not passed its gate is
  reported as `inferred`, and its confidence is capped at **0.3**, whatever the
  program's author wrote. The cap lifts per construct when its gate passes and
  the pass is recorded in this document.
- The uncertainty report (`runtime/uncertainty.py`) carries confidence along
  the chain that produced a value: a protein's level in the mitochondrion is no
  better grounded than the weakest of its gene, its `produces` rule, the mRNA
  export, the import transport and the compartment it sits in. The output
  names that weakest link.
- The Evidence explorer (`genomeos evidence`) counts the new blocks like any
  other; the share of this layer that is measured rather than inferred is
  quoted in every report on it.
- A stage that does not reproduce its measurement is reported as failed and
  its constructs stay optional. None is shipped because it runs.

## 3. Licence boundary

Everything in §4 to §8 that is language or runtime (the parser, the IR types,
the located runtime, the regime and allocation code, the homeostat
controller, the standard-library modules of constants) lives in
`genomeos/lang`, `genomeos/ir`, `genomeos/runtime` and `genomeos/std`, and is
**Apache 2.0** (LICENSING.md). The engine never imports the application:
validation data (MitoCarta, UniProt locations, PaxDb, turnover tables) is read
by **AGPL** scripts under `scripts/` and application modules, which generate
`.bio` programs or pass numbers in. Measured results go to `data/results/`
under their sources' terms. `tests/test_engine_boundary.py` enforces the
direction; this document marks each item with **[engine]** or **[app]**.

## 4. Structure (stage 1)

### 4.1 `compartment` [engine]

```
compartment Extracellular { membrane: no; evidence: curated "Alberts MBoC 6e"; confidence: 0.9 }
compartment PlasmaMembrane { parent: Extracellular; membrane: yes }
compartment Cytosol       { parent: PlasmaMembrane; volume: 0.54; translation: yes
                            evidence: experimental "Alberts MBoC 6e Table 12-1 (hepatocyte)"; confidence: 0.7 }
compartment Nucleus       { parent: Cytosol; volume: 0.06; genome: nuclear }
compartment Mitochondrion { parent: Cytosol; volume: 0.22; genome: chrM; translation: yes; copies: 1 }
compartment ER            { parent: Cytosol; volume: 0.12 }
```

| property | form | meaning |
|---|---|---|
| `parent` | `Id` | the compartment that contains this one; exactly one root |
| `membrane` | `yes \| no` | a membrane faces its parent and its children at once |
| `volume` | `number` | fraction of the cell's volume (recorded in stage 1; used from stage 2) |
| `genome` | `chrM, ... \| nuclear` | chromosomes read here; `nuclear` = every chromosome except chrM/MT |
| `translation` | `yes \| no` | ribosomes are present (cytosolic or mitochondrial) |
| `copies` | `integer` | copies per cell (recorded in stage 1; heteroplasmy in stage 3) |

**Semantics.** Compartments form a containment tree. Two compartments are
*adjacent* when one is the other's parent. This is the same shape as a
process-bigraph place graph, so a located program maps onto a Composite whose
stores nest the same way (`runtime/compose.py`).

**BioIR.** `Compartment(Entity)`: `parent`, `membrane`, `volume`, `genome`,
`translation`, `copies`; `kind = "compartment"`.

### 4.2 Locations [engine]

```
gene MT-ND1  { locus: chrM:3306-4262(+); location: Mitochondrion; produces: MT-ND1p }
gene SDHB    { location: Nucleus; produces: SDHBp }
protein SDHBp { location: Mitochondrion; signals: presequence
                evidence: curated "UniProt P21912 transit peptide 1-28"; confidence: 0.9 }
```

- In a located program **every gene and protein has a location**. A missing
  one is a compile error, never a default.
- A gene with a locus is located by the compartment whose `genome` holds its
  chromosome; a stated `location` that disagrees is a compile error (a chrM
  gene declared nuclear is a mislocalised gene).
- A protein's `location` is the set of compartments it occupies when it
  works (`location: Nucleus, Cytosol`). `signals` names the targeting signals
  it carries (presequence, signal_peptide, NLS, ...), which is what transport
  machinery recognises. A protein is never moved because of where it is
  declared to be; it is moved because a transport recognises it.
- A species' state is per compartment: `SDHB.mRNA@Nucleus`,
  `SDHB.mRNA@Cytosol`, `SDHBp@Cytosol`, `SDHBp@Mitochondrion` are four
  variables.

**The central dogma, located.** A gene is transcribed in its own
compartment. Its mRNA is translated where it meets ribosomes: in its own
compartment if that one has `translation: yes` (the mitochondrion), otherwise
in the first compartment an mRNA transport carries it to (the cytosol, through
the nuclear pore). The protein is born there. If no route from the gene to a
ribosome exists, `produces` is a compile error.

**Rule legality.** A rule acts at a *site*, the compartment of its target
(a gene's compartment for `activates` and `inhibits`, the product's for
`binds`). The source is readable at the site if it is located there, or is
located in a membrane adjacent to the site. Otherwise the rule is a compile
error: a nuclear factor cannot regulate a mitochondrial gene by being
declared, and a cytosolic kinase acts on a nuclear gene only once a transport
carries it into the nucleus and it lists the nucleus among its locations.

**Mislocalisation is visible, never silent.** After a run the runtime reports
two lists: *stranded* species (declared in a compartment they never reached)
and *ectopic* species (present in a compartment their declaration does not
list). Rules read levels at their site, so a stranded subunit changes every
prediction downstream of it.

**BioIR.** `Gene.location: list[str]`, `Protein.location: list[str]`,
`Protein.signals: list[str]` (empty lists in unlocated programs).

### 4.3 `transport` [engine]

```
transport NuclearPore { from: Nucleus; to: Cytosol; cargo: mRNA; capacity: 1e4; affinity: 100
                        evidence: inferred "order of magnitude"; confidence: 0.3 }
transport TOM_TIM23   { from: Cytosol; to: Mitochondrion; cargo: signal = presequence
                        capacity: 5e3; affinity: 50; via: TOMM40p
                        evidence: curated "Wiedemann & Pfanner 2017, Annu Rev Biochem 86:685"; confidence: 0.6 }
```

| property | form | meaning |
|---|---|---|
| `from`, `to` | `Id` | adjacent compartments (anything else is a compile error) |
| `cargo` | `Id, Id \| mRNA \| signal = S` | what it carries: named species, every mRNA, or proteins carrying signal S |
| `capacity` | `amount per hour` | maximal total flux |
| `affinity` | `amount` | cargo level giving half-maximal flux |
| `via` | `Id` | protein whose presence gates the capacity (knock it out and the route closes) |

**Semantics.** One transporter shared by all its cargos, with competitive
saturation: flux of cargo *i* is `capacity · (x_i/K) / (1 + Σ_j x_j/K)`, with
`x` read in the `from` compartment. This is the standard form for one carrier
with competing substrates, and it already makes transport *finite*: a flood of
one cargo slows every other. With `via`, the capacity is multiplied by the
`via` protein's Hill occupancy in either adjacent compartment.

**BioIR.** `Transport(Entity)`: `from_compartment`, `to_compartment`, `cargo`,
`capacity`, `affinity`, `via`; `kind = "transport"`.

### 4.4 Stage 1 gate (what would falsify it) [app data, engine run]

Program generated from MitoCarta3.0 (Rath et al. 2021, NAR 49:D1541) and the
verified chrM translation, with UniProt signals:

1. The 13 chrM proteins are born and stay in the mitochondrion with **no
   transport** used, at zero anywhere else.
2. Every nuclear-encoded MitoCarta protein reaches the mitochondrion only
   through a transport: remove the import transports and the count reaching it
   must fall to zero.
3. The published rho0 phenotype (King & Attardi 1989, Science 246:500):
   with the mitochondrial genome removed, complexes I, III, IV and V have no
   activity and complex II, entirely nuclear-encoded, keeps it. Blocking
   import removes all five. Stripping one subunit's presequence strands it in
   the cytosol and removes only its own complex, and the stranded list names
   it.
4. The red blood cell: a located program with no nucleus, no mitochondrion and
   no gene at all compiles and runs on pre-loaded proteins. A runtime that
   assumes a genome fails this gate.

Falsified if any of the four does not hold on the generated program.

## 5. Economy (stage 2, not implemented)

### 5.1 `pool` [engine]

```
pool Ribosomes { location: Cytosol; size: 5e6; regenerates: 1e5 /h
                 evidence: experimental "HeLa, BioNumbers (to verify against the primary source)"; confidence: 0.5 }
pool RNAPII    { location: Nucleus; size: 6.5e4
                 evidence: experimental "Kimura et al. 1999, Mol Cell Biol 19:5383 (HeLa)"; confidence: 0.6 }
pool ATP       { location: Cytosol; size: 3e9; regenerates: from Glycolysis, OXPHOS }
```

A pool is a finite, counted resource with a size (capacity), a regeneration
rate or the processes that regenerate it, and a location. Pools are drawn on
by `cost` clauses and are returned by release (a ribosome is freed when a
protein is finished; ATP is converted to ADP, not destroyed).

**BioIR.** `Pool(Entity)`: `location`, `size`, `regeneration`, `returns_as`.

### 5.2 `cost` [engine]

```
gene SDHB { ...; cost: RNAPII 1 per transcript; ATP 2 per nucleotide }
protein SDHBp { ...; cost: Ribosomes 1 per chain; ATP 4 per residue }
event divide { ...; cost: ATP 1e10 }
```

A cost is paid from a named pool per unit of work. Defaults for the
per-residue and per-nucleotide costs live in `bio.std.energetics`, each with
its citation: about 4 ATP equivalents per peptide bond (two for aminoacylation,
two GTP for elongation) and about 2 per transcribed nucleotide in polymerisation
alone, with the full biosynthetic costs of Lynch & Marinov 2015 (PNAS
112:15690) kept separate so a program states which accounting it uses.

### 5.3 `allocation` [engine]

The allocation policy is a **scientific claim**, so it is named in the
program and never hidden in the engine:

```
allocation ribosome_share { pool: Ribosomes; policy: proportional }      # share ∝ demand
allocation stress         { pool: ATP; policy: priority; order: Na_K_ATPase, translation, transcription }
allocation growth_law     { pool: Ribosomes; policy: optimise; objective: growth
                            evidence: inferred "modelling device (Scott et al. 2010)"; confidence: 0.2 }
```

Policies: `proportional` (demand-weighted shares), `priority` (an ordered
list), `competitive` (saturable, as for transports), `optimise` (a declared
objective, always labelled a modelling device, never the default; see §7.3).

**Gates, both required.**

- (a) *Burden.* Expressing a costly protein must reduce the output of
  unrelated genes sharing the pool, in the documented direction and order of
  magnitude: Ceroni et al. 2015 (Nat Methods 12:415) in *E. coli* and Frei et
  al. 2020 (Nat Commun 11:4641) in mammalian cells. Falsified if no reduction
  appears, or if it appears with the pool switched off.
- (b) *Abundance.* Allocation must improve the prediction of absolute protein
  copy numbers from transcript abundance, against PaxDb (Wang et al. 2015,
  Proteomics 15:3163) or another open absolute set, relative to the current
  one-to-one mapping. **Stated in advance:** a shared pool with proportional
  allocation rescales every gene by the same factor, so it cannot change a rank
  correlation. The test is therefore on absolute error in log copy numbers and
  on the total protein per cell (Milo 2013, BioEssays 35:1050), and a rank
  improvement is only credited to terms that differ per gene (length cost,
  half-life, translation efficiency). If (b) fails, the pool layer stays
  optional and this document says so.

## 6. Energy and mitochondria (stage 3, not implemented)

A reduced metabolic module (glycolysis, TCA, oxidative phosphorylation) in
`bio.std.core_metabolism`, producing the ATP the pools spend, with
stoichiometries from the textbook and yields with their spread (about 2.5 ATP
per NADH and 1.5 per FADH2, Hinkle 2005, BBA 1706:1; about 30 to 32 ATP per
glucose aerobically against 2 anaerobically).

Mitochondria are compartments with `genome: chrM` and `copies: N`. Each copy
may carry its own genome variant, which makes heteroplasmy expressible:
`heteroplasmy: m.3243A>G = 0.7`, with the biochemical threshold effect
(Rossignol et al. 2003, Biochem J 370:751) as the check.

**Gate.** (i) The ATP budget of a modelled cell type lands within the
published order of magnitude, and its split across consumers matches the
hierarchy of Buttgereit & Brand 1995 (Biochem J 312:163: protein synthesis
and the Na⁺/K⁺-ATPase first). (ii) Lowering oxygen shifts production to
glycolysis and raises glucose consumption (the Pasteur effect); removing
glucose collapses a glycolysis-only cell. (iii) The red blood cell (no
mitochondria) runs on glycolysis alone and keeps its ATP and 2,3-BPG within
the ranges of Mulquiney & Kuchel 1999 (Biochem J 342:581).

## 7. Control (specified here, implemented after stage 2)

### 7.1 `homeostat` [engine]

Homeostasis is **feedback**: a sensor reads a variable, a controller compares
it with a set point, actuators act with a gain. It is never an objective the
runtime optimises.

```
homeostat cytosolic_pH {
  variable: pH@Cytosol
  set_point: 7.2; band: 7.0..7.4
  evidence: experimental "Casey, Grinstein & Orlowski 2010, Nat Rev Mol Cell Biol 11:50"; confidence: 0.8
  sensor: SLC9A1p                   # NHE1 is its own proton sensor
  actuators: SLC9A1p exports H+, SLC4A7p imports HCO3-
  gain: 2.0 /h                      # inferred until the recovery rate is fitted to a published trace
  cost: ATP via Na_K_ATPase         # stage 2: the Na+ gradient it spends
}
```

| property | form | meaning |
|---|---|---|
| `variable` | `species@compartment` | the regulated quantity |
| `set_point`, `band` | `number`, `lo..hi` | the published value and viable range |
| `sensor` | `Id` | what reads the variable |
| `actuators` | `rule or transport ids` | what the controller drives |
| `gain` | `number /time` | response strength (proportional; `integral:` optional) |

Each number carries its own evidence. **BioIR.** `Homeostat`: `variable`,
`set_point`, `band`, `sensor`, `actuators`, `gain`, `integral`, evidence and
confidence.

Candidates with published set points, in the order they would be built:

| Homeostat | Set point | Actuator whose removal is the test | Documented failure |
|---|---|---|---|
| Na⁺ and K⁺ (pump-leak) | Na⁺ ~10–15 mM in, K⁺ ~140 mM in | Na⁺/K⁺-ATPase (ouabain) | Na⁺ rises, K⁺ falls, the cell swells (Tosteson & Hoffman 1960, J Gen Physiol 44:169) |
| Cytosolic pH | ~7.2 | NHE1 (amiloride, knockout) | no recovery after an acid load (Casey et al. 2010) |
| Energy charge | 0.8–0.9 (Atkinson 1968, Biochemistry 7:4030) | AMPK sensing | charge falls further under the same demand |
| Cytosolic Ca²⁺ | ~100 nM | PMCA and SERCA | sustained Ca²⁺ rise |
| Redox (glutathione) | GSH:GSSG >30:1 cytosol, ~1–3:1 ER (Hwang et al. 1992, Science 257:1496) | glutathione reductase | ratio collapses under oxidant load |
| Folding load (ER) | BiP free fraction | IRE1, PERK, ATF6 | unresolved UPR, then apoptosis |
| DNA damage | lesion count | p53 pathway (compiled) | arrest lost, damage carried into division |
| Osmolarity and volume | ~290 mOsm/kg | regulatory volume decrease channels | no recovery after swelling |

**Gate.** At least two homeostats hold their variable inside the published
band under a perturbation, and each fails **in the documented direction** when
its actuator is removed. A controller that cannot be broken correctly is not a
controller. First pair: Na⁺/K⁺ and cytosolic pH.

### 7.2 Specialisation is state, not inheritance

A neuron does not extend a generic cell with extra code. Every cell carries
the same genome and differs in which parts are readable: the project's reader
measures 42% to 77% of coding genes read per cell type, from the same
sequence. So:

- **There is no `extends` and there will not be.** A cell type is a *state*
  that gates one shared rule set, as `cell_type` context and `when` guards
  already do in BioIR. Subtype inheritance would encode a falsehood and force
  rules to be copied per type.
- Specialised behaviour is selection from shared code: the same `homeostat`,
  `transport` and `pool` blocks apply in every cell whose state makes their
  species present.
- Behaviour we cannot yet derive from mechanism is written as a `role`:

```
role contraction {
  when: cell_type = Cardiomyocyte
  performs: contraction; consumes: ATP 1e8 /s
  basis: phenomenological           # required; `mechanistic` only once its rules exist
  evidence: inferred "stated behaviour, no mechanism"; confidence: 0.2
}
```

A `role` with `basis: phenomenological` is counted as `inferred` in every
report, whatever evidence is written, and its confidence is capped at 0.3, so
the Evidence explorer shows exactly how much of a cell's specialised behaviour
is stated rather than derived.

**Fate and commitment (open, to be folded in).** Fate is a function of
position, signals, state, genome, epigenome and developmental time, with
commitment as a ratchet that makes identity progressively harder to reverse.
genomeos-d1 (area E) is sending the semantics it needs for `commitment` and
`competence` blocks and the runtime support for real cell-to-cell contacts, so
that lateral inhibition emerges rather than being stated. They will be
specified in this section, not in a second document.

### 7.3 No central loop

Nothing in a cell runs `while alive:`. The chemistry is parallel and
asynchronous, and the order in which a simulator updates it is a modelling
choice that changes results: in `runtime/boolean.py` synchronous and
asynchronous updating give different attractors for the Fauré et al. 2006
cell cycle (Bioinformatics 22:e124). So the engine bakes in no loop order.

- Temperature is not regulated by a human cell. It is set by the organism
  and the environment and belongs in `organism { environment: temperature = 37 }`.
- Apoptosis is not a runtime exception. It is a regulated program with its
  own machinery, written as a `decision` with `action: die` and guards
  (caspase activity, the p53 pathway), which v0.3 already supports.
- An objective function appears only as a declared `allocation` or flux
  balance with `policy: optimise`, labelled a modelling device, with low
  confidence. Fitting behaviour to an assumed purpose is forbidden by the
  project's evidence rules.

## 8. Execution regime (declared per run)

```
regime default {
  treatment: auto          # continuous | stochastic | auto
  threshold: 50            # auto: species under this many copies are treated stochastically
  units: copies            # stochastic treatment is refused for arbitrary units
  update: continuous       # continuous (ODE) | synchronous | asynchronous | event
  allocation: competitive  # default policy for pools and transports with none named
  seed: 0
}
```

Low copy numbers are stochastic and high ones are not, and transcription is
bursty at single-cell level (Suter et al. 2011, Science 332:472). The runtime
is told, or chooses under `auto`, a continuous or stochastic treatment **per
species**, and every result records the regime it ran under: the treatment
of each species, the update scheme, the allocation policy and the seed. A
result without its regime is not a result.

`stochastic` uses tau-leaping on the same fluxes as the continuous
treatment, so the two differ only in noise. Bursting, when it is added, is a
telegraph promoter (`gene X { bursting: on 0.5 /h; off 5 /h }`), falsified
against the transcriptome-wide burst frequencies and sizes of Larsson et al.
2019 (Nature 565:251).

**Stage 4, division and dilution (not implemented).** A `divide` event
partitions contents (binomially for low copy numbers) instead of duplicating
them, guarded by a size or resource checkpoint, so growth costs something.
Gate: over several divisions dilution must separate a stable protein from a
short-lived one. `bio.std.human_turnover` holds cell lifespans (Sender & Milo
2021), not protein half-lives; the gate needs a protein turnover set (for
example Mathieson et al. 2018, Nat Commun 9:689), which is **open** below.

## 9. Gates and status

| Stage | Constructs | Gate | Status |
|---|---|---|---|
| 1 | compartment, location, signals, transport, regime record | chrM, MitoCarta import, rho0, red blood cell | in progress |
| 2 | pool, cost, allocation | burden; absolute abundance vs PaxDb | not started (waits for Albert) |
| 3 | core metabolism, mitochondrial copies, heteroplasmy | ATP budget; oxygen and glucose dependence; red blood cell glycolysis | not started |
| 4 | partitioning division, checkpoint | dilution vs protein turnover | not started |
| control | homeostat, role | two homeostats hold and break correctly | specified only |

## 10. Open decisions for Albert

1. **Opt-in or mandatory.** Located rules apply to programs that declare a
   compartment. Should v0.5 make locations mandatory everywhere?
2. **Amounts or concentrations.** Stage 1 keeps amounts per compartment, so
   Hill thresholds are amounts. Concentration needs volumes, which vary per cell
   type; the choice changes what every threshold in the project means.
3. **Membranes as compartments.** A membrane is a node of the tree facing both
   sides. The alternative (a membrane as the edge between two spaces) is closer
   to a bigraph link but makes receptors harder to place.
4. **Whose cost accounting.** Polymerisation only, or full biosynthetic cost
   of precursors (Lynch & Marinov). They differ several-fold.
5. **Default allocation.** `competitive` for transports is a choice of
   mechanism; pools may deserve no default at all, forcing every program to name one.
6. **The protein turnover set** for stage 4, and whether a human set is
   required or a mouse one (Schwanhäusser et al. 2011, Nature 473:337) is acceptable.
7. **Which homeostats first.** The proposal is Na⁺/K⁺ and pH because their
   failure modes are the best documented.
