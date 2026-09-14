# BioLang v0.4: structure, economy and control (specification)

**Status: DRAFT for Albert's review, 2026-09-14.** Nothing past stage 1 is to
be implemented until he has steered the constructs below. Stage 1
(compartments and transport) is implemented and passed its four gates; the
numbers are in §9.1. Sections marked *open* are
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
compartment ER            { parent: Cytosol; volume: 0.09 }
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
*adjacent* when one is the other's parent, or when both sit either side of one
membrane compartment, which a transport then crosses. This is the same shape as a
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
  works (`location: Nucleus, Cytosol`); `initial` gives the amount it holds in
  each of them at time zero, which is how a cell with no genome is written. `signals` names the targeting signals
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
| `from`, `to` | `Id` | adjacent compartments, or the two sides of one membrane (anything else is a compile error) |
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

### 7.2a Fate: commitment, competence and how a factor is read [engine]

Fate is a function of position, signals, state, genome, epigenome and
developmental time. Two constructs make the *irreversibility* of that function
expressible, proposed by area E (genomeos-d1) after testing it against the
worm's factor atlas:

```
commitment intestine {
  programme: Intestine                 # what the cell is committing to
  establish: ELT-2.exposure(lineage) >= 0.8   # a level, or a time-integrated read
  locks: cell_type                     # what can no longer change
  maintain: ELT-2 >= 0.2               # what must hold for the lock to stand
  excludes: HLH-1, ELT-1               # factors masked in a committed cell
  hysteresis: enter 0.8, leave 0.3     # different thresholds in and out
  release: never                       # never | after 60 min below leave | experiment
  inherit: daughters                    # the state is inherited as mechanism
  evidence: inferred "Fukushige & Krause 2005; Yuzyuk et al. 2009"; confidence: 0.3
}

competence early_muscle {
  when: stage = Cleavage
  allows: Muscle, Intestine, Hypodermis
  closes: on commitment                # at 350 min | on commitment | after generation 8
  evidence: inferred "Fukushige & Krause 2005 (HLH-1 makes muscle only in early embryos)"; confidence: 0.3
}
```

**Semantics.** A commitment is checked before `differentiate`: once established,
a fate outside the locked programme is refused (and counted, never silently
applied), signals can no longer change `cell_type`, and excluded factors are
masked in that cell's context. A fate outside an open `competence` is refused
and reported as **"outside competence"**, which is a different answer from
UNKNOWN: the program said no, rather than saying nothing. Commitment state is
inherited by daughters when `inherit: daughters`, so it is mechanism rather
than bookkeeping.

**The evidence, stated plainly, because it does not support the picture.** Area
E tested the ratchet on the Ma et al. 2021 atlas and it came out **negative,
and in the wrong direction**: factor states get *less* clustered with
developmental time rather than more (p 0.005), competing programmes are
mutually exclusive from the start rather than becoming so, and mother-to-
daughter persistence does not change with time. The atlas is the wrong
instrument — no perturbation, and reporter protein carries across divisions —
so the justification for these constructs is published perturbation work, not
our own measurement:

- Fukushige & Krause 2005: ectopic HLH-1 makes muscle only in early embryos.
- Yuzyuk et al. 2009: MES-2 ends that plasticity around gastrulation.

**Falsifying measurement (a perturbation, not an atlas).** A program with
`commitment` and `competence` must reproduce a published perturbation series:
forcing a master factor early changes fate, forcing the same factor after the
competence window closes does not, and removing the closing machinery
(MES-2) restores the late response. Until such a series is reproduced, every
number in a `commitment` or `competence` block is `inferred` at 0.3 or below,
and anything derived from one is reported the same way. **Nothing in this
specification claims our own data show a ratchet.**

**The series, run (2026-09-14), and the gate passed.**
`data/organisms/celegans/plasticity.bio` is that series as a program, beside
area E's worm and importing it unchanged. It says three things and nothing
else: HLH-1 is sufficient for a muscle fate (one `decision`); it is sufficient
only while a window is open, and the window is closed by MES-2 (one
`competence`); and it is never sufficient in a cell that has already
differentiated (one `commitment`). The transgene is its own factor, `hsHLH-1`,
as the published heat-shock construct is its own copy of the gene, so the wild
type carries none of it and the worm's own development is untouched — checked:
the 610 cells at 800 min are identical to `embryo_factors.bio`'s, type by type.
Each experiment forces the factor at a stated time, which is new
(`add: hsHLH-1 at 350 min`; before this, `add:` was a zygote load and the late
arm had no form at all). `bio test` runs all five.

| Arm | Published outcome | What the program does |
|---|---|---|
| forced at 60 min | almost every cell becomes muscle, whatever its lineage (F&K 2005) | **610 of 610** are muscle; no neurons, no intestine; 1,340 endogenous fate assignments refused because the cell is already committed to muscle |
| forced at 350 min | after the window the same factor does nothing (F&K 2005) | the embryo keeps its own fates exactly: 117 muscle, 233 neurons, 34 intestine, every type unchanged; 1,180 refusals counted as *outside competence* |
| forced at 350 min, no MES-2 | plasticity is prolonged (Yuzyuk et al. 2009) | the window never closes and **594 of 610** become muscle |
| forced at 350 min, no PHA-4 | losing a regulator that is not the closing machinery does **not** prolong plasticity (their own control) | unchanged from the wild type (117 muscle), so the closing is attributed to MES-2 and not to any knockout |
| forced at 700 min, no MES-2 | lost completely in terminally differentiated cells (F&K 2005) | of the 84 cells that convert, **not one had differentiated**: all are precursors or untyped cells. Every cell already differentiated at 700 min is still itself at 800 (233 neurons, 26 intestine, 99 hypodermis) |

**What each construct is worth**, measured by removing it
(`scripts/plasticity_gate.py`, `data/results/plasticity_gate.json`); this is
the part that matters, because a program that passes with the construct and
also without it has not tested the construct:

| Arm | as written | `competence` removed | `commitment` removed |
|---|---|---|---|
| forced at 60 min | 610 muscle, passes | 610 muscle, passes | **174 muscle, fails** |
| forced at 350 min | 117 muscle, passes | **594 muscle, fails** | 117 muscle, passes |
| 350 min, no MES-2 | 594 muscle, passes | 594 muscle, passes | 610 muscle, passes |
| 350 min, no PHA-4 | 117 muscle, passes | **594 muscle, fails** | 117 muscle, passes |
| 700 min, no MES-2 | 201 muscle, passes | 201 muscle, passes | **610 muscle, fails** |

Both are load-bearing, and in different arms: the window is what makes a late
factor inert, and the lock is what makes a cell that has already differentiated
refuse — and also what keeps a *converted* cell converted, which is why the
early arm collapses to 174 without it (the endogenous fate rules simply
overwrite the ectopic fate later). So the two constructs are not one
construct wearing two names, and the series separates them.

**What the runtime does.** A `differentiate` decision may name the window it
needs (`competence: early_plasticity`); outside the window that decision is
refused and counted as **outside competence**, never silently applied and never
silently dropped. A `commitment` locks `cell_type` when its `establish` clause
matches, refuses any later fate outside the programme (counted separately), and
passes the lock to the daughters with `inherit: daughters`. Both counts are in
every run's summary beside `ambiguous_fates` and `revised_fates`.

**Decided while implementing, for Albert to confirm or overturn.**

1. **A competence governs the decisions that name it**, not every decision that
   happens to target a fate it allows. The alternative (a window that refuses
   any fate in `allows` after it closes) would refuse the worm's own terminal
   fates, which are assigned from 200 min onwards: the window would abolish the
   embryo it was meant to describe. `allows` is still checked — a decision whose
   target the window does not allow is a compile error — and a window that
   governs no decision is a compile error too, since it would be a claim the
   runtime never checks.
2. **`closed_by` names the machinery, so removing it keeps the window open.**
   This is the only way the third arm of the series can exist; it is an addition
   to §7.2a as first written, where a window closed unconditionally.
3. **Commitment is established at differentiation**, so a cell *born after* the
   perturbation is not protected: 12 of 610 cells in the last arm (8 of them in
   the E lineage) take the forced fate at birth, because the program has no
   determination before differentiation — a precursor is committed to nothing.
   That is a real limitation, not a rounding error: in the embryo those cells
   descend from a committed precursor. Stating it would need a `commitment` on
   the founder lineages, which the atlas-driven program does not have.
4. **Clauses specified here but not gated are compile errors, not no-ops:**
   `maintain`, `excludes`, `hysteresis`, any `release` other than `never`, and
   the `.exposure(window)` / `.mean(window)` reads in an `establish` clause.
   The integrated reads in particular are **not implemented**: area E's
   `_integrated` factors are a per-cell lookup generated from the atlas, not a
   quantity the runtime integrates, so a program cannot yet ask for a window it
   did not precompute. §7.2a's example `establish: ELT-2.exposure(lineage) >= 0.8`
   therefore does not compile today.
5. **The confidence cap lifts for these two constructs** (§2), because the gate
   above passed: the program declares 0.6 for the window and the lock and 0.8
   for the rule and the arms, each with its citation. What is *not* earned by
   this run: the window's 180 min is F&K's "first 3 hours" read onto our
   timeline, and the cell counts are our worm program's, not theirs.

**How a factor is read is itself a measurement, and it has one.** Area E
measured three readings of the same rules on 555 terminal cells: the
instantaneous peak makes 62 errors at threshold 0.2, exposure summed along the
lineage 38, and the mean over the cell's own life 29 — and the instantaneous
read is the worst at every threshold from 0.1 to 0.5. So the language gets
numeric factor levels and both integrated forms, with the window **declared
rather than assumed**, because the data do not settle which window a cell
integrates over:

```
when: ELT-2 >= 0.3                     # the level now
when: ELT-2.exposure(lineage) >= 0.8   # summed along the lineage path
when: ELT-2.mean(cell) >= 0.25         # mean over this cell's own life
when: ELT-2.mean(lineage) >= 0.25      # the same two reads over the other window
when: ratio(HLH-1, ELT-1) >= 2
noise: 0.05                            # on a threshold, so a decision is not a knife edge
```

An `exposure` or `mean` read that does not name its window is a compile
error, not a default: the three readings give different answers, and
the difference is the finding. Presence tests (`F = present`) keep their v0.3
meaning.

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

**Fate precedence is declared, not implied by line order (implemented 2026-09-14).**
Until now the Body fired every matching `differentiate` at a decision point in
module order, so the *last* match won, while BIOLANG-v0.3.md said the first
did. Rule order in a file was semantically load-bearing and undocumented, and
area E's `fates.bio` had to list its rules bottom to top to work. The fix is a
regime property plus an explicit precedence:

```
regime worm { fates: first }            # first | last (legacy, still the default)
decision factor_fate_00 { action: differentiate; priority: 24; ... }
```

- `fates: first`: one fate per decision point — the matching decision with the
  highest `priority`, ties to the first in module order. After a change of type
  the chain continues only through a decision the new type *enables*; a
  decision that already applied before the change is a competitor, not a
  successor, so it cannot overwrite the fate.
- `fates: last`: the legacy behaviour, kept and documented rather than removed,
  because it is now a stated mode and not an accident; `priority` with
  `fates: last` is refused, so the two cannot be mixed by mistake.
- **`first` is the default since 2026-09-14** (§10 decision 8, resolved below).
  A program that declares no regime therefore takes one fate per decision point.
  **What an unmigrated program sees:** where two rules matched the same cell, the
  fate is now chosen by `priority` (all zero unless the program says otherwise,
  so by module order) rather than by which rule sat lowest in the file, and
  `ambiguous_fates` in the run's summary says at how many decision points that
  choice was made. A program whose rules are written bottom-to-top for the old
  behaviour should either add priorities or declare `regime { fates: last }`.
- Both modes report `ambiguous_fates` (decision points where equal-precedence
  rules disagreed about the target) and `revised_fates` (a terminal fate changed
  again at a later decision point, which `commitment` will refuse).

Measured on area E's `embryo_factors.bio` to 800 min (1,439 cells): as
committed then, under the legacy default, nothing changed and the runtime reported **45 ambiguous decision
points** — exactly the cells whose fate depended on line order; under
`fates: first` without priorities 45 cells change (40 hypodermis and 5 muscle
become neurons); with priorities emitted by the generator (the file's last rule
highest) **all 1,439 fates are identical, 0 ambiguous, and 620 fewer decisions
fire** (3,707 against 4,327), because the overwritten firings are gone.

Area E landed that migration (`0e93b5b`): the generator writes the rules in
natural order with the first highest (24 down to 1), so every factor rule
outranks the lineage lookup at 0, and `embryo_factors.bio` declares
`regime worm { fates: first }`; their run reproduces the measurement exactly
(1,439 cells, 0 fates differing, 45 ambiguous points to 0, 4,327 firings to
3,707, terminal score unchanged at 496 of 555 and 902 of 961 to the adult).
The default was then flipped, and **every other committed organism program
gives identical results under both modes**: the worm embryo, its designs and
mutants (2,183 cells each, every knockout experiment identical), the human
body, haematopoiesis with its designs and mutants, the lineage demo and the
French-flag demo — 0 cells differing, 0 ambiguous points anywhere. Population
shares are not competitors, so a split is never counted as ambiguous.

### 7.4 Contacts, neighbours and noise (Body runtime, for emergent lateral inhibition) [engine, implemented]

Today a contact `signal` names its sender (`from: cell = P2`), so lateral
inhibition can only be stated. For it to emerge, area E needs the Body to
provide, in this order:

1. **Neighbours per cell**, from grid adjacency when the organism has a space,
   or from an imported time-resolved contact table
   (`contacts: celegans_contacts.tsv`, rows of time, cell, cell, area).
2. **Contact signals against current neighbours**, not named senders: a
   receiver reads the *amount* of ligand summed over its neighbours at that
   moment (`mode: contact; reads: amount`), so the same rule works for every
   pair that happens to touch.
3. **A per-cell network stepped in sync across neighbours between events**
   (the network runtime, one instance per cell, coupled through the contact
   amounts), under the declared `regime` update scheme.
4. **Seeded noise**, recorded in the result. Not a convenience: the anchor-cell
   choice between Z1.ppp and Z4.aaa is close to 50:50, and area E's
   implementation of Collier et al. 1996 shows two equal cells never diverge
   without noise, while with noise 198 of 200 pairs diverge and the first cell
   wins 51% of the time. The noise is the mechanism.
5. **Daughters placed along the division axis their names imply** (a/p, l/r,
   d/v), so neighbour relations follow the lineage.
6. **Replicate runs with a diff that scores alternatives**: an assert form such
   as `exactly one of Z1.ppp, Z4.aaa is AnchorCell` evaluated across seeds, with
   the split reported (for the anchor cell, near 50:50 is the pass).

**Falsifying measurement.** Two equivalent cells with the Collier circuit and no
stated winner must (a) not diverge with noise off, (b) diverge in nearly every
replicate with noise on, and (c) split about evenly across replicates; a
runtime that picks the same winner every time has smuggled in an order.

**Implemented 2026-09-14, and the gate passed.** What the Body now does:

```
organism AcVu {
  space: 2 x 1; placement: names          # daughters along the axis their names imply
  contacts: contacts.tsv                   # or: neighbours from a time-resolved table instead of a grid
  cell_network: 6 min                      # every cell's own network, stepped together
  replicates: 100                          # outcomes decided by noise are scored over seeds
  seed: 0
  assert: exactly one of Z1.ppp, Z4.aaa is AnchorCell at 40 h in >= 90%
  assert: Z1.ppp is AnchorCell at 40 h in 35..65%
}
signal DeltaAtContact { mode: contact; ligand: Dp; reads: amount; sets: Dext }
param noise = 0.05
```

- **Neighbours** (`Body.neighbours`) come from the latest snapshot of a contact
  table at or before now, weighted by contact area, or from the four grid
  neighbours when there is no table. Nothing reads node geometry, so a change
  to where nodes are called (area B's orientation-aware boundaries) moves
  nothing here.
- **Contact amounts.** A `reads: amount` signal gives the receiver the ligand
  summed over the neighbours that match the sender condition — the ligand's
  level in each neighbour's own network, or a count of touching senders when
  there are no networks — so one rule serves every pair that touches.
- **Per-cell networks, synchronous.** Every `cell_network` minutes all inputs
  are read from every cell first, then every cell's copy of the module's
  network is advanced, then every cell decides again on its new levels, which
  its `when` clauses can read (`Dp = >=0.6`). No cell sees another's update
  from the same step.
- **Seeded noise as mechanism.** Each cell draws from its own stream, seeded
  by the run's seed and the cell's name, never by the order cells are stored
  in; a noisy network without a seed is refused. The seed is in the summary.
- **Replicate asserts** (`exactly one of A, B is T in >= P%`, `A is T in
  lo..hi%`) are scored across `replicates` seeds by `bio test`.
- **A revised fate** now drops its old terminal name as well as being counted,
  so a cell that stopped being the anchor cell is not still called AC.

The gate (`scripts/body_contacts_gate.py`, `data/results/body_contacts_gate.json`,
200 seeds, `data/demo/lateral_inhibition.bio`):

| Condition | Result |
|---|---|
| noise off | 0 of 200 runs diverge; the two cells' Delta identical to the last digit; both stay progenitors |
| noise on | **200 of 200** runs give exactly one anchor cell; Z1.ppp wins **96 (48%)**, Z4.aaa 104 |
| daughters created in the opposite order | the **same winner in 200 of 200 seeds** — no order smuggled in |
| a program that names the winner | the share assert fails at 100%, as it must |

Area E's direct reference gives 198 of 200 and 51%. The comparison is
qualitative on purpose: the Body's noise is multiplicative on every network
species, Collier's reference adds it to Delta only. Every committed organism
program and all 17 of their knockout experiments give identical cells on the
Body before and after this change. The circuit's two Hill terms are Collier's
constants (curated); the noise size, time scale and translation rate are
inferred at 0.3, and the program says so.

## 8. Execution regime (declared per run)

```
regime default {
  treatment: auto          # continuous | stochastic | auto
  threshold: 50            # auto: species under this many copies are treated stochastically
  units: copies            # stochastic treatment is refused for arbitrary units
  update: continuous       # continuous (ODE) | synchronous | asynchronous | event
  allocation: competitive  # default policy for pools and transports with none named
  fates: first            # Body: one fate per decision point by precedence (last = legacy)
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
| 1 | compartment, location, signals, transport, regime record | chrM, MitoCarta import, rho0, red blood cell | **passed 2026-09-14** (below) |
| 2 | pool, cost, allocation | burden; absolute abundance vs PaxDb | not started (waits for Albert) |
| 3 | core metabolism, mitochondrial copies, heteroplasmy | ATP budget; oxygen and glucose dependence; red blood cell glycolysis | not started |
| 4 | partitioning division, checkpoint | dilution vs protein turnover | not started |
| control | homeostat, role | two homeostats hold and break correctly | specified only |
| fate | commitment, competence | a published perturbation series (Fukushige & Krause 2005; Yuzyuk et al. 2009) | **passed 2026-09-14** (§7.2a): five arms reproduced, and each construct fails an arm when removed |
| fate | integrated reads (`.exposure(window)`, `.mean(window)`) | area E's three readings on 555 terminal cells | specified only, **not implemented**: the `_integrated` factors are a precomputed lookup |
| fate precedence | `regime fates`, decision `priority` | area E's worm fates identical under explicit priorities | **implemented and now the default**: 1,439 of 1,439 identical, 45 ambiguous points to 0 |
| contacts | neighbours, contact amounts, per-cell networks, seeded noise, replicate asserts | Collier 1996 equivalence group splits about evenly | **passed 2026-09-14**: 0 of 200 diverge without noise, 200 of 200 with it, first cell 48%, same winner per seed in either creation order |

### 9.1 Stage 1 as measured (2026-09-14)

`scripts/located_mitochondrion.py` generates a located program for all 1,136
MitoCarta3.0 genes and runs it on `runtime/located.py`; the result is
`data/results/located_mitochondrion.json`, and the committed programs
`data/organisms/human/oxphos.bio` (generated) and
`data/organisms/human/erythrocyte.bio` test themselves under `bio test`.

| Gate | Result |
|---|---|
| 1. the 13 chrM proteins | 13 of 13 are made inside the mitochondrion, stay there, use no transport, and are unaffected when import is closed |
| 2. nuclear-encoded import | 1,123 of 1,123 reach the mitochondrion with the routes open, **0** with them closed, and all 1,123 are then listed as stranded; closing only the presequence route leaves 418 (the internal-signal proteins), closing only the internal route leaves 705 |
| 3. rho0 (King & Attardi 1989) | without mtDNA complexes I, III, IV and V fall to zero activity and complex II, entirely nuclear-encoded, keeps it; closing import removes all five; removing SDHB's presequence removes complex II alone and names SDHB as the one stranded protein |
| 4. no genome (red blood cell) | a program with 0 genes, no nucleus and no mitochondrion compiles and runs; haemoglobin assembles from the chains the cell was born with and mass balance holds |

Six modelling errors are refused at compile time rather than run: a gene or a
protein with no location, a chrM gene declared nuclear, a nuclear gene whose
mRNA reaches no ribosome, a rule across compartments, a transport between
compartments that are not adjacent, and a transport gated by an undeclared
protein.

**How much of the layer is guesswork.** In the OXPHOS program, 419 facts:
90.2% curated or experimental, 8.8% inferred, 1.0% predicted. Of the 12
structure constructs in the two committed programs, 5 are experimental
(compartment volumes), 4 curated and **3 inferred — the three transport
capacities, which are not measured**; two time-scale parameters are inferred
as well. Over the whole MitoCarta program the import signal is curated from
UniProt for 549 proteins, predicted by TargetP for 156, and **inferred for
418** ("no presequence annotated, so an internal signal"), which is the
honest cost of covering the whole inventory. The runtime carries this to the
output: every located species reports the weakest link on the chain that put
it there, and for an imported protein that link is the transport at 0.3.

**Decided while implementing, for Albert to confirm or overturn.**

1. **A transport may cross one membrane.** Band 3 moves bicarbonate from
   outside to the cytosol, and the plasma membrane is a compartment between
   them, so "adjacent" means parent-child *or* both neighbours of one membrane.
2. **Assembly consumes its subunits, one of each.** Without it the red cell
   made forty times more haemoglobin than it had chains — an economy layer
   that produces matter from nothing is exactly the failure this stage exists
   to prevent. Real stoichiometry (α2β2) waits for stage 2's costs, so a
   complex count is currently the chain count.
3. **`initial` on a protein** is the bootstrap state S0 the architecture
   already requires; the red blood cell cannot be written without it.
4. **Amounts, not concentrations** (open decision 2 stands): volumes are
   declared and recorded but no rate reads them yet.
5. **What stage 1 deliberately does not model:** mitochondrial translation
   still works with import closed, because ribosomes are not yet a resource;
   in a cell the mitoribosome is imported. That dependence is stage 2's, and
   the program says so in its own evidence field.

## 10. Open decisions for Albert

1. **Opt-in or mandatory.** Located rules apply to programs that declare a
   compartment. Should v0.5 make locations mandatory everywhere?
2. **Amounts or concentrations** (stage 2 is held on this one: pools and costs would bake it in). Stage 1 keeps amounts per compartment, so
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
8. **~~When `fates: first` becomes the default.~~ Resolved 2026-09-14:** area E
   migrated (`0e93b5b`) and the default is now `first`. The measurement behind it:
   the worm's 1,439 fates and its terminal score are unchanged, its 45 ambiguous
   decision points fall to 0 and 620 fewer decisions fire, and every other
   committed organism program is identical under both modes, experiments included.
   `fates: last` remains available as the documented legacy mode.
9. **Commitment and competence as refusals or as rates.** As specified, a
   locked cell refuses a fate outright; the published plasticity data could
   also be read as a declining probability. The perturbation series decides —
   **and it has now been run (§7.2a), without settling this one.** A refusal
   reproduces all five arms at the resolution we score them (fates at 800 min),
   so nothing in the series demands a rate. But Fukushige & Krause describe the
   response as *declining rapidly over the subsequent hour* rather than stopping
   at an edge, and a step function cannot be that. What would decide it is their
   per-stage conversion frequencies, scored against the same program with the
   window as a probability instead of a gate; that is a run, not an opinion, and
   it needs the paper's own tables rather than its abstract. Until then the
   implementation is a refusal and this document says the edge is sharper than
   the biology.

10. **Order as a construct: position, direction and time.** The language can
    say that a gene is read, not in which order a cluster is read. The
    known-locus benchmark made this concrete: HOXD defeats the node model and
    its colinear activation order is untestable with any layer we hold,
    because nothing in BioLang can state an order. Proposed shape, for Albert:
    an `order` clause on a `domain` naming the axis (genomic position,
    direction of opening, time of activation) and the measured sequence, with
    the falsifier a published activation series (for HOXD, the limb and trunk
    colinearity data). Nothing is implemented; it should not bake node
    identity in, since area B's boundary caller may move every node.

**Which of these can be settled by measurement rather than preference.**

| Decision | Settled by | How |
|---|---|---|
| 1. opt-in or mandatory locations | preference | a language-design choice; no experiment decides it |
| 2. amounts or concentrations | **partly measurable** | run stage 2's burden gate both ways on cells of different volume; if only one reproduces the published burden scaling, it decides. Held for Albert until then |
| 3. membranes as nodes or edges | preference, lightly constrained | both passed stage 1; the band-3 transport forced the "across one membrane" rule, which either shape can express |
| 4. whose cost accounting | **measurable** | the ATP budget gate (Buttgereit & Brand's hierarchy, Lynch & Marinov's totals): the accounting that lands in the published order of magnitude wins |
| 5. default allocation | **measurable** | gate (b), absolute abundance against PaxDb, run under each policy |
| 6. protein turnover set | **measurable** | the dilution gate run with the human and the mouse sets; if they disagree beyond their own spread, the human set is required |
| 7. which homeostats first | data availability, not preference | the pair whose failure direction is best documented; Na⁺/K⁺ and pH stand |
| 8. `fates: first` default | **settled by measurement** (resolved above) | |
| 9. commitment as refusal or rate | **measured, and still open** | the series ran and passed as a refusal (§7.2a); the published decline over the fourth hour is what a rate would fit, and separating them needs the per-stage frequencies |
| 10. order as a construct | measurable once built | a published activation series |

Decision 9's instrument is built and run (§7.2a, 2026-09-14): the series is
`data/organisms/celegans/plasticity.bio`, `bio test` runs it, and removing
either construct breaks an arm. The next engine work that needs no preference
from Albert is therefore decision 10's construct sketch and the two pieces
§7.2a leaves missing — the integrated reads as a language read rather than a
precomputed lookup, and a `commitment` that a precursor can hold before it
differentiates; 4, 5 and 6 need stage 2 or 4 first, and 2 holds stage 2.
