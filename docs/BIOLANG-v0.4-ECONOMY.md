# BioLang v0.4: structure, economy and control (specification)

**Status on 2026-10-10.** Stage 1 (structure) is built and passed its gates on 2026-09-14 (§9.1).
Stage 2 (economy) was implemented on 2026-09-21: gate (a) passed on three clauses of four, gate (b)
failed, and the pool layer is optional (§9.2). Stage 3 (energy and mitochondria) was implemented on
2026-10-10: all ten registered clauses are met, and of the gate as §6 words it the magnitude, the split
and the red cell's ranges are not askable and the hierarchy is reproduced as an input (§9.5, §9.6).
Stage 4 (division and dilution) was implemented on 2026-09-28: the mechanism passes all eight
registered clauses and the gate as §8 words it is refuted (§9.4). `competence`, `commitment`, `order`
and `regime` are built; `homeostat` and `role` are not, and no engine code applies §2's 0.3 cap. The
current grammar, [BIOLANG-GRAMMAR.md](BIOLANG-GRAMMAR.md), has 27 block kinds plus `module` and
`import`: v0.4's 26 and v0.5's `reader` (2026-10-10, [BIOLANG-v0.5-READER.md](BIOLANG-v0.5-READER.md)),
which gates rules by chromatin state at rule selection; BioIR is still versioned 0.4.

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
  *(2026-10-10: no code in `genomeos/lang`, `genomeos/ir` or `genomeos/runtime` applies this cap or
  the `inferred` relabelling; it holds only where a program's author writes the confidence.)*
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
compartment RedCytosol    { parent: PlasmaMembrane; volume: 1.0; absolute_volume: 94 fL
                            evidence: experimental "Diez-Silva et al. 2010, MRS Bull 35:382"; confidence: 0.7 }
```

| property | form | meaning |
|---|---|---|
| `parent` | `Id` | the compartment that contains this one; exactly one root |
| `membrane` | `yes \| no` | a membrane faces its parent and its children at once |
| `volume` | `number` | fraction of the cell's volume (recorded in stage 1; used from stage 2) |
| `absolute_volume` | `number fL\|pL\|nL\|uL\|mL\|L\|um3 \| unknown` | this compartment's own volume, beside the fraction; what a concentration divides by |
| `genome` | `chrM, ... \| nuclear` | chromosomes read here; `nuclear` = every chromosome except chrM/MT |
| `translation` | `yes \| no` | ribosomes are present (cytosolic or mitochondrial) |
| `copies` | `integer` | copies per cell (recorded in stage 1; heteroplasmy in stage 3) |
| `heteroplasmy` | `m.POS REF>ALT = fraction, ...` | stage 3: the fraction of the copies carrying each mtDNA variant; needs `genome: chrM` (§6) |

**Absolute volumes and concentrations** (§10 decision 2, resolved 2026-09-19; implemented
2026-09-19). `volume` keeps its meaning and its name: a fraction of the cell, which is what the
containment tree needs and what the committed programs state. `absolute_volume` is a second,
independent fact — this compartment's own volume — and it is the one a threshold can be divided by,
because a fraction of a cell cannot tell a big cell from a small one. It is stated per compartment
and never derived: nothing multiplies a fraction by a cell volume on a program's behalf.

The accepted units are litres and its SI prefixes down to `aL`, plus `um3`, which is exactly `fL`.
Both spellings are accepted because the literature quotes a red cell in femtolitres and most cells in
cubic micrometres and means the same unit, so a citation is transcribed rather than converted by
hand. The unit is **required**: a bare number already means a fraction on this block, so
`absolute_volume: 94` is a compile error rather than an ambiguity. `absolute_volume: unknown` is a
statement a program may make, the same one `param x = unknown` makes, and it is what
`human.oxphos` says — its fractions are a hepatocyte's but it declares no cell type, so no absolute
volume can be cited for it.

What the field buys is that a rule's `threshold` may now be a concentration: `threshold: 100 nM`
instead of `threshold: 5661`. A molar unit on the threshold is what says which it is; the located
runtime converts it once, at each compartment the rule acts in, as `c · V · N_A` (Avogadro is exact
by the 2019 SI definition of the mole), and everything downstream stays an amount. Three refusals
come with it, all of them the same refusal — a concentration without a volume is not a number:

- a rule with a concentration threshold acting in a compartment with no `absolute_volume` is a
  **compile error** (the seventh modelling error refused before a run, §9.1);
- the same rule in the unlocated network runtime, which has no places at all, is refused;
- a concentration threshold under `units: au` is refused: arbitrary units have no copy number, so
  the regime must declare `units: copies`.

**What it does not do.** No rate reads a volume yet, and stage 2 is not built here (§5)
*(2026-10-10: stage 2 was implemented on 2026-09-21, §9.2)*. The field
makes the comparison *expressible*, which is what §9's registered gate needed to be runnable at all;
it does not decide which unit a threshold should be stated in. The demonstration that the two units
are the same statement in one cell and different statements across two is
`data/demo/concentration_threshold.bio` and `tests/test_concentration.py`.

**Semantics.** Compartments form a containment tree. Two compartments are
*adjacent* when one is the other's parent, or when both sit either side of one
membrane compartment, which a transport then crosses. This is the same shape as a
process-bigraph place graph, so a located program maps onto a Composite whose
stores nest the same way (`runtime/compose.py`).

**BioIR.** `Compartment(Entity)`: `parent`, `membrane`, `volume`, `genome`,
`translation`, `copies`, `absolute_volume_fl`, `heteroplasmy` (stage 3); `kind = "compartment"`. `Rule.threshold_unit` is empty
for an amount and names a molar unit for a concentration.

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

## 5. Economy (stage 2, implemented 2026-09-21; gate (a) passed, gate (b) failed)

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
*(2026-10-10: the type also has `regenerates_from`, which keeps `regenerates: from Glycolysis,
OXPHOS`; no runtime reads it, since the processes it names belong to stage 3, which is not
implemented.)* *(Later on 2026-10-10: stage 3 reads it. `runtime/metabolism.py` runs Glycolysis, TCA
and OXPHOS for the one pool that names them, and refuses any other name, §6.)*

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
*(2026-10-10: `optimise` compiles, and the economy runtime refuses to run it with an error that the
runtime solves no objective, `runtime/economy.py`; the other three run.)*

*(2026-10-10, SUPERSEDING the line above the same day, lane-deadkeys, ROADMAP section 5 item 13 A:
`optimise` no longer compiles. The CHECK refuses it, `lang/parser.py`'s `_OPTIMISE_REFUSAL`, in an
`allocation` block and in a `regime`'s `allocation:`, naming the policy and the four that can be run,
because a program that cannot run must not pass `bio check`. The economy runtime's refusal stays for
an IR built without this parser. The `growth_law` example above is therefore no longer a program; it
is kept as what the key would mean if an engine ever honoured it. Whether `optimise` stays in the
language is Albert's decision, so it is still listed and still parsed as far as its objective.)*

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

**Both gates have now run. (a) passed on its three askable clauses; (b) failed,
and the trap stated above turned out to be understated.** §9.2 has the numbers.
The short form: the sentence "cannot change a rank correlation" is true of the
absolute-error test as well, because `proportional` — and, as implemented,
`competitive` too — hands every demander of a pool the same scalar, and a
baseline that fits a constant absorbs it exactly. **So the pool layer stays
optional, as this paragraph agreed in advance.** It is expressible, it is
charged, it reproduces a burden of the documented shape, and at proteome scale
it does not improve a single prediction.

*(2026-10-10: a fourth runnable policy, `saturable`, gives each demander its own factor
`K / (K + d)` in its own demand `d`, with K the program's `half_saturation`; §9.7 registered gate (b)
for it and §9.8 has the run: on held-out genes its log10 MAE is 0.949823 against the baseline's
0.949969, no improvement detected at the 0.05 bar, and the pool layer stays optional.)*

## 6. Energy and mitochondria (stage 3, implemented 2026-10-10; ten registered clauses met, the gate as worded partly not askable)

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

**As built (2026-10-10).** The paragraphs above are left as written. The language gains one key,
`heteroplasmy: m.3243A>G = 0.7` on a compartment whose genome reads chrM (a fraction of its `copies`;
fractions of two variants add and may not pass 1, since nothing here knows their linkage). The
module is `bio.std.core_metabolism` (`genomeos/std/core_metabolism.bio`): thirteen `param`s, the
textbook stoichiometry of glycolysis and of pyruvate dehydrogenase with the TCA cycle, Hinkle's P/O
ratios at both ends (2.5 and 1.5; 2.3 and 1.4), and half an O₂ per electron pair. The engine holds no
such number. `genomeos/runtime/metabolism.py` refills the one pool that names `from Glycolysis, TCA,
OXPHOS` per hour, reading the cell's own facts from the program: the `Glucose` and `O2` pools'
`regenerates` rates, `param oxphos_capacity_per_copy`, and `param bpg_shunt` for the red cell. Supply
follows demand, pyruvate is oxidised first and the rest becomes lactate; every reading is computed at
both ends of the spread (32 and 28 ATP per glucose); heteroplasmy removes the capacity of the copies
carrying the variant. What the cell cannot make is divided among its consumers by the program's own
`allocation`, through stage 2's `Economy`, which now takes a per-pool capacity in place of the
standing stock. The registration is §9.5 and the measurement §9.6.

## 7. Control (specified here, implemented after stage 2)

*(2026-10-10: stage 2 was implemented on 2026-09-21. `homeostat` (§7.1) and `role` (§7.2) are still
not built; §7.2a, §7.4, §7.5 and §7.6 record what was.)*

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
   perturbation had nothing to inherit: 12 of 610 cells in the last arm (8 of
   them in the E lineage) took the forced fate at birth, because the program had
   no determination before differentiation — a precursor was committed to
   nothing. **Closed 2026-09-15, and not by a new construct.** Once area E's fate
   rules read the integral directly, those parents reach a terminal fate *before*
   they divide, commit, and the lock is inherited: **0 of the 12 now take the
   forced fate**, under either `recheck` mode, and the five published arms are
   unchanged.

   *Corrected 2026-09-15, after area E checked the attribution.* This document
   first credited their new `commitment terminal_fate`. It is not what does the
   work, and a one-block ablation shows why not — **either block alone is
   enough**, because both establish on a terminal type and the daughters inherit
   whichever locked first:

   | ablation | of the 12 born after the induction, forced to muscle |
   |---|---|
   | as written | 0 |
   | area E's `terminal_fate` stripped | 0 |
   | this program's own `differentiated` stripped | 0 |
   | **both stripped** | **12** |

   So what closed the hole is the parents reaching a terminal fate before
   dividing; the lock that then protects the daughters was already in
   `plasticity.bio` and had been all along. The coarse ablation that cleared
   every commitment at once could not see that, which is why the first version of
   this paragraph named the wrong cause.

   The *shape* of the hole is still real, and worth stating because a future
   program will meet it: a cell that runs a programme for a long time **without
   differentiating** has nothing to pass on. A program can now say so itself,
   with no new construct, because `establish` on an integrated read is this
   section's own example and compiles since the reads landed:

   ```
   commitment gut_programme {
     programme: Intestine                      # what it commits to, before it becomes it
     establish: ELT-2.exposure(lineage) >= 30  # a programme running, which is what the runtime can see
     inherit: daughters; release: never
   }
   ```

   So the two questions this raised answer themselves. *What establishes a
   commitment in a cell that has not differentiated* is a **programme running**,
   and the only thing the runtime can observe of that is a sustained read — the
   same quantity the fate rule itself reads, at a lower threshold. And it is
   **the same construct with an earlier trigger**, not a weaker one: a
   determined precursor differs from a terminally differentiated cell in *when*
   it locks, not in how firmly, which is why an early forcing beats it (no
   exposure has accumulated yet) and a late one does not. A separate construct
   with weaker claims would have been a second way of saying the same thing, and
   it is not built. `tests/test_plasticity.py` pins both: the 12 cells, and a
   precursor that never differentiates protecting its daughters.
4. **Clauses specified here but not gated are compile errors, not no-ops:**
   `maintain`, `excludes`, `hysteresis` and any `release` other than `never`.
   The `.exposure(window)` / `.mean(window)` reads were in that list when the
   series was written and were implemented on 2026-09-15 (above), so
   `establish: ELT-2.exposure(lineage) >= 0.8` compiles and runs; a read that
   does not name its window is still a compile error.
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

**Implemented 2026-09-15.** The Body integrates what a cell actually carried:
a factor counts 1 while the cell carries it, a species with a per-cell network
level counts its level, and the integral is banked at every decision point (which
is where factors change) and clamped when the cell divides or dies. `exposure` is
in minutes; `mean` divides by the window, so a plain factor reads as the fraction
of the window it was carried (0 to 1) and a network species as its mean level.
`cell` is the cell's own life, `lineage` the path from the first cell, carried to
the daughters at each division. A program that asks for no integrated read pays
nothing: the reads a program names are collected once at start-up, and every
committed organism program and every experiment arm is identical before and after
this change.

Two things it deliberately does not do. **It is absolute, not normalised:** area
E's generated `_integrated` factors call a factor present at "20% of its largest
path exposure over the terminal cells", and that maximum is a statistic of the
finished run, which no cell can read while it is deciding. A program states a
threshold in minutes (or a fraction, with `mean`). **It integrates carriage, not
concentration:** without a `cell_network` there is no level to integrate, so
`exposure` counts time carried; with one, the level is used and the answer changes
meaning. A program should say which it means in its evidence.

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

**A later decision point is not an appeal (2026-09-15).** Area E found the hole
by giving the worm a `cell_network` cadence: every step re-decides every cell,
and a rule that had lost the precedence contest simply won it later, so the
Sulston fate score fell with the cadence and for no other reason. Reproduced on
`embryo_factors.bio` to 800 min: **484 of 555 with a 6-minute network against
501 without it**, and the 17 cells are exactly rules taking a fate back. The
rule the code already applied inside one decision point — *a decision that
already applied before the change is a competitor, not a successor* — now
carries across decision points: a decision of **lower** precedence that **could
already have applied when the fate was settled** is refused and counted as
`overruled_fates`. A decision whose guard was false then is a new reading and
still applies, which is how a contact signal changes a cell's mind; equal
precedence still applies, so nothing is frozen. With the fix the score is 501
of 555 with the cadence and without it, every committed organism program and
every experiment arm is byte-identical, and 68 attempts are counted instead of
being taken silently.

**Two decisions may share an id, but not silently (2026-09-15).** Candidate
lists were indexed by decision id, so a decision *without* a `cell` clause that
shared an id with one that had it was dropped from every candidate list and
never fired — nothing said so, and it cost area E an hour. Indexing is by
identity now. Shared ids stay legal, because the worm states mechanism before
the generated lookup and both are called `div_EMS`, but they are not free: a
cell's `fired` list records ids, so two decisions with one name cannot be told
apart in it. Every run therefore reports `duplicate_decision_ids` beside
`ambiguous_fates` (the worm's five founder divisions show up there), and a new
rule should be given a new name.

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
- **A death re-reads contacts, as a birth does (2026-09-15).** `_add` makes a
  newborn's neighbours decide again; a death only freed the grid site, so a
  receiver went on reading a ligand from a cell that was no longer there. The
  worm kills 110 cells, so this is not a corner case. The neighbours are taken
  before the site is freed and each re-decides, which is the same path a birth
  takes; with per-cell networks the next step already re-read everyone, and that
  is still the condition under which this is skipped.
- **`asymmetric` segregates a factor and never creates one (2026-09-15).**
  `asymmetric: X -> D` used to write X into the keeper even when the mother had
  none, which is a fate appearing from nothing — the failure the whole evidence
  discipline exists to prevent. A keeper now keeps only what its mother had.

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

### 7.5 Looking again: when a cell decides, and why not on a clock [engine]

A cell decides when it is born, when a signal reaches it, when a stage starts,
when a factor is forced on it, and when it finishes dividing. Nothing else makes
a cell look again — so a rule whose guard becomes true *between* those moments
never fires, and the language has no way to say "look again later".

Area E hit this the moment their fate rules read a growing integral
(`ELT-2.exposure(lineage) >= 15`, §7.2a): the threshold is crossed at some
instant of the cell's own life, and unless an unrelated event happens then, the
cell never notices. The only way to be re-read today is to declare
`cell_network: 6 min`, which exists to step per-cell networks (§7.4) and
re-decides every cell as a side effect. Used as a clock it is measurably wrong:
the worm scored **386 of 555 against 522** at a 6-minute cadence, because a
clock gives *every* rule extra chances, not only the one whose quantity arrived.
Area E recovered it with `commitment`. Stepping a network must not change which
fates are taken.

**Two shapes, and the one chosen.**

- *A declared interval* (`recheck: 6 min` on a decision or a cell type). Simple
  and obviously deterministic, but it is a clock with a free parameter, and the
  measurement above is the argument against it: the answer would depend on an
  interval nobody can derive, and every rule pays for one rule's need. A cell
  does not consult a clock.
- *A crossing* — the cell decides again at the instant the quantity a guard
  names reaches the threshold that guard states. **Chosen.** It has no free
  parameter, it is what the biology does (a cell responds when a quantity
  arrives), and it is *exactly* computable rather than merely approximable: the
  runtime already integrates each read as a rate held constant between decision
  points (§7.2a), so for `F.exposure(w) >= T` the crossing is
  `now + (T - E)/r` for the current exposure `E` and rate `r`, and for
  `F.mean(w) >= T` it is `now + (T·span - E)/(r - T)` when `r > T`. Both are
  closed forms, so the re-decision is scheduled, not searched for. Determinism
  is not weakened: the crossing goes into the same event queue as every other
  event and draws no randomness, so a replicate under one seed is unchanged.

**Semantics.** Declared in the regime, because it is an execution scheme and
"a result without its regime is not a result" (§8):

```
regime worm { recheck: crossings }      # crossings (default) | none
```

- Under `crossings`, whenever a cell finishes deciding, the runtime computes the
  earliest future instant at which any *unmet* threshold of an integrated read
  named by a decision that could apply to it would be reached at the rate now in
  force, and schedules one re-decision there. A rate of zero schedules nothing;
  a factor that arrives or leaves already re-decides the cell, which recomputes
  the crossing, so a cell holds at most one pending re-decision.
- A crossing is a *re-reading*, not a new right: precedence (§7.3) applies
  unchanged, so a rule that lost when the fate was settled cannot win at a
  crossing, and `commitment` still refuses what it refuses.
- Under `none` nothing is scheduled and a program must arrange its own
  re-decisions, which is the behaviour before this section existed.
- Every run reports `rechecks` (crossings taken) beside `ambiguous_fates`.

**BioIR.** `Regime.recheck: str`; no new entity, and no change to `Decision` —
the guard already names the quantity and the threshold, which is the whole
input. A program that names no integrated read schedules nothing and is
provably unchanged.

**Falsifying measurement.** *If re-decision is a first-class thing, stepping a
network must stop changing fates as a side effect.* The test is area E's worm
with fate rules on integrated reads, scored against Sulston: the terminal fates
must be **identical with no network and with a 6-minute network**, and identical
to what the program scores today with `commitment` doing the work. Falsified if
the cadence still moves a fate, or if `recheck: crossings` and `recheck: none`
differ on any program that names no integrated read.

**Implemented and measured 2026-09-15; it passes, and it refutes part of the
argument above.** The worm scores **522 of 555 in all four combinations** of
cadence (none, 6 min) and `recheck` (none, crossings), with **not one cell
differing** in type or terminal name, and all 14 committed organism programs are
identical under `crossings` and `none` — the only difference anywhere is 250
further decisions firing on the three worm programs that name a read, which
change nothing. So the falsifier holds and the default is `crossings`.

But the sentence "a lock compensating for a scheduler" was wrong, and the
measurement says so. Strip `commitment terminal_fate` from the worm and the
score depends on when the cell looks, whatever the scheduler:

| `embryo_factors.bio` | no network | 6-minute network |
|---|---|---|
| as written (`commitment` declared) | 522 | 522 |
| `commitment` stripped, `recheck: none` | 522 | **386** |
| `commitment` stripped, `recheck: crossings` | **412** | 386 |

A crossing is *better* than a clock (412 against 386) and it is principled where
a clock has a free parameter, but it does not restore 522: with a monotonically
growing integral **every later look is a genuinely new reading**, so any extra
look moves fates unless the fate is locked. Area E's conclusion stands exactly as
they wrote it — *a fate written on a growing integral is only a fate if the cell
cannot take it back* — and `commitment` is doing real work, not covering for the
scheduler. What §7.5 buys is that the *when* is now declared, parameter-free and
exact, instead of being smuggled in by a construct that means something else.

The second lesson is in the 412 itself: these rules were learned against the
decision points the program already had (births, signals, stages), so looking
more often is not automatically truer. A program that changes `recheck` is
changing what its rules were fitted to, and should say so.

**Decided while implementing.**

1. **A crossing that has effectively happened has happened.** A wake-up computed
   to land exactly on its threshold can miss it by a float's width and reschedule
   itself for ever — the worm hung on the first run. A threshold within 1e-9 is
   reached, and a wake-up sooner than 1e-6 minutes is now.
2. **Revising nothing is not a revision.** A decision that lands on the type the
   cell already holds was counted in `revised_fates` and took the cell's terminal
   name away, so 250 of 555 terminal cells kept the right type and stopped
   answering to their names — the §7.4 anchor-cell bug, reached from the other
   side. Fixed, and it changes nothing under `recheck: none`, where
   `revised_fates` was already 0 on every committed program.
3. **The mean over a window of no length is the value in force.** At the instant
   a cell is born `mean(cell)` is 0/0; its limit is the rate now in force, so a
   cell reading it at birth is taking the instantaneous read, which is the
   honest answer and needs no wake-up.

### 7.6 `order`: a sequence the runtime can be wrong about [engine] (specified, §10 decision 10)

*(2026-10-10: `order` was built on 2026-09-15 in `144f399`; its gate results are in §9.)*

The language can say *when* a thing happens (`timer`, `stage`, a crossing) and
*where* it is (`locus`, `compartment`), but not that A happens before B. The
known-locus benchmark hit this on HOXD, where position on the chromosome *is*
the order of execution — the cluster is read 3′ to 5′, and the benchmark's
verdict was that colinearity is untestable with any layer the project holds,
because nothing can write it down.

**First, whether it earns its place, because mostly it does not.** Wherever a
program knows the times, order is already expressible and needs no construct:

```
assert: lineage AB at 20 min >= 1        # AB exists by 20 minutes
assert: lineage EMS at 20 min = 0        # and EMS does not: AB is born first
```

That is the worm's founder order, written with v0.3 and checked today. So an
`order` block earns its place in exactly one case: **a sequence that is known
without its times.** There the alternatives are all worse — invent times and
assert them (fabricating the facts §2 exists to prevent), or write a relay of
rules (asserting that each member switches the next on, which for a HOX cluster
is false: colinearity is progressive chromatin opening, not a cascade), or write
nothing, which is what the benchmark did. Two published cases have that shape:
HOXD limb and trunk colinearity, and the β-globin switch, where HBE1, HBG1/2 and
HBB come on in the order of their distance from the control region. Both are *the
order is in the layout*, not *the order is in a clock*.

**The construct.**

```
order hoxd_colinear {
  members: HOXD1, HOXD4, HOXD3, HOXD8, HOXD9, HOXD10, HOXD11, HOXD12, HOXD13
  axis: position                   # position | time
  direction: decreasing            # position only: along the chromosome
  evidence: experimental "..."; confidence: 0.8
}
order founders { members: AB|P1, EMS|P2, MS|E|C|P3, D|P4; axis: time; observe: birth }
```

**A step may name members a source does not order** (`B|C`). Both falsifiers
needed it, which is why it is here rather than in a later version: the worm's
reference puts MS, E, C and P3 at one instant, and a construct that took only a
total sequence would have made the author invent an order between them, which is
what §2 forbids. An inversion therefore needs both a stated order and a real gap
in the run.

**Semantics, and the one rule that keeps it from being a second `stage`.** An
`order` **never drives anything**. It sets no time, fires no rule and changes no
level; it is a claim *about* a run, and its whole content is that it can be
contradicted. `stage` and `timer` cause things to happen; `order` states what
must be true of the things that happen, and is checked.

- **`axis: position` is checked at compile time**, against the coordinates the
  program already holds: the declared sequence must be monotonic in the members'
  loci, in `direction`. A program whose sequence contradicts its own genome is a
  compile error, not a warning — the same treatment as a chrM gene declared
  nuclear (§9.1). This is the check with teeth today, because the loci are
  measured and we hold them.
- **`axis: time` is checked against a run.** The runtime records when each member
  first crossed `threshold` (for a network) or when its cell was born or took its
  fate (for a Body), and reports `order_taken`, the number of **inversions**
  against the declared sequence, and the members that never happened. An `order`
  with inversions does not stop the run: it is reported, as `ambiguous_fates` is,
  because a program that predicts the wrong sequence is a result, not a crash.
- Members that a run never activates are reported separately from members that
  activate out of order, since "did not happen" and "happened too early" are
  different failures and a single score would hide both.
- `direction: opening` appears in §10's sketch — the direction chromatin opens —
  and is **refused**: we hold no measurement of it, and a clause the runtime
  cannot check would be the comment this construct exists not to be.
- It must not bake node identity in (§10): `members` are genes, never nodes, so
  area B's boundary caller can move every node without touching an `order`.

**BioIR.** `Order(Entity)`: `members`, `axis`, `direction`, `observe`,
`threshold`; `kind = "order"`. No change to any existing type. **[engine]** —
language and runtime are Apache 2.0; the loci and any activation series that test
it are read by **[app]** scripts, as §3 requires.

**Falsifying measurement, and the part we cannot do.** The honest position, stated
before implementing:

- *What we can falsify today:* the **position** claim. We hold the HOXD cluster's
  coordinates (`data/cache/gencode_genes.tsv`: HOXD13 at chr2:176,092,689 through
  HOXD1 at chr2:176,188,658), so a program can state the colinear sequence and the
  runtime can refuse it if it contradicts the genome. Falsified if a correct
  sequence is refused, or a scrambled one accepted.
- *What we can falsify today in a run:* the **worm's founder birth order**, where
  the reference gives every cell's birth time independently of the program. An
  `order` over AB, EMS, P2, MS, E, C, D, P4 must report zero inversions on the
  wild type, and must report the inversion when a timer is perturbed to reverse a
  pair. Falsified if a reversed run is reported as in order.
- **What we cannot do, and what it would take.** The motivating case is not
  testable: *the project holds no expression time course for a human HOX cluster*,
  so the claim "HOXD activates in the order it is laid out" can be checked for
  consistency with the genome but never **reproduced**. What it would need is a
  staged time course over a differentiating limb bud or an embryonic stem cell
  line — per-gene expression for the nine members at three or more time points,
  which is an external dataset, not a run. Until that exists, the `time` axis of
  this construct is gated only by the worm, and this document says so rather than
  claiming HOXD as a pass.

**Built and measured 2026-09-15**, in exactly those two places.
`data/demo/hoxd_order.bio` is the benchmark's untestable sentence written down:
the cluster from GENCODE v50, refused if two members are swapped. The worm's
founders report zero inversions and none missing, a reversed sequence is reported
with the order the run actually took, and a member that never happens is counted
apart from one that happens early. `bio test` checks every order in a committed
program, and a test pins the rule that keeps this from being a second `stage`:
adding an order to the worm changes no field of the run's summary.

**Two things writing it found, both of which changed the specification above.**

1. **On GENCODE v50 the HOXD3 gene record spans HOXD4**
   (chr2:176,136,077-176,174,483 against 176,151,544-176,156,078), so the
   canonical textbook sequence HOXD1, HOXD3, HOXD4 **is not monotonic in
   transcription start** and the first version of the check refused it. That is a
   question about HOXD3's annotated start rather than about colinearity, so the
   refusal now names the overlap and suggests the step form — and the demo states
   the order this annotation supports, with the discrepancy recorded in the file.
2. **The steps came from the falsifiers, not from the design.** Neither of the two
   things we can actually check is a total order, and pretending otherwise would
   have required inventing facts in both axes at once: the worm's four tied
   founders, and HOXD3 against HOXD4.

*Still not implemented:* `observe: activation` (a network's first threshold
crossing), refused rather than accepted, because the dataset that would gate it
is the HOX time course we do not hold.

## 8. Execution regime (declared per run)

```
regime default {
  treatment: auto          # continuous | stochastic | auto
  threshold: 50            # auto: species under this many copies are treated stochastically
  units: copies            # stochastic treatment is refused for arbitrary units
  update: continuous       # continuous (ODE) | synchronous | asynchronous | event
  allocation: competitive  # default policy for pools and transports with none named
  fates: first            # Body: one fate per decision point by precedence (last = legacy)
  recheck: crossings      # Body: a cell decides again when a read it names reaches a threshold (§7.5)
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

**Stage 4, division and dilution (implemented 2026-09-28; the mechanism passes
its eight registered clauses and the gate as worded below is refuted).** A
`divide` event partitions contents (binomially for low copy numbers) instead of
duplicating them, guarded by a size or resource checkpoint, so growth costs
something. Gate: over several divisions dilution must separate a stable protein
from a short-lived one. `bio.std.human_turnover` holds cell lifespans (Sender &
Milo 2021), not protein half-lives; the gate needs a protein turnover set (for
example Mathieson et al. 2018, Nat Commun 9:689), which is **open** below.

The paragraph above is left as it was written, the way §9.1 item 4 leaves its
own superseded decision standing. What it asks for is built — `partition:
duplicate | contents | binomial` on an `event`, `genomeos/runtime/division.py`,
and the checkpoint as the event's own `when` clause — and the sentence that
states its gate is **wrong**, for a reason written down before the instrument
existed. Partitioning multiplies every species in the cell by the same one-half,
so the division term cancels out of the stable:short ratio and **dilution
separates nothing**; in a steady state it compresses the separation 8.80-fold.
§9.3 registers that and §9.4 measures it. What dilution does do on its own is
remove a protein the program gives no half-life to, which turnover cannot touch,
and that is the clause the stage passes. The turnover set stays open, and no
clause here needs it.

## 9. Gates and status

| Stage | Constructs | Gate | Status |
|---|---|---|---|
| 1 | compartment, location, signals, transport, regime record | chrM, MitoCarta import, rho0, red blood cell | **passed 2026-09-14** (below) |
| 2 | pool, cost, allocation | burden; absolute abundance vs PaxDb | not started; **no longer blocked** — decision 2 was resolved 2026-09-19 and the absolute volume it needed is implemented (§4.1), so the gate's concentration arm is now expressible |
| 2 (superseding the row above, 2026-09-21) | pool, cost, allocation | burden; absolute abundance vs PaxDb | **implemented; gate (a) passed on three clauses of four, gate (b) failed** (§9.2). The burden reproduces the documented shape and vanishes when the pool is switched off; the abundance gate improves nothing, for a reason derived from the runtime and committed before the data were fetched. §5.3's consequence stands: the pool layer is **optional**. The row above is left standing rather than rewritten, the way §9.1 item 4 leaves its own superseded decision standing |
| 3 | core metabolism, mitochondrial copies, heteroplasmy | ATP budget; oxygen and glucose dependence; red blood cell glycolysis | not started |
| 3 (superseding the row above, 2026-10-10) | core metabolism (`bio.std.core_metabolism`), mitochondrial copies, `heteroplasmy` | ATP budget; oxygen and glucose dependence; red blood cell glycolysis | **implemented; all ten registered clauses met, and the gate as §6 words it is partly not askable** (§9.5 registered, §9.6 measured). Not askable: the budget's magnitude (so §10 decision 4 stays open), its split, the red cell's ATP and 2,3-BPG ranges. Reproduced as an input, not predicted: Buttgereit & Brand's hierarchy, which appears under the `priority` the program writes and cannot appear under `proportional` or `competitive`. By construction: the Pasteur shift. Measured: the yields 32 and 28, 6 O₂ per glucose, the red cell's 1.62, the biochemical threshold at the declared reserve, and stage 4's division paid once ATP is regenerated over a cycle |
| 4 | partitioning division, checkpoint | dilution vs protein turnover | **implemented 2026-09-28; the mechanism passes all eight registered clauses and the gate as §8 words it is refuted rather than passed** (§9.3 registered, §9.4 measured). Partitioning is gene-blind, exactly as stage 2's `_share` was: dilution contributes 0.0 of the 9.392136 decades that separate a stable protein from a short-lived one, and in a steady state it compresses the separation 8.80-fold. The stage's own claim — conservation, and the removal of a protein nothing degrades — holds: 2,000 of 2,000 binomial splits conserve exactly, and a protein with no declared half-life falls by exactly 2^-8 over eight divisions and not at all when contents are duplicated. `partition` defaults to `duplicate`, so no committed program changes. The protein turnover set (§10 decision 6) stays open and no clause needed it |
| control | homeostat, role | two homeostats hold and break correctly | specified only |
| fate | commitment, competence | a published perturbation series (Fukushige & Krause 2005; Yuzyuk et al. 2009) | **passed 2026-09-14** (§7.2a): five arms reproduced, and each construct fails an arm when removed |
| fate | integrated reads (`.exposure(window)`, `.mean(window)`) | area E's three readings on 555 terminal cells | **implemented 2026-09-15**, in absolute units; the gate is area E's to run, by writing their fate rules as declared reads instead of a precomputed lookup |
| fate precedence | `regime fates`, decision `priority` | area E's worm fates identical under explicit priorities | **implemented and now the default**: 1,439 of 1,439 identical, 45 ambiguous points to 0 |
| re-decision | `regime recheck`, crossings of an integrated read | stepping a network must stop changing fates as a side effect | **passed 2026-09-15** (§7.5): 522 of 555 in all four combinations of cadence and recheck, 0 cells differing, all 14 programs identical |
| order | `order` along position or time | a sequence that contradicts the genome is refused; the worm's founder births are reproduced and a reversed one reported | **position passed 2026-09-15** (§7.6, HOXD from GENCODE v50); **time passed on the worm**; the motivating claim, that HOXD *activates* in that order, is **not testable**: no human HOX time course |
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

Seven modelling errors are refused at compile time rather than run: a gene or a
protein with no location, a chrM gene declared nuclear, a nuclear gene whose
mRNA reaches no ribosome, a rule across compartments, a transport between
compartments that are not adjacent, a transport gated by an undeclared
protein, and (since 2026-09-19) a rule whose threshold is a concentration acting
in a compartment that declares no absolute volume.

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
4. **~~Amounts, not concentrations~~ (open decision 2 stands): volumes are
   declared and recorded but no rate reads them yet.** *Superseded 2026-09-19:* decision 2 is
   resolved and the absolute volume is implemented (§4.1). A threshold may now be stated as a
   concentration and is converted at the compartment where its rule acts. No *rate* reads a volume
   yet, which is still stage 2's, so the sentence above is true of rates and no longer true of
   thresholds. `human.erythrocyte` states `absolute_volume: 94 fL` (Diez-Silva et al. 2010, MRS
   Bulletin 35:382, "Healthy Red Blood Cell Structure": the cytosol "averages 94 um3 at 300
   mOsmol/kg"), and `human.oxphos` states `absolute_volume: unknown`, because its fractions are a
   hepatocyte's while the program declares no cell type — one of the two citations the decision
   priced exists, the other has no cell to be about. The two units are shown to agree in one cell and
   to disagree across two of different volume in `tests/test_concentration.py`: 20,000 molecules of a
   repressor with a threshold of 100 nM and a threshold of 5,661 molecules (the same statement at 94
   fL) give target levels 1.041 and 1.041 at 94 fL, and 1000.757 against 1.041 at 940 fL, a factor of
   961. That is the arm of the gate the language could not express.
5. **What stage 1 deliberately does not model:** mitochondrial translation
   still works with import closed, because ribosomes are not yet a resource;
   in a cell the mitoribosome is imported. That dependence is stage 2's, and
   the program says so in its own evidence field.
   *(2026-10-10: stage 2 made a ribosome pool expressible, but `human.oxphos` declares no pool, so
   this still holds for that program.)*

### 9.2 Stage 2 as measured (2026-09-21)

The language surface is `genomeos/lang/parser.py` and `genomeos/ir/model.py`; the
arithmetic is `genomeos/runtime/economy.py`; the two gates are
`genomeos/runtime/burden_gate.py` and `genomeos/runtime/abundance_gate.py`, each
of which is a **pre-registration committed before the instrument that runs it**.
The results are `data/results/burden_gate.json` and
`data/results/abundance_gate.json`.

**Gate (a), burden: passed on three clauses of four, and the fourth was
registered as unaskable before it was run.** Holding an unrelated gene's own
demand fixed at 100 and expressing the costly protein at seven registered levels
from 1e3 to 5e7, the unrelated gene's factor reads 1.0, 1.0, 1.0, 1.0, 0.50,
0.25, 0.05. It falls, it never rises, and the pools that go short are named at
every level where it moves. **The falsifier holds:** with every pool multiplied
by 1e6 and the same demands, every factor stays 1.0, so the reduction is not the
bug that looks exactly like the finding. The **magnitude** clause — agreement
with Ceroni 2015 or Frei 2020 in order of magnitude — is *not askable here*: this
project holds neither measured curve, the verdict prints "not askable" rather
than an agreement, and a test asserts that clause cannot contribute to a pass.
So the honest reading is narrow: the implementation produces a burden of the
documented **shape** where a shared pool is short, and the levels at which the
curve bends are a property of the demo program's round HeLa-scale capacities
rather than of a cell.

**Gate (b), absolute abundance: failed, and the failure is algebraic rather than
empirical.** PaxDb's integrated human whole-organism set (19,483 proteins) against
GTEx v8 median TPM in *Cells - Cultured fibroblasts* — the tissue declared before
the fetch — joined on 18,420 symbols, 15,159 analysed after dropping 2,423 with
zero TPM, 516 below the 10-copy floor and 322 with no residue count.

| Reading | Number |
|---|---|
| distinct factors the runtime gives 15,159 genes | **1**, under `proportional` and `competitive` alike, with the pools short and with them not |
| log10 MAE, one-to-one baseline with a fitted constant | 0.952005 |
| log10 MAE, allocation arm with a fitted constant | 0.952005 |
| improvement | **0.000000000** (bar to pass: 0.05, bootstrap 95% [0.0, 0.0]) |
| Spearman, both arms | 0.529265, identical — reported and **excluded** from the verdict |
| allocation arm with its constant *derived* rather than fitted | 1.551887, i.e. **worse** by 0.599881 |
| total protein per cell from the declared ribosome capacity, nothing fitted | **1.36e10**, against Milo 2013's 4e9–1.2e10 |
| OLS slope of log10 copies on log10 TPM | **0.515** [0.499, 0.531] |

The reason was derived from `Economy._share` and **committed before the transcript
arm was fetched** (`genomeos/runtime/abundance_gate.py`, `AMENDMENT`). `_share`
returns `capacity / wanted` to every demander of an oversubscribed pool — one
scalar with no gene index — and an entity's factor is the min over the pools it
draws on. Every protein-coding gene draws the same two, so the factor is one
global constant F and `log10 P_alloc = log10 P_base + log10 F`. The baseline fits
a constant, so the two arms are the same arm. In the run F was 1.754e-4, which
moves every prediction by 3.8 decades, and the MAE is unchanged to six decimals.
**§5.3's trap was understated:** it said the pool layer cannot move a *rank*
correlation, and the truth is that it cannot move the absolute-error test either.

Three findings the gate bought, none of them about allocation:

1. **At the demand a real transcriptome implies, §5.1's capacities are not even
   binding.** The pools had to be driven 1e4 times harder before anything was
   short. At proteome scale the pool layer currently does nothing at all.
2. **`competitive` is as gene-blind as `proportional`.** It returns
   `capacity / (capacity + wanted)` where `wanted` is the pool's *total* demand,
   so it differs from `proportional` only in the value of the global constant.
   Its own comment says each demander's share falls off; every demander's share
   is the same number. Three policies of the four are gene-blind; `priority` is
   the one that is not, and it needs an ordering over 19,000 genes that no
   measurement supplies. **This is §10 decision 5's answer, and it is a negative:
   gate (b) cannot choose a default allocation policy, because the policies it
   was meant to choose between make identical predictions.**
3. **The secondary is the one number the pool layer supplies and nothing fits,
   and it lands.** 5e6 ribosomes at 5.6 aa/s over an expression-weighted mean of
   491 residues is 5.7e4 chains/s; a 46 h median protein half-life
   (Schwanhäusser et al. 2011, Nature 473:337) makes that 1.36e10 proteins per
   cell, against Milo 2013's 4e9–1.2e10. It could have been wrong by decades and
   is high by 1.1×. That is a statement about §5.1's capacity, not about
   allocation.

**What would make the primary able to move**, named in the registration before
the run: a term that differs per gene. The measured slope sizes it — protein
abundance rises 0.515 decades per decade of transcript, so **0.485 decades per
decade of per-gene compression** is what a future policy would have to supply,
and allocation as implemented supplies none of it. The form that would do it is a
*per-demander* saturable share, `capacity / (capacity + demand_i)` rather than the
per-pool form `competitive` computes today; a per-gene half-life is stage 4's. A
test asserts the per-demander form does beat the one-to-one arm on a compressed
truth, so the instrument is not merely asserted to be able to work.

**The coverage trap, registered in advance and then measured.** PaxDb holds a
protein because somebody could quantify it. The covered set's median transcript
level is 11 TPM; the 10,505 symbols expressed in GTEx and absent from PaxDb have
a median of 0.37 TPM, about thirty times lower. Every figure above is about the
covered set, which carries 86% of the assumed protein count per cell. The protein
and transcript arms are bought; the residue counts are free, from the packaged
proteome.

**Population mismatch, stated and not resolved:** GTEx is post-mortem bulk tissue
and PaxDb is an average over experiments, so even the closest pairing is two
populations. The gate asks whether a shared-capacity model beats a
proportionality constant at predicting a steady-state distribution, not whether
either is right about one cell.

**Consequence, as §5.3 agreed before any of this: the pool layer stays optional.**
It is expressible, it is charged, it refuses six modelling errors, it reproduces
a burden of the documented shape, and it improves no prediction of protein
abundance. A program that wants a burden should declare pools; a program that
wants abundance gains nothing by declaring them.

### 9.3 Stage 4 registered (2026-09-28), before the instrument existed

The registration is `genomeos/runtime/division_gate.py`, committed before
`genomeos/runtime/division.py` and `scripts/division_gate.py` were written and
before any number existed. It registers a **negative about the gate's own
wording**, which is the part worth reading before the result.

**§8 asks for something partitioning cannot give.** "Over several divisions
dilution must separate a stable protein from a short-lived one." Partitioning
multiplies every species in the cell by the same one-half; it carries no index
over proteins, exactly as stage 2's `Economy._share` carried none over genes
(§9.2). In a chase, after `n` divisions of interval `T`, a protein of half-life
`t` stands at `A0 · 2^-n · 2^(-nT/t)`, so

    log10(A_stable / A_short) = n · T · log10(2) · (1/t_short − 1/t_stable)

and the division term has cancelled. **Dilution separates nothing.** It removes
both proteins by the same factor, and turnover — which does the separating —
separates them just as well in a cell that never divides. In a steady state
with synthesis on it is worse than neutral: a level is
`s / (k_degradation + k_division)`, so dividing adds the same constant to both
denominators and **compresses** the separation. With the registered constants
(240 h and 6 h, dividing every 24 h) the stable:short ratio falls from 40.0 to
4.545, a compression of **8.80×**. The instrument must reproduce that, and a run
that showed dilution separating the two would mean the partitioning is indexed
by protein, which is a bug of exactly the shape that looks like the finding.

**What dilution does do, and what the gate therefore asks instead.** A protein
the program gives no half-life is never degraded, so in a cell that does not
divide it never leaves. Division is its only removal, by exactly `2^-n`. That is
the one separation dilution makes on its own — between what turnover can remove
and what only division can — and it is the registered positive clause.

**The eight clauses**, each able to fail: the closed form in both arms; the
dilution-only protein at exactly `2^-8` partitioned and exactly 1.0 duplicating;
the two arms' separation identical; the steady-state compression; exact
conservation at every division, every species and every seed; the binomial arm's
mean and variance against their analytic values with a whole-number daughter and
zero variance above the regime threshold; and two falsifiers — equal half-lives
must make the two proteins identical, and an unmet checkpoint must fire zero
divisions and leave the decay-only amounts.

**Registered as not askable rather than attempted.** Whether these half-lives
are any real protein's: `bio.std.human_turnover` holds *cell* lifespans (Sender
& Milo 2021), no protein turnover table is under `data/knowledge`, and §10
decision 6 stays open. Whether real partitioning is binomial: Huh & Paulsson
2011 show clustering makes it noisier, and this project holds no
partitioning-error dataset, so the binomial arm is checked against its own
analytic moments and never against a measurement. And whether a per-gene
half-life would move stage 2's gate (b): §9.2 named a per-gene term as the thing
that could, and a half-life is one, but that is a proteome-scale run needing the
turnover set this paragraph has just said the project does not hold. **Nothing
in stage 4 is evidence about gate (b).**

**Figures that must not move**, registered: area E's 1,439 of 1,439 fates,
terminal 496 of 555, 45 ambiguous decision points to 0 and 522 of 555 under
§7.5 — `runtime/body.py` is not touched and `partition` defaults to `duplicate`,
so no program that does not ask for partitioning behaves differently; stage 2's
burden factors and its 0.952005 MAE; stage 1's 13 of 13 and 1,123 of 1,123; and
the engine package's 12 of 12.

### 9.4 Stage 4 as measured (2026-09-28)

The language surface is one key, `partition` on an `event`
(`genomeos/lang/grammar.py`, `genomeos/ir/model.py`, `genomeos/lang/parser.py`);
the arithmetic is `genomeos/runtime/division.py`; the gate is
`genomeos/runtime/division_gate.py`, a pre-registration committed in `4ff3cd4`
before the instrument existed, run by `scripts/division_gate.py` into
`data/results/division_gate.json`.

**The negative first, and it was registered before the run.** Over eight
divisions 24 h apart, from 1e6 of each, a protein declaring a 240 h half-life and
one declaring 6 h end 9.392136 decades apart. They end **the same 9.392136
decades apart when contents are duplicated instead of partitioned**. Dilution
contributes 0.0 decades of it, to the last bit of a double. In a steady state
with synthesis on the stable:short ratio is 40.0 in a cell that never divides and
4.5455 in one that divides daily, so dividing **compresses** the separation
8.80-fold. §8's sentence — "over several divisions dilution must separate a
stable protein from a short-lived one" — is therefore refuted, not passed, and
for the same reason §9.2 gave about allocation: the operation carries no index
over the things it is asked to tell apart.

**What dilution does do, which is the clause the stage passes.** A protein the
program gives no half-life is not given one by the engine. It is never degraded,
so in a cell that does not divide it never leaves; division is its only removal.
Over eight divisions it falls by exactly 2^-8 = 0.00390625 when contents are
partitioned and by exactly nothing when they are duplicated. That is the one
separation dilution makes on its own, and it is the reason the stage is worth
having.

| Reading | Number |
|---|---|
| stable:short separation, contents partitioned | 9.392136 decades |
| stable:short separation, contents duplicated | 9.392136 decades |
| dilution's contribution to it | **0.0** |
| steady-state stable:short ratio, not dividing / dividing | 40.0000 / 4.5455 (**compression 8.80**) |
| a protein with no declared half-life, after 8 divisions | exactly 2^-8 partitioned, exactly 1.0 duplicated |
| binomial split of 20 copies over 2,000 seeds | mean 9.8895 (analytic 10), variance 4.8993 (analytic 5) |
| splits conserving the parent exactly | **2,000 of 2,000** at 20 copies and 2,000 of 2,000 at 1e6 |
| variance at 1e6 copies, above the regime threshold | **0.0**, i.e. deterministic halves |
| divisions fired with the checkpoint unmet | **0**, and the amounts are the decay-only values |
| registered clauses met | 8 of 8 |

**A second finding, measured and deliberately not scored.** The demo program
declares §5.2's own example, `event divide { cost: ATP 1e10 }`, against §5.1's
own ATP pool of 3e9, and the resource checkpoint therefore refuses **every**
division: one division costs more ATP than the cell holds. Neither number is
wrong. A pool is a standing stock and a division cost is a draw over a whole
cycle, and nothing in the engine converts between them — `regenerates` takes a
rate, this pool declares none (`from Glycolysis, OXPHOS`), and `Economy` never
integrates a rate over an interval. So the **resource** form of §8's checkpoint
is only askable against a standing stock today; the size form, which is the
event's own `when` clause through `ir.model.matches`, is what the registered
clause 8 tests and what passes. Making the resource form askable is stage 3's
job, since stage 3 is what supplies ATP over time. It is reported and excluded
from the verdict because the registration did not name it, exactly as gate (a)
excluded its unasked magnitude clause.

**What did not move.** `runtime/body.py` is untouched and `partition` defaults to
`duplicate`, so every committed program parses and runs as before: area E's
1,439 fates, 496 of 555 terminal and 522 of 555 under §7.5 stand, as do stage 2's
burden factors and its 0.952005, and stage 1's 13 of 13 and 1,123 of 1,123. The
engine package holds 12 of 12 with the new module and the regenerated grammar
inside it: rebuilding it gives 36 engine files to 38, 177 packaged test functions
to 192, and its isolated suite 285 passed to 300, with the same 3 skipped and 3
xfailed. `data/results/engine_package.json` itself is **not** updated by this
commit: overwriting it would delete lines three earlier commits by other lanes
put there, which is a force this lane does not take on somebody else's behalf.
The numbers above are the run; the file follows when whoever owns that result
rebuilds it.

**And what stage 4 does not license anybody to say.** Not that these half-lives
are any protein's: they are the program's declarations, this project holds no
protein turnover table, and §10 decision 6 is still open. Not that real
partitioning is binomial: Huh & Paulsson 2011 make it noisier than that, and the
binomial arm here is checked against its own analytic moments, never against a
measurement. And not that stage 2's gate (b) would now move: a per-gene half-life
is the per-gene term §9.2 asked for, but whether it supplies the 0.485 decades
per decade is a proteome-scale run needing the turnover set that does not exist
here. Nothing measured above is evidence about it.

### 9.5 Stage 3 registered (2026-10-10), before the instrument existed

The registration is `genomeos/runtime/metabolism_gate.py`, committed before
`genomeos/runtime/metabolism.py`, `bio.std.core_metabolism`, the two demo programs
(`data/demo/stage3_energy.bio`, `data/demo/stage3_erythrocyte.bio`) and
`scripts/metabolism_gate.py` were written, and before any number existed.

**Whether stage 3 may be built at all, checked first.** This document was a draft
for Albert's review, and nothing past stage 1 was to be built until he had steered
it. Of §10's decisions the one stage 3 touches is 4, whose cost accounting, and §10's
own table records it as settled *by* this stage, not as a precondition of it: "the
ATP budget gate (Buttgereit & Brand's hierarchy, Lynch & Marinov's totals): the
accounting that lands in the published order of magnitude wins". ROADMAP-HISTORY
item 5 reads "No decision blocks code any more", and ROADMAP area A lists stage 3 as
its next item 2 with no steer pending. So stage 3 is free to build, and decision 4
stays open unless the magnitude clause can be asked, which, below, it cannot.

**What the sources say, read before anything was built** (PubMed abstracts; the
bodies of the 1995 and 1999 papers are scans this project does not hold).
Buttgereit & Brand's hierarchy is an *order of sensitivity to energy supply*, not a
split of the budget: "pathways of macromolecule biosynthesis (protein synthesis and
RNA/DNA synthesis) were most sensitive to energy supply, followed by sodium cycling
and then calcium cycling across the plasma membrane. Mitochondrial proton leak was
the least sensitive". §6's parenthesis, "protein synthesis and the Na⁺/K⁺-ATPase
first", reads it as a split, and the abstract does not support that reading: it
states no shares and puts sodium cycling *after* macromolecule synthesis. Hinkle
2005 gives P/O "about 2.5 with NADH-linked substrates and 1.5 with succinate", and
2.3 and 1.4 as "also possible"; that is the spread the yields carry. Rossignol et
al. 2003 define the biochemical threshold effect as inhibiting a respiratory complex
"up to a critical value, without affecting the rate of mitochondrial respiration or
ATP synthesis". Mulquiney, Bubb & Kuchel 1999 (Biochem J 342:567, the first of the
three papers whose second §6 cites) estimate that "19%" of glycolytic carbon flux
passes through the 2,3-BPG shunt; none of the three abstracts states the ATP or
2,3-BPG ranges §6 names.

**The design.** One energy currency, refilled per hour by the processes a pool names
in `regenerates: from` (Glycolysis, TCA, OXPHOS), with every stoichiometric number
read from the program through `bio.std.core_metabolism`; the engine holds none.
Supply follows demand: pyruvate is oxidised as far as oxygen, mitochondrial capacity
and glucose allow, and the rest becomes lactate. That routing is the module's one
structural claim, and it makes the cell Pasteur-competent by construction, so the
Warburg effect is not expressible in stage 3. Every yield-dependent number is
computed at both ends of the spread, the high arm (malate–aspartate shuttle, P/O 2.5
and 1.5) and the low arm (glycerol 3-phosphate shuttle, 2.3 and 1.4), and a clause
holds only in both. Heteroplasmy removes the OXPHOS capacity of the copies carrying
the declared variant; the engine holds no table of which mtDNA variants are
pathogenic.

**The negatives, registered before the run.** Under `proportional` and
`competitive` every ATP consumer gets the same scalar (§9.2), so *no* hierarchy can
appear; `priority` carries one because the program writes it. The verdict is
committed to in advance: **the hierarchy is reproduced as an input, not predicted.**
The Pasteur clause holds by construction of the routing. The threshold sits at
h* = 1 − 12·D/(Y·C), set by the program's declared reserve, and with no reserve there
is none (the falsifier). Every cell here collapses without glucose, since the module
has no other fuel.

**The ten clauses**, each able to fail, with the numbers worked by hand:
1. yields per glucose 32.0 (high) and 28.0 (low), 30.0 for P/O 2.5 and 1.5 through the
   glycerol 3-phosphate shuttle, 2.0 anaerobic, and exactly 6 O₂ for 12 electron pairs
   in both arms (the combustion stoichiometry checks the module's constants);
2. with full oxygen and glucose, supply equals the declared demand of 1e10 ATP/h,
   glucose uptake is D/Y, O₂ six times that, lactate 0, OXPHOS share 28/32 and 24/28;
3. lowering oxygen (1, 0.75, 0.5, 0.25, 0.1, 0 of supply) never raises the OXPHOS
   share and ends it at 0, never lowers glucose uptake, and hits the registered table
   at every step (supply 1e10, 1e10, 1e10, 7.75e9, 5.5e9, 4e9 high);
4. with uptake not limiting, anoxia multiplies glucose uptake by exactly 16.0 and 14.0;
5. the red cell: TCA and OXPHOS make exactly 0, the yield is 2 − 2·0.19 = 1.62, and
   every reading is identical to the bit at every oxygen level;
6. with no glucose the red cell supplies 0 and every consumer is cut to 0;
7. under the priority order the first cut falls on macromolecule synthesis, then
   sodium, then calcium, with the anoxic factors exactly 1.0, 0.6, 0, 0; under
   `proportional` and `competitive` the four factors are equal at every step;
8. along a heteroplasmy sweep, OXPHOS ATP is unchanged to the bit below h* (0.8125
   high, 0.785714 low) and strictly lower above it; with no reserve it falls at once;
9. stage 4's demo division (1e10 ATP) is refused against the standing stock, as §9.4
   found, paid when the pool is regenerated over 24 h (2.19e11 and 1.71e11 available),
   and refused again with no glucose;
10. the refusals: heteroplasmy without chrM, a fraction outside 0..1 or summing past 1,
    and, at run time, an unknown process, OXPHOS without a mitochondrion, TCA without
    OXPHOS, a missing substrate pool, capacity or stoichiometric constant.

**Registered as not askable rather than attempted.** The budget's magnitude: no
measured per-cell ATP turnover is held, so the demo's rates are declared, and §10
decision 4 stays open. The budget's split: no shares are held. The red cell's ATP
and 2,3-BPG ranges: not held, and a steady-state concentration is a kinetic quantity
a module of stoichiometry and capacities cannot produce. Where Rossignol's measured
thresholds sit; how far real cells raise glucose uptake (they also lower demand);
aerobic glycolysis; other fuels; and proton leak, which spends the gradient rather
than ATP.

**Figures that must not move**: the 47 files and 288 checks of `bio test data/demo
genomeos/std data/organisms` at `d32fe44`; stage 2's burden factors and 0.952005
(the economy is unchanged when no capacity override is passed); stage 4's 8 of 8
(`runtime/division.py` untouched); area E's figures (`runtime/body.py` untouched).
Committed programs that write `regenerates: from Glycolysis, OXPHOS` with no
mitochondrion still compile; the runtime refuses them only when asked to run them.

### 9.6 Stage 3 as measured (2026-10-10)

The registration is `genomeos/runtime/metabolism_gate.py` (`e09edf7`), committed before anything it
measures; the language surface is `ac89b41`, the runtime, the standard-library module and the two demo
programs `fa5a4f1`, and the gate script `4668ee2`, which wrote `data/results/metabolism_gate.json` from
the committed tree (its manifest is complete and records its own code as committed). The script
refuses to run when a demo program declares a number the registration does not hold. Every
prediction below was worked by hand into the registration; nothing was tuned after a run.

**All ten registered clauses are met. The gate as §6 words it is not passed as worded**, because
three of its asks are not askable here and one is met only as an input:

- **(i) The ATP budget.** Its magnitude is **not askable**: this project holds no measured per-cell ATP
  turnover, so the demo's 1e10 ATP per hour is declared, and §10 decision 4, which this magnitude was
  to settle, stays open. Its split is **not askable**: Buttgereit & Brand's abstract states a hierarchy
  of *sensitivity to supply*, not shares, and §6's parenthesis "protein synthesis and the Na⁺/K⁺-ATPase
  first" reads it as a split the abstract does not support (it puts sodium cycling after
  macromolecule synthesis). The hierarchy itself is **reproduced as an input, not predicted**: it
  appears under the `priority` order the program writes and cannot appear under `proportional` or
  `competitive`, which hand every consumer one scalar — §9.2's finding again.
- **(ii) Oxygen and glucose.** Lowering oxygen shifts production to glycolysis and raises glucose
  uptake, **by construction** of the module's one routing claim (oxidise first); removing glucose
  collapses the glycolysis-only cell, and every other cell too, since the module has no other fuel.
- **(iii) The red cell** runs on glycolysis alone, with TCA and OXPHOS at exactly 0. Its ATP and
  2,3-BPG ranges are **not askable**: not held, and a concentration is a kinetic quantity a module of
  stoichiometry and capacities cannot produce.

| Clause | Reading | High arm | Low arm |
|---|---|---|---|
| 1 | ATP per glucose, aerobic / anaerobic | 32.0 / 2.0 | 28.0 / 2.0 |
| 1 | P/O 2.5 and 1.5 through the glycerol 3-phosphate shuttle (§6's "30") | 30.0 | — |
| 1 | O₂ per glucose, electron pairs | 6.0, 12.0 | 6.0, 12.0 |
| 2 | full oxygen and glucose: supply, glucose, O₂, lactate | 1e10, 3.125e8, 1.875e9, 0 | 1e10, 3.5714e8, 2.1429e9, 0 |
| 2 | OXPHOS share of the ATP made | 0.875 | 0.857143 |
| 3 | supply along oxygen 1, 0.75, 0.5, 0.25, 0.1, 0 | 1e10, 1e10, 1e10, 7.75e9, 5.5e9, 4e9 | 1e10, 1e10, 1e10, 7.25e9, 5.3e9, 4e9 |
| 3 | OXPHOS share along the same steps | 0.875 … 0.4516, 0.2545, **0** | 0.857 … 0.4138, 0.2264, **0** |
| 3 | glucose uptake, full oxygen to none | 3.125e8 → 2e9 (the ceiling) | 3.5714e8 → 2e9 |
| 4 | anoxia with uptake not limiting: glucose ×, supply | **×16.0**, 1e10 | **×14.0**, 1e10 |
| 5 | red cell: yield, glucose, lactate, 2,3-BPG made per hour | 1.62, 1.2346e8, 2.4691e8, 4.6914e7 | the same |
| 5 | red cell, every reading at oxygen 1 … 0 | identical to the bit | identical |
| 5 | red cell glucose without the shunt / with it | 0.81 | 0.81 |
| 6 | no glucose: red cell supply, its pumps' factors | 0, 0 and 0 | 0, 0 and 0 |
| 7 | priority: step of first cut, protein / RNA-DNA / Na⁺ / Ca²⁺ | 3 / 4 / 5 / 6 | 3 / 3 / 5 / 6 |
| 7 | priority factors at anoxia, Ca²⁺ / Na⁺ / RNA-DNA / protein | 1.0 / 0.6 / 0 / 0 | 1.0 / 0.6 / 0 / 0 |
| 7 | proportional, competitive: step of first cut, all four | 3, together | 3, together |
| 8 | heteroplasmy: last unchanged / first lower point (registered h*) | 0.80 / 0.82 (0.8125) | 0.78 / 0.80 (0.785714) |
| 8 | OXPHOS ATP at h = 0 … 0.80, 0.82, 0.90, 1.0 | 8.75e9, 8.4e9, 4.67e9, 0 | 8.57e9 … 0 |
| 8 | no reserve: OXPHOS ATP at h = 0.2 (of h = 0) | 7.0e9 (8.75e9) | 6.86e9 (8.57e9) |
| 9 | stage 4's division (1e10 ATP) against the stock 3e9 | refused, factor 0.3 | refused |
| 9 | ATP available over 24 h, regenerated / without glucose | 2.19e11 paid / 3e9 refused | 1.71e11 paid / 3e9 refused |
| 10 | refusals asked / refused; the unbroken control runs | 13 / 13; yes | — |

**Two findings, measured and deliberately not scored**, because the registration did not name them.
First, under `competitive` every consumer of an oversubscribed pool gets `capacity / (capacity +
demand)`: at oxygen 0.25 that is 0.4366 each, so the cell spends 4.37e9 of the 7.75e9 ATP it makes and
leaves the rest unspent. That is stage 2's policy as written ("nobody is served in full"), and it is
the regime's default for a pool with no allocation named. Second, stage 4's resource checkpoint, which
§9.4 found refused every division because a standing stock cannot pay a draw over a cycle, is now
askable: the same 1e10 ATP is paid when the pool is regenerated over one 24 h interval and refused
again when glucose is removed. That is clause 9, and it is scored; what is not is whether 24 h is any
real cell's cycle.

**What did not move.** `bio test data/demo genomeos/std data/organisms`: 47 files and 288 of 288 checks
at `d32fe44`; 54 files and 354 of 354 after, the four files other lanes added in between carrying 28
and stage 3's three files 38, so the 47 files present before still pass with 288 checks between them. Stage 2's burden and abundance
tests, stage 4's division and gate tests and the grammar, language, economy, engine-boundary and
licence-header tests pass; `runtime/body.py` and `runtime/division.py` are untouched, and `Economy`
behaves as before when no capacity is passed. Stage 2's and stage 4's demos, which write §5.1's
`regenerates: from Glycolysis, OXPHOS` with no mitochondrion, compile as before; the runtime refuses
them when asked to run them, since OXPHOS there has nowhere to happen.

**And what stage 3 does not license anybody to say.** Not that any cell spends 1e10 ATP an hour, or
spends it in equal quarters: those are the demo's declarations. Not that the hierarchy was predicted:
it was written in. Not that real cells raise glucose uptake 16-fold without oxygen: that is the
module's ceiling, and real cells also cut demand. Not that the threshold sits where Rossignol's
measured thresholds sit: its position is the declared reserve. Not that a cell ferments only without
oxygen: the routing assumes it, and the Warburg effect cannot be written in stage 3.

### 9.7 Gate (b) for a per-demander policy, registered (2026-10-10), before the run

§9.2 named the one change that could make gate (b) answerable — a share that saturates in each
demander's own demand — and sized it: 0.485 decades per decade of compression. `91e2a85` built it.
`allocation { pool: P; policy: saturable; half_saturation: K }` serves each demander `d K / (K + d)`
of its own demand `d`, then enforces the pool on what is served (one cut, the pool used exactly). K is
the program's number, kept in the allocation's `attrs`; a missing, unknown, zero, negative or
non-finite K is refused by the parser and by the runtime, `half_saturation` under any other policy is
refused, and a regime's default cannot be `saturable`. A demander is cut by its own saturation even
when the pool is not short; that is the claim, and it is not a burden, because no other demander's
demand enters its factor until the pool is short. `proportional`, `priority` and `competitive` are
unchanged: `scale()` and `report()` were compared with `==` against the previous `economy.py` on 3,000
random programs.

The registration is `data/results/saturable_gate_registration.json`, written by
`scripts/saturable_gate.py register` (code frozen in `504c536`) with no protein value read, and its
binding words are `REGISTRATION` in that script; `run` refuses if either differs or if a pinned file
changed. In short:

- **Form.** One pool, Ribosomes, 5e6 (§5.1, as gate (b)); each gene charges 1 per chain at a work
  rate equal to its TPM (gate (b)'s convention); the allocation is parsed from `.bio` text, and the
  predictions are the runtime's factors, checked against the closed form to 1e-9.
- **Parameters.** K (on a grid, log10 K from −1 to 7 in 0.01 steps, TPM units) and the constant c,
  both **fitted on training genes only**; no independent measurement supplies a per-gene
  half-saturation in TPM units, so none is cited. The capacity is §5.1's and is not fitted.
- **Split.** Held out when the first 32 bits of sha256 of the gene symbol are odd, fixed before any
  protein value was read. **Circularity guard: no protein value of a held-out gene enters any fit**;
  a test replaces every held-out value and asserts that no fitted number moves.
- **Metric and bar.** Gate (b)'s: log10 MAE on held-out genes against the one-to-one baseline fitted on
  the same training genes; PASS needs an improvement of at least 0.05 with the 95% interval of 2,000
  bootstrap resamples of held-out genes above 0, and gate (b)'s secondary in band.
- **Reported beside the verdict.** The slope of the policy's held-out predictions on log10 TPM against
  the registered 0.515 [0.499, 0.531]; a two-parameter straight-line comparator, which says whether an
  improvement is specific to the saturable form; the gene-blind arms through the same runtime; the
  K = capacity form §9.2 named; and whether the capacity binds at real demand.
- **Planted contrasts, run at registration.** On a planted saturable truth the policy passes and
  recovers log10 K 1.48 against 1.477, while `proportional` and `competitive` move the held-out error
  by 0; on a null truth it does not pass; on a straight compressed truth (slope 0.5, no saturation) it
  **passes the bar and the comparator line does better**.
- **Stated before the run.** That last contrast is the expectation for real data too: the bar is
  against a slope-1 baseline and the measured slope is 0.515, so almost any fitted compression is
  expected to clear it. A PASS would therefore be weak evidence for the saturable form, and the
  comparator is the reading that can tell the form from compression in general. Nor can the gate tell
  a saturable allocation from noise in the transcript arm attenuating the slope.

**Amendment 1 (2026-10-10), before any number was seen.** Run 1 computed every reading and was then
refused at the write by the result contract, because the manifest's `partitions` was a string; the
payload went to `data/quarantine/results/saturable_gate_result.json`. The run prints nothing before
the write and the file was not opened. The fix (`411b8ff`, `4f71f08`, insertions only) makes
`partitions` a dict and has `run` compare its payload with run 1's key for key; the new pins are in
`data/results/saturable_gate_registration_amendment1.json`, whose binding words are the registration's.
The reading is run 2's only if the two payloads are identical; otherwise run 1's is the registered run.

### 9.8 Gate (b) for the per-demander policy as measured (2026-10-10)

**The registered reading, verbatim:** *On held-out genes, the saturable policy's log10 MAE is 0.949823
against the fitted baseline's 0.949969: an improvement of 0.000146 [0.000023, 0.000287]; no improvement
detected at gate (b)'s registered bar of 0.05.*

| Reading | Number |
|---|---|
| genes analysed (gate (b)'s join, rebuilt) | 15,159: 7,562 training, 7,597 held out |
| K fitted on training | log10 K 4.39, about 24,547 TPM; 1 gene of 15,159 has a TPM above it |
| held-out log10 MAE, baseline / saturable | 0.949969 / 0.949823 |
| improvement | **0.000146** [0.000023, 0.000287], bar 0.05: **NO_IMPROVEMENT_DETECTED** |
| comparator line (b = 0.608 fitted on training), held-out MAE | 0.879610, better than the saturable arm by 0.070212 [0.061179, 0.079399] |
| slope of the policy's held-out predictions on log10 TPM | 0.999, outside the registered 0.515 [0.499, 0.531] |
| K = capacity arm (the form §9.2 named) | 0.949968 |
| `proportional`, `competitive` (pool short) | one factor each; held-out improvement 0.0 and 1.1e-16 |
| runtime against closed form | largest gap 0.0 in log10 |
| secondary, total protein per cell (no allocation in it) | 1.36e10, within the band, as in gate (b) |

The comparator reading, verbatim: *A two-parameter straight line does better than the saturable
policy on the same held-out genes (0.879610 against 0.949823, difference -0.070212 [-0.079399,
-0.061179]), so the improvement is not specific to the saturable form.*

**Whether the capacity binds at real demand, as a stated fact:** at the fitted K the ribosome pool
serves 7.859e5 against a capacity of 5e6, out of a demand of 8.712e5, so **the capacity does not
bind**; every factor below 1 comes from K, not from the pool. As in gate (b), §5.1's capacity is not
reached at the demand a real transcriptome implies.

**The expectation registered before the run is refuted.** It said almost any fitted compression would
clear a slope-1 baseline. The saturable form did not: it has slope 1 below K and flattens above it, and
the K that minimises training error sits where it touches one gene, so it supplies almost none of the
0.485 decades per decade. A straight line in log-log, with the same number of fitted parameters, does
lower the held-out error. That is a statement about which functional form fits this join, not about
how a cell allocates ribosomes.

**Amendment 1's comparison, and a defect in it.** Run 2 compared its payload with run 1's key for key
and found one difference: `registration_amendment`, the key that names the amendment file written after
run 1, which run 1 therefore could not carry. The comparison should have left that key out. Every
other key, the verdict, the reading and every number among them, is identical. By the registered rule
run 1's payload is the registered run, and its reading is the one above word for word. Run 1's file
was deleted when run 2 wrote the result (`save_result` removes a superseded quarantine copy,
genomeos/results.py); its sha256 is kept among run 2's declared inputs.

**Consequence.** Gate (b) is still failed for the pool layer, now with a per-demander policy as well
as the three gene-blind ones, and §5.3's consequence stands: the pool layer stays optional. Results:
`data/results/saturable_gate_registration.json`, `saturable_gate_registration_amendment1.json` and
`saturable_gate_result.json`; code `scripts/saturable_gate.py`; tests `tests/test_saturable_gate.py`.

## 10. Open decisions for Albert

1. **Opt-in or mandatory.** Located rules apply to programs that declare a
   compartment. Should v0.5 make locations mandatory everywhere?
   **Priced 2026-09-15** (`scripts/mandatory_locations.py`,
   `data/results/mandatory_locations.json`). This was filed as a preference, but
   what it needs is a price, and a price is measurable. Over every committed
   program: **27 of 36 are unaffected**, having no genes or proteins at all (the
   organism layer is about cells, not molecules), and **2 of the remaining 9 are
   located already**. The cost falls entirely on the other 7, and it splits in
   two:

   - **Where the program is about human molecules, the facts exist.** 193 genes
     in `noncoding_chr21.bio` carry no locus, but 163 of them are symbols the
     packaged proteome already holds, so their compartment is a lookup (a gene
     is placed by the compartment whose genome holds its chromosome), not an
     invention. UniProt places the proteins: SDHB in the mitochondrial inner
     membrane, TP53 in cytoplasm and nucleus.
   - **Where it is not, there is no source at all.** All **19** proteins that
     would need a compartment — the repressilator's TetR, LacI and CI, the
     worm's factors, the demo circuits — are not human symbols, so no human
     source can ever place them. Their author would write a compartment with no
     citation, which under §2 is `inferred` at 0.3 and **caps the confidence of
     everything computed from it**.

   So mandatory is free for the human molecular programs and converts a
   language-design preference into a confidence penalty on every program that is
   not about human molecules. That is an argument for the rule we have (located
   when you declare a compartment), or for a narrower mandatory — required only
   where the program names a genome whose molecules a source can place. It is
   Albert's call, but it is no longer a call made without a number.
2. **Amounts or concentrations** (stage 2 is held on this one: pools and costs would bake it in). Stage 1 keeps amounts per compartment, so
   Hill thresholds are amounts. Concentration needs volumes, which vary per cell
   type; the choice changes what every threshold in the project means.
   **Priced 2026-09-17** (`scripts/amounts_or_concentrations.py`,
   `data/results/amounts_or_concentrations.json`), and the first finding is that §9's
   proposed test cannot be run. It says to run the burden gate "both ways on cells of
   different volume". `volume` on a compartment is documented in two places — the
   compartment table above and BIOLANG-GRAMMAR.md — as **a fraction of the cell's
   volume**, and a fraction cannot express cells of different volume. A concentration
   is molecules over an *absolute* volume and no v0.4 construct carries one, so the
   concentration arm of the gate has nothing to divide by. **The blocker is a language
   field, not data.**
   What that field would cost is countable, and it is smaller than the objection above
   suggests. Of **39 committed programs, 2 declare a compartment** — `human.erythrocyte`
   (Cytosol at 1.0) and `human.oxphos` (four fractions from Alberts Table 12-1). Both
   stand for cell types whose **absolute** volume is published, so neither author would
   invent one and §2's penalty for an uncited fact does not apply. "Volumes vary per
   cell type" is true in general and costs **two citations** at today's scale.
   **And the price only rises.** It is this small *because* stage 2 is not built: **0
   pools and 0 costs** are committed, so almost nothing is denominated in either unit
   yet. Every program stage 2 adds raises the cost of changing the answer. That is an
   argument about *when* to decide, not about which unit is right — and this decision
   does not say which unit is right, because the gate still cannot separate them until
   the language can express an absolute volume.
   *(2026-10-10: stage 2 was implemented on 2026-09-21. Of 45 committed `.bio` files, 5
   declare a compartment — the two above and three demos, `concentration_threshold`,
   `stage2_pools` and `stage4_division` — and 4 pools are committed, all in the last two.)*

   **Resolved 2026-09-19 by Albert: add the absolute volume, unblock stage 2.** A
   `compartment` gains an absolute volume beside the fraction it already carries, so a
   concentration becomes expressible and §9's gate becomes runnable — it was not runnable
   before, since a fraction of a cell cannot express cells of different volume. The price
   is the one that was measured: one language field, and one citation per cell type that
   declares a compartment, which is two today and both published. The fraction stays,
   because it is what the existing programs state and what the tree needs; the absolute
   volume is what a threshold can be divided by. This resolves the decision by making the
   arm expressible rather than by choosing the unit: stage 2 is no longer held, and the
   gate can now be run both ways and report which unit a threshold should be stated in.
3. **Membranes as compartments.** A membrane is a node of the tree facing both
   sides. The alternative (a membrane as the edge between two spaces) is closer
   to a bigraph link but makes receptors harder to place.
4. **Whose cost accounting.** Polymerisation only, or full biosynthetic cost
   of precursors (Lynch & Marinov). They differ several-fold.
5. **Default allocation.** `competitive` for transports is a choice of
   mechanism; pools may deserve no default at all, forcing every program to name one.
   **Status 2026-10-10 (§9.7, §9.8).** A fourth runnable policy, the per-demander `saturable`, was
   put to gate (b) on held-out genes: no improvement detected at the registered bar (0.000146
   [0.000023, 0.000287] against 0.05). Gate (b) still cannot choose a default; the decision stays open.
6. **The protein turnover set** for stage 4, and whether a human set is
   required or a mouse one (Schwanhäusser et al. 2011, Nature 473:337) is acceptable.
   **Still open, and stage 4 shipped without it (2026-09-28, §9.4).** Every clause
   the stage could ask is a property of its own arithmetic — conservation, the
   closed form, the gene-blindness of dilution — so none of them needs a measured
   half-life, and the registration recorded the comparison as *not askable* rather
   than making one up. What a turnover set would buy is a different question from
   the one stage 4 answers: whether a per-gene half-life supplies the 0.485 decades
   per decade of per-gene compression §9.2 named, which is a proteome-scale run.
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
    colinearity data). It should not bake node identity in, since area B's
    boundary caller may move every node. **Specified 2026-09-15 as §7.6**, with
    the part that cannot be done stated there: the position claim is checkable
    against coordinates we hold and the time claim is gated only by the worm,
    because the project holds no expression time course for a human HOX cluster.

**Which of these can be settled by measurement rather than preference.**
*What is actually decidable today, 2026-09-15, without the API and without stage
2:* **1** is priced (below, done); **8** is resolved and **9** has been run and
stayed open, needing Fukushige & Krause's per-stage frequencies, which is a paper
to read rather than a run to do; **4, 5 and 6 are blocked behind decision 2**,
not behind any missing work of ours, because each needs the pools or the
division layer that decision 2 holds; **3 and 7** are genuinely preference and
data availability, and both are already argued in place; **10** must be built
before it can be measured. So of the ten, exactly one was decidable and cheap,
and it is done. *(2026-10-10: decision 2 was resolved on 2026-09-19 and stage 2 implemented on
2026-09-21; 10 was built on 2026-09-15, `144f399`.)*

| Decision | Settled by | How |
|---|---|---|
| 1. opt-in or mandatory locations | **priced, 2026-09-15** | not a truth but a cost, and the cost is countable: 27 of 36 programs unaffected, 2 located already, 163 gene facts a lookup, and **19 proteins with no possible source** whose invented compartments would cap every confidence downstream |
| 2. amounts or concentrations | **resolved by Albert, 2026-09-19** | a compartment states an absolute volume beside its fraction, so a threshold may be a concentration and is converted where its rule acts (§4.1, §9.1 item 4). The test this row proposed — the burden gate both ways on cells of different volume — was **not** what settled it and has still not been run: gate (a) ran on a demo program that declares no cell (§9.2), so the volume arm of that gate remains open |
| 3. membranes as nodes or edges | preference, lightly constrained | both passed stage 1; the band-3 transport forced the "across one membrane" rule, which either shape can express |
| 4. whose cost accounting | **measurable** | the ATP budget gate (Buttgereit & Brand's hierarchy, Lynch & Marinov's totals): the accounting that lands in the published order of magnitude wins |
| 5. default allocation | **asked and answered with a negative, 2026-09-21** | gate (b) ran under each policy and they make identical predictions: `proportional` and `competitive` both hand every demander of a pool one scalar, so no abundance measurement can choose between them (§9.2). `priority` is the only per-gene policy and needs an ordering over 19,000 genes that no measurement supplies. The decision is not Albert's to make on this evidence, because the evidence cannot see the difference |
| 6. protein turnover set | **measurable** | the dilution gate run with the human and the mouse sets; if they disagree beyond their own spread, the human set is required |
| 7. which homeostats first | data availability, not preference | the pair whose failure direction is best documented; Na⁺/K⁺ and pH stand |
| 8. `fates: first` default | **settled by measurement** (resolved above) | |
| 9. commitment as refusal or rate | **measured, and still open** | the series ran and passed as a refusal (§7.2a); the published decline over the fourth hour is what a rate would fit, and separating them needs the per-stage frequencies |
| 10. order as a construct | measurable once built *(2026-10-10: built 2026-09-15, `144f399`; position passed, time passed on the worm, activation not testable, §9)* | a published activation series |

Decision 9's instrument is built and run (§7.2a, 2026-09-14): the series is
`data/organisms/celegans/plasticity.bio`, `bio test` runs it, and removing
either construct breaks an arm. The integrated reads landed on 2026-09-15 and
are area E's to gate; re-decision is §7.5; and the precursor-level `commitment`
turned out to need no construct at all (§7.2a decision 3). The next engine work
that needs no preference from Albert is therefore decision 10's construct
sketch. 4 and 6 still need stage 2 or 4 first; 2 no longer holds stage 2, which
shipped and was gated on 2026-09-21 (§9.2); and 5 is now closed by a negative
rather than by a preference. *(2026-10-10: decision 10's construct was built on 2026-09-15,
`144f399`. Stage 4 was implemented on 2026-09-28, so 6 no longer waits on a stage; it stays open
and no stage 4 clause needed it (§9.4). 4 waits on stage 3's ATP budget, which is not
implemented.)* *(Later on 2026-10-10: stage 3 is implemented, and the budget's magnitude, which was
to settle 4, is registered not askable: no measured per-cell ATP turnover is held (§9.5, §9.6). So 4
stays open; the demo states polymerisation-only costs and says so.)*
