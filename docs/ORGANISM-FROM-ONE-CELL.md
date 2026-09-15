# From one cell to an organism: research, design, plan

Status: proposed 2026-09-10; all six steps built 2026-09-11 (see the end of the document and docs/PROGRESS.md).

The ask: a computational model that starts with one cell and, by reading the
genome under timers and signals, grows, maintains and kills cells until an
organism exists; not an aesthetic body, but a language, code and libraries
that replicate as faithfully as today's knowledge allows how an organism
grows from its genetic load, controlling each cell, group of cells and node.

This document has three parts: what was found (research), what to build
(design), and how to build it (plan and proposal).

---

## Part 1. Research

### 1.1 What the premise gets right, and the one correction

A cell does decide to divide, differentiate, migrate, rest or die by reading
the part of the genome that is open to it, under the control of timers
(cell-cycle length, developmental clocks) and signals (contact, gradients).
That loop is codable, and the invariant lineage of C. elegans proves it can
be checked cell by cell against reality.

The correction: the seed is not DNA alone. The zygote carries maternally
deposited RNA and proteins (in C. elegans: SKN-1, PIE-1, PAL-1, POP-1,
PAR proteins), an epigenetic starting state and mitochondria, and the early
embryo is "dominated by the maternal deposit" with few new transcription
factors induced before the 15-cell stage (Cole et al. 2024, eLife). Every
later decision also uses signals from neighbours. GenomeOS already separates
genome from state, so the honest model is: one genome, one bootstrap cell
state, one environment, then everything else emergent. This keeps the
"only the genetic load" spirit without pretending the biology is simpler.

### 1.2 The ground truth that exists (C. elegans)

C. elegans is the only organism with a complete, named, timed, cell-by-cell
ground truth, which makes it the only place where "faithful" can be measured
rather than claimed.

| fact | value | source |
|---|---|---|
| somatic cells, adult hermaphrodite | 959 | Sulston et al. 1983 |
| nuclei generated in the embryo | 671, of which 113 die | Sulston 1983 |
| cells at hatching | 558 | Sulston 1983 |
| total programmed deaths (embryo + larva) | 131 | Sulston 1983 / 1977 |
| embryogenesis at 20 °C | ~14 h | WormAtlas |
| cell-cycle length variability across 18 embryos | median CV 0.044 | Richards et al. 2013 |
| C and D lineages more variable | median CV 0.061 | Richards 2013 |
| founder lineages present | by the 28-cell stage | Cole 2024 |
| most cells' progeny single-fated | by the 102-cell stage | Cole 2024 |

Datasets, checked on 2026-09-10 for availability and licence:

1. **Complete timed lineage tree** (WormWeb, Nikhil Bhatla, CC BY 1.0).
   One JSON file, 320 KB, fetched and parsed: 2,183 nodes, 1,092 terminal
   cells, 131 deaths, per-cell birth and division times in minutes from P0
   to adult (max 5,763 min), terminal names (AINL, ILshL, "mu bod" ...) and
   tissue types (neuron 306, hypodermis 195, muscle 152, reproductive 151,
   intestine 34, ...). Embryonic and post-embryonic. This is the reference
   the simulator is scored against.
2. **Richards et al. 2013** (Dev Biol): 18 wild-type embryos traced to the
   >550-cell stage; Supplementary Table 2 gives mean division times,
   cycle lengths, SDs, migration distances per cell. The AceTree model zip
   link on the Murray lab site is dead (404); the supplementary Excel via
   PMC3548946 remains.
3. **Packer et al. 2019** (Science, GEO GSE126954): 86,024 embryonic single
   cells with `lineage`, `cell.type`, `cell.subtype`, `embryo.time` columns
   in `GSE126954_cell_annotation.csv.gz` (header verified). Gives
   lineage → terminal cell type → time; 502 terminal and pre-terminal types.
4. **Ma et al. 2021** (Nat Methods, Zenodo 10.5281/zenodo.4737593, CC BY 4.0):
   4D protein atlas of 266 transcription factors in 1,204 lineaged cells to
   the bean stage (~600 cells), 1.25-min resolution, one CSV per TF reporter
   (88 MB zip). 84 lineage-specific, 101 tissue-specific, 58 time-specific
   TFs; ten TFs distinguish a cell from 93% of others. This is the measured
   "reader configuration" per cell per minute.
5. **Cole et al. 2024** (eLife, GEO GSE83523): 840 cells, 1- to 102-cell
   stages, 119 cell states; founder patterning codes (AB: ceh-43, ceh-13,
   ceh-36, unc-30; MS: ceh-32 / ceh-51; C: vab-7 / hmbx-1; E: GATA factors).
6. **Boolean model of early fate** (Genetics 2025, "Quantitative resolution
   of cell fate in the early embryogenesis of C. elegans"): 17 nodes
   (CDC-25.1, CDK-1, WEE-1.1, PAR-6, LGL-1, PLK-1, GLP-1, POP-1, SKN-1,
   PIE-1, PAL-1, TBX-35, TBX-37, PHA-4, HLH-1, END-3, LIN-26), 1- to
   16-cell stage, literature-derived rules, 58–82% consistency with
   transcriptomes; weakest for maternal factors (post-transcriptional). No
   model file published; the rules are in the paper and can be transcribed
   as cited BioLang rules.
7. **Digital Development** (Du et al. 2014 Cell; NAR 2016; online, HTTP 200):
   204 genes perturbed, 1,368 embryos, ~600,000 digitised cells, 820
   predicted cell-specific fate changes, a gene → cell → landscape model.
   This is the knockout validation set: what a lineage does when a gene is
   lost.
8. **Post-embryonic lineage** (Sulston & Horvitz 1977, WormAtlas): ~50
   blast cells resume division after hatching through L1–L4; included in
   dataset 1 as timed nodes.

WormBase REST is behind a Cloudflare challenge for non-browser clients, so
it is not a reliable programmatic source; dataset 1 plus WormAtlas cover
the same lineage.

### 1.3 The human numbers (for the top of the hierarchy)

| fact | value | source |
|---|---|---|
| cells, reference adult male 70 kg | ~36 × 10^12 | Hatton et al. 2023, PNAS |
| adult female 60 kg / child 32 kg | ~28 / ~17 × 10^12 | Hatton 2023 |
| nucleated cells, male | ~7 × 10^12 | Hatton 2023 |
| cell groups × tissues | 1,264 groups (~400 types) × 60 tissues | Hatton 2023; humancelltreemap.mis.mpg.de (online) |
| cells replaced per day | 0.33 × 10^12 (~90% blood) | Sender & Milo 2021, Nat Med; github.com/milo-lab/cellular_turnover |
| biomass replaced per day | 80 ± 20 g | Sender & Milo 2021 |
| gastrulation | days 14–21, Carnegie stages 7–9 | 2024 CS7/CS8/CS9 spatial atlases |
| segmentation clock period human / mouse | 5–6 h / 2–3 h | Matsuda et al. 2020 |
| tempo cause | ~2× slower protein degradation and cycle in human cells | Rayon et al. 2020; 2024 SILAC follow-up |

The tempo result matters for the design: the same program runs at a
species-specific speed set largely by protein stability, so timers need a
tempo scalar, not a rewrite, between species.

### 1.4 State of the art in growing embryos in software

- **Whole-embryo mechanistic models** exist only for small windows: mouse
  blastocyst inner-cell-mass fate (spatial-stochastic, PLOS Comp Biol
  2024), early lineage specification with division (2017), C. elegans
  morphogenetic domain signalling (CSBJ 2022). A 2026 review
  ("Multiscale modeling of embryonic morphogenesis from fertilization to
  organogenesis") confirms there is no end-to-end embryo simulator.
- **Agent frameworks**: CompuCell3D and Morpheus (Cellular Potts, easy,
  do not scale past ~10^5 cells); PhysiCell (off-lattice, ~10^6 cells) and
  PhysiBoSS 2.0 (December 2024; a MaBoSS Boolean network inside every
  PhysiCell agent, which is exactly the per-cell reader/decider idiom);
  Tissue Forge (particle physics); BioDynaMo (1.7 × 10^9 agents on one
  server, 5 × 10^11 distributed on 84k cores, Apache 2.0, C++). None ships
  cell cycle, apoptosis or differentiation rules out of the box; those are
  always model-specific, which is what BioLang is for.
- **Coarse-graining** of populations is standard in tumour modelling
  (hybrid agent + continuum, coarse-grained ABMs): agents where decisions
  happen, counts where they do not.
- **Virtual cell (AI)**: Arc's Virtual Cell Challenge 2025 (over 1,200
  teams; winners BioMap and Altos Labs) predicts perturbation responses in
  one cell type from data; it does not grow anything and has no lineage or
  time. Useful later as a predictor for "what does this cell do if gene X
  is lost", not as an engine.

Conclusion: nobody has an end-to-end "one cell to organism" engine; the
parts exist (lineage data, per-cell TF atlases, Boolean fate logic, agent
frameworks, coarse-graining). The gap is a language and a runtime that
tie them together with evidence and a measurable score. That is the
GenomeOS position.

### 1.5 What GenomeOS already has, and the concrete holes

Inventory from the codebase (2026-09-10, 125 tests passing, 4 skipped):

| exists | where | limit |
|---|---|---|
| lineage engine, P0 → ~600 named cells, division queue | genomeos/organism/celegans.py | 9 scalar params from a .bio file; no gene state; every cell divides forever except P4; EMS mislabelled as lineage P; no truth comparison |
| cell runtime (ageing, years) | genomeos/runtime/cell.py | parameters hardcoded in Python, not read from BioLang; division is a counter, never creates a cell; nothing dies |
| gene network ODE with external clamp | genomeos/runtime/grn.py | fine as the continuous decider |
| Boolean network (.bnet), attractors | genomeos/runtime/boolean.py | fine as the fast decider |
| 2-D fields, sources, fate rules | genomeos/runtime/spatial.py | fixed grid; no division, movement, add or remove |
| clock-and-wavefront segmentation | genomeos/runtime/segmentation.py | the timer + signal gating precedent |
| one network per cell with a clamped signal | genomeos/runtime/gastrulation.py | the per-cell agent idiom to generalise |
| process-bigraph composite | genomeos/runtime/compose.py | where a body process joins the shared clock |
| uncertainty per level | genomeos/runtime/uncertainty.py | the `organism` level is defined and never populated |
| BioIR: Gene, Transcript, Protein, CellType, Rule, Event, Effect, Parameter, RegulatoryElement | genomeos/ir/model.py | no Timer, Stage, Organism, Lineage, Signal, Population; Event is a rate (1/yr), not a delay |
| BioLang blocks: gene, protein, region, element, rule, param, cell_type, event, transcript | genomeos/lang/parser.py | no organism / stage / timer / signal; v0.3 `experiment` and `organism` reserved in docs/ARCHITECTURE.md |
| std library | genomeos/std | cell_types.bio (6 types), ageing.bio; no development or timer modules |
| domains and regulatory targets | genomeos/genome/domains.py, regulation.py | node model exists at genome level, not yet a BioIR block |
| BioLib timer.* (9 libraries) | genomeos/lib/catalog.py | catalogue only |

The single largest hole: nothing interprets `Module.events`. Event and
Effect exist in the IR and in bio.std.ageing, `Event.applies(context)`
exists, and no runtime reads them. One → many (a cell creating a cell) does
not exist in any tick loop. Terminal differentiation exists only as the
hardcoded P4.

---

## Part 2. Design

### 2.1 Principles

1. **One genome, one bootstrap state, one environment.** The organism
   block names all three. Nothing else is input.
2. **Hierarchy of control, as in the body.** Agents where a decision is
   pending; populations (counts) where a lineage is committed. This is the
   node idea applied to cells, and it is what makes 10^13 tractable.
3. **Every rule is cited; UNKNOWN is legal; the organism level of the
   uncertainty report is populated** from the rules that actually fired.
4. **Faithfulness is a number.** Every run of a reference organism prints
   its distance to the ground truth: fates, timings, counts, deaths.
5. **Language first.** All biology lives in BioLang modules; the runtime
   is a small interpreter; zero dependencies in the core; process-bigraph
   remains the optional scheduler; PhysiCell/BioDynaMo are not adopted now
   (C++ and heavy), but the Body runtime keeps a clean agent interface so
   an external physics engine can be plugged in later.
6. **Stream, distil, discard.** Raw atlases are downloaded once, reduced to
   small summaries in data/results or data/knowledge, and deleted.

### 2.2 The cell loop

Per cell (agent) or per population, at each of its events:

```
read     active nodes = genome ∩ epigenetic mask(cell_type, maternal factors) ∩ signals(neighbours, fields)
update   network state (Boolean or ODE) over the active genes
decide   divide | differentiate | migrate | quiesce | die     (event whose `when:` matches)
write    daughters inherit the mask, plus asymmetric inheritance where a rule says so
```

Timers gate which events are eligible: cell-cycle timers per lineage, stage
timers for the whole organism, tempo scalar per species.

### 2.3 BioIR v0.3 additions

| entity | fields | note |
|---|---|---|
| `Organism` | species, genome ref, bootstrap (cell_type, factors), environment, stages, tempo | the program root |
| `Stage` | name, starts (min or trigger), ends, evidence | gates `when: stage = X` |
| `Timer` | name, duration, unit, per (lineage / cell_type / stage), tempo-scaled, evidence | delay/deadline primitive; distinct from `Event.rate` |
| `Signal` | name, kind (contact / gradient / systemic), source, range, evidence | read by cells; contact needs neighbours, gradient needs a field |
| `Decision` (Event subtype) | action ∈ {divide, differentiate, migrate, quiesce, die}, `when:`, daughters (names or rule), to (cell_type), asymmetric factors | the five verbs of development |
| `Population` | lineage root, cell_type, count, timer, mask | runtime-only coarse node, also serialisable |
| `Domain` | chrom, start, end, genes, boundaries | promised in NODES-READER-WRITER.md, promoted into the IR so the mask can be per node |

`CellType.expresses` remains the reader mask; masks compose with maternal
factors and with Domain-level on/off.

### 2.4 BioLang v0.3 syntax (proposed)

```
module bio.organism.celegans

import bio.std.development
import celegans.lineage_rules

organism Celegans {
  species: Caenorhabditis elegans
  genome: WBcel235
  tempo: 1.0
  bootstrap: cell_type = Zygote; factors = SKN-1, PIE-1, PAL-1, POP-1, PAR-3, PAR-2
  environment: temperature = 20 C
  evidence: experimental "Sulston et al. 1983"
}

stage Cleavage   { from: 0 min;  to: 100 min }
stage Gastrulation { from: 100 min; to: 350 min }

timer ab_cycle { duration: 20 min; per: lineage = AB; lengthening: 1.12; evidence: experimental "Richards et al. 2013" }

signal Notch1 { kind: contact; ligand: APX-1; receptor: GLP-1; evidence: experimental "Priess 2005" }

decision divide {
  when: cell_type = Blastomere, timer = elapsed
  daughters: anterior, posterior
  asymmetric: POP-1 high -> anterior
  evidence: experimental "Lin et al. 1998"
}

decision differentiate {
  when: lineage = E, stage = Gastrulation, END-3 = present
  to: IntestinalPrecursor
  evidence: experimental "Maduro 2010"; confidence: 0.8
}

decision die {
  when: cell = ABalaaaalar
  evidence: experimental "Sulston 1983"
}

observe count, fates, deaths every 10 min
assert count at 100 min in 24..34
compare lineage against data/results/celegans_lineage_cells.json
```

The `experiment` block from the v0.3 note keeps its planned role (fork,
perturb, compare); `organism` is the thing an experiment runs.

### 2.5 Body runtime

`genomeos/runtime/body.py`:

- `Body`: a tree whose nodes are `CellAgent` or `Population`; a
  discrete-event queue for decisions (division, death, stage change),
  plus optional continuous fields from `spatial.py` stepped between events.
- `CellAgent`: name, lineage, cell_type, mask, factors (dict), network
  state (Boolean network by default, `NetworkRuntime` optional), position
  (optional), timers, birth time.
- `Population`: created when a branch is committed (no eligible decision
  other than its own cycle timer); holds count, cycle timer, cell_type;
  expands back to agents when a decision becomes eligible for it.
- `Reader`: computes the active gene set for a cell (mask ∩ domains ∩
  signals) and feeds the decider.
- `LineageDiff`: scores a run against a reference tree: fraction of
  reference cells born, fate accuracy on terminal cells, median timing
  error, deaths matched, extra or missing cells.
- Uncertainty: every decision that fires contributes its evidence to the
  `organism` level; timers to `cellular`; rules to `molecular`.

Existing engines plug in: Boolean or GRN as deciders, spatial fields as
signal carriers, segmentation as a timer generator, gastrulation logic as a
decision set, ageing runtime as the maintenance phase for adults.

### 2.6 Human scale by hierarchy

The same runtime, different resolution: Carnegie stages as `Stage`s,
lineage populations (ectoderm → neural crest → ...) as `Population`s,
adult tissues seeded from Hatton's 1,264 cell groups, turnover from Sender
& Milo, tempo scalar ~2 relative to mouse. This yields a body of ~36 × 10^12
cells as a few thousand populations with cited counts and turnover rates,
grown from one population over developmental time, with the uncertainty
report saying honestly that nearly all of it is `curated` counts, not
mechanism. Agents are used only where a question needs them (a niche, a
clonal expansion, a tumour).

### 2.7 What "reads the genome" means here, honestly

The gene identities are real (WBcel235 and GRCh38 annotations; the repo
already compiles C. elegans chromosome III). The reader mask is real
(measured TF presence per cell from the Ma atlas; ENCODE domains for
human). The decision rules are literature-derived and cited (`curated`,
`experimental`), not derived from sequence. The step from sequence to rule
stays `inferred` or UNKNOWN until predictive models (AlphaGenome, Virtual
Cell) earn a place, and the report will say so per level.

---

## Part 3. Plan and proposal

Each step ends with tests and lint green, a docs/PROGRESS.md entry, and a
result file where real data was touched. Steps are ordered so each is
useful alone.

### Step 0. Truth data, distilled (small)

- Fetch the WormWeb lineage JSON; convert to a compact
  `data/results/celegans_lineage_cells.json` (name, parent, born, divides,
  dies, fate) with licence and source recorded; ~150 KB; committed.
- Distil Packer 2019 annotation to lineage → cell type → time (small
  table); distil Richards 2013 Table 2 to per-cell cycle length mean/SD;
  delete raw files. `genomeos data status` lists them.
- Tests: tree invariants (959 somatic leaves alive in the adult, 131
  deaths, 558 at hatching within tolerance).

### Step 1. BioIR and BioLang v0.3 core

- IR: `Organism`, `Stage`, `Timer`, `Signal`, `Decision`, `Population`,
  `Domain` (promoted), with `to_dict`/`from_dict` round-trip.
- Parser: new block kinds `organism`, `stage`, `timer`, `signal`,
  `decision`; `observe`, `assert`, `compare` lines inside `organism`.
- `genomeos/std/development.bio`: the five decision verbs, generic cell
  types (Zygote, Blastomere, Progenitor, PostMitotic), cited.
- Docs: docs/BIOLANG-v0.3.md, docs/BIOIR-v0.3.md.
- Tests: parse, round-trip, reference checks, error messages.

### Step 2. Body runtime and the event interpreter

- `genomeos/runtime/body.py`: Body, CellAgent, Population, Reader,
  discrete-event loop, coarsening and expansion, LineageDiff, organism-level
  uncertainty.
- The interpreter for `Module.events` and `Decision`s (the missing piece);
  `cell.py` ageing keeps working and gains a `from_module` path.
- Replace the hardcoded logic in `celegans.py` with a BioLang module run
  by the Body; fix the EMS lineage label; keep the current test passing.
- CLI: `genomeos grow FILE --until MIN --observe --compare`.
- Tests: one → two → four with correct names and times; coarsening
  preserves counts; diff against truth is deterministic.

### Step 3. Milestone A: C. elegans to the 28-cell founders from factors

- `celegans/early_fate.bio`: the 17-node Boolean logic (Genetics 2025)
  plus Notch (GLP-1/APX-1) and Wnt/POP-1 asymmetry as cited rules;
  bootstrap factors in the organism block; cycle timers from Richards.
- Score: founder identities and division order match Sulston; timing
  error against Richards; consistency with Cole 2024 founder codes.
- Web: Organism tab (count vs time, lineage tree, fate map, diff panel),
  modelled on the existing develop endpoint.

### Step 4. Milestone B: embryo to hatching (558 cells, 113 deaths)

- Lineage-level timers per founder from Richards; fates per terminal cell
  from Packer and Ma (reader masks as `cell_type expresses` lists derived
  from the TF atlas, cited per TF); deaths as cited decisions.
- Populations used for committed branches; agents where fates split.
- Score: fate accuracy over 1,092 terminal cells, deaths matched, count
  curve error, uncertainty report per level.
- Knockout check against Digital Development for a handful of the 204
  genes (fate changes predicted vs observed).

### Step 5. Milestone C: larva to adult, then maintenance

- Post-embryonic blast cells, moult timers L1–L4, adult 959; germline as a
  Population; adult somatic cells post-mitotic; ageing runtime as the
  maintenance phase driven by bio.std.ageing through the same interpreter.

### Step 6. Milestone D: human body as populations

- `human/body.bio`: Carnegie stages, germ layers, tissue populations from
  Hatton (distilled), turnover from Sender & Milo (distilled), tempo scalar
  with citation; grow to adult and run maintenance for years; report says
  what is mechanism and what is counts.
- HG002 remains the test human: variants that hit rule genes flow into
  the decisions through the existing variant-effect path.

### Estimated size

| step | new code | tests | data |
|---|---|---|---|
| 0 | ~200 lines | 3 | ~200 KB committed |
| 1 | ~500 lines | 10 | none |
| 2 | ~700 lines | 10 | none |
| 3 | ~300 lines + .bio | 6 | none |
| 4 | ~400 lines + .bio | 6 | ~300 KB distilled |
| 5 | ~200 lines + .bio | 4 | none |
| 6 | ~300 lines + .bio | 5 | ~200 KB distilled |

Steps 0–3 are one continuous run of work and produce the first measurable
result (founders from factors, scored). Steps 4–6 each stand alone.

### Decisions needed before starting

1. Names: `organism`, `stage`, `timer`, `signal`, `decision` as the v0.3
   block kinds, and `genomeos grow` as the command. Alternative: fold
   `decision` into `event` with an `action:` field.
2. Order: C. elegans first (steps 0–5), human populations after (step 6).
   Alternative: run step 6 in parallel after step 2, since it needs only
   populations.
3. Scope of step 4 fates: start with tissue-level fates (neuron, muscle,
   hypodermis, intestine, pharynx, reproductive, death) rather than the
   502 subtypes; subtypes later from Packer.

Sources are listed in section 1; datasets and licences in 1.2 and 1.3.

---

## Built (2026-09-11)

| step | state | where |
|---|---|---|
| 0 truth data | done | `genomeos data distil --only celegans_lineage` and `celegans_packer2019`; data/results/celegans_lineage_cells.json, data/results/celegans_*.json |
| 1 BioIR and BioLang v0.3 | done | genomeos/ir/model.py, genomeos/lang/parser.py, genomeos/std/development.bio, docs/BIOLANG-v0.3.md |
| 2 Body runtime | done | genomeos/runtime/body.py, genomeos/organism/diff.py, `genomeos grow` |
| 3 founders from factors | done | data/organisms/celegans/founders.bio (PAR, SKN-1, PIE-1, PAL-1, POP-1, Wnt, Notch) |
| 4 embryo to hatching | done | data/organisms/celegans/embryo.bio: 1,439 cells, 555 fates, 110 deaths reproduced |
| 5 larva to adult | done | lineage_larva.bio: 961 cells alive in the adult, 131 deaths |
| 6 human body as populations | done | data/organisms/human/body.bio + generated tissues.bio from Sender & Milo 2021 (`genomeos data distil --only human_cell_turnover`); 1.8e12 cells at birth, 2.84e13 at 20 years, 3.19e11 replaced per day; low confidence by construction |

Measured against the reference (`genomeos grow data/organisms/celegans/embryo.bio --until 6000 --compare`):
cells born 2,183/2,183, parent mismatches 0, terminal fates 961/961, deaths 131/131, division timing median 12 min
(embryo median 8 min, p90 25 min). The alive-cell curve runs under the reference between 350 and 500 min because
per-generation mean timers smooth the real spread; the program's own `assert` flags it.

Time axis (`genomeos data distil --only celegans_time_axis`): one reference minute is about 0.6 to 0.7 real minutes
(Ma 2021 frames: births 0.58, divisions 0.61; Packer 2019 times: 0.70), so the WormWeb axis runs slow by about 1.5;
topology, fates and deaths are unaffected.

Knockouts against the published set (`genomeos data distil --only celegans_digital_development`): 7 of 11 modelled
founder transformations of Du et al. 2014 reproduced; the misses were the PAR polarity genes, whose rules were not written.
(11 of 11 since 2026-09-14, when they were: see the PAR section at the end.)

Reader (`genomeos data distil --only celegans_tf_atlas`): the Ma 2021 transcription-factor atlas as one `express`
decision per cell, imported by the worm program; textbook factors predict their tissues (ELT-2 intestine 0.96, HLH-1
body-wall muscle 0.94, PHA-4 pharynx 0.85), which is the atlas and the lineage checking each other.

Design (`genomeos grow FILE --design`): BioForge over organism programs; the C. elegans and haematopoiesis design
files answer their questions with the expected factor (POP-1, PIE-1, GFI1), each returned as a predicted experiment.

Haematopoiesis (`genomeos grow data/organisms/human/haematopoiesis.bio --until "3 yr"`): the first human
mechanism module; every compartment gated by the factors the genetics showed necessary, flows solved from the
measured outputs, pools within 15% of design after three years, 2.8e11 cells replaced per day, nine knockouts
reproducing their published phenotypes.

Human (`genomeos grow data/organisms/human/body.bio --until "20 yr"`): 2-cell at day 1, 8-cell at day 3,
1.8e12 cells at birth, 2.87e13 at 18 years, 2.84e13 at 20 years, 3.19e11 cells replaced per day against the
published 3.3e11; erythrocytes 2.5e13, neurons 1.2e11, cardiomyocytes 3.2e9. The uncertainty report reads
"low" at every level, which is the point: the human program is counts and shares with citations, and the
language now has a place to put mechanism when it is known.

---

## Fate as context, tested (2026-09-14, genomeos-d1)

The working model: a cell's fate is a function of position, signals, current state, genome, epigenome and
developmental time, read by a regulatory network that then reinforces itself until the identity is hard to
reverse. Four parts of that model were tested against the worm, where every cell is named; the results are
in `data/results/celegans_fate_rules.json` and `data/results/celegans_commitment.json`
(`scripts/celegans_fates.py`, `scripts/celegans_commitment.py`; modules `organism/atlas_levels.py`,
`fate_rules.py`, `commitment.py`, `lateral.py`). The Ma 2021 archive was streamed again and kept as levels
over time (peak, mean over the cell's life, exposure in minutes at the factor's maximum, onset frame), a
9.5 MB local table under `data/knowledge/celegans`.

### 1. Terminal fates from measured factors

Scored on the 555 embryonic terminal cells (the atlas ends at the bean stage; the 406 larval cells keep the
lookup). Before: the observed lookup, 555/555 by construction. The rules take precedence over it in
`embryo_factors.bio` (`fates.bio`, each rule with a `priority` above the lookup's 0 under
`regime { fates: first }`, gated on a terminal cell type and the embryonic stages, so the lineage still says which cells stop dividing and when, and the factors say what
they become).

| rule set, read | decided by factors | right | wrong | score with lookup fallback | lookup-free (factors, else sublineage majority) |
|---|---|---|---|---|---|
| sublineage majority, no factors | 555 | 364 | 191 | 65.6% | 364 |
| textbook, instantaneous threshold | 225 | 164 | 61 | 89.0% | |
| textbook, integrated along the lineage | 216 | 178 | 38 | 93.2% | 405 |
| learned, instantaneous, sublineage held out | 378 | 280 | 98 | 82.3% | |
| learned, integrated, sublineage held out | 488 | 345 | 143 | 74.2% | |
| textbook then learned neuron rules, integrated, held out | 449 | 365 | 84 | 84.9% | 408 |
| the program (`embryo_factors.bio`, neuron rules fitted on all cells), through the Body and LineageDiff | 437 | 378 | 59 | 496/555 = 89.4% | |

To the adult the program scores 902/961 (the larval cells are the lookup's). Per tissue, in the program:
intestine 34/34 and hypodermis 69/69 kept, and every intestinal cell is called by ELT-2 or ELT-7 exposure
alone; neuron 212/226 (200 called by factors), body muscle 113/122 (78 called; HLH-1 and UNC-120 are
body-wall muscle, pharyngeal muscle is not called and 9 muscle cells go to neuron rules); the factors do
worse than the lookup for sheath glia 6/22 and socket glia 9/18 (lost to hypodermis and neuron: glia are
sisters of sensory neurons and carry ELT-1, LIN-26 and NHR-25 and the neuronal factors of their
lineage), coelomocytes 0/4 (HLH-1-positive MS cells, called muscle), valve 6/8, rectal 1/2, excretory 4/5.
Pharyngeal marginal (8/9), epithelium (12/13), gland (5/5) and the rest have no rule of their own and lose
at most a cell to another tissue's rule.

Where factors decide, they beat the lineage alone: on the 216 cells the textbook rules decide, the rules are
right 178 times and the sublineage majority 137; on the 449 the held-out hybrid decides, 365 against 321.
The factors carry per-cell information beyond composition: with tissues shuffled within each sublineage the
held-out learned rules reach 18 to 29% of cells, against 50% (instantaneous) and 62% (integrated) on the
real labels. Learned rules fail completely for a tissue only one founder makes: holding E out leaves no
intestine to learn from, and the integrated list calls all 34 intestinal cells coelomocytes.

Time integration, tested at every threshold (textbook rules; fraction of the factor's maximum):

| fraction | peak in own life: right / wrong | mean over own life | exposure integrated along the lineage |
|---|---|---|---|
| 0.1 | 152 / 191 | 167 / 41 | 179 / 80 |
| 0.2 | 163 / 62 | 159 / 29 | 178 / 38 |
| 0.3 | 147 / 35 | 142 / 16 | 161 / 20 |
| 0.5 | 118 / 9 | 119 / 4 | 136 / 8 |

The instantaneous peak is the worst read at every threshold. The two time-averaged reads beat it; the
lineage-integrated read makes the most correct calls at every threshold and the mean over the cell's own
life makes the fewest errors. Cells behave as if they read sustained exposure, not a moment above a line;
which window is integrated is not settled by these data.

### 2. Commitment as a ratchet: not visible in the atlas (negative)

Over 755 cells with complete lifetimes, in five birth-time bins (50 min each), against 200 shuffled-time
permutations:

- **States do not become fewer and more separated.** Factors per cell rise from 10 to 79. Nearest-neighbour
  distance relative to the spread (discreteness) rises, 0.38 to 0.71 (p 0.005 in the direction opposite to
  commitment). Discreteness in excess of a curveball randomisation that keeps every cell's count and every
  factor's frequency is flat, -0.14 to -0.12 (p 0.49 and 0.52): cells are more clustered than their margins
  imply at every time, and no more so later. Grouped by the tissue their descendants mostly make, the
  silhouette is negative throughout (-0.15 to -0.12, p 0.25): until the bean stage the factor state tracks
  lineage history more than eventual fate.
- **Competing programmes are exclusive from the start, and do not become more so.** Among cells with an
  intestine, body-muscle, hypodermis or neuron programme factor on, the share with two programmes goes 0, 0,
  2.5%, 11% (relative to the curveball expectation 0, 0.10, 0.37); a declining trend has p 0.995.
  Programme factors switch on in cells already restricted (Cole 2024: most cells single-fated by the
  102-cell stage), so the data cannot see a programme suppressing its competitors.
- **The ratchet on mother to daughter pairs is not measurable.** 88 pairs with a programme in the mother:
  persistence 0.71 (reporter perdurance makes this an upper bound), switching 9%, neither trends with time
  (p 0.27, 0.48).

The atlas is the wrong instrument for hysteresis: it has no perturbation, and reporter protein perdures.
The published evidence for commitment in the worm is perturbational (below), and that is what the construct
proposal cites.

### 3. Competence windows: not separable from lineage in an invariant lineage (negative as a test)

For every de novo onset of a factor (present in a cell, absent in its tracked mother), does the onset time
bin say something about the fate beyond the factor and the founder sublineage? Over all factors, 15,353
onsets in 931 strata: 0.078 bits beyond the within-stratum permutation null (p 0.005); with the generation
also fixed, 1,081 onsets: 0.041 bits (p 0.005) on descendant composition and nothing on the dominant tissue
(p 0.14); for the textbook programme factors only, 124 onsets: 0.02 bits (p 0.16). Time carries a little
information at fixed factor, sublineage and depth, but in the worm the time at fixed depth is branch
identity (branches cycle at different speeds), so this does not separate competence from lineage history.
The evidence for competence windows in the worm is experimental and not in these data: ectopic HLH-1
converts blastomeres to muscle only in early embryos (Fukushige & Krause 2005, Development 132:1795), and
MES-2 (Polycomb) ends that plasticity around the start of gastrulation (Yuzyuk et al. 2009, Dev Cell
16:699). That is enough to justify a construct; it is not a result of this analysis.

### 4. Lateral inhibition and noise

Of 327 sister pairs with complete lives, 314 are anterior/posterior: the anterior sister carries more POP-1
in 255 of 271 measurable pairs (median 1.4-fold, sign test z 14.5; Lin et al. 1998), which is instruction,
not noise; the asymmetry is no larger in the 124 pairs whose descendants' tissues differ (p 0.20). Nine
left/right sister pairs and 107 ABpl/ABpr bilateral homologues are the equivalent case: 5 and 18 of them
differ in dominant tissue, with a larger factor distance than the concordant ones (0.32 against 0.17; 0.59
against 0.51). But the measurement floor swamps it: the 24 factors followed by two reporter strains
disagree on presence in the same cell 68% of the time, more than sisters disagree on the same factors
(45%) or homologues (53%). The Notch target REF-1 (Neves & Priess 2005) is present in 21 cells; sisters
discordant on it are no more often fate-discordant than for factors of the same frequency (rank 17 of 43).
Negative: no lateral-inhibition signature is detectable at this resolution, and the embryonic Notch
decisions of the worm are inductions by named neighbours, which the program already writes.

Where equivalent cells do diverge in the worm, it is noise plus selection: Z1.ppp and Z4.aaa each become the
anchor cell about half the time (Kimble & Hirsh 1979; Seydoux & Greenwald 1989), while the reference names
Z1.ppp and scores only that. `organism/lateral.py` holds the reference behaviour a contact runtime must meet
(Collier et al. 1996): two equal cells with a deterministic rule never diverge; with noise 198 of 200 pairs
diverge and the first cell wins 51%. A fate rule deterministic in position and signal is wrong exactly
there.

### Constructs proposed to the engine (genomeos-c1)

Sent for docs/BIOLANG-v0.4 (the IR and runtime are genomeos-c1's):

- `commitment NAME { programme: factors; establish: condition (level or integrated exposure); locks:
  cell_type; maintain: condition; excludes: factors or programmes; hysteresis: enter X, leave Y; release:
  never | after N min below leave | experiment; inherit: daughters | none }`. Semantics: evaluated before
  `differentiate`; once established, a differentiate decision to a type outside the locked type's
  descendants is refused (the Body let the last applying decision win, against the v0.3 specification, so
  `fates.bio` first listed its rules bottom to top; since genomeos-c1's `regime { fates: first }` and
  decision `priority` (a20c922) the rules carry explicit priorities and the order is not load-bearing, but a
  later decision point can still revise a fate until `commitment` refuses it); signal-driven re-decisions no longer change the fate; excluded factors are masked in the
  cell; the state is inherited as stated mechanism.
- `competence NAME { when: context; allows: cell types; closes: at time | on commitment | after generation
  }`. A differentiate decision whose target is outside the cell's open competence is refused and reported as
  outside competence, not UNKNOWN.
- Numeric factor levels and time reads in `when:`: `ELT-2.exposure >= 30` (minutes at the maximum level,
  accumulated along the lineage between events), `ELT-2.mean >= 0.2`, and `ratio(A, B)`; `noise:` on a
  threshold, drawn per cell from the seeded generator.
- For contacts: neighbours per cell from grid adjacency or an imported time-resolved contact table; contact
  signals evaluated against current neighbours (not named senders) with an amount (ligand level times
  contact); a per-cell network stepped synchronously across coupled neighbours between events with seeded
  noise; daughters placed along the division axis their names say; replicate runs and a diff that scores an
  equivalence group ("exactly one of Z1.ppp, Z4.aaa is the anchor cell") rather than one fixed outcome.

---

## The founders decided by contact, not by name (2026-09-14, genomeos-d2)

The constructs above landed in the engine on 2026-09-14 (contacts 21d3dbc; `commitment` and
`competence` 126e65b), so the worm's founder program was rewritten to use them.
`data/organisms/celegans/founders_contacts.bio` replaces founders.bio's three named signals
(`from: cell = P2; to: cell = ABp`) with contact signals that name no cell: each receiver reads the
ligand summed over the cells touching it (`mode: contact; reads: amount`), who touches whom comes
from `contacts_embryo.tsv`, which cells present a ligand is an `express` decision carrying its
measurement, and which cells can hear a Delta is the maternal GLP-1 inherited by AB.
`embryo_contacts.bio` is `embryo_factors.bio` with that layer swapped in and nothing else changed.
Measured by `scripts/celegans_contacts.py` into `data/results/celegans_contacts.json`:

| test | result |
|---|---|
| the eight founder identities (ABa, ABp, EMS, MS, E, C, D, P4) | **8 of 8**, from contact and inheritance alone |
| **the rearrangement**: ABa put where ABp is, in the contact table only | the fates **swap**: ABa becomes ABpPrecursor and ABp ABaPrecursor |
| control: the same program with no contact table, so nothing touches anything | nothing is induced: ABp is ABa-like, E takes the MS fate |
| the eight published knockouts (apx-1, glp-1, mom-2, mom-5, pie-1, skn-1, pop-1, pal-1) | **8 of 8** reproduced, as with named senders |
| **the whole embryo against Sulston to 800 min** | **496 of 555 terminal fates, 89.4% — exactly the number before**; 1,439/1,439 cells born, 0 parent mismatches, deaths 110/110, timing median 14 min, the same five confusions |
| cells differing between the two programs, by fate, terminal name and birth time | **1 of 1,439**: EMS itself, which is now typed EMSPolarised rather than EMSPrecursor |

So the answer to "what do contacts change about the lineage score" is **nothing, and that is the
result**: 1,438 of 1,439 cells are identical to the last digit, and the founder layer that produces
them no longer contains a statement of the embryo's geometry. What changes is what the program can
now be wrong about. Named senders make ABp's fate unfalsifiable — the rule says ABp, so ABp gets it.
Read from contact, the same rule predicts the blastomere rearrangement of Priess & Thomson 1987
(Cell 48:241) and Hutter & Schnabel 1994 (Development 120:2051): move ABa into ABp's place and it
takes ABp's fate. That experiment is a three-line edit of a data file here and cannot be written at
all in founders.bio.

**What the table is, honestly.** `contacts_embryo.tsv` is curated from the published geometry of the
cleavage-stage embryo (the 4-cell rhomboid in which ABp touches P2 and ABa does not; MS reaching
ABalp and ABara at the 12-cell stage), not measured contact areas. No time-resolved contact table for
C. elegans is held locally and the Ma 2021 atlas carries expression per cell per minute, not
positions, so every area in it is 1.0 or 0.5 and means "touching" or "a minor face". A measured table
would replace the file without a line changing in any program that reads it. That is the next
measurement area E should fetch.

**Two runtime changes area E would like (reported, not made; the engine is another lane's).**

1. **A contact amount is re-read when a neighbour is born, but not when one leaves.** *Corrected the
   same night.* This section first said the Body never re-reads, and the program carried two
   presence signals to force it; both were wrong and both are gone. `Body._add` already ends with
   `for name in self.neighbours(c): self._resolve(...)` whenever there are amount signals and no
   per-cell network, so a receiver born before its sender (ABp at 28 min, P2 at 36 to 42 min) does
   hear it — provided the contact table already lists the pair at the moment the sender is born.
   What is missing is the other half: a cell dividing or dying does not make its old neighbours read
   again, and with a `cell_network` declared the birth hook is skipped entirely in favour of the
   cadence. Requested: re-resolve `neighbours(...)` in `_divide` and `_cull` as well, and keep the
   birth hook when a network is present.

   The way this was found is the more useful half. The contact table's snapshots were first keyed to
   the division times of one seeded run, and `run_experiment` runs every timer at its mean: P1
   divides at 42 min with the seed and 36 min at the mean, so under the mean the 4-cell snapshot had
   not started when P2 was born, P2 had no neighbours, and **three knockouts that had been
   reproduced were silently missed** (apx-1, glp-1 and pie-1, all through ABp). The fix is in the
   table, not the program: each snapshot is now taken at an AB division, which has no measured
   spread, and lists the union of the configurations holding until the next one, with aliveness
   filtering the rest. A contact table keyed to a time a timer can move is a trap, and nothing in
   the language warns about it — a `contacts:` table whose rows never match a living pair should be
   an error, not silence.
2. **Two decisions may share an id, and the second becomes unreachable.** `Body.__init__` builds
   `named = {x.id for v in self._by_cell.values() for x in v}` and then `_general = [d for d in
   module.decisions if d.id not in named]`, so a decision keyed on `cell_type` is dropped from every
   candidate list when *another* decision somewhere in the program happens to have the same id and a
   `cell` clause. This cost an hour: `div_EMS` here collided with `div_EMS` in the generated
   `lineage_embryo.bio` and never fired, silently and with no warning. Requested: key that set on
   object identity, and warn (or refuse) when two decisions share an id.

A third, smaller one: `cell_network` re-resolves every live cell on its cadence, and re-resolving a
cell lets the next-highest-precedence `differentiate` fire, because `_pick_fate` skips decisions
already in `c.fired`. With `regime { fates: first }` that overwrites the intended winner — measured
here as 496 → 482 of 555 terminal fates with a 2-minute cadence and no other change, with
Hypodermis → Neuron replacing Neuron → Hypodermis in the confusions. A cadence should not be able to
change a fate the program already settled.

## AC/VU as an equivalence group: 48% Z1.ppp, and the ways an order gets smuggled in (2026-09-14, genomeos-d2)

`data/organisms/celegans/acvu.bio` is the one decision in the worm the invariant lineage cannot
predict. Z1.ppp and Z4.aaa are interchangeable: in any one animal exactly one becomes the anchor
cell, and across animals it is about half and half (Kimble & Hirsh 1979, Dev Biol 70:396; Kimble
1981, Dev Biol 87:286; Seydoux & Greenwald 1989, Cell 57:1237). The two cells in the program are
identical — same rules, same initial levels, same circuit (Collier et al. 1996 on LIN-12/LAG-2);
they touch (`contacts_acvu.tsv`), each reads the other's Delta as an amount at the contact, and only
the seed differs. Nothing names the winner. 200 seeds, `scripts/celegans_acvu.py`,
`data/results/celegans_acvu.json`:

| condition | result |
|---|---|
| the program as written | **199 of 200** runs give exactly one anchor cell; **Z1.ppp 96 (48.2%)**, Z4.aaa 103 |
| noise set to zero | **0 of 200** diverge; the two cells' Delta identical to the last digit (max gap 0.0) |
| each division's daughters created in the opposite order | the **same winner in 200 of 200 seeds** |
| against the reference lineage, which records Z1.ppp as `gon herm anch` | **96 of 200 = 48%**, and that is the ceiling for an honest program |

The engine's own two-cell gate gives 104/96 on the same circuit; this is the same result inside the
worm's lineage, at its measured division times, with the four cells of the group present.

**The last row is the point.** The worm program (`lineage_larva.bio`) names Z1.ppp the anchor cell
and scores 100% on that cell, because the reference recorded one animal. A program that reproduces
the biology scores 48%. The AC/VU cell is the one place among the 1,092 terminal cells where the
lineage score is not the right instrument, and a fate rule deterministic in position and signal is
wrong exactly there.

**Two ways the coin flip was silently loaded, both found by running the arms rather than by
reading the code.**

1. **A daughter inherits its mother's whole network state.** Z1.pp divides 56 minutes before Z4.aa.
   List Z1.ppp as touching Z4.aa — which is true geometry — and Z4.aa spends those 56 minutes
   receiving Delta, so Z4.aaa is born with LIN-12 already active and loses: **Z1.ppp wins 100 of 100
   seeds**. The contact table therefore lists the surface between the two members of the group only
   from the moment both exist, and says so. The general problem is that there is no way in the
   language to say a species resets at division, and no warning when a signal received by a mother
   decides a contest between her daughter and someone else's.
2. **An instantaneous level read at an arbitrary time is not a fate.** Before they touch, both cells
   sit at maximum Delta, so `Dp >= 0.6` fires for both at birth; during the symmetric transient both
   sit at maximum Notch, so `Np >= 0.6` fires for both. The fate reads are gated on the L3 stage and
   the run goes to 4,200 min because the circuit takes 10 to 17 hours of model time to resolve.
   What is missing is the sustained read of §7.2a (`Dp.exposure`, `Dp.mean`), which did not land:
   with it the rule would say "high Delta for long enough", which is what the cell does.

**The experiment that defines the group, run.** Kimble 1981's ablations are what make this an
equivalence group rather than a coin flip between two cells that were going to differ anyway, and
they are now four `experiment` blocks in the program and four arms in the script, scored over 100
seeds each because "every time" is a claim about every run:

| ablation | result |
|---|---|
| Z1.ppp | Z4.aaa is the anchor cell in **100 of 100** runs |
| Z4.aaa | Z1.ppp is the anchor cell in **100 of 100** runs |
| both | **no anchor cell in any of 100** runs; the two flanking cells stay ventral uterine |
| Z1.ppa (a flanking cell: the control) | the decision is untouched — 99 of 100 with one anchor, split 43 to 56 |

Nothing about the ablations is stated in the fate rules: the survivor wins because there is no
longer a neighbour presenting Delta at its contact, which is the same mechanism that makes the
wild-type decision a coin flip. `ablate` is an environment value the experiments set, so the
wild-type program never reads it.

One more control, on the membership: the equivalence group is two cells, not four. Z1.ppa and
Z4.aap flank the pair and are ventral uterine precursors in every animal (Kimble 1981), so only
Z1.ppp and Z4.aaa present LAG-2 here. Let all four present it and the four-cell chain settles into
the alternating pattern lateral inhibition gives on a line — the two outer cells keep Delta, both
central cells lose — and **0 of 100 runs produce an anchor cell**. The restriction is load-bearing
and is stated, with its ablation citation, rather than assumed.

One more language note: which two cells are in the group is written as `cell = Z1ppp|Z4aaa` on the
anchor decision. `competence` cannot say it. A window governs the cells its `when` matches and
leaves every other cell competent, so "only these two may take this fate" has to be written inside
out, as a window closed at time zero for everyone else. A `competence` with an `only:` sense, or a
`restricts:` clause, would say it directly.

## Glia from factors: sheath partly, socket not at all (2026-09-14, genomeos-d2)

Glia were the clearest failure of the factor rules: of 555 embryonic terminal cells, 22 are sheath
glia and 18 socket glia, and the rules called 6 and 9 of them. The rest went to the hypodermis and
neuron rules, because a glial cell is the sister of a sensory neuron, carries that neuron's factors,
and carries the ectodermal factors (ELT-1, LIN-26, NHR-25) the hypodermis rules read.
`scripts/celegans_glia.py` asks whether anything in the Ma 2021 atlas separates them;
`data/results/celegans_glia.json` holds the answer.

**The question is answerable.** All 40 glia are tracked in the atlas and so are all 40 sisters, and
under every read of the atlas **no glial cell has the same factor set as a non-glial cell**. Nothing
rules separability out in principle, so a negative here is about the factors, not about coverage.

**Sheath glia: yes, partly, and by the gene the literature names.** PROS-1/Prospero is the published
sheath-glia regulator (Wallace et al. 2016, Development 143:3016). Alone it is worth nothing
(precision 0.27: it claims 22 cells, 8 of them neurons and 4 coelomocytes). Conjoined with the
factors that say which class the cell is in, it sharpens:

| rule, on the read the program's reader gives a cell | claims | sheath | precision |
|---|---|---|---|
| PROS-1 | 22 | 6 | 0.27 |
| PROS-1 ∧ NHR-25 (the non-neuronal ectodermal class) | 10 | 6 | 0.60 |
| PROS-1 ∧ SOX-2 (the neural class) | 6 | 5 | 0.83 |
| **PROS-1 ∧ NHR-25 ∧ SOX-2** | **5** | **5** | **1.00** |

The third rule is now in the program (`fate_rules.CITED_GLIA_RULES` → `fates.bio`, priority between
the muscle and hypodermis rules, read instantaneously because the lineage-integrated read carries no
sheath signal at all). Measured through the Body and LineageDiff: **sheath glia 6/22 → 11/22**, with
hypodermis still 69/69, neuron 212/226 and muscle 113/122 — nothing is lost elsewhere — and the
program as a whole **496 → 501 of 555 (89.4% → 90.3%)**, 902 → 907 of 961 to the adult. The two
other published sheath factors in the atlas do worse: MLS-2 (Yoshimura et al. 2008) precision 0.27,
HLH-17 (McMiller & Johnson 2005) precision 0.045 — it is present in 264 of 555 cells.

**Socket glia: no. These cells are not separable by the factors we have.** EGL-13/SoxD, the
published socket factor (Feng et al. 2013), claims 34 cells and gets 1 right. The best socket rule
that exists at all, SOX-2 ∧ PAG-3, claims 22 cells and gets 9 — it takes 12 neurons to win 9 socket
glia, and putting it in the program **costs 4 terminal fates net**. So socket glia have no rule, and
that is measured rather than forgotten.

**The ceiling, exhaustively.** Over all 250 factors with support, all 25,867 pairs and every triple
— the whole hypothesis class a BioLang `when:` clause can express on measured factor presence —
scored by what the rule would be worth to the program:

| target | best single | best pair | best triple |
|---|---|---|---|
| sheath | 0 fates | +4 (PROS-1 ∧ SOX-2) | **+6** (ALY-1 ∧ EOR-1 ∧ PROS-1) |
| socket | 0 fates | +2 (MLS-2 ∧ MML-1) | **+3** (CEH-32 ∧ HAM-2 ∧ HLH-4) |

Every one of those is chosen on the same 555 cells it is scored on, so they are upper bounds, and
the best sheath rules all contain PROS-1 while the best socket rules contain nothing anyone has
named. The rule adopted is worth +5 of the +6 ceiling and is cited rather than fitted to the maximum.

**And the fitted rules do not generalise.** Rules fitted on seven founder sublineages at a median
in-sample precision of 1.0 score, on the eighth: glia precision **0.345** (peak read, 10 right of 29
claimed), 0.30 on the mean read, 0.25 on the integrated read; sheath 0.214, socket 0.286. Against a
base rate of 7.2% that is four to nine times chance — the signal is real — and it is nowhere near
the 0.8 precision the program requires of a rule.

**Why, as far as the data can say.** A glial cell differs from its sister by a median of **39
factors**, two sisters of the same tissue by **35**, and two sisters of different tissues by **40**:
the sister difference is almost all measurement and lineage, and carries little about the fate. The
atlas's own floor is in the same place — the 23 factors followed by two reporter strains disagree on
presence in the same cell 20.7% of the time (median per factor). A rule that has to separate two
cells born from one mother a few minutes apart is working below that floor. What would separate them
is not another transcription factor at this resolution: it is the two things this atlas does not
have — the ligand a glial cell receives from its neighbours, and the terminal genes (the glial
secretome PROS-1 drives) that come after the bean stage where the atlas stops.

One consequence in another lane's file, recorded here. `plasticity.bio`'s terminal arm asserts that
no cell already differentiated at 700 min is lost when HLH-1 arrives, with the bound for each type
being the number of such cells — and the hypodermis bound was 99 because five sheath glia were
counted as hypodermis. With the sheath rule the number is 94, so the bound is 94 and a fourth
assert (`type GliaSheath at 800 min >= 11`) now carries the five cells that moved. The arm is
stronger, not weaker: it protects the same cells under their right names.

## PAR polarity: the last four misses against the published knockout set (2026-09-14, genomeos-d2)

The PAR genes were the only misses the worm had against Digital Development (Du et al. 2014, Cell
156:359): 7 of 11 modelled founder transformations reproduced, and all four misses were par-2 and
par-3. The reason was structural rather than biological — the first two divisions segregated the
maternal factors unconditionally, so a program in which PAR-2 or PAR-3 is absent divided exactly as
the wild type did. The segregation now depends on the domain that does it, which is one extra
decision per division in each founder layer:

```
div_P0            when: cell = P0, PAR-3 = present   asymmetric: PIE-1, SKN-1, PAL-1, PAR-2 -> P1 ...
div_P0_no_par3    when: cell = P0                    asymmetric: PIE-1, PAR-2 -> P1 ...
div_P1            when: cell = P1, PAR-2 = present   asymmetric: PIE-1 -> P2
div_P1_no_par2    when: cell = P1                    asymmetric: PIE-1 -> none
```

with two more statements of the same kind: the EMS fate rule reads `cell = AB|EMS|P2` rather than
`EMS|P2` (a cell with SKN-1 and no PIE-1 takes the EMS fate, whichever cell it is), and the MOM-2
ligand is presented by P2 only while P2 has PIE-1, as APX-1 already was.

| knockout | observed (Du et al. 2014) | program |
|---|---|---|
| par-3 | AB adopts EMS | AB inherits SKN-1 and PAL-1 and becomes EMSPrecursor; the P lineage is untouched, so ABp is still induced and E is still endoderm |
| par-2 | ABp adopts ABa; P2 adopts EMS; E adopts MS | neither daughter of P1 keeps PIE-1, so P2 becomes EMSPrecursor, presents neither APX-1 nor MOM-2, and ABp and E lose their inductions in turn — **one lost identity, three transformations** |

**Measured (`scripts/celegans_par.py`, `data/results/celegans_digital_development.json`, no network
call — the table is the one already distilled): 7 of 11 → 11 of 11 modelled transformations
reproduced, 0 missed, in both the named-sender and the contact programs.** The 14 changes among the
8-cell AB granddaughters stay out of scope: the program has no type for an ABalp identity.

Nothing about the wild type moves: PAR-2 and PAR-3 are maternal factors of the zygote, so the
conditional decision is the one that applies in every wild-type run — 1,439 of 1,439 cells, 501 of
555 terminal fates, deaths 110/110, unchanged. `mutants.bio` gains the two experiments, so the
claims are asserts the program checks on itself.

One honest extra: gating MOM-2 on PIE-1 also makes the program predict E adopting MS in *pie-1*,
which the published table does not list for that gene. It follows from the same mechanism (a P2
that has lost its germline identity signals with neither ligand), it is a prediction rather than a
miss, and it is the obvious thing to check against a pie-1 lineage.

A runtime wart found on the way, reported rather than patched: `asymmetric: X -> D` **creates** X in
the keeper when the mother did not have it (`factors[factor] = factors.get(factor, "present")`), so
in the par-2 run P3 and P4 regain the PIE-1 their mother lost. Nothing scored here depends on it,
but an asymmetric clause should divide what exists, not conjure it.

## The fate rules rewritten on the read the runtime computes (2026-09-15, genomeos-d3)

The engine landed the integrated reads of §7.2a on 2026-09-15 (8042a57): `F.exposure(cell)`,
`F.exposure(lineage)`, `F.mean(cell)` and `F.mean(lineage)` compile and run, a factor counts 1
while carried, the integral is banked at every decision point and clamped when the cell divides
or dies, and a read that does not name its window is a compile error. That is what area E's fate
rules had been waiting for, because until now they did not read a factor at all: `fates.bio`
tested `ELT-2_integrated = present`, and `ELT-2_integrated` was a per-cell lookup generated into
`exposure.bio` by `scripts/celegans_fates.py`, which called a factor present at **20% of its
largest path exposure over the finished run**. That is a statistic of the answer. No cell deciding
at 480 min can hold the maximum over 555 terminal cells, several of which are not yet born.

`exposure.bio` is deleted. `fates.bio` now reads `ELT-2.exposure(lineage) >= 15`, and the Body
integrates the reader's own presence calls. `scripts/celegans_fate_reads.py` is the generator and
the measurement; `data/results/celegans_fate_reads.json` holds it.

### The threshold, and why it is not a knob

**15 minutes, one AB cell cycle at 20 °C** (Sulston et al. 1983; Richards et al. 2013 measure 15
to 16 min for the first AB rounds). The rule says: *the path to this cell carried this factor
across at least one division* — not that a reporter crossed a line in one frame. The number is in
minutes, which is what the read is in, and a deciding cell holds it.

It is not fitted, and this is measured rather than asserted. On the cited textbook rules alone,
with nothing learned anywhere, every threshold from 1 to 60 min gives **the identical 164 cells
claimed and the identical 25 errors**; for the whole program through the Body, 1 to 30 min give
the identical program and the identical score. Above the plateau the score rises — 120 min scores
531 in sample against 522 — but only by claiming 62 fewer cells. That is buying a number with
abstention, and the plateau is where the threshold was taken.

`mean(lineage) >= 0.25` scores better than the shipped rule (533 in sample, 515 held out) and is
**not** shipped, for exactly that reason: it is not on a plateau. 0.20, 0.22 and 0.25 all give
different answers, so 0.25 is a point chosen because it scored well on the cells it is scored on.
It is in the result file as the alternative it is.

### The three scores, all through the Body and the existing LineageDiff

Same rule shape throughout — the cited textbook rules for intestine and body-wall muscle, the
cited sheath-glia rule, the cited hypodermis rules, then neuron rules learned from the atlas.
Only the neuron rules are fitted, so only they are held out: eight folds, one founder sublineage
held out at a time, each fold scored only on the lineage its rules never saw. The baseline is the
precomputed lookup **regenerated and run through the same Body in the same eight folds**, so the
comparison is like for like rather than against a number in an old table.

| read | embryonic fates | cells the factors decide | to the adult | held out | decided |
|---|---|---|---|---|---|
| instantaneous, `F = present` | 483 / 555 | 412 | 889 / 961 | 467 / 555 | 393 |
| **baseline**: precomputed `_integrated` lookup | **501 / 555** | 437 | **907 / 961** | **476 / 555** | 449 |
| **rewrite**: `F.exposure(lineage) >= 15` | **522 / 555** | 332 | **928 / 961** | **483 / 555** | 340 |
| `F.mean(lineage) >= 0.25` (not shipped) | 533 / 555 | 285 | 939 / 961 | 515 / 555 | 251 |
| `F.exposure(cell)`, any threshold | 555 / 555 | **5** | 961 / 961 | 555 / 555 | 5 |
| `F.mean(cell)`, any threshold | 555 / 555 | **5** | 961 / 961 | 555 / 555 | 5 |

So the rewrite is ahead of the lookup it replaces: **501 → 522 of 555 in sample (90.3% → 94.1%),
907 → 928 of 961 to the adult, and 476 → 483 of 555 held out (85.8% → 87.0%)**. Cells born
1,439/1,439, parent mismatches 0, deaths 110/110, unchanged.

**Three negatives in that table, and they are the interesting part.**

1. **The gain held out is a third of the gain in sample** (+7 against +21). The rewrite fits its
   own training set harder than the lookup did: its in-sample-to-held-out gap is 39 fates against
   the baseline's 25. A continuous read gives the neuron learner a larger effective hypothesis
   class than a binary lookup does, and it overfits accordingly. The honest number is 483.
2. **The rewrite claims 105 fewer cells** (332 against 437; 340 against 449 held out). It is more
   accurate on what it claims — 89% against 88% in sample, 88% against 82% held out — and it
   abstains more, and part of the score is the abstention. Scored without the lineage lookup as
   the fallback (factors, else the founder sublineage's majority tissue), the baseline reaches
   414 of 555 and the rewrite 377. **On that measure the rewrite is worse**, and it is worse for
   the reason this rewrite exists: the precomputed statistic normalised every factor by its own
   maximum over the run, so a factor with a small total exposure still called cells, and an
   absolute threshold in minutes cannot do that. That normalisation was carrying real
   information, and it was information no cell has.
3. **The cited sheath-glia rule stopped being load-bearing.** Under the precomputed read it was
   worth +5 terminal fates (sheath 6/22 → 11/22, genomeos-d2). Under the integrated read the
   program scores 522 with it and 522 without it, and sheath glia reach 20/22 either way — the
   hypodermis and neuron rules that used to steal them no longer fire on those cells, so there is
   nothing left for PROS-1 ∧ NHR-25 ∧ SOX-2 to win back. The rule is kept, because it is a cited
   claim about mechanism that the runtime now checks, and it is recorded here as worth nothing.

Per tissue, the rewrite against the baseline: muscle 113 → 121 of 122, sheath glia 11 → 20 of 22,
socket glia 9 → 16 of 18, intestine 34/34 and hypodermis 69/69 kept, coelomocytes still 0 of 4
(HLH-1-positive MS cells, called muscle), and neurons 212 → 208 of 226, which is the one tissue
that loses. The confusions are the same ones: 14 neurons to hypodermis, 4 to muscle.

### Does the ordering survive a real runtime? Yes. Does the cell window win anywhere? No

The measurement this replaces said the instantaneous peak makes 62 errors, exposure along the
lineage 38, and the mean over the cell's own life 29 — the instantaneous read worst at every
threshold from 0.1 to 0.5, with the mean over the cell's own life making the fewest errors. All
three were computed from the atlas's record of each cell's whole life.

**The instantaneous-versus-integrated ordering holds, and by more.** On the cited rules alone,
with nothing fitted, read by the Body at each cell's own decision point: the instantaneous read
claims 225 cells and gets **61 wrong**; `exposure(lineage)` claims 164 and gets **25 wrong**.
Through the whole program the same ordering: 483 of 555 instantaneous against 522 integrated in
sample, 467 against 483 held out. A cell that reads sustained carriage beats a cell that reads a
moment, in a runtime that integrates it honestly.

**The cell window wins nothing, because it is empty when the cell decides.** Every one of the 555
embryonic terminal cells gets **exactly one decision point, at its birth**, so at the instant its
fate is settled `t - born` is zero: `F.exposure(cell)` is 0 and `F.mean(cell)` is 0 for every
factor and every cell, and a rule written on either fires for nothing. The five cells still
decided in those two arms are the sheath-glia rule, which reads instantaneous presence.

That is the second statistic of the future that came out of this rewrite, and it is a bigger one
than the normalisation. "The mean over the cell's own life" was the best-scoring read in the
earlier measurement (29 errors), and it cannot be a fate rule: it is a summary of what the cell
went on to carry *after* it decided what to be. A cell can only integrate its own window if it
decides again later, which in this program it never does. So the honest answer to "which window
does a cell integrate over" is that this program can only ask about the lineage window, and the
question the earlier table appeared to settle was never posed to a deciding cell.

What would pose it: a re-decision cadence (a `cell_network` step re-resolves every live cell, so a
terminal cell born at 455 min would read its own window at 461, 467, …), or a `differentiate` that
is allowed to name a delay after birth. The first exists in the engine today and changes when
every fate is taken, which is a separate experiment and not this one.

### What the read means here, stated as §7.2a asks

It integrates **carriage, not concentration**: this program declares no `cell_network`, so there
is no level to integrate and a factor counts 1 for every minute the cell carried it. What is
carried is the reader's presence call (Ma et al. 2021, max adjusted expression ≥ 20% of the
factor's maximum), so `ELT-2.exposure(lineage) >= 15` means "the path to this cell spent at least
15 minutes above the atlas's presence threshold for ELT-2", and the 20% in that sentence is the
reader's, measured per factor over the atlas and not over this run. The rules' own threshold adds
no second normalisation, which is the whole change.

One consequence in another lane's file, recorded here as genomeos-d2 recorded the last one.
`plasticity.bio`'s terminal arm bounds each type by the number of cells of that type already
differentiated at 700 min, so those numbers move whenever the fate rules move: Neuron 233 → 216,
Hypodermis 94 → 81, GliaSheath 11 → 20, Intestine 26 unchanged. Nothing about the series, the
`competence` or the `commitment` changed — all five arms pass, and the two arms that assert the
wild type is untouched still read 130 muscle cells inside their 110..130 bound.

### A fate written on a growing integral is only a fate if the cell cannot take it back

The engine fixed a defect tonight and pinned it: a `cell_network` cadence used to let a rule that had
lost the precedence contest take a fate back, costing the worm 17 terminal fates, and
`test_a_network_cadence_does_not_change_which_fates_are_taken` asserts that the worm's Sulston score
is the same with a 6-minute cadence and without it. **The rewrite broke that invariant, and the
runtime was not at fault — the read was.**

`F.exposure(lineage)` grows for as long as the cell carries the factor, so a rule whose threshold was
not met when the fate was settled meets it later for no reason but the clock. `_pick_fate` refuses a
lower-precedence decision at a later decision point only when it *could already have applied* then
(`d.applies(c.fate_ctx)`), and a threshold on a growing integral is never in that set — so with enough
decision points every rule whose factor the cell ever carried fires, and the highest-priority one wins
whatever the threshold says. Measured, at a 6-minute cadence: **522 → 386 of 555**, 473 fates revised.

The obvious engine patch is **not** the fix, and that was measured rather than assumed: dropping
`d.applies(c.fate_ctx)`, so a strictly lower-precedence decision can never revise a settled fate,
recovers 31 of the 136 (386 → 417). Most of the revisions come from rules of *higher* priority whose
guard was genuinely false at birth, which the runtime is right to allow.

The fix is the construct §7.2a already has, and the worm's own program now declares it — the first
use of `commitment` outside the plasticity series:

```
commitment terminal_fate {
  establish: cell_type = Coelomocyte|Epithelium|...|Valve    # every terminal type
  locks: cell_type
  inherit: daughters
  release: never
  evidence: experimental "Fukushige & Krause 2005, Development 132:1795"
}
```

Its ablation, which is the only justification for it, is the same program with the block stripped
(`stability_under_a_re_decision_cadence` in the result file; these rows run the timers from the
program's declared seed rather than at their mean, as the engine's regression test does, which is why
the `mean` row reads 503 where the table above reads 533):

| program | cadence 0 | cadence 6 min |
|---|---|---|
| the lookup this replaces (`_integrated`, a guard that never changes) | 501 / 555 | 501 / 555 |
| `F.exposure(lineage) >= 15`, commitment removed | 522 / 555 | **386 / 555** |
| `F.exposure(lineage) >= 15`, **as shipped** | 522 / 555 | **522 / 555** |
| `F.mean(lineage) >= 0.25`, commitment removed | 503 / 555 | **327 / 555** |
| `F.mean(lineage) >= 0.25`, as shipped | 503 / 555 | **503 / 555** |

At cadence 0 the two are identical to the cell — 1,439 cells, 522 of 555, deaths 110/110, the same
in `embryo_contacts.bio` — so the block costs nothing in the program as it ships and is what makes
the fate a fate the moment anything asks the cell twice. `bio test` stays 202/202, and the five arms
of `plasticity.bio`, which declares a `commitment` of its own over a subset of the same types, pass
unchanged.

What this says about the language, and it is the general point rather than a worm detail: `exposure`
and `mean` are absolute and monotone, and a fate is not. A cell decides what to be from what it has
carried **so far**, and nothing in `when:` can say "as of now, once". So the cadence invariant is an
invariant **for time-invariant guards**, and any program whose fate guards are integrated reads needs
a `commitment` to have a stable fate at all. That belongs in §7.2a beside the two limits the engine
already states there.
