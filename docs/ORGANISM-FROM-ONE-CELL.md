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
founder transformations of Du et al. 2014 reproduced; the misses are the PAR polarity genes, whose rules are not written.

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
