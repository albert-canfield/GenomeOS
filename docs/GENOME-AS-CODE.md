# The genome as code: a reverse-engineering guide

GenomeOS's purpose is to *unbuild* the genome: understand it, organise it,
and produce a clear base for genetic building and construction. This document
is the map. It answers the four questions the project asks of DNA, and it
names the "libraries" that already exist in the code. The machine-readable
version is `genomeos/lib/catalog.py` (`genomeos libs -v`).

## 0. What kind of code is it

Before reverse engineering anything, be precise about what the artefact is.

| Software concept | Genome reality | Consequence for GenomeOS |
|---|---|---|
| Source file | 3.1 billion bases, two copies (maternal, paternal), plus 16.6 kb of mtDNA | keep it lossless; diploid from day one |
| Functions | ~19,400 protein-coding genes, ~20,000 non-coding RNA genes | `Gene` entities from GENCODE |
| Code vs comments/history | 1.5% coding; ~45% transposable elements and repeats; ~8% under constraint | most bases are not instructions; mark `UNKNOWN`, do not invent |
| Config / feature flags | ~1 million cis-regulatory elements (promoters, enhancers, silencers, insulators) | `Region` entities; context-gated `Rule`s |
| Runtime that reads the code | the cell it is already in (ribosomes, polymerases, mitochondria) | genome + bootstrap state, never genome alone |
| Program counter | none. Execution is massively parallel and concentration-driven | continuous and stochastic engines, not an instruction pointer |
| Namespaces / modules | conserved gene families and pathways reused across tissues | BioLib libraries with `when:` contexts |
| Version control | evolution; the pangenome records the population's branches | pangenome-graph design target |
| Runtime code generation | V(D)J recombination in lymphocytes rewrites DNA to make antibodies | events that modify `G` locally |
| Self-modifying config | methylation and chromatin marks change with age and environment | epigenome lives in `S`, read by clocks |

The most useful single idea: **the genome is a library of parts and rules,
not a script.** The body plan is not written as a sequence of steps; it
emerges from rules that fire when their conditions are met. That is why
BioLang is declarative.

## 1. Time and ageing: the timers

The first consideration. There is no single "age gene" or master countdown.
Ageing is the sum of several timers plus accumulated damage, and construction
(development) uses different timers from maintenance (adulthood).

### Construction timers

| Timer | Mechanism | Library | Period / scale |
|---|---|---|---|
| Segmentation clock | HES7/Notch oscillation lays down one somite per cycle | `timer.segmentation_clock` | ~5 h in human, ~2 h in mouse |
| Heterochronic switch | LIN28 falls, let-7 rises, gating the transition from early to late programmes | `timer.developmental_timing` | days to years |
| Puberty | KISS1/GnRH axis released when MKRN3 brake lifts; ESR1 then closes growth plates | `timer.puberty_reproduction` | years |
| Growth axis | GH/IGF-1/mTOR set the rate of building | `timer.growth_longevity_axis` | continuous |

Construction timing is mostly **relative**, not absolute. Stages are triggered
by state (concentration thresholds, cell counts, completed structures), which
is why development is robust to temperature and nutrition within limits.

### Maintenance timers and counters

| Timer | Mechanism | Library | What it measures |
|---|---|---|---|
| Telomere | each division loses 50–100 bp unless TERT is active | `timer.telomere` | divisions |
| Epigenetic clock | methylation drift at thousands of CpG sites, read by DNMT/TET machinery | `core.chromatin_epigenome` | time, modified by pace |
| Somatic mutation burden | ~20–50 mutations per cell per year, tissue dependent | `core.dna_repair` | time × exposure |
| Senescence register | p16INK4a (CDKN2A) rises; SASP secretome spreads the signal | `timer.senescence_checkpoint` | damage |
| Circadian | 24 h cycle in every cell; degrades with age | `timer.circadian` | phase, not age |
| Depletion | ovarian reserve, stem-cell pools, naive T-cell pool | `timer.puberty_reproduction`, `timer.stem_cell_niches` | remaining capacity |

The maintenance schedule is best located through failure: progeroid genes
(`timer.progeroid_maintenance`: LMNA, WRN, ERCC6/8, TERT) are the modules
whose loss accelerates ageing. Each one points to a maintenance library.

**What v0.1 already runs:** `genomeos age` simulates telomere, mutation,
epigenetic and senescence clocks per cell with cited parameters.

## 2. Structure: the blueprint

The second consideration. The blueprint is encoded as **positional code +
reusable signalling + master switches**, not as a drawing.

1. **Coordinates.** Axis genes establish a coordinate system in the embryo.
   The HOX clusters are the clearest case: 39 genes in four clusters,
   expressed in the same order along the chromosome as along the body.
   Left–right comes from NODAL/LEFTY/PITX2; dorsal–ventral from BMP versus
   its antagonists. (`blueprint.axes`)
2. **Segmentation.** The trunk is built as repeated units by an oscillator
   coupled to a moving gradient. (`blueprint.segmentation`)
3. **Germ layers.** Gastrulation partitions cells into ectoderm, mesoderm,
   endoderm, each with its own master factors. (`blueprint.germ_layers`)
4. **The shared API.** About a dozen signalling pathways (WNT, Hedgehog,
   BMP/TGF-β, FGF, Notch, RTK/RAS, PI3K/mTOR, JAK/STAT, Hippo, retinoic
   acid, Ephrin, Semaphorin) are reused by every organ. The same pathway
   yields different outcomes depending on context. This is exactly the
   `when:` clause in BioLang. (`blueprint.signalling_toolkit`)
5. **Master switches per organ.** PAX6 for eye, NKX2-5/GATA4/TBX5 for
   heart, PDX1 for pancreas, MYOD1 for muscle, TP63 for skin, SRY→SOX9 versus
   FOXL2 for gonads. Switching one on in the wrong place can build the organ
   there (ectopic eyes from PAX6 are the classic demonstration).
   (`blueprint.organ_*`)

So the blueprint is compressed: coordinates × shared signals × context ×
master switches. GenomeOS must represent those four things explicitly; a
gene list alone cannot reproduce structure.

## 3. The parts list

The third consideration: what needs to exist. The genome does not contain a
parts list as such; the list is what the blueprint produces. But curated
ontologies already enumerate it, and GenomeOS imports them rather than
retyping them.

| Level | Source | Size | GenomeOS entity (planned v0.2) |
|---|---|---|---|
| Organ systems | Uberon | 11 | `Tissue` / `Organ` |
| Organs | Uberon | ~78 | `Organ` |
| Anatomical terms | Uberon | ~15,000 | `Region`-like anatomy nodes |
| Cell types | Cell Ontology | ~2,900 | `CellType` |
| Cell states | Human Cell Atlas / CELLxGENE | tens of millions of profiled cells | `CellState` templates |
| Molecular parts | GENCODE, UniProt | ~19,400 proteins, ~20,000 ncRNAs | `Gene`, `Transcript`, `Protein` |
| Interactions | Reactome, CellPhoneDB | ~2,900 pathways, ~2,000 ligand–receptor pairs | `Rule` |

The inventory question then becomes computable: **for each part, which
libraries does it need?** A cardiomyocyte needs `core.*`, `core.energy` at
high rate, `blueprint.organ_heart`, `systems.nervous` (ion channels) and
`systems.circulation_respiration`. Membership tables are the Phase 1
deliverable.

## 4. How systems connect

The fourth consideration. Integration in the genome is written as
**communication protocols**, and they are all in the code:

| Channel | Speed | Reach | Library |
|---|---|---|---|
| Ligand–receptor (paracrine) | minutes | neighbours | `systems.ligand_receptor_protocol` |
| Hormones (endocrine) | minutes to hours | whole body | `systems.endocrine` |
| Synapses | milliseconds | addressed | `systems.nervous` |
| Immune cytokines and MHC presentation | hours | whole body, addressed by antigen | `systems.immune` |
| Extracellular matrix and adhesion | structural | mechanical, local | `systems.extracellular_matrix_adhesion` |
| Gap junctions | milliseconds | direct neighbours | `systems.extracellular_matrix_adhesion` (GJA1) |
| Metabolic exchange | hours | via blood | `systems.metabolic_homeostasis` |

Every channel is a ligand encoded by one gene reaching a receptor encoded by
another. That means the dependency graph between tissues can be *computed*
from expression data plus the ligand–receptor table, and represented as
`Rule`s with `when: tissue = X`. Building that graph is Phase 2.

## 5. The reverse-engineering method

Software reverse engineering has a standard sequence. It applies here.

| Step | Software | Genome | GenomeOS deliverable |
|---|---|---|---|
| 1. Inventory | list symbols and files | GENCODE genes, cCREs, repeats | `Gene`/`Region` entities from GFF3 |
| 2. Group by module | find libraries and namespaces | GO, Reactome, gene families | library membership tables |
| 3. Find entry points | main(), init | pluripotency factors, gastrulation | `blueprint.pluripotency` bootstrap |
| 4. Trace calls | call graph | ligand–receptor, TF–target networks | `Rule` graph |
| 5. Find the timers | schedulers, counters | telomere, clocks, LIN28/let-7 | `timer.*` engines |
| 6. Find what is essential | delete and see what breaks | DepMap essential genes (~2,000), OMIM, knockout mice | confidence per rule |
| 7. Recompile and run | rebuild from source | simulate cell → tissue → development | BioVM |
| 8. Modify and compare | patch and test | perturbation in silico, twin versus twin | BioForge |

Steps 1 to 2 are mechanical and can start now with public data. Steps 4 to 6
are where AI models (AlphaGenome, STATE) supply `predicted` evidence for what
experiments have not yet measured. Steps 7 to 8 are the platform's reason to
exist.

## 6. What a "clear base for building" looks like

The end state is a library tree a bio-engineer can read the way a developer
reads a standard library:

```
bio/
  core/        replication transcription splicing translation folding
               degradation dna_repair cell_cycle apoptosis energy
               cytoskeleton membrane_transport chromatin_epigenome
  blueprint/   pluripotency germ_layers axes segmentation signalling_toolkit
               organ_heart organ_nervous organ_eye organ_endoderm_gut
               organ_kidney_urogenital organ_musculoskeletal
               organ_blood_immune organ_skin
  timer/       telomere circadian segmentation_clock developmental_timing
               puberty_reproduction growth_longevity_axis
               senescence_checkpoint progeroid_maintenance stem_cell_niches
  systems/     endocrine nervous immune circulation_respiration
               extracellular_matrix_adhesion ligand_receptor_protocol
               metabolic_homeostasis
  parts/       cell_types tissues_organs genome_census
```

Each library will hold: its genes (from GENCODE), its rules (from Reactome,
GO and literature), its parameters (with evidence), its contexts (which
cell types and stages), and its tests (known phenotypes the module must
reproduce). At that point "genetic building" means composing libraries under
constraints and running the result, with the platform reporting where the
prediction is grounded and where it is a guess.
