# The protein layer: a federated compiler, not a database

The human proteome is far better catalogued than whole-organism behaviour,
so GenomeOS does not build its own protein database. It compiles what the
public databases already know into one stable object per protein, keyed by
its **UniProt accession**, with every section carrying its source, its
evidence kind and a confidence. Simulations read the compiled object from
the local knowledge cache; they never call the databases.

```
                Ensembl
                   ↓
Genome ─────→ Gene / Transcript(s)
                   ↓
                UniProt  ← accession is the primary identifier
                   │
     ┌──────┬──────┼──────┬──────┐
     ↓      ↓      ↓      ↓      ↓
 InterPro Reactome STRING  HPA   PDB / AlphaFold
     └──────┴──────┼──────┴──────┘
                   ↓
       ProteinDefinition (data/knowledge/proteins/SYMBOL.json)
```

`genomeos protein TP53 --compile` runs it; the Molecules tab shows the
result; `genomeos proteome --chrom chr21` compiles a whole chromosome and
keeps only the coverage table as a result.

## What "mapped" means, per question

| question | source | evidence kind | confidence | section |
|---|---|---|---|---|
| which gene, which transcripts, how many protein products | Ensembl REST `lookup/symbol?expand=1` | curated | 0.95 | `genomic_origin` |
| amino-acid sequence, name, existence level | UniProtKB/Swiss-Prot (reviewed) | curated | 0.95 | `identity` |
| what it does, where in the cell | UniProt comments | curated | 0.95 | `function` |
| isoforms | UniProt alternative products + Ensembl cross-references (which transcript makes which isoform) | curated | 0.95 | `isoforms` |
| domains | InterPro via UniProt cross-references, plus UniProt features | curated | 0.9 | `domains` |
| modifications, processing | UniProt features (modified residue, glycosylation, disulfide, lipidation; signal, propeptide, chain) | curated | 0.95 | `modifications`, `processing` |
| experimental 3-D structures | PDB via UniProt cross-references, with **method** (X_RAY, CRYO_EM, NMR) and resolution | experimental | 0.95 | `structures_experimental` |
| predicted 3-D structure | AlphaFold DB (pLDDT per residue) | predicted | 0.7 | `structures_predicted` |
| pathways | Reactome via UniProt cross-references | curated | 0.85 | `pathways` |
| associations | STRING network API, combined score ≥ 0.7, channel scores kept; `physical_evidence` when the experimental channel ≥ 0.4 | predicted | 0.5 | `interactions` |
| where it is expressed | Human Protein Atlas entry JSON: tissue and single-cell specificity, nTPM per tissue, main subcellular location | experimental | 0.8 | `expression` |
| disease associations | UniProt disease comments with MIM ids | curated | 0.95 | `diseases` |

Predicted is never treated as observed: experimental and predicted
structures live in different sections, and STRING associations are
`predicted` with the experimental channel exposed so a physical interaction
can be told from a text-mining one.

## The model, and the distinction that matters

```
Gene
 ↓
Transcript(s)             Ensembl: TP53 has 40 transcripts, 37 protein products
 ↓
Protein isoform(s)        UniProt: P04637-1 … P04637-9
 ↓
modified protein states   30 annotated modification sites on p53
```

Never `Gene → Protein`. In BioIR the `Protein` entity now carries
`accession`, `isoforms`, `domains`, `structures`, `pathways` and
`interactions`; the compiled JSON is the full record behind it.

A **ProteinDefinition** is what the protein *is*. A **ProteinState**
(`genomeos.molecules.ProteinState`) is one protein in one place at one time:
tissue, cell type, level, localisation, modifications, isoform, with its own
evidence. `states_from_definition()` derives the first states from the HPA
section (protein × tissue × nTPM). Runtime rules act on states, not on
definitions.

## Measured, not assumed: coverage

`genomeos proteome --chrom chr21` compiles every protein-coding gene of a
chromosome and records, per question above, how many proteins have an
answer. The result is `data/results/proteome_chr21.json`; the compiled
definitions stay in the local cache (not committed, a few tens of KB each).
This is the measurable form of "how much of the proteome do we know".
Chromosome 21, 216 protein-coding genes, compiled in one pass
(3 read-through genes have no reviewed UniProt entry:
CFAP298-TCP10L, GET1-SH3BGR, IFNAR2-IL10RB):

| question | proteins | fraction |
|---|---|---|
| genomic origin | 216 | 100.0% |
| sequence | 213 | 98.6% |
| name | 213 | 98.6% |
| function | 194 | 89.8% |
| domains | 208 | 96.3% |
| pathways | 134 | 62.0% |
| interactions | 178 | 82.4% |
| expression | 216 | 100.0% |
| structure experimental | 98 | 45.4% |
| structure predicted | 212 | 98.1% |
| disease | 55 | 25.5% |

Sequence, name, domains and a predicted structure are close to complete;
function is known for nine in ten; pathways for six in ten; an experimental
structure for under half; a recorded disease association for a quarter.
That is the shape Albert's table predicted, now measured on real data. It
is the first milestone of the protein layer:

**Human genome → complete protein knowledge graph**: click any protein and
see everything science currently records about it, with every fact carrying
its source and confidence. From there the reaction and pathway parts of the
graph (Reactome also exports SBML, which the runtime already executes) turn
knowledge into execution.

## From knowledge to execution: running a pathway

Reactome exports every pathway as SBML (`ContentService/exporter/event/ID.sbml`),
with species annotated by the UniProt parts they contain and reactions with
inputs, outputs, catalysts and inhibitors, but **no rate laws**. So it is
not integrated like a BioModels file; it is run as a reachability graph
(`genomeos.molecules.reactome.PathwayModel`): a reaction fires when every
input and catalyst is present, its outputs become present, repeat to a
fixpoint. A knockout removes every entity containing the protein and
reports the reactions and products that are lost.

```
genomeos pathway R-HSA-69541 --knockout TP53
genomeos pathway --gene TP53 --knockout TP53      # every pathway TP53 belongs to, ranked by loss
```

Stabilization of p53 (Reactome v97): 28 entities, 15 reactions, 13 reachable
from the pathway's inputs. Without TP53, 6 of the 13 are lost (MDM2 binds
TP53, CHEK2 and ATM phosphorylation, ubiquitination and translocation) and
the phosphorylated and ubiquitinated p53 tetramers can no longer be made.
The Reactome content is curated; the reachability logic is inferred at 0.6,
because it ignores kinetics and treats inhibitors as reported, not applied.
The Molecules tab runs it: click a pathway chip with a gene loaded.

Across all 46 pathways TP53 belongs to (`data/results/knockout_TP53.json`),
144 reachable reactions are lost without it. The most dependent are the
regulation of TP53 activity through phosphorylation (18 of 21 reactions,
86%), association with co-factors (80%), methylation (67%) and acetylation
(64%); pathways where p53 is one input among many, such as ALK fusion
signalling, lose 2 to 3%. The ranking is what a reachability model can say;
how much each loss matters in a given cell needs the rates it does not have.

### Kinetics where a curated model exists (2026-09-11)

Reachability says which reactions survive a knockout; it cannot say how much
or when. BioModels holds hand-curated ODE models (BIOMD… ids, each one
reproducing its paper) for many of the same pathways, with rate laws,
parameters and initial conditions. `genomeos pathway R-HSA-5673001
--kinetic` takes the Reactome pathway's name ("RAF/MAP kinase cascade"),
searches BioModels with it (backing off to fewer words until something
curated matches), keeps the SBML under `data/knowledge/biomodels` and runs
it on the in-house engine (`genomeos/molecules/biomodels.py`); `--model
BIOMD…` picks a model directly, `--knockout SYMBOL` holds the matching
species at zero. A gene symbol lands on a model whose species are called
x1, x2, x3 through the model's own MIRIAM annotations: the species annotated
with the UniProt accession of the symbol is the one knocked out.

| model | knockout | what the run says |
|---|---|---|
| Hornberg2005 ERK cascade (BIOMD0000000084) | RAF1 (x1, x1p held at 0) | the MEK-P and ERK-P transients vanish (peaks −100%); steady states barely move (−4%), which is why the comparison counts peaks as well as final levels |
| Kholodenko2000 MAPK oscillator (BIOMD0000000010) | MEK (MKK, MKK-P, MKK-PP) | the ERK-PP oscillation is gone (peak 299 → 10, −97%); Mos-P steady state rises 62% because the negative feedback through ERK is cut |

The engine gained what curated models need: SBML functionDefinitions
(lambda calls), rateRules, species display names and annotations, and the
amount/concentration semantics (a species given as an amount in a
compartment of size ≠ 1 is a concentration in the maths, and a reaction's
substance-per-time rate changes it by rate/volume). Evidence: the model is
curated (BioModels), the run is derived (our RK4 integration), the knockout
effect is inferred from the run; time units are the model's own. Known
limit, recorded rather than hidden: Schoeberl2002's 100-species EGF-MAPK
model (BIOMD0000000019) does not reproduce its published transient on this
engine (the receptor binds, the cascade stays flat), while Sarma2012,
Huang1996, Kholodenko2000, Hornberg2005 and McClean2007 run; the
libRoadRunner adapter (roadmap item 12) is the reference for such cases.
Results land in `data/results/kinetic_<model>.json` without the time series.
The Molecules tab has the same thing under the pathway card: a kinetic model
found by the loaded pathway's name or given as a BIOMD id, a knockout, the
changed species as a table and the top species' time courses as a chart,
knockout curves dashed.

## Post-translational state (2026-09-12)

**Writers as rules.** The BioLang importer (`genomeos protein X --bio`,
`genomeos protein X --lib --bio`) now emits, after the `protein` block, one
rule per writer: `rule ATM modifies TP53 { strength: 1.0; evidence: curated
"UniProt: 5 modified residues written by ATM"; confidence: 0.8 }`. A program
that imports a protein therefore carries who acts on it, as rules the
runtime can apply, and the packaged proteome table stores the writers per
protein so the offline block has them too.


The model separated what a protein *is* from one protein in one place at one
time (ProteinState) and nothing populated a state's modifications. Now the
UniProt features already in every compiled definition are read as the
modifiable sites of the protein, with their chemistry and, where UniProt
names it, their writer ("Phosphoserine; by CDK5, PRPK, AMPK, NUAK1 and ATM"):
`genomeos ptm TP53` lists them, `genomeos ptm --writer PKA` lists a writer's
substrates, `genomeos ptm` gives the genome (`genomeos/molecules/ptm.py`,
`data/results/ptm_genome_wide.json`; the index is local under the protein
cache). The knowledge graph gains a `modifies` edge from writer to substrate,
curated from the same features.

| | |
|---|---|
| compiled proteins with at least one site | 13,083 of 19,452 |
| modifiable sites | 96,362: phospho 41,661, disulfide 18,329, glyco 17,603, acetyl 7,041, methyl 2,879, hydroxy 1,979, lipid 1,207, other 4,895 |
| named writers | 352, with 2,985 writer → substrate edges |
| writers with the most substrates | PKA 166, PKC 164, CDK1 95, CK2 94, FAM20C 94, PKB 68, PRMT5 59, AMPK 57, AKT1 56, ATM 53 |

Read with its limit: a site is one that *can* be modified, recorded from the
literature UniProt cites; whether it *is* modified in a given cell is a
measurement (phosphoproteomics) the project does not hold. Writers are
given as UniProt writes them, so "PKA" and "PKC" are families, not genes,
and land in the graph as nodes of their own. Ubiquitin cross-links are not
in the compiled `modifications` section yet, which is why no ubiquitin class
appears.

## Phosphosite observation (2026-09-27)

**What the field means:** `observed_in_cell_types_or_tissues` counts the
distinct cell lines or tissues in whose public mass-spectrometry data a
phosphosite was identified (Ochoa et al. 2020 reanalysis, 1% site-level FDR);
it does **not** say what fraction of the protein is phosphorylated there
(occupancy), whether the site is phosphorylated in a given cell, or that the
phosphorylation does anything. A curated site with no entry is absent from
that reference, which is not evidence that it is never phosphorylated.

Area C recorded measured modification state as blocked: the layer above says
only that a site *can* be modified. No open per-site, per-tissue occupancy
table exists, and none is built here. Two partial routes were probed.

**Route 1: the Ochoa et al. 2020 reference phosphoproteome.** Ochoa, Jarnuczak
et al., "The functional landscape of the human phosphoproteome", Nat.
Biotechnol. 38, 365-373 (2020), doi:10.1038/s41587-019-0344-3; author
manuscript PMC7100915. 112 PRIDE datasets from 104 cell types or tissues,
6,801 raw files reanalysed jointly with MaxQuant (PRIDE PXD012174); 119,809
sites pass 1% site-level FDR, 116,258 of them on reviewed UniProt proteins.

- Where. The article's Supplementary Tables 2 and 3 (Springer Nature
  `41587_2019_344_MOESM4_ESM.xlsx`, 54,465,047 bytes, and `MOESM5`,
  3,556,910 bytes; Content-Length checked 2026-09-27) and the authors' R
  package funscoR, https://github.com/evocellnet/funscoR, whose `data/`
  holds the same reference as R data files: `phosphoproteome.rda` (197,440
  bytes), `feature_spectral_counts.rda` (367,201), `feature_ms_pride.rda`
  (1,427,549).
- Terms. The article is under exclusive licence to Springer Nature; the
  PMC manuscript permits viewing and text and data mining for academic
  research under Nature's conditions, and the supplementary spreadsheets
  carry no licence of their own. The funscoR package declares
  `License: LGPL` in its DESCRIPTION (no version, no LICENSE file), which
  permits copying and redistribution of the package, data included. This
  project takes the funscoR files, not the spreadsheets, and keeps only
  per-site counts joined to its own curated sites, with the citation.
  funscoR's `psp.rda` is parsed from PhosphoSitePlus and is **not** used:
  PhosphoSitePlus terms are non-commercial.
- One row. `feature_spectral_counts`: `acc` (UniProt accession), `residue`,
  `position`, `Biological_samples`, `Spectral_Counts`; 116,258 rows,
  `Biological_samples` 1 to 83 (median 3), `Spectral_Counts` median 87.
  The paper's text names this feature "the number of different cell lines
  or tissues in which the site had been identified", which is what the
  field here is called. It is **not** a count of the 6,801 experiments: no
  per-site experiment count is published, so "observed in N experiments"
  cannot be delivered from this source.
- Size. 2.0 MB for the three files; well under the 2 GB stop.

**Route 2: CPTAC phosphoproteomics through cBioPortal.** The public API
(https://www.cbioportal.org/api, the client in `genomeos/cancer/cbioportal.py`)
lists 12 phosphoprotein profiles as `GENERIC_ASSAY` / `LIMIT-VALUE`
(profile ids `brca_cptac_2020_phosphoproteome`,
`luad_cptac_2020_phosphoproteome`, `lusc_cptac_2021_phosphoproteome`,
`ucec_cptac_2020_phosphoproteome`, `gbm_cptac_2021_phosphoproteome`,
`paad_cptac_2021_phosphoproteome`, `brain_cptac_2020_phosphoprotein`,
`coad_cptac_2019_phosphoprotein_quantification`, and the CPTAC
quantifications in `brca_tcga`, `ov_tcga` and their PanCan Atlas studies). Entities per profile, from `/api/generic-assay-meta`:
38,751 (breast 2020), 41,188 (lung adenocarcinoma), 18,806 (breast, TCGA
PanCan). Values are log2 abundance ratios to a pooled reference across
tumours, so the only honest name is `relative_abundance_in_tumours`.
Terms: the cBioPortal FAQ says data are under the ODC Open Database
License unless a study says otherwise (attribution, share-alike on a
derived database). The join is the obstacle: sites are named by gene symbol
and a RefSeq protein position (`NP_000010.1_1_1_69_69`, `A2M_S710s`,
`AAAS_pS462`), not by UniProt accession and residue, so each needs a RefSeq
to UniProt residue mapping first. Route 1 is smaller and joins directly, so
it is the one built; route 2 is recorded as reachable and usable, not built.

**Registration (written before the join).**

- Field: `observed_in_cell_types_or_tissues` (integer, 1 to 104) and
  `spectral_count` (peptide-spectrum matches), per curated site; source
  `Ochoa et al. 2020 / funscoR`. Never `occupancy`, `modified` or `active`.
- Join: UniProt accession of the compiled definition and the site's
  `start`, curated sites of class `phospho` only (the reference holds only
  phosphosites), and the reference residue must equal the residue at that
  position in the current UniProt sequence. A pair whose residue disagrees
  (sequence changed since the 2017 proteome) is counted and dropped, not
  joined.
- Expectation: 45% to 70% of the 41,661 curated phospho sites gain an
  observation, which is 19% to 30% of all 96,362 curated sites; the other
  classes gain none by construction.
- Checks that would show the join wrong: NPM1 S125 (P06748, the
  constitutive CK2 site) must be observed; no joined site may sit on a
  residue other than S, T or Y in the current sequence (TP53 M1 and every
  disulfide cysteine must not appear); residue disagreements must stay under
  2% of accession-position matches, or the numbering is off.

**Result, against the registration.** `scripts/phosphosite_observation.py`
(run with `uv run --with rdata`, the R data read in a temporary directory and
deleted) writes `data/results/phosphosite_observation.json` (0.55 MB, the
joined sites only); `ptm.observed_in(accession, position, observation)`
returns one site's record with the meaning attached.

| | |
|---|---|
| curated phospho sites | 41,661 of 96,362 |
| accession and position found in the reference | 29,299 |
| residue disagrees with the current sequence, dropped | 20 (0.07%) |
| joined: `observed_in_cell_types_or_tissues` set | 29,279 on 6,222 proteins |
| share of curated phospho sites | 70.3% (registered 45% to 70%) |
| share of all curated sites | 30.4% (registered 19% to 30%) |
| cell types or tissues per joined site | median 12; 2,238 in one; 16,699 in ten or more; max 83 |
| reference sites not joined to a curated site | 86,979 |

The coverage lands 0.3 points above the registered range, on the high side:
the expectation was too low and is reported as missed, not moved. All three
checks pass: NPM1 S125 is observed, TP53 M1 is not, no joined site sits on a
residue other than S, T or Y, and residue disagreements are 0.07% against a
2% bound. Curated sites of other classes that fall on a reference position
(acetyl 212, glyco 81, ADP-ribosyl 41, nitro 9, sulfo 3, other 2, methyl 1)
are counted and not joined: the same serine or threonine can carry either
chemistry, and the reference says only that it was seen phosphorylated.

Read with its limits. The reference and UniProt's curation are not
independent: many UniProt phosphosites were annotated from the same
large-scale studies Ochoa et al. reanalysed, so the 70% is partly
agreement with itself, not an external confirmation rate. The 86,979
reference sites without a curated counterpart are observations the curated
layer lacks; they are not added as sites here. Occupancy stays blocked.

## Relative abundance in tumours (2026-09-28)

**What the field means:** `relative_abundance_in_tumours` gives, per CPTAC
study on cBioPortal, the number of tumours with a value for a phosphosite and
the median of the per-tumour log2 ratio of that site's abundance to a pooled
reference of tumours from the same study. It is **not** occupancy (the
fraction of the protein phosphorylated there), not an absolute level
comparable across studies, and never per patient in the repository. A
curated site with no entry is absent from these profiles, which is not
evidence that it is never phosphorylated.

This builds route 2 of the section above. The obstacle was the site names:
gene symbol plus a position on a RefSeq protein, not UniProt accession and
residue. `scripts/cptac_phospho.py map` builds that mapping;
`ptm.parse_cptac_entity` and `ptm.transfer_position` are its two rules.

**Mapping, before any join.**

- Entities. Twelve phospho profiles, 430,832 entities. Four name their
  RefSeq protein (lung adenocarcinoma 2020, glioblastoma 2021, pancreatic
  2021, paediatric brain 2020); seven name the gene only (breast 2020,
  endometrial 2020, colon 2019, the two TCGA breast and two TCGA ovarian
  quantifications); the lung squamous 2021 profile holds 7,729 gene-level
  aggregates whose ids end in `acetylprotein` and is excluded whole.
- Single sites only. An entity counts only when it names one localised
  S/T/Y site. Dropped at parse time: several sites in one entity (lung
  5,919; glioblastoma 17,202; brain 464; breast 2020 5,911; TCGA breast
  13,100 and PanCan 1,880; TCGA ovarian 456 and PanCan 31), site not
  localised (lung 7,537; breast 2020 6,926), a second glioblastoma entity
  for the same site (the `.1` suffix, 7,620), not a RefSeq protein (11
  lung smORF or YP entities).
- RefSeq to UniProt. 12,025 RefSeq protein accessions, versioned as CPTAC
  used them. UniProt REST ID mapping (`RefSeq_Protein` to `UniProtKB`,
  release 2026_03, CC BY 4.0): the reviewed entries among the project's
  compiled definitions. A superseded version (`NP_x.1` when `.2` is
  current) has no cross-reference, so the unversioned accession is asked
  instead (3,893 proteins); the sequence compared is still the exact
  version, fetched from NCBI E-utilities. Result: 7,521 identical to the
  UniProt canonical sequence, 3,123 different, 1,074 with no UniProt
  entry, 305 with no reviewed entry among the compiled definitions, 2 with
  two compiled entries.
- Position transfer. Identical sequences carry the position (87,207 site
  entities); otherwise the 15-residue window around the site must occur
  exactly once in the UniProt sequence (27,029 carried). The residue letter
  must match at the RefSeq end (27 disagree, 0.02%) and at the UniProt end.
  Dropped: window not found (lung 494, glioblastoma 401, pancreatic 925,
  brain 63), window occurs more than once (12), RefSeq protein not mapped
  to one compiled entry (lung 2,682, glioblastoma 4,747, pancreatic 4,769,
  brain 426), several entities on one UniProt site in one profile (all
  dropped as ambiguous: 51, 0, 161, 96).
- Mapped, RefSeq-keyed: lung 24,491, glioblastoma 40,357, pancreatic
  45,586, brain 3,494; 66,427 distinct UniProt sites on 8,791 proteins
  (S 84%, T 14%, Y 2% of site entities).
- Gene-keyed profiles. The position is on an unnamed RefSeq isoform. Taken
  on the UniProt canonical of the one compiled entry for the gene symbol,
  the residue letter disagrees for 31,424 of 219,179 site entities
  (14.3%; 4.3% colon, 13.5% breast 2020, 14.6% endometrial, 16.2% to 17.6% TCGA
  breast, 18.4% to 19.2% TCGA ovarian). That rate says many positions are on
  another isoform, and a letter that matches by chance on a wrong isoform
  cannot be told apart, so these profiles are counted and **not
  committed**.

**Registration (written after the mapping, before any curated site was
joined or any value fetched; constants `ptm.CPTAC_REGISTRATION`).**

- Join: UniProt accession + position + residue letter must equal a curated
  site of class `phospho` (S, T or Y). Committed tier: RefSeq-keyed
  profiles only; a gene-keyed tier would need a residue mismatch under 2%,
  the same bound as the RefSeq tier, and at 14.3% it is counted only.
- Expected coverage of the 41,661 curated phospho sites: **35% to 55%**
  from the RefSeq-keyed profiles, 45% to 65% with the gene-keyed profiles
  counted. Why: the Ochoa reference, 116,258 sites from 104 cell types or
  tissues, reached 70.3%; these profiles give 66,427 mapped sites (109,474
  with gene-keyed ones) from four tumour types, single localised sites
  only, and tumour tissue lacks the cell-line studies much of UniProt's
  curation came from.
- Checks that would show the join wrong: EGFR Y1092 (P00533; legacy Y1068)
  present in lung adenocarcinoma; NPM1 S125 (P06748) present; TP53 M1 absent;
  no joined site off S, T or Y; RefSeq-end residue mismatches under 2%; the
  median of the per-site, per-study median log2 ratios within ±0.5 (pooled
  reference ratios centre near zero). AKT1 S473 is reported without an
  expectation: its tryptic peptide is poorly seen in global
  phosphoproteomes.
- Field: `relative_abundance_in_tumours`, per site a list of (study,
  tumours with a value, median log2 ratio). Never `occupancy`, `modified`
  or `active`; no per-tumour value is written to the repository (values are
  summarised in memory and never cached).
- Licences. cBioPortal data are under the ODC Open Database License 1.0
  unless a study says otherwise. Obligations and how the result meets
  them: attribution (the result's `licence` field and this section name
  CPTAC, cBioPortal and every profile used); share-alike (the per-site
  summary is a derived database and the result file declares it is offered
  under ODbL 1.0, separately from the code licences); keep open (it is a
  plain JSON file in the public repository). UniProt CC BY 4.0 and NCBI
  RefSeq are attributed the same way.

**Result, against the registration, negatives first.**
`scripts/cptac_phospho.py join` writes
`data/results/phosphosite_tumour_abundance.json` (3.6 MB, per-site
per-study summaries only, manifest complete); 133 MB fetched from
cBioPortal (gzip), no AlphaGenome request. `ptm.relative_abundance_in_tumours(accession,
position, result)` returns one site's list with the meaning attached.

- **The pancreatic 2021 profile is not in log2-ratio units.** The median of
  its per-site medians is 18.8 (lung 0.04, glioblastoma -0.02, brain -0.29):
  those are log2 intensities, whatever the profile's name says. It is
  excluded from the committed field. The registered pooled check (median of
  the RefSeq-keyed site medians within ±0.5) **passed at 0.11 with the
  pancreatic values inside it**, so that check was too weak to catch a
  whole profile in the wrong units; the per-profile medians found it. The
  gene-keyed TCGA breast PanCan profile sits at -0.52, just outside the
  bound; it is not committed either way.
- **Coverage missed high on both registered ranges.** RefSeq-keyed
  profiles, as registered (pancreatic included): 55.3% of the 41,661
  curated phospho sites counted per definition file (55.7% of 41,244
  distinct accession-position sites; 59 accessions sit in more than one
  definition file), against 35% to 55%. All profiles: 68.6% against 45% to
  65%. The expectation was too low, as it was for the Ochoa join, and is
  reported as missed, not moved. Committed after the units exclusion
  (lung adenocarcinoma, glioblastoma, paediatric brain): 19,104 sites on
  5,105 proteins, 45.9%.
- Not reached: 22,140 distinct curated phospho sites have no committed
  value. Mapped CPTAC sites that fall on a curated site of another class
  (glyco 205, acetyl 162, ADP-ribosyl 65, sulfo 25, nitro 7, other 4, lipid
  1, methyl 1) are counted and not joined.
- Checks: EGFR Y1092 (legacy Y1068) is present in lung adenocarcinoma
  (32 tumours, median log2 ratio -0.23) and glioblastoma; NPM1 S125 present
  (brain, 217 tumours); TP53 M1 absent; no joined site off S, T or Y;
  RefSeq-end residue mismatches 0.02% against 2%; committed median of
  site medians -0.02. AKT1 S473 is in no profile, as the registration
  allowed.

| committed profile | curated sites with a value | tumours per site (median) |
|---|---|---|
| lung adenocarcinoma 2020 | 12,330 | 102 |
| glioblastoma 2021 | 16,119 | 73 |
| paediatric brain 2020 | 2,594 | 217 |

Sites in all three: 1,745; in one only: 8,910. Read with its limits: a
median log2 ratio near zero says the site sits near the pooled reference of
that study, not that it is unphosphorylated; ratios from different studies
share no reference and are not compared. The gene-keyed profiles (breast,
endometrial, colon, TCGA breast and ovarian) stay uncommitted until their
RefSeq isoforms are known; `mapping_drops_by_profile` in the result keeps
every drop count per profile.

## Isoform-level expression (2026-09-12)

The model is Gene → Transcript(s) → Protein isoform(s), and expression had
been recorded per gene only. `genomeos rna X --isoforms` reads GTEx v8's
median TPM per transcript per tissue (RSEM isoform quantification, paged
from the portal API, cached under `data/knowledge/expression`) and labels
each transcript with the compiled definition's name, canonical flag and
protein length; per tissue it names the dominant isoform and its share, and
overall how often the canonical transcript is the one the tissue actually
makes (`gtex_isoforms`, `summarise_isoforms` in `genomeos/molecules/rna.py`).

Two worked cases. TP53: 28 transcripts, 54 tissues, the canonical TP53-201
(393 aa) dominates in all 54, the rest are near zero; the canonical choice
is the expression. APP: 17 transcripts, the canonical APP-201 (770 aa)
dominates in 33 tissues, APP-204 (751 aa) in 14, and APP-202, the 695-residue
neuronal isoform, in the cerebellum (79% of the gene's TPM) and cerebellar
hemisphere, which is the textbook: APP695 is the neuron's isoform, APP751
and APP770 the periphery's. The layer says which protein a tissue makes,
which is what `ProteinState.isoform` was for.

## The knowledge graph itself

The compiled definitions are one graph (`genomeos.molecules.graph`, no
network): proteins linked by STRING associations (score kept, physical
channel flagged), pathways they belong to (Reactome), InterPro domains and
HPA tissues, every edge carrying the evidence of the section it came from.
`genomeos graph` summarises it, `genomeos graph APP` prints one
neighbourhood, and the Molecules tab lays the neighbourhood out as a small
force graph. After chromosome 21 and the mitochondrial genome:

| measure | value |
|---|---|
| nodes | 4,175 (2,919 proteins of which 260 compiled, 458 pathways, 794 domains) |
| edges | 32,520 (30,404 associations, 15,139 with experimental support) |
| components among compiled proteins | 133, largest 49 |
| hubs | MRPL39 1271, PWP2 1265, MRPS6 1210, NDUFV3 1101, LTN1 1048 |

The hubs are mitochondrial ribosomal and respiratory proteins: STRING's
co-expression channel links every member of a large complex to every
other, so degree measures complex size, not importance. Read it with the
evidence kind in view. Tissue edges are few because HPA records per-tissue
values only for tissue-enhanced genes; GTEx (`genomeos rna`) covers the
rest. The mitochondrial proteome compiles at 100% on every question, the
13 best-studied proteins in the cell.

## The milestone: the whole human proteome compiled (2026-09-11)

The genome-wide job finished: 19,478 protein-coding genes on
25 chromosomes, each compiled from UniProt, InterPro, PDB, AlphaFold,
Reactome, STRING and the Human Protein Atlas, with the genomic origin from
the local GENCODE models. This is the measured answer to "how much of the
proteome do we know":

| question | proteins | fraction |
|---|---|---|
| genomic origin | 19,353 | 99.4% |
| sequence | 19,310 | 99.1% |
| name | 19,310 | 99.1% |
| function | 16,804 | 86.3% |
| domains | 19,280 | 99.0% |
| pathways | 11,333 | 58.2% |
| interactions | 15,882 | 81.5% |
| expression | 19,373 | 99.5% |
| structure experimental | 8,799 | 45.2% |
| structure predicted | 19,121 | 98.2% |
| disease | 4,960 | 25.5% |

Sequence, name and domains are essentially complete; a predicted structure
exists for almost every protein; function is known for six in seven;
pathways for six in ten; an experimental structure for under half; a
recorded disease association for a quarter. The mitochondrial proteins and
chromosome 19's zinc fingers sit at the two ends of that range.

The translation check ran on every chromosome: of 19,310 compiled
genes, 90.9% translate from the canonical transcript to exactly the
reviewed UniProt sequence and 97.8% match some annotated isoform
exactly. The 96 that disagree beyond boundary differences were read
against the gene models and fall into five groups, each detected by a
machine-checkable signature rather than by hand:

| why the reference does not give the curated protein | genes | examples |
|---|---|---|
| same length, different sequence: frame or isoform choice | 38 | AGAP9, ARL9, ASPRV1, C16orf82, CALML4, DEFB112 … |
| different length, partial similarity: exon boundary or isoform choice | 20 | CAPS, CCDC28A, CPEB2, CXXC4, DRC8, DUSP13B … |
| frameshift allele in hg38: the canonical CDS is not a whole number of codons (untagged) | 18 | CYP2D7, FAM246C, GPATCH4, IFNL4, KIR2DS4, OR10AC1 … |
| nonsense allele in hg38: translation stops before 70% of the CDS and of the curated protein | 13 | AKR7L, CASP12, DEFB109D, FCGR2C, MUC19, OR1P1 … |
| different protein (UniProt entry or gene model does not describe the same product) | 7 | C10orf95, C12orf76, EPM2A, FAM174C, GAGE12B, RTL8C … |

Two of those groups are facts about hg38, not about the compiler. The
reference genome is one haplotype, and at 31 loci it carries a
loss-of-function allele where UniProt describes the working protein: a
frameshift (the canonical CDS is not a whole number of codons and GENCODE
does not tag it incomplete; OR2B8 matches UniProt for 30 residues, then the
frame shifts and a stop follows at codon 32) or a nonsense allele (our
translation stops before 70% of the CDS while the curated protein is more
than 30% longer; CASP12, FCGR2C, IFNL4 and CYP2D7 are the textbook cases).
The frameshift signature fires on 29 canonical transcripts genome-wide; 7 of
them translate exactly anyway, six of which are the mitochondrial genes
whose CDS ends mid-codon by design and is completed by polyadenylation
(MT-CO3, MT-CYB, MT-ND1 to MT-ND4), which is the check that the signature
detects a property of the annotation and not a fault of the engine.
A negative result worth keeping: CDS length not divisible by three on its
own is useless as a signature, since it holds for 5.2% of coding
transcripts across 8,602 genes through 5'-incomplete models; the
untagged-and-canonical condition is what makes it specific.

Reading the disagreements also found a second compiler bug after the
synonym match: UniProt's query parser treats some symbols as stop words
(`gene_exact:WAS` returned every reviewed human protein), so the compiler
falls back to a plain `gene:` query when no entry names the symbol as
primary, and marks the identity section `symbol_match: false` when even
that fails (168 definitions, all renamed genes or identical paralogs
sharing one entry).

The knowledge graph built from the compiled definitions
(`graph_genome.json`): 41,982 nodes (19,797 proteins,
2,302 pathways, 19,846 domains, 37 tissues) and
327,024 edges, 95,098 of the associations with experimental
support; 15,165 of 19,283 compiled proteins sit in one connected
component. The compiled definitions occupy about 300 MB locally and are
not committed; the summaries are.

### The disagreements, read with the isoform the body makes (2026-09-12)

The 96 genes whose canonical transcript does not translate to the curated
protein were triaged by mechanism; the isoform layer adds the question the
triage could not ask: is the canonical transcript even the one the body
makes? `scripts/disagreements_isoforms.py` takes GTEx's dominant transcript
for each gene, translates it locally and compares
(`data/results/translation_disagreements_isoforms.json`, 151 s):

| verdict | genes |
|---|---|
| the body makes the canonical transcript and it still differs: the reference allele is the disagreement | 38 |
| the body makes another isoform and that one differs too | 18 |
| the dominant transcript is non-coding in the local models | 16 |
| no expressed isoform in GTEx, or the gene is not in GTEx (olfactory receptors mostly) | 24 |
| the isoform the body makes matches UniProt (a canonical-choice artefact) | **0** |

Not one disagreement is the annotation picking the wrong isoform: where
GTEx sees expression, the transcript it sees is the disagreeing one, and
where it sees another, that one disagrees as well. The frameshift and
nonsense alleles of hg38 stay what the triage said (the reference carries
an allele the curated protein does not), and the "same length, different
sequence" group (18 of its 38 measured cases make the canonical transcript)
is now most likely the same thing at single-residue scale: reference
alleles, not isoforms. What the isoform check cannot settle is the 40 genes
GTEx does not resolve, and it does not need to: their triage already names
the mechanism.

### The disagreements, held against three real genotypes (2026-09-12)

If the reference carries a minor allele where UniProt describes the common
protein, a real person should carry the curated version.
`scripts/disagreements_genotype.py` takes each disagreement gene, applies the
SNVs a local individual carries inside the canonical CDS (alt allele, every
site), translates and compares again, for HG002 and both parents
(`data/results/translation_disagreements_genotype.json`):

| HG002 | genes |
|---|---|
| no variant of this person inside the CDS | 40 |
| the person's SNVs leave the disagreement as it is (or make it worse) | 52 |
| only indels inside the CDS, not applied (a frame change is possible) | 1 |
| the person's alleles give the curated protein | **0** |

HG003 and HG004 give the same picture (46 / 46 / 1 and 38 / 54 / 1, none
restored). So the 96 are not common alleles this Ashkenazi trio happens
to carry the other way: where the trio does vary inside these CDS, the
variation is unrelated to the curated protein and often makes the
translation shorter still (olfactory receptors that are pseudogenes in
most people, OR12D1 among them). Together with the isoform check above,
the disagreements are now bounded from three sides: not isoform choice,
not this trio's alleles, and by mechanism either a reference allele the
population does not share widely, or a gene model and a UniProt entry that
describe different products. What would settle each one is population
allele frequency at the exact site (gnomAD), which the lookup layer already
fetches one variant at a time.

### The residue-level disagreements settled by allele frequency (2026-09-12)

For the 38 genes whose canonical translation has the curated protein's
length but not its sequence, each differing residue names a codon and the
single-base change that would give UniProt's residue is a variant with a
genomic position. `scripts/disagreements_frequency.py` asks gnomAD and
dbSNP (through Ensembl VEP, cached by the lookup layer) about each such
allele, up to six per gene
(`data/results/translation_disagreements_frequency.json`, 112 s):

| verdict | genes |
|---|---|
| the curated allele is the common one: the reference carries a minor allele (OR9H1 rs7555046 at 0.99; FCGBP, four sites at 1.0) | 2 |
| the curated allele is a known polymorphism, not the majority (HLA-DQA1, MICA, GSTT2, SERPINA2, PRB4, SAMD1) | 6 |
| the curated allele is in dbSNP but rare or unmeasured | 10 |
| the curated allele is unknown to dbSNP and gnomAD: the gene model and the UniProt entry describe different things | 16 |
| no residue-level difference found (alignment rather than substitution) | 4 |

That is the last side of the bound. Eight of the 38 are population
variation, of which two are the reference carrying the minor allele (the
classic case, FCGBP's four sites all at frequency 1.0) and six are
polymorphic loci where UniProt chose one haplotype (the HLA and MICA
entries, as expected of the MHC). The 16 with alleles no population
resource has seen are not alleles at all: the UniProt sequence there comes
from a different gene model, a different assembly or a corrected entry.
Across all four checks, the 96 disagreements now read as: 31 hg38
frameshift and nonsense alleles (by mechanism), 8 residue-level population
variants, 7 different products, and the rest gene-model differences that
no isoform, no genotype of the trio and no population allele explains.

## The proteome as a library

The compiled proteome is now a BioLib layer that ships with GenomeOS:
`genomeos/lib/data/proteome.json.gz`, 19,283 human proteins in 2.3 MB,
distilled from the local cache by `genomeos.lib.proteome.distil()`. Per
protein it keeps the accession, length, existence level, a one-line
function, location, up to ten InterPro domains, twelve Reactome pathways,
ten physically supported partners, the tissue and cell-type pattern, up to
five diseases, the count of experimental structures and whether AlphaFold
has a model, plus the evidence per field and a `symbol_match` flag. It
answers offline:

```
genomeos protein BRCA1 --lib          # the protein from the packaged table, no network
genomeos protein BRCA1 --lib --bio    # the same as a BioLang protein block
genomeos libs --proteome              # counts and evidence for the whole library
```

Storage rule respected: the full definitions (290 MB, rebuildable with
the proteome job) stay local; the library is what travels. Regenerate it
after a new compile with `genomeos.lib.proteome.distil()`.

## Verified: the engine reads genes the way the curators do

`genomeos verify --chrom C` translates the canonical transcript of every
compiled gene from the reference sequence and compares it with the reviewed
UniProt sequence, then tries every coding isoform:

| chromosome | genes | canonical identical | some isoform identical | same protein, other boundaries | real disagreements |
|---|---|---|---|---|---|
| chr21 | 213 | 91.5% | 99.1% | 1 | 1 (SH3BGR) |
| chr22 | 426 | 91.1% | 96.5% | 7 | 4 (CYP2D7, FAM246C, KLHDC7B, GSTT2) |
| chrM | 13 | 100% | 100% | 0 | 0 |

A canonical-choice difference (16 on chr21) is GENCODE and UniProt naming
different isoforms as canonical; the protein is still produced exactly by
another transcript. The two real disagreements are recorded with their
identities in `translation_vs_uniprot_chr21.json` for a look.

The pass found two things the engine did not know. Selenoproteins (25 human
genes, GENCODE tag `seleno`) read an in-frame UGA as selenocysteine, so
SELENOM stopped at residue 47 instead of 145; the engine now has a
selenocysteine table chosen by the tag, and the three selenoproteins on
chr22 match exactly. And a gene can own two reviewed UniProt entries (MIEF1:
the 463-residue protein and a 70-residue microprotein from an upstream ORF);
the compiler now takes the longest entry with the matching primary name.

This pass also caught a compiler bug: UniProt's `gene_exact` query matches
synonyms, so the first hit for MIF was a 560-residue protein that carries
MIF as an alias. The compiler now takes the entry whose primary gene name is
the symbol; 16 cached definitions were recompiled.

## BioLang importer

`genomeos protein TP53 --bio` writes the compiled definition as a BioLang
`protein` block (accession, sequence, isoforms, domains, pathways,
interactions with physical evidence, structures), so a program can carry
real protein knowledge with its evidence and confidence. The parser accepts
these properties on any `protein` block.

## Storage

Stream, distil, discard applies: one compiled definition is small enough
to keep locally, the coverage summary is the committed result, and a
`--refresh` recompiles from the sources. A transient failure of one source
is recorded as `evidence: none` with the error and is refetched on the
next load, so a bad hour at Ensembl never becomes a permanent gap.
