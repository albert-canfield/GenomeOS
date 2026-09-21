# Therapeutic targeting and mechanism reasoning

GenomeOS reads a tumour's alterations and asks the question a target biologist
asks:

> What is molecularly different about these cells, can that difference be
> reached selectively from outside, and what therapeutic mechanism could
> exploit it without damaging normal tissue?

It produces research hypotheses, rankings, evidence and mechanistic reasoning.
It is not a drug designer, and the outputs are not clinical recommendations.

```
MUTATION
   -> MOLECULAR CONSEQUENCE
      -> CELLULAR CONSEQUENCE
         -> EXTERNALLY ACCESSIBLE DIFFERENCE
            -> NORMAL-vs-TUMOUR DIFFERENTIAL
               -> THERAPEUTIC VULNERABILITY
                  -> COMPATIBLE MECHANISM
                     -> EVIDENCE
                        -> UNCERTAINTY
```

## The three target spaces, searched in parallel

A nuclear mutation does not get `THERAPEUTIC TARGET = NO`. It gets three
separate answers:

```
                    CANCER ALTERATION
                          |
          +---------------+---------------+
          |               |               |
   DIRECT SURFACE     NEOANTIGEN      INDIRECT EFFECT
      TARGET            / HLA            TARGET
          |               |               |
   surface binder    TCR / TCR-mimic   altered pathway gives
                                       an exploitable surface
                                       phenotype
```

### The target the tumour never altered

CD19 and BCMA are not mutated, amplified or rearranged. They are lineage
antigens the malignant cell carries in quantity and most healthy tissue does
not, and a pipeline that reaches a gene only through an alteration cannot
propose either — a hole in the middle of the target space rather than a missing
refinement. `--scan N` closes it, using the patient's own tumour RNA against
the packaged healthy-tissue atlas over the membrane and CD-marker universe:

```bash
genomeos therapeutic --tumour t.vcf --rna tumour_tpm.tsv --scan 6
```

On the demo B-cell lymphoma (`data/demo/expression/`, synthetic values in a
plausible range) the scan proposes FCRL5 at 5.9x and CD19 at 4.3x, both direct
surface targets ranking above the TP53 driver, and neither reachable by any
other route in this pipeline.

**The rejections are the more instructive half.** MS4A1 (CD20) is turned down
at 2.2x, CD79A at 1.4x, CD22 at 1.0x — the targets of the most successful
antibodies in oncology, every one failing a selectivity screen because healthy
spleen is full of the same lineage: CD20 sits at 247 nTPM there. They are
reported as turned down with the number rather than dropped, because a
threshold here is a sort order and never a verdict. The screen measures
selectivity against healthy tissue, and an antigen shared with a healthy
lineage can still be the right target when losing that lineage is survivable.

Three refusals hold. The scan runs only on patient RNA — a cohort describes a
cancer type and cannot say what this tumour displays, so without patient RNA
the answer is that no expression-driven target can be proposed. Raised
transcript is never called protein, still less protein on the surface, and the
missing evidence is recorded as such. And every candidate it produces carries
the sentence the CD19 story actually teaches: hitting an unaltered antigen hits
the healthy lineage that shares it, and B-cell aplasia is not a side effect of
CD19 therapy but the same event seen from the other side.

## Five concepts, kept apart

Collapsing these is what makes a pipeline assume every target needs a toxin.

| Concept | Question it answers | Example |
|---|---|---|
| **Target** | what distinguishes the cancer cell | HER2 |
| **Binder** | what recognises the target | antibody-like |
| **Mechanism** | how therapeutic action follows | ADCC |
| **Cargo** | optional material carried | none |
| **Effector** | what ultimately acts on the cell | NK cells |

`payload = none` is a complete, active mechanism. The model says so: an ADCC
fit carries `payload_required: false`, `cargo: none`, `effector: endogenous
immune system`, `immune_components: [NK cells, macrophages]`.

## The pipeline

| Stage | Module | What it establishes |
|---|---|---|
| consequence | `genomeos.cancer.tumour` | gene, consequence, protein change, driver status |
| localisation | `therapeutics/localisation.py` | compartment, topology, which residues face outwards |
| accessibility | `therapeutics/pipeline.py` | the target class, from localisation and never from the mutation |
| expression and selectivity | `therapeutics/expression.py` | healthy-tissue levels, tumour levels when supplied, cohort behaviour, the differential |
| trafficking | `therapeutics/trafficking.py` | internalisation, endosomal and lysosomal routing, recycling, shedding |
| neoantigen / HLA | `therapeutics/neoantigen.py` | whether the variant yields a novel peptide, and the layers of presentation |
| pathway-induced | `therapeutics/pipeline.py` | surface proteins associated with a disrupted driver |
| normal-tissue exclusion | `therapeutics/safety.py` | tissues a binder must avoid, and how firmly that is known |
| structure and epitope | `therapeutics/structure.py` | resolved regions, and mutation-created surface epitopes |
| mechanism selection | `therapeutics/mechanisms.py` | 15 mechanisms with declarative requirements |
| scoring | `therapeutics/scoring.py` | component scores and one inspectable overall |
| design dataset | `therapeutics/design.py` | the specification a binder-design system consumes |

Every stage takes a provider, never an API. `therapeutics/providers.py` defines
the protocols (`ProteinAnnotationProvider`, `ExpressionProvider`,
`TumourExpressionProvider`, `CohortExpressionProvider`, `StructureProvider`,
`TraffickingProvider`, `TranscriptSequenceProvider`, `NeoantigenProvider`,
`CancerEvidenceProvider`, `HomologyProvider`) so a source can be replaced
without touching the reasoning, and so the whole pipeline runs offline against
the caches.

## Three questions about expression, kept apart

Expression is where a target pipeline is most tempted to cheat, so the three
questions are answered separately and never substituted for one another.

| Question | Source | What it can and cannot say |
|---|---|---|
| Does healthy tissue carry it? | Human Protein Atlas consensus nTPM | population-level, measured, the basis of the safety score |
| Does this cancer type raise it above normal? | a cBioPortal cohort, scored against the study's own normal samples | a property of the cancer type, never of this patient |
| Does *this* tumour carry it? | the patient's own RNA-seq | the only answer that is about this patient |

The cohort is the one worth explaining. `--cohort brca_tcga_pan_can_atlas_2018`
fetches, per gene, the distribution of tumour-versus-normal z-scores across the
study, and with it the fraction of tumours in which the gene is raised above
normal tissue at all. That fraction is a far better selectivity prior than a
tissue-specificity class, because it is measured against normal samples instead
of inferred from how widely a gene is expressed. In breast cancer it separates
MKI67 (raised in 66% of 1,082 tumours) from ERBB2 (11%) from RUNX1 (0.8%).

It is still a prior. A cohort describes a cancer type, so its contribution is
capped at 0.6, below what a patient measurement can reach, and the report says
so in the score's own basis line.

Patient RNA, when supplied, overrides everything. Two comparisons are then
possible and only one of them is safe:

- the gene's rank within the patient's own transcriptome, which needs no unit
  and is therefore valid;
- the raw value against healthy-tissue nTPM, which is an approximation because
  TPM and nTPM are different normalisations of the same idea.

A patient TPM is never compared against a cohort's RSEM or z-score. Those are
different scales, and a percentile computed across them would be meaningless.

Copy number, supplied with `--cnv`, is treated as a third thing again: it bounds
how much protein a cell could display and never shows that it does, so it scores
at most 0.5 and the basis line says why.

### Healthy tissue, shipped rather than fetched

Every candidate needs the same question answered — how much of this protein a
healthy person already carries, and where — and it used to cost one Human
Protein Atlas request per gene. Offline that returned "no cached HPA record",
normal-tissue safety was unknown, and every candidate was capped for missing
safety data; and scanning the membrane universe for a target was impossible,
because five thousand requests is not a scan.

The Atlas answers a whole protein class in one request, so
`scripts/normal_tissue.py` fetches the two that matter — its 5,573 predicted
membrane proteins and its 384 CD markers, which is where the antigens of the
approved cell therapies live — keeps the twenty tissues GenomeOS weighs, and
ships the result gzipped inside the package
(`genomeos/therapeutics/data/normal_tissue.json.gz`, 5,588 genes, 396 KB, two
requests, seven seconds). The raw download is discarded and the record of the
run is committed as `normal_tissue_atlas`. Rows keep the Atlas's own column
names, so a local row and a freshly fetched one go through `normal_profile`
unchanged.

**Two of the twenty tissues had never been read.** The Atlas returns the
`t_RNA_skin_1` and `t_RNA_stomach_1` fields under the titles "skin 1" and
"stomach 1"; the lookup asked for "skin" and "stomach" and silently got
nothing, so both were recorded as an explicit absence on every gene ever
scored. Skin is where an EGFR antibody's classic toxicity shows: EGFR carries
48.4 nTPM there and it now appears among EGFR's tissues at risk, where it
belongs. Column titles are derived from the field ids now rather than written
down, and a test asserts all twenty are read. The benchmark does not move,
because both surface cases already sit on the poor-safety cap — EGFR's raw
score is 0.572 and the cap holds it at 0.5 — which is why a benchmark is not
the only control a change needs.

## Accessibility is not a membrane word

The single most important correction in this module: UniProt annotates SRC as
`Cell membrane`, and SRC sits on the *inner* leaflet on a myristoyl anchor.
A binder in the blood will never touch it. So GenomeOS claims external
accessibility only with positive evidence of an outward-facing part:

- a curated extracellular topological domain, or
- a transmembrane segment with a plasma-membrane location, or
- a GPI anchor, or
- a curated cell-surface annotation, or
- a signal peptide.

Without one of those, a membrane annotation is held below the reachability
threshold and the orientation is reported as *not established*. A lipid anchor
with no transmembrane segment is reported as *cytoplasmic face*.

The same distinction runs one level deeper: a surface protein whose *mutation*
sits in the cytoplasmic tail has no mutation-specific extracellular epitope,
however accessible the protein is. The demo tumour's APP N770K is exactly that
case, and the report says so.

## Reconstructing the altered protein

A frameshift's novel peptide stretch cannot be read off a VCF record. The indel
has to be applied to the transcript the consequence was called on, and the
result translated. So the transcript id is carried down from Ensembl VEP, the
coding sequence is fetched, the HGVS coding change is applied, and both proteins
are translated and compared.

That turns a rule about the variant class into a measurement of what the protein
becomes:

| Variant | What the sequence shows |
|---|---|
| frameshift | the novel C-terminal stretch, which exists in no healthy protein |
| in-frame indel | no new residue, but a junction that exists in no healthy cell |
| substitution | one changed residue, with its wild-type counterpart kept beside it |
| premature stop | no residue differs from the reference; the protein is simply shorter |

Two details matter. Alignment is decided by the reading frame, not by protein
length: a frameshift with no downstream stop can produce a protein of exactly the
same length that shares no sequence with the reference. And a frameshift scores
above a substitution on neoantigen strength, because its residues have no
wild-type counterpart at all, which is the easiest discrimination problem in the
pipeline rather than the hardest.

The reconstruction also settles isoform disagreements. The demo's RUNX1 Y480* is
called on a 481-residue transcript while the canonical UniProt entry is 453
residues, so the canonical route could only report a mismatch. Rebuilding from
ENST00000675419 confirms the truncation from sequence instead.

Changes outside the coding sequence, and forms this module does not reconstruct,
are declined and recorded as missing rather than approximated.

## Evidence, and the refusal to fill gaps

```python
class Evidence:
    source: str
    source_type: "database" | "publication" | "clinical_trial" | "patient_data" | "prediction" | "derived"
    claim: str
    level: "clinical" | "human" | "preclinical" | "in_vitro" | "computational" | "inferred"
    confidence: float
```

Clinical observation and computational prediction never merge: the report
groups evidence by level and the scores carry the level that produced them.

Where a fact cannot be established the pipeline records a `Missing` with the
reason and what would resolve it. Three rules follow from that:

1. A DNA mutation never implies the protein is made. Without RNA or protein
   data, tumour expression is `unavailable` and the report states the chain:
   *DNA evidence present, RNA evidence unavailable, protein evidence
   unavailable, so surface expression cannot be established.*
2. An unknown dimension is dropped from the overall score, never defaulted.
   Component coverage is reported beside the score.
3. Safety only ever lowers a score. Missing normal-tissue data caps the overall
   at 0.6 and design readiness at 0.5; a target no mechanism can address caps
   design readiness at 0.15.

## Mechanism compatibility

Fifteen mechanisms, each declaring hard **gates** and weighted **factors**:

| Mechanism | Cargo | Effector | Needs internalisation |
|---|---|---|---|
| blocking_antibody, agonist_antibody | none | the target's own signalling | no |
| adcc, adcp, complement_recruitment | none | NK cells, macrophages, complement | no (penalised by it) |
| t_cell_engager, nk_cell_engager | none | T cell, NK cell | no |
| adc | cytotoxic payload | released payload | yes |
| targeted_radionuclide | radionuclide | emitted radiation | no |
| immune_marker_delivery | immunogenic marker | endogenous immune system | yes |
| rna_delivery | mRNA | translated program | yes |
| tumour_suppressor_restoration, genome_editing | coding RNA, editing machinery | restored protein, edited genome | yes |
| tcr_based, tcr_mimic | none | T cell | not applicable |

A failed gate sets compatibility to zero and prints the reason. Weighted
factors whose input is unknown are dropped and listed as blocking unknowns, and
the score is multiplied by an input-coverage factor so a mechanism cannot rank
highly on ignorance. Mechanism maturity is a separate, visible multiplier.

### Not refused is not established

A gate has three answers, not two: met, refused, and *unanswered*. A mechanism
whose hard requirement went unanswered is **provisional**: capped at 0.25,
listed with the requirement named, and never the candidate's preferred
mechanism. `best_mechanism` ranges only over mechanisms whose requirements were
answered and returns None when there are none; `best_provisional_mechanism`
gives the nearest open one with its open requirement, and the report prints
both.

This was the benchmark's last standing complaint and it closed two rows at
once. BRAF's plasma-membrane compartment is curated at 0.45, below the 0.6
reachability threshold, so the surface requirement was neither met nor refused;
every surface mechanism was scored provisionally at 0.25 and still headed the
list, because nothing else scored at all. The same for PIK3CA. Both now report
no preferred mechanism, and name `surface_accessible` as the reason.

The control that says this is a fix and not a silencing: give the same BRAF
case an HLA genotype and the preferred mechanism appears — `tcr_based` at 0.44,
the peptide/HLA route, which is the right answer for a protein no binder
reaches. EGFR and ERBB2 keep `blocking_antibody` unchanged, because their
surface requirement is answered rather than merely unrefused. It is the same
mistake as the other three the benchmark found: treating an unanswered question
as permission.

Everything marked experimental stays experimental. `immune_marker_delivery` is
the "mark the cell for the immune system" concept, and it is classed
experimental with the note that GenomeOS can cite no clinical precedent for the
full chain.

#### The same mistake, one layer down: a fusion is not the whole gene

The surface gate read the gene's curated localisation, which describes the
full-length protein. A fusion keeps one side of a junction. For EML4-ALK the
pipeline therefore offered a blocking antibody at **0.75** as the preferred
mechanism against what is, in the tumour, a cytoplasmic kinase: ALK is the 3'
partner and its extracellular domain is not in the product. There is no
approved antibody against it and there could not be — crizotinib, alectinib and
lorlatinib are small molecules that work inside the cell.

GenomeOS reconstructs neither the junction nor the partner orientation, so it
cannot say the ectodomain survives. For a candidate whose every origin is a
fusion the surface requirement is now **unanswered** rather than met: every
surface mechanism is provisional, and none heads the list. ALK stays at rank 1
with the fusion, its recurrent partner and its cohort frequency intact.

Half of the defect is left, and is pinned in
`test_the_rest_of_the_fusion_defect_is_recorded_rather_than_argued_away` rather
than argued away. The candidate is still classed `direct_surface` with an
accessibility of 1.0, because the class and the score also come from the gene
rather than from the product. Closing that needs a measurement the project does
not hold: the junction, the 5'/3' orientation, or transcript evidence for the
retained domains. The `.sv` format records a gene and a partner and nothing
else, so there is nowhere to read it from today.

### Repair is not elimination

Genome editing and tumour-suppressor restoration carry a mandatory caveat:

> Correcting or replacing one gene does not restore a normal cell. An
> established malignancy carries several cooperating drivers, passenger
> mutations, copy-number change, aneuploidy and epigenetic dysregulation;
> reverting one of them removes one dependency, not the disease.

Nothing hard-codes "killing is better". Elimination mechanisms tend to rank
higher because they have clinical precedent and fewer unmet requirements, and
both of those are visible in the contributions.

## Patient data levels

The report states which of nine levels the analysis actually had, so a
conclusion is never read above its data:

| Level | Input |
|---|---|
| 1 | tumour VCF |
| 2 | + matched normal genome |
| 3 | + tumour RNA-seq |
| 4 | + copy number and structural variation |
| 5 | + tumour proteomics |
| 6 | + surface proteomics |
| 7 | + HLA genotype |
| 8 | + immunopeptidomics |
| 9 | + single-cell tumour data |

### Alterations that are not point mutations

Level 4 is an input, not an annotation. A gene enters the candidate list
through a copy-number or structural call exactly as it does through a
variant, with the same origin record, the same cohort frequency and a score on
the same scale (`genomeos/cancer/alterations.py`, distilled knowledge in
`cancer_alterations_msk_impact_2017`). Three consequences, each a control in
`tests/test_cancer_alterations.py`:

- **An amplified oncogene with no mutation of its own is a candidate.** With a
  12-copy ERBB2 and no ERBB2 variant anywhere in the VCF, ERBB2 ranks first as
  a direct surface target. Without this path, the same inputs produced one
  candidate and it was not ERBB2.
- **A homozygously deleted gene is never a target.** It is classed
  `unsuitable`, dropped from the surface targets, and a hard requirement —
  `gene_product_present` — fails every mechanism in the ontology that has to
  recognise a product, because the tumour makes none of it. Only restoration
  and editing are exempt, and they fail in turn on delivery, which is the
  honest answer. A deep deletion is among the most actionable findings in a
  tumour and what it points at is the dependency the loss creates, which this
  pipeline does not model.
- **Copy number is still not expression.** A count of copies scores
  `(copies - 2) / 8` capped at 0.5; a discrete call with no count scores 0.2
  for an amplification and 0.1 for a gain, below what the weakest count
  producing that call would reach, because a call gives a direction and no
  amount. A fusion is recorded as a rearrangement between two genes and its
  junction sequence is not reconstructed, so no novel peptide is claimed
  from it — and, since 2026-09-21, no outward-facing epitope either.

These are unit controls, and for a long time they were the only evidence for
three of the four routes: every scored case in the therapeutic benchmark
reached its target through a point mutation, so the benchmark measured the
variant path end to end and the other three not at all. The seventh case, CD19,
closes that for expression. Copy number and structural variants are still
covered by controls alone, and the benchmark result says so in its note rather
than leaving the gap to be inferred from the case list.

## The Therapeutic Design Dataset

The bridge from cancer genomics to molecular design. It says what must be
recognised and, with equal weight, what must not be.

```
therapeutic_analysis.json        interpretation, ranked, with evidence
therapeutic_design_dataset.json  the complete machine-readable specification
target_positive_set.json         what a future binder SHOULD recognise
target_negative_set.json         what it MUST avoid
evidence_graph.json              variant -> transcript -> protein -> structure
                                 -> expression -> target -> mechanism
therapeutic_report.txt           the human-readable analysis
```

The negative set is the part that makes selectivity a solvable problem:

- the wild-type form of the same protein, when recognition is mutation-specific
- the wild-type peptide/HLA complex and similar self complexes, for a peptide target
- the same protein on every healthy tissue that expresses it, with the level
- the closest human paralogues with their sequence identity (Ensembl Compara)
- common polymorphic variants, flagged as not screened

Every specification keeps its genomic provenance: chromosome, position,
reference, alternate, transcript, protein change, variant type, somatic status,
VAF, copy number and clonality. The evidence graph makes the traversal work in
both directions.

`encoding_strategy_compatibility` is a modality classification only: whether a
candidate is conceptually compatible with an administered protein, an encoded
binder, a cell-expressed binder and so on. **No therapeutic nucleotide or
protein sequence, construct, vector, expression cassette, delivery formulation
or laboratory protocol is produced anywhere in this module.** Deciding what
molecule to build is a separate downstream stage.

## What the loss makes indispensable

A homozygously deleted gene is never a target — the section above says so and
`genomeos/cancer/alterations.py` enforces it — and the sentence that followed was
an admission: *what a deletion points at is the dependency the loss creates,
which this pipeline does not model.* It does now.

`genomeos/therapeutics/synthetic_lethality.py` asks one question over public
cell-line data:

```
cell lines that lost A  ->  do they need B more than the lines that kept A?
```

DepMap 24Q4 Public is the source (`genomeos/therapeutics/depmap.py`): Chronos
CRISPR gene effect, DepMap's own likely-loss-of-function matrix, PureCN absolute
copy number, and the MSIsensor2 score per line. The release is tens of gigabytes
and none of it is kept — each file is streamed over HTTP, the header gives the
column indices of the genes actually asked for, every row is discarded after
those few numbers are taken out of it, and what lands on disk is a gzipped
per-gene cache under `data/knowledge/therapeutics/depmap/`. About 1.3 GB is read
over the wire for this sweep and 2 MB of it is kept — gene effect 1.5 MB,
absolute copy number 428 KB, damaging mutations 60 KB, the signatures and the
model table 92 KB — and `--offline` then reproduces the run from that cache
alone. The portal API is behind a bot check, so files are addressed by their
figshare ids, which is the citable form anyway.

### Declared before the run

| | |
|---|---|
| test | Mann-Whitney rank sum, one-sided, tie- and continuity-corrected |
| why a rank test | gene effect is not normal, the lost group is small, and one hypersensitive line would carry a mean |
| multiple testing | Benjamini-Hochberg across every tested pair; `q <= 0.10` |
| group sizes | at least 5 lost and 20 intact lines, otherwise `untestable` rather than negative |
| loss | a damaging mutation (DepMap `LikelyLoF`) **or** absolute copy number below 0.5; *kept* needs both calls present, because a missing copy-number value cannot rule out a deletion |
| MSI-high | MSIsensor2 >= 20 |
| candidate space | 206 DNA-repair and replication-stress genes in 12 pathways, all 42,230 ordered pairs, plus 318 paralogue pairs from the project's own Ensembl Compara cache (206/206 genes resolved, identity >= 20%) |

The seven pre-registered controls are in `PREREGISTERED`, each with the reason it
is there — including two that were expected to **fail**: BRCA1 and BRCA2 loss
against PARP1 dependency, the textbook synthetic-lethal pair clinically and
famously weak in CRISPR knockout screens, because a PARP inhibitor traps PARP1 on
DNA and deleting PARP1 does not reproduce that.

### The run

19,749 of the 42,548 declared pairs were testable (171 loss genes reached 5 lost
lines; 6,830 pairs are reported as untestable with the counts that made them so),
over 1,178 lines with gene effect, 1,929 with mutation calls, 1,607 with absolute
copy number and 1,929 with an MSI score, 112 of them MSI-high. **44 pairs at
q <= 0.10.**

**The control for the sweep itself: 20 permutations of the loss labels, group
sizes preserved, produced 0 hits every single time** (mean 0.0, max 0). The null
is flat, so the 44 are read against zero and not against a floor.

| Control | Expected | Lost lines | Cliff's delta | median gene-effect difference | p | Verdict |
|---|---|---|---|---|---|---|
| ARID1A damaging -> ARID1B | recovered | 118 | -0.421 | -0.135 | 1e-14 | **recovered** |
| SMARCA4 damaging -> SMARCA2 | recovered | 59 | -0.424 | -0.114 | 2e-08 | **recovered** |
| MTAP deletion -> PRMT5 | recovered | 175 | -0.367 | -0.162 | 4e-14 | **recovered** |
| MTAP deletion -> MAT2A | recovered | 175 | -0.188 | -0.144 | 5e-05 | **recovered** |
| MSI-high -> WRN | recovered | 74 | -0.467 | -0.370 | <1e-8 | **recovered** |
| BRCA1 damaging -> PARP1 | weak | 21 | -0.271 | -0.047 | 0.017 | recovered, against the pre-registration |
| BRCA2 damaging -> PARP1 | weak | 36 | -0.213 | -0.028 | 0.015 | recovered, against the pre-registration |

All five positive controls come back, the strongest of them the MSI-high line
state against WRN — the only control whose "loss" is not a gene at all — where
MSI-high lines sit at a median WRN gene effect of -0.42 against -0.05 for the
microsatellite-stable lines. So 5 of 7 controls agree with what was written down,
and the two that do not are the two the pre-registration told us to read
carefully. They are worth
the space, because they are a lesson about verdict rules rather than about PARP1.
The rank shift is real and it is nothing: BRCA1-mutant lines sit at a median
PARP1 gene effect of **-0.215** against **-0.168** for the rest, on a scale where
0 is no effect and -1 is the median common essential. Neither group depends on
PARP1. With 1,178 lines a 0.047 shift clears a significance threshold comfortably,
and a verdict rule made only of p-values calls that "recovered" — which is why
the effect size sits in the same table and why the pre-registration asked for an
honest report rather than a pass mark. Under the sweep-wide correction across
19,749 tests both pairs fall back to `nominal` (q = 0.457).

### The 44, by the twelve strongest

| Loss | Dependency | Lost lines | Cliff's delta | q |
|---|---|---|---|---|
| MTAP | WDR77 | 177 | -0.441 | 0.0 |
| ARID1A | ARID1B | 121 | -0.428 | 0.0 |
| MTAP | PRMT5 | 177 | -0.360 | 0.0 |
| ARID1A | WRN | 121 | -0.339 | 0.0 |
| CDKN2B | WDR77 | 274 | -0.286 | 0.0 |
| CDKN2A | WDR77 | 427 | -0.245 | 0.0 |
| SMARCA4 | SMARCA2 | 63 | -0.443 | 1e-05 |
| EP300 | CREBBP | 81 | -0.393 | 1e-05 |
| KMT2D | RAD50 | 103 | -0.324 | 8e-05 |
| KMT2D | WRN | 103 | -0.323 | 8e-05 |
| MSH3 | WRN | 32 | -0.496 | 0.0015 |
| RECQL | WRN | 10 | -0.823 | 0.0052 |

Three readings, and two of them are warnings.

**WDR77 above PRMT5 is the sweep checking itself.** WDR77 is MEP50, PRMT5's
obligate partner in the same methylosome, and it was in the candidate list as a
methionine-salvage gene rather than as a known answer. A sweep that recovers a
complex's two halves with the same loss, in the right order of effect size, is
measuring the complex and not the annotation.

**CDKN2A and CDKN2B against WDR77 are the same event as MTAP, not three
findings.** MTAP sits beside CDKN2A and CDKN2B at 9p21 and is deleted with them:
of the 287 lines this sweep calls MTAP-lost, **273 are also called CDKN2A-lost
and 260 CDKN2B-lost**, while CDKN2A-loss covers 653 lines in all. A pair test
cannot say which of two genes lost together carries the dependency, so the result
lists them under `loss_genes_that_are_the_same_event` with the overlap, and
`confounders` names the three other ways a hit here can be right about the
statistics and wrong about the mechanism — co-deletion, lineage composition, and
a p-value standing in for an effect size.

**RECQL -> WRN rests on 10 lines.** It has the largest effect in the table and
the least evidence behind it, which is what a group size of 10 looks like even
after a rank test and an FDR. EP300 -> CREBBP and the ARID1A/SMARCA4 rows are
the paralogue pattern the candidate space was built to find; WRN appearing under
four different losses (ARID1A, KMT2D, MSH3, RECQL) is consistent with WRN's known
dependence on replication-stress and mismatch-repair state and is the part of
this list most worth an independent look.

One more thing the controls earned. On the first complete run the MSI-high
control came back `untestable`, and the reason was not biology:
`OmicsSignatures.csv` leaves its index column unnamed, the table was read by
`ModelID`, every row was dropped, and an empty MSI table was cached in silence.
A pre-registered control that *should* be the cleanest hit in the sweep is what
made a silent parse failure visible. The key now falls back to the first column
and an empty table raises instead of being cached, with a test for both.

Every hit is emitted as `prediction` evidence at `computational` level with
confidence capped at 0.7, because a knockout in culture is evidence about a
target and not about a patient. **A dependency is a hypothesis about a target. No
molecule, construct, formulation or protocol follows from it here.**

## Marker logic: A, A and B, A and not C

`--scan` scores one antigen at a time, and the section above shows what that
costs: CD20 turned down at 2.2x because healthy spleen carries the same lineage.
The refusal is right about the single marker and wrong about the target space,
because a cell can be addressed by a *combination* — two antigens that must both
be present, or one that must be present while another, carried by the healthy
tissue at risk, must be absent. Bispecifics, logic-gated CARs and dual-antigen
engagers are all built on that, and a pipeline that scores one gene at a time
cannot propose any of them.

`genomeos/therapeutics/selectivity.py` scores gates instead of genes, over the
5,588-gene packaged membrane and CD-marker atlas, with DepMap 24Q4 cell-line RNA
as the tumour side (`scripts/marker_selectivity.py`, result
`marker_selectivity`; the 507 MB expression matrix is streamed the same way and
distils to a 44 MB cache of those 5,588 genes). Two numbers are kept apart and
never merged:

| | What it is | Source |
|---|---|---|
| tumour coverage | fraction of one lineage's cell lines where every required marker is present at log2(TPM+1) >= 3 and no excluded marker is | DepMap `OmicsExpressionProteinCodingGenesTPMLogp1` |
| healthy burden | the worst of the eight essential organs the gate still reaches, in nTPM | the packaged HPA atlas |

`score = coverage / (1 + burden / 10 nTPM)`. Two rules inside the burden do the
real work. An AND gate falls towards its lower marker in each tissue but **never
below a quarter of the higher one**, because bulk tissue RNA cannot say whether
two markers sit on the same healthy cell; without that floor, pairing any target
with a transcript the atlas reports as absent drove the burden to zero and
thousands of gates tied at the top on an absence twenty bulk tissues cannot
establish. A NOT marker, by contrast, vetoes its tissue outright, because a
measured presence is evidence in a way an absence is not. The 137 lines DepMap
itself calls non-cancerous are dropped before anything is scored; an earlier pass
kept them and let 42 fibroblast lines supply "tumour coverage", and that pass was
discarded without its benchmark verdict being read.

154,991 gates over 1,588 lines in 23 lineages. **37 of the top 200 gates beat
their own best single marker**, and every one of them is an `A and not C` gate:
the veto is what buys selectivity here, not the second required antigen, which
the co-expression floor caps at a four-fold gain.

### The pre-registered claim, and its verdict: not supported

Declared in code before the ranking was read: gates containing one of 16
approved or clinical surface targets (ERBB2, EGFR, CD19, MS4A1, CD22, TNFRSF17,
CD38, CD33, TACSTD2, FOLH1, MSLN, CEACAM5, DLL3, CLDN18, GPC3, CD70) must
outrank random gates from the same pool, at p <= 0.05, Cliff's delta >= 0.10, and
above the 95th percentile of a placebo null.

**It fails, and not narrowly.** The 9,778 gates containing an approved target sit
at median percentile 25.1 against 51.5 for 2,000 random gates: Cliff's delta
**-0.417**, one-sided p = 1.0. Read against zero, the screen ranks the targets of
oncology *below* arbitrary pairs of membrane proteins.

### The control is not neutral, and that changes the reading

The placebo null draws 200 gene sets of 16 from markers forced into the pool the
same way the benchmark genes are, and it is **not flat**: its mean delta is
**-0.496**. Half a Cliff's delta of the benchmark's deficit belongs to the
comparison, not to the targets — a marker carried into the pool is not
comparable to one that scored its way in. Against that null the benchmark set
sits at the **87th percentile**, +0.078 above the null mean, still short of the
declared 95th-percentile bar.

A post-hoc sensitivity analysis says where the rest of the deficit lives. Raising
the presence threshold — a marker at 7 TPM is not the same object as a marker at
60 TPM, and an approved target is always the second kind — widens the gap
monotonically and clears the null's own bar at every stricter threshold:

| present at log2(TPM+1) >= | benchmark delta | null mean | gap | above null p95 |
|---|---|---|---|---|
| 3.0 (declared) | -0.417 | -0.496 | +0.078 | no |
| 4.0 | -0.359 | -0.525 | +0.166 | yes |
| 5.0 | -0.276 | -0.584 | +0.308 | yes |
| 6.0 | -0.191 | -0.593 | +0.402 | yes |

The top of the ranking says the same thing from the other end. ZACN, a
ligand-gated ion channel subunit, heads the single markers on 96% of liver lines
at log2(TPM+1) of 3.0 to 5.3 — about 8 to 40 TPM — with 0.0 nTPM across the eight
essential organs. It is a real measurement and a useless target: the score
rewards *absence from eight organs plus barely-detectable tumour presence*, which
is a description of a transcript near the detection floor rather than of an
antigen a binder could see.

### What the approved targets actually do in the screen

Each gene's *best* gate ranks far better than its average one, which is the
pre-registration's own mistake made visible: a target programme picks a target
with its best partner, not its mean partner.

| Target | Best gate | Lineage | Coverage | Burden | Percentile |
|---|---|---|---|---|---|
| EGFR | EGFR and not C | Pleura | 0.91 | 0.0 | 99.9 |
| ERBB2 | ERBB2 and not C | Pleura | 0.91 | 1.3 | 99.5 |
| CD70 | CD70 and not C | Lymphoid | 0.81 | 0.3 | 99.3 |
| MSLN | MSLN and not C | Pleura | 0.86 | 3.0 | 96.5 |
| CEACAM5 | CEACAM5 and not C | Bowel | 0.65 | 0.6 | 93.9 |
| CD33 | CD33 and not C | Myeloid | 0.89 | 5.1 | 92.1 |
| CD19 | CD19 and not C | Lymphoid | 0.54 | 0.3 | 83.9 |
| CLDN18 | CLDN18 and not C | Esophagus/Stomach | 0.16 | 0.2 | 31.6 |
| FOLH1 | FOLH1 and not C | Skin | 0.14 | 0.6 | 29.3 |

Every one of them is an `A and not C` gate, and the excluded marker of each is
named in `best_per_benchmark_gene` rather than summarised here, because which
tissue marker does the vetoing is the whole content of the proposal.

On that statistic the 16 benchmark markers reach median percentile **84.3**
against **71.4** for the 64 placebo markers — delta 0.133, p = 0.21 with 16
against 64, which is a direction and not a result. It is reported as post-hoc
because it was written after the pre-registered verdict was read.

Two further honest failures are visible in that table. CLDN18 and FOLH1 rank low
because cell lines of the relevant lineage mostly do not express them — a
property of the model system, since claudin-18.2 and PSMA are expressed in the
tumours the approved drugs treat. And EGFR, ERBB2 and MSLN all place their best
gate in **Pleura**, not in the lineage their drug is approved for, because the
screen takes each gate's best lineage and pleural mesothelioma lines happen to
carry all three. The screen measures a marker's reach across cell lines, and it
does not know what disease anyone is treating.

### What this says about the method

The negative result is specific and useful: a coverage-versus-essential-organs
score over cell-line RNA is not a target ranking. Three things are missing, and
naming them is the deliverable rather than a ranked list nobody should act on —
abundance (a threshold crossing is not a surface density), the other twelve
healthy tissues (a lineage antigen's cost is in spleen, colon and skin, none of
them essential organs, which is exactly why CD19's own cost is B-cell aplasia),
and protein (both sides of this screen are transcripts). Every candidate is
emitted as `prediction` evidence at `computational` level with confidence capped
at 0.5 for that reason, and naming the marker logic is the whole output: no
binder, no format, no cargo, no protocol.

## Use

```bash
# the full analysis, the report and the five dataset layers
genomeos therapeutic --tumour data/demo/cancer_tumour.vcf \
  --hla "HLA-A*02:01,HLA-B*07:02" --out out/ --report

# with a cohort of the same cancer type
genomeos cancer expression --study brca_tcga_pan_can_atlas_2018
genomeos therapeutic --tumour tumour.vcf --cohort brca_tcga_pan_can_atlas_2018

# one target's specification
genomeos therapeutic --tumour data/demo/cancer_tumour.vcf --spec RUNX1

# richer patient input
genomeos therapeutic --tumour tumour.vcf --normal normal.vcf \
  --rna tumour_tpm.tsv --cnv copy_number.tsv --purity 0.7 --hla "HLA-A*02:01"

# the driver is an amplification, a deletion or a fusion, not a mutation
genomeos cancer alterations --study msk_impact_2017      # the cohort, once
genomeos therapeutic --tumour t.vcf --cnv gistic.tsv --cna-format gistic --sv fusions.tsv
genomeos cancer tumour --vcf t.vcf --cnv copy_number.tsv --sv fusions.tsv

# offline, against the caches only
genomeos therapeutic --tumour data/demo/cancer_tumour.vcf --offline

# from the tumour-only cancer pipeline
genomeos cancer tumour --vcf data/demo/cancer_tumour.vcf --therapeutic

# the dependency sweep and the marker-logic screen: DepMap is streamed and discarded,
# only the requested gene columns are cached, and --offline runs on that cache alone
uv run python scripts/synthetic_lethality.py
uv run python scripts/marker_selectivity.py
```

```python
from genomeos.therapeutics import analyse_vcf, text_report, write_outputs

analysis = analyse_vcf("data/demo/cancer_tumour.vcf", hla=["HLA-A*02:01"])
print(text_report(analysis))
write_outputs(analysis, "out/")
```

`POST /api/therapeutics {"vcf": ..., "hla": ..., "rna": ..., "top": 8}` returns
the same analysis, and the **Targets** tab in `genomeos serve` renders it.

## Sources

| Source | Used for | Licence |
|---|---|---|
| Ensembl VEP | consequence, gnomAD frequency, COSMIC ids | Apache 2.0 service, EMBL-EBI terms |
| UniProtKB/Swiss-Prot | sequence, subcellular location, topology, keywords, processing | CC BY 4.0 |
| InterPro, PDB, AlphaFold, Reactome, STRING | domains, structures, pathways, associations | open, per-source terms |
| Human Protein Atlas | consensus healthy-tissue RNA, immunofluorescence location | CC BY-SA 4.0 |
| Open Targets Platform | tractability, approved drugs and trials, safety liabilities | CC0 1.0 |
| Ensembl Compara | human paralogues with sequence identity | Apache 2.0 service |
| Ensembl sequence | transcript coding sequences, for reconstructing an altered protein | Apache 2.0 service |
| cBioPortal expression profiles | cohort tumour-versus-normal behaviour per gene | study terms, portal ODbL |
| cBioPortal (distilled) | driver frequencies and hotspots | ODbL, distilled summary committed |
| DepMap 24Q4 Public (figshare) | CRISPR gene effect, damaging mutations, absolute copy number, MSI score, cell-line RNA | CC BY 4.0, doi:10.25452/figshare.plus.27993248.v1 |
| MHCflurry 2 (optional) | peptide/HLA class I binding, percentile rank | Apache 2.0; models downloaded by the user |
| NetMHCpan 4.1 (optional) | peptide/HLA class I binding, percentile rank | academic licence from DTU; never shipped or downloaded by GenomeOS |

Every response is cached under `data/knowledge/therapeutics/`, so an analysis
is reproducible without the network. No page is scraped where an API exists.

Two peptide/HLA predictors can be enabled, and neither is assumed
(`genomeos/therapeutics/binding.py`). **MHCflurry 2** (Apache 2.0) installs with
`uv sync --extra hla` plus `mhcflurry-downloads fetch models_class1_presentation`
and is the default when present. **NetMHCpan 4.1** is the field's reference
predictor under an academic licence from DTU: GenomeOS never ships, downloads or
redistributes it, and calls a binary the user installed, found on `PATH` or at
`$NETMHCPAN`. Choose with `--hla-predictor mhcflurry|netmhcpan|none`.

With neither installed, `NoNeoantigenPredictor` reports binding as unavailable
rather than inventing affinities; the analysis still runs and says what is
missing. With one installed, every mutation-spanning peptide is scored against
the patient's class I alleles and reported as a percentile rank (binder at
<= 2.0, strong binder at <= 0.5) with `computational` evidence. Binding is still
not presentation: proteasomal cleavage, TAP transport and immunopeptidomics stay
unestablished, so presentation confidence remains capped at 0.25 on binding
alone.

## What this is not

GenomeOS does not say a therapy will work. It says an alteration creates a
candidate vulnerability worth further investigation, and then says exactly how
much of that is established, predicted, or missing.
