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
  from it.

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
