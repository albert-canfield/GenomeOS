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
This is the measurable form of "how much of the proteome do we know", and
the first milestone of the protein layer:

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
