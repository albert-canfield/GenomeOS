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

The translation check ran on every chromosome as it compiled: of
19,249 compiled genes, 90.9% translate from the
canonical transcript to exactly the reviewed UniProt sequence, 97.8%
match some annotated isoform exactly, and 99 disagree beyond
boundary differences (`translation_vs_uniprot_<chrom>.json` lists each).

The knowledge graph built from the compiled definitions
(`graph_genome.json`): 41,982 nodes (19,797 proteins,
2,302 pathways, 19,846 domains, 37 tissues) and
327,024 edges, 95,098 of the associations with experimental
support; 15,165 of 19,283 compiled proteins sit in one connected
component. The compiled definitions occupy about 300 MB locally and are
not committed; the summaries are.

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
