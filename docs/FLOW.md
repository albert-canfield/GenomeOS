# The upward flow, and the types that carry it

**Status on 2026-10-10.** Current. The layer table and the four commands hold:
`flow`, `lookup`, `rna` and `regulation` are all parsers in `genomeos/cli.py`, and
`genomeos/flow/trace.py` is still the single home of the DNA–RNA–protein relationship
(`trace`, `residue_of`, `genomic_of_residue`, `substitute`). The regulation confidences
quoted below are the ones the code uses (`genomeos/genome/regulation.py`: promoter of, 0.8;
nearest TSS in domain, 0.4; same domain, 0.25), and node calls are still inferred at 0.4
(`genomeos/genome/domains.py`). One line has moved on: "Hi-C would make it curated" now has
a route — `genomeos/genome/hic.py` streams measured 4D Nucleome boundary calls and loads
disabled without an account key (docs/DATA.md, "Optional: 4D Nucleome boundary calls").

Genetic information flows upward:

```
DNA sequence → regulation → RNA → proteins → cell behaviour → tissues → organs → organism
```

Every layer has a home in GenomeOS, and `genomeos flow GENE --chrom chrN`
walks one gene through all of them with the evidence at each step.

| layer | what carries it | where in GenomeOS | evidence |
|---|---|---|---|
| DNA sequence | bases, loci, chromosomes, an individual's alleles | `Sequence`, `Locus`, `Chromosome`, `Genome`, `IndexedGenome`, `Variant` (VCF) | curated (assembly), measured (a person's VCF) |
| regulation | promoters, enhancers, CTCF boundaries, CpG islands, motifs that name required factors | `RegulatoryElement` with targets bounded by the CTCF domain (`genomeos regulation`), ENCODE cCREs, pattern file, learned signals | curated (ENCODE), inferred (targets), predicted (patterns) |
| RNA | transcripts and isoforms, exons spliced to mRNA, non-coding RNA classes, where the gene is expressed | `Transcript` (GENCODE), `splice()`, `transcribe()`, `genomeos rna GENE` (isoforms with lengths and tags; GTEx median TPM in 54 tissues), anatomy counts of ncRNA | curated; measured (GTEx) |
| protein | translation with the right genetic code and initiators, UniProt record, AlphaFold structure | `translate_transcript()`, `Protein`, molecules layer | curated (UniProt), predicted (AlphaFold) |
| cell behaviour | rules with context, libraries the gene belongs to, cell-type expression | `Rule`, `CellType`, `Module.active_rules()`, `KnowledgeBase.libraries_of()`, network runtime | curated (GO/Reactome), inferred |
| tissues, organs | blueprint and systems libraries, spatial development, ligand–receptor protocol | `blueprint.*`, `systems.*`, spatial/segmentation/gastrulation engines | curated, inferred |
| organism | lineage from a zygote, timers, twins, whole-genome inventory | lineage engine, ageing runtime, BioTwin, anatomy | measured (twin), inferred |

## Your type list against BioIR

| you | GenomeOS today | note |
|---|---|---|
| Base, BasePair | characters of `Sequence`; `complement()` gives the pair | no separate class: a base is a position, not an object, at this scale |
| Sequence | `Sequence` | validated IUPAC, reverse complement, GC, repeats |
| Locus | `Locus` | 0-based, half-open, stranded |
| Chromosome, Genome | `Chromosome`, `Genome`, `IndexedGenome` | linear now; pangenome graph is the design target |
| Allele | `Variant.alts` + haplotypes from `apply_variants` | worth promoting to its own BioIR entity when twins carry both haplotypes as first-class objects |
| Gene, Transcript, Protein | `Gene`, `Transcript`, `Protein` | from GENCODE/Ensembl; protein sequence by our translation, checked against UniProt |
| RegulatoryElement | `RegulatoryElement` (BioIR) with class, locus, the node that bounds it and `targets` with a basis each; BioLang `element` block; `genomeos regulation --gene G` | ENCODE gives the elements (curated); the targets are inferred: a promoter reaches the TSS within 1 kb (0.8), an enhancer reaches the coding genes of its CTCF domain, nearest TSS first (0.4, others 0.25) |
| Domain (node) | `genomeos.genome.domains.Domain`: intervals between CTCF-only boundaries, with the genes, promoters and enhancers inside; a lane in the block map | inferred, 0.4; Hi-C would make it curated |
| CellState | `CellState` (ageing runtime) and `CellType` (identity) | state and identity are kept separate on purpose |

## Where the DNA–RNA–protein relationship lives

The relationship has one home: `genomeos/flow/trace.py`. `trace(genome,
transcript)` builds a `CentralDogmaTrace` that maps every genomic base to
its position in the spliced mRNA, its codon and its residue, and back
(`residue_of(pos)`, `genomic_of_residue(n)`), and can trace one base change
upward (`substitute(pos, ref, alt)` → synonymous / missense / nonsense /
start lost, with HGVS c. and p. names). Everything above uses it:

| consumer | what it does with the trace |
|---|---|
| `genomeos flow GENE` | prints the walk with sizes and evidence |
| `/api/flow`, the **Flow** tab | three lanes to scale (DNA with exons, mRNA with UTRs and CDS, residues by chemistry), linked on hover; a variant box shows what one base does to the protein |
| Molecules tab | checks our translation against UniProt, then shows the compiled protein definition (see docs/PROTEIN.md) |
| variant classifier | the coding consequence is the same computation |

The regulation layer is drawn above the DNA lane in the Flow tab. APP sits
alone in its node (chr21:D91, 175 kb): one promoter element, 98 enhancers
can reach it, 89 of them inside the gene, two CTCF insulators bound the node.

Checked on APP: chr21:25,897,620 C>T on the plus strand of a minus-strand
gene reads G>A in the mRNA, codon 673 GCA→ACA, p.Ala673Thr, the protective
Icelandic variant.

## One variant through every layer

`genomeos lookup chr17:7675088 C>T` (any chromosome, 1-based, plus strand),
`/api/lookup`, and the "Look a variant up" card in the Flow tab put the
layers together for one change:

| layer | what comes back | evidence |
|---|---|---|
| consequence | Ensembl VEP on the canonical transcript: effect, gene, HGVS c. and p., SIFT, PolyPhen | curated / predicted |
| known | ClinVar significance, dbSNP and COSMIC ids, gnomAD frequency, PubMed citations; or "novel" | curated |
| local trace | the GenomeOS trace of the same change when the chromosome is local, checked against VEP | derived |
| protein | UniProt features at the residue (domains, regions, sites, cleavage products, modifications), structures available | curated |
| pathways | for a truncating change, the Reactome reactions lost with the protein absent | curated + inferred |

TP53 R175H comes back pathogenic in ClinVar with 65 citations and 311
experimental structures; APP A673T comes back protective, at the residue
where the amyloid-beta chains begin (672), which is the mechanism: the
change weakens beta-secretase cleavage. The local trace agrees with VEP on
both.

## From ATGCG… to a protein: the executable piece

This is the first piece of biology GenomeOS could run, and it is fully
tested: `transcribe` (strand-aware), `splice` (exons in transcript order),
`coding_sequence` (phase-trimmed CDS), `translate` (standard and vertebrate
mitochondrial codes, initiator codons AUU/AUA/AUC/GUG read as Met). Proof:
all 13 mitochondrial proteins at their published lengths from chrM DNA;
2,703 chromosome-21 transcripts translate with no internal stop; APP-201
translated from chromosome 21 is 100% identical to UniProt P05067.

Run it on any gene:

```
genomeos gene MT-CO1 --gff3 data/results/gencode_v50_chr21_chrM.gff3.gz --genome data/reference/chrM.fa.gz --chrom chrM --table mito
genomeos flow APP --chrom chr21
```

## Individuals on the flow (2026-09-11)

Below the DNA lane the Flow tab draws one row of ticks per local individual
(the test human and every genome imported with `genomeos individual import`):
every variant of that person inside the gene as a grey tick, protein-changing
coding SNVs in red, and the summary lists the changes with their genotype;
clicking one traces it through DNA → RNA → protein like any variant typed
into the box. Genotypes are measured (each person's file); consequences are
derived by the local trace on the canonical transcript. Served in `/api/flow`
as `individuals`.

## The trace on the tissue's own isoform (2026-10-10)

The consequence trace used to follow the canonical transcript only, so a
variant's predicted consequence "in the cerebellum" could describe APP770
when the cerebellum makes APP695. `trace_gene_in_tissue` (genomeos/flow/tissue.py)
now takes a tissue and traces that tissue's dominant isoform beside the
canonical one; `TissueTrace.substitute` returns the consequence on the
tissue's isoform, the canonical consequence beside it, and how they compare
(same, class_changed, residue_changed, change_differs).

The tissue's dominant isoform comes from the cached GTEx v8 median transcript
TPM (`data/knowledge/expression/gtex_isoforms_*.json`, read only, never
fetched). It is ESTABLISHED when GTEx measured the tissue, the top isoform
reaches 1 TPM, the next isoform is below 0.9 of the top, and the top isoform
codes a protein in the local GENCODE v50 models. Otherwise the call is
UNKNOWN with its reason (no_isoform_data, tissue_not_measured, low_expression,
tie, dominant_not_in_local_models, dominant_non_coding) and the consequence is
UNKNOWN: the trace never falls back to the canonical transcript. The isoform is
GTEx's median across donors, not one person's own expression, and both
consequences are predictions of the local trace.

**Registered census (`data/results/flow_tissue_census_registration.json`,
written before any count).** Over HG002's own coding SNVs in the genes the
project holds GTEx isoform tables for (99 tables; 96 were fetched for the
translation-disagreement genes, plus APP, FUT2 and TP53, so this is the
project's set and not a sample of the genome), how many change consequence
class or residue number between the canonical transcript and each tissue's
dominant isoform, per tissue, with denominators. A census of disagreement
between two predictions, not an accuracy claim; no prediction is registered.
The script (`scripts/flow_tissue_census.py`), tissue.py and trace.py are
pinned by sha256 and the tables by digest; the run refuses if any moved. The
committed result carries counts only: no position, allele, genotype or protein
change of the person.

**Census result (`data/results/flow_tissue_census_result.json`).** Over 299
coding SNVs of HG002 in 46 genes with cached GTEx isoform tables, the
consequence predicted on a tissue's dominant isoform differs in class or
residue number from the one predicted on the canonical transcript in 551 of
5,705 variant-tissue pairs where the tissue's isoform is established (238
class changes, 313 residue-number changes, 0 other-amino-acid-only); 10,441
of 16,146 pairs are UNKNOWN (low expression 9,554, dominant non-coding 649,
tie 215, dominant not in the local models 23). Every one of the 54 tissues
has at least 3 disagreeing variants (median 10 of a median 127 established;
testis 32 of 156, cerebellum 13 of 37, whole blood 3 of 72); 11 of the 46
genes carry all 551. All 299 are coding on the canonical transcript; none is
coding only on a tissue's isoform. Both consequences are predictions of the
local trace on GTEx's population-median isoform, not measurements and not
this person's own expression; the genes are the project's cached set (mostly
translation-disagreement genes), not a sample of the genome.
