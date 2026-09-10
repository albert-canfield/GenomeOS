# The upward flow, and the types that carry it

Genetic information flows upward:

```
DNA sequence → regulation → RNA → proteins → cell behaviour → tissues → organs → organism
```

Every layer has a home in GenomeOS, and `genomeos flow GENE --chrom chrN`
walks one gene through all of them with the evidence at each step.

| layer | what carries it | where in GenomeOS | evidence |
|---|---|---|---|
| DNA sequence | bases, loci, chromosomes, an individual's alleles | `Sequence`, `Locus`, `Chromosome`, `Genome`, `IndexedGenome`, `Variant` (VCF) | curated (assembly), measured (a person's VCF) |
| regulation | promoters, enhancers, CTCF boundaries, CpG islands, motifs that name required factors | `Region` / UNKNOWN blocks, ENCODE cCREs, pattern file, learned signals | curated (ENCODE), predicted (patterns) |
| RNA | transcripts and isoforms, exons spliced to mRNA, non-coding RNA classes | `Transcript` (GENCODE), `splice()`, `transcribe()`, anatomy counts of ncRNA | curated |
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
| RegulatoryElement | `Region` (role may be UNKNOWN) and ENCODE cCRE blocks | the next BioIR entity to add explicitly, with class (promoter, enhancer, insulator) and the genes it reaches; then `Domain` above it |
| CellState | `CellState` (ageing runtime) and `CellType` (identity) | state and identity are kept separate on purpose |

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
