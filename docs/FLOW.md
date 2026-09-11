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
