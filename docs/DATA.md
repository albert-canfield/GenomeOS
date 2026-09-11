# Open genome data for testing

All of the following are open access and were verified reachable on
2026-09-10. Fetch with `scripts/fetch_reference.sh <target>`; files land in
`data/reference/` (gitignored).

## A real individual: HG002 (recommended test subject)

HG002 (NA24385) is the son of the Genome in a Bottle Ashkenazi trio and a
Personal Genome Project participant with open consent, so the data are usable
in research and commercially. He is the best-characterised human genome in
existence.

| Target | What | Size | URL |
|---|---|---|---|
| `hg002` | T2T diploid assembly v1.1 (both haplotypes, telomere to telomere) | ~1 GB | https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz |
| `hg002-vcf` | GIAB v4.2.1 benchmark variants against GRCh38 | ~100 MB | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz |
| `hg001-vcf` | NA12878 / HG001 benchmark variants (the classic CEPH genome) | ~100 MB | https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz |

Raw reads (Illumina, PacBio HiFi, ONT) for the same individuals are under
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/ when read-level
work (telomere length estimation, methylation from ONT) begins.

## Reference assemblies

| Target | What | Size | URL |
|---|---|---|---|
| `chrM` | hg38 mitochondrial genome (16,569 bp) | 6 KB | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz |
| `chr21` | hg38 chromosome 21 (46.7 Mb) | 13 MB | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz |
| `chrN` | any hg38 chromosome | 10–80 MB | same pattern |
| `hg38` | GRCh38 full assembly | ~940 MB | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz |
| `hs1` | T2T-CHM13 v2.0, gapless | ~900 MB | https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz |

## Annotation and knowledge

| Target | What | URL |
|---|---|---|
| `gencode` | GENCODE Release 50 GFF3 (genes, transcripts, exons, CDS) | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gff3.gz |
| — | Gene Ontology basic | https://current.geneontology.org/ontology/go-basic.obo |
| — | Reactome V96 (SBML, BioPAX) | https://reactome.org/download-data/ |
| — | Uberon anatomy, Cell Ontology | http://purl.obolibrary.org/obo/uberon/uberon-ext.obo , http://purl.obolibrary.org/obo/cl.obo |
| — | ClinVar VCF (weekly) | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz |
| — | dbSNP b157 | https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz |

## Populations

| Target | What | URL |
|---|---|---|
| `1kg-chr21` | 1000 Genomes, 3,202 samples, chr21 genotypes | https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/ |
| — | HPRC Release 2 pangenome (GFA / GBZ) | https://github.com/human-pangenomics/hpp_pangenome_resources |
| — | CELLxGENE Census (single-cell expression, ~33M cells) | `pip install cellxgene-census` |

## Verified today

```
genomeos info data/reference/chrM.fa.gz        # 16,569 bp
genomeos orfs data/reference/chrM.fa.gz --table mito --min-aa 200
    -> MT-CO1 513 aa, MT-CO2 227 aa, MT-ATP6 226 aa (match published lengths)
genomeos info data/reference/chr21.fa.gz       # 46,709,983 bp in 0.5 s
```
Known limitation: the ORF finder only starts at AUG. Human mitochondrial
genes ND1, ND2 and ND5 use AUA/AUU starts, so their ORFs are found at the
first internal AUG. Alternative start codons are on the Phase 1 list.

## Optional: AlphaGenome (predicted regulatory effects)

Google DeepMind's AlphaGenome predicts what a sequence change does to
expression, splicing and chromatin in each tissue. It is optional: the core
never depends on it, and its output enters GenomeOS only as `predicted`
evidence capped at confidence 0.7. Without it the features below are loaded
but disabled, and `genomeos predict --status` and the Progress tab say why
and how to enable them.

| | Feature | What the key enables | State |
|---|---|---|---|
| a | variant effect | predicted expression change per tissue for any variant: `genomeos predict chr21:25897620 C>T`, `/api/predict` | built |
| b | regulatory blocks | predicted target gene and tissue for UNKNOWN regulatory elements, tightening the CTCF-domain inference | planned |
| c | sequence grammar | splice-site and chromatin signals the learned matrices cannot capture, entering as predicted evidence | planned |
| d | twin | predicted expression differences between an individual's haplotypes and the reference | planned |

What each feature is worth and how it fits the architecture:
[ALPHAGENOME.md](ALPHAGENOME.md).

To enable:

1. Request a key (free for non-commercial use) at
   https://deepmind.google.com/science/alphagenome/account/terms.
2. Install the client: `uv sync --extra predict`.
3. Put the key in a file named `.env` at the project root (git-ignored,
   never committed) or export it in your shell:

```
ALPHAGENOME_API_KEY=your-key-here
```

4. Check: `uv run genomeos predict --status` prints `enabled`.

The key is read from the environment first and from `.env` second; it is
never written to a result, a log or a job file.

## HG002 on every autosome (2026-09-11)

`scripts/twin_genome_wide.py` applied the GIAB v4.2.1 benchmark calls to
all 22 local autosomes and kept only the statistics
(`hg002_twin_by_chromosome.json`): 4,048,342 PASS variants
(3,460,251 SNVs), 1,619,891 applied on haplotype 1 and
4,032,208 on haplotype 2 (unphased heterozygotes go to
haplotype 2 by policy), 0 reference mismatches and
840 skipped overlaps. Zero mismatches over four million
calls is the check that the fetched hg38 sequences and the benchmark agree
base for base. chrX and chrY are not in the benchmark set. No haplotype
FASTA is kept; `genomeos twin build` makes one on demand.

## Your own genome (`genomeos individual import`, 2026-09-11)

A geneticist's own data enters the way HG002 does. One VCF against GRCh38,
plain or gzipped, from any caller:

    genomeos individual import /path/genome.vcf.gz --name ME --note "GATK 4.5, 2026-08, own consent"
    genomeos individual list
    genomeos individual check --name ME               # apply to GRCh38, count reference mismatches: the assembly check
    genomeos individual genes --name ME --chrom chr21
    genomeos individual predict --name ME --gene APP --chrom chr21   # AlphaGenome feature d, needs the key
    genomeos lookup chr17:7675088 C>T          # now says whether ME carries it
    genomeos report TP53 --chrom chr17 ...      # lists ME's variants inside the gene
    genomeos twin build --vcf data/individuals/ME/ME_chr21.vcf --sample ME --chrom chr21 --genome data/reference/chr21.fa

The file is streamed once and split into one small PASS file per chromosome
under `data/individuals/<name>/` with a manifest (source file, sample column,
counts per chromosome, phased genotypes, rows skipped as filtered or on
other contigs). Chromosome names `21`, `chr21`, `MT` are normalised to hg38's;
contigs and decoys are dropped and counted. The directory is git-ignored and
no layer sends a byte of it anywhere: every lookup, dossier and gene walk is
local. Genotypes are measured evidence (the caller's); consequences are
derived by the local trace on the canonical transcript. The Twin tab has the
same import, list and gene-by-gene table, and the predicted effect of the
person's regulatory variants on one gene (docs/ALPHAGENOME.md, feature d).
Per-person results never land under `data/results`: the imported files, the
gene walks and the predictions stay under git-ignored directories, so a
named person's genome cannot be committed by accident.

`individual check` is the first thing to run after an import: it applies the
person's variants to the local GRCh38 sequence, haplotype by haplotype, keeps
the statistics only (no FASTA) and counts the calls whose reference allele
disagrees with the reference base. Up to 0.5% mismatches is a match; a file
against GRCh37/hg19 disagrees at most positions and is called out as such,
before anyone reads a consequence off the wrong coordinates. HG002's own
benchmark gives 0 mismatches over 4 million calls, which is what the check
is calibrated against.

Checked on HG002's chr21 rows re-imported under another name: 55,210 PASS
variants; 194 of 221 coding genes carry a variant, 269 coding SNVs (135
missense, 133 synonymous, 1 nonsense), the KRTAP10 cluster on top as
expected for a highly polymorphic family.

## Optional: peptide/HLA binding predictors

The therapeutic pipeline enumerates the peptides a mutation creates; whether
one of them binds the patient's HLA needs a predictor, and GenomeOS ships
none. Two can be enabled, and the analysis runs without either (binding is
reported as unavailable, never invented).

| Predictor | Licence | How to enable |
|---|---|---|
| MHCflurry 2 | Apache 2.0, usable by anyone | `uv sync --extra hla` then `uv run mhcflurry-downloads fetch models_class1_presentation` |
| NetMHCpan 4.1 | academic licence from DTU Health Tech | request the package at https://services.healthtech.dtu.dk/services/NetMHCpan-4.1/, install it, then put `netMHCpan` on `PATH` or set `NETMHCPAN=/path/to/netMHCpan` |

`genomeos therapeutic --hla-predictor mhcflurry|netmhcpan|none` chooses one;
the default is the first installed. `/api/features` and the Progress tab list
both with their licence and their state. GenomeOS never downloads or
redistributes NetMHCpan; a commercial user needs a licence from DTU.
