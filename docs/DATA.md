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
