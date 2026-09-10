#!/usr/bin/env bash
# Download a real human chromosome (or whole assembly) into data/reference/.
#   scripts/fetch_reference.sh chr21          # hg38 chromosome 21, ~13 MB gz
#   scripts/fetch_reference.sh chrM           # mitochondrial genome, tiny
#   scripts/fetch_reference.sh hs1            # T2T-CHM13 v2.0 complete genome, ~900 MB
#   scripts/fetch_reference.sh gencode        # GENCODE 50 GFF3 annotation, ~60 MB
#   scripts/fetch_reference.sh hg002          # HG002 T2T diploid assembly v1.1 (a real individual), ~1 GB
#   scripts/fetch_reference.sh hg002-vcf      # HG002 GIAB v4.2.1 benchmark variants vs GRCh38
#   scripts/fetch_reference.sh hg001-vcf      # NA12878 / HG001 GIAB benchmark variants vs GRCh38
#   scripts/fetch_reference.sh 1kg-chr21      # 1000 Genomes 3,202 samples, chr21 genotypes
# All of these are open-access. HG001/HG002 are Personal Genome Project
# participants with open consent (usable in research and commercially).
set -euo pipefail
target="${1:-chr21}"
dest="$(dirname "$0")/../data/reference"
mkdir -p "$dest"
case "$target" in
  hs1)     url="https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz" ;;
  hg38)    url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz" ;;
  gencode) url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gff3.gz" ;;
  hg002)   url="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.1.fasta.gz" ;;
  hg002-vcf) url="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" ;;
  hg001-vcf) url="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/NISTv4.2.1/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" ;;
  1kg-chr21) url="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_raw_GT_with_annot/20201028_CCDG_14151_B01_GRM_WGS_2020-08-05_chr21.recalibrated_variants.vcf.gz" ;;
  chr*)    url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/${target}.fa.gz" ;;
  *) echo "unknown target: $target" >&2; exit 1 ;;
esac
echo "fetching $url"
curl -L --fail --progress-bar -o "$dest/$(basename "$url")" "$url"
echo "saved to $dest/$(basename "$url")"
