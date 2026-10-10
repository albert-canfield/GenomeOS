A demo tumour above data level 1: one public tumour with its matched normal genome.

The demo tumour `data/demo/cancer_tumour.vcf` is synthetic and runs at data level 1. This
pair runs at data level 2. Level 2 means the matched normal's own variant calls: the profile
holds them and the pipeline subtracts them. A somatic call set made against a normal elsewhere
is level-1 data. TCGA and CPTAC keep the normal's calls under controlled access, so neither can
reach level 2 from open data. docs/THERAPEUTICS.md, "A demo case above data level 1", records
that decision and the choice of this case. Both were committed before these files were built.

| file | what it is |
|---|---|
| `HG008-T.vcf` | NYGC's high-confidence somatic calls in HG008-T, the tumour column only, with AD, DP and an allele fraction from AD |
| `HG008-N-D.vcf` | every HaplotypeCaller call in the matched normal HG008-N-D, common variants included, with GT only |
| `provenance.json` | sources, NIST's md5s (checked), the region rule and its size, and what each rewrite dropped |
| `vep_cache.jsonl` | the Ensembl VEP annotations of the tumour records, so the demo also runs offline |
| `level_comparison.json` | the same tumour file at level 2 (with the normal) and at level 1 (without it) |
| `therapeutics/` | the level-2 analysis: the report and the five dataset layers |

**The case.** GIAB HG008 is a pancreatic ductal adenocarcinoma cell line, HG008-T, from the Liss
lab's batch 0823p23 (passage 23), with the same patient's normal duodenal tissue, HG008-N-D.
NIST says the patient "explicitly consented for fully public dissemination of genomic data and
cell lines". The tumour is a cell line, not a biopsy.

**The source.** The NYGC Somatic Pipeline v7 analysis of 2024-04-12, on GRCh38-GIABv3:
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data_somatic/HG008/Liss_lab/analysis/NYGC-somatic-pipeline_20240412/

**What was changed.**

- **The region.** Both files are cut to the same region set: the union of every GENCODE v50 CDS
  feature, padded by 2 bp on each side. That is 234,760 intervals and 39,601,470 bases. A record
  is kept when any base of its reference allele falls inside. This keeps 107 of the tumour's
  11,142 calls and 33,714 of the normal's 5,023,572.
- **The tumour.** INFO is dropped. The source records carry no GT, so none is written.
- **The normal.** It is reduced to GT. FILTER is kept as written.

The data are the National Institute of Standards and Technology's, under the NIST data use
policy: not subject to copyright in the United States, with a royalty-free worldwide grant
elsewhere. GenomeOS changed them on 2026-10-10 as described above.

**What level 2 does not give.** Tumour RNA (level 3) is not supplied, so expression in this
tumour, and with it surface expression, is reported as unavailable. HG008's RNA-seq exists only
as raw reads, from NIST passages 21 and 103 rather than this batch. No HLA is supplied either,
so no peptide is prioritised.

Rebuild and re-run:

```bash
uv run python scripts/demo_hg008_level2.py build
uv run python scripts/demo_hg008_level2.py run
genomeos therapeutic --tumour data/demo/hg008_pdac/HG008-T.vcf --normal data/demo/hg008_pdac/HG008-N-D.vcf
```

`build` streams the 1.5 GB normal file once, hashing it as it arrives, and stores none of it
outside this directory.
