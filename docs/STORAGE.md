# Storage economics: stream, distil, discard

Real biology data is large (a 30x genome is ~100 GB of reads) and this project
runs on an ordinary laptop. The rule is therefore:

1. **Stream when possible.** A test that needs reads pulls them over HTTP into
   the estimator and never writes them to disk (`genomeos telomere <url>`).
2. **Distil every real-data run into a small summary** under `data/results/`:
   what was measured, from which input, with which method, and the numbers.
   Summaries are committed to git and are the durable outcome.
3. **Discard the raw input** once its summary exists (`genomeos data clean --yes`).
   Tests read the summary when the raw file is gone, so nothing is lost.
4. **Embed what was learned.** Calibrations and rules extracted from summaries
   move into the engine (parameters with evidence, library membership tables),
   so the software carries the knowledge, not the data.

```
genomeos data status          disk use, largest files, what can be distilled
genomeos data distil          run every distiller whose inputs are present
genomeos data clean           dry run: what would be deleted and how much freed
genomeos data clean --yes     delete raw inputs that have summaries
genomeos data results         list summaries
```

## Working on another chromosome

`genomeos data fetch --chrom chr22` brings one more chromosome to the same
footing as chromosome 21, within the same economics: the sequence
(UCSC, 13-80 MB gz, kept under data/reference), its GENCODE rows (the
genome-wide 60 MB file is streamed once and only this chromosome's rows
are written, 2-10 MB), the ENCODE elements (streamed, rows kept) and the
RepeatMasker annotation (20 MB fetched, 1 MB kept). After that Blocks, Flow,
regulation, domains, unknown, lookup, rna and proteome all work for it;
`genomeos data status` shows the footprint and `data clean` removes what
has a summary.

`genomeos data fetch --individual --chrom chr22` adds the test human: the
156 MB GIAB HG002 benchmark is streamed once and every chromosome's PASS
rows are kept as a file under data/reference (158 MB for the 22
autosomes, not committed; chr21's subset stays the committed distilled
copy). After that `genomeos twin build` and `genomeos lookup` work on any
fetched chromosome; the Progress tab has it as the `fetch_hg002` job.

## Summaries produced so far

| Summary | Raw input it replaces | Size before → after |
|---|---|---|
| `hg002_telomere_stream` | 6 GB FASTQ, streamed (2 M reads sampled) | 0 bytes on disk at any time |
| `clock_GSE41169` | GEO methylation matrix (184 MB) | ~15 KB |
| `clinvar_chr21_agreement` | genome-wide ClinVar VCF (185 MB) | ~2 KB (chr21 subset 22 MB kept for the variant test) |
| `library_members` | Reactome mapping (175 MB) | ~300 KB |
| `hg002_chr21` | HG002 benchmark VCF (149 MB) + haplotype FASTA (91 MB) | ~2 MB chr21-only VCF |
| `gencode_chr21_chrM` | genome-wide GENCODE 50 (153 MB) | 1.9 MB subset |
| `clinvar_chr21_coding` | ClinVar chr21 (22 MB subset) | 1.1 MB coding records |
| `hg002_methylation_stream` | two nanopore bedMethyl files (1.2 GB), streamed | 0 bytes on disk; ~30 KB of betas at 418 clock CpGs |
| `clock_probes_hg38.json` (packaged) | Illumina 450k manifest (29 MB) | 13 KB |

What must stay for the code to run: the chr21 and chrM reference (58 MB), the
GENCODE annotation (153 MB, used by many tests), GO (31 MB + 15 MB), and the
Cell Ontology (3 MB). Everything else can be re-fetched with
`scripts/fetch_reference.sh` if a distiller has to be rerun.
