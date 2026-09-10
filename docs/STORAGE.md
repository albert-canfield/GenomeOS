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

## Summaries produced so far

| Summary | Raw input it replaces | Size before → after |
|---|---|---|
| `hg002_telomere_stream` | 6 GB FASTQ, streamed (2 M reads sampled) | 0 bytes on disk at any time |
| `clock_GSE41169` | GEO methylation matrix (184 MB) | ~15 KB |
| `clinvar_chr21_agreement` | genome-wide ClinVar VCF (185 MB) | ~2 KB (chr21 subset 22 MB kept for the variant test) |
| `library_members` | Reactome mapping (175 MB) | ~300 KB |
| `hg002_chr21` | HG002 benchmark VCF (149 MB) + haplotype FASTA (91 MB) | ~5 MB chr21-only VCF |

What must stay for the code to run: the chr21 and chrM reference (58 MB), the
GENCODE annotation (153 MB, used by many tests), GO (31 MB + 15 MB), and the
Cell Ontology (3 MB). Everything else can be re-fetched with
`scripts/fetch_reference.sh` if a distiller has to be rerun.
