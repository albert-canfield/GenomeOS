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

The compiled proteome follows the same rule: 290 MB of definitions stay in
data/knowledge (rebuildable), and a 2.3 MB distilled table travels inside the
package as the proteome library (docs/PROTEIN.md).

## The reference cache: one file per chromosome, not two

Rule 3 above says discard the raw input once a summary exists. The reference
sequence is the exception: it is not raw input to be distilled, it *is* the
subject, and every lane reads it. So it had to be kept -- and until now it was
kept twice. `Genome.from_fasta` streams `chrN.fa.gz`; `IndexedGenome` seeks
inside `chrN.fa` through its `.fai`, because plain gzip cannot be seeked. Both
forms were therefore live, and the cache carried 3.1 GB of flat FASTA beside
960 MB of gzip for the same 3.1 Gb of sequence.

Blocked gzip removes the duplication rather than trading it away. A BGZF file
is an ordinary gzip file written as a chain of independent ~64 KiB members;
`gzip.open` reads it unchanged, and a companion `.gzi` index turns any
uncompressed offset into a seek plus one inflate. It is the format BAM has
always used, and the project already decodes it from a remote BAM
(`genome/bam_range.py`), so this is one reader more, not one dependency more.
`genomeos/genome/bgzf.py` writes and seeks it with the standard library alone.

```
data/reference/chrN.fa.gz        blocked gzip, level 9
data/reference/chrN.fa.gz.gzi    block offsets (10 KB per chromosome)
data/reference/chrN.fa.gz.fai    the same .fai as before: offsets into the decompressed bytes
```

`IndexedGenome` takes either form and callers name neither: `reference_fasta(chrom)`
returns whichever is cached, preferring the blocked one, and a freshly fetched
UCSC `.fa.gz` (plain gzip, not blocked) still falls back to being decompressed
once, so nothing breaks before conversion.

**What it costs.** Measured on chr21 (46.7 Mb), warm cache, against the flat
`.fa` it replaces:

| Access pattern | Flat `.fa` | Blocked `.fa.gz` |
|---|---|---|
| whole chromosome in one fetch (the lexicon's pattern) | 0.035 s | 0.067 s |
| 50,000 sorted loci of 200-3,000 bp (the element sweep) | 0.21 s | 0.17 s |
| 5,000 unsorted random loci of 100-2,000 bp | 0.052 s | 0.310 s |
| 2,000 unsorted random loci of 10-50 kb | 0.055 s | 0.263 s |

A sorted sweep is *faster* blocked, because consecutive loci land in the same
cached block and cost no syscall at all. The worst case, unsorted random loci,
costs 52 microseconds more per fetch: a job doing 100,000 scattered fetches
pays five seconds. Nothing in the project is a lane that becomes slower in any
way that can be noticed, and the size is 12.8 MB against 60.3 MB for the pair.

**The proof that nothing was lost.** `scripts/compact_reference.py` writes the
blocked file beside the old two, checks that its decompressed bytes hash to the
SHA-256 of the flat `.fa` *and* to the SHA-256 of the plain `.fa.gz` -- both
files it is about to delete are proved redundant, not assumed to be -- records
the hash and the length in `data/results/reference_bgzf.json`, and only then
removes them. It skips any chromosome a running process has open (`lsof`), so a
scoring run never loses the file underneath it. The manifest is committed, so
`tests/test_bgzf.py` re-proves the conversion on every run long after the
originals are gone; the same file also round-trips synthetic FASTA with
soft-masking, N runs and IUPAC codes, and asserts that a file with one base
changed is *rejected*, so the proof can fail.

What was not done: 2-bit packing (four bases per byte, the UCSC `.2bit` idea)
would be smaller still, around 800 MB for the whole cache, but it cannot
represent soft-masking or any non-ACGT character without side tables. That is a
loss of information, so it is not on the table under the rule above.

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
