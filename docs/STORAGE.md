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

**What it costs.** Measured on chr2 (242 Mb), the same bytes read both ways on
the same machine, warm cache:

| Access pattern | Flat `.fa` | Blocked `.fa.gz` | Per fetch |
|---|---|---|---|
| whole chromosome in one fetch (the lexicon's pattern) | 0.41 s | 0.47 s | +60 ms, once |
| 200,000 sorted loci of 200-3,000 bp (the element sweep) | 0.53 s | 0.82 s | +1.5 µs |
| 20,000 unsorted random loci of 100-2,000 bp | 0.28 s | 1.56 s | +64 µs |
| 247 MB on disk | | 80 MB | |

The whole-chromosome read, which is what the heavy lanes actually do, costs 15%
more and is still under half a second for the largest chromosome. A sorted
sweep costs 1.5 microseconds more per locus, because consecutive loci land in
the same cached block and cost no syscall at all; on chr21 the sorted sweep was
in fact *faster* blocked than flat. The worst case is scattered random access,
at 64 microseconds more per fetch: a job doing 100,000 unsorted fetches pays
six seconds. No lane in the project becomes slower in a way anyone can notice.

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

The split HG002 variant files went the same way. `iter_vcf` already took either
form, but several readers opened the file directly as text; they now go through
`individuals.open_variants`, and the conversion proves equality twice over --
the decompressed bytes hash to the original, and both forms are parsed through
`iter_vcf` and compared variant by variant in order. 141 MB to 22 MB, 3,993,132
variants. `data/individuals/` is never touched: a person's imported genome is
not cache and is not compressed, pruned or moved.

**What it came to.** `data/reference` fell from 4.4 GB to 1.4 GB, and the only
flat `.fa` left is the chromosome that was being scored while this ran.

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

## The audit at the end of the sweeps (2026-09-16)

Both genome-wide sweeps finished on this day — every chromosome scored by deletion, every chromosome
read by the human panel — so this is the high-water mark of what the project holds, and the moment to
ask whether any of it can be distilled further. **`genomeos data clean` frees 0 B**: every raw input
with a summary has already been discarded, and what is left is either irreplaceable or a working cache.

**7.9 GB in `data/`, of which 357 MB is committed** (`data/results`, the summaries) and the rest is
git-ignored and local.

| what | size | cost to rebuild |
|---|---|---|
| `knowledge/epigenome` | 2.5 GB | hours of ENCODE streaming; 1,320 float32 arrays, already gzipped |
| `reference` | 1.4 GB | minutes from UCSC; mostly BGZF since the conversion, and read constantly |
| `knowledge/alphagenome` | 1.3 GB | **778,780 model requests** — the sweep itself, bounded by a daily quota |
| `knowledge/human_panel` | 1.1 GB | **267 GB of Cactus alignment streamed**, about 13 hours |
| `results` | 357 MB | the durable outcome; committed, and the only part that survives a clean checkout |
| `individuals` + `twins` | 427 MB | derived from the GIAB HG002 benchmark, minutes to hours |

**Why nothing more is compressed.** The two large caches that look compressible are not free to touch.
The epigenome signal is already `float32` inside gzip, at about 4 MB per cell type, mark and
chromosome. The per-element deletion tables (594 MB of plain JSON under `alphagenome/all_elements`)
would gzip to roughly a fifth, but three modules read them with `read_text()` — `attribution/
element_types.py`, `attribution/crispri.py` and the chain's own reader — so transparent `.gz` support
belongs in one shared loader, changed with their owners rather than under them. It is worth about
480 MB and it is not urgent at 7.9 GB.

**The disk pressure was never this repository.** Four separate runs were stopped by the disk floor on
2026-09-15 and two more on 2026-09-16, on a 460 GB disk holding a 7.9 GB project. The largest single
consumer measured on this machine was **80 GB of stale Chrome code-sign clones** in the system temp
tree (`/private/var/folders/.../X/com.google.Chrome.code_sign_clone`, 58 of them dated weeks earlier),
which macOS recreates on launch and never collects. That is outside the repository and outside what
this project should delete on its own; it is named here because six runs died of it and the next
session to lose a chromosome to free space should look there first, not at `data/`.

**The rule that earned its place.** Streaming and discarding is why 267 GB of alignment costs 1.1 GB
and 778,780 model requests cost 1.3 GB. The one place it was not followed — keeping a per-element
table rather than a summary — is also the one place a lane later had to pack caches by hand to free
3.16 GB, losslessly and with SHA-256 proofs.
