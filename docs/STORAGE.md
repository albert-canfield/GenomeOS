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

## Against the GitHub limits (2026-09-27)

The daily merge of `dev` into `main` starts now, so the question is no longer how much disk the
project uses but how much of it GitHub is asked to keep. Those are different numbers by a factor of
146: **9.85 GiB on disk, 69.2 MiB in a clone.** Every limit below was read from `docs.github.com`
on this day, and each measurement names the command that produced it.

### What the repository actually is

| what | size | command |
|---|---|---|
| whole working tree | 9.85 GiB | `du -sk .` |
| `data/` | 8.53 GiB | `du -sk data` |
| `.venv/` (ignored) | 1.15 GiB | `du -sk .venv` |
| `.git/` | 147 MiB | `du -sk .git` |
| `genomeos/` | 14.5 MiB on disk, **7.2 MiB tracked** | `du -sk genomeos`; `git ls-files genomeos \| xargs stat -f %z` |
| `tests/` | 8.1 MiB on disk, **1.2 MiB tracked** | same, the difference is `__pycache__` |
| `scripts/` | 2.2 MiB on disk, **1.1 MiB tracked** | same |
| `docs/` | 2.0 MiB, all tracked | same |

**Of the 8.53 GiB in `data/`, 256.5 MiB is tracked and 8.28 GiB is git-ignored — 2.9 per cent
committed, 97.1 per cent local.** That ratio is the stream-distil-discard rule, measured. The whole
tracked tree is 1,365 files and 268.4 MiB (`git ls-tree -r -l HEAD`), of which `data/` is 761 files.

### The history, and a finding that has nothing to do with GitHub

`git count-objects -vH` reports **8,029 loose objects, 145.84 MiB, and `packs: 0`**. This repository
has never been packed. Across all 728 commits there are 3,600 distinct blobs totalling **492.4 MiB
uncompressed**, held as 113.9 MiB of individually zlib'd loose files.

Packed, the same history is far smaller. Measured without touching the repository, by piping every
reachable object through `git pack-objects --stdout` into `wc -c`:

```
git rev-list --all --objects | awk '{print $1}' | git pack-objects --stdout | wc -c
72564509
```

**69.2 MiB is what a full clone transfers**, and it is the number every GitHub repository limit
should be compared against. `git gc` was **not** run: three other lanes hold this checkout and one of
them moved the `dev` tip from `36f36d6` to `a1f8418` while the measurement was in progress, which is
exactly the "if in doubt, do not" case. It is worth about 77 MiB of local `.git` when the lanes are idle.

### The largest files, and the largest that ever existed

The biggest tracked file at `HEAD` is `data/results/origin_genome_wide.json` at **6.10 MiB**
(`git ls-tree -r -l HEAD | sort -nr`); then `eqtl_targets.json` 5.87, `human_panel_chr1.json` 5.53,
`human_panel_chr2.json` 4.38, `human_panel_chr11.json` 3.38, and a further fifteen between 2.49 and
3.38 MiB, every one of them a `human_panel_chr*`, `motifs_chr*`, `panel_union_five_arms.json` or
`noncoding_chr21.bio`.

The largest blob that ever existed is a **superseded** version of the same file at **6.12 MiB**
(`git rev-list --objects --all | git cat-file --batch-check`). Three versions of `eqtl_targets.json`
(5.87, 3.75, 3.08 MiB) and three of `noncoding_chr21.bio` survive in history; nothing was ever
committed and deleted that is larger than what is tracked now. **The high-water mark of this
repository's per-file size, over its whole history, is 6.12 MiB.**

### What compresses, and the thing that does not

| extension | tracked bytes | in the pack | ratio |
|---|---|---|---|
| `.json` (654 files) | 234.4 MiB | **29.8 MiB** | 12.7 % |
| `.gz` (49 files) | 18.6 MiB | **18.5 MiB** | 99.4 % |
| `.py` (532 files) | 6.3 MiB | 2.3 MiB | 36 % |
| `.bio`, `.md`, `.vcf` | 8.0 MiB | 1.7 MiB | 21 % |

**The 49 already-gzipped files are 6.9 per cent of the tracked bytes and 35 per cent of the clone.**
Pretty-printing, by contrast, is nearly free to the repository: compacting
`origin_genome_wide.json` cuts it from 6.10 to 3.57 MiB on disk (41 per cent) but from 0.37 to
0.32 MiB gzipped (14 per cent). Compact JSON would be a token decision, not a storage one.

### Growth

Tracked snapshot, one commit per day (`git rev-list -1 --before=<day> dev`, then `git ls-tree -r -l`):

```
2026-09-12  293 commits  882 files  162.0 MiB
2026-09-14  412 commits 1080 files  214.6 MiB
2026-09-17  610 commits 1271 files  257.4 MiB
2026-09-22  728 commits 1365 files  268.4 MiB
2026-09-27  728 commits 1365 files  268.4 MiB   (idle since 09-22)
```

Clone size, packed at three tips: **45.3 MiB on 09-13, 63.3 MiB on 09-17, 69.2 MiB today.** Over the
fourteen days that is **+1.71 MiB per day**; over the busiest four it was **+4.51 MiB per day**. Eleven
of the last sixteen days carried commits, between 5 and 104 of them.

| target | at 1.71 MiB/day | at the peak 4.51 MiB/day |
|---|---|---|
| 1 GiB (GitHub's "ideally") | 558 days, 2028-04-07 | 212 days, 2027-04-26 |
| 5 GiB ("strongly recommended") | 2,953 days, 2034-10-28 | 1,121 days, 2029-10-22 |
| 10 GiB (on-disk recommendation) | 5,947 days | 2,258 days |

**The daily merge itself costs nothing.** `dev` is four commits ahead of `main` and the objects unique
to it are 23, totalling **0.55 MiB** (`git rev-list --objects dev --not origin/main`). A merge commit
introduces no blobs; the data was already pushed on `dev`. The heaviest single day measured, 09-17
with 95 commits, added **12.01 MiB** of new compressed objects — that is the largest push this
project has ever needed.

### The limits, with headroom

Verified from `docs.github.com` on 2026-09-27. Nothing below is remembered.

| limit | GitHub's figure | kind | ours | headroom |
|---|---|---|---|---|
| per-file push block | **100 MiB** | enforced, outright reject | 6.12 MiB ever | 16.3× |
| per-file warning | **50 MiB** | warning from Git | 6.12 MiB ever | 8.2× |
| browser upload | 25 MiB | enforced | n/a, all commits from the CLI | — |
| repository, recommended | **"ideally less than 1 GB"** | recommendation | 69.2 MiB | 14.8×, 558 days |
| repository, soft | **"less than 5 GB is strongly recommended"** | recommendation | 69.2 MiB | 74× |
| repository on-disk | **10 GB** | recommendation | 69.2 MiB | 148× |
| per push | **2 GiB** | **enforced** | 12.01 MiB worst day | 170× |
| push rate | 6 per minute | enforced | one per commit cycle | — |
| directory width | **3,000 entries** | recommendation | **691 in `data/results`** | **165 days** |
| directory depth | 50 | recommendation | 3 | — |
| branches | 5,000 | recommendation | 2 | — |
| files per repository | **none stated** | — | 1,365 | — |
| files per commit | **none stated** | — | 106 worst | — |
| files per rendered diff | **300** (25 renderable) | enforced on display | **106 in one commit** | — |
| lines per rendered diff | 20,000 lines or 1 MB | enforced on display | exceeded by one day's work | **none** |

**What GitHub does when a repository is too large:** nothing automatic. The wording is "If your
repository excessively impacts our infrastructure, you might receive an email from GitHub Support
asking you to take corrective action." There is no size at which pushes start failing, other than the
2 GiB per push, which is a property of one pack and not of the repository.

**There is no stated limit on the number of files in a repository or in a single commit.** The only
file counts GitHub publishes are display caps on diffs and a 3,000-entry-per-directory recommendation
made for the sake of Git performance, not storage.

### Git LFS would make this project worse

The quotas, verified: **10 GiB storage and 10 GiB bandwidth free** on GitHub Free, Pro and Free for
organizations; 250 GiB of each on Team and Enterprise Cloud. Per-file ceiling 2 GB on Free and Pro.
Bandwidth is "billed for each GiB of data downloaded", and every clone and every CI checkout counts
against the repository owner's allowance.

It is not needed, and it would not help:

1. **Nothing is near the block.** The largest object this project has ever committed is 6.12 MiB,
   sixteen times below the 100 MiB reject and eight times below the 50 MiB warning.
2. **LFS gives up the compression this project lives on.** LFS objects are stored as uploaded — no
   zlib, no deltas between versions. These are text JSON files that compress to 12.7 per cent. The
   492.4 MiB of blobs that git holds in a 69.2 MiB pack would become roughly 492 MiB of LFS storage,
   and each new version of `human_panel_chr1.json` would be stored whole instead of as a delta.
3. **And it would put the transfer on a meter.** A checkout would pull about 268 MiB of LFS objects
   where git transfers 69.2 MiB — about **38 clones or CI runs to exhaust the free monthly
   bandwidth**, against plain git transfer which is unmetered.

A 3.9× larger download, on a budget, to solve a problem that is sixteen times away. No.

### The recommendation

**Safe to keep committing, as it is, for years.** At 1.71 MiB per day the 1 GB recommendation is 558
days out and the 5 GB one is eight years out; at the fastest rate the project has ever sustained,
seven months and three years. No enforced limit is within reach of any projection: the per-file block
is 16× away and the per-push limit 170× away.

**The first thing that will break is not a size limit at all — it is the daily merge's diff.** GitHub
renders at most 300 files and 20,000 lines in one diff, and a single commit here has already touched
106 files while a working day runs to 104 commits. The pull request will merge correctly and refuse to
show itself. **Second, in about 165 days, `data/results` crosses the 3,000-entry directory-width
recommendation**: 691 tracked entries today, growing 14 a day over the last fortnight and 39 a day
while the sweeps ran. That one is real but gradual, and the fix is a shallow split of the directory,
which is a path change across modules and belongs in a plan of its own.

**One thing is tracked that arguably should not be.** `data/results/ccres_chr*.bed.gz` — 25 files
besides chr21 — is **11.46 MiB, and because it is already gzipped it is 16.6 per cent of the clone**,
the largest single saving available. These are streamed ENCODE SCREEN track rows, re-fetchable in
minutes, and the project's own `.gitignore` already keeps the analogous `rmsk`, `superdups` and
`dnase` tracks to chr21 as the worked example. The tests do not need them: they write their own
gzipped fixtures into `tmp_path`. But the all-enhancer scoring chain reads them chromosome by
chromosome, so this is a recommendation with an owner and a plan, not a deletion.

**`.gitignore` is doing its job, and it was verified rather than assumed.** The five largest ignored
paths are `data/knowledge` 6.31 GiB, `data/reference` 1.41 GiB, `.venv` 1.15 GiB,
`data/individuals` 321 MiB and `data/twins/HG002_chr22.fa` 98.5 MiB, and every one is covered by a
rule that states its reason — fetched datasets, a virtual environment, imported personal genomes,
a named person's biological state. The sharper test: `git ls-files | git check-ignore --no-index
--stdin -v` finds **30 tracked files that an ignore pattern would otherwise exclude, and all 30 are
re-admitted by an explicit `!` negation** — the chr21 `dnase_*` and `reader_*` set, `rmsk_chr21`,
`superdups_chr21`, `gencode_v50_chr21_chrM`, `clinvar_chr21_coding.vcf.gz`. Nothing is tracked by
accident.

*Correction, same day: `data/individuals` is **329 MiB**, not 321 MiB as written above — `du -sk` reports 336,864 KiB. No other figure depended on it.*

**The one change with the best size-to-effort ratio is `git gc --prune=now`, run when the other lanes
are idle.** It is one command, it is not needed for any GitHub limit, and it takes `.git` from 147 MiB
to about 70 MiB — but the reason to run it is the 8,029 loose objects with no pack, which every
`git status` and every `git add` in this shared checkout has to stat one at a time.

**Can someone on an ordinary connection clone this repository comfortably today? Yes.** 69.2 MiB
transferred, about 29 seconds on a 20 Mbit line, checking out to 268 MiB — smaller than a great many
ordinary application repositories. They will not get `data/knowledge` or `data/reference`, and they
are not meant to: that is the 8.28 GiB the rule at the top of this document exists to keep out.
