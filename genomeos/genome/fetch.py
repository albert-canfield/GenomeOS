"""Fetch everything GenomeOS needs to work on one more chromosome.

`genomeos data fetch --chrom chr22` brings a chromosome to the same footing
as chromosome 21: the sequence (UCSC hg38, ~13-80 MB gz), its GENCODE gene
models (the genome-wide file streamed once and only this chromosome's rows
kept, ~2-10 MB gz), the ENCODE regulatory elements (streamed, rows kept) and
the RepeatMasker annotation (UCSC API, distilled). After that every view
works for it: Blocks, Flow, regulation, domains, UNKNOWN, lookup, rna.
Nothing else is written; `genomeos data status` shows the footprint.
"""

from __future__ import annotations

import gzip
import io
import shutil
import urllib.request
from pathlib import Path
from typing import Any

UCSC_FASTA = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz"
GENCODE = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gff3.gz"
)
REFERENCE = Path("data/reference")
RESULTS = Path("data/results")


def gencode_chrom_path(chrom: str) -> Path:
    return RESULTS / f"gencode_v50_{chrom}.gff3.gz"


def _download(url: str, dest: Path, progress=None) -> None:
    req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (fetch)"})
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    with urllib.request.urlopen(req, timeout=600) as r, open(tmp, "wb") as fh:  # noqa: S310
        shutil.copyfileobj(r, fh, 1 << 20)
    tmp.rename(dest)
    if progress:
        progress(f"{dest.name}: {dest.stat().st_size / 1e6:.1f} MB")


def fetch_sequence(chrom: str, progress=None) -> Path:
    dest = REFERENCE / f"{chrom}.fa.gz"
    if not dest.exists():
        _download(UCSC_FASTA.format(chrom=chrom), dest, progress)
    return dest


def fetch_gencode_chrom(chrom: str, progress=None) -> Path:
    """Stream the genome-wide GENCODE GFF3 and keep only this chromosome's rows (gzip)."""
    dest = gencode_chrom_path(chrom)
    if dest.exists():
        return dest
    req = urllib.request.Request(GENCODE, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
    dest.parent.mkdir(parents=True, exist_ok=True)
    kept = 0
    with (
        urllib.request.urlopen(req, timeout=600) as resp,  # noqa: S310
        gzip.open(io.BufferedReader(resp, 1 << 20), "rt") as fh,
        gzip.open(dest, "wt") as out,
    ):
        for line in fh:
            if line.startswith("#"):
                if line.startswith("##gff-version"):
                    out.write(line)
                continue
            if line.split("\t", 1)[0] == chrom:
                out.write(line)
                kept += 1
    if progress:
        progress(f"GENCODE rows kept for {chrom}: {kept:,} ({dest.stat().st_size / 1e6:.1f} MB)")
    return dest


def fetch_chromosome(
    chrom: str, progress=None, elements: bool = True, repeats: bool = True
) -> dict[str, Any]:
    out: dict[str, Any] = {"chrom": chrom}
    out["sequence"] = str(fetch_sequence(chrom, progress))
    out["gencode"] = str(fetch_gencode_chrom(chrom, progress))
    if elements:
        from genomeos.genome.regulatory import load_ccres, save_ccres, stream_ccres

        if not load_ccres(chrom):
            els = stream_ccres({chrom}, progress=progress)
            save_ccres(chrom, els)
            out["ccres"] = len(els)
    if repeats:
        from genomeos.genome.repeats import fetch_repeats, load_repeats, save_repeats, summarise
        from genomeos.results import save_result

        reps = load_repeats(chrom)
        if not reps:
            reps = fetch_repeats(chrom)
            save_repeats(chrom, reps)
            save_result(f"rmsk_{chrom}", summarise(chrom, reps))
        out["repeats"] = len(reps)
    return out


GIAB_HG002 = (
    "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/"
    "NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
)


def individual_vcf_path(chrom: str, sample: str = "HG002") -> Path:
    """The per-chromosome variant file of the test human: chr21 is the committed distilled subset,
    the rest are split once from the 156 MB GIAB file into data/reference (not committed)."""
    committed = RESULTS / f"{sample}_{chrom}.vcf"
    return committed if committed.exists() else REFERENCE / f"{sample}_{chrom}.vcf"


def fetch_individual(sample: str = "HG002", progress=None, url: str = GIAB_HG002) -> dict[str, int]:
    """Stream the GIAB benchmark VCF once and keep every chromosome's PASS rows as a small file
    (about 60 MB in total for 22 autosomes), so twins and lookups work on any fetched chromosome."""
    header = (
        "##fileformat=VCFv4.2\n##source=GIAB {sample} v4.2.1 benchmark, {chrom} PASS subset via GenomeOS\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample}\n"
    )
    REFERENCE.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
    counts: dict[str, int] = {}
    handles: dict[str, Any] = {}
    try:
        with (
            urllib.request.urlopen(req, timeout=900) as resp,  # noqa: S310
            gzip.open(io.BufferedReader(resp, 1 << 20), "rt") as fh,
        ):
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 10 or f[6] not in ("PASS", "."):
                    continue
                chrom = f[0]
                if chrom not in handles:
                    if (RESULTS / f"{sample}_{chrom}.vcf").exists():
                        handles[chrom] = None  # already distilled and committed
                    else:
                        handles[chrom] = open(REFERENCE / f"{sample}_{chrom}.vcf", "w")  # noqa: SIM115
                        handles[chrom].write(header.format(sample=sample, chrom=chrom))
                        if progress:
                            progress(f"{chrom}: splitting {sample} variants")
                h = handles[chrom]
                if h is None:
                    continue
                gt = f[9].split(":")[0]
                h.write(f"{chrom}\t{f[1]}\t{f[2]}\t{f[3]}\t{f[4]}\t.\tPASS\t.\tGT\t{gt}\n")
                counts[chrom] = counts.get(chrom, 0) + 1
    finally:
        for h in handles.values():
            if h is not None:
                h.close()
    return counts
