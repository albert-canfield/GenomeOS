"""Investigate the UNKNOWN blocks of every human chromosome: stream, distil, discard.

Once: stream the ENCODE cCRE registry (64 MB) and keep a small per-chromosome
subset; download the genome-wide GENCODE annotation and split it per
chromosome in memory. Then per chromosome: download the FASTA, run
`investigate`, save `unknown_<chrom>.json`, delete the FASTA. Peak disk use is
one chromosome.

    uv run python scripts/unknown_genome_wide.py [--chroms chr1 chr2 ...]
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import tempfile
import time
import urllib.request
from pathlib import Path

from genomeos.genome import Annotation, Genome
from genomeos.genome.regulatory import RESULTS, save_ccres, stream_ccres
from genomeos.genome.unknown import investigate
from genomeos.jobs import heartbeat
from genomeos.results import load_result, save_result

UCSC = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/{chrom}.fa.gz"
GENCODE = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_50/gencode.v50.annotation.gff3.gz"
)
ALL = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM"]


def fetch(url: str, dest: Path) -> None:
    with urllib.request.urlopen(url, timeout=600) as r, open(dest, "wb") as fh:  # noqa: S310
        shutil.copyfileobj(r, fh, 1 << 20)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", nargs="*", default=ALL)
    args = ap.parse_args()
    # a chromosome is done only if it was classified with the curated repeat annotation
    todo = [c for c in args.chroms if not (load_result(f"unknown_{c}") or {}).get("curated_repeats")]
    if not todo:
        print("all chromosomes already investigated")
        return
    missing_ccre = [c for c in todo if not (RESULTS / f"ccres_{c}.bed.gz").exists()]
    if missing_ccre:
        print(f"streaming ENCODE cCREs for {len(missing_ccre)} chromosomes…", flush=True)
        by: dict[str, list] = {}
        for c in stream_ccres(set(missing_ccre)):
            by.setdefault(c.chrom, []).append(c)
        for chrom, els in by.items():
            save_ccres(chrom, els)
    tmp = Path(tempfile.mkdtemp(prefix="genomeos-unknown-"))
    try:
        gff = tmp / "gencode.gff3.gz"
        print("fetching genome-wide GENCODE annotation (split per chromosome, then deleted)…", flush=True)
        fetch(GENCODE, gff)
        per: dict[str, list[str]] = {}
        with gzip.open(gff, "rt") as fh:
            for line in fh:
                if not line.startswith("#"):
                    per.setdefault(line.split("\t", 1)[0], []).append(line)
        gff.unlink()
        summary = load_result("unknown_genome_wide") or {"chromosomes": {}}
        for chrom in todo:
            t0 = time.time()
            fa = tmp / f"{chrom}.fa.gz"
            fetch(UCSC.format(chrom=chrom), fa)
            gpath = tmp / f"{chrom}.gff3"
            gpath.write_text("".join(per.get(chrom, [])))
            # curated repeats for this chromosome: fetched, distilled to a summary, BED deleted after use
            from genomeos.genome.repeats import RESULTS as RMSK_DIR
            from genomeos.genome.repeats import fetch_repeats, save_repeats
            from genomeos.genome.repeats import summarise as summarise_repeats

            try:
                reps = fetch_repeats(chrom)
                save_repeats(chrom, reps)
                save_result(f"rmsk_{chrom}", summarise_repeats(chrom, reps))
            except Exception as e:  # noqa: BLE001
                print(f"  {chrom}: RepeatMasker unavailable ({str(e)[:60]}); sequence rules only", flush=True)
            seq = Genome.from_fasta(fa).chromosomes[chrom].sequence
            ann = Annotation.from_gff3(gpath, {chrom})
            r = investigate(seq, ann, chrom, progress=lambda m, c=chrom: print(f"  {c}: {m}", flush=True))
            r["seconds"] = round(time.time() - t0)
            save_result(f"unknown_{chrom}", r)
            summary["chromosomes"][chrom] = {
                "unknown_blocks": r["unknown_blocks"],
                "unknown_bp": r["unknown_bp"],
                "classified_fraction": r["classified_fraction"],
                "curated_repeats": r.get("curated_repeats", False),
                "by_class": r["by_class"],
                "seconds": r["seconds"],
            }
            save_result("unknown_genome_wide", summary)
            heartbeat("unknown_genome_wide")
            fa.unlink()
            gpath.unlink()
            if chrom != "chr21":
                (RMSK_DIR / f"rmsk_{chrom}.bed.gz").unlink(missing_ok=True)
            mb = r["unknown_bp"] / 1e6
            print(
                f"{chrom}: {mb:.1f} Mb unknown, {r['classified_fraction']:.0%} classified, {r['seconds']} s",
                flush=True,
            )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    # never give up: a failed pass (a broken download, a source outage) is retried with a growing pause
    pause = 60
    for _attempt in range(1000):
        try:
            main()
            break
        except Exception as e:  # noqa: BLE001
            print(f"pass failed ({str(e)[:100]}); next pass in {pause} s", flush=True)
            time.sleep(pause)
            pause = min(pause * 2, 1800)
