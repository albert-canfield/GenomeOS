"""Whole-genome anatomy by stream, distil, discard.

For each human chromosome: download the chromosome FASTA and the rows of the
genome-wide GENCODE annotation for it, count the blocks, save the inventory,
delete the files. Peak disk use is one chromosome (<250 MB); the output is a
single JSON with every chromosome's anatomy.

    uv run python scripts/anatomy_genome_wide.py [--chroms chr1 chr2 ...]
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import tempfile
import time
import urllib.request
from pathlib import Path

from genomeos.genome import Annotation, Genome, anatomy_of, design_lessons
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
    ap.add_argument("--result", default="anatomy_hg38_by_chromosome")
    args = ap.parse_args()
    existing = load_result(args.result) or {"chromosomes": {}}
    done = existing["chromosomes"]
    tmp = Path(tempfile.mkdtemp(prefix="genomeos-anatomy-"))
    try:
        gff = tmp / "gencode.gff3.gz"
        print("fetching genome-wide GENCODE annotation (streamed once, split per chromosome)…", flush=True)
        fetch(GENCODE, gff)
        # split annotation per chromosome so each is small
        per: dict[str, list[str]] = {}
        with gzip.open(gff, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                per.setdefault(line.split("\t", 1)[0], []).append(line)
        gff.unlink()
        for chrom in args.chroms:
            if chrom in done:
                print(f"{chrom}: already done", flush=True)
                continue
            t0 = time.time()
            fa = tmp / f"{chrom}.fa.gz"
            fetch(UCSC.format(chrom=chrom), fa)
            gpath = tmp / f"{chrom}.gff3"
            gpath.write_text("".join(per.get(chrom, [])))
            seq = Genome.from_fasta(fa).chromosomes[chrom].sequence
            ann = Annotation.from_gff3(gpath, {chrom})
            a = anatomy_of(f"human {chrom}", seq, ann, chrom)
            d = a.to_dict()
            d["lessons"] = design_lessons([a])
            d["seconds"] = round(time.time() - t0)
            done[chrom] = d
            fa.unlink()
            gpath.unlink()
            existing["chromosomes"] = done
            existing["assembly"] = "GRCh38 (UCSC hg38) + GENCODE 50"
            existing["disk_policy"] = (
                "each chromosome downloaded, counted and deleted; nothing kept but this summary"
            )
            save_result(args.result, existing)
            f = d["composition_fraction"]
            print(
                f"{chrom}: {d['length'] / 1e6:.1f} Mb, {d['genes']['protein_coding']} coding genes, "
                f"CDS {f['cds']:.2%}, intron {f['intron']:.1%}, {d['seconds']} s",
                flush=True,
            )
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    print("done:", len(done), "chromosomes")


if __name__ == "__main__":
    main()
