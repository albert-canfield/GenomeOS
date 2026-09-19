# SPDX-License-Identifier: AGPL-3.0-or-later
"""Fold finished chromosomes' per-element deletion answers into one archive each.

    uv run python scripts/pack_element_cache.py                 # every chromosome whose job is complete
    uv run python scripts/pack_element_cache.py --chrom chr22    # one chromosome
    uv run python scripts/pack_element_cache.py --keep-files     # write the archive, remove nothing

The cache holds one small JSON per element (about 14 kB each, 20,000 per chromosome, heading for 13 GB and
a million files genome-wide). A chromosome that is finished never needs its files touched individually
again: packed, it is one gzipped file about eight times smaller, and `load_cached` reads it transparently.
Only chromosomes whose all-elements job says `complete` are packed by default, and the files are removed
only after every one of them is in the archive.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from genomeos.predict.enhancer_target import CACHE, pack

RESULTS = Path("data/results")


def complete_chromosomes(results: Path = RESULTS) -> list[str]:
    out = []
    for f in sorted(results.glob("enhancer_targets_all_chr*.json")):
        try:
            d = json.loads(f.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if d.get("complete"):
            out.append(d["chrom"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", nargs="*", help="chromosomes to pack (default: those with a complete job)")
    ap.add_argument(
        "--keep-files", action="store_true", help="write the archive but keep the per-element files"
    )
    ap.add_argument("--cache", default=str(CACHE))
    args = ap.parse_args()
    chroms = args.chrom or complete_chromosomes()
    if not chroms:
        print("nothing to pack: no chromosome has a complete all-elements job")
        return 0
    saved = 0.0
    for chrom in chroms:
        d = Path(args.cache) / chrom
        before = sum(f.stat().st_size for f in d.glob("*.json")) if d.exists() else 0
        r = pack(chrom, Path(args.cache), remove=not args.keep_files)
        if not r["elements"]:
            print(f"{chrom}: no cached answers")
            continue
        after = r.get("archive_mb", 0) * 1e6
        if before:
            saved += before - after
        print(
            f"{chrom}: {r['elements']:,} answers in {r['archive_mb']} MB"
            + (
                f" (from {before / 1e6:.0f} MB in {r['packed']:,} files, {r['removed']:,} removed)"
                if before
                else ""
            )
        )
    if saved:
        print(f"freed {saved / 1e9:.2f} GB")
    return 0


if __name__ == "__main__":
    sys.exit(main())
