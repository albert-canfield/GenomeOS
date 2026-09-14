# SPDX-License-Identifier: AGPL-3.0-or-later
"""Code chromosomes under generic models and the project's knowledge, and measure bits per base.

Each argument is CHROM:TRAIN, the chromosome to code and the chromosomes (comma-separated) its
static models are fitted on.
No model API is called; the only request is UCSC's CpG islands, once per chromosome, cached under
data/knowledge/compress. The summary is committed as data/results/compress_<chrom>.json.

    uv run python scripts/compress.py chr21:chr22 chr22:chr21
"""

import argparse
import time
from pathlib import Path

from genomeos.attribution.compress import run_and_save


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pairs", nargs="*", default=["chr21:chr22", "chr22:chr21"])
    ap.add_argument("--no-markov", action="store_true", help="skip the order-0..16 baseline table")
    ap.add_argument("--no-compressors", action="store_true", help="skip xz and bzip2")
    ap.add_argument("--spill", help="directory to keep model codes in files instead of memory")
    ap.add_argument("--name", help="result name for a sensitivity run (default compress_<chrom>)")
    args = ap.parse_args()
    for pair in args.pairs:
        chrom, train = pair.split(":")
        path = run_and_save(
            chrom,
            train,
            progress=lambda m: print(m, flush=True),
            markov=not args.no_markov,
            compressors=not args.no_compressors,
            name=args.name,
            spill=Path(args.spill) if args.spill else None,
        )
        print(f"{time.strftime('%H:%M:%S')} {chrom}: {path}", flush=True)


if __name__ == "__main__":
    main()
