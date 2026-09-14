# SPDX-License-Identifier: AGPL-3.0-or-later
"""Code chromosomes under generic models and the project's knowledge, and measure bits per base.

Each argument is CHROM:TRAIN, the chromosome to code and the chromosomes (comma-separated) its
static models are fitted on. No model API is called; the only request is UCSC's CpG islands, once
per chromosome, cached under data/knowledge/compress.

    uv run python scripts/compress.py chr21:chr22 chr22:chr21   # full stack, compress_<chrom>
    uv run python scripts/compress.py --pass chr1:chr22         # verdict pass, compress_pass_<chrom>
    uv run python scripts/compress.py --rollup                  # compress_genome_wide from the passes

The full stack adds the order 0-16 table, xz and bzip2, the mixing grid and the partner-primed
control; the verdict pass produces every stack and layer verdict in bounded memory and disk.
"""

import argparse
import signal
import sys
import time
from pathlib import Path

from genomeos.attribution.compress import (
    DISK_FLOOR_GB,
    KARYOTYPE,
    SEGMENT_BASES,
    DiskFloorError,
    rollup,
    run_and_save,
    run_pass,
)
from genomeos.results import save_result


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("pairs", nargs="*", default=["chr21:chr22", "chr22:chr21"])
    ap.add_argument("--pass", dest="verdict", action="store_true", help="the bounded verdict pass")
    ap.add_argument("--rollup", action="store_true", help="sum the verdict passes into the genome")
    ap.add_argument("--no-markov", action="store_true", help="skip the order-0..16 baseline table")
    ap.add_argument("--no-compressors", action="store_true", help="skip xz and bzip2")
    ap.add_argument("--spill", help="directory for temporary files (default data/cache/compress)")
    ap.add_argument("--segment", type=int, default=SEGMENT_BASES, help="bases mixed at a time (pass)")
    ap.add_argument("--floor-gb", type=float, default=DISK_FLOOR_GB, help="free disk never to cross")
    ap.add_argument("--name", help="result name for a sensitivity run")
    ap.add_argument("--heartbeat", help="job name to beat for (pass)")
    args = ap.parse_args()

    def say(m: str) -> None:
        print(m, flush=True)

    if args.rollup:
        path = save_result(args.name or "compress_genome_wide", rollup(list(KARYOTYPE)))
        say(f"{time.strftime('%H:%M:%S')} rollup: {path}")
        return
    for pair in args.pairs:
        chrom, train = pair.split(":")
        if args.verdict:
            # a pass started in the background inherits SIGINT ignored: say explicitly that both
            # signals stop it through its cleanup, which deletes the temporaries
            def stop(signum, frame):
                raise DiskFloorError(f"signal {signum}")

            signal.signal(signal.SIGTERM, stop)
            signal.signal(signal.SIGINT, stop)
            beat = None
            if args.heartbeat:
                from genomeos.jobs import heartbeat

                beat = lambda: heartbeat(args.heartbeat)  # noqa: E731
            try:
                result = run_pass(
                    chrom,
                    train,
                    progress=say,
                    segment=args.segment,
                    spill=Path(args.spill) if args.spill else Path("data/cache/compress"),
                    floor_gb=args.floor_gb,
                    heartbeat=beat,
                )
            except DiskFloorError as e:
                say(f"{time.strftime('%H:%M:%S')} {chrom}: stopped: {e}")
                sys.exit(3)
            path = save_result(args.name or f"compress_pass_{chrom}", result)
        else:
            path = run_and_save(
                chrom,
                train,
                progress=say,
                markov=not args.no_markov,
                compressors=not args.no_compressors,
                name=args.name,
                spill=Path(args.spill) if args.spill else None,
            )
        say(f"{time.strftime('%H:%M:%S')} {chrom}: {path}")


if __name__ == "__main__":
    main()
