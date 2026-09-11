"""Compile the proteome of every chromosome, smallest first; keep one coverage
summary per chromosome (data/results/proteome_<chrom>.json) and the compiled
definitions in the local knowledge cache. Resumable: a chromosome with a
result is skipped. Run through the jobs registry (Progress tab) or directly."""

from __future__ import annotations

import sys
import time

from genomeos.molecules.proteome import CHROM_LENGTHS, compile_chromosome
from genomeos.results import load_result, save_result

ORDER = sorted(CHROM_LENGTHS, key=lambda c: CHROM_LENGTHS[c])


def main() -> None:
    chroms = [c for c in ORDER if not load_result(f"proteome_{c}")]
    print(f"{len(ORDER) - len(chroms)} chromosomes already compiled; {len(chroms)} to go: {' '.join(chroms)}")
    t0 = time.time()
    for c in chroms:
        r = compile_chromosome(c, log=sys.stdout)
        save_result(f"proteome_{c}", r)
        cov = r["coverage_fraction"]
        print(
            f"{c}: {r['coding_genes']} genes in {r['seconds']} s; sequence {cov['sequence']:.1%}, "
            f"function {cov['function']:.1%}, pathways {cov['pathways']:.1%}, "
            f"experimental structure {cov['structure_experimental']:.1%}  (total {time.time() - t0:.0f} s)",
            flush=True,
        )
    print("done: every chromosome compiled")


if __name__ == "__main__":
    main()
