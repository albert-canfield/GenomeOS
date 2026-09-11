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
        r = None
        for attempt in range(3):
            try:
                r = compile_chromosome(c, log=sys.stdout)
                break
            except Exception as e:  # noqa: BLE001  (a source outage must not end the whole run)
                print(f"{c}: attempt {attempt + 1} failed ({str(e)[:80]}); waiting 60 s", flush=True)
                time.sleep(60)
        if r is None:
            print(f"{c}: skipped after 3 attempts", flush=True)
            continue
        save_result(f"proteome_{c}", r)
        # with local models and sequence, also verify our translation against UniProt for this chromosome
        from pathlib import Path

        from genomeos.genome import Annotation, IndexedGenome, default_gencode
        from genomeos.molecules.verify import verify_chromosome

        fa = Path(f"data/reference/{c}.fa.gz")
        gff = default_gencode({c})
        if gff and fa.exists():
            ann = Annotation.from_gff3(gff, {c})
            genome = IndexedGenome(fa)
            try:
                v = verify_chromosome(c, ann, genome)
            finally:
                genome.close()
            save_result(f"translation_vs_uniprot_{c}", v)
            print(
                f"{c}: translation vs UniProt: {v['exact_some_isoform_fraction']:.1%} identical for some "
                f"isoform, {v['disagreement_count']} real disagreements",
                flush=True,
            )
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
