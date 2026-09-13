# SPDX-License-Identifier: AGPL-3.0-or-later
"""GC and interspersed-repeat share of every canonical promoter (TSS ± 1 kb), all chromosomes.

The operator test's null (genome/motifs.py `distil`): expectations are taken within strata of GC
and repeat content, so GC-rich profiles co-hitting and transposon-matching profiles co-hitting are
explained by composition rather than reported as regulatory logic.

    uv run python scripts/promoter_composition_genome_wide.py

Resumable per chromosome through a partial file; the result is `promoter_composition_genome_wide.json`.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

from genomeos.genome.motifs import promoter_composition
from genomeos.jobs import heartbeat
from genomeos.results import save_result

JOB = "promoter_composition_genome_wide"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
PARTIAL = Path("data/knowledge/jaspar/promoter_composition_partial.json")


def main() -> None:
    done: dict[str, dict] = json.loads(PARTIAL.read_text()) if PARTIAL.exists() else {}
    for chrom in CHROMS:
        if chrom in done:
            continue
        t0 = time.time()
        done[chrom] = promoter_composition(chrom)
        PARTIAL.parent.mkdir(parents=True, exist_ok=True)
        PARTIAL.write_text(json.dumps(done))
        heartbeat(JOB)
        print(f"{chrom}: {len(done[chrom])} promoters, {time.time() - t0:.0f} s", flush=True)
    genes: dict[str, dict] = {}
    for per in done.values():
        genes.update(per)
    gcs = sorted(v["gc"] for v in genes.values())
    reps = sorted(v["repeat"] for v in genes.values() if v.get("repeat") is not None)
    save_result(
        JOB,
        {
            "promoters": len(genes),
            "flank": 1000,
            "gc_median": gcs[len(gcs) // 2] if gcs else None,
            "repeat_median": reps[len(reps) // 2] if reps else None,
            "repeat_free_share": round(sum(1 for r in reps if r == 0) / len(reps), 4) if reps else None,
            "genes": genes,
            "evidence": (
                "curated: GRCh38 sequence around the canonical TSS of GENCODE 50 coding genes; "
                "interspersed repeats (SINE, LINE, LTR, DNA, retroposon) from RepeatMasker"
            ),
        },
    )
    print(f"genome: {len(genes)} promoters, median GC {gcs[len(gcs) // 2] if gcs else None}", flush=True)


if __name__ == "__main__":
    main()
