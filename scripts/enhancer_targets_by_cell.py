# SPDX-License-Identifier: AGPL-3.0-or-later
"""Add each element's deletion effect on the cell lines' own tracks to a chromosome's committed results.

    uv run python scripts/enhancer_targets_by_cell.py --chrom chr21

The enhancer_targets_<chrom> and constrained_targets_<chrom> results get, per element, `predicted_by_cell`
and `predicted_coding_by_cell`: {cell: log2 fold change of the named target on that cell line's RNA-seq
track} for K562, HepG2, GM12878 and IMR-90. Elements whose cached answer predates the per-cell field are
scored again (one request each); everything else is read from the cache. The closure test then takes the
cell's own magnitude instead of the tissue that moved most.
"""

from __future__ import annotations

import argparse
import sys
import time

from genomeos.jobs import heartbeat
from genomeos.predict import AlphaGenomeAdapter, status
from genomeos.predict.enhancer_target import Context, cache_path, has_cells, load_cached
from genomeos.results import load_result, save_result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument("--results", default="enhancer_targets,constrained_targets")
    args = ap.parse_args()
    chrom = args.chrom
    ctx = Context(chrom)
    adapter = AlphaGenomeAdapter()
    scorer = None
    rescored = 0
    for prefix in args.results.split(","):
        name = f"{prefix}_{chrom}"
        res = load_result(name)
        if not res:
            print(f"{name}: no result")
            continue
        for i, e in enumerate(res["elements"], 1):
            heartbeat("enhancer_targets_by_cell")
            el = ctx.element_by_id.get(e["id"])
            if el is None:
                continue
            hit = load_cached(chrom, e["id"])
            if not has_cells(hit):
                st = status()
                if not st["enabled"]:
                    print(f"AlphaGenome is disabled: {st['reason']}")
                    return 2
                cache_path(chrom, e["id"]).unlink(missing_ok=True)
                while True:
                    try:
                        if scorer is None:
                            scorer = adapter._live_scorer(threshold=0.0)  # noqa: SLF001
                        r = ctx.score(scorer, el)
                        break
                    except Exception as ex:  # noqa: BLE001 - retry, never give up
                        print(f"  {e['id']}: {type(ex).__name__}: {str(ex)[:100]}; retry in 30s", flush=True)
                        time.sleep(30)
                        scorer = None
                        adapter = AlphaGenomeAdapter()
                rescored += 1
            else:
                r = ctx.score(lambda *a: [], el)
            e["predicted_by_cell"] = r.get("predicted_by_cell")
            e["predicted_coding_by_cell"] = r.get("predicted_coding_by_cell")
            if i % 20 == 0:
                print(f"{name}: {i}/{len(res['elements'])} ({rescored} scored again)", flush=True)
        res["cells"] = ["K562", "HepG2", "GM12878", "IMR-90"]
        save_result(name, {k: v for k, v in res.items() if k not in ("result", "date")})
        with_cells = sum(1 for e in res["elements"] if e.get("predicted_by_cell"))
        print(f"done {name}: {with_cells} of {len(res['elements'])} named elements carry per-cell effects")
    return 0


if __name__ == "__main__":
    sys.exit(main())
