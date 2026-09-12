# SPDX-License-Identifier: AGPL-3.0-or-later
"""Hold the attribution against VISTA's measured enhancers on one chromosome.

    uv run python scripts/vista_score.py --chrom chr21              # registry, node, constraint, deletion
    uv run python scripts/vista_score.py --chrom chr21 --no-model   # without AlphaGenome

Saves data/results/vista_<chrom>.json: positives against negatives on what GenomeOS says blind.
"""

from __future__ import annotations

import argparse
import sys
import time

from genomeos.attribution import vista
from genomeos.attribution.constraint import phylop_over_blocks
from genomeos.jobs import heartbeat
from genomeos.predict import AlphaGenomeAdapter, status
from genomeos.predict.enhancer_target import Context, load_cached
from genomeos.results import save_result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument("--no-model", action="store_true", help="skip the AlphaGenome deletions")
    args = ap.parse_args()
    chrom, name = args.chrom, f"vista_{args.chrom}"
    elements = vista.load_loci(chrom)
    if not elements:
        print(f"{chrom}: no VISTA elements")
        return 1
    ctx = Context(chrom)
    intervals, keep = [], []
    last_end = -1
    for e in elements:  # the bigWig summariser wants non-overlapping intervals
        if e.start >= last_end:
            intervals.append((e.start, e.end))
            keep.append(e)
            last_end = e.end
    t0 = time.time()
    stats, cost = phylop_over_blocks(chrom, intervals)
    print(f"{chrom}: phyloP over {len(keep)} VISTA elements in {time.time() - t0:.0f} s ({cost})", flush=True)
    rows = vista.annotate(ctx, keep, stats)
    if not args.no_model:
        st = status()
        if not st["enabled"]:
            print(f"AlphaGenome is disabled: {st['reason']}; scoring without the model")
        else:
            adapter = AlphaGenomeAdapter()
            scorer = None
            for i, row in enumerate(rows, 1):
                heartbeat(name)
                while True:
                    try:
                        if scorer is None and load_cached(chrom, row["id"]) is None:
                            scorer = adapter._live_scorer(threshold=0.0)  # noqa: SLF001
                        vista.score(ctx, scorer or (lambda *a: []), row)
                        break
                    except Exception as ex:  # noqa: BLE001 - retry, never give up
                        print(
                            f"  {row['id']}: {type(ex).__name__}: {str(ex)[:100]}; retry in 30s", flush=True
                        )
                        time.sleep(30)
                        scorer = None
                        adapter = AlphaGenomeAdapter()
                p = row.get("predicted")
                tgt = (
                    f"{p['gene']} ({p['log2_fold_change']:+.2f} in {p['tissue'][:24]})" if p else "no effect"
                )
                print(
                    f"{chrom}: {i}/{len(rows)}  {row['id']} {row['status']} "
                    f"{';'.join(row['tissues']) or '-'} -> {tgt}; "
                    f"group {row.get('predicted_group')} agrees {row.get('tissue_agrees')}",
                    flush=True,
                )
    vista.save_knowledge(chrom, rows)
    summary = vista.summarise(rows)
    out = {
        "chrom": chrom,
        "elements": len(elements),
        "annotated": len(rows),
        "summary": summary,
        "rows": [
            {
                k: r.get(k)
                for k in (
                    "id",
                    "start",
                    "end",
                    "status",
                    "tissues",
                    "ccre_classes",
                    "domain",
                    "inferred",
                    "constrained_fraction",
                    "predicted",
                    "verdict_coding",
                    "predicted_group",
                    "tissue_agrees",
                    "skipped",
                )
            }
            for r in rows
        ],
        "phylop_cost": cost,
        "evidence": f"{vista.EVIDENCE}; measured: Zoonomia phyloP; curated: ENCODE cCREs; "
        "inferred: CTCF-only nodes; predicted: AlphaGenome deletion effect",
    }
    save_result(name, out)
    pos, neg = summary["positive"], summary["negative"]
    print(
        f"done: {pos['elements']} positive / {neg['elements']} negative; "
        f"enhancer-like {pos['enhancer_like']} / {neg['enhancer_like']}; "
        f"constrained {pos['constrained']} / {neg['constrained']}; "
        f"named {pos['fraction_with_target']} / {neg['fraction_with_target']}; "
        f"tissue agrees {pos['tissue_agrees']} of {pos['tissue_judged']} judged"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
