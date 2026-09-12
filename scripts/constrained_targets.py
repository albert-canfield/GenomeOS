# SPDX-License-Identifier: AGPL-3.0-or-later
"""Targets for the constrained enhancers: where the composition budget says sequence is conserved and
the registry says it is regulatory, ask the deletion tool which gene it reaches.

The budget's regulatory tier holds constrained bases with no target named. This script ranks a
chromosome's distal enhancer-like elements by their phyloP-constrained fraction (Zoonomia, through the
bigWig range reader), scores the most constrained ones by deleting each in AlphaGenome and reading which
gene moves (genomeos/predict/enhancer_target.py, cached per element), and compares them with the uniform
sample of 200 the genome-wide job scored: do constrained elements name a gene more often, more strongly,
and inside the node more often? Predicted evidence on top of measured constraint.

    uv run python scripts/constrained_targets.py --chrom chr21 [--top 100]
"""

from __future__ import annotations

import argparse
import sys
import time

from genomeos.attribution.constraint import phylop_over_blocks
from genomeos.jobs import heartbeat
from genomeos.predict import AlphaGenomeAdapter, status
from genomeos.predict.enhancer_target import Context, load_cached, summarise
from genomeos.results import load_result, save_result


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument("--top", type=int, default=100, help="most constrained elements to score")
    ap.add_argument(
        "--min-constrained", type=float, default=0.2, help="constrained fraction an element needs"
    )
    args = ap.parse_args()
    chrom, name = args.chrom, f"constrained_targets_{args.chrom}"
    st = status()
    if not st["enabled"]:
        print(f"AlphaGenome is disabled: {st['reason']}. {st['how']}")
        return 2
    ctx = Context(chrom)
    distal = ctx.distal_enhancers()
    distal.sort(key=lambda e: e.locus.start)
    intervals, keep = [], []
    last_end = -1
    for e in distal:  # the bigWig summariser wants non-overlapping intervals
        if e.locus.start >= last_end:
            intervals.append((e.locus.start, e.locus.end))
            keep.append(e)
            last_end = e.locus.end
    t0 = time.time()
    stats, cost = phylop_over_blocks(chrom, intervals, progress=None)
    print(
        f"{chrom}: phyloP over {len(keep):,} distal enhancers in {time.time() - t0:.0f} s ({cost})",
        flush=True,
    )
    ranked = sorted(zip(keep, stats, strict=False), key=lambda es: -(es[1].fraction_above or 0.0))
    chosen = [(e, s) for e, s in ranked if (s.fraction_above or 0.0) >= args.min_constrained][: args.top]
    print(
        f"{chrom}: {len(chosen)} elements with ≥ {args.min_constrained:.0%} constrained bases to score",
        flush=True,
    )
    adapter = AlphaGenomeAdapter()
    scorer = None
    rows = []
    for i, (e, s) in enumerate(chosen, 1):
        heartbeat(name)
        while True:
            try:
                if scorer is None and load_cached(chrom, e.id) is None:
                    scorer = adapter._live_scorer(threshold=0.0)  # noqa: SLF001
                r = ctx.score(scorer or (lambda *a: []), e)
                break
            except Exception as ex:  # noqa: BLE001 - retry, never give up
                print(f"  {e.id}: {type(ex).__name__}: {str(ex)[:100]}; retry in 30s", flush=True)
                time.sleep(30)
                scorer = None
                adapter = AlphaGenomeAdapter()
        r["constrained_fraction"] = round(s.fraction_above or 0.0, 3)
        r["phylop_mean"] = round(s.mean, 3) if s.bases else None
        rows.append(r)
        p = r["predicted"]
        tgt = f"{p['gene']} ({p['log2_fold_change']:+.2f} in {p['tissue'][:24]})" if p else "no effect"
        print(
            f"{chrom}: {i}/{len(chosen)}  {e.id} constrained {r['constrained_fraction']:.0%} -> {tgt}; "
            f"{r['verdict_coding']}",
            flush=True,
        )
    cons = summarise(rows)
    uniform = (load_result(f"enhancer_targets_{chrom}") or {}).get("summary", {})
    strong_named = [abs(r["predicted"]["log2_fold_change"]) for r in rows if r.get("predicted")]
    out = {
        "chrom": chrom,
        "distal_enhancers": len(keep),
        "min_constrained_fraction": args.min_constrained,
        "constrained_elements_available": sum(
            1 for _, s in ranked if (s.fraction_above or 0.0) >= args.min_constrained
        ),
        "scored": len(rows),
        "summary": cons,
        "mean_abs_log2fc_when_named": round(sum(strong_named) / len(strong_named), 3)
        if strong_named
        else None,
        "uniform_sample": {
            k: uniform.get(k)
            for k in (
                "elements_scored",
                "fraction_with_target",
                "coding_target_agrees_with_nearest",
                "fraction_inside_domain",
                "strong",
                "silencer_like",
            )
        },
        "elements": [
            {
                k: r.get(k)
                for k in (
                    "id",
                    "start",
                    "end",
                    "constrained_fraction",
                    "phylop_mean",
                    "inferred",
                    "predicted",
                    "predicted_coding",
                    "verdict_coding",
                    "domain",
                )
            }
            for r in rows
        ],
        "phylop_cost": cost,
        "evidence": "measured: Zoonomia phyloP constraint; predicted: AlphaGenome deletion effect "
        "per element",
    }
    save_result(name, out)
    print(
        f"done: constrained {cons['fraction_with_target']:.0%} name a gene "
        f"(uniform {uniform.get('fraction_with_target')}); "
        f"inside node {cons['fraction_inside_domain']} (uniform {uniform.get('fraction_inside_domain')}); "
        f"strong {cons['strong']} of {cons['with_predicted_target']}",
        flush=True,
    )
    ctx.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
