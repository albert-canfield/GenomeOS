# SPDX-License-Identifier: AGPL-3.0-or-later
"""Predicted enhancer targets for one chromosome (AlphaGenome feature b): a sampled, resumable job.

Takes every distal enhancer-like element (dELS) whose reach the domain inference
already names, samples them evenly along the chromosome, deletes each one in
AlphaGenome and reads which gene moves. Per-element answers are cached under
data/knowledge/alphagenome/elements, so a rerun costs nothing for what is done;
the summary lands in data/results/enhancer_targets_<chrom>.json after every ten
elements and at the end. Loops until the sample is complete.

    uv run python scripts/enhancer_targets.py --chrom chr21 [--sample 200]
"""

from __future__ import annotations

import argparse
import sys
import time

from genomeos.jobs import heartbeat
from genomeos.predict import AlphaGenomeAdapter, status
from genomeos.predict.enhancer_target import MIN_EFFECT, Context, load_cached, summarise
from genomeos.results import load_result, save_result

MIN_DISTANCE = 2_000  # promoters are another question; keep elements clear of the TSS they reach
MAX_DISTANCE = 400_000  # stay well inside the 1 Mb window AlphaGenome sees


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--sample", type=int, default=200, help="elements to score, spread along the chromosome")
    ap.add_argument("--min-effect", type=float, default=MIN_EFFECT)
    args = ap.parse_args()
    chrom = args.chrom
    name = f"enhancer_targets_{chrom}"
    st = status()
    if not st["enabled"]:
        print(f"AlphaGenome is disabled: {st['reason']}. {st['how']}")
        return 2

    try:
        ctx = Context(chrom)
    except FileNotFoundError as ex:
        print(ex)
        return 1
    distal = ctx.distal_enhancers(MIN_DISTANCE, MAX_DISTANCE)
    n = min(args.sample, len(distal))
    step = len(distal) / n if n else 1
    picked = [distal[int(i * step)] for i in range(n)]
    print(f"{chrom}: {len(ctx.elements):,} elements, {len(distal):,} distal enhancers; {n} to score")

    adapter = AlphaGenomeAdapter()
    scorer = None
    rows: list[dict] = []
    done = 0
    failures = 0
    for i, e in enumerate(picked, 1):
        heartbeat(name)
        while True:
            try:
                if scorer is None and load_cached(chrom, e.id) is None:
                    scorer = adapter._live_scorer(threshold=0.0)  # noqa: SLF001
                r = ctx.score(scorer or (lambda *a: []), e, min_effect=args.min_effect)
                break
            except Exception as ex:  # noqa: BLE001 - the API drops calls; wait and retry, never give up
                failures += 1
                wait = min(300, 10 * failures)
                print(f"  {e.id}: {type(ex).__name__}: {str(ex)[:120]}; retry in {wait}s", flush=True)
                time.sleep(wait)
                scorer = None
                adapter = AlphaGenomeAdapter()
        rows.append(r)
        done += 1
        p = r["predicted"]
        tgt = (
            f"{p['gene']} ({p['action']}, {p['log2_fold_change']:+.2f} in {p['tissue']})"
            if p
            else "no effect"
        )
        pc = r.get("predicted_coding")
        cod = (
            f"; coding {pc['gene']} {pc['log2_fold_change']:+.2f}"
            if pc and pc["gene"] != (p or {}).get("gene")
            else ""
        )
        inferred = e.targets[0]["gene"]
        print(
            f"{chrom}: {i}/{n}  {e.id} -> {tgt}{cod}; inferred {inferred}; {r['verdict_coding']}", flush=True
        )
        if i % 10 == 0 or i == n:
            write(chrom, rows, len(distal), len(ctx.elements), st["model"], args.min_effect)
    ctx.close()
    print(f"{chrom}: done, {done} elements scored, {failures} retries")
    return 0


def write(chrom: str, rows: list[dict], distal: int, elements: int, model: str, min_effect: float) -> None:
    prev = load_result(f"enhancer_targets_{chrom}") or {}
    payload = {
        "chrom": chrom,
        "model": model,
        "evidence": "predicted",
        "method": "delete the element as a variant in its 1 Mb window; RNA-seq gene scorer over all tracks",
        "min_effect_log2fc": min_effect,
        "elements_on_chromosome": elements,
        "distal_enhancers_with_domain_target": distal,
        "summary": summarise(rows),
        "elements": [
            {
                "id": r["id"],
                "start": r["start"],
                "end": r["end"],
                "domain": r.get("domain", ""),
                "inferred": r.get("inferred"),
                "predicted": r.get("predicted"),
                "predicted_coding": r.get("predicted_coding"),
                "verdict": r.get("verdict"),
                "verdict_coding": r.get("verdict_coding"),
                "top_genes": [
                    {"gene": g["gene"], "drop": g["max_drop_log2fc"], "tissue": g["max_drop_tissue"]}
                    for g in r.get("top_genes", [])[:3]
                ],
            }
            for r in rows
        ],
    }
    if prev.get("summary", {}).get("elements_scored", 0) > len(rows):
        return  # never overwrite a larger finished sample with a smaller one
    save_result(f"enhancer_targets_{chrom}", payload)


if __name__ == "__main__":
    sys.exit(main())
