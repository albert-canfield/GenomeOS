# SPDX-License-Identifier: AGPL-3.0-or-later
"""Every enhancer inside a node of one chromosome deleted in AlphaGenome, with the effect per cell line.

    uv run python scripts/enhancer_targets_all.py --chrom chr21

The sampled jobs scored 200 uniform and 100 constrained elements per chromosome; the gene-level closure
needs every element that reaches a gene scored on the cell lines' own tracks (K562, HepG2, GM12878,
IMR-90). This is that run: about 12,000 elements on chr21, one request each, cached per element, the
summary rewritten every 200 elements so the job is resumable and the Progress tab sees it move; the
per-element table (megabytes per chromosome) stays local under data/knowledge/alphagenome. The daily
quota is respected, not fought: a RESOURCE_EXHAUSTED answer is waited out for the seconds it names; when
it keeps answering that for QUOTA_PATIENCE seconds the run sleeps an hour before asking again, and says so.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

from genomeos.jobs import heartbeat
from genomeos.predict import AlphaGenomeAdapter, status
from genomeos.predict.enhancer_target import Context, load_cached, summarise
from genomeos.results import save_result

SAVE_EVERY = 200
ELEMENTS_DIR = Path("data/knowledge/alphagenome/all_elements")  # the per-element table, local
QUOTA_PATIENCE = 600  # seconds of continuous quota answers before the hour's sleep
QUOTA_SLEEP = 3600


def retry_seconds(message: str) -> int | None:
    """The seconds a quota answer asks for ("Quota exceeded; retry in 60s"), or None for other errors."""
    if "RESOURCE_EXHAUSTED" not in message and "Quota" not in message:
        return None
    m = re.search(r"retry in (\d+)s", message)
    return int(m.group(1)) if m else 60


def compact(r: dict) -> dict:
    keys = (
        "id",
        "start",
        "end",
        "domain",
        "inferred",
        "predicted",
        "predicted_coding",
        "predicted_by_cell",
        "predicted_coding_by_cell",
        "verdict_coding",
    )
    return {k: r.get(k) for k in keys}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument("--limit", type=int, help="stop after this many elements (default: all)")
    args = ap.parse_args()
    chrom, name = args.chrom, f"enhancer_targets_all_{args.chrom}"
    st = status()
    if not st["enabled"]:
        print(f"AlphaGenome is disabled: {st['reason']}. {st['how']}")
        return 2
    ctx = Context(chrom)
    distal = {e.id for e in ctx.distal_enhancers()}
    elements = [e for e in ctx.elements if e.cls == "enhancer" and e.domain]
    elements.sort(key=lambda e: (e.id not in distal, e.locus.start))  # the distal ones first
    if args.limit:
        elements = elements[: args.limit]
    adapter = AlphaGenomeAdapter()
    scorer = None
    rows: list[dict] = []
    requests = 0
    quota_since: float | None = None
    t0 = time.time()

    def save(final: bool = False) -> None:
        """The committed result is the summary; the element table (megabytes per chromosome) stays local."""
        s = summarise(rows) if rows else {}
        table = ELEMENTS_DIR / f"{chrom}.json"
        table.parent.mkdir(parents=True, exist_ok=True)
        table.write_text(json.dumps([compact(r) for r in rows]))
        save_result(
            name,
            {
                "chrom": chrom,
                "elements_total": len(elements),
                "scored": len(rows),
                "requests_this_run": requests,
                "cells": ["K562", "HepG2", "GM12878", "IMR-90"],
                "complete": final and len(rows) >= len(elements),
                "summary": s,
                "elements_where": str(table),
                "evidence": "predicted: AlphaGenome deletion effect per element, per gene and per cell line "
                "(K562, HepG2, GM12878, IMR-90 RNA-seq tracks); inferred: CTCF-only nodes",
            },
        )

    for i, e in enumerate(elements, 1):
        heartbeat(name)
        cached = load_cached(chrom, e.id)
        while True:
            try:
                if scorer is None and (
                    cached is None or not all("by_cell" in g for g in cached["genes"][:1])
                ):
                    scorer = adapter._live_scorer(threshold=0.0)  # noqa: SLF001
                if cached is not None and not all("by_cell" in g for g in cached["genes"][:1]):
                    from genomeos.predict.enhancer_target import cache_path

                    cache_path(chrom, e.id).unlink(
                        missing_ok=True
                    )  # an old cache without the per-cell fields
                    cached = None
                if cached is None:
                    requests += 1
                r = ctx.score(scorer or (lambda *a: []), e)
                quota_since = None
                break
            except Exception as ex:  # noqa: BLE001 - wait for the quota, retry the rest, never give up
                wait = retry_seconds(str(ex))
                if wait is not None:
                    quota_since = quota_since or time.time()
                    if time.time() - quota_since > QUOTA_PATIENCE:
                        print(
                            f"  quota exhausted for {QUOTA_PATIENCE // 60} min: sleeping an hour", flush=True
                        )
                        heartbeat(name)
                        time.sleep(QUOTA_SLEEP)
                        quota_since = time.time()
                    else:
                        time.sleep(min(120, wait))
                else:
                    print(f"  {e.id}: {type(ex).__name__}: {str(ex)[:100]}; retry in 30s", flush=True)
                    time.sleep(30)
                scorer = None
                adapter = AlphaGenomeAdapter()
        rows.append(r)
        if i % SAVE_EVERY == 0 or i == len(elements):
            save(final=i == len(elements))
            p = r.get("predicted")
            print(
                f"{chrom}: {i:,}/{len(elements):,} elements ({requests:,} requests, "
                f"{(time.time() - t0) / 60:.0f} min); last {e.id} -> {p['gene'] if p else 'no effect'}",
                flush=True,
            )
    save(final=True)
    s = summarise(rows)
    print(
        f"done: {len(rows):,} elements, {s['fraction_with_target']:.1%} name a gene, "
        f"{s['coding_target_agrees_with_nearest']} agree with the nearest TSS; {requests:,} requests"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
