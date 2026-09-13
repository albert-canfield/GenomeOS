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
A request that never returns (chr19 stopped that way once, silently) is cut off by the client's own
timeout and retried, so the run cannot hang on one element. The run is paced by round trips rather than
by the quota, so WORKERS requests are in flight at once and a quota answer to any of them holds all of
them; --workers 1 is the old sequential behaviour.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from genomeos.jobs import heartbeat
from genomeos.predict import AlphaGenomeAdapter, status
from genomeos.predict.enhancer_target import Context, cache_path, has_cells, load_cached, summarise
from genomeos.results import save_result

SAVE_EVERY = 200
ELEMENTS_DIR = Path("data/knowledge/alphagenome/all_elements")  # the per-element table, local
QUOTA_PATIENCE = 600  # seconds of continuous quota answers before the hour's sleep
QUOTA_SLEEP = 3600
CALL_TIMEOUT = 300  # the client's own timeout: a request that has not answered by then raises
WORKERS = 8  # requests in flight at once; 8 gave 400 requests a minute against 51 sequential, no waits


class Pacer:
    """One back-off shared by the workers: a quota answer stops all of them, not just its own."""

    def __init__(self, patience: float = QUOTA_PATIENCE, long_sleep: float = QUOTA_SLEEP) -> None:
        self.lock = threading.Lock()
        self.until = 0.0
        self.since: float | None = None
        self.patience = patience
        self.long_sleep = long_sleep
        self.waits = 0

    def wait(self) -> None:
        while True:
            with self.lock:
                left = self.until - time.time()
            if left <= 0:
                return
            time.sleep(min(left, 5))

    def hold(self, seconds: float) -> float:
        """Record a quota answer; returns how long the workers will now wait."""
        with self.lock:
            now = time.time()
            self.since = self.since or now
            self.waits += 1
            if now - self.since > self.patience:
                seconds = self.long_sleep
                self.since = now
            self.until = max(self.until, now + seconds)
            return self.until - now

    def clear(self) -> None:
        with self.lock:
            self.since = None


def _locked(lock: threading.Lock, fn, *a):
    with lock:
        return fn(*a)


def worker_scorer(local: threading.local):
    """One live scorer per worker thread, its client given the call timeout."""
    if getattr(local, "scorer", None) is None:
        from alphagenome.models import dna_client  # type: ignore[import-not-found]

        a = AlphaGenomeAdapter()
        a._client = dna_client.create(a.api_key, timeout=CALL_TIMEOUT)  # noqa: SLF001
        local.scorer = a._live_scorer(threshold=0.0)  # noqa: SLF001
    return local.scorer


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
    ap.add_argument("--workers", type=int, default=WORKERS, help="requests in flight at once")
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
    # the sequence reader is one file handle: the workers take turns at it, the requests overlap
    fetch_lock = threading.Lock()
    raw_fetch = ctx.genome.fetch
    ctx.genome.fetch = lambda locus: _locked(fetch_lock, raw_fetch, locus)
    pacer = Pacer()
    local = threading.local()
    done: dict[int, dict] = {}
    requests = 0
    counted = threading.Lock()
    t0 = time.time()

    def save(final: bool = False) -> None:
        """The committed result is the summary; the element table (megabytes per chromosome) stays local."""
        rows = [done[i] for i in sorted(done)]
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
                "workers": args.workers,
                "cells": ["K562", "HepG2", "GM12878", "IMR-90"],
                "complete": final and len(rows) >= len(elements),
                "summary": s,
                "elements_where": str(table),
                "evidence": "predicted: AlphaGenome deletion effect per element, per gene and per cell line "
                "(K562, HepG2, GM12878, IMR-90 RNA-seq tracks); inferred: CTCF-only nodes",
            },
        )

    def score_one(index: int, e) -> tuple[int, dict]:
        """One element, retried until it answers: the cache first, then a request through this worker."""
        nonlocal requests
        while True:
            pacer.wait()
            cached = load_cached(chrom, e.id)
            if cached is not None and not has_cells(cached):
                cache_path(chrom, e.id).unlink(missing_ok=True)  # cached before the per-cell fields
                cached = None
            try:
                if cached is None:
                    with counted:
                        requests += 1
                    r = ctx.score(worker_scorer(local), e)
                else:
                    r = ctx.score(lambda *a: [], e)
                pacer.clear()
                return index, r
            except Exception as ex:  # noqa: BLE001 - wait for the quota, retry the rest, never give up
                wait = retry_seconds(str(ex))
                if wait is not None:
                    held = pacer.hold(min(120, wait))
                    if held > 300:
                        print(f"  quota exhausted: holding {held / 60:.0f} min", flush=True)
                else:
                    print(f"  {e.id}: {type(ex).__name__}: {str(ex)[:100]}; retrying", flush=True)
                    time.sleep(5)
                local.scorer = None

    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        pending = {pool.submit(score_one, i, e) for i, e in enumerate(elements)}
        for fut in as_completed(pending):
            index, r = fut.result()
            done[index] = r
            n = len(done)
            heartbeat(name)
            if n % SAVE_EVERY == 0 or n == len(elements):
                save(final=n == len(elements))
                p = r.get("predicted")
                print(
                    f"{chrom}: {n:,}/{len(elements):,} elements ({requests:,} requests, "
                    f"{(time.time() - t0) / 60:.0f} min, {args.workers} workers); "
                    f"last {r.get('id')} -> {p['gene'] if p else 'no effect'}",
                    flush=True,
                )
    save(final=True)
    rows = [done[i] for i in sorted(done)]
    s = summarise(rows)
    print(
        f"done: {len(rows):,} elements, {s['fraction_with_target']:.1%} name a gene, "
        f"{s['coding_target_agrees_with_nearest']} agree with the nearest TSS; {requests:,} requests, "
        f"{pacer.waits} quota waits"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
