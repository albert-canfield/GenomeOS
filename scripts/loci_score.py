# SPDX-License-Identifier: AGPL-3.0-or-later
"""Delete every enhancer inside the known-locus panel's windows, and nothing else.

    uv run python scripts/loci_score.py                     # every locus that still needs asking
    uv run python scripts/loci_score.py --locus HBB_LCR     # one locus
    uv run python scripts/loci_score.py --dry-run           # what it would ask, per locus

The chromosome sweep costs about 770 hours of model time genome-wide; the twelve loci of
docs/LOCI-BENCHMARK.md need 2,222 elements. This driver spends the quota on those windows only.
It drives the existing scorer (`genomeos/predict/enhancer_target.py`) exactly as
`scripts/enhancer_targets_all.py` does - the same `Context`, the same per-element cache under
data/knowledge/alphagenome/elements, the same shared pacer for quota answers - and rewrites
nothing of it, so when the sweep resumes a chromosome it finds these answers in its cache.

Order: the elements overlapping the published element of each locus first, so the element-level
verdict (does the ZRS's own deletion name SHH?) lands in the first seconds, then the rest of that
locus's window, smallest window first. A chromosome the sweep has already finished is skipped,
because `attribution/targets.py` already reads it.

The result (`loci_deletions`) is the per-element table itself rather than a summary: 2,222 rows are
small enough to commit, and the benchmark on another machine then needs no local cache.
"""

from __future__ import annotations

import argparse
import json
import sys
import threading
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from enhancer_targets_all import (  # noqa: E402 - the worked example's machinery, imported not copied
    CHECK_EVERY,
    IN_FLIGHT,
    STALL_EXIT,
    WORKERS,
    Pacer,
    _locked,
    compact,
    retry_seconds,
    worker_scorer,
)

from genomeos.benchmark.loci import PANEL  # noqa: E402
from genomeos.predict import status  # noqa: E402
from genomeos.predict.enhancer_target import Context, cache_path, has_cells, load_cached  # noqa: E402
from genomeos.results import RESULTS_DIR, load_result, save_result  # noqa: E402

NAME = "loci_deletions"
SAVE_EVERY = 100


def windows(only: str | None = None) -> list:
    """The panel's loci, smallest window first: the order early answers land in."""
    loci = [e for e in PANEL if only is None or e.locus == only]
    return sorted(loci, key=lambda e: e.window[1] - e.window[0])


def already_swept(chrom: str, results_dir: Path = RESULTS_DIR) -> bool:
    """Whether the chromosome sweep has finished this chromosome (then nothing is asked for it)."""
    r = load_result(f"enhancer_targets_all_{chrom}", results_dir) or {}
    return bool(r.get("complete"))


def to_ask(ctx: Context, expect) -> list:
    """The enhancer elements of one locus window, those over the published element first."""
    w0, w1 = expect.window
    e0, e1 = expect.element
    els = [
        x
        for x in ctx.elements
        if x.cls == "enhancer" and x.domain and x.locus.start < w1 and x.locus.end > w0
    ]
    els.sort(key=lambda x: (not (x.locus.start < e1 and x.locus.end > e0), x.locus.start))
    return els


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--locus", help="one locus of the panel (default: all that need asking)")
    ap.add_argument("--workers", type=int, default=WORKERS, help="requests in flight at once")
    ap.add_argument("--dry-run", action="store_true", help="count the requests, ask nothing")
    args = ap.parse_args()
    st = status()
    if not args.dry_run and not st["enabled"]:
        print(f"AlphaGenome is disabled: {st['reason']}. {st['how']}")
        return 2

    kept: dict[str, dict] = {}  # element id -> row, over every locus
    per_locus: dict[str, dict] = {}
    old = load_result(NAME) or {}
    for row in old.get("elements", []):
        kept[row["id"]] = row
    pacer = Pacer()
    local = threading.local()
    requests = 0
    counted = threading.Lock()
    t0 = time.time()

    def save() -> None:
        rows = sorted(kept.values(), key=lambda r: (r["chrom"], r["start"]))
        save_result(
            NAME,
            {
                "note": (
                    "AlphaGenome deletions inside the known-locus panel's windows only"
                    " (docs/LOCI-BENCHMARK.md); the chromosome sweep is not run for these"
                ),
                "elements": rows,
                "loci": per_locus,
                "requests_this_run": requests,
                "workers": args.workers,
                "cells": ["K562", "HepG2", "GM12878", "IMR-90"],
                "evidence": (
                    "predicted: AlphaGenome deletion effect per element, per gene and per cell line"
                    " (K562, HepG2, GM12878, IMR-90 RNA-seq tracks); inferred: CTCF-only nodes"
                ),
            },
        )

    for expect in windows(args.locus):
        chrom = expect.chrom
        if already_swept(chrom):
            per_locus[expect.locus] = {"chrom": chrom, "skipped": "the chromosome sweep finished it"}
            print(f"{expect.locus}: {chrom} is swept; nothing to ask", flush=True)
            continue
        ctx = Context(chrom)
        fetch_lock = threading.Lock()
        raw_fetch = ctx.genome.fetch
        ctx.genome.fetch = lambda locus, _r=raw_fetch, _l=fetch_lock: _locked(_l, _r, locus)
        elements = to_ask(ctx, expect)
        mine = [e for e in elements if e.id not in kept]
        print(
            f"{expect.locus}: {len(elements)} elements in {(expect.window[1] - expect.window[0]) / 1e6:.2f}"
            f" Mb, {len(mine)} to score",
            flush=True,
        )
        if args.dry_run:
            per_locus[expect.locus] = {"chrom": chrom, "elements": len(elements), "to_score": len(mine)}
            ctx.close()
            continue
        spent0, done = requests, 0
        t1 = time.time()

        def score_one(e, _ctx=ctx):
            nonlocal requests
            while True:
                pacer.wait()
                cached = load_cached(_ctx.chrom, e.id)
                if cached is not None and not has_cells(cached):
                    cache_path(_ctx.chrom, e.id).unlink(missing_ok=True)
                    cached = None
                try:
                    if cached is None:
                        with counted:
                            requests += 1
                        r = _ctx.score(worker_scorer(local), e)
                    else:
                        r = _ctx.score(lambda *a: [], e)
                    pacer.clear()
                    return r
                except Exception as ex:  # noqa: BLE001 - wait out the quota, retry the rest
                    held = retry_seconds(str(ex))
                    if held is not None:
                        waited = pacer.hold(min(120, held))
                        if waited > 300:
                            print(f"  quota exhausted: holding {waited / 60:.0f} min", flush=True)
                    else:
                        print(f"  {e.id}: {type(ex).__name__}: {str(ex)[:100]}; retrying", flush=True)
                        time.sleep(5)
                    local.scorer = None

        todo = iter(mine)
        pool = ThreadPoolExecutor(max_workers=max(1, args.workers))
        pending = set()
        for _ in range(max(1, args.workers) * IN_FLIGHT):
            nxt = next(todo, None)
            if nxt is None:
                break
            pending.add(pool.submit(score_one, nxt))
        last = time.time()
        while pending:
            finished, pending = wait(pending, timeout=CHECK_EVERY, return_when=FIRST_COMPLETED)
            if not finished:
                idle = time.time() - last
                if idle > STALL_EXIT:
                    save()
                    print(f"{expect.locus}: nothing scored for {idle / 60:.0f} min; saving and stopping")
                    return 3
                print(f"{expect.locus}: waiting ({idle / 60:.0f} min since the last element)", flush=True)
                continue
            for fut in finished:
                r = fut.result()
                kept[r["id"]] = {**compact(r), "chrom": chrom}
                done += 1
                last = time.time()
                if done % SAVE_EVERY == 0:
                    save()
                    p = r.get("predicted")
                    print(
                        f"{expect.locus}: {done}/{len(mine)} ({requests - spent0} requests,"
                        f" {(time.time() - t1) / 60:.1f} min); last {r['id']} ->"
                        f" {p['gene'] if p else 'no effect'}",
                        flush=True,
                    )
                nxt = next(todo, None)
                if nxt is not None:
                    pending.add(pool.submit(score_one, nxt))
        pool.shutdown(wait=False)
        ctx.close()
        per_locus[expect.locus] = {
            "chrom": chrom,
            "elements": len(elements),
            "scored_this_run": done,
            "requests": requests - spent0,
            "minutes": round((time.time() - t1) / 60, 1),
        }
        save()
        print(
            f"{expect.locus}: done, {requests - spent0} requests in {(time.time() - t1) / 60:.1f} min",
            flush=True,
        )

    if not args.dry_run:
        save()
    print(
        f"{len(kept)} elements in the table, {requests} requests this run, "
        f"{pacer.waits} quota waits, {(time.time() - t0) / 60:.1f} min"
    )
    print(json.dumps(per_locus, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
