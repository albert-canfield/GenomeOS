# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run the registered in-silico mutagenesis over the VISTA panel (attribution/satmut_vista.py).

    uv run python scripts/satmut_vista.py                       # plan and constraint: no model request
    uv run python scripts/satmut_vista.py --run --quota-handed  # spend the budget

The design is `attribution.satmut_vista.PRE_REGISTRATION`, committed before this file existed
(0ae59ef) and amended once before any request (1d25598), when the panel turned out to match at group
level rather than pair by pair. Nothing here may change it: this assembles the windows the
registration describes, scores them through the sweep's own deletion path and cache, and reads the
result by the registered bars.

Part one tiles the panel and reads Zoonomia phyloP over every window — a bigWig range read, no model
quota — and writes `satmut_vista_plan`. Part two refuses to start without --quota-handed, holds the
shared model key, spends one request per window in a fixed order, looks once at half the budget for
futility only, and writes `satmut_vista`.
"""

from __future__ import annotations

import argparse
import json
import math
import random
import time
from contextlib import contextmanager
from pathlib import Path
from statistics import mean

from genomeos import jobs
from genomeos.attribution import satmut_vista as reg
from genomeos.predict import AlphaGenomeAdapter
from genomeos.predict.enhancer_target import Context, score_element
from genomeos.results import load_result, save_result

HOLDER = "genomeos-9c"
PLAN = Path("data/knowledge/alphagenome/satmut_vista_windows.json")
POSITIVE_GROUPS = ("limb", "neural")


def say(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


@contextmanager
def the_key(what: str):
    jobs.take_key(HOLDER, what)
    try:
        yield
    finally:
        jobs.drop_key(HOLDER)


def windows_of(record: dict) -> list[tuple[int, int]]:
    """`WINDOWS_PER_ELEMENT` windows of `WINDOW` bp, evenly spaced across the element."""
    span = record["end"] - record["start"]
    n = reg.WINDOWS_PER_ELEMENT
    if span < reg.WINDOW * 2:
        return []
    step = (span - reg.WINDOW) / max(1, n - 1)
    return [(record["start"] + int(i * step), record["start"] + int(i * step) + reg.WINDOW) for i in range(n)]


def spearman(xs: list[float], ys: list[float]) -> float | None:
    """Rank correlation, ties averaged. None when either side is constant: no order to correlate."""
    n = len(xs)
    if n < 4:
        return None

    def ranks(v: list[float]) -> list[float]:
        order = sorted(range(n), key=lambda i: v[i])
        out = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            r = (i + j) / 2 + 1
            for k in range(i, j + 1):
                out[order[k]] = r
            i = j + 1
        return out

    rx, ry = ranks(xs), ranks(ys)
    mx, my = mean(rx), mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    den = math.sqrt(sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry))
    return (num / den) if den else None


def assemble() -> dict:
    """Every window of both arms with its phyloP, and no model request."""
    from genomeos.attribution.constraint import phylop_over_blocks

    t0 = time.time()
    panel = load_result("across_panel_vista")
    if not panel:
        raise SystemExit("across_panel_vista has not been computed")
    rows: list[dict] = []
    by_chrom: dict[str, list[dict]] = {}
    for rec in panel["records"]:
        arm = "positive" if rec["group"] in POSITIVE_GROUPS else "negative"
        for start, end in windows_of(rec):
            row = {
                "element": rec["id"],
                "arm": arm,
                "group": rec["group"],
                "chrom": rec["chrom"],
                "start": start,
                "end": end,
            }
            rows.append(row)
            by_chrom.setdefault(rec["chrom"], []).append(row)

    cost = {"requests": 0, "mb": 0.0}
    for chrom, group in sorted(by_chrom.items()):
        stats, spent = phylop_over_blocks(chrom, [(r["start"], r["end"]) for r in group])
        for row, st in zip(group, stats, strict=True):
            row["phylop_mean"] = round(st.mean, 4) if st.bases else None
        cost["requests"] += spent.get("requests", 0)
        cost["mb"] += spent.get("mb_fetched", 0.0)
        say(f"{chrom}: phyloP over {len(group)} windows")

    PLAN.parent.mkdir(parents=True, exist_ok=True)
    PLAN.write_text(json.dumps(rows))
    arms = {a: sum(1 for r in rows if r["arm"] == a) for a in ("positive", "negative")}
    return {
        "result": "satmut_vista_plan",
        "pre_registration": reg.PRE_REGISTRATION,
        "amendment": reg.AMENDMENT,
        "windows": len(rows),
        "windows_by_arm": arms,
        "elements": len({r["element"] for r in rows}),
        "windows_without_phylop": sum(1 for r in rows if r["phylop_mean"] is None),
        "constraint_cost": {"requests": cost["requests"], "mb": round(cost["mb"], 2)},
        "seconds": round(time.time() - t0, 1),
    }


def read(rows: list[dict]) -> dict:
    """The registered reading: mean within-element Spearman per arm, and the difference."""
    by_element: dict[str, list[dict]] = {}
    for r in rows:
        if r.get("effect") is not None and r.get("phylop_mean") is not None:
            by_element.setdefault(r["element"], []).append(r)
    per_element = {}
    for eid, got in by_element.items():
        rho = spearman([x["effect"] for x in got], [x["phylop_mean"] for x in got])
        if rho is not None:
            per_element[eid] = {"arm": got[0]["arm"], "windows": len(got), "spearman": round(rho, 4)}
    arms = {
        a: [v["spearman"] for v in per_element.values() if v["arm"] == a] for a in ("positive", "negative")
    }
    means = {a: (round(mean(v), 4) if v else None) for a, v in arms.items()}
    diff = (
        round(means["positive"] - means["negative"], 4)
        if means["positive"] is not None and means["negative"] is not None
        else None
    )
    return {
        "elements_read": len(per_element),
        "elements_by_arm": {a: len(v) for a, v in arms.items()},
        "mean_spearman": means,
        "difference_positive_minus_negative": diff,
        "per_element": per_element,
    }


def shuffled_null(rows: list[dict], observed: float | None, draws: int = reg.DRAWS) -> dict:
    """The registered null: shuffle the constraint values WITHIN each element, keep both distributions."""
    rng = random.Random(23)
    by_element: dict[str, list[dict]] = {}
    for r in rows:
        if r.get("effect") is not None and r.get("phylop_mean") is not None:
            by_element.setdefault(r["element"], []).append(r)
    diffs = []
    for _ in range(draws):
        drawn = []
        for got in by_element.values():
            phy = [x["phylop_mean"] for x in got]
            rng.shuffle(phy)
            drawn.extend({**x, "phylop_mean": p} for x, p in zip(got, phy, strict=True))
        d = read(drawn)["difference_positive_minus_negative"]
        if d is not None:
            diffs.append(d)
    above = sum(1 for d in diffs if observed is not None and d >= observed)
    return {
        "draws": len(diffs),
        "mean": round(mean(diffs), 4) if diffs else None,
        "at_or_above_observed": above,
        "p_one_sided": round((above + 1) / (len(diffs) + 1), 4) if diffs else None,
    }


def verdict(diff: float | None, above: int, draws: int) -> str:
    """The registered bars, applied without interpretation."""
    if diff is None:
        return "not readable: an arm has no element with both readings"
    p = (above + 1) / (draws + 1) if draws else 1.0
    if diff >= 0.10 and p <= 0.01:
        return f"success: +{diff} at p {p:.4f}"
    if diff >= 0.05:
        return f"weak: +{diff} at p {p:.4f}, claims nothing"
    return f"failure: {diff:+} is below +0.05 - the arms are related the same way, as section 15 found"


def run(args) -> dict:
    t0 = time.time()
    rows = json.loads(PLAN.read_text())
    order = sorted(rows, key=lambda r: (r["chrom"], r["start"]))
    budget = args.max_requests
    adapter = AlphaGenomeAdapter()
    scorer = adapter._live_scorer(threshold=0.0)  # noqa: SLF001  (as scripts/enhancer_targets_all.py does)
    done: list[dict] = []
    stopped = ""
    with the_key(f"VISTA in-silico mutagenesis, budget {budget}"):
        ctx_of: dict[str, Context] = {}
        for row in order:
            if len(done) >= budget or stopped:
                break
            chrom = row["chrom"]
            if chrom not in ctx_of:
                for c in list(ctx_of):
                    ctx_of.pop(c).close()
                try:
                    ctx_of[chrom] = Context(chrom)
                except FileNotFoundError:
                    say(f"{chrom}: not fetched, its windows are skipped and counted")
                    continue
            ctx = ctx_of[chrom]
            wid = f"satmut_{chrom}_{row['start']}_{row['end']}"
            hit = score_element(
                scorer, ctx.genome.fetch, chrom, wid, row["start"], row["end"], coding=ctx.coding
            )
            pred = hit.get("predicted") or {}
            log2 = pred.get("log2_fold_change")
            done.append({**row, "effect": abs(log2) if log2 is not None else 0.0, "id": wid})
            if len(done) % 100 == 0:
                say(f"{len(done)}/{budget} windows scored")
            if len(done) == budget // 2:
                mid = read(done)["difference_positive_minus_negative"]
                say(f"futility look at {len(done)}: difference {mid}")
                if mid is not None and mid < 0.0:
                    stopped = f"futility at {len(done)} windows: difference {mid} below zero"
        for c in list(ctx_of):
            ctx_of.pop(c).close()

    reading = read(done)
    diff = reading["difference_positive_minus_negative"]
    null = shuffled_null(done, diff)
    above = null.get("at_or_above_observed") or 0
    return {
        "result": "satmut_vista",
        "pre_registration": reg.PRE_REGISTRATION,
        "amendment": reg.AMENDMENT,
        "windows_scored": len(done),
        "stopped": stopped,
        **{k: v for k, v in reading.items() if k != "per_element"},
        "null": null,
        "verdict": verdict(diff, above, null.get("draws") or 0),
        "per_element": reading["per_element"],
        "seconds": round(time.time() - t0, 1),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--quota-handed", action="store_true")
    ap.add_argument("--max-requests", type=int, default=2 * reg.ELEMENTS_PER_ARM * reg.WINDOWS_PER_ELEMENT)
    args = ap.parse_args(argv)

    if not args.run:
        out = assemble()
        say(f"saved {save_result(out['result'], out)}")
        say(f"{out['windows']} windows over {out['elements']} elements: {out['windows_by_arm']}")
        return 0
    if not args.quota_handed:
        raise SystemExit("part two spends model requests: pass --quota-handed once the key is yours")
    out = run(args)
    say(f"saved {save_result(out['result'], out)}")
    say(out["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
