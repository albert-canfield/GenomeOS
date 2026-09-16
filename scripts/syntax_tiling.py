# SPDX-License-Identifier: AGPL-3.0-or-later
"""The 69 syntax blocks tiled and deleted against two matched arms, as pre-registered on 2026-09-17.

    uv run python scripts/syntax_tiling.py                       # assemble and match: no model request
    uv run python scripts/syntax_tiling.py --run --quota-handed  # spend the budget

The design is `attribution/syntax_tiling.PRE_REGISTRATION`, committed before this script existed
(9b2a4ea). Nothing here may change it: this file assembles the windows the registration describes,
scores them through the sweep's own scorer and cache, and reads the result by the registered bars.

Part one writes `syntax_tiling_plan` (how many windows each arm has, how many syntax windows found a
match in each arm, and the requests the run would cost). Part two refuses to start without
--quota-handed, holds the shared model key for the length of the run, spends one request per window in
a fixed order, looks once at half the budget for futility only, and writes `syntax_tiling`.
"""

from __future__ import annotations

import argparse
import json
import time
from contextlib import contextmanager
from pathlib import Path

from genomeos import jobs
from genomeos.attribution import syntax_tiling as st
from genomeos.predict import AlphaGenomeAdapter
from genomeos.predict.enhancer_target import Context, score_element
from genomeos.results import save_result

HOLDER = "genomeos-9c"
PLAN = Path("data/knowledge/alphagenome/syntax_tiling_windows.json")


def say(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


@contextmanager
def the_key(what: str):
    jobs.take_key(HOLDER, what)
    try:
        yield
    finally:
        jobs.drop_key(HOLDER)


def assemble(chroms: list[str]) -> dict:
    """Every arm's windows, matched, with no request made."""
    t0 = time.time()
    rows: list[dict] = []
    for chrom in chroms:
        try:
            ctx = Context(chrom)
        except FileNotFoundError:
            say(f"{chrom}: not fetched, skipped")
            continue
        coding_tss = sorted(
            (g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start)
            for g in ctx.annotation.protein_coding()
            if g.locus.chrom == chrom
        )
        spans = st.scored_starts(chrom)
        got = st.windows_of(chrom, ctx, coding_tss, spans)
        kept = st.match_windows(got)
        ctx.close()
        rows.extend(kept)
        n_syntax = sum(1 for r in kept if r["arm"] == "syntax")
        if n_syntax:
            say(f"{chrom}: {n_syntax} syntax windows kept, {len(kept)} rows with their matches")
    PLAN.parent.mkdir(parents=True, exist_ok=True)
    PLAN.write_text(json.dumps(rows))
    syntax = [r for r in rows if r["arm"] == "syntax"]
    return {
        "result": "syntax_tiling_plan",
        "pre_registration": st.PRE_REGISTRATION,
        "windows_by_arm": {a: sum(1 for r in rows if r["arm"] == a) for a in st.ARMS},
        "syntax_windows": len(syntax),
        "syntax_windows_matched_in_both": sum(1 for r in syntax if len(r.get("matched", [])) == 2),
        "syntax_windows_matched_in_neither": sum(1 for r in syntax if not r.get("matched")),
        "syntax_windows_already_scored": sum(1 for r in syntax if r["already_scored"]),
        "requests_if_run": len(rows),
        "blocks_by_arm": {a: len({r["block"] for r in rows if r["arm"] == a}) for a in st.ARMS},
        "seconds": round(time.time() - t0, 1),
    }


def read(rows: list[dict]) -> dict:
    """The registered reading: the two primaries, the secondary, and the per-arm counts behind them."""
    by_arm = {a: [r for r in rows if r["arm"] == a] for a in st.ARMS}
    hits = {a: sum(1 for r in v if st.moved(r)) for a, v in by_arm.items()}
    n = {a: len(v) for a, v in by_arm.items()}
    syntax = by_arm["syntax"]
    old = [r for r in syntax if r["already_scored"]]
    new = [r for r in syntax if not r["already_scored"]]
    return {
        "by_arm": {
            a: {"windows": n[a], "moved": hits[a], "rate": round(hits[a] / n[a], 4) if n[a] else None}
            for a in st.ARMS
        },
        "syntax_minus_neutral": st.difference(hits["syntax"], n["syntax"], hits["neutral"], n["neutral"]),
        "syntax_minus_relaxed": st.difference(hits["syntax"], n["syntax"], hits["relaxed"], n["relaxed"]),
        "secondary_already_scored_against_new": st.difference(
            sum(1 for r in old if st.moved(r)), len(old), sum(1 for r in new if st.moved(r)), len(new)
        ),
    }


def verdict(reading: dict) -> str:
    """The registered bars, applied without interpretation."""
    a, b = reading["syntax_minus_neutral"], reading["syntax_minus_relaxed"]
    if a["difference"] is None or b["difference"] is None:
        return "not readable: an arm has no windows"
    strong = [x for x in (a, b) if x["difference"] >= 0.10 and x["p_one_sided"] <= 0.01]
    if len(strong) == 2:
        return "success: at least +0.10 against both arms at p 0.01"
    if a["difference"] < 0.05:
        return "failure: below +0.05 against the neutral arm"
    if len(strong) == 1 or (a["difference"] >= 0.05 and b["difference"] >= 0.05):
        return "weak: the observation survives, and nothing is claimed for the blocks"
    return "failure: below +0.05 against the neutral arm"


def run(args) -> dict:
    t0 = time.time()
    rows = json.loads(PLAN.read_text())
    order = sorted(rows, key=lambda r: (r["chrom"], r["start"], r["arm"]))
    budget = args.max_requests
    look_at = budget // 2
    adapter = AlphaGenomeAdapter()
    scorer = adapter._live_scorer(threshold=0.0)  # noqa: SLF001  (as scripts/enhancer_targets_all.py does)
    done: list[dict] = []
    stopped = ""
    with the_key(f"the 69 syntax blocks tiled, budget {budget}"):
        ctx_of: dict[str, Context] = {}
        for row in order:
            if len(done) >= budget or stopped:
                break
            chrom = row["chrom"]
            if chrom not in ctx_of:
                for c in list(ctx_of):
                    ctx_of.pop(c).close()
                ctx_of[chrom] = Context(chrom)
            ctx = ctx_of[chrom]
            wid = f"tile_{chrom}_{row['start']}_{row['end']}"
            hit = score_element(
                scorer,
                ctx.genome.fetch,
                chrom,
                wid,
                row["start"],
                row["end"],
                coding=ctx.coding,
            )
            done.append({**row, "predicted": hit.get("predicted"), "id": wid})
            if len(done) % 25 == 0:
                say(f"{len(done)}/{budget} windows scored")
            if len(done) == look_at:
                mid = read(done)
                upper = mid["syntax_minus_neutral"]["upper_95"]
                say(f"futility look at {look_at}: syntax - neutral upper 95% = {upper}")
                if upper is not None and upper < 0.05:
                    stopped = f"futility at {look_at} windows: upper 95% {upper} below 0.05"
        for c in list(ctx_of):
            ctx_of.pop(c).close()
    reading = read(done)
    return {
        "result": "syntax_tiling",
        "pre_registration": st.PRE_REGISTRATION,
        "windows_scored": len(done),
        "stopped": stopped,
        **reading,
        "verdict": verdict(reading),
        "seconds": round(time.time() - t0, 1),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="", help="comma-separated; default every chromosome")
    ap.add_argument("--run", action="store_true", help="spend the budget (part two)")
    ap.add_argument("--quota-handed", action="store_true", help="confirm the model key is yours")
    ap.add_argument("--max-requests", type=int, default=2400, help="the registered budget")
    args = ap.parse_args(argv)

    if not args.run:
        chroms = args.chroms.split(",") if args.chroms else [f"chr{c}" for c in [*range(1, 23), "X"]]
        out = assemble(chroms)
        say(f"saved {save_result(out['result'], out)}")
        say(f"{out['syntax_windows']} syntax windows, {out['requests_if_run']} requests if run")
        return 0
    if not args.quota_handed:
        raise SystemExit("part two spends model requests: pass --quota-handed once the key is yours")
    out = run(args)
    say(f"saved {save_result(out['result'], out)}")
    say(out["verdict"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
