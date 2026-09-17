# SPDX-License-Identifier: AGPL-3.0-or-later
"""Re-run the variation lane on the chromosomes still carrying the old base counts, one at a time.

    uv run python scripts/variation_rerun_chain.py                 # every chromosome but chr21
    uv run python scripts/variation_rerun_chain.py --chrom chr22    # one of them

`measured_bp` and `human_constrained_bp` used to count every kilobase a Gnocchi bin gave an
interval; since 2026-09-17 they count the interval's own bases, and the old quantity is kept as
`touched_bin_bp`. chr21 was re-run with the fix and its deltas published (GRAMMAR-BY-COMPARISON.md,
"The kilobase that was called a base"). This chain brings the other 23 into the same meaning, so no
two files disagree about what the field is.

The stopping rule is fixed here before the chain runs, at genomeos-9c's request. On chr21 the
constrained fraction moved by at most 0.0073 absolute in a tier of any size, which is edge-bin
weighting: the bins at an interval's edges are now weighted by how much of them the interval covers.
A larger move would mean something other than the crediting changed, so the chain STOPS on the first
chromosome where a tier holding at least MIN_TIER_BP of measured bases moves its fraction by more
than FRACTION_TOLERANCE absolute, and says which tier, rather than finishing the sweep and burying
it. Tiers below that size are reported and never stop the chain, because a handful of blocks can
move a ratio by any amount (chr21's structural tier is three blocks: 0.5 to 0.5839).

Every re-run file keeps `touched_bin_bp` populated, so an old file and a new one can be compared
without either quantity going missing. Saves data/results/variation_rerun_chain.json.
"""

from __future__ import annotations

import argparse
import json
import shutil
import time
from pathlib import Path
from typing import Any

from genomeos.attribution import variation
from genomeos.results import RESULTS_DIR, save_result

FRACTION_TOLERANCE = 0.02  # absolute move in a tier's constrained fraction that stops the chain
MIN_TIER_BP = 100_000  # a tier smaller than this cannot stop the chain; its ratio is noise
CHROMS = tuple(f"chr{c}" for c in list(range(1, 23)) + ["X", "Y"])
FIELDS = ("measured_bp", "human_constrained_bp", "touched_bin_bp", "human_constrained_fraction")


def walk_cost(cost: Any) -> tuple[int, float]:
    """Requests and megabytes summed over a nested cost dict (one entry per read the lane made)."""
    requests, mb = 0, 0.0
    if isinstance(cost, dict):
        requests += int(cost.get("requests") or 0)
        mb += float(cost.get("mb_fetched") or 0)
        for v in cost.values():
            if isinstance(v, dict):
                r, m = walk_cost(v)
                requests += r
                mb += m
    return requests, mb


def before_after(before: dict[str, Any], after: dict[str, Any]) -> dict[str, Any]:
    """Per tier, the two absolute counts and the fraction, before and after, with the verdict."""
    out: dict[str, Any] = {"tiers": {}, "stopped_by": None}
    for tier, b in (before.get("by_tier") or {}).items():
        a = (after.get("by_tier") or {}).get(tier)
        if not a:
            out["tiers"][tier] = {"missing_after": True}
            continue
        fb, fa = b.get("human_constrained_fraction"), a.get("human_constrained_fraction")
        move = None if fb is None or fa is None else round(abs(fa - fb), 4)
        row = {
            "measured_bp": [b.get("measured_bp"), a.get("measured_bp")],
            "human_constrained_bp": [b.get("human_constrained_bp"), a.get("human_constrained_bp")],
            "touched_bin_bp": a.get("touched_bin_bp"),
            "fraction": [fb, fa],
            "fraction_move": move,
            "big_enough_to_judge": (a.get("measured_bp") or 0) >= MIN_TIER_BP,
        }
        if move is not None and row["big_enough_to_judge"] and move > FRACTION_TOLERANCE:
            out["stopped_by"] = tier
        out["tiers"][tier] = row
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", help="one chromosome instead of the whole chain")
    ap.add_argument("--results", default=str(RESULTS_DIR))
    args = ap.parse_args()
    results = Path(args.results)
    todo = [args.chrom] if args.chrom else [c for c in CHROMS if c != "chr21"]
    done: dict[str, Any] = {}
    requests = mb = 0.0
    t0 = time.time()
    stopped = None
    for chrom in todo:
        path = results / f"variation_{chrom}.json"
        if not path.exists():
            done[chrom] = {"refused": "no variation result for this chromosome yet"}
            print(f"{chrom}: refused, no earlier result", flush=True)
            continue
        before = json.loads(path.read_text())
        keep = path.with_suffix(".before.json")
        shutil.copyfile(path, keep)
        t1 = time.time()
        after = variation.run_and_save(chrom)
        cost = (after.get("cost") or {}) if isinstance(after, dict) else {}
        r, m = walk_cost(cost)  # the lane reports a cost per read, not one total
        requests += r
        mb += m
        cmp = before_after(before, after)
        cmp["seconds"] = round(time.time() - t1, 1)
        cmp["cost"] = cost
        cmp["range_reads"] = {"requests": r, "mb_fetched": round(m, 2)}
        done[chrom] = cmp
        keep.unlink(missing_ok=True)
        worst = max((r.get("fraction_move") or 0) for r in cmp["tiers"].values() if isinstance(r, dict))
        print(
            f"{chrom}: re-run in {cmp['seconds']}s, largest fraction move {worst}, "
            f"stopped_by {cmp['stopped_by']}",
            flush=True,
        )
        if cmp["stopped_by"]:
            stopped = {"chrom": chrom, "tier": cmp["stopped_by"], "tiers": cmp["tiers"]}
            print(
                f"STOP: {chrom}'s {cmp['stopped_by']} tier moved its constrained fraction by more "
                f"than {FRACTION_TOLERANCE}; the chain stops here as pre-registered",
                flush=True,
            )
            break
    payload = {
        "question": (
            "bring every chromosome's measured_bp and human_constrained_bp onto the interval's own "
            "bases, as chr21 already is, without any two files disagreeing about the field"
        ),
        "stopping_rule": (
            f"stop on the first chromosome where a tier with at least {MIN_TIER_BP} measured bases "
            f"moves its constrained fraction by more than {FRACTION_TOLERANCE} absolute; fixed "
            "before the chain ran"
        ),
        "chromosomes_attempted": todo,
        "stopped": stopped,
        "range_reads": {"requests": int(requests), "mb_fetched": round(mb, 2)},
        "seconds": round(time.time() - t0, 1),
        "per_chromosome": done,
    }
    path = save_result("variation_rerun_chain", payload)
    print(
        f"\n{len(done)} chromosomes, {int(requests)} requests, {round(mb, 2)} MB, "
        f"{payload['seconds']}s; stopped: {stopped['chrom'] if stopped else 'no'}\n-> {path}"
    )
    return 1 if stopped else 0


if __name__ == "__main__":
    raise SystemExit(main())
