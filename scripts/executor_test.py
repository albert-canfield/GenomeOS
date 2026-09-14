# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does a storage unit execute its value? Part one assembles; part two spends the model quota.

    uv run python scripts/executor_test.py    # shortlist, controls, criterion, plan: no model call
    uv run python scripts/executor_test.py --run --quota-handed [--max-requests N]
    uv run python scripts/executor_test.py --wide             # the widened E1's plan: no model call
    uv run python scripts/executor_test.py --wide --run --quota-handed --max-pairs 1000

Part one writes data/results/executor_shortlist.json and the full tables under
data/knowledge/human_panel/executor. Part two refuses to start without --quota-handed, spends requests
in the pre-registered order through the existing AlphaGenome variant path, and writes
data/results/executor_test.json with requests spent against units answered.

--wide is the same two parts for the widened E1: every significant MPRAVarDB allele pair genome-wide
against a measured null from the same library (executor_mpra_wide.json, then
executor_mpra_wide_run.json).
"""

from __future__ import annotations

import argparse
import gzip
import json
import time

from genomeos.attribution import executor as ex
from genomeos.results import save_result


def wide(args, say, t0: float) -> None:
    """The widened E1, in the same two parts: assemble without a request, then spend the quota."""
    table = ex.LANE / "wide_assembled.json.gz"
    if not args.run:
        rows = ex.mpra_genome_wide(progress=say)
        say(f"{len(rows)} MPRAVarDB rows; building tests, read-out genes and measured-null controls")
        built = ex.wide_pairs(rows, progress=say)
        ex.LANE.mkdir(parents=True, exist_ok=True)
        with gzip.open(table, "wt") as fh:
            json.dump(built, fh, default=str)
        say("checking the direction convention against GTEx lymphocyte eQTLs: no model request")
        summary = ex.wide_summary(built, rows, ex.wide_sign_check(built["pairs"], progress=say))
        summary["seconds"] = round(time.time() - t0, 1)
        say(f"saved {save_result('executor_mpra_wide', summary)}")
        return
    if not args.quota_handed:
        raise SystemExit("part two spends AlphaGenome requests: pass --quota-handed once the quota is yours")
    with gzip.open(table, "rt") as fh:
        pairs = json.load(fh)["pairs"]
    wanted = set(args.endpoints.split(",")) & set(ex.WIDE_ENDPOINTS) or {"E1W_mpra_wide"}
    pairs = [p for p in pairs if p["endpoint"] in wanted]
    out = ex.run_pairs(
        pairs,
        ex.live_scorer(),
        max_requests=args.max_requests,
        progress=say,
        max_pairs=args.max_pairs,
        looks_at=tuple(sorted(wanted)),
        endpoints=ex.WIDE_ENDPOINTS,
        split="cell",
    )
    rows = out.pop("rows")
    out["no_call_grid"] = {e: [ex.recall(rows, e, t) for t in (0.0, 0.001, 0.01)] for e in sorted(wanted)}
    with gzip.open(ex.LANE / f"wide_rows_{'-'.join(sorted(wanted))}.json.gz", "wt") as fh:
        json.dump(rows, fh)
    primary = ex.LANE / "run_rows_E1_mpra-E2_eqtl.json.gz"
    if primary.exists():
        with gzip.open(primary, "rt") as fh:
            out["against_e2"] = ex.compare_endpoints(json.load(fh), rows, b=sorted(wanted)[0])
    out["pre_registration"] = ex.E1_WIDE
    out["seconds"] = round(time.time() - t0, 1)
    say(f"saved {save_result('executor_mpra_wide_run', out)}")


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="chr21,chr22")
    ap.add_argument("--run", action="store_true", help="spend model requests (part two)")
    ap.add_argument("--quota-handed", action="store_true", help="confirm the quota holder has handed it over")
    ap.add_argument("--max-requests", type=int, default=None)
    ap.add_argument("--max-pairs", type=int, default=None, help="stop at this many pairs as well")
    ap.add_argument(
        "--endpoints", default="E1_mpra,E2_eqtl", help="E3 runs only when the quota holder agrees"
    )
    ap.add_argument(
        "--wide",
        action="store_true",
        help="the widened E1: every significant MPRAVarDB allele pair genome-wide, against measured nulls",
    )
    args = ap.parse_args(argv)
    t0 = time.time()

    def say(*a) -> None:
        print(f"[{time.time() - t0:6.0f} s]", *a, flush=True)

    if args.wide:
        wide(args, say, t0)
        return
    table = ex.LANE / "assembled.json.gz"
    if not args.run:
        assembled = ex.assemble(args.chroms.split(","), progress=say)
        ex.LANE.mkdir(parents=True, exist_ok=True)
        with gzip.open(table, "wt") as fh:
            json.dump(assembled, fh, default=str)
        summary = ex.summarise(assembled)
        summary["seconds"] = round(time.time() - t0, 1)
        say(f"saved {save_result('executor_shortlist', summary)}")
        return
    if not args.quota_handed:
        raise SystemExit("part two spends AlphaGenome requests: pass --quota-handed once the quota is yours")
    with gzip.open(table, "rt") as fh:
        pairs = json.load(fh)["matched"]["pairs"]
    wanted = set(args.endpoints.split(","))
    pairs = [p for p in pairs if p["endpoint"] in wanted]
    out = ex.run_pairs(
        pairs,
        ex.live_scorer(),
        max_requests=args.max_requests,
        progress=say,
        max_pairs=args.max_pairs,
    )
    rows = out.pop("rows")
    out["no_call_grid"] = {e: [ex.recall(rows, e, t) for t in (0.0, 0.001, 0.01)] for e in sorted(wanted)}
    with gzip.open(ex.LANE / f"run_rows_{'-'.join(sorted(wanted))}.json.gz", "wt") as fh:
        json.dump(rows, fh)
    primary = ex.LANE / "run_rows_E1_mpra-E2_eqtl.json.gz"
    if "E3_eqtl_linked" in wanted and primary.exists():
        with gzip.open(primary, "rt") as fh:
            out["hold_out_read"] = ex.compare_endpoints(json.load(fh), rows)
    out["criterion"] = ex.CRITERION
    out["amendments"] = list(ex.AMENDMENTS)
    out["hold_out"] = ex.HOLD_OUT
    name = "executor_test" if wanted == {"E1_mpra", "E2_eqtl"} else "executor_holdout"
    out["seconds"] = round(time.time() - t0, 1)
    say(f"saved {save_result(name, out)}")


if __name__ == "__main__":
    main()
