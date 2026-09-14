# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does a storage unit execute its value? Part one assembles; part two spends the model quota.

    uv run python scripts/executor_test.py    # shortlist, controls, criterion, plan: no model call
    uv run python scripts/executor_test.py --run --quota-handed [--max-requests N]

Part one writes data/results/executor_shortlist.json and the full tables under
data/knowledge/human_panel/executor. Part two refuses to start without --quota-handed, spends requests
in the pre-registered order through the existing AlphaGenome variant path, and writes
data/results/executor_test.json with requests spent against units answered.
"""

from __future__ import annotations

import argparse
import gzip
import json
import time

from genomeos.attribution import executor as ex
from genomeos.results import save_result


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="chr21,chr22")
    ap.add_argument("--run", action="store_true", help="spend model requests (part two)")
    ap.add_argument("--quota-handed", action="store_true", help="confirm the quota holder has handed it over")
    ap.add_argument("--max-requests", type=int, default=None)
    args = ap.parse_args(argv)
    t0 = time.time()

    def say(*a) -> None:
        print(f"[{time.time() - t0:6.0f} s]", *a, flush=True)

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
    out = ex.run_pairs(pairs, ex.live_scorer(), max_requests=args.max_requests, progress=say)
    rows = out.pop("rows")
    with gzip.open(ex.LANE / "run_rows.json.gz", "wt") as fh:
        json.dump(rows, fh)
    out["criterion"] = ex.CRITERION
    out["seconds"] = round(time.time() - t0, 1)
    say(f"saved {save_result('executor_test', out)}")


if __name__ == "__main__":
    main()
