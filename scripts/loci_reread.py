# SPDX-License-Identifier: AGPL-3.0-or-later
"""Re-read every scored frame through the repaired deletion reader. No request, ever.

    uv run python scripts/loci_reread.py --registration   # what was committed before any re-read
    uv run python scripts/loci_reread.py --frame NAME     # one frame, before and after
    uv run python scripts/loci_reread.py --all            # every frame, and save the result

`genomeos/benchmark/loci_reread.py` carries the registration: the reading, the direction every
affected rate was expected to move in, and the falsifiers, all committed before the one-line change
to `loci.read_deletion` was made. The re-read spends nothing - it ranks rows the finished
all-element sweep had already bought - and it refuses to report a number unless the pre-2026-09-21
loop reproduces every stored reading first.
"""

from __future__ import annotations

import argparse
import json
import sys

from genomeos.benchmark import loci_reread as rr
from genomeos.results import save_result


def say(msg: str) -> None:
    print(msg, file=sys.stderr)


def show(frame: dict) -> None:
    b, a = frame["before"], frame["after"]
    print(f"\n{frame['frame']}  (docs/LOCI-BENCHMARK.md section {frame['section']})")
    print(
        f"  fidelity: {frame['fidelity']['loci_checked']} loci and"
        f" {frame['fidelity']['negatives_checked']} negatives replayed, 0 mismatches"
    )
    for key in (
        "target_derived",
        "target_derived_where_the_model_could_answer",
        "target_heuristic",
        "target_looked_up",
    ):
        x, y = b.get(key) or {}, a.get(key) or {}
        arrow = "" if x.get("k") == y.get("k") else "   <-- moved"
        print(
            f"  {key:52s} {x.get('k')}/{x.get('n')} ({x.get('rate')})"
            f" -> {y.get('k')}/{y.get('n')} ({y.get('rate')}){arrow}"
        )
    print(f"  deletion layer rank-1 changed at {frame['rank_one_changed']} of {frame['loci']} loci")
    for m in frame["movements"]:
        flip = (
            "hit->miss"
            if m["hit_before"] and not m["hit_after"]
            else ("miss->hit" if m["hit_after"] and not m["hit_before"] else "no flip")
        )
        print(
            f"    {m['locus']:38s} {m['was']} ({m['was_log2']}) -> {m['now']} ({m['now_log2']})"
            f"  published {m['published']}  {flip}"
        )
    neg = frame["negatives"]
    print(
        f"  negatives: {neg['checked']} checked, {len(neg['claims_changed'])} claim(s) moved,"
        f" {len(neg['element_counts_changed'])} element count(s) moved"
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--registration", action="store_true")
    ap.add_argument("--frame")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--save", action="store_true", help="with --all, write data/results/loci_reread.json")
    args = ap.parse_args()

    if args.registration:
        print(json.dumps(rr.PREREGISTRATION, indent=1))
        return 0
    if args.frame:
        show(rr.reread_frame(args.frame))
        return 0
    if args.all:
        out = rr.reread_all()
        for name in out["frames"]:
            show(out["frames"][name])
        print(f"\nrequests spent: {out['requests_spent']}")
        if args.save:
            save_result("loci_reread", out)
            say("wrote data/results/loci_reread.json")
        return 0
    ap.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
