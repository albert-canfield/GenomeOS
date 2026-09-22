# SPDX-License-Identifier: AGPL-3.0-or-later
"""Sort every miss in every scored frame into the classes registered before they were counted.

    uv run python scripts/loci_miss.py --registration      # what was committed first
    uv run python scripts/loci_miss.py --frame NAME        # one frame's distribution
    uv run python scripts/loci_miss.py --all [--save]      # every frame, pooled, and the verdict

`genomeos/benchmark/loci_miss.py` carries the registration: the classes, the order they are tested
in, the threshold that would make the case for a second overlap-tolerant rate reported beside the
strict one, the threshold that would close the question, and the trap - a benchmark that widens what
counts as a hit after seeing its rate is a benchmark tuning itself. It spends nothing: every reading
is on disk and GENCODE is a local annotation.
"""

from __future__ import annotations

import argparse
import json
import sys

from genomeos.benchmark import loci_miss as lm
from genomeos.results import save_result


def say(msg: str) -> None:
    print(msg, file=sys.stderr)


def show(frame: dict) -> None:
    h = frame["headline"]
    print(f"\n{frame['frame']}  (docs/LOCI-BENCHMARK.md section {frame['section']})")
    print(f"  {frame['loci']} loci, {frame['hits']} strict hits, {frame['misses']} strict misses")
    print(f"  headline over {h['n']} ({h['denominator']}):")
    print(f"    strict                   {h['strict_k']}/{h['n']} ({h['strict_rate']})")
    print(f"    overlap-tolerant, A DESCRIPTION  {h['tolerant_k']}/{h['n']} ({h['tolerant_rate']})")
    if h["loci_it_would_add"]:
        print(f"    it would add: {', '.join(h['loci_it_would_add'])}")
    print(
        f"    the nearest-TSS-in-node baseline, both readings: {h['baseline_strict_k']}/{h['n']}"
        f" -> {h['baseline_tolerant_k']}/{h['n']}"
    )
    print("  the misses, by class (any derived layer / the deletion layer alone):")
    for name, n in frame["miss_classes"].items():
        print(f"    {name:26s} {n:3d}   {frame['miss_classes_deletion_layer_only'][name]:3d}")
    for o in frame["overlapping_misses"]:
        print(f"    overlap: {o['locus']} -> {o['named']} ({o['layer']}) over {o['targets']}")
    if frame["control_hits_that_are_not_exact"]:
        print(f"  CONTROL FAILED: {frame['control_hits_that_are_not_exact']}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--registration", action="store_true")
    ap.add_argument("--frame")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--save", action="store_true", help="with --all, write data/results/loci_miss.json")
    args = ap.parse_args()

    if args.registration:
        print(json.dumps({"classes": lm.CLASSES, "preregistration": lm.PREREGISTRATION}, indent=1))
        return 0
    if args.frame:
        show(lm.classify_frame(args.frame, lm.rows_for(args.frame, lm.RESULTS_DIR)))
        return 0
    if args.all:
        out = lm.classify_all()
        for v in out["frames"].values():
            show(v)
        p = out["pooled"]
        print(f"\npooled over {p['misses']} misses:")
        for name, n in p["classes"].items():
            print(f"  {name:26s} {n:3d}   {p['classes_deletion_layer_only'][name]:3d}")
        print(f"  overlap share of misses: {p['overlap_share_of_misses']}")
        print(f"  with one overlap fewer it would read: {p['one_fewer_overlap_would_read']}")
        for f, s in p["overlap_share_dropping_one_fourth_frame_reading"].items():
            print(f"  the same share with {f} dropped from the pool: {s}")
        gap = p["largest_tolerant_rate_gap_at_n_over_50"]
        print(f"  largest tolerant-minus-strict rate gap at n >= 50: {gap}")
        print(f"  tolerant minus strict, per frame: {p['tolerant_minus_strict_per_frame']}")
        v = out["verdict"]
        print(f"  case for a second rate: {v['case_for_a_second_rate']}")
        print(f"  question closed: {v['question_closed']}")
        print(f"requests spent: {out['requests_spent']}")
        if args.save:
            save_result("loci_miss", out)
            say("wrote data/results/loci_miss.json")
        return 0
    ap.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
