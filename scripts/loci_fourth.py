# SPDX-License-Identifier: AGPL-3.0-or-later
"""Draw and grade a fourth frame of published enhancer-gene loci, by rule rather than by eye.

    uv run python scripts/loci_fourth.py --draw                # the rule, free, no request
    uv run python scripts/loci_fourth.py --plan                # reach and cost, free, no request
    uv run python scripts/loci_fourth.py --intervals --run --quota-handed
    uv run python scripts/loci_fourth.py --score               # read and score, no request

The roadmap item is "the known-locus benchmark with more loci, not more readings of the one it has",
and section 21 says which loci are worth having: the three existing frames agree on the derived rate
to 0.007 while the nearest-gene rule over them moves 0.445, because the frames differ in how far each
element sits from its target. So the fourth frame fixes that geometry with a rule stated before the
loci are known - every held-out ENCODE CRISPR positive whose target is NOT its own nearest coding TSS
- and `genomeos/benchmark/loci_fourth.py` holds the registration and the rule. This file drives them.

Four steps, the free ones first:

  --draw       apply the rule to the held-out arm. Prints what it kept, what it rejected and why,
               before anything is scored and before a request exists.
  --plan       `loci.read_reach` per drawn locus (gene BODY against the scorer's 1 Mb input, never
               TSS distance), then whether the finished all-chromosome sweep has already deleted an
               element inside the perturbed interval - in which case the answer costs nothing.
  --intervals  deletes the perturbed elements no annotation covers, one request each, through the
               same `Context.score_region` the other three frames' stated intervals use.
  --score      `loci.build` over the graded loci and their matched negatives. No request.
"""

from __future__ import annotations

import argparse
import json
import time
from contextlib import contextmanager
from typing import Any

from genomeos import jobs
from genomeos.benchmark import loci_fourth as fourth
from genomeos.predict import AlphaGenomeAdapter, status
from genomeos.predict.enhancer_target import Context
from genomeos.results import save_result

HOLDER = "lane-loci4"


def say(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


@contextmanager
def the_key(what: str):
    jobs.take_key(HOLDER, what)
    try:
        yield
    finally:
        jobs.drop_key(HOLDER)


def show_draw(draw: dict[str, Any]) -> None:
    say(
        f"{draw['perturbed_elements_with_a_regulated_gene']} perturbed element(s) with a regulated gene"
        f" in the held-out arm; the rule keeps {draw['kept']}"
    )
    for reason, n in sorted(draw["rejected_by_reason"].items(), key=lambda kv: -kv[1]):
        say(f"  {n:>4}  {reason}")
    for e in draw["panel"]:
        say(
            f"{e.locus}: {e.chrom}:{e.element[0]}-{e.element[1]}, target {'/'.join(e.targets)} at"
            f" {e.distance} bp, nearest coding TSS {e.nearest_gene_trap}"
        )
    if draw["control"]:
        c = draw["control"]
        say(f"control: {c.locus}, target {'/'.join(c.targets)} at {c.distance} bp (the shortcut is right)")
    if draw["over_cap"]:
        say(f"over the cap of {fourth.MAX_LOCI} and not drawn: {', '.join(draw['over_cap'])}")


def show_plan(rows: list[dict[str, Any]]) -> None:
    for r in rows:
        if not r["askable"]:
            far = r["targets_out_of_reach"][0] if r["targets_out_of_reach"] else {}
            say(
                f"{r['locus']}: UNASKABLE - {far.get('target')} is {far.get('distance_bp', 0) // 1000} kb"
                " from the element and the model's input reaches 524 kb each way. No request"
            )
            continue
        say(
            f"{r['locus']}: graded, {r['length']} bp, target at {r['distance']} bp,"
            f" {r['ccres_over_element']} cCRE(s) over the element,"
            f" {r['already_scored_over_element']} already deleted, {r['requests']} request(s)"
        )
    say(
        f"{sum(1 for r in rows if r['graded'])} graded of {len(rows)} loci (one is the control);"
        f" {sum(1 for r in rows if not r['askable'])} died at the reach filter;"
        f" {sum(r['requests'] for r in rows)} request(s) to spend"
    )


def score_intervals(rows: list[dict[str, Any]], draw: dict[str, Any]) -> dict[str, Any]:
    """Delete the perturbed elements the sweep has not already covered, one request each."""
    t0 = time.time()
    asked = [r for r in rows if r["requests"]]
    adapter = AlphaGenomeAdapter()
    scorer = adapter._live_scorer(threshold=0.0)  # noqa: SLF001  (as scripts/loci_third.py does)
    elements: list[dict[str, Any]] = []
    read: list[dict[str, Any]] = []
    by_name = {e.locus: e for e in [*draw["panel"]] + ([draw["control"]] if draw["control"] else [])}
    with the_key(f"known-locus fourth frame, {len(asked)} perturbed intervals"):
        for row in asked:
            ctx = Context(row["chrom"])
            try:
                say(f"{row['locus']}: deleting {row['chrom']}:{row['element'][0]}-{row['element'][1]}")
                hit = ctx.score_region(scorer, row["element"][0], row["element"][1])
            finally:
                ctx.close()
            # score_region substitutes an overlapping ENCODE element when there is one, which would
            # answer about a different interval than the perturbed one. Recorded, never hidden.
            want = f"{row['chrom']}_{row['element'][0]}_{row['element'][1]}"
            substituted = hit.get("id") != want
            elements.append(
                {
                    **hit,
                    "chrom": row["chrom"],
                    "locus": row["locus"],
                    "substituted_annotated_element": substituted,
                }
            )
            pred = hit.get("predicted_coding") or hit.get("predicted") or {}
            any_gene = hit.get("predicted") or {}
            targets = list(by_name[row["locus"]].targets)
            read.append(
                {
                    "locus": row["locus"],
                    "element": f"{row['chrom']}:{row['element'][0]}-{row['element'][1]}",
                    "substituted_annotated_element": substituted,
                    "named_first_coding": pred.get("gene"),
                    "named_first_any": any_gene.get("gene"),
                    "log2_fold_change": pred.get("log2_fold_change"),
                    "action": pred.get("action"),
                    "nearest_coding_tss": by_name[row["locus"]].nearest_gene_trap,
                    "tissue": pred.get("tissue"),
                    "top_genes": [g["gene"] for g in (hit.get("top_genes") or []) if g.get("gene")][:6],
                    "hit": bool(pred.get("gene") in targets),
                    "targets": targets,
                }
            )
            say(f"{row['locus']}: named {pred.get('gene')} first among coding genes")
    return {
        "result": fourth.INTERVALS,
        "elements": elements,
        "read": read,
        "note": (
            "AlphaGenome deletions of the ENCODE CRISPR benchmark's own perturbed intervals, drawn"
            " into the fourth frame by genomeos/benchmark/loci_fourth.py and not already covered by"
            " the finished sweep. Read by loci._deletion_rows and labelled `stated_interval`, so these"
            " rows are never pooled with registry ones. loci.STATED_INTERVAL_RESULTS must list this"
            " result's name for that to happen: see loci_fourth.NEEDED_LOCI_HUNK and"
            " loci_fourth.register_stated_intervals"
        ),
        "requests_this_run": len(read),
        "seconds": round(time.time() - t0, 1),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--draw", action="store_true", help="apply the rule and print the frame; no request")
    ap.add_argument("--plan", action="store_true", help="reach and cost per drawn locus; no request")
    ap.add_argument("--intervals", action="store_true", help="delete the elements the sweep missed")
    ap.add_argument("--score", action="store_true", help="read and score the graded loci")
    ap.add_argument("--run", action="store_true", help="with --intervals: actually spend the requests")
    ap.add_argument("--quota-handed", action="store_true")
    ap.add_argument("--no-network", action="store_true", help="local layers only, for a quick check")
    args = ap.parse_args(argv)
    if not (args.draw or args.plan or args.intervals or args.score):
        args.draw = True

    draw = fourth.select(progress=say)
    if args.draw or args.plan or args.intervals:
        show_draw(draw)
    if len(draw["panel"]) < 5:
        say(f"the rule returned {len(draw['panel'])} loci; the registration says no rate is quoted below 5")

    rows = fourth.plan(draw) if (args.plan or args.intervals) else []
    if args.plan or args.intervals:
        show_plan(rows)

    if args.intervals:
        if not args.run:
            say("pass --run to spend the requests above")
            return 0
        if not args.quota_handed:
            raise SystemExit("this spends model requests: pass --quota-handed once the key is yours")
        st = status()
        if not st["enabled"]:
            raise SystemExit(f"AlphaGenome is disabled: {st['reason']}. {st['how']}")
        out = score_intervals(rows, draw)
        say(f"saved {save_result(fourth.INTERVALS, out)}")

    if args.score:
        out = fourth.run(network=not args.no_network, progress=say, draw=draw)
        a = out["aggregate"]
        say(
            f"target derived {a['target_derived']['k']}/{a['target_derived']['n']},"
            f" heuristic {a['target_heuristic']['k']}/{a['target_heuristic']['n']},"
            f" chance floor {a['target_by_chance']['expected']}"
        )
        print(json.dumps(out["reach"], indent=1))
        print(json.dumps({k: v for k, v in out["baselines"].items() if k != "per_locus"}, indent=1))
        print(json.dumps(out["positive_control"], indent=1))
        print(json.dumps(out["beside_the_other_three"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
