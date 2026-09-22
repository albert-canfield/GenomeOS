# SPDX-License-Identifier: AGPL-3.0-or-later
"""Grade eleven more published enhancer-gene loci against the machinery.

    uv run python scripts/loci_candidates.py --plan            # reach and cost, free, no request
    uv run python scripts/loci_candidates.py --intervals --run --quota-handed
    uv run python scripts/loci_candidates.py --score           # read and score, no request

The roadmap item is "the known-locus benchmark with more loci, not more readings of the one it has".
The candidates, the hit rules, the denominators, the controls and the predictions are registered in
`genomeos/benchmark/loci_candidates.py` before anything is scored; this file only drives them.

Three steps, in an order chosen so the free ones come first:

  --plan       `loci.read_reach` per candidate (gene BODY against the scorer's 1 Mb input, never TSS
               distance), then whether any annotation already drew an interval over the published
               element - in which case the finished sweep has deleted it already and the element-level
               answer costs nothing. Prints the request count before a request exists.
  --intervals  deletes the published elements no annotation covers, one request each, through the
               same `Context.score_region` the panel's stated intervals use. Skips anything the reach
               filter calls unaskable: a request there buys a certain negative.
  --score      `loci.build` over the graded candidates and their matched negatives. No request.
"""

from __future__ import annotations

import argparse
import json
import time
from contextlib import contextmanager
from typing import Any

from genomeos import jobs
from genomeos.benchmark import loci_candidates as cands
from genomeos.predict import AlphaGenomeAdapter, status
from genomeos.predict.enhancer_target import Context
from genomeos.results import save_result

HOLDER = "genomeos-ae"


def say(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


@contextmanager
def the_key(what: str):
    jobs.take_key(HOLDER, what)
    try:
        yield
    finally:
        jobs.drop_key(HOLDER)


def show_plan(rows: list[dict[str, Any]]) -> None:
    for r in rows:
        if not r["askable"]:
            far = r["targets_out_of_reach"][0]
            say(
                f"{r['locus']}: UNASKABLE - {far['target']} is {far.get('distance_bp', 0) // 1000} kb"
                " from the element and the model's input reaches 524 kb each way. No request"
            )
            continue
        tag = "exhibit, not graded" if r["exhibit_only"] else "graded"
        say(
            f"{r['locus']}: {tag}, targets in reach {r['targets_in_reach']},"
            f" {r['ccres_over_element']} cCRE(s) over the element,"
            f" {r['already_scored_over_element']} already deleted, {r['requests']} request(s)"
        )
    say(
        f"{sum(1 for r in rows if r['graded'])} graded of {len(rows)} candidates;"
        f" {sum(1 for r in rows if not r['askable'])} died at the reach filter;"
        f" {sum(r['requests'] for r in rows)} request(s) to spend"
    )


def score_intervals(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Delete the published elements no annotation drew, one request each."""
    t0 = time.time()
    asked = [r for r in rows if r["requests"]]
    adapter = AlphaGenomeAdapter()
    scorer = adapter._live_scorer(threshold=0.0)  # noqa: SLF001  (as scripts/loci_score.py does)
    elements: list[dict[str, Any]] = []
    read: list[dict[str, Any]] = []
    by_name = cands.by_name()
    with the_key(f"known-locus candidates, {len(asked)} stated intervals"):
        for row in asked:
            ctx = Context(row["chrom"])
            try:
                say(f"{row['locus']}: deleting {row['chrom']}:{row['element'][0]}-{row['element'][1]}")
                hit = ctx.score_region(scorer, row["element"][0], row["element"][1])
            finally:
                ctx.close()
            # score_region substitutes an overlapping ENCODE element when there is one, which would
            # answer about a different interval than the published one. Recorded, never hidden.
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
                    "tissue": pred.get("tissue"),
                    "top_genes": [g["gene"] for g in (hit.get("top_genes") or []) if g.get("gene")][:6],
                    "hit": bool(pred.get("gene") in targets),
                    "targets": targets,
                }
            )
            say(f"{row['locus']}: named {pred.get('gene')} first among coding genes")
    return {
        "result": cands.INTERVALS,
        "elements": elements,
        "read": read,
        "note": (
            "AlphaGenome deletions of intervals genomeos/benchmark/loci_candidates.py STATED from a"
            " publication, which no registry drew an element over. Read by loci._deletion_rows and"
            " labelled `stated_interval`, so these rows are never pooled with registry ones"
        ),
        "requests_this_run": len(read),
        "seconds": round(time.time() - t0, 1),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--plan", action="store_true", help="reach and cost per candidate; no request")
    ap.add_argument("--intervals", action="store_true", help="delete the elements no annotation covers")
    ap.add_argument("--score", action="store_true", help="read and score the graded candidates")
    ap.add_argument("--run", action="store_true", help="with --intervals: actually spend the requests")
    ap.add_argument("--quota-handed", action="store_true")
    ap.add_argument("--no-network", action="store_true", help="local layers only, for a quick check")
    args = ap.parse_args(argv)
    if not (args.plan or args.intervals or args.score):
        args.plan = True

    rows = cands.plan()
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
        out = score_intervals(rows)
        say(f"saved {save_result(cands.INTERVALS, out)}")

    if args.score:
        out = cands.run(network=not args.no_network, progress=say)
        a = out["aggregate"]
        say(
            f"target derived {a['target_derived']['k']}/{a['target_derived']['n']},"
            f" heuristic {a['target_heuristic']['k']}/{a['target_heuristic']['n']},"
            f" chance floor {a['target_by_chance']['expected']}"
        )
        print(json.dumps(out["reach"], indent=1))
        print(json.dumps(out["beside_the_panel"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
