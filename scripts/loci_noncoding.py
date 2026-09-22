# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fourth frame's non-coding drops, scored as a frame of their own. No request, ever.

    uv run python scripts/loci_noncoding.py --frame     # the branch, the Ensembl join, n; free
    uv run python scripts/loci_noncoding.py --plan      # reach and cost; free, and must read 0
    uv run python scripts/loci_noncoding.py --score     # read and score through loci.build

`loci_fourth.assess` stops fourteen held-out elements at its first branch, for having a regulated
target that is not protein coding in GENCODE. `--frame` joins the benchmark file's own
`measuredGeneEnsemblId` to GENCODE and shows that twelve of them are stale gene symbols, so the frame
this scores is n = 2. `--score` never spends a request: both elements are in the finished sweep, and
`loci_noncoding.plan` raises rather than let the run buy one.
"""

from __future__ import annotations

import argparse
import json
import sys

from genomeos.benchmark import loci_noncoding as nc


def say(msg: str) -> None:
    print(msg, file=sys.stderr)


def show_frame() -> dict:
    branch = nc.dropped_for_non_coding()
    resolved = nc.resolve(branch)
    split = nc.classify(resolved)
    print(f"the branch stops at {len(branch)} element(s); the Ensembl join splits them:")
    for r in resolved:
        for x in r["resolved"]:
            kind = "STALE SYMBOL" if x["stale_symbol"] else "non-coding"
            print(
                f"  {r['chrom']}:{r['start']}-{r['end']} {r['cell']:8s}"
                f" {x['published_symbol']:10s} {x['ensembl_id']} -> {x['current_symbol']}"
                f" ({x['gene_type']})  {kind}"
            )
    print(f"\nn = {len(split['non_coding'])} genuinely non-coding; {len(split['stale_symbol'])} stale")
    rows = nc.with_nearest_coding(split["stale_symbol"])
    corrected = nc.corrected_rule(rows)
    print("\nwhat the fourth frame's rule does to the stale ones once the symbol is current:")
    for k, v in corrected["by_verdict"].items():
        print(f"  {v:3d}  {k}")
    return {"branch": len(branch), "non_coding": len(split["non_coding"]), "corrected": corrected}


def show_plan() -> list[dict]:
    branch = nc.classify(nc.resolve(nc.dropped_for_non_coding()))["non_coding"]
    panel = nc.frame(nc.with_nearest_coding(branch))
    rows = nc.plan(panel)
    for r in rows:
        print(
            f"{r['locus']:32s} askable={r['askable']!s:5s}"
            f" already_scored={r['already_scored_over_element']} requests={r['requests']}"
        )
    print(f"\n{sum(r['requests'] for r in rows)} request(s) to spend (budget {nc.REQUEST_BUDGET})")
    return rows


def show_scores(out: dict) -> None:
    agg = out["aggregate"]
    print(json.dumps({k: agg[k] for k in sorted(agg) if "target" in k or k == "n"}, indent=1))
    print("\nper layer:")
    for name, v in out["layers"]["by_layer"].items():
        print(f"  {name:11s} {v['provenance']:10s} {v['hit']}/{v['n']}  pending={v['pending']}")
    print("\nwhat it named instead:")
    for layer, tally in out["layers"]["what_it_named_instead"].items():
        print(f"  {layer:11s} {tally}")
    b = out["baseline"]
    print(f"\nregistered baseline (any-gene, no coding-first break): {b['k']}/{b['n']}")
    for p in b["per_locus"]:
        print(f"  {p['locus']:32s} would name {p['would_name']} at {p['effect']} -> hit={p['hit']}")
    print(f"\nrequests spent: {out['requests_spent']}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--frame", action="store_true", help="the branch and the Ensembl join; free")
    ap.add_argument("--plan", action="store_true", help="reach and cost per locus; free")
    ap.add_argument("--score", action="store_true", help="read and score; no request")
    ap.add_argument("--no-network", action="store_true", help="local layers only")
    args = ap.parse_args(argv)
    if not (args.frame or args.plan or args.score):
        ap.error("pass --frame, --plan or --score")
    if args.frame:
        show_frame()
    if args.plan:
        show_plan()
    if args.score:
        show_scores(nc.run(network=not args.no_network, progress=say))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
