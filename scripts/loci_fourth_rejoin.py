# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fourth frame re-drawn with the benchmark file's own Ensembl ids joined. No request, ever.

    uv run python scripts/loci_fourth_rejoin.py --draw    # both draws, and the difference; free
    uv run python scripts/loci_fourth_rejoin.py --cost    # what scoring the difference would spend
    uv run python scripts/loci_fourth_rejoin.py --score   # score the corrected draw

Step 4 of `loci_fourth`'s rule asks whether a regulated target is protein coding, and it asked it of
a symbol column as old as the screens that filled it. `loci_fourth.current_symbols` joins the file's
own `measuredGeneEnsemblId` to GENCODE and `select(join_ensembl_ids=True)` applies the rule to the
names GENCODE uses now - the rule as written rather than as implemented. `--draw` prints both draws
side by side, including the loci the correction REMOVES, because a stale symbol can have let an
element in as easily as it kept one out.

The corrected frame is scored as `loci_fourth_rejoined`, beside `loci_fourth` and never over it:
section 22's numbers stay reproducible and both readings stay in the record.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from genomeos.benchmark import loci_fourth as lf

REJOINED = "loci_fourth_rejoined"
#: this draw's own distilled GTEx hits. `loci.stream_gtex` clears any directory whose window
#: manifest changed, and the corrected draw's windows are not section 22's, so sharing
#: `data/knowledge/loci_fourth` would delete the cache the published frame is reproducible from.
GTEX_DIR = Path("data/knowledge/loci_fourth_rejoined")


def say(msg: str) -> None:
    print(msg, file=sys.stderr)


#: "every element the rule passes", which is what `loci_fourth.CAP_RAISE` set MAX_LOCI to mean when
#: it moved the cap from 24 to 60: 60 was the number of held-out elements that passed the rule as it
#: was then implemented, chosen so the cap could not act as a filter and nothing was left over for a
#: later hand to pick from. The corrected rule passes more than 60, so under that same principle the
#: cap is again the number that pass - not a new number anybody chose.
UNCAPPED = 10**9


def draws(cap: int = lf.MAX_LOCI) -> tuple[dict, dict]:
    as_implemented = lf.select(join_ensembl_ids=False, progress=say)
    as_written = lf.select(join_ensembl_ids=True, cap=cap, progress=say)
    return as_implemented, as_written


def show_draw(as_implemented: dict, as_written: dict) -> None:
    old = {e.locus: e for e in as_implemented["panel"]}
    new = {e.locus: e for e in as_written["panel"]}
    old_at = {(e.chrom, *e.element): e for e in as_implemented["panel"]}
    new_at = {(e.chrom, *e.element): e for e in as_written["panel"]}
    print(f"renamed symbols in the held-out arm: {len(as_written['ensembl_join']['renamed_symbols'])}")
    for published, now in as_written["ensembl_join"]["renamed_symbols"].items():
        print(f"  {published} -> {now}")
    print(
        f"\nelements whose regulated target was renamed:"
        f" {as_written['ensembl_join']['elements_with_a_renamed_target']}"
    )
    print(f"\ndrawn as implemented: {len(old)}   drawn as written: {len(new)}")
    added = [k for k in new_at if k not in old_at]
    removed = [k for k in old_at if k not in new_at]
    print(f"  ADDED by the join:   {len(added)}")
    for k in sorted(added):
        e = new_at[k]
        print(f"    {e.locus:38s} {e.chrom}:{e.element[0]}-{e.element[1]}  trap {e.nearest_gene_trap}")
    print(f"  REMOVED by the join: {len(removed)}")
    for k in sorted(removed):
        e = old_at[k]
        print(f"    {e.locus:38s} {e.chrom}:{e.element[0]}-{e.element[1]}  trap {e.nearest_gene_trap}")
    renamed_in_place = [
        (old_at[k].locus, new_at[k].locus)
        for k in old_at
        if k in new_at and old_at[k].locus != new_at[k].locus
    ]
    print(f"  same element, new name: {len(renamed_in_place)}")
    for was, now in sorted(renamed_in_place):
        print(f"    {was} -> {now}")
    print("\nrejections, as implemented:")
    for reason, n in sorted(as_implemented["rejected_by_reason"].items()):
        print(f"  {n:4d}  {reason}")
    print("rejections, as written:")
    for reason, n in sorted(as_written["rejected_by_reason"].items()):
        print(f"  {n:4d}  {reason}")


def show_cost(as_written: dict) -> int:
    lf.register_stated_intervals()
    rows = lf.plan(as_written, lf.RESULTS_DIR)
    askable = [r for r in rows if r["askable"]]
    need = [r for r in rows if r["requests"]]
    spend = sum(r["requests"] for r in rows)
    print(
        f"drawn {len(rows)}, askable {len(askable)}, already covered by the sweep"
        f" {sum(1 for r in askable if r['already_scored_over_element'])}, would spend {spend}"
    )
    for r in need:
        print(f"  {r['locus']}: {r['requests']} request(s), no scored deletion over the element")
    return spend


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--draw", action="store_true")
    ap.add_argument("--cost", action="store_true")
    ap.add_argument("--score", action="store_true")
    ap.add_argument("--budget", type=int, default=10, help="refuse to score if the draw costs more")
    ap.add_argument(
        "--uncapped",
        action="store_true",
        help="cap = every element the corrected rule passes, which is what CAP_RAISE set 60 to mean",
    )
    args = ap.parse_args()
    cap = UNCAPPED if args.uncapped else lf.MAX_LOCI

    if args.draw:
        show_draw(*draws(cap))
        return 0
    if args.cost:
        show_cost(lf.select(join_ensembl_ids=True, cap=cap, progress=say))
        return 0
    if args.score:
        as_written = lf.select(join_ensembl_ids=True, cap=cap, progress=say)
        spend = show_cost(as_written)
        if spend > args.budget:
            say(
                f"refusing: the corrected draw would spend {spend} requests, over the budget of {args.budget}"
            )
            return 1
        out = lf.run(draw=as_written, progress=say, name=REJOINED, gtex_dir=GTEX_DIR)
        agg = out["aggregate"]
        for key in ("target_derived", "target_derived_where_the_model_could_answer", "target_heuristic"):
            x = agg[key]
            print(f"  {key:52s} {x['k']}/{x['n']} ({x['rate']})")
        return 0
    ap.print_help()
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
