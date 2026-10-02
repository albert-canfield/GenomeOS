# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registered counts of direction rule v2 against v1, from the deletion cache and nothing else.

0 model requests, no network, no money: every value read here is already on disk, in
`data/results/enhancer_targets_*` and `data/knowledge/alphagenome/elements/<chrom>.json.gz`.

What is reported is fixed in `attribution/direction_v2.COUNTS_REGISTERED` and was written and
committed before any count was taken. Nothing here validates anything in either direction
(`direction_v2.COUNTS_VALIDATE_NOTHING`), and 0 flips is arithmetic rather than agreement
(`direction_v2.FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE`).

Usage: python3 scripts/direction_v2.py [--chrom chr21 ...] [--all]
"""

from __future__ import annotations

import argparse
from collections import Counter
from typing import Any

from genomeos.attribution import compile as cp
from genomeos.attribution import direction_v2 as dv

#: the whole assembly, in the order the project names it
ALL_CHROMOSOMES = tuple(f"chr{i}" for i in range(1, 23)) + ("chrX", "chrY")


def chromosome_counts(chrom: str) -> dict[str, Any]:
    """Every predicted-layer rule of one chromosome under both versions of the rule.

    One rule per attributed element, which is the population `compile.compile_chromosome` writes a
    `rule ... activates|inhibits` line for. The measured layer is excluded: a measured rule carries a
    screen's effect and no AlphaGenome direction of its own.
    """
    elements = cp._attributed(chrom)
    tally: Counter[str] = Counter()
    reasons: Counter[str] = Counter()
    for e in elements:
        pc = e["predicted_coding"]
        v1 = dv.direction_call(pc, dv.V1)
        axis = "represses_target" if v1.action == dv.INHIBITS else "activates_target"
        tally[f"v1_{axis}"] += 1
        row = dv.cached_row(chrom, e["id"], pc.get("gene"))
        v2 = dv.direction_call(pc, dv.V2, row=row) if row is not None else dv.direction_call(pc, dv.V2)
        if not v2.resolved:
            tally[f"unresolved_from_{axis}"] += 1
            reasons[f"{axis}/{v2.reason}"] += 1
        elif v2.action != v1.action:
            tally[f"flip_from_{axis}"] += 1
        else:
            tally[f"agrees_from_{axis}"] += 1
        # SELECTION_EXCLUDED_DIAGNOSTIC: a property of the cache, not a v2 output and not a count of
        # wrong rules. by_cell[cell] against the sign of the value the maximum selected.
        if row is not None:
            other = [v for src, v in dv.retained_values(row, v1.cell) if src == "by_cell"]
            selected = float(pc["log2_fold_change"])
            if other and other[0] != 0.0 and (other[0] > 0) != (selected > 0):
                tally["by_cell_value_opposite_to_the_selected_extreme"] += 1
    return {"chrom": chrom, "rules": len(elements), "tally": dict(tally), "reasons": dict(reasons)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="the registered counts of direction rule v2 against v1")
    ap.add_argument("--chrom", action="append", default=None, help="repeatable; default chr21")
    ap.add_argument("--all", action="store_true", help="every chromosome of the assembly")
    args = ap.parse_args(argv)
    chroms = list(args.chrom or (ALL_CHROMOSOMES if args.all else ["chr21"]))
    total: Counter[str] = Counter()
    reasons: Counter[str] = Counter()
    rules = 0
    for chrom in chroms:
        r = chromosome_counts(chrom)
        rules += r["rules"]
        total.update(r["tally"])
        reasons.update(r["reasons"])
        print(f"{chrom}: {r['rules']} predicted rules")
        for k in sorted(r["tally"]):
            print(f"    {k}: {r['tally'][k]}")
    print()
    print(f"population: {rules} predicted-layer rules over {len(chroms)} chromosome(s)")
    for k in sorted(total):
        print(f"  {k}: {total[k]}")
    print("why v2 did not resolve, by the v1 axis it came from:")
    for k in sorted(reasons):
        print(f"  {k}: {reasons[k]}")
    print()
    print(dv.COUNTS_VALIDATE_NOTHING)
    print(dv.NOT_AN_IMPROVEMENT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
