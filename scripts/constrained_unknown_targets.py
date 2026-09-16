# SPDX-License-Identifier: AGPL-3.0-or-later
"""What the completed sweep says about the real unknown: the constrained-unknown blocks that are not copies.

    uv run python scripts/constrained_unknown_targets.py [--chroms chr21,chr22]

Milestone 1.3 asks for the constrained-unknown blocks attributed to a gene and a tissue. The sweep
finished on 2026-09-16, so every ENCODE element inside a node has been deleted in AlphaGenome on every
chromosome, and the answer is now arithmetic over tables already on disk: no model request.

For every UNKNOWN block (the organiser's join of tier, human axis and copy flag) this counts the
elements that lie inside it, how many of them move a gene, and what they name. It reports the
constrained-unknown tier with the copies taken out -- the organiser's "real unknown" -- split by the
case the two axes make, and it reports it **against the other tiers**, because a rate with no control
is what this project has had to withdraw three times:

- per element, the share that moves a gene, which does not depend on how long a block is;
- per block, the share carrying at least one element that moves a gene, inside length deciles, since a
  longer block holds more elements for no biological reason.

The result is `constrained_unknown_targets`. It names genes but claims nothing about them: the same
model names a target at 87% of matched random windows, so a named target here is a lead, not a finding.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any

from genomeos.attribution import organise
from genomeos.results import save_result

ELEMENTS = Path("data/knowledge/alphagenome/all_elements")
MIN_LOG2 = 0.1  # the sweep's own threshold for "this deletion moved a gene"
DECILES = 10


def elements_of(chrom: str) -> list[dict[str, Any]]:
    """The chromosome's scored elements, or nothing if the sweep never wrote it."""
    p = ELEMENTS / f"{chrom}.json"
    if not p.exists():
        return []
    return json.loads(p.read_text())


def moves(element: dict[str, Any]) -> bool:
    """Did deleting this element move a gene, by the sweep's own threshold?"""
    pred = element.get("predicted") or {}
    log2 = pred.get("log2_fold_change")
    return bool(pred.get("gene")) and log2 is not None and abs(log2) >= MIN_LOG2


def inside(block: dict[str, Any], els: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Elements whose midpoint falls in the block: an element is not split between two blocks."""
    lo, hi = block["start"], block["end"]
    return [e for e in els if lo <= (e["start"] + e["end"]) // 2 < hi]


def read_chromosome(chrom: str) -> list[dict[str, Any]]:
    """Every UNKNOWN block of one chromosome with the sweep's answers attached."""
    els = elements_of(chrom)
    if not els:
        return []
    els.sort(key=lambda e: e["start"])
    rows = []
    for block in organise.blocks(chrom):
        got = inside(block, els)
        moving = [e for e in got if moves(e)]
        rows.append(
            {
                "chrom": chrom,
                "block": f"{chrom}:{block['start']}-{block['end']}",
                "length": block["length"],
                "tier": block.get("tier"),
                "case": block.get("case"),
                "copy": bool(block.get("copy")),
                "elements": len(got),
                "moving": len(moving),
                "targets": sorted({(e["predicted"] or {}).get("gene") for e in moving} - {None}),
                "cells": sorted({(e["predicted"] or {}).get("tissue") for e in moving} - {None}),
                "strongest": max(
                    (abs((e["predicted"] or {}).get("log2_fold_change", 0)) for e in moving), default=0
                ),
            }
        )
    return rows


def rate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Per-element and per-block rates for a set of blocks, with the counts they came from."""
    elements = sum(r["elements"] for r in rows)
    moving = sum(r["moving"] for r in rows)
    with_any = [r for r in rows if r["elements"]]
    return {
        "blocks": len(rows),
        "blocks_with_an_element": len(with_any),
        "elements": elements,
        "moving": moving,
        "moves_per_element": round(moving / elements, 4) if elements else None,
        "blocks_with_a_moving_element": sum(1 for r in rows if r["moving"]),
        "share_of_blocks_with_an_element_that_move": (
            round(sum(1 for r in with_any if r["moving"]) / len(with_any), 4) if with_any else None
        ),
        "median_length_kb": round(median([r["length"] for r in rows]) / 1000, 1) if rows else None,
    }


def by_length_decile(target: list[dict[str, Any]], control: list[dict[str, Any]]) -> dict[str, Any]:
    """The per-block comparison inside length deciles, so block size cannot produce the difference."""
    pool = sorted(r["length"] for r in control + target)
    if not pool:
        return {}
    cuts = [pool[int(len(pool) * i / DECILES)] for i in range(1, DECILES)]

    def decile(length: int) -> int:
        return sum(length >= c for c in cuts)

    buckets: dict[int, dict[str, list]] = defaultdict(lambda: {"target": [], "control": []})
    for r in target:
        buckets[decile(r["length"])]["target"].append(r)
    for r in control:
        buckets[decile(r["length"])]["control"].append(r)
    rows, weighted, weight = {}, 0.0, 0
    for d in sorted(buckets):
        t, c = buckets[d]["target"], buckets[d]["control"]
        if not t or not c:
            continue
        t_rate = sum(1 for r in t if r["moving"]) / len(t)
        c_rate = sum(1 for r in c if r["moving"]) / len(c)
        rows[str(d)] = {
            "target_blocks": len(t),
            "control_blocks": len(c),
            "target": round(t_rate, 4),
            "control": round(c_rate, 4),
            "difference": round(t_rate - c_rate, 4),
        }
        weighted += (t_rate - c_rate) * len(t)
        weight += len(t)
    return {
        "per_decile": rows,
        "deciles_compared": len(rows),
        "difference_weighted_by_target_blocks": round(weighted / weight, 4) if weight else None,
    }


def collect(chroms: list[str]) -> dict[str, Any]:
    t0 = time.time()
    rows: list[dict[str, Any]] = []
    for chrom in chroms:
        got = read_chromosome(chrom)
        rows.extend(got)
        print(f"{chrom}: {len(got)} blocks, {sum(r['elements'] for r in got)} elements inside", flush=True)

    real_unknown = [r for r in rows if r["tier"] == "constrained_unknown" and not r["copy"]]
    copies = [r for r in rows if r["tier"] == "constrained_unknown" and r["copy"]]
    others = {t: [r for r in rows if r["tier"] == t] for t in ("neutral", "fossil", "regulatory")}
    by_case = {
        case: rate([r for r in real_unknown if r["case"] == case])
        for case in sorted({r["case"] for r in real_unknown if r["case"]})
    }
    named = [r for r in real_unknown if r["targets"]]
    return {
        "result": "constrained_unknown_targets",
        "chromosomes": chroms,
        "min_log2": MIN_LOG2,
        "real_unknown": rate(real_unknown),
        "constrained_unknown_copies": rate(copies),
        "other_tiers": {t: rate(v) for t, v in others.items()},
        "by_case": by_case,
        "against_neutral_in_length_deciles": by_length_decile(real_unknown, others["neutral"]),
        "blocks_with_a_named_target": len(named),
        "named": sorted(
            (
                {
                    "block": r["block"],
                    "case": r["case"],
                    "elements": r["elements"],
                    "moving": r["moving"],
                    "targets": r["targets"][:6],
                    "cells": r["cells"][:4],
                    "strongest_abs_log2": round(r["strongest"], 3),
                }
                for r in named
            ),
            key=lambda r: -r["strongest_abs_log2"],
        ),
        "reading": (
            "a named target here is a lead and not a finding: the same model names a target at 87% of "
            "matched random windows, which is why the rates are reported against the other tiers and "
            "inside length deciles rather than on their own"
        ),
        "seconds": round(time.time() - t0, 1),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="", help="comma-separated; default every chromosome with a table")
    ap.add_argument(
        "--no-save",
        action="store_true",
        help="print without writing the result: a run on a few chromosomes must not overwrite the genome's",
    )
    args = ap.parse_args(argv)
    chroms = (
        args.chroms.split(",")
        if args.chroms
        else sorted((p.stem for p in ELEMENTS.glob("chr*.json")), key=lambda c: (len(c), c))
    )
    out = collect(chroms)
    if args.no_save:
        print("\nnot saved (--no-save)")
    else:
        print(f"\nsaved {save_result(out['result'], out)}")
    ru, ne = out["real_unknown"], out["other_tiers"]["neutral"]
    print(f"real unknown: {ru['blocks']} blocks, {ru['elements']} elements, moves {ru['moves_per_element']}")
    print(f"neutral tier: {ne['blocks']} blocks, {ne['elements']} elements, moves {ne['moves_per_element']}")
    deciles = out["against_neutral_in_length_deciles"]
    print(f"per block inside length deciles: {deciles.get('difference_weighted_by_target_blocks')}")
    print(f"blocks with a named target: {out['blocks_with_a_named_target']} of {ru['blocks']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
