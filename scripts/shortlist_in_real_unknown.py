# SPDX-License-Identifier: AGPL-3.0-or-later
"""Do the calibration's reliable targets fall inside the real unknown? The library's internal bridge.

    uv run python scripts/shortlist_in_real_unknown.py

The CRISPRi screens called 105 of 105 pairs regulated where the predicted drop clears 0.2
(`target_calibration`), and 24,114 of the sweep's coding targets reach that drop genome-wide. The
882 blocks of the real unknown are the sequence the MPRA library of `data/knowledge/library` tiles.
Where the two overlap, an oligo carries a prior from a measurement that is not the model's own
effect size, and the library gets a direct test of whether the drop band holds off the population
the screens chose. Where they do not overlap, that is the result: the library then carries no
internal bridge to the CRISPRi evidence and the shortlist has to be tested somewhere else.

Reads committed results and the local element tables only; no network, no model request. Saves
data/results/shortlist_in_real_unknown.json.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from typing import Any

from genomeos.attribution import organise, target_calibration
from genomeos.results import save_result

SHORTLIST_DROP = 0.2  # the band where the screens called 105 of 105 (target_calibration)
CASES = ("syntax", "relaxed", "recent", "tolerant", "unmeasured")
OLIGO = 300  # the library's tile length, for the oligo count an overlap implies


def real_unknown_blocks(chrom: str) -> list[dict[str, Any]]:
    """The blocks of one chromosome's real unknown: constrained-unknown tier, copies removed."""
    return [r for r in organise.blocks(chrom) if r["tier"] == "constrained_unknown" and not r["copy"]]


def shortlist_elements(chrom: str) -> tuple[list[dict[str, Any]], Counter[str]]:
    """Deleted elements of one chromosome whose coding target clears the band, with what was refused."""
    path = target_calibration.ELEMENTS / f"{chrom}.json"
    if not path.exists():
        return [], Counter({"no_elements_cached": 1})
    els = json.loads(path.read_text())
    starts = target_calibration.gene_starts(chrom, target_calibration.REFERENCE)
    classes = (
        target_calibration.registry_classes(chrom, target_calibration.RESULTS)
        if (target_calibration.RESULTS / f"ccres_{chrom}.bed.gz").exists()
        else {}
    )
    out, refused = [], Counter()
    for el in els:
        row = target_calibration.element_row(el, "predicted_coding", starts, classes)
        if row is None:
            refused["no_coding_target"] += 1
            continue
        if "missing" in row:
            refused[row["missing"]] += 1
            continue
        below = row["drop"] <= SHORTLIST_DROP
        if below:
            refused["below_the_band"] += 1
        out.append(
            {
                "id": el.get("id", ""),
                "start": el["start"],
                "end": el["end"],
                "gene": el["predicted_coding"]["gene"],
                "drop": round(row["drop"], 4),
                "tissue": el["predicted_coding"].get("tissue", ""),
                "in_the_band": not below,
            }
        )
    return out, refused


def intersect(chrom: str) -> dict[str, Any]:
    """The shortlist elements of one chromosome that lie inside one of its real-unknown blocks."""
    blocks = real_unknown_blocks(chrom)
    elements, refused = shortlist_elements(chrom)
    hits: dict[str, list[dict[str, Any]]] = defaultdict(list)
    scored: dict[str, dict[str, Any]] = {}
    cases: dict[str, str] = {}
    for b in blocks:
        bid = f"{chrom}:{b['start']}-{b['end']}"
        cases[bid] = b["case"] or "unmeasured"
        for e in elements:
            if e["start"] < b["end"] and e["end"] > b["start"]:
                hits[bid].append(e)
        inside = hits.get(bid, [])
        if inside:
            scored[bid] = {
                "elements": len(inside),
                "max_drop": round(max(e["drop"] for e in inside), 4),
                "in_the_band": sum(e["in_the_band"] for e in inside),
            }
    return {
        "blocks": len(blocks),
        "blocks_bp": sum(b["length"] for b in blocks),
        "scored_elements": len(elements),
        "shortlist_elements": sum(e["in_the_band"] for e in elements),
        "refused": dict(sorted(refused.items())),
        "cases": cases,
        "scored_blocks": scored,
        "hits": {
            k: [e for e in v if e["in_the_band"]] for k, v in hits.items() if any(e["in_the_band"] for e in v)
        },
    }


def main() -> int:
    per_chrom, by_case, hit_blocks, elements_in = {}, Counter(), {}, 0
    scored_by_case: Counter[str] = Counter()
    scored_blocks: dict[str, dict[str, Any]] = {}
    blocks_total = shortlist_total = 0
    refused_total: Counter[str] = Counter()
    for chrom in target_calibration.CHROMS:
        r = intersect(chrom)
        blocks_total += r["blocks"]
        shortlist_total += r["shortlist_elements"]
        refused_total.update(r["refused"])
        for bid, info in r["scored_blocks"].items():
            scored_by_case[r["cases"][bid]] += 1
            scored_blocks[bid] = {"case": r["cases"][bid], **info}
        for bid, els in r["hits"].items():
            case = r["cases"][bid]
            by_case[case] += 1
            elements_in += len(els)
            hit_blocks[bid] = {
                "case": case,
                "elements": [e["id"] for e in els],
                "genes": sorted({e["gene"] for e in els}),
                "max_drop": max(e["drop"] for e in els),
            }
        per_chrom[chrom] = {
            "real_unknown_blocks": r["blocks"],
            "shortlist_elements": r["shortlist_elements"],
            "blocks_with_a_shortlist_element": len(r["hits"]),
        }
        print(
            f"{chrom}: {r['blocks']} real-unknown blocks, {r['shortlist_elements']} shortlist "
            f"elements, {len(r['hits'])} blocks hit",
            flush=True,
        )
    blocks_by_case: Counter[str] = Counter()
    for chrom in target_calibration.CHROMS:
        for b in real_unknown_blocks(chrom):
            blocks_by_case[b["case"] or "unmeasured"] += 1
    payload = {
        "question": (
            "do the calibration's reliable targets (coding drop above "
            f"{SHORTLIST_DROP}, the band where the CRISPRi screens called 105 of 105) fall inside "
            "the 882 blocks of the real unknown, so the MPRA library carries its own bridge to the "
            "screens' evidence"
        ),
        "shortlist_drop": SHORTLIST_DROP,
        "real_unknown_blocks": blocks_total,
        "real_unknown_blocks_by_case": dict(blocks_by_case),
        "shortlist_elements_genome_wide": shortlist_total,
        "blocks_with_any_scored_element": len(scored_blocks),
        "blocks_with_any_scored_element_by_case": dict(scored_by_case),
        "max_drop_of_the_best_scored_element_per_block": {
            k: v["max_drop"] for k, v in sorted(scored_blocks.items(), key=lambda kv: -kv[1]["max_drop"])[:20]
        },
        "blocks_with_a_shortlist_element": sum(by_case.values()),
        "blocks_with_a_shortlist_element_by_case": dict(by_case),
        "shortlist_elements_inside_a_block": elements_in,
        "oligos_an_overlap_implies": elements_in * 2 if elements_in else 0,
        "oligo_length": OLIGO,
        "blocks": hit_blocks,
        "per_chromosome": per_chrom,
        "refused": dict(sorted(refused_total.items())),
        "evidence": (
            "predicted: AlphaGenome deletion per element, banded by the CRISPRi-fitted calibration; "
            "inferred: the real-unknown tier from constraint on two axes with copies removed"
        ),
        "reading": (
            "empty means the shortlist and the real unknown do not overlap, so the library holds no "
            "internal positive set tied to the screens and the drop band has to be tested elsewhere"
        ),
    }
    path = save_result("shortlist_in_real_unknown", payload)
    print(
        f"\n{blocks_total} real-unknown blocks, {shortlist_total} shortlist elements genome-wide, "
        f"{sum(by_case.values())} blocks hit, {elements_in} elements inside a block"
    )
    print(f"by case: {dict(by_case)} of {dict(blocks_by_case)}")
    print(f"blocks holding any scored element: {len(scored_blocks)} ({dict(scored_by_case)})")
    print(f"-> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
