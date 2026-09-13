# SPDX-License-Identifier: AGPL-3.0-or-later
"""The block organiser: every UNKNOWN block read with all the evidence there is, copies first.

The budget gave each block a tier from its sequence class and mammalian constraint.
Two more readings have since been measured per block: the human axis (area J,
gnomAD Gnocchi: is the block depleted of variation among people, and which of the
four cases does the pair of axes make), and the copy flag (curated segmental
duplications: is the block a copy of sequence elsewhere, where a mammalian
alignment is paralogous and human variants do not map). This module joins the
three per block and re-reads the tier in that order:

1. a copy (half or more of the block duplicated) is read as a copy before anything
   else is said about it, whatever its tier;
2. a constrained_unknown block that is not a copy is read by its case: constrained
   on both axes is the sharpest candidate for something unannotated; held across
   mammals but variable among people is a frame whose value varies, or lost
   function; measured on neither is unmeasured, not neutral;
3. the other tiers keep their label, with the copy flag and the case attached.

Nothing is recomputed: the three inputs are `budget_<chrom>`, `variation_<chrom>` and
`duplication_<chrom>`. The result keeps the per-tier tallies, the copies and the
candidates, not a second copy of every block; `blocks()` returns the full join on
demand. The real unknown of a chromosome is what remains of the constrained_unknown
tier after the copies, split by case.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Any

from genomeos.attribution.budget import TIERS
from genomeos.results import RESULTS_DIR, load_result, save_result

COPY_MIN = 0.5  # duplicated fraction from which a block is read as a copy first
CASE_READING = {
    "syntax": "constrained across mammals and among people, unannotated: the sharpest candidate",
    "relaxed": "held across mammals, variable among people: a frame whose value varies, or lost function",
    "recent": "free across mammals, constrained among people: recent function or selection",
    "tolerant": "free on both scales at kilobase resolution: check the alignment before attributing",
    "unmeasured": "no human-axis coverage (the block is under a kilobase or on chrY): mammals only",
}


def blocks(chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """Every UNKNOWN block of the chromosome with its tier, both constraint axes, the copy flag
    and the attributed elements inside it."""
    budget = load_result(f"budget_{chrom}", results_dir)
    if not budget:
        raise FileNotFoundError(f"no budget_{chrom} result; run genomeos budget --chrom {chrom}")
    var = {
        (b["start"], b["end"]): b
        for b in (load_result(f"variation_{chrom}", results_dir) or {}).get("blocks", [])
    }
    dup = {
        (b["start"], b["end"]): b
        for b in (load_result(f"duplication_{chrom}", results_dir) or {}).get("blocks", [])
    }
    elements = _attributed(chrom, results_dir)
    out = []
    for b in sorted(budget["blocks"], key=lambda x: x["start"]):
        key = (b["start"], b["end"])
        v, d = var.get(key), dup.get(key)
        case = (v or {}).get("case") or {}
        gn = (v or {}).get("gnocchi") or {}
        dup_frac = (d or {}).get("duplicated_fraction")
        inside = [e for e in elements if b["start"] <= e["start"] and e["end"] <= b["end"]]
        row = {
            "start": b["start"],
            "end": b["end"],
            "length": b["length"],
            "class": b["class"],
            "tier": b["guess"]["tier"],
            "confidence": b["guess"]["confidence"],
            "mammal_fraction": (b.get("phylop") or {}).get("fraction_above"),
            "conserved_elements": (b.get("elements") or {}).get("n"),
            "human_fraction": gn.get("fraction_above"),
            "case": case.get("case") or ("unmeasured" if v is not None else None),
            "duplicated_fraction": dup_frac,
            "partners": (d or {}).get("pairs"),
            "copy": bool(dup_frac is not None and dup_frac >= COPY_MIN),
            "attributed_elements": len(inside),
            "targets": sorted({e["target"] for e in inside}),
        }
        row["reading"] = reading(row)
        out.append(row)
    return out


def _attributed(chrom: str, results_dir: Path) -> list[dict[str, Any]]:
    from genomeos.attribution.targets import attributed

    return [
        {"id": e["id"], "start": e["start"], "end": e["end"], "target": e["predicted_coding"]["gene"]}
        for e in attributed(chrom, results_dir)
    ]


def reading(row: dict[str, Any]) -> str:
    """One line per block, copies first."""
    if row["copy"]:
        n = row.get("partners") or 0
        return (
            f"copy: {row['duplicated_fraction'] * 100:.0f}% of the block is a curated segmental duplication "
            f"({n} partner{'s' if n != 1 else ''}); read the source first, attribute second"
        )
    if row["tier"] == "constrained_unknown":
        return CASE_READING.get(row["case"] or "unmeasured", CASE_READING["unmeasured"])
    base = row["tier"]
    if row["case"] and row["case"] != "unmeasured":
        base += f", {row['case']} on the two axes"
    if row["attributed_elements"]:
        n = row["attributed_elements"]
        base += f", {n} attributed element{'s' if n != 1 else ''}"
    return base


def organise(chrom: str, results_dir: Path = RESULTS_DIR, top: int = 25) -> dict[str, Any]:
    t0 = time.time()
    rows = blocks(chrom, results_dir)
    by_tier: dict[str, dict[str, Any]] = {}
    for t in TIERS:
        sel = [r for r in rows if r["tier"] == t]
        copies = [r for r in sel if r["copy"]]
        cases: dict[str, dict[str, int]] = {}
        for r in sel:
            if r["copy"]:
                continue
            c = r["case"] or "unmeasured"
            cases.setdefault(c, {"blocks": 0, "bp": 0})
            cases[c]["blocks"] += 1
            cases[c]["bp"] += r["length"]
        by_tier[t] = {
            "blocks": len(sel),
            "bp": sum(r["length"] for r in sel),
            "copies": len(copies),
            "copies_bp": sum(r["length"] for r in copies),
            "after_copies": len(sel) - len(copies),
            "after_copies_bp": sum(r["length"] for r in sel if not r["copy"]),
            "cases_after_copies": cases,
        }
    cu = [r for r in rows if r["tier"] == "constrained_unknown" and not r["copy"]]
    rank = {"syntax": 0, "relaxed": 1, "recent": 2, "unmeasured": 3, "tolerant": 4, None: 3}
    cu.sort(key=lambda r: (rank.get(r["case"], 3), -(r["mammal_fraction"] or 0)))
    copies = sorted((r for r in rows if r["copy"]), key=lambda r: -r["length"])
    return {
        "chrom": chrom,
        "blocks": len(rows),
        "with_human_axis": sum(1 for r in rows if r["case"] and r["case"] != "unmeasured"),
        "with_copy_flag": sum(1 for r in rows if r["duplicated_fraction"] is not None),
        "copies": len(copies),
        "copies_bp": sum(r["length"] for r in copies),
        "copy_min_fraction": COPY_MIN,
        "by_tier": by_tier,
        "real_unknown": {
            "blocks": len(cu),
            "bp": sum(r["length"] for r in cu),
            "by_case": by_tier["constrained_unknown"]["cases_after_copies"],
        },
        "candidates": cu[:top],
        "largest_copies": copies[:top],
        "inputs": {
            "budget": f"budget_{chrom}",
            "variation": f"variation_{chrom}",
            "duplication": f"duplication_{chrom}",
            "elements": [f"constrained_targets_{chrom}", f"enhancer_targets_{chrom}"],
        },
        "evidence": {
            "tier": "inferred: sequence class plus Zoonomia constraint (the budget)",
            "case": "inferred: Zoonomia phyloP against gnomAD Gnocchi (area J); unmeasured below a kilobase",
            "copy": "curated: UCSC genomicSuperDups, >= 1 kb at >= 90% identity; copy at >= 50% of the block",
        },
        "cost": {"seconds": round(time.time() - t0, 2)},
    }


def distil(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    per: dict[str, dict[str, Any]] = {}
    by_tier: dict[str, dict[str, Any]] = {
        t: {"blocks": 0, "bp": 0, "copies": 0, "copies_bp": 0, "cases_after_copies": {}} for t in TIERS
    }
    real = {"blocks": 0, "bp": 0, "by_case": {}}
    for c in chroms:
        r = load_result(f"organised_{c}", results_dir)
        if not r:
            continue
        per[c] = {
            "blocks": r["blocks"],
            "copies": r["copies"],
            "real_unknown_blocks": r["real_unknown"]["blocks"],
            "real_unknown_bp": r["real_unknown"]["bp"],
            "syntax_candidates": (r["real_unknown"]["by_case"].get("syntax") or {}).get("blocks", 0),
        }
        for t, v in r["by_tier"].items():
            d = by_tier[t]
            for k in ("blocks", "bp", "copies", "copies_bp"):
                d[k] += v[k]
            for case, cv in v["cases_after_copies"].items():
                e = d["cases_after_copies"].setdefault(case, {"blocks": 0, "bp": 0})
                e["blocks"] += cv["blocks"]
                e["bp"] += cv["bp"]
        real["blocks"] += r["real_unknown"]["blocks"]
        real["bp"] += r["real_unknown"]["bp"]
        for case, cv in r["real_unknown"]["by_case"].items():
            e = real["by_case"].setdefault(case, {"blocks": 0, "bp": 0})
            e["blocks"] += cv["blocks"]
            e["bp"] += cv["bp"]
    return {
        "chromosomes": len(per),
        "per_chromosome": per,
        "by_tier": by_tier,
        "real_unknown": real,
        "copy_min_fraction": COPY_MIN,
        "note": (
            "The budget's constrained_unknown tier re-read with copies set apart and the human axis "
            "attached: the real unknown is what remains, split by the case the two axes make."
        ),
    }


def run_and_save(chrom: str, results_dir: Path = RESULTS_DIR, **kw) -> dict[str, Any]:
    out = organise(chrom, results_dir, **kw)
    save_result(f"organised_{chrom}", out, results_dir)
    return out
