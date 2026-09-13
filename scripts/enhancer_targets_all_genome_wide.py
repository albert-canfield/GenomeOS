# SPDX-License-Identifier: AGPL-3.0-or-later
"""One table over every chromosome's all-elements summary, as the chain lands them.

    uv run python scripts/enhancer_targets_all_genome_wide.py

Reads data/results/enhancer_targets_all_chr*.json (summaries only; the element tables stay local) and
writes enhancer_targets_all_genome_wide.json: per chromosome the counts and rates, and the totals over the
chromosomes complete so far, with the requests spent.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

from genomeos.results import save_result


def aggregate(results_dir: Path = Path("data/results")) -> dict:
    rows = {}
    tot = {"elements": 0, "scored": 0, "named": 0, "coding_named": 0, "strong": 0, "silencer_like": 0}
    verdicts: dict[str, int] = {}
    requests = 0
    for f in sorted(results_dir.glob("enhancer_targets_all_chr*.json")):
        d = json.loads(f.read_text())
        s = d.get("summary") or {}
        vc = s.get("verdicts_coding") or {}
        coding_named = sum(n for k, n in vc.items() if k != "no predicted effect")
        rows[d["chrom"]] = {
            "elements": d.get("elements_total"),
            "scored": d.get("scored", 0),
            "complete": bool(d.get("complete")),
            "fraction_with_target": s.get("fraction_with_target"),
            "coding_target_agrees_with_nearest": s.get("coding_target_agrees_with_nearest"),
            "fraction_inside_domain": s.get("fraction_inside_domain"),
            "strong": s.get("strong"),
            "silencer_like": s.get("silencer_like"),
            "requests_this_run": d.get("requests_this_run"),
        }
        requests += d.get("requests_this_run") or 0
        if d.get("complete"):
            tot["elements"] += d.get("elements_total") or 0
            tot["scored"] += d.get("scored") or 0
            tot["named"] += s.get("with_predicted_target") or 0
            tot["coding_named"] += coding_named
            tot["strong"] += s.get("strong") or 0
            tot["silencer_like"] += s.get("silencer_like") or 0
            for k, n in vc.items():
                verdicts[k] = verdicts.get(k, 0) + n
    agree = verdicts.get("agrees with nearest TSS in domain", 0)
    inside = agree + verdicts.get("another gene in the same domain", 0)
    return {
        "chromosomes": rows,
        "complete_chromosomes": sum(1 for r in rows.values() if r["complete"]),
        "genome": {
            **tot,
            "fraction_with_target": round(tot["named"] / tot["scored"], 3) if tot["scored"] else None,
            "coding_target_agrees_with_nearest": round(agree / tot["coding_named"], 3)
            if tot["coding_named"]
            else None,
            "coding_target_inside_domain": round(inside / tot["coding_named"], 3)
            if tot["coding_named"]
            else None,
            "verdicts_coding": dict(sorted(verdicts.items(), key=lambda kv: -kv[1])),
        },
        "requests_total": requests,
        "cells": ["K562", "HepG2", "GM12878", "IMR-90"],
        "elements_where": "data/knowledge/alphagenome/all_elements/<chrom>.json (local)",
        "evidence": "predicted: AlphaGenome deletion effect per element, per gene and per cell line; "
        "inferred: CTCF-only nodes",
    }


def main() -> int:
    out = aggregate()
    save_result("enhancer_targets_all_genome_wide", out)
    g = out["genome"]
    print(
        f"{out['complete_chromosomes']} chromosomes complete: {g['scored']:,} elements, "
        f"{g['fraction_with_target']} name a gene, nearest TSS {g['coding_target_agrees_with_nearest']}, "
        f"inside the node {g['coding_target_inside_domain']}; {out['requests_total']:,} requests so far"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
