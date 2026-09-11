# SPDX-License-Identifier: AGPL-3.0-or-later
"""Predicted enhancer targets on every chromosome (AlphaGenome feature b), one after another.

Runs scripts/enhancer_targets.py per chromosome, smallest first, skipping the ones whose sample is
complete; a chromosome that fails (quota, network) is retried after a wait and the loop only ends when
every chromosome has its sample. Per-element answers are cached, so nothing is ever scored twice.

    uv run python scripts/enhancer_targets_genome_wide.py [--sample 200]
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

from genomeos.jobs import heartbeat
from genomeos.results import load_result, save_result

# smallest first, so the first results arrive early; chrM has no enhancers
ORDER = [
    "chr21",
    "chr22",
    "chrY",
    "chr19",
    "chr20",
    "chr18",
    "chr17",
    "chr16",
    "chr15",
    "chr14",
    "chr13",
    "chr12",
]
ORDER += ["chr11", "chr10", "chr9", "chr8", "chrX", "chr7", "chr6", "chr5", "chr4", "chr3", "chr2", "chr1"]


def scored(chrom: str) -> int:
    r = load_result(f"enhancer_targets_{chrom}")
    return int((r or {}).get("summary", {}).get("elements_scored", 0))


def aggregate(chroms: list[str]) -> dict:
    """One summary over every chromosome scored: the node-model test at genome scale."""
    rows = {}
    tot = {"elements_scored": 0, "with_predicted_target": 0, "strong": 0, "silencer_like": 0}
    verdicts: dict[str, int] = {}
    tissues: dict[str, int] = {}
    for c in chroms:
        r = load_result(f"enhancer_targets_{c}") or {}
        s = r.get("summary", {})
        if not s:
            continue
        rows[c] = {
            "elements_scored": s["elements_scored"],
            "fraction_with_target": s["fraction_with_target"],
            "coding_target_agrees_with_nearest": s.get("coding_target_agrees_with_nearest"),
            "fraction_inside_domain": s.get("fraction_inside_domain"),
            "silencer_like": s.get("silencer_like", 0),
            "distal_enhancers": r.get("distal_enhancers_with_domain_target"),
        }
        for k in tot:
            tot[k] += s.get(k, 0)
        for k, n in s.get("verdicts_coding", {}).items():
            verdicts[k] = verdicts.get(k, 0) + n
        for k, n in s.get("top_tissues", {}).items():
            tissues[k] = tissues.get(k, 0) + n
    named_coding = sum(n for k, n in verdicts.items() if k != "no predicted effect")
    agree = verdicts.get("agrees with nearest TSS in domain", 0)
    inside = agree + verdicts.get("another gene in the same domain", 0)
    return {
        "chromosomes": rows,
        "genome": {
            "chromosomes": len(rows),
            **tot,
            "fraction_with_target": round(tot["with_predicted_target"] / max(1, tot["elements_scored"]), 3),
            "coding_targets_named": named_coding,
            "coding_target_agrees_with_nearest": round(agree / max(1, named_coding), 3),
            "coding_target_inside_domain": round(inside / max(1, named_coding), 3),
            "verdicts_coding": dict(sorted(verdicts.items(), key=lambda kv: -kv[1])),
            "top_tissues": dict(sorted(tissues.items(), key=lambda kv: -kv[1])[:15]),
        },
        "evidence": "predicted: AlphaGenome RNA-seq gene scorer on element deletions; inferred: CTCF-domain "
        "targets",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", type=int, default=200)
    args = ap.parse_args()
    waits = 0
    while True:
        todo = [c for c in ORDER if scored(c) < args.sample and (Path("data/reference") / f"{c}.fa").exists()]
        if not todo:
            done = [c for c in ORDER if scored(c) >= args.sample]
            summary = aggregate(done)
            save_result("enhancer_targets_genome_wide", summary)
            print(f"done: {len(done)} chromosomes with {args.sample} elements each", flush=True)
            print(json.dumps(summary["genome"], indent=1), flush=True)
            return 0
        chrom = todo[0]
        heartbeat("enhancer_targets_genome_wide")
        print(f"{chrom}: {scored(chrom)}/{args.sample} scored; {len(todo)} chromosomes to go", flush=True)
        rc = subprocess.call(  # noqa: S603
            [sys.executable, "scripts/enhancer_targets.py", "--chrom", chrom, "--sample", str(args.sample)]
        )
        if rc != 0 or scored(chrom) < args.sample:
            waits += 1
            wait = min(900, 60 * waits)
            print(f"{chrom}: not complete (exit {rc}); waiting {wait}s before the next attempt", flush=True)
            time.sleep(wait)
        else:
            waits = 0


if __name__ == "__main__":
    sys.exit(main())
