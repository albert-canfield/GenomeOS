# SPDX-License-Identifier: AGPL-3.0-or-later
"""One summary over every chromosome's constrained-enhancer targets, next to the uniform samples.

uv run python scripts/constrained_targets_genome_wide.py
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

from genomeos.jobs import CHROMOSOMES, heartbeat
from genomeos.results import save_result

ORDER = ["chr21", "chr22", "chrY", "chr19", "chr20", "chr18", "chr17", "chr16", "chr15", "chr14", "chr13"]
ORDER += [c for c in CHROMOSOMES if c not in ORDER and c != "chrM"][::-1]


def run_missing(top: int) -> None:
    """Score every chromosome without a result, one subprocess each, heartbeat per chromosome."""
    for i, chrom in enumerate(ORDER, 1):
        if Path(f"data/results/constrained_targets_{chrom}.json").exists():
            continue
        heartbeat(f"{chrom}: scoring the {top} most constrained distal enhancers ({i}/{len(ORDER)})")
        cmd = [sys.executable, "scripts/constrained_targets.py", "--chrom", chrom, "--top", str(top)]
        r = subprocess.run(cmd, check=False)  # noqa: S603 - our own script
        if r.returncode != 0:
            heartbeat(f"{chrom}: exit {r.returncode}, continuing with the next chromosome")


def aggregate(results_dir: Path = Path("data/results")) -> dict:
    """Every constrained_targets_chr*.json under results_dir folded into one genome-wide summary."""
    rows = {}
    tot = {"scored": 0, "named": 0, "strong": 0, "silencer_like": 0, "available": 0, "distal": 0}
    verdicts: dict[str, int] = {}
    uni = {"scored": 0, "named": 0, "strong": 0, "silencer_like": 0, "agree": 0.0, "inside": 0.0, "n": 0}
    for f in sorted(results_dir.glob("constrained_targets_chr*.json")):
        d = json.loads(f.read_text())
        s = d["summary"]
        c = d["chrom"]
        rows[c] = {
            "distal_enhancers": d["distal_enhancers"],
            "constrained_available": d["constrained_elements_available"],
            "scored": s["elements_scored"],
            "fraction_with_target": s["fraction_with_target"],
            "strong": s["strong"],
            "silencer_like": s["silencer_like"],
            "coding_target_agrees_with_nearest": s.get("coding_target_agrees_with_nearest"),
            "fraction_inside_domain": s.get("fraction_inside_domain"),
            "uniform_fraction_with_target": (d.get("uniform_sample") or {}).get("fraction_with_target"),
            "uniform_fraction_inside_domain": (d.get("uniform_sample") or {}).get("fraction_inside_domain"),
        }
        tot["scored"] += s["elements_scored"]
        tot["named"] += s["with_predicted_target"]
        tot["strong"] += s["strong"]
        tot["silencer_like"] += s["silencer_like"]
        tot["available"] += d["constrained_elements_available"]
        tot["distal"] += d["distal_enhancers"]
        for k, n in s.get("verdicts_coding", {}).items():
            verdicts[k] = verdicts.get(k, 0) + n
        u = d.get("uniform_sample") or {}
        if u.get("elements_scored"):
            uni["scored"] += u["elements_scored"]
            uni["named"] += round(u["elements_scored"] * (u.get("fraction_with_target") or 0))
            uni["strong"] += u.get("strong") or 0
            uni["silencer_like"] += u.get("silencer_like") or 0
            uni["agree"] += u.get("coding_target_agrees_with_nearest") or 0
            uni["inside"] += u.get("fraction_inside_domain") or 0
            uni["n"] += 1
    named_coding = sum(n for k, n in verdicts.items() if k != "no predicted effect")
    agree = verdicts.get("agrees with nearest TSS in domain", 0)
    inside = agree + verdicts.get("another gene in the same domain", 0)
    out = {
        "chromosomes": rows,
        "genome": {
            "chromosomes": len(rows),
            "distal_enhancers": tot["distal"],
            "constrained_elements_available": tot["available"],
            "scored": tot["scored"],
            "with_predicted_target": tot["named"],
            "fraction_with_target": round(tot["named"] / max(1, tot["scored"]), 3),
            "strong": tot["strong"],
            "silencer_like": tot["silencer_like"],
            "coding_targets_named": named_coding,
            "coding_target_agrees_with_nearest": round(agree / max(1, named_coding), 3),
            "coding_target_inside_domain": round(inside / max(1, named_coding), 3),
            "verdicts_coding": dict(sorted(verdicts.items(), key=lambda kv: -kv[1])),
        },
        "uniform_samples": {
            "scored": uni["scored"],
            "fraction_with_target": round(uni["named"] / max(1, uni["scored"]), 3),
            "strong": uni["strong"],
            "silencer_like": uni["silencer_like"],
            "coding_target_agrees_with_nearest": round(uni["agree"] / max(1, uni["n"]), 3),
            "coding_target_inside_domain": round(uni["inside"] / max(1, uni["n"]), 3),
        },
        "evidence": "measured: Zoonomia phyloP constraint; "
        "predicted: AlphaGenome deletion effect per element; inferred: CTCF-only nodes",
    }
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="score the chromosomes without a result first")
    ap.add_argument("--top", type=int, default=100)
    args = ap.parse_args()
    if args.run:
        run_missing(args.top)
    out = aggregate()
    save_result("constrained_targets_genome_wide", out)
    print(json.dumps(out["genome"], indent=1))
    print("uniform:", json.dumps(out["uniform_samples"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
