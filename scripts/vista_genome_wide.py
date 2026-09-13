# SPDX-License-Identifier: AGPL-3.0-or-later
"""VISTA's measured enhancers against the attribution, every chromosome, then one summary.

    uv run python scripts/vista_genome_wide.py          # fold the vista_chr*.json that exist
    uv run python scripts/vista_genome_wide.py --run    # score the missing chromosomes first

With --run this is the `vista_genome_wide` job: smallest chromosome first, one result per
chromosome as it lands, resumable, then `vista_genome_wide.json` with the positives against the
negatives over the whole set and per chromosome.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from genomeos.attribution import vista
from genomeos.jobs import CHROMOSOMES, heartbeat
from genomeos.results import save_result

ORDER = ["chr21", "chr22", "chr19", "chr20", "chr18", "chr17", "chr16", "chr15", "chr14", "chr13"]
ORDER += [c for c in CHROMOSOMES if c not in ORDER and c not in ("chrM", "chrY")][::-1]  # VISTA has no chrY


def run_missing(passes: int = 3) -> None:
    """Score every chromosome without a result; a chromosome that fails (a dropped connection) is retried
    on the next pass, up to `passes`."""
    for n in range(1, passes + 1):
        missing = [c for c in ORDER if not Path(f"data/results/vista_{c}.json").exists()]
        if not missing:
            return
        for i, chrom in enumerate(missing, 1):
            if Path(f"data/results/vista_{chrom}.json").exists():
                continue  # landed meanwhile (another worker)
            heartbeat("vista_genome_wide")
            print(f"{chrom}: VISTA elements scored ({i}/{len(missing)}, pass {n})", flush=True)
            r = subprocess.run(
                [sys.executable, "scripts/vista_score.py", "--chrom", chrom],
                check=False,
                env={**os.environ, "GENOMEOS_JOB": "vista_genome_wide"},
            )  # noqa: S603
            if r.returncode != 0:
                print(f"{chrom}: exit {r.returncode}, continuing with the next chromosome", flush=True)


def aggregate(results_dir: Path = Path("data/results")) -> dict:
    """Every vista_chr*.json folded: one summary over all rows, and the per-chromosome headline."""
    rows: list[dict] = []
    per: dict[str, dict] = {}
    cost = {"requests": 0, "mb_fetched": 0.0}
    for f in sorted(results_dir.glob("vista_chr*.json")):
        d = json.loads(f.read_text())
        for r in d["rows"]:
            r = dict(r)
            r["enhancer_like"] = any(c in ("dELS", "pELS") for c in r.get("ccre_classes") or [])
            r["any_ccre"] = bool(r.get("ccre_classes"))
            r["groups"] = sorted(vista.tissue_groups(r.get("tissues") or []))
            rows.append(r)
        s = d["summary"]
        per[d["chrom"]] = {
            "positive": s["positive"]["elements"],
            "negative": s["negative"]["elements"],
            "enhancer_like": [s["positive"]["enhancer_like"], s["negative"]["enhancer_like"]],
            "constrained": [s["positive"]["constrained"], s["negative"]["constrained"]],
            "fraction_with_target": [
                s["positive"]["fraction_with_target"],
                s["negative"]["fraction_with_target"],
            ],
            "tissue_agrees": s["positive"]["tissue_agrees"],
            "tissue_judged": s["positive"]["tissue_judged"],
        }
        c = d.get("phylop_cost") or {}
        cost["requests"] += c.get("requests", 0)
        cost["mb_fetched"] = round(cost["mb_fetched"] + c.get("mb_fetched", 0.0), 1)
    return {
        "chromosomes": per,
        "elements": len(rows),
        "genome": vista.summarise(rows),
        "phylop_cost": cost,
        "evidence": f"{vista.EVIDENCE}; measured: Zoonomia phyloP; curated: ENCODE cCREs; "
        "inferred: CTCF-only nodes; predicted: AlphaGenome deletion effect",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true", help="score the chromosomes without a result first")
    args = ap.parse_args()
    if args.run:
        run_missing()
    out = aggregate()
    save_result("vista_genome_wide", out)
    g = out["genome"]
    print(json.dumps({k: g[k] for k in ("positive", "negative", "separation")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
