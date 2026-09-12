# SPDX-License-Identifier: AGPL-3.0-or-later
"""lentiMPRA activity against GenomeOS, every chromosome, then one summary over all elements.

    uv run python scripts/mpra_genome_wide.py          # fold the mpra_chr*.json that exist
    uv run python scripts/mpra_genome_wide.py --run    # score the missing chromosomes first

With --run this is the `mpra_genome_wide` job: smallest chromosome first, one result per chromosome
as it lands, resumable, then `mpra_genome_wide.json` (rows are not repeated there; the per-chromosome
files carry them).
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

from genomeos.attribution import mpra
from genomeos.jobs import CHROMOSOMES, heartbeat
from genomeos.results import save_result

ORDER = ["chr21", "chr22", "chrY", "chr19", "chr20", "chr18", "chr17", "chr16", "chr15", "chr14", "chr13"]
ORDER += [c for c in CHROMOSOMES if c not in ORDER and c != "chrM"][::-1]


def run_missing(passes: int = 3) -> None:
    """Score every chromosome without a result; a chromosome that fails (a dropped connection) is retried
    on the next pass, up to `passes`."""
    for n in range(1, passes + 1):
        missing = [c for c in ORDER if not Path(f"data/results/mpra_{c}.json").exists()]
        if not missing:
            return
        for i, chrom in enumerate(missing, 1):
            if Path(f"data/results/mpra_{chrom}.json").exists():
                continue  # landed meanwhile (another worker)
            heartbeat("mpra_genome_wide")
            print(f"{chrom}: lentiMPRA elements scored ({i}/{len(missing)}, pass {n})", flush=True)
            r = subprocess.run(
                [sys.executable, "scripts/mpra_score.py", "--chrom", chrom],
                check=False,
                env={**os.environ, "GENOMEOS_JOB": "mpra_genome_wide"},
            )  # noqa: S603
            if r.returncode != 0:
                print(f"{chrom}: exit {r.returncode}, continuing with the next chromosome", flush=True)


def aggregate(results_dir: Path = Path("data/results")) -> dict:
    rows: list[dict] = []
    per: dict[str, dict] = {}
    requests = 0
    for f in sorted(results_dir.glob("mpra_chr*.json")):
        d = json.loads(f.read_text())
        rows.extend(d["rows"])
        s = d["summary"]
        per[d["chrom"]] = {
            "elements": d["annotated"],
            "active": {c: v["active"] for c, v in s["cells"].items()},
            "predicted_rho_same_cell": {
                c: ((v.get("predicted_dnase_vs_activity") or {}).get("same_cell") or {}).get("spearman")
                for c, v in s["cells"].items()
                if c in mpra.MODELLED
            },
        }
        requests += (d.get("model_cost") or {}).get("requests", 0)
    return {
        "chromosomes": per,
        "elements": len(rows),
        "genome": mpra.summarise(rows),
        "model_requests": requests,
        "evidence": f"{mpra.EVIDENCE}; experimental: ENCODE DNase peaks (reader); measured: Zoonomia phyloP; "
        "curated: ENCODE cCREs; predicted: AlphaGenome DNase per cell line",
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="store_true")
    args = ap.parse_args()
    if args.run:
        run_missing()
    out = aggregate()
    save_result("mpra_genome_wide", out)
    g = out["genome"]
    for cell, s in g["cells"].items():
        print(cell, json.dumps({k: s.get(k) for k in ("elements", "active", "predicted_dnase_vs_activity")}))
    print("specific:", json.dumps(g.get("specific_elements")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
