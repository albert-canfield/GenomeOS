# SPDX-License-Identifier: AGPL-3.0-or-later
"""VISTA enhancers across species: which factor families do limb and neural enhancers keep in place
more often than conserved elements that drove no expression?

    uv run python scripts/across_panel.py [--per-group 40]

Each locus is aligned to mouse, opossum, chicken, frog and zebrafish (Ensembl Compara LASTZ), scanned
with JASPAR, and its factor-strict held sites recorded (knowledge/across.py). Per-locus records are
cached, so the run resumes. The result is `across_panel_vista.json`.
"""

from __future__ import annotations

import argparse
import time

from genomeos.genome.motifs import load_motifs
from genomeos.jobs import heartbeat
from genomeos.knowledge.across import panel_record, panel_summary, select_panel, vista_loci
from genomeos.results import save_result

JOB = "across_panel_vista"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-group", type=int, default=40)
    ap.add_argument("--seed", type=int, default=5)
    args = ap.parse_args()
    loci = select_panel(vista_loci(), per_group=args.per_group, seed=args.seed)
    motifs = load_motifs()
    print(
        f"{len(loci)} loci: "
        + ", ".join(f"{g} {sum(1 for x in loci if x['group'] == g)}" for g in ("limb", "neural", "negative")),
        flush=True,
    )
    records = []
    t0 = time.time()
    for i, loc in enumerate(loci, 1):
        try:
            records.append(panel_record(loc, motifs))
        except (OSError, ValueError, KeyError) as e:
            print(f"  {loc['id']}: skipped ({str(e)[:80]})", flush=True)
            continue
        heartbeat(JOB)
        if i % 10 == 0:
            print(f"  {i}/{len(loci)} loci, {time.time() - t0:.0f} s", flush=True)
    s = panel_summary(records)
    s["records"] = records
    s["seed"] = args.seed
    s["evidence"] = (
        "curated: VISTA transgenic outcomes, Ensembl Compara LASTZ alignments; predicted: JASPAR hits; "
        "inferred: families held at factor-strict sites, positives against constraint-matched negatives, "
        "one-sided Fisher with Benjamini–Hochberg"
    )
    save_result(JOB, s)
    for grp, v in s["groups"].items():
        print(
            f"{grp}: {v['positives']} positives vs {v['negatives']} negatives; held sites per kb median "
            f"{v['density_median_positive']:.1f} vs {v['density_median_negative']:.1f} "
            f"(p {v['density_p_greater']}); "
            f"families at q<=0.05: {', '.join(r['family'] for r in v['families_q05']) or 'none'}",
            flush=True,
        )


if __name__ == "__main__":
    main()
