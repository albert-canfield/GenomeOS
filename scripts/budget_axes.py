# SPDX-License-Identifier: AGPL-3.0-or-later
"""The stored composition budget restated under review R7, one chromosome at a time, no requests.

    uv run python scripts/budget_axes.py [--chroms chr21 ...]

Each `budget_<chrom>` keeps its values. `budget_axes_<chrom>` holds the same blocks, constraint and
scores with the classifier's RepeatMasker origin read for every block (the regulatory ones included),
the constraint reading as `evidence_status`, and tiers and labels that read no function into missing
constraint; `budget_axes_genome_wide` sums them.
"""

from __future__ import annotations

import argparse

from genomeos import manifest as mf
from genomeos.attribution.budget import AXES_NAME, AXES_TIERS, LEGACY_TIER, distil, restate_and_save
from genomeos.results import RESULTS_DIR, save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", nargs="*", default=CHROMS)
    args = ap.parse_args()
    for c in args.chroms:
        o = restate_and_save(c)
        print(c, len(o["blocks"]), {t: v["blocks"] for t, v in o["by_tier"].items()}, flush=True)
    s = distil(name=AXES_NAME, tiers=AXES_TIERS)
    s["note"] = (
        "The R7 restatement of budget_genome_wide: the same blocks, constraint and scores; "
        "repeat_unconstrained "
        "was fossil and unconstrained_unknown was neutral, and constraint is evidence of selection only. "
        + s["note"]
    )
    paths = [RESULTS_DIR / f"{AXES_NAME}_{c}.json" for c in CHROMS]
    save_result(
        f"{AXES_NAME}_genome_wide",
        s,
        manifest={
            "sources": [{"accession": "this repository, budget_axes_<chrom>", "version": "sha256"}],
            "inputs": [mf.input_entry(p, partition=None) for p in paths if p.exists()],
            "assembly": "GRCh38",
            "coordinates": {"base": 0, "interval": "half-open"},
            "parameters": {"tiers": list(AXES_TIERS), "legacy_tier": LEGACY_TIER, "chromosomes": CHROMS},
            "exclusions": [],
            "partitions": "n/a: a sum of per-chromosome rules, no evaluation",
        },
    )
    print({t: v["blocks"] for t, v in s["by_tier"].items()}, "guessed", s["guessed_fraction"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
