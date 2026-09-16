# SPDX-License-Identifier: AGPL-3.0-or-later
"""Is the human panel's matched background a neutral baseline? Read its own results and ask.

    uv run python scripts/panel_tier_baseline.py

Every tier ratio the panel reports is a tier's fixed rate over the rate of GC- and
replication-timing-matched background windows, so a ratio of 1 is meant to say "this tier behaves like
background". Seven of 23 chromosomes fail a control, and on all six autosomes and X the failing check
is the neutral tier missing that background. This asks whether those are six special chromosomes or
the tail of something that applies to every chromosome, using only the committed per-chromosome
results: no alignment is streamed and no model is called.

It reports, for each tier, on how many chromosomes the ratio sits above 1 with a sign test, the median
ratio, and every tier ratio divided by the same chromosome's neutral ratio -- the reading the panel
lane arrived at independently, that a tier means something beside the other unknown sequence rather
than beside the background. chrY is excluded, as everywhere else in the panel: it is a gene-poor
chromosome whose coding checks fail by construction.

Writes data/results/panel_tier_baseline.json.
"""

from __future__ import annotations

import glob
import json
from math import comb
from pathlib import Path
from statistics import median

from genomeos.results import save_result

TIERS = ("cds", "fossil", "regulatory", "constrained_unknown", "neutral")
EXCLUDED = ("chrY",)
RESULTS = Path("data/results")


def sign_test_above_one(values: list[float]) -> dict[str, float | int]:
    """How many of these sit above 1, and how surprising that is if each were a coin toss."""
    n = len(values)
    above = sum(v > 1 for v in values)
    p = sum(comb(n, k) for k in range(above, n + 1)) / 2**n
    return {"chromosomes": n, "above_one": above, "one_sided_p": round(p, 6)}


def read_chromosomes() -> dict[str, dict[str, float]]:
    """Every committed per-chromosome panel result, as tier -> ratio against the matched background."""
    out = {}
    for path in sorted(glob.glob(str(RESULTS / "human_panel_chr*.json"))):
        chrom = Path(path).stem.split("_")[-1]
        if chrom in EXCLUDED:
            continue
        pooled = json.loads(Path(path).read_text())["pooled"]
        out[chrom] = {t: pooled[t]["ratio_gc_rt"] for t in TIERS if t in pooled}
    return out


def collect() -> dict[str, object]:
    by_chrom = read_chromosomes()
    chroms = sorted(by_chrom)
    against_background, against_neutral = {}, {}
    for tier in TIERS:
        ratios = [by_chrom[c][tier] for c in chroms if tier in by_chrom[c]]
        if not ratios:  # a tier no result carries: say nothing about it rather than divide by nothing
            continue
        against_background[tier] = {"median": round(median(ratios), 3), **sign_test_above_one(ratios)}
        if tier == "neutral":
            continue
        relative = [by_chrom[c][tier] / by_chrom[c]["neutral"] for c in chroms if tier in by_chrom[c]]
        if not relative:
            continue
        against_neutral[tier] = {
            "median": round(median(relative), 3),
            "below_neutral": sum(v < 1 for v in relative),
            "chromosomes": len(relative),
        }
    return {
        "result": "panel_tier_baseline",
        "chromosomes": chroms,
        "excluded": list(EXCLUDED),
        "against_matched_background": against_background,
        "against_the_neutral_tier": against_neutral,
        "per_chromosome": by_chrom,
        "reading": (
            "every non-coding tier reads above the matched background on almost every chromosome, so "
            "the background is not a neutral baseline: it fixes more than the sequence it is the "
            "control for. A tier's distance from 1 therefore carries that offset, and the comparison "
            "that does not is a tier against the neutral tier of the same chromosome."
        ),
    }


def main() -> int:
    result = collect()
    save_result(result["result"], result)
    bg = result["against_matched_background"]
    assert isinstance(bg, dict)
    for tier, row in bg.items():
        print(
            f"{tier:20} above 1 on {row['above_one']:2} of {row['chromosomes']} "
            f"(p {row['one_sided_p']}), median {row['median']}"
        )
    print()
    rel = result["against_the_neutral_tier"]
    assert isinstance(rel, dict)
    for tier, row in rel.items():
        print(
            f"{tier:20} against the neutral tier: median {row['median']}, "
            f"below it on {row['below_neutral']} of {row['chromosomes']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
