# SPDX-License-Identifier: AGPL-3.0-or-later
"""What is the panel's matched background made of, and how much of its offset does that explain?

    uv run python scripts/panel_background_composition.py                 # every chromosome with a store
    uv run python scripts/panel_background_composition.py --chrom chr21   # one, for a quick read
    uv run python scripts/panel_background_composition.py --phylop 400    # plus sampled phyloP

`panel_tier_baseline.py` established on 2026-09-16 that the GC- and replication-timing-matched
background is not a neutral baseline: every non-coding tier reads above it, a systematic +6 to +11
points that belongs to the control. It did not measure what the background holds. This does.

Each chromosome's background is rebuilt exactly as `human_panel.build_background` builds it, every
kilobase is annotated from local data (GENCODE coding bases, phastCons conserved elements, segmental
duplication, RepeatMasker classes, ENCODE cCREs, distance to a coding TSS, GC and replication
timing), and the offset is decomposed by rebuilding the background with each candidate removed. The
covariate shares, the rates and the ratios are reported side by side, with the denominator beside
every one: a covariate that was measured and explains nothing and a covariate that could not be
measured are different rows, never the same zero.

chrY is measured and reported but kept out of the medians, as everywhere else in the panel.

Writes data/results/panel_background_composition.json.
"""

from __future__ import annotations

import argparse
import glob
import json
from math import comb
from pathlib import Path
from statistics import median
from typing import Any

from genomeos.attribution.human_panel import CACHE, TIERS
from genomeos.attribution.panel_background import (
    CANDIDATES,
    CLAIM_WITH_COVERAGE,
    CLAIM_WITHOUT_COVERAGE,
    COVERAGE_KEY,
    COVERAGE_MEASURE,
    GROUPS,
    SOURCES,
    analyse,
    reads_explained,
)
from genomeos.results import RESULTS_DIR, save_result

EXCLUDED = ("chrY",)
KEYS = (
    "as_the_panel_builds_it",
    *(f"without_{c}" for c in CANDIDATES),
    "without_any_of_them",
    COVERAGE_KEY,
)
STRATIFICATIONS = ("gc_and_timing", "gc_timing_and_tss", "gc_timing_tss_and_coverage")


def chromosomes(results_dir: Path = RESULTS_DIR, cache: Path = CACHE) -> list[str]:
    """Every chromosome with both a committed panel result and a local store to rebuild from."""
    out = []
    for p in sorted(glob.glob(str(results_dir / "human_panel_chr*.json"))):
        c = Path(p).stem.split("_")[-1]
        if (cache / c / "meta.json").exists():
            out.append(c)
    # smallest store first: a long run that has to be stopped has read the cheap chromosomes already
    return sorted(out, key=lambda c: (cache / c / "sites.tsv.gz").stat().st_size)


def sign_test_above_one(values: list[float]) -> dict[str, Any]:
    n = len(values)
    above = sum(v > 1 for v in values)
    return {
        "chromosomes": n,
        "above_one": above,
        "one_sided_p": round(sum(comb(n, k) for k in range(above, n + 1)) / 2**n, 6) if n else None,
    }


def across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Every tier's ratio under every exclusion, over the chromosomes where that exclusion was made.

    The explained share is computed on the *same* chromosomes as the ratio it is compared with, so a
    covariate assessed on fewer chromosomes cannot borrow the others' offset.
    """
    out: dict[str, Any] = {}
    for tier in TIERS:
        row: dict[str, Any] = {}
        for key in KEYS:
            pairs = [
                (r["offset"]["as_the_panel_builds_it"].get(tier), _ratio(r["offset"][key], tier))
                for r in runs
                if tier in r["offset"]["as_the_panel_builds_it"] and _assessed(r["offset"][key])
            ]
            pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
            if not pairs:
                row[key] = {"chromosomes_assessed": 0, "median": None, "explains_of_the_offset": None}
                continue
            base, here = median(a for a, _ in pairs), median(b for _, b in pairs)
            share = round((base - here) / (base - 1), 4) if base != 1 else None
            row[key] = {
                "chromosomes_assessed": len(pairs),
                "median": round(here, 4),
                "median_as_built_on_the_same_chromosomes": round(base, 4),
                "explains_of_the_offset": share,
                "reads": reads_explained(share) if key != "as_the_panel_builds_it" else "the offset itself",
                **sign_test_above_one([b for _, b in pairs]),
            }
        out[tier] = row
    return out


def _assessed(entry: Any) -> bool:
    return isinstance(entry, dict) and entry.get("assessed", True)


def _ratio(entry: Any, tier: str) -> float | None:
    return entry.get(tier) if isinstance(entry, dict) else None


def composition_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Median feature share per group, with the number of chromosomes each share was assessed on."""
    out: dict[str, Any] = {}
    groups = ["background", "background_outside_every_block", *GROUPS]
    for g in groups:
        rows = [r["composition"][g] for r in runs if r["composition"].get(g, {}).get("bases")]
        if not rows:
            out[g] = {"chromosomes_assessed": 0}
            continue
        entry: dict[str, Any] = {
            "chromosomes_assessed": len(rows),
            "bases_total": sum(r["bases"] for r in rows),
            "gc_median": round(median(r["gc_median"] for r in rows if r.get("gc_median")), 4),
            "tss_median": int(median(r["tss_median"] for r in rows)),
            "tss_assessed_chromosomes": sum(1 for r in rows if r.get("tss_assessed")),
            # the background's unit is the kilobase the panel bins by; a tier's is its block
            "length_median": int(median(r.get("length_median", 1000) for r in rows)),
            "coverage_measure": COVERAGE_MEASURE,
            "coverage_median": round(
                median(r["coverage_median"] for r in rows if r.get("coverage_median") is not None), 4
            ),
            "coverage_assessed_chromosomes": sum(1 for r in rows if r.get("coverage_median") is not None),
            "recurring_per_kb_median": round(median(r["recurring_per_kb"] for r in rows), 3),
            "features": {},
        }
        for f in rows[0]["features"]:
            vals = [r["features"][f]["share"] for r in rows if r["features"][f]["share"] is not None]
            entry["features"][f] = {
                "share_median": round(median(vals), 4) if vals else None,
                "chromosomes_assessed": len(vals),
                "reads": "not assessed on any chromosome"
                if not vals
                else f"measured on {len(vals)}, median {round(100 * median(vals), 2)}% of the bases",
            }
        out[g] = entry
    return out


def standardised_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """The tier-against-background window comparison, pooled as medians of the per-chromosome runs."""
    out: dict[str, Any] = {}
    for tier in TIERS:
        rows = [r["standardised"][tier] for r in runs if tier in r["standardised"]]
        if not rows:
            out[tier] = {"chromosomes_assessed": 0}
            continue
        last = STRATIFICATIONS[-1]
        out[tier] = {
            "chromosomes_assessed": len(rows),
            "targets": sum(r[last]["targets"] for r in rows),
            "coverage_measure": COVERAGE_MEASURE,
            "matched_difference": {
                k: round(median(r[k]["matched"]["difference"] for r in rows), 4) for k in STRATIFICATIONS
            },
            "dropped_for_want_of_a_control": {
                k: sum(r[k]["dropped_for_want_of_a_control"] for r in rows) for k in STRATIFICATIONS
            },
            "targets_matched": {k: sum(r[k]["targets_matched"] for r in rows) for k in STRATIFICATIONS},
            "imbalance_medians": {
                axis: {
                    "target": round(median(r[last]["imbalance"][axis]["target_median"] for r in rows), 4),
                    "control": round(median(r[last]["imbalance"][axis]["control_median"] for r in rows), 4),
                }
                for axis in ("length", "gc", "tss", "coverage")
            },
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chrom", action="append", help="one chromosome; repeatable (default: all with a store)")
    ap.add_argument("--phylop", type=int, default=0, help="kilobases of background to sample for phyloP")
    ap.add_argument("--name", default="panel_background_composition")
    args = ap.parse_args(argv)
    chroms = args.chrom or chromosomes()
    runs = []
    for c in chroms:
        print(f"{c} ...", flush=True)
        r = analyse(c, phylop=args.phylop)
        b = r["composition"]["background"]
        print(
            f"  {b['kilobases']:>7} kb measured, {b['recurring_per_kb']}/kb, "
            f"exons {b['features']['exon_any']['share']}, conserved "
            f"{b['features']['conserved_elements']['share']}, "
            f"segdup {b['features']['segmental_duplication']['share']}",
            flush=True,
        )
        runs.append(r)
    kept = [r for r in runs if r["chrom"] not in EXCLUDED]
    conditioned = all(r["claim_available"] == CLAIM_WITH_COVERAGE for r in kept) if kept else False
    rebuilt_arm = sum(1 for r in kept if r["coverage_conditioned"]["by_rebuilding_the_background"])
    result = {
        "chromosomes": [r["chrom"] for r in runs],
        "excluded_from_the_medians": list(EXCLUDED),
        "evidence": SOURCES,
        "coverage_measure": COVERAGE_MEASURE,
        # one of two sentences, never a third, and carried as data so that no reader can upgrade it
        "claim_available": CLAIM_WITH_COVERAGE if conditioned else CLAIM_WITHOUT_COVERAGE,
        "chromosomes_with_coverage_conditioned": {
            "in_the_window_comparison": sum(
                1 for r in kept if r["coverage_conditioned"]["in_the_window_comparison"]
            ),
            "by_rebuilding_the_background": rebuilt_arm,
            "chromosomes": len(kept),
        },
        "reproduces_the_panel": {
            "chromosomes_checked": sum(1 for r in runs if r["reproduces_the_panel"]["tiers_checked"]),
            "largest_absolute_difference": max(
                (
                    r["reproduces_the_panel"]["largest_absolute_difference"]
                    for r in runs
                    if r["reproduces_the_panel"]["largest_absolute_difference"] is not None
                ),
                default=None,
            ),
        },
        "across_chromosomes": across(kept),
        "composition": composition_across(kept),
        "standardised": standardised_across(kept),
        "per_chromosome": {
            r["chrom"]: {
                k: r[k]
                for k in (
                    "reproduces_the_panel",
                    "composition",
                    "feature_splits",
                    "offset",
                    "explained",
                    "constrained_bases",
                    "standardised",
                    "coverage",
                )
            }
            for r in runs
        },
    }
    save_result(args.name, result)
    for tier in TIERS:
        row = result["across_chromosomes"][tier]
        base = row["as_the_panel_builds_it"]
        print(f"\n{tier}: median {base['median']} above 1 on {base['above_one']} of {base['chromosomes']}")
        for key in KEYS[1:]:
            r = row[key]
            print(
                f"  {key:34} median {r['median']} (assessed on {r['chromosomes_assessed']}) {r.get('reads')}"
            )
    print("\n" + result["claim_available"])
    print(json.dumps({"result": args.name, "chromosomes": len(runs)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
