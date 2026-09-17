# SPDX-License-Identifier: AGPL-3.0-or-later
"""All four covariates out of the panel's background at once: the union arm the two lanes left undone.

    uv run python scripts/panel_union_arms.py                 # every chromosome with a store
    uv run python scripts/panel_union_arms.py --chrom chr21   # one, for a quick read
    uv run python scripts/panel_union_arms.py --from-result   # re-distil, measuring nothing

`panel_background_composition.py` took exons, conserved elements and promoter proximity out of the
panel's matched background -- 41% (fossil), 56% (regulatory), 25% (neutral). `panel_leftover_classes.py`
took coding-gene introns out -- 56% / 66% / 30%, more than the other three together. The two sets
overlap and are not additive, and their union was never measured. This measures it: the background
rebuilt with all four gone, beside each single arm and beside the sum of the singles, on the same
chromosomes and by the same method.

The pairwise overlap between the covariates is reported first, because it is cheap and it is what
decides whether the union lands near the largest single or near the sum. Every arm carries the share
of the remaining background that is block sequence and the `degenerate` flag the earlier lanes wrote,
and a degenerate arm is excluded from every reported share: an exclusion that leaves a background made
of blocks has turned the comparison into a tier average against a tier.

chrY is measured and reported but kept out of the medians, as everywhere else in the panel.

Writes data/results/panel_union_arms.json.
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
from genomeos.attribution.panel_background import reads_explained
from genomeos.attribution.panel_union import (
    ARMS,
    CLAIM_WITH_COVERAGE,
    CLAIM_WITHOUT_COVERAGE,
    COVARIATES,
    COVERAGE_MEASURE,
    FROM_LANE,
    NEAR,
    SOURCES,
    STRATIFICATIONS,
    UNION,
    analyse,
    where_it_sits,
)
from genomeos.results import RESULTS_DIR, load_result, save_result

EXCLUDED = ("chrY",)
LAST = STRATIFICATIONS[-1]
UNION_ARM = f"without_{UNION}"
COMPARISONS = (
    "the_covariates_against_the_covariate_free_background",
    *(f"{t}_against_the_covariate_free_background" for t in TIERS),
)
PER_CHROMOSOME = (
    "chromosome_bases",
    "claim_available",
    "reproduces_the_panel",
    "eligible",
    "overlap",
    "groups",
    "feature_splits",
    "offset",
    "explained",
    "where_the_union_sits",
    "comparisons",
    "coverage",
)


def chromosomes(results_dir: Path = RESULTS_DIR, cache: Path = CACHE) -> list[str]:
    """Every chromosome with a committed panel result and a local store to rebuild from."""
    out = []
    for p in sorted(glob.glob(str(results_dir / "human_panel_chr*.json"))):
        c = Path(p).stem.split("_")[-1]
        if (cache / c / "meta.json").exists():
            out.append(c)
    # smallest store first: a long run that has to be stopped has read the cheap chromosomes already
    return sorted(out, key=lambda c: (cache / c / "sites.tsv.gz").stat().st_size)


def sign_test(values: list[float], above: float = 0.0) -> dict[str, Any]:
    n = len(values)
    k = sum(v > above for v in values)
    return {
        "chromosomes": n,
        "above": k,
        "one_sided_p": round(sum(comb(n, i) for i in range(k, n + 1)) / 2**n, 6) if n else None,
    }


def _median(values: list[Any]) -> float | None:
    """The median of the values that exist, or None: a covariate absent everywhere is not a zero."""
    kept = [v for v in values if v is not None]
    return round(median(kept), 4) if kept else None


def overlap_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """The pairwise overlaps summed over the chromosomes, with both denominators beside each."""
    pairs = sorted({k for r in runs for k in r["overlap"]["pairs"]})
    out: dict[str, Any] = {}
    for key in pairs:
        rows = [
            r["overlap"]["pairs"][key] for r in runs if r["overlap"]["pairs"][key]["shared_bases"] is not None
        ]
        if not rows:
            out[key] = {"chromosomes_assessed": 0, "reads": "not assessed on any chromosome"}
            continue
        shared = sum(x["shared_bases"] for x in rows)
        a, b = sum(x["bases_a"] for x in rows), sum(x["bases_b"] for x in rows)
        out[key] = {
            "chromosomes_assessed": len(rows),
            "shared_bases": shared,
            "bases_a": a,
            "bases_b": b,
            "share_of_a": round(shared / a, 4) if a else None,
            "share_of_b": round(shared / b, 4) if b else None,
            "share_of_the_smaller": round(shared / min(a, b), 4) if min(a, b) else None,
        }
    union = sum(r["overlap"]["union_bases"] for r in runs)
    summed = sum(r["overlap"]["sum_of_the_singles_bases"] for r in runs)
    return {
        "chromosomes": len(runs),
        "pairs": out,
        "union_bases": union,
        "sum_of_the_singles_bases": summed,
        "bases_counted_more_than_once_by_the_sum": summed - union,
        "the_sum_over_the_union": round(summed / union, 4) if union else None,
        "chromosome_bases": sum(r["overlap"]["chromosome_bases"] for r in runs),
    }


def groups_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Median covariates per group. Length, GC, distance to a coding TSS and coverage, for every one."""
    out: dict[str, Any] = {}
    for name in sorted({k for r in runs for k in r["groups"]}):
        rows = [r["groups"][name] for r in runs if r["groups"].get(name, {}).get("bases")]
        rows = [x for x in rows if x.get("gc_median") is not None]
        if not rows:
            out[name] = {"chromosomes_assessed": 0, "why": "no chromosome carries this group"}
            continue
        out[name] = {
            "chromosomes_assessed": len(rows),
            "units": sum(x.get("kilobases", x.get("blocks", 0)) for x in rows),
            "bases": sum(x["bases"] for x in rows),
            # the background's unit is the kilobase the panel bins by; a tier's is its block
            "length_median": int(median(x.get("length_median", 1000) for x in rows)),
            "gc_median": round(median(x["gc_median"] for x in rows), 4),
            "tss_median": int(median(x["tss_median"] for x in rows)),
            "tss_assessed_chromosomes": sum(1 for x in rows if x.get("tss_assessed")),
            "coverage_measure": COVERAGE_MEASURE,
            "coverage_median": _median([x.get("coverage_median") for x in rows]),
            "coverage_assessed_chromosomes": sum(1 for x in rows if x.get("coverage_median") is not None),
            "recurring_per_kb_median": round(median(x["recurring_per_kb"] for x in rows), 3),
        }
    return out


def offset_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Every tier's ratio under every arm, on the chromosomes where that arm was made and assessed.

    The explained share is computed on the SAME chromosomes as the ratio it is compared with, so an
    arm assessed on fewer chromosomes cannot borrow the others' offset. A degenerate arm keeps its
    numbers and its flag and is never given an explained share.
    """
    out: dict[str, Any] = {}
    for tier in TIERS:
        row: dict[str, Any] = {}
        built = [
            r["offset"]["as_the_panel_builds_it"][tier]
            for r in runs
            if r["offset"]["as_the_panel_builds_it"].get(tier) is not None
        ]
        row["as_the_panel_builds_it"] = {
            "chromosomes_assessed": len(built),
            "median": round(median(built), 4) if built else None,
            **sign_test(built, 1.0),
        }
        for key in ARMS:
            arm = [r["offset"].get(key, {}) for r in runs]
            degenerate = sum(1 for a in arm if a.get("degenerate"))
            block_shares = [a.get("block_share_of_the_remaining_background") for a in arm]
            pairs = [
                (r["offset"]["as_the_panel_builds_it"].get(tier), r["offset"][key].get(tier))
                for r in runs
                if r["offset"].get(key, {}).get("assessed")
            ]
            pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
            entry: dict[str, Any] = {
                "chromosomes_assessed": len(pairs),
                "chromosomes_degenerate": degenerate,
                "block_share_of_the_remaining_background_median": _median(block_shares),
                "block_share_of_the_remaining_background_max": (
                    max([v for v in block_shares if v is not None], default=None)
                ),
            }
            if not pairs:
                entry.update(
                    {
                        "median": None,
                        "explains_of_the_offset": None,
                        "reads": "not assessed",
                        "why": "the arm is degenerate on every chromosome"
                        if degenerate
                        else "the covariate has no bases on any chromosome",
                    }
                )
                row[key] = entry
                continue
            base, here = median(a for a, _ in pairs), median(b for _, b in pairs)
            share = round((base - here) / (base - 1), 4) if base != 1 else None
            entry.update(
                {
                    "median": round(here, 4),
                    "median_as_built_on_the_same_chromosomes": round(base, 4),
                    "explains_of_the_offset": share,
                    "reads": reads_explained(share),
                    **sign_test([b for _, b in pairs], 1.0),
                }
            )
            row[key] = entry
        out[tier] = row
    return out


def union_across(offsets: dict[str, Any]) -> dict[str, Any]:
    """Where the union sits between the largest single and the sum, from the pooled medians.

    `panel_union.where_it_sits` is given the pooled shares in the shape it reads per chromosome, so
    the same rule decides the same sentence at both scales and there is no second implementation.
    """
    out: dict[str, Any] = {}
    for tier, row in offsets.items():
        shaped = {
            key: {
                "explains": row[key]["explains_of_the_offset"],
                "assessed": bool(row[key]["chromosomes_assessed"]),
                "degenerate": bool(row[key]["chromosomes_degenerate"]),
            }
            for key in ARMS
        }
        out[tier] = where_it_sits(shaped)
        out[tier]["chromosomes_assessed"] = {k: row[k]["chromosomes_assessed"] for k in ARMS}
    return out


def comparisons_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """The window comparisons pooled as medians, with the imbalance, the presence and the dropped count."""
    out: dict[str, Any] = {}
    for key in COMPARISONS:
        rows = [r["comparisons"][key] for r in runs if r["comparisons"].get(key, {}).get("assessed")]
        if not rows:
            out[key] = {"chromosomes_assessed": 0, "reads": "not assessed on any chromosome"}
            continue
        out[key] = {
            "chromosomes_assessed": len(rows),
            "targets": sum(x["targets"] for x in rows),
            "controls": sum(x["controls"] for x in rows),
            "coverage_measure": COVERAGE_MEASURE,
            "input_presence": {
                name: {
                    "kind": rows[0]["input_presence"][name]["kind"],
                    "chromosomes_where_it_is_free": sum(
                        1 for x in rows if x["input_presence"][name]["universal"]
                    ),
                }
                for name in rows[0]["input_presence"]
            },
            "targets_below_full_coverage": sum(x["targets_below_full_coverage"] for x in rows),
            "controls_below_full_coverage": sum(x["controls_below_full_coverage"] for x in rows),
            "matched_difference": {
                k: round(median(x[k]["matched"]["difference"] for x in rows), 4) for k in STRATIFICATIONS
            },
            "targets_matched": {k: sum(x[k]["targets_matched"] for x in rows) for k in STRATIFICATIONS},
            "dropped_for_want_of_a_control": {
                k: sum(x[k]["dropped_for_want_of_a_control"] for x in rows) for k in STRATIFICATIONS
            },
            "imbalance_medians": {
                axis: {
                    side: _median([x[LAST]["imbalance"][axis][f"{side}_median"] for x in rows])
                    for side in ("target", "control")
                }
                for axis in ("length", "gc", "tss", "coverage")
            },
            **sign_test([x[LAST]["matched"]["difference"] for x in rows]),
        }
    return out


def measure(chroms: list[str]) -> list[dict[str, Any]]:
    """Read every chromosome from its store, which is where the hours go."""
    runs = []
    for c in chroms:
        print(f"{c} ...", flush=True)
        r = analyse(c)
        if not r.get("assessed"):
            print(f"  skipped: {r.get('why')}", flush=True)
            continue
        arm = r["offset"][UNION_ARM]
        print(
            f"  union takes {arm.get('excluded_bases'):,} bases, leaves a background "
            f"{arm.get('block_share_of_the_remaining_background')} blocks, "
            f"degenerate {arm.get('degenerate')}",
            flush=True,
        )
        runs.append(r)
    return runs


def from_result(name: str, results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """The per-chromosome runs a committed result already carries, so the distillation can be redone.

    Reading the stores is the expensive half and distilling them is not, so a wording or a reading
    that has to be corrected should not cost an hour of rereading sequence that has not changed.
    Nothing is measured here. The one thing recomputed is each chromosome's `where_the_union_sits`,
    from that chromosome's own stored `explained`: it is a reading of numbers rather than a number,
    and leaving the stale one beside a freshly distilled top level is how a file comes to hold two
    answers to one question.
    """
    prev = load_result(name, results_dir) or {}
    rows = prev.get("per_chromosome") or {}
    if not rows:
        raise FileNotFoundError(f"no per-chromosome runs in {name}; run without --from-result first")
    print(f"re-distilling {len(rows)} chromosomes from {name}, measuring nothing", flush=True)
    runs = [{"chrom": c, "assessed": True, **r} for c, r in rows.items()]
    for r in runs:
        r["where_the_union_sits"] = {t: where_it_sits(r["explained"][t]) for t in TIERS}
    return runs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chrom", action="append", help="one chromosome; repeatable (default: all with a store)")
    ap.add_argument("--name", default="panel_union_arms")
    ap.add_argument(
        "--from-result",
        action="store_true",
        help="re-distil the top level from the committed per-chromosome runs, measuring nothing",
    )
    args = ap.parse_args(argv)
    runs = from_result(args.name) if args.from_result else measure(args.chrom or chromosomes())
    kept = [r for r in runs if r["chrom"] not in EXCLUDED]
    if not kept:
        print(json.dumps({"result": None, "chromosomes": 0}))
        return 1
    offsets = offset_across(kept)
    conditioned = all(r["claim_available"] == CLAIM_WITH_COVERAGE for r in kept)
    result = {
        "chromosomes": [r["chrom"] for r in runs],
        "excluded_from_the_medians": list(EXCLUDED),
        "evidence": SOURCES,
        "coverage_measure": COVERAGE_MEASURE,
        "covariates": list(COVARIATES),
        # from the module, not from a run: a re-distillation has no run to read it off, and a
        # definition that differed between the two would be the drift this lane exists to avoid
        "covariate_definitions": FROM_LANE,
        "near_threshold": NEAR,
        # one of two sentences, never a third, and carried as data so that no reader can upgrade it
        "claim_available": CLAIM_WITH_COVERAGE if conditioned else CLAIM_WITHOUT_COVERAGE,
        "reproduces_the_panel": {
            "chromosomes_checked": sum(1 for r in runs if r["reproduces_the_panel"]["tiers_checked"]),
            "tiers_checked": sum(r["reproduces_the_panel"]["tiers_checked"] for r in runs),
            "largest_absolute_difference": max(
                (
                    r["reproduces_the_panel"]["largest_absolute_difference"]
                    for r in runs
                    if r["reproduces_the_panel"]["largest_absolute_difference"] is not None
                ),
                default=None,
            ),
        },
        "eligible": {
            "chromosome_bases": sum(r["chromosome_bases"] for r in kept),
            "kilobases_in_the_panel_alignment": sum(
                r["eligible"]["kilobases_in_the_panel_alignment"] for r in kept
            ),
            "windows_in_the_panel_alignment": sum(
                r["eligible"]["windows_in_the_panel_alignment"] for r in kept
            ),
            "covariates_declared": list(COVARIATES),
            "covariates_with_bases_somewhere": sorted(
                {n for r in kept for n in r["eligible"]["covariates_with_bases_on_this_chromosome"]}
            ),
            "tiers_declared": list(TIERS),
        },
        "overlap": overlap_across(kept),
        "groups": groups_across(kept),
        "offset": offsets,
        "where_the_union_sits": union_across(offsets),
        "comparisons": comparisons_across(kept),
        "per_chromosome": {r["chrom"]: {k: r[k] for k in PER_CHROMOSOME} for r in runs},
    }
    save_result(args.name, result)
    print("\nthe covariates, pairwise:")
    for key, row in result["overlap"]["pairs"].items():
        print(f"  {key:44} {row.get('shared_bases'):>13,} shared  {row.get('share_of_the_smaller')}")
    print(
        f"\nunion {result['overlap']['union_bases']:,} bases, "
        f"sum of the singles {result['overlap']['sum_of_the_singles_bases']:,}"
    )
    for tier in TIERS:
        row = offsets[tier]
        base = row["as_the_panel_builds_it"]
        print(f"\n{tier}: median {base['median']} above 1 on {base['above']} of {base['chromosomes']}")
        for key in ARMS:
            r = row[key]
            print(
                f"  {key:34} median {r['median']} ({r['chromosomes_assessed']}) "
                f"blocks {r['block_share_of_the_remaining_background_median']} {r.get('reads')}"
            )
        print(f"  -> {result['where_the_union_sits'][tier]['reads']}")
    print("\n" + result["claim_available"])
    print(json.dumps({"result": args.name, "chromosomes": len(runs)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
