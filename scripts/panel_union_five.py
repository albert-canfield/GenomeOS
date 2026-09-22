# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fifth covariate, the one that runs the other way, put into the panel's union of exclusions.

    uv run python scripts/panel_union_five.py                  # every chromosome with a store
    uv run python scripts/panel_union_five.py --chrom chr21    # one, for a quick read
    uv run python scripts/panel_union_five.py --from-result    # re-distil, measuring nothing

`panel_union_arms.py` took four covariates out of the panel's matched background at once -- fossil
98.7%, regulatory 141.0%, neutral 55.5% of the offset -- and named the omission in its own "left
undone": segmental duplication, measured by the first lane of all, *widens* the offset when it is
removed (29% / 34% / 33%) and was not in the union. This adds it, and rebuilds the four-way arm in
the same run so the two unions are a difference between arms rather than between runs.

Every share is carried with its sign. The narrowing singles and the widening ones are summed
separately as well as together, and there is no mean anywhere that mixes them: a union that contains
a widening term can explain LESS than the same union without it, and if that is what the arms do it
is the result and not an anomaly to smooth.

The segmental-duplication track is per chromosome under `data/results`; `--tracks-dir` points the run
at a checkout that holds the ones this one does not.

chrY is measured and reported but kept out of the medians, as everywhere else in the panel.

Writes data/results/panel_union_five_arms.json.
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
from genomeos.attribution.panel_background import TRACKS, reads_explained
from genomeos.attribution.panel_union_five import (
    AGREE,
    ARMS,
    CLAIM_WITH_COVERAGE,
    CLAIM_WITHOUT_COVERAGE,
    COVARIATES,
    COVERAGE_MEASURE,
    FIFTH,
    FROM_LANE,
    NEAR,
    SIGN_RULE,
    SOURCES,
    STRATIFICATIONS,
    UNION,
    UNION_FOUR,
    WIDENING_WHEN_MEASURED,
    WITH_A_WIDENING_TERM,
    analyse,
    where_it_sits,
    where_the_five_sit,
    which_way_each_covariate_runs,
)
from genomeos.results import RESULTS_DIR, load_result, save_result

EXCLUDED = ("chrY",)
LAST = STRATIFICATIONS[-1]
FIVE_ARM = f"without_{UNION}"
FOUR_ARM = f"without_{UNION_FOUR}"
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
    "which_way_each_covariate_runs",
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
    four_union = sum(r["overlap"]["the_four_way"]["union_bases"] for r in runs)
    four_sum = sum(r["overlap"]["the_four_way"]["sum_of_the_singles_bases"] for r in runs)
    return {
        "chromosomes": len(runs),
        "pairs": out,
        "union_bases": union,
        "sum_of_the_singles_bases": summed,
        "bases_counted_more_than_once_by_the_sum": summed - union,
        "the_sum_over_the_union": round(summed / union, 4) if union else None,
        "the_four_way": {
            "union_bases": four_union,
            "sum_of_the_singles_bases": four_sum,
            "bases_counted_more_than_once_by_the_sum": four_sum - four_union,
            "the_sum_over_the_union": round(four_sum / four_union, 4) if four_union else None,
        },
        "bases_the_fifth_adds_to_the_union": union - four_union,
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


def which_way_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Which way each covariate runs, pooled: the two rates, and how many chromosomes read each way."""
    out: dict[str, Any] = {}
    for name in COVARIATES:
        rows = [r["which_way_each_covariate_runs"][name] for r in runs]
        readings: dict[str, int] = {}
        for x in rows:
            readings[x["reads"]] = readings.get(x["reads"], 0) + 1
        out[name] = {
            "chromosomes": len(rows),
            "declared_as_widening": name in WIDENING_WHEN_MEASURED,
            "recurring_per_kb_where_it_is": _median([x["recurring_per_kb_where_it_is"] for x in rows]),
            "recurring_per_kb_where_it_is_not": _median(
                [x["recurring_per_kb_where_it_is_not"] for x in rows]
            ),
            "share_of_background_bases_median": _median([x["share_of_background_bases"] for x in rows]),
            "chromosomes_reading": readings,
            "chromosomes_where_the_rate_and_the_arm_agree": sum(
                1 for x in rows if x["the_rate_and_the_arm"] == AGREE
            ),
        }
    return out


def offset_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Every tier's ratio under every arm, on the chromosomes where that arm was made and assessed.

    The explained share is computed on the SAME chromosomes as the ratio it is compared with, so an
    arm assessed on fewer chromosomes cannot borrow the others' offset. A degenerate arm keeps its
    numbers and its flag and is never given an explained share. The share keeps its sign: an arm that
    widens the offset reads negative here and is never folded into a mean with one that narrows it.
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
                "bases_taken_out_of_the_background_median": _median(
                    [a.get("excluded_bases") for a in arm if a.get("assessed")]
                ),
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
                    "widens_the_offset": bool(share is not None and share < 0),
                    "reads": reads_explained(share),
                    **sign_test([b for _, b in pairs], 1.0),
                }
            )
            row[key] = entry
        out[tier] = row
    return out


def union_across(offsets: dict[str, Any]) -> dict[str, Any]:
    """Both readings from the pooled medians: the four-way as its own lane reads it, and the five-way.

    The shapes the two functions read are built once and given to both, so the four-way row here is
    `panel_union.where_it_sits` on the four-way arm -- that lane's rule, that lane's words -- and the
    five-way row is the sign-aware reading beside it, not a second opinion about the same arm.
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
        out[tier] = {
            "the_four_way_as_it_stands": where_it_sits(shaped),
            "the_five_way": where_the_five_sit(shaped),
            "chromosomes_assessed": {k: row[k]["chromosomes_assessed"] for k in ARMS},
        }
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


def measure(chroms: list[str], tracks_dir: Path = TRACKS) -> list[dict[str, Any]]:
    """Read every chromosome from its store, which is where the hours go."""
    runs = []
    for c in chroms:
        print(f"{c} ...", flush=True)
        r = analyse(c, tracks_dir=tracks_dir)
        if not r.get("assessed"):
            print(f"  skipped: {r.get('why')}", flush=True)
            continue
        five, four = r["offset"][FIVE_ARM], r["offset"][FOUR_ARM]
        print(
            f"  four takes {four.get('excluded_bases'):,} bases, five takes "
            f"{five.get('excluded_bases'):,}, leaving a background "
            f"{five.get('block_share_of_the_remaining_background')} blocks, "
            f"degenerate {five.get('degenerate')}",
            flush=True,
        )
        for tier in ("fossil", "regulatory", "neutral"):
            sit = r["where_the_union_sits"].get(tier, {}).get("the_five_way", {})
            print(
                f"    {tier:12} four {sit.get('four_way_union_explains')} "
                f"five {sit.get('union_explains')} "
                f"delta {sit.get('the_fifth_covariate_moves_the_union_by')}",
                flush=True,
            )
        runs.append(r)
    return runs


def from_result(name: str, results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """The per-chromosome runs a committed result already carries, so the distillation can be redone.

    Reading the stores is the expensive half and distilling them is not, so a wording or a reading
    that has to be corrected should not cost hours of rereading sequence that has not changed.
    Nothing is measured here. What is recomputed is each chromosome's `where_the_union_sits` and its
    `which_way_each_covariate_runs`, from that chromosome's own stored `explained` and
    `feature_splits`: both are readings of numbers rather than numbers, and leaving a stale one beside
    a freshly distilled top level is how a file comes to hold two answers to one question.
    """
    prev = load_result(name, results_dir) or {}
    rows = prev.get("per_chromosome") or {}
    if not rows:
        raise FileNotFoundError(f"no per-chromosome runs in {name}; run without --from-result first")
    print(f"re-distilling {len(rows)} chromosomes from {name}, measuring nothing", flush=True)
    runs = [{"chrom": c, "assessed": True, **r} for c, r in rows.items()]
    for r in runs:
        r["where_the_union_sits"] = {
            t: {
                "the_four_way_as_it_stands": where_it_sits(r["explained"][t]),
                "the_five_way": where_the_five_sit(r["explained"][t]),
            }
            for t in TIERS
        }
        r["which_way_each_covariate_runs"] = which_way_each_covariate_runs(
            r["feature_splits"], r["explained"]
        )
    return runs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chrom", action="append", help="one chromosome; repeatable (default: all with a store)")
    ap.add_argument("--name", default="panel_union_five_arms")
    ap.add_argument(
        "--tracks-dir",
        default=str(TRACKS),
        help="where the per-chromosome segmental-duplication and cCRE tracks are (default: data/results)",
    )
    ap.add_argument(
        "--from-result",
        action="store_true",
        help="re-distil the top level from the committed per-chromosome runs, measuring nothing",
    )
    args = ap.parse_args(argv)
    runs = (
        from_result(args.name)
        if args.from_result
        else measure(args.chrom or chromosomes(), Path(args.tracks_dir))
    )
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
        "the_covariate_this_lane_adds": FIFTH,
        "covariates_declared_widening": list(WIDENING_WHEN_MEASURED),
        # from the module, not from a run: a re-distillation has no run to read it off, and a
        # definition that differed between the two would be the drift this lane exists to avoid
        "covariate_definitions": FROM_LANE,
        "sign_rule": SIGN_RULE,
        "what_a_union_with_a_widening_term_means": WITH_A_WIDENING_TERM,
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
            "covariates_without_bases_anywhere": sorted(
                {
                    n
                    for n in COVARIATES
                    if all(n not in r["eligible"]["covariates_with_bases_on_this_chromosome"] for r in kept)
                }
            ),
            "tiers_declared": list(TIERS),
        },
        "overlap": overlap_across(kept),
        "groups": groups_across(kept),
        "which_way_each_covariate_runs": which_way_across(kept),
        "offset": offsets,
        "where_the_union_sits": union_across(offsets),
        "comparisons": comparisons_across(kept),
        "per_chromosome": {r["chrom"]: {k: r[k] for k in PER_CHROMOSOME} for r in runs},
    }
    save_result(args.name, result)
    print("\nwhich way each covariate runs, from the background's own kilobases:")
    for name, row in result["which_way_each_covariate_runs"].items():
        print(
            f"  {name:24} {row['recurring_per_kb_where_it_is']} per kb where it is, "
            f"{row['recurring_per_kb_where_it_is_not']} where it is not"
        )
    print("\nthe covariates, pairwise:")
    for key, row in result["overlap"]["pairs"].items():
        print(f"  {key:48} {row.get('shared_bases'):>13,} shared  {row.get('share_of_the_smaller')}")
    print(
        f"\nfive-way union {result['overlap']['union_bases']:,} bases, "
        f"four-way {result['overlap']['the_four_way']['union_bases']:,}, "
        f"sum of the five {result['overlap']['sum_of_the_singles_bases']:,}"
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
        sit = result["where_the_union_sits"][tier]["the_five_way"]
        print(f"  -> {sit['reads']}")
        print(f"  -> {sit['the_fifth_covariate_reads']}")
    print("\n" + result["claim_available"])
    print(json.dumps({"result": args.name, "chromosomes": len(runs)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
