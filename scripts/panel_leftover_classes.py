# SPDX-License-Identifier: AGPL-3.0-or-later
"""What is the 2,006 Mb the budget never cut into a block, and does the offset live there?

    uv run python scripts/panel_leftover_classes.py                 # every chromosome with a store
    uv run python scripts/panel_leftover_classes.py --chrom chr21   # one, for a quick read

`panel_background_composition.py` measured what the panel's matched background holds and found that
composition explains at most half of the offset every non-coding tier reads above it. It then named
its own leftover: 71% of the background lies outside every unknown block and nobody had measured what
that sequence is. This does, and it asks the question the leftover poses -- is the offset a property
of a biological class, or of how the budget cuts blocks?

`genome/unknown.unknown_blocks` is replayed with the bookkeeping it discards, its blocks are checked
against the committed budget one for one, and the never-cut sequence is split into the seven classes
that partition it. The background is then rebuilt with each class taken out through the panel's own
`build_background`, and the weight is carried by `compare.standardised` over 100-base windows, where
both arms are groups of sequence and neither is a tier average.

chrY is measured and reported but kept out of the medians, as everywhere else in the panel.

Writes data/results/panel_leftover_classes.json.
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
from genomeos.attribution.panel_leftover import (
    CARRIES_MOST,
    CLAIM_WITH_COVERAGE,
    CLAIM_WITHOUT_COVERAGE,
    CLASSES,
    COVERAGE_MEASURE,
    INTERGENIC,
    LENGTH_RULE_ONLY,
    SEPARABLE_MIN,
    SOURCES,
    UNACCOUNTED,
    analyse,
    verdict,
)
from genomeos.results import RESULTS_DIR, load_result, save_result

EXCLUDED = ("chrY",)
NAMES = (*CLASSES, UNACCOUNTED)
ARMS = (
    *(f"without_{c}" for c in NAMES),
    "without_every_gene_body",
    "without_the_intergenic_never_cut",
    "without_every_never_cut_base",
)
STRATIFICATIONS = ("gc_and_timing", "gc_timing_and_tss", "gc_timing_tss_and_coverage")
LAST = STRATIFICATIONS[-1]
BLOCK_ARM = "every_block_against_the_never_cut"
# what a per-chromosome run keeps in the result: everything the aggregators below read, so that the
# distillation can be redone from the committed file without touching a store again
PER_CHROMOSOME = (
    "chromosome_bases",
    "claim_available",
    "reproduces_the_budget",
    "eligible",
    "cut",
    "background",
    "leftover",
    "split",
    "by_dominant_class",
    "offset",
    "explained",
    "comparisons",
    "coverage",
)


def chromosomes(results_dir: Path = RESULTS_DIR, cache: Path = CACHE) -> list[str]:
    """Every chromosome with a committed panel result, a budget and a local store to rebuild from."""
    out = []
    for p in sorted(glob.glob(str(results_dir / "human_panel_chr*.json"))):
        c = Path(p).stem.split("_")[-1]
        if (cache / c / "meta.json").exists() and (results_dir / f"budget_{c}.json").exists():
            out.append(c)
    # smallest store first: a long run that has to be stopped has read the cheap chromosomes already
    return sorted(out, key=lambda c: (cache / c / "sites.tsv.gz").stat().st_size)


def _median(values: list[Any]) -> float | None:
    """The median of the values that exist, or None: a covariate absent on every chromosome is not 0."""
    kept = [v for v in values if v is not None]
    return round(median(kept), 4) if kept else None


def sign_test(values: list[float], above: float = 0.0) -> dict[str, Any]:
    n = len(values)
    k = sum(v > above for v in values)
    return {
        "chromosomes": n,
        "above": k,
        "one_sided_p": round(sum(comb(n, i) for i in range(k, n + 1)) / 2**n, 6) if n else None,
    }


def split_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """The never-cut bases summed over the chromosomes, by class, with the reason beside each."""
    by: dict[str, Any] = {}
    total = sum(r["split"]["by_class"][n]["bases"] for r in runs for n in NAMES)
    for n in NAMES:
        b = sum(r["split"]["by_class"][n]["bases"] for r in runs)
        by[n] = {
            "bases": b,
            "share_of_the_never_cut": round(b / total, 6) if total else None,
            "chromosomes_with_bases": sum(1 for r in runs if r["split"]["by_class"][n]["bases"]),
            "why_the_budget_did_not_cut_it": runs[0]["split"]["by_class"][n]["why_the_budget_did_not_cut_it"],
        }
    length_rule = sum(by[n]["bases"] for n in INTERGENIC)
    share = length_rule / total if total else None
    return {
        "chromosomes": len(runs),
        "never_cut_bases": total,
        "kilobases_outside_every_block": sum(r["split"]["kilobases"] for r in runs),
        "bases_kept_by_the_background": sum(r["split"]["bases_kept"] for r in runs),
        "bases_in_no_class": sum(r["split"]["bases_in_no_class"] for r in runs),
        "by_class": by,
        "separability": {
            "bases_declined_for_being_inside_a_gene_span": total - length_rule,
            "bases_declined_for_the_length_rule_alone": length_rule,
            "share_that_can_separate_the_two_readings": round(share, 6) if share is not None else None,
            "threshold": SEPARABLE_MIN,
            "separable": bool(share is not None and share >= SEPARABLE_MIN),
        },
    }


def groups_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Median covariates per group: the background, the leftover, and the leftover by dominant class.

    Length, GC, distance to a coding TSS and coverage are printed for every group that is compared,
    and a group assessed on fewer chromosomes says so rather than borrowing the others'.
    """
    out: dict[str, Any] = {}
    wanted = ["background", "leftover"]
    labels = sorted({k for r in runs for k in r["by_dominant_class"]})
    for name in wanted + labels:
        rows = [
            (r[name] if name in wanted else r["by_dominant_class"].get(name, {}))
            for r in runs
            if (r[name] if name in wanted else r["by_dominant_class"].get(name, {})).get("bases")
        ]
        rows = [x for x in rows if x.get("gc_median") is not None]
        if not rows:
            out[name] = {"chromosomes_assessed": 0, "why": "no chromosome carries this group"}
            continue
        out[name] = {
            "chromosomes_assessed": len(rows),
            "kilobases": sum(x["kilobases"] for x in rows),
            "bases": sum(x["bases"] for x in rows),
            "unit_bases": 1000,
            "gc_median": round(median(x["gc_median"] for x in rows), 4),
            "tss_median": int(median(x["tss_median"] for x in rows)),
            "tss_assessed_chromosomes": sum(1 for x in rows if x.get("tss_assessed")),
            "coverage_measure": COVERAGE_MEASURE,
            "coverage_median": round(
                median(x["coverage_median"] for x in rows if x.get("coverage_median") is not None), 4
            ),
            "coverage_assessed_chromosomes": sum(1 for x in rows if x.get("coverage_median") is not None),
            "recurring_per_kb_median": round(median(x["recurring_per_kb"] for x in rows), 3),
        }
    return out


def offset_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """Every tier's ratio under every exclusion, on the chromosomes where that exclusion was made.

    The explained share is computed on the SAME chromosomes as the ratio it is compared with, so a
    class assessed on fewer chromosomes cannot borrow the others' offset. A degenerate arm keeps its
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
            pairs = [
                (r["offset"]["as_the_panel_builds_it"].get(tier), r["offset"][key].get(tier))
                for r in runs
                if r["offset"].get(key, {}).get("assessed")
            ]
            pairs = [(a, b) for a, b in pairs if a is not None and b is not None]
            degenerate = sum(1 for r in runs if r["offset"].get(key, {}).get("degenerate"))
            if not pairs:
                row[key] = {
                    "chromosomes_assessed": 0,
                    "chromosomes_degenerate": degenerate,
                    "median": None,
                    "explains_of_the_offset": None,
                    "reads": "not assessed",
                    "why": "the arm is degenerate or the class has no bases"
                    if degenerate
                    else "the class has no bases on any chromosome",
                }
                continue
            base, here = median(a for a, _ in pairs), median(b for _, b in pairs)
            share = round((base - here) / (base - 1), 4) if base != 1 else None
            row[key] = {
                "chromosomes_assessed": len(pairs),
                "chromosomes_degenerate": degenerate,
                "median": round(here, 4),
                "median_as_built_on_the_same_chromosomes": round(base, 4),
                "explains_of_the_offset": share,
                "reads": reads_explained(share),
                **sign_test([b for _, b in pairs], 1.0),
            }
        out[tier] = row
    return out


def comparisons_across(runs: list[dict[str, Any]]) -> dict[str, Any]:
    """The window comparisons pooled as medians, with the imbalance and the dropped count kept."""
    out: dict[str, Any] = {}
    keys = [BLOCK_ARM] + [f"{n}_against_the_rest_of_the_never_cut" for n in NAMES]
    for key in keys:
        rows = [r["comparisons"][key] for r in runs if r["comparisons"].get(key, {}).get("assessed")]
        if not rows:
            out[key] = {"chromosomes_assessed": 0, "reads": "not assessed on any chromosome"}
            continue
        diffs = [x[LAST]["matched"]["difference"] for x in rows]
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
            **sign_test(diffs),
        }
    return out


def carried_by(offsets: dict[str, Any]) -> dict[str, Any]:
    """Which classes move at least `CARRIES_MOST` of a tier's offset, and which widen it."""
    out: dict[str, Any] = {}
    for tier, row in offsets.items():
        carries = []
        widens = []
        for key in ARMS:
            share = row[key].get("explains_of_the_offset")
            if share is None:
                continue
            if share >= CARRIES_MOST:
                carries.append((key, share))
            elif share < 0:
                widens.append((key, share))
        out[tier] = {
            "threshold": CARRIES_MOST,
            "classes_carrying_most_of_the_offset": [k for k, _ in sorted(carries, key=lambda kv: -kv[1])],
            "classes_widening_it": [k for k, _ in sorted(widens, key=lambda kv: kv[1])],
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
        s = r["split"]
        print(
            f"  {s['kilobases']:>7} kb never cut, "
            f"{s['separability']['share_that_can_separate_the_two_readings']} separable, "
            f"blocks vs never cut "
            f"{r['comparisons'][BLOCK_ARM].get(LAST, {}).get('matched', {}).get('difference')}",
            flush=True,
        )
        runs.append(r)
    return runs


def from_result(name: str, results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """The per-chromosome runs a committed result already carries, so the distillation can be redone.

    Reading the stores is the expensive half and distilling them is not, and a wording or an
    aggregation that has to be corrected should not cost an hour of rereading sequence that has not
    changed. Nothing is measured here: the rows are exactly what the committed file holds.
    """
    prev = load_result(name, results_dir) or {}
    rows = prev.get("per_chromosome") or {}
    if not rows:
        raise FileNotFoundError(f"no per-chromosome runs in {name}; run without --from-result first")
    print(f"re-distilling {len(rows)} chromosomes from {name}, measuring nothing", flush=True)
    return [{"chrom": c, "assessed": True, **r} for c, r in rows.items()]


def _gap(comps: dict[str, Any], key: str) -> float | None:
    """The matched difference with coverage held, or None where the arm was never assessed."""
    row = comps.get(key) or {}
    return row["matched_difference"][LAST] if row.get("chromosomes_assessed") else None


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--chrom", action="append", help="one chromosome; repeatable (default: all with a store)")
    ap.add_argument("--name", default="panel_leftover_classes")
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
    split = split_across(kept)
    offsets = offset_across(kept)
    comps = comparisons_across(kept)
    rule_key = f"{LENGTH_RULE_ONLY}_against_the_rest_of_the_never_cut"
    block_gap = _gap(comps, BLOCK_ARM)
    rule_gap = _gap(comps, rule_key)
    conditioned = all(r["claim_available"] == CLAIM_WITH_COVERAGE for r in kept)
    result = {
        "chromosomes": [r["chrom"] for r in runs],
        "excluded_from_the_medians": list(EXCLUDED),
        "evidence": SOURCES,
        "coverage_measure": COVERAGE_MEASURE,
        # one of two sentences, never a third, and carried as data so that no reader can upgrade it
        "claim_available": CLAIM_WITH_COVERAGE if conditioned else CLAIM_WITHOUT_COVERAGE,
        **verdict(
            split["separability"]["separable"],
            block_gap,
            rule_gap,
            comps[rule_key].get("targets") or 0,
        ),
        "reproduces_the_budget": {
            "chromosomes_checked": len(runs),
            "chromosomes_identical": sum(1 for r in runs if r["reproduces_the_budget"]["identical"]),
            "blocks_replayed": sum(r["reproduces_the_budget"]["blocks_replayed"] for r in runs),
            "blocks_only_in_the_replay": sum(r["reproduces_the_budget"]["only_in_the_replay"] for r in runs),
            "blocks_only_in_the_committed_budget": sum(
                r["reproduces_the_budget"]["only_in_the_committed_budget"] for r in runs
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
            "classes_declared": list(NAMES),
            "classes_with_bases_somewhere": sorted(
                {n for r in kept for n in r["eligible"]["classes_with_bases_on_this_chromosome"]}
            ),
            "tiers_declared": list(TIERS),
        },
        "split": split,
        "groups": groups_across(kept),
        "offset": offsets,
        "carried_by": carried_by(offsets),
        "comparisons": comps,
        "per_chromosome": {r["chrom"]: {k: r[k] for k in PER_CHROMOSOME} for r in runs},
    }
    save_result(args.name, result)
    print("\nthe never-cut sequence, by class:")
    for n in NAMES:
        row = split["by_class"][n]
        print(f"  {n:34} {row['bases']:>12,} bases  {row['share_of_the_never_cut']}")
    print(f"\nseparable: {split['separability']['share_that_can_separate_the_two_readings']}")
    print(f"blocks against the never cut: {block_gap}   {LENGTH_RULE_ONLY}: {rule_gap}")
    for tier in TIERS:
        row = offsets[tier]
        base = row["as_the_panel_builds_it"]
        print(f"\n{tier}: median {base['median']} above 1 on {base['above']} of {base['chromosomes']}")
        for key in ARMS:
            r = row[key]
            print(f"  {key:44} median {r['median']} ({r['chromosomes_assessed']}) {r.get('reads')}")
    print("\n" + result["claim_available"])
    print(result["verdict"])
    print(result["mechanism"])
    print(json.dumps({"result": args.name, "chromosomes": len(runs)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
