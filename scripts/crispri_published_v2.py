# SPDX-License-Identifier: AGPL-3.0-or-later
"""The CRISPRi published-benchmark result rebuilt under a new name: 2,000 draws, and no deletion gain
where the deletion feature was never available.

    uv run python scripts/crispri_published_v2.py

Writes `data/results/crispri_published_v2.json`. The 2026-09-27 result,
`data/results/crispri_published.json`, is left exactly as it is: it stays the record of what was
computed and quoted then, and this file is the corrected run beside it, not a replacement of it.

Two defects of that result are what this run exists to remove, both already fixed in
`genomeos/attribution/crispri.py` at `bc28f79`:

1. every interval there rests on 200 chromosome resamples and records nothing about how it was made.
   `crispri.BOOTSTRAPS` is now 2,000 and `crispri._interval_provenance` states, beside each interval,
   the number of resampling clusters, the draws requested, the draws dropped and whether the kept
   draws reach `crispri.MIN_RESAMPLES`;
2. three held-out strata there report a deletion gain and an interval for cell types that carry no
   deletion value at all (HCT116 +0.0004, Jurkat +0.0038, WTC11 -0.0049). `crispri.score_published`
   now routes every per-cell gain through `crispri.gain_where_available`, so such a stratum reports
   `gain: null` with `crispri.UNAVAILABLE_GAIN` as the reason and the gain is not computed at all.

Nothing else about the measurement changes. The model, the features, the frozen weights, the
estimator, the seed and the pair sets are those of the 2026-09-27 run, and no AlphaGenome request is
made: every deletion value comes from the sweep's existing cache, as `alphagenome_requests: 0` says.

The result carries two blocks this script adds for the reader:

* `beside_the_2026_09_27_result`: each interval with its population named, the 2026-09-27 figure at
  200 draws quoted (transcribed in `BEFORE`, never recomputed) and this run's figure beside it. The
  populations are not interchangeable: "all 4,378 held-out pairs" and "the 1,918 K562 held-out
  pairs" are different sets of pairs and are never compared with each other here.
* `what_did_not_change`: the point AUPRCs, the bands against the published figures, the
  coverage-matched arm and the request count, each checked against the 2026-09-27 value.

`result_manifest.import_closure` answers the shared checkout: several lanes hold uncommitted files
in this working tree, so the closure of this script is computed by `import_closure` and set against
the uncommitted code the stamp names, rather than asserted by hand.
"""

from __future__ import annotations

import ast
import json
import sys
import time
from collections import deque
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from scripts.crispri_published import manifest as published_manifest  # noqa: E402

NAME = "crispri_published_v2"
ENTRY = "scripts/crispri_published_v2.py"
OLD = "data/results/crispri_published.json"
OLD_DATE = "2026-09-27"
OLD_DRAWS = 200

#: What this run changes about the 2026-09-27 run, as the manifest's exclusions state it.
WHY_V2 = (
    "rebuilt under a new name on 2026-10-02: 2,000 chromosome resamples per interval in place of "
    "200, each interval stating its clusters and draws, and no deletion gain where the deletion "
    f"feature is unavailable; the {OLD_DATE} result at {OLD} is kept unchanged beside this one"
)

#: The three held-out cell types the sweep never scored a deletion value in. They are not dropped:
#: their model AUPRCs are still reported, because a model AUPRC is a real output of a fitted model.
#: Only the deletion gain is refused, because the two models differ there by their fitted weights
#: alone, which is a property of the model form and not a measurement of the deletion feature.
WITHOUT_A_DELETION_VALUE = ("HCT116", "Jurkat", "WTC11")


# --- the import closure of this script, computed rather than listed ------------------------------


def _module_files(module: str, root: Path) -> list[str]:
    """The repository files a dotted module name reads when it is imported: the module itself and the
    `__init__.py` of every package above it. Empty for a module that is not in this repository."""
    parts = module.split(".")
    found: list[str] = [
        pkg.relative_to(root).as_posix()
        for i in range(1, len(parts))
        if (pkg := root.joinpath(*parts[:i], "__init__.py")).is_file()
    ]
    for candidate in (root.joinpath(*parts[:-1], f"{parts[-1]}.py"), root.joinpath(*parts, "__init__.py")):
        if candidate.is_file():
            found.append(candidate.relative_to(root).as_posix())
            break
    return found


def _imported_modules(path: Path, root: Path) -> list[str]:
    """Every dotted module name one file imports, with relative imports resolved against its package.

    A `from X import a, b` contributes `X` and also `X.a` and `X.b`, because the name after `import`
    may itself be a submodule; `_module_files` keeps only the ones that are files in this repository.
    """
    package = list(path.relative_to(root).with_suffix("").parts[:-1])
    out: list[str] = []
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            out.extend(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package[: len(package) - (node.level - 1)]
                prefix = ".".join([*base, *([node.module] if node.module else [])])
            else:
                prefix = node.module or ""
            if prefix:
                out.append(prefix)
                out.extend(f"{prefix}.{a.name}" for a in node.names)
    return out


def import_closure(entry: str = ENTRY, root: Path = ROOT) -> list[str]:
    """Every file in this repository that `entry` imports, directly or at any remove.

    Found by parsing the import statements, not by importing: the answer is then the same in any
    checkout of the same revision, which is what lets a rebuild in a clean worktree reproduce it.
    Modules outside the repository (the standard library, the dependencies) are not on the list.
    """
    seen = {entry}
    queue = deque([entry])
    while queue:
        rel = queue.popleft()
        path = root / rel
        if not path.is_file():
            continue
        for module in _imported_modules(path, root):
            for f in _module_files(module, root):
                if f not in seen:
                    seen.add(f)
                    queue.append(f)
    return sorted(seen)


def closure_block(root: Path = ROOT) -> dict[str, Any]:
    """The closure, and the uncommitted code of this shared checkout set against it.

    Only facts a clean rebuild reaches too are recorded here. Which files were uncommitted when this
    ran is the run's property, not the result's, and the stamp records it under
    `result_manifest.code.dirty_code_paths`; what is recorded here is the intersection, which is
    empty in this checkout exactly when no uncommitted file is on the closure, and empty in a clean
    worktree because there is nothing to intersect. A foreign file that did reach the closure would
    therefore show up as a difference in `scripts/manifest_rebuild.py`, not pass unnoticed.
    """
    closure = import_closure(root=root)
    dirty = mf.code_revision(root).get("dirty_code_paths") or []
    return {
        "entry": ENTRY,
        "how": "scripts/crispri_published_v2.py import_closure: the import statements parsed and "
        "followed transitively, repository files only, sorted; not a list written by hand",
        "files": closure,
        "count": len(closure),
        "uncommitted_code_at_write": "named in result_manifest.code.dirty_code_paths, which a clean "
        "rebuild leaves empty; the intersection below is what both runs can state",
        "uncommitted_code_on_the_closure": sorted(set(dirty) & set(closure)),
    }


# --- the 2026-09-27 figures, transcribed ---------------------------------------------------------

#: The intervals of the 2026-09-27 result, each with the population it was measured on.
#: `path` reads the same field of this run. The figures are transcribed from
#: data/results/crispri_published.json as committed (written at 0126346, crispri.BOOTSTRAPS = 200);
#: tests/test_crispri_published_v2.py checks every one of them against that file, so a mistyped
#: figure fails rather than being quoted. They are never recomputed here.
BEFORE: tuple[tuple[str, str, str, dict[str, Any]], ...] = (
    (
        "training_loco",
        "training_published_split/deletion_gain",
        "all 10,356 K562 training pairs, 471 regulated, hold-one-chromosome-out, unweighted",
        {"gain": 0.2173, "ci95": [0.1768, 0.2606], "resamples": 200},
    ),
    (
        "heldout_pooled",
        "heldout_published_pairs/deletion_gain",
        "all 4,378 held-out pairs pooled over five cell types, 190 regulated, weighted by "
        "direct-effect probability",
        {"gain": 0.1095, "ci95": [0.0607, 0.1741], "resamples": 200},
    ),
    (
        "heldout_GM12878",
        "heldout_published_pairs/per_cell_type_weighted/GM12878/deletion_gain",
        "the 68 GM12878 held-out pairs, 16 regulated, weighted; GM12878 has a deletion value",
        {"gain": 0.0146, "ci95": [-0.0941, 0.2052], "resamples": 200},
    ),
    (
        "heldout_K562",
        "heldout_published_pairs/per_cell_type_weighted/K562/deletion_gain",
        "the 1,918 K562 held-out pairs, 118 regulated, weighted; K562 has a deletion value",
        {"gain": 0.1361, "ci95": [0.0768, 0.2277], "resamples": 200},
    ),
    (
        "heldout_HCT116",
        "heldout_published_pairs/per_cell_type_weighted/HCT116/deletion_gain",
        "the 396 HCT116 held-out pairs, 34 regulated, weighted; HCT116 has no deletion value",
        {"gain": 0.0004, "ci95": [-0.0002, 0.004], "resamples": 200},
    ),
    (
        "heldout_Jurkat",
        "heldout_published_pairs/per_cell_type_weighted/Jurkat/deletion_gain",
        "the 75 Jurkat held-out pairs, 7 regulated, weighted; Jurkat has no deletion value",
        {"gain": 0.0038, "ci95": [0.0, 0.0059], "resamples": 200},
    ),
    (
        "heldout_WTC11",
        "heldout_published_pairs/per_cell_type_weighted/WTC11/deletion_gain",
        "the 1,921 WTC11 held-out pairs, 15 regulated, weighted; WTC11 has no deletion value",
        {"gain": -0.0049, "ci95": [-0.0174, 0.0], "resamples": 200},
    ),
    (
        "arm1_k562_all_pairs",
        "coverage_arms_k562_heldout/arm1_all_pairs",
        "the 1,918 K562 held-out pairs, 118 regulated, 174 uncovered, unweighted average precision",
        {"gain": 0.1361, "ci95": [0.081, 0.2248], "resamples": 200},
    ),
    (
        "arm3_k562_coverage_indicator_control",
        "coverage_arms_k562_heldout/arm3_coverage_indicator_control",
        "the same 1,918 K562 held-out pairs, the coverage indicator in place of the deletion "
        "features, unweighted average precision",
        {"gain": 0.0008, "ci95": [-0.0009, 0.0031], "resamples": 200},
    ),
    (
        "post_hoc_dnase_training_loco",
        "post_hoc_positive_filter/dnase_only_diagnostic/training_published_split/deletion_gain",
        "all 10,356 K562 training pairs, 471 regulated, hold-one-chromosome-out, unweighted, DNase "
        "in place of sqrt(DNase x H3K27ac); post hoc, not registered",
        {"gain": 0.2244, "ci95": [0.1819, 0.2676], "resamples": 200},
    ),
    (
        "post_hoc_dnase_heldout_pooled",
        "post_hoc_positive_filter/dnase_only_diagnostic/heldout_pooled_weighted/deletion_gain",
        "all 4,378 held-out pairs pooled over five cell types, 190 regulated, weighted, DNase only; "
        "post hoc, not registered",
        {"gain": 0.1636, "ci95": [0.1015, 0.2368], "resamples": 200},
    ),
    (
        "post_hoc_dnase_heldout_k562_covered",
        "post_hoc_positive_filter/dnase_only_diagnostic/heldout_k562_covered_headline_estimator/"
        "deletion_gain",
        "the 1,744 covered K562 held-out pairs, 114 regulated, unweighted average precision, DNase "
        "only; post hoc, not registered",
        {"gain": 0.1792, "ci95": [0.1023, 0.2867], "resamples": 200},
    ),
)

#: What the rebuild is expected to leave alone, transcribed from the same committed file: the point
#: AUPRCs of every model in every stratum, the bands against the published figures, the
#: coverage-matched arm (whose 1,000 draws are `crispri.MATCH_DRAWS`, which this run does not
#: change), the estimator check and the request count. The draws govern the intervals, not the
#: points, so a point that moved would be a finding; the block records each comparison either way.
UNCHANGED_BEFORE: dict[str, Any] = {
    "estimator_check_raw_distance/training_unweighted": 0.4359,
    "estimator_check_raw_distance/heldout_weighted": 0.3631,
    "training_published_split/models/distance/auprc": 0.4234,
    "training_published_split/models/activity + distance/auprc": 0.5068,
    "training_published_split/models/activity + distance + deletion/auprc": 0.7241,
    "training_published_split/against_encode_re2g": "above the published interval",
    "training_published_split/against_encode_re2g_extended": "inside the published interval",
    "training_published_split/baseline_against_abc": "below the published interval",
    "heldout_published_pairs/models/distance/auprc": 0.3631,
    "heldout_published_pairs/models/activity + distance/auprc": 0.5674,
    "heldout_published_pairs/models/activity + distance + deletion/auprc": 0.677,
    "heldout_published_pairs/against_encode_re2g": "above the published interval",
    "heldout_published_pairs/baseline_against_abc": "above the published interval",
    "heldout_published_pairs/pairs_without_a_deletion_value": 2392,
    "heldout_published_pairs/per_cell_type_weighted/GM12878/models/distance/auprc": 0.8203,
    "heldout_published_pairs/per_cell_type_weighted/GM12878/models/activity + distance/auprc": 0.7848,
    "heldout_published_pairs/per_cell_type_weighted/GM12878/models/"
    "activity + distance + deletion/auprc": 0.7994,
    "heldout_published_pairs/per_cell_type_weighted/HCT116/models/distance/auprc": 0.3167,
    "heldout_published_pairs/per_cell_type_weighted/HCT116/models/activity + distance/auprc": 0.4946,
    "heldout_published_pairs/per_cell_type_weighted/HCT116/models/"
    "activity + distance + deletion/auprc": 0.495,
    "heldout_published_pairs/per_cell_type_weighted/Jurkat/models/distance/auprc": 0.4465,
    "heldout_published_pairs/per_cell_type_weighted/Jurkat/models/activity + distance/auprc": 0.5748,
    "heldout_published_pairs/per_cell_type_weighted/Jurkat/models/"
    "activity + distance + deletion/auprc": 0.5786,
    "heldout_published_pairs/per_cell_type_weighted/K562/models/distance/auprc": 0.4083,
    "heldout_published_pairs/per_cell_type_weighted/K562/models/activity + distance/auprc": 0.5911,
    "heldout_published_pairs/per_cell_type_weighted/K562/models/activity + distance + deletion/auprc": 0.7272,
    "heldout_published_pairs/per_cell_type_weighted/WTC11/models/distance/auprc": 0.4188,
    "heldout_published_pairs/per_cell_type_weighted/WTC11/models/activity + distance/auprc": 0.5197,
    "heldout_published_pairs/per_cell_type_weighted/WTC11/models/"
    "activity + distance + deletion/auprc": 0.5148,
    "coverage_arms_k562_heldout/arm2_coverage_matched/regulated_coverage_before": 0.9661,
    "coverage_arms_k562_heldout/arm2_coverage_matched/non_regulated_coverage": 0.9056,
    "coverage_arms_k562_heldout/arm2_coverage_matched/regulated_kept": 107,
    "coverage_arms_k562_heldout/arm2_coverage_matched/regulated_covered": 114,
    "coverage_arms_k562_heldout/arm2_coverage_matched/draws": 1000,
    "coverage_arms_k562_heldout/arm2_coverage_matched/median_gain": 0.1449,
    "coverage_arms_k562_heldout/arm2_coverage_matched/range": [0.124, 0.1601],
    "coverage_arms_k562_heldout/arm2_coverage_matched/share_above_zero": 1.0,
    "post_hoc_positive_filter/dnase_only_diagnostic/training_published_split/models/"
    "dnase + distance/auprc": 0.4958,
    "post_hoc_positive_filter/dnase_only_diagnostic/training_published_split/models/"
    "dnase + distance + deletion/auprc": 0.7201,
    "post_hoc_positive_filter/dnase_only_diagnostic/heldout_pooled_weighted/models/"
    "dnase + distance/auprc": 0.4757,
    "post_hoc_positive_filter/dnase_only_diagnostic/heldout_pooled_weighted/models/"
    "dnase + distance + deletion/auprc": 0.6393,
    "post_hoc_positive_filter/dnase_only_diagnostic/heldout_pooled_weighted/"
    "against_encode_re2g": "above the published interval",
    "post_hoc_positive_filter/dnase_only_diagnostic/heldout_pooled_weighted/"
    "baseline_against_abc": "inside the published interval",
    "post_hoc_positive_filter/dnase_only_diagnostic/heldout_k562_covered_headline_estimator/models/"
    "dnase + distance/auprc": 0.5021,
    "post_hoc_positive_filter/dnase_only_diagnostic/heldout_k562_covered_headline_estimator/models/"
    "dnase + distance + deletion/auprc": 0.6813,
    "alphagenome_requests": 0,
}

MISSING = object()


def at(payload: Any, path: str) -> Any:
    """The value a slash-separated path names, or `MISSING` when the path is not there."""
    node = payload
    for key in path.split("/"):
        if not isinstance(node, dict) or key not in node:
            return MISSING
        node = node[key]
    return node


def _width(ci: Any) -> float | None:
    return round(ci[1] - ci[0], 4) if isinstance(ci, list) and len(ci) == 2 else None


def beside_the_old_result(result: dict[str, Any]) -> dict[str, Any]:
    """Every interval of this run beside the 2026-09-27 one, each with its own population named.

    One row is one population. The rows are not comparable with each other: the held-out pooled row
    covers all 4,378 held-out pairs and the K562 row covers the 1,918 K562 held-out pairs within
    them, so neither figure is evidence about the other's population.
    """
    rows: dict[str, Any] = {}
    for key, path, population, before in BEFORE:
        now = at(result, path)
        now = {} if now is MISSING or not isinstance(now, dict) else now
        kept = {
            k: now.get(k)
            for k in (
                "gain",
                "ci95",
                "resamples",
                "clusters",
                "draws_requested",
                "draws_dropped",
                "met_minimum",
                "unavailable",
            )
            if k in now
        }
        rows[key] = {
            "population": population,
            f"result_of_{OLD_DATE}": {**before, "draws": OLD_DRAWS, "file": OLD},
            "this_run": kept,
            "found_in_this_run": bool(kept),
            "gain_now_refused": bool(kept) and kept.get("gain") is None,
            "gain_moved": (
                None
                if kept.get("gain") is None
                else round(kept["gain"] - before["gain"], 4)
                if "gain" in kept
                else None
            ),
            "ci95_width": {
                f"result_of_{OLD_DATE}": _width(before.get("ci95")),
                "this_run": _width(kept.get("ci95")),
            },
        }
    return {
        "reading": "each row is one population; rows are never compared with one another. The "
        f"{OLD_DATE} figures are quoted from {OLD} at {OLD_DRAWS} draws, transcribed in this "
        "script's BEFORE and not recomputed. This run's figures rest on "
        f"{crispri.BOOTSTRAPS} requested draws.",
        "rows": rows,
        "gains_now_refused": sorted(k for k, v in rows.items() if v["gain_now_refused"]),
        "rows_not_found_in_this_run": sorted(k for k, v in rows.items() if not v["found_in_this_run"]),
    }


def what_did_not_change(result: dict[str, Any]) -> dict[str, Any]:
    """Each point AUPRC, band and count of the 2026-09-27 result set against this run's."""
    checked: dict[str, Any] = {}
    for path, before in UNCHANGED_BEFORE.items():
        now = at(result, path)
        checked[path] = {
            f"result_of_{OLD_DATE}": before,
            "this_run": None if now is MISSING else now,
            "present": now is not MISSING,
            "same": now is not MISSING and now == before,
        }
    moved = sorted(p for p, v in checked.items() if not v["same"])
    return {
        "reading": "the number of draws governs the intervals, not the points; these are the "
        f"{OLD_DATE} points, bands and counts checked one by one against this run",
        "fields_checked": len(checked),
        "fields_that_moved": moved,
        "all_unchanged": not moved,
        "match_draws_unchanged": crispri.MATCH_DRAWS,
        "fields": checked,
    }


def manifest(training: list[crispri.Pair], heldout: list[crispri.Pair]) -> dict[str, Any]:
    """The 2026-09-27 script's provenance contract, with this run's two additions.

    The sources, inputs, assembly, coordinates, parameters and partitions are the same contract, read
    from `scripts/crispri_published.py` rather than copied, so the two runs cannot drift apart in
    what they claim to have read. What is added is why this run exists, what a stratum without a
    deletion value now reports, and the import closure of this script.
    """
    m = dict(published_manifest(training, heldout))
    m["exclusions"] = [
        *m.get("exclusions", []),
        WHY_V2,
        "held-out "
        + ", ".join(WITHOUT_A_DELETION_VALUE)
        + " report no deletion gain at all: "
        + crispri.UNAVAILABLE_GAIN
        + ". Their model AUPRCs are still reported; it is the gain that is refused",
    ]
    m["parameters"] = {
        **m.get("parameters", {}),
        "min_resamples": crispri.MIN_RESAMPLES,
        "interval_provenance": "every interval states clusters, draws_requested, draws_dropped and "
        "met_minimum (crispri._interval_provenance)",
        "gain_where_unavailable": crispri.UNAVAILABLE_GAIN,
    }
    m["import_closure"] = closure_block()
    m["supersedes"] = {
        "file": OLD,
        "date": OLD_DATE,
        "kept": "unchanged; this run is written beside it under a new name, not over it",
        "why": WHY_V2,
    }
    return m


def main() -> int:
    t0 = time.time()
    for name in (crispri.TRAINING, crispri.HELDOUT):
        crispri.fetch(name)
    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    result = crispri.score_published(training, heldout, crispri.DeletionTable(), crispri.ElementCache())
    result["beside_the_2026_09_27_result"] = beside_the_old_result(result)
    result["what_did_not_change"] = what_did_not_change(result)
    result[mf.KEY] = manifest(training, heldout)  # save_result takes it from the payload
    path = save_result(NAME, result)

    tr, ho = result["training_published_split"], result["heldout_published_pairs"]
    print(f"training LOCO, all 10,356 training pairs: gain {json.dumps(tr['deletion_gain'])}")
    print(f"held out, all 4,378 pairs pooled weighted: gain {json.dumps(ho['deletion_gain'])}")
    for cell, c in ho["per_cell_type_weighted"].items():
        g = c["deletion_gain"]
        print(
            f"  {cell:8s} available {str(c['deletion_available']):5s} "
            f"gain {g['gain']} ci {g['ci95']} resamples {g['resamples']}"
        )
    beside = result["beside_the_2026_09_27_result"]
    print(f"gains now refused: {beside['gains_now_refused']}")
    for key, row in beside["rows"].items():
        print(
            f"  {key:40s} {row['population'][:52]:52s} "
            f"{OLD_DATE} {row[f'result_of_{OLD_DATE}']['gain']} {row[f'result_of_{OLD_DATE}']['ci95']} "
            f"-> now {row['this_run'].get('gain')} {row['this_run'].get('ci95')}"
        )
    did_not = result["what_did_not_change"]
    print(f"unchanged: {did_not['fields_checked']} fields checked, moved {did_not['fields_that_moved']}")
    closure = result[mf.KEY]["import_closure"]
    print(
        f"import closure: {closure['count']} files, uncommitted on it "
        f"{closure['uncommitted_code_on_the_closure']}"
    )
    print(f"({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
