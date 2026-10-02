# SPDX-License-Identifier: AGPL-3.0-or-later
"""The CRISPRi benchmark result rebuilt under a new name: 2,000 draws, every interval stating the
clusters it rests on, and no interval where the clusters are too few to read one.

    uv run python scripts/crispri_benchmark_v2.py

Writes `data/results/crispri_benchmark_v2.json`. The committed result,
`data/results/crispri_benchmark.json` (recorded date 2026-09-22), is left exactly as it is: it stays
the record of what was computed and quoted then, and this file is the corrected run beside it, not a
replacement of it.

This is the file README's one CRISPRi sentence rests on, through docs/CRISPRI-RESULT.md, and it was
the last result still carrying intervals made under the rules of before 2026-10-02. It holds five
`deletion_gain` intervals, every one at `resamples: 200` and none recording a cluster count, so no
reader could tell whether any of them met the minimum number of resampling units an interval needs.
The defect class is the one `data/results/crispri_published_v2.json` was rebuilt for; this file was
left out of that work.

What this run changes, all of it already in `genomeos/attribution/crispri.py`:

1. `crispri.BOOTSTRAPS` is 2,000, not 200, and `crispri._interval_provenance` states beside every
   interval its `clusters`, `draws_requested`, `draws_dropped`, `met_minimum` and `enough_clusters`;
2. an interval resting on fewer than `crispri.MIN_CLUSTERS_FOR_AN_INTERVAL` chromosomes is not
   reported at all: `ci95` is `null` and `interval_unreliable` says why. The point estimate stands,
   because it does not depend on the resampling, and no replacement interval is written by hand;
3. every stratum's gain goes through `crispri.stratum_gain`, which is `crispri.gain_where_available`:
   a stratum whose pairs carry no deletion value reports `gain: null` with `crispri.UNAVAILABLE_GAIN`,
   and an empty stratum reports `gain: null` with `crispri.NO_PAIRS_GAIN`. Neither is ever a number.

Nothing else about the measurement changes. The features, the fit, the leave-one-chromosome-out
scheme, the frozen held-out weights, the estimator, the seed and the pair sets are those of the
committed run, and no AlphaGenome request is made: every deletion value is read from the sweep's
existing per-element cache, as `alphagenome_requests: 0` says.

Two blocks are added for the reader:

* `beside_the_committed_result`: each of the five intervals with its population named, the committed
  figure at 200 draws quoted (transcribed in `BEFORE`, never recomputed) and this run's figure beside
  it. A row holds one population and the rows are never compared with one another: "the 1,744 covered
  K562 held-out pairs" and "the 62 covered GM12878 held-out pairs" are different sets of pairs, and
  `population_unchanged` says, per row, whether the committed row and this run's row are even about
  the same pairs.
* `what_did_not_change`: the points, counts, weights, refusals and precisions of the committed result
  set one by one against this run's, so what the rebuild left alone is on the record beside what it
  moved.

`result_manifest.code_cleanliness` answers the shared checkout: several lanes hold uncommitted files
in this working tree, so the counting path of this script is computed by `import_closure` and set
against the uncommitted code git reports, rather than asserted by hand. `own_code_is_committed` and
`foreign_uncommitted_code_on_the_counting_path` are the two `scripts/manifest_rebuild.py` requires to
hold on both sides.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from scripts.crispri_published import invalid_connections  # noqa: E402
from scripts.crispri_published_v2 import MISSING, at, import_closure  # noqa: E402

NAME = "crispri_benchmark_v2"
ENTRY = "scripts/crispri_benchmark_v2.py"
TESTS = "tests/test_crispri_benchmark_v2.py"
OLD = "data/results/crispri_benchmark.json"
OLD_DATE = "2026-09-22"  # the `date` the committed file itself records
OLD_DRAWS = 200
WRITTEN = "2026-10-02"

#: What this run changes about the committed run, as the manifest's exclusions state it.
WHY_V2 = (
    f"rebuilt under a new name on {WRITTEN}: {crispri.BOOTSTRAPS} chromosome resamples per interval in "
    f"place of {OLD_DRAWS}, each interval stating the clusters and draws it rests on, no interval at all "
    f"below {crispri.MIN_CLUSTERS_FOR_AN_INTERVAL} clusters, and no deletion gain where a stratum carries "
    f"no deletion value; the committed result at {OLD} (date {OLD_DATE}) is kept unchanged beside this one"
)

#: The files this lane holds. Everything else uncommitted in this checkout belongs to another lane, and
#: `code_cleanliness` reports the two apart rather than together.
OWN_CODE = frozenset({ENTRY, TESTS, "genomeos/attribution/crispri.py"})


# --- the committed figures, transcribed ----------------------------------------------------------

#: The five `deletion_gain` intervals of the committed `data/results/crispri_benchmark.json`, each with the
#: population it was measured on. `path` reads the same field of this run; `pairs_path` and
#: `positives_path` read the pair and positive counts that say which pairs the row is about, so the row can
#: state whether the two sides share a population at all. The figures and the counts are transcribed from
#: the committed file, which carries no `result_manifest` and was written when `crispri.BOOTSTRAPS` was
#: 200; `tests/test_crispri_benchmark_v2.py` checks every one of them against that file, so a mistyped
#: figure fails rather than being quoted. None of them is recomputed here.
BEFORE: tuple[dict[str, Any], ...] = (
    {
        "key": "training_leave_chromosome_out",
        "path": "training_leave_chromosome_out/deletion_gain",
        "pairs_path": "training_leave_chromosome_out/models/distance/pairs",
        "positives_path": "training_leave_chromosome_out/models/distance/positives",
        "population": "the 9,237 covered K562 training pairs, 451 regulated, "
        "hold-one-chromosome-out, unweighted average precision",
        "pairs": 9237,
        "positives": 451,
        "before": {"gain": 0.2282, "ci95": [0.1854, 0.2772], "resamples": 200},
    },
    {
        "key": "heldout_K562",
        "path": "heldout/K562/deletion_gain",
        "pairs_path": "heldout/K562/models/distance/pairs",
        "positives_path": "heldout/K562/models/distance/positives",
        "population": "the 1,744 covered K562 held-out pairs, 114 regulated, frozen weights, "
        "unweighted average precision. This is the figure nearest README",
        "pairs": 1744,
        "positives": 114,
        "before": {"gain": 0.1407, "ci95": [0.082, 0.2313], "resamples": 200},
    },
    {
        "key": "heldout_GM12878",
        "path": "heldout/GM12878/deletion_gain",
        "pairs_path": "heldout/GM12878/models/distance/pairs",
        "positives_path": "heldout/GM12878/models/distance/positives",
        "population": "the 62 covered GM12878 held-out pairs, 14 regulated, frozen weights, "
        "unweighted average precision",
        "pairs": 62,
        "positives": 14,
        "before": {"gain": 0.0352, "ci95": [-0.0261, 0.1686], "resamples": 200},
    },
    {
        "key": "heldout_elements_not_in_training_K562",
        "path": "heldout_elements_not_in_training/K562/deletion_gain",
        "pairs_path": "heldout_elements_not_in_training/K562/models/distance/pairs",
        "positives_path": "heldout_elements_not_in_training/K562/models/distance/positives",
        "population": "the 1,580 covered K562 held-out pairs on elements that overlap no training "
        "element, 91 regulated, frozen weights, unweighted average precision",
        "pairs": 1580,
        "positives": 91,
        "before": {"gain": 0.1489, "ci95": [0.0761, 0.2484], "resamples": 200},
    },
    {
        "key": "heldout_elements_not_in_training_GM12878",
        "path": "heldout_elements_not_in_training/GM12878/deletion_gain",
        "pairs_path": "heldout_elements_not_in_training/GM12878/models/distance/pairs",
        "positives_path": "heldout_elements_not_in_training/GM12878/models/distance/positives",
        "population": "the 62 covered GM12878 held-out pairs on elements that overlap no training "
        "element, 14 regulated, frozen weights, unweighted average precision",
        "pairs": 62,
        "positives": 14,
        "before": {"gain": 0.0352, "ci95": [-0.0261, 0.1686], "resamples": 200},
    },
)

#: What the rebuild is expected to leave alone, transcribed from the same committed file: every model's
#: point AUPRC and the AUROCs of the leave-chromosome-out arm, the coverage counts, the fitted weights of
#: the deletion model, the three one-call-per-element thresholds, the element-level arms, the deletion
#: census and the three refusals. The number of draws governs the intervals, not the points, so a point
#: that moved would be a finding; `what_did_not_change` records each comparison either way.
UNCHANGED_BEFORE: dict[str, Any] = {
    "verdict": "passed",
    "deletion_reader": "the sweep's per-element response cache (data/knowledge/alphagenome/elements), "
    "every gene in the scorer's 1 Mb window on the cell's own track",
    "coverage/training_pairs": 10356,
    "coverage/training_pairs_on_a_deleted_element": 9237,
    "coverage/training_positives_on_a_deleted_element": 451,
    "coverage/heldout_pairs": 4378,
    "coverage/heldout_pairs_on_a_deleted_element": 3849,
    "training_leave_chromosome_out/models/distance/auprc": 0.4304,
    "training_leave_chromosome_out/models/distance/auroc": 0.8925,
    "training_leave_chromosome_out/models/activity + distance/auprc": 0.5113,
    "training_leave_chromosome_out/models/activity + distance/auroc": 0.9215,
    "training_leave_chromosome_out/models/activity + distance + deletion/auprc": 0.7395,
    "training_leave_chromosome_out/models/activity + distance + deletion/auroc": 0.9428,
    "training_single_predictors/distance/auprc": 0.4413,
    "training_single_predictors/activity over distance/auprc": 0.519,
    "training_single_predictors/deletion (top target x drop)/auprc": 0.4647,
    "heldout/GM12878/models/distance/auprc": 0.8207,
    "heldout/GM12878/models/activity + distance/auprc": 0.8648,
    "heldout/GM12878/models/activity + distance + deletion/auprc": 0.9001,
    "heldout/GM12878/models/distance/pairs": 62,
    "heldout/GM12878/models/distance/positives": 14,
    "heldout/GM12878/passes": True,
    "heldout/K562/models/distance/auprc": 0.3753,
    "heldout/K562/models/activity + distance/auprc": 0.5501,
    "heldout/K562/models/activity + distance + deletion/auprc": 0.6909,
    "heldout/K562/models/distance/pairs": 1744,
    "heldout/K562/models/distance/positives": 114,
    "heldout/K562/passes": True,
    "heldout/HCT116/refused": "no AlphaGenome line for this cell type in the deletion table",
    "heldout/Jurkat/refused": "no AlphaGenome line for this cell type in the deletion table",
    "heldout/WTC11/refused": "no AlphaGenome line for this cell type in the deletion table",
    "heldout_elements_not_in_training/GM12878/models/distance/auprc": 0.8207,
    "heldout_elements_not_in_training/GM12878/models/activity + distance/auprc": 0.8648,
    "heldout_elements_not_in_training/GM12878/models/activity + distance + deletion/auprc": 0.9001,
    "heldout_elements_not_in_training/GM12878/models/distance/pairs": 62,
    "heldout_elements_not_in_training/GM12878/models/distance/positives": 14,
    "heldout_elements_not_in_training/K562/models/distance/auprc": 0.3365,
    "heldout_elements_not_in_training/K562/models/activity + distance/auprc": 0.4971,
    "heldout_elements_not_in_training/K562/models/activity + distance + deletion/auprc": 0.646,
    "heldout_elements_not_in_training/K562/models/distance/pairs": 1580,
    "heldout_elements_not_in_training/K562/models/distance/positives": 91,
    "weights/activity + distance + deletion/intercept": 5.2275,
    "weights/activity + distance + deletion/log_distance": -0.3231,
    "weights/activity + distance + deletion/log_activity": 0.2465,
    "weights/activity + distance + deletion/activity_over_distance": 0.5696,
    "weights/activity + distance + deletion/top_target": 1.3287,
    "weights/activity + distance + deletion/deletion_drop": 25.2989,
    "one_call_per_element_k562/0.0/all_elements/deletion/calls": 225,
    "one_call_per_element_k562/0.0/all_elements/deletion/right": 179,
    "one_call_per_element_k562/0.0/all_elements/deletion/precision": 0.7956,
    "one_call_per_element_k562/0.1/all_elements/deletion/calls": 128,
    "one_call_per_element_k562/0.1/all_elements/deletion/right": 123,
    "one_call_per_element_k562/0.1/all_elements/deletion/precision": 0.9609,
    "one_call_per_element_k562/0.2/all_elements/deletion/calls": 86,
    "one_call_per_element_k562/0.2/all_elements/deletion/right": 86,
    "one_call_per_element_k562/0.2/all_elements/deletion/precision": 1.0,
    "one_call_per_element_k562/0.1/elements": 3472,
    "one_call_per_element_k562/0.1/elements_with_a_regulated_gene": 400,
    "element_level/K562 training/activity/auprc": 0.2798,
    "element_level/K562 training/closest distance/auprc": 0.5114,
    "element_level/K562 training/deletion drop/auprc": 0.7526,
    "element_level/GM12878 held out/activity/auprc": 0.6693,
    "element_level/GM12878 held out/closest distance/auprc": 0.8207,
    "element_level/GM12878 held out/deletion drop/auprc": 0.6465,
    "element_level/K562 held out/activity/auprc": 0.4077,
    "element_level/K562 held out/closest distance/auprc": 0.4359,
    "element_level/K562 held out/deletion drop/auprc": 0.7608,
    "deletion_census/training/K562/covered": 9237,
    "deletion_census/heldout/K562/covered": 1744,
    "deletion_census/training/K562/regulated": 451,
    "deletion_census/heldout/K562/regulated": 114,
    "deletion_census/training/K562/gene_is_the_top_target": 246,
    "deletion_census/heldout/K562/gene_is_the_top_target": 40,
    "deletion_census/training/K562/structural_zero_answered": 5581,
    "deletion_census/heldout/K562/structural_zero_answered": 1112,
    "deletion_census/training/K562/share_answered": 0.6207,
    "deletion_census/heldout/K562/share_answered": 0.6526,
    "coverage_by_arm/training/arms/regulated/covered": 451,
    "coverage_by_arm/training/arms/regulated/fraction": 0.9575,
    "coverage_by_arm/training/elements_fully_covered": 3472,
    "coverage_by_arm/heldout/arms/regulated/covered": 182,
    "coverage_by_arm/heldout/arms/regulated/fraction": 0.9579,
    "coverage_by_arm/heldout/elements_fully_covered": 1489,
}

#: The interval fields `_interval_provenance` writes, plus the refusal field, in the order a reader wants
#: them. Only the ones a row actually carries are quoted back.
INTERVAL_FIELDS = (
    "gain",
    "ci95",
    "resamples",
    "clusters",
    "draws_requested",
    "draws_dropped",
    "met_minimum",
    "enough_clusters",
    "interval_unreliable",
    "unavailable",
)


# --- the shared checkout this run happened in ----------------------------------------------------


def code_cleanliness(root: Path = ROOT) -> dict[str, Any]:
    """Which uncommitted code this shared checkout held, split into this lane's and other lanes', and
    whether any of it is on the counting path. Read from git and from the imports, never asserted.

    The field names are the convention `scripts/cell2_eligibility.py` set and
    `scripts/manifest_rebuild.py` matches, so `own_code_is_committed` and
    `foreign_uncommitted_code_on_the_counting_path` must hold on both sides of a rebuild while the lists
    that merely describe the tree a run happened in may take their clean-worktree values.
    """
    rev = mf.code_revision(root)
    path = import_closure(ENTRY, root)
    dirty = list(rev.get("dirty_code_paths") or [])
    own = [p for p in dirty if p in OWN_CODE]
    foreign = [p for p in dirty if p not in OWN_CODE]
    return {
        "git_sha": rev.get("git_sha"),
        "dirty": rev.get("dirty"),
        "own_uncommitted_code": own,
        "own_code_is_committed": not own,
        "foreign_uncommitted_code": foreign,
        "foreign_uncommitted_code_on_the_counting_path": [p for p in foreign if p in path],
        "counting_path": path,
        "counting_path_count": len(path),
        "counting_path_is_computed": (
            "the transitive import closure of this script, computed from the files' import statements at "
            "write time (scripts.crispri_published_v2.import_closure); not a hand list"
        ),
        "note": (
            "several sessions work in this one checkout. A file under foreign_uncommitted_code belongs to "
            "another lane; this lane did not write it and did not commit it. A foreign file outside the "
            "counting path cannot have entered a number here, and "
            "foreign_uncommitted_code_on_the_counting_path names any that could"
        ),
    }


# --- the committed result beside this one --------------------------------------------------------


def _width(ci: Any) -> float | None:
    return round(ci[1] - ci[0], 4) if isinstance(ci, list) and len(ci) == 2 else None


def beside_the_committed_result(result: dict[str, Any]) -> dict[str, Any]:
    """The five committed intervals beside this run's, each row naming its own population.

    One row is one population, and no row is compared with another: the leave-chromosome-out row covers
    the 9,237 covered K562 training pairs and the held-out K562 row covers 1,744 quite different pairs,
    so neither figure is evidence about the other's population. Within a row the two sides are set against
    each other only when `population_unchanged` holds, which is read from the pair and positive counts
    this run reports rather than assumed; where it does not hold, no difference is taken.
    """
    rows: dict[str, Any] = {}
    for spec in BEFORE:
        now = at(result, spec["path"])
        now = {} if now is MISSING or not isinstance(now, dict) else now
        kept = {k: now[k] for k in INTERVAL_FIELDS if k in now}
        pairs_now, positives_now = at(result, spec["pairs_path"]), at(result, spec["positives_path"])
        same_population = pairs_now == spec["pairs"] and positives_now == spec["positives"]
        row = {
            "population": spec["population"],
            "population_unchanged": same_population,
            "population_counts": {
                f"result_of_{OLD_DATE}": {"pairs": spec["pairs"], "positives": spec["positives"]},
                "this_run": {
                    "pairs": None if pairs_now is MISSING else pairs_now,
                    "positives": None if positives_now is MISSING else positives_now,
                },
            },
            f"result_of_{OLD_DATE}": {
                **spec["before"],
                "draws": OLD_DRAWS,
                "file": OLD,
                "clusters_recorded": False,
                "readable_against_the_cluster_rule": (
                    "no: the committed row records no cluster count, so whether it reached "
                    f"crispri.MIN_CLUSTERS_FOR_AN_INTERVAL = {crispri.MIN_CLUSTERS_FOR_AN_INTERVAL} "
                    "cannot be told from the file"
                ),
            },
            "this_run": kept,
            "found_in_this_run": bool(kept),
            "gain_now_refused": bool(kept) and kept.get("gain") is None,
            "interval_now_withdrawn": bool(kept) and kept.get("gain") is not None and not kept.get("ci95"),
            "ci95_width": {
                f"result_of_{OLD_DATE}": _width(spec["before"].get("ci95")),
                "this_run": _width(kept.get("ci95")),
            },
        }
        if same_population and isinstance(kept.get("gain"), int | float):
            row["gain_moved"] = round(kept["gain"] - spec["before"]["gain"], 4)
        else:
            row["gain_not_compared"] = (
                "the two rows are not about the same pairs, so their gains are not set against each other"
                if not same_population
                else "this run reports no gain for this stratum"
            )
        rows[spec["key"]] = row
    return {
        "reading": (
            f"each row is one population; rows are never compared with one another. The {OLD_DATE} "
            f"figures are quoted from {OLD} at {OLD_DRAWS} draws, transcribed in this script's BEFORE and "
            f"not recomputed. This run's figures rest on {crispri.BOOTSTRAPS} requested draws, and an "
            f"interval below {crispri.MIN_CLUSTERS_FOR_AN_INTERVAL} clusters is withdrawn rather than "
            "replaced: the point estimate stands, the interval is not reported"
        ),
        "rows": rows,
        "intervals_in_the_committed_result": len(BEFORE),
        "intervals_now_withdrawn": sorted(k for k, v in rows.items() if v["interval_now_withdrawn"]),
        "intervals_that_survive": sorted(k for k, v in rows.items() if v["this_run"].get("ci95") is not None),
        "gains_now_refused": sorted(k for k, v in rows.items() if v["gain_now_refused"]),
        "populations_that_changed": sorted(k for k, v in rows.items() if not v["population_unchanged"]),
        "rows_not_found_in_this_run": sorted(k for k, v in rows.items() if not v["found_in_this_run"]),
    }


def what_did_not_change(result: dict[str, Any]) -> dict[str, Any]:
    """Each point, count, weight, refusal and precision of the committed result set against this run's."""
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
        "reading": (
            "the number of draws governs the intervals, not the points; these are the committed result's "
            "points, counts, weights, refusals and precisions checked one by one against this run"
        ),
        "fields_checked": len(checked),
        "fields_that_moved": moved,
        "all_unchanged": not moved,
        "fields": checked,
    }


# --- provenance ----------------------------------------------------------------------------------


@mf.depends_on_models("alphagenome")  # the sweep's model, as far as the disk says (R9)
def manifest(training: list[crispri.Pair], heldout: list[crispri.Pair]) -> dict[str, Any]:
    """What this result read, which bytes, and under which rules its intervals were made.

    The pair tables and the two deletion stores are pinned by sha256, so a rebuild in a clean worktree
    either reads the same bytes or says which input it could not find.
    """
    inputs = [
        mf.input_entry(crispri.KNOWLEDGE / name, partition=CRISPRI_SPLIT_OF[name], pairs=len(pairs))
        for name, pairs in ((crispri.TRAINING, training), (crispri.HELDOUT, heldout))
    ]
    inputs.append(mf.input_entry(crispri.ELEMENTS, partition=None, role="DeletionTable: elements by overlap"))
    inputs.append(
        mf.input_entry(
            crispri.ELEMENT_CACHE, partition=None, role="ElementCache: every gene's deletion value"
        )
    )
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison resources/crispr_data, EPCrisprBenchmark "
                "training_K562 and heldout_5_cell_types (Gschwind et al.)",
                "version": "main branch, unpinned upstream; fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion (all-element table and "
                "per-element response cache)",
                "version": "AlphaGenome as served during the 2026-09 all-element sweep (unpinned); "
                "pinned here by sha256",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "features": {k: list(v) for k, v in crispri.FEATURES.items()},
            "fit": "crispri.logistic_fit, lam 1e-3, 25 rounds, on the covered training pairs",
            "training_scoring": "hold-one-chromosome-out over the covered training pairs",
            "heldout_scoring": "frozen weights, per cell type, covered pairs only, unweighted",
            "estimator": "crispri.average_precision (step-wise, unweighted); the published-benchmark "
            "estimator crispri.benchmark_auprc belongs to crispri_published_v2, not to this result",
            "model_cells": list(crispri.MODEL_CELLS),
            "min_distance": crispri.MIN_DISTANCE,
            "reach": crispri.REACH,
            "bootstraps": crispri.BOOTSTRAPS,
            "bootstrap_seed": 0,
            "min_resamples": crispri.MIN_RESAMPLES,
            "min_clusters_for_an_interval": crispri.MIN_CLUSTERS_FOR_AN_INTERVAL,
            "resampling_unit": "whole chromosomes, percentile bootstrap",
            "interval_provenance": "every interval states clusters, draws_requested, draws_dropped, "
            "met_minimum and enough_clusters (crispri._interval_provenance)",
            "gain_where_unavailable": crispri.UNAVAILABLE_GAIN,
            "gain_where_no_pairs": crispri.NO_PAIRS_GAIN,
            "one_call_thresholds": [0.0, 0.1, 0.2],
            "alphagenome_requests": 0,
        },
        "exclusions": [
            {
                "file": name,
                "dropped_by_valid_connection": invalid_connections(name),
                "why": "the benchmark marks these pairs as no test of an enhancer-gene link",
            }
            for name in (crispri.TRAINING, crispri.HELDOUT)
        ]
        + [
            WHY_V2,
            "a pair no deleted element overlaps is not scored in any arm here: every arm is the covered "
            "pairs of its stratum, and the coverage is reported beside it",
            "held-out HCT116, Jurkat and WTC11 are refused outright rather than scored with deletion "
            "features of zero: the sweep kept four cell lines and these three carry no deletion value, so "
            + crispri.UNAVAILABLE_GAIN,
            "an interval resting on fewer than "
            f"{crispri.MIN_CLUSTERS_FOR_AN_INTERVAL} chromosomes is not reported at all; the point "
            "estimate is, because it does not depend on the resampling, and no replacement interval is "
            "written by hand",
            "per-cell held-out blocks for a cell with no regulated pair are omitted",
        ],
        "partitions": {
            crispri.TRAINING: CRISPRI_SPLIT_OF[crispri.TRAINING],
            crispri.HELDOUT: CRISPRI_SPLIT_OF[crispri.HELDOUT],
            "training": "fitted, and scored by hold-one-chromosome-out "
            "(training_leave_chromosome_out, training_single_predictors)",
            "heldout": "evaluation only, frozen weights (heldout, heldout_elements_not_in_training)",
        },
        "code_cleanliness": code_cleanliness(),
        "supersedes": {
            "file": OLD,
            "date": OLD_DATE,
            "kept": "unchanged; this run is written beside it under a new name, not over it",
            "why": WHY_V2,
            "readers_of_the_old_file": "README's one CRISPRi sentence, through docs/CRISPRI-RESULT.md",
        },
    }


def _line(label: str, g: dict[str, Any]) -> str:
    return (
        f"{label:48s} gain {g.get('gain')} ci {g.get('ci95')} clusters {g.get('clusters')} "
        f"draws {g.get('resamples')}/{g.get('draws_requested')} enough {g.get('enough_clusters')}"
    )


def main() -> int:
    t0 = time.time()
    for name in (crispri.TRAINING, crispri.HELDOUT):
        crispri.fetch(name)
    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    result = crispri.score(training, heldout, crispri.DeletionTable(), crispri.ElementCache())
    result["beside_the_committed_result"] = beside_the_committed_result(result)
    result["what_did_not_change"] = what_did_not_change(result)
    result[mf.KEY] = manifest(training, heldout)  # save_result takes it from the payload
    path = save_result(NAME, result)

    print(f"coverage: {result['coverage']}")
    print(
        _line("K562 training, leave-chromosome-out", result["training_leave_chromosome_out"]["deletion_gain"])
    )
    for cell, h in result["heldout"].items():
        if "refused" in h:
            print(f"  held out {cell}: refused ({h['refused']})")
            continue
        print("  " + _line(f"held out {cell}", h["deletion_gain"]) + f" passes {h['passes']}")
    for cell, h in result["heldout_elements_not_in_training"].items():
        print("  " + _line(f"held out {cell}, elements not in training", h["deletion_gain"]))
    beside = result["beside_the_committed_result"]
    print(f"intervals withdrawn: {beside['intervals_now_withdrawn']}")
    print(f"intervals that survive: {beside['intervals_that_survive']}")
    print(f"gains refused: {beside['gains_now_refused']}")
    print(f"populations that changed: {beside['populations_that_changed']}")
    did_not = result["what_did_not_change"]
    print(f"unchanged: {did_not['fields_checked']} fields checked, moved {did_not['fields_that_moved']}")
    clean = result[mf.KEY]["code_cleanliness"]
    print(
        f"counting path: {clean['counting_path_count']} files; own code committed "
        f"{clean['own_code_is_committed']}; foreign uncommitted on the path "
        f"{clean['foreign_uncommitted_code_on_the_counting_path']}"
    )
    print(f"verdict: {result['verdict']}  ({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
