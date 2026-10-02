# SPDX-License-Identifier: AGPL-3.0-or-later
"""Gate first, then the paired comparison against ENCODE-rE2G on identical held-out pairs
(data/results/re2g_paired.json).

    uv run --frozen python scripts/re2g_paired.py

Applies the specification registered in `genomeos/attribution/re2g.py` and committed before any
comparator score was read. Nothing here changes it. The order is the registered order and is not
negotiable inside this script:

1. join the five ENCODE rE2G prediction files to the 4,378 held-out pairs by the benchmark's own rule,
   unpredicted pairs at 0, scores aggregated by `sum`;
2. run the gate on the pooled population: reproduce the published 0.556151 within 0.02 or stop;
3. only if the gate passes, score our own frozen model and report the paired delta per population.

A failed gate writes its miss and exits non-zero without computing any delta, so a failure cannot be
mistaken for a comparison. 0 model requests: our own scores come from the deletion table and element
cache already on disk.
"""

from __future__ import annotations

import ast
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, re2g  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "re2g_paired"
PREDICTIONS = Path("data/cache/re2g/predictions")
PORTAL = "https://www.encodeproject.org/files/{a}/@@download/{a}.bed.gz"
HEADLINE = Path("data/results/crispri_published.json")
ROOT = Path(__file__).resolve().parents[1]
PACKAGE = "genomeos"
#: This lane's own files, so uncommitted code can be split into ours and other sessions'.
OWN_CODE = (
    "genomeos/attribution/re2g.py",
    "scripts/re2g_register.py",
    "scripts/re2g_paired.py",
    "scripts/re2g_k562_only_range.py",
    "scripts/re2g_identity.py",
    "scripts/re2g_likeforlike.py",
    "tests/test_re2g.py",
    "tests/test_re2g_paired.py",
)
#: The commit that fixed the specification, before any comparator score was read.
REGISTRATION_COMMIT = "d2ca94d"
REGISTRATION_FILE = "data/results/re2g_registration.json"


def counting_path(entry: Path | None = None) -> list[str]:
    """Every module of this repository the entry script can reach by import, as the transitive closure
    of its import statements, read from the files' syntax at write time rather than listed by hand."""
    entry = entry or Path(__file__).resolve()
    seen: dict[str, Path] = {}
    queue = [entry]
    while queue:
        path = queue.pop()
        rel = str(path.relative_to(ROOT))
        if rel in seen:
            continue
        seen[rel] = path
        try:
            tree = ast.parse(path.read_text())
        except (OSError, SyntaxError):
            continue
        names: set[str] = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names.update(a.name for a in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module and not node.level:
                names.add(node.module)
                names.update(f"{node.module}.{a.name}" for a in node.names)
        for name in names:
            if not (name == PACKAGE or name.startswith(PACKAGE + ".")):
                continue
            parts = name.split(".")
            for candidate in ([*parts[:-1], f"{parts[-1]}.py"], [*parts, "__init__.py"]):
                if _is_file_exactly(ROOT, candidate):
                    queue.append(ROOT.joinpath(*candidate))
    return sorted(seen)


def _is_file_exactly(root: Path, parts: list[str]) -> bool:
    """A file at `parts` below `root`, spelled as the directories spell it.

    Carries lane-rebuild's fix (`cd263bc`). `Path.is_file` on a case-insensitive filesystem answers yes
    for `genomeos/genome/Genome.py` when only `genome.py` is there, so `from genomeos.genome import
    Genome` -- a class, not a module -- would otherwise put a file that does not exist on the closure,
    and the closure would differ between this machine and a case-sensitive one.
    """
    node = root
    for part in parts:
        try:
            if part not in {q.name for q in node.iterdir()}:
                return False
        except OSError:
            return False
        node = node / part
    return node.is_file()


def code_cleanliness() -> dict[str, Any]:
    """Which uncommitted code the stamp names, split into this lane's and other sessions', and whether
    any of it lies on the computed counting path. Read from git and from the imports, not asserted."""
    rev = mf.code_revision()
    path = counting_path()
    dirty = list(rev.get("dirty_code_paths") or [])
    own = [p for p in dirty if p in OWN_CODE]
    foreign = [p for p in dirty if p not in OWN_CODE]
    return {
        "git_sha": rev["git_sha"],
        "dirty": rev["dirty"],
        "own_uncommitted_code": own,
        "own_code_is_committed": not own,
        "foreign_uncommitted_code": foreign,
        "foreign_uncommitted_code_on_the_counting_path": [p for p in foreign if p in path],
        "counting_path": path,
        "counting_path_is_computed": (
            "the transitive import closure of this script over the repository's own package, computed "
            "from the files' import statements at write time, not a hand-written list"
        ),
        "the_registrations_own_dirt": {
            "result": REGISTRATION_FILE,
            "named": ["genomeos/attribution/experiment_value.py", "tests/test_experiment_value.py"],
            "whose": "lane-s8, untracked in this shared checkout when the registration was stamped",
            "since": "committed at bee8f20",
            "why_it_could_not_enter": (
                "neither path is on the counting path computed above, and the registration reads no "
                "value from either; the registration writes no measurement at all"
            ),
        },
    }


def registration_unchanged() -> dict[str, Any]:
    """Compare every registered element against the committed registration, rather than asserting it.

    Reads `REGISTRATION_FILE` as it stands at `REGISTRATION_COMMIT` -- the commit made before any
    comparator score was read -- and checks the live module's constants against it field by field. The
    additions since that commit are listed by name, and the diffstat is read from git.
    """
    blob = subprocess.run(
        ["git", "show", f"{REGISTRATION_COMMIT}:{REGISTRATION_FILE}"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    ).stdout
    was = json.loads(blob)
    checks = {
        "gate_target": (was["gate_target"]["weighted_auprc"], re2g.GATE_TARGET),
        "gate_tolerance": (was["gate_target"]["tolerance"], re2g.GATE_TOLERANCE),
        "gate_interval": (was["gate_target"]["interval"], list(re2g.GATE_INTERVAL)),
        "join_rule": (was["join_rule"], dict(re2g.JOIN_RULE)),
        "populations": (was["populations"], {k: dict(v) for k, v in re2g.POPULATIONS.items()}),
        "draws_requested": (was["statistic_parameters"]["draws_requested"], re2g.DRAWS),
        "seed": (was["statistic_parameters"]["seed"], re2g.SEED),
        "cluster": (was["statistic_parameters"]["cluster"], re2g.CLUSTER),
        "aggregate_function": (was["result_manifest"]["parameters"]["aggregate_function"], "sum"),
        "fill_value": (was["result_manifest"]["parameters"]["fill_value"], 0),
        "reading_better": (was["reading"]["lower_bound_above_zero"], re2g.READS_BETTER),
        "reading_none": (was["reading"]["interval_covers_zero"], re2g.NO_DIFFERENCE),
        "reading_worse": (was["reading"]["upper_bound_below_zero"], re2g.READS_WORSE),
        "feature_unavailable": (was["feature_unavailable"]["strata"], list(re2g.FEATURE_UNAVAILABLE)),
        "benchmark_commit": (was["join_rule_source"]["commit"], re2g.BENCHMARK_COMMIT),
    }
    differs = sorted(k for k, (a, b) in checks.items() if a != b)
    numstat = subprocess.run(
        ["git", "diff", "--numstat", REGISTRATION_COMMIT, "--", "genomeos/attribution/re2g.py"],
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    ).stdout.split()
    return {
        "compared_against": f"{REGISTRATION_FILE} at {REGISTRATION_COMMIT}",
        "elements_compared": sorted(checks),
        "elements_that_differ": differs,
        "registration_unchanged": not differs,
        "module_diff_since_registration": (
            {"insertions": int(numstat[0]), "deletions": int(numstat[1])} if len(numstat) >= 2 else None
        ),
        "additions_outside_the_registered_logic": ["PREDICTION_COLUMNS", "_overlaps", "join_predictions"],
        "why_this_is_not_an_amendment": (
            "the additions implement the join rule the registration already fixed in words; no "
            "registered value, threshold, population, wording or convention was altered, which the "
            "comparison above checks rather than claims. Insertions only, no deletions"
        ),
    }


def prediction_path(cell: str) -> Path:
    return PREDICTIONS / f"{re2g.COMPARATOR_FILES[cell]['file']}.bed.gz"


def comparator_scores(pairs: list[crispri.Pair], aggregate: str) -> tuple[list[float], dict[str, Any]]:
    """rE2G score per pair: joined where a prediction overlaps, 0 where none does, per cell type."""
    scores = [0.0] * len(pairs)
    joined: set[int] = set()
    report: dict[str, Any] = {}
    for cell in sorted(re2g.COMPARATOR_FILES):
        path = prediction_path(cell)
        got, hits = re2g.join_predictions(pairs, path, cell, aggregate=aggregate)
        of_cell = [i for i, p in enumerate(pairs) if p.cell == cell]
        for i, v in got.items():
            scores[i] = v
            joined.add(i)
        report[cell] = {
            "file": re2g.COMPARATOR_FILES[cell]["file"],
            "pairs": len(of_cell),
            "pairs_joined": len(got),
            "pairs_unpredicted_scored_zero": len(of_cell) - len(got),
            "overlapping_prediction_rows": hits,
            "positives": sum(1 for i in of_cell if pairs[i].regulated),
            "positives_joined": sum(1 for i in of_cell if pairs[i].regulated and i in got),
        }
    report["_total"] = {
        "pairs": len(pairs),
        "pairs_joined": len(joined),
        "pairs_unpredicted_scored_zero": len(pairs) - len(joined),
        "aggregate_function": aggregate,
    }
    return scores, report


def our_scores(
    training: list[crispri.Pair], heldout: list[crispri.Pair]
) -> tuple[dict[str, list[float]], dict[str, Any]]:
    """The frozen models' per-pair held-out scores, by the recipe `crispri.score_published` uses.

    Weights are fitted on the covered training pairs and frozen before the held-out pairs are scored,
    which is independence statement (b). Nothing is refitted on held-out data.
    """
    table, cache = crispri.DeletionTable(), crispri.ElementCache()
    crispri.annotate(training + heldout, table, cache)
    for p in training + heldout:
        p.features["covered"] = float(p.covered)
    covered = [p for p in training if p.covered]
    weights = {
        name: crispri.logistic_fit(crispri.matrix(covered, cols), [p.regulated for p in covered])
        for name, cols in crispri.FEATURES.items()
    }
    scores = {
        name: crispri.logistic_score(weights[name], crispri.matrix(heldout, cols))
        for name, cols in crispri.FEATURES.items()
    }
    answered = [p for p in heldout if p.features.get("deletion_answered") == 1.0]
    both_zero = [p for p in heldout if p.features["deletion_drop"] == 0 and p.features["top_target"] == 0]
    provenance = {
        "fitted_on": f"{len(covered)} covered training pairs of {len(training)}",
        "frozen_before_scoring_heldout": True,
        "features": {name: list(cols) for name, cols in crispri.FEATURES.items()},
        # Named carefully, because the first version of this block called element coverage "pairs with
        # a deletion value" and so overstated the feature's reach by more than threefold. Coverage
        # means a registry element overlaps the pair; a value means the sweep actually scored this
        # gene in this pair's own cell line.
        "heldout_pairs_with_an_overlapping_registry_element": sum(1 for p in heldout if p.covered),
        "heldout_pairs_with_a_deletion_value_answered": len(answered),
        "heldout_pairs_with_a_predicted_drop_above_zero": sum(
            1 for p in heldout if p.features["deletion_drop"] > 0
        ),
        "heldout_pairs_where_both_deletion_columns_are_exactly_zero": len(both_zero),
        "what_that_means": (
            f"on {len(both_zero)} of {len(heldout)} pooled pairs the frozen "
            "'activity + distance + deletion' model is 'activity + distance' with its two deletion "
            "columns at zero, so the pooled delta is mostly a comparison of activity and distance "
            "against ENCODE-rE2G rather than of the deletion feature. The deletion feature acts on the "
            "K562 and GM12878 pairs, which is why K562 is the registered primary population"
        ),
    }
    return scores, provenance


def zero_fill_split(rows: list[crispri.Pair]) -> dict[str, Any]:
    """Three kinds of zero kept apart on one population, because a value never looked at is not a
    measured value -- the rule `crispri.UNAVAILABLE_GAIN` enforces one level up.

    `answered_nonzero` is the sweep answering with a drop; `answered_zero` is the sweep answering and
    predicting no drop; `unanswered_filled` is nothing ever answered, split into pairs with no
    overlapping registry element at all and pairs whose element was in the registry but whose gene the
    sweep was silent about. The headline's own coverage record counts the first of those two.
    """

    def count(kept: list[crispri.Pair]) -> dict[str, Any]:
        return {
            "pairs": len(kept),
            "positives": sum(1 for p in kept if p.regulated),
            "weighted_positives": round(sum(p.weight for p in kept if p.regulated), 2),
        }

    answered = [p for p in rows if p.features.get("deletion_answered") == 1.0]
    unanswered = [p for p in rows if p.features.get("deletion_answered") != 1.0]
    return {
        "answered_nonzero": count([p for p in answered if p.features["deletion_drop"] > 0]),
        "answered_zero": count([p for p in answered if p.features["deletion_drop"] == 0]),
        "unanswered_filled": count(unanswered),
        "unanswered_filled_split": {
            "no_overlapping_registry_element": count([p for p in unanswered if not p.covered]),
            "element_in_the_registry_but_the_sweep_was_silent_about_this_gene": count(
                [p for p in unanswered if p.covered]
            ),
        },
        "why_kept_apart": (
            "a zero because the sweep answered and predicted no drop, and a zero because nothing was "
            "ever answered, are different outputs and are not pooled here. The frozen model sees the "
            "same 0.0 in both cases, which is a property of the model, stated rather than hidden"
        ),
    }


def coverage_reconciliation(k562: list[crispri.Pair], headline: dict[str, Any]) -> dict[str, Any]:
    """Our K562 coverage against the headline's own record, and why 'unanswered' is a different count."""
    arms = headline["coverage_arms_k562_heldout"]  # a top-level key, not under heldout_published_pairs
    recorded = arms["arm1_all_pairs"]
    uncovered = sum(1 for p in k562 if not p.covered)
    unanswered = sum(1 for p in k562 if p.features.get("deletion_answered") != 1.0)
    return {
        "headline_record": "crispri_published.json coverage_arms_k562_heldout.arm1_all_pairs",
        "headline_uncovered_pairs": recorded["uncovered_pairs"],
        "ours_uncovered_pairs": uncovered,
        "uncovered_agrees": uncovered == recorded["uncovered_pairs"],
        "ours_unanswered_pairs": unanswered,
        "headline_arm1_deletion_gain": recorded.get("gain"),
        "these_are_two_different_counts": (
            "the headline's 174 is the COVERAGE count: pairs with no overlapping registry element. Ours "
            f"is {uncovered}, which agrees exactly. The ANSWERED count is a different question -- "
            f"whether the sweep returned a value for this pair's gene in this cell -- and gives "
            f"{unanswered} unanswered, of which {uncovered} are the uncovered ones and the rest are "
            "pairs whose element is in the registry but whose gene the sweep was silent about. The two "
            "numbers differing is not a join difference; comparing them to each other would be the error"
        ),
        "the_headlines_coverage_argument_stands_beside_this": {
            "arm2_coverage_matched_median_gain": arms["arm2_coverage_matched"]["median_gain"],
            "arm3_coverage_indicator_control_gain": arms["arm3_coverage_indicator_control"]["gain"],
            "note": (
                "the headline already argued coverage-invariance on its own terms; the split above "
                "describes where the feature fires and does not contradict it"
            ),
        },
    }


#: What our activity term is built from, and what it shares with the comparator. Facts read from this
#: project's own module and committed result, cited rather than re-derived.
SHARED_INPUTS = {
    "our_activity_term": (
        "sqrt(DNase x H3K27ac), the geometric mean of the benchmark table's own DHS.RPM and H3K27ac.RPM "
        "columns at the tested element (crispri.py, _annotate_one and the contact annotator)"
    ),
    "read_from_the_same_file_as_the_labels": (
        "crispri.py's POSITIVE_FILTER already records it: the activity columns come from the same "
        "benchmark table that carries the Regulated labels, and the held-out positives were selected "
        "partly on chromatin at the tested element"
    ),
    "dnase_is_a_shared_input": (
        "ENCODE-rE2G's published held-out model is DNase-only, so DNase enters both models. The "
        "comparison is therefore not between disjoint evidence"
    ),
    "h3k27ac_is_selection_correlated": (
        "H3K27ac enters our activity term and also entered the selection of the positives: every one of "
        "the 190 held-out positives sits in an H3K27ac element while 1,438 of the 4,188 negatives do "
        "not (docs/ROADMAP.md, 2026-09-27)"
    ),
    "the_honest_statement": (
        "the comparison advantages our model through a shared input (DNase) and a selection-correlated "
        "one (H3K27ac). This does not invalidate the registered reading, which is about ranking on "
        "these pairs, but the reading may not be quoted as evidence of a better model of enhancer-gene "
        "regulation in general, and leaving this out would misrepresent it"
    ),
    "the_existing_diagnostic": (
        "crispri.dnase_only_diagnostic and crispri_published.json's post_hoc_positive_filter block "
        "already address this; they are cited here rather than replaced. Nothing is refitted to "
        "equalise the inputs: that would be a different study and is not authorised"
    ),
}


def main() -> None:
    missing = [c for c in re2g.COMPARATOR_FILES if not prediction_path(c).exists()]
    if missing:
        raise SystemExit(f"prediction files not in the cache for: {', '.join(missing)}")

    training = crispri.load(crispri.TRAINING)
    heldout = crispri.load(crispri.HELDOUT)
    labels = [p.regulated for p in heldout]
    weights = [p.weight for p in heldout]

    print("joining the five prediction files by the registered rule (aggregate=sum) ...")
    theirs, join_report = comparator_scores(heldout, "sum")
    print(
        f"  joined {join_report['_total']['pairs_joined']} of {len(heldout)} pairs; "
        f"{join_report['_total']['pairs_unpredicted_scored_zero']} scored 0"
    )

    reproduced = crispri.benchmark_auprc(theirs, labels, weights)
    verdict = re2g.gate(reproduced)
    print(
        f"GATE: reproduced {verdict['reproduced']} against {re2g.GATE_TARGET} "
        f"(tolerance {re2g.GATE_TOLERANCE}) -> {'PASS' if verdict['passed'] else 'STOP'}"
    )

    payload: dict[str, Any] = {
        # First key on purpose. The positive figure below is withdrawn by the like-for-like comparison,
        # and no reader should meet the positive before the thing that withdraws it.
        "read_this_first": (
            "THE BINDING RESULT OF THIS LANE IS NOT IN THIS FILE. The like-for-like comparison "
            "(data/results/re2g_likeforlike.json, registered in advance as re2g.SECOND_REGISTRATION) "
            "matches our feature set to the comparator's input and reads NO DIFFERENCE DETECTED on the "
            "K562 primary; that file carries the figure and its interval, which are deliberately not "
            "copied here so the two cannot drift apart. By its own pre-registered falsifier, "
            "README MAY NOT SAY THE DELETION MODEL RANKS BETTER THAN ENCODE-rE2G, and the comparison in "
            "this file MAY NOT BE CITED ALONE. The reason is that this file's activity term reads "
            "H3K27ac, which all 190 held-out positives carry and which the comparator's published "
            "held-out model does not read, so the benchmark's own positive selection favours our side "
            "here. Everything below is reported as registered and is true of these pairs; it is not "
            "evidence that this project's model is the better model of enhancer-gene regulation"
        ),
        "lane": "lane-re2g",
        "registration": "data/results/re2g_registration.json",
        "claim": re2g.CLAIM,
        "endpoint": re2g.ENDPOINT,
        "join_rule": dict(re2g.JOIN_RULE),
        "join_rule_commit": {
            "repository": "EngreitzLab/CRISPR_comparison",
            "tag": re2g.BENCHMARK_TAG,
            "commit": re2g.BENCHMARK_COMMIT,
        },
        "comparator_files": {
            cell: {
                **dict(v),
                "url": PORTAL.format(a=v["file"]),
                **dict(zip(("sha256", "bytes_on_disk"), _hash(prediction_path(cell)), strict=True)),
            }
            for cell, v in re2g.COMPARATOR_FILES.items()
        },
        "join": join_report,
        "gate": verdict,
        "independence_and_exposure": dict(re2g.INDEPENDENCE),
        "falsifier": re2g.FALSIFIER,
        "registration_unchanged": registration_unchanged(),
        "code_cleanliness": code_cleanliness(),
        "alphagenome_requests": 0,
    }

    if not verdict["passed"]:
        payload["status"] = (
            "STOPPED at the gate: the comparator scores do not reproduce the published held-out figure "
            "within the registered tolerance, so no comparison is made"
        )
        payload["paired_delta"] = None
        payload["reading"] = None
        _write(payload, heldout)
        print("the gate failed; no delta computed. " + verdict["reason"])
        raise SystemExit(1)

    print("gate passed; scoring our own frozen model (0 model requests) ...")
    ours_all, ours_provenance = our_scores(training, heldout)
    ours = ours_all["activity + distance + deletion"]

    populations = {
        "primary_k562": [i for i, p in enumerate(heldout) if p.cell == "K562"],
        "secondary_pooled": list(range(len(heldout))),
        "gm12878": [i for i, p in enumerate(heldout) if p.cell == "GM12878"],
    }
    deltas: dict[str, Any] = {}
    for name, idx in populations.items():
        rows = [heldout[i] for i in idx]
        out = re2g.paired_delta([ours[i] for i in idx], [theirs[i] for i in idx], rows)
        out["population"] = {
            "name": name,
            "pairs": len(rows),
            "positives": sum(1 for p in rows if p.regulated),
            "weighted_positives": round(sum(p.weight for p in rows if p.regulated), 2),
        }
        out["ours_auprc"] = crispri.bench_metrics([ours[i] for i in idx], rows, weighted=True)["auprc"]
        out["re2g_auprc"] = crispri.bench_metrics([theirs[i] for i in idx], rows, weighted=True)["auprc"]
        out["deletion_reach"] = {
            "pairs_with_a_deletion_value_answered": sum(
                1 for p in rows if p.features.get("deletion_answered") == 1.0
            ),
            "pairs_with_a_predicted_drop_above_zero": sum(1 for p in rows if p.features["deletion_drop"] > 0),
            "pairs_where_both_deletion_columns_are_zero": sum(
                1 for p in rows if p.features["deletion_drop"] == 0 and p.features["top_target"] == 0
            ),
            "of_pairs": len(rows),
        }
        deltas[name] = out
        print(
            f"  {name}: delta {out['delta_auprc']} ci {out['ci95']} "
            f"clusters {out['clusters']} draws {out['resamples']}/{out['draws_requested']} "
            f"-> {out['reading']}"
        )

    for cell in re2g.FEATURE_UNAVAILABLE:
        deltas[cell.lower()] = re2g.delta_where_available(cell, lambda: {})
        print(f"  {cell}: feature unavailable, no number reported")

    payload["our_model"] = ours_provenance
    payload["paired_delta"] = deltas
    payload["reading"] = {
        name: deltas[name]["reading"] for name in populations if deltas[name].get("reading")
    }
    payload["reading_must_travel_with"] = {
        "rule": (
            "the registered wording is used word for word and NEVER alone. Each of the three "
            "qualifiers below is part of the same sentence wherever the reading is quoted, including in "
            "README, docs/ROADMAP.md, docs/ATTRIBUTION.md and any reply to the owner"
        ),
        "lower_bound": (
            f"the K562 primary interval's lower bound is {deltas['primary_k562']['ci95'][0]}, unrounded"
        ),
        "the_comparator_is_a_reconstruction": (
            f"the comparator is ENCODE-rE2G as reconstructed from the ENCODE portal prediction files, "
            f"which reproduce its published pooled weighted AUPRC {verdict['miss']} low "
            f"({verdict['reproduced']} against {re2g.GATE_TARGET:.4f}). The primary's lower bound is "
            "about an eighth of that shortfall. No adjustment was made to the comparator's scores"
        ),
        "the_positive_selection_favours_an_h3k27ac_model": (
            "all 190 held-out positives carry H3K27ac at the tested element (52 'H3K27ac', 138 'High "
            "H3K27ac', and none in 'No H3K27ac', 'CTCF element' or 'H3K27me3 element'), while 1,438 of "
            "the non-regulated pairs lack it (589 No H3K27ac, 552 CTCF, 297 H3K27me3). Our activity "
            "term reads H3K27ac; ENCODE-rE2G's published held-out model reads DNase only. The "
            "benchmark's own positive selection therefore hands an H3K27ac-reading model a separation "
            "that the comparator cannot use, and the primary's lower bound sits well inside it. Counts "
            "quoted from crispri_published.json post_hoc_positive_filter"
        ),
        "what_may_not_be_said": (
            "this reading may not be cited alone as evidence about ENCODE-rE2G. A like-for-like "
            "comparison that matches our feature set to the comparator's inputs is registered "
            "separately in genomeos/attribution/re2g.py (SECOND_REGISTRATION) and is the one README may "
            "cite about ENCODE-rE2G"
        ),
    }
    # Descriptive decomposition: the two-feature model against the same comparator on the same pairs,
    # so the deletion feature's share of each delta is visible as the difference between the two.
    two = ours_all["activity + distance"]
    decomposition: dict[str, Any] = {}
    for name, idx in populations.items():
        rows = [heldout[i] for i in idx]
        out = re2g.paired_delta([two[i] for i in idx], [theirs[i] for i in idx], rows)
        decomposition[name] = {
            "activity_plus_distance_auprc": crispri.bench_metrics([two[i] for i in idx], rows, weighted=True)[
                "auprc"
            ],
            "re2g_auprc": deltas[name]["re2g_auprc"],
            "delta_auprc": out["delta_auprc"],
            "ci95": out["ci95"],
            "reading": out["reading"],
            "full_model_delta": deltas[name]["delta_auprc"],
            "deletion_features_share": round(deltas[name]["delta_auprc"] - (out["delta_auprc"] or 0), 4),
        }
        print(f"  [decomp] {name}: activity+distance vs rE2G {out['delta_auprc']} ci {out['ci95']}")
    headline_arm1 = json.loads(HEADLINE.read_text())["coverage_arms_k562_heldout"]["arm1_all_pairs"]
    decomposition["what_this_says"] = (
        "the two-feature 'activity + distance' model, scored against the same reconstructed "
        "ENCODE-rE2G on the identical pairs. Descriptive, not a registered endpoint. Asked plainly: "
        "does a two-feature activity-and-distance model rank at or above ENCODE-rE2G as reconstructed? "
        "On the K562 primary, NO -- the point estimate is below zero. On the pooled population, NOT "
        "DETECTABLY -- the interval covers zero. So the advantage the registered reading reports is "
        "carried by the deletion features on both populations, which is the opposite of what the count "
        "of zero-valued columns suggests and is why it was measured. Every figure here still inherits "
        "the shared-input qualifications in shared_inputs"
    )
    decomposition["a_count_of_zeroes_does_not_measure_a_contribution"] = (
        "Both the coordinator and this lane first argued from the 3,528 pooled pairs whose deletion "
        "columns are zero to the conclusion that the pooled delta was largely not a test of the deletion "
        "feature. The inference does not hold, and measuring it gives the opposite answer: with the "
        "deletion columns removed, the remaining model does not clear zero against the same comparator "
        "on either population. A column that is zero on most pairs can still carry the ranking, because "
        "what it does on the pairs where it fires is to lift them past the rest. The lesson kept here is "
        "that the share of pairs a feature is non-zero on is not a measurement of what the feature "
        "contributes, and only the comparison with the feature removed is"
    )
    decomposition["cross_check_against_the_headlines_own_deletion_gain"] = {
        "ours_k562_deletion_share": decomposition["primary_k562"]["deletion_features_share"],
        "headline_arm1_gain_k562_heldout": headline_arm1["gain"],
        "headline_arm1_ci95": headline_arm1["ci95"],
        "note": (
            "the share is the full model's delta minus the two-feature model's against a common "
            "comparator, so the comparator cancels and it is the deletion gain on the same pairs. It is "
            "compared here against the committed headline's own K562 held-out deletion gain, derived "
            "independently and at 200 draws, as a check that this lane scored the frozen model"
        ),
    }
    payload["decomposition_descriptive"] = decomposition

    k_two = decomposition["primary_k562"]["delta_auprc"]
    p_two = decomposition["secondary_pooled"]["delta_auprc"]
    payload["status"] = (
        "gate passed; paired comparison of two frozen models on a reused benchmark. HEADLINE, measured "
        "rather than argued from a count: on 3,528 of the 4,378 pooled pairs both deletion columns are "
        "zero, which invites the conclusion that the pooled delta is not a test of the deletion "
        "feature. THE DECOMPOSITION DOES NOT SUPPORT THAT CONCLUSION. Scored against the same "
        f"reconstructed ENCODE-rE2G on the identical pairs, 'activity + distance' alone gives {k_two} "
        f"on K562 and {p_two} pooled, and neither interval clears zero, while the full model gives "
        f"{deltas['primary_k562']['delta_auprc']} and {deltas['secondary_pooled']['delta_auprc']}. So "
        "the deletion features carry the advantage on BOTH populations, and the count of zero columns "
        "does not measure the feature's contribution. What the zero columns do mean is that the pooled "
        "advantage comes from lifting the 850 pairs where the feature is non-zero relative to the other "
        "strata, not from ranking within the strata that have no value"
    )

    # Post hoc, labelled: the primary delta on the pairs the sweep actually answered.
    answered_idx = [
        i for i in populations["primary_k562"] if heldout[i].features.get("deletion_answered") == 1.0
    ]
    answered_rows = [heldout[i] for i in answered_idx]
    post_hoc = re2g.paired_delta(
        [ours[i] for i in answered_idx], [theirs[i] for i in answered_idx], answered_rows
    )
    payload["post_hoc_answered_pairs_only"] = {
        "label": "POST HOC. Not the registered estimate and does not replace it",
        "population": "the held-out K562 pairs the sweep answered",
        "pairs": len(answered_rows),
        "positives": sum(1 for p in answered_rows if p.regulated),
        "weighted_positives": round(sum(p.weight for p in answered_rows if p.regulated), 2),
        "delta_auprc": post_hoc["delta_auprc"],
        "ci95": post_hoc["ci95"],
        "clusters": post_hoc["clusters"],
        "resamples": post_hoc["resamples"],
        "reading_if_this_were_registered": post_hoc["reading"],
        "agrees_with_the_registered_primary": (post_hoc["reading"] == deltas["primary_k562"]["reading"]),
        "note": (
            "the registered primary is all 1,918 K562 pairs with the benchmark's fill rule. If this "
            "sensitivity and the registered estimate pointed in different directions the reading would "
            "have to carry both; whether they do is recorded in agrees_with_the_registered_primary"
        ),
    }
    payload["zero_fill_split_primary_k562"] = zero_fill_split(
        [heldout[i] for i in populations["primary_k562"]]
    )
    payload["coverage_reconciliation_primary_k562"] = coverage_reconciliation(
        [heldout[i] for i in populations["primary_k562"]], json.loads(HEADLINE.read_text())
    )
    payload["shared_inputs"] = dict(SHARED_INPUTS)
    payload["sensitivities"] = _sensitivities(heldout, ours, labels, weights)
    payload["stated_beside_the_numbers"] = {
        "the_comparator_reproduction_falls_short_and_that_favours_us": (
            f"the gate reproduced {verdict['reproduced']} against the published "
            f"{re2g.GATE_TARGET:.4f}, a shortfall of {verdict['miss']}. The paired delta is measured "
            "against this reproduction, not against the published figure, so the shortfall inflates "
            "the pooled delta by roughly that much: our pooled 0.677 stands "
            f"{round(0.677 - re2g.GATE_TARGET, 4)} above the published figure unpaired, against "
            f"{deltas['secondary_pooled']['delta_auprc']} paired against the reproduction. The paired "
            "interval is the registered statistic; this is the size of the gap it does not account for"
        ),
        "the_activity_term_is_not_the_comparator_s": (
            "already on this project's record (docs/ROADMAP.md, 2026-09-27): the held-out comparison "
            "flatters this project. Every one of the 190 held-out positives sits in an H3K27ac element "
            "while 1,438 of the 4,188 negatives do not, and the activity term here reads H3K27ac where "
            "the published held-out ENCODE-rE2G reads DNase only. This is not one of the four "
            "registered independence statements; it is a fifth exposure and it is not quantified here"
        ),
        "what_the_deletion_feature_is_on_a_non_model_cell_types_pair": (
            "nothing, and this was checked in code rather than assumed. `crispri.annotate` computes the "
            "two deletion columns only when the pair's own cell line is one the sweep scored "
            "(MODEL_CELLS: K562, HepG2, GM12878, IMR-90); for any other cell it takes the branch "
            "`(0.0, [])`. Counted on the held-out pairs: HCT116 0 of 396, Jurkat 0 of 75 and WTC11 0 of "
            "1,921 carry a deletion value answered, a predicted drop above zero, or a top-target flag. "
            "So NO feature derived from other cell types is applied to those strata: both columns are "
            "exactly zero there. What is non-zero on 1,616 WTC11, 363 HCT116 and 64 Jurkat pairs is "
            "`covered`, which only records that a registry element overlaps the tested element and is "
            "not a model column. An earlier draft of this result mislabelled that coverage count as "
            "'pairs with a deletion value' and so overstated the feature's reach more than threefold"
        ),
        "what_the_pooled_figure_therefore_rests_on": (
            "on 3,528 of the 4,378 pooled pairs both deletion columns are exactly zero, so the frozen "
            "'activity + distance + deletion' model is 'activity + distance' plus two constant-zero "
            "terms there. Within each of those strata the zero columns cannot change the ranking at "
            "all; across the pooled set they shift those pairs by a constant relative to the K562 and "
            "GM12878 pairs that do carry a value. The pooled delta is therefore mostly a comparison of "
            "activity and distance against ENCODE-rE2G, and it is a property of this comparison, NOT a "
            "cross-cell-type validation of the deletion feature. The registered primary is K562, where "
            "the feature is that cell line's own"
        ),
        "unavailable_refers_to_the_per_cell_gain": (
            "the per-stratum 'feature unavailable' wording for HCT116, Jurkat and WTC11 is about the "
            "per-cell deletion gain: there is no value in those cells, so a difference between the two "
            "model forms there would be an artefact of the form. On the evidence above it is also true "
            "of the model's input in those strata, which are exactly zero. Neither statement is "
            "weakened: no number is reported for those three strata in either sense"
        ),
        "the_gm12878_interval_rests_on_seven_clusters": (
            f"GM12878 has {deltas['gm12878']['clusters']} chromosome clusters over 68 pairs and 16 "
            "positives, so its interval is as coarse as that allows and its width should be read as "
            "such; clusters, requested and dropped draws are reported with it"
        ),
        "the_primary_lower_bound_is_close_to_zero": (
            f"the K562 interval's lower bound is {deltas['primary_k562']['ci95'][0]}, not rounded. It "
            "clears zero, and only just. The registered reading holds and is not softened here, and it "
            "is not to travel without this bound beside it"
        ),
    }
    _write(payload, heldout)


def _hash(path: Path) -> tuple[str, int]:
    sha, size, _ = mf.sha256_of(path)
    return sha, size


def _sensitivities(
    heldout: list[crispri.Pair], ours: list[float], labels: list[bool], weights: list[float]
) -> dict[str, Any]:
    """The two sensitivities registered in advance, reported beside the primary and never in place of it."""
    out: dict[str, Any] = {}
    by_max, _ = comparator_scores(heldout, "max")
    out["aggregate_max"] = {
        "what": "several overlapping predictions aggregated by max instead of the registered sum",
        "re2g_pooled_auprc": round(crispri.benchmark_auprc(by_max, labels, weights) or 0, 4),
        "gate_on_this_variant": re2g.gate(crispri.benchmark_auprc(by_max, labels, weights)),
        "note": "pre-declared; the gate is judged on sum alone, so this is information, not a verdict",
    }
    joined = [i for i, v in enumerate(by_max) if v != 0.0]
    kept = [heldout[i] for i in joined]
    if any(p.regulated for p in kept) and not all(p.regulated for p in kept):
        drop = re2g.paired_delta([ours[i] for i in joined], [by_max[i] for i in joined], kept)
        out["drop_unpredicted_pairs"] = {
            "what": "the unpredicted pairs dropped instead of scored 0, on the pairs that joined",
            "pairs": len(kept),
            "positives": sum(1 for p in kept if p.regulated),
            "delta_auprc": drop["delta_auprc"],
            "ci95": drop["ci95"],
            "reading": drop["reading"],
            "note": "pre-declared sensitivity; the primary keeps the benchmark's own fill rule of 0",
        }
    return out


def _write(payload: dict[str, Any], heldout: list[crispri.Pair]) -> None:
    inputs = [
        mf.input_entry(
            crispri.KNOWLEDGE / crispri.HELDOUT,
            partition="held-out",
            url=crispri.BASE_URL + crispri.HELDOUT,
        ),
        mf.input_entry(
            crispri.KNOWLEDGE / crispri.TRAINING,
            partition="training",
            url=crispri.BASE_URL + crispri.TRAINING,
        ),
    ] + [
        mf.input_entry(
            prediction_path(cell),
            partition=cell,
            url=PORTAL.format(a=re2g.COMPARATOR_FILES[cell]["file"]),
            accession=re2g.COMPARATOR_FILES[cell]["file"],
            annotation=re2g.COMPARATOR_FILES[cell]["annotation"],
        )
        for cell in sorted(re2g.COMPARATOR_FILES)
    ]
    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "Gschwind et al. 2026, Nature, doi:10.1038/s41586-026-10781-4, "
                "Supplementary Table 3 (the gate figure) and Supplementary Table 12 (the accessions)",
                "version": "version of record",
            },
            {
                "accession": "EngreitzLab/CRISPR_comparison held-out benchmark and merge rule",
                "version": f"{re2g.BENCHMARK_TAG} ({re2g.BENCHMARK_COMMIT})",
            },
            {
                "accession": "ENCODE ENCODE-rE2G element gene links, one per benchmark biosample "
                + ", ".join(f"{c}:{v['file']}" for c, v in sorted(re2g.COMPARATOR_FILES.items())),
                "version": "released 2024-06-06, software distal-regulation-encode_re2g",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "gate_target": re2g.GATE_TARGET,
            "gate_tolerance": re2g.GATE_TOLERANCE,
            "draws_requested": re2g.DRAWS,
            "seed": re2g.SEED,
            "cluster": re2g.CLUSTER,
            "aggregate_function": "sum",
            "fill_value": 0,
            "overlap": "the pipeline's inclusive findOverlaps on both sides' raw BED numbers",
        },
        "exclusions": [
            "no pair is excluded: the benchmark's published pair set and labels are used whole",
            f"no number is reported for {', '.join(re2g.FEATURE_UNAVAILABLE)}: the deletion feature "
            "does not exist there",
        ],
        "partitions": {
            "primary_k562": f"{sum(1 for p in heldout if p.cell == 'K562')} held-out K562 pairs",
            "secondary_pooled": f"all {len(heldout)} held-out pairs",
            "gm12878": f"{sum(1 for p in heldout if p.cell == 'GM12878')} held-out GM12878 pairs",
        },
    }
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(json.dumps(payload["gate"], indent=1))


if __name__ == "__main__":
    main()
