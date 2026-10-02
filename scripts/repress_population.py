# SPDX-License-Identifier: AGPL-3.0-or-later
"""Gate 1 on the repression population, then the registered prediction and the internal check
(data/results/repress2_population.json).

    uv run --frozen python scripts/repress_population.py
    uv run --frozen python scripts/repress_population.py --chroms chr21 --result repress2_population_chr21

In the order `data/results/repress2_registration.json` fixed, which was committed before this script
ran:

1. **Gate 1, blind.** How many measured pairs with a significant signed effect attach to a
   `represses_target` rule, and to an `activates_target` rule, and over how many independent loci.
   The predicate is `repress2.significant`, which reads whether the pair's outcome is one of the two
   labels a significant pair takes and never which of the two, so no direction enters the count.
   Both floors are `repress2`'s imports. Below either floor the count goes on the record as a no-go.
2. **The registered prediction**, tested descriptively: the repression share by effect band and
   reader state, with each registered comparison reported by name as a strict rise or not.
3. **The 0-request internal check**: sign concordance between the compiled call and the element's own
   `predicted_coding_by_cell` across the four retained cell lines, for repress and activate rules
   separately, with the calls that disagree in every retained cell counted as outliers of that reading
   rather than as a mechanism.

Nothing is refitted, no threshold is moved, no verdict is moved, no compiled label is changed and
nothing is deleted. No model is called and nothing is downloaded. The rules are
`genomeos.attribution.not_open_profile.rules` unchanged and the state per rule is
`genomeos.attribution.context_evidence.state_for` unchanged, so this run and the profile it builds on
cannot describe different rules.
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter, defaultdict
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import compile as cp  # noqa: E402
from genomeos.attribution import context_evidence as ce  # noqa: E402
from genomeos.attribution import crispri  # noqa: E402
from genomeos.attribution import measured as ms  # noqa: E402
from genomeos.attribution import not_open_profile as nop  # noqa: E402
from genomeos.attribution import repress2 as rp  # noqa: E402
from genomeos.attribution.targets import attributed  # noqa: E402
from genomeos.genome import reader  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

RESULT = "repress2_population"
REGISTRATION = "repress2_registration"
PROFILE = "not_open_profile"
CHROMS = tuple(f"chr{c}" for c in [*range(1, 23), "X", "Y"])

#: The profile figures this run must reproduce, from data/results/not_open_profile.json, so that a
#: discrepancy in the rule enumeration shows up here rather than being carried into a share.
ASSESSABLE_ON_THE_RECORD = 81_635
REPRESS_ASSESSABLE_ON_THE_RECORD = 30_480
ACTIVATE_ASSESSABLE_ON_THE_RECORD = 51_103

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/repress2.py",
    "scripts/repress_register.py",
    "scripts/repress_population.py",
    "tests/test_repress2.py",
)

#: The same per-chromosome results `scripts/not_open_profile.py` reads, each digested by its own path.
INPUT_GLOBS = (
    "budget_chr*.json",
    "budget_axes_chr*.json",
    "variation_chr*.json",
    "duplication_chr*.json",
    "domains_chr*.json",
    "unknown_chr*.json",
    "reader_*_chr*.json",
    "ccres_chr*.bed.gz",
    "enhancer_targets_chr*.json",
    "enhancer_targets_all_chr*.json",
    "constrained_targets_chr*.json",
)


def reproduces(
    per_state: Counter[str], per_axis_state: dict[str, Counter[str]], whole_genome: bool
) -> dict[str, Any]:
    """Whether this run's rule enumeration gives the profile's own figures, where it covers the same
    population. A run over fewer chromosomes says so rather than comparing against a genome-wide one."""
    if not whole_genome:
        return {
            "checked": False,
            "why": "the profile's figures are over all 24 chromosomes and this run covers fewer",
        }
    assessable = per_state[ce.STATE_OPEN] + per_state[ce.STATE_NOT_OPEN]
    repress = sum(per_axis_state[rp.REPRESSES].values())
    activate = sum(per_axis_state[rp.ACTIVATES].values())
    return {
        "checked": True,
        "assessable_rules": assessable,
        "assessable_rules_on_the_record": ASSESSABLE_ON_THE_RECORD,
        "represses_target_assessable": repress,
        "represses_target_assessable_on_the_record": REPRESS_ASSESSABLE_ON_THE_RECORD,
        "activates_target_assessable": activate,
        "activates_target_assessable_on_the_record": ACTIVATE_ASSESSABLE_ON_THE_RECORD,
        "agrees": (
            assessable == ASSESSABLE_ON_THE_RECORD
            and repress == REPRESS_ASSESSABLE_ON_THE_RECORD
            and activate == ACTIVATE_ASSESSABLE_ON_THE_RECORD
        ),
        "source": f"data/results/{PROFILE}.json",
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chroms", default=",".join(CHROMS))
    ap.add_argument("--result", default=RESULT, help="a run over fewer chromosomes writes its own name")
    args = ap.parse_args()
    chroms = [c for c in args.chroms.split(",") if c]
    if args.result == RESULT and tuple(chroms) != CHROMS:
        raise SystemExit(f"{RESULT} is the genome-wide reading; a run over {len(chroms)} needs --result")
    whole_genome = tuple(chroms) == CHROMS

    started = time.time()
    table = ce.mapping()
    readers = ce.Readers(RESULTS_DIR)

    # ---- gate 1, blind: nothing below reads a direction -----------------------------------------
    keys: dict[tuple[str, bool], list] = defaultdict(list)
    ladders: dict[str, Counter[str]] = defaultdict(Counter)
    # ---- the descriptive counts, (effect band, reader state, activity axis) ----------------------
    cells: Counter[tuple[str, str, str]] = Counter()
    per_state: Counter[str] = Counter()
    per_axis_state: dict[str, Counter[str]] = defaultdict(Counter)
    other_axis: Counter[str] = Counter()
    measured_source_axis: Counter[str] = Counter()
    # ---- the internal check ----------------------------------------------------------------------
    signs: dict[str, list[dict[str, str]]] = defaultdict(list)
    rules_seen = 0

    for chrom in chroms:
        elements = attributed(chrom, RESULTS_DIR)
        by_cell_of = {e["id"]: e.get(rp.BY_CELL_FIELD) for e in elements}
        call_of = {e["id"]: float(e["predicted_coding"]["log2_fold_change"]) for e in elements}
        _layer, measured_rows = cp._measured_rows(chrom, elements, RESULTS_DIR, None)

        rules = nop.rules(chrom, RESULTS_DIR)
        rules_seen += len(rules)
        for activity in rp.ACTIVITIES:
            for match_cell in (True, False):
                keys[(activity, match_cell)].extend(
                    rp.eligible_links(rules, measured_rows, activity, match_cell)
                )
            ladders[activity].update(rp.ladder(rules, measured_rows, activity))

        for r in rules:
            state, _rest = ce.parse_value(ce.state_for(r.cell, chrom, r.start, r.end, readers, table))
            per_state[state] += 1
            if state == ce.STATE_NOT_ASSESSABLE:
                continue
            per_axis_state[r.activity_axis][state] += 1
            if r.source == nop.SOURCE_MEASURED:
                measured_source_axis[r.activity_axis] += 1
                continue  # a measured rule carries no AlphaGenome effect, so it has no band
            if r.activity_axis not in rp.ACTIVITIES:
                other_axis[r.activity_axis] += 1
                continue
            cells[(r.effect_band or nop.WEAK, state, r.activity_axis)] += 1
            signs[r.activity_axis].append(rp.cell_signs(call_of[r.element], by_cell_of.get(r.element)))

    # ---- gate 1 ---------------------------------------------------------------------------------
    gates: dict[str, Any] = {}
    for activity in rp.ACTIVITIES:
        for match_cell in (True, False):
            name = f"{activity}|{'cell_matched' if match_cell else 'cell_ignored'}"
            population = (
                "measured CRISPRi pairs with a significant signed effect on the same element and gene "
                f"as a rule whose activity axis is exactly {activity!r}"
                + (
                    ", and in the cell that rule's `when: cell_type` names"
                    if match_cell
                    else ", in any cell the screen measured; reported beside the cell-matched count and "
                    "never pooled with it"
                )
            )
            gates[name] = rp.gate(keys[(activity, match_cell)], population)

    passed = [k for k, v in gates.items() if v["meets_both_floors"]]
    repress_matched = f"{rp.REPRESSES}|cell_matched"
    gate_1: dict[str, Any] = {
        "order": "taken before any direction or sign of an outcome was read",
        "call": rp.GATE_CALL,
        "blinded_fields": list(rp.GATE_FIELDS),
        "floors": rp.FLOORS_CALL,
        "per_population": gates,
        "populations_at_or_above_both_floors": passed,
        "reading": rp.GATE_PASS if repress_matched in passed else rp.GATE_NO_GO,
        "what_follows": (
            "the repression population carries a direction reading"
            if repress_matched in passed
            else "no direction is read on the repression population and no rate is reported for it; "
            "the lane continues only as the descriptive part below"
        ),
        "where_a_zero_falls": {
            "call": rp.LADDER_CALL,
            "steps": list(rp.LADDER_STEPS),
            "per_activity_axis": {
                a: {step: ladders[a].get(step, 0) for step in rp.LADDER_STEPS} for a in rp.ACTIVITIES
            },
        },
        "for_scale": (
            "the record says 212 measured rule links exist in all (element x gene x cell, "
            "data/results/response_map_coverage.json), so a repression population too small to decide "
            "is an expected outcome of this count and not a fault in it"
        ),
    }

    # ---- the registered prediction ---------------------------------------------------------------
    cells_table = rp.table(cells)
    compared = rp.comparisons(cells_table)
    reading = rp.prediction_reading(compared)
    prediction: dict[str, Any] = {
        "mechanism": rp.MECHANISM,
        "prediction": rp.PREDICTION,
        "refutation_as_registered": rp.REFUTATION,
        "registered_in": f"data/results/{REGISTRATION}.json, committed before this run",
        "population": (
            "the assessable rules of the predicted layer whose activity axis is exactly "
            f"{rp.REPRESSES!r} or exactly {rp.ACTIVATES!r}: measured-layer rules carry the screen's "
            "effect size and no AlphaGenome effect, so they have no band and are counted apart"
        ),
        "by_effect_band_and_reader_state": cells_table,
        "comparisons": compared,
        "verdict": reading,
        "rules_outside_this_population": {
            "measured_layer_rules_by_activity_axis": dict(measured_source_axis.most_common()),
            "predicted_layer_rules_on_a_combined_activity_axis": dict(other_axis.most_common()),
            "not_assessable": per_state[ce.STATE_NOT_ASSESSABLE],
            "why": rp.ACTIVITY_CALL,
        },
    }

    # ---- the internal check ----------------------------------------------------------------------
    internal = {
        "call": rp.CONCORDANCE_CALL,
        "retained_cells": list(rp.RETAINED_CELLS),
        "outlier": rp.OUTLIER_CALL,
        "not_validation": rp.CONCORDANCE_IS_NOT_VALIDATION,
        "requests": 0,
        "per_activity_axis": {a: rp.concordance(signs[a]) for a in rp.ACTIVITIES},
    }

    payload: dict[str, Any] = {
        "result": args.result,
        "date": date.today().isoformat(),
        "lane": "lane-repress2",
        "registered_before_this_run": f"data/results/{REGISTRATION}.json",
        "describes": f"data/results/{PROFILE}.json",
        "chromosomes": chroms,
        "rules_enumerated": rules_seen,
        "per_state": dict(per_state),
        "assessable_rules": per_state[ce.STATE_OPEN] + per_state[ce.STATE_NOT_OPEN],
        "reproduces_the_profile": reproduces(per_state, per_axis_state, whole_genome),
        "by_activity_axis_and_state": {
            a: dict(c) for a, c in sorted(per_axis_state.items(), key=lambda kv: -sum(kv[1].values()))
        },
        "gate_1": gate_1,
        "registered_prediction": prediction,
        "internal_check": internal,
        # nested, never spread: `registration()` carries its own `gate_1` and `internal_check`
        # descriptions, and spreading it over this payload replaced the counts above with them.
        "registration_as_committed": rp.registration(),
        "status": "the result of the registered run; the registration itself is the other file",
        "reading": {
            "moves": (
                "the repression population is counted against both imported floors before any "
                "direction is read, and the registered mechanism is tested descriptively against the "
                "refutation condition committed with it"
            ),
            "does_not_move": (
                "no verdict moved, no rule was deleted, no compiled label changed, no threshold or "
                "registered definition changed, and nothing was refitted. No rule is established to "
                "be wrong by anything here"
            ),
            "never": (
                "not_open_in_reader is never reported as the element being closed or the rule being "
                "contradicted, open_in_reader is never reported as validation, and the sign "
                "concordance across the four retained cells is never reported as validation either"
            ),
        },
        "code_cleanliness": ce.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
        "money": "none: every input was already on disk",
        "seconds": round(time.time() - started, 1),
    }

    # Every input by its own path. `manifest.files_entry` records a group under a label and
    # `scripts/manifest_rebuild.py` resolves an input by its `path`, so a grouped input is reported
    # absent without any file being checked; one entry per file keeps the rebuild able to compare.
    paths = {ce.TRACK_METADATA, Path(f"data/results/{REGISTRATION}.json")}
    for name in ms.CRISPRI_FILES:
        p = crispri.KNOWLEDGE / name
        if p.exists():
            paths.add(p)
    for cell in ce.READER_TERMS:
        for chrom in chroms:
            p = RESULTS_DIR / reader.peaks_path(cell, chrom).name
            if p.exists():
                paths.add(p)
    for glob in INPUT_GLOBS:
        for chrom in chroms:
            paths.update(RESULTS_DIR.glob(glob.replace("chr*", chrom)))
    inputs = [mf.input_entry(p, partition=None) for p in sorted(paths, key=lambda q: q.as_posix())]

    payload["result_manifest"] = {
        "sources": [
            {
                "accession": "ENCODE DNase-seq narrowPeak, GRCh38, released, 13 biosamples",
                "version": "as reader v1 cached them under data/results/dnase_*_chr*.bed.gz",
            },
            {"accession": "AlphaGenome output track metadata", "version": str(ce.TRACK_METADATA)},
            {
                "accession": "ENCODE CRISPRi enhancer-gene benchmark "
                "(EngreitzLab/CRISPR_comparison, Gschwind et al.)",
                "version": "the two benchmark tables cached under data/knowledge, digested in inputs",
            },
            {
                "accession": "the compiled non-coding programs of the 24 chromosomes",
                "version": "re-enumerated in this run from the results on disk, not read from a file",
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "openness_call": ce.OPENNESS_CALL,
            "measured_span": ce.SPAN_CALL,
            "effect_band_call": nop.EFFECT_BAND_CALL,
            "gate_call": rp.GATE_CALL,
            "independent_locus_rule": rp.LOCUS_RULE,
            "independent_locus_span": rp.LOCUS_SPAN,
            "locus_floor": rp.LOCUS_FLOOR,
            "link_floor": rp.POSITIVE_FLOOR,
            "retained_cells": list(rp.RETAINED_CELLS),
            "chromosomes": chroms,
            "aggregate_function": "count",
            "fill_value": 0,
        },
        "exclusions": [
            "no rule is excluded: every rule the compiler emits is enumerated and each rule outside a "
            "reported population is counted in its own row under rules_outside_this_population",
            "a rule whose activity axis holds more than one term is in neither the repression nor the "
            "activation population",
        ],
        "partitions": {
            "gate_1_cell_matched": "the measured pair's cell equals the rule's own `when: cell_type`",
            "gate_1_cell_ignored": "element and gene only, never pooled with the cell-matched count",
            "descriptive": "the assessable rules of the predicted layer, by effect band and reader state",
        },
    }

    path = save_result(args.result, payload)
    print(f"{args.result}: {path}")
    for name, g in gates.items():
        print(
            f"  gate 1 {name}: {g['eligible_measured_links']} links, "
            f"{g['independent_loci']} independent loci -> "
            f"{'at or above both floors' if g['meets_both_floors'] else 'no-go'}"
        )
    for a in rp.ACTIVITIES:
        print(f"  ladder {a}: " + ", ".join(f"{s}={ladders[a].get(s, 0)}" for s in rp.LADDER_STEPS))
    for key, cell in cells_table.items():
        print(f"  {key}: repression share {cell['repression_share']} of {cell['assessable_rules']}")
    print(
        f"  prediction: {'supported' if reading['supported'] else 'refuted'}; "
        f"failed {reading['which_failed']}"
    )
    for activity, c in internal["per_activity_axis"].items():
        print(
            f"  {activity}: disagrees in all four retained cells "
            f"{c['disagrees_in_all_four']} of {c['rules_scored_in_all_four_retained_cells']}"
        )


if __name__ == "__main__":
    main()
