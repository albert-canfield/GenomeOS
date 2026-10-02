# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the per-assertion direction-v2 emitter before it reads a single direction
(data/results/respmap_direction_v2_registration.json).

    uv run --frozen python scripts/respmap_v2_register.py

Why an emitter had to be built at all: `scripts/direction_v2.py` accumulates two
`collections.Counter`s and writes no file, and no committed result carries the v2 class of any one
assertion. Searched across `data/results`, each of the nine reason tokens of
`attribution.direction_v2.UNRESOLVED_REASONS` appears in 0 files. A class that cannot be read
cannot be reported, and a per-assertion class may not be re-derived from a Counter total.

This registration names the population, the denominator, the exhaustive class set, the per-row
fields and every count that the run will report, and it names them before the run. The question,
the population, the class rule and the readings all come from `genomeos.attribution.respmap_v2`,
so the registration and the code that applies it cannot drift apart. The rule itself is
`attribution.direction_v2`, imported unchanged: no threshold is set here and none is moved.

No direction is read, no class is computed, no count is taken, no rule is edited and no published
figure is touched. 0 model requests, no network, no money.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import direction_v2 as dv  # noqa: E402
from genomeos.attribution import respmap_v2 as rv  # noqa: E402
from genomeos.predict import enhancer_target as et  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "respmap_direction_v2_registration"
WILL_WRITE = "data/results/respmap_direction_v2.json"
SOURCE = "data/results/response_map_increment2.json"
CACHE_DIR = Path("data/knowledge/alphagenome/elements")

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/respmap_v2.py",
    "scripts/respmap_v2_register.py",
    "scripts/respmap_v2.py",
    "tests/test_respmap_v2.py",
)

#: Every per-row field the run will write, named before it writes one.
ROW_FIELDS = (
    "assertion_id",
    "chrom",
    "element",
    "gene",
    "cell",
    "published_action",
    "published_axis",
    "v1_action",
    "published_value",
    "v2_class",
    "v2_action",
    "v2_unresolved_reason",
    "cached_row_found",
    "retained_values",
    "published_value_is_among_retained",
    "only_retained_is_by_cell_opposite_to_published",
    "by_cell_value_opposite_to_published",
    "cell_is_a_cached_track",
)

#: Every count the run will report, named before it is taken.
COUNTS_NAMED = (
    "denominator: the published predicted assertions classified, one row each",
    "by_class: one count per member of respmap_v2.CLASSES, exhaustive and mutually exclusive",
    "unresolved_by_reason: one count per reason of direction_v2.UNRESOLVED_REASONS that any row "
    "carries, with the reasons no row carries listed apart so that an absent reason is visibly "
    "absent and not missing",
    "by_published_axis: the same class counts split by the axis the map published, reported apart "
    "and never pooled",
    "no_cached_row: rows for which the deletion cache holds no row for the named target",
    "published_value_among_retained: rows whose published value is numerically among the cell's "
    "retained values, which is what FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE rests on",
    "unresolved_one_value_only_by_cell_opposite_to_published: the escape branch of "
    "ABSENCE_GLOSS_CORRECTED, counted exactly over this population",
    "selection_excluded_diagnostic_hits: the diagnostic of "
    "direction_v2.SELECTION_EXCLUDED_DIAGNOSTIC over this population, a property of the cache and "
    "not a v2 output",
    "cells_that_are_cached_tracks: rows whose assigned cell is one of predict.enhancer_target.CELLS",
)


def inputs() -> list[dict[str, Any]]:
    """The named files the run will read, digested before it runs."""
    out: list[dict[str, Any]] = []
    p = Path(SOURCE)
    if p.exists():
        out.append(
            mf.input_entry(p, partition="the published response map whose predicted assertions are counted")
        )
    archives = sorted(q.as_posix() for q in CACHE_DIR.glob("chr*.json.gz"))
    if archives:
        out.append(
            mf.files_entry(
                "per_element_response_cache",
                archives,
                partition="the finished deletion sweep's own per-element cache, read for a sign only",
            )
        )
    return out


def payload() -> dict[str, Any]:
    """The registration, every binding word of it imported from the module that applies it."""
    entries = inputs()
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-propb",
        "question": rv.QUESTION,
        "why_an_emitter_had_to_be_built": (
            "scripts/direction_v2.py accumulates two collections.Counter objects and writes no file: "
            "it calls no json.dump, opens no file for writing, and records no element id, locus, "
            "gene or cell for any call. Each of the nine reason tokens of "
            "direction_v2.UNRESOLVED_REASONS was searched across data/results and appears in 0 "
            "files. So no committed artefact carried the v2 class of any one published assertion, "
            "and the count could not be taken by reading one. A class is therefore emitted per "
            "assertion and read back from committed rows; it is never re-derived from a Counter total"
        ),
        "population": rv.POPULATION,
        "denominator": {
            "value": 166,
            "rule": rv.POPULATION,
            "fixed_by": "the published payload's own counts.by_status, committed before this lane ran",
            "the_other_361_are_observed": (
                "361 of the 527 published assertions are `observed` - 195 CRISPRi perturbations and "
                "166 ENCODE DNase readings - and `observed` is reserved for what an experiment "
                "measured, so they carry no AlphaGenome direction for v2 to leave unresolved"
            ),
        },
        "increment_3_is_not_this_population": rv.INCREMENT_3_IS_NOT_THIS_POPULATION,
        "classes": list(rv.CLASSES),
        "class_rule": rv.CLASS_RULE,
        "unresolved_reasons": list(dv.UNRESOLVED_REASONS),
        "sign_rule_v2": dv.CONSISTENCY_RULE,
        "magnitude_floor": {"value": dv.MAGNITUDE_FLOOR, "imported": dv.MAGNITUDE_FLOOR_IS_IMPORTED},
        "amendment_1": dv.AMENDMENT_1,
        "unresolved_means": rv.UNRESOLVED_MEANS,
        "unresolved_is_first_class": dv.UNRESOLVED_IS_FIRST_CLASS,
        "flip_is_reported_as_measured": rv.FLIP_IS_REPORTED_AS_MEASURED,
        "absence_gloss_corrected": rv.ABSENCE_GLOSS_CORRECTED,
        "earlier_tallies_not_touched": rv.EARLIER_TALLIES_NOT_TOUCHED,
        "selection_excluded_diagnostic": rv.DIAGNOSTIC_IS_NOT_A_V2_OUTPUT,
        "row_fields": list(ROW_FIELDS),
        "counts_named_in_advance": list(COUNTS_NAMED),
        "nothing_rescored": rv.NOTHING_RESCORED,
        "validates_nothing": rv.VALIDATES_NOTHING,
        "not_an_improvement": dv.NOT_AN_IMPROVEMENT,
        "wording": dv.WORDING,
        "no_recommendation": rv.NO_RECOMMENDATION,
        "cannot_establish": list(rv.CANNOT_ESTABLISH),
        "requests": {"alphagenome": 0, "network": "none", "money": "none: every input is on disk"},
        "what_the_run_will_write": WILL_WRITE,
        "frozen_inputs": [
            {k: e[k] for k in ("path", "label", "sha256", "bytes", "files") if k in e} for e in entries
        ],
        "result_manifest": {
            "sources": [
                {
                    "accession": "the published increment-2 response map, as "
                    "scripts/response_map_increment2.py wrote it",
                    "version": "the committed data/results/response_map_increment2.json on this "
                    "machine, digested in inputs",
                },
                {
                    "accession": "the finished genome-wide AlphaGenome deletion sweep's per-element "
                    "response cache",
                    "version": "the per-chromosome archives on this machine, digested in inputs as a group",
                },
            ],
            "inputs": entries,
            "assembly": "GRCh38",
            "coordinates": {"base": 0, "interval": "half-open"},
            "parameters": {
                "population": rv.POPULATION,
                "denominator": 166,
                "classes": list(rv.CLASSES),
                "class_rule": rv.CLASS_RULE,
                "sign_rule": dv.CONSISTENCY_RULE,
                "magnitude_floor": dv.MAGNITUDE_FLOOR,
                "magnitude_floor_imported_from": "genomeos.predict.enhancer_target.MIN_EFFECT",
                "cached_cell_tracks": list(et.CELLS),
                "unresolved_reasons": list(dv.UNRESOLVED_REASONS),
                "aggregate_function": "a count of rows per class; no mean, no rate, no interval",
                "fill_value": "none: an assertion with no cached row is classified unresolved under "
                "its own reason no_row_for_target and is never dropped or imputed",
            },
            "exclusions": [
                rv.POPULATION,
                "the 361 `observed` assertions are outside the question by their status: an "
                "experiment's measurement carries no AlphaGenome direction for v2 to withhold",
                "no assertion is excluded for its magnitude, its cell, its chromosome or its "
                "class: all 166 are classified and every class is reported",
                rv.EARLIER_TALLIES_NOT_TOUCHED,
            ],
            "partitions": {
                "by_published_axis": "the axis the map published for the assertion, "
                "activates_target or represses_target, reported apart and never pooled",
                "by_class": "respmap_v2.CLASSES, exhaustive and mutually exclusive",
                "unresolved_by_reason": "direction_v2.UNRESOLVED_REASONS, one bucket per refusal",
            },
            "code_cleanliness": mf.code_cleanliness(Path(__file__).resolve(), OWN_CODE),
        },
    }


def main() -> int:
    path = save_result(RESULT, payload())
    print(f"{RESULT}: {path}")
    print(f"  population: {166} published `predicted` assertions of {SOURCE}")
    print(f"  classes: {', '.join(rv.CLASSES)}")
    print(f"  unresolved reasons: {len(dv.UNRESOLVED_REASONS)}, imported from direction_v2")
    print(f"  floor: {dv.MAGNITUDE_FLOOR}, imported; consistency is unanimity, nothing to tune")
    print(f"  row fields named in advance: {len(ROW_FIELDS)}")
    print(f"  counts named in advance: {len(COUNTS_NAMED)}")
    print(f"  the run will write: {WILL_WRITE}")
    print("  no direction read, no class computed, no count taken, 0 model requests, no money")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
