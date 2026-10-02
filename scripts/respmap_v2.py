# SPDX-License-Identifier: AGPL-3.0-or-later
"""Emit direction rule v2's class for every published `predicted` assertion of the response map
(data/results/respmap_direction_v2.json).

    uv run --frozen python scripts/respmap_v2.py

One row per assertion - element, gene, cell, the action the map published, the v2 class and, when
v2 did not resolve, its reason - so the count is read off committed rows and never off a
`collections.Counter` total. Registered first, in
`data/results/respmap_direction_v2_registration.json`, which names the population, the denominator,
the class set, every per-row field and every count before any of them was taken.

The rule is `attribution.direction_v2`, imported and applied unchanged. Nothing is re-scored, no
rule or threshold is edited, no published figure is touched, and the response-map payload is read
and never written. 0 model requests, no network, no money: the two inputs are the committed
response map and the finished deletion sweep's own per-element cache, both on disk.
"""

from __future__ import annotations

import json
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import direction_v2 as dv  # noqa: E402
from genomeos.attribution import respmap_v2 as rv  # noqa: E402
from genomeos.predict import enhancer_target as et  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = rv.RESULT
REGISTRATION = "data/results/respmap_direction_v2_registration.json"
SOURCE = "data/results/response_map_increment2.json"
CACHE_DIR = Path("data/knowledge/alphagenome/elements")

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to another lane.
OWN_CODE = (
    "genomeos/attribution/respmap_v2.py",
    "scripts/respmap_v2_register.py",
    "scripts/respmap_v2.py",
    "tests/test_respmap_v2.py",
)


def inputs() -> list[dict[str, Any]]:
    """The files this run read, by their own paths and digests."""
    out: list[dict[str, Any]] = []
    for name, part in (
        (REGISTRATION, "this run's own registration, committed before it ran"),
        (SOURCE, "the published response map whose predicted assertions are counted"),
    ):
        p = Path(name)
        if p.exists():
            out.append(mf.input_entry(p, partition=part))
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


def main() -> int:
    started = time.time()
    source = json.loads(Path(SOURCE).read_text())
    out = rv.emit(source)
    rows, counts = out["rows"], out["counts"]
    entries = inputs()
    payload: dict[str, Any] = {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-propb",
        "question": rv.QUESTION,
        "registration": REGISTRATION,
        "population": rv.POPULATION,
        "increment_3_is_not_this_population": rv.INCREMENT_3_IS_NOT_THIS_POPULATION,
        "classes": list(rv.CLASSES),
        "class_rule": rv.CLASS_RULE,
        "counts": counts,
        "flip_is_reported_as_measured": rv.FLIP_IS_REPORTED_AS_MEASURED,
        "absence_gloss_corrected": rv.ABSENCE_GLOSS_CORRECTED,
        "earlier_tallies_not_touched": rv.EARLIER_TALLIES_NOT_TOUCHED,
        "selection_excluded_diagnostic": rv.DIAGNOSTIC_IS_NOT_A_V2_OUTPUT,
        "unresolved_means": rv.UNRESOLVED_MEANS,
        "nothing_rescored": rv.NOTHING_RESCORED,
        "validates_nothing": rv.VALIDATES_NOTHING,
        "not_an_improvement": dv.NOT_AN_IMPROVEMENT,
        "wording": dv.WORDING,
        "no_recommendation": rv.NO_RECOMMENDATION,
        "cannot_establish": list(rv.CANNOT_ESTABLISH),
        "assertions": rows,
        "alphagenome_requests": 0,
        "money": "none: every input read was already on disk",
        "seconds": round(time.time() - started, 1),
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
                "denominator": counts["denominator"],
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
    path = save_result(RESULT, payload)
    print(f"{RESULT}: {path}")
    print(f"  denominator: {counts['denominator']} published `predicted` assertions")
    for k in rv.CLASSES:
        print(f"    {k}: {counts['by_class'][k]}")
    print("  unresolved by reason:")
    for k, v in sorted(counts["unresolved_by_reason"].items(), key=lambda kv: -kv[1]):
        print(f"    {k}: {v}")
    print(f"  reasons no assertion carries: {', '.join(counts['unresolved_reasons_with_no_assertion'])}")
    for axis, block in counts["by_published_axis"].items():
        print(f"  {axis}: {block}")
    print(f"  no cached row: {counts['no_cached_row']}")
    print(f"  published value among the cell's retained values: {counts['published_value_among_retained']}")
    print(
        "  one_value_only whose single retained value is by_cell and opposite in sign: "
        f"{counts['unresolved_one_value_only_by_cell_opposite_to_published']}"
    )
    print(
        "  selection-excluded diagnostic hits (not a v2 output): "
        f"{counts['selection_excluded_diagnostic_hits']}"
    )
    print(f"  cells that are cached tracks: {counts['cells_that_are_cached_tracks']}")
    print()
    print(rv.FLIP_IS_REPORTED_AS_MEASURED)
    print(rv.VALIDATES_NOTHING)
    print(rv.NO_RECOMMENDATION)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
