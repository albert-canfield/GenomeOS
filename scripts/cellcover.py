# SPDX-License-Identifier: AGPL-3.0-or-later
"""Price a second value for the cells direction rule v2 has only one of (data/results/cellcover.json).

    uv run --frozen python scripts/cellcover.py

Reads the committed data/results/respmap_direction_v2.json rows, the saved copy of the model
client's track metadata, and ONE loose per-element answer as the track-axis witness. It opens no
chromosome archive: `load_cached` inflates whole archives and the population's elements are
archive-only (cellcover.AD_HOC_READ_INCIDENT). Reports, per assigned
cell, how many RNA-seq tracks carry it and how that was established; the number of requests a
second value for the rows that can have one would take; and the rows that cannot have one at any
price. The 19 `one_track_seen_twice` rows are reported apart.

Registered first in data/results/cellcover_registration.json, which is committed before this
script runs. Every binding word is imported from `genomeos.attribution.cellcover`.

0 model requests, no network, no money. No client is imported and no key is read.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import cellcover as cc  # noqa: E402
from genomeos.predict.alphagenome_adapter import CELL_TRACKS  # noqa: E402
from genomeos.predict.enhancer_target import CELLS  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "cellcover"

OWN_CODE = (
    "genomeos/attribution/cellcover.py",
    "scripts/cellcover_register.py",
    "scripts/cellcover.py",
    "tests/test_cellcover.py",
)

#: The population's two counts, as lane-propb committed them. Read and checked, never re-derived.
EXPECTED_ROWS = {cc.ONE_VALUE_ONLY: 93, cc.ONE_TRACK_SEEN_TWICE: 19}

CONTEXT_UNKNOWN_TOKENS = ("", "not recorded", "context_unknown", "CONTEXT_UNKNOWN")


def rows_of(reason: str, result: dict[str, Any]) -> list[dict[str, Any]]:
    """The committed assertion rows carrying one unresolved reason, read from the file's own rows."""
    out = [a for a in result["assertions"] if a.get("v2_unresolved_reason") == reason]
    want = EXPECTED_ROWS[reason]
    if len(out) != want:
        raise SystemExit(
            f"refusing: the committed v2 result carries {len(out)} rows under {reason!r}, not the "
            f"{want} this run was registered over. The population is lane-propb's and is read, not "
            "re-derived; a different count means a different population and the registration no "
            "longer describes it."
        )
    for r in out:
        if str(r.get("cell") or "").strip() in CONTEXT_UNKNOWN_TOKENS:
            raise SystemExit(f"refusing: row {r['assertion_id']!r} names no cell to price")
    return out


def retained_fields(rows: list[dict[str, Any]]) -> dict[str, int]:
    """Each row's retained-field pattern, from the committed v2 rows and not from any cache.

    The committed rows carry `retained_values` with each value's field name, so the pattern that
    shows the assigned cell is v1's own selected extreme is read from a 160 KB JSON file. No element
    answer is opened: the population's elements are archive-only and `load_cached` inflates whole
    chromosome archives (cellcover.AD_HOC_READ_INCIDENT).
    """
    out: Counter[str] = Counter()
    for r in rows:
        out["+".join(v["field"] for v in r["retained_values"])] += 1
    return dict(out)


def report(reason: str, rows: list[dict[str, Any]], axis: cc.Axis, witness: int) -> dict[str, Any]:
    """One partition's figures: per-cell coverage, the split, the price."""
    cover = cc.coverage(Counter(str(r["cell"]) for r in rows), axis)
    parts = cc.split(rows, cover)
    fields = retained_fields(rows)
    rep = cc.Report(
        reason=reason,
        rows=len(rows),
        cells=cover,
        purchasable=cc.price(parts["purchasable"]),
        unpurchasable=cc.price(parts["unpurchasable_at_any_price"]),
        axis_width=axis.width,
        witness_width=witness,
    )
    return {
        "reason": reason,
        "counts": rep.counts(),
        "retained_field_patterns": fields,
        "per_cell": [c.to_dict() for c in cover],
        "purchasable": rep.purchasable.to_dict(),
        "unpurchasable_at_any_price": {
            **rep.unpurchasable.to_dict(),
            "requests": 0,
            "request_unit": "no request buys a second value for these rows: one track is one column "
            "of every response, now and in future",
            "elements": [],
            "elements_named_but_not_priced": list(rep.unpurchasable.elements),
        },
    }


def payload() -> dict[str, Any]:
    result = json.loads(cc.V2_RESULT.read_text())
    axis = cc.rna_axis()
    if axis.width != cc.EXPECTED_AXIS:
        raise SystemExit(
            f"refusing: the saved track metadata holds {axis.width} distinct RNA-seq track names, "
            f"not {cc.EXPECTED_AXIS}. AXIS_IDENTIFIED rests on that number equalling the cached "
            "rows' own n_tracks, and the per-cell counts are not reported without it."
        )
    witness = cc.axis_witness()
    if witness["tracks"] != cc.EXPECTED_AXIS:
        raise SystemExit(
            f"refusing: the axis witness records {witness['tracks']} tracks, not {cc.EXPECTED_AXIS}. "
            "AXIS_IDENTIFIED rests on that equalling the metadata's distinct RNA-seq track-name "
            "count, and the per-cell counts are not reported without it."
        )
    main_rows = rows_of(cc.ONE_VALUE_ONLY, result)
    nineteen = rows_of(cc.ONE_TRACK_SEEN_TWICE, result)
    entries = [
        mf.input_entry(cc.V2_RESULT, partition="the committed v2 classification, the population"),
        mf.input_entry(cc.TRACK_METADATA, partition="the saved copy of the client's track metadata"),
        mf.input_entry(
            cc.AXIS_WITNESS,
            partition="one loose per-element answer, read for its response record's axis width; it "
            "is not a row of the population",
        ),
    ]
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-cellcover",
        "registration": "data/results/cellcover_registration.json, committed before this run",
        "question": cc.QUESTION,
        "population": cc.POPULATION,
        "premise_corrected": cc.PREMISE_CORRECTED,
        "request_unit": cc.REQUEST_UNIT,
        "purchasable_rule": cc.PURCHASABLE_RULE,
        "unpurchasable": cc.UNPURCHASABLE,
        "what_it_does_not_buy": cc.WHAT_IT_DOES_NOT_BUY,
        "how_established": cc.HOW_ESTABLISHED,
        "track_label_rule": cc.TRACK_LABEL_RULE,
        "axis_identified": cc.AXIS_IDENTIFIED,
        "metadata_copy_limitation": cc.METADATA_COPY_LIMITATION,
        "the_nineteen": cc.THE_NINETEEN,
        "no_recommendation": cc.NO_RECOMMENDATION,
        "cell2_group_not_used": cc.CELL2_GROUP_NOT_USED,
        "readings_carried_unstrengthened": cc.READINGS_CARRIED,
        "cannot_establish": list(cc.CANNOT_ESTABLISH),
        "retained_cells_today": list(CELLS),
        "cells_whose_columns_a_fresh_response_already_records": list(CELL_TRACKS),
        "axis": {
            "rna_seq_metadata_rows": axis.rows_read,
            "distinct_track_names": axis.width,
            "source": axis.source,
            "note": cc.AXIS_IDENTIFIED,
            "response_record_axis_width": witness["tracks"],
            "response_record_axis_sha256": witness["tracks_sha256"],
        },
        "axis_witness": witness,
        "axis_witness_is_not_the_population": cc.AXIS_WITNESS_IS_NOT_THE_POPULATION,
        "ad_hoc_read_incident": cc.AD_HOC_READ_INCIDENT,
        "measured_by_a_read_now_forbidden": cc.MEASURED_BY_A_READ_NOW_FORBIDDEN,
        "one_value_only": report(cc.ONE_VALUE_ONLY, main_rows, axis, witness["tracks"]),
        "one_track_seen_twice": report(cc.ONE_TRACK_SEEN_TWICE, nineteen, axis, witness["tracks"]),
        "alphagenome_requests": 0,
        "money": "none: no request was sent and no key was read",
        "result_manifest": {
            "sources": [
                {
                    "accession": "the committed direction-rule-v2 classification of the published "
                    "response map (230509a)",
                    "version": "the committed data/results/respmap_direction_v2.json, digested in inputs",
                },
                {
                    "accession": "the model client's track metadata as a peer session saved it",
                    "version": "the saved copy on this machine, digested in inputs; whether it is the "
                    "client's current table is not established",
                },
                {
                    "accession": "one answer of the deletion sweep, as the axis witness",
                    "version": "the loose per-element JSON digested in inputs",
                },
            ],
            "inputs": entries,
            "input_count": len(entries),
            "assembly": "GRCh38: the assembly the cached sweep and the response map were written on. "
            "No coordinate of this result's own is in GRCh38 or any build",
            "coordinates": "n/a: no interval is read, written or compared; a row is an element id and "
            "a gene symbol, a cell is a label",
            "parameters": {
                "values_needed_by_clause_1": cc.VALUES_NEEDED,
                "expected_axis_width": cc.EXPECTED_AXIS,
                "rna_output": cc.RNA_OUTPUT,
                "reasons_priced": [cc.ONE_VALUE_ONLY, cc.ONE_TRACK_SEEN_TWICE],
                "retained_cells_today": list(CELLS),
                "no_threshold_set_here": "no magnitude, consistency or confidence threshold is set or "
                "moved: no value is read at all, only track counts and row counts",
            },
            "exclusions": [
                "the 54 resolved rows of the published 166: not unresolved, so not priced",
                "the six reason tokens of direction_v2.UNRESOLVED_REASONS no row of this population "
                "carries: excluded by being absent, not by a choice",
                "every metadata row whose output is not rna_seq: the sweep's scorer is the RNA_SEQ "
                "gene scorer, so no other output is a column a cell's value could come from",
            ],
            "partitions": {
                "one_value_only": "the 93 rows, priced per cell and per element",
                "one_track_seen_twice": "the 19 rows, reported apart and never pooled with the 93",
            },
            "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        },
    }


def main() -> int:
    p = payload()
    path = save_result(RESULT, p)
    for key in ("one_value_only", "one_track_seen_twice"):
        c = p[key]["counts"]
        print(f"{key}: {json.dumps(c)}")
    print(f"written: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
