# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the cell-coverage costing before it is run (data/results/cellcover_registration.json).

    uv run --frozen python scripts/cellcover_register.py

This registration names the population, the rule by which a cell is priced, the request unit, the
per-cell fields and every count the run will report, and it names them before the run writes a
figure. Every binding word is imported from `genomeos.attribution.cellcover`, the module that
applies them, so the registration and the code cannot drift apart.

Honesty about the order, recorded here rather than claimed away
---------------------------------------------------------------
`EXPLORATION_PRECEDED_THIS` states plainly that this lane read the committed v2 rows and the saved
track metadata interactively before writing this file, and names the figures it saw. This
registration is therefore NOT a prediction of those figures and may not be read as one. What it
binds is the rule, the fields and the count list, so that the committed program's output can be
held against a written rule rather than against this lane's word. The separate commit this file
gets before the run is what makes the order checkable in git; it does not make the figures
unseen, and saying otherwise would be the error it is meant to prevent.

No request is sent, no key is read, no figure is computed and no published number is touched.
0 model requests, no network, no money.
"""

from __future__ import annotations

import sys
from datetime import date
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import cellcover as cc  # noqa: E402
from genomeos.predict.enhancer_target import CELLS  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "cellcover_registration"
WILL_WRITE = "data/results/cellcover.json"

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to a peer.
OWN_CODE = (
    "genomeos/attribution/cellcover.py",
    "scripts/cellcover_register.py",
    "scripts/cellcover.py",
    "tests/test_cellcover.py",
)

EXPLORATION_PRECEDED_THIS = (
    "stated against this lane's own interest. This lane read the committed "
    "data/results/respmap_direction_v2.json rows and the saved track metadata interactively BEFORE "
    "writing this registration, and had already seen: 93 rows under `one_value_only` over 59 "
    "distinct cells, 19 under `one_track_seen_twice`, 371 distinct RNA-seq track names, and the "
    "split of the 59 into 44 one-track cells and 15 two-track cells. So this file does not predict "
    "those figures and no reader may treat it as having done so. What it fixes in advance of the "
    "COMMITTED PROGRAM is the rule by which a cell is priced (PURCHASABLE_RULE), how a track count "
    "may be established (HOW_ESTABLISHED), the label rule (TRACK_LABEL_RULE), the request unit "
    "(REQUEST_UNIT), the per-cell fields and the count list - so that the program's output is held "
    "against written rules rather than against this lane's recollection, and so that a figure "
    "produced by a rule invented after seeing the data would be visible as a departure from this "
    "file. The separate commit is the audit of the ORDER OF THE RULE, not of the order of the "
    "figures"
)

PER_CELL_FIELDS = (
    "cell",
    "rows",
    "rna_seq_tracks",
    "track_exists",
    "established_how",
    "second_value_exists_in_any_response",
    "ontology_terms",
    "assay_titles",
    "data_sources",
    "biosample_types",
    "label_via",
)

REFUSALS = (
    "the program refuses rather than reporting, if the saved metadata's distinct RNA-seq track-name "
    f"count is not {cc.EXPECTED_AXIS}, or if any cached gene row of the population records a "
    f"different n_tracks: the per-cell counts are grounded on that equality (AXIS_IDENTIFIED) and "
    "must not be reported without it",
    "the program refuses if the committed v2 result's two reason counts are not 93 and 19: the "
    "population is lane-propb's and is read, never re-derived",
    "the program refuses if any row's assigned cell is the empty string or the compiler's "
    "context-unknown token, because such a row names no cell to price",
)

NOT_PRICED = (
    "the 54 resolved rows of the 166 are not priced, no count here is added to lane-propb's 112, 93 "
    "or 19, and the earlier 14 and 10 of other lanes are not touched, not added and not mentioned "
    "as comparable: they overlap each other and belong to other populations"
)


def inputs() -> list[dict[str, Any]]:
    """The named files the run will read, digested before it runs."""
    out: list[dict[str, Any]] = []
    if cc.V2_RESULT.exists():
        out.append(
            mf.input_entry(
                cc.V2_RESULT,
                partition="the committed v2 classification whose unresolved rows are the population",
            )
        )
    if cc.TRACK_METADATA.exists():
        out.append(
            mf.input_entry(
                cc.TRACK_METADATA,
                partition="the saved copy of the model client's own track metadata, read for its "
                "RNA-seq rows only; re-fetching it is a request this lane may not make",
            )
        )
    if cc.AXIS_WITNESS.exists():
        out.append(
            mf.input_entry(
                cc.AXIS_WITNESS,
                partition="one loose per-element answer, read for its response record's own track-axis "
                "width; it is not a row of the population",
            )
        )
    return out


def payload() -> dict[str, Any]:
    """The registration, every binding word imported from the module that applies it."""
    entries = inputs()
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-cellcover",
        "question": cc.QUESTION,
        "population": cc.POPULATION,
        "exploration_preceded_this": EXPLORATION_PRECEDED_THIS,
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
        "axis_witness": str(cc.AXIS_WITNESS),
        "axis_witness_is_not_the_population": cc.AXIS_WITNESS_IS_NOT_THE_POPULATION,
        "ad_hoc_read_incident": cc.AD_HOC_READ_INCIDENT,
        "measured_by_a_read_now_forbidden": cc.MEASURED_BY_A_READ_NOW_FORBIDDEN,
        "no_recommendation": cc.NO_RECOMMENDATION,
        "cell2_group_not_used": cc.CELL2_GROUP_NOT_USED,
        "readings_carried_unstrengthened": cc.READINGS_CARRIED,
        "cannot_establish": list(cc.CANNOT_ESTABLISH),
        "not_priced": NOT_PRICED,
        "retained_cells_today": list(CELLS),
        "values_needed_by_clause_1": cc.VALUES_NEEDED,
        "expected_axis_width": cc.EXPECTED_AXIS,
        "per_cell_fields": list(PER_CELL_FIELDS),
        "counts_named_in_advance": list(cc.COUNTS_NAMED),
        "refusals": list(REFUSALS),
        "what_the_run_will_write": WILL_WRITE,
        "requests": {
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
            "network": "none: every input is already on this machine",
        },
        "result_manifest": {
            "sources": [
                {
                    "accession": "the committed direction-rule-v2 classification of the published "
                    "response map, as scripts/respmap_v2.py wrote it at 230509a",
                    "version": "the committed data/results/respmap_direction_v2.json on this machine, "
                    "digested in inputs",
                },
                {
                    "accession": "the model client's track metadata (output_metadata) as a peer "
                    "session saved it to disk",
                    "version": "the saved copy on this machine, digested in inputs; whether it is the "
                    "client's current table is not established (metadata_copy_limitation)",
                },
                {
                    "accession": "one answer of the deletion sweep, read as the track-axis witness",
                    "version": "the loose per-element JSON digested in inputs; no chromosome archive "
                    "is opened, because the population's elements are archive-only and the loader "
                    "inflates whole archives (cellcover.AD_HOC_READ_INCIDENT)",
                },
            ],
            "inputs": entries,
            "input_count": len(entries),
            "assembly": "GRCh38: the assembly the cached deletion sweep and the response map were "
            "written on. No coordinate of this registration's own is in GRCh38 or any build",
            "coordinates": "n/a: no interval is read, written or compared. A row is identified by its "
            "element id and target gene symbol, and a cell by its label",
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
                "the 54 resolved rows of the published 166: they are not unresolved and so are not priced",
                "the six reason tokens of direction_v2.UNRESOLVED_REASONS that no row of this "
                "population carries: they are excluded by being absent, not by a choice",
                "every metadata row whose output is not rna_seq (chip_histone, cage, dnase, atac, "
                "contact_maps, procap): the deletion sweep's scorer is the RNA_SEQ gene scorer, so no "
                "other output is a column of the response a cell's value could come from",
            ],
            "partitions": {
                "one_value_only": "the 93 rows, priced per cell and per element",
                "one_track_seen_twice": "the 19 rows, reported apart and never pooled with the 93",
            },
            "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        },
    }


def main() -> int:
    path = save_result(RESULT, payload())
    print(f"registered: {path}")
    print(f"the run will write: {WILL_WRITE}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
