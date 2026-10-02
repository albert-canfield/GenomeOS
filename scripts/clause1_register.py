# SPDX-License-Identifier: AGPL-3.0-or-later
"""Register the clause-(1) count before it is taken (data/results/clause1_registration.json).

    uv run --frozen python scripts/clause1_register.py

Every binding word is imported from `genomeos.attribution.clause1`, the module that applies them,
so the registration and the code cannot drift apart. This file names the population, the five track
classes and the rule that decides them, the rule by which a row satisfies clause (1), the rule by
which the assigned cell counts as an argmax, the way numerical distinctness is established, every
count the run will report and every refusal that stops it reporting.

`EXPLORATION_PRECEDED_THIS` in the module states plainly that this lane read the committed v2 rows,
the committed cellcover prose and the metadata copy's header before this file was written, and names
the figures it had seen. This registration is therefore NOT a prediction of those figures. What it
binds is the class set, the rules, the count list and the refusals, so the committed program's
output can be held against written rules rather than against this lane's recollection.

This lane edits no rule: `attribution/direction_v2.py` is imported and read, and no clause,
amendment or threshold of it is changed here or anywhere by this lane.

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
from genomeos.attribution import clause1 as c1  # noqa: E402
from genomeos.attribution import direction_v2 as dv2  # noqa: E402
from genomeos.predict.enhancer_target import CELLS  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "clause1_registration"
WILL_WRITE = "data/results/clause1.json"

#: This lane's own files. Everything else uncommitted in this shared checkout belongs to a peer.
OWN_CODE = (
    "genomeos/attribution/clause1.py",
    "scripts/clause1_register.py",
    "scripts/clause1_count.py",
    "tests/test_clause1.py",
)


def inputs() -> list[dict[str, Any]]:
    """The named files the run will read, digested before it runs."""
    out: list[dict[str, Any]] = []
    if c1.V2_RESULT.exists():
        out.append(
            mf.input_entry(
                c1.V2_RESULT,
                partition="the committed v2 classification of the published response map; its 166 "
                "predicted rows are this lane's whole population and are read, never re-derived",
            )
        )
    if c1.TRACK_METADATA.exists():
        out.append(
            mf.input_entry(
                c1.TRACK_METADATA,
                partition="the saved copy of the model client's own track metadata, read for its "
                "RNA-seq rows only; re-fetching it is a request this lane may not make",
            )
        )
    return out


def payload() -> dict[str, Any]:
    """The registration, every binding word imported from the module that applies it."""
    entries = inputs()
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-clause1",
        "question": c1.QUESTION,
        "population": c1.POPULATION,
        "denominator_is_separate": c1.DENOMINATOR_IS_SEPARATE,
        "exploration_preceded_this": c1.EXPLORATION_PRECEDED_THIS,
        "edits_no_rule": (
            "this lane changes no clause, no amendment and no threshold of "
            "genomeos/attribution/direction_v2.py or of any other rule. The module is imported and "
            "read. If the answer is that clause (1) does not secure its stated reason, that is a "
            "finding for Albert and the supervisor and not a licence to change the rule: changing "
            "it would invalidate every count already taken under it"
        ),
        "clause_1_quoted": (
            "SIGN_RULE_V2 clause (1), quoted from direction_v2: `at least two values are retained "
            "for that cell name, so that no direction rests on one track - which is the defect v1 "
            "has`. The words after the comma are the clause's stated reason"
        ),
        "sign_rule_v2": dv2.SIGN_RULE_V2,
        "amendment_1": dv2.AMENDMENT_1,
        "cell_evidence": dv2.CELL_EVIDENCE,
        "assigned_cell": dv2.ASSIGNED_CELL,
        "track_label_rule": c1.TRACK_LABEL_RULE,
        "name_match_is_not_a_confirmation": c1.NAME_MATCH_IS_NOT_A_CONFIRMATION,
        "metadata_copy_limitation": c1.METADATA_COPY_LIMITATION,
        "track_classes": list(c1.TRACK_CLASSES),
        "class_rule": c1.CLASS_RULE,
        "per_cell_fields": list(c1.PER_CELL_FIELDS),
        "established_per_cell": c1.ESTABLISHED_PER_CELL,
        "clause_1_satisfied_rule": c1.CLAUSE_1_SATISFIED_RULE,
        "argmax_rule": c1.ARGMAX_RULE,
        "argmax_is_a_selection": c1.ARGMAX_IS_A_SELECTION,
        "outside_retained_cells_cannot_resolve": c1.OUTSIDE_RETAINED_CELLS_CANNOT_RESOLVE,
        "distinctness_rule": c1.DISTINCTNESS_RULE,
        "readings_carried_unstrengthened": c1.READINGS_CARRIED,
        "no_recommendation": c1.NO_RECOMMENDATION,
        "validates_nothing": c1.VALIDATES_NOTHING,
        "cannot_establish": list(c1.CANNOT_ESTABLISH),
        "counts_named_in_advance": list(c1.COUNTS_NAMED),
        "refusals": list(c1.REFUSALS),
        "committed_v2_figures_checked_first": c1.COMMITTED_V2_FIGURES,
        "retained_cells_today": list(CELLS),
        "values_needed_by_clause_1": cc.VALUES_NEEDED,
        "expected_axis_width": cc.EXPECTED_AXIS,
        "what_the_run_will_write": WILL_WRITE,
        "requests": {
            "alphagenome_requests": 0,
            "money": "none: no request is sent and no key is read",
            "network": "none: every input is already on this machine",
        },
        "no_archive_is_opened": (
            "no element answer is read through load_cached or any archive loader. A peer's inline "
            "loop over 93 elements reached 23 GB RSS tonight because that loader falls back to "
            "whole chromosome archives. Everything this lane needs is in the two committed inputs "
            "digested above, and the run opens nothing else"
        ),
        "result_manifest": {
            "sources": [
                {
                    "accession": "the committed direction-rule-v2 classification of the published "
                    "increment-2 response map",
                    "version": "the committed data/results/respmap_direction_v2.json on this "
                    "machine, digested in inputs",
                },
                {
                    "accession": "the model client's track metadata (output_metadata) as a peer "
                    "session saved it to disk",
                    "version": "the saved copy on this machine, digested in inputs; whether it is "
                    "the client's current table is not established (metadata_copy_limitation)",
                },
            ],
            "inputs": entries,
            "input_count": len(entries),
            "assembly": "GRCh38: the assembly the cached deletion sweep and the response map were "
            "written on. No coordinate of this registration's own is in GRCh38 or any build",
            "coordinates": "n/a: no interval is read, written or compared. A row is identified by "
            "its assertion id and a cell by its label",
            "parameters": {
                "values_needed_by_clause_1": cc.VALUES_NEEDED,
                "expected_axis_width": cc.EXPECTED_AXIS,
                "rna_output": cc.RNA_OUTPUT,
                "track_classes": list(c1.TRACK_CLASSES),
                "retained_cells_today": list(CELLS),
                "forbidden_column_fragment": c1.FORBIDDEN_COLUMN_FRAGMENT,
                "no_threshold_set_here": "no magnitude, consistency or confidence threshold is set "
                "or moved. MAGNITUDE_FLOOR stays the imported MIN_EFFECT and consistency stays "
                "unanimity; this lane reads values only to ask whether two of them are equal",
            },
            "exclusions": [
                "every assertion of the published response map whose status is not `predicted`: "
                "they are outside the v2 result this lane reads and so outside this population",
                "every metadata row whose output is not rna_seq: the deletion sweep's scorer is the "
                "RNA_SEQ gene scorer, so no other output is a column a cell's value could come from",
                "every cached element answer and every chromosome archive: none is opened "
                "(no_archive_is_opened)",
            ],
            "partitions": {
                "the_166_rows": "this lane's row denominator, read from the committed v2 result",
                "the_distinct_assigned_cells": "this lane's cell denominator, the distinct cells of "
                "those 166 rows",
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
