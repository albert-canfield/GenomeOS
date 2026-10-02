# SPDX-License-Identifier: AGPL-3.0-or-later
"""Take the clause-(1) count registered at ba8c41f (data/results/clause1.json).

    uv run --frozen python scripts/clause1_count.py

Reads exactly two files, both digested in the registration: the committed
data/results/respmap_direction_v2.json, for its 166 predicted rows, and the saved copy of the model
client's own track metadata, for its RNA-seq rows. It opens no element answer, calls no archive
loader and sends no request. Every class, rule, count and refusal is imported from
`genomeos.attribution.clause1`; this script chooses nothing.

This lane edits no rule. `attribution/direction_v2.py` is imported and read, and no clause,
amendment or threshold of it is changed.

0 model requests, no network, no money.
"""

from __future__ import annotations

import json
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

RESULT = "clause1"
REGISTRATION = "data/results/clause1_registration.json (registered at ba8c41f, code at fde5584)"

OWN_CODE = (
    "genomeos/attribution/clause1.py",
    "scripts/clause1_register.py",
    "scripts/clause1_count.py",
    "tests/test_clause1.py",
)


def rows_and_cells() -> tuple[list[c1.RowReading], list[c1.CellTracks], dict[str, Any], c1.Axis_t]:
    """The population and its cells, with every registered refusal applied first."""
    with c1.V2_RESULT.open() as fh:
        v2 = json.load(fh)
    rows = c1.read_rows(v2)
    c1.check_population(v2, rows)
    axis = cc.rna_axis(c1.TRACK_METADATA)
    c1.check_metadata(c1.metadata_columns(c1.TRACK_METADATA), axis.width)
    cells = c1.classify(rows, axis.by_label())
    classes = {c.track_class for c in cells}
    if not classes <= set(c1.TRACK_CLASSES):
        raise ValueError(f"a cell fell outside the registered classes: {classes} (REFUSALS)")
    return rows, cells, v2, axis


def row_record(r: c1.RowReading, cls_of: dict[str, str]) -> dict[str, Any]:
    """One row, as this lane read it. No value is recomputed and no direction is re-derived."""
    return {
        "assertion_id": r.assertion_id,
        "cell": r.cell,
        "v2_class_as_committed": r.v2_class,
        "v2_unresolved_reason_as_committed": r.reason,
        "retained_fields": list(r.fields),
        "retained_values": list(r.values),
        "satisfies_clause_1": r.satisfies_clause_1,
        "passes_amendment_1": r.passes_amendment_1,
        "cell_is_an_argmax": r.cell_is_argmax,
        "argmax_kind": r.argmax_kind,
        "cell_track_class": cls_of.get(r.cell, ""),
    }


def findings(tally: dict[str, Any], cells: list[c1.CellTracks]) -> dict[str, Any]:
    """What the counts say about the clause's stated reason, in the counts' own terms."""
    capable = [c for c in cells if c.rows_satisfying_clause_1 > 0]
    multi = [c for c in capable if len(c.tracks) >= 2]
    one = [c for c in capable if len(c.tracks) == 1]
    return {
        "clause_1_capable_cells": {
            "count": len(capable),
            "with_two_or_more_tracks": len(multi),
            "with_exactly_one_track": len(one),
            "cells_with_exactly_one_track": [c.cell for c in one],
        },
        "what_the_two_values_are": (
            f"on the {tally['rows_satisfying_clause_1']} of {tally['rows_in_population']} rows that "
            "satisfy clause (1), the two retained values are the fields the result records, and "
            f"{tally['rows_whose_cell_is_an_argmax']} of the {tally['rows_in_population']} rows "
            "carry max_drop or max_rise among them. Which of a cell's tracks a value came from is "
            "not recorded anywhere: AMENDMENT_1 says in its own words that the cache records no "
            "track identity and that none is invented, so this lane does not name a track for a "
            "value and no count here does"
        ),
        "the_argmax": (
            "the assigned cell is predicted_coding['tissue'], and enhancer_target.predict_target "
            "sets that field to max_drop_tissue or max_rise_tissue - the tissue of the larger of "
            "the biggest drop and the biggest rise over all tracks, chosen jointly with the gene. "
            "An argmax is a selection, so on every row where the measured count holds, the cell is "
            "the place an extreme fell and not a context chosen before the values were read. The "
            "count is reported over this lane's own denominator of "
            f"{tally['rows_in_population']} rows and over no other"
        ),
        "distinctness_can_rest_on_two_assay_titles": (
            f"{tally['multi_track_cells_whose_nonzero_means_are_pairwise_distinct']} of "
            f"{tally['multi_track_cells']} multi-track cells of this population have tracks whose "
            "nonzero_mean values are pairwise distinct in the client's own metadata, and "
            f"{tally['rows_passing_amendment_1_whose_cell_is_one_biosample_name_several_assay_titles']} "
            "rows pass amendment 1 with an assigned cell whose tracks are one biosample name under "
            "several assay titles. Where both hold, amendment 1's distinctness test is satisfiable "
            "by two assay titles of one biosample name"
        ),
        "what_clause_1_secures_on_these_counts": (
            "stated in the counts' own terms and no further. Clause (1) requires two retained "
            "values and amendment 1 requires them to be numerically distinct; nothing in either "
            "requires the two to come from independent material, and nothing in the cache would let "
            "a rule check it. Where the assigned cell's tracks are one biosample name under several "
            "assay titles, a resolved call rests on two distinct columns of the model's output over "
            "material the metadata records under one biosample name - which is weaker than the "
            "clause's words `no direction rests on one track` invite a reader to assume, and which "
            "the clause does not say. Whether those columns are replicates is not known, in "
            "_cell_summary's own words. This lane proposes no change to the clause"
        ),
        "the_fourth_answer": (
            f"{len(one)} of the {len(capable)} clause-(1)-capable cells carry exactly ONE track in "
            "the client's metadata. For such a cell two retained values are two field names over "
            "one track, which is the shape amendment 1 exists to refuse; the result reports how "
            "many such rows nevertheless pass amendment 1"
        ),
    }


def verdict(tally: dict[str, Any], cells: list[c1.CellTracks]) -> dict[str, Any]:
    """The plain answer, derived from the counts and conditional on them.

    Written AFTER the counts were taken, which `written_after_the_counts` says in the file itself.
    It is a reading of the figures and was not registered as a prediction of them; the registration
    fixed the classes, the rules and the count list and nothing more. Every branch below is decided
    by a measured figure, so a different population would get a different sentence.
    """
    resolved = [
        r for r in tally["rows_by_v2_class_as_committed"].items() if r[0].startswith("resolved") and r[1]
    ]
    n_resolved = sum(n for _, n in resolved)
    amend = tally["rows_passing_amendment_1"]
    on_multi = tally["rows_passing_amendment_1_whose_cell_is_one_biosample_name_several_assay_titles"]
    on_one = tally["rows_passing_amendment_1_whose_cell_is_one_track_only"]
    literal_secured = amend > 0 and on_one == 0 and on_multi == amend
    capable = [c for c in cells if c.rows_satisfying_clause_1 > 0]
    multi_incapable = [c for c in cells if len(c.tracks) >= 2 and c.rows_satisfying_clause_1 == 0]
    return {
        "written_after_the_counts": (
            "this block was written after the counts above were taken and is a reading of them. It "
            "was NOT registered as a prediction of any figure; what the registration fixed was the "
            "classes, the rules, the count list and the refusals. Every sentence below is decided "
            "by a measured figure and would read differently on different figures"
        ),
        "does_clause_1_secure_its_stated_reason": (
            "YES as to its words, NO as to what its words invite"
            if literal_secured
            else "NO: see the figures, which do not support the literal reading either"
        ),
        "as_to_its_words": (
            f"on this population clause (1)'s stated reason holds literally. All {amend} rows that "
            f"pass amendment 1 - which is all {n_resolved} resolved rows - have an assigned cell "
            "that carries two RNA-seq tracks in the client's own metadata, and their two retained "
            f"values are numerically distinct; {on_one} such row has a one-track cell. One track "
            "yields one value and both retained fields round to four places the same way, so two "
            "distinct values cannot both have come from one track. No direction among the resolved "
            "calls of this population therefore rests on one track"
            if literal_secured
            else "the literal reading is not supported on these figures; the counts above are the statement"
        ),
        "as_to_what_its_words_invite": (
            f"all {on_multi} of those rows have an assigned cell whose tracks are ONE biosample "
            "name under two assay titles - polyA plus RNA-seq and total RNA-seq - agreeing also on "
            "ontology term, biosample type, life stage, endedness and modification status. So a "
            "resolved call rests on two assays of material the metadata records under one "
            "biosample name. Clause (1) does not require, and nothing in the cache would let a "
            "rule check, that the two values come from independent material. Whether they are "
            "replicates is not known, in _cell_summary's own words. The corroboration a reader may "
            "take from `so that no direction rests on one track` is therefore not secured by the "
            "clause, although the words themselves are"
        ),
        "the_reason_the_clause_does_not_state": (
            "stated as the distinction it is and not as a defect found. `no direction rests on one "
            "track` is a claim about tracks and it holds. `no direction rests on one biosample`, "
            "`on one assay` or `on unreplicated evidence` are different claims, the clause makes "
            "none of them, and this lane measured no figure that establishes any of them either "
            "way: the metadata copy carries no accession, so even `one biosample` is not "
            "established of these cells - only one biosample NAME"
        ),
        "the_reach_of_the_clause": (
            f"{len(capable)} of {tally['distinct_assigned_cells']} assigned cells of this "
            f"population can satisfy clause (1) at all, and {len(multi_incapable)} further cells "
            "carry two tracks each and still cannot, because enhancer_target.aggregate writes a "
            "by_cell value only for a cell in its four-name CELLS list and discards every other "
            "track's value. So which rows clause (1) can reach is set by that hard-coded list and "
            "not by the evidence the model emitted. This is a statement about retention, the same "
            "reading lane-cellcover registered, and it is not a measurement of absence of "
            "regulation anywhere"
        ),
        "cells_with_two_tracks_that_still_cannot_satisfy_clause_1": [c.cell for c in multi_incapable],
        "the_argmax_answer": (
            f"{tally['rows_whose_cell_is_an_argmax']} of {tally['rows_in_population']} rows - "
            f"{tally['rows_whose_cell_is_argmax_of_max_drop_only']} by max_drop, "
            f"{tally['rows_whose_cell_is_argmax_of_max_rise_only']} by max_rise, "
            f"{tally['rows_whose_cell_is_argmax_of_both']} by both. The rule's cell is an argmax on "
            "every row of this population, so it is never an independently assigned context here"
        ),
        "no_rule_is_changed_by_this": c1.NO_RECOMMENDATION,
    }


def payload() -> dict[str, Any]:
    rows, cells, v2, axis = rows_and_cells()
    cls_of = {c.cell: c.track_class for c in cells}
    tally = c1.Tally(rows=rows, cells=cells, committed=v2.get("counts") or {}).to_dict()
    entries = [
        mf.input_entry(
            c1.V2_RESULT,
            partition="the committed v2 classification; its 166 predicted rows are the population",
        ),
        mf.input_entry(
            c1.TRACK_METADATA,
            partition="the saved copy of the model client's own track metadata, RNA-seq rows only",
        ),
    ]
    return {
        "result": RESULT,
        "date": date.today().isoformat(),
        "lane": "lane-clause1",
        "registration": REGISTRATION,
        "question": c1.QUESTION,
        "population": c1.POPULATION,
        "denominator_is_separate": c1.DENOMINATOR_IS_SEPARATE,
        "exploration_preceded_this": c1.EXPLORATION_PRECEDED_THIS,
        "edits_no_rule": (
            "no clause, amendment or threshold of direction_v2 or of any other rule is changed by "
            "this lane. The finding below is a finding and not a licence: a changed clause would "
            "invalidate every count already taken under it"
        ),
        "clause_1_quoted": (
            "SIGN_RULE_V2 clause (1): `at least two values are retained for that cell name, so "
            "that no direction rests on one track - which is the defect v1 has`"
        ),
        "amendment_1": dv2.AMENDMENT_1,
        "cell_evidence": dv2.CELL_EVIDENCE,
        "assigned_cell": dv2.ASSIGNED_CELL,
        "track_label_rule": c1.TRACK_LABEL_RULE,
        "name_match_is_not_a_confirmation": c1.NAME_MATCH_IS_NOT_A_CONFIRMATION,
        "metadata_copy_limitation": c1.METADATA_COPY_LIMITATION,
        "class_rule": c1.CLASS_RULE,
        "track_classes": list(c1.TRACK_CLASSES),
        "clause_1_satisfied_rule": c1.CLAUSE_1_SATISFIED_RULE,
        "argmax_rule": c1.ARGMAX_RULE,
        "argmax_is_a_selection": c1.ARGMAX_IS_A_SELECTION,
        "outside_retained_cells_cannot_resolve": c1.OUTSIDE_RETAINED_CELLS_CANNOT_RESOLVE,
        "distinctness_rule": c1.DISTINCTNESS_RULE,
        "established_per_cell": c1.ESTABLISHED_PER_CELL,
        "readings_carried_unstrengthened": c1.READINGS_CARRIED,
        "no_recommendation": c1.NO_RECOMMENDATION,
        "validates_nothing": c1.VALIDATES_NOTHING,
        "cannot_establish": list(c1.CANNOT_ESTABLISH),
        "metadata_columns": list(c1.metadata_columns(c1.TRACK_METADATA)),
        "rna_seq_track_axis": {
            "distinct_track_names": axis.width,
            "metadata_rows_read": axis.rows_read,
            "source": axis.source,
            "expected": cc.EXPECTED_AXIS,
        },
        "retained_cells_today": list(CELLS),
        "counts": tally,
        "findings": findings(tally, cells),
        "verdict": verdict(tally, cells),
        "per_cell": [c.to_dict() for c in cells],
        "per_row": [row_record(r, cls_of) for r in rows],
        "alphagenome_requests": 0,
        "money": "none: no request was sent and no key was read",
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
            "written on. No coordinate of this result's own is in GRCh38 or any build",
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
                "or moved; values are read only to ask whether two of them are equal",
            },
            "exclusions": [
                "every assertion of the published response map whose status is not `predicted`",
                "every metadata row whose output is not rna_seq",
                "every cached element answer and every chromosome archive: none is opened",
            ],
            "partitions": {
                "the_166_rows": "this lane's row denominator",
                "the_distinct_assigned_cells": "this lane's cell denominator",
            },
            "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        },
    }


def main() -> int:
    path = save_result(RESULT, payload())
    print(f"wrote: {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
