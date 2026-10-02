# SPDX-License-Identifier: AGPL-3.0-or-later
"""Does direction rule v2's clause (1) secure its own stated reason? The classes and the counts.

**This file is the registration.** Every class, every rule by which a cell is classified and every
count the run will report is written here, and the first commit of this module carries nothing but
these constants and the pure functions that apply them: no count is taken, no file is written and
no figure of this lane's own exists yet.

**This lane edits no rule.** `attribution/direction_v2.py` is read and imported; not one of its
constants, clauses, amendments or thresholds is changed, and nothing here may be read as a proposal
to change one. Changing clause (1) would invalidate every count already taken under it, so the
answer to the question - either way - is a finding and not a licence.

What clause (1) says, and what it gives as its reason
------------------------------------------------------
`direction_v2.SIGN_RULE_V2` clause (1): *"at least two values are retained for that cell name, so
that no direction rests on one track - which is the defect v1 has"*. The words after the comma are
the clause's **stated reason**, and this lane asks whether the clause secures it.

`direction_v2.AMENDMENT_1` then required the two retained values to be **numerically distinct**,
because on a legacy row `by_cell[cell]` can be the very track the maximum selected. The amendment
can only move a call from resolved to unresolved; it moves no threshold.

Where a second value can come from at all
-----------------------------------------
`direction_v2.retained_values` retains at most three per gene row: `max_drop` when
`max_drop_tissue` is the cell, `max_rise` when `max_rise_tissue` is the cell, and `by_cell[cell]`
when the cell is one of `predict.enhancer_target.CELLS`. `OUTSIDE_RETAINED_CELLS_CANNOT_RESOLVE`
records what follows from that for a cell outside `CELLS`, as arithmetic and before any count.
"""

from __future__ import annotations

import csv
import math
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.attribution import cellcover as cc
from genomeos.attribution import direction_v2 as dv2
from genomeos.predict.enhancer_target import CELLS

# ---- the question, and the population this lane will not mix with another ----------------------

QUESTION = (
    "does direction rule v2's clause (1) secure its own stated reason? The clause requires two "
    "retained values for the rule's assigned cell `so that no direction rests on one track`, and "
    "amendment 1 requires the two to be numerically distinct. This lane counts three things and "
    "argues none of them: (a) for every cell that can satisfy clause (1) in this population, "
    "whether the tracks carrying its label are one biosample name under several assay titles, "
    "several biosample names, or undeterminable from the client's own saved track metadata; (b) how "
    "often the rule's assigned cell is an argmax over the track axis rather than a context chosen "
    "independently of the values read; (c) whether two assay titles of one biosample name differ "
    "numerically, so that amendment 1's distinctness test can be satisfied by them"
)

POPULATION = (
    "the 166 assertion rows of the committed data/results/respmap_direction_v2.json - every "
    "assertion of the published increment-2 response map whose status is `predicted`, as that "
    "committed result records them. The rows are READ and never re-derived: no element answer is "
    "opened, no archive is loaded and no direction is recomputed. The cell denominator is the "
    "distinct assigned cells among those 166 rows"
)

DENOMINATOR_IS_SEPARATE = (
    "this lane's counts are not added to lane-propb's 112, 93 or 19, nor to the 26, 67, 14 or 10 of "
    "other lanes. Several of those overlap and may not be summed. Every count below names the "
    "denominator it is taken over - 166 rows, or the distinct cells of those rows - and the two "
    "reason counts this lane reads from the committed v2 result (93 and 19) are quoted as that "
    "result's figures and are not re-counted as this lane's own"
)

EXPLORATION_PRECEDED_THIS = (
    "stated against this lane's own interest. Before writing this file, this lane read the "
    "committed data/results/respmap_direction_v2.json header fields and one sample row, the "
    "committed data/results/cellcover.json prose fields, the column header of the saved track "
    "metadata copy, and the source of direction_v2.retained_values, direction_v2.SIGN_RULE_V2, "
    "direction_v2.AMENDMENT_1, enhancer_target.aggregate, enhancer_target.predict_target, "
    "enhancer_target._cell_summary and cellcover.rna_axis. It had therefore already seen these "
    "figures: 166 rows, 54 resolved, 0 opposite, 112 unresolved, 93 one_value_only, 19 "
    "one_track_seen_twice, 73 `cells_that_are_cached_tracks`, 371 distinct RNA-seq track names, "
    "lane-cellcover's 59 cells split 44 one-track and 15 two-track, and amendment 1's chr21 figure "
    "of 271 of 619. So this registration is NOT a prediction of any figure and no reader may treat "
    "it as one. What it binds is the class set, the rule by which each class is established, the "
    "count list and the refusals, so that the committed program's output can be held against "
    "written rules rather than against this lane's recollection. The separate commit this file gets "
    "before the run audits the order of the RULE, not the order of the figures"
)

# ---- what a cell's label is, and what the metadata copy can and cannot establish ---------------

TRACK_LABEL_RULE = (
    "a track's label is cellcover.track_label, imported and not restated: the gtex_tissue when the "
    "metadata row carries one that is not empty and not the string `nan`, else the biosample_name. "
    "This is predict.alphagenome_adapter.tissue_names' own rule, so a cell label on a compiled rule "
    "is a track label by construction. A cell is matched to tracks by EQUALITY of that label, "
    "character for character. No ontology lookup, no synonym table and no judgement resolves a row"
)

NAME_MATCH_IS_NOT_A_CONFIRMATION = (
    "the saved metadata copy carries no biosample accession and no experiment accession - its "
    "columns are name, strand, Assay title, ontology_curie, biosample_name, biosample_type, "
    "biosample_life_stage, gtex_tissue, data_source, endedness, genetically_modified, nonzero_mean, "
    "output, histone_mark - so two rows that agree on biosample_name CANNOT be established from it "
    "to be one physical specimen. They could be two specimens of one cell line or two donors of one "
    "tissue. This lane therefore never writes `same biosample`: the class is "
    "`one_biosample_name_several_assay_titles`, the per-cell record names every field that was "
    "compared, and `specimen_identity_established` is false on every multi-track cell. The run "
    "REFUSES if the copy ever gains a column whose name contains `accession`, because then this "
    "limitation would no longer be the honest one"
)

METADATA_COPY_LIMITATION = cc.METADATA_COPY_LIMITATION

# ---- the classes, fixed here -------------------------------------------------------------------

#: No RNA-seq track of the saved copy carries the cell's label. Clause (1) can rest on no track of
#: this cell at all, and why the label is absent is not established here.
ABSENT = "absent_from_the_metadata_copy"
#: Exactly one RNA-seq track carries the label. Two retained values for such a cell are two FIELD
#: NAMES over one track, which is the shape amendment 1 exists to refuse.
ONE_TRACK = "one_track_only"
#: Two or more tracks, all agreeing on biosample_name, carrying two or more distinct `Assay title`
#: values. The tracks are distinct columns of the model's output; whether they are one specimen is
#: NOT established (NAME_MATCH_IS_NOT_A_CONFIRMATION), and whether they are replicates is not
#: known (REPLICATION_READING_CARRIED).
ONE_BIOSAMPLE_NAME_SEVERAL_ASSAY_TITLES = "one_biosample_name_several_assay_titles"
#: Two or more tracks, all agreeing on biosample_name, all carrying one `Assay title`. They differ
#: in some other recorded column, which the per-cell record names.
ONE_BIOSAMPLE_NAME_ONE_ASSAY_TITLE = "one_biosample_name_one_assay_title"
#: Two or more tracks whose biosample_name values differ. The distinct names are recorded.
SEVERAL_BIOSAMPLE_NAMES = "several_biosample_names"

#: Every class, in the order the result reports them. A cell falls in exactly one.
TRACK_CLASSES = (
    ABSENT,
    ONE_TRACK,
    ONE_BIOSAMPLE_NAME_SEVERAL_ASSAY_TITLES,
    ONE_BIOSAMPLE_NAME_ONE_ASSAY_TITLE,
    SEVERAL_BIOSAMPLE_NAMES,
)

CLASS_RULE = (
    "a cell's class is decided by the saved copy's RNA-seq rows whose label equals the cell label, "
    "and by nothing else. 0 such rows is `absent_from_the_metadata_copy`; exactly 1 is "
    "`one_track_only`; 2 or more with one distinct biosample_name and 2 or more distinct `Assay "
    "title` values is `one_biosample_name_several_assay_titles`; 2 or more with one biosample_name "
    "and one `Assay title` is `one_biosample_name_one_assay_title`; 2 or more with 2 or more "
    "distinct biosample_names is `several_biosample_names`. The classes are exhaustive and "
    "mutually exclusive by construction and a test asserts both. `one_track_only` is a FOURTH "
    "answer that the three the question offered - same biosample, different biosamples, "
    "undeterminable - does not contain, and it is reported as its own class rather than folded "
    "into any of them"
)

PER_CELL_FIELDS = (
    "cell",
    "rows_in_population",
    "rna_seq_tracks",
    "track_class",
    "established_how",
    "established",
    "specimen_identity_established",
    "track_names",
    "assay_titles",
    "biosample_names",
    "ontology_curies",
    "biosample_types",
    "biosample_life_stages",
    "endedness",
    "genetically_modified",
    "nonzero_means",
    "nonzero_means_pairwise_distinct",
    "in_retained_cells",
    "rows_satisfying_clause_1",
    "rows_resolved",
    "rows_one_track_seen_twice",
)

ESTABLISHED_PER_CELL = (
    "every per-cell row carries `established_how` - `client_track_metadata` when at least one "
    "RNA-seq row of the copy carries the label by equality, `absent` when none does - and a prose "
    "`established` naming the fields that were compared on that cell's tracks and what the "
    "comparison did and did not settle. The difference between a label match and a specimen "
    "identity is visible per row, in `specimen_identity_established`, which is false on every "
    "multi-track cell for the reason NAME_MATCH_IS_NOT_A_CONFIRMATION gives"
)

# ---- clause (1) and the argmax, as this lane will measure them ---------------------------------

CLAUSE_1_SATISFIED_RULE = (
    "a ROW satisfies clause (1) when the committed result records two or more entries in its "
    "`retained_values`, which is exactly the test direction_v2._v2 applies (`len(values) < 2` is "
    "`one_value_only`). A CELL can satisfy clause (1) in this population when at least one row "
    "assigned to it satisfies it. Satisfying clause (1) is not resolving: amendment 1, the zero "
    "test, the sign test and the magnitude floor all come after it, and a row can satisfy clause "
    "(1) and still be unresolved"
)

ARGMAX_RULE = (
    "the rule's assigned cell is an ARGMAX on a row when the committed row's `retained_values` "
    "carries the field `max_drop` or the field `max_rise`, because direction_v2.retained_values "
    "retains `max_drop` exactly when the gene row's `max_drop_tissue` equals the cell and "
    "`max_rise` exactly when `max_rise_tissue` equals it. Those two fields are the tissue of the "
    "minimum and of the maximum over the whole track axis, as enhancer_target.aggregate writes "
    "them. The measurement is taken from the committed rows; it is not re-derived from any element "
    "answer, and nothing is loaded"
)

ARGMAX_IS_A_SELECTION = (
    "stated before any count, because it follows from the code and not from the data. "
    "direction_v2.ASSIGNED_CELL says the assigned cell is `predicted_coding['tissue']`, and "
    "enhancer_target.predict_target sets that field to `max_drop_tissue` or `max_rise_tissue` - "
    "the tissue of the larger of the biggest drop and the biggest rise over all tracks, jointly "
    "with the gene. So on every rule the compiler wrote from predict_target, the cell is the place "
    "an extreme fell and not a context chosen before the values were read. The count below tests "
    "whether the committed rows agree with that reading of the code; a figure below 166 would mean "
    "some row's cell is not an argmax and the code's own account is incomplete"
)

OUTSIDE_RETAINED_CELLS_CANNOT_RESOLVE = (
    "arithmetic, registered before the count. For a cell outside predict.enhancer_target.CELLS no "
    "`by_cell` value is ever written (aggregate writes one only `if tissue in cells`), so such a "
    "cell can reach two retained values only by being BOTH `max_drop_tissue` and `max_rise_tissue`. "
    "aggregate starts min and max at 0.0 and lowers min only on a negative value and raises max "
    "only on a positive one, so both fields naming one cell means one retained value is below zero "
    "and the other above it - which direction_v2 records as `signs_disagree_in_cell`. Therefore no "
    "cell outside CELLS can ever RESOLVE under v2, however many tracks it has, and every resolved "
    "call's cell is one of "
    f"{list(CELLS)}. The run reports the measured count of resolved rows whose cell is outside "
    "CELLS and REFUSES if it is not 0"
)

DISTINCTNESS_RULE = (
    "whether two assay titles of one biosample name differ numerically is established two ways and "
    "both are reported. (a) From the client's own metadata: the `nonzero_mean` column is a number "
    "per track, and the run reports, per multi-track cell, whether its tracks' nonzero_mean values "
    "are pairwise distinct. Distinct nonzero_mean values establish that the tracks are different "
    "columns of the model's output, not one column seen twice. (b) From the committed rows: among "
    "rows whose assigned cell is a multi-track cell, how many carry two numerically distinct "
    "retained values - that is, pass amendment 1 - and how many are recorded "
    "`one_track_seen_twice`. The two together say whether amendment 1's distinctness test CAN be "
    "satisfied by two assay titles of one biosample name. Neither is a measurement of replication "
    "and neither establishes a specimen identity"
)

# ---- the readings this lane carries forward, unstrengthened -------------------------------------

READINGS_CARRIED = {
    "withholds_never_reverses": (
        "under SIGN_RULE_V2 v2 can WITHHOLD a direction and can never reverse one: v1's own "
        "selected value is itself one of the assigned cell's retained values and v2 resolves only "
        "when every retained value shares one sign, so a resolved v2 call always carries v1's sign. "
        "Carried word for word from the committed v2 result and not strengthened here"
    ),
    "unresolved_is_not_absence": dv2.UNRESOLVED_MEANS,
    "not_known_to_be_replication": (
        "enhancer_target._cell_summary's own words about a cell's several emitted values: `the "
        "tracks are not known to be biological replicates, so it is not a validated cell-level "
        "effect`. Nothing in this lane's counts makes them replicates, and no count here is "
        "evidence that they are or are not"
    ),
    "cell2_group_not_biological_independence": (
        "attribution.cell2's grouping is an operational grouping, NOT established biological "
        "independence. It is neither used nor needed here: every cell is counted on its own label "
        "and no cell is merged with another"
    ),
    "cache_is_legacy": dv2.CACHE_IS_LEGACY,
}

REPLICATION_READING_CARRIED = READINGS_CARRIED["not_known_to_be_replication"]

NO_RECOMMENDATION = (
    "this lane makes NO recommendation about clause (1), amendment 1 or any threshold, and proposes "
    "no change to any of them. It reports what the clause does and does not secure and stops there. "
    "No published number is touched, no direction is recomputed, no result of another lane is "
    "rewritten, and no count here supersedes one"
)

VALIDATES_NOTHING = (
    "no count here says a v2 call is right or wrong. Nothing here is a validation of v1 or of v2, "
    "nothing here is an accuracy figure, and the three classes are a description of the model's "
    "track axis as the saved metadata records it - not a measurement of biology"
)

CANNOT_ESTABLISH = (
    "whether two tracks that agree on biosample_name are one physical specimen: the saved metadata "
    "copy carries no biosample or experiment accession (NAME_MATCH_IS_NOT_A_CONFIRMATION)",
    "whether any two tracks of one cell label are biological replicates: _cell_summary says in its "
    "own words that they are not known to be",
    "whether the saved track metadata copy is the client's current output_metadata: re-fetching it "
    "is a request this lane may not make (METADATA_COPY_LIMITATION)",
    "which of a cell's tracks a given retained value came from: the cache records no track "
    "identity, and AMENDMENT_1 says in its own words that none is invented",
    "whether a different assigned cell would resolve any row: that needs values this lane does not "
    "have and requests it may not make",
    "whether clause (1) would be better or worse if it were written differently: this lane edits no "
    "rule and tests no alternative",
)

COUNTS_NAMED = (
    "rows_in_population",
    "rows_satisfying_clause_1",
    "rows_not_satisfying_clause_1",
    "rows_by_v2_class_as_committed",
    "rows_by_unresolved_reason_as_committed",
    "distinct_assigned_cells",
    "clause_1_capable_cells",
    "cells_by_track_class",
    "clause_1_capable_cells_by_track_class",
    "rows_whose_cell_is_an_argmax",
    "rows_whose_cell_is_argmax_of_max_drop_only",
    "rows_whose_cell_is_argmax_of_max_rise_only",
    "rows_whose_cell_is_argmax_of_both",
    "rows_whose_cell_is_not_an_argmax",
    "resolved_rows_whose_cell_is_outside_retained_cells",
    "multi_track_cells_whose_nonzero_means_are_pairwise_distinct",
    "rows_passing_amendment_1",
    "rows_passing_amendment_1_whose_cell_is_one_biosample_name_several_assay_titles",
    "rows_passing_amendment_1_whose_cell_is_one_track_only",
    "rows_recorded_one_track_seen_twice_by_track_class",
)

REFUSALS = (
    "the program refuses rather than reporting if the committed v2 result's counts.denominator is "
    "not 166, or its by_class and unresolved_by_reason counts are not the committed 54 / 0 / 112 "
    "and 93 / 19: the population is lane-propb's and is read, never re-derived",
    "the program refuses if the committed result's assertion list is not 166 rows long, or if any "
    "row lacks `cell`, `retained_values` or `v2_class`",
    "the program refuses if the saved metadata copy's distinct RNA-seq track-name count is not "
    f"{cc.EXPECTED_AXIS}: every per-cell class is grounded on that axis and must not be reported "
    "without it",
    "the program refuses if the saved metadata copy carries any column whose name contains "
    "`accession`, because NAME_MATCH_IS_NOT_A_CONFIRMATION would then be the wrong limitation and "
    "the `undeterminable` reading would be understated",
    "the program refuses if any row's assigned cell is the empty string or the compiler's "
    "context-unknown token: such a row names no cell to classify",
    "the program refuses if the measured count of resolved rows whose cell is outside "
    "predict.enhancer_target.CELLS is not 0, because OUTSIDE_RETAINED_CELLS_CANNOT_RESOLVE says it "
    "is 0 by arithmetic and a non-zero count would mean this lane has misread the code",
    "the program refuses if the track classes are not exhaustive and mutually exclusive over the "
    "cells it classified",
)

V2_RESULT = cc.V2_RESULT
TRACK_METADATA = cc.TRACK_METADATA

#: The committed figures of the v2 result that the run checks before it counts anything.
COMMITTED_V2_FIGURES = {
    "denominator": 166,
    "by_class": {
        "resolved_agrees_with_published": 54,
        "resolved_opposite_to_published": 0,
        "unresolved": 112,
    },
    "unresolved_by_reason": {"one_value_only": 93, "one_track_seen_twice": 19},
}

#: A column name fragment whose presence would make NAME_MATCH_IS_NOT_A_CONFIRMATION wrong.
FORBIDDEN_COLUMN_FRAGMENT = "accession"

#: The compiler's token for a rule whose cell is unknown; such a row names no cell to classify.
CONTEXT_UNKNOWN_TOKENS = ("", "CONTEXT_UNKNOWN", "context_unknown")


# ---- the pure functions that apply the rules above ----------------------------------------------


@dataclass(frozen=True)
class CellTracks:
    """One assigned cell, its tracks in the saved metadata copy, and its class."""

    cell: str
    rows_in_population: int
    tracks: tuple[dict[str, Any], ...]
    rows_satisfying_clause_1: int = 0
    rows_resolved: int = 0
    rows_one_track_seen_twice: int = 0

    def _distinct(self, column: str) -> tuple[str, ...]:
        return tuple(sorted({str(t.get(column) or "") for t in self.tracks}))

    @property
    def assay_titles(self) -> tuple[str, ...]:
        return self._distinct("Assay title")

    @property
    def biosample_names(self) -> tuple[str, ...]:
        return self._distinct("biosample_name")

    @property
    def track_class(self) -> str:
        """`CLASS_RULE`, word for word."""
        if not self.tracks:
            return ABSENT
        if len(self.tracks) == 1:
            return ONE_TRACK
        if len(self.biosample_names) > 1:
            return SEVERAL_BIOSAMPLE_NAMES
        if len(self.assay_titles) > 1:
            return ONE_BIOSAMPLE_NAME_SEVERAL_ASSAY_TITLES
        return ONE_BIOSAMPLE_NAME_ONE_ASSAY_TITLE

    @property
    def nonzero_means(self) -> tuple[float | None, ...]:
        out: list[float | None] = []
        for t in self.tracks:
            raw = t.get("nonzero_mean")
            try:
                val = float(raw)  # type: ignore[arg-type]
            except (TypeError, ValueError):
                out.append(None)
                continue
            # `nan` is how the copy records a value it does not carry, and nan != nan, so a
            # non-finite reading is recorded as unreadable rather than as a distinct number.
            out.append(round(val, 10) if math.isfinite(val) else None)
        return tuple(out)

    @property
    def nonzero_means_pairwise_distinct(self) -> bool | None:
        """True when every track's nonzero_mean differs, None when fewer than two tracks or any is
        unreadable - never asserted when a value could not be read."""
        vals = self.nonzero_means
        if len(vals) < 2 or any(v is None for v in vals):
            return None
        return len(set(vals)) == len(vals)

    @property
    def specimen_identity_established(self) -> bool:
        """False on every multi-track cell: NAME_MATCH_IS_NOT_A_CONFIRMATION."""
        return False

    @property
    def established(self) -> str:
        """What the comparison on this cell's tracks settled, and what it did not."""
        if not self.tracks:
            return (
                "no RNA-seq row of the saved metadata copy carries this label by equality, so no "
                "track of this cell was found; why the label is absent is not established here"
            )
        fields = (
            "Assay title, biosample_name, ontology_curie, biosample_type, biosample_life_stage, "
            "endedness, genetically_modified and nonzero_mean"
        )
        if len(self.tracks) == 1:
            return (
                "exactly one RNA-seq row of the copy carries this label by equality; its "
                f"{fields} are recorded. Two retained values for this cell are two field names over "
                "one track, not two tracks"
            )
        return (
            f"{len(self.tracks)} RNA-seq rows of the copy carry this label by equality, compared on "
            f"{fields}. The distinct values of each are recorded. The copy carries no biosample or "
            "experiment accession, so agreement on biosample_name does NOT establish one physical "
            "specimen and `same biosample` is not claimed; whether the tracks are replicates is not "
            "known (_cell_summary's own words)"
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "cell": self.cell,
            "rows_in_population": self.rows_in_population,
            "rna_seq_tracks": len(self.tracks),
            "track_class": self.track_class,
            "established_how": "client_track_metadata" if self.tracks else "absent",
            "established": self.established,
            "specimen_identity_established": self.specimen_identity_established,
            "track_names": [str(t.get("name") or "") for t in self.tracks],
            "assay_titles": list(self.assay_titles),
            "biosample_names": list(self.biosample_names),
            "ontology_curies": list(self._distinct("ontology_curie")),
            "biosample_types": list(self._distinct("biosample_type")),
            "biosample_life_stages": list(self._distinct("biosample_life_stage")),
            "endedness": list(self._distinct("endedness")),
            "genetically_modified": list(self._distinct("genetically_modified")),
            "nonzero_means": list(self.nonzero_means),
            "nonzero_means_pairwise_distinct": self.nonzero_means_pairwise_distinct,
            "in_retained_cells": self.cell in CELLS,
            "rows_satisfying_clause_1": self.rows_satisfying_clause_1,
            "rows_resolved": self.rows_resolved,
            "rows_one_track_seen_twice": self.rows_one_track_seen_twice,
        }


@dataclass
class RowReading:
    """One committed assertion row, read for the three questions and for nothing else."""

    assertion_id: str
    cell: str
    v2_class: str
    reason: str | None
    fields: tuple[str, ...]
    values: tuple[float, ...]

    @property
    def satisfies_clause_1(self) -> bool:
        """`CLAUSE_1_SATISFIED_RULE`: two or more retained values."""
        return len(self.values) >= 2

    @property
    def passes_amendment_1(self) -> bool:
        """Two or more retained values, numerically distinct - amendment 1's own test."""
        return self.satisfies_clause_1 and len(set(self.values)) >= 2

    @property
    def cell_is_argmax(self) -> bool:
        """`ARGMAX_RULE`: a `max_drop` or `max_rise` field was retained for this cell."""
        return "max_drop" in self.fields or "max_rise" in self.fields

    @property
    def argmax_kind(self) -> str:
        drop, rise = "max_drop" in self.fields, "max_rise" in self.fields
        if drop and rise:
            return "both"
        if drop:
            return "max_drop_only"
        if rise:
            return "max_rise_only"
        return "none"


def read_rows(result: dict[str, Any]) -> list[RowReading]:
    """The committed v2 result's assertion rows, read and not re-derived."""
    rows: list[RowReading] = []
    for a in result["assertions"]:
        for key in ("cell", "retained_values", "v2_class"):
            if key not in a:
                raise ValueError(f"row {a.get('assertion_id')!r} lacks {key!r} (REFUSALS)")
        cell = str(a["cell"] or "")
        if cell.strip() in CONTEXT_UNKNOWN_TOKENS:
            raise ValueError(f"row {a.get('assertion_id')!r} names no cell to classify (REFUSALS)")
        retained = tuple(a["retained_values"] or ())
        rows.append(
            RowReading(
                assertion_id=str(a.get("assertion_id") or ""),
                cell=cell,
                v2_class=str(a["v2_class"]),
                reason=a.get("v2_unresolved_reason"),
                fields=tuple(str(r["field"]) for r in retained),
                values=tuple(float(r["value"]) for r in retained),
            )
        )
    return rows


def check_population(result: dict[str, Any], rows: list[RowReading]) -> None:
    """`REFUSALS`, on the committed v2 result."""
    want = COMMITTED_V2_FIGURES
    counts = result.get("counts") or {}
    if counts.get("denominator") != want["denominator"]:
        raise ValueError(f"committed v2 denominator is not {want['denominator']} (REFUSALS)")
    if counts.get("by_class") != want["by_class"]:
        raise ValueError("committed v2 by_class counts differ from the registered ones (REFUSALS)")
    if counts.get("unresolved_by_reason") != want["unresolved_by_reason"]:
        raise ValueError("committed v2 unresolved_by_reason differs from the registered (REFUSALS)")
    if len(rows) != want["denominator"]:
        raise ValueError(f"the committed assertion list is not {want['denominator']} rows (REFUSALS)")


def metadata_columns(path: Path = TRACK_METADATA) -> tuple[str, ...]:
    """The saved copy's column names, read so the accession refusal is checked and not asserted."""
    with path.open(newline="") as fh:
        header = next(csv.reader(fh))
    return tuple(header)


def check_metadata(columns: tuple[str, ...], axis_width: int) -> None:
    """`REFUSALS`, on the saved metadata copy."""
    bad = [c for c in columns if FORBIDDEN_COLUMN_FRAGMENT in c.lower()]
    if bad:
        raise ValueError(
            f"the saved metadata copy carries {bad!r}: NAME_MATCH_IS_NOT_A_CONFIRMATION would be "
            "the wrong limitation, so this lane refuses to report the undeterminable reading "
            "(REFUSALS)"
        )
    if axis_width != cc.EXPECTED_AXIS:
        raise ValueError(f"the RNA-seq track axis is {axis_width}, not {cc.EXPECTED_AXIS} (REFUSALS)")


def classify(rows: list[RowReading], by_label: dict[str, list[dict[str, Any]]]) -> list[CellTracks]:
    """One `CellTracks` per distinct assigned cell of the population, ordered by rows then name."""
    per_cell: dict[str, list[RowReading]] = {}
    for r in rows:
        per_cell.setdefault(r.cell, []).append(r)
    out = [
        CellTracks(
            cell=cell,
            rows_in_population=len(rs),
            tracks=tuple(by_label.get(cell, ())),
            rows_satisfying_clause_1=sum(1 for r in rs if r.satisfies_clause_1),
            rows_resolved=sum(1 for r in rs if r.v2_class.startswith("resolved")),
            rows_one_track_seen_twice=sum(1 for r in rs if r.reason == "one_track_seen_twice"),
        )
        for cell, rs in per_cell.items()
    ]
    out.sort(key=lambda c: (-c.rows_in_population, c.cell))
    return out


@dataclass
class Tally:
    """Every count named in `COUNTS_NAMED`, and nothing else."""

    rows: list[RowReading]
    cells: list[CellTracks]
    committed: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        rows, cells = self.rows, self.cells
        cls_of = {c.cell: c.track_class for c in cells}
        multi = [c for c in cells if len(c.tracks) >= 2]
        amend = [r for r in rows if r.passes_amendment_1]
        twice = [r for r in rows if r.reason == "one_track_seen_twice"]
        capable = [c for c in cells if c.rows_satisfying_clause_1 > 0]
        return {
            "rows_in_population": len(rows),
            "rows_satisfying_clause_1": sum(1 for r in rows if r.satisfies_clause_1),
            "rows_not_satisfying_clause_1": sum(1 for r in rows if not r.satisfies_clause_1),
            "rows_by_v2_class_as_committed": self.committed.get("by_class"),
            "rows_by_unresolved_reason_as_committed": self.committed.get("unresolved_by_reason"),
            "distinct_assigned_cells": len(cells),
            "clause_1_capable_cells": len(capable),
            "cells_by_track_class": {k: sum(1 for c in cells if c.track_class == k) for k in TRACK_CLASSES},
            "clause_1_capable_cells_by_track_class": {
                k: sum(1 for c in capable if c.track_class == k) for k in TRACK_CLASSES
            },
            "rows_whose_cell_is_an_argmax": sum(1 for r in rows if r.cell_is_argmax),
            "rows_whose_cell_is_argmax_of_max_drop_only": sum(
                1 for r in rows if r.argmax_kind == "max_drop_only"
            ),
            "rows_whose_cell_is_argmax_of_max_rise_only": sum(
                1 for r in rows if r.argmax_kind == "max_rise_only"
            ),
            "rows_whose_cell_is_argmax_of_both": sum(1 for r in rows if r.argmax_kind == "both"),
            "rows_whose_cell_is_not_an_argmax": sum(1 for r in rows if not r.cell_is_argmax),
            "resolved_rows_whose_cell_is_outside_retained_cells": sum(
                1 for r in rows if r.v2_class.startswith("resolved") and r.cell not in CELLS
            ),
            "multi_track_cells": len(multi),
            "multi_track_cells_whose_nonzero_means_are_pairwise_distinct": sum(
                1 for c in multi if c.nonzero_means_pairwise_distinct is True
            ),
            "rows_passing_amendment_1": len(amend),
            "rows_passing_amendment_1_whose_cell_is_one_biosample_name_several_assay_titles": sum(
                1 for r in amend if cls_of.get(r.cell) == ONE_BIOSAMPLE_NAME_SEVERAL_ASSAY_TITLES
            ),
            "rows_passing_amendment_1_whose_cell_is_one_track_only": sum(
                1 for r in amend if cls_of.get(r.cell) == ONE_TRACK
            ),
            "rows_recorded_one_track_seen_twice_by_track_class": dict(
                Counter(cls_of.get(r.cell, "") for r in twice)
            ),
        }
