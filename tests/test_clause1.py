# SPDX-License-Identifier: AGPL-3.0-or-later
"""The clause-(1) classification: the classes, the rules that decide them, and the refusals.

Every test here is on `genomeos.attribution.clause1`'s own functions and on stub rows. No committed
result is recomputed, no element answer is opened and no archive loader is called.
"""

from __future__ import annotations

import csv
from pathlib import Path

import pytest

from genomeos.attribution import clause1 as c1
from genomeos.attribution import direction_v2 as dv2
from genomeos.predict.enhancer_target import CELLS


def track(name: str, assay: str, biosample: str, **extra: object) -> dict[str, object]:
    row: dict[str, object] = {
        "name": name,
        "Assay title": assay,
        "biosample_name": biosample,
        "ontology_curie": "CL:0000000",
        "biosample_type": "cell line",
        "biosample_life_stage": "adult",
        "endedness": "paired",
        "genetically_modified": "False",
        "nonzero_mean": "0.5",
        "gtex_tissue": "",
        "output": "rna_seq",
    }
    row.update(extra)
    return row


def row(cell: str, fields: tuple[str, ...], values: tuple[float, ...], **extra: object) -> dict[str, object]:
    out: dict[str, object] = {
        "assertion_id": f"r|{cell}|{'+'.join(fields)}",
        "cell": cell,
        "v2_class": "unresolved",
        "v2_unresolved_reason": None,
        "retained_values": [{"field": f, "value": v} for f, v in zip(fields, values, strict=True)],
    }
    out.update(extra)
    return out


# ---- this lane edits no rule -------------------------------------------------------------------


def test_clause_1_and_amendment_1_are_quoted_not_restated():
    """The clause's own words and its reason, as direction_v2 holds them."""
    assert "at least two values are retained for that cell name" in dv2.SIGN_RULE_V2
    assert "so that no direction rests on one track" in dv2.SIGN_RULE_V2
    assert "numerically distinct" in dv2.AMENDMENT_1


def test_the_magnitude_floor_and_the_default_are_untouched():
    """Nothing this lane imports has moved a threshold or the default version."""
    from genomeos.predict.enhancer_target import MIN_EFFECT

    assert dv2.MAGNITUDE_FLOOR is MIN_EFFECT
    assert dv2.DIRECTION_RULE_DEFAULT == dv2.V1


def test_no_recommendation_and_validates_nothing_are_registered():
    assert "NO recommendation" in c1.NO_RECOMMENDATION
    assert "no count here says a v2 call is right or wrong" in c1.VALIDATES_NOTHING


def test_the_replication_reading_is_carried_word_for_word():
    """_cell_summary's own wording, not a stronger one."""
    import inspect

    from genomeos.predict import enhancer_target as et

    doc = " ".join((inspect.getdoc(et._cell_summary) or "").split())
    quoted = "the tracks are not known to be biological replicates"
    assert quoted in doc
    assert quoted in " ".join(c1.REPLICATION_READING_CARRIED.split())


def test_cell2_grouping_reading_is_carried_unstrengthened():
    assert (
        "an operational grouping, NOT established biological independence"
        in c1.READINGS_CARRIED["cell2_group_not_biological_independence"]
    )


def test_withhold_never_reverse_is_carried():
    carried = c1.READINGS_CARRIED["withholds_never_reverses"]
    assert "WITHHOLD" in carried and "never reverse" in carried


def test_unresolved_is_not_absence():
    assert c1.READINGS_CARRIED["unresolved_is_not_absence"] is dv2.UNRESOLVED_MEANS
    assert "never means the element has no action" in dv2.UNRESOLVED_MEANS


# ---- the classes ------------------------------------------------------------------------------


def test_the_five_classes_are_exhaustive_and_exclusive():
    cases = {
        c1.ABSENT: (),
        c1.ONE_TRACK: (track("a", "polyA plus RNA-seq", "K562"),),
        c1.ONE_BIOSAMPLE_NAME_SEVERAL_ASSAY_TITLES: (
            track("a", "polyA plus RNA-seq", "K562"),
            track("b", "total RNA-seq", "K562"),
        ),
        c1.ONE_BIOSAMPLE_NAME_ONE_ASSAY_TITLE: (
            track("a", "total RNA-seq", "K562"),
            track("b", "total RNA-seq", "K562", endedness="single"),
        ),
        c1.SEVERAL_BIOSAMPLE_NAMES: (
            track("a", "total RNA-seq", "donor one"),
            track("b", "total RNA-seq", "donor two"),
        ),
    }
    seen = set()
    for want, tracks in cases.items():
        got = c1.CellTracks(cell="X", rows_in_population=1, tracks=tracks).track_class
        assert got == want, (want, got)
        seen.add(got)
    assert seen == set(c1.TRACK_CLASSES)


def test_several_biosample_names_wins_over_several_assay_titles():
    """A cell whose tracks differ in BOTH is reported by the stronger difference."""
    cell = c1.CellTracks(
        cell="Whole Blood",
        rows_in_population=1,
        tracks=(
            track("a", "polyA plus RNA-seq", "donor one"),
            track("b", "total RNA-seq", "donor two"),
        ),
    )
    assert cell.track_class == c1.SEVERAL_BIOSAMPLE_NAMES


def test_specimen_identity_is_never_established_and_the_row_says_so():
    cell = c1.CellTracks(
        cell="K562",
        rows_in_population=3,
        tracks=(track("a", "polyA plus RNA-seq", "K562"), track("b", "total RNA-seq", "K562")),
    )
    d = cell.to_dict()
    assert d["specimen_identity_established"] is False
    assert "does NOT establish one physical specimen" in d["established"]
    assert "not known" in d["established"]
    assert d["established_how"] == "client_track_metadata"


def test_an_absent_cell_says_absent_and_claims_nothing():
    d = c1.CellTracks(cell="nowhere", rows_in_population=1, tracks=()).to_dict()
    assert d["track_class"] == c1.ABSENT
    assert d["established_how"] == "absent"
    assert "not established here" in d["established"]


def test_one_track_cell_says_two_field_names_over_one_track():
    d = c1.CellTracks(
        cell="K562", rows_in_population=1, tracks=(track("a", "total RNA-seq", "K562"),)
    ).to_dict()
    assert d["track_class"] == c1.ONE_TRACK
    assert "two field names over one track" in d["established"]


def test_nonzero_means_pairwise_distinct_is_none_when_a_value_cannot_be_read():
    cell = c1.CellTracks(
        cell="K562",
        rows_in_population=1,
        tracks=(
            track("a", "polyA plus RNA-seq", "K562", nonzero_mean="nan"),
            track("b", "total RNA-seq", "K562", nonzero_mean="0.9"),
        ),
    )
    assert cell.nonzero_means_pairwise_distinct is None


def test_nonzero_means_distinct_and_identical_are_told_apart():
    distinct = c1.CellTracks(
        cell="K562",
        rows_in_population=1,
        tracks=(
            track("a", "polyA plus RNA-seq", "K562", nonzero_mean="0.1"),
            track("b", "total RNA-seq", "K562", nonzero_mean="0.2"),
        ),
    )
    same = c1.CellTracks(
        cell="K562",
        rows_in_population=1,
        tracks=(
            track("a", "polyA plus RNA-seq", "K562", nonzero_mean="0.1"),
            track("b", "total RNA-seq", "K562", nonzero_mean="0.1"),
        ),
    )
    assert distinct.nonzero_means_pairwise_distinct is True
    assert same.nonzero_means_pairwise_distinct is False


# ---- clause (1), amendment 1 and the argmax, as this lane reads them ---------------------------


def test_clause_1_is_two_retained_values_and_nothing_more():
    one = c1.read_rows({"assertions": [row("K562", ("max_drop",), (-0.8,))]})[0]
    two = c1.read_rows({"assertions": [row("K562", ("max_drop", "by_cell"), (-0.8, -0.7))]})[0]
    assert one.satisfies_clause_1 is False
    assert two.satisfies_clause_1 is True


def test_amendment_1_refuses_two_identical_values_and_accepts_two_distinct():
    same = c1.read_rows({"assertions": [row("K562", ("max_drop", "by_cell"), (-0.8, -0.8))]})[0]
    diff = c1.read_rows({"assertions": [row("K562", ("max_drop", "by_cell"), (-0.8, -0.7))]})[0]
    assert same.satisfies_clause_1 is True and same.passes_amendment_1 is False
    assert diff.passes_amendment_1 is True


def test_the_argmax_reading_matches_direction_v2s_own_retention():
    """`max_drop` is retained exactly when max_drop_tissue is the cell, so the field names it."""
    gene_row = {
        "gene": "G",
        "max_drop_log2fc": -0.8,
        "max_drop_tissue": "K562",
        "max_rise_log2fc": 0.0,
        "max_rise_tissue": "",
        "by_cell": {"K562": -0.7},
    }
    found = dv2.retained_values(gene_row, "K562")
    assert [f for f, _ in found] == ["max_drop", "by_cell"]
    reading = c1.read_rows(
        {"assertions": [row("K562", tuple(f for f, _ in found), tuple(v for _, v in found))]}
    )[0]
    assert reading.cell_is_argmax is True
    assert reading.argmax_kind == "max_drop_only"

    elsewhere = dv2.retained_values({**gene_row, "max_drop_tissue": "HepG2"}, "K562")
    assert [f for f, _ in elsewhere] == ["by_cell"]
    not_argmax = c1.read_rows({"assertions": [row("K562", ("by_cell",), (-0.7,))]})[0]
    assert not_argmax.cell_is_argmax is False
    assert not_argmax.argmax_kind == "none"


def test_a_cell_that_is_both_argmax_fields_is_reported_as_both():
    r = c1.read_rows({"assertions": [row("spleen", ("max_drop", "max_rise"), (-0.8, 0.5))]})[0]
    assert r.argmax_kind == "both" and r.cell_is_argmax is True


def test_outside_retained_cells_two_values_always_straddle_zero():
    """OUTSIDE_RETAINED_CELLS_CANNOT_RESOLVE, held against direction_v2 itself."""
    gene_row = {
        "gene": "G",
        "max_drop_log2fc": -0.8,
        "max_drop_tissue": "spleen",
        "max_rise_log2fc": 0.5,
        "max_rise_tissue": "spleen",
        "by_cell": {},
    }
    assert "spleen" not in CELLS
    call = dv2.direction_call({"tissue": "spleen", "action": "activates"}, dv2.V2, row=gene_row)
    assert call.token is dv2.UNRESOLVED
    assert call.resolved is False
    assert call.reason == "signs_disagree_in_cell"


# ---- the refusals ------------------------------------------------------------------------------


def test_the_population_check_refuses_a_moved_figure():
    good = {
        "counts": c1.COMMITTED_V2_FIGURES,
        "assertions": [row("K562", ("max_drop",), (-0.8,))] * 166,
    }
    rows = c1.read_rows(good)
    c1.check_population(good, rows)
    moved = {**good, "counts": {**c1.COMMITTED_V2_FIGURES, "denominator": 165}}
    with pytest.raises(ValueError, match="denominator"):
        c1.check_population(moved, rows)
    with pytest.raises(ValueError, match="166 rows"):
        c1.check_population(good, rows[:-1])


def test_a_row_naming_no_cell_is_refused():
    for token in c1.CONTEXT_UNKNOWN_TOKENS:
        with pytest.raises(ValueError, match="names no cell"):
            c1.read_rows({"assertions": [row(token, ("max_drop",), (-0.8,))]})


def test_a_row_missing_a_field_is_refused():
    bad = row("K562", ("max_drop",), (-0.8,))
    del bad["retained_values"]
    with pytest.raises(ValueError, match="retained_values"):
        c1.read_rows({"assertions": [bad]})


def test_an_accession_column_refuses_the_undeterminable_reading():
    with pytest.raises(ValueError, match="accession"):
        c1.check_metadata(("name", "biosample_accession"), c1.cc.EXPECTED_AXIS)


def test_a_moved_axis_width_refuses():
    with pytest.raises(ValueError, match="track axis"):
        c1.check_metadata(("name", "Assay title"), c1.cc.EXPECTED_AXIS - 1)
    c1.check_metadata(("name", "Assay title"), c1.cc.EXPECTED_AXIS)


def test_the_saved_metadata_copy_carries_no_accession_column_today():
    """The limitation is checked against the file, not asserted."""
    if not c1.TRACK_METADATA.exists():
        pytest.skip("the saved track metadata copy is not on this machine")
    cols = c1.metadata_columns()
    assert not [c for c in cols if "accession" in c.lower()]
    assert "biosample_name" in cols and "Assay title" in cols and "nonzero_mean" in cols


# ---- the tally --------------------------------------------------------------------------------


def test_the_tally_reports_every_count_named_in_advance():
    rows = c1.read_rows(
        {
            "assertions": [
                row("K562", ("max_drop", "by_cell"), (-0.8, -0.7), v2_class="resolved_agrees_with_published"),
                row(
                    "HepG2",
                    ("max_drop", "by_cell"),
                    (-0.8, -0.8),
                    v2_unresolved_reason="one_track_seen_twice",
                ),
                row("spleen", ("max_drop",), (-0.8,), v2_unresolved_reason="one_value_only"),
            ]
        }
    )
    by_label = {
        "K562": [
            track("a", "polyA plus RNA-seq", "K562"),
            track("b", "total RNA-seq", "K562", nonzero_mean="0.9"),
        ],
        "HepG2": [track("c", "total RNA-seq", "HepG2")],
        "spleen": [track("d", "total RNA-seq", "spleen")],
    }
    cells = c1.classify(rows, by_label)
    tally = c1.Tally(rows=rows, cells=cells, committed=c1.COMMITTED_V2_FIGURES).to_dict()
    for name in c1.COUNTS_NAMED:
        assert name in tally, name
    assert tally["rows_in_population"] == 3
    assert tally["rows_satisfying_clause_1"] == 2
    assert tally["rows_passing_amendment_1"] == 1
    assert tally["rows_whose_cell_is_an_argmax"] == 3
    assert tally["rows_whose_cell_is_not_an_argmax"] == 0
    assert tally["distinct_assigned_cells"] == 3
    assert tally["clause_1_capable_cells"] == 2
    assert tally["cells_by_track_class"][c1.ONE_BIOSAMPLE_NAME_SEVERAL_ASSAY_TITLES] == 1
    assert tally["cells_by_track_class"][c1.ONE_TRACK] == 2
    assert tally["multi_track_cells_whose_nonzero_means_are_pairwise_distinct"] == 1
    assert tally["rows_recorded_one_track_seen_twice_by_track_class"] == {c1.ONE_TRACK: 1}
    assert tally["resolved_rows_whose_cell_is_outside_retained_cells"] == 0


def test_classify_counts_each_cell_on_its_own_label():
    rows = c1.read_rows(
        {"assertions": [row("K562", ("max_drop",), (-0.8,)), row("K562", ("by_cell",), (-0.7,))]}
    )
    cells = c1.classify(rows, {"K562": [track("a", "total RNA-seq", "K562")]})
    assert len(cells) == 1 and cells[0].rows_in_population == 2


def test_the_class_rule_names_one_track_only_as_a_fourth_answer():
    assert "FOURTH answer" in c1.CLASS_RULE


def test_the_denominator_is_kept_apart():
    for n in ("112", "93", "19", "26", "67", "14", "10"):
        assert n in c1.DENOMINATOR_IS_SEPARATE
    assert "not added" in c1.DENOMINATOR_IS_SEPARATE


def test_metadata_columns_reads_only_the_header(tmp_path: Path):
    p = tmp_path / "copy.csv"
    with p.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["name", "Assay title"])
        w.writerow(["x", "total RNA-seq"])
    assert c1.metadata_columns(p) == ("name", "Assay title")


# ---- the counting script's own helpers ---------------------------------------------------------


def _count_module():
    import importlib.util

    path = Path(__file__).resolve().parents[1] / "scripts" / "clause1_count.py"
    spec = importlib.util.spec_from_file_location("clause1_count", path)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_axis_type_is_cellcovers_own():
    assert c1.Axis_t is c1.cc.Axis


def test_row_record_reports_the_committed_class_as_committed():
    mod = _count_module()
    r = c1.read_rows(
        {
            "assertions": [
                row(
                    "K562",
                    ("max_drop", "by_cell"),
                    (-0.8, -0.7),
                    v2_class="resolved_agrees_with_published",
                )
            ]
        }
    )[0]
    rec = mod.row_record(r, {"K562": c1.ONE_BIOSAMPLE_NAME_SEVERAL_ASSAY_TITLES})
    assert rec["v2_class_as_committed"] == "resolved_agrees_with_published"
    assert rec["satisfies_clause_1"] is True and rec["passes_amendment_1"] is True
    assert rec["cell_is_an_argmax"] is True and rec["argmax_kind"] == "max_drop_only"
    assert rec["cell_track_class"] == c1.ONE_BIOSAMPLE_NAME_SEVERAL_ASSAY_TITLES


def test_findings_name_the_denominator_and_make_no_recommendation():
    mod = _count_module()
    rows = c1.read_rows({"assertions": [row("K562", ("max_drop", "by_cell"), (-0.8, -0.7))]})
    cells = c1.classify(
        rows,
        {
            "K562": [
                track("a", "polyA plus RNA-seq", "K562"),
                track("b", "total RNA-seq", "K562", nonzero_mean="0.9"),
            ]
        },
    )
    tally = c1.Tally(rows=rows, cells=cells, committed={}).to_dict()
    f = mod.findings(tally, cells)
    assert "proposes no change to the clause" in f["what_clause_1_secures_on_these_counts"]
    assert "1 rows" in f["the_argmax"] or "denominator of 1" in f["the_argmax"]
    assert f["clause_1_capable_cells"]["with_two_or_more_tracks"] == 1
    assert "not known" in f["what_clause_1_secures_on_these_counts"]
