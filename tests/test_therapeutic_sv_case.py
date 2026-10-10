# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The structural-variant case: the record read as its caller exports it.

Two things are pinned here, and they are different kinds of claim.

The first is the reader. `read_sv_table` took four columns in a fixed order,
and a real structural-variant export does not have that shape: cBioPortal's
first column is the sample identifier and its fourth is a chromosome, so the
positional reading invents a fusion on a gene named after the sample and reads
"4" as the end of the product. A table whose first row names its columns is now
read by those names, and the tests below hold both readings against each other —
the hand-written fixture of the ten-case benchmark must parse to exactly what it
parsed to before, because the ten rows are a regression guard and not only a
result.

The second is the case's registration, `data/results/therapeutic_sv_case_registration.json`,
and the result beside it. The registration was committed before the case was
run; these tests assert that the published row answers the three quantities the
registration named, and they skip rather than pass when the result is absent,
because a missing file is not a measurement.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from genomeos.cancer.alterations import read_sv_table
from genomeos.results import load_result

FIXTURES = Path("data/demo/alterations")
CASE_SV = FIXTURES / "fgfr3_tacc3_urothelial.sv"
ALK_SV = FIXTURES / "alk_eml4_fusion.sv"


# ---- the reader -------------------------------------------------------------------


def test_the_hand_written_fixture_parses_exactly_as_before():
    """The benchmark's ALK fixture, which the ten rows depend on, is untouched.

    `gene partner class end` names its columns, so it now goes through the
    named reading; what it must produce is what the positional reading produced,
    field for field, including the orientation that makes ALK the 3' partner.
    """
    alterations = read_sv_table(str(ALK_SV))
    assert len(alterations) == 1
    a = alterations[0]
    assert (a.gene, a.kind, a.partner, a.detail, a.orientation) == (
        "ALK",
        "fusion",
        "EML4",
        "INVERSION",
        "3'",
    )
    assert a.in_frame == ""
    assert a.junction == ""
    assert a.source == "patient structural-variant table (alk_eml4_fusion.sv)"
    assert a.evidence == ["ALK-EML4 fusion in this tumour's structural-variant table"]


def test_a_gene_and_partner_header_does_not_infer_an_orientation(tmp_path):
    """The gene order in a fusion's name is a naming convention, not a measurement."""
    p = tmp_path / "named.sv"
    p.write_text("gene\tpartner\tclass\nALK\tEML4\tINVERSION\n")
    (a,) = read_sv_table(str(p))
    assert a.orientation == "", "a gene/partner pair says nothing about which end is which"
    assert (a.gene, a.partner, a.detail) == ("ALK", "EML4", "INVERSION")


def test_a_site_pair_carries_the_orientation_in_its_column_names():
    """site1 is the 5' partner and site2 the 3', so one row yields both ends."""
    alterations = read_sv_table(str(CASE_SV))
    assert [(a.gene, a.orientation) for a in alterations] == [("FGFR3", "5'"), ("TACC3", "3'")]
    fgfr3, tacc3 = alterations
    assert fgfr3.partner == "TACC3" and tacc3.partner == "FGFR3"
    assert all(a.kind == "fusion" for a in alterations)
    assert all(a.in_frame == "in-frame" for a in alterations)
    assert all(a.junction == "4:1806934 to 4:1739702 (GRCh37)" for a in alterations)
    assert all(a.detail == "FGFR3-TACC3" for a in alterations)


def test_the_sample_identifier_is_not_read_as_a_gene():
    """The defect the named reading exists for: a candidate named after the sample."""
    genes = {a.gene for a in read_sv_table(str(CASE_SV))}
    assert genes == {"FGFR3", "TACC3"}
    assert not any(g.startswith("TCGA") for g in genes), "the sample identifier is not a gene"


def test_a_table_with_no_header_is_still_read_positionally(tmp_path):
    p = tmp_path / "bare.sv"
    p.write_text("ALK\tEML4\tINVERSION\t3\n")
    (a,) = read_sv_table(str(p))
    assert (a.gene, a.partner, a.detail, a.orientation) == ("ALK", "EML4", "INVERSION", "3'")


def test_absent_values_are_absent_and_not_text(tmp_path):
    """A caller writes NA where it has no value; NA is not an event class."""
    p = tmp_path / "na.sv"
    p.write_text(
        "Site1_Hugo_Symbol\tSite2_Hugo_Symbol\tVariant_Class\tSite2_Effect_On_Frame\tAnnotation\n"
        "FGFR3\tTACC3\tNA\tNA\tNA\n"
    )
    first, second = read_sv_table(str(p))
    assert first.detail == "" and first.in_frame == "" and first.junction == ""
    assert (first.gene, second.gene) == ("FGFR3", "TACC3")


def test_the_annotation_is_kept_verbatim_when_the_source_gives_one(tmp_path):
    p = tmp_path / "annotated.sv"
    p.write_text(
        "Site1_Hugo_Symbol\tSite2_Hugo_Symbol\tAnnotation\tSite1_Position\tSite2_Position\n"
        "FGFR3\tTACC3\tFGFR3 exon 17 to TACC3 exon 11\t1806934\t1739702\n"
    )
    first, _ = read_sv_table(str(p))
    assert first.junction == "FGFR3 exon 17 to TACC3 exon 11", "the source's own words, not ours"


def test_a_header_needs_a_gene_column_and_one_more(tmp_path):
    """One recognised word is not a header; the row stays data."""
    p = tmp_path / "oneword.sv"
    p.write_text("GENE\nALK\n")
    out = read_sv_table(str(p))
    assert [a.gene for a in out] == ["ALK"], "the old header skip still applies to a bare GENE row"


def test_a_missing_file_is_empty_and_not_an_error(tmp_path):
    assert read_sv_table(str(tmp_path / "nothing.sv")) == []


# ---- the registration and the run ------------------------------------------------

REGISTRATION = load_result("therapeutic_sv_case_registration")
RESULT = load_result("therapeutic_sv_case_result")


@pytest.mark.skipif(not REGISTRATION, reason="registration not present")
def test_the_registration_names_the_three_quantities_and_the_guard():
    r = REGISTRATION
    assert r["case"]["gene"] == "FGFR3"
    assert r["case"]["sample"] == "TCGA-CF-A3MG-01"
    assert r["case"]["expect"] == "out_of_scope_expected"
    assert set(r["registered"]) >= {"recovered", "route", "top_mechanism", "regression_guard"}
    assert r["eleventh_row"]["added"] == "whatever the outcome"


@pytest.mark.skipif(not RESULT, reason="the case has not been run")
def test_the_ten_rows_did_not_move():
    guard = RESULT["regression_guard"]
    assert guard["rows_compared"] == 10
    assert guard["fields_that_differ"] == [], guard["fields_that_differ"]


@pytest.mark.skipif(not RESULT, reason="the case has not been run")
def test_the_case_row_answers_the_three_quantities():
    row = RESULT["row"]
    assert row["gene"] == "FGFR3"
    assert row["driver_call"] == "structural variant"
    assert row["recovered"] is True, "FGFR3 is reachable only through the structural call"
    assert row["expected"] == "out_of_scope_expected"
    # The route question for this case: no surface route may be claimed for a
    # product whose extent the record does not establish.
    assert row["pass"] is True
    acc = row["surface_accessibility"]
    assert acc is None or acc <= 0.5, f"a surface route claimed at accessibility {acc}"
    assert row["mechanism_sane"] is True, row.get("mechanism_concerns")
