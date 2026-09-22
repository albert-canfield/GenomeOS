# SPDX-License-Identifier: AGPL-3.0-or-later
"""The Ensembl-id join that the fourth frame's step 4 should have been using all along."""

from __future__ import annotations

from genomeos.benchmark import loci_fourth as lf


def test_the_join_is_registered_before_it_is_applied():
    for key in (
        "landed",
        "the_defect",
        "what_changed",
        "it_cuts_both_ways",
        "why_it_is_a_flag_and_not_a_replacement",
        "the_cap_binds_again_and_what_is_done_about_it",
    ):
        assert lf.ENSEMBL_JOIN.get(key), f"the join says nothing about {key}"


def test_with_current_symbols_rewrites_targets_and_keeps_the_published_ones():
    rows = [
        {
            "chrom": "chr2",
            "start": 10,
            "end": 20,
            "cell": "HCT116",
            "targets": ["SSFA2"],
            "distances": {"SSFA2": 5000},
        }
    ]
    out = lf.with_current_symbols(rows, {"SSFA2": "ITPRID2"})
    assert out[0]["targets"] == ["ITPRID2"]
    assert out[0]["published_targets"] == ["SSFA2"]
    assert out[0]["distances"] == {"ITPRID2": 5000}


def test_a_target_nobody_renamed_is_left_exactly_alone():
    rows = [
        {"chrom": "chr1", "start": 1, "end": 2, "cell": "K562", "targets": ["MYC"], "distances": {"MYC": 7}}
    ]
    out = lf.with_current_symbols(rows, {"SSFA2": "ITPRID2"})
    assert out[0]["targets"] == ["MYC"] and out[0]["published_targets"] == ["MYC"]


def test_assess_is_untouched_and_it_is_the_names_it_is_given_that_changed():
    """The rule's own steps must behave identically; only the symbols reaching them are corrected."""
    row = {
        "chrom": "chr2",
        "start": 100,
        "end": 200,
        "mid": 150,
        "cell": "HCT116",
        "targets": ["ITPRID2"],
        "distance": 5000,
    }
    stale = {**row, "targets": ["SSFA2"]}
    coding = {"ITPRID2", "NEUROD1"}
    assert lf.assess(stale, ("NEUROD1", 900), coding, [])["reason"] == (
        "a regulated target is not protein coding in GENCODE"
    )
    assert lf.assess(row, ("NEUROD1", 900), coding, [])["reason"] is None
