# SPDX-License-Identifier: AGPL-3.0-or-later
"""Amendment 1 of the second-targets registration and its run (scripts/secondtargets_run_readout.py), on
synthetic tables: the join, the amended fold column, and that nothing else of the hit rule moved."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import secondtargets_run_readout as rr  # noqa: E402

from genomeos.attribution import secondtargets as s2  # noqa: E402

HEAD = (
    "pairs4merge\tbeta\tintercept\tfold_change.transcript_remaining\tpvalue.raw\tpvalue.empirical\t"
    "pvalue.empirical.adjusted\tENSG\n"
)


def row(pair: str, fold: str, padj: str) -> str:
    return f"{pair}\t0\t0\t{fold}\t0.5\t0.5\t{padj}\tE\n"


def test_the_join_upper_cases_the_group_and_drops_the_gene_version():
    assert rr.norm("chr1.10703_top_two:ENSG00000143862.4") == "CHR1.10703_TOP_TWO:ENSG00000143862"
    assert rr.norm('"CHR1.10703_TOP_TWO:ENSG00000143862"') == "CHR1.10703_TOP_TWO:ENSG00000143862"
    assert rr.norm("ACTB_TSS:ENSG00000075624") == "ACTB_TSS:ENSG00000075624"
    assert rr.wanted_of(["chr1.1_top_two:ENSG1", "CHR1.1_TOP_TWO:ENSG1"]) == {
        "CHR1.1_TOP_TWO:ENSG1": ["CHR1.1_TOP_TWO:ENSG1", "chr1.1_top_two:ENSG1"]
    }


def test_only_the_fold_column_and_the_pair_column_differ_from_the_registered_rule():
    changed = {k for k in s2.HIT if s2.HIT[k] != rr.HIT_AMENDED[k]}
    assert changed == {"fold_column", "pair_column"}
    assert rr.HIT_AMENDED["fold_column"] == "fold_change.transcript_remaining"
    assert rr.HIT_AMENDED["pair_column"] in s2.HIT["pair_column"]  # the registered first choice
    assert (rr.HIT_AMENDED["p_below"], rr.HIT_AMENDED["fold_below"]) == (0.1, 1.0)


def test_the_amended_reader_keys_its_answers_by_the_registered_ids():
    lines = [
        HEAD,
        row("CHR1.1_TOP_TWO:ENSG1", "0.7", "0.01"),  # a decrease, called: hit
        row("CHR1.1_TOP_TWO:ENSG2", "1.3", "0.01"),  # an increase, called: not a hit
        row("CHR1.1_TOP_TWO:ENSG3", "0.9", "0.40"),  # not called
        row("CHR1.1_TOP_TWO:ENSG4", "NA", "NA"),  # not measured
        row("CHR9.9_TOP_TWO:ENSG9", "0.1", "0.001"),  # not wanted
    ]
    keys = ["chr1.1_top_two:ENSG1", "chr1.1_top_two:ENSG2", "chr1.1_top_two:ENSG3", "chr1.1_top_two:ENSG4"]
    got = rr.read_amended(lines, rr.wanted_of(keys + ["chr1.1_top_two:ENSG5"]))
    assert got == {
        "chr1.1_top_two:ENSG1": True,
        "chr1.1_top_two:ENSG2": False,
        "chr1.1_top_two:ENSG3": False,
        "chr1.1_top_two:ENSG4": None,
    }


def test_the_amended_reader_stops_before_any_value_when_a_column_is_absent():
    with pytest.raises(KeyError):
        rr.read_amended(["pairs4merge\tfold_change.transcript_remap\tpvalue.empirical.adjusted\n"], {})


def test_the_identifier_pass_reads_the_pair_column_only():
    lines = [HEAD, row("CHR1.1_TOP_TWO:ENSG1", "x", "y"), row("ACTB_TSS:ENSG7", "x", "y"), "\n"]
    assert rr.identifiers(lines, "pairs4merge") == {"CHR1.1_TOP_TWO:ENSG1", "ACTB_TSS:ENSG7"}
    assert rr.missing_columns(HEAD.rstrip("\n").split("\t"), rr.HIT_AMENDED) == []
    assert rr.missing_columns(HEAD.rstrip("\n").split("\t"), {**s2.HIT, "pair_column": "pairs4merge"}) == [
        "fold_change.transcript_remap"
    ]


def test_the_run_writes_a_new_result_and_never_the_not_run_record():
    assert rr.RESULT == "secondtargets_result_run" and rr.SUPERSEDES.endswith("secondtargets_result.json")
    assert rr.AMENDED_CODE == ("scripts/secondtargets_run_readout.py",)
    assert rr.FETCH["bytes"] == 37256050 and len(rr.FETCH["sha256"]) == 64
    assert rr.FETCH["url"].endswith("GSE120861_all_deg_results.at_scale.txt.gz")
