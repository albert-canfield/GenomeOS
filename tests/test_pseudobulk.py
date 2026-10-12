# SPDX-License-Identifier: AGPL-3.0-or-later
"""`genomeos/attribution/pseudobulk.py`: the rule, the minimum, and what the minimum forbids.

The minimum is the whole protection in this lane. It is stated before anything is counted, so these
tests pin the two numbers and the ORDER in which they are applied: a gene count is never read to
rescue a cell-type count that did not clear, and NO_SOURCE_QUALIFIES is reachable and complete.
"""

from __future__ import annotations

import pytest

from genomeos.attribution import pseudobulk as pb


def test_minimum_is_four_cell_types_and_a_thousand_genes():
    """The two numbers, pinned. A later edit that loosens either one fails here by name."""
    assert pb.MINIMUM["cell_types"] == 4
    assert pb.MINIMUM["genes"] == 1000


def test_four_is_the_smallest_n_whose_permutations_can_reach_the_control_alpha():
    """The reason for 4, as arithmetic and not as prose.

    The registered control is M2h's within-gene cell permutation at alpha 0.05. A gene with n cells
    has n! distinct permutations, and the control p is (1 + hits) / (1 + count), which cannot fall to
    0.05 or below unless at least 19 other arrangements exist. 3! = 6 cannot; 4! = 24 can.
    """
    import math

    alpha = 0.05
    needed = math.ceil(1 / alpha) - 1
    smallest = next(n for n in range(2, 12) if math.factorial(n) - 1 >= needed)
    assert smallest == pb.MINIMUM["cell_types"] == 4
    assert math.factorial(3) - 1 < needed <= math.factorial(4) - 1


@pytest.mark.parametrize(
    ("cells", "genes", "expected"),
    [
        (13, 8328, "QUALIFIES"),
        (4, 1000, "QUALIFIES"),
        (3, 8328, "NO_SOURCE_QUALIFIES"),
        (0, None, "NO_SOURCE_QUALIFIES"),
        (13, 999, "NO_SOURCE_QUALIFIES"),
        (4, None, "NO_SOURCE_QUALIFIES"),
    ],
)
def test_verdict_needs_both_numbers_and_is_not_rescued_by_either(cells, genes, expected):
    """Both clauses bind. Thirteen cell types with 999 genes is NOT a qualifying inventory."""
    assert pb.verdict(cells, genes) == expected


def test_a_count_at_the_minimum_is_in_and_one_below_is_out():
    """The boundary is inclusive on the minimum, exclusive below it, in both clauses."""
    assert pb.verdict(4, 1000) == "QUALIFIES"
    assert pb.verdict(3, 1000) == "NO_SOURCE_QUALIFIES"
    assert pb.verdict(4, 999) == "NO_SOURCE_QUALIFIES"


def test_covered_set_is_the_thirteen_m2h_and_m2i_read_and_is_not_widened():
    """The covered set is read from the committed records, and every entry names its bulk files.

    `data/results/freshrna_result.json` and `freshrna_sources.json` are tracked, so
    `pb.covered_cells()` is readable in every checkout of this commit; its failure is a broken
    checkout, which `tracked_paths.must_be_present` is the helper for.
    """
    import tracked_paths as tp

    tp.must_be_present(pb.FRESHRNA_RESULT, pb.FRESHRNA_SOURCES)
    covered = pb.covered_cells()
    assert len(covered) == 13
    assert [c["cell"] for c in covered] == [
        "K562",
        "HepG2",
        "GM12878",
        "IMR_90",
        "H1",
        "SK_N_SH",
        "astrocyte",
        "hepatocyte",
        "keratinocyte",
        "cardiac_muscle_cell",
        "CD14_positive_monocyte",
        "ovary",
        "testis",
    ]
    for entry in covered:
        assert entry["bulk_assay"] in pb.BULK_ASSAYS
        assert entry["bulk_experiment"].startswith("ENCSR")
        assert entry["bulk_files"], entry["cell"]
        for f in entry["bulk_files"]:
            assert f["accession"].startswith("ENCFF")
            assert len(f["sha256"]) == 64


def test_the_scale_is_m2h_s_constant_and_nothing_in_this_module_fits_it():
    """+0.2146 from `ac5580b`, to the precision M2h committed, and no fitter in this module."""
    assert round(pb.SCALE, 4) == 0.2146
    source = (pb.ROOT / "genomeos" / "attribution" / "pseudobulk.py").read_text()
    for forbidden in ("curve_fit", "minimize", "lstsq", "polyfit", "LinearRegression"):
        assert forbidden not in source


def test_the_wording_limit_is_carried_and_never_strengthened():
    """The registered wording. Agreement licenses the modality swap for THIS endpoint only."""
    assert "FOR THIS ENDPOINT ONLY" in pb.WORDING_LIMIT
    assert "never" in pb.WORDING_LIMIT
    assert "single-nucleus RNA equals bulk RNA" in pb.WORDING_LIMIT
    assert pb.SELECTION_RULE["wording_limit"] == pb.WORDING_LIMIT


def test_no_source_qualifies_is_declared_an_acceptable_and_complete_answer():
    """If the covered cell types do not offer both modalities, saying so is the lane's whole value."""
    assert "COMPLETE" in pb.SELECTION_RULE["acceptable_answer"]
    assert "NO_SOURCE_QUALIFIES" in pb.SELECTION_RULE["acceptable_answer"]
    assert "NO_SOURCE_QUALIFIES" in pb.MINIMUM["below_either"]


def test_the_rule_adds_no_biosample_and_makes_no_model_request():
    """The stop list of 2026-09-27, named in the rule itself."""
    clause = pb.SELECTION_RULE["no_biosample"]
    assert "adds NO biosample" in clause
    assert "AlphaGenome request" in clause
    assert "area J" in clause
