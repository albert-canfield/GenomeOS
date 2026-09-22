# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration that fixes the miss classes first, and the classifier that sorts misses into them."""

from __future__ import annotations

from dataclasses import dataclass

import pytest

from genomeos.benchmark import loci_miss as lm


@dataclass
class _Locus:
    chrom: str
    start: int
    end: int


@dataclass
class _Gene:
    id: str
    symbol: str
    locus: _Locus


class _Ann:
    """The two attributes `classify_gene` uses of a GENCODE annotation, and nothing else."""

    def __init__(self, genes: list[_Gene]) -> None:
        self.genes = {g.id: g for g in genes}


TARGET = _Gene("HMGA1", "HMGA1", _Locus("chr6", 34_236_000, 34_246_000))
ANN = _Ann(
    [
        TARGET,
        # inside the target's body: the control's promoted gene, the case the whole question is about
        _Gene("ENSG00000288879", "ENSG00000288879", _Locus("chr6", 34_240_000, 34_242_000)),
        _Gene("NEAR", "NEAR", _Locus("chr6", 34_250_000, 34_252_000)),  # 4 kb
        _Gene("MIDDLING", "MIDDLING", _Locus("chr6", 34_300_000, 34_310_000)),  # 54 kb
        _Gene("FAR", "FAR", _Locus("chr6", 40_000_000, 40_010_000)),  # 5.7 Mb
        _Gene("TWIN.1", "TWIN", _Locus("chr6", 34_260_000, 34_261_000)),
        _Gene("TWIN.2", "TWIN", _Locus("chr6", 39_000_000, 39_001_000)),
    ]
)


def _cls(named: str | None, node: list[str] | None = None) -> str:
    return lm.classify_gene(ANN, named, ["HMGA1"], node or [])["class"]


def test_the_registration_states_what_it_has_to_state_before_a_single_miss_is_classified():
    """A registration missing any of these is not one, and no distribution below it may be read."""
    p = lm.PREREGISTRATION
    for key in (
        "written",
        "the_question",
        "what_this_is_not",
        "what_is_classified",
        "which_gene_is_classified",
        "the_classes",
        "the_classes_are_fixed_first",
        "only_overlap_can_be_tolerated",
        "what_would_make_the_case_for_a_second_rate",
        "what_would_close_the_question",
        "in_between",
        "the_trap",
        "the_classifier_control",
        "falsifiers",
        "cost",
    ):
        assert p.get(key), f"the registration says nothing about {key}"
    assert p["the_classes"] == [c["name"] for c in lm.CLASSES]
    assert len(p["falsifiers"]) >= 5


def test_the_trap_is_named_and_the_tolerant_rate_is_bound_to_be_reported_beside_the_strict_one():
    """The one sentence this work exists to not fall foul of, and the constraint that follows it."""
    trap = lm.PREREGISTRATION["the_trap"]
    assert "widens what counts as a hit after seeing its rate" in trap
    assert "BESIDE" in trap and "never in place of it" in trap
    assert "0 AlphaGenome requests" in lm.PREREGISTRATION["cost"]


def test_only_the_overlap_class_can_ever_feed_a_tolerant_rate():
    """Admitting a distance class would score the nearest-gene baseline as a hit, and that baseline
    is the control the fourth frame was drawn to beat."""
    assert lm.TOLERATED == ("exact", "overlaps_the_target_body")
    for name in ("within_10_kb", "within_100_kb", "same_ctcf_node", "elsewhere"):
        assert name not in lm.TOLERATED


def test_a_published_symbol_is_an_exact_match_and_nothing_else():
    assert _cls("HMGA1") == "exact"


def test_a_gene_whose_body_lies_inside_the_targets_is_an_overlap_not_a_distance():
    """The control's case: ENSG00000288879 sits inside HMGA1, so it is the right place under
    another name rather than a gene 2 kb away."""
    out = lm.classify_gene(ANN, "ENSG00000288879", ["HMGA1"], [])
    assert out["class"] == "overlaps_the_target_body"
    assert out["gap_bp"] == 0


def test_the_distance_classes_cut_where_the_registration_says_they_cut():
    assert _cls("NEAR") == "within_10_kb"
    assert _cls("MIDDLING") == "within_100_kb"
    assert _cls("FAR") == "elsewhere"


def test_a_far_gene_in_the_elements_own_node_is_a_node_member_and_a_near_one_is_still_near():
    """The order matters: a gene 4 kb away is reported as near even when both sit in the node."""
    assert _cls("FAR", ["FAR", "HMGA1"]) == "same_ctcf_node"
    assert _cls("NEAR", ["NEAR", "HMGA1"]) == "within_10_kb"
    assert _cls("FAR", ["FAR"]) == "elsewhere"


def test_a_symbol_gencode_uses_more_than_once_is_not_locatable_rather_than_placed():
    """Y_RNA names 57 loci on chr2; a distance to whichever copy comes first is the wrong gene."""
    assert _cls("TWIN") == "not_locatable"
    assert _cls("NOT_IN_GENCODE") == "not_locatable"
    assert _cls(None) == "named_nothing"


def test_a_published_target_with_no_single_body_leaves_the_named_gene_unplaced():
    assert lm.classify_gene(ANN, "NEAR", ["TWIN"], [])["class"] == "not_locatable"


def test_the_locus_takes_the_best_class_over_its_derived_layers_and_only_derived_ones():
    """The headline hit is a union over derived layers, so the locus class has to be one too, and a
    heuristic layer that names the target must not turn a miss into an overlap."""
    row = {
        "locus": "CONTROL_HMGA1_WTC11_chr6_34235k",
        "expected": {"chrom": "chr6", "targets": ["HMGA1"]},
        "readings": {"node": {"genes_in_node": ["HMGA1", "FAR"]}},
        "score": {
            "target_hit_derived": False,
            "scored": {
                "target": {
                    "by_layer": {
                        "deletion": {"provenance": "derived", "named": ["FAR"], "hit": False},
                        "eqtl": {"provenance": "derived", "named": ["ENSG00000288879"], "hit": False},
                        "node": {"provenance": "heuristic", "named": ["HMGA1"], "hit": True},
                    }
                }
            },
        },
    }
    out = lm.classify_locus(row, {"chr6": ANN})
    assert out["class"] == "overlaps_the_target_body"
    assert out["class_from_layer"] == "eqtl"
    assert out["deletion_class"] == "same_ctcf_node"
    assert "node" not in out["layers"]


def test_a_strict_hit_that_does_not_classify_as_exact_stops_everything():
    """The classifier's control. A hit that lands anywhere else means the lookup or the join is
    wrong, and the registration says nothing else in the output may be read."""
    row = {
        "locus": "BROKEN",
        "expected": {"chrom": "chr6", "targets": ["HMGA1"]},
        "readings": {},
        "score": {
            "target_hit_derived": True,
            "scored": {
                "target": {"by_layer": {"deletion": {"provenance": "derived", "named": ["FAR"], "hit": True}}}
            },
        },
    }
    frame = lm.classify_frame("broken", [row], {"chr6": ANN})
    assert frame["control_hits_that_are_not_exact"] == ["BROKEN"]


def test_the_tolerant_rate_is_computed_beside_the_strict_one_and_names_what_it_would_add():
    rows = [
        {
            "locus": f"L{i}",
            "expected": {"chrom": "chr6", "targets": ["HMGA1"]},
            "readings": {},
            "score": {
                "target_hit_derived": named == "HMGA1",
                "scored": {
                    "target": {
                        "by_layer": {
                            "deletion": {
                                "provenance": "derived",
                                "named": [named],
                                "hit": named == "HMGA1",
                            }
                        }
                    }
                },
            },
        }
        for i, named in enumerate(["HMGA1", "ENSG00000288879", "NEAR", "FAR"])
    ]
    frame = lm.classify_frame("toy", rows, {"chr6": ANN})
    h = frame["headline"]
    assert (h["strict_k"], h["n"]) == (1, 4)
    assert (h["tolerant_k"], h["loci_it_would_add"]) == (2, ["L1"])
    assert frame["miss_classes"]["overlaps_the_target_body"] == 1
    assert frame["miss_classes"]["within_10_kb"] == 1
    assert frame["miss_classes"]["elsewhere"] == 1
    assert "DESCRIPTION" in h["reading"]


def test_every_frame_the_repaired_reader_scored_has_a_published_rate_to_be_checked_against():
    """Falsifier 4: a frame whose strict rate is not pinned here could drift without anyone seeing."""
    assert set(lm.PUBLISHED_STRICT) == set(lm.FRAMES)


def test_the_classifier_refuses_a_frame_that_is_not_on_disk():
    with pytest.raises(FileNotFoundError):
        lm.rows_for("loci_not_a_frame", lm.RESULTS_DIR)
