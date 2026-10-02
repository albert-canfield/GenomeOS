"""The registration for `is a rule's cell the element's or the target gene's`, and its instrument.

These tests guard the registration's words, not a number: that the control is the label base rate
and never the track count, that the band rule is relative, that reading (c) says "the data cannot
tell", that the clause saying what this population cannot answer survives, and that the RSS check
RAISES instead of warning.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from genomeos.attribution import label_gene as lg

REGISTRATION = Path("data/results/label_gene_registration.json")


@pytest.mark.skipif(
    not REGISTRATION.exists(), reason="the registration is written by its own committed script"
)
def test_registration_on_disk_carries_every_registered_field_unchanged() -> None:
    """The committed file must still say what the module says: save_result adds fields, never
    rewrites one."""
    on_disk = json.loads(REGISTRATION.read_text(encoding="utf-8"))
    for key, value in lg.registration().items():
        assert on_disk[key] == value, f"{key} on disk differs from the module's registration"


def test_control_one_is_the_base_rate_and_not_the_track_count() -> None:
    text = lg.CONTROL_ONE
    assert "sum of p squared" in text
    assert "NOT one over the track count" in text
    assert "base rate and not the track count is the control" in text


def test_the_clause_this_population_cannot_answer_is_registered() -> None:
    clause = lg.registration()["cannot_be_answered"]
    assert "CANNOT BE RUN on this population at all" in clause
    assert 'It does NOT mean "an element property"' in clause


def test_reading_c_says_the_data_cannot_tell_and_not_no_information() -> None:
    assert "THE DATA CANNOT TELL" in lg.READINGS["c"]
    assert "Never 'no information'" in lg.READINGS["c"]


def test_the_unit_of_analysis_is_genes() -> None:
    assert "THE UNIT OF ANALYSIS IS GENES, NOT PAIRS" in lg.UNIT_OF_ANALYSIS
    assert lg.MIN_CLUSTERS == 10


def test_the_training_exposure_caveat_is_registered_in_advance() -> None:
    caveat = lg.TRAINING_EXPOSURE_CAVEAT
    assert "RECALL OF SEEN DATA rather than biology" in caveat
    assert "must not be softened" in caveat


def test_the_model_consequence_is_a_proposal_and_recommends_nothing() -> None:
    assert "PROPOSAL ONLY" in lg.MODEL_CONSEQUENCE
    assert "NO recommendation about adopting" in lg.MODEL_CONSEQUENCE


def test_control_one_on_a_known_distribution() -> None:
    rows = [lg.RuleRow(f"e{i}", "G", "A" if i < 2 else "B", "chr21", i, i + 1) for i in range(4)]
    # two of four A, two of four B: 0.5^2 + 0.5^2
    assert lg.control_one(rows) == pytest.approx(0.5)


def test_control_one_is_not_one_over_the_number_of_labels_when_skewed() -> None:
    labels = ["A"] * 90 + ["B"] * 5 + ["C"] * 5
    rows = [lg.RuleRow(f"e{i}", "G", lab, "chr21", i, i + 1) for i, lab in enumerate(labels)]
    assert lg.control_one(rows) == pytest.approx(0.9**2 + 0.05**2 + 0.05**2)
    assert lg.control_one(rows) > 1 / 3


def test_distance_bins_cover_the_registered_edges() -> None:
    assert lg.distance_bin(0) == 0
    assert lg.distance_bin(1_000) == 0
    assert lg.distance_bin(1_001) == 1
    assert lg.distance_bin(10_000) == 1
    assert lg.distance_bin(100_000) == 2
    assert lg.distance_bin(1_000_000) == 3
    assert lg.distance_bin(10_000_000) == 4
    assert lg.distance_bin(40_000_000) == 5
    assert len(lg.BIN_LABELS) == len(lg.DISTANCE_BIN_EDGES_BP) - 1
    assert lg.DISTANCE_BIN_EDGES_BP[-1] == math.inf


def test_bands_are_relative_and_read_as_registered() -> None:
    assert lg.band_of([0.95, 1.10]) == "at_control"
    assert lg.band_of([1.60, 2.40]) == "far_above"
    assert lg.band_of([1.20, 1.90]) == "crosses_a_band_edge"
    assert lg.band_of([1.40, 1.45]) == "crosses_a_band_edge"
    assert lg.band_of(None) == "no_interval"


def test_the_four_readings_are_selected_by_the_registered_rule_alone() -> None:
    which, _ = lg.read_the_comparison("far_above", "far_above", "at_control")
    assert which == "a"
    which, words = lg.read_the_comparison("far_above", "far_above", "far_above")
    assert which == "a" and "(a-ii)" in words
    which, _ = lg.read_the_comparison("at_control", "far_above", "far_above")
    assert which == "a_prime"
    which, _ = lg.read_the_comparison("at_control", "at_control", "at_control")
    assert which == "b"
    which, words = lg.read_the_comparison("crosses_a_band_edge", "far_above", "at_control")
    assert which == "c" and "DATA CANNOT TELL" in words
    which, words = lg.read_the_comparison("no_interval", "far_above", "at_control")
    assert which == "c" and "INCONCLUSIVE BY RULE" in words


def test_a_pure_gene_property_is_not_denied_by_the_chains_second_link() -> None:
    """A pure gene property puts the different-gene arm AT control 1; that must still read (a)."""
    which, words = lg.read_the_comparison("far_above", "far_above", "at_control")
    assert which == "a"
    assert "(a-i)" in words


def test_identical_resample_share_is_reported_for_every_interval() -> None:
    flat = lg.percentile_interval([0.4] * 500)
    assert flat["identical_share"] == 1.0 and flat["degenerate"] is True
    assert flat["ci95"] == [0.4, 0.4]
    varied = lg.percentile_interval([i / 1000 for i in range(1000)])
    assert varied["identical_share"] == pytest.approx(0.001)
    assert varied["degenerate"] is False
    assert lg.percentile_interval([])["ci95"] is None


def test_degenerate_bootstrap_is_named_and_routed() -> None:
    assert lg.bootstrap_degenerate(0, 10) is True
    assert lg.bootstrap_degenerate(10, 10) is True
    assert lg.bootstrap_degenerate(3, 10) is False
    assert lg.bootstrap_degenerate(0, 0) is False
    assert "DEGENERATE" in lg.DEGENERACY_RULE
    assert "WILSON bound on an EFFECTIVE sample size" in lg.DEGENERACY_RULE


def test_wilson_bound_is_available_for_the_degenerate_route() -> None:
    assert lg.wilson(0, 0) is None
    low = lg.wilson(0, 20)
    assert low is not None and low[0] == 0.0 and 0.0 < low[1] < 0.2


def test_the_rss_check_raises_and_does_not_warn() -> None:
    with pytest.raises(MemoryError):
        lg.check_rss(ceiling=1)
    assert lg.check_rss() > 0
    assert "RAISES - it does not warn" in lg.RSS_RULE
    assert lg.RSS_CEILING_BYTES == 4 * 1024**3


def test_the_label_to_gtex_column_rule_is_the_registered_one() -> None:
    assert lg.normalise_tissue("Brain_Cerebellar_Hemisphere") == "brain_cerebellar_hemisphere"
    assert lg.normalise_tissue("Brain - Cerebellar Hemisphere") == "brain_cerebellar_hemisphere"
    assert lg.normalise_tissue("Adipose - Visceral (Omentum)") == "adipose_visceral_omentum"
    assert lg.normalise_tissue("CD8-positive, alpha-beta T cell") == "cd8_positive_alpha_beta_t_cell"


def test_cluster_multiplicities_draw_whole_genes_with_replacement() -> None:
    import random

    multiplicity = lg.cluster_multiplicities(50, random.Random(7))
    assert len(multiplicity) == 50
    assert sum(multiplicity) == 50
    assert max(multiplicity) > 1


def test_genes_needed_is_reported_and_never_licenses_a_reading() -> None:
    assert lg.genes_needed(0.0, 2.0) is None
    assert lg.genes_needed(1.0, 0.0) is None
    few = lg.genes_needed(1.0, 1.0)
    many = lg.genes_needed(1.0, 9.0)
    assert few is not None and many is not None and many > few
    assert "never used to license a reading" in (lg.genes_needed.__doc__ or "")


def test_declared_inputs_are_hashed_in_the_registration() -> None:
    inputs = lg.registration()["declared_inputs"]
    assert "data/organisms/human/noncoding_chr21.bio" in inputs
    assert all(len(digest) == 64 for digest in inputs.values())


def test_the_registration_promises_no_spend_and_no_weakening() -> None:
    assert "0 model requests" in lg.NO_SPEND
    assert "No test, threshold, pin or registration is weakened" in lg.NEVER_WEAKENED
