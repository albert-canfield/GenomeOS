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


def test_amendment_one_is_additive_and_discloses_that_it_appends() -> None:
    text = lg.AMENDMENT_1
    assert "DISCLOSED HONESTLY" in text
    assert "47835d1" in text
    assert "purely ADDITIVE" in text
    assert "scripts/label_gene_count.py does not exist in the tree" in text
    body = lg.amendment_1_payload()
    assert body["amendment_1"] == text
    assert "47835d1" in body["amends"]
    # the amendment is its own file: the landed registration is not rewritten
    assert "amendment_1" not in lg.registration()
    # nothing the registration already fixed may move
    assert body["band_at_control"] == [0.80, 1.25]
    assert body["band_far_above"] == 1.50
    assert body["min_clusters_genes"] == 10
    assert body["reading_c_words"] == "THE DATA CANNOT TELL"


def test_the_binding_sentence_on_disagreement_is_registered_before_either_number() -> None:
    assert (
        "IF THE PAIR-POOLED PRIMARY AND THE GENE-EQUAL-WEIGHT SENSITIVITY DISAGREE, THE RESULT "
        "MUST SAY SO IN THOSE WORDS" in lg.AMENDMENT_1
    )
    assert "DECIDES NOTHING" in lg.AMENDMENT_1


def test_disagreement_fires_on_either_registered_condition() -> None:
    pooled = {"ci95": [1.60, 1.90], "band": "far_above"}
    same = lg.disagreement(pooled, {"band": "far_above"}, 1.70)
    assert same["disagree"] is False
    bands = lg.disagreement(pooled, {"band": "at_control"}, 1.70)
    assert bands["disagree"] is True
    assert "DISAGREE." in bands["statement"]
    outside = lg.disagreement(pooled, {"band": "far_above"}, 2.40)
    assert outside["disagree"] is True
    assert outside["equal_weight_point_outside_the_pooled_interval"] is True
    assert "selects the reading on the PAIR-POOLED primary" in outside["the_reading_is_still_the_pooled_one"]


# --------------------------------------------------------------------------- #
# The counting path, on a synthetic population small enough to check by hand.
# A @1000 B @1500 C @11000 on gene G1; D @2000 E @1001000 on gene G2.
# Labels: A=L1 B=L1 C=L2 D=L1 E=L2.
# --------------------------------------------------------------------------- #


def _block(element: str, gene: str, label: str, start: int, verb: str = "activates") -> str:
    """One element block and its predicted rule line, built so no source line runs long."""
    evidence = f'predicted "AlphaGenome deletion, {label}" effect -0.1'
    return (
        f"element {element} {{\n  locus: chr21:{start}-{start}\n}}\n"
        f"rule {element} {verb} {gene} {{ strength: 0.1; when: cell_type = {label}; "
        f"evidence: {evidence} }}\n"
    )


_MEASURED = (
    "element F_measured {\n  locus: chr21:1000-1000\n}\n"
    "rule F_measured activates G1 { strength: 0.2; when: cell_type = K562; "
    'evidence: experimental "a screen"; confidence: 0.9 }\n'
)

SYNTHETIC = (
    "module t\n\n"
    + _block("A", "G1", "L1", 1000)
    + _block("B", "G1", "L1", 1500, verb="inhibits")
    + _block("C", "G1", "L2", 11000)
    + _block("D", "G2", "L1", 2000)
    + _block("E", "G2", "L2", 1001000)
    + _MEASURED
)


@pytest.fixture
def synthetic(tmp_path):
    path = tmp_path / "synthetic.bio"
    path.write_text(SYNTHETIC, encoding="utf-8")
    rows, census = lg.parse_population(str(path))
    return rows, census


def test_the_population_excludes_the_experimental_rule(synthetic) -> None:
    rows, census = synthetic
    assert census["rule_lines"] == 6
    assert census["predicted_rules"] == 5
    assert census["experimental_rules"] == 1
    assert sorted(r.element for r in rows) == ["A", "B", "C", "D", "E"]
    assert {r.gene for r in rows} == {"G1", "G2"}
    assert next(r for r in rows if r.element == "C").midpoint == 11000


def test_one_gene_per_element_is_verified_and_raises_when_it_fails(synthetic) -> None:
    rows, _ = synthetic
    assert lg.verify_one_gene_per_element(rows) == {"distinct_elements": 5, "rules": 5}
    doubled = rows + [lg.RuleRow("A", "G9", "L1", "chr21", 1000, 1000)]
    with pytest.raises(AssertionError, match="more than one target gene"):
        lg.verify_one_gene_per_element(doubled)


def test_control_one_on_the_synthetic_population(synthetic) -> None:
    rows, _ = synthetic
    assert lg.control_one(rows) == pytest.approx(0.6**2 + 0.4**2)


def test_vector_bins_agrees_elementwise_with_the_scalar_rule() -> None:
    np = pytest.importorskip("numpy")
    probes = [
        0,
        1,
        999,
        1000,
        1001,
        9999,
        10000,
        10001,
        100000,
        999999,
        1000000,
        5_000_000,
        10_000_000,
        46_000_000,
    ]
    assert list(lg.vector_bins(np.array(probes, dtype=float))) == [lg.distance_bin(float(p)) for p in probes]


def test_aggregate_pairs_counts_the_synthetic_pairs_by_hand(synthetic) -> None:
    pytest.importorskip("numpy")
    rows, _ = synthetic
    counts = lg.aggregate_pairs(rows)
    assert counts.genes == ["G1", "G2"]
    within_t = counts.within_total.sum(axis=0)
    within_c = counts.within_concordant.sum(axis=0)
    assert list(within_t) == [1, 2, 0, 1, 0, 0]
    assert list(within_c) == [1, 0, 0, 0, 0, 0]
    diff_t = counts.diff_total.sum(axis=1)
    diff_c = counts.diff_concordant.sum(axis=1)
    assert list(diff_t) == [2, 1, 0, 3, 0, 0]
    assert list(diff_c) == [2, 0, 0, 1, 0, 0]
    # every unordered pair of five elements is counted exactly once
    assert int(within_t.sum() + diff_t.sum()) == 10


def test_the_buffer_flush_does_not_change_a_single_count(synthetic) -> None:
    pytest.importorskip("numpy")
    rows, _ = synthetic
    whole = lg.aggregate_pairs(rows)
    original = lg.FLUSH_AT
    try:
        lg.FLUSH_AT = 1
        flushed = lg.aggregate_pairs(rows)
    finally:
        lg.FLUSH_AT = original
    assert list(flushed.within_total.sum(axis=0)) == list(whole.within_total.sum(axis=0))
    assert list(flushed.diff_concordant.sum(axis=1)) == list(whole.diff_concordant.sum(axis=1))


def test_contributing_bins_drops_a_thin_bin_and_reports_the_dropped_weight(synthetic) -> None:
    pytest.importorskip("numpy")
    rows, _ = synthetic
    bins = lg.contributing_bins(lg.aggregate_pairs(rows))
    # every synthetic bin is far below the minimum of 10 pairs, so none contributes
    assert bins["contributing"] == []
    assert bins["dropped_within_gene_weight"] == 1.0
    assert bins["inconclusive_by_rule"] is True


def test_point_estimates_match_the_hand_count(synthetic) -> None:
    pytest.importorskip("numpy")
    rows, _ = synthetic
    counts = lg.aggregate_pairs(rows)
    point = lg.point_estimates(counts, keep=[0, 1, 3])
    assert point["within_gene_pairs"] == 4
    assert point["within_gene_concordant_pairs"] == 1
    assert point["within_gene_concordance_primary"] == pytest.approx(0.25)
    assert point["different_gene_pairs"] == 6
    assert point["different_gene_concordant_pairs"] == 3
    # weights over the contributing bins 0, 1 and 3 are 1/4, 2/4 and 1/4 of the within-gene pairs
    assert point["within_gene_concordance_matched"] == pytest.approx(0.25 * 1.0)
    assert point["different_gene_concordance_matched"] == pytest.approx(
        0.25 * 1.0 + 0.5 * 0.0 + 0.25 * (1 / 3)
    )


def test_a_bin_with_no_pairs_in_an_arm_takes_weight_zero_and_the_rest_renormalise() -> None:
    within_c = [1, 0, 0, 0, 0, 0]
    within_t = [2, 2, 0, 0, 0, 0]
    diff_c = [1, 0, 0, 0, 0, 0]
    diff_t = [2, 0, 0, 0, 0, 0]
    within, different, used = lg._matched(within_c, within_t, diff_c, diff_t, keep=[0, 1])
    assert within == pytest.approx(0.5)
    assert different == pytest.approx(0.5)
    assert used == pytest.approx(0.5)


def test_the_dyadic_resample_weights_a_different_gene_pair_by_both_its_genes(synthetic) -> None:
    pytest.importorskip("numpy")
    rows, _ = synthetic
    counts = lg.aggregate_pairs(rows)
    draws = lg.resample(counts, keep=[0, 1, 3], draws=50, seed=3)
    assert draws["clusters"] == 2
    assert draws["met_minimum"] is False  # 2 genes is below the floor of 10
    assert draws["resamples_used"] <= 50
    assert "PRODUCT of its two genes'" in draws["cluster_rule"]


def test_concentration_is_printed_with_no_threshold(synthetic) -> None:
    pytest.importorskip("numpy")
    rows, _ = synthetic
    conc = lg.concentration(lg.aggregate_pairs(rows), top=1)
    assert conc["within_gene_pairs"] == 4
    assert conc["top_genes"] == ["G1"]
    assert conc["share_of_pairs_from_the_top_genes"] == pytest.approx(0.75)
    assert conc["genes_with_at_least_one_pair"] == 2
    assert "no reading turns on it" in conc["no_threshold"]


def test_the_equal_weight_sensitivity_excludes_units_with_no_pairs(synthetic) -> None:
    pytest.importorskip("numpy")
    rows, _ = synthetic
    equal = lg.gene_equal_weight(lg.aggregate_pairs(rows), keep=[0, 1, 3])
    assert equal["labelled"].startswith("SENSITIVITY")
    assert equal["genes_contributing"] == 2
    assert equal["gene_pairs_contributing"] == 1
    assert equal["genes_excluded_no_pairs_in_a_contributing_bin"] == 0
    # G1 has pairs in bins 0 and 1 only, G2 in bin 3 only: each reweights over its own bins
    assert equal["within_gene_equal_weight"] == pytest.approx(((0.25 * 1.0) / 0.75 + 0.0) / 2)


def test_a_ratio_interval_carries_its_band_assumptions_and_identical_share() -> None:
    block = lg.ratio_interval([0.5] * 100, 0.25, "flat / control")
    assert block["ci95"] == [2.0, 2.0]
    assert block["band"] == "far_above"
    assert block["identical_resample_share"] == 1.0
    assert block["degenerate"] is True
    assert "GENES are the clusters" in block["instrument_assumptions"]
    assert "held fixed through every bootstrap resample" in block["instrument_assumptions"]
    with pytest.raises(ValueError):
        lg.ratio_interval([0.5], 0.0, "bad")
    with pytest.raises(ValueError, match="same resamples"):
        lg.ratio_interval([0.5, 0.6], [0.5], "mismatched")


def test_a_ratio_of_two_arms_uses_the_same_resamples() -> None:
    block = lg.ratio_interval([0.6] * 50, [0.3] * 50, "arm / arm")
    assert block["ci95"] == [2.0, 2.0]
    assert "same gene resample" in block["instrument_assumptions"].lower()


def test_the_gtex_normalisation_maps_a_compiled_label_onto_a_column_name() -> None:
    assert lg.normalise_tissue("Small_Intestine_Terminal_Ileum") == lg.normalise_tissue(
        "Small Intestine - Terminal Ileum"
    )
