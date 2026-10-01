# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registered parts of the ENCODE-rE2G paired comparison, on synthetic inputs only.

No comparator score and no benchmark table is read here. What is checked is that the registration says
what it must say, that the gate refuses as registered, that the reading cannot be worded any other way,
and that the paired statistic is actually paired.
"""

from __future__ import annotations

import random

import pytest

from genomeos.attribution import crispri, re2g


def pair(chrom: str, start: int, regulated: bool, weight: float = 1.0) -> crispri.Pair:
    return crispri.Pair(
        chrom=chrom,
        start=start,
        end=start + 500,
        gene=f"G{start}",
        cell="K562",
        dataset="synthetic",
        distance=10_000.0,
        dhs=1.0,
        h3k27ac=1.0,
        regulated=regulated,
        weight=weight,
    )


def synthetic(n_chrom: int = 8, per_chrom: int = 12) -> list[crispri.Pair]:
    rng = random.Random(7)
    out = []
    for c in range(n_chrom):
        for i in range(per_chrom):
            out.append(
                pair(f"chr{c + 1}", 1000 + i * 2000, rng.random() < 0.3, weight=0.5 + rng.random() / 2)
            )
    return out


# --- the gate ---------------------------------------------------------------------------------------


def test_gate_target_is_the_published_figure_and_the_tolerance_is_the_registered_one():
    assert pytest.approx(0.556151) == re2g.GATE_TARGET
    assert re2g.GATE_INTERVAL == (0.467852, 0.631224)
    assert re2g.GATE_TOLERANCE == 0.02


def test_gate_passes_only_inside_the_tolerance():
    assert re2g.gate(re2g.GATE_TARGET)["passed"] is True
    assert re2g.gate(re2g.GATE_TARGET + 0.02)["passed"] is True
    assert re2g.gate(re2g.GATE_TARGET - 0.02)["passed"] is True


def test_gate_stops_on_a_miss_larger_than_the_tolerance_in_either_direction():
    for value in (re2g.GATE_TARGET + 0.0201, re2g.GATE_TARGET - 0.0201, 0.0, 1.0):
        verdict = re2g.gate(value)
        assert verdict["passed"] is False, value
        assert verdict["may_compare"] is False, value
        assert "not the published ones" in verdict["reason"]


def test_a_gate_that_could_not_be_run_is_not_a_pass():
    verdict = re2g.gate(None)
    assert verdict["passed"] is False
    assert verdict["may_compare"] is False
    assert verdict["reproduced"] is None and verdict["miss"] is None


def test_the_gate_records_the_miss_it_refuses_on():
    assert re2g.gate(0.3)["miss"] == pytest.approx(round(abs(0.3 - re2g.GATE_TARGET), 4))


# --- the reading ------------------------------------------------------------------------------------


def test_the_reading_is_the_three_registered_strings_and_nothing_else():
    assert re2g.reading([0.01, 0.2]) == "ranks better than ENCODE-rE2G on these pairs"
    assert re2g.reading([-0.05, 0.2]) == "no difference detected"
    assert re2g.reading([-0.3, -0.01]) == "ranks worse"


def test_an_interval_touching_zero_reads_as_no_difference_detected():
    assert re2g.reading([0.0, 0.2]) == re2g.NO_DIFFERENCE
    assert re2g.reading([-0.2, 0.0]) == re2g.NO_DIFFERENCE


def test_no_interval_means_no_reading_rather_than_a_borrowed_one():
    assert re2g.reading(None) is None
    assert re2g.reading([]) is None
    assert re2g.reading([0.1]) is None


# --- the statistic ----------------------------------------------------------------------------------


def test_paired_delta_reports_clusters_draws_requested_and_draws_dropped():
    pairs = synthetic()
    rng = random.Random(1)
    ours = [rng.random() for _ in pairs]
    theirs = [rng.random() for _ in pairs]
    out = re2g.paired_delta(ours, theirs, pairs, draws=50)
    for field in ("clusters", "draws_requested", "draws_dropped", "resamples", "met_minimum"):
        assert field in out, field
    assert out["clusters"] == len({p.chrom for p in pairs})
    assert out["draws_requested"] == 50
    assert out["draws_dropped"] == 50 - out["resamples"]


def test_paired_delta_names_its_seed_and_carries_the_registered_reading():
    pairs = synthetic()
    rng = random.Random(2)
    ours = [rng.random() for _ in pairs]
    theirs = [rng.random() for _ in pairs]
    out = re2g.paired_delta(ours, theirs, pairs, draws=40)
    assert out["seed"] == re2g.SEED
    assert "delta_auprc" in out and "gain" not in out
    assert out["reading"] in (re2g.READS_BETTER, re2g.NO_DIFFERENCE, re2g.READS_WORSE, None)


def test_the_delta_is_paired_the_same_draws_serve_both_models():
    """A model scored against itself must give a delta of exactly zero in every draw.

    That can only happen if each draw resamples one index set and scores both models on it. Two
    independent resamples would give a non-degenerate interval here, so this fails on an unpaired
    bootstrap.
    """
    pairs = synthetic()
    rng = random.Random(3)
    same = [rng.random() for _ in pairs]
    out = re2g.paired_delta(same, list(same), pairs, draws=60)
    assert out["delta_auprc"] == 0.0
    assert out["ci95"] == [0.0, 0.0]
    assert out["reading"] == re2g.NO_DIFFERENCE


def test_swapping_the_two_models_negates_the_delta_on_the_same_seed():
    pairs = synthetic()
    rng = random.Random(4)
    ours = [rng.random() for _ in pairs]
    theirs = [rng.random() for _ in pairs]
    a = re2g.paired_delta(ours, theirs, pairs, draws=60)
    b = re2g.paired_delta(theirs, ours, pairs, draws=60)
    assert a["delta_auprc"] == pytest.approx(-b["delta_auprc"], abs=1e-4)
    assert a["ci95"][0] == pytest.approx(-b["ci95"][1], abs=1e-4)
    assert a["ci95"][1] == pytest.approx(-b["ci95"][0], abs=1e-4)


def test_the_registered_draw_count_is_at_least_the_minimum_the_project_requires():
    assert re2g.DRAWS >= 2000
    assert re2g.DRAWS >= crispri.MIN_RESAMPLES


# --- strata without the feature ---------------------------------------------------------------------


def test_the_three_strata_without_the_deletion_feature_report_null_and_the_reason():
    for cell in ("HCT116", "Jurkat", "WTC11"):
        out = re2g.delta_where_available(cell, lambda: {"gain": 0.123, "ci95": [0.1, 0.2]})
        assert out["gain"] is None, cell
        assert out["ci95"] is None, cell
        assert out["unavailable"] == crispri.UNAVAILABLE_GAIN, cell


def test_a_stratum_with_the_feature_gets_its_number():
    out = re2g.delta_where_available("K562", lambda: {"gain": 0.123, "ci95": [0.1, 0.2]})
    assert out == {"gain": 0.123, "ci95": [0.1, 0.2]}


def test_the_unavailable_strata_are_exactly_the_three_registered_ones():
    assert re2g.FEATURE_UNAVAILABLE == ("HCT116", "Jurkat", "WTC11")


# --- the registration itself ------------------------------------------------------------------------


def test_the_populations_carry_a_count_beside_every_name():
    for key, pop in re2g.POPULATIONS.items():
        assert pop["name"] and isinstance(pop["pairs"], int), key
        assert isinstance(pop["positives"], int) and pop["positives"] > 0, key
        assert pop["weighted_positives"] > 0, key
    assert re2g.POPULATIONS["primary"]["pairs"] == 1918
    assert re2g.POPULATIONS["secondary"]["pairs"] == 4378
    assert re2g.POPULATIONS["gm12878"]["pairs"] == 68


def test_the_pooled_population_matches_the_published_pair_and_positive_counts():
    pooled = re2g.POPULATIONS["secondary"]
    assert (pooled["pairs"], pooled["positives"]) == (4378, 190)
    assert pooled["weighted_positives"] == pytest.approx(157.39, abs=0.01)


def test_all_four_independence_statements_are_present_and_non_empty():
    assert len(re2g.INDEPENDENCE) == 4
    assert set(re2g.INDEPENDENCE) == {
        "a_comparator_trained_on_this_compendium",
        "b_our_weights_frozen",
        "c_reused_benchmark",
        "d_deletion_feature_exposure",
    }
    for key, text in re2g.INDEPENDENCE.items():
        assert len(text) > 40, key


def test_the_deletion_feature_exposure_is_stated_as_not_established():
    assert "not established" in re2g.INDEPENDENCE["d_deletion_feature_exposure"]


def test_the_join_rule_states_every_part_it_must_state():
    assert set(re2g.JOIN_RULE) == {"key", "overlap", "several_overlaps", "unpredicted", "biosample"}
    assert "sum" in re2g.JOIN_RULE["several_overlaps"]
    assert "fill_value, 0" in re2g.JOIN_RULE["unpredicted"]


def test_the_comparator_preference_order_is_kept_and_preference_one_is_recorded_absent():
    assert len(re2g.COMPARATOR_PREFERENCE) == 2
    assert re2g.COMPARATOR_PREFERENCE[0].startswith("1.")
    assert re2g.COMPARATOR_PREFERENCE[1].startswith("2.")
    assert re2g.COMPARATOR_IN_TABLE.startswith("absent:")


def test_the_falsifier_withdraws_the_old_wording_whatever_the_outcome():
    assert "either way" in re2g.FALSIFIER


# --- the comparator and the budget ------------------------------------------------------------------


def test_a_comparator_file_is_named_for_every_benchmark_cell_type():
    assert set(re2g.COMPARATOR_FILES) == {"K562", "GM12878", "HCT116", "Jurkat", "WTC11"}
    for cell, entry in re2g.COMPARATOR_FILES.items():
        assert entry["annotation"].startswith("ENCSR"), cell
        assert entry["file"].startswith("ENCFF"), cell
        assert entry["bytes"] > 0, cell


def test_the_bytes_the_pooled_gate_needs_are_summed_from_the_named_files():
    assert sum(v["bytes"] for v in re2g.COMPARATOR_FILES.values()) == re2g.COMPARATOR_BYTES_REQUIRED


def test_the_pooled_gate_costs_more_than_the_authorised_budget():
    """The reason the comparison is registered and not run, kept as an assertion rather than prose.

    If a later lane is given a larger budget it must change DOWNLOAD_BUDGET_BYTES deliberately, and this
    test is where that decision becomes visible.
    """
    assert re2g.COMPARATOR_BYTES_REQUIRED > re2g.DOWNLOAD_BUDGET_BYTES
    assert re2g.COMPARATOR_FILES["K562"]["bytes"] < re2g.DOWNLOAD_BUDGET_BYTES


def test_the_blocker_states_that_no_file_is_behind_a_login_or_a_payment():
    assert "not_a_paywall" in re2g.BLOCKER
    assert "login" in re2g.BLOCKER["not_a_paywall"]


def test_the_blocker_names_the_pooled_only_gate_and_the_ambiguous_hct116_sample():
    assert "no per-cell-type held-out AUPRC" in re2g.BLOCKER["the_gate_is_pooled_only"]
    assert "16 HCT116" in re2g.BLOCKER["hct116_is_ambiguous"]


def test_the_join_rule_is_pinned_to_the_commit_the_paper_names_not_to_a_branch():
    assert re2g.BENCHMARK_TAG == "v1.0.0"
    assert len(re2g.BENCHMARK_COMMIT) == 40
    assert re2g.BENCHMARK_COMMIT == "50587422e6b11259ead6fbc6f867681c788f39b7"


def test_a_k562_only_file_could_pass_the_gate_for_the_wrong_reason():
    """The bracket is the real reason the one affordable file cannot serve the gate.

    The first version of this blocker said a K562-only file must miss the pooled figure by more than the
    tolerance. That was wrong: with 36.0% of the positive weight forced into one tie at the bottom, the
    attainable pooled figure runs from no skill inside K562 to perfect separation inside K562, and the
    published target sits inside that span. So the gate could be passed without the scores being the
    published ones, which is worse than failing, and this test is what keeps that straight.
    """
    low, high = re2g.K562_ONLY_RANGE
    assert low < high
    assert low <= re2g.GATE_TARGET <= high
    assert "could pass it for the wrong reason" in re2g.BLOCKER["the_pooled_gate_needs_all_five_cell_types"]


def test_the_blocker_does_not_claim_the_k562_only_file_must_fail_the_gate():
    text = re2g.BLOCKER["the_pooled_gate_needs_all_five_cell_types"]
    assert "missed by far more than" not in text
    assert "whatever either model does" not in text


# --- the join, on synthetic prediction files ---------------------------------------------------------


def _prediction_file(tmp_path, rows: list[tuple[str, int, int, str, str, float]]):
    """A minimal ENCODE rE2G `element gene links` file: the real header, only the read columns filled."""
    import gzip

    header = [
        "#chr",
        "start",
        "end",
        "name",
        "class",
        "TargetGene",
        "TargetGeneEnsemblID",
        "TargetGeneTSS",
        "isSelfPromoter",
        "CellType",
        "distanceToTSS.Feature",
        "normalizedDNase_prom.Feature",
        "3DContact.Feature",
        "ABC.Score.Feature",
        "numCandidateEnhGene.Feature",
        "numTSSEnhGene.Feature",
        "sumNearbyEnhancers.Feature",
        "ubiquitousExpressedGene.Feature",
        "Score",
    ]
    path = tmp_path / "pred.bed.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("\t".join(header) + "\n")
        for chrom, start, end, gene, cell, score in rows:
            row = [""] * len(header)
            row[0], row[1], row[2] = chrom, str(start), str(end)
            row[5], row[9], row[18] = gene, cell, str(score)
            fh.write("\t".join(row) + "\n")
    return path


def test_the_overlap_test_is_the_pipelines_inclusive_one():
    assert re2g._overlaps(100, 200, 200, 300) is True  # touching at one base
    assert re2g._overlaps(100, 200, 201, 300) is False
    assert re2g._overlaps(100, 200, 50, 100) is True
    assert re2g._overlaps(100, 200, 50, 99) is False
    assert re2g._overlaps(100, 200, 120, 130) is True  # contained


def test_a_pair_joins_only_its_own_gene_and_its_own_cell_type(tmp_path):
    pairs = [pair("chr1", 1000, True)]
    pairs[0].gene = "GENE"
    rows = [
        ("chr1", 1000, 1500, "GENE", "K562", 0.9),
        ("chr1", 1000, 1500, "OTHER", "K562", 0.8),  # wrong gene
        ("chr2", 1000, 1500, "GENE", "K562", 0.7),  # wrong chromosome
    ]
    got, hits = re2g.join_predictions(pairs, _prediction_file(tmp_path, rows), "K562")
    assert got == {0: 0.9}
    assert hits == 1


def test_a_pair_of_another_cell_type_is_not_served_by_this_file(tmp_path):
    pairs = [pair("chr1", 1000, True)]
    pairs[0].gene = "GENE"
    pairs[0].cell = "WTC11"
    rows = [("chr1", 1000, 1500, "GENE", "K562", 0.9)]
    got, hits = re2g.join_predictions(pairs, _prediction_file(tmp_path, rows), "K562")
    assert got == {} and hits == 0


def test_several_overlapping_predictions_are_summed_by_default_and_maxed_on_request(tmp_path):
    pairs = [pair("chr1", 1000, True)]
    pairs[0].gene = "GENE"
    rows = [
        ("chr1", 1000, 1100, "GENE", "K562", 0.3),
        ("chr1", 1200, 1400, "GENE", "K562", 0.4),
    ]
    path = _prediction_file(tmp_path, rows)
    summed, hits = re2g.join_predictions(pairs, path, "K562", aggregate="sum")
    assert summed[0] == pytest.approx(0.7)
    assert hits == 2
    maxed, _ = re2g.join_predictions(pairs, path, "K562", aggregate="max")
    assert maxed[0] == pytest.approx(0.4)


def test_an_unpredicted_pair_is_absent_from_the_mapping_so_the_caller_fills_it_with_zero(tmp_path):
    pairs = [pair("chr1", 1000, True), pair("chr1", 90000, False)]
    for p in pairs:
        p.gene = "GENE"
    rows = [("chr1", 1000, 1500, "GENE", "K562", 0.9)]
    got, _ = re2g.join_predictions(pairs, _prediction_file(tmp_path, rows), "K562")
    assert 0 in got and 1 not in got


def test_an_unknown_aggregate_function_is_refused_rather_than_guessed(tmp_path):
    pairs = [pair("chr1", 1000, True)]
    path = _prediction_file(tmp_path, [("chr1", 1000, 1500, "GENE", "K562", 0.9)])
    with pytest.raises(ValueError, match="sum or max"):
        re2g.join_predictions(pairs, path, "K562", aggregate="mean")


# --- the second registration ------------------------------------------------------------------------


def test_the_second_registration_states_its_prior_exposure_rather_than_claiming_blindness():
    second = re2g.SECOND_REGISTRATION
    exposure = second["prior_exposure"]
    assert "NOT blind" in exposure["stated_because_it_is_real"]
    seen = exposure["heldout_pooled_weighted_already_seen"]
    assert seen["dnase + distance"] == 0.4757
    assert seen["dnase + distance + deletion"] == 0.6393
    assert exposure["heldout_k562_already_seen"]["pairs"] == 1744


def test_the_second_registration_reuses_the_three_registered_readings_unchanged():
    reading = re2g.SECOND_REGISTRATION["reading"]
    assert reading["lower_bound_above_zero"] == re2g.READS_BETTER
    assert reading["interval_covers_zero"] == re2g.NO_DIFFERENCE
    assert reading["upper_bound_below_zero"] == re2g.READS_WORSE


def test_the_second_registration_changes_our_features_and_never_the_comparators():
    models = re2g.SECOND_REGISTRATION["models"]
    assert "DNASE_FEATURES" in models["ours"]
    assert "unchanged" in models["comparator"]
    assert (
        "never theirs" in re2g.SECOND_REGISTRATION["why_it_exists"]
        or "Nothing of the comparator" in (models["comparator"])
    )


def test_the_second_registration_names_the_k562_pairs_as_its_only_registered_population():
    assert "1,918 held-out K562 pairs" in re2g.SECOND_REGISTRATION["population"]
    assert "descriptive" in re2g.SECOND_REGISTRATION["population"]


def test_the_second_registrations_falsifier_forbids_citing_the_h3k27ac_comparison_alone():
    assert "may not be cited alone" in re2g.SECOND_REGISTRATION["falsifier"]
