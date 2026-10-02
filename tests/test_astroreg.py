# SPDX-License-Identifier: AGPL-3.0-or-later
"""The astrocyte CRISPRi registration: labels, intervals, readings, power classes and the sizing."""

from __future__ import annotations

import pytest

from genomeos.attribution import astroreg, crispri


class TestSizing:
    def test_the_test_is_sized_on_the_positive_class_not_on_all_hits(self):
        """74 loci over 133 decreases, not 86 over 158 hits: the narrower rule governs."""
        assert astroreg.REGISTERED_LOCI == 74
        assert astroreg.REGISTERED_POSITIVES == 133
        assert astroreg.ELIGIBILITY_LOCI == 86
        assert astroreg.ELIGIBILITY_HITS == 158
        assert astroreg.REGISTERED_LOCI < astroreg.ELIGIBILITY_LOCI

    def test_the_hits_reconcile_with_the_positives_and_the_increases(self):
        assert astroreg.REGISTERED_POSITIVES + astroreg.REGISTERED_INCREASES_HELD_APART == (
            astroreg.ELIGIBILITY_HITS
        )

    def test_both_floors_still_clear_on_the_narrower_class(self):
        from genomeos.attribution import fresh

        assert astroreg.REGISTERED_POSITIVES >= fresh.POSITIVE_FLOOR
        assert astroreg.REGISTERED_LOCI >= fresh.LOCUS_FLOOR

    def test_the_reconciliation_says_the_wider_figure_may_not_be_quoted(self):
        assert "may not be quoted for this" in astroreg.LOCUS_RECONCILIATION
        assert "86" in astroreg.LOCUS_RECONCILIATION and "74" in astroreg.LOCUS_RECONCILIATION

    def test_the_interval_and_power_rules_are_sized_on_the_registered_loci(self):
        assert "74" in astroreg.INTERVAL_RULE
        assert "74" in astroreg.POWER_RULE
        assert "86" not in astroreg.INTERVAL_RULE
        assert "86" not in astroreg.POWER_RULE


class TestLabels:
    def test_a_significant_decrease_is_a_positive(self):
        assert astroreg.label_of(True, True, True) == "positive"

    def test_a_significant_increase_is_held_apart_and_is_never_a_negative(self):
        assert astroreg.label_of(True, False, True) == "increase_held_apart"
        assert astroreg.label_of(True, False, False) == "increase_held_apart"

    def test_a_well_powered_non_hit_is_the_primary_negative(self):
        assert astroreg.label_of(False, False, True) == "negative"

    def test_an_underpowered_non_hit_is_excluded_not_counted_as_a_negative(self):
        assert astroreg.label_of(False, False, False) == "excluded_underpowered"

    def test_the_increase_rule_is_registered_in_words(self):
        assert "never folded into the negatives" in astroreg.LABELS["increases_held_apart"]

    def test_the_negative_rule_names_the_power_convention_of_the_k562_benchmark(self):
        assert "0.25" in astroreg.LABELS["primary_negative"]
        assert "0.15" in astroreg.LABELS["sensitivity_negative"]


class TestClusterBootstrap:
    def _data(self, n_clusters, per=6):
        a, b, lab, cl = [], [], [], []
        for c in range(n_clusters):
            for i in range(per):
                pos = i % 3 == 0
                a.append(0.9 if pos else 0.1)
                b.append(0.5)
                lab.append(pos)
                cl.append(f"c{c}")
        return a, b, lab, cl

    def test_enough_clusters_gives_an_interval(self):
        out = astroreg.cluster_bootstrap(*self._data(20), n=200)
        assert "ci95" in out
        assert out["clusters"] == 20
        assert out["enough_clusters"]

    def test_too_few_clusters_gives_no_ci95_at_all(self):
        """Below the minimum there is no interval to misread as a 95% CI."""
        out = astroreg.cluster_bootstrap(*self._data(4), n=200)
        assert "ci95" not in out
        assert out["unreliable"] == "interval unreliable: 4 clusters"
        assert not out["enough_clusters"]

    def test_the_minimum_is_the_committed_one(self):
        assert crispri.MIN_CLUSTERS_FOR_AN_INTERVAL == 10
        out = astroreg.cluster_bootstrap(*self._data(9), n=100)
        assert "ci95" not in out
        assert "ci95" in astroreg.cluster_bootstrap(*self._data(10), n=200)

    def test_the_point_estimate_does_not_depend_on_the_draws(self):
        a, b, lab, cl = self._data(15)
        one = astroreg.cluster_bootstrap(a, b, lab, cl, n=50)
        two = astroreg.cluster_bootstrap(a, b, lab, cl, n=500)
        assert one["point"] == two["point"]

    def test_the_same_seed_reproduces_the_interval(self):
        a, b, lab, cl = self._data(15)
        one = astroreg.cluster_bootstrap(a, b, lab, cl, seed=7, n=200)
        two = astroreg.cluster_bootstrap(a, b, lab, cl, seed=7, n=200)
        assert one["ci95"] == two["ci95"]

    def test_draws_with_no_positive_are_dropped_and_counted(self):
        out = astroreg.cluster_bootstrap(*self._data(12), n=200)
        assert out["kept_draws"] + out["dropped_draws_with_no_positive"] == 200

    def test_a_clustering_of_one_cluster_is_refused(self):
        out = astroreg.cluster_bootstrap(*self._data(1), n=100)
        assert "ci95" not in out


class TestWider:
    def test_the_wider_interval_wins(self):
        narrow = {"ci95": [-0.01, 0.01], "clusters": 74}
        wide = {"ci95": [-0.2, 0.3], "clusters": 22}
        assert astroreg.wider(narrow, wide)["ci95"] == [-0.2, 0.3]
        assert astroreg.wider(wide, narrow)["ci95"] == [-0.2, 0.3]

    def test_an_interval_without_a_ci_cannot_win_on_width(self):
        refused = {"unreliable": "interval unreliable: 4 clusters"}
        real = {"ci95": [-0.1, 0.1]}
        assert astroreg.wider(refused, real)["ci95"] == [-0.1, 0.1]

    def test_when_none_has_a_ci_the_refusal_travels(self):
        refused = {"unreliable": "interval unreliable: 4 clusters"}
        assert "unreliable" in astroreg.wider(refused, {"unreliable": "x"})


class TestReadings:
    def test_a_lower_bound_above_zero_reads_as_detected(self):
        assert (
            astroreg.reading({"ci95": [0.01, 0.2]})
            == "the deletion gain is detected in cultured fetal astrocytes"
        )

    def test_an_interval_covering_zero_reads_as_no_gain_detected(self):
        assert astroreg.reading({"ci95": [-0.05, 0.2]}) == "no gain detected"

    def test_an_interval_covering_zero_is_never_no_effect(self):
        assert "no effect" not in astroreg.reading({"ci95": [-0.05, 0.2]})
        assert "never" in astroreg.READINGS["never"]

    def test_an_upper_bound_below_zero_reads_as_lowering_ranking(self):
        assert astroreg.reading({"ci95": [-0.3, -0.02]}) == "the deletion feature lowers ranking here"

    def test_a_refused_interval_reads_as_unreliable_and_not_as_an_outcome(self):
        out = astroreg.reading({"unreliable": "interval unreliable: 5 clusters"})
        assert out == "interval unreliable: 5 clusters"

    @pytest.mark.parametrize("lo,hi", [(0.0, 0.2), (-0.2, 0.0)])
    def test_a_bound_exactly_on_zero_is_not_read_as_detected(self, lo, hi):
        assert astroreg.reading({"ci95": [lo, hi]}) == "no gain detected"


class TestPowerClasses:
    def test_at_least_point_eight_is_confirmatory(self):
        assert astroreg.power_class(0.8)["class"] == "confirmatory"
        assert astroreg.power_class(0.95)["run"]

    def test_between_point_five_and_point_eight_is_exploratory_and_must_be_labelled(self):
        out = astroreg.power_class(0.62)
        assert out["class"] == "exploratory"
        assert out["run"]
        assert "every quote" in out["label_required"]

    def test_below_point_five_is_refused_and_not_run(self):
        out = astroreg.power_class(0.42)
        assert out["class"] == "refused"
        assert out["run"] is False

    @pytest.mark.parametrize(
        "share,klass",
        [
            (0.0, "refused"),
            (0.4999, "refused"),
            (0.5, "exploratory"),
            (0.7999, "exploratory"),
            (1.0, "confirmatory"),
        ],
    )
    def test_the_class_boundaries_are_where_they_were_registered(self, share, klass):
        assert astroreg.power_class(share)["class"] == klass


class TestRegisteredTerms:
    def test_the_cache_reuse_refusal_gives_the_selection_reason(self):
        assert "most extreme of the 371 tracks" in astroreg.CACHE_REUSE_REFUSED
        assert "selection on a correlate of the outcome" in astroreg.CACHE_REUSE_REFUSED

    def test_the_track_is_recorded_as_exposed_not_as_missing(self):
        assert "returned an astrocyte output" in astroreg.ASTROCYTE_TRACK_EXPOSED
        assert "224" in astroreg.ASTROCYTE_TRACK_EXPOSED

    def test_the_same_code_path_term_names_the_membership_gate(self):
        assert "p.cell in cells" in astroreg.SAME_CODE_PATH
        assert "membership gate" in astroreg.SAME_CODE_PATH

    def test_the_primary_claim_refits_nothing(self):
        assert "Nothing is refitted" in astroreg.PRIMARY_CLAIM

    def test_the_not_claimed_paragraph_is_carried_word_for_word(self):
        joined = " ".join(astroreg.NOT_CLAIMED)
        assert "not a replication of the K562 result" in joined
        assert "not independence of the genomic regions" in joined
        assert "not primary brain tissue" in joined
        assert "no claim at all until the comparison is registered" in joined

    def test_the_exposure_list_names_the_psychencode_selection_and_its_effect(self):
        joined = " ".join(astroreg.EXPOSURE)
        assert "flatters the activity term in both arms equally" in joined
        assert "leaves the gain comparison less affected" in joined

    def test_the_exposure_list_admits_the_chromatin_already_read(self):
        joined = " ".join(astroreg.EXPOSURE)
        assert "astrocyte CHROMATIN features have already been read here" in joined

    def test_the_build_is_settled_and_needs_no_liftover(self):
        assert "No liftover is needed" in astroreg.BUILD

    def test_the_cost_is_framed_per_positive_against_the_precedent(self):
        assert astroreg.HCT116_PRECEDENT == {"requests": 705, "covered_pairs": 363, "regulated": 34}
        assert "expected requests per" in astroreg.COST_FRAMING

    def test_terms_carries_every_registered_section(self):
        t = astroreg.terms()
        for key in (
            "feasibility_gate",
            "build",
            "labels",
            "primary_claim",
            "interval_rule",
            "readings",
            "power_rule",
            "same_code_path",
            "exposure",
            "not_claimed",
            "cost_framing",
        ):
            assert key in t, key
        assert t["requests_sent_by_this_lane"] == 0
        assert t["money_spent_by_this_lane"] == 0

    def test_terms_carries_the_locus_reconciliation(self):
        assert "locus_reconciliation" in astroreg.terms()["eligibility"]
