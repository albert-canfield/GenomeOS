# SPDX-License-Identifier: AGPL-3.0-or-later
"""The coverage analysis of N1's amendment-2 result: pure arithmetic, source locations, no measurement.

Every test here runs on committed bytes. None of them opens a study file, and the one test that touches
`data/results/` reads the committed result to check that the figures this analysis explains are still the
figures it was written against -- so a later rewrite of that result cannot leave the analysis asserting
numbers nothing holds any more.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from genomeos.attribution import n1_coverage as nc

MODULE = Path("genomeos/attribution/n1_perturb_response.py")
RESULT = Path("data/results/n1_result_amendment_2.json")

#: The responder counts the committed result records, in the order the result lists them is irrelevant.
COMMITTED_COUNTS = [27, 11, 3, 3, 3, 2, 2, 1, 1, 1, 1, 1] + [0] * 44


def _module_text() -> str:
    if not MODULE.exists():
        pytest.skip(f"{MODULE} is absent")
    return MODULE.read_text()


def _result() -> dict:
    if not RESULT.exists():
        pytest.skip(f"{RESULT} is absent")
    return json.loads(RESULT.read_text())


# --- the concentration finding ----------------------------------------------------------------------


def test_concentration_on_the_committed_counts_refutes_an_even_spread():
    c = nc.concentration(COMMITTED_COUNTS)
    assert c["factors"] == 56
    assert c["responder_cells"] == 56
    assert c["zero_factors_observed"] == 44
    assert c["defined_observed"] == 12
    assert c["per_factor_mean"] == pytest.approx(1.0)
    assert c["per_factor_variance_ddof1"] == pytest.approx(15.1636, abs=1e-4)
    assert c["variance_to_mean_ratio"] == pytest.approx(15.1636, abs=1e-4)
    # The headline: an even spread at the SAME mean clears the floor, so the shortfall is not density.
    assert c["even_spread"]["expected_defined"] == pytest.approx(35.40, abs=0.01)
    assert c["even_spread"]["clears_the_floor"] is True
    assert c["even_spread"]["expected_defined"] > c["floor"]
    # And the observed zero count is far from that model, so the committed counts refute it.
    assert c["even_spread"]["observed_zeros_in_sd"] == pytest.approx(6.5, abs=0.1)


def test_concentration_of_a_perfectly_even_spread_has_no_overdispersion():
    c = nc.concentration([1] * 56)
    assert c["per_factor_mean"] == pytest.approx(1.0)
    assert c["per_factor_variance_ddof1"] == pytest.approx(0.0)
    assert c["variance_to_mean_ratio"] == pytest.approx(0.0)
    assert c["zero_factors_observed"] == 0
    # Same mean as the committed counts, so the even-spread expectation is identical: the mean alone
    # does not distinguish the two, which is exactly why the variance is the finding.
    assert c["even_spread"]["expected_defined"] == pytest.approx(35.40, abs=0.01)


# --- the gamma-mixed extrapolation, and the bracket it may not be quoted without --------------------


def test_gamma_mixed_shape_is_none_without_overdispersion():
    assert nc.gamma_mixed_shape(1.0, 1.0) is None
    assert nc.gamma_mixed_shape(1.0, 0.5) is None
    assert nc.gamma_mixed_shape(1.0, 2.0) == pytest.approx(1.0)


def test_gamma_mixed_zero_probability_matches_the_geometric_case():
    # k = 1 is the geometric distribution: P(0) = 1 / (1 + mu).
    for mu in (0.5, 1.0, 4.0):
        assert nc.gamma_mixed_zero_probability(mu, 1.0) == pytest.approx(1 / (1 + mu))


def test_fitted_extrapolation_falls_short_of_the_floor_at_every_universe_the_study_holds():
    out = nc.gamma_mixed_universe(COMMITTED_COUNTS, 496, [496, 860, 8248])
    assert out["fits"] is True
    assert out["shape_k"] == pytest.approx(0.070604, abs=1e-6)
    assert out["mixing_coefficient_of_variation"] == pytest.approx(3.763, abs=1e-3)
    got = {r["universe_genes"]: r["expected_defined"] for r in out["at_universes"]}
    assert got[496] == pytest.approx(9.8, abs=0.05)
    assert got[860] == pytest.approx(11.5, abs=0.05)
    assert got[8248] == pytest.approx(17.9, abs=0.05)
    # The claim the result file rests on: short of the floor even at every gene the pseudobulk holds.
    assert all(v < nc.FLOOR for v in got.values())
    # And the fit is honest about the universe it was fitted on: 9.8 against the 12 observed.
    assert got[496] < 12


def test_the_universe_the_fit_would_need_for_the_floor_exceeds_any_genome():
    out = nc.gamma_mixed_universe_for_defined(COMMITTED_COUNTS, 496, nc.FLOOR)
    assert out["solvable"] is True
    assert out["universe_genes_needed"] > 1e6


def test_exact_upper_rate_is_the_binomial_bound_and_not_an_approximation():
    # One gene, no responder: any rate up to 0.95 gives no responder with probability at least 0.05.
    assert nc.exact_upper_rate(1, 0, 0.05) == pytest.approx(0.95, abs=1e-9)
    r = nc.exact_upper_rate(496, 0, 0.05)
    assert (1 - r) ** 496 == pytest.approx(0.05, abs=1e-6)
    assert r == pytest.approx(0.006022, abs=1e-6)
    with pytest.raises(NotImplementedError):
        nc.exact_upper_rate(496, 1, 0.05)


def test_the_bracket_carries_both_ends_including_the_one_that_is_not_excluded():
    b = nc.distribution_free_bracket(COMMITTED_COUNTS, 496, 364)
    assert b["low_end"]["defined_at_any_universe"] == 12
    assert b["high_end"]["upper_per_gene_rate"] == pytest.approx(0.006022, abs=1e-6)
    assert b["high_end"]["responders_per_observed_universe_at_that_rate"] == pytest.approx(2.99, abs=0.01)
    assert b["high_end"]["probability_at_least_one_responder_among_the_added_genes"] == pytest.approx(
        0.889, abs=1e-3
    )
    assert "unlucky" in b["what_separates_them"]


# --- the endpoint, answered two ways --------------------------------------------------------------


def test_the_endpoint_is_named_by_the_result_s_own_keys():
    assert nc.endpoint_from_result_keys({"genes_missing": 0})["endpoint"] == "published_anderson_darling"
    assert nc.endpoint_from_result_keys({"genes_nonfinite_response": 0})["endpoint"] == "derived_t_statistic"
    both = nc.endpoint_from_result_keys({"genes_missing": 0, "genes_nonfinite_response": 0})
    assert both["endpoint"] == "ambiguous"
    assert nc.endpoint_from_result_keys({"factor": "NRF1"})["endpoint"] == "ambiguous"


def test_the_committed_result_s_keys_name_the_published_p_value_path():
    result = _result()
    record = result["factors"][0]
    out = nc.endpoint_from_result_keys(record)
    assert out["endpoint"] == "published_anderson_darling"
    assert out["keys_present"] == ["genes_missing"]
    assert out["keys_absent"] == ["genes_nonfinite_response"]
    # The same result's own read log: no normalized pseudobulk and so no perturbed row of it.
    assert not any("normalized" in s for s in result["datasets_read"])


def test_every_source_site_the_analysis_cites_occurs_exactly_once():
    sites = nc.source_sites(_module_text())
    for name, lines in sites.items():
        assert len(lines) == 1, f"{name}: expected one occurrence, found {lines}"
    # The predicate forms the label from a published adjusted p-value against the producers' level.
    assert "AD_LEVEL" in nc.SITES["responder_predicate"]
    assert nc.SITES["ad_level_constant"].endswith("0.05")


def test_the_derived_statistic_is_absent_from_run_v3_and_the_check_can_see_it_when_present():
    text = _module_text()
    absent = nc.t_is_absent_from(text, "run_v3", ["RESPONDS_T", "pooled_t("])
    assert absent["absent"] is True
    assert absent["span"] is not None
    assert absent["marker_lines_inside_the_span"] == {}
    # Every marker occurs somewhere, so the absence above is a fact about run_v3's span and not about
    # a needle that never matches anything.
    assert all(absent["marker_lines"][m] for m in ("RESPONDS_T", "pooled_t("))
    # Positive control: the amendment-1 analysis DOES pool the statistic, and the same check says so.
    present = nc.t_is_absent_from(text, "_analyse_v2", ["pooled_t("])
    assert present["absent"] is False
    assert present["marker_lines_inside_the_span"]["pooled_t("]


def test_function_span_returns_none_for_a_name_the_module_does_not_define():
    assert nc.function_span(_module_text(), "no_such_function") is None


# --- what the label conflates, and the field that would decompose it -------------------------------


def test_the_missing_field_is_named_with_its_cost():
    out = nc.missing_decomposition_field()
    assert out["field"] == "responders_outside_universe"
    assert len(out["would_separate"]) == 2
    assert "outcome blindness" in out["cost_to_form"]


def test_the_committed_result_still_holds_the_figures_this_analysis_explains():
    result = _result()
    analysis = result["analysis"]
    cov = analysis["coverage"]
    assert result["status"] == "insufficient_coverage"
    assert cov == {
        "factors_analysed": 56,
        "defined": 12,
        "undefined": 44,
        "undefined_reasons": {"no_responders": 44},
        "floor": 30,
    }
    assert cov["floor"] == nc.FLOOR
    # The three sibling causes the one label is NOT hiding.
    assert cov["undefined_reasons"].get("all_responders", 0) == 0
    assert cov["undefined_reasons"].get("no_stratum_with_both_classes", 0) == 0
    assert {r["genes_missing"] for r in result["factors"]} == {0}
    assert {r["genes"] for r in result["factors"]} == {496}
    assert result["coverage"]["universe"] == 860
    assert result["coverage"]["in_file_rows"] == 496
    assert result["coverage"]["absent_from_file"] == 364
    assert result["coverage"]["listed_twice"] == 0
    # The counts the arithmetic above is run on are the result's own.
    assert sorted((r["responders"] for r in result["factors"]), reverse=True) == COMMITTED_COUNTS
    # The floor is not met and no estimate is reported: nothing here recovers one.
    assert analysis["estimate"] is None
    assert analysis["criteria"] is None
    assert analysis["secondary_unstratified"]["factors"] == 12
    assert math.isclose(analysis["secondary_unstratified"]["estimate"], -0.06105263899814837)
