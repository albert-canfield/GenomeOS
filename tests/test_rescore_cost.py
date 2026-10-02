# SPDX-License-Identifier: AGPL-3.0-or-later
"""The costing's rules, pinned: the window boundary, the request unit, the ceiling and the no-go wording."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import cell2, fresh
from genomeos.attribution import direction_link as dl
from genomeos.attribution import increase_links as il
from genomeos.attribution import reptest as rt
from genomeos.attribution import rescore_cost as rc
from tests.committed_data import must_be_committed

RESULT_RELATIVE = "data/results/rescore_cost.json"
RESULT = Path(RESULT_RELATIVE)


# ---- the window, which is the binding constraint --------------------------------------------------


def test_the_sweep_already_asks_at_the_services_maximum():
    assert rc.MAX_SUPPORTED_WINDOW == max(rc.SUPPORTED_SEQUENCE_LENGTHS) == 1_048_576
    assert rc.WINDOW_USED_BY_THE_SWEEP == rc.MAX_SUPPORTED_WINDOW


def test_window_bounds_centre_the_width_the_way_the_committed_chain_does():
    # the variant's reference interval is [start - 1, end), and resize centres the width on its centre
    start, end, width = 1_000_001, 1_000_501, 1_048_576
    centre = (start - 1 + end) // 2
    assert rc.window_centre(start, end) == centre
    lo, hi = rc.window_bounds(start, end, width)
    assert (lo, hi) == (centre - width // 2, centre + width // 2)
    assert hi - lo == width


def test_window_needed_is_the_smallest_width_that_would_contain_the_tss():
    start, end = 1_000_001, 1_000_501
    centre = rc.window_centre(start, end)
    for offset in (-700_000, -524_288, -1, 0, 1, 524_287, 900_000):
        needed = rc.window_needed(start, end, centre + offset)
        lo, hi = rc.window_bounds(start, end, needed)
        assert lo <= centre + offset < hi, offset
        smaller = needed - 2
        if smaller > 0:
            lo2, hi2 = rc.window_bounds(start, end, smaller)
            assert not (lo2 <= centre + offset < hi2), offset


def test_a_width_over_the_maximum_has_no_supported_window_and_that_is_the_refusal():
    assert rc.smallest_supported_window(1) == 16_384
    assert rc.smallest_supported_window(1_048_576) == 1_048_576
    assert rc.smallest_supported_window(1_048_578) is None
    assert rc.smallest_supported_window(3_000_000) is None


def test_the_tss_inside_the_scored_window_needs_no_wider_window_than_the_sweep_used():
    start, end = 1_000_001, 1_000_501
    centre = rc.window_centre(start, end)
    needed = rc.window_needed(start, end, centre - 374_345)
    assert needed <= rc.MAX_SUPPORTED_WINDOW
    assert rc.smallest_supported_window(needed) == rc.WINDOW_USED_BY_THE_SWEEP


# ---- purchasability: exhaustive, disjoint, and in the order that makes it so ----------------------


def test_purchasability_is_one_of_three_and_the_window_refusal_comes_first():
    over = 3_000_000
    within = 100_000
    assert rc.purchasability(True, over, False) == rc.UNPURCHASABLE_BEYOND_MAX_WINDOW
    assert rc.purchasability(False, over, False) == rc.UNPURCHASABLE_BEYOND_MAX_WINDOW
    assert rc.purchasability(True, within, False) == rc.PURCHASABLE_CELL_TRACK
    assert rc.purchasability(False, within, False) == rc.NOT_A_PURCHASE_NAMING
    # a cell the sweep already keeps is not a cell purchase, whatever else holds
    assert rc.purchasability(True, within, True) == rc.NOT_A_PURCHASE_NAMING
    for n in (over, within):
        for named in (True, False):
            for cached in (True, False):
                assert rc.purchasability(named, n, cached) in rc.PURCHASABILITY


# ---- the request unit ----------------------------------------------------------------------------


def test_the_request_unit_is_per_element_so_two_links_on_one_element_are_one_request():
    assert rc.requests_for(["E1", "E1", "E2"]) == 2
    assert rc.requests_for([]) == 0
    assert "PER ELEMENT" in rc.REQUEST_UNIT
    assert "NOT per element per cell" in rc.REQUEST_UNIT


def test_the_request_unit_cites_where_each_clause_was_read():
    for clause, citation in rc.REQUEST_UNIT_READ_FROM.items():
        assert citation.startswith(("genomeos/", "scripts/")), clause
    assert any("crispri_hct116" in c for c in rc.REQUEST_UNIT_READ_FROM.values())
    assert "AstroREG" in rc.REQUEST_UNIT_IS_CHECKED_NOT_CARRIED_OVER


# ---- the population rules, which may not be the shortfall -----------------------------------------


def test_the_cells_not_covered_are_computed_and_include_one_the_constant_misses():
    seen = ["K562", "GM12878", "HCT116", "WTC11", "Jurkat"]
    assert rc.cells_not_covered(seen) == ["HCT116", "Jurkat", "WTC11"]
    assert "Jurkat" not in dl.CELLS_NOT_CACHED
    assert "Jurkat" in rc.cells_not_covered(seen)


def test_the_population_rules_name_no_count_no_arm_and_no_floor():
    for rule in (rc.POPULATION_RULE_CELL_ARM, rc.POPULATION_RULE_WINDOW_ARM):
        low = rule.lower()
        for forbidden in ("27", "18", "9 links", "short of", "floor of"):
            assert forbidden not in low, rule
    assert "never read" in rc.FOLLOWS_A_NO_GO_WHOSE_OUTCOME_WAS_NEVER_READ.lower()
    assert "shortfall is not the population" in rc.THE_SHORTFALL_MAY_NOT_BE_BOUGHT


# ---- the floors are imported objects and are not moved -------------------------------------------


def test_the_floors_are_the_imported_objects_themselves():
    assert rc.POSITIVE_FLOOR is rt.POSITIVE_FLOOR is fresh.POSITIVE_FLOOR
    assert rc.LOCUS_FLOOR is rt.LOCUS_FLOOR is fresh.LOCUS_FLOOR is cell2.POOLED_LOCUS_FLOOR


def test_the_module_defines_no_floor_of_its_own():
    text = Path("genomeos/attribution/rescore_cost.py").read_text()
    assert "POSITIVE_FLOOR = rt.POSITIVE_FLOOR" in text
    assert "LOCUS_FLOOR = rt.LOCUS_FLOOR" in text
    assert "FLOOR = 30" not in text and "FLOOR = 20" not in text


# ---- the ceiling ---------------------------------------------------------------------------------


def test_a_ceiling_short_of_either_floor_reads_as_a_no_go_with_the_margin():
    v = rc.ceiling(21, 7, 18)
    assert v["links_at_the_ceiling"] == 28
    assert v["meets_both_floors_at_the_ceiling"] is False
    assert v["reading"] == rc.CEILING_MISSES_A_FLOOR
    assert v["short_by"] == [
        "links at the ceiling 28, short of the imported floor 30 by 2",
        "independent loci at the ceiling 18, short of the imported floor 20 by 2",
    ]
    assert "nothing to authorise" in v["what_follows"]


def test_a_ceiling_at_both_floors_reads_as_a_could_and_never_as_an_authorisation():
    v = rc.ceiling(21, 9, 20)
    assert v["meets_both_floors_at_the_ceiling"] is True
    assert v["reading"] == rc.CEILING_CLEARS_BOTH_FLOORS
    assert "Nothing here authorises" in v["reading"]
    assert v["short_by"] == []


def test_one_floor_cleared_and_one_missed_is_still_a_no_go():
    v = rc.ceiling(21, 20, 13)
    assert v["links_at_the_ceiling"] == 41
    assert v["meets_both_floors_at_the_ceiling"] is False
    assert len(v["short_by"]) == 1


def test_no_reading_softens_a_short_ceiling():
    for word in rc.FORBIDDEN_OF_A_SHORT_CEILING:
        assert word not in rc.CEILING_MISSES_A_FLOOR.split("NEVER")[0]
    assert rc.FORBIDDEN_OF_A_SHORT_CEILING == rt.FORBIDDEN_OF_A_SHORT_COUNT
    assert len(rc.CEILING_READINGS) == 2


def test_the_ceiling_says_in_words_that_it_is_not_a_count():
    assert "not a count of answerable links" in rc.CEILING_IS_NOT_A_COUNT
    assert rc.ceiling(21, 7, 18)["ceiling_is_not_a_count"] == rc.CEILING_IS_NOT_A_COUNT


# ---- requests per answerable positive ------------------------------------------------------------


def test_requests_per_answerable_positive_is_given_per_arm_and_pooled():
    f = rc.per_request_figures(24, 7, 16, 5, 3)
    ratios = f["requests_per_answerable_positive_at_the_ceiling"]
    assert ratios["increase_arm_only"] == round(24 / 7, 4)
    assert ratios["decrease_arm_only"] == 1.5
    assert ratios["both_arms_pooled"] == round(24 / 23, 4)
    assert "increase arm's" in ratios["which_one_bears_on_the_decision"]
    assert f["requests_serving_no_link_the_test_could_score"] == 3
    assert "CEILING" in f["what_these_are_not"]


def test_a_zero_denominator_gives_no_ratio_rather_than_a_made_up_one():
    f = rc.per_request_figures(24, 0, 0, 0, 24)
    ratios = f["requests_per_answerable_positive_at_the_ceiling"]
    assert ratios["increase_arm_only"] is None
    assert ratios["both_arms_pooled"] is None
    assert f["requests_per_locus_at_the_ceiling"] is None


# ---- R2, not weakened and not shadowed -----------------------------------------------------------


def test_r2_is_the_shared_wording_and_the_check_is_the_shared_one():
    assert rc.R2_WORDING is il.R2_WORDING
    assert "increase_links.check_no_mechanism_claim" in rc.PROHIBITION_IS_ENFORCED_BY
    text = Path("genomeos/attribution/rescore_cost.py").read_text()
    for word in il.FORBIDDEN_OF_AN_INCREASE_LINK:
        assert f'"{word}"' not in text, word  # no word list of its own


def test_the_registration_states_it_registers_no_purchase_and_no_test():
    reg = rc.registration()
    assert reg["lane"] == "lane-rescore"
    assert "not a purchase" in reg["what_this_is"]
    assert "BEFORE any request is sent" in reg["the_registration_that_would_spend_this_is_a_new_one"]
    assert reg["there_is_no_third"] is True
    assert reg["floors_are_imported"] == rc.FLOORS_ARE_IMPORTED
    assert "0 model requests" in reg["budget"]


# ---- the committed result's own invariants -------------------------------------------------------


@pytest.fixture(scope="module")
def result() -> dict:
    # The per-element cache this lane COUNTED is machine-local; the result it wrote is not. git
    # tracks data/results/rescore_cost.json, so the old reason named the wrong file: absence here is
    # a deleted or corrupted published result, and it FAILS naming the path.
    return json.loads(must_be_committed(RESULT_RELATIVE).read_text())


def test_the_result_spends_nothing(result):
    assert result["result_manifest"]["parameters"]["model_requests"] == 0
    assert result["result_manifest"]["parameters"]["alphagenome_requests"] == 0
    assert result["cost"]["window_population_requests"] == 0
    assert result["window_population"]["requests_available"] == 0


def test_the_results_headline_is_one_of_the_two_registered_readings(result):
    assert result["headline"] in rc.CEILING_READINGS
    assert result["headline"] == result["ceiling"]["reading"]


def test_the_results_purchasability_sums_to_each_population(result):
    assert result["cell_population"]["purchasability_sums_to_the_population"] is True
    assert sum(result["cell_population"]["purchasability"].values()) == (result["cell_population"]["links"])
    by = result["cell_population"]["by_direction"]
    assert (
        by["increase_derived"]["links"] + by["decrease_derived"]["links"]
        == (result["cell_population"]["links"])
    )


def test_the_results_request_count_is_the_distinct_element_count(result):
    detail = result["cell_population"]["links_detail"]
    assert result["cell_population"]["requests_under_the_unit"] == len({r["element"] for r in detail})
    assert result["cost"]["requests_total"] == result["cell_population"]["requests_under_the_unit"]


def test_the_result_agrees_with_the_count_it_is_costed_against(result):
    assert result["inherited"]["figures_agree"] is True
    assert (
        result["inherited"]["increase_links_re_enumerated_here"]
        == (result["inherited"]["figures_as_read"]["increase_derived_links"])
    )


def test_the_result_carries_the_admissibility_statement_and_not_in_a_footnote(result):
    assert result["follows_a_no_go_whose_outcome_was_never_read"] == (
        rc.FOLLOWS_A_NO_GO_WHOSE_OUTCOME_WAS_NEVER_READ
    )
    assert result["the_shortfall_may_not_be_bought"] == rc.THE_SHORTFALL_MAY_NOT_BE_BOUGHT


def test_the_results_precedents_were_read_and_not_typed_in(result):
    p = result["precedents"]
    assert p["HCT116 second cell type, delivered"]["requests"] == 705
    assert p["HCT116 second cell type, delivered"]["covered_pairs"] == 363
    assert p["AstroREG request plan, authorised and not sent"]["requests"] == 1322
    assert p["AstroREG request plan, authorised and not sent"]["requests_serving_no_scored_pair"] == 90
    declared = {i["path"] for i in result["result_manifest"]["inputs"]}
    assert "data/results/crispri_published_v2.json" in declared
    assert "data/results/astroreg_request_plan.json" in declared


def test_the_result_declares_the_inherited_count_and_the_benchmark_tables_as_inputs(result):
    declared = {i["path"] for i in result["result_manifest"]["inputs"]}
    assert "data/results/reptest_answerable.json" in declared
    assert any(p.startswith("data/knowledge/crispri/") for p in declared)
    assert "per_element_response_cache" in declared
    assert all(not Path(p).is_absolute() for p in declared)


def test_the_most_generous_ceiling_is_reported_and_names_the_assumption_it_makes(result):
    g = result["most_generous_ceiling"]
    assert g["rule"] == rc.MOST_GENEROUS_CEILING_RULE
    assert "amendment to a registered predicate, which this lane does not make" in g["rule"]
    assert g["links_at_the_ceiling"] >= result["ceiling"]["links_at_the_ceiling"]
    assert g["reading"] in rc.CEILING_READINGS
    # it is an upper bound on an upper bound and may not be read as a count either
    assert g["ceiling_is_not_a_count"] == rc.CEILING_IS_NOT_A_COUNT


def test_the_naming_gap_population_costs_nothing_and_is_not_counted_as_answerable(result):
    n = result["naming_gap_population"]
    assert n["requests_available"] == 0
    assert "NO REQUEST WOULD HELP" in n["cost"]
    assert "none of these is counted as answerable here" in n["it_is_not_an_answerable_link"]
    assert n["links"] == sum(n["by_direction"].values())


def test_no_link_of_the_window_population_could_be_bought_at_any_width(result):
    w = result["window_population"]
    assert w["requests_available"] == 0
    for r in w["links_detail"]:
        assert r["smallest_supported_window_that_would_contain_the_tss"] is None
        assert r["window_width_needed"] > rc.MAX_SUPPORTED_WINDOW
