# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration of the significant-increase feasibility lane, and the structural cause it rests on.

Nothing here opens a benchmark table or a result: these are assertions about the registration's own
content, about the committed extractor's behaviour on synthetic rows, and about the arithmetic of the
gate. The point of them is that a floor, a predicate or a reading cannot be moved after a count without
a test failing.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos import manifest as mf
from genomeos.attribution import cell2, fresh
from genomeos.attribution import increases as inc
from genomeos.attribution import measured as ms
from genomeos.attribution import repress2 as rp

# ---- the floors and the locus convention are imports, not choices --------------------------------


def test_both_floors_are_the_imported_ones() -> None:
    assert inc.POSITIVE_FLOOR is fresh.POSITIVE_FLOOR
    assert inc.LOCUS_FLOOR is fresh.LOCUS_FLOOR
    assert fresh.LOCUS_FLOOR == cell2.POOLED_LOCUS_FLOOR
    assert (inc.POSITIVE_FLOOR, inc.LOCUS_FLOOR) == (30, 20)


def test_the_floors_are_the_ones_lane_repress2_was_gated_against() -> None:
    assert (inc.POSITIVE_FLOOR, inc.LOCUS_FLOOR) == (rp.POSITIVE_FLOOR, rp.LOCUS_FLOOR)


def test_the_locus_convention_is_cell2s_and_carries_the_not_independence_wording() -> None:
    assert inc.LOCUS_RULE is cell2.INDEPENDENT_LOCUS_RULE
    assert inc.LOCUS_SPAN is cell2.INDEPENDENT_LOCUS_SPAN
    assert "not established biological independence" in inc.LOCUS_RULE
    assert inc.NOT_BIOLOGICAL_INDEPENDENCE is rp.NOT_BIOLOGICAL_INDEPENDENCE


def test_the_registration_names_where_each_floor_came_from() -> None:
    floors = inc.registration()["floors"]
    assert floors["links_imported_from"] == "genomeos.attribution.fresh.POSITIVE_FLOOR"
    assert "cell2.POOLED_LOCUS_FLOOR" in floors["independent_loci_imported_from"]
    assert floors["neither_chosen_here"] is True


# ---- the eligibility predicate reads the cached outcome field and nothing else -------------------


def test_is_increase_reads_the_cached_outcome_field() -> None:
    assert inc.is_increase({"outcome": ms.INCREASE}) is True
    assert inc.is_increase({"outcome": ms.DECREASE}) is False
    assert inc.is_increase({"outcome": ms.NULL_INFORMATIVE}) is False
    assert inc.is_increase({}) is False
    # an effect size of either sign cannot make a pair eligible on its own
    assert inc.is_increase({"effect_size": 0.9, "regulated": False}) is False


def test_significant_is_repress2s_predicate_itself() -> None:
    for outcome in ms.OUTCOMES:
        assert inc.significant({"outcome": outcome}) is rp.significant({"outcome": outcome})
    assert inc.SIGNIFICANT_OUTCOMES is rp.SIGNIFICANT_OUTCOMES
    assert set(inc.SIGNIFICANT_OUTCOMES) == {ms.DECREASE, ms.INCREASE}


def test_the_breakdown_is_exhaustive_at_a_denominator() -> None:
    from collections import Counter

    seen = Counter({ms.DECREASE: 3, ms.INCREASE: 1})
    out = inc.breakdown(seen)
    assert set(out) == set(ms.OUTCOMES)
    assert out[ms.NULL_INFORMATIVE] == 0
    assert sum(out.values()) == 4


# ---- the structural cause, asserted against the committed extractor ------------------------------


def _row(pairs: list[dict[str, object]]) -> dict[str, object]:
    return {"id": "e1", "measured": {"crispri": {"pairs": pairs}}}


def _pair(**over: object) -> dict[str, object]:
    base: dict[str, object] = {
        "gene": "GENE1",
        "cell": "K562",
        "regulated": False,
        "effect_size": 0.4,
        "outcome": ms.INCREASE,
        "split": ms.TRAINING,
    }
    base.update(over)
    return base


def test_a_significant_increase_raises_no_measured_rule() -> None:
    """The structural cause this lane counts around: `rule_links` emits nothing from an increase."""
    assert ms.rule_links(_row([_pair()])) == []


def test_a_significant_decrease_does_raise_one_and_its_action_is_activates() -> None:
    links = ms.rule_links(_row([_pair(regulated=True, effect_size=-0.4, outcome=ms.DECREASE)]))
    assert [(g, a) for g, a, _s, _c, _sp in links] == [("GENE1", "activates")]


def test_the_inhibits_branch_is_unreachable_through_the_committed_filter() -> None:
    """A positive effect size cannot reach the `inhibits` branch, because `regulated` gates it."""
    actions = {a for _g, a, _s, _c, _sp in ms.rule_links(_row([_pair(effect_size=0.9)]))}
    assert actions == set()


def test_the_increase_is_recorded_on_the_block_even_though_no_rule_is_raised() -> None:
    """`genes_increased` keeps it, which is why this lane can count it without changing an extractor."""
    layer = ms.Layer(
        chrom="chr1",
        crispri=[
            ms.CrispriPair(
                chrom="chr1",
                start=1000,
                end=2000,
                gene="GENE1",
                cell="K562",
                dataset="d",
                reference="r",
                regulated=False,
                significant=True,
                effect_size=0.4,
                p_adjusted=0.001,
            )
        ],
    )
    block = layer.for_element(1000, 2000)["crispri"]
    assert block["genes_increased"] == ["GENE1"]
    assert block["genes_regulated"] == []
    assert block["genes_not_regulated"] == []  # an increase is never listed as a measured null (R2)


# ---- the element-id convention cannot drift from lane-repress2's ---------------------------------


@pytest.mark.parametrize("element", ["e1", "chr21:100-200", "e1_measured", "x_measured_measured"])
def test_base_element_id_agrees_with_repress2s(element: str) -> None:
    class _R:
        element = ""

    r = _R()
    r.element = element
    assert inc.base_element_id(element) == rp._element_id(r)


# ---- the gate's arithmetic ------------------------------------------------------------------------


def _keys(n: int, span: int) -> list[cell2.LocusKey]:
    """`n` links, each `span` bases from the last, so the locus count is controllable."""
    return [cell2.LocusKey("K562", "chr1", i * span, i * span + 100, f"G{i}") for i in range(n)]


def test_a_population_below_either_floor_reads_as_the_registered_no_go() -> None:
    g = inc.gate(_keys(5, 10_000_000), "p", "call")
    assert g["links"] == 5
    assert g["independent_loci"] == 5
    assert g["meets_both_floors"] is False
    assert g["reading"] == inc.NO_GO
    assert g["short_by"] == [
        f"links 5, short of {inc.POSITIVE_FLOOR} by {inc.POSITIVE_FLOOR - 5}",
        f"independent loci 5, short of {inc.LOCUS_FLOOR} by {inc.LOCUS_FLOOR - 5}",
    ]


def test_enough_links_in_too_few_loci_is_still_a_no_go() -> None:
    # 40 links all within 1 Mb of the next: the chained rule makes them one locus
    keys = [cell2.LocusKey("K562", "chr1", i * 1000, i * 1000 + 100, f"G{i}") for i in range(40)]
    g = inc.gate(keys, "p", "call")
    assert (g["links"], g["independent_loci"]) == (40, 1)
    assert g["meets_both_floors"] is False
    assert g["short_by"] == [f"independent loci 1, short of {inc.LOCUS_FLOOR} by {inc.LOCUS_FLOOR - 1}"]


def test_both_floors_met_reads_as_the_registered_go_and_the_gate_stops_there() -> None:
    g = inc.gate(_keys(40, 10_000_000), "p", "call")
    assert (g["links"], g["independent_loci"]) == (40, 40)
    assert g["meets_both_floors"] is True
    assert g["reading"] == inc.GO
    assert g["short_by"] == []
    # the gate reports no direction, no rate and no effect of any kind
    assert not {k for k in g if "rate" in k or "direction" in k or "effect" in k}


def test_the_gate_counts_a_repeated_link_once() -> None:
    k = cell2.LocusKey("K562", "chr1", 0, 100, "G")
    assert inc.gate([k, k, k], "p", "call")["links"] == 1


def test_an_empty_population_is_a_no_go_and_not_an_error() -> None:
    g = inc.gate([], "p", "call")
    assert (g["links"], g["independent_loci"], g["meets_both_floors"]) == (0, 0, False)
    assert g["reading"] == inc.NO_GO


# ---- the readings were written before the outcome and admit no encouraging short count -----------


def test_there_are_exactly_two_readings() -> None:
    assert inc.READINGS == (inc.GO, inc.NO_GO)
    assert inc.registration()["readings_before_the_outcome_was_seen"]["there_is_no_third"] is True


@pytest.mark.parametrize(
    "word", ["close", "promising", "nearly enough", "a good start", "enough for a pilot"]
)
def test_the_no_go_forbids_each_encouraging_reading_by_name(word: str) -> None:
    assert word in inc.NO_GO
    assert "never" in inc.NO_GO


def test_the_go_reading_stops_the_lane_rather_than_starting_the_test() -> None:
    assert "stop" in inc.GO.lower()
    assert "registered as its own lane" in inc.GO
    for text in (inc.GO, inc.WHAT_FOLLOWS["go"]):
        assert "proposal" in text


def test_the_no_go_makes_the_inherited_finding_structural_and_carries_its_wording() -> None:
    assert inc.INHERITED_NO_GO == "the measured layer holds no repression call to test"
    assert inc.INHERITED_NO_GO in inc.WHAT_FOLLOWS["no_go"]
    assert "structural" in inc.WHAT_FOLLOWS["no_go"]


# ---- the data-or-code attribution is fixed before the count --------------------------------------


def test_every_ladder_step_after_the_first_has_a_registered_attribution() -> None:
    for step in inc.LADDER_STEPS[1:]:
        assert step in inc.ATTRIBUTION_OF_STEPS
        assert inc.attribution_of(step)["data_or_code"] in (inc.DATA, inc.CODE)


def test_only_the_significance_step_is_attributed_to_the_data() -> None:
    to_data = [s for s, (kind, _why) in inc.ATTRIBUTION_OF_STEPS.items() if kind == inc.DATA]
    assert to_data == ["and_significant"]


def test_the_extractor_itself_is_attributed_to_code_by_construction() -> None:
    assert inc.ATTRIBUTION_OF_STEPS["the_extractor"][0] == inc.CODE
    assert "not the finding of this count" in inc.ATTRIBUTION_OF_STEPS["the_extractor"][1]


# ---- the ladder, its mirror and the two populations ---------------------------------------------


def test_every_ladder_step_says_how_it_mirrors_lane_repress2s() -> None:
    assert set(inc.LADDER_MIRRORS) == set(inc.LADDER_STEPS)
    for step in inc.LADDER_STEPS[1:]:
        assert inc.LADDER_MIRRORS[step]


def test_the_inherited_ladder_is_carried_with_its_file_and_its_scope() -> None:
    assert inc.INHERITED_LADDER["file"] == "data/results/repress2_population.json"
    assert inc.INHERITED_LADDER["scope"] == "genome-wide, all 24 chromosomes"
    assert inc.INHERITED_LADDER["rules_of_this_axis"] == 156_925
    assert inc.INHERITED_LADDER["element_carries_a_crispri_pair"] == 357
    assert inc.INHERITED_LADDER["and_a_pair_on_the_rules_own_gene"] == 11
    assert inc.INHERITED_LADDER["and_in_the_rules_own_cell"] == 0
    assert inc.INHERITED_LADDER["and_that_pair_is_significant"] == 0
    assert inc.INHERITED_LADDER["gate_1_reading"] is rp.GATE_NO_GO


def test_the_inherited_measured_layer_count_is_carried_verbatim() -> None:
    m = inc.INHERITED_MEASURED_LAYER
    assert (m["assessable_measured_layer_rules"], m["on_a_pure_activates_target_axis"]) == (195, 143)
    assert m["on_a_pure_represses_target_axis"] == 0


def test_the_structural_cause_is_imported_and_not_restated() -> None:
    assert inc.STRUCTURAL_CAUSE is rp.MEASURED_AXIS_IS_NOT_THE_LINK_DIRECTION


def test_the_two_populations_are_never_pooled_and_p2_reconciles_with_their_step() -> None:
    assert inc.POPULATIONS == (inc.P1, inc.P2)
    assert "never pooled" in inc.NEVER_POOLED
    assert "at or below" in inc.P2_RECONCILES_WITH
    assert "and_in_the_rules_own_cell" in inc.P2_RECONCILES_WITH


def test_the_contest_rule_excludes_rather_than_relaxes() -> None:
    assert "is in neither population" in inc.P1_CONTEST_RULE
    assert "not relaxed if P1 falls short" in inc.P1_CONTEST_RULE


def test_the_overlap_rule_is_the_committed_one_and_is_quoted_with_its_value() -> None:
    assert str(ms.RECIPROCAL_OVERLAP) in inc.OVERLAP_CALL
    assert str(ms.REACH) in inc.OVERLAP_CALL


# ---- what the lane does not do, and what it cannot establish -------------------------------------


def test_the_registration_says_no_extractor_is_changed() -> None:
    assert "no extractor is changed by this lane" in inc.NO_EXTRACTOR_CHANGED
    assert "a written proposal and not an edit" in inc.NO_EXTRACTOR_CHANGED
    assert inc.registration()["no_extractor_changed"] is inc.NO_EXTRACTOR_CHANGED


def test_the_cannot_establish_section_refuses_both_directions_of_overreading() -> None:
    assert "never that a test of it would come out one way" in inc.CANNOT_ESTABLISH
    assert "never that the repression calls are right" in inc.CANNOT_ESTABLISH
    assert "never established biological independence" in inc.CANNOT_ESTABLISH


def test_the_registration_takes_no_count_and_makes_no_request() -> None:
    payload = inc.registration()
    assert payload["requests"] == 0
    assert payload["money"].startswith("none")
    # no count, share, rate or verdict about a population is present before the run
    assert not {k for k in payload if k in ("ladder_counts", "per_population", "verdict")}


def test_the_registration_is_a_json_safe_payload() -> None:
    import json

    assert json.loads(json.dumps(inc.registration()))["lane"] == "lane-increase"


# ---- the committed registration result ------------------------------------------------------------

REGISTRATION = Path("data/results/increase_registration.json")


def _registration() -> dict:
    return json.loads(REGISTRATION.read_text())


def test_the_registration_result_is_committed_and_carries_both_imported_floors() -> None:
    floors = _registration()["floors"]
    assert (floors["links"], floors["independent_loci"]) == (inc.POSITIVE_FLOOR, inc.LOCUS_FLOOR)
    assert (floors["links"], floors["independent_loci"]) == (30, 20)


def test_the_committed_registration_took_no_count() -> None:
    payload = _registration()
    assert payload["requests"] == 0
    assert "ladder" in payload and "steps" in payload["ladder"]
    # the ladder is registered as a list of step names, never as a step name mapped to a number
    assert payload["ladder"]["steps"] == list(inc.LADDER_STEPS)
    assert "verdict" not in payload


def test_the_committed_registration_names_what_the_run_will_write() -> None:
    assert _registration()["what_the_run_will_write"] == "data/results/increase_population.json"


def test_the_committed_registration_carries_lane_repress2s_wording_unchanged() -> None:
    follows = _registration()["follows_from"]
    assert follows["registered_no_go_carried_verbatim"] == inc.INHERITED_NO_GO
    assert follows["their_blind_ladder"]["and_in_the_rules_own_cell"] == 0
    assert follows["commits"] == ["7a5dd8b", "fd07258", "37a256c", "e9cc033"]


def test_the_committed_registrations_manifest_is_complete_with_its_cleanliness_block() -> None:
    m = _registration()["result_manifest"]
    assert m["complete"] is True
    assert not m.get("problems")
    assert set(mf.CLEANLINESS_KEYS) <= set(m["code_cleanliness"])
    assert m["code_cleanliness"]["own_code_is_committed"] is True
