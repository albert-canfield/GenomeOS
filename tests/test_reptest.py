# SPDX-License-Identifier: AGPL-3.0-or-later
"""The answerability count's registration, pinned: the floors are imported, the causes are exhaustive
and disjoint, a short count is a no-go with its margin, the collapse identity's premise is checked on
this record shape rather than assumed, and R2 is refused in code and not in prose."""

from __future__ import annotations

import pytest

from genomeos.attribution import cell2, fresh
from genomeos.attribution import direction_link as dl
from genomeos.attribution import increase_links as il
from genomeos.attribution import measured as ms
from genomeos.attribution import reptest as rt

# ---- the floors are imported, never restated and never chosen here -------------------------------


def test_floors_are_the_imported_objects_and_not_copies() -> None:
    assert rt.POSITIVE_FLOOR is fresh.POSITIVE_FLOOR
    assert rt.LOCUS_FLOOR is fresh.LOCUS_FLOOR
    assert fresh.LOCUS_FLOOR is cell2.POOLED_LOCUS_FLOOR
    assert (rt.POSITIVE_FLOOR, rt.LOCUS_FLOOR) == (30, 20)


def test_the_floors_name_where_they_came_from() -> None:
    f = rt.FLOORS
    assert f["answerable_links_at_least"] == 30
    assert f["independent_loci_at_least"] == 20
    assert f["answerable_links_floor_imported_from"] == "genomeos.attribution.fresh.POSITIVE_FLOOR"
    assert "cell2.POOLED_LOCUS_FLOOR" in f["independent_loci_floor_imported_from"]
    assert "ANSWERABLE" in f["applied_to"]


def test_the_locus_rule_and_its_disclaimer_are_cell2s_own_words() -> None:
    assert rt.LOCUS_RULE is cell2.INDEPENDENT_LOCUS_RULE
    assert rt.LOCUS_SPAN == cell2.INDEPENDENT_LOCUS_SPAN == 1_000_000
    assert "not established biological independence" in rt.LOCUS_RULE
    assert "not established biological independence" in rt.NOT_BIOLOGICAL_INDEPENDENCE


# ---- the population is the extractor's output, not a join rebuilt beside it ----------------------


def test_the_population_is_read_out_of_the_extractor_under_v2() -> None:
    assert rt.EXTRACTOR is ms.EXTRACTOR_V2 == "v2"
    call = rt.POPULATION_IS_THE_EXTRACTORS_OUTPUT
    assert "rule_links_detail" in call
    assert "measured.rows" in call
    assert "is not a join rebuilt beside it" in call


def test_the_by_construction_property_is_carried_not_restated() -> None:
    assert rt.POPULATION_IS_NOT_COVERAGE is il.BY_CONSTRUCTION
    assert rt.DECREASE_ARM_IS_NOT_MODEL_COVERAGE is dl.DECREASE_ARM_IS_NOT_MODEL_COVERAGE
    assert "NOT evidence that the model covered" in rt.POPULATION_IS_NOT_COVERAGE
    assert "is a coverage figure" in rt.COUNT_IS_NOT_COVERAGE
    assert rt.COUNT_IS_NOT_COVERAGE.startswith("neither the 48 nor the answerable subset")


def test_the_predicate_is_lane_directions_own_and_not_a_second_wording() -> None:
    assert rt.PREDICTED_CALL is dl.PREDICTED_CALL
    assert rt.SIGN_CALL is dl.SIGN_CALL
    assert rt.CACHED_CELLS is dl.CACHED_CELLS
    assert rt.CELLS_NOT_CACHED is dl.CELLS_NOT_CACHED
    assert "elements_hct116" in rt.ANSWERABLE_CALL and "NOT read" in rt.ANSWERABLE_CALL


# ---- the causes: exhaustive, disjoint, and one per link ------------------------------------------


def _record(genes: list[dict[str, object]]) -> dict[str, object]:
    return {"id": "e1", "genes": genes}


def test_a_cell_the_cache_carries_no_track_of_is_its_own_cause() -> None:
    for cell in rt.CELLS_NOT_CACHED:
        cause, value = rt.cause_of(_record([{"gene": "G", "by_cell": {"K562": -1.0}}]), "G", cell)
        assert cause == rt.CELL_NOT_CACHED
        assert value is None


def test_an_element_the_archive_does_not_hold_is_its_own_cause() -> None:
    assert rt.cause_of(None, "G", "K562") == (rt.ELEMENT_ABSENT, None)


def test_a_gene_outside_the_window_is_counted_apart_from_a_gene_off_the_track() -> None:
    not_in_window = rt.cause_of(_record([{"gene": "OTHER", "by_cell": {"K562": -1.0}}]), "G", "K562")
    not_on_track = rt.cause_of(_record([{"gene": "G", "by_cell": {"HepG2": -1.0}}]), "G", "K562")
    assert not_in_window == (rt.GENE_NOT_IN_WINDOW, None)
    assert not_on_track == (rt.GENE_NOT_ON_TRACK, None)
    assert rt.GENE_NOT_IN_WINDOW != rt.GENE_NOT_ON_TRACK


def test_an_exactly_zero_value_carries_no_direction_and_is_not_answerable() -> None:
    cause, value = rt.cause_of(_record([{"gene": "G", "by_cell": {"K562": 0.0}}]), "G", "K562")
    assert cause == rt.VALUE_EXACTLY_ZERO
    assert value == 0.0
    assert dl.sign_of(0.0) is None


@pytest.mark.parametrize("v", [-2.5, -1e-9, 1e-9, 3.0])
def test_a_signed_value_on_the_links_own_cell_track_is_answerable(v: float) -> None:
    cause, value = rt.cause_of(_record([{"gene": "G", "by_cell": {"K562": v}}]), "G", "K562")
    assert cause == rt.ANSWERABLE
    assert value == v


def test_every_cause_cause_of_can_return_is_in_CAUSES_and_each_has_its_own_words() -> None:
    returned = {
        rt.cause_of(*args)[0]
        for args in (
            (_record([{"gene": "G", "by_cell": {"K562": -1.0}}]), "G", "WTC11"),
            (None, "G", "K562"),
            (_record([{"gene": "OTHER", "by_cell": {"K562": -1.0}}]), "G", "K562"),
            (_record([{"gene": "G", "by_cell": {"HepG2": -1.0}}]), "G", "K562"),
            (_record([{"gene": "G", "by_cell": {"K562": 0.0}}]), "G", "K562"),
            (_record([{"gene": "G", "by_cell": {"K562": -1.0}}]), "G", "K562"),
        )
    }
    assert returned == set(rt.CAUSES)
    assert set(rt.CAUSE_TEXT) == set(rt.CAUSES)
    assert len(set(rt.CAUSE_TEXT.values())) == len(rt.CAUSES)


def test_the_causes_must_sum_to_the_population_or_the_check_refuses() -> None:
    assert rt.causes_reconcile({rt.ANSWERABLE: 21, rt.CELL_NOT_CACHED: 9, rt.GENE_NOT_IN_WINDOW: 18}, 48)
    assert not rt.causes_reconcile({rt.ANSWERABLE: 21, rt.CELL_NOT_CACHED: 9}, 48)
    assert not rt.causes_reconcile({"invented_cause": 48}, 48)


# ---- the verdict: both floors or a no-go, with the margin ----------------------------------------


def _links(n: int, loci: int, cell: str = "K562") -> list[dict[str, object]]:
    """`n` links spread over `loci` distinct cell2 loci: distinct genes, 10 Mb apart, so no chaining."""
    out = []
    for i in range(n):
        g = i % loci
        out.append(
            {
                "cell": cell,
                "chrom": "chr1",
                "start": g * 10_000_000,
                "end": g * 10_000_000 + 500,
                "gene": f"G{g}",
            }
        )
    return out


def test_the_answerable_count_lane_direction_measured_is_a_no_go_with_both_margins() -> None:
    v = rt.verdict(_links(21, 13), population=48)
    assert v["answerable_links"] == 21
    assert v["population_measured_links"] == 48
    assert v["independent_loci_of_the_answerable_links"] == 13
    assert v["meets_both_floors"] is False
    assert v["short_by"] == [
        "answerable links 21, short of the imported floor 30 by 9",
        "independent loci 13, short of the imported floor 20 by 7",
    ]
    assert v["reading"] is rt.GATE_NO_GO
    assert v["what_follows"] == "nothing further is registered and no direction is read"


def test_clearing_one_floor_and_missing_the_other_is_still_a_no_go() -> None:
    many_links_few_loci = rt.verdict(_links(40, 11), population=48)
    assert many_links_few_loci["meets_both_floors"] is False
    assert len(many_links_few_loci["short_by"]) == 1
    assert "independent loci 11, short of the imported floor 20 by 9" in many_links_few_loci["short_by"]


def test_both_floors_clear_only_at_or_above_both_and_then_a_test_may_be_registered() -> None:
    assert rt.verdict(_links(29, 29), population=48)["meets_both_floors"] is False
    go = rt.verdict(_links(30, 20), population=48)
    assert go["meets_both_floors"] is True
    assert go["reading"] is rt.GATE_GO
    assert go["what_follows"].startswith("register a direction test")


def test_the_no_go_wording_forbids_the_softening_words_and_names_them() -> None:
    assert rt.FORBIDDEN_OF_A_SHORT_COUNT is dl.FORBIDDEN_OF_A_SHORT_OR_UNDETECTED_OUTCOME
    assert rt.FORBIDDEN_OF_A_SHORT_COUNT == (
        "close",
        "promising",
        "nearly enough",
        "a good start",
        "enough for a pilot",
    )
    for word in rt.FORBIDDEN_OF_A_SHORT_COUNT:
        assert word in rt.GATE_NO_GO  # named as forbidden, in the registered sentence that forbids them
    assert "NEVER as close" in rt.GATE_NO_GO
    assert rt.READINGS == (rt.GATE_GO, rt.GATE_NO_GO) and rt.THERE_IS_NO_THIRD is True


# ---- the statistic: the collapse identity's premise, checked on this record shape ----------------


def test_the_collapse_identity_is_carried_verbatim_and_not_rewritten() -> None:
    assert rt.COLLAPSE is dl.COLLAPSE
    assert rt.TEST_STATISTIC is dl.TEST_STATISTIC
    assert rt.CHANCE == dl.CHANCE == 0.5
    assert "2 * (balanced_accuracy - 0.5) * (n_u - n_d) / (n_d + n_u)" in rt.COLLAPSE


@pytest.mark.parametrize(
    ("n_down", "n_up", "expect"),
    [
        (174, 21, "decreases arm is the larger one"),
        (5, 40, "increases arm is the larger one"),
        (30, 30, "IDENTICALLY ZERO"),
    ],
)
def test_collapse_applies_reports_the_branch_this_populations_own_counts_put_it_in(
    n_down: int, n_up: int, expect: str
) -> None:
    got = rt.collapse_applies(n_down, n_up)
    assert got["answered_decrease_links"] == n_down
    assert got["answered_increase_links"] == n_up
    assert got["arms_are_equal"] is (n_down == n_up)
    assert expect in got["implication"]
    assert got["premise_recomputed_on_this_population"] is True
    assert got["difference_is_identifiable_as_a_skill_statistic"] is False
    assert got["statistic_a_test_would_be_registered_on"] is dl.TEST_STATISTIC


def test_with_no_answered_link_in_either_arm_nothing_is_identifiable() -> None:
    got = rt.collapse_applies(0, 0)
    assert got["arm_size_factor_n_u_minus_n_d_over_n"] is None
    assert "no statistic is identifiable at all" in got["implication"]


@pytest.mark.parametrize(
    ("n_down", "a", "n_up", "b"),
    [(174, 157, 21, 9), (174, 100, 21, 21), (30, 30, 30, 0), (30, 15, 30, 15), (12, 3, 47, 40)],
)
def test_the_identity_holds_on_populations_of_this_lanes_own_denominators(
    n_down: int, a: int, n_up: int, b: int
) -> None:
    """The algebra is not assumed to transfer: it is recomputed through lane-direction's own functions
    on link records of the shape this lane's arms have, at this lane's own denominators."""
    answered = [
        {"arm": dl.DECREASES, "measured_sign": -1, "predicted_sign": -1 if i < a else 1, "agrees": i < a}
        for i in range(n_down)
    ] + [
        {"arm": dl.INCREASES, "measured_sign": 1, "predicted_sign": 1 if i < b else -1, "agrees": i < b}
        for i in range(n_up)
    ]
    by_arm = {arm: [r for r in answered if r["arm"] == arm] for arm in dl.ARMS}
    ba = dl.balanced_accuracy(answered)
    got = dl.excess(answered)
    predicted = 2 * (ba - 0.5) * (n_up - n_down) / (n_down + n_up)
    assert got == pytest.approx(predicted, abs=1e-12)
    assert dl.difference(by_arm) == pytest.approx(a / n_down - b / n_up, abs=1e-12)
    if n_down == n_up:
        assert got == pytest.approx(0.0, abs=1e-12)
    collapse = rt.collapse_applies(n_down, n_up)
    assert collapse["arm_size_factor_n_u_minus_n_d_over_n"] == pytest.approx(
        (n_up - n_down) / (n_down + n_up), abs=1e-6
    )


# ---- R2, refused in code -------------------------------------------------------------------------


def test_r2_is_enforced_by_the_committed_function_and_not_reimplemented_here() -> None:
    assert rt.registration()["binding_wording_r2"] is il.R2_WORDING
    assert not hasattr(rt, "FORBIDDEN_OF_AN_INCREASE_LINK")  # not copied, not weakened, not shadowed
    assert not hasattr(rt, "check_no_mechanism_claim")


@pytest.mark.parametrize(
    "planted",
    [
        {"molecular_role": "silencer"},
        {"molecular_role": "repressor"},
        {"note": "this element represses the gene"},
        {"note": "evidence of repression"},
    ],
)
def test_a_planted_mechanism_claim_on_an_increase_derived_link_is_refused(planted: dict) -> None:
    record = {
        "gene": "G",
        "cell": "K562",
        "derived_from": ms.INCREASE,
        "outcome": il.OUTCOME_TEXT,
        "molecular_role": il.MOLECULAR_ROLE,
        **planted,
    }
    with pytest.raises(ValueError):
        il.check_no_mechanism_claim([record])


def test_an_increase_derived_link_may_say_the_gene_went_up_on_knockdown() -> None:
    il.check_no_mechanism_claim(
        [
            {
                "gene": "G",
                "cell": "K562",
                "derived_from": ms.INCREASE,
                "outcome": il.OUTCOME_TEXT,
                "molecular_role": il.MOLECULAR_ROLE,
            }
        ]
    )
    assert il.OUTCOME_TEXT == "increase on knockdown"
    assert il.MOLECULAR_ROLE is None


# ---- the registration is complete and states what it does not do ---------------------------------


def test_the_registration_carries_the_inherited_no_go_word_for_word() -> None:
    r = rt.registration()
    assert r["acts_on"]["lane_direction_no_go_carried_verbatim"] is dl.GATE_NO_GO
    assert r["acts_on"]["figures"]["lane_direction_increase_arm_measured"] == 48
    assert r["acts_on"]["figures"]["lane_direction_increase_arm_answered"] == 21
    assert r["acts_on"]["figures"]["v2_inhibits_links_added"] == 48
    assert r["acts_on"]["figures"]["conflicts_today"] == 0


def test_the_registration_names_every_piece_the_count_needs_and_registers_no_test() -> None:
    r = rt.registration()
    for key in (
        "population",
        "answerable_predicate",
        "sign_rule_zero_and_absent",
        "causes",
        "causes_are_exhaustive",
        "floors",
        "gate_call",
        "readings",
        "statistic",
        "statistic_is_checked_not_assumed",
        "binding_wording_r2",
        "by_construction",
        "count_is_not_coverage",
        "conflict_is_not_a_finding",
        "independent_locus_rule",
        "not_biological_independence",
        "budget",
        "order",
    ):
        assert r[key], key
    assert "No direction test is registered" in r["what_is_registered"]
    assert "0 model requests" in r["budget"]
    assert set(r["causes"]) == set(rt.CAUSES)


def test_a_conflict_count_of_zero_is_not_a_finding_that_the_measurements_agree() -> None:
    text = rt.registration()["conflict_is_not_a_finding"]
    assert text is il.CONFLICT_IS_NOT_A_FINDING
    assert "not a finding that the measurements agree" in text
