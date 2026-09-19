# SPDX-License-Identifier: AGPL-3.0-or-later
"""Standardisation with a memory: the four mistakes of 2026-09-16 and 2026-09-17, as tests."""

from __future__ import annotations

from genomeos.compare import (
    Strata,
    difference,
    imbalance,
    input_presence,
    spearman,
    standardised,
    stratum_rates,
)

STRATA = Strata(gc=(0.35, 0.45, 0.55), tss=(1_000, 5_000, 20_000, 100_000))


def row(hit: bool, gc: float = 0.5, tss: int = 10_000) -> dict:
    return {"moved": hit, "gc": gc, "tss": tss}


def test_a_large_stratum_cannot_outvote_a_small_one() -> None:
    """The published defect: pooling every control once per target let one big stratum decide."""
    targets = [row(True, tss=500), row(True, tss=1_000_000)]
    controls = [row(False, tss=500) for _ in range(2)] + [row(True, tss=1_000_000) for _ in range(200)]

    out = standardised(targets, controls, STRATA)

    assert out["matched"]["b"] == 0.5  # mean of the two strata's rates
    assert out["raw"]["b"] > 0.98  # what pooling would have said
    assert out["targets_matched"] == 2


def test_the_control_arm_is_not_given_more_precision_than_it_has() -> None:
    """A standardised control rate has the n of the matched targets, not of the control rows.

    With the pooled n of 202 the p-value was 0.0001-scale; with the honest n of 2 it cannot be, and
    0.079 on two targets is what an arm this thin is allowed to say.
    """
    targets = [row(True, tss=500), row(True, tss=1_000_000)]
    controls = [row(False, tss=500) for _ in range(2)] + [row(True, tss=1_000_000) for _ in range(200)]

    out = standardised(targets, controls, STRATA)

    assert out["matched"]["p_one_sided"] > 0.05


def test_the_imbalance_is_in_the_result_and_not_left_to_the_reader() -> None:
    """The 882-block reading died of an arm imbalance nothing reported: a tenfold gap in TSS distance."""
    targets = [row(True, tss=80_000) for _ in range(10)]
    controls = [row(True, tss=420_000) for _ in range(10)]

    out = standardised(targets, controls, STRATA)

    assert out["imbalance"]["tss"]["target_median"] == 80_000
    assert out["imbalance"]["tss"]["control_median"] == 420_000
    assert out["imbalance"]["tss"]["ratio"] < 0.2


def test_targets_with_no_control_are_counted_not_hidden() -> None:
    """A comparison that keeps a tenth of its targets is not the comparison it looks like."""
    targets = [row(True, tss=500)] + [row(True, tss=50_000_000) for _ in range(9)]
    controls = [row(False, tss=500)]

    out = standardised(targets, controls, STRATA)

    assert out["targets_matched"] == 1
    assert out["dropped_for_want_of_a_control"] == 9


def test_stratum_rates_are_computed_once_over_the_pool() -> None:
    rates = stratum_rates([row(True), row(False), row(True)], STRATA)

    assert list(rates.values()) == [2 / 3]


def test_an_empty_arm_gives_no_difference_rather_than_a_crash() -> None:
    assert difference(0, 0, 3, 10)["difference"] is None
    out = standardised([], [row(True)], STRATA)
    assert out["raw"]["difference"] is None
    assert out["matched"]["difference"] is None


def test_equal_groups_read_as_no_difference() -> None:
    same = [row(True), row(False), row(True), row(False)]

    out = standardised(list(same), list(same), STRATA)

    assert out["raw"]["difference"] == 0.0
    assert out["matched"]["difference"] == 0.0
    assert out["imbalance"]["gc"]["ratio"] == 1.0


def test_strata_split_on_every_covariate_given() -> None:
    strata = Strata(gc=(0.5,), tss=(1_000,))

    assert strata.key({"gc": 0.4, "tss": 500}) == (0, 0)
    assert strata.key({"gc": 0.6, "tss": 5_000}) == (1, 1)
    assert strata.covariates == ("gc", "tss")


def test_imbalance_tolerates_a_missing_covariate() -> None:
    out = imbalance([{"gc": 0.5, "tss": None}], [{"gc": 0.5, "tss": 10}], Strata(gc=(0.5,), tss=(10,)))

    assert out["tss"]["target_median"] is None
    assert out["gc"]["ratio"] == 1.0


# ------------------------------------------------------- bought or free, before the stratum is run
def test_an_input_present_on_every_row_is_free_and_a_stratum_on_it_cannot_bite() -> None:
    """The 0.3529 -> 0.3578 result: no movement, because the input was never missing anywhere.

    A value on constrained sequence is read from a genome-wide phyloP track and a trio's variants, so
    every window has it. Reporting that as "conditioning on coverage changed nothing" would read as a
    claim passing a test when no test was administered.
    """
    targets = [{"phylop": 1.0} for _ in range(17)]
    controls = [{"phylop": 0.5} for _ in range(85)]

    out = input_presence(targets, controls, {"constraint": "phylop"})

    assert out["constraint"]["kind"] == "free"
    assert out["constraint"]["universal"]
    assert "arithmetic rather than evidence" in out["constraint"]["reading"]


def test_an_input_missing_on_some_rows_is_bought_and_names_the_shortfall() -> None:
    """60 of 85 control windows had a deletion spent on them; that gap is where the artefact lived."""
    targets = [{"scored": True} for _ in range(17)]
    controls = [{"scored": True} for _ in range(60)] + [{"scored": None} for _ in range(25)]

    out = input_presence(targets, controls, {"deletion": "scored"})

    assert out["deletion"]["kind"] == "bought"
    assert out["deletion"]["controls_with_the_input"] == 60
    assert "60/85" in out["deletion"]["reading"]


def test_the_arms_are_reported_apart_because_a_one_sided_gap_is_its_own_failure() -> None:
    """An input universal in the controls and missing in the targets is not the symmetric case.

    Pooling would read 95 of 102 and look like mild patchiness; apart, it reads 10 of 17 against
    85 of 85, which is a different and worse thing.
    """
    targets = [{"x": 1} for _ in range(10)] + [{"x": None} for _ in range(7)]
    controls = [{"x": 1} for _ in range(85)]

    out = input_presence(targets, controls, {"thing": "x"})["thing"]

    assert (out["targets_with_the_input"], out["controls_with_the_input"]) == (10, 85)
    assert out["kind"] == "bought"


def test_naming_the_outcome_instead_of_the_input_is_the_same_substitution_one_level_down() -> None:
    """`deletion_scored` is 60 of 85 and `deletion_target` is 54: the input, and what it produced.

    Passing the outcome measures how often the claim succeeded rather than how often it could be
    attempted. Both read "bought" here, so the function cannot catch it — the docstring has to, and
    this test pins the two numbers so the distinction stays visible.
    """
    controls_input = [{"k": True} for _ in range(60)] + [{"k": None} for _ in range(25)]
    controls_outcome = [{"k": True} for _ in range(54)] + [{"k": None} for _ in range(31)]
    targets = [{"k": True} for _ in range(17)]

    got_input = input_presence(targets, controls_input, {"d": "k"})["d"]
    got_outcome = input_presence(targets, controls_outcome, {"d": "k"})["d"]

    assert got_input["controls_with_the_input"] == 60
    assert got_outcome["controls_with_the_input"] == 54
    assert got_input["kind"] == got_outcome["kind"] == "bought"


def test_the_upper_bound_is_one_sided_because_the_p_beside_it_is() -> None:
    """A one-sided p with a two-sided z is two confidence levels under one name.

    Until 2026-09-18 `upper_95` used 1.96, the two-sided 95% constant, which makes it a one-sided
    97.5% bound. `attribution/executor.py`'s own `difference_upper` had always used 1.645 and had
    always been right, so the shared helper was the odd one out while nine lanes imported it. No
    verdict moved: no gate read this field, executor gated on its own, and the only other consumer
    was cancelled before it ran.
    """
    import math

    out = difference(60, 100, 40, 100)
    pa, pb = 0.6, 0.4
    se = math.sqrt(pa * (1 - pa) / 100 + pb * (1 - pb) / 100)

    assert out["upper_95"] == round(pa - pb + 1.645 * se, 4)
    assert out["upper_95"] < round(pa - pb + 1.96 * se, 4), "one-sided is tighter than two-sided"


def test_the_shared_helper_and_the_executor_lane_now_agree() -> None:
    """Two implementations of one bound must not disagree on which 95% they mean."""
    import importlib

    ex = importlib.import_module("genomeos.attribution.executor")

    mine = difference(60, 100, 40, 100)["upper_95"]
    theirs = ex.difference_upper(60, 100, 40, 100)

    assert theirs is not None
    assert abs(mine - theirs) < 1e-4


# -------------------------------------------------------- one Spearman, with the bar at the call site
def test_the_shared_spearman_matches_the_de_facto_canonical_one() -> None:
    """Eight copies existed and the algorithm agreed in all of them; only the bar diverged.

    `attribution/closure.py` held the version other lanes delegated to. This asserts the shared one
    reproduces it exactly over random inputs, including the None cases, so the two cannot drift while
    both exist.
    """
    import random

    from genomeos.attribution.closure import spearman as canonical

    rng = random.Random(7)
    for _ in range(2000):
        n = rng.randint(1, 12)
        xs = [round(rng.uniform(-3, 3), 3) for _ in range(n)]
        ys = [round(rng.uniform(-3, 3), 3) for _ in range(n)]
        mine, theirs = spearman(xs, ys), canonical(xs, ys)
        assert (mine is None) == (theirs is None)
        if mine is not None:
            assert abs(mine - theirs) < 1e-12


def test_the_minimum_n_is_an_argument_so_the_bar_is_visible_where_it_is_chosen() -> None:
    """The divergence was 3 / 4 / 10 hidden in eight module constants, deciding admissibility."""
    three = ([1, 2, 3], [3, 2, 1])

    assert spearman(*three) == -1.0
    assert spearman(*three, min_n=4) is None
    assert spearman([1, 2, 3, 4], [1, 2, 3, 4], min_n=4) == 1.0
    # 3 is a floor, not a default worth overriding downwards: below it there is no order to correlate
    assert spearman([1, 2], [2, 1], min_n=1) is None


def test_a_constant_side_is_none_rather_than_zero() -> None:
    """ "No ranking" and "no relationship" are different statements and must not share an output."""
    assert spearman([1, 1, 1, 1], [1, 2, 3, 4]) is None
    assert spearman([1, 2, 3, 4], [2, 2, 2, 2]) is None
