# SPDX-License-Identifier: AGPL-3.0-or-later
"""The activator combination rules (lane-combine, 2026-09-28): what each can and cannot produce."""

from __future__ import annotations

import itertools

import pytest

from genomeos.runtime import grn

TERMS = [0.1, 0.35, 0.8]


def test_the_runtime_still_combines_by_the_mean():
    assert grn.ACTIVATOR_COMBINATION == "mean"
    assert grn.combine_activators(TERMS) == pytest.approx(sum(TERMS) / 3)


def test_every_rule_agrees_on_one_activator_and_on_none():
    for rule in grn.ACTIVATOR_RULES:
        assert grn.combine_activators([0.37], rule) == pytest.approx(0.37)
        assert grn.combine_activators([], rule) == 1.0  # a gene with no activator runs at full rate


def test_the_four_rules():
    assert grn.combine_activators(TERMS, "sum_capped") == 1.0
    assert grn.combine_activators([0.1, 0.2], "sum_capped") == pytest.approx(0.3)
    assert grn.combine_activators(TERMS, "max") == 0.8
    assert grn.combine_activators([0.5, 0.5], "or") == pytest.approx(0.75)
    with pytest.raises(ValueError):
        grn.combine_activators(TERMS, "median")


def test_under_the_mean_single_removal_changes_sum_to_zero():
    """The parameter-free identity the registration states: sum_i (A_-i - A) = 0 under the mean."""
    a = grn.combine_activators(TERMS, "mean")
    changes = [grn.combine_activators(TERMS[:i] + TERMS[i + 1 :], "mean") - a for i in range(3)]
    assert sum(changes) == pytest.approx(0.0, abs=1e-12)
    assert any(c > 0 for c in changes)  # removing a below-average activator raises the drive


@pytest.mark.parametrize("rule", ["sum_capped", "max", "or"])
def test_the_other_rules_never_fall_when_an_activator_is_added(rule):
    for k in range(1, 4):
        for subset in itertools.combinations(TERMS, k):
            for extra in (0.0, 0.05, 0.5, 1.0):
                assert (
                    grn.combine_activators([*subset, extra], rule)
                    >= grn.combine_activators(list(subset), rule) - 1e-12
                )


def test_max_cannot_give_two_single_removal_decreases():
    a = grn.combine_activators(TERMS, "max")
    drops = [grn.combine_activators(TERMS[:i] + TERMS[i + 1 :], "max") < a for i in range(3)]
    assert sum(drops) == 1


def test_the_choice_test_is_registered_before_its_run():
    t = grn.COMBINATION_TEST
    assert t["pass_rule"].startswith("switch to the single surviving rule")
    assert t["exclude_cases_contradicting_every_rule"] is True
    assert set(t["cases"]["bothma2015_kni_early"]["contradicts"]) == set(grn.ACTIVATOR_RULES)
    assert "gastrulation_census" in t["not_judged"]
