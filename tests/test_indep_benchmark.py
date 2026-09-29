# SPDX-License-Identifier: AGPL-3.0-or-later
"""The independent benchmark's metric and registered reading (genomeos/attribution/indep.py)."""

from genomeos.attribution import indep


def test_average_precision_perfect_and_tied():
    assert indep.average_precision([3, 2, 1], [True, False, False]) == 1.0
    # all tied: one threshold, precision = base rate
    assert indep.average_precision([1, 1, 1, 1], [True, False, False, False]) == 0.25
    assert indep.average_precision([1, 2], [False, False]) is None


def test_average_precision_matches_the_step_definition():
    # ranks: T F T F -> recall 0.5 at P=1, recall 1 at P=2/3
    ap = indep.average_precision([4, 3, 2, 1], [True, False, True, False])
    assert abs(ap - (0.5 * 1 + 0.5 * 2 / 3)) < 1e-12


def test_outcome_keeps_increase_apart():
    assert indep.outcome({"fdr": 0.01, "effect": -0.3}) == "decrease"
    assert indep.outcome({"fdr": 0.01, "effect": 0.3}) == "increase"
    assert indep.outcome({"fdr": 0.2, "effect": -0.9}) == "not_detected"


def test_near_tss():
    tss = [1_000, 50_000]
    assert indep.near_tss(tss, 1_500, 1_800)
    assert not indep.near_tss(tss, 3_000, 3_500)
    assert indep.near_tss(tss, 48_500, 49_100)


def test_verdict_readings():
    assert indep.verdict([0.01, 0.1], [0.2, 0.4], 0.05, True) == "PASS"
    assert indep.verdict([-0.1, -0.01], [0.2, 0.4], 0.05, True) == "FAIL"
    assert indep.verdict([-0.05, 0.1], [0.01, 0.04], 0.05, True) == "FAIL"
    assert indep.verdict([-0.05, 0.1], [0.1, 0.3], 0.05, True) == "NOT ESTABLISHED"
    assert indep.verdict([0.01, 0.1], [0.2, 0.4], 0.05, False) == "NOT READABLE"


def test_element_bootstrap_moves_clusters_together():
    pairs = [
        {
            "element": f"e{i // 2}",
            "gene": f"g{i % 3}",
            "label": "decrease" if i % 4 == 0 else "not_detected",
            "a": float(i % 4 == 0),
            "b": float(i),
        }
        for i in range(40)
    ]
    draws = indep.bootstrap(pairs, ["a", "b"], n=50)
    assert len(draws["a"]) == len(draws["a-b"]) > 0
    assert all(x == 1.0 for x in draws["a"])


def test_registration_names_its_rules():
    text = indep.PREREGISTERED
    for phrase in (
        "PASS.",
        "FALSIFIER.",
        "NOT ESTABLISHED",
        "COVERAGE.",
        "never imputed",
        "500 kb",
        "ELEMENTS",
    ):
        assert phrase in text
