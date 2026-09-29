# SPDX-License-Identifier: AGPL-3.0-or-later
"""Stage 4's gate against its registration, and the ways it is allowed to fail.

The instrument is only worth its run if it can return "failed". These build the eight runs the gate
reads, break one thing at a time, and check the verdict notices — including the one break that would
look exactly like a finding: a partitioning that separated the two proteins.
"""

from __future__ import annotations

import copy
import importlib.util
import math
from pathlib import Path

import pytest

from genomeos.runtime.division_gate import (
    LOW_COPIES,
    PREDICTED_SEPARATION_DECADES,
    PREDICTED_STEADY_RATIO_DIVIDING,
    PREDICTED_STEADY_RATIO_STATIC,
    SEEDS,
    START_AMOUNT,
)

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("division_gate", ROOT / "scripts" / "division_gate.py")
assert spec and spec.loader
dg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dg)

STABLE, SHORT, UNDEG = dg.STABLE, dg.SHORT, dg.UNDEGRADED


def chase(stable: float, short: float, undeg: float, fired: int = 8) -> dict:
    final = {STABLE: stable, SHORT: short, UNDEG: undeg}
    return {
        "final": final,
        "closed_form": dict(final),
        "matches_closed_form": True,
        "divisions_fired": fired,
        "conserved": True,
        "separation_decades": math.log10(stable) - math.log10(short),
    }


def passing() -> dict:
    """The runs the real instrument produced, as data, so a break is one edited field."""
    part = chase(2243.5514746035838, 9.094947017729282e-07, START_AMOUNT * 2.0**-8)
    dup = chase(574349.1774985174, 0.00023283064365386963, START_AMOUNT)
    return {
        "partitioned": part,
        "duplicating": dup,
        "equal_half_lives": chase(2243.5, 2243.5, START_AMOUNT * 2.0**-8),
        "checkpoint_unmet": chase(574349.1774985174, 0.00023283064365386963, START_AMOUNT, fired=0),
        "binomial_low": {
            "mean": 9.8895,
            "variance": 4.8993,
            "analytic_mean": LOW_COPIES / 2,
            "analytic_variance": LOW_COPIES / 4,
            "conserved": SEEDS,
            "whole_daughters": SEEDS,
        },
        "binomial_high": {"variance": 0.0, "conserved": SEEDS},
        "steady_state": {
            "ratio_static": PREDICTED_STEADY_RATIO_STATIC,
            "ratio_dividing": PREDICTED_STEADY_RATIO_DIVIDING,
        },
        "resource_checkpoint": {"finding": "reported, not scored"},
    }


def test_the_measured_run_passes_all_eight_clauses() -> None:
    got = dg.verdict(passing())

    assert got["passed"] and all(got["clauses"].values())
    assert got["separation_partitioned_decades"] == got["separation_duplicating_decades"]
    assert got["dilution_contributes_decades"] == 0.0
    assert got["separation_partitioned_decades"] == pytest.approx(PREDICTED_SEPARATION_DECADES, rel=1e-12)


def test_a_partitioning_that_separated_the_two_would_fail() -> None:
    """The break that looks like the finding: dilution is gene-blind, so this cannot happen."""
    runs = copy.deepcopy(passing())
    runs["partitioned"]["separation_decades"] *= 1.5

    got = dg.verdict(runs)

    assert not got["passed"] and not got["clauses"]["3_dilution_separates_nothing"]


def test_a_division_that_created_matter_would_fail() -> None:
    runs = copy.deepcopy(passing())
    runs["binomial_low"]["conserved"] = SEEDS - 1

    got = dg.verdict(runs)

    assert not got["passed"] and not got["clauses"]["5_conservation"]


def test_degrading_a_protein_with_no_declared_half_life_would_fail() -> None:
    runs = copy.deepcopy(passing())
    runs["partitioned"]["final"][UNDEG] *= 0.9

    got = dg.verdict(runs)

    assert not got["passed"] and not got["clauses"]["2_dilution_only"]


def test_dividing_with_the_checkpoint_unmet_would_fail() -> None:
    runs = copy.deepcopy(passing())
    runs["checkpoint_unmet"]["divisions_fired"] = 3

    got = dg.verdict(runs)

    assert not got["passed"] and not got["clauses"]["8_checkpoint"]


def test_two_equal_half_lives_that_differed_would_fail() -> None:
    runs = copy.deepcopy(passing())
    runs["equal_half_lives"]["final"][SHORT] *= 1 + 1e-12

    got = dg.verdict(runs)

    assert not got["passed"] and not got["clauses"]["7_equal_half_lives"]


def test_a_binomial_split_of_the_wrong_width_would_fail() -> None:
    runs = copy.deepcopy(passing())
    runs["binomial_low"]["variance"] = LOW_COPIES / 8.0

    got = dg.verdict(runs)

    assert not got["passed"] and not got["clauses"]["6_binomial"]


def test_the_resource_checkpoint_cannot_contribute_to_a_pass() -> None:
    """It was not registered, so it is reported and excluded — as gate (a)'s magnitude clause was."""
    runs = copy.deepcopy(passing())
    runs["resource_checkpoint"] = {"finding": "anything at all", "discriminates": False}

    got = dg.verdict(runs)

    assert got["passed"]
    assert "reported" not in got["clauses"]
    assert got["resource_checkpoint_excluded"] == "anything at all"


def test_the_verdict_never_claims_the_gate_as_worded_passed() -> None:
    got = dg.verdict(passing())

    assert got["gate_as_worded"].startswith("REFUTED")
    assert "not askable" in got["measured_turnover"]
    assert "cannot" in got["verdict"]
