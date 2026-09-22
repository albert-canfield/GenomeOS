# SPDX-License-Identifier: AGPL-3.0-or-later
"""Stage 2's gate (a) against its registration, and the two ways it is allowed to fail."""

from __future__ import annotations

import importlib.util
from pathlib import Path

from genomeos.runtime.burden_gate import LEVELS, PRE_REGISTRATION

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("burden_gate", ROOT / "scripts" / "burden_gate.py")
assert spec and spec.loader
bg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bg)


def rows(factors: list[float]) -> list[dict]:
    return [{"level": lvl, "unrelated_factor": f} for lvl, f in zip(LEVELS, factors, strict=False)]


FREE = rows([1.0] * len(LEVELS))


def test_a_falling_monotonic_curve_with_a_clean_falsifier_passes() -> None:
    got = bg.verdict(rows([1.0, 1.0, 1.0, 1.0, 0.5, 0.25, 0.05]), FREE)

    assert got["passed"] and got["direction"] and got["monotonic"] and got["falsifier_holds"]


def test_no_reduction_fails_the_gate() -> None:
    """§5.3's first failure: the burden simply does not appear."""
    got = bg.verdict(rows([1.0] * len(LEVELS)), FREE)

    assert not got["passed"]
    assert "direction" in got["verdict"]


def test_a_reduction_with_the_pool_switched_off_fails_the_gate() -> None:
    """§5.3's second failure, and the one a bug passes.

    A model that cuts an unrelated gene whatever the capacity is not reproducing a burden. Only the
    unconstrained arm can tell that apart from the real thing, which is why the registration names it.
    """
    still_cut = rows([1.0, 1.0, 1.0, 1.0, 0.5, 0.25, 0.05])

    got = bg.verdict(still_cut, still_cut)

    assert not got["passed"]
    assert "falsifier" in got["verdict"]
    assert "§5.3 calls falsified" in got["falsifier_reads"]


def test_a_curve_that_recovers_under_more_demand_fails_monotonicity() -> None:
    """A resource that is competed for does not become more available under more demand."""
    got = bg.verdict(rows([1.0, 1.0, 1.0, 0.5, 0.9, 0.25, 0.05]), FREE)

    assert not got["passed"] and not got["monotonic"]


def test_the_magnitude_clause_cannot_contribute_to_a_pass() -> None:
    """An unasked question is not evidence; the registration says so and the verdict must honour it."""
    got = bg.verdict(rows([1.0, 1.0, 1.0, 1.0, 0.5, 0.25, 0.05]), FREE)

    assert "not askable" in got["magnitude"]
    assert "NOT ASKABLE" in PRE_REGISTRATION["outcomes"]["magnitude"]
    assert "unasked question is not evidence" in PRE_REGISTRATION["verdict_rule"]
