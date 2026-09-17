# SPDX-License-Identifier: AGPL-3.0-or-later
"""The VISTA in-silico mutagenesis instrument: its windows, its rank correlation, its bars."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("satmut_vista", ROOT / "scripts" / "satmut_vista.py")
assert spec and spec.loader
sv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sv)

from genomeos.attribution import satmut_vista as reg  # noqa: E402


def test_the_windows_are_evenly_spaced_and_stay_inside_the_element() -> None:
    """Evenly spaced, not the first twenty: a long element sampled at one end is not the element."""
    got = sv.windows_of({"start": 1_000, "end": 3_000})

    assert len(got) == reg.WINDOWS_PER_ELEMENT
    assert got[0][0] == 1_000
    assert all(e - s == reg.WINDOW for s, e in got)
    assert got[-1][1] <= 3_000
    assert got[0][1] <= got[1][0] or got[1][0] >= got[0][0]  # ordered along the element


def test_an_element_shorter_than_two_windows_yields_none() -> None:
    assert sv.windows_of({"start": 0, "end": reg.WINDOW}) == []


def test_spearman_is_a_rank_correlation_with_ties_averaged() -> None:
    assert sv.spearman([1, 2, 3, 4], [1, 2, 3, 4]) == 1.0
    assert sv.spearman([1, 2, 3, 4], [4, 3, 2, 1]) == -1.0
    assert sv.spearman([1, 2, 3, 4], [1, 1, 2, 2]) is not None  # ties are allowed


def test_a_constant_side_has_no_order_to_correlate() -> None:
    """None rather than 0: 'no ranking' and 'no relationship' are different statements."""
    assert sv.spearman([1, 1, 1, 1], [1, 2, 3, 4]) is None
    assert sv.spearman([1, 2], [2, 1]) is None  # too few windows to rank


def test_the_verdict_applies_the_registered_bars_without_interpretation() -> None:
    assert sv.verdict(0.15, 0, 200).startswith("success")
    assert sv.verdict(0.15, 100, 200).startswith("weak")  # size without the p is weak, as registered
    assert sv.verdict(0.07, 0, 200).startswith("weak")
    assert sv.verdict(0.01, 0, 200).startswith("failure")
    assert sv.verdict(-0.2, 0, 200).startswith("failure")
    assert "not readable" in sv.verdict(None, 0, 200)


def test_the_reading_needs_both_a_scored_effect_and_a_constraint() -> None:
    """A window with one of the two carries no correlation and must not be counted as zero."""
    rows = [
        {"element": "e1", "arm": "positive", "effect": 1.0, "phylop_mean": 2.0},
        {"element": "e1", "arm": "positive", "effect": 2.0, "phylop_mean": None},
    ]

    out = sv.read(rows)

    assert out["elements_read"] == 0  # one usable window cannot be ranked


def test_an_element_the_effect_floor_flattened_drops_out_rather_than_scoring_zero() -> None:
    """The run's own surprise: 58% of windows named no target, so their effect is exactly 0.0.

    `MIN_EFFECT` is 0.1, a floor built for deleting a whole element, and a 25 bp substitution often
    falls under it. An element whose twenty windows all fall under it has no order to correlate. It
    must leave the reading, not enter it as a zero correlation: 'the model said nothing here' and 'the
    model's sensitivity is unrelated to constraint here' are different statements, and averaging the
    second over elements that only did the first would be reading the floor rather than the panel.
    """
    flat = [
        {"element": "e1", "arm": "positive", "effect": 0.0, "phylop_mean": p}
        for p in (0.5, 1.0, 1.5, 2.0, 2.5)
    ]

    out = sv.read(flat)

    assert out["elements_read"] == 0
    assert out["mean_spearman"]["positive"] is None


def test_the_run_of_2026_09_17_reads_as_failure_by_the_registered_bars() -> None:
    """The bars were set before the number; this holds them to it after.

    +0.0453 at one-sided p 0.1741 over 3,200 windows. A later hand that widened the weak band would
    turn this run's recorded difference into a claim it never earned, so the recorded difference is
    what the bars are tested against.
    """
    got = sv.verdict(0.0453, 34, 200)

    assert got.startswith("failure")
    assert "section 15" in got
