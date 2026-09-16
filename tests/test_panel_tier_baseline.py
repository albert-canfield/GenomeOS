# SPDX-License-Identifier: AGPL-3.0-or-later
"""The check that asks whether the panel's matched background is a neutral baseline."""

from __future__ import annotations

import importlib.util
from pathlib import Path

SPEC = importlib.util.spec_from_file_location(
    "panel_tier_baseline", Path(__file__).resolve().parents[1] / "scripts" / "panel_tier_baseline.py"
)
assert SPEC and SPEC.loader
baseline = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(baseline)


def test_a_tier_above_one_everywhere_is_not_a_coin_toss() -> None:
    row = baseline.sign_test_above_one([1.1] * 22)
    assert row == {"chromosomes": 22, "above_one": 22, "one_sided_p": round(1 / 2**22, 6)}


def test_a_tier_at_the_background_is_unsurprising() -> None:
    row = baseline.sign_test_above_one([1.1] * 11 + [0.9] * 11)
    assert row["above_one"] == 11
    assert row["one_sided_p"] > 0.4


def test_exactly_one_is_not_above_one() -> None:
    """A ratio of exactly the background is the null, not evidence for it."""
    assert baseline.sign_test_above_one([1.0] * 5)["above_one"] == 0


def test_tiers_are_reported_against_both_baselines(monkeypatch) -> None:
    """The point of the check: the same tier read against the background and against the neutral tier."""
    monkeypatch.setattr(
        baseline,
        "read_chromosomes",
        lambda: {
            "chr1": {"cds": 0.4, "constrained_unknown": 0.9, "neutral": 1.2},
            "chr2": {"cds": 0.4, "constrained_unknown": 1.1, "neutral": 1.2},
        },
    )
    out = baseline.collect()

    background = out["against_matched_background"]
    assert background["constrained_unknown"]["above_one"] == 1
    assert background["neutral"]["above_one"] == 2
    # 0.9/1.2 and 1.1/1.2 are both below the neutral tier, where against the background one was above
    assert out["against_the_neutral_tier"]["constrained_unknown"]["below_neutral"] == 2
    assert "neutral" not in out["against_the_neutral_tier"]
    assert out["chromosomes"] == ["chr1", "chr2"]


def test_chry_is_excluded_like_everywhere_else_in_the_panel() -> None:
    assert "chrY" in baseline.EXCLUDED
