# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Gate (b) asked of the per-demander `saturable` policy: the planted contrasts and the guards.

The planted transcriptomes are the registration's own (scripts/saturable_gate.py, REGISTRATION
["planted_contrasts"]): where the truth is per-gene saturation, `saturable` must find it and pass while
`proportional` and `competitive` cannot move the error at all; where nothing saturates it must not
manufacture a pass; and where the truth is a straight compressed line the comparator must say so.
The circularity guard is tested by changing every held-out protein value and asking whether any fitted
number moved.
"""

from __future__ import annotations

import hashlib
import json
import math

import numpy as np
import pytest

from scripts.saturable_gate import (
    GENE_BLIND_TOLERANCE_LOG10,
    LOG10_K_GRID,
    PLANTED_K,
    REGISTRATION,
    is_held_out,
    planted,
    run_gate,
    verdict,
)


@pytest.fixture(scope="module")
def saturable_truth() -> dict:
    return run_gate(*planted("saturable_truth"))


def test_the_split_rule_is_the_registered_hash_and_halves_the_genes() -> None:
    symbols = planted("null_truth")[0]
    held = sum(is_held_out(s) for s in symbols)
    assert 0.45 < held / len(symbols) < 0.55
    for s in ("ACTB", "GAPDH", "TP53", "MYC"):
        assert is_held_out(s) == (int(hashlib.sha256(s.encode()).hexdigest()[:8], 16) % 2 == 1)


def test_planted_saturation_is_found_and_passes(saturable_truth: dict) -> None:
    key, reading = verdict(saturable_truth, secondary_in_band=True)
    assert key == "PASS", reading
    assert abs(saturable_truth["fitted_on_training"]["saturable_log10_k"] - math.log10(PLANTED_K)) <= 0.3
    assert saturable_truth["instrument"]["runtime_vs_closed_form_max_log10_gap"] <= 1e-12


def test_the_gene_blind_policies_cannot_move_the_error_where_the_answer_is_saturation(
    saturable_truth: dict,
) -> None:
    """The contrast: the old policies fail where the new one passes, through the same runtime."""
    for policy, row in saturable_truth["gene_blind_arms"].items():
        assert row["distinct_factors"] == 1, policy
        assert row["log10_factor"] < 0.0, f"{policy}: the pool must be short or this proves nothing"
        assert abs(row["improvement_over_baseline"]) <= GENE_BLIND_TOLERANCE_LOG10, policy


def test_a_null_truth_does_not_manufacture_a_pass() -> None:
    g = run_gate(*planted("null_truth"))
    key, _ = verdict(g, secondary_in_band=True)
    assert key == "NO_IMPROVEMENT_DETECTED"
    assert g["held_out"]["improvement"] < 0.05


def test_a_straight_compressed_truth_is_named_by_the_comparator() -> None:
    """An improvement over the baseline is not, alone, evidence for the saturable form."""
    g = run_gate(*planted("straight_truth"))
    assert g["comparator_key"] == "comparator_line_better"
    assert g["held_out"]["improvement"] > 0.0, "the baseline is worse still, which is the point"


def test_no_held_out_protein_value_enters_any_fit() -> None:
    """The circularity guard: replace every held-out outcome and nothing fitted may move."""
    symbols, tpm, y = planted("saturable_truth")
    g1 = run_gate(symbols, tpm, y)
    held = np.array([is_held_out(s) for s in symbols])
    scrambled = y.copy()
    scrambled[held] = np.random.default_rng(1).normal(9.0, 4.0, size=int(held.sum()))
    g2 = run_gate(symbols, tpm, scrambled)
    assert g1["fitted_on_training"] == g2["fitted_on_training"]
    assert g1["held_out"] != g2["held_out"], "the held-out readings must see the change"


def test_the_grid_and_the_registered_readings_are_what_the_registration_says() -> None:
    assert (LOG10_K_GRID[0], LOG10_K_GRID[-1], len(LOG10_K_GRID)) == (-1.0, 7.0, 801)
    for key in ("PASS", "NO_IMPROVEMENT_DETECTED", "INSTRUMENT_FAILED"):
        assert key in REGISTRATION["readings"]
    assert "no protein value of a held-out gene enters any fit" in REGISTRATION["circularity_guard"]
    assert "explains" not in json.dumps(REGISTRATION["readings"])
