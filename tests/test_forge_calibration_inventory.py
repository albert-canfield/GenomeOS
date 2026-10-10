# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area H, 2026-10-10: the registration that decides whether BioForge's design answers have a calibration
record, pinned before the census that applies it (genomeos/forge/calibration_inventory.py)."""

from __future__ import annotations

import pytest

from genomeos.forge.calibration_inventory import (
    MAX_INTERVAL_WIDTH,
    MINIMUM_INDEPENDENT_TRIPLES,
    REGISTRATION,
    clopper_pearson,
    minimum_triples,
    worst_width,
)


def test_the_exact_interval_reproduces_the_published_numbers() -> None:
    """The same [0.3976, 1] and [0.1581, 1] that docs/BIOFORGE-CONFIDENCE.md §4(c) printed at 4/4 and 2/2."""
    assert clopper_pearson(4, 4) == pytest.approx((0.3976, 1.0), abs=5e-5)
    assert clopper_pearson(2, 2) == pytest.approx((0.1581, 1.0), abs=5e-5)
    lo, hi = clopper_pearson(19, 38)
    assert lo == pytest.approx(1 - hi, abs=1e-9) and 1 / 3 < lo < hi < 2 / 3
    with pytest.raises(ValueError):
        clopper_pearson(3, 2)


def test_the_minimum_is_the_smallest_n_whose_worst_interval_is_narrower_than_a_third() -> None:
    assert MAX_INTERVAL_WIDTH == 1 / 3
    assert minimum_triples() == MINIMUM_INDEPENDENT_TRIPLES == 38
    assert worst_width(38) < 1 / 3 <= worst_width(37)
    assert REGISTRATION["the_minimum"]["n"] == MINIMUM_INDEPENDENT_TRIPLES


def test_the_registration_names_every_circularity_rule_and_refuses_a_rate_on_circular_triples() -> None:
    assert set(REGISTRATION["circular_if_any"]) == {
        "same_source",
        "acceptance_check",
        "authored_frame",
        "evidence_is_the_model",
    }
    assert "compute no rate" in REGISTRATION["decision"]
    assert "never counted as an independent outcome" in REGISTRATION["never"]
