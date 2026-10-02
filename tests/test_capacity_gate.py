# SPDX-License-Identifier: AGPL-3.0-or-later
"""scripts/capacity_gate.py: refuses while the machine pages, and never passes on a reading it lacks.

The decision is a pure function so both directions can be planted. Three wrong signals preceded it on
2026-10-02 -- free pages alone, free plus inactive, and the compressor -- so the tests pin the signal
that was finally chosen (the RATE of swapouts) and pin that an unreadable signal is not a yes, which
is the mistake a gate makes when it is written to be convenient.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import capacity_gate as cg  # noqa: E402

QUIET = cg.Reading(swapouts=100, pressure_free_percent=69, pressure_level=None, free_disk_gb=50.0)


def test_a_quiet_machine_may_start_a_heavy_job():
    ok, reasons = cg.decide(100, 100, QUIET)
    assert ok, reasons
    assert reasons == []


def test_a_machine_that_paged_between_the_samples_is_refused():
    ok, reasons = cg.decide(100, 137, QUIET)
    assert not ok
    assert any("paging" in r and "100 -> 137" in r for r in reasons), reasons


def test_memory_pressure_at_warn_or_critical_is_refused():
    for level in ("warn", "critical"):
        ok, reasons = cg.decide(100, 100, cg.Reading(100, 5, level, 50.0))
        assert not ok, level
        assert any(level in r for r in reasons), reasons


def test_free_disk_below_the_floor_is_refused_because_swap_grows_on_disk():
    ok, reasons = cg.decide(100, 100, cg.Reading(100, 69, None, cg.MIN_FREE_DISK_GB - 0.1))
    assert not ok
    assert any("free disk" in r and "swap files on disk" in r for r in reasons), reasons


@pytest.mark.parametrize(
    "before,after,reading,missing",
    [
        (None, 100, QUIET, "swapouts"),
        (100, None, QUIET, "swapouts"),
        (100, 100, cg.Reading(100, 69, None, None), "free disk"),
        (100, 100, cg.Reading(100, None, None, 50.0), "memory_pressure"),
    ],
)
def test_an_unreadable_signal_is_never_a_yes(before, after, reading, missing):
    ok, reasons = cg.decide(before, after, reading)
    assert not ok, (reasons, missing)
    assert any(missing in r for r in reasons), (reasons, missing)


def test_an_absent_pressure_level_alone_is_normal_on_this_machine():
    # memory_pressure here prints a free percentage and no level line; that must not read as a failure
    ok, reasons = cg.decide(100, 100, cg.Reading(100, 69, None, 50.0))
    assert ok, reasons


def test_the_swapouts_field_is_parsed_from_real_vm_stat_text():
    text = "Pages free:   100.\nSwapins:    389256.\nSwapouts:     833625.\n"
    assert cg.vm_swapouts(text) == 833625
    assert cg.vm_swapouts("Pages free: 1.\n") is None


def test_the_refusal_prints_every_figure_so_the_number_can_be_argued_with():
    line = cg.figures(100, 137, cg.Reading(137, 12, "warn", 3.5))
    for fragment in ("swapouts 100 -> 137", "rose 37", "warn", "12%", "3.5 GB"):
        assert fragment in line, (fragment, line)
