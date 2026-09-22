# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration that licenses the reader fix, and the historical loop that checks the harness."""

from __future__ import annotations

from genomeos.benchmark import loci_reread as rr

CCDC26_ROW = {
    "predicted": {"gene": "CCDC26", "log2_fold_change": -1.691, "action": "activates", "tissue": "K562"},
    "predicted_coding": {
        "gene": "GSDMC",
        "log2_fold_change": -0.1433,
        "action": "activates",
        "tissue": "K562",
    },
    "run": "enhancer_targets_all",
}


def test_the_registration_states_what_it_has_to_state_before_anything_is_re_read():
    """A registration missing any of these is not one, and the fix is not licensed by it."""
    p = rr.PREREGISTRATION
    for key in (
        "written",
        "the_two_changes",
        "registered_together_measured_apart",
        "what_the_reader_will_do",
        "why_this_reading_and_not_the_favourable_one",
        "scope_what_is_not_being_changed",
        "expected_direction_per_section",
        "known_before_the_registration",
        "falsifiers",
        "if_the_fourth_frames_0.200_RISES",
        "if_the_fix_is_wrong",
        "budget",
    ):
        assert p.get(key), f"the registration says nothing about {key}"
    for name, section in rr.FRAMES.items():
        number = section.split(" - ")[0]
        assert any(k.startswith(f"{number} ") for k in p["expected_direction_per_section"]), (
            f"{name} is re-read and no direction was registered for section {number}"
        )
    for key in ("the harness", "the mechanism", "double counting", "the claim tables", "unfalsifiability"):
        assert p["falsifiers"].get(key)


def test_the_registered_budget_is_zero_and_says_so_without_an_estimate():
    assert "0 AlphaGenome requests" in rr.PREREGISTRATION["budget"]


def test_the_historical_loop_reproduces_the_defect_it_is_kept_to_check():
    """The pre-2026-09-21 loop must still hide CCDC26 behind GSDMC, whatever the shipped reader does.

    It is the harness's fidelity check: if this loop reproduces a saved run's stored readings, the
    re-read is reading the rows that run read. It is therefore frozen and must never be repaired.
    """
    out = rr.historical_deletion_reading([CCDC26_ROW])
    assert out["target"] == "GSDMC"
    assert [g["gene"] for g in out["targets"]] == ["GSDMC"]
    assert out["elements_scored"] == 1


def test_the_historical_loop_still_falls_through_when_no_coding_gene_is_named():
    out = rr.historical_deletion_reading([{**CCDC26_ROW, "predicted_coding": None}])
    assert out["target"] == "CCDC26"
