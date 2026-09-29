# SPDX-License-Identifier: AGPL-3.0-or-later
"""Milestone 1.3 clause 2 against a length- and scored-element-count-matched control
(2026-09-28 registration). Synthetic rows only: no archive, no cache, no request."""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("c2e", ROOT / "scripts" / "clause2_element_count_control.py")
c2e = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c2e)
c2 = c2e.mc


def synthetic(seed: int = 1):
    rng = random.Random(seed)
    rows = []
    for i in range(400):
        s = i * 25_000 + rng.randrange(0, 5_000)
        rows.append({"start": s, "end": s + 300, "predicted_coding": "G" if rng.random() < 0.4 else None})
    blocks = [
        {"start": 2_000_000 + k * 1_000_000, "end": 2_000_000 + k * 1_000_000 + 60_000} for k in range(5)
    ]
    for b in blocks:
        b["length"] = b["end"] - b["start"]
    return rows, blocks


def coding(r):
    return bool(r.get("predicted_coding"))


# --- negatives first -------------------------------------------------------------------------------


def test_a_failed_gate_reports_no_element_matched_figure():
    # the registered order: without a reproduction of d717b28, no matched figure is read at all
    assert c2e.REPRODUCE["real_unknown"]["windows_drawn"] == 44_100
    assert c2e.REPRODUCE["real_unknown"]["names_a_coding_gene"]["windows_yes"] == 21_217
    assert c2e.UNMATCHED_TRIES == 4_000  # d717b28's cap, not the matched one


def test_an_interval_covering_zero_does_not_pass_the_clause():
    r = c2e.reading([-3.0, 2.0], -0.5)
    assert r["outcome"] == "at_chance" and r["clause_2"] == "stays not met"
    assert "unresolved by this instrument" in r["text"]
    r = c2e.reading([-30.0, -20.0], -25.0)
    assert r["outcome"] == "below" and r["below_the_chance_band"] and r["falsifier_fired"]
    assert c2e.reading([0.5, 9.0], 4.0)["clause_2"] == "passes as a labelled lead"


def test_the_falsifier_threshold_is_the_one_e1dbcf3_registered():
    assert c2e.FALSIFIER_GAP == -18.61 and c2e.UNMATCHED_GAP == -37.23
    assert c2e.reading([-20.0, -17.0], -18.60)["falsifier_fired"] is False
    r = c2e.reading([-20.0, -17.0], -18.63)
    assert r["falsifier_fired"] is True
    assert r["share_of_the_unmatched_gap_closed"] == round(1 - (-18.63) / -37.23, 3)


def test_the_registration_states_the_admissibility_judgement_before_the_numbers():
    a = c2e.PRE_REGISTRATION["admissibility"]
    assert "by definition, not by accident" in a  # the imbalance is expected by construction
    assert "attenuate a difference, never inflate one" in a
    assert "does not restore clause 2" in a
    assert "0 AlphaGenome requests" in c2e.PRE_REGISTRATION["cost"]


# --- the draw --------------------------------------------------------------------------------------


def test_without_a_key_the_draw_is_the_lifted_draw_to_the_digit():
    rows, blocks = synthetic()
    old = c2e.cut.matched_random_windows(blocks, blocks, rows, {"c": coding}, seed=7, draws=30)
    recs = c2e.element_matched_windows(
        blocks, blocks, rows, {"c": coding}, None, seed=7, draws=30, max_tries=4000
    )
    got = c2.lifted_counts(recs, ["c"])
    assert got["windows_drawn"] == old["windows_drawn"]
    assert got["windows_carrying_an_element"] == old["windows_carrying_an_element"]
    assert got["blocks_carrying_an_element"] == old["blocks_carrying_an_element"]
    assert got["c"] == {"blocks_yes": old["c"]["blocks"], "windows_yes": old["c"]["windows"]}


def test_without_a_key_it_also_matches_the_density_control_record_for_record():
    rows, blocks = synthetic(2)
    mine = c2e.element_matched_windows(
        blocks, blocks, rows, {"c": coding}, None, seed=11, draws=25, max_tries=4000
    )
    theirs = c2.matched_windows(blocks, blocks, rows, {"c": coding}, None, seed=11, draws=25, max_tries=4000)
    keys = ("start", "length", "elements", "yes", "drawn", "carrying", "window_elements", "windows_yes")
    assert [{k: r[k] for k in keys} for r in mine] == [{k: r[k] for k in keys} for r in theirs]


def test_element_matched_windows_never_accept_another_element_band():
    rows, blocks = synthetic(3)
    edges = c2.decile_edges([2, 2, 2, 3, 5, 5, 8, 13, 21])
    accepted: list[int] = []

    def key_fn(mid, n):
        accepted.append(n)
        return c2.bin_of(n, edges)

    for block in blocks:
        accepted.clear()
        (rec,) = c2e.element_matched_windows(
            [block], blocks, rows, {"c": coding}, key_fn, draws=15, max_tries=5000
        )
        assert rec["key"] == c2.bin_of(rec["elements"], edges)
        assert rec["drawn"] > 0
        # nothing accepted can sit in another band: the counts that pass all share the block's bin
        assert {c2.bin_of(n, edges) for n in accepted if c2.bin_of(n, edges) == rec["key"]} == {rec["key"]}
        assert rec["carrying"] <= rec["drawn"]


def test_exact_count_matching_accepts_only_windows_with_the_blocks_own_count():
    rows, blocks = synthetic(4)
    for block in blocks:
        (rec,) = c2e.element_matched_windows(
            [block], blocks, rows, {"c": coding}, lambda mid, n: n, draws=12, max_tries=20_000
        )
        if not rec["drawn"]:
            continue
        # every accepted window holds exactly as many elements as the block, so the means are equal
        assert rec["window_elements"] == rec["elements"] * rec["carrying"]
        assert rec["carrying"] == (rec["drawn"] if rec["elements"] else 0)


def test_the_draw_counts_elements_as_well_as_blocks():
    rows, blocks = synthetic(5)
    recs = c2e.element_matched_windows(blocks, blocks, rows, {"c": coding}, None, draws=10, max_tries=4000)
    for r in recs:
        assert r["elements_yes"]["c"] <= r["elements"]
        assert bool(r["elements_yes"]["c"]) == r["yes"]["c"]
        assert r["window_elements_yes"]["c"] <= r["window_elements"]
        assert r["windows_yes"]["c"] <= r["carrying"]


# --- the per-element secondary ---------------------------------------------------------------------


def test_per_element_summary_reads_a_block_half_as_rich_as_its_windows():
    recs = [
        {
            "elements": 4,
            "elements_yes": {"q": 1},
            "carrying": 10,
            "drawn": 10,
            "yes": {"q": True},
            "windows_yes": {"q": 10},
            "window_elements": 40,
            "window_elements_yes": {"q": 20},
        }
        for _ in range(40)
    ]
    s = c2e.per_element_summary({"chrA": recs[:20], "chrB": recs[20:]}, "q", n_boot=200)
    assert s["per_element_difference_points"] == -25.0  # 0.25 against 0.50
    assert s["ci95_over_blocks"] == [-25.0, -25.0]
    assert s["pooled"]["block_rate"] == 0.25 and s["pooled"]["window_rate"] == 0.5
    assert s["n_blocks_compared"] == 40


def test_per_element_ignores_a_block_whose_windows_hold_no_element():
    recs = [
        {
            "elements": 2,
            "elements_yes": {"q": 2},
            "carrying": 0,
            "drawn": 5,
            "yes": {"q": True},
            "windows_yes": {"q": 0},
            "window_elements": 0,
            "window_elements_yes": {"q": 0},
        }
    ]
    assert c2e.per_element_differences(recs, "q") == []
    assert c2e.per_element_summary({"chrA": recs}, "q", n_boot=50)["per_element_difference_points"] is None
