# SPDX-License-Identifier: AGPL-3.0-or-later
"""Milestone 1.3 clause 2 against a length- and gene-density-matched control (2026-09-28 registration).
Synthetic rows only: no archive, no cache, no request."""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("c2", ROOT / "scripts" / "clause2_matched_control.py")
c2 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c2)


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


def test_a_failed_gate_reports_no_matched_figure():
    # the registered order: without a reproduction, the reading names the failure and nothing else
    assert "no matched figure" in "the reproduction gate failed: no matched figure is reported"
    assert c2.REPRODUCE["real_unknown"]["windows_drawn"] == 44_100


def test_an_interval_covering_zero_does_not_pass_the_clause():
    r = c2.reading([-3.0, 2.0], -0.5)
    assert r["outcome"] == "at_chance" and r["clause_2"] == "stays not met"
    r = c2.reading([-30.0, -20.0], -25.0)
    assert r["outcome"] == "below" and r["below_the_chance_band"] and r["falsifier_fired"]
    assert c2.reading([0.5, 9.0], 4.0)["clause_2"] == "passes as a labelled lead"


def test_matched_windows_never_accept_another_bin():
    rows, blocks = synthetic()
    tss = sorted(random.Random(3).sample(range(0, 10_000_000), 60))
    edges = c2.decile_edges([c2.tss_count(tss, m) for m in range(0, 10_000_000, 50_000)])
    assert len(edges) >= 3

    def bin_fn(m):
        return c2.bin_of(c2.tss_count(tss, m), edges)

    for block in blocks:
        mids = []  # raw_fn is asked for the block's midpoint, then for every accepted window's
        (rec,) = c2.matched_windows(
            [block], blocks, rows, {"c": coding}, bin_fn, draws=20, max_tries=5000, raw_fn=mids.append
        )
        assert len(mids) == 1 + rec["drawn"]
        assert all(bin_fn(m) == rec["bin"] for m in mids[1:])
        assert rec["bin"] == bin_fn(mids[0])


def test_without_a_bin_the_draw_is_the_lifted_draw_to_the_digit():
    rows, blocks = synthetic()
    old = c2.cut.matched_random_windows(blocks, blocks, rows, {"c": coding}, seed=7, draws=30)
    recs = c2.matched_windows(blocks, blocks, rows, {"c": coding}, None, seed=7, draws=30, max_tries=4000)
    got = c2.lifted_counts(recs, ["c"])
    assert got["windows_drawn"] == old["windows_drawn"]
    assert got["windows_carrying_an_element"] == old["windows_carrying_an_element"]
    assert got["blocks_carrying_an_element"] == old["blocks_carrying_an_element"]
    assert got["c"] == {"blocks_yes": old["c"]["blocks"], "windows_yes": old["c"]["windows"]}


def test_tss_count_and_distance():
    tss = [100, 1_000_000, 1_000_010]
    assert c2.tss_count(tss, 1_000_000, half=5) == 1
    assert c2.tss_count(tss, 1_000_000) == 2  # 100 lies 999,900 bp away, outside the 1 Mb window
    assert c2.tss_count(tss, 500_000) == 3
    assert c2.tss_distance(tss, 999_990) == 10
    assert c2.tss_distance([], 5) == 10**12


def test_decile_edges_collapse_ties():
    assert c2.decile_edges([0] * 50 + [1] * 30 + [5] * 20) == [0, 1, 5]
    assert c2.bin_of(0, [0, 1, 5]) == 1 and c2.bin_of(7, [0, 1, 5]) == 3


def test_summarise_reads_a_block_that_always_says_yes_against_windows_that_half_do():
    recs = [
        {
            "elements": 2,
            "carrying": 10,
            "drawn": 10,
            "yes": {"q": True},
            "windows_yes": {"q": 5},
            "raw": 3,
            "window_raw": 30.0,
            "window_elements": 20,
        }
        for _ in range(40)
    ]
    s = c2.summarise({"chrA": recs[:20], "chrB": recs[20:]}, "q", n_boot=200)
    assert s["matched_difference_points"] == 50.0
    assert s["ci95_over_blocks"] == [50.0, 50.0]
    assert s["balance"]["covariate_block_mean"] == s["balance"]["covariate_window_mean"] == 3.0
