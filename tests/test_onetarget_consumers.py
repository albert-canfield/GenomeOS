# SPDX-License-Identifier: AGPL-3.0-or-later
"""The one-target consumers' window reading (2026-09-27, third registration in docs/ATTRIBUTION.md)."""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

from genomeos.attribution import syntax_tiling as st
from genomeos.attribution import unknown_scoring as us
from genomeos.attribution.targets import genes_at_bar
from genomeos.predict.enhancer_target import MIN_EFFECT, predict_target

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("cut", ROOT / "scripts" / "constrained_unknown_targets.py")
cut = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cut)


def gene(name: str, drop: float, rise: float) -> dict:
    return {
        "gene": name,
        "max_drop_log2fc": drop,
        "max_rise_log2fc": rise,
        "max_drop_tissue": "t",
        "max_rise_tissue": "t",
    }


def random_record(rng: random.Random) -> dict:
    return {
        "genes": [
            gene(f"G{i}", -abs(rng.gauss(0, 0.12)), abs(rng.gauss(0, 0.08)))
            for i in range(rng.randrange(0, 9))
        ]
    }


def test_the_window_head_is_predict_targets_gene_and_the_yes_or_no_never_differs():
    rng = random.Random(7)
    for _ in range(2000):
        rec = random_record(rng)
        head = predict_target(rec["genes"])
        got = genes_at_bar(rec, MIN_EFFECT)
        assert bool(got) == (head is not None)
        if head:
            assert got[0][0] == head["gene"]
            assert got[0][1] == head["log2_fold_change"]


def test_an_uncached_element_is_none_not_an_empty_list():
    assert genes_at_bar(None) is None
    assert st.genes_at_bar_of(None) is None
    assert genes_at_bar({"genes": []}) == []


def test_syntax_tiling_moved_agrees_with_some_gene_at_the_bar():
    # the registered falsifier for syntax_tiling: a cached element where `moved` and "some gene at the
    # bar" disagree
    rng = random.Random(11)
    for _ in range(2000):
        rec = random_record(rng)
        row = {"predicted": predict_target(rec["genes"])}
        assert st.moved(row) == bool(st.genes_at_bar_of(rec))


def test_the_lifted_control_draws_the_windows_unknown_scoring_draws(monkeypatch):
    rng = random.Random(3)
    blocks = []
    pos = 0
    for _ in range(40):
        pos += rng.randrange(5_000, 60_000)
        length = rng.randrange(500, 20_000)
        blocks.append({"start": pos, "end": pos + length, "length": length})
        pos += length
    rows = []
    for i in range(600):
        s = rng.randrange(0, pos)
        rows.append({"id": f"e{i}", "start": s, "end": s + 300, "names_coding": rng.random() < 0.4})
    monkeypatch.setattr(us, "unknown_blocks", lambda chrom, results_dir=None: blocks)
    want = us.matched_random_windows("chrT", rows)
    got = cut.matched_random_windows(blocks, blocks, rows, {"q": lambda r: r["names_coding"]}, seed=us.SEED)
    assert got["windows_drawn"] == want["windows_drawn"]
    assert got["windows_carrying_an_element"] == want["windows_carrying_an_element"]
    assert got["blocks_carrying_an_element"] == want["unknown_blocks_carrying_an_element"]
    assert got["undrawable_blocks"] == want["undrawable_blocks"]
    assert got["q"]["blocks"] == want["from_the_named_gene"]["unknown_blocks_naming_a_coding_gene"]
    assert got["q"]["windows"] == want["from_the_named_gene"]["random_windows_naming_a_coding_gene"]


def test_pooled_reads_the_registered_band():
    per = {
        "a": {
            "blocks": 1,
            "windows_drawn": 50,
            "undrawable_blocks": 0,
            "windows_carrying_an_element": 100,
            "blocks_carrying_an_element": 10,
            "q": {"blocks": 4, "windows": 60},
        }
    }
    got = cut.pooled(per, ["q"])["q"]
    assert got["difference_in_points"] == -20.0 and got["reading"] == "below_chance"
    per["a"]["q"] = {"blocks": 6, "windows": 64}
    assert cut.pooled(per, ["q"])["q"]["reading"] == "at_chance"
