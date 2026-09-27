# SPDX-License-Identifier: AGPL-3.0-or-later
"""What the whole-chromosome scoring bought the 98%: the joins, the strata and the ablation's rebuild."""

import json

import pytest

from genomeos.attribution.unknown_scoring import (
    _rebuild,
    add_strata,
    block_level,
    coverage,
    gc_fraction,
    rho_p,
    stratified,
    unknown_blocks,
)


def test_gc_fraction_ignores_ambiguous_bases():
    assert gc_fraction("GGCCAATT") == 0.5
    assert gc_fraction("NNNN") is None


def test_unknown_blocks_join_tier_constraint_and_case(tmp_path):
    (tmp_path / "unknown_chrT.json").write_text(
        json.dumps({"blocks": [{"start": 100, "end": 600, "length": 500, "class": "unique_intergenic"}]})
    )
    (tmp_path / "budget_chrT.json").write_text(
        json.dumps(
            {
                "blocks": [
                    {
                        "start": 100,
                        "end": 600,
                        "phylop": {"fraction_above": 0.07},
                        "guess": {"tier": "constrained_unknown"},
                    }
                ]
            }
        )
    )
    (tmp_path / "variation_chrT.json").write_text(
        json.dumps(
            {
                "blocks": [
                    {"start": 100, "end": 600, "case": {"case": "syntax"}, "gnocchi": {"fraction_above": 0.5}}
                ]
            }
        )
    )
    b = unknown_blocks("chrT", tmp_path)[0]
    assert (b["tier"], b["case"], b["mammal_fraction"], b["human_fraction"]) == (
        "constrained_unknown",
        "syntax",
        0.07,
        0.5,
    )


def _row(flag, metric, stratum=(0, 0, 0), **kw):
    return {"in_unknown": flag, "metric": metric, "stratum": stratum, **kw}


def test_stratified_reports_no_difference_when_there_is_none():
    rows = [_row(i % 2 == 0, 1.0 if i % 3 else 0.0) for i in range(60)]
    out = stratified(rows, "metric")
    assert out["unknown_n"] == 30 and out["rest_n"] == 30
    assert abs(out["difference_matched"]) < 0.3 and out["p_permuted_within_strata"] > 0.2


def test_stratified_finds_a_difference_inside_shared_strata_and_drops_lonely_strata():
    rows = [_row(True, 1.0) for _ in range(30)] + [_row(False, 0.0) for _ in range(30)]
    rows.append(_row(True, 1.0, stratum=(9, 9, 9)))  # no control in this stratum: left out
    out = stratified(rows, "metric")
    assert out["difference_matched"] == 1.0
    assert out["p_permuted_within_strata"] < 0.01
    assert out["unknown_in_shared_strata"] == 30 and out["strata_shared"] == 1


def test_add_strata_bins_by_length_gc_and_distance():
    rows = [
        {"length": n, "gc": 0.4 + n / 1000, "nearest_coding_tss": 1000 * n, "in_unknown": n % 2 == 0}
        for n in range(1, 41)
    ]
    add_strata(rows)
    assert len({tuple(r["stratum"]) for r in rows}) > 5
    assert all(len(r["stratum"]) == 3 for r in rows)


def test_block_level_uses_blocks_as_the_unit():
    # one constrained block whose elements all move, ten others whose elements never do
    blocks = {(0, 100): [{"_mammal": True, "moves_gene": 1.0} for _ in range(50)]}
    for i in range(10):
        blocks[(1000 * (i + 1), 1000 * (i + 1) + 100)] = [{"_mammal": False, "moves_gene": 0.0}]
    out = block_level(blocks, "moves_gene")
    assert out["constrained_blocks"] == 1 and out["other_blocks"] == 10
    assert out["difference"] == 1.0
    # fifty elements inside one block are still one unit, so the p cannot be small
    assert out["p_permuted_over_blocks"] > 0.05


def test_rho_p_is_small_only_for_a_real_correlation():
    x = list(range(20))
    assert rho_p(x, x) == pytest.approx(1 / 2001, abs=1e-4)
    assert rho_p(x, [0] * 20) is None
    assert rho_p(x, [3, 1, 2] * 6 + [2, 1]) > 0.05


class _Peaks:
    """A stand-in for the reader's peak index: everything overlapping the given window is open."""

    def __init__(self, window):
        self.window = window

    def overlapping(self, start, end):
        a, b = self.window
        return [(a, b, 1.0)] if start < b and end > a else []


def test_rebuild_drops_the_elements_the_keep_rule_refuses():
    closure = {
        "cells": ["K562"],
        "genes": [
            {
                "gene": "G",
                "cells": {
                    "K562": {
                        "expression": 0.9,
                        "expressed": True,
                        "promoter_open": True,
                        "promoter_percentile": 0.5,
                    }
                },
            }
        ],
    }
    by_gene = {
        "G": [
            {"id": "in", "start": 10, "end": 20, "effect": 0.5, "by_cell": {"K562": 0.4}},
            {"id": "out", "start": 10, "end": 20, "effect": 0.5, "by_cell": {"K562": 0.2}},
        ]
    }
    indexes = {"K562": _Peaks((0, 100))}
    full = _rebuild(closure, by_gene, indexes, lambda e: True)[0]["cells"]["K562"]
    assert full["input_both"] == pytest.approx(0.6) and full["scored_elements"] == 2
    one = _rebuild(closure, by_gene, indexes, lambda e: e["id"] == "in")[0]["cells"]["K562"]
    assert one["input_both"] == pytest.approx(0.4) and one["active_elements"] == 1
    none = _rebuild(closure, by_gene, indexes, lambda e: False)[0]
    assert none["elements"] == 0 and none["cells"]["K562"]["input_both"] is None
    # an element outside the cell's peaks is scored but not active, so it does not enter the input
    closed = {"K562": _Peaks((500, 600))}
    shut = _rebuild(closure, by_gene, closed, lambda e: True)[0]["cells"]["K562"]
    assert shut["input_both"] == 0 and shut["active_elements"] == 0


def test_coverage_counts_blocks_before_and_after_the_sweep(tmp_path):
    (tmp_path / "unknown_chrT.json").write_text(
        json.dumps(
            {
                "blocks": [
                    {"start": 0, "end": 1000, "length": 1000, "class": "unique_intergenic"},
                    {"start": 2000, "end": 3000, "length": 1000, "class": "unique_intergenic"},
                ]
            }
        )
    )
    (tmp_path / "budget_chrT.json").write_text(
        json.dumps(
            {
                "blocks": [
                    {"start": 0, "end": 1000, "guess": {"tier": "regulatory"}, "phylop": {}},
                    {"start": 2000, "end": 3000, "guess": {"tier": "regulatory"}, "phylop": {}},
                ]
            }
        )
    )
    for name in ("constrained_targets_chrT", "enhancer_targets_chrT"):
        (tmp_path / f"{name}.json").write_text(
            json.dumps(
                {
                    "elements": [
                        {"id": "old", "start": 100, "end": 200, "predicted_coding": {"gene": "A"}},
                    ]
                }
            )
        )
    rows = [
        {
            "id": "old",
            "block": [0, 1000],
            "unknown_bp": 100,
            "names_coding": True,
            "in_unknown": True,
            "moves_gene": True,
            "cells_acting": 0,
            "tier": "regulatory",
        },
        {
            "id": "new",
            "block": [2000, 3000],
            "unknown_bp": 150,
            "names_coding": True,
            "in_unknown": True,
            "moves_gene": True,
            "cells_acting": 2,
            "tier": "regulatory",
        },
    ]
    out = coverage("chrT", rows, tmp_path)
    assert out["blocks_with_a_named_target_before"] == 1
    assert out["blocks_with_a_named_target_after"] == 2
    assert out["blocks_with_a_target_and_a_cell"] == 1
    assert out["element_bp_with_a_named_target_inside_unknown"] == 250
    assert out["element_bp_with_a_named_target_inside_unknown_before"] == 100


# ---- the window reading, 2026-09-27 -------------------------------------------------------


class _Responses:
    """A stand-in for `ElementResponses` holding one element's window."""

    def __init__(self, windows):
        self._w = windows

    def genes(self, chrom, element_id):
        return self._w.get((chrom, element_id))


def test_window_reading_names_the_silence_it_cannot_answer():
    from genomeos.attribution.unknown_scoring import NOT_READ, window_reading

    unread = window_reading(None, "chrT", "E1", set())
    assert unread["cache_silence"] == NOT_READ
    assert unread["head_abs_log2"] is None and unread["cells_acting_window"] is None
    missing = window_reading(_Responses({}), "chrT", "E1", set())
    assert missing["cache_silence"] and missing["cache_silence"] != NOT_READ
    # a silence is never a zero: nothing here may be read as "the model predicts no effect"
    assert missing["head_abs_log2"] is None


def test_window_reading_takes_the_head_over_every_gene_and_the_cells_over_every_gene():
    from genomeos.attribution.unknown_scoring import window_reading

    window = {
        ("chrT", "E1"): {
            # the head of the window on the max-across-tracks rule, but under the bar on every line
            "TOP": {
                "gene": "TOP",
                "max_drop_log2fc": -0.9,
                "max_rise_log2fc": 0.1,
                "by_cell": {"K562": -0.05, "HepG2": -0.02, "GM12878": None, "IMR-90": 0.0},
            },
            # not the head, but it reaches the bar on two lines' own tracks
            "OTHER": {
                "gene": "OTHER",
                "max_drop_log2fc": -0.3,
                "max_rise_log2fc": 0.0,
                "by_cell": {"K562": -0.25, "HepG2": 0.4, "GM12878": -0.01, "IMR-90": None},
            },
        }
    }
    got = window_reading(_Responses(window), "chrT", "E1", {"OTHER"})
    assert got["cache_silence"] is None
    assert got["window_genes"] == 2
    assert got["head_gene"] == "TOP" and got["head_abs_log2"] == 0.9
    assert got["head_signed_log2"] == -0.9
    # the coding head is the coding gene's, not the window's
    assert got["head_coding_gene"] == "OTHER" and got["head_coding_abs_log2"] == 0.3
    # two lines reach the bar through OTHER, which the compact table's top gene could not show
    assert got["cells_acting_window"] == 2
    assert got["window_by_cell"]["HepG2"] == 0.4
    # a gene not scored on a line's track does not hide another gene's answer on it
    assert got["window_by_cell"]["GM12878"] == -0.01
    # a scored 0.0 is an answer and stays one, where a None is a silence and stays None
    assert got["window_by_cell"]["IMR-90"] == 0.0


def test_the_censoring_check_reports_a_disagreement_when_there_is_one():
    from genomeos.attribution.unknown_scoring import the_table_was_not_censoring

    rows = [
        # the table named a coding gene and the window agrees at the same threshold
        {
            "id": "a",
            "cache_silence": None,
            "names_coding": True,
            "head_coding_abs_log2": 0.4,
            "head_coding_gene": "A",
        },
        # the table named none and the window's coding head is under the bar: they agree
        {
            "id": "b",
            "cache_silence": None,
            "names_coding": False,
            "head_coding_abs_log2": 0.05,
            "head_coding_gene": "B",
        },
        # the table named none and the window's coding head clears the bar: censoring
        {
            "id": "c",
            "cache_silence": None,
            "names_coding": False,
            "head_coding_abs_log2": 0.5,
            "head_coding_gene": "C",
        },
    ]
    got = the_table_was_not_censoring(rows)
    assert got["names_a_coding_gene_compact_table"] == 1
    assert got["names_a_coding_gene_window_at_the_same_threshold"] == 2
    assert got["disagreements"] == 1 and got["first_disagreements"] == ["c"]
    assert got["names_a_coding_gene_window_ungated"] == 3


def test_matched_random_windows_reject_windows_that_overlap_the_unknown_space(tmp_path):
    from genomeos.attribution.unknown_scoring import matched_random_windows

    (tmp_path / "unknown_chrT.json").write_text(
        json.dumps({"blocks": [{"start": 0, "end": 100, "length": 100, "class": "unique_intergenic"}]})
    )
    rows = [
        # inside the unknown block, and it names a coding gene
        {
            "id": "u",
            "start": 10,
            "end": 30,
            "block": [0, 100],
            "names_coding": True,
            "cells_acting": 1,
            "cells_acting_window": 1,
        },
        # outside it, naming nothing, so the random windows that can be drawn all name nothing
        {
            "id": "r",
            "start": 500,
            "end": 520,
            "block": None,
            "names_coding": False,
            "cells_acting": None,
            "cells_acting_window": 0,
        },
        {
            "id": "s",
            "start": 900,
            "end": 920,
            "block": None,
            "names_coding": False,
            "cells_acting": None,
            "cells_acting_window": 0,
        },
    ]
    got = matched_random_windows("chrT", rows, tmp_path)
    assert got["blocks"] == 1 and got["undrawable_blocks"] == 0
    assert got["windows_drawn"] > 0
    named = got["from_the_named_gene"]
    assert named["unknown_blocks_naming_a_coding_gene"] == 1
    assert named["unknown_block_rate"] == 1.0
    # no random window can hold the element inside the unknown block, so the control names nothing
    assert named["random_windows_naming_a_coding_gene"] == 0
    assert named["difference_in_points"] == 100.0


def test_against_vista_reports_the_head_beside_the_structural_zero(tmp_path, monkeypatch):
    import genomeos.attribution.unknown_scoring as us

    d = tmp_path / "data" / "knowledge" / "vista"
    d.mkdir(parents=True)
    (d / "rows_chrT.json").write_text(
        json.dumps(
            [
                {"id": "v1", "start": 0, "end": 100, "status": "positive"},
                {"id": "v2", "start": 200, "end": 300, "status": "negative"},
            ]
        )
    )
    monkeypatch.chdir(tmp_path)
    rows = [
        # covers the positive and clears the gate
        {"id": "a", "start": 10, "end": 20, "moves_gene": True, "abs_log2": 0.5, "head_abs_log2": 0.5},
        # covers the negative and does not: the old field had to invent a zero for it
        {"id": "b", "start": 210, "end": 220, "moves_gene": False, "abs_log2": None, "head_abs_log2": 0.08},
    ]
    got = us.against_vista("chrT", rows)
    assert got["mean_max_abs_log2_negative"] == 0.0
    assert got["negatives_scored_as_a_structural_zero"] == 1
    assert got["mean_max_head_abs_log2_negative"] == 0.08
    assert got["mean_max_abs_log2_positive"] == got["mean_max_head_abs_log2_positive"] == 0.5


def test_window_reading_leaves_a_line_no_gene_was_scored_on_unanswered():
    from genomeos.attribution.unknown_scoring import window_reading

    window = {
        ("chrT", "E1"): {
            "ONE": {
                "gene": "ONE",
                "max_drop_log2fc": -0.4,
                "max_rise_log2fc": 0.0,
                "by_cell": {"K562": -0.3, "HepG2": None, "GM12878": None, "IMR-90": None},
            }
        }
    }
    got = window_reading(_Responses(window), "chrT", "E1", set())
    assert got["window_by_cell"]["K562"] == -0.3
    # three lines are unanswered, and an unanswered line is None, never 0.0
    assert got["window_by_cell"]["HepG2"] is None
    assert got["cells_acting_window"] == 1
