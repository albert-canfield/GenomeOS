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
