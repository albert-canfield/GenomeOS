# SPDX-License-Identifier: AGPL-3.0-or-later
"""Discovered element types: the behaviour vector, the clustering, the controls' building blocks."""

import math

import pytest

from genomeos.attribution import element_types as et


def _element(genes, nearest="A", distances=None):
    return {
        "cells": et.summarise_cells(genes),
        "nearest_gene": nearest,
        "nearest_distance": 1000,
        "target_distance": distances or {c: 5000 for c in et.CELLS},
    }


def test_summarise_cells_takes_the_strongest_signed_effect_per_line():
    genes = [
        {"gene": "A", "by_cell": {"K562": -0.3, "HepG2": 0.02, "GM12878": -0.01, "IMR-90": 0.0}},
        {"gene": "B", "by_cell": {"K562": 0.1, "HepG2": 0.25, "GM12878": 0.0, "IMR-90": -0.05}},
    ]
    per = et.summarise_cells(genes)
    assert per["K562"]["top"] == -0.3 and per["K562"]["gene"] == "A"
    assert per["HepG2"]["top"] == 0.25 and per["HepG2"]["gene"] == "B"
    assert per["K562"]["moving"] == ("A", "B")
    assert per["GM12878"]["moving"] == ()
    assert per["IMR-90"]["n"] == 2


def test_behaviour_vector_reads_strength_sign_breadth_and_disagreement():
    genes = [
        {"gene": "A", "by_cell": {"K562": -0.3, "HepG2": 0.2, "GM12878": -0.01, "IMR-90": 0.0}},
        {"gene": "B", "by_cell": {"K562": 0.0, "HepG2": 0.0, "GM12878": 0.0, "IMR-90": 0.0}},
    ]
    v = dict(zip(et.FEATURES, et.behaviour(_element(genes)), strict=True))
    assert v["log_strongest"] == pytest.approx(math.log10(0.3 + et.FLOOR))
    assert v["sign_strongest"] == -1.0
    assert v["lines_acting"] == 0.5
    assert v["direction_differs"] == 1.0  # K562 drops, HepG2 rises
    assert v["acts_on_nearest"] == 1.0
    assert v["share_window_moving"] == 0.5
    assert v["log_target_distance"] == pytest.approx(math.log10(5001))


def test_leaving_a_line_out_removes_its_effect_from_the_vector():
    genes = [{"gene": "A", "by_cell": {"K562": -0.5, "HepG2": 0.01, "GM12878": 0.01, "IMR-90": 0.01}}]
    el = _element(genes)
    full = dict(zip(et.FEATURES, et.behaviour(el), strict=True))
    loo = dict(zip(et.FEATURES, et.behaviour(el, ("HepG2", "GM12878", "IMR-90")), strict=True))
    assert full["lines_acting"] == 0.25 and loo["lines_acting"] == 0.0
    assert loo["log_strongest"] < full["log_strongest"]
    assert et.behaviour({"cells": {}}, et.CELLS) is None


def test_kmeans_separates_two_blobs_and_assign_transfers_them():
    pytest.importorskip("numpy")
    import random

    rng = random.Random(1)
    a = [[rng.gauss(0, 0.1), rng.gauss(0, 0.1)] for _ in range(200)]
    b = [[rng.gauss(5, 0.1), rng.gauss(5, 0.1)] for _ in range(200)]
    fit = et.kmeans(a + b, 2, seed=3)
    assert len(set(fit["labels"][:200])) == 1 and len(set(fit["labels"][200:])) == 1
    assert fit["labels"][0] != fit["labels"][-1]
    new = et.assign([[0.05, -0.02], [4.9, 5.1]], fit["centroids"])
    assert new == [fit["labels"][0], fit["labels"][-1]]
    s = et.silhouette(a + b, fit["labels"], sample=300)
    assert s > 0.9


def test_choose_k_prefers_the_real_number_of_groups():
    pytest.importorskip("numpy")
    import random

    rng = random.Random(2)
    centres = [(0, 0), (6, 0), (0, 6)]
    x = [[cx + rng.gauss(0, 0.3), cy + rng.gauss(0, 0.3)] for cx, cy in centres for _ in range(150)]
    choice = et.choose_k(x, ks=(2, 3, 4, 5))
    assert choice["k"] == 3
    assert set(choice["curve"]) == {2, 3, 4, 5}


def test_rand_index_agreement_and_cramers_v():
    assert et.adjusted_rand([0, 0, 1, 1], [1, 1, 0, 0]) == pytest.approx(1.0)
    assert et.matched_agreement([0, 0, 1, 1], [1, 1, 0, 0]) == 1.0
    assert et.adjusted_rand([0, 0, 1, 1, 2, 2], [0, 1, 0, 1, 0, 1]) < 0.1
    assert et.cramers_v([0, 0, 1, 1], ["a", "a", "b", "b"]) == pytest.approx(1.0)
    assert et.cramers_v([0, 1, 0, 1], ["a", "a", "b", "b"]) == pytest.approx(0.0)
    assert et.cramers_v([0, 0], ["a", "b"]) is None


def test_controls_keep_marginals_and_strata():
    x = [[i, 10 * i] for i in range(50)]
    shuffled = et.column_shuffle(x, seed=4)
    assert sorted(r[0] for r in shuffled) == list(range(50))
    assert sorted(r[1] for r in shuffled) == [10 * i for i in range(50)]
    assert any(r[1] != 10 * r[0] for r in shuffled)
    keys = ["a"] * 10 + ["b"] * 10
    donors = et.stratified_donors(keys, seed=5)
    assert sorted(donors) == list(range(20))
    assert all(keys[d] == keys[i] for i, d in enumerate(donors))


def test_scale_and_fit_predict_use_the_training_scale():
    means, sds = et.scaler([[1.0, 5.0], [3.0, 5.0]])
    assert means == [2.0, 5.0] and sds == [1.0, 1.0]
    assert et.scale([[3.0, 6.0]], means, sds) == [[1.0, 1.0]]
    train = [[float(i)] for i in range(20)]
    pred = et.fit_predict(train, [2.0 * i + 1 for i in range(20)], [[30.0]], lam=1e-6)
    assert pred[0] == pytest.approx(61.0, rel=1e-3)
    assert et.one_hot(2, 4) == [0.0, 1.0, 0.0]


def test_tasks_follow_the_epigenome_layers_definitions():
    units = [{"lfc": -0.2}, {"lfc": 0.15}, {"lfc": 0.01}]
    t = et.tasks_of(units)
    assert t["acts"] == ([1.0, 1.0, 0.0], [0, 1, 2])
    assert t["rise_among_acting"] == ([0.0, 1.0], [0, 1])
    assert t["magnitude"][0][0] == pytest.approx(math.log10(0.2 + et.FLOOR))
