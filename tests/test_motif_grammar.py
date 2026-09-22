"""Motif grammar: site collapse, the three nested feature sets, the shuffle control and the ridge."""

import random

from genomeos.attribution.motif_grammar import (
    GAP_BINS,
    PAIR_FEATURES,
    assemble,
    block_widths,
    bootstrap_models,
    collapse,
    composition,
    count_features,
    fit,
    fold_assignment,
    gap_bin,
    grammar_features,
    in_phase,
    interval,
    pair_index,
    predict,
    ranks,
    shuffle_sites,
    spearman,
    totals,
    verdict,
)
from genomeos.attribution.motif_grammar import accumulate as accumulate_rows


def test_collapse_keeps_one_best_site_per_family_and_position():
    hits = [
        ("KLF", 10, 20, 0, 0.96),
        ("KLF", 12, 22, 1, 0.99),  # overlaps the first on the other strand: the same site
        ("KLF", 40, 50, 0, 0.97),  # separate
        ("FOX", 14, 24, 0, 0.95),  # another family over the same bases: kept
    ]
    sites = collapse(hits)
    assert sites == [("KLF", 12, 22, 1), ("FOX", 14, 24, 0), ("KLF", 40, 50, 0)]


def test_composition_reads_gc_and_cpg():
    gc, gc2, oe = composition("GCGCGCGCGC")
    assert gc == 1.0 and gc2 == 1.0 and oe > 1.0
    assert composition("AAAATTTT")[0] == 0.0
    assert composition("NNNN") == [0.0, 0.0, 0.0]


def test_gap_bins_and_helical_phase():
    assert [gap_bin(g) for g in (0, 9, 10, 19, 20, 49, 50, 5000)] == [0, 0, 1, 1, 2, 2, 3, 3]
    assert gap_bin(-1) is None
    assert in_phase(0.0) and in_phase(10.5) and in_phase(21.0) and in_phase(9.0)
    assert not in_phase(5.0) and not in_phase(15.5)


def test_count_features_are_log1p_per_unit_plus_totals():
    units = ["KLF", "FOX", "GATA"]
    sites = [("KLF", 0, 8, 0), ("KLF", 20, 28, 1), ("FOX", 40, 48, 0), ("ETS", 60, 68, 0)]
    f = count_features(sites, units)
    assert round(f[0], 4) == round(__import__("math").log1p(2), 4)  # two KLF sites
    assert 1 in f and 2 not in f  # FOX present, GATA absent
    assert f[len(units)] > 0 and f[len(units) + 1] > 0  # total sites, distinct families
    assert count_features([], units) == {}


def test_grammar_features_count_pairs_by_gap_orientation_and_phase():
    pairs = pair_index(["A", "B"])
    assert pairs == {("A", "A"): 0, ("A", "B"): 1, ("B", "B"): 2}
    sites = [("A", 0, 10, 0), ("B", 15, 25, 0)]  # gap 5, same strand, centres 15 apart: out of phase
    f = grammar_features(sites, pairs)
    base = pairs[("A", "B")] * PAIR_FEATURES
    assert f[base + 0] > 0  # gap bin 0 (under 10 bp)
    assert f[base + len(GAP_BINS)] > 0  # same strand
    assert base + len(GAP_BINS) + 1 not in f  # not in helical phase
    turned = [("A", 0, 10, 0), ("B", 21, 31, 1)]  # centres 21 apart: two turns, opposite strand
    g = grammar_features(turned, pairs)
    assert g[base + len(GAP_BINS) + 1] > 0
    assert base + len(GAP_BINS) not in g
    overlapping = grammar_features([("A", 0, 10, 0), ("B", 5, 15, 0)], pairs)
    assert overlapping == {}  # sites that overlap cannot both be bound, and are not a pair
    assert grammar_features([("C", 0, 10, 0), ("D", 20, 30, 0)], pairs) == {}  # families outside the list


def test_shuffle_keeps_counts_and_positions_but_moves_the_grammar():
    rng = random.Random(1)
    sites = [("A", 0, 10, 0), ("B", 20, 30, 0), ("A", 40, 50, 1), ("B", 60, 70, 1)]
    units = ["A", "B"]
    pairs = pair_index(units)
    moved = 0
    for _ in range(40):
        s = shuffle_sites(sites, rng)
        assert sorted(x[0] for x in s) == sorted(x[0] for x in sites)
        assert sorted((x[1], x[2]) for x in s) == sorted((x[1], x[2]) for x in sites)
        assert count_features(s, units) == count_features(sites, units)
        if grammar_features(s, pairs) != grammar_features(sites, pairs):
            moved += 1
    assert moved > 0


def test_nested_feature_blocks_are_prefixes_of_each_other():
    units = ["A", "B"]
    pairs = pair_index(units)
    w = block_widths(units, pairs)
    assert w["a"] == 4 and w["b"] == 4 + len(units) + 2
    assert w["c"] == w["b"] + len(pairs) * PAIR_FEATURES
    row = assemble([0.5, 0.25, 0.9], [("A", 0, 10, 0), ("B", 20, 30, 0)], units, pairs)
    assert row[0] == 1.0 and row[1] == 0.5
    assert all(i < w["c"] for i in row)
    assert any(i >= w["b"] for i in row)  # the grammar block is populated


def test_ridge_recovers_a_planted_signal_and_prefers_the_right_columns():
    rng = random.Random(3)
    rows, y = [], []
    for _ in range(300):
        a, b = rng.random(), rng.random()
        rows.append({0: 1.0, 1: a, 2: b})
        y.append(3.0 * a - 2.0 * b + rng.gauss(0, 0.05))
    xtx, xty = accumulate_rows(rows, {"cell": y}, 3)
    blocks = {"chr1": (xtx, xty)}
    tx, ty = totals(blocks, ["chr1"], 3, ["cell"])
    w = fit(tx, ty["cell"], 3, 3, 1.0)
    assert w is not None and 2.5 < w[1] < 3.5 and -2.5 < w[2] < -1.5
    rho = spearman(predict(w, rows), y)
    assert rho is not None and rho > 0.95
    small = fit(tx, ty["cell"], 3, 2, 1.0)  # a nested prefix: the intercept and a only
    assert small is not None and len(small) == 2


def test_ranks_spearman_and_interval():
    assert ranks([10, 20, 20, 30]) == [1.0, 2.5, 2.5, 4.0]
    assert round(spearman([1, 2, 3], [1, 2, 3]), 9) == 1.0
    assert round(spearman([1, 2, 3], [3, 2, 1]), 9) == -1.0
    assert spearman([1, 2], [1, 2]) is None
    lo, hi = interval(list(range(101)))
    assert (lo, hi) == (2, 98)


def test_bootstrap_and_verdict_read_the_intervals():
    rng = random.Random(5)
    truth = [rng.gauss(0, 1) for _ in range(200)]
    preds = {
        "a": [rng.gauss(0, 1) for _ in truth],  # noise
        "b": [t + rng.gauss(0, 1.2) for t in truth],
        "c": [t + rng.gauss(0, 0.2) for t in truth],  # much better
        "shuffled": [t + rng.gauss(0, 1.2) for t in truth],
    }
    out = bootstrap_models(truth, preds, draws=120, seed=7)
    assert out["c"]["point"] > out["b"]["point"] > out["a"]["point"]
    assert out["c-b"]["ci95"][0] > 0 and out["c-shuffled"]["ci95"][0] > 0
    v = verdict({"K562": out, "HepG2": out})
    assert v["preregistered_claim_holds"] and v["arrangement_reading_holds"]
    flat = {
        "K562": {"c-b": {"ci95": [-0.01, 0.01]}, "c-shuffled": {"ci95": [-0.01, 0.01]}},
    }
    n = verdict(flat)
    assert not n["preregistered_claim_holds"]
    assert n["statement"].startswith("grammar features add nothing")


def test_fold_assignment_balances_chromosomes():
    sizes = {"chr1": 5000, "chr2": 4000, "chr3": 1000, "chr4": 900}
    folds = fold_assignment(sizes, folds=2)
    load = [0, 0]
    for c, f in folds.items():
        load[f] += sizes[c]
    assert set(folds.values()) == {0, 1} and abs(load[0] - load[1]) < 1000
