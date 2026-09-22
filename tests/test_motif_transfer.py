"""Motif transfer: the rank statistics, the chunked scan, the nested VISTA models and the block readings.

Everything here is synthetic. No network, no reference genome, no cached scan.
"""

import math
import random

from genomeos.attribution.motif_transfer import (
    MODELS,
    auprc,
    auroc,
    block_coverage,
    blocks_verdict,
    bootstrap_auc,
    chunk_sequence,
    count_columns,
    count_row,
    covariates,
    fit_and_score,
    mann_whitney,
    match_by_length,
    measured_reading,
    merge_sites,
    model_columns,
    model_row,
    permutation_null,
    predicted_level,
    sweep_agreement,
    unit_counts,
    vista_feature_groups,
    vista_verdict,
    windows_of,
)


def test_auroc_is_the_rank_sum_and_handles_ties():
    assert auroc([1.0, 2.0, 3.0, 4.0], [0, 0, 1, 1]) == 1.0
    assert auroc([4.0, 3.0, 2.0, 1.0], [0, 0, 1, 1]) == 0.0
    assert auroc([1.0, 1.0, 1.0, 1.0], [0, 0, 1, 1]) == 0.5  # every score tied: chance
    assert auroc([1.0, 2.0, 2.0, 3.0], [0, 0, 1, 1]) == 0.875  # the tie splits between the classes
    assert auroc([1.0, 2.0], [1, 1]) is None  # one class only
    assert auroc([], []) is None


def test_auprc_is_average_precision():
    assert auprc([3.0, 2.0, 1.0], [1, 1, 0]) == 1.0
    assert auprc([3.0, 2.0, 1.0], [0, 1, 1]) == (0.5 + 2 / 3) / 2
    assert auprc([1.0, 1.0], [1, 0]) == 0.5  # tied: precision at the end of the group
    assert auprc([1.0, 2.0], [0, 0]) is None


def test_mann_whitney_reads_the_level_and_its_direction():
    high = mann_whitney([2.0, 3.0, 4.0, 5.0], [0.0, 0.5, 1.0, 1.5])
    assert high["auroc"] == 1.0 and high["z"] > 0 and high["n"] == [4, 4]
    low = mann_whitney([0.0, 0.5], [2.0, 3.0])
    assert low["auroc"] == 0.0 and low["z"] < 0
    same = mann_whitney([1.0, 2.0], [1.0, 2.0])
    assert same["auroc"] == 0.5 and same["z"] == 0.0


def test_bootstrap_reports_intervals_and_differences_between_models():
    rng = random.Random(7)
    labels = [i % 2 for i in range(200)]
    good = [y + rng.gauss(0, 0.3) for y in labels]
    useless = [rng.gauss(0, 1) for _ in labels]
    out = bootstrap_auc(labels, {"good": good, "useless": useless}, (("good", "useless"),), draws=120)
    assert out["good"]["auroc"] > 0.9 and out["good"]["ci95"][0] > 0.5
    assert out["good"]["auprc"] > 0.8
    assert 0.3 < out["useless"]["auroc"] < 0.7
    assert out["good - useless"]["ci95"][0] > 0
    assert out["good"]["above_chance"] == 1.0


def test_permutation_null_sits_at_chance():
    rng = random.Random(3)
    labels = [i % 2 for i in range(120)]
    scores = [y + rng.gauss(0, 0.2) for y in labels]
    null = permutation_null(scores, labels, draws=200)
    assert null["draws"] == 200
    assert 0.44 < null["mean"] < 0.56  # the permuted label carries no information
    assert null["ci95"][0] < 0.5 < null["ci95"][1]
    assert null["share_at_or_above_observed"] < 0.05


def test_chunking_covers_a_long_sequence_with_overlap():
    seq = "".join(random.Random(1).choice("ACGT") for _ in range(5000))
    pieces = chunk_sequence(seq, size=1000, overlap=50)
    assert pieces[0][0] == 0
    assert all(len(s) <= 1000 for _, s in pieces)
    assert pieces[-1][0] + len(pieces[-1][1]) == len(seq)  # the tail is covered
    for off, sub in pieces:
        assert seq[off : off + len(sub)] == sub
    assert chunk_sequence("ACGT", size=1000) == [(0, "ACGT")]


def test_merge_sites_removes_the_duplicates_two_chunks_produce():
    doubled = [("KLF", 100, 110, 0), ("KLF", 100, 110, 0), ("FOX", 104, 114, 1)]
    assert merge_sites(doubled) == [("KLF", 100, 110, 0), ("FOX", 104, 114, 1)]
    # two overlapping sites of one unit, seen in different chunks: one site
    assert merge_sites([("KLF", 100, 110, 0), ("KLF", 105, 115, 1)]) == [("KLF", 100, 110, 0)]
    assert unit_counts([("A", 0, 8, 0), ("A", 9, 17, 1), ("B", 20, 28, 0)]) == {"A": 2, "B": 1}


def test_count_row_is_section_17_b_with_composition_as_its_leading_block():
    units = ["KLF", "FOX"]
    row = count_row([0.5, 0.25, 0.9], {"KLF": 2, "ETS": 1}, units)
    assert count_columns(units) == 4 + len(units) + 2
    assert row[0] == 1.0 and row[1] == 0.5 and row[2] == 0.25 and row[3] == 0.9
    assert row[4] == math.pow(math.log1p(2), 1)  # two KLF sites
    assert 5 not in row  # no FOX site
    assert row[4 + len(units)] == math.log1p(3)  # three sites in all, ETS included
    assert row[5 + len(units)] == math.log1p(2)  # two distinct units
    empty = count_row([0.0, 0.0, 0.0], {}, units)
    assert empty == {0: 1.0}


def test_windows_of_tiles_and_drops_the_n_runs():
    seq = "ACGT" * 100 + "N" * 200 + "ACGT" * 50
    kept, dropped = windows_of(seq, size=200)
    assert dropped == 1 and len(kept) == 3
    assert all(len(w) == 200 and "N" not in w for w in kept)
    assert windows_of("ACGT" * 10, size=200) == ([], 0)  # shorter than one window


def test_length_matching_takes_each_control_once_and_keeps_the_distribution():
    target = [{"length": n} for n in (1000, 5000, 20000, 100000)]
    pool = [{"length": n} for n in (900, 1100, 4800, 21000, 98000, 500000)]
    got = match_by_length(target, pool)
    assert len(got) == len(target)
    assert len({id(g) for g in got}) == len(got)  # no control used twice
    assert sorted(g["length"] for g in got) == [1100, 4800, 21000, 98000]
    assert len(match_by_length(target, pool[:2])) == 2  # a short pool matches what it can


def test_vista_features_are_densities_and_fall_back_when_conservation_is_missing():
    units = ["KLF", "FOX"]
    element = {"length": 2000, "gc": [0.5, 0.25, 0.8], "counts": {"KLF": 4}, "counts_shuffled": {"FOX": 2}}
    fallback = {"phylop_mean": 1.0, "constrained_fraction": 0.3}
    g = vista_feature_groups({"phylop_mean": 2.5, "constrained_fraction": 0.7}, element, units, fallback)
    assert g["conservation"] == [2.5, 0.7]
    assert g["counts"][0] == math.log1p(4 * 1000 / 2000)  # two sites per kb, not four
    assert g["counts_shuffled"][0] == 0.0 and g["counts_shuffled"][1] > 0
    missing = vista_feature_groups({}, element, units, fallback)
    assert missing["conservation"] == [1.0, 0.3]
    assert model_columns(g, MODELS["length"]) == 2
    assert model_columns(g, MODELS["both"]) == 2 + 2 + 3 + len(units) + 2
    row = model_row(g, MODELS["conservation"])
    assert row[0] == 1.0 and row[2] == 2.5 and row[3] == 0.7


def test_nested_models_recover_a_planted_label_and_the_useless_block_adds_nothing():
    rng = random.Random(11)
    rows = {}
    labels = {}
    for chrom in ("chr1", "chr2", "chr3", "chr4", "chr5"):
        rows[chrom] = []
        labels[chrom] = []
        for _ in range(120):
            signal = rng.gauss(0, 1)
            noise = rng.gauss(0, 1)
            rows[chrom].append({0: 1.0, 1: signal, 2: noise})
            labels[chrom].append(1.0 if signal + rng.gauss(0, 0.4) > 0 else 0.0)
    train = ["chr1", "chr2", "chr3", "chr4"]
    chosen, preds = fit_and_score(rows, labels, train, ["chr5"], 3)
    truth = [int(v) for v in labels["chr5"]]
    assert chosen["inner_auroc"] > 0.8
    assert auroc(preds, truth) > 0.8
    noise_only = {c: [{0: 1.0, 1: r[2]} for r in rows[c]] for c in rows}
    _, weak = fit_and_score(noise_only, labels, train, ["chr5"], 2)
    assert abs(auroc(weak, truth) - 0.5) < 0.15


def test_covariates_state_what_is_unmatched_between_the_classes():
    meta = [
        {
            "status": "positive",
            "length": 2000,
            "gc": 0.5,
            "cpg_oe": 0.8,
            "phylop_mean": 2.0,
            "constrained_fraction": 0.5,
            "sites_per_kb": 40.0,
        },
        {
            "status": "positive",
            "length": 2200,
            "gc": 0.52,
            "cpg_oe": 0.9,
            "phylop_mean": 2.2,
            "constrained_fraction": 0.6,
            "sites_per_kb": 44.0,
        },
        {
            "status": "negative",
            "length": 1000,
            "gc": 0.4,
            "cpg_oe": 0.5,
            "phylop_mean": 1.0,
            "constrained_fraction": 0.2,
            "sites_per_kb": 20.0,
        },
        {
            "status": "negative",
            "length": 1100,
            "gc": 0.42,
            "cpg_oe": 0.6,
            "phylop_mean": 1.2,
            "constrained_fraction": 0.3,
            "sites_per_kb": 22.0,
        },
    ]
    cov = covariates(meta)
    assert cov["positive"]["elements"] == 2 and cov["negative"]["elements"] == 2
    assert cov["positive"]["median_length"] > cov["negative"]["median_length"]
    assert cov["positives_above_negatives_auroc"]["length"] == 1.0
    assert cov["positives_above_negatives_auroc"]["phylop_mean"] == 1.0
    assert "not a background matched" in cov["reading"]


def test_block_coverage_counts_the_missing_as_named_categories():
    tiers = {
        "real_unknown": [
            {"length": 10000, "mpra_keys": ["a"], "tested_elements": 2, "moving_elements": 1},
            {"length": 20000, "mpra_keys": [], "tested_elements": 1, "moving_elements": 0},
            {"length": 30000, "mpra_keys": [], "tested_elements": 0, "moving_elements": 0},
        ],
        "neutral": [{"length": 5000, "mpra_keys": ["b", "c"], "tested_elements": 1, "moving_elements": 1}],
    }
    cov = block_coverage(tiers)
    assert cov["real_unknown"]["blocks"] == 3
    assert cov["real_unknown"]["blocks_with_a_measured_element"] == 1
    assert cov["real_unknown"]["blocks_without_any_measured_element"] == 2
    assert cov["real_unknown"]["blocks_with_an_element_the_sweep_tested"] == 2
    assert cov["real_unknown"]["blocks_with_a_sweep_lead"] == 1
    assert cov["neutral"]["measured_elements"] == 2


def _synthetic_model(rng: random.Random):
    """A tiny lentiMPRA stand-in: two chromosomes, one block each.

    Column 1 stands for the counts and column 2 for composition; activity leans mostly on the counts, so
    the (b) weights, which read both, must beat the (a) weights, which read composition alone.
    """
    rows = {"chr1": [], "chr2": []}
    activity = {"chr1": {}, "chr2": {}}
    meta = {"chr1": [], "chr2": []}
    for c, lift in (("chr1", 1.0), ("chr2", 0.0)):
        vals = []
        for i in range(40):
            x, g = rng.random(), rng.random()
            rows[c].append({0: 1.0, 1: x, 2: g})
            vals.append(2.0 * x + 0.3 * g + lift + rng.gauss(0, 0.1))
            meta[c].append(
                {
                    "key": f"{c}:{i}",
                    "chrom": c,
                    "start": i * 1000,
                    "end": i * 1000 + 200,
                    "gc": 0.5,
                    "sites": 30,
                }
            )
        activity[c] = {"K562": vals, "HepG2": list(vals), "WTC11": list(vals)}
    return {
        "rows": rows,
        "activity": activity,
        "meta": meta,
        "weights": {
            "a": {cell: [0.0, 0.0, 0.3] for cell in ("K562", "HepG2", "WTC11")},
            "b": {cell: [0.0, 2.0, 0.3] for cell in ("K562", "HepG2", "WTC11")},
        },
        "elements_excluded_by_location": 80,
        "training_elements": 0,
    }


def test_measured_reading_compares_the_tiers_on_measured_activity():
    rng = random.Random(5)
    model = _synthetic_model(rng)
    tiers = {
        "real_unknown": [
            {
                "length": 40000,
                "mpra_keys": [m["key"] for m in model["meta"]["chr1"]],
                "tested_elements": 3,
                "moving_elements": 2,
            }
        ],
        "neutral": [
            {
                "length": 40000,
                "mpra_keys": [m["key"] for m in model["meta"]["chr2"]],
                "tested_elements": 3,
                "moving_elements": 1,
            }
        ],
    }
    out = measured_reading(model, tiers)
    assert out["coverage"]["real_unknown"]["measured_elements"] == 40
    assert out["measured_elements_with_activity_in_all_three_lines"]["neutral"] == 40
    # chr1 carries the +1.0 lift, so the real unknown reads ABOVE the neutral tier here
    assert out["level"]["K562"]["auroc"] > 0.8 and out["level"]["K562"]["z"] > 2
    # (b) reads column 1, (a) is an intercept only: the count model beats composition
    assert out["transfer"]["K562"]["b-a"]["ci95"][0] > 0
    assert out["transfer"]["K562"]["elements"] == 40


def test_predicted_level_and_sweep_agreement_read_their_own_columns():
    real = [
        {
            "length": 10000,
            "case": "syntax",
            "tested_elements": 2,
            "moving_elements": 1,
            "windows": 50,
            "mpra_keys": [],
            "mean_K562": 1.0,
            "max_K562": 2.0,
            "shuffled_mean_K562": 0.2,
        },
        {
            "length": 12000,
            "case": "syntax",
            "tested_elements": 2,
            "moving_elements": 0,
            "windows": 60,
            "mpra_keys": [],
            "mean_K562": 0.4,
            "max_K562": 1.0,
            "shuffled_mean_K562": 0.1,
        },
    ]
    control = [
        {
            "length": 10500,
            "tested_elements": 0,
            "moving_elements": 0,
            "windows": 52,
            "mpra_keys": [],
            "mean_K562": 0.3,
            "max_K562": 0.9,
            "shuffled_mean_K562": 0.1,
        },
        {
            "length": 11500,
            "tested_elements": 0,
            "moving_elements": 0,
            "windows": 57,
            "mpra_keys": [],
            "mean_K562": 0.2,
            "max_K562": 0.8,
            "shuffled_mean_K562": 0.05,
        },
    ]
    level = predicted_level(real, control, cells=("K562",))
    assert level["K562"]["real_unknown_vs_matched_neutral"]["auroc"] == 1.0
    assert level["K562"]["real_unknown_vs_its_shuffle"]["auroc"] == 1.0
    assert "not a measurement" in level["reading"]
    agree = sweep_agreement(real, cell="K562")
    assert agree["coverage"]["with_a_lead"] == 1
    assert "two model readings" in agree["reading"]
    assert agree["tested_blocks_only"]["count_model_mean"]["auroc"] == 1.0


def test_verdicts_name_the_negative_and_the_contradiction():
    fails = vista_verdict(
        {
            "held_out": {
                "counts": {"auroc": 0.52, "ci95": [0.47, 0.57]},
                "counts - composition": {"auroc": 0.01, "ci95": [-0.03, 0.05]},
                "both - conservation": {"auroc": 0.0, "ci95": [-0.04, 0.04]},
                "counts - counts_shuffled": {"auroc": 0.0, "ci95": [-0.05, 0.05]},
            },
            "tissue": {"by_group": {}},
        }
    )
    assert not fails["preregistered_claim_holds"]
    assert not fails["counts_separate_positives_from_negatives"]
    assert "does not transfer" in fails["statement"]
    holds = vista_verdict(
        {
            "held_out": {
                "counts": {"auroc": 0.7, "ci95": [0.65, 0.75]},
                "counts - composition": {"auroc": 0.08, "ci95": [0.04, 0.12]},
                "both - conservation": {"auroc": 0.05, "ci95": [0.01, 0.09]},
                "counts - counts_shuffled": {"auroc": 0.06, "ci95": [0.02, 0.1]},
            },
            "tissue": {"by_group": {"neural": {"auroc": 0.6, "ci95": [0.52, 0.68]}}},
        }
    )
    assert holds["preregistered_claim_holds"] and holds["a_tissue_group_is_predicted"]
    measured = {
        "coverage": {"real_unknown": {"measured_elements": 227}},
        "level": {"K562": {"auroc": 0.61, "z": 3.4}},
        "transfer": {"K562": {"b-a": {"point": 0.2, "ci95": [0.05, 0.35]}}},
    }
    level = {"K562": {"real_unknown_vs_matched_neutral": {"auroc": 0.55, "z": 2.5}}}
    agreement = {"tested_blocks_only": {"count_model_mean": {"auroc": 0.52, "ci95": [0.47, 0.58]}}}
    v = blocks_verdict(measured, level, agreement)
    assert v["measured_elements"] == 227
    assert v["count_model_transfers_inside_the_real_unknown"]
    assert "contradicts" in v["measured_statement"]
    below = blocks_verdict({**measured, "level": {"K562": {"auroc": 0.4, "z": -3.1}}}, level, agreement)
    assert "LESS active" in below["measured_statement"]
