# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, stage 2: the worm's constraint tiers (genomeos/attribution/simpler_tiers.py).

The cut is checked against Benjamini-Hochberg computed directly on sorted p values; the matched
difference is checked on synthetic elements where the answer is known, including one where the raw
difference is all confounding and the matched one must be near zero. The population is recomputed
from the git-ignored local copies when they are present (only that test skips without them).
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
import tracked_paths as tp

from genomeos.attribution import simpler
from genomeos.attribution import simpler_tiers as st

ROOT = Path(__file__).resolve().parents[1]
LOCAL = (ROOT / simpler.TABLE_S2, ROOT / simpler.GFF3)


def _bh_direct(p: np.ndarray, q: float = 0.05) -> float:
    """Benjamini-Hochberg on sorted p values: the score of the largest rejected p."""
    ps = np.sort(p)
    m = len(ps)
    ok = ps <= q * np.arange(1, m + 1) / m
    k = int(np.nonzero(ok)[0].max())
    return float(-np.log10(ps[k]))


def test_bh_cut_matches_benjamini_hochberg_on_sorted_p_values() -> None:
    rng = np.random.default_rng(1)
    p = np.concatenate([rng.random(95_000), rng.random(5_000) * 1e-5])
    s = (-np.log10(p)).astype(np.float32)
    hist = np.bincount(st.abs_bins(s), minlength=st.N_BINS)
    cut = st.bh_cut(hist)
    assert cut is not None
    assert abs(cut - _bh_direct(p)) <= 0.05
    # the rejected count at the cut keeps the FDR bound: m * 10^-cut / R <= q
    rejected = int((s >= np.float32(cut)).sum())
    assert len(p) * 10 ** (-cut) / rejected <= 0.05 + 1e-9


def test_bh_cut_none_without_signal_and_batch_agrees() -> None:
    rng = np.random.default_rng(2)
    null = (-np.log10(rng.random(50_000))).astype(np.float32)
    h0 = np.bincount(st.abs_bins(null), minlength=st.N_BINS)
    p = np.concatenate([rng.random(45_000), rng.random(5_000) * 1e-6])
    h1 = np.bincount(st.abs_bins((-np.log10(p)).astype(np.float32)), minlength=st.N_BINS)
    batch = st.bh_cut_batch(np.stack([h0, h1]))
    single = [st.bh_cut(h0), st.bh_cut(h1)]
    for b, one in zip(batch, single, strict=True):
        assert (math.isnan(b) and one is None) or b == one
    assert single[1] is not None


def test_float32_scores_land_in_the_bin_their_comparison_uses() -> None:
    v = np.array([2.27, 2.269, -2.27, 0.0, 25.0], dtype=np.float32)
    assert st.abs_bins(v).tolist() == [227, 226, 227, 0, st.N_BINS - 1]
    assert (v >= np.float32(2.27)).tolist() == [True, False, False, False, True]


def test_intervals_merge_subtract_inside() -> None:
    s, e = st.merge([(10, 20), (15, 30), (40, 50), (50, 55), (60, 60)])
    assert s.tolist() == [10, 40] and e.tolist() == [30, 55]
    assert st.subtract(0, 100, s, e) == [(0, 10), (30, 40), (55, 100)]
    assert st.subtract(12, 28, s, e) == []
    assert st.subtract(25, 45, s, e) == [(30, 40)]
    assert st.inside(np.array([9, 10, 29, 30, 54, 55]), s, e).tolist() == [
        False,
        True,
        True,
        False,
        True,
        False,
    ]


def test_tier_is_the_human_regulatory_rule() -> None:
    assert st.CONSTRAINED_MIN == 0.05
    assert st.tier(100, 5) == "constrained"
    assert st.tier(100, 4) == "not_constrained"
    assert st.tier(0, 0) == "not_measured"


def _rows(n: int, seed: int, effect: float, confound: bool) -> list[dict]:
    """Synthetic elements: tissue sets the base concordance; the constrained tier adds `effect`.
    With `confound`, constrained elements are mostly in the high-concordance tissue."""
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n):
        tissue = simpler.TISSUES[i % 2]  # Germline high, Neurons low
        base = 0.7 if tissue == "Germline" else 0.3
        p_con = (0.8 if tissue == "Germline" else 0.2) if confound else 0.5
        constrained = rng.random() < p_con
        y = int(rng.random() < base + (effect if constrained else 0.0))
        rows.append(
            {
                "chrom": "chrI",
                "mid": i * 20_000,
                "tissue": tissue,
                "domain_case": "host",
                "domain_dist": 500,
                "length_measured": 300,
                "tier": "constrained" if constrained else "not_constrained",
                "c_domain": y,
            }
        )
    return rows


def test_matched_difference_finds_a_true_effect() -> None:
    res = st.matched_difference(_rows(4_000, 3, 0.15, False), [200.0, 400.0], reps=300, perms=200)
    assert 10 < res["delta_points"] < 20
    assert res["delta_ci95_points"][0] > 0
    assert res["permutation"]["p_one_sided"] < 0.01


def test_matched_difference_removes_a_confounded_raw_difference() -> None:
    rows = _rows(4_000, 4, 0.0, True)
    res = st.matched_difference(rows, [200.0, 400.0], reps=300, perms=200)
    u = res["unmatched"]
    raw = 100 * (
        u["constrained_concordant"] / u["constrained"]
        - u["not_constrained_concordant"] / u["not_constrained"]
    )
    assert raw > 15  # the tissue makes the raw difference large
    assert abs(res["delta_points"]) < 5  # matched on tissue it is not there
    assert res["delta_ci95_points"][0] < 0 < res["delta_ci95_points"][1]


def test_permuted_labels_keep_each_stratum_count() -> None:
    rng = np.random.default_rng(5)
    s = rng.integers(0, 7, 1_000)
    t = rng.random(1_000) < 0.3
    n_per = np.bincount(s[t], minlength=7)
    tp = st._permuted_labels(s, n_per, rng)
    assert np.bincount(s[tp], minlength=7).tolist() == n_per.tolist()
    assert (tp != t).any()


def test_outcome_labels() -> None:
    assert st.outcome_label(0.1, 2.0) == "CONSTRAINT_ADDS"
    assert st.outcome_label(-2.0, -0.1) == "CONSTRAINT_BELOW"
    assert st.outcome_label(-1.0, 1.0) == "NO_DIFFERENCE_DETECTED"


def test_reading_fills_every_placeholder_and_keeps_the_fixed_sentences() -> None:
    primary = {
        "n": 5000,
        "constrained_in_support": 1200,
        "not_constrained_in_support": 3500,
        "constrained_concordant_in_support": 700,
        "constrained_rate_pct": 58.333,
        "matched_not_constrained_rate_pct": 54.1,
        "delta_points": 4.233,
        "delta_ci95_points": [1.2, 7.0],
        "permutation": {"p_one_sided": 0.002},
    }
    res = {
        "outcome": "CONSTRAINT_ADDS",
        "falsifier": "NOT_FALSIFIED",
        "primary": primary,
        "strata_beside": {"between_genes": {"delta_ci95_points": [-1.0, 5.0]}},
    }
    text = st.reading(res)
    assert "{" not in text and "}" not in text
    assert "700 of 1,200 constrained elements (58.3%)" in text
    assert text.endswith(st.REGISTRATION["outcome_wording"]["always"])
    assert st.REGISTRATION["outcome_wording"]["between_genes_caveat"] in text
    res["strata_beside"]["between_genes"]["delta_ci95_points"] = [0.5, 5.0]
    assert st.REGISTRATION["outcome_wording"]["between_genes_caveat"] not in st.reading(res)
    assert st.reading({"outcome": "UNDERPOWERED"}) == st.REGISTRATION["outcome_wording"]["UNDERPOWERED"]


def test_sampling_draws_a_tenth_of_the_groups_with_the_seed() -> None:
    from genomeos.attribution.bigwig import LeafItem

    leaves = [LeafItem(0, i * 1000, 0, i * 1000 + 1000, i * 100, 100) for i in range(1_000)]
    a = st.sample_groups(leaves)
    b = st.sample_groups(leaves)
    assert [g[0].offset for g in a] == [g[0].offset for g in b]
    assert len(a) == math.ceil(st.SAMPLE_FRACTION * math.ceil(1_000 / st.SAMPLE_GROUP))
    assert all(len(g) == st.SAMPLE_GROUP for g in a[:-1])


@pytest.mark.skipif(not all(p.exists() for p in LOCAL), reason="git-ignored Table S2 or WBcel235 GFF3 absent")
def test_population_is_stage_one_population() -> None:
    chroms = simpler.read_genes(ROOT / simpler.GFF3)
    rows = st.population(
        chroms, simpler.read_elements(ROOT / simpler.TABLE_S2, kinds=(simpler.ELEMENT_TYPE,))
    )
    reg = json.loads((ROOT / "data/results/simpler_worm_registration.json").read_text())
    assert len(rows) == reg["census"]["elements"]
    assert sum(r["domain_case"] == "host" for r in rows) == reg["census"]["host_elements"]


TIERS = ROOT / "data/results/simpler_worm_tiers.json"


def test_stored_cut_is_benjamini_hochberg_on_the_stored_histogram() -> None:
    tp.must_be_present(TIERS, was="tiers not yet written")
    d = json.loads(TIERS.read_text())
    assert d["rna_sheet_read"] is False and d["bytes_stored"] == 0
    assert d["element_columns"] == list(st.ELEMENT_COLUMNS)
    for t in d["tracks"].values():
        hists = t["abs_hist_by_class"]
        total = sum(np.array(h["pos"]) + np.array(h["neg"]) for h in hists.values())
        assert st.bh_cut(total) == t["cut"]
        i = round(t["cut"] / st.HIST_BIN)
        above = sum(sum(h["pos"][i:]) for h in hists.values())
        assert round(above / total.sum(), 5) == t["constrained_share_all_scored"]
        assert t["cut_interval_95"][0] <= t["cut"] <= t["cut_interval_95"][1]
        assert t["bytes_moved_total"] <= 200_000_000
        cds, intron = (t["by_class"][c]["constrained_share"] for c in ("cds", "intron"))
        assert t["instrument_check"]["passed"] is (cds > intron)
    for row in d["elements"]:
        assert 0 <= row[4] <= row[3] and 0 <= row[6] <= row[5]


def test_fraction_tertiles_split_top_from_bottom_and_drop_the_middle() -> None:
    rows = [{"measured": 100, "above": a, "tier": "x"} for a in range(0, 30)] + [{"measured": 0, "above": 0}]
    fedges = st.fraction_edges(rows)
    assert len(fedges) == 2 and fedges[0] < fedges[1]
    out = st.by_fraction_tertile(rows, fedges)
    top = [r["above"] for r in out if r["tier"] == "constrained"]
    bottom = [r["above"] for r in out if r["tier"] == "not_constrained"]
    assert min(top) > max(bottom)
    assert len(top) + len(bottom) < 30 and abs(len(top) - len(bottom)) <= 2


REGISTRATION = ROOT / "data/results/simpler_worm_tiers_registration.json"


def test_committed_registration_is_the_module_registration() -> None:
    tp.must_be_present(REGISTRATION, was="registration not yet written")
    reg = json.loads(REGISTRATION.read_text())
    assert reg["written_before_tiers_were_read_against_rna"] is True
    assert reg["registration"] == json.loads(json.dumps(st.REGISTRATION))
    cen = reg["census"]
    assert cen["constrained"] + cen["not_constrained"] + cen["not_measured"] == cen["elements"]
    assert sum(sum(v.values()) for v in cen["by_tissue"].values()) == cen["elements"] - cen["not_measured"]
    assert cen["stop_rule_met"] is (
        min(cen["constrained_in_support"], cen["not_constrained_in_support"]) < 200
    )
    tiers = json.loads(TIERS.read_text())
    assert reg["registered_inputs"]["tiers"]["sha256"] == simpler.sha256(TIERS)
    for name, cut in reg["cuts"].items():
        assert cut["cut"] == tiers["tracks"][name]["cut"]
