# SPDX-License-Identifier: AGPL-3.0-or-later
"""The tip-noise pair null for area J: the model, the calibration, and the outcome rule."""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pytest

from genomeos.knowledge import pair_null as pn
from genomeos.knowledge import species_null as sn

# ------------------------------------------------------------------------------------ helpers


def balanced(depth: int) -> str:
    """A balanced binary tree with 2**depth tips t0, t1, ..."""
    names = iter(range(2**depth))

    def build(k: int) -> str:
        return f"t{next(names)}" if k == 0 else f"({build(k - 1)},{build(k - 1)})"

    return build(depth) + ";"


def grafen_tree(depth: int = 6) -> tuple[sn.TreeArrays, np.ndarray]:
    root = sn.parse_newick(balanced(depth))
    t = sn.arrays(root, [f"t{i}" for i in range(2**depth)])
    return t, sn.covariance(t, sn.grafen_lengths(t))


def draw(cov: np.ndarray, lam: float, sigma: float, rng: np.random.Generator) -> np.ndarray:
    v = lam * cov + (1 - lam) * np.eye(len(cov))
    return np.linalg.cholesky(v) @ rng.standard_normal(len(cov)) * sigma


def synthetic(n_sp: int = 12, n_genes: int = 600, seed: int = 8) -> sn.PairInputs:
    """Presence for three strata; the Mammalia stratum is absent from the first three species."""
    rng = np.random.default_rng(seed)
    bits = rng.random((n_genes, n_sp)) < np.linspace(0.2, 0.9, n_sp)
    strata = ("Life", "Metazoa", "Mammalia")
    for i in range(n_genes):
        if strata[i % 3] == "Mammalia":
            bits[i, :3] = False
    genes = [f"g{i:04d}" for i in range(n_genes)]
    width = (n_sp + 3) // 4
    hexes = {
        g: format(int("".join("1" if b else "0" for b in row), 2), f"0{width}x")
        for g, row in zip(genes, bits, strict=True)
    }
    presence = {"species": [f"s{i}" for i in range(n_sp)], "genes": hexes}
    origin = {"genes": {g: {"origin": strata[i % 3]} for i, g in enumerate(genes)}}
    return sn.pair_inputs(presence, origin)


TREE12 = "(((s0,s1,s2),(s3,(s4,s5))),((s6,s7),(s8,(s9,(s10,s11)))));"


# -------------------------------------------------------------------------------------- model


def test_at_lambda_one_the_rate_is_the_brownian_reml_rate() -> None:
    _, cov = grafen_tree(5)
    rng = np.random.default_rng(1)
    y = rng.normal(size=len(cov)) + 0.4
    evals, evecs = np.linalg.eigh(cov)
    _, s2, mu = pn.reml_profile(y, evals, evecs, np.array([1.0]))
    bm_mu, bm_rate = sn.Brownian(cov, sims=10, seed=2).rate(y)
    assert s2[0] == pytest.approx(bm_rate, rel=1e-9)
    assert mu[0] == pytest.approx(bm_mu, rel=1e-9)
    # at lambda = 0 the model is white noise: the plain variance about the plain mean
    _, s2_0, mu_0 = pn.reml_profile(y, evals, evecs, np.array([0.0]))
    assert mu_0[0] == pytest.approx(y.mean()) and s2_0[0] == pytest.approx(y.var(ddof=1))


def test_reml_tells_tip_noise_from_motion_along_the_tree() -> None:
    _, cov = grafen_tree(6)
    evals, evecs = np.linalg.eigh(cov)
    rng = np.random.default_rng(3)
    white = [pn.reml_lambda(draw(cov, 0.0, 0.5, rng), evals, evecs) for _ in range(30)]
    tree = [pn.reml_lambda(draw(cov, 0.9, 0.5, rng), evals, evecs) for _ in range(30)]
    assert np.mean([f.lam for f in white]) < 0.25
    assert np.mean([f.lam for f in tree]) > 0.6
    assert np.median([f.sigma2 for f in white]) == pytest.approx(0.25, rel=0.25)
    assert all(f.columns == 64 for f in white)


def test_a_constant_deviation_fits_no_rate() -> None:
    _, cov = grafen_tree(3)
    evals, evecs = np.linalg.eigh(cov)
    fit = pn.reml_lambda(np.full(8, 0.2), evals, evecs)
    assert fit.sigma2 == pytest.approx(0.0, abs=1e-20) and fit.mu == pytest.approx(0.2)


def test_the_banks_hold_the_tree_and_the_tip_noise() -> None:
    _, cov = grafen_tree(3)
    tn = pn.TipNoise(cov, sims=40_000, seed=4)
    assert np.allclose(np.cov(tn.b1, rowvar=False), cov, atol=0.03)
    assert np.allclose(np.cov(tn.w2, rowvar=False), np.eye(8), atol=0.03)
    assert abs(float(np.corrcoef(tn.b1[:, 0], tn.b2[:, 0])[0, 1])) < 0.03
    assert abs(float(np.corrcoef(tn.b1[:, 0], tn.w1[:, 0])[0, 1])) < 0.03


def test_the_null_pearson_is_the_row_by_row_pearson() -> None:
    _, cov = grafen_tree(3)
    tn = pn.TipNoise(cov, sims=50, seed=5)
    rng = np.random.default_rng(6)
    h = np.r_[0.0, 0.0, rng.uniform(0.2, 0.5, 6)]
    a = tn.side(rng.normal(size=8), rng.normal(size=8) * h, h)
    b = tn.side(rng.normal(size=8), rng.normal(size=8), np.ones(8))
    fast = tn.null_r(a, b, chunk=7)
    slow = [
        sn.pearson(
            a.fixed + a.p * tn.b1[i].astype(float) + a.q * tn.w1[i].astype(float),
            b.fixed + b.p * tn.b2[i].astype(float) + b.q * tn.w2[i].astype(float),
        )
        for i in range(50)
    ]
    assert np.allclose(fast, slow, atol=1e-6)


def test_a_species_with_no_room_to_move_draws_no_deviation() -> None:
    d = synthetic()
    mammalia = np.flatnonzero(d.origin_of_row == "Mammalia")[:20]
    h = pn.binomial_scale(d, mammalia)
    assert np.all(h[:3] == 0) and np.all(h[3:] > 0)
    m = d.stratum_mean["Mammalia"]
    assert np.allclose(h, np.sqrt(m * (1 - m)))
    life = np.flatnonzero(d.origin_of_row == "Life")[:10]
    mixed = np.r_[mammalia, life]
    want = np.sqrt((20 * m * (1 - m) + 10 * d.stratum_mean["Life"] * (1 - d.stratum_mean["Life"])) / 30)
    assert np.allclose(pn.binomial_scale(d, mixed), want)
    t = sn.arrays(sn.parse_newick(TREE12), d.species)
    tn = pn.TipNoise(sn.covariance(t, sn.grafen_lengths(t)), sims=20, seed=7)
    f = np.arange(12.0)
    dev = sn.profile(d, mammalia) - sn.age_expectation(d, mammalia)
    assert np.all(dev[:3] == 0)
    side = tn.side(f, dev, h)
    assert np.all(side.p[:3] == 0) and np.all(side.q[:3] == 0) and np.all(side.fixed[:3] == f[:3])
    assert side.fit.columns == 9
    assert np.all(pn.scale_of(d, mammalia, "lambda") == 1)
    with pytest.raises(ValueError):
        pn.scale_of(d, mammalia, "brownian")


def p_value(null: np.ndarray, r_obs: float) -> float:
    return (int((null >= r_obs).sum()) + 1) / (len(null) + 1)


def test_tip_noise_fixes_the_failure_section_22_recorded() -> None:
    """White-noise deviations around an alike age expectation: §22's Brownian null finds the random
    pair alike far too often, the tip-noise null about as often as it should."""
    _, cov = grafen_tree(6)
    n = len(cov)
    rng = np.random.default_rng(8)
    split = np.r_[np.ones(n // 2), -np.ones(n // 2)]
    fa = split + rng.normal(scale=0.3, size=n)
    fb = split + rng.normal(scale=0.3, size=n)
    bm = sn.Brownian(cov, sims=2000, seed=9)
    tn = pn.TipNoise(cov, sims=2000, seed=10)
    ones = np.ones(n)
    p_bm, p_tn = [], []
    for _ in range(200):
        da, db = rng.normal(scale=0.4, size=n), rng.normal(scale=0.4, size=n)
        r_obs = sn.pearson(fa + da, fb + db)
        va, vb = bm.rate(da)[1], bm.rate(db)[1]
        p_bm.append(p_value(bm.null_r(fa, fb, va**0.5, vb**0.5), r_obs))
        p_tn.append(p_value(tn.null_r(tn.side(fa, da, ones), tn.side(fb, db, ones)), r_obs))
    share_bm = np.mean(np.array(p_bm) <= 0.05)
    share_tn = np.mean(np.array(p_tn) <= 0.05)
    assert share_bm > 0.5
    assert share_tn <= pn.CAL_MAX_SHARE
    assert pn.summarise(p_tn)["passes"] and not pn.summarise(p_bm)["passes"]


# -------------------------------------------------------------------------------- calibration


def test_the_first_calibration_set_is_section_22_s_random_pairs() -> None:
    d = synthetic()
    t = sn.arrays(sn.parse_newick(TREE12), d.species)
    bm = sn.Brownian(sn.covariance(t, sn.grafen_lengths(t)), sims=300, seed=11)
    genes = d.genes
    sets = {"x": set(genes[0:40]), "y": set(genes[30:90]), "z": set(genes[200:260])}
    pairs = sn.testable_pairs(sets)
    old = sn.null_a(d, bm, sets, halves=False)["calibration"]
    rand = pn.random_pairs(d, pairs, pn.CAL1_SEED)
    assert pn.brownian_hits(d, bm, rand) == old["at_or_below_005"]
    for (_, _, oa, ob), (_, _, ra, rb) in zip(pairs, rand, strict=True):
        assert sorted(d.origin_of_row[ra].tolist()) == sorted(d.origin_of_row[sn.rows_of(d, oa)].tolist())
        assert sorted(d.origin_of_row[rb].tolist()) == sorted(d.origin_of_row[sn.rows_of(d, ob)].tolist())
    again = pn.random_pairs(d, pairs, pn.CAL2_SEED)
    assert any(not np.array_equal(a[2], b[2]) for a, b in zip(rand, again, strict=True))


def test_summarise_applies_both_ceilings() -> None:
    ok = pn.summarise([0.5] * 95 + [0.04] * 5)
    assert ok["passes"] and ok["at_or_below_005"] == 5 and ok["q_at_or_below_005"] == 0
    wide = pn.summarise([0.5] * 89 + [0.04] * 11)
    assert not wide["passes"] and wide["share_at_or_below_005"] == pytest.approx(0.11)
    # exactly 10% passes; a heavy tail fails even under 10%
    assert pn.summarise([0.5] * 90 + [0.04] * 10)["passes"]
    tail = pn.summarise([0.5] * 92 + [1e-5] * 2 + [0.04] * 6)
    assert tail["q_at_or_below_005"] == 2 and not tail["passes"]


def test_every_step_runs_on_a_small_world() -> None:
    d = synthetic()
    t = sn.arrays(sn.parse_newick(TREE12), d.species)
    tn = pn.TipNoise(sn.covariance(t, sn.grafen_lengths(t)), sims=400, seed=12)
    genes = d.genes
    sets = {"x": set(genes[0:40]), "y": set(genes[30:90]), "z": set(genes[200:260])}
    pairs = sn.testable_pairs(sets)
    for model in pn.MODELS:
        cal = pn.calibrate(d, tn, pn.random_pairs(d, pairs, 1), model)
        assert cal["pairs"] == 3 and len(cal["rows"]) == 3
        assert all(0 <= x <= 1 for r in cal["rows"] for x in r["lambda"])
    read = pn.read_pairs(d, tn, "lambda_binomial", pairs)
    assert read["tested"] == 3
    for row in read["pairs"]:
        assert 1 / 401 <= row["p_pair"] <= 1 and row["q_pair"] >= row["p_pair"]
        a, b = row["libraries"]
        want = sn.pearson(
            sn.profile(d, sn.rows_of(d, sets[a] - sets[b])) - d.genome,
            sn.profile(d, sn.rows_of(d, sets[b] - sets[a])) - d.genome,
        )
        assert row["r_exclusive"] == pytest.approx(want)
    hv = pn.halves(d, tn, "lambda", sets)
    assert hv["libraries"] == 3
    stored = {
        "pairs": [
            {"libraries": list(p[:2]), "r_exclusive": r["r_exclusive"]}
            for p, r in zip(pairs, read["pairs"], strict=True)
        ]
    }
    want_ok = {"n_tested": 3, "n_tested_14": 3, "n_mismatch": 0, "same_pairs": True}
    assert pn.reproduce(d, pairs, stored) == want_ok
    stored["pairs"][0]["r_exclusive"] += 0.01
    assert pn.reproduce(d, pairs, stored)["n_mismatch"] == 1


# -------------------------------------------------------------------------------------- outcomes


def placeholders(text: str) -> set[str]:
    return set(re.findall(r"{(\w+)}", text))


def stored_14(n: int = 4, called: int = 2) -> dict:
    pairs = [{"libraries": [f"l{i}", f"m{i}"], "r_exclusive": 0.1 * i} for i in range(n)]
    return {"pairs": pairs, "significant_q05": pairs[:called]}


def cal(hits: int, q: int = 0, n: int = 100) -> dict:
    return {
        "pairs": n,
        "at_or_below_005": hits,
        "share_at_or_below_005": hits / n,
        "q_at_or_below_005": q,
        "q_share_at_or_below_005": q / n,
        "passes": hits / n <= pn.CAL_MAX_SHARE and q / n <= pn.CAL_MAX_Q_SHARE,
    }


def read(qs: list[float]) -> dict:
    return {"pairs": [{"libraries": [f"l{i}", f"m{i}"], "p_pair": q, "q_pair": q} for i, q in enumerate(qs)]}


def halves_of(share: float) -> dict:
    return {"libraries": 10, "at_or_below_005": round(share * 10), "share_at_or_below_005": share}


REPRO = {"n_tested": 4, "n_tested_14": 4, "n_mismatch": 0, "same_pairs": True}


def test_the_rule_is_applied_in_its_registered_order() -> None:
    s = stored_14()
    assert pn.select({**REPRO, "n_mismatch": 1}, s)[0] == "not_reproduced"
    assert pn.select({**REPRO, "same_pairs": False}, s)[0] == "not_reproduced"
    both_fail = {"lambda_binomial": cal(20), "lambda": cal(12)}
    key, text, f = pn.select(REPRO, s, both_fail, 205, 205, None)
    assert key == "null_not_calibrated" and "{" not in text and f["selected"] is None
    # the scaled null first; the plain one only when the scaled one fails
    assert pn.choose_null({"lambda_binomial": cal(5), "lambda": cal(5)}) == "lambda_binomial"
    assert pn.choose_null({"lambda_binomial": cal(5, q=2), "lambda": cal(6)}) == "lambda"
    with pytest.raises(ValueError):  # the run must select what the rule selects
        pn.select(REPRO, s, {"lambda_binomial": cal(5), "lambda": cal(5)}, 205, 205, "lambda", cal(5))
    ok = {"lambda_binomial": cal(20), "lambda": cal(6)}
    key, text, f = pn.select(REPRO, s, ok, 205, 205, "lambda", cal(11))
    assert key == "not_confirmed" and "Pagel's lambda on the deviations as they stand" in text
    key, _, _ = pn.select(REPRO, s, ok, 205, 205, "lambda", cal(4), read([0.5] * 4), halves_of(0.8))
    assert key == "none_survive"
    key, _, _ = pn.select(REPRO, s, ok, 205, 205, "lambda", cal(4), read([0.5] * 4), halves_of(0.4))
    assert key == "none_survive_low_power"
    some = read([0.01, 0.5, 0.02, 0.5])
    key, text, f = pn.select(REPRO, s, ok, 205, 205, "lambda", cal(4), some, halves_of(0.8))
    assert key == "some_survive" and "{" not in text
    assert (f["n_surv"], f["k_14"], f["k_new"], f["lost_14"]) == (2, 1, 1, 1)
    assert f["survivors_among_14"] == [["l0", "m0"]]


def test_every_wording_states_only_figures_the_rule_fills() -> None:
    early = {"n_tested", "n_tested_14", "n_mismatch", "n_14"}
    cal1 = {
        "n_cal", "bm_c1_hits", "bm_c1_hits_22", "cal_bin_hits", "cal_bin_share", "cal_bin_q",
        "cal_lam_hits", "cal_lam_share", "cal_lam_q",
    }  # fmt: skip
    chosen = {"null_name", "cal1_hits", "cal1_q", "n_cal2", "cal2_hits", "cal2_share", "cal2_q"}
    late = {"n_surv", "k_14", "k_new", "lost_14", "halves_hits", "n_halves"}
    assert placeholders(pn.P_OUTCOMES["not_reproduced"]) <= early
    assert placeholders(pn.P_OUTCOMES["null_not_calibrated"]) <= early | cal1
    assert placeholders(pn.P_OUTCOMES["not_confirmed"]) <= early | cal1 | chosen
    for key in ("none_survive", "none_survive_low_power", "some_survive"):
        assert placeholders(pn.P_OUTCOMES[key]) <= early | cal1 | chosen | late
    assert placeholders(pn.CAL_SENTENCE) <= cal1
    assert set(pn.P_OUTCOMES) == {
        "not_reproduced", "null_not_calibrated", "not_confirmed", "none_survive",
        "none_survive_low_power", "some_survive",
    }  # fmt: skip


def test_no_wording_calls_a_pair_a_dependency() -> None:
    for text in pn.P_OUTCOMES.values():
        if "dependency" in text:
            assert re.search(r"no library pair is called a dependency", text, re.IGNORECASE), text
            assert "is a dependency" not in text.replace("is not a dependency", "")
    assert "not a dependency" in pn.P_OUTCOMES["some_survive"]


# --------------------------------------------------------------------------- the registration

REGISTRATION = Path("data/results/pair_null_registration.json")
ROOT = Path(__file__).resolve().parents[1]


def register_module():
    spec = importlib.util.spec_from_file_location(
        "pair_null_register", ROOT / "scripts" / "pair_null_register.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_registration_as_written_holds_no_figure_and_the_module_s_words() -> None:
    import tracked_paths as tp

    tp.must_be_present(REGISTRATION, was="the registration is not written yet")
    written = json.loads(REGISTRATION.read_text())
    assert written["result"] == "pair_null_registration"
    reg = register_module()
    hits = [
        k
        for k in reg.forbidden_at_any_depth(written)
        if not k.startswith(("result_manifest.", "null.outcome_wordings", "null.decision_rule"))
    ]
    assert hits == []
    # the run fills these words and no others; the module that fills them must hold the same
    assert written["null"]["outcome_wordings"] == pn.P_OUTCOMES
    assert written["null"]["calibration_sentence"] == pn.CAL_SENTENCE
    assert written["null"]["decision_rule"] == pn.P_RULE
    assert set(written["frozen_code"]["blobs"]) == set(reg.FROZEN)
    assert written["tree_as_built"]["tips"] == 355 and written["tree_as_built"]["grafen_diagonal_is_one"]
    assert written["null"]["parameters"]["MODELS"] == list(pn.MODELS)
    assert written["null"]["parameters"]["CAL1_SEED"] == sn.A_CAL_SEED
