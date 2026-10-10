# SPDX-License-Identifier: AGPL-3.0-or-later
"""The mixed calibration of area J's pair null: the draw, the reproduction gate and the outcome rule."""

from __future__ import annotations

import importlib.util
import json
import re
from pathlib import Path

import numpy as np
import pytest

from genomeos.knowledge import pair_mix as pm
from genomeos.knowledge import pair_null as pn
from genomeos.knowledge import species_null as sn

# ------------------------------------------------------------------------------------ helpers


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


def world(sims: int = 400) -> tuple[sn.PairInputs, pn.TipNoise, dict[str, set[str]], list[tuple]]:
    d = synthetic()
    t = sn.arrays(sn.parse_newick(TREE12), d.species)
    tn = pn.TipNoise(sn.covariance(t, sn.grafen_lengths(t)), sims=sims, seed=12)
    genes = d.genes
    sets = {"x": set(genes[0:40]), "y": set(genes[30:90]), "z": set(genes[200:260])}
    return d, tn, sets, sn.testable_pairs(sets)


def makeup(d: sn.PairInputs, rows: np.ndarray) -> list[str]:
    return sorted(s for s in d.origin_of_row[rows].tolist() if s)


def balanced(depth: int) -> str:
    names = iter(range(2**depth))

    def build(k: int) -> str:
        return f"t{next(names)}" if k == 0 else f"({build(k - 1)},{build(k - 1)})"

    return build(depth) + ";"


# ------------------------------------------------------------------------------------- the draw


def test_the_draw_matches_the_partner_s_make_up_and_shares_no_gene() -> None:
    d = synthetic()
    lib = np.arange(0, 40)
    partner = np.arange(40, 100)
    rng = np.random.default_rng(1)
    for _ in range(20):
        rows, short = pm.matched_rows_excluding(d, partner, lib, rng)
        assert short == 0
        assert makeup(d, rows) == makeup(d, partner)
        assert not set(rows.tolist()) & set(lib.tolist())


def test_a_pool_too_small_falls_short_and_says_so() -> None:
    d = synthetic()
    mammalia = d.pools["Mammalia"]
    lib = mammalia[2:]  # leaves two Mammalia genes in the pool
    asked = np.r_[mammalia[:5], d.pools["Life"][:5]]  # five Mammalia, five Life
    rows, short = pm.matched_rows_excluding(d, asked, lib, np.random.default_rng(2))
    assert short == 3
    assert makeup(d, rows).count("Mammalia") == 2 and makeup(d, rows).count("Life") == 5
    assert not set(rows.tolist()) & set(lib.tolist())


def test_mixed_pairs_keep_each_library_on_its_side_and_its_partner_s_age() -> None:
    d, _, _, pairs = world()
    mixed = pm.mixed_pairs(d, pairs, seed=3)
    assert len(mixed) == len(pairs) == 3
    for m, (a, b, oa, ob) in zip(mixed, pairs, strict=True):
        ra, rb = sn.rows_of(d, oa), sn.rows_of(d, ob)
        assert m["libraries"] == [a, b]
        s1, s2, short = m["first_library_real"]
        assert np.array_equal(s1, ra) and short == 0
        assert makeup(d, s2) == makeup(d, rb) and not set(s2.tolist()) & set(ra.tolist())
        assert np.allclose(sn.age_expectation(d, s2), sn.age_expectation(d, rb))
        s1, s2, short = m["second_library_real"]
        assert np.array_equal(s2, rb) and short == 0
        assert makeup(d, s1) == makeup(d, ra) and not set(s1.tolist()) & set(rb.tolist())
        assert np.allclose(sn.age_expectation(d, s1), sn.age_expectation(d, ra))


def test_the_draw_order_is_the_registered_one() -> None:
    d, _, _, pairs = world()
    mixed = pm.mixed_pairs(d, pairs, seed=pm.MIX_SEED)
    rng = np.random.default_rng(pm.MIX_SEED)
    for m, (_, _, oa, ob) in zip(mixed, pairs, strict=True):
        ra, rb = sn.rows_of(d, oa), sn.rows_of(d, ob)
        beside_a, _ = pm.matched_rows_excluding(d, rb, ra, rng)
        beside_b, _ = pm.matched_rows_excluding(d, ra, rb, rng)
        assert np.array_equal(m["first_library_real"][1], beside_a)
        assert np.array_equal(m["second_library_real"][0], beside_b)
    again = pm.mixed_pairs(d, pairs, seed=pm.MIX_SEED + 1)
    assert any(
        not np.array_equal(x["first_library_real"][1], y["first_library_real"][1])
        for x, y in zip(mixed, again, strict=True)
    )


# ------------------------------------------------------------------------------- the calibration


def test_every_step_runs_on_a_small_world() -> None:
    d, tn, _, pairs = world()
    read = pn.read_pairs(d, tn, pm.MODEL, pairs)
    ok = pm.reproduce_read(d, tn, pairs, read)
    assert ok["n_mismatch_read"] == 0 and ok["same_pairs"] and ok["same_survivors"]
    assert ok["n_tested"] == ok["n_stored"] == 3
    read["pairs"][0]["draws_at_or_above"] += 1
    assert pm.reproduce_read(d, tn, pairs, read)["n_mismatch_read"] == 1
    read["pairs"][0]["draws_at_or_above"] -= 1
    read["pairs"][1]["q_pair"] = 0.0 if read["pairs"][1]["q_pair"] > pn.Q_ALPHA else 1.0
    assert not pm.reproduce_read(d, tn, pairs, read)["same_survivors"]

    mix = pm.calibrate_mixed(d, tn, pm.mixed_pairs(d, pairs, seed=4))
    for k in pm.DIRECTIONS:
        assert mix[k]["pairs"] == 3 and len(mix["rows"][k]) == 3
        assert all(r["q_pair"] >= r["p_pair"] and 1 / 401 <= r["p_pair"] <= 1 for r in mix["rows"][k])
    assert mix["pooled"]["pairs"] == 6
    assert mix["lambda_library_sides"]["n"] == mix["lambda_random_sides"]["n"] == 6
    assert mix["random_genes_short"] == {"mixed_pairs_short": 0, "genes_short": 0}
    assert pm.passes(mix) == (mix["first_library_real"]["passes"] and mix["second_library_real"]["passes"])
    lib = pm.per_library(mix)
    assert sum(x["mixed_pairs"] for x in lib) == 6
    survivors = {tuple(read["pairs"][0]["libraries"])}
    assert pm.among_read(mix, survivors)["pairs"] == 1


def test_a_library_side_fits_what_it_fitted_in_the_read() -> None:
    """A side's fit is a function of its own deviation, so a library keeps §23's lambda in its mixed
    pairs, and the random side beside it is the only thing that changes."""
    d, tn, _, pairs = world()
    read = pn.read_pairs(d, tn, pm.MODEL, pairs)
    mix = pm.calibrate_mixed(d, tn, pm.mixed_pairs(d, pairs, seed=5))
    for lib_row, first, second in zip(
        read["pairs"], mix["rows"]["first_library_real"], mix["rows"]["second_library_real"], strict=True
    ):
        assert first["lambda"][0] == lib_row["lambda"][0]
        assert second["lambda"][1] == lib_row["lambda"][1]


def test_the_mixed_check_has_teeth() -> None:
    """A library deviation that lines up with its partner's age expectation reads as alike to random
    sets far too often; one of the same size that does not line up does so less often."""
    root = sn.parse_newick(balanced(6))
    t = sn.arrays(root, [f"t{i}" for i in range(64)])
    cov = sn.covariance(t, sn.grafen_lengths(t))
    n = len(cov)
    tn = pn.TipNoise(cov, sims=2000, seed=13)
    rng = np.random.default_rng(14)
    split = np.r_[np.ones(n // 2), -np.ones(n // 2)]
    fa = split + rng.normal(scale=0.3, size=n)
    fb = -split + rng.normal(scale=0.3, size=n)
    ones = np.ones(n)
    centred = (fb - fb.mean()) / fb.std()
    aligned_lib = 1.5 * centred
    plain_lib = 1.5 * rng.permutation(centred)
    p_aligned, p_plain = [], []
    for _ in range(100):
        db = rng.normal(scale=0.4, size=n)
        for lib, ps in ((aligned_lib, p_aligned), (plain_lib, p_plain)):
            r_obs = sn.pearson(fa + lib, fb + db)
            null = tn.null_r(tn.side(fa, lib, ones), tn.side(fb, db, ones))
            ps.append((int((null >= r_obs).sum()) + 1) / (len(null) + 1))
    assert not pn.summarise(p_aligned)["passes"]
    assert np.mean(np.array(p_plain) <= 0.05) < np.mean(np.array(p_aligned) <= 0.05)


# -------------------------------------------------------------------------------------- outcomes


def placeholders(text: str) -> set[str]:
    return set(re.findall(r"{(\w+)}", text))


def stored_14(n: int = 4, called: int = 2) -> dict:
    pairs = [{"libraries": [f"l{i}", f"m{i}"], "r_exclusive": 0.1 * i} for i in range(n)]
    return {"pairs": pairs, "significant_q05": pairs[:called]}


FIG23 = {"n_surv": 154, "k_14": 6, "lost_14": 29}
R14 = {"n_tested": 860, "n_tested_14": 860, "n_mismatch": 0, "same_pairs": True}
R23 = {
    "n_tested": 860,
    "n_stored": 860,
    "n_mismatch_read": 0,
    "same_pairs": True,
    "n_survivors_recomputed": 154,
    "n_survivors_stored": 154,
    "same_survivors": True,
}


def cal(hits: int, q: int = 0, n: int = 860) -> dict:
    return {
        "pairs": n,
        "at_or_below_005": hits,
        "share_at_or_below_005": hits / n,
        "q_at_or_below_005": q,
        "q_share_at_or_below_005": q / n,
        "passes": hits / n <= pn.CAL_MAX_SHARE and q / n <= pn.CAL_MAX_Q_SHARE,
    }


def mix_of(a: dict, b: dict) -> dict:
    lam = {"n": 1720, "min": 0.1, "q1": 0.6, "median": 0.7, "q3": 0.8, "max": 1.0}
    hits = a["at_or_below_005"] + b["at_or_below_005"]
    pooled = cal(hits, a["q_at_or_below_005"] + b["q_at_or_below_005"], 1720)
    return {
        "first_library_real": a,
        "second_library_real": b,
        "pooled": pooled,
        "lambda_library_sides": lam,
        "lambda_random_sides": lam,
    }


def test_the_rule_is_applied_in_its_registered_order() -> None:
    s = stored_14(35 + 10, 35)
    for bad14, bad23 in (
        ({"n_mismatch": 1}, {}),
        ({"same_pairs": False}, {}),
        ({}, {"n_mismatch_read": 2}),
        ({}, {"same_pairs": False}),
        ({}, {"same_survivors": False, "n_survivors_recomputed": 153}),
    ):
        key, text, _ = pm.select_mix({**R14, **bad14}, {**R23, **bad23}, s, FIG23)
        assert key == "not_reproduced" and "{" not in text
        with pytest.raises(ValueError):  # no mixed pair may be computed once rule 1 fails
            pm.select_mix({**R14, **bad14}, {**R23, **bad23}, s, FIG23, mix_of(cal(40), cal(40)))
    with pytest.raises(ValueError):  # the read reproduced, so the mixed pairs must be there
        pm.select_mix(R14, R23, s, FIG23)
    key, text, f = pm.select_mix(R14, R23, s, FIG23, mix_of(cal(40), cal(50)))
    assert key == "calibrated" and "{" not in text and "stands as read" in text
    assert "154 of 860" in text and "6 of §14's 35" in text and "decides nothing" in text
    # exactly 10% and exactly 1% pass; one direction over either ceiling fails the whole
    assert pm.select_mix(R14, R23, s, FIG23, mix_of(cal(86, 8), cal(10)))[0] == "calibrated"
    assert pm.select_mix(R14, R23, s, FIG23, mix_of(cal(87), cal(10)))[0] == "not_calibrated"
    assert pm.select_mix(R14, R23, s, FIG23, mix_of(cal(10), cal(20, 9)))[0] == "not_calibrated"
    key, text, f = pm.select_mix(R14, R23, s, FIG23, mix_of(cal(10), cal(120)))
    assert key == "not_calibrated" and "withdrawn to no reading" in text and "{" not in text
    assert "29 of §14's 35" in text and f["b_share"] == "14.0%"


def test_every_wording_states_only_figures_the_rule_fills() -> None:
    early = {
        "n_tested", "n_mismatch_14", "same_pairs_14", "n_mismatch_read", "same_pairs_23",
        "n_23_recomputed", "same_survivors", "n_23", "k_23_14", "lost_23_14", "n_14",
    }  # fmt: skip
    late = {
        "n_a", "a_hits", "a_share", "a_q", "n_b", "b_hits", "b_share", "b_q", "lam_lib", "lam_rand",
        "n_pool", "pool_hits", "pool_share", "pool_q",
    }  # fmt: skip
    assert placeholders(pm.MIX_OUTCOMES["not_reproduced"]) <= early
    for key in ("calibrated", "not_calibrated"):
        assert placeholders(pm.MIX_OUTCOMES[key]) <= early | late
    assert placeholders(pm.POOLED_SENTENCE) <= late
    assert set(pm.MIX_OUTCOMES) == {"not_reproduced", "calibrated", "not_calibrated"}


def test_no_wording_calls_a_pair_a_dependency() -> None:
    for text in pm.MIX_OUTCOMES.values():
        if "dependency" in text:
            assert re.search(r"no library pair is called a dependency", text, re.IGNORECASE), text
            assert "is a dependency" not in text.replace("is not a dependency", "")
    assert "not a dependency" in pm.MIX_OUTCOMES["calibrated"]


def test_the_null_calibrated_is_the_one_section_23_read() -> None:
    assert pn.MODELS[0] == pm.MODEL
    assert pm.MIX_SEED not in {
        pn.P_SEED, pn.CAL1_SEED, pn.CAL2_SEED, pn.HALVES_SEED, sn.A_SEED, sn.B_SEED, sn.B_SIM_SEED,
    }  # fmt: skip


# --------------------------------------------------------------------------- the registration

ROOT = Path(__file__).resolve().parents[1]
REGISTRATION = Path("data/results/pair_mix_registration.json")


def register_module():
    spec = importlib.util.spec_from_file_location(
        "pair_mix_register", ROOT / "scripts" / "pair_mix_register.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_registration_as_written_holds_no_figure_and_the_module_s_words() -> None:
    import tracked_paths as tp

    tp.must_be_present(REGISTRATION, was="the registration is not written yet")
    written = json.loads(REGISTRATION.read_text())
    assert written["result"] == "pair_mix_registration"
    reg = register_module()
    hits = [
        k
        for k in reg.forbidden_at_any_depth(written)
        if not k.startswith(
            ("result_manifest.", "mixed_calibration.outcome_wordings", "mixed_calibration.decision_rule")
        )
    ]
    assert hits == []
    # the run fills these words and no others; the module that fills them must hold the same
    mix = written["mixed_calibration"]
    assert mix["outcome_wordings"] == pm.MIX_OUTCOMES
    assert mix["pooled_sentence"] == pm.POOLED_SENTENCE
    assert mix["decision_rule"] == pm.MIX_RULE
    assert set(mix["what_happens_to_the_154"]) == set(pm.MIX_OUTCOMES)
    assert mix["parameters"]["MODEL"] == pm.MODEL and mix["parameters"]["MIX_SEED"] == pm.MIX_SEED
    assert set(written["frozen_code"]["blobs"]) == set(reg.FROZEN)
    assert written["null_parameters"]["equal_to_section_23_registration"]
    assert written["tree_as_built"]["tips"] == 355 and written["tree_as_built"]["grafen_diagonal_is_one"]
    inputs = {Path(e["path"]).name for e in written["result_manifest"]["inputs"]}
    assert {"pair_null_result.json", "pair_null_registration.json", "profiling_genome_wide.json"} <= inputs


RESULT = Path("data/results/pair_mix_result.json")


def test_the_committed_reading_is_the_rule_applied_to_its_own_figures() -> None:
    import hashlib

    import tracked_paths as tp

    tp.must_be_present(
        RESULT,
        REGISTRATION,
        Path("data/results/profiling_genome_wide.json"),
        Path("data/results/pair_null_result.json"),
    )
    out = json.loads(RESULT.read_text())
    assert out["registration"]["registration_sha256"] == hashlib.sha256(REGISTRATION.read_bytes()).hexdigest()
    stored_14 = json.loads(Path("data/results/profiling_genome_wide.json").read_text())["null"]
    s23 = json.loads(Path("data/results/pair_null_result.json").read_text())
    mix = out.get("mixed_calibration")
    key, text, _ = pm.select_mix(
        out["reproduction"]["section_14"], out["reproduction"]["section_23"], stored_14, s23["figures"], mix
    )
    assert (key, text) == (out["outcome"], out["reading"])
    registered = json.loads(REGISTRATION.read_text())["mixed_calibration"]["what_happens_to_the_154"]
    assert out["section_23_consequence_registered"] == registered[key]
    if mix is not None:
        # each direction's figures are the registered summary of its own rows
        for k in pm.DIRECTIONS:
            ps = [r["p_pair"] for r in mix["rows"][k] if r["p_pair"] is not None]
            assert pn.summarise(ps) == mix[k]
            assert len(ps) == out["reproduction"]["section_14"]["n_tested"]
