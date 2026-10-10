# SPDX-License-Identifier: AGPL-3.0-or-later
"""M2h (lane-magcal): one scale per form, fitted on training genes and read on held-out ones.

The arithmetic of the scale and of Delta MSE, the split, the four-way decision table of the
supervisor's ruling, the control's refit, and that the registration says in advance everything the
brief required it to say: prior exposure, one number per form, the falsifier, the inconclusive case,
and the two vacuity readings (a scale interval including 0, and a negative scale).
"""

from __future__ import annotations

import hashlib
import importlib.util
import math
from pathlib import Path

import pytest

from genomeos.attribution import combine2 as c2
from genomeos.attribution import magcal as mc
from genomeos.attribution import multrna as mr

ROOT = Path(__file__).resolve().parents[1]


def _script(name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / f"{name}.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def rows(*pd: tuple[str, list[float], list[float]]) -> list[dict]:
    return [
        {"gene": s, "symbol": s, "chrom": "chr1", "n_rules": 2, "n_inhibiting": 0, "P": p, "D": d}
        for s, p, d in pd
    ]


# ---- the scale and the metric --------------------------------------------------------------------
def test_delta_mse_at_scale_1_is_combine2s_delta_mse():
    rs = rows(("A", [1.0, -1.0], [0.5, -2.0]), ("B", [2.0, 0.0], [1.0, 3.0]))
    assert mc.delta_mse(rs, 1.0) == pytest.approx(c2.dmse(rs))


def test_delta_mse_at_scale_0_is_exactly_the_baseline():
    rs = rows(("A", [3.0, -4.0], [0.5, -2.0]))
    assert mc.delta_mse(rs, 0.0) == 0.0


def test_fitted_scale_is_the_closed_form_and_maximises_training_delta_mse():
    rs = rows(("A", [1.0, -1.0, 2.0], [0.4, -0.2, 0.9]), ("B", [0.5, 1.5, -1.0], [0.1, 0.8, -0.3]))
    a = mc.fit_scale(rs)
    dp, pp, _dd, _n = mc._sums(rs)
    assert a == pytest.approx(dp / pp)
    best = mc.delta_mse(rs, a)
    for other in (a - 0.3, a - 0.01, a + 0.01, a + 0.3, 0.0, 1.0):
        assert mc.delta_mse(rs, other) <= best + 1e-12


def test_a_perfectly_scaled_prediction_recovers_its_scale():
    rs = rows(("A", [1.0, -2.0, 3.0], [0.25, -0.5, 0.75]), ("B", [0.5, -0.5, 1.0], [0.125, -0.125, 0.25]))
    assert mc.fit_scale(rs) == pytest.approx(0.25)


def test_scale_is_negative_when_the_prediction_points_the_wrong_way():
    rs = rows(("A", [1.0, -2.0], [-0.5, 1.0]))
    assert mc.fit_scale(rs) == pytest.approx(-0.5)


def test_fit_scale_is_none_with_no_prediction_variance():
    assert mc.fit_scale(rows(("A", [0.0, 0.0], [1.0, -1.0]))) is None


def test_training_delta_mse_at_the_fitted_scale_is_never_negative():
    rs = rows(("A", [1.0, -1.0, 0.5], [-3.0, 2.0, 1.0]), ("B", [2.0, 0.5, -1.5], [0.2, -0.1, 0.4]))
    assert mc.delta_mse(rs, mc.fit_scale(rs)) >= -1e-12


# ---- the split -----------------------------------------------------------------------------------
def test_split_is_saturables_convention():
    sat = _script("saturable_gate")
    for s in ("OCA2", "HERC2", "TP53", "ADAR", "AGO3", "A3GALT2", "MYC"):
        assert mc.is_held_out(s) == sat.is_held_out(s)
        want = int(hashlib.sha256(s.encode("utf-8")).hexdigest()[:8], 16) % 2 == 1
        assert mc.is_held_out(s) is want


def test_split_partitions_every_gene_once_and_whole():
    rs = rows(*[(f"G{i}", [float(i)] * 13, [1.0] * 13) for i in range(200)])
    train, held = mc.split_rows(rs)
    assert len(train) + len(held) == len(rs)
    assert not ({r["symbol"] for r in train} & {r["symbol"] for r in held})
    assert all(len(r["D"]) == 13 for r in train + held)
    assert 0 < len(held) < len(rs)


def test_the_fit_sees_no_held_out_gene():
    rs = rows(*[(f"G{i}", [1.0, -1.0], [0.3, -0.3]) for i in range(60)])
    train, held = mc.split_rows(rs)
    poisoned = [{**r, "D": [50.0, -50.0]} for r in held]
    assert mc.fit_scale(train) != mc.fit_scale(train + poisoned)


# ---- the four-way table, the supervisor's ruling --------------------------------------------------
def test_held_out_interval_decides_on_the_rules_side():
    assert mc.decide(0.1, 0.4, -0.5, 0.5) == "BETTER"
    assert mc.decide(0.1, 0.4, 0.2, 0.4) == "BETTER"


def test_held_out_interval_decides_on_the_baselines_side_and_that_is_the_falsifier():
    assert mc.decide(-0.4, -0.1, -0.5, 0.5) == "WORSE"
    assert mc.decide(-0.4, -0.1, 0.2, 0.4) == "WORSE"


def test_a_straddle_with_a_dead_scale_is_the_negative_not_a_null():
    assert mc.decide(-0.2, 0.3, -0.05, 0.05) == "CALIBRATION_REMOVES"
    assert "no magnitude of the compiled rules improves on no rule" in mc.OUTCOMES["CALIBRATION_REMOVES"]
    assert "negative result" in mc.OUTCOMES["CALIBRATION_REMOVES"]
    assert "is not 'no difference detected'" in mc.OUTCOMES["CALIBRATION_REMOVES"]


def test_a_straddle_with_a_live_scale_is_inconclusive():
    assert mc.decide(-0.2, 0.3, 0.1, 0.4) == "INCONCLUSIVE"
    assert mc.OUTCOMES["INCONCLUSIVE"].startswith("{form}: inconclusive")


def test_a_negative_scale_is_reported_beside_whichever_reading_applies():
    assert mc.scale_is_negative(-0.4, -0.1)
    assert not mc.scale_is_negative(-0.4, 0.1)
    assert "the rules point the wrong way on average" in mc.SCALE_CLAUSES["negative"]


def test_every_verdict_sentence_carries_the_scale_and_the_interval():
    for kind in ("BETTER", "WORSE", "CALIBRATION_REMOVES", "INCONCLUSIVE", "NOT_TESTABLE"):
        assert "{scale}" in mc.OUTCOMES[kind]
    assert "{a}" in mc.SCALE_CLAUSES["fitted"] and "{alo}" in mc.SCALE_CLAUSES["fitted"]


def test_combined_reading_always_ends_with_the_prior_exposure_sentence():
    per = {f: {"reading": f"{f} clause", "outcome": "INCONCLUSIVE", "scale": 0.1} for f in mc.FORMS}
    out = mc.combined_reading(per)
    assert out.endswith(mc.PRIOR_EXPOSURE_CLAUSE)
    assert out.index("mean clause") < out.index("competition clause") < out.index("multiplicative clause")


# ---- the control ---------------------------------------------------------------------------------
def test_the_control_refits_the_scale_and_keeps_each_genes_own_values():
    train = rows(*[(f"T{i}", [float(i % 5) - 2, 1.0, -1.0], [0.3, -0.2, 0.5]) for i in range(40)])
    held = rows(*[(f"H{i}", [float(i % 3) - 1, 0.5, -0.5], [0.1, 0.4, -0.2]) for i in range(40)])
    perm = mc.permuted_control(train, held, 25, 7)
    assert len(perm) == 25
    assert all(math.isfinite(x) for x in perm)
    assert mc.permuted_control(train, held, 25, 7) == perm
    assert mc.permuted_control(train, held, 25, 8) != perm


def test_the_control_leaves_an_uninformative_form_at_zero():
    flat = rows(*[(f"G{i}", [0.0] * 3, [1.0, -1.0, 0.0]) for i in range(20)])
    assert all(abs(x) < 1e-12 for x in mc.permuted_control(flat, flat, 5, 1))


# ---- the registration says in advance what the brief required -------------------------------------
def test_prior_exposure_names_every_earlier_reader_of_this_rna():
    t = mc.PRIOR_EXPOSURE
    for s in ("item 7", "7d", "7h", "capped-sum census", "not fresh", "does not eliminate"):
        assert s in t, s
    assert "item 7" in mc.PRIOR_EXPOSURE_CLAUSE or "items 7" in mc.PRIOR_EXPOSURE_CLAUSE


def test_one_number_per_form_is_stated():
    assert "one real number" in mc.PREDICTIONS["calibrated"]
    assert "No per-gene parameter, no per-cell parameter" in mc.PREDICTIONS["calibrated"]
    assert "training pairs only" in mc.FIT


def test_the_falsifier_and_the_inconclusive_case_are_registered():
    assert "excludes 0 on the BASELINE's side" in mc.FALSIFIER
    assert "straddling 0 falsifies nothing" in mc.FALSIFIER
    assert "MANUFACTURE A STRADDLE but it cannot manufacture a win" in mc.NOT_VACUOUS


def test_which_governs_carries_the_four_way_ruling():
    t = mc.WHICH_GOVERNS
    assert t.startswith("THE HELD-OUT INTERVAL DECIDES")
    for s in ('"better"', '"worse"', '"inconclusive"', "the rules point the wrong way on average"):
        assert s in t, s


def test_endpoint_baseline_and_vocabulary_are_7ds_and_7hs_verbatim():
    assert mc.MEASURED == c2.MEASURED
    assert c2.PREDICTIONS["no_rule_baseline"] == mc.NO_RULE_BASELINE
    assert mc.MIN_GENES == c2.MIN_GENES == 100
    assert mc.CELLS == c2.CELLS and len(mc.CELLS) == 13
    assert mc.CARRIED["seven_h_control"] == mr.CONTROL
    assert mc.BOOTSTRAP["seed"] == c2.BOOTSTRAP["seed"]
    assert mc.BOOTSTRAP["unit"] == c2.BOOTSTRAP["unit"]
    for kind in ("BETTER", "WORSE"):
        assert "measured log2(TPM+1) across the thirteen reader cell types" in mc.OUTCOMES[kind]
    assert "closer to measured RNA than" in mc.OUTCOMES["BETTER"]
    assert "further from measured RNA than" in mc.OUTCOMES["WORSE"]
    assert "explains expression" not in " ".join(mc.OUTCOMES.values())


def test_the_three_forms_and_their_sources_are_named():
    assert mc.FORMS == ("mean", "competition", "multiplicative")
    assert mc.SOURCES["multiplicative"] == "data/results/multrna_result.json"
    assert set(mc.FORM_WORDS) >= set(mc.FORMS)


def test_bound_fields_all_exist_and_registered_fields_covers_them():
    fields = mc.registered_fields()
    assert set(fields) == set(mc.BOUND)
    for name in mc.BOUND.values():
        assert getattr(mc, name) is not None


# ---- the run's harness on the committed tables ----------------------------------------------------
def test_the_harness_reproduces_every_committed_delta_mse():
    run = _script("magcal_run")
    rs = run.load_rows()
    hc = run.harness(rs)
    assert hc["problems"] == []
    for form in mc.FORMS:
        per = hc["per_form"][form]
        assert per["abs_difference"] < mc.TOLERANCE
        assert per["recomputed_delta_mse_at_scale_1"] < 0  # all three read worse uncalibrated


def test_the_harness_catches_reversed_measurements():
    run = _script("magcal_run")
    rs = run.load_rows()
    rs["mean"] = [{**r, "D": [-x for x in r["D"]]} for r in rs["mean"]]
    assert any("Delta MSE" in p for p in run.harness(rs)["problems"])


def test_a_within_gene_shift_of_D_cannot_move_delta_mse_because_both_are_centred():
    rs = rows(("A", [1.0, -1.0, 0.0], [0.4, -0.2, -0.2]))
    assert mc.delta_mse(rs, 1.0) == pytest.approx(
        mc.delta_mse(rows(("A", [1.0, -1.0, 0.0], [3.4, 2.8, 2.8])), 1.0)
    )


def test_the_harness_catches_a_dropped_gene():
    run = _script("magcal_run")
    rs = run.load_rows()
    rs["multiplicative"] = rs["multiplicative"][:-1]
    assert any("genes" in p for p in run.harness(rs)["problems"])
