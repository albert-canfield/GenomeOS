# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-maxform: the algebra of `max` and `or`, and that the score machinery is lane-magcal's."""

from __future__ import annotations

import json
import math
import random

import pytest

from genomeos.attribution import combine2 as c2
from genomeos.attribution import magcal as mc
from genomeos.attribution import maxform as xf
from genomeos.runtime import grn

H = xf.H


def _rho(k_act: int, k_inh: int, rng: random.Random) -> dict[str, float]:
    rho = {f"A{i}": rng.uniform(0.4, 0.99) for i in range(k_act)}
    rho.update({f"I{j}": rng.uniform(1.01, 1.3) for j in range(k_inh)})
    return rho


# ---- the two forms are the engine's, not a reimplementation --------------------------------------


def test_both_forms_are_in_the_engines_rule_list() -> None:
    assert set(xf.FORMS) <= set(grn.ACTIVATOR_RULES)


def test_drive_is_the_engines_combine_activators() -> None:
    rng = random.Random(7)
    for _ in range(200):
        rho = _rho(rng.randint(1, 4), rng.randint(0, 1), rng)
        for form in xf.FORMS:
            f = xf.fit(form, rho)
            if not f.ok:
                continue
            open_ = {el: rng.random() < 0.5 for el in rho}
            terms = [t if open_.get(el, False) else 0.0 for el, t in f.terms.items()]
            assert xf.drive(f, open_) == pytest.approx(grn.combine_activators(terms, form), rel=0, abs=1e-15)


# ---- `max`: no fit with two activators, the mean exactly with one ---------------------------------


def test_max_has_no_fit_with_two_or_more_activating_elements() -> None:
    rng = random.Random(11)
    for _ in range(300):
        rho = _rho(rng.randint(2, 5), rng.randint(0, 2), rng)
        f = xf.fit("max", rho)
        assert not f.ok
        assert f.reason == "max_needs_one_activating_element"


def test_max_with_one_activator_is_the_runtime_mean_element_for_element() -> None:
    rng = random.Random(13)
    checked = 0
    for _ in range(400):
        rho = _rho(1, rng.randint(1, 3), rng)
        f, mean = xf.fit("max", rho), c2.fit("mean", rho)
        assert f.ok == mean.ok, rho
        if not f.ok:
            continue
        checked += 1
        assert f.q == pytest.approx(mean.q, abs=1e-15)
        assert f.vmax_frac == pytest.approx(mean.vmax_frac, abs=1e-12)
        for el in rho:
            assert f.s[el] == pytest.approx(mean.s[el], abs=1e-15)
        for _ in range(8):
            open_ = {el: rng.random() < 0.5 for el in rho}
            assert xf.level(f, open_) == pytest.approx(c2.level(mean, open_), abs=1e-15)
    assert checked > 50


def test_the_census_question_for_max_is_zero_by_construction() -> None:
    """PREMISE: on 7d's population `max` either has no fit or gives the mean's prediction exactly."""
    checks = xf.max_identity_checks(draws=3000)
    assert checks["max_fits_two_activators"] == 0.0
    assert checks["max_differs_from_mean_at_one"] == 0.0


# ---- `or`: the closed form, its uniqueness, and that it is not the mean ---------------------------


def test_or_reproduces_every_removal_observation() -> None:
    rng = random.Random(17)
    checked = 0
    for _ in range(400):
        rho = _rho(rng.randint(2, 4), rng.randint(0, 2), rng)
        f = xf.fit("or", rho)
        if not f.ok:
            continue
        checked += 1
        allopen = dict.fromkeys(rho, True)
        base = xf.level(f, allopen)
        for el in rho:
            held = {**allopen, el: False}
            assert xf.level(f, held) / base == pytest.approx(rho[el], abs=1e-12)
    assert checked > 50


def test_or_strengths_and_declared_strength_convention() -> None:
    rng = random.Random(19)
    for _ in range(300):
        rho = _rho(rng.randint(2, 5), 0, rng)
        f = xf.fit("or", rho)
        if not f.ok:
            continue
        top = max(f.weights.values())
        for el, d in f.weights.items():
            assert 0.0 < f.s[el] <= 1.0
            assert f.terms[el] == pytest.approx(d / (top + d), abs=1e-15)
        assert max(f.s[el] for el in f.weights) == pytest.approx(1.0, abs=1e-15)
        assert max(f.weights, key=lambda e: f.weights[e]) == max(f.s, key=lambda e: f.s[e])


def test_or_q_is_the_product_form_and_at_least_the_means() -> None:
    rng = random.Random(23)
    for _ in range(300):
        k = rng.randint(1, 4)
        rho = _rho(k, 0, rng)
        f, mean = xf.fit("or", rho), c2.fit("mean", rho)
        d = {el: 1.0 - r for el, r in rho.items()}
        top = max(d.values())
        want = top * (math.prod(1.0 + x / top for x in d.values()) - 1.0)
        if f.ok:
            assert f.q == pytest.approx(want, rel=1e-12)
            if mean.ok:
                assert f.q >= mean.q - 1e-12
                if k == 1:
                    assert f.q == pytest.approx(mean.q, abs=1e-15)
        else:
            assert f.reason in ("basal_not_positive", "inhibitor_strength_above_1")


def test_every_form_agrees_on_a_single_removal_and_or_differs_once_two_are_closed() -> None:
    """FIT['where_the_forms_can_differ']: both forms are fitted to the single-removal observations, so
    they agree there by construction; `or` separates from the mean only once two or more of a gene's
    activating elements are closed in the same cell."""
    rho = {"A": 0.6, "B": 0.85}
    f, mean = xf.fit("or", rho), c2.fit("mean", rho)
    assert f.ok and mean.ok
    for one in ({"A": True, "B": False}, {"A": False, "B": True}, {"A": True, "B": True}):
        assert xf.level(f, one) == pytest.approx(c2.level(mean, one), abs=1e-12)
    both_closed = {"A": False, "B": False}
    assert abs(math.log2(xf.level(f, both_closed)) - math.log2(c2.level(mean, both_closed))) == pytest.approx(
        0.5849625007211562, abs=1e-9
    )
    three = {"A": 0.6, "B": 0.85, "C": 0.9}
    g, m3 = xf.fit("or", three), c2.fit("mean", three)
    assert g.ok and m3.ok
    two_closed = {"A": True, "B": False, "C": False}
    assert abs(math.log2(xf.level(g, two_closed)) - math.log2(c2.level(m3, two_closed))) > 1e-3


def test_fit_refuses_a_form_it_does_not_carry() -> None:
    with pytest.raises(ValueError):
        xf.fit("sum_capped", {"A": 0.5, "B": 0.6})


# ---- the registration -----------------------------------------------------------------------------


def test_form_interval_is_98_75_percent_over_a_family_of_four() -> None:
    lo, hi = xf.FORM_INTERVAL
    assert (lo, hi) == (0.625, 99.375)
    assert hi - lo == pytest.approx(98.75)
    assert "four" in xf.FAMILY


def test_registered_fields_cover_bound_and_serialise() -> None:
    fields = xf.registered_fields()
    assert set(fields) == set(xf.BOUND)
    json.dumps(fields)
    assert fields["prior_exposure"] == xf.PRIOR_EXPOSURE
    for item in ("item 7", "7d", "7h", "M2h"):
        assert item in xf.PRIOR_EXPOSURE
    assert xf.POPULATION == c2.POPULATION
    assert xf.SPLIT == mc.SPLIT
    assert xf.OUTCOMES == mc.OUTCOMES


def test_inconclusive_and_falsifier_are_stated() -> None:
    assert "INCONCLUSIVE" in xf.INCONCLUSIVE
    assert "CALIBRATION_REMOVES" in xf.INCONCLUSIVE
    assert "NOT_TESTABLE" in xf.INCONCLUSIVE
    assert "98.75%" in xf.FALSIFIER
    assert "mean(D^2)" in xf.HONEST_SIZE


def test_census_decision_is_7d_min_genes() -> None:
    few = {"or": {"differs_from_the_mean": {"informative_symbol_on_one_gene_id": 99}}}
    many = {"or": {"differs_from_the_mean": {"informative_symbol_on_one_gene_id": 100}}}
    assert xf.MIN_GENES == 100
    assert xf.decide_from_census(few, "or") == "NOT_TESTABLE"
    assert xf.decide_from_census(many, "or") == "register_and_run"


# ---- the score machinery is lane-magcal's ---------------------------------------------------------


def _rows(n: int, seed: int = 5) -> list[dict]:
    rng = random.Random(seed)
    out = []
    for i in range(n):
        p = [rng.gauss(0, 1) for _ in range(13)]
        out.append(
            {
                "gene": f"G{i}",
                "symbol": f"S{i}",
                "chrom": "chr1",
                "start": i * 200_000,
                "P": p,
                "D": [0.3 * x + rng.gauss(0, 1) for x in p],
            }
        )
    return out


def test_bootstrap_reproduces_magcals_at_magcals_interval(monkeypatch: pytest.MonkeyPatch) -> None:
    """BOOTSTRAP['interval']: the same machinery; only the percentile pair is this lane's."""
    rows = _rows(120)
    monkeypatch.setattr(xf, "FORM_INTERVAL", mc.FORM_INTERVAL)
    assert xf.bootstrap_delta(rows, 0.5, 200, 1) == mc.bootstrap_delta(rows, 0.5, 200, 1)
    assert xf.bootstrap_scale(rows, 200, 2) == mc.bootstrap_scale(rows, 200, 2)


def test_this_lanes_interval_is_wider_than_magcals() -> None:
    rows = _rows(120)
    lo, hi = xf.bootstrap_delta(rows, 0.5, 400, 3)
    mlo, mhi = mc.bootstrap_delta(rows, 0.5, 400, 3)
    assert lo <= mlo and hi >= mhi


def test_cluster_bootstrap_is_labelled_post_hoc() -> None:
    rows = _rows(120)
    got = xf.cluster_bootstrap_delta(rows, 0.5, 200, 4)
    assert got["post_hoc"] is True
    assert got["decides"] == "nothing"
    assert got["span_bp"] == 1_000_000
    assert got["clusters"] == len({r["start"] // 1_000_000 for r in rows})


def test_head_to_head_is_descriptive_and_refuses_under_min_genes() -> None:
    rows = {r["gene"]: r for r in _rows(120)}
    few = {k: rows[k] for k in list(rows)[:10]}
    assert xf.head_to_head(rows, few, 100, 5, (1.0, 1.0))["compared"] is False
    got = xf.head_to_head(rows, rows, 200, 6, (1.0, 1.0))
    assert got["compared"] is True
    assert got["descriptive"] is True and got["decides"] == "nothing"
    assert got["mse_other_minus_mse_own"] == pytest.approx(0.0, abs=1e-12)


def test_the_four_way_table_is_magcals() -> None:
    assert xf.WHICH_GOVERNS == mc.WHICH_GOVERNS
    assert mc.decide(0.1, 0.2, 0.1, 0.3) == "BETTER"
    assert mc.decide(-0.2, -0.1, 0.1, 0.3) == "WORSE"
    assert mc.decide(-0.1, 0.2, -0.1, 0.3) == "CALIBRATION_REMOVES"
    assert mc.decide(-0.1, 0.2, 0.1, 0.3) == "INCONCLUSIVE"
