# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-regionscale: the registered design, planted in every direction.

The tests that matter most here are the ones that would catch THE error this lane exists to avoid: a
region's scale fitted on the genes it is then scored on. `test_PLANTED_the_region_scale_does_not_see_a
_single_held_out_value` and `test_PLANTED_perturbing_a_held_out_D_cannot_move_any_region_scale` hold
it directly, and `test_PLANTED_both_arms_are_scored_on_the_identical_rows` holds the other half of the
design. `decide_region` is planted in every one of its branches, and NOT_VACUOUS is planted in both
of its directions on synthetic tables.

The test-local resample and permutation counts are deliberately smaller than the REGISTERED ones; the
registered values live in `regionscale.BOOTSTRAP` and `regionscale.PERMUTATIONS` and are what the run
uses. A test asserts that those registered values are magcal's and 1000 respectively, so a reduced
count here cannot silently become the registered one.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from genomeos import manifest as mf
from genomeos.attribution import hetblock as hb
from genomeos.attribution import magcal as mc
from genomeos.attribution import regionscale as rs
from genomeos.attribution import scalehet as sh

ROOT = Path(__file__).resolve().parents[1]
CELLS = 13


# ---- synthetic tables ----------------------------------------------------------------------------
def _rows(slopes: dict[str, float], per_chrom: int, noise: float, seed: int) -> list[dict]:
    """A synthetic genes table in 7h's shape: one slope per chromosome, D = slope * P + noise."""
    rng = np.random.default_rng(seed)
    out: list[dict] = []
    for ci, (chrom, a) in enumerate(sorted(slopes.items())):
        for g in range(per_chrom):
            p = rng.normal(0.0, 1.0, CELLS)
            d = a * p + rng.normal(0.0, noise, CELLS)
            out.append(
                {
                    "gene": f"G{ci:02d}_{g:05d}",
                    "symbol": f"SYM{ci:02d}_{g:05d}",
                    "chrom": chrom,
                    "P": [float(x) for x in p],
                    "D": [float(x) for x in d],
                }
            )
    return out


def _pipeline(rows: list[dict], resamples: int = 200, perms: int = 100, seed: int = 7) -> dict:
    """The registered design, end to end, at test-local counts. The ONLY place the test fits anything.

    It is written here, in the test, exactly as scripts/regionscale_run.py writes it, so a divergence
    between the two shows up as a failing planted direction rather than as a silent difference.
    """
    counts = rs.split_counts(rows)
    adm, _fall = rs.admitted_regions(counts)
    train = [r for r in rows if not mc.is_held_out(r["symbol"])]
    held = [r for r in rows if mc.is_held_out(r["symbol"])]
    t_dp, t_pp, _t_n = rs.gene_sums(train)
    h_dp, h_pp, h_n = rs.gene_sums(held)
    a_g = rs.scale_from_sums(sum(t_dp), sum(t_pp))
    assert a_g is not None
    per_region, bad = rs.region_scales(t_dp, t_pp, [r["chrom"] for r in train], adm)
    s_r = rs.scales_for([r["chrom"] for r in held], per_region, a_g)
    s_g = [a_g] * len(held)
    dm_r = rs.delta_mse_rowwise(h_dp, h_pp, h_n, s_r)
    dm_g = rs.delta_mse_rowwise(h_dp, h_pp, h_n, s_g)
    boot = rs.bootstrap_paired(h_dp, h_pp, h_n, s_r, s_g, resamples, seed)
    null = rs.permuted_delta(
        t_dp,
        t_pp,
        [r["chrom"] for r in train],
        h_dp,
        h_pp,
        h_n,
        [r["chrom"] for r in held],
        adm,
        a_g,
        perms,
        seed,
    )
    delta = dm_r - dm_g
    p_perm = rs.permutation_p(delta, null)
    outcome = rs.decide_region(
        boot["delta_ci"][0], boot["delta_ci"][1], boot["delta_ci_nominal"][0], p_perm, len(adm), bad
    )
    return {
        "admitted": adm,
        "a_global": a_g,
        "per_region": per_region,
        "delta": delta,
        "dm_region": dm_r,
        "dm_global": dm_g,
        "boot": boot,
        "p_perm": p_perm,
        "outcome": outcome,
        "null": null,
        "held": held,
        "train": train,
    }


# ---- the registration is bound, complete and non-empty -------------------------------------------
def test_registered_fields_cover_BOUND_and_none_is_empty():
    f = rs.registered_fields()
    assert set(f) == set(rs.BOUND)
    for k, v in f.items():
        assert v not in (None, "", [], {}), k


def test_the_registered_bootstrap_and_permutation_counts_are_magcals_and_1000():
    assert rs.BOOTSTRAP["resamples"] == mc.BOOTSTRAP["resamples"]
    assert rs.BOOTSTRAP["seed"] == mc.BOOTSTRAP["seed"]
    assert rs.PERMUTATIONS["count"] == mc.PERMUTATIONS["count"] == 1000
    assert rs.FORM_INTERVAL == mc.FORM_INTERVAL
    assert rs.CONTROL_ALPHA == mc.CONTROL_ALPHA == 0.05
    assert rs.MIN_GENES == mc.MIN_GENES == 100
    assert rs.MIN_REGIONS == sh.MIN_CHROMOSOMES == 10


def test_the_exposure_is_stated_TWICE_and_names_every_reader_of_this_RNA():
    both = rs.PRIOR_EXPOSURE[0] + " " + rs.PRIOR_EXPOSURE_CLAUSE
    for who in (
        "item 7",
        "7d",
        "7h",
        "capped-sum census",
        "M2h",
        "M2i",
        "scalesplit",
        "recovered",
        "scalehet",
        "hetblock2",
    ):
        assert who in both, who
    assert "LOOKED AT, TESTED, AND TESTED AGAIN" in rs.PRIOR_EXPOSURE[0]
    assert "LOOKED AT, TESTED, AND TESTED AGAIN" in rs.PRIOR_EXPOSURE_CLAUSE
    assert len(rs.PRIOR_EXPOSURE) == 2


def test_the_honest_other_half_of_hetblock2s_reading_is_carried_not_strengthened():
    # hetblock2's own registered caveat: at 0.0167 the primary Q clears and the secondary R does NOT
    for text in (
        rs.OUTCOMES["NO_IMPROVEMENT"],
        rs.OUTCOMES["WORSE"],
        rs.OUTCOMES["IMPROVES"],
        rs.OUTCOMES["INCONCLUSIVE_NOT_TESTABLE"],
    ):
        assert "0.0167" in text
        assert "the secondary R does NOT" in text
    assert str(hb.LANE)  # hetblock's module is imported read-only and not edited by this lane


def test_the_regional_unit_is_the_chromosome_and_says_why():
    assert "THE CHROMOSOME" in rs.REGIONAL_UNIT
    assert "the unit the heterogeneity was MEASURED in" in rs.REGIONAL_UNIT


def test_the_held_out_design_states_the_impossibility_and_how_it_is_resolved():
    assert "leaving the region out" in rs.HELD_OUT_DESIGN
    assert "WITHIN-REGION SPLIT" in rs.HELD_OUT_DESIGN
    assert "DISJOINT, PER REGION, BY CONSTRUCTION" in rs.HELD_OUT_DESIGN
    assert "cannot be tested out of sample at this gene count" in rs.HELD_OUT_DESIGN
    low = rs.OUTCOMES["INCONCLUSIVE_NOT_TESTABLE"].lower()
    assert "cannot be tested out of sample at this gene count" in low


def test_both_outcomes_are_pre_written_and_say_what_no_improvement_means():
    t = rs.OUTCOMES["NO_IMPROVEMENT"]
    assert "REAL BUT NOT USABLE" in t
    assert "SINGLE GLOBAL" in t and "WHOLE CORRECTION IN PRACTICE" in t
    assert "COMPLETE ANSWER AND NOT AN ABSENCE OF ONE" in t
    assert "STAND UNTOUCHED" in t
    assert "NO PER-REGION MODEL IS" in t


def test_no_further_model_is_proposed_in_any_outcome():
    for kind, text in rs.OUTCOMES.items():
        if kind == "HARNESS_DEFECT":
            continue
        assert "NO PER-REGION MODEL IS" in text or "PROPOSES NO FURTHER MODEL" in text, kind
    assert "NO FURTHER MODEL IS PROPOSED" in rs.STOP_RULE


def test_the_fraction_of_variance_is_carried_beside_every_better_and_worse():
    assert "OF THE VARIANCE" in rs.OUTCOMES["IMPROVES"]
    assert "{pct_delta:.2f}" in rs.OUTCOMES["IMPROVES"]
    assert "{pct_delta_abs:.2f}" in rs.OUTCOMES["WORSE"]
    for kind in ("IMPROVES", "NO_IMPROVEMENT", "WORSE", "MIXED"):
        assert "{pct_region:.2f}" in rs.OUTCOMES[kind], kind
        assert "{pct_global:.2f}" in rs.OUTCOMES[kind], kind


def test_the_falsifier_names_what_inconclusive_is_and_what_it_is_not():
    f = rs.FALSIFIER
    assert "NO_IMPROVEMENT and WORSE" in f
    assert "A small DELTA is NOT inconclusive" in f
    assert "NO_IMPROVEMENT DOES NOT FALSIFY scalehet OR hetblock2" in f


def test_the_interval_and_its_correction_are_decided_in_the_registration():
    assert "DECIDED NOW AND NOT AFTER" in rs.INTERVAL
    assert "98.33%" in rs.INTERVAL and "95%" in rs.INTERVAL
    assert "CORRECTION, DECIDED NOW" in rs.ALPHA
    assert "NO Bonferroni correction is applied" in rs.ALPHA
    assert rs.NOMINAL_INTERVAL == (2.5, 97.5)


# ---- the split, the admission and the feasibility -------------------------------------------------
def test_split_counts_are_counts_and_sum_to_the_table():
    rows = _rows({"chrA": 0.3, "chrB": 0.3}, 120, 0.5, 1)
    c = rs.split_counts(rows)
    assert sum(d["train"] + d["held"] for d in c.values()) == len(rows)
    assert set(c) == {"chrA", "chrB"}


def test_admission_needs_MIN_GENES_in_BOTH_halves_and_a_short_half_falls_back():
    counts = {
        "chrBig": {"train": 400, "held": 400},
        "chrShortTrain": {"train": 40, "held": 400},
        "chrShortHeld": {"train": 400, "held": 40},
        "chrTiny": {"train": 5, "held": 5},
    }
    adm, fall = rs.admitted_regions(counts)
    assert adm == ["chrBig"]
    assert fall == ["chrShortHeld", "chrShortTrain", "chrTiny"]


def test_is_testable_is_the_registered_floor_on_the_number_of_regions():
    assert rs.is_testable([f"c{i}" for i in range(rs.MIN_REGIONS)])
    assert not rs.is_testable([f"c{i}" for i in range(rs.MIN_REGIONS - 1)])


def test_a_non_admitted_region_is_scored_by_the_global_scale_and_is_NOT_dropped():
    labels = ["chrA", "chrB", "chrA", "chrZ"]
    s = rs.scales_for(labels, {"chrA": 0.9, "chrB": 0.1}, 0.5)
    assert s == [0.9, 0.1, 0.9, 0.5]
    assert len(s) == len(labels)  # nothing dropped


# ---- the metric and the interval ARE magcal's ----------------------------------------------------
def test_the_rowwise_metric_at_a_constant_scale_IS_magcal_delta_mse():
    rows = _rows({"chrA": 0.4, "chrB": 0.1}, 60, 0.5, 2)
    dps, pps, ns = rs.gene_sums(rows)
    for a in (-0.7, 0.0, 0.2145911381525326, 1.0, 3.3):
        mine = rs.delta_mse_rowwise(dps, pps, ns, [a] * len(rows))
        assert abs(mine - mc.delta_mse(rows, a)) <= rs.TOLERANCE


def test_the_group_fit_IS_magcal_fit_scale():
    rows = _rows({"chrA": 0.4, "chrB": 0.1}, 60, 0.5, 3)
    for chrom in ("chrA", "chrB"):
        sub = [r for r in rows if r["chrom"] == chrom]
        dps, pps, _ = rs.gene_sums(sub)
        mine = rs.scale_from_sums(sum(dps), sum(pps))
        theirs = mc.fit_scale(sub)
        assert mine is not None and theirs is not None
        assert abs(mine - theirs) <= rs.TOLERANCE


def test_the_per_arm_bootstrap_at_a_constant_scale_IS_magcal_bootstrap_delta_within_TOLERANCE():
    rows = _rows({"chrA": 0.4, "chrB": 0.1}, 60, 0.5, 4)
    dps, pps, ns = rs.gene_sums(rows)
    a = 0.2145911381525326
    const = [a] * len(rows)
    boot = rs.bootstrap_paired(dps, pps, ns, const, const, 300, 20261010)
    want = mc.bootstrap_delta(rows, a, 300, 20261010)
    for mine in (boot["region_arm_ci"], boot["global_arm_ci"]):
        assert max(abs(mine[0] - want[0]), abs(mine[1] - want[1])) <= rs.TOLERANCE
    assert boot["region_arm_ci"] == boot["global_arm_ci"]


def test_PLANTED_the_residual_against_magcal_is_a_SUMMATION_ORDER_artefact_not_an_algebraic_one():
    """HARNESS (4)'s recorded reason, pinned: accumulate each gene's term in magcal's own per-cell
    order and this lane's draw sequence reproduces magcal.bootstrap_delta EXACTLY; accumulate
    sum(P D) and sum(P^2) first, as this lane must for its permutation null, and the same draws give
    the same interval to about 1e-17. The difference is in the order of the additions, not the
    algebra."""
    rows = _rows({"chrA": 0.4, "chrB": 0.1}, 60, 0.5, 4)
    a = 0.2145911381525326
    resamples, seed = 300, 20261010
    want = mc.bootstrap_delta(rows, a, resamples, seed)

    # magcal's own per-cell accumulation, then THIS LANE's draw sequence
    percell = []
    for r in rows:
        t = 0.0
        for p, d in zip(r["P"], r["D"], strict=True):
            t += 2.0 * a * p * d - a * a * p * p
        percell.append(t)
    ta = np.array(percell)
    na = np.array([float(len(r["D"])) for r in rows])
    rng = np.random.default_rng(seed)
    stats = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, len(rows), len(rows))
        stats[i] = ta[idx].sum() / na[idx].sum()
    lo, hi = np.percentile(stats, list(rs.FORM_INTERVAL))
    assert [float(lo), float(hi)] == [want[0], want[1]]  # EXACTLY, so the draws are magcal's

    # this lane's group-sums order: the same interval, to a floating-point residual only
    dps, pps, ns = rs.gene_sums(rows)
    boot = rs.bootstrap_paired(dps, pps, ns, [a] * len(rows), [a] * len(rows), resamples, seed)
    resid = max(abs(boot["region_arm_ci"][0] - want[0]), abs(boot["region_arm_ci"][1] - want[1]))
    assert resid <= 1e-12, resid
    assert resid <= rs.TOLERANCE


def test_two_identical_arms_give_a_DELTA_interval_of_exactly_zero():
    rows = _rows({"chrA": 0.4, "chrB": 0.1}, 60, 0.5, 5)
    dps, pps, ns = rs.gene_sums(rows)
    const = [0.3] * len(rows)
    boot = rs.bootstrap_paired(dps, pps, ns, const, const, 200, 11)
    assert boot["delta_ci"] == [0.0, 0.0]
    assert boot["delta_ci_nominal"] == [0.0, 0.0]


def test_the_paired_bootstrap_uses_the_SAME_draws_for_both_arms():
    # if the arms were drawn independently, the DELTA interval of two arms differing by a constant
    # offset would be wide; paired, it is exactly that offset.
    rows = _rows({"chrA": 0.4}, 80, 0.5, 6)
    dps, pps, ns = rs.gene_sums(rows)
    boot = rs.bootstrap_paired(dps, pps, ns, [0.4] * len(rows), [0.2] * len(rows), 200, 12)
    lo, hi = boot["delta_ci"]
    rlo, rhi = boot["region_arm_ci"]
    glo, ghi = boot["global_arm_ci"]
    assert (hi - lo) < (rhi - rlo) + (ghi - glo)


# ---- THE ERROR THIS LANE EXISTS TO AVOID ---------------------------------------------------------
def test_PLANTED_the_region_scale_does_not_see_a_single_held_out_value():
    rows = _rows({f"chr{i:02d}": 0.3 for i in range(12)}, 260, 0.4, 20)
    out = _pipeline(rows, resamples=60, perms=20)
    train = out["train"]
    for chrom, a in out["per_region"].items():
        sub = [r for r in train if r["chrom"] == chrom]
        theirs = mc.fit_scale(sub)
        assert theirs is not None
        assert abs(a - theirs) <= rs.TOLERANCE, chrom


def test_PLANTED_perturbing_a_held_out_D_cannot_move_any_region_scale():
    rows = _rows({f"chr{i:02d}": 0.3 for i in range(12)}, 260, 0.4, 21)
    before = _pipeline(rows, resamples=40, perms=10)["per_region"]
    poisoned = [dict(r) for r in rows]
    for r in poisoned:
        if mc.is_held_out(r["symbol"]):
            r["D"] = [d * 100.0 + 7.0 for d in r["D"]]
    after = _pipeline(poisoned, resamples=40, perms=10)["per_region"]
    assert set(before) == set(after)
    for chrom in before:
        assert before[chrom] == after[chrom], chrom


def test_PLANTED_perturbing_a_TRAINING_D_does_move_the_region_scale_so_the_test_above_is_not_vacuous():
    rows = _rows({f"chr{i:02d}": 0.3 for i in range(12)}, 260, 0.4, 22)
    before = _pipeline(rows, resamples=40, perms=10)["per_region"]
    poisoned = [dict(r) for r in rows]
    for r in poisoned:
        if not mc.is_held_out(r["symbol"]) and r["chrom"] == "chr00":
            r["D"] = [d * 5.0 for d in r["D"]]
    after = _pipeline(poisoned, resamples=40, perms=10)["per_region"]
    assert before["chr00"] != after["chr00"]


def test_PLANTED_both_arms_are_scored_on_the_identical_rows():
    rows = _rows({f"chr{i:02d}": 0.3 for i in range(11)} | {"chrSmall": 0.3}, 240, 0.4, 23)
    small = [r for r in rows if r["chrom"] == "chrSmall"][:40]
    rows = [r for r in rows if r["chrom"] != "chrSmall"] + small
    out = _pipeline(rows, resamples=40, perms=10)
    assert "chrSmall" not in out["admitted"]
    # the fallback chromosome's held-out genes are still in the scored set
    held_chroms = {r["chrom"] for r in out["held"]}
    assert "chrSmall" in held_chroms
    # and both arms used the same denominator, so their Delta MSEs are comparable row for row
    _dp, _pp, ns = rs.gene_sums(out["held"])
    assert sum(ns) == len(out["held"]) * CELLS


def test_PLANTED_the_global_arm_IS_a_single_scale_fitted_on_all_training_genes():
    rows = _rows({f"chr{i:02d}": 0.3 for i in range(12)}, 240, 0.4, 24)
    out = _pipeline(rows, resamples=40, perms=10)
    theirs = mc.fit_scale(out["train"])
    assert theirs is not None
    assert abs(out["a_global"] - theirs) <= rs.TOLERANCE


# ---- the null ------------------------------------------------------------------------------------
def test_the_null_holds_every_groups_two_counts_exactly_and_leaves_the_admitted_set_alone():
    assert "held EXACTLY at their observed values" in rs.NULL
    assert "the admitted set is" in rs.NULL and "unchanged" in rs.NULL
    assert "invariant under any permutation of labels" in rs.NULL


def test_PLANTED_on_a_homogeneous_table_the_null_delta_sits_at_or_below_zero():
    rows = _rows({f"chr{i:02d}": 0.3 for i in range(12)}, 240, 0.4, 25)
    out = _pipeline(rows, resamples=40, perms=120)
    null = out["null"]
    assert len(null) == 120
    assert float(np.median(null)) <= 0.0


def test_the_permutation_p_is_magcals_one_sided_form_and_never_zero():
    assert rs.permutation_p(1.0, [0.0] * 10) == pytest.approx(1 / 11)
    assert rs.permutation_p(-1.0, [0.0] * 10) == pytest.approx(11 / 11)
    assert rs.permutation_p(5.0, [0.0] * 999) > 0.0


# ---- decide_region, planted in EVERY direction ----------------------------------------------------
def test_decide_region_IMPROVES_needs_the_interval_AND_the_permutation():
    assert rs.decide_region(0.01, 0.05, 0.015, 0.001, 22) == "IMPROVES"


def test_decide_region_MIXED_when_the_interval_clears_but_the_permutation_does_not():
    assert rs.decide_region(0.01, 0.05, 0.015, 0.40, 22) == "MIXED"
    assert rs.why_mixed(0.01, 0.40) == rs.MIXED_CLAUSES["interval_without_permutation"]


def test_decide_region_MIXED_when_the_permutation_clears_but_the_interval_does_not():
    assert rs.decide_region(-0.01, 0.05, -0.005, 0.001, 22) == "MIXED"
    assert rs.why_mixed(-0.01, 0.001) == rs.MIXED_CLAUSES["permutation_without_interval"]


def test_decide_region_MIXED_when_the_CHOSEN_WIDTH_and_not_the_evidence_would_decide():
    # the registered 98.33% interval includes 0, the nominal 95% interval does not
    assert rs.decide_region(-0.001, 0.05, 0.002, 0.40, 22) == "MIXED"
    assert rs.why_mixed(-0.001, 0.40) == rs.MIXED_CLAUSES["width_decided_it"]


def test_decide_region_NO_IMPROVEMENT_when_nothing_clears():
    assert rs.decide_region(-0.01, 0.01, -0.008, 0.40, 22) == "NO_IMPROVEMENT"


def test_decide_region_WORSE_when_the_interval_lies_wholly_below_zero():
    assert rs.decide_region(-0.05, -0.01, -0.045, 0.90, 22) == "WORSE"


def test_decide_region_INCONCLUSIVE_NOT_TESTABLE_below_the_region_floor():
    assert rs.decide_region(0.01, 0.05, 0.015, 0.001, rs.MIN_REGIONS - 1) == ("INCONCLUSIVE_NOT_TESTABLE")


def test_decide_region_INCONCLUSIVE_DEGENERATE_when_a_fit_is_undefined():
    assert rs.decide_region(0.01, 0.05, 0.015, 0.001, 22, ["chr7"]) == "INCONCLUSIVE_DEGENERATE"


def test_the_width_can_only_turn_a_NEGATIVE_verdict_into_MIXED_never_into_a_positive_one():
    # a nominal interval clearing 0 never yields IMPROVES on its own
    for p in (0.001, 0.4):
        assert rs.decide_region(-0.001, 0.05, 0.002, p, 22) == "MIXED"


def test_an_undefined_group_fit_is_reported_and_not_substituted():
    dps = [0.0, 0.0]
    pps = [0.0, 0.0]
    per, bad = rs.region_scales(dps, pps, ["chrA", "chrA"], ["chrA"])
    assert per == {} and bad == ["chrA"]


# ---- NOT_VACUOUS, planted in both directions ------------------------------------------------------
def test_PLANTED_a_real_per_region_slope_difference_reads_IMPROVES():
    slopes = {f"chr{i:02d}": (0.05 if i % 2 == 0 else 0.55) for i in range(12)}
    out = _pipeline(_rows(slopes, 300, 0.3, 30), resamples=300, perms=200)
    assert len(out["admitted"]) >= rs.MIN_REGIONS
    assert out["delta"] > 0.0
    assert out["outcome"] == "IMPROVES", (out["outcome"], out["delta"], out["boot"]["delta_ci"])


def test_PLANTED_one_slope_everywhere_does_NOT_read_IMPROVES():
    out = _pipeline(_rows({f"chr{i:02d}": 0.3 for i in range(12)}, 300, 0.3, 31), 300, 200)
    assert len(out["admitted"]) >= rs.MIN_REGIONS
    assert out["outcome"] != "IMPROVES", (out["outcome"], out["delta"], out["boot"]["delta_ci"])


# ---- the readings ---------------------------------------------------------------------------------
_NUMBERS = {
    "held_genes": 8328,
    "held_pairs": 108264,
    "n_admitted": 22,
    "n_fallback": 2,
    "a_global": 0.2145911381525326,
    "delta": 0.001,
    "fw": rs.FORM_INTERVAL_WORDS,
    "delta_lo": -0.001,
    "delta_hi": 0.003,
    "nw": rs.NOMINAL_INTERVAL_WORDS,
    "nom_lo": -0.0008,
    "nom_hi": 0.0028,
    "p_perm": 0.2,
    "n_perm": 1000,
    "dm_region": 0.0443,
    "dm_global": 0.0434,
    "pct_region": 2.05,
    "pct_global": 2.01,
    "pct_delta": 0.04,
    "pct_delta_abs": 0.04,
    "why": rs.MIXED_CLAUSES["width_decided_it"],
    "min_regions": rs.MIN_REGIONS,
    "what": "chr7's training fit is undefined",
}


@pytest.mark.parametrize("kind", sorted(set(rs.OUTCOMES) - {"HARNESS_DEFECT"}))
def test_every_registered_outcome_formats_with_the_numbers_the_run_supplies(kind):
    text = rs.reading(kind, **_NUMBERS)
    assert "{" not in text and "}" not in text


@pytest.mark.parametrize("kind", sorted(set(rs.OUTCOMES) - {"HARNESS_DEFECT"}))
def test_every_reading_carries_the_prior_exposure_clause_ALWAYS(kind):
    assert rs.full_reading(kind, **_NUMBERS).endswith(rs.PRIOR_EXPOSURE_CLAUSE)


def test_the_harness_defect_reading_names_what_went_wrong():
    assert "boom" in rs.reading("HARNESS_DEFECT", what="boom")


def test_quantiles_and_percentile_range_are_the_committed_lanes_own():
    xs = [float(i) for i in range(100)]
    assert rs.quantiles(xs) == mc.quantiles(xs)
    assert rs.percentile_range(xs) == sh.percentile_range(xs)


# ---- this lane edits nobody else's file ------------------------------------------------------------
def test_this_lane_adds_only_its_own_files():
    for p in (
        "genomeos/attribution/regionscale.py",
        "scripts/regionscale_harness.py",
        "tests/test_regionscale.py",
    ):
        assert (ROOT / p).exists(), p


def test_the_stop_rule_forbids_a_second_run_and_any_edit_of_the_prior_lanes():
    assert "ONE RUN" in rs.STOP_RULE
    assert "No file of magcal's, scalehet's," in rs.STOP_RULE
    assert "read-only" in rs.STOP_RULE


def test_no_nan_leaks_from_an_empty_scored_set():
    assert math.isnan(rs.delta_mse_rowwise([], [], [], []))


# ---- the committed harness record ----------------------------------------------------------------
# These tests read data/results/regionscale_harness.json and are committed WITH that record (22e12de).
# No pytest.skip is keyed on it (4f2652c): if the record is absent the tests FAIL, which is correct,
# because the registration and the run both refuse to exist without it.
HARNESS_RESULT = ROOT / "data/results/regionscale_harness.json"


def _harness() -> dict:
    return json.loads(HARNESS_RESULT.read_text())


def test_the_committed_harness_reproduces_M2h_FIRST_and_has_no_problems():
    h = _harness()
    assert h["outcome"] == "HARNESS_REPRODUCES"
    assert h["problems"] == []
    rep = h["m2h_reproduction_FIRST"]
    assert rep["reproduces"] is True
    got, want = rep["recomputed"], rep["committed"]
    assert got["training_genes"] == want["training_genes"] == 8181
    assert got["held_out_genes"] == want["held_out_genes"] == 8328
    assert got["outcome"] == want["outcome"] == "BETTER"
    for k in ("scale", "held_out_delta_mse", "control_p"):
        assert abs(float(got[k]) - float(want[k])) <= rs.M2H_TOLERANCE, k
    assert abs(float(got["scale"]) - rs.GLOBAL_SCALE_COMMITTED) <= rs.M2H_TOLERANCE


def test_the_committed_harness_reading_does_NOT_overstate_the_agreement_with_magcal():
    # the registered claim is agreement WITHIN TOLERANCE; "bit-for-bit" was the overstatement 027d414
    # removed, and it must not come back in the stored sentence.
    h = _harness()
    assert "bit-for-bit" not in h["reading"]
    assert "within TOLERANCE" in h["reading"]
    assert "within TOLERANCE" in rs.HARNESS
    assert "WITHIN TOLERANCE AND NOT BIT-FOR-BIT" in rs.HARNESS


def test_the_committed_harness_shows_this_lanes_metric_IS_magcals():
    p = _harness()["one_path"]
    assert p["metric_abs_difference"] <= rs.TOLERANCE
    assert p["arm_interval_max_abs_difference"] <= rs.TOLERANCE
    assert p["delta_interval_of_two_identical_arms"] == [0.0, 0.0]
    assert p["delta_interval_nominal_of_two_identical_arms"] == [0.0, 0.0]
    assert p["constant_scale_used"] == rs.GLOBAL_SCALE_COMMITTED


def test_the_committed_harness_says_the_held_out_design_CAN_be_made_honest():
    f = _harness()["design_feasibility"]
    assert f["the_held_out_design_is_testable"] is True
    assert f["admitted_count"] >= rs.MIN_REGIONS
    assert f["min_genes_each_half"] == rs.MIN_GENES
    for c in f["admitted"]:
        d = f["genes_per_chromosome_per_half"][c]
        assert d["train"] >= rs.MIN_GENES and d["held"] >= rs.MIN_GENES, c
    for c, d in f["fallback_to_global_both_arms"].items():
        assert d["train"] < rs.MIN_GENES or d["held"] < rs.MIN_GENES, c
    # nothing is dropped: every held-out gene is in the scored set, admitted or fallback
    assert (
        f["held_out_genes_in_admitted"] + f["held_out_genes_in_fallback"]
        == rs.M2H_COMMITTED["held_out_genes"]
    )
    assert f["training_genes_in_admitted"] <= rs.M2H_COMMITTED["training_genes"]


def test_the_committed_harness_confirms_M2hs_about_2_percent_of_the_variance():
    v = _harness()["variance_denominator"]
    assert v["confirms"] is True
    assert 1.0 <= v["m2h_percent_of_variance"] <= 3.0
    assert v["held_out_pairs"] == rs.M2H_COMMITTED["held_out_genes"] * CELLS


def test_the_committed_harness_computed_NO_value_of_this_lanes_design():
    h = _harness()
    assert "no per-chromosome scale was fitted" in h["nothing_of_this_lanes_design_was_read"]
    for forbidden in ("per_chromosome", "delta", "per_region_scales", "null"):
        assert forbidden not in h, forbidden
    assert h["alphagenome_requests"] == 0


def test_the_committed_harness_manifest_is_complete_declares_its_readers_and_uses_relative_paths():
    m = _harness()["result_manifest"]
    assert m["complete"] is True, m.get("problems")
    assert len(m["readers"]) == len(rs.CELLS) == 13
    for e in m["inputs"]:
        assert not e["path"].startswith("/"), e["path"]  # d72f959: input_entry paths RELATIVE
    assert m["code_cleanliness"]["own_code_is_committed"] is True


# ---- the registration script's contract ----------------------------------------------------------
def test_the_register_script_refuses_without_a_reproducing_harness(tmp_path, monkeypatch):
    """The registration cannot be written on a harness that did not reproduce M2h's numbers."""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "regionscale_register", ROOT / "scripts/regionscale_register.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    bad = tmp_path / "h.json"
    bad.write_text(json.dumps({"outcome": "HARNESS_DEFECT", "m2h_reproduction_FIRST": {}}))
    monkeypatch.setattr(mod, "HARNESS_RESULT", str(bad))
    monkeypatch.setattr(mod, "OUT", str(tmp_path / "nope.json"))
    with pytest.raises(SystemExit, match="did not reproduce"):
        mod.refuse_unless_harness_reproduced()
    # and it refuses to overwrite a registration that already exists
    there = tmp_path / "already.json"
    there.write_text("{}")
    monkeypatch.setattr(mod, "OUT", str(there))
    with pytest.raises(SystemExit, match="already exists"):
        mod.refuse_unless_harness_reproduced()


# ---- the committed registration ------------------------------------------------------------------
# Read data/results/regionscale_registration.json, committed WITH these tests (22e12de). No skip is
# keyed on it (4f2652c): the run refuses to start without it, so an absent record must FAIL here.
REGISTRATION = ROOT / "data/results/regionscale_registration.json"


def _registration() -> dict:
    return json.loads(REGISTRATION.read_text())


def test_the_committed_registration_binds_every_field_of_the_module_WORD_FOR_WORD():
    reg = _registration()
    fields = rs.registered_fields()
    drift = [f for f, v in fields.items() if reg.get(f) != v]
    assert drift == [], drift
    assert len(fields) == len(rs.BOUND) == 45


def test_the_committed_registration_pins_every_input_the_run_will_read():
    reg = _registration()
    pinned = reg["pinned_sha256"]
    assert "data/results/multrna_result.json" in pinned
    assert "data/results/regionscale_harness.json" in pinned
    for p, want in pinned.items():
        assert not p.startswith("/"), p
        assert mf.input_entry(p)["sha256"] == want, p


def test_the_committed_registration_carries_the_harness_that_reproduced_M2h():
    # "harness_record", NOT "harness". The module has a BOUND field called "harness" -- the
    # registered text of what the harness must check -- and the writer's provenance block was
    # originally called the same thing, so `**registered_fields()` silently overwrote it and the
    # registration lost the PROOF that M2h reproduced. Reading the wrong key here is what caught
    # it; this test asserts on the evidence, and the one below asserts the description is intact.
    reg = _registration()
    h = reg["harness_record"]
    assert h["outcome"] == "HARNESS_REPRODUCES"
    assert h["m2h_reproduces"] is True
    assert "REPRODUCE" in h["m2h_reading"]
    assert reg["harness"] == rs.registered_fields()["harness"]
    assert isinstance(reg["harness"], str)


def test_the_committed_registration_declares_the_design_feasible_on_COUNTS_alone():
    d = _registration()["the_held_out_design_is_feasible"]
    assert d["testable"] is True
    assert d["admitted_count"] >= rs.MIN_REGIONS
    assert d["min_genes_each_half"] == rs.MIN_GENES
    assert d["min_regions"] == rs.MIN_REGIONS
    assert set(d["admitted"]) & set(d["fallback_to_global_both_arms"]) == set()
    assert "depend on no measured value" in d["these_are_counts"]


def test_the_committed_registration_holds_NO_held_out_value_of_this_design():
    # The names of the result quantities do appear in this file, and they MUST: both outcomes are
    # pre-written as TEMPLATES before any value exists, so `dm_region`, `dm_global` and `p_perm`
    # are present as UNFILLED FORMAT PLACEHOLDERS inside the outcome sentences. A test that merely
    # banned the substring failed on those placeholders; it was MIS-SPECIFIED, not the registration.
    # The correct and stronger assertion is the one below: the placeholders are STILL UNFILLED, and
    # the quantity names occur NOWHERE outside the templates.
    reg = _registration()
    blob = json.dumps(reg)
    for forbidden in ("delta_ci", "per_region_scales", "null_q", "delta_mse_region", "delta_mse_global"):
        assert forbidden not in blob, forbidden
    templated = {"outcomes", "mixed_clauses"}
    elsewhere = json.dumps({k: v for k, v in reg.items() if k not in templated})
    for name, placeholder in (
        ("dm_region", "{dm_region:+.6g}"),
        ("dm_global", "{dm_global:+.6g}"),
        ("p_perm", "{p_perm:.6g}"),
    ):
        assert name not in elsewhere, name
        carriers = [k for k, text in reg["outcomes"].items() if name in text]
        assert carriers, name
        for k in carriers:
            assert placeholder in reg["outcomes"][k], (k, placeholder)
    assert "no per-chromosome scale" in reg["what_is_not_in_this_file"]
    assert reg["alphagenome_requests"] == 0


def test_the_committed_registration_states_the_exposure_twice_and_the_interval_decided_now():
    reg = _registration()
    assert len(reg["prior_exposure"]) == 2
    assert "LOOKED AT, TESTED, AND TESTED AGAIN" in reg["prior_exposure"][0]
    assert "LOOKED AT, TESTED, AND TESTED AGAIN" in reg["prior_exposure_clause"]
    assert "DECIDED NOW AND NOT AFTER" in reg["interval"]
    assert "CORRECTION, DECIDED NOW" in reg["alpha"]


def test_the_committed_registration_pre_writes_BOTH_outcomes_and_proposes_no_model():
    o = _registration()["outcomes"]
    assert "REAL BUT NOT USABLE" in o["NO_IMPROVEMENT"]
    assert "WHOLE CORRECTION IN PRACTICE" in o["NO_IMPROVEMENT"]
    assert "BETTER THAN THE SINGLE GLOBAL SCALE" in o["IMPROVES"]
    assert "cannot be tested out of sample at this gene count" in o["INCONCLUSIVE_NOT_TESTABLE"].lower()
    for kind, text in o.items():
        if kind != "HARNESS_DEFECT":
            assert "NO PER-REGION MODEL IS" in text or "PROPOSES NO FURTHER MODEL" in text, kind


def test_the_committed_registration_manifest_is_complete_and_declares_its_readers():
    m = _registration()["result_manifest"]
    assert m["complete"] is True, m.get("problems")
    assert len(m["readers"]) == 13
    assert m["code_cleanliness"]["own_code_is_committed"] is True


# ---- the run script's contract, and the TRIAL RUN it discloses ------------------------------------
def _run_module():
    import importlib.util

    spec = importlib.util.spec_from_file_location("regionscale_run", ROOT / "scripts/regionscale_run.py")
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_run_script_refuses_to_run_twice_and_refuses_without_a_registration(tmp_path, monkeypatch):
    """STOP_RULE's ONE RUN, and the registration that must exist: both driven through the real script."""
    mod = _run_module()
    there = tmp_path / "already.json"
    there.write_text("{}")
    monkeypatch.setattr(mod, "OUT", str(there))
    with pytest.raises(SystemExit, match="ONE RUN"):
        mod.refuse_unless_registered()
    monkeypatch.setattr(mod, "OUT", str(tmp_path / "nope.json"))
    monkeypatch.setattr(mod, "REGISTRATION", str(tmp_path / "absent.json"))
    with pytest.raises(SystemExit, match="does not exist"):
        mod.refuse_unless_registered()


def test_the_run_script_refuses_a_registration_that_DRIFTED_from_the_module(tmp_path, monkeypatch):
    """A registration that no longer binds a BOUND field word for word is a different experiment."""
    mod = _run_module()
    drifted = dict(rs.registered_fields())
    drifted["question"] = drifted["question"] + " AND SOMETHING ELSE"
    drifted["pinned_sha256"] = {}
    p = tmp_path / "reg.json"
    p.write_text(json.dumps(drifted))
    monkeypatch.setattr(mod, "OUT", str(tmp_path / "nope.json"))
    monkeypatch.setattr(mod, "REGISTRATION", str(p))
    with pytest.raises(SystemExit, match="word for word"):
        mod.refuse_unless_registered()


def test_the_run_script_refuses_when_a_PINNED_INPUT_MOVED(tmp_path, monkeypatch):
    """A moved input makes the registered numbers unreproducible, so the run does not start."""
    mod = _run_module()
    body = dict(rs.registered_fields())
    body["pinned_sha256"] = {"data/results/multrna_result.json": "0" * 64}
    p = tmp_path / "reg.json"
    p.write_text(json.dumps(body))
    monkeypatch.setattr(mod, "OUT", str(tmp_path / "nope.json"))
    monkeypatch.setattr(mod, "REGISTRATION", str(p))
    with pytest.raises(SystemExit, match="pinned input moved"):
        mod.refuse_unless_registered()


def test_the_run_DISCLOSES_the_trial_run_that_preceded_the_registration_commit():
    """The disclosure is a module constant, so it cannot be dropped from one code path only.

    A trial run computed DELTA before the registration was committed. It wrote nothing into
    data/results -- save_result refused on an untracked entry script and quarantined it -- and no value
    of it was read. What makes the registration honest regardless is a TIMESTAMP and not an assurance,
    and the constant has to say so, name the quarantined file's sha256, and state plainly which claim is
    no longer available.
    """
    d = _run_module().TRIAL_RUN_DISCLOSURE
    assert "TRIAL RUN HAPPENED BEFORE THE REGISTRATION WAS COMMITTED" in d
    assert "37fd95a18592207ea70eb2e1d19ea4faaa67959ae6dd7c2f591e19da0fa2f1ea" in d
    assert "NOTHING REACHED data/results" in d
    # the honesty check must be GIT and not a file mtime: the module holding every registered field is
    # unchanged since 027d414, hours before the trial, and a committed test asserts the copy word for word
    assert "WHAT SETTLES THE REGISTRATION'S HONESTY IS GIT" in d
    assert "NOT A FILE MTIME" in d
    assert "027d414" in d and "16:52:54" in d
    assert "no longer claimable" in d.lower()
    assert "none was READ" in d


# ---- the committed result -------------------------------------------------------------------------
# Read data/results/regionscale_result.json, committed BEFORE these tests (22e12de; and the record
# could not be committed with them, because scripts/check_staged.py refuses a result whose manifest says
# the lane's own code -- which declares this file -- was uncommitted when the numbers were computed).
# No skip is keyed on it (4f2652c): an absent record must FAIL here.
RESULT = ROOT / "data/results/regionscale_result.json"


def _result() -> dict:
    return json.loads(RESULT.read_text())


def test_the_committed_result_outcome_is_what_the_REGISTERED_DECISION_RULE_returns():
    """The verdict is recomputed from the record's own numbers, so the rule and not the writer decides."""
    r = _result()
    assert r["outcome"] in rs.OUTCOMES
    # No branch on the record's own content: the committed record HAS measured numbers, so a record
    # that lost them must fail here rather than quietly take a shorter path (the spirit of 4f2652c).
    m = r["measured"]
    lo, hi = m["delta_ci"]
    nom_lo, _nom_hi = m["delta_ci_nominal"]
    assert rs.decide_region(lo, hi, nom_lo, m["p_perm"], r["design"]["admitted_count"]) == r["outcome"]
    # and the disagreement that produced MIXED is the registered clause for THAT disagreement, not
    # another one: the permutation reached ALPHA while the interval did not clear 0
    if r["outcome"] == "MIXED":
        assert rs.why_mixed(lo, m["p_perm"]) in r["reading"]
        assert rs.MIXED_CLAUSES["permutation_without_interval"] in r["reading"]
        assert lo <= 0.0 and m["p_perm"] <= rs.CONTROL_ALPHA


def test_the_committed_result_reading_IS_the_registered_sentence_and_nothing_else():
    """The reading is the registered template with numbers filled, with the exposure clause ALWAYS."""
    r = _result()
    assert r["reading"].endswith(rs.PRIOR_EXPOSURE_CLAUSE)
    template = rs.OUTCOMES[r["outcome"]]
    prefix = template.split("{")[0]
    assert r["reading"].startswith(prefix), prefix[:60]
    assert "NO PER-REGION MODEL IS" in r["reading"] or "PROPOSES NO FURTHER MODEL" in r["reading"]


def test_the_committed_result_carries_the_FRACTION_OF_VARIANCE_beside_better_or_worse():
    """DECISION (6): the fraction of the held-out variance is in the sentence, so a summary cannot drop it."""
    r = _result()
    assert "% of the held-out variance" in r["reading"] or "% OF THE VARIANCE" in r["reading"]
    d = r["descriptives"]
    var = r["held_out_variance"]
    assert (
        abs(d["arms"]["region_arm_percent_of_variance"] - 100.0 * r["measured"]["region_arm_delta_mse"] / var)
        < 1e-9
    )
    assert (
        abs(d["arms"]["global_arm_percent_of_variance"] - 100.0 * r["measured"]["global_arm_delta_mse"] / var)
        < 1e-9
    )


def test_the_committed_result_DELTA_is_exactly_the_difference_of_the_two_arms():
    r = _result()
    m = r["measured"]
    assert m["delta"] == m["region_arm_delta_mse"] - m["global_arm_delta_mse"]
    assert m["resamples"] == rs.BOOTSTRAP["resamples"] and m["bootstrap_seed"] == rs.BOOTSTRAP["seed"]
    assert m["permutations"] == rs.PERMUTATIONS["count"]
    assert len(m["per_region_scales"]) == r["design"]["admitted_count"]


def test_the_committed_result_global_scale_IS_M2hs_committed_one():
    """A check and not a value: the same fit on the same training genes must return +0.2145911381525326."""
    g = _result()["global_scale_check"]
    assert g["agrees"] is True
    assert g["m2h_committed"] == rs.GLOBAL_SCALE_COMMITTED
    assert g["abs_difference"] <= rs.M2H_TOLERANCE


def test_the_committed_result_per_chromosome_table_is_named_a_DESCRIPTIVE_and_decides_nothing():
    r = _result()
    d = r["descriptives"]
    assert len(d["per_admitted_chromosome"]) == r["design"]["admitted_count"]
    assert "EACH DECIDING NOTHING" in d["per_admitted_chromosome_decides_nothing"]
    assert d["p_against_the_stricter_divisor"]["threshold"] == 0.0167
    assert "THE VERDICT USES 0.05 AND ONLY 0.05" in d["p_against_the_stricter_divisor"]["what_it_is"]


def test_the_committed_result_PROPOSES_NO_MODEL_and_discloses_the_trial_run():
    r = _result()
    assert "NO FURTHER MODEL IS PROPOSED" in r["no_further_model_is_proposed"]
    assert r["trial_run_before_the_registration_was_committed"] == _run_module().TRIAL_RUN_DISCLOSURE
    assert r["alphagenome_requests"] == 0


def test_the_committed_result_manifest_is_complete_declares_its_readers_and_uses_relative_paths():
    m = _result()["result_manifest"]
    assert m["complete"] is True, m.get("problems")
    assert len(m["readers"]) == len(rs.CELLS) == 13
    for e in m["inputs"]:
        assert not e["path"].startswith("/"), e["path"]  # d72f959: input_entry paths RELATIVE
    assert m["code_cleanliness"]["own_code_is_committed"] is True
    assert m["code_cleanliness"]["foreign_uncommitted_code_on_the_counting_path"] == []
