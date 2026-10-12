# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalehet: is M2h's magnitude one global scale, or does it vary by region?

What is guarded here: that the grouping is by an address and not by a measured value; that every word
and function that fits, decides or reads a calibrated-form outcome is M2h's BY IDENTITY and not by
copy, so a word cannot drift; that the permutation null holds each group's gene count exactly and
moves a gene as one block; that the heterogeneity statistic is the registered Cochran form and the p
magcal's own form; that `decide_heterogeneity` is planted in ALL FOUR of its directions, including the
MIXED case that claims neither answer; that a planted regional slope reads HETEROGENEOUS and a planted
single slope reads CONSISTENT_WITH_ONE_SCALE, so the test is not vacuous; that no reading strengthens
the gain or calls a scale an explanation; and that the harness refuses when M2h's committed numbers do
not come back.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos.attribution import magcal as mc
from genomeos.attribution import scalehet as sh
from genomeos.attribution import scalesplit as ss

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "genomeos/attribution/scalehet.py"
HARNESS = ROOT / "scripts/scalehet_harness.py"
REGISTER = ROOT / "scripts/scalehet_register.py"
RUN = ROOT / "scripts/scalehet_run.py"
HARNESS_RESULT = ROOT / "data/results/scalehet_harness.json"
REGISTRATION = ROOT / "data/results/scalehet_registration.json"
RESULT = ROOT / "data/results/scalehet_result.json"
sys.path.insert(0, str(ROOT / "scripts"))


def rows_on(chrom: str, n: int, slope: float, start: int = 0) -> list[dict]:
    """n synthetic genes on one chromosome with D = slope * P, thirteen cells each."""
    out = []
    for i in range(n):
        p = [((i + j + start) % 7) - 3.0 for j in range(13)]
        out.append(
            {
                "gene": f"g{chrom}{i}",
                "symbol": f"S{chrom}{i}",
                "chrom": chrom,
                "n_rules": 1,
                "n_inhibiting": 0,
                "P": p,
                "D": [slope * x for x in p],
            }
        )
    return out


def statistics_of(rows: list[dict]) -> dict:
    """The observed statistics on rows, admitted by the registered gene floor."""
    dps, pps = sh.gene_sums(rows)
    counts = sh.gene_counts(rows)
    adm, _exc = sh.admitted_chromosomes(counts)
    return sh.observed_statistics(dps, pps, [r["chrom"] for r in rows], adm)


def run_on(rows: list[dict], count: int = 200, seed: int = 1) -> tuple[str, float, float]:
    """The whole registered procedure on synthetic rows: (outcome, p_Q, p_R)."""
    dps, pps = sh.gene_sums(rows)
    counts = sh.gene_counts(rows)
    adm, _exc = sh.admitted_chromosomes(counts)
    obs = statistics_of(rows)
    keep = [i for i, r in enumerate(rows) if r["chrom"] in adm]
    null = sh.permutation_null(
        [dps[i] for i in keep], [pps[i] for i in keep], obs["gene_counts"], count, seed
    )
    p_q = sh.permutation_p(obs["q"], null["q"])
    p_r = sh.permutation_p(obs["r"], null["r"])
    return sh.decide_heterogeneity(p_q, p_r, len(adm), null["degenerate"]), p_q, p_r


# ---- the grouping depends on no measured value ---------------------------------------------------
def test_the_grouping_is_the_chromosome_name_and_the_gene_counts_are_a_count():
    rows = rows_on("chr1", 150, 0.3) + rows_on("chr2", 120, 0.3) + rows_on("chrY", 9, 0.3)
    counts = sh.gene_counts(rows)
    assert counts == {"chr1": 150, "chr2": 120, "chrY": 9}
    # the count does not change when every measured value is replaced
    for r in rows:
        r["P"] = [0.0] * 13
        r["D"] = [99.0] * 13
    assert sh.gene_counts(rows) == counts


def test_the_gene_floor_admits_at_min_genes_and_excludes_below_it():
    counts = {"chrA": sh.MIN_GENES, "chrB": sh.MIN_GENES - 1, "chrC": sh.MIN_GENES + 1}
    adm, exc = sh.admitted_chromosomes(counts)
    assert adm == ["chrA", "chrC"]
    assert exc == ["chrB"]


def test_the_registered_floor_is_magcals_own_and_is_100():
    assert sh.MIN_GENES is mc.MIN_GENES
    assert sh.MIN_GENES == 100
    assert sh.MIN_CHROMOSOMES == 10


def test_a_chromosome_below_the_floor_enters_no_statistic():
    rows = rows_on("chr1", 150, 0.3) + rows_on("chr2", 120, 0.3) + rows_on("chrY", 9, 9.0)
    obs = statistics_of(rows)
    assert obs["groups"] == ["chr1", "chr2"]
    assert sum(obs["gene_counts"]) == 270  # chrY's nine genes are not in it
    # and chrY's wildly different slope has moved nothing
    same = statistics_of(rows_on("chr1", 150, 0.3) + rows_on("chr2", 120, 0.3))
    assert obs["q"] == pytest.approx(same["q"])
    assert obs["global_scale"] == pytest.approx(same["global_scale"])


# ---- every calibrated-form word and function is M2h's, by identity -------------------------------
def test_the_fit_the_interval_the_alpha_and_the_cells_are_magcals_objects_not_copies():
    assert sh.FORM_INTERVAL is mc.FORM_INTERVAL
    assert sh.FORM_INTERVAL_WORDS is mc.FORM_INTERVAL_WORDS
    assert sh.CONTROL_ALPHA is mc.CONTROL_ALPHA
    assert sh.CELLS is mc.CELLS
    assert sh.MEASURED is mc.MEASURED
    assert sh.TOLERANCE is mc.TOLERANCE
    assert sh.FORM_WORDS is mc.FORM_WORDS[sh.FORM]
    assert sh.OUTCOMES["HARNESS_DEFECT"] is mc.OUTCOMES["HARNESS_DEFECT"]
    assert sh.M2H_TOLERANCE is ss.M2H_TOLERANCE
    assert sh.quantiles([1.0, 2.0, 3.0]) == mc.quantiles([1.0, 2.0, 3.0])


def test_m2hs_committed_numbers_are_scalesplits_and_are_not_retyped():
    assert dict(ss.M2H_COMMITTED) == sh.M2H_COMMITTED
    assert dict(ss.SEVEN_H_COMMITTED) == sh.SEVEN_H_COMMITTED
    assert sh.M2H_COMMITTED["scale"] == pytest.approx(0.2145911381525326)
    assert sh.M2H_COMMITTED["held_out_delta_mse"] == pytest.approx(0.0433912108101358)


def test_this_lane_edits_no_peers_file():
    src = MODULE.read_text()
    assert "def fit_scale" not in src  # magcal's fit is imported, never reimplemented
    assert "def delta_mse" not in src
    assert "def bootstrap_" not in src
    assert "def permuted_control" not in src
    for owned in ("magcal.py", "scalesplit.py", "multrna.py", "combine2.py", "recovered.py"):
        assert f"open('genomeos/attribution/{owned}'" not in src


def test_the_per_chromosome_scale_path_is_magcals_fit_on_the_same_rows():
    rows = rows_on("chr1", 150, 0.37) + rows_on("chr2", 130, 0.11, start=3)
    dps, pps = sh.gene_sums(rows)
    for chrom in ("chr1", "chr2"):
        idx = [i for i, r in enumerate(rows) if r["chrom"] == chrom]
        mine = sh.scale_from_sums(sum(dps[i] for i in idx), sum(pps[i] for i in idx))
        assert mine == pytest.approx(mc.fit_scale([rows[i] for i in idx]), abs=1e-12)
    assert sh.scale_from_sums(sum(dps), sum(pps)) == pytest.approx(mc.fit_scale(rows), abs=1e-12)


def test_an_undefined_fit_is_none_and_is_not_zero():
    assert sh.scale_from_sums(0.0, 0.0) is None
    assert sh.scale_from_sums(1.0, 0.0) is None


# ---- the statistic is the registered form -------------------------------------------------------
def test_q_is_the_precision_weighted_sum_of_squared_deviations():
    q = sh.heterogeneity_q([0.1, 0.3], [2.0, 5.0], 0.2)
    assert q == pytest.approx(2.0 * 0.01 + 5.0 * 0.01)


def test_q_is_zero_exactly_when_every_group_sits_at_the_global_scale():
    assert sh.heterogeneity_q([0.2, 0.2, 0.2], [1.0, 9.0, 4.0], 0.2) == 0.0


def test_q_weights_a_precise_group_more_than_a_noisy_one_at_the_same_deviation():
    precise = sh.heterogeneity_q([0.3], [100.0], 0.2)
    noisy = sh.heterogeneity_q([0.3], [1.0], 0.2)
    assert precise > noisy


def test_r_is_the_max_minus_min_of_the_group_scales():
    assert sh.range_stat([0.1232, 0.2176, 0.2874]) == pytest.approx(0.2874 - 0.1232)


def test_the_p_is_magcals_one_sided_form_and_never_zero():
    assert sh.permutation_p(5.0, [1.0] * 1000) == pytest.approx(1 / 1001)
    assert sh.permutation_p(5.0, [9.0] * 1000) == pytest.approx(1.0)
    assert sh.permutation_p(5.0, [5.0] * 1000) == pytest.approx(1.0)  # >= is a tie, counted
    assert sh.permutation_p(1.0, [0.0] * 10) > 0.0


def test_the_registered_permutation_count_is_magcals_and_the_seed_is_this_lanes():
    assert sh.PERMUTATIONS["count"] == mc.PERMUTATIONS["count"] == 1000
    assert sh.PERMUTATIONS["seed"] == 20261012
    assert sh.PERMUTATIONS["seed"] != mc.PERMUTATIONS["seed"]


# ---- the null holds the gene counts and moves a gene as one block -------------------------------
def test_the_null_holds_every_groups_gene_count_exactly():
    dps = [float(i + 1) for i in range(60)]
    pps = [1.0 + (i % 3) for i in range(60)]
    sizes = [30, 20, 10]
    null = sh.permutation_null(dps, pps, sizes, 5, seed=7)
    assert null["degenerate"] is None
    assert len(null["per_group"]) == 3
    assert all(len(g) == 5 for g in null["per_group"])
    # every draw is a partition of the same 60 genes into blocks of 30, 20 and 10, so the TOTAL sums
    # -- and therefore the global scale -- are invariant under the null, which is the point of it
    assert len({round(x, 12) for x in null["global_scales"]}) == 1
    assert null["global_scales"][0] == pytest.approx(sum(dps) / sum(pps))


def test_the_null_refuses_sizes_that_do_not_sum_to_the_genes_given():
    with pytest.raises(ValueError, match="group sizes sum to"):
        sh.permutation_null([1.0, 2.0], [1.0, 1.0], [1, 2], 2, seed=1)


def test_the_null_is_reproducible_at_its_seed_and_differs_at_another():
    args = ([float(i) for i in range(40)], [1.0] * 40, [20, 20], 20)
    a = sh.permutation_null(*args, seed=11)["q"]
    b = sh.permutation_null(*args, seed=11)["q"]
    c = sh.permutation_null(*args, seed=12)["q"]
    assert a == b
    assert a != c


def test_a_gene_travels_as_one_block_so_its_two_sums_are_never_separated():
    rows = rows_on("chr1", 4, 0.5)
    dps, pps = sh.gene_sums(rows)
    for r, dp, pp in zip(rows, dps, pps, strict=True):
        assert dp == pytest.approx(sum(p * d for p, d in zip(r["P"], r["D"], strict=True)))
        assert pp == pytest.approx(sum(p * p for p in r["P"]))
    # the per-gene sums are all the permutation sees, so a gene's thirteen pairs cannot be split
    assert len(dps) == len(rows) == 4


def test_a_group_whose_sum_of_p_squared_is_zero_is_degenerate_not_a_scale_of_zero():
    obs = sh.observed_statistics([0.0, 1.0], [0.0, 2.0], ["chrA", "chrB"], ["chrA", "chrB"])
    assert obs["degenerate"]
    assert "undefined" in obs["degenerate"]


# ---- the decision table, planted in all four directions -----------------------------------------
def test_both_above_alpha_reads_consistent_with_one_scale():
    assert sh.decide_heterogeneity(0.9, 0.7, 23) == "CONSISTENT_WITH_ONE_SCALE"
    assert sh.decide_heterogeneity(0.0501, 0.0501, 23) == "CONSISTENT_WITH_ONE_SCALE"


def test_both_at_or_below_alpha_reads_heterogeneous():
    assert sh.decide_heterogeneity(0.000999, 0.001, 23) == "HETEROGENEOUS"
    assert sh.decide_heterogeneity(0.05, 0.05, 23) == "HETEROGENEOUS"  # <= is a rejection


def test_the_two_statistics_on_opposite_sides_read_mixed_and_claim_neither():
    assert sh.decide_heterogeneity(0.01, 0.9, 23) == "MIXED"
    assert sh.decide_heterogeneity(0.9, 0.01, 23) == "MIXED"
    assert "NEITHER ANSWER IS CLAIMED" in sh.OUTCOMES["MIXED"]


def test_too_few_admitted_chromosomes_reads_inconclusive_whatever_the_p_values():
    assert sh.decide_heterogeneity(0.000999, 0.000999, sh.MIN_CHROMOSOMES - 1) == "INCONCLUSIVE"
    assert sh.decide_heterogeneity(0.9, 0.9, 0) == "INCONCLUSIVE"
    assert sh.decide_heterogeneity(0.000999, 0.000999, sh.MIN_CHROMOSOMES) == "HETEROGENEOUS"


def test_a_degenerate_null_reads_inconclusive_and_cannot_be_overridden_by_a_small_p():
    assert sh.decide_heterogeneity(0.000999, 0.000999, 23, degenerate="a group's P^2 sum is 0") == (
        "INCONCLUSIVE"
    )


def test_r_alone_can_never_declare_either_answer():
    # whatever R says, with Q above alpha the verdict is never HETEROGENEOUS
    for p_r in (0.000999, 0.01, 0.049, 0.05):
        assert sh.decide_heterogeneity(0.5, p_r, 23) == "MIXED"
    # and with Q below alpha it is never CONSISTENT
    for p_r in (0.06, 0.5, 1.0):
        assert sh.decide_heterogeneity(0.001, p_r, 23) == "MIXED"


def test_the_verdict_uses_alpha_0_05_and_not_the_stricter_bonferroni_divisor():
    assert sh.CONTROL_ALPHA == 0.05
    # a p between 0.0167 and 0.05 must REJECT, which is the lax direction the registration argues for
    assert sh.decide_heterogeneity(0.03, 0.03, 23) == "HETEROGENEOUS"


# ---- the test is not vacuous: both outcomes are reachable on planted data ------------------------
def test_PLANTED_one_slope_everywhere_reads_consistent_with_one_scale():
    rows = (
        rows_on("chr1", 400, 0.21)
        + rows_on("chr2", 300, 0.21, start=1)
        + rows_on("chr3", 250, 0.21, start=2)
        + rows_on("chr4", 200, 0.21, start=3)
        + rows_on("chr5", 150, 0.21, start=4)
        + rows_on("chr6", 120, 0.21, start=5)
        + rows_on("chr7", 110, 0.21, start=6)
        + rows_on("chr8", 105, 0.21, start=7)
        + rows_on("chr9", 102, 0.21, start=8)
        + rows_on("chr10", 101, 0.21, start=9)
    )
    outcome, p_q, p_r = run_on(rows)
    assert outcome == "CONSISTENT_WITH_ONE_SCALE"
    assert p_q > sh.CONTROL_ALPHA and p_r > sh.CONTROL_ALPHA


def test_PLANTED_a_different_slope_on_some_chromosomes_reads_heterogeneous():
    rows = (
        rows_on("chr1", 400, 0.21)
        + rows_on("chr2", 300, 0.21, start=1)
        + rows_on("chr3", 250, 0.21, start=2)
        + rows_on("chr4", 200, 0.21, start=3)
        + rows_on("chr5", 150, 0.21, start=4)
        + rows_on("chr6", 120, 0.60, start=5)
        + rows_on("chr7", 110, 0.60, start=6)
        + rows_on("chr8", 105, 0.60, start=7)
        + rows_on("chr9", 102, 0.60, start=8)
        + rows_on("chr10", 101, 0.60, start=9)
    )
    outcome, p_q, _p_r = run_on(rows)
    assert outcome == "HETEROGENEOUS"
    assert p_q <= sh.CONTROL_ALPHA


def test_PLANTED_noise_at_a_small_gene_count_does_not_by_itself_read_heterogeneous():
    """The whole point: a WIDE SPREAD at small counts is what the null itself produces.

    Ten groups, one slope, but the per-group gene counts are small, so the fitted per-group scales
    spread visibly. The spread must not be read as regional variation.
    """
    import random

    rng = random.Random(4)
    rows: list[dict] = []
    for k in range(10):
        chrom = f"chr{k + 1}"
        for i in range(110):
            p = [rng.gauss(0, 1) for _ in range(13)]
            rows.append(
                {
                    "gene": f"g{chrom}{i}",
                    "symbol": f"S{chrom}{i}",
                    "chrom": chrom,
                    "n_rules": 1,
                    "n_inhibiting": 0,
                    "P": p,
                    "D": [0.21 * x + rng.gauss(0, 1.5) for x in p],
                }
            )
    obs = statistics_of(rows)
    assert max(obs["scales"]) - min(obs["scales"]) > 0.05  # the spread IS visible
    outcome, p_q, _p_r = run_on(rows, count=300, seed=5)
    assert outcome == "CONSISTENT_WITH_ONE_SCALE"
    assert p_q > sh.CONTROL_ALPHA


# ---- the readings: every registered sentence fills, and none strengthens -------------------------
NUMBERS = {
    "n_admitted": 23,
    "a_global": 0.2146,
    "q_obs": 1.234,
    "q_med": 1.1,
    "q_lo": 0.4,
    "q_hi": 2.9,
    "p_q": 0.37,
    "r_obs": 0.1642,
    "a_min": 0.1232,
    "a_max": 0.2874,
    "chrom_min": "chr13",
    "chrom_max": "chr14",
    "r_med": 0.16,
    "r_lo": 0.09,
    "r_hi": 0.26,
    "p_r": 0.41,
    "m2h_scale": 0.2146,
    "count": 1000,
    "seed": 20261012,
}


@pytest.mark.parametrize("kind", ["CONSISTENT_WITH_ONE_SCALE", "HETEROGENEOUS", "MIXED"])
def test_every_registered_reading_fills_from_the_run_numbers(kind):
    text = sh.reading(kind, **NUMBERS)
    assert "{" not in text and "}" not in text
    assert len(text) > 400


def test_the_inconclusive_reading_names_its_reason():
    text = sh.reading("INCONCLUSIVE", what="only 4 chromosomes cleared the gene floor")
    assert "only 4 chromosomes cleared the gene floor" in text
    assert "{" not in text


def test_the_full_reading_always_carries_the_prior_exposure_clause():
    for kind in ("CONSISTENT_WITH_ONE_SCALE", "HETEROGENEOUS", "MIXED"):
        assert sh.full_reading(kind, **NUMBERS).endswith(sh.PRIOR_EXPOSURE_CLAUSE)
    assert sh.full_reading("INCONCLUSIVE", what="x").endswith(sh.PRIOR_EXPOSURE_CLAUSE)


def test_the_prior_exposure_clause_says_the_pattern_was_already_looked_at():
    c = sh.PRIOR_EXPOSURE_CLAUSE
    assert "ALREADY" in c and "POST HOC" in c
    assert "weaker than a test of one that has not" in c
    for lane in ("items 7, 7d and", "capped-sum census", "M2h", "M2i", "scalesplit", "recovered"):
        assert lane in c


def test_the_gain_is_never_strengthened_in_any_reading():
    for kind in ("CONSISTENT_WITH_ONE_SCALE", "HETEROGENEOUS"):
        text = sh.reading(kind, **NUMBERS)
        assert "about 2% of the variance" in text
        for forbidden in ("explains expression", "explains the expression", "most of the variance"):
            assert forbidden not in text.lower()


def test_no_reading_calls_a_scale_an_explanation_of_expression():
    blob = " ".join(sh.OUTCOMES.values()) + " " + " ".join(sh.LIMITATIONS)
    assert "MAGNITUDE CORRECTION" in blob
    assert not re.search(r"scale\s+explains", blob, re.I)


def test_the_consistent_reading_states_the_answer_plainly_and_as_a_full_answer():
    text = sh.reading("CONSISTENT_WITH_ONE_SCALE", **NUMBERS)
    assert "CONSISTENT WITH A SINGLE GLOBAL SCALE" in text
    assert "THE MAGNITUDE IS ONE CONSTANT" in text
    assert "failure to reject is not a proof" in text
    assert "FULL ANSWER" in sh.QUESTION


def test_the_heterogeneous_reading_proposes_no_regional_model():
    text = sh.reading("HETEROGENEOUS", **NUMBERS)
    assert "no regional model is proposed" in text
    assert "no per-region scale is fitted" in text
    assert "separate registration and a separate decision" in text
    assert "no cause is named" in text


# ---- the registration binds the run -------------------------------------------------------------
def test_every_bound_field_names_a_real_constant_and_registered_fields_returns_them_all():
    fields = sh.registered_fields()
    assert set(fields) == set(sh.BOUND)
    for field, name in sh.BOUND.items():
        assert hasattr(sh, name), f"{field} names a missing constant {name}"
        assert fields[field] is not None


def test_the_binding_set_covers_every_decision_the_run_could_otherwise_move():
    for field in (
        "question",
        "min_genes",
        "min_genes_rule",
        "min_chromosomes",
        "statistic",
        "null",
        "permutations",
        "p_value",
        "alpha",
        "interval",
        "prior_exposure",
        "falsifier",
        "decision",
        "outcomes",
        "stop_rule",
    ):
        assert field in sh.BOUND


def test_the_registration_states_the_exposure_the_floor_the_interval_and_the_falsifier():
    reg = sh.registered_fields()
    assert "MIN_GENES = 100" in reg["min_genes_rule"]
    assert "13 genes" in reg["min_genes_rule"]  # chrY, named with its count
    assert "98.33%" in reg["interval"]
    assert "IS FALSIFIED by p_Q <= 0.05" in reg["falsifier"]
    assert "INCONCLUSIVE IS" in reg["falsifier"]
    assert any("ALREADY SEEN" in x for x in reg["prior_exposure"])
    assert "0.05" in reg["alpha"] and "0.0167" in reg["alpha"]


def test_the_alpha_argues_the_lax_direction_rather_than_the_flattering_one():
    assert "LAXER" in sh.ALPHA
    assert "HARDER TO OBTAIN" in sh.ALPHA
    assert "NO Bonferroni correction is applied" in sh.ALPHA


def test_the_null_is_described_as_holding_the_gene_counts_and_as_exchangeability():
    blob = " ".join(sh.NULL)
    assert "EACH CHROMOSOME'S GENE COUNT EXACTLY AT ITS OBSERVED VALUE" in blob
    assert "EXCHANGEABLE" in blob
    assert "20261012" in blob


def test_the_limitations_name_the_coarse_unit_the_exchangeability_and_the_shared_elements():
    blob = " ".join(sh.LIMITATIONS)
    assert "THE REGION UNIT IS THE CHROMOSOME" in blob
    assert "EXCHANGEABLE BETWEEN CHROMOSOMES" in blob
    assert "SHARE A COMPILED ELEMENT" in blob
    assert "anti-conservative" in blob


# ---- the harness refuses, and the run cannot precede the registration ---------------------------
def test_the_harness_reproduces_m2hs_committed_numbers():
    tp.must_be_present(HARNESS_RESULT, was="the harness has not been run yet")
    d = json.loads(HARNESS_RESULT.read_text())
    assert d["outcome"] == "HARNESS_REPRODUCES"
    rep = d["m2h_reproduction"]
    assert rep["reproduces"] is True
    for k, v in rep["committed"].items():
        got = rep["recomputed"][k]
        if isinstance(v, str):
            assert got == v
        elif isinstance(v, list):
            for g, w in zip(got, v, strict=True):
                assert abs(g - w) <= sh.M2H_TOLERANCE
        else:
            assert abs(float(got) - float(v)) <= sh.M2H_TOLERANCE


def test_the_harness_reads_no_statistic_of_this_lane():
    tp.must_be_present(HARNESS_RESULT, was="the harness has not been run yet")
    d = json.loads(HARNESS_RESULT.read_text())
    blob = json.dumps(d)
    for forbidden in ("heterogeneity_q", '"q_obs"', '"p_q"', '"per_chromosome_scales"'):
        assert forbidden not in blob
    assert "no Q, no range statistic" in d["nothing_of_this_lanes_statistic_was_read"]


def test_the_harness_refuses_when_m2hs_numbers_do_not_come_back():
    import scalehet_harness as hn

    rows = hn.load_rows()
    rows[0]["D"] = [x + 7.0 for x in rows[0]["D"]]  # one gene's measurement moved
    _out, problems = hn.reproduce_m2h(rows)
    assert problems
    assert any("M2h" in p for p in problems)


def test_the_harness_refuses_a_table_that_is_not_7hs():
    import scalehet_harness as hn

    rows = hn.load_rows()[:-1]
    _out, problems = hn.table_checks(rows)
    assert problems
    assert any("genes" in p for p in problems)


def test_the_harness_declares_its_thirteen_readers():
    import scalehet_harness as hn

    assert list(sh.CELLS) == hn.READERS
    assert len(hn.READERS) == 13


def test_the_harness_declares_its_inputs_relative_to_the_repository():
    import scalehet_harness as hn

    for p in hn.SOURCE_FILES:
        assert not p.startswith("/")
        assert p.startswith("data/results/")
        assert (ROOT / p).exists()


def test_the_registration_is_committed_and_binds_the_module():
    tp.must_be_present(REGISTRATION, was="the registration has not been written yet")
    reg = json.loads(REGISTRATION.read_text())
    want = sh.registered_fields()
    for field, v in want.items():
        assert reg[field] == v, f"the registration's {field} differs from the module"
    assert reg["alphagenome_requests"] == 0


def test_the_result_is_the_registered_decision_recomputed_and_not_a_sentence_typed_by_hand():
    tp.must_be_present(RESULT, was="the run has not been done yet")
    d = json.loads(RESULT.read_text())
    outcome = sh.decide_heterogeneity(d["p_q"], d["p_r"], d["admitted_count"], d.get("degenerate"))
    assert d["outcome"] == outcome
    assert d["reading"] == sh.full_reading(outcome, **d["reading_numbers"])
    assert d["registration"] == "data/results/scalehet_registration.json"
    assert d["permutations"] == sh.PERMUTATIONS
    assert d["alphagenome_requests"] == 0


def test_the_results_p_values_are_the_registered_form_recomputed_from_its_own_nulls():
    tp.must_be_present(RESULT, was="the run has not been done yet")
    d = json.loads(RESULT.read_text())
    assert d["p_q"] == pytest.approx(sh.permutation_p(d["q_observed"], d["null_q"]))
    assert d["p_r"] == pytest.approx(sh.permutation_p(d["r_observed"], d["null_r"]))
    assert len(d["null_q"]) == len(d["null_r"]) == sh.PERMUTATIONS["count"]


def test_the_result_reproduces_scalesplits_post_hoc_spread_rather_than_quoting_it():
    tp.must_be_present(RESULT, was="the run has not been done yet")
    d = json.loads(RESULT.read_text())
    per = d["per_chromosome"]
    scales = [c["scale"] for c in per if c["admitted"]]
    assert d["q_observed"] == pytest.approx(
        sh.heterogeneity_q(scales, [c["pp"] for c in per if c["admitted"]], d["global_scale"])
    )
    assert d["r_observed"] == pytest.approx(sh.range_stat(scales))


def test_the_scripts_exist_and_the_module_carries_the_licence_header():
    for p in (MODULE, HARNESS, REGISTER, RUN):
        assert p.exists(), f"{p} is missing"
        assert p.read_text().startswith("# SPDX-License-Identifier: AGPL-3.0-or-later")
