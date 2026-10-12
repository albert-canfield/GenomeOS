# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-scalesplit: does M2h's one scale survive a CHROMOSOME-held-out split?

What is guarded here: that the chromosome partition is a partition and depends on no measured value;
that every word and every function that fits, bootstraps, decides or reads an outcome is M2h's, by
IDENTITY and not by copy, so a word cannot drift; that the transport verdict is the registered function
of the two folds' outcomes, planted in all three directions, and that it cannot bury a failure; that
the magnitude comparison against M2h never says 'stronger' or 'weaker' on overlapping intervals; and
that the harness refuses when M2h's committed numbers do not come back.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest

from genomeos.attribution import magcal as mc
from genomeos.attribution import scalesplit as ss

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "scripts/scalesplit_run.py"
HARNESS = ROOT / "scripts/scalesplit_harness.py"
REGISTER = ROOT / "scripts/scalesplit_register.py"
MODULE = ROOT / "genomeos/attribution/scalesplit.py"
MR_RESULT = ROOT / "data/results/multrna_result.json"
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


# ---- the split is a partition, and depends on no measured value ---------------------------------
def test_the_two_chromosome_sets_are_disjoint_and_cover_24_chromosomes():
    odd, even = set(ss.ODD_CHROMOSOMES), set(ss.EVEN_AND_SEX_CHROMOSOMES)
    assert not odd & even
    assert len(odd) == 11 and len(even) == 13
    assert len(odd | even) == 24
    assert odd == {f"chr{i}" for i in range(1, 22, 2)}
    assert even == {f"chr{i}" for i in range(2, 23, 2)} | {"chrX", "chrY"}


def test_the_two_folds_are_complementary_and_every_gene_is_held_out_exactly_once():
    rows = rows_on("chr1", 3, 0.5) + rows_on("chr2", 4, 0.5) + rows_on("chrY", 2, 0.5)
    ta, ha = ss.split_rows_by_chromosome(rows, "odd_held")
    tb, hb = ss.split_rows_by_chromosome(rows, "even_and_sex_held")
    assert [r["symbol"] for r in ha] == [r["symbol"] for r in tb]
    assert [r["symbol"] for r in ta] == [r["symbol"] for r in hb]
    assert len(ta) + len(ha) == len(rows)
    held_once = [r["symbol"] for r in ha] + [r["symbol"] for r in hb]
    assert sorted(held_once) == sorted(r["symbol"] for r in rows)


def test_the_partition_depends_on_no_measured_value():
    rows = rows_on("chr1", 5, 0.5) + rows_on("chr2", 5, 0.5)
    before = [r["symbol"] for r in ss.split_rows_by_chromosome(rows, "odd_held")[1]]
    for r in rows:  # destroy every measured value; the partition must not move
        r["P"] = [0.0] * 13
        r["D"] = [99.0] * 13
    after = [r["symbol"] for r in ss.split_rows_by_chromosome(rows, "odd_held")[1]]
    assert before == after


def test_a_chromosome_named_by_neither_set_is_reported_unnamed_and_not_silently_trained_on():
    assert ss.chromosome_set("chr1") == "odd"
    assert ss.chromosome_set("chrX") == "even_and_sex"
    assert ss.chromosome_set("chrM") == "unnamed"
    rows = rows_on("chrM", 3, 0.5)
    for fold in ss.FOLDS:
        train, held = ss.split_rows_by_chromosome(rows, fold)
        assert not train and not held


def test_PLANTED_an_unnamed_chromosome_in_the_table_is_a_harness_defect():
    from scalesplit_harness import chromosome_checks

    _, problems = chromosome_checks(rows_on("chr1", 2, 0.5) + rows_on("chrM", 2, 0.5))
    assert any("named by neither set" in p for p in problems)


# ---- M2h's words and functions are carried by IDENTITY, so they cannot drift ---------------------
@pytest.mark.parametrize(
    "name",
    [
        "WHICH_GOVERNS",
        "NOT_VACUOUS",
        "OUTCOMES",
        "SCALE_CLAUSES",
        "MIN_GENES",
        "CONTROL_ALPHA",
        "PERMUTATIONS",
        "FORM_INTERVAL",
        "FORM_INTERVAL_WORDS",
        "TOLERANCE",
        "CELLS",
        "FIT",
    ],
)
def test_the_registered_vocabulary_is_magcals_own_object(name):
    assert getattr(ss, name) is getattr(mc, name)


def test_this_lane_reimplements_none_of_magcals_machinery():
    src = MODULE.read_text() + RUN.read_text() + HARNESS.read_text()
    for fn in (
        "fit_scale",
        "bootstrap_scale",
        "bootstrap_delta",
        "permuted_control",
        "delta_mse",
        "decide",
        "scale_is_negative",
        "reading",
        "scale_clause",
        "quantiles",
    ):
        assert not re.search(rf"^\s*def {fn}\b", src, re.M), f"{fn} is reimplemented, not imported"
        assert f"mc.{fn}" in src, f"{fn} is not called through magcal"


def test_this_lane_writes_no_file_of_magcals_multrnas_or_combine2s():
    src = MODULE.read_text() + RUN.read_text() + HARNESS.read_text() + REGISTER.read_text()
    for path in ("magcal.py", "multrna.py", "combine2.py", "magcal_result", "multrna_result"):
        assert not re.search(rf"write_text\([^)]*{path}", src)
    assert 'RESULT = "scalesplit_result"' in RUN.read_text()
    assert 'RESULT = "scalesplit_harness"' in HARNESS.read_text()


def test_the_scale_is_refitted_on_the_training_chromosomes_and_not_taken_from_M2h():
    assert "REFITTED ON THE TRAINING CHROMOSOMES" in ss.REFITTED
    assert "NOT taken as M2h's constant" in ss.REFITTED
    run = RUN.read_text()
    assert "mc.fit_scale(train)" in run
    assert "mc.bootstrap_scale(train," in run
    # M2h's committed scale is quoted as context only, never used as a fold's scale
    assert not re.search(r"scale\s*=\s*ss\.M2H_COMMITTED", run)


# ---- the four-way table is M2h's, applied as written --------------------------------------------
@pytest.mark.parametrize(
    "lo,hi,alo,ahi,want",
    [
        (0.01, 0.05, 0.2, 0.23, "BETTER"),
        (-0.05, -0.01, 0.2, 0.23, "WORSE"),
        (-0.01, 0.01, -0.05, 0.05, "CALIBRATION_REMOVES"),
        (-0.01, 0.01, 0.2, 0.23, "INCONCLUSIVE"),
        (0.0, 0.05, 0.2, 0.23, "INCONCLUSIVE"),
    ],
)
def test_the_held_out_interval_decides_and_the_scale_interval_only_splits_the_straddle(
    lo, hi, alo, ahi, want
):
    assert mc.decide(lo, hi, alo, ahi) == want


def test_calibration_removes_is_registered_as_a_negative_and_not_as_no_difference():
    assert "negative result and is not 'no difference detected'" in ss.OUTCOMES["CALIBRATION_REMOVES"]
    assert "calibration removes the rules" in ss.FALSIFIER


# ---- the transport verdict: planted in all three directions -------------------------------------
@pytest.mark.parametrize(
    "a,b,want",
    [
        ("BETTER", "BETTER", "SURVIVES_BOTH_FOLDS"),
        ("BETTER", "WORSE", "SURVIVES_ONE_FOLD"),
        ("WORSE", "BETTER", "SURVIVES_ONE_FOLD"),
        ("BETTER", "INCONCLUSIVE", "SURVIVES_ONE_FOLD"),
        ("BETTER", "CALIBRATION_REMOVES", "SURVIVES_ONE_FOLD"),
        ("WORSE", "WORSE", "DOES_NOT_SURVIVE"),
        ("INCONCLUSIVE", "CALIBRATION_REMOVES", "DOES_NOT_SURVIVE"),
        ("NOT_TESTABLE", "NOT_FITTABLE", "DOES_NOT_SURVIVE"),
    ],
)
def test_the_transport_verdict_is_the_registered_function_of_the_two_folds(a, b, want):
    assert ss.transport_verdict({"odd_held": a, "even_and_sex_held": b}) == want


def test_a_failure_to_transport_is_the_FIRST_thing_the_reading_says():
    for verdict in ("SURVIVES_ONE_FOLD", "DOES_NOT_SURVIVE"):
        text = ss.registered_reading(verdict, ["fold one reads", "fold two reads"], ["a finding"])
        assert text.startswith(ss.TRANSPORT[verdict])
        assert "DOES NOT" in ss.TRANSPORT[verdict]
        assert "fold one reads; fold two reads" in text
        assert text.endswith(ss.PRIOR_EXPOSURE_CLAUSE)


def test_the_two_negative_transport_sentences_are_not_softened():
    assert "DOES NOT SURVIVE" in ss.TRANSPORT["DOES_NOT_SURVIVE"]
    assert "LEAK" in ss.TRANSPORT["DOES_NOT_SURVIVE"]
    assert "not softened" in ss.TRANSPORT["DOES_NOT_SURVIVE"]
    assert "DOES NOT TRANSPORT CLEANLY" in ss.TRANSPORT["SURVIVES_ONE_FOLD"]
    assert "stated as one" in ss.TRANSPORT["SURVIVES_ONE_FOLD"]


# ---- the magnitude comparison decides nothing and never strengthens -----------------------------
def test_overlapping_intervals_read_consistent_with_the_same_effect_and_never_stronger():
    s = ss.magnitude_finding(0.030, 0.050, 0.0379, 0.04894, "x")
    assert "CONSISTENT WITH THE SAME EFFECT" in s
    rest = s.replace("neither 'stronger' nor 'weaker' is said", "")
    for word in ("stronger", "weaker", "strengthen", "weaken"):
        assert word not in rest


def test_a_fold_interval_wholly_below_M2hs_reads_smaller_and_wholly_above_reads_larger():
    assert "SMALLER" in ss.magnitude_finding(0.001, 0.010, 0.0379, 0.04894, "x")
    assert "LARGER" in ss.magnitude_finding(0.10, 0.20, 0.0379, 0.04894, "x")
    assert ss.overlaps((0.0, 1.0), (1.0, 2.0)) is True
    assert ss.overlaps((0.0, 1.0), (1.1, 2.0)) is False


def test_no_test_of_the_difference_against_M2h_is_registered():
    assert "no test of the difference is registered and none will be computed" in ss.MAGNITUDE_COMPARISON
    assert "same 16,509 genes repartitioned" in ss.MAGNITUDE_COMPARISON


# ---- the harness: M2h's committed numbers must come back ----------------------------------------
def test_PLANTED_a_wrong_committed_number_makes_the_m2h_reproduction_fail(monkeypatch):
    import scalesplit_harness as h

    rows = rows_on("chr1", 200, 0.5) + rows_on("chr2", 200, 0.5)
    monkeypatch.setattr(ss, "M2H_COMMITTED", {**ss.M2H_COMMITTED, "scale": 999.0}, raising=True)
    _, problems = h.reproduce_m2h(rows)
    assert any("M2h scale" in p for p in problems)


def test_the_harness_states_that_it_reproduces_M2hs_numbers_before_anything_is_split():
    assert "M2h'S COMMITTED GENE-LEVEL NUMBERS REPRODUCE" in ss.HARNESS
    assert "stops the run (HARNESS_DEFECT) and nothing is split by chromosome" in ss.HARNESS
    assert "M2h's split is NOT re-run as this lane's endpoint" in ss.WHICH_SPLIT_IS_WHOSE


# ---- what the registration must say in advance --------------------------------------------------
def test_every_registered_field_round_trips_through_BOUND():
    fields = ss.registered_fields()
    assert set(fields) == set(ss.BOUND)
    for field, name in ss.BOUND.items():
        want = getattr(ss, name)
        assert fields[field] == (list(want) if isinstance(want, tuple) else want)
    for key in (
        "split",
        "refitted",
        "interval",
        "prior_exposure",
        "falsifier",
        "transport",
        "transport_rule",
        "magnitude_comparison",
        "post_hoc",
    ):
        assert key in fields


def test_the_prior_exposure_names_every_item_that_read_this_RNA():
    for item in ("item 7", "7d", "7h", "capped-sum", "M2h", "M2i"):
        assert item in ss.PRIOR_EXPOSURE
    assert "RE-ANALYSIS of a measurement already read, not a fresh test" in ss.PRIOR_EXPOSURE
    assert "M2i" in ss.PRIOR_EXPOSURE_CLAUSE and "RE-ANALYSIS" in ss.PRIOR_EXPOSURE_CLAUSE


def test_the_interval_and_its_correction_are_decided_in_the_registration():
    assert ss.FORM_INTERVAL_WORDS == "98.33% CI"
    assert "DECIDED NOW AND NOT AFTER" in ss.INTERVAL
    assert "97.5%" in ss.INTERVAL and "WIDER" in ss.INTERVAL
    assert "never as a Bonferroni over two" in ss.INTERVAL


def test_the_falsifier_and_inconclusive_are_the_supervisors_four_way_table():
    assert "upper bound < 0, outcome WORSE" in ss.FALSIFIER
    assert "INCONCLUSIVE when the scale interval excludes 0" in ss.FALSIFIER
    assert "a live scale with a straddling held-out interval" in ss.FALSIFIER


def test_the_post_hoc_decides_nothing():
    assert "DECIDING NOTHING" in ss.POST_HOC
    assert "no interval is placed on a per-chromosome scale" in ss.POST_HOC
    assert "cannot change a verdict" in ss.STOP_RULE


def test_the_wording_rules_of_this_series_are_kept():
    joined = " ".join(ss.OUTCOMES.values()) + " ".join(ss.TRANSPORT.values())
    assert "explains expression" not in joined
    assert "closer to measured RNA than" in joined
    assert "further from measured RNA than" in joined
    assert "about 2% of the variance" in ss.TRANSPORT["SURVIVES_BOTH_FOLDS"]
    assert "about 2% of the variance" in ss.QUESTION


# ---- the real committed table, structurally ------------------------------------------------------
def test_the_committed_table_carries_a_chromosome_for_every_gene_and_both_sets_clear_MIN_GENES():
    # no skip: 7h's result is a tracked artefact, so its absence is a failure and not a skip
    assert MR_RESULT.exists(), f"{MR_RESULT} is tracked and must be present"
    tab = json.loads(MR_RESULT.read_text())["genes_table"]
    ix = {c: i for i, c in enumerate(tab["columns"])}
    assert "chrom" in ix
    chroms = [r[ix["chrom"]] for r in tab["rows"]]
    assert len(chroms) == ss.SEVEN_H_COMMITTED["genes"]
    assert all(c for c in chroms)
    n_odd = sum(1 for c in chroms if c in ss.ODD_CHROMOSOMES)
    n_even = sum(1 for c in chroms if c in ss.EVEN_AND_SEX_CHROMOSOMES)
    assert n_odd + n_even == len(chroms)
    assert n_odd >= ss.MIN_GENES and n_even >= ss.MIN_GENES
