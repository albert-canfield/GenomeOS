# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fold over the all-elements summaries: complete chromosomes counted, partial ones listed."""

import importlib.util
import json
from pathlib import Path

spec = importlib.util.spec_from_file_location("agw", Path("scripts/enhancer_targets_all_genome_wide.py"))
agw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(agw)


def _res(chrom, total, scored, complete, named, agree, other, outside, strong):
    return {
        "chrom": chrom,
        "elements_total": total,
        "scored": scored,
        "complete": complete,
        "requests_this_run": 10,
        "summary": {
            "with_predicted_target": named,
            "fraction_with_target": named / scored,
            "strong": strong,
            "silencer_like": 1,
            "coding_target_agrees_with_nearest": agree / (agree + other + outside),
            "fraction_inside_domain": (agree + other) / (agree + other + outside),
            "verdicts_coding": {
                "agrees with nearest TSS in domain": agree,
                "another gene in the same domain": other,
                "gene outside the domain": outside,
                "no predicted effect": scored - (agree + other + outside),
            },
        },
    }


def test_aggregate_counts_complete_chromosomes_only(tmp_path):
    (tmp_path / "enhancer_targets_all_chr21.json").write_text(
        json.dumps(_res("chr21", 100, 100, True, 60, 30, 10, 10, 20))
    )
    (tmp_path / "enhancer_targets_all_chr22.json").write_text(
        json.dumps(_res("chr22", 200, 50, False, 30, 15, 5, 5, 10))
    )
    out = agw.aggregate(tmp_path)
    assert out["complete_chromosomes"] == 1 and set(out["chromosomes"]) == {"chr21", "chr22"}
    g = out["genome"]
    assert g["scored"] == 100 and g["named"] == 60 and g["coding_named"] == 50
    assert g["fraction_with_target"] == 0.6 and g["coding_target_agrees_with_nearest"] == 0.6
    assert g["coding_target_inside_domain"] == 0.8 and out["requests_total"] == 20
    assert out["chromosomes"]["chr22"]["complete"] is False


def test_spread_is_reported_per_chromosome_not_only_pooled(tmp_path):
    (tmp_path / "enhancer_targets_all_chr21.json").write_text(
        json.dumps(_res("chr21", 100, 100, True, 60, 30, 10, 10, 20))
    )
    (tmp_path / "enhancer_targets_all_chr1.json").write_text(
        json.dumps(_res("chr1", 100, 100, True, 90, 10, 10, 30, 20))
    )
    s = agw.aggregate(tmp_path)["spread"]["fraction_with_target"]
    assert s["n"] == 2 and s["min"] == 0.6 and s["max"] == 0.9 and s["range"] == 0.3
    assert s["lowest"] == "chr21" and s["highest"] == "chr1"
    assert s["per_chromosome"] == {"chr21": 0.6, "chr1": 0.9}


def test_the_acrocentric_chromosomes_are_split_out_and_the_controls_travel(tmp_path):
    (tmp_path / "enhancer_targets_all_chr21.json").write_text(
        json.dumps(_res("chr21", 100, 100, True, 60, 30, 10, 10, 20))
    )
    (tmp_path / "enhancer_targets_all_chr1.json").write_text(
        json.dumps(_res("chr1", 100, 100, True, 90, 10, 10, 30, 20))
    )
    out = agw.aggregate(tmp_path)
    t = out["size_trend"]
    assert t["acrocentric_and_chrY_only"]["chroms"] == ["chr21"]
    assert t["acrocentric_and_chrY_only"]["pooled"]["fraction_with_target"] == 0.6
    assert t["excluding_acrocentric_and_chrY"]["pooled"]["fraction_with_target"] == 0.9
    assert t["fraction_with_target"] is None  # two chromosomes is not a trend
    c = out["controls"]
    assert c["fraction_with_target"]["control"] == 0.87
    assert c["coding_target_inside_domain"]["control"] == 0.791
    assert c["coding_target_inside_domain"]["excess"] == 0.026


def test_inside_counts_pairs_the_partition_keeps_together():
    starts = [0, 100, 200]  # three nodes: [0,100), [100,200), [200,...)
    # element and target in the same node, then split by the boundary at 100
    assert agw._inside(starts, [(10, 50), (120, 150), (250, 260)]) == 3
    assert agw._inside(starts, [(10, 150), (90, 210)]) == 0
    assert agw._inside(starts, [(10, 50), (10, 150)]) == 1


def test_the_node_control_skips_a_chromosome_with_no_local_table(tmp_path, monkeypatch):
    monkeypatch.setattr(agw, "ARCHIVE", tmp_path)
    out = agw.node_control(["chr21", "chrNotAChromosome"])
    assert out["chromosomes"] == 0 and out["inside"] is None
    assert out["skipped"]["chr21"] == "no local element table"
    assert out["skipped"]["chrNotAChromosome"] == "no length"


def test_the_sign_rule_is_frozen_so_a_failure_cannot_be_refitted_away():
    """The rule was fixed on 2026-09-15 from twelve chromosomes and every later one is a held-out test
    of it. Moving the threshold or widening the refusal band after a chromosome lands would turn a
    falsifiable claim into a description of whatever arrived, so it has to break a test to do it."""
    assert agw.SIGN_RULE["threshold_per_mb"] == 5.8
    assert agw.SIGN_RULE["too_close_per_mb"] == 0.3
    assert agw.SIGN_RULE["fixed_at"] == "2026-09-15"
    assert agw.SIGN_RULE["fixed_on"] == [
        "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
        "chr18", "chr19", "chr20", "chr21", "chr22", "chrY",
    ]  # fmt: skip
    assert "chr11" not in agw.SIGN_RULE["fixed_on"]  # the first held-out test, and it held
    assert "chrX" not in agw.SIGN_RULE["fixed_on"]  # the one the rule can lose on


def test_the_chrX_failure_cannot_be_recomputed_away():
    """chrX measured +0.0007 under the null in force when the rule was fixed and -0.0007 under the
    better null adopted afterwards. The quantity is zero and the sign belongs to the seed, so a verdict
    re-derived from the current control would turn a failure into a pass every time the instrument is
    improved. The verdict as landed is frozen, and this asserts it stays frozen and stays a failure."""
    assert agw.HELD_OUT_LOG["chrX"]["verdict"] == "WRONG"
    assert agw.HELD_OUT_LOG["chrX"]["excess"] == 0.0007  # as measured at 20 draws, when it landed
    assert "NOT being widened" in agw.HELD_OUT_LOG["chrX"]["note"]

    # the current control says the opposite; the frozen verdict must win, and must say so
    control = {"per_chromosome": {"chrX": {"excess": -0.0007}}}
    p = agw.sign_prediction(control, {"chrX": 5.26})
    assert p["per_chromosome"]["chrX"]["status"] == "held out: WRONG"
    assert "recomputes_differently_now" in p["per_chromosome"]["chrX"]
    assert p["verdict"] == "FAILED on chrX"
    assert p["held_out_wrong"] == 1


def test_the_sign_rule_refuses_the_chromosomes_too_near_its_line():
    assert agw.sign_call(7.0) == "positive"
    assert agw.sign_call(4.8) == "negative"
    assert agw.sign_call(5.87) == "too close to call"  # chr9, and the rule says so rather than guessing
    assert agw.sign_call(5.26) == "negative"  # chrX, the one autosome-or-X predicted to lose
    assert agw.sign_call(None) is None


def test_the_sign_rule_scores_only_the_chromosomes_it_was_not_fitted_on():
    control = {
        "per_chromosome": {
            "chr21": {"excess": -0.0208},  # in sample, fitted on
            "chr11": {"excess": 0.02},  # held out, predicted positive at 7.0/Mb
            "chr9": {"excess": -0.001},  # held out but inside the refusal band
        }
    }
    density = {"chr21": 4.86, "chr11": 7.0, "chr9": 5.87, "chrX": 5.26}
    p = agw.sign_prediction(control, density)
    assert p["held_out_chromosomes_scored"] == 1 and p["held_out_right"] == 1
    assert p["verdict"] == "holding"
    assert p["per_chromosome"]["chr21"]["status"] == "fitted on, not a test"
    assert p["per_chromosome"]["chr11"]["status"] == "held out: right"
    assert p["per_chromosome"]["chr9"]["status"] == "held out: refused"
    assert p["per_chromosome"]["chrX"]["status"] == "still to land"
    assert p["standing_predictions"] == {"chrX": "negative"}
    assert p["chromosomes_by_status"]["held out: right"] == ["chr11"]
    assert p["scoreboard"] == (
        "held out: 1 right, 0 wrong, 1 refused; 1 fitted on and not tests; 1 still to land"
    )


def test_the_sign_rule_records_a_failure_rather_than_hiding_it():
    # chr8 has not landed and is not in the frozen log, so its verdict is derived live
    control = {"per_chromosome": {"chr8": {"excess": -0.05}}}  # predicted positive, comes back negative
    p = agw.sign_prediction(control, {"chr8": 7.65})
    assert p["held_out_wrong"] == 1 and p["held_out_right"] == 0
    assert p["per_chromosome"]["chr8"]["status"] == "held out: WRONG"
    assert p["verdict"] == "FAILED on chr8"  # it says so in the committed result, in capitals


def test_history_keeps_one_entry_per_chromosome_set(tmp_path):
    (tmp_path / "enhancer_targets_all_chr21.json").write_text(
        json.dumps(_res("chr21", 100, 100, True, 60, 30, 10, 10, 20))
    )
    first = agw.with_history(agw.aggregate(tmp_path), None)
    assert [h["chromosomes"] for h in first["history"]] == [1]
    assert first["history"][0]["fraction_with_target"] == 0.6

    # the same set re-folded replaces its entry rather than adding one
    again = agw.with_history(agw.aggregate(tmp_path), first)
    assert len(again["history"]) == 1

    # a chromosome lands and the earlier reading is kept beside it
    (tmp_path / "enhancer_targets_all_chr1.json").write_text(
        json.dumps(_res("chr1", 100, 100, True, 90, 10, 10, 30, 20))
    )
    third = agw.with_history(agw.aggregate(tmp_path), again)
    assert [h["chromosomes"] for h in third["history"]] == [1, 2]
    assert third["history"][0]["fraction_with_target"] == 0.6
    assert third["history"][1]["fraction_with_target"] == 0.75


def test_the_length_trend_needs_four_chromosomes(tmp_path):
    for chrom, named in (("chr1", 90), ("chr2", 80), ("chr21", 60), ("chr22", 62)):
        (tmp_path / f"enhancer_targets_all_{chrom}.json").write_text(
            json.dumps(_res(chrom, 100, 100, True, named, 30, 10, 10, 20))
        )
    t = agw.aggregate(tmp_path)["size_trend"]
    assert t["chromosomes"] == 4
    assert t["fraction_with_target"] == 1.0  # the rate rises with length, perfectly ranked


def test_the_module_level_controls_are_never_mutated_by_a_fold(tmp_path):
    """`aggregate` hands out the CONTROLS dict itself; a fold that wrote into it would edit the
    constant for every later fold in the same process. Both branches of main copy before writing."""
    before = json.dumps(agw.CONTROLS, sort_keys=True)
    (tmp_path / "enhancer_targets_all_chr21.json").write_text(
        json.dumps(_res("chr21", 100, 100, True, 60, 30, 10, 10, 20))
    )
    out = agw.aggregate(tmp_path)
    out["controls"] = {
        **out["controls"],
        "coding_target_inside_domain": {
            **out["controls"]["coding_target_inside_domain"],
            "measured_on_this_set": {"whatever": 1},
        },
    }
    assert json.dumps(agw.CONTROLS, sort_keys=True) == before
    assert "measured_on_this_set" not in agw.CONTROLS["coding_target_inside_domain"]


def test_the_measurability_rule_is_frozen_and_owns_its_in_sample_misses():
    """A second rule, fixed after the first failed. It must state its own fitted accuracy honestly and
    name both chromosomes it already gets wrong, or it is a description dressed as a prediction."""
    assert agw.MEASURABILITY_RULE["threshold_per_mb"] == 6.55
    assert agw.MEASURABILITY_RULE["fixed_on_chromosomes"] == 18
    assert len(agw.MEASURABILITY_RULE["known_exceptions"]) == 2
    assert any("chr20" in x for x in agw.MEASURABILITY_RULE["known_exceptions"])
    assert any("chrY" in x for x in agw.MEASURABILITY_RULE["known_exceptions"])
    for chrom in ("chr3", "chr7", "chr4", "chr6", "chr5", "chr1"):
        assert chrom not in agw.MEASURABILITY_RULE_FIXED_ON  # the six it can lose on


def test_the_measurability_rule_scores_only_chromosomes_it_was_not_fitted_on():
    control = {
        "per_chromosome": {
            "chr20": {"boundaries_per_mb": 7.40, "distinguishable_from_random": False, "p": 0.70},
            "chr3": {"boundaries_per_mb": 7.26, "distinguishable_from_random": True, "p": 0.01},
            "chr7": {"boundaries_per_mb": 7.10, "distinguishable_from_random": False, "p": 0.40},
        }
    }
    m = agw.score_measurability(control)
    assert m["per_chromosome"]["chr20"]["held_out"] is False  # fitted on, and a known miss
    assert m["per_chromosome"]["chr20"]["held"] is False
    assert m["held_out_right"] == 1 and m["held_out_wrong"] == 1
    assert m["verdict"] == "FAILED on chr7"
