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
