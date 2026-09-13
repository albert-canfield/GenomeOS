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
