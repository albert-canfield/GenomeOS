# SPDX-License-Identifier: AGPL-3.0-or-later
"""The genome-wide fold of the constrained-enhancer results, on two synthetic chromosomes."""

import importlib.util
import json
from pathlib import Path

spec = importlib.util.spec_from_file_location("ctgw", Path("scripts/constrained_targets_genome_wide.py"))
ctgw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ctgw)


def _chrom(chrom, scored, named, agree, other, outside, strong, uni_named):
    return {
        "chrom": chrom,
        "distal_enhancers": 1000,
        "constrained_elements_available": 50,
        "summary": {
            "elements_scored": scored,
            "with_predicted_target": named,
            "fraction_with_target": named / scored,
            "strong": strong,
            "silencer_like": 3,
            "coding_target_agrees_with_nearest": agree / (agree + other + outside),
            "fraction_inside_domain": (agree + other) / (agree + other + outside),
            "verdicts_coding": {
                "agrees with nearest TSS in domain": agree,
                "another gene in the same domain": other,
                "gene outside the domain": outside,
                "no predicted effect": scored - named,
            },
        },
        "uniform_sample": {
            "elements_scored": 200,
            "fraction_with_target": uni_named / 200,
            "strong": 40,
            "silencer_like": 30,
            "coding_target_agrees_with_nearest": 0.65,
            "fraction_inside_domain": 0.9,
        },
    }


def test_aggregate_totals_and_ratios(tmp_path):
    (tmp_path / "constrained_targets_chr21.json").write_text(
        json.dumps(_chrom("chr21", 100, 80, 40, 20, 10, 30, 120))
    )
    (tmp_path / "constrained_targets_chr22.json").write_text(
        json.dumps(_chrom("chr22", 50, 40, 20, 10, 5, 15, 130))
    )
    (tmp_path / "constrained_targets_genome_wide.json").write_text("{}")  # must be ignored by the glob
    out = ctgw.aggregate(tmp_path)
    g = out["genome"]
    assert g["chromosomes"] == 2
    assert g["scored"] == 150 and g["with_predicted_target"] == 120
    assert g["fraction_with_target"] == 0.8
    assert g["strong"] == 45
    assert g["coding_targets_named"] == 105
    assert g["coding_target_agrees_with_nearest"] == round(60 / 105, 3)
    assert g["coding_target_inside_domain"] == round(90 / 105, 3)
    u = out["uniform_samples"]
    assert u["scored"] == 400 and u["fraction_with_target"] == 0.625
    assert out["chromosomes"]["chr22"]["uniform_fraction_with_target"] == 0.65
    assert set(out["chromosomes"]) == {"chr21", "chr22"}


def test_order_is_smallest_first_and_skips_chrM():
    assert ctgw.ORDER[0] == "chr21"
    assert "chrM" not in ctgw.ORDER
    assert len(ctgw.ORDER) == 24 and ctgw.ORDER[-1] == "chr1"
