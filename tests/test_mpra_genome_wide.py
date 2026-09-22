# SPDX-License-Identifier: AGPL-3.0-or-later
"""The genome-wide fold of the lentiMPRA results, on two synthetic chromosomes."""

import importlib.util
import json
from pathlib import Path

from genomeos.attribution import mpra

spec = importlib.util.spec_from_file_location("mgw", Path("scripts/mpra_genome_wide.py"))
mgw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mgw)


def _row(key, act_k, act_h, pk, ph):
    return {
        "key": key,
        "ccre_class": "dELS",
        "activity": {"K562": act_k, "HepG2": act_h},
        "active": {"K562": act_k >= mpra.ACTIVE, "HepG2": act_h >= mpra.ACTIVE},
        "open": {"K562": pk > 0.5, "HepG2": ph > 0.5},
        "constrained_fraction": 0.1,
        "predicted_dnase": {"K562": pk, "HepG2": ph},
    }


def _result(chrom, rows):
    return {
        "chrom": chrom,
        "annotated": len(rows),
        "rows": rows,
        "summary": mpra.summarise(rows),
        "model_cost": {"requests": 5, "seconds": 9.0},
    }


def test_aggregate_folds_rows(tmp_path):
    a = [_row("a", 2.0, -1.0, 0.9, 0.1), _row("b", -0.5, 1.5, 0.2, 0.8), _row("c", -1.0, -1.0, 0.1, 0.1)]
    b = [_row("d", 1.2, 1.3, 0.7, 0.7), _row("e", -0.2, -0.3, 0.3, 0.2), _row("f", 1.1, -0.9, 0.6, 0.2)]
    (tmp_path / "mpra_chr21.json").write_text(json.dumps(_result("chr21", a)))
    (tmp_path / "mpra_chr22.json").write_text(json.dumps(_result("chr22", b)))
    (tmp_path / "mpra_genome_wide.json").write_text("{}")
    out = mgw.aggregate(tmp_path)
    assert out["elements"] == 6 and set(out["chromosomes"]) == {"chr21", "chr22"}
    assert out["model_requests"] == 10
    g = out["genome"]
    assert g["cells"]["K562"]["elements"] == 6 and g["cells"]["K562"]["active"] == 0.5
    assert g["cells"]["K562"]["predicted_dnase_vs_activity"]["same_cell"]["spearman"] > 0.8
    assert g["specific_elements"]["elements"] == 3
    assert g["specific_elements"]["predicted_dnase_higher_in_the_active_cell"] == 1.0
    assert out["chromosomes"]["chr21"]["active"] == {"K562": 0.333, "HepG2": 0.333}


def test_order_is_smallest_first():
    assert mgw.ORDER[0] == "chr21" and "chrM" not in mgw.ORDER and len(mgw.ORDER) == 24
