# SPDX-License-Identifier: AGPL-3.0-or-later
"""The genome-wide fold of the VISTA results, on two synthetic chromosomes."""

import importlib.util
import json
from pathlib import Path

spec = importlib.util.spec_from_file_location("vgw", Path("scripts/vista_genome_wide.py"))
vgw = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vgw)


def _result(chrom, rows):
    from genomeos.attribution import vista

    full = []
    for r in rows:
        r = dict(r)
        r["enhancer_like"] = any(c in ("dELS", "pELS") for c in r["ccre_classes"])
        r["any_ccre"] = bool(r["ccre_classes"])
        r["groups"] = sorted(vista.tissue_groups(r["tissues"]))
        full.append(r)
    return {
        "chrom": chrom,
        "rows": rows,
        "summary": vista.summarise(full),
        "phylop_cost": {"requests": 3, "mb_fetched": 0.2},
    }


def _row(vid, status, tissues, classes, cf, predicted, group=None, agrees=None):
    return {
        "id": vid,
        "status": status,
        "tissues": tissues,
        "ccre_classes": classes,
        "domain": "d",
        "constrained_fraction": cf,
        "predicted": predicted,
        "verdict_coding": "agrees with nearest TSS in domain" if predicted else "no predicted effect",
        "predicted_group": group,
        "tissue_agrees": agrees,
    }


def test_aggregate_folds_rows_and_keeps_per_chromosome_headline(tmp_path):
    strong = {"strength": "strong", "tissue": "Brain_Cerebellum"}
    weak = {"strength": "weak", "tissue": "K562"}
    a = _result(
        "chr21",
        [
            _row("hs1", "positive", ["fb"], ["dELS"], 0.5, strong, "neural", True),
            _row("hs2", "negative", [], ["dELS"], 0.1, None),
        ],
    )
    b = _result(
        "chr22",
        [
            _row("hs3", "positive", ["ht"], [], 0.3, weak),
            _row("hs4", "negative", [], ["CTCF-only"], 0.25, weak),
        ],
    )
    (tmp_path / "vista_chr21.json").write_text(json.dumps(a))
    (tmp_path / "vista_chr22.json").write_text(json.dumps(b))
    (tmp_path / "vista_genome_wide.json").write_text("{}")
    out = vgw.aggregate(tmp_path)
    assert out["elements"] == 4 and set(out["chromosomes"]) == {"chr21", "chr22"}
    g = out["genome"]
    assert g["positive"]["elements"] == 2 and g["negative"]["elements"] == 2
    assert g["positive"]["enhancer_like"] == 0.5 and g["negative"]["enhancer_like"] == 0.5
    assert g["positive"]["constrained"] == 1.0 and g["negative"]["constrained"] == 0.5
    assert g["positive"]["fraction_with_target"] == 1.0 and g["negative"]["fraction_with_target"] == 0.5
    assert g["positive"]["tissue_judged"] == 1 and g["positive"]["tissue_agrees"] == 1.0
    assert out["chromosomes"]["chr21"]["fraction_with_target"] == [1.0, 0.0]
    assert out["phylop_cost"] == {"requests": 6, "mb_fetched": 0.4}


def test_order_has_no_chrY_or_chrM():
    assert vgw.ORDER[0] == "chr21" and "chrY" not in vgw.ORDER and "chrM" not in vgw.ORDER
    assert len(vgw.ORDER) == 23
