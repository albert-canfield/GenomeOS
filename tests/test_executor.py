# SPDX-License-Identifier: AGPL-3.0-or-later
"""The executor test: the read-out, the matching, the statistics and the stopping rule."""

import pytest

from genomeos.attribution import executor as ex


def test_readout_prefers_the_named_tissue_then_the_mean():
    effects = [["ABO", "Whole_Blood", 0.4], ["ABO", "Liver", -0.2], ["OTHER", "Whole_Blood", 1.0]]
    assert ex.readout(effects, "ABO", "Whole Blood", None)["value"] == 0.4
    mean = ex.readout(effects, "ABO", "Spleen", None)
    assert mean["value"] == pytest.approx(0.1) and "no tissue match" in mean["how"]
    assert ex.readout(effects, "MISSING", None, None)["value"] is None


def test_readout_in_a_cell_picks_the_strongest_gene_when_none_is_named():
    effects = [["A", "GM12878 RNA", 0.2], ["B", "GM12878 RNA", -0.9], ["B", "K562 RNA", 2.0]]
    r = ex.readout(effects, None, None, "GM12878")
    assert r["gene"] == "B" and r["value"] == pytest.approx(-0.9)
    assert ex.readout(effects, None, None, "HepG2")["value"] is None


def test_panel_r2_reads_linkage_and_ignores_the_uninformative():
    a = ["=", "=", "G", "G", None, "*"]
    b = ["=", "=", "T", "T", "=", "="]
    assert ex.r2(a * 3, b * 3) == pytest.approx(1.0)
    c = ["=", "T", "=", "T", "=", "T"]
    assert ex.r2(a * 3, c * 3) < 0.3
    assert ex.r2(["="] * 12, ["="] * 12) is None


def test_value_bucket_and_order_put_mpra_first():
    assert (ex.value_bucket(2), ex.value_bucket(3), ex.value_bucket(7)) == ("2", "3", "4+")
    unit = {"rt_stratum": 2, "recurring_values": 2}
    assert ex.order_key(unit, {"endpoint": "E1_mpra", "rank_key": -1.0}) < ex.order_key(
        unit, {"endpoint": "E2_eqtl", "rank_key": 1e-20}
    )
    late = {"rt_stratum": 0, "recurring_values": 5}
    assert ex.order_key(unit, {"endpoint": "E2_eqtl", "rank_key": 0.5}) < ex.order_key(
        late, {"endpoint": "E2_eqtl", "rank_key": 1e-9}
    )


def test_tested_variant_prefers_a_significant_mpra_row():
    unit = {
        "measured": [
            {
                "pos": 100,
                "ref": "A",
                "alt": "G",
                "rsid": "rs1",
                "eqtl": [{"gene": "GENE1", "tissue": "Liver", "slope": -0.3, "pval": 1e-9}],
                "mpra": [{"cell": "GM12878", "log2fc": 0.8, "fdr": 0.01, "pvalue": 0.001, "study": "s"}],
                "dapg": [],
            },
            {
                "pos": 150,
                "ref": "C",
                "alt": "T",
                "rsid": "rs2",
                "eqtl": [{"gene": "GENE2", "tissue": "Lung", "slope": 0.5, "pval": 1e-12}],
                "dapg": [{"gene": "GENE2", "tissue": "Lung", "pip": 0.9}],
                "mpra": [],
            },
        ]
    }
    t = ex.tested_variant(unit)
    assert t["endpoint"] == "E1_mpra" and t["pos"] == 100 and t["measured_sign"] == 1 and t["gene"] == "GENE1"
    only_eqtl = {"measured": [unit["measured"][1]]}
    t2 = ex.tested_variant(only_eqtl)
    assert t2["endpoint"] == "E2_eqtl" and t2["measured_sign"] == 1 and t2["tissue"] == "Lung"
    linked = {"measured": [dict(unit["measured"][1], dapg=[])]}
    assert ex.tested_variant(linked)["endpoint"] == "E3_eqtl_linked"
    utr = {
        "measured": [
            dict(unit["measured"][0], mpra=[{**unit["measured"][0]["mpra"][0], "study": "3'UTR screen"}])
        ]
    }
    assert ex.tested_variant(utr)["endpoint"] == "E3_eqtl_linked"
    nofdr = [{"cell": "GM12878", "log2fc": -0.4, "fdr": None, "pvalue": 0.004, "study": "s"}]
    assert ex.mpra_significant(nofdr) and not ex.mpra_significant([{**nofdr[0], "pvalue": 0.2}])
    assert (
        ex.tested_variant(
            {"measured": [{"pos": 1, "eqtl": [], "mpra": [], "ref": "A", "alt": "T", "rsid": None}]}
        )
        is None
    )


def test_statistics():
    assert ex.fisher_greater(10, 10, 0, 10) < 1e-4
    assert ex.fisher_greater(5, 10, 5, 10) > 0.5
    assert ex.difference_upper(8, 10, 5, 10) > 0.3
    assert ex.binom_two_sided(10, 10) == pytest.approx(0.001953, abs=1e-5)
    assert ex.binom_two_sided(5, 10) == 1.0


def _pair(endpoint, agree_unit, agree_control, chrom="chr21"):
    sign = 1
    return {
        "endpoint": endpoint,
        "chrom": chrom,
        "unit": f"{chrom}:1-2",
        "borrowed": {"measured_sign": sign},
        "test_result": {"agrees": agree_unit},
        "control_result": {"agrees": agree_control},
    }


def test_tally_and_decide():
    done = [_pair("E1_mpra", True, False) for _ in range(20)] + [
        _pair("E1_mpra", True, False, "chr22") for _ in range(20)
    ]
    t = ex.tally(done, "E1_mpra")
    assert t["difference"] == 1.0 and t["p_one_sided"] < 1e-9
    assert ex.decide(done, "E1_mpra", 0.0005) == "success"
    flat = [_pair("E2_eqtl", i % 2 == 0, i % 2 == 0) for i in range(1000)]
    assert ex.decide(flat, "E2_eqtl", 0.002) == "futility"
    assert ex.decide([], "E2_eqtl", 0.01, final=True) == "negative"


def test_run_pairs_spends_one_request_a_side_and_stops(tmp_path, monkeypatch):
    calls = []

    def scorer(chrom, pos, ref, alt):
        calls.append(pos)
        return [("GENE1", "Liver", 0.5 if pos % 2 else -0.5)]

    pairs = [
        {
            "endpoint": "E2_eqtl",
            "unit": f"chr21:{i}-{i + 200}",
            "test": {"pos": 1001 + 2 * i, "ref": "A", "alt": "G"},
            "control": {"pos": 2000 + 2 * i, "ref": "C", "alt": "T"},
            "borrowed": {"gene": "GENE1", "tissue": "Liver", "cell": None, "measured_sign": 1},
            "order": i,
        }
        for i in range(4)
    ]
    out = ex.run_pairs(pairs, scorer, cache=tmp_path, max_requests=4)
    assert out["requests_spent"] == 4 and out["pairs_done"] == 2 and len(calls) == 4
    assert out["E2_eqtl"]["units_agree"] == 2 and out["E2_eqtl"]["controls_agree"] == 0
    assert out["E2_eqtl"]["verdict"] in ("negative", "success", "futility")
    again = ex.run_pairs(pairs[:1], scorer, cache=tmp_path)
    assert again["requests_spent"] == 0  # the cache answers both sides
    _ = monkeypatch


def test_criterion_is_complete_and_pre_registered():
    for key in ("claim_tested", "prediction", "success", "stopping", "order", "alpha_spending", "caveat"):
        assert ex.CRITERION[key]
    assert "before any model request" in ex.CRITERION["written"]
    assert ex.plan([{"endpoint": "E1_mpra"}, {"endpoint": "E2_eqtl"}])["requests_if_all_run"] == 4
