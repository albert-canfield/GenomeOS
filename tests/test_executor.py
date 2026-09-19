# SPDX-License-Identifier: AGPL-3.0-or-later
"""The executor test: the read-out, the matching, the statistics and the stopping rule."""

import itertools

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


_ids = itertools.count()


def _pair(endpoint, agree_unit, agree_control, chrom="chr21"):
    sign = 1
    return {
        "endpoint": endpoint,
        "chrom": chrom,
        "unit": f"{chrom}:{next(_ids)}",
        "control": {"unit": f"{chrom}:c{next(_ids)}"},
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
            "control": {"unit": f"chr21:c{i}", "pos": 2000 + 2 * i, "ref": "C", "alt": "T"},
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


def test_a_unit_with_two_controls_counts_once():
    a = _pair("E2_eqtl", True, False)
    b = dict(a, control={"unit": "chr21:other"}, control_result={"agrees": True})
    t = ex.tally([a, b], "E2_eqtl")
    assert (t["units_called"], t["controls_called"], t["controls_agree"]) == (1, 2, 1)


def test_run_pairs_stops_at_a_pair_count(tmp_path):
    def scorer(chrom, pos, ref, alt):
        return [("GENE1", "Liver", 0.5)]

    pairs = [
        {
            "endpoint": "E3_eqtl_linked",
            "unit": f"chr21:{i}-{i + 200}",
            "test": {"pos": 1001 + 2 * i, "ref": "A", "alt": "G"},
            "control": {"unit": f"chr21:c{i}", "pos": 2000 + 2 * i, "ref": "C", "alt": "T"},
            "borrowed": {"gene": "GENE1", "tissue": "Liver", "cell": None, "measured_sign": 1},
            "order": i,
        }
        for i in range(6)
    ]
    assert ex.run_pairs(pairs, scorer, cache=tmp_path, max_pairs=3)["pairs_done"] == 3


def test_hold_out_comparison_reads_the_gap_the_way_it_was_written():
    strong = [_pair("E2_eqtl", i < 34, i < 14) for i in range(40)]  # 0.85 against 0.35
    diluted = [_pair("E3_eqtl_linked", i < 17, i < 15) for i in range(40)]  # 0.425 against 0.375
    out = ex.compare_endpoints(strong, diluted)
    assert out["gap_lower_95"] > 0 and "diluted" in out["reading"]
    same = [_pair("E3_eqtl_linked", i < 34, i < 14) for i in range(40)]
    out = ex.compare_endpoints(strong, same)
    assert out["gap"] == 0 and out["second_endpoint_meets_success"] and "neighbourhood" in out["reading"]
    thin = [_pair("E3_eqtl_linked", i < 4, i < 2) for i in range(6)]  # the same way, too few to say
    assert "neither separated" in ex.compare_endpoints(strong, thin)["reading"]
    assert ex.compare_endpoints(strong, [])["gap"] is None


def _mpra_row(pos, alt, cell, log2fc, fdr, study="A study", chrom="chr1"):
    return {
        "chrom": chrom,
        "chromStart": pos - 1,
        "ref": "A",
        "alt": alt,
        "rsid": f"rs{pos}",
        "cellLine": cell,
        "log2FC": str(log2fc),
        "fdr": str(fdr),
        "pvalue": str(fdr),
        "mpraStudy": study,
    }


def test_wide_population_and_nulls():
    rows = [
        _mpra_row(1000, "G", "GM12878", 1.2, 0.01),  # significant
        _mpra_row(1000, "G", "HepG2", 0.4, 0.2),  # the same variant, not significant: never a control
        _mpra_row(1400, "C", "GM12878", 0.02, 0.9),  # a measured null, far enough away
        _mpra_row(1050, "T", "GM12878", 0.03, 0.9),  # a null inside the tested 200-bp window
        _mpra_row(2000, "GG", "GM12878", 2.0, 0.001),  # not a single-base swap
        _mpra_row(3000, "G", "HEK293T", 2.0, 0.001),  # no track for the cell
        _mpra_row(4000, "G", "GM12878", 2.0, 0.001, study="Saturation mutagenesis of twenty"),
        _mpra_row(4500, "C", "GM12878", 0.01, 0.8, study="Saturation mutagenesis of twenty"),
    ]
    tests, nulls = ex.wide_split(rows)
    assert {t["pos"] for t in tests} == {1000, 4000}
    assert {t["endpoint"] for t in tests} == {"E1W_mpra_wide", "E1S_mpra_saturation"}
    assert {n["pos"] for n in nulls} == {1400, 1050, 4500}  # the significant variant is never a null
    assert [t["measured_sign"] for t in tests if t["pos"] == 1000] == [1]


def test_wide_readout_never_lets_the_model_choose_the_gene():
    effects = [["FAR", "GM12878", 0.5], ["MINE", "GM12878", 0.002]]
    strict = ex.readout(effects, "ABSENT", None, "GM12878", fallback=False)
    assert strict["value"] is None and "gene not in the window" in strict["how"]
    assert ex.readout(effects, "ABSENT", None, "GM12878")["value"] == 0.5  # the old E1 rule
    assert ex.readout(effects, "MINE", None, "GM12878", fallback=False)["value"] == 0.002


def test_decide_can_split_on_the_cell_line():
    done = [dict(_pair("E1W_mpra_wide", True, False), borrowed={"measured_sign": 1, "cell": "GM12878"})]
    done += [dict(_pair("E1W_mpra_wide", True, False), borrowed={"measured_sign": 1, "cell": "Jurkat"})]
    done += [
        dict(_pair("E1W_mpra_wide", True, False, "chr9"), borrowed={"measured_sign": 1, "cell": "Jurkat"})
    ]
    assert set(ex.by_group(done, "E1W_mpra_wide", "cell")) == {"GM12878", "Jurkat"}
    assert ex.decide(done, "E1W_mpra_wide", 0.5, split="cell") == "success"


def test_pre_registration_of_the_widened_e1_is_complete():
    for key in (
        "claim_tested",
        "population",
        "strata",
        "readout",
        "controls",
        "order",
        "budget",
        "success",
        "confirms_e2_e3",
        "contradicts_e2_e3",
        "fails_to_resolve",
        "convention",
        "caveat",
    ):
        assert ex.E1_WIDE[key]
    assert "before any model request" in ex.E1_WIDE["written"]


def test_criterion_is_complete_and_pre_registered():
    for key in ("claim_tested", "prediction", "success", "stopping", "order", "alpha_spending", "caveat"):
        assert ex.CRITERION[key]
    assert "before any model request" in ex.CRITERION["written"]
    assert ex.plan([{"endpoint": "E1_mpra"}, {"endpoint": "E2_eqtl"}])["requests_if_all_run"] == 4


def test_round_robin_puts_breadth_before_depth():
    def pair(chrom, unit, endpoint="E2R_eqtl_replication"):
        return {"endpoint": endpoint, "unit": f"{chrom}:{unit}", "borrowed": {"measured_sign": 1}}

    pairs = [pair("chr6", i) for i in range(3) for _ in range(2)]
    pairs += [pair("chr19", i) for i in range(2) for _ in range(2)]
    pairs += [pair("chr7", 0, "E3R_eqtl_linked_replication")]
    out = ex.round_robin(pairs)
    assert [p["unit"] for p in out[:4]] == ["chr6:0", "chr6:0", "chr19:0", "chr19:0"]
    assert out[-1]["endpoint"] == "E3R_eqtl_linked_replication"  # E2R is spent first
    assert [p["order"] for p in out] == list(range(1, len(out) + 1))


def test_leave_one_out_finds_the_chromosome_that_carries_a_difference():
    rows = [_pair("E2R_eqtl_replication", True, False, "chr6") for _ in range(20)]
    rows += [_pair("E2R_eqtl_replication", False, False, "chr7") for _ in range(20)]
    out = ex.leave_one_out(rows, "E2R_eqtl_replication")
    assert out["lowest"]["chromosome_removed"] == "chr6" and out["lowest"]["difference"] == 0.0
    assert out["holds_above_0.05"] is False


def test_the_replication_pre_registration_names_every_outcome():
    for key in (
        "population",
        "unchanged",
        "order",
        "replicates",
        "fails_to_replicate",
        "resolves_nothing",
        "dilution_again",
        "what_a_success_would_still_not_mean",
        "caveat",
    ):
        assert ex.E2_REPLICATION[key]
    assert set(ex.REPLICATION_OF.values()) == set(ex.REPLICATION_ENDPOINTS)
    assert "chr21" in ex.DISCOVERY_CHROMS and "chr22" in ex.DISCOVERY_CHROMS
