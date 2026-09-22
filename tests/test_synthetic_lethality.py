"""Synthetic lethality: the statistics, the loss calls and the pre-registration.

Every number here is made up on purpose. The point of these tests is that the
rank test, the loss definitions and the verdict rules behave as declared before
DepMap is ever contacted, so a real run cannot quietly redefine what "recovered"
means. No test touches the network.
"""

from __future__ import annotations

import gzip
import io
import json

import pytest

from genomeos.therapeutics import depmap
from genomeos.therapeutics import synthetic_lethality as sl


def test_ranks_average_ties_and_report_the_tie_correction():
    rank_of, ties = sl.ranks({"a": 1.0, "b": 2.0, "c": 2.0, "d": 5.0})
    assert rank_of == {"a": 1.0, "b": 2.5, "c": 2.5, "d": 4.0}
    assert ties == 6.0  # one run of two: 2^3 - 2
    assert sl.ranks({"only": 1.0}) == ({"only": 1.0}, 0.0)


def test_rank_sum_is_one_sided_and_directional():
    low = {f"L{i}": float(i) for i in range(6)}
    high = {f"H{i}": 10.0 + i for i in range(20)}
    rank_of, ties = sl.ranks(low | high)
    p_low, delta_low, _u = sl.rank_sum_p(rank_of, ties, list(low))
    p_high, delta_high, _u2 = sl.rank_sum_p(rank_of, ties, list(high))
    assert p_low < 0.001 and delta_low == -1.0
    assert p_high > 0.99 and delta_high == 1.0
    assert sl.rank_sum_p(rank_of, ties, []) == (1.0, 0.0, 0.0)


def test_benjamini_hochberg_is_monotone_and_never_below_the_p_value():
    ps = [0.001, 0.02, 0.5, 0.9]
    qs = sl.benjamini_hochberg(ps)
    assert qs == sorted(qs)
    assert all(q >= p for p, q in zip(ps, qs, strict=True))
    assert sl.benjamini_hochberg([]) == []
    assert sl.benjamini_hochberg([0.5]) == [0.5]


def test_loss_definitions_are_kept_apart():
    damaging = {"m1": {"A": 1.0}, "m2": {"A": 0.0}, "m3": {"A": 0.0}}
    cn = {"m1": {"A": 2.1}, "m2": {"A": 0.0}, "m4": {"A": 1.9}}
    msi = {"m1": 35.0, "m2": 3.0}
    assert sl.loss_status("A", "damaging_mutation", damaging, cn, msi) == {
        "m1": True,
        "m2": False,
        "m3": False,
    }
    assert sl.loss_status("A", "deletion", damaging, cn, msi) == {"m1": False, "m2": True, "m4": False}
    # the union: m3 has no copy-number call, so it cannot be called intact
    assert sl.loss_status("A", "loss", damaging, cn, msi) == {"m1": True, "m2": True}
    assert sl.loss_status("A", "msi_high", damaging, cn, msi) == {"m1": True, "m2": False}


def _planted(dependency_effect: float = -1.2, lost: int = 10, intact: int = 90):
    effect = {}
    for i in range(lost):
        effect[f"L{i}"] = {"B": dependency_effect - 0.01 * i, "C": 0.01 * ((i % 5) - 2)}
    for i in range(intact):
        effect[f"K{i}"] = {"B": -0.1 + 0.001 * i, "C": 0.01 * ((i % 7) - 3)}
    status = {m: m.startswith("L") for m in effect}
    return effect, status


def test_a_planted_pair_is_recovered_and_an_unrelated_gene_is_not():
    effect, status = _planted()
    store = sl._RankStore(effect)
    hit = sl.test_pair("A", "B", "loss", status, store)
    miss = sl.test_pair("A", "C", "loss", status, store)
    assert hit["cliffs_delta"] == -1.0
    assert hit["p"] < 1e-5
    assert hit["median_difference"] < -1.0
    assert miss["p"] > 0.05
    qs = sl.benjamini_hochberg([hit["p"], miss["p"]])
    assert qs[0] <= sl.FDR < qs[1] or qs[1] > sl.FDR


def test_an_underpowered_pair_is_untestable_rather_than_negative():
    effect, _status = _planted(lost=2)
    status = {m: m.startswith("L") for m in effect}
    store = sl._RankStore(effect)
    row = sl.test_pair("A", "B", "loss", status, store)
    assert row["verdict"] == "untestable"
    assert "2 lost" in row["reason"]
    missing = sl.test_pair("A", "NOSUCHGENE", "loss", status, store)
    assert missing["verdict"] == "untestable"


def test_the_shuffled_null_destroys_the_planted_pair():
    effect, status = _planted()
    store = sl._RankStore(effect)
    null = sl.shuffled_null([("A", "B"), ("A", "C")], {"A": status}, store, permutations=10)
    assert null["permutations"] == 10
    assert null["mean_hits"] <= 1.0
    assert len(null["hits_per_permutation"]) == 10


def test_the_candidate_space_is_the_panel_plus_paralogue_pairs_both_ways():
    pairs, space = sl.candidate_pairs({"SMARCA4": [{"symbol": "SMARCA2", "identity": 74.0}]})
    assert space["panel_genes"] == len(sl.PANEL)
    assert ("SMARCA4", "SMARCA2") in pairs and ("SMARCA2", "SMARCA4") in pairs
    assert ("BRCA1", "PARP1") in pairs
    assert all(a != b for a, b in pairs)
    # a distant paralogue is below the declared identity floor and is dropped
    _pairs, small = sl.candidate_pairs(
        {"MTAP": [{"symbol": "ADI1", "identity": 5.0}]}, panel=("MTAP", "PRMT5")
    )
    assert small["paralogue_pairs"] == 0


def test_the_preregistration_is_explicit_about_the_control_expected_to_fail():
    parp = [c for c in sl.PREREGISTERED if c["dependency"] == "PARP1"]
    assert len(parp) == 2
    assert all(c["expect"] == "weak" for c in parp)
    assert {c["dependency"] for c in sl.PREREGISTERED if c["expect"] == "recovered"} == {
        "ARID1B",
        "SMARCA2",
        "PRMT5",
        "MAT2A",
        "WRN",
    }
    assert sl.FDR == 0.10 and sl.MIN_LOST == 5 and sl.MSI_HIGH == 20.0


def test_control_verdicts_follow_the_declared_rules():
    effect = {}
    for i in range(12):
        effect[f"L{i}"] = {"ARID1B": -1.1 - 0.01 * i, "PARP1": -0.2, "SMARCA2": -0.1}
    for i in range(60):
        effect[f"K{i}"] = {"ARID1B": -0.1 + 0.002 * i, "PARP1": -0.2, "SMARCA2": -0.1}
    damaging = {m: {"ARID1A": 1.0 if m.startswith("L") else 0.0} for m in effect}
    rows = sl.control_rows(damaging, {}, {}, sl._RankStore(effect), [])
    arid = next(r for r in rows if r["dependency"] == "ARID1B")
    assert arid["verdict"] == "recovered" and arid["agrees_with_preregistration"]
    parp = next(r for r in rows if r["dependency"] == "PARP1")
    assert parp["verdict"] in ("not_recovered", "untestable")
    assert parp["agrees_with_preregistration"]  # "weak" is satisfied by not recovering it
    smarca = next(r for r in rows if r["dependency"] == "SMARCA2")
    assert smarca["verdict"] == "untestable"  # no SMARCA4 calls in this synthetic table
    assert not smarca["agrees_with_preregistration"]


def test_evidence_from_a_hit_is_a_prediction_and_never_clinical():
    ev = sl.evidence_for(
        {
            "lost": "MTAP",
            "loss": "deletion",
            "dependency": "PRMT5",
            "lines_lost": 40,
            "lines_intact": 900,
            "cliffs_delta": -0.6,
            "q": 0.0,
        }
    )
    assert ev.source_type == "prediction" and ev.level == "computational"
    assert ev.confidence <= 0.7
    assert "PRMT5" in ev.claim and "MTAP" in ev.claim


def test_depmap_column_symbols_and_the_release_are_pinned():
    assert depmap.symbol("PRMT5 (10419)") == "PRMT5"
    assert depmap.symbol('"MTAP (4507)"') == "MTAP"
    assert depmap.RELEASE == "DepMap 24Q4 Public"
    assert depmap.provenance()["files"]["gene_effect"]["name"] == "CRISPRGeneEffect.csv"
    assert any("cell lines, not patient tumours" in limit for limit in depmap.LIMITS)


def test_an_unnamed_index_column_falls_back_and_an_empty_table_is_loud(monkeypatch, tmp_path):
    # OmicsSignatures.csv leaves its index column unnamed; asking for ModelID used to
    # produce an empty table in silence, and the MSI control came back untestable.
    csv_text = ",MSIScore,Ploidy\nACH-000001,3.68,2.1\nACH-000002,66.5,3.0\n"

    def _response(text):
        return io.BytesIO(text.encode())

    monkeypatch.setattr(depmap, "_stream", lambda file_id, timeout=180: _response(csv_text))
    got = depmap.stream_narrow("signatures", "ModelID", ("MSIScore",))
    assert got["key"] == "" and got["key_requested"] == "ModelID"
    assert got["values"]["ACH-000002"]["MSIScore"] == "66.5"
    monkeypatch.setattr(depmap, "_stream", lambda file_id, timeout=180: _response(",MSIScore\n"))
    with pytest.raises(depmap.DepMapUnavailableError, match="returned no rows"):
        depmap.stream_narrow("signatures", "ModelID", ("MSIScore",))


def test_the_depmap_cache_is_used_offline_and_an_uncovered_request_is_refused(tmp_path):
    payload = {
        "dataset": "gene_effect",
        "release": depmap.RELEASE,
        "requested": ["PRMT5", "WRN"],
        "genes": ["PRMT5", "WRN"],
        "missing": [],
        "columns": [],
        "values": {"ACH-000001": {"PRMT5": -1.1, "WRN": -0.2}},
    }
    path = tmp_path / "gene_effect.json.gz"
    with gzip.open(path, "wt") as fh:
        json.dump(payload, fh)
    got = depmap.matrix("gene_effect", {"PRMT5"}, net=False, cache_dir=tmp_path)
    assert got["from_cache"] and got["values"]["ACH-000001"]["PRMT5"] == -1.1
    with pytest.raises(depmap.DepMapUnavailableError, match="network disabled"):
        depmap.matrix("gene_effect", {"PRMT5", "NEWGENE"}, net=False, cache_dir=tmp_path)


def test_co_deleted_loss_genes_are_reported_as_one_event():
    lines = [f"ACH-{i}" for i in range(20)]
    mtap = {m: i < 8 for i, m in enumerate(lines)}
    cdkn2a = {m: i < 12 for i, m in enumerate(lines)}  # the same 8 lines plus 4 more
    unrelated = {m: i >= 15 for i, m in enumerate(lines)}
    rows = sl.co_loss({"MTAP": mtap, "CDKN2A": cdkn2a, "ARID1A": unrelated}, ["ARID1A", "CDKN2A", "MTAP"])
    assert [r["genes"] for r in rows] == [["CDKN2A", "MTAP"]]
    assert rows[0]["shared_lines"] == 8
    assert rows[0]["overlap_of_the_smaller_set"] == 1.0
    assert "cannot say" in rows[0]["reading"]
