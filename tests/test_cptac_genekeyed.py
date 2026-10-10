# SPDX-License-Identifier: AGPL-3.0-or-later
"""Gene-keyed CPTAC sites: the registered gene-key rule, the isoform position rule, the result's gates."""

from __future__ import annotations

import json

from genomeos.molecules import cptac_keys as ck
from genomeos.molecules import ptm


def _hgnc() -> ck.Hgnc:
    rows = [
        {
            "hgnc_id": "HGNC:1",
            "symbol": "CNOT12",
            "status": "Approved",
            "prev_symbol": "TNKS1BP1",
            "alias_symbol": "TAB182",
            "uniprot_ids": "Q9C0C2",
            "ensembl_gene_id": "ENSG01",
        },
        {
            "hgnc_id": "HGNC:2",
            "symbol": "AFDN",
            "status": "Approved",
            "prev_symbol": "MLLT4",
            "alias_symbol": "AF6",
            "uniprot_ids": "P55196",
            "ensembl_gene_id": "ENSG02",
        },
        {
            "hgnc_id": "HGNC:3",
            "symbol": "GENEA",
            "status": "Approved",
            "prev_symbol": "OLD1",
            "alias_symbol": "SHARED",
            "uniprot_ids": "",
            "ensembl_gene_id": "ENSG03",
        },
        {
            "hgnc_id": "HGNC:4",
            "symbol": "GENEB",
            "status": "Approved",
            "prev_symbol": '"OLD1|GENEA"',
            "alias_symbol": "SHARED",
            "uniprot_ids": '"P1|P2"',
            "ensembl_gene_id": "ENSG04",
        },
        {
            "hgnc_id": "HGNC:5",
            "symbol": "GONE",
            "status": "Entry Withdrawn",
            "prev_symbol": "",
            "alias_symbol": "",
            "uniprot_ids": "P9",
            "ensembl_gene_id": "",
        },
    ]
    return ck.Hgnc.from_rows(rows)


def test_gene_key_tiers_in_order_and_ambiguity_refused():
    h = _hgnc()
    assert h.record("AFDN") == ("resolved", "HGNC:2", "approved symbol")
    assert h.record("MLLT4") == ("resolved", "HGNC:2", "previous symbol")
    assert h.record("AF6") == ("resolved", "HGNC:2", "alias symbol")
    # previous symbol of two records: refused, never the first one found
    assert h.record("OLD1")[0] == "ambiguous"
    # approved for one record and previous symbol of another: refused
    assert h.record("GENEA")[0] == "ambiguous"
    assert h.record("SHARED")[0] == "ambiguous"
    assert h.record("GONE") == ("unresolved", None, "symbol unknown to HGNC")
    assert h.record("NOPE")[0] == "unresolved"


def test_gene_key_needs_exactly_one_compiled_entry():
    h = _hgnc()
    compiled = {"Q9C0C2", "P55196", "P1", "P2", "X3"}
    by_ens = {"ENSG03": ["X3"]}
    assert ck.resolve_gene_key("TNKS1BP1", h, compiled, by_ens)[:2] == ("resolved", "Q9C0C2")
    # no UniProt accession on the record: the Ensembl gene among the compiled entries decides
    st, acc, how = ck.resolve_gene_key("GENEB", h, compiled, by_ens)
    assert st == "ambiguous" and acc is None and "2 compiled UniProt entries" in how
    h.approved["SOLO"] = "HGNC:3"
    h.previous.pop("GENEA", None)
    assert ck.resolve_gene_key("SOLO", h, compiled, by_ens)[:2] == ("resolved", "X3")
    assert ck.resolve_gene_key("SOLO", h, compiled, {})[0] == "unresolved"


CANON = "MAAAKLSPEDRRTYWGGHNVCQPLKSTREDQ"


def test_site_on_the_canonical_resolves_unchanged():
    pos = CANON.index("S") + 1
    assert ck.resolve_site("S", pos, CANON, []) == ("resolved", pos, "canonical")


def test_site_on_an_isoform_with_an_n_terminal_extension_is_carried_by_its_window():
    iso = "MDEFGHIKL" + CANON
    site = CANON.index("KST") + 2  # S26 on the canonical
    st, pos, how = ck.resolve_site("S", site + 9, CANON, [iso])
    # the canonical is too short to carry S at the isoform's position, so only the isoform carries it
    assert (st, pos, how) == ("resolved", site, "isoform window")


def test_carriers_that_disagree_are_refused():
    # the isoform lacks six residues; T20 on it is T26 on the canonical, which carries T20 by itself
    pre, tail = "MKWLPRHGDE", "NQFTCVMIKTLGERHWSTPDNQ"
    canon, iso = pre + "AAAAAA" + tail, pre + tail
    assert canon[19] == iso[19] == canon[25] == "T"
    assert ck.resolve_site("T", 20, canon, []) == ("resolved", 20, "canonical")
    st, pos, how = ck.resolve_site("T", 20, canon, [iso])
    assert (st, pos) == ("ambiguous", None) and "disagree" in how


def test_no_carrier_is_unresolved_and_isoform_only_sequence_is_off_canonical():
    assert ck.resolve_site("Y", 2, CANON, [])[0] == "unresolved"
    assert ck.resolve_site("S", 999, CANON, [])[0] == "unresolved"
    iso = CANON + "WWWWWWWWSWWWWWWWW"
    st, pos, _ = ck.resolve_site("S", len(CANON) + 9, CANON, [iso])
    assert (st, pos) == ("off_canonical", None)


def test_estimated_wrong_scales_the_unresolved_by_the_decoy_odds():
    assert ck.estimated_wrong(100, 10, 90) == 100 * 10 / 90
    assert ck.estimated_wrong(5, 3, 0) is None


def test_the_registered_bounds_are_the_2026_09_28_bounds_not_loosened():
    g = ck.REGISTRATION["gates"]
    assert g["G2_unresolved_rate_below"] == ptm.CPTAC_REGISTRATION["gene_keyed_mismatch_bound"] == 0.02
    assert tuple(g["G1_median_of_site_medians_within"]) == tuple(
        ptm.CPTAC_REGISTRATION["median_log2_ratio_within"]
    )
    assert g["G0_rule_wrong_rate_on_refseq_truth_below"] == 0.02
    assert g["G3_estimated_wrong_rate_below"] == 0.02
    assert set(ck.GENE_KEYED_PROFILES).isdisjoint(ck.REFSEQ_KEYED_PROFILES)
    assert len(ck.GENE_KEYED_PROFILES) == 7


def test_the_registration_result_matches_the_code_constants():
    if not ck.REGISTRATION_RESULT.exists():
        return
    r = json.loads(ck.REGISTRATION_RESULT.read_text())
    assert r["outcome_seen"] is False
    reg = r["registration"]
    for k, v in ck.REGISTRATION["gates"].items():
        assert reg["gates"][k] == (list(v) if isinstance(v, tuple) else v)
    assert reg["position_rule"] == ck.REGISTRATION["position_rule"]
    assert reg["gene_key_rule"] == ck.REGISTRATION["gene_key_rule"]


def test_the_result_commits_only_profiles_that_passed_every_gate():
    if not ck.RESULT.exists():
        return
    r = json.loads(ck.RESULT.read_text())
    g = ck.REGISTRATION["gates"]
    assert sorted(r["committed_profiles"] + r["retired_profiles"]) == sorted(ck.GENE_KEYED_PROFILES)
    assert r["studies"] == r["committed_profiles"]
    for p, st in r["per_profile"].items():
        assert (st["outcome"] == "committed") == (p in r["committed_profiles"])
        assert (st["gates_failed"] == []) == (st["outcome"] == "committed")
        if st["outcome"] == "committed":
            assert r["calibration_on_refseq_keyed_truth"]["G0_pass"] is True
            assert st["unresolved_rate"] < g["G2_unresolved_rate_below"]
            assert st["estimated_wrong_rate"] < g["G3_estimated_wrong_rate_below"]
            lo, hi = g["G1_median_of_site_medians_within"]
            assert lo <= st["median_of_site_medians"] <= hi
        # every unresolved or ambiguous gene key is listed with its entity count
        keys = st["gene_keys"]
        assert len(keys["unresolved_listed"]) == keys["unresolved"]
        assert len(keys["ambiguous_listed"]) == keys["ambiguous"]
        assert sum(v[1] for v in keys["unresolved_listed"].values()) == keys["entities_with_unresolved_key"]
        assert len(st["unresolved_entities"]) == st["sites"]["unresolved"]
    assert r["off_acceptor_joined"] == []
    for acc, rows in list(r["sites"].items())[:200]:
        for _pos, residue, per_study in rows:
            assert residue in "STY"
            for i, n, med in per_study:
                assert isinstance(n, int) and n >= 1 and 0 <= i < len(r["studies"])
                assert isinstance(med, float)
        got = ptm.relative_abundance_in_tumours(acc, rows[0][0], r)
        assert got and got["meaning"] == ptm.CPTAC_MEANING


def test_the_apply_path_on_a_synthetic_profile(tmp_path, monkeypatch):
    import sys
    from pathlib import Path

    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
    import cptac_genekeyed as g

    prof = "ucec_cptac_2020_phosphoproteome"
    ents = [
        {"stableId": sid, "genericEntityMetaProperties": {"GENE_SYMBOL": sym}}
        for sid, sym in (
            ("MLLT4_S7", "MLLT4"),
            ("MLLT4_S35", "MLLT4"),
            ("MLLT4_Y2", "MLLT4"),
            ("NOPE_S7", "NOPE"),
        )
    ]
    (tmp_path / f"meta_{prof}.json").write_text(json.dumps(ents))
    monkeypatch.setattr(g, "CACHE", tmp_path)
    compiled = {"P55196": {"gene": "AFDN", "sequence": CANON, "sites": [(7, "phospho")], "ensg": "ENSG02"}}
    res = g.Resolver(compiled, _hgnc(), {"ENSG02": ["MDEFGHIKL" + CANON]})
    st = g.run_profile(prof, res, {"MLLT4_S7": ("P55196", 7)})
    keep = st.pop("_kept")
    assert sorted(x[2] for x in keep) == [7, 26]
    assert st["gene_keys"]["unresolved_listed"] == {"NOPE": ["symbol unknown to HGNC", 1]}
    assert st["gene_keys"]["resolved_how"] == {"previous symbol, UniProt accession": 1}
    assert st["sites"]["resolved"] == 2 and st["sites"]["unresolved"] == 1
    assert st["unresolved_entities"] == ["MLLT4_Y2"]
    assert st["sites_with_resolved_key"] == 3 and st["unresolved_rate"] == round(1 / 3, 4)
    assert st["against_the_canonical_only_mapping_of_2026_09_28"] == {"same canonical site": 1}


def _result(studies: list[str], sites: dict) -> dict:
    return {"studies": studies, "sites": sites}


def test_the_tiers_reader_labels_each_study_row_and_refuses_a_residue_disagreement():
    a = _result(["luad"], {"P1": [[5, "S", [[0, 30, -0.1]]]]})
    b = _result(
        ["brca", "coad"], {"P1": [[5, "S", [[0, 12, 0.2], [1, 40, 0.0]]]], "P2": [[9, "T", [[1, 3, 1.0]]]]}
    )
    got = ptm.tumour_abundance_tiers("P1", 5, {"refseq-keyed": a, "gene-keyed": b})
    assert got["residue"] == "S" and got["meaning"] == ptm.CPTAC_MEANING
    rows = got[ptm.CPTAC_FIELD]
    assert [(r["tier"], r["study"], r["tumours_with_value"]) for r in rows] == [
        ("refseq-keyed", "luad", 30),
        ("gene-keyed", "brca", 12),
        ("gene-keyed", "coad", 40),
    ]
    assert ptm.tumour_abundance_tiers("P2", 9, {"refseq-keyed": a, "gene-keyed": b})[ptm.CPTAC_FIELD][0][
        "tier"
    ] == ("gene-keyed")
    assert ptm.tumour_abundance_tiers("P3", 1, {"refseq-keyed": a, "gene-keyed": None}) is None
    bad = _result(["x"], {"P1": [[5, "T", [[0, 1, 0.0]]]]})
    try:
        ptm.tumour_abundance_tiers("P1", 5, {"refseq-keyed": a, "gene-keyed": bad})
    except ValueError:
        pass
    else:
        raise AssertionError("a residue disagreement between tiers must not be merged")


def test_the_committed_tiers_agree_and_share_no_study():
    a, b = ptm.load_tumour_abundance(), ptm.load_tumour_abundance(ptm.CPTAC_GENEKEYED_RESULT)
    if not a or not b:
        return
    assert not set(a["studies"]) & set(b["studies"])
    ka = {(acc, r[0]): r[1] for acc, rows in a["sites"].items() for r in rows}
    kb = {(acc, r[0]): r[1] for acc, rows in b["sites"].items() for r in rows}
    shared = set(ka) & set(kb)
    assert shared and all(ka[k] == kb[k] for k in shared)
    acc, pos = sorted(shared)[0]
    got = ptm.tumour_abundance_tiers(acc, pos)
    assert {r["tier"] for r in got[ptm.CPTAC_FIELD]} == {"refseq-keyed", "gene-keyed"}


def test_census_of_the_committed_gene_keyed_result():
    # pinned 2026-10-10 from data/results/cptac_genekeyed_result.json (lane-cptac); a re-run that moves
    # any of these must say why
    if not ck.RESULT.exists():
        return
    r = json.loads(ck.RESULT.read_text())
    assert r["retired_profiles"] == ["brca_tcga_pan_can_atlas_2018_phosphoprotein_quantification"]
    assert len(r["committed_profiles"]) == 6
    assert r["joined_committed"] == 19813 and r["proteins_with_joined_sites"] == 5223
    cal = r["calibration_on_refseq_keyed_truth"]
    assert (cal["pooled"]["correct"], cal["pooled"]["wrong"]) == (79929, 256)
