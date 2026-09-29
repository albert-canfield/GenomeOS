# SPDX-License-Identifier: AGPL-3.0-or-later
"""CPTAC phosphosite entities: parsing, RefSeq-to-UniProt position transfer, the committed summary."""

from __future__ import annotations

from genomeos.molecules import ptm


def _e(sid: str, **props: str) -> dict:
    return {"stableId": sid, "genericEntityMetaProperties": props}


def test_luad_single_localised_site_keeps_the_refseq_accession():
    r = ptm.parse_cptac_entity(
        "luad_cptac_2020_phosphoproteome", _e("NP_000010.1_1_1_69_69", GENE_SYMBOL="ACAT1", NAME="ACAT1_S69s")
    )
    assert r == {"refseq": "NP_000010.1", "gene": "ACAT1", "residue": "S", "position": 69}


def test_luad_unlocalised_and_multi_site_entities_are_dropped_with_reasons():
    p = "luad_cptac_2020_phosphoproteome"
    assert ptm.parse_cptac_entity(p, _e("NP_113584.3_1_0_3446_3447", NAME="HUWE1_S3446s"))["drop"] == (
        "site not localised"
    )
    assert ptm.parse_cptac_entity(p, _e("NP_057269.1_3_3_87_91", NAME="JPT1_S87s_S88s_S91s"))["drop"] == (
        "several sites in one entity"
    )
    assert ptm.parse_cptac_entity(p, _e("smORF_G086960_1_1_20_20", NAME="G086960_T20t"))["drop"] == (
        "not a RefSeq protein"
    )


def test_gbm_suffix_duplicate_and_brain_and_gene_keyed_formats():
    g = "gbm_cptac_2021_phosphoproteome"
    assert (
        ptm.parse_cptac_entity(g, _e("AAAS_S495:NP_056480.1", GENE_SYMBOL="AAAS"))["refseq"] == "NP_056480.1"
    )
    assert "drop" in ptm.parse_cptac_entity(g, _e("AAAS_S495.1:NP_056480.1", GENE_SYMBOL="AAAS"))
    assert "drop" in ptm.parse_cptac_entity(g, _e("APC2_S1764_S1767:NP_005874.1", GENE_SYMBOL="APC2"))
    b = ptm.parse_cptac_entity(
        "brain_cptac_2020_phosphoprotein",
        _e("AAK1_14_34_1_1_S21", DESCRIPTION="NP_055726.3", NAME="AAK1 S21 14-34 1_1"),
    )
    assert b == {"refseq": "NP_055726.3", "gene": "AAK1", "residue": "S", "position": 21}
    u = ptm.parse_cptac_entity("ucec_cptac_2020_phosphoproteome", _e("PRR5_ARHGAP8_S555", GENE_SYMBOL="PRR5"))
    assert u == {"refseq": None, "gene": "PRR5", "residue": "S", "position": 555}
    t = ptm.parse_cptac_entity(
        "brca_tcga_phosphoprotein_quantification", _e("A2M_pS710", GENE_SYMBOL="A2M", PHOSPHOSITE="pS710")
    )
    assert t["position"] == 710 and t["refseq"] is None
    assert "drop" in ptm.parse_cptac_entity(
        "ov_tcga_phosphoprotein_quantification", _e("A2M_pT118_T119", PHOSPHOSITE="pT118_T119")
    )
    assert ptm.parse_cptac_entity("lusc_cptac_2021_phosphoproteome", _e("ZZZ3_acetylprotein"))["drop"] == (
        "gene-level aggregate, not a site"
    )


def test_transfer_carries_identical_positions_and_refuses_ambiguous_windows():
    s = "MAAAKLSPEDRRTYWGGHNVCQ"
    assert ptm.transfer_position(s, s, 7) == (7, "identical")
    # the target carries four extra residues well before the site; the window still finds it
    long = "MDEFGHIKLMNPAAAKLSPEDRRTYWGGHNVCQ"
    site = long.index("S") + 1
    t = "MDEFWWWWGHIKLMNPAAAKLSPEDRRTYWGGHNVCQ"
    pos, how = ptm.transfer_position(long, t, site)
    assert how == "window" and pos == site + 4 and t[pos - 1] == "S"
    assert ptm.transfer_position(s, "MXXXXXXXXXXXXX", 7)[0] is None
    rep = "KLSPEDRRTY" * 3
    assert ptm.transfer_position("KLSPEDRRTY", rep, 3, flank=2)[1].endswith(
        "more than once in the UniProt sequence"
    )
    assert ptm.transfer_position(s, s, 99)[0] is None


def test_the_committed_summary_carries_no_per_patient_value_and_reads_back():
    r = ptm.load_tumour_abundance()
    if r is None:
        return
    assert r["field"] == ptm.CPTAC_FIELD
    for rows in list(r["sites"].values())[:200]:
        for _pos, residue, per_study in rows:
            assert residue in "STY"
            for i, n, med in per_study:
                assert isinstance(n, int) and n >= 1 and 0 <= i < len(r["studies"])
                assert isinstance(med, float)
    acc, rows = next(iter(r["sites"].items()))
    got = ptm.relative_abundance_in_tumours(acc, rows[0][0], r)
    assert got and got["meaning"] == ptm.CPTAC_MEANING
