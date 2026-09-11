# SPDX-License-Identifier: AGPL-3.0-or-later
"""Post-translational state read out of compiled definitions (no network)."""

import json

from genomeos.molecules.ptm import build_index, classify, sites, states, summary, writers_of


def test_writers_and_classes():
    assert writers_of("Phosphoserine; by CDK5, PRPK, AMPK, NUAK1 and ATM") == [
        "CDK5",
        "PRPK",
        "AMPK",
        "NUAK1",
        "ATM",
    ]
    assert writers_of("Phosphothreonine; by autocatalysis") == ["(self)"]
    assert writers_of("N6-acetyllysine") == []
    assert classify("Phosphoserine; by ATM", "Modified residue") == "phospho"
    assert classify("N6-acetyllysine; by KAT5", "Modified residue") == "acetyl"
    assert (
        classify("Glycyl lysine isopeptide (Lys-Gly) (interchain with G-Cter in ubiquitin)", "Cross-link")
        == "ubiquitin"
    )
    assert classify("N-linked (GlcNAc...) asparagine", "Glycosylation") == "glyco"


def test_sites_states_and_index(tmp_path):
    defn = {
        "gene": "TP53",
        "id": "UniProt:P04637",
        "sections": {
            "identity": {"items": {"accession": "P04637"}},
            "modifications": {
                "items": [
                    {
                        "type": "Modified residue",
                        "description": "Phosphoserine; by HIPK4",
                        "start": 9,
                        "end": 9,
                    },
                    {
                        "type": "Modified residue",
                        "description": "N6-acetyllysine; by KAT5 and KAT8",
                        "start": 120,
                        "end": 120,
                    },
                ]
            },
        },
    }
    s = sites(defn)
    assert [x["class"] for x in s] == ["phospho", "acetyl"] and s[1]["writers"] == ["KAT5", "KAT8"]
    st = states(defn)
    assert st[0]["modifications"] == ["Phosphoserine at 9"] and st[0]["written_by"] == ["HIPK4"]
    (tmp_path / "TP53.json").write_text(json.dumps(defn))
    (tmp_path / "EMPTY.json").write_text(
        json.dumps({"gene": "EMPTY", "sections": {"identity": {"items": {}}}})
    )
    idx = build_index(tmp_path, tmp_path / "_ptm_index.json")
    assert idx["proteins"] == 2 and idx["proteins_with_sites"] == 1 and idx["sites"] == 2
    assert idx["writers"]["KAT5"] == {"TP53": 1} and idx["by_class"] == {"phospho": 1, "acetyl": 1}
    sm = summary(idx)
    assert sm["writers"] == 3 and sm["writer_edges"] == 3 and sm["top_writers"][0]["substrates"] == 1
