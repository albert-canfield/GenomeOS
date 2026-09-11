import gzip
import json

from genomeos.lib import proteome


def _defn(gene, acc, length=100):
    return {
        "gene": gene,
        "sections": {
            "identity": {
                "items": {
                    "accession": acc,
                    "name": f"{gene} protein",
                    "length": length,
                    "existence": "1: Evidence at protein level",
                    "symbol_match": True,
                }
            },
            "function": {"items": {"summary": ["Does X. " * 60], "location": ["Nucleus", "Cytoplasm"]}},
            "domains": {"items": {"interpro": [{"id": "IPR1", "name": "Dom1"}], "features": []}},
            "pathways": {"items": [{"id": "R-HSA-1", "name": "P"}]},
            "interactions": {
                "items": [
                    {"partner": "B", "score": 0.9, "physical_evidence": True},
                    {"partner": "C", "score": 0.8, "physical_evidence": False},
                ]
            },
            "expression": {"items": {"tissue_specificity": "Tissue enhanced"}},
            "diseases": {"items": [{"name": "Syndrome Y"}]},
            "structures_experimental": {"items": [{"id": "1ABC"}]},
            "structures_predicted": {"items": [{"id": acc}]},
        },
    }


def test_distil_load_block(tmp_path, monkeypatch):
    cache = tmp_path / "cache"
    cache.mkdir()
    (cache / "A.json").write_text(json.dumps(_defn("A", "P1")))
    (cache / "Z.json").write_text(json.dumps({"gene": "Z", "sections": {"identity": {"items": None}}}))
    out = tmp_path / "proteome.json.gz"
    monkeypatch.chdir(tmp_path)  # no data/results here: chrom stays None
    s = proteome.distil(cache, out)
    assert s["proteins"] == 1 and out.exists()
    monkeypatch.setattr(proteome, "_LOADED", {})
    monkeypatch.setattr(proteome, "PACKAGED", out)
    t = proteome.load(out)
    r = t["A"]
    assert r["accession"] == "P1" and r["partners"] == ["B"] and r["function"].endswith("…")
    assert r["structures_experimental"] == 1 and r["alphafold"] and r["existence"] == "1"
    b = proteome.block("a")
    assert b.startswith("# A protein (100 aa)") and "protein A {" in b and "accession: P1;" in b
    assert "interactions: B;" in b and "confidence: 0.9" in b
    assert "writers" in r
    assert proteome.block("nope").startswith("# NOPE: not in")
    with gzip.open(out, "rt") as fh:
        assert json.load(fh)["evidence"]["pathways"].startswith("curated")


def test_packaged_table_ships_with_the_code():
    t = proteome.load()
    assert len(t) > 19000 and t["TP53"]["accession"] == "P04637"
    assert proteome.summary()["proteins"] == len(t)
