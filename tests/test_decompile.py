"""The decompiled locus assembles saved layers and names the absent ones (area J step 6)."""

import json

from genomeos.decompile import LAYERS, _elements, _node, _reader, render


def test_layers_from_results_and_absence(tmp_path):
    node = {
        "id": "chrT:D1",
        "start": 0,
        "end": 5000,
        "length": 5000,
        "genes": ["G1"],
        "coding_genes": 1,
        "promoters": 1,
        "enhancers": 3,
        "evidence": "inferred",
    }
    (tmp_path / "domains_chrT.json").write_text(json.dumps({"domains": [node]}))
    (tmp_path / "reader_K562_chrT.json").write_text(json.dumps({"cell_type": "K562", "silent_genes": ["G2"]}))
    (tmp_path / "reader_HepG2_chrT.json").write_text(
        json.dumps({"cell_type": "HepG2", "silent_genes": ["G1"]})
    )
    (tmp_path / "reader_K562_vs_HepG2_chrT.json").write_text(json.dumps({"silent_genes": []}))
    pc = {"gene": "G1", "action": "activates", "log2_fold_change": -0.4, "tissue": "liver", "confidence": 0.4}
    el = {"id": "E1", "start": 100, "end": 300, "predicted_coding": pc}
    other = {**el, "id": "E2", "predicted_coding": {"gene": "G9"}}
    (tmp_path / "enhancer_targets_chrT.json").write_text(json.dumps({"elements": [el, other]}))
    var = {"id": "E1", "mammal_fraction": 0.3, "gnocchi": {"fraction_above": 0.5}, "case": {"case": "syntax"}}
    (tmp_path / "variation_chrT.json").write_text(json.dumps({"elements": [var]}))
    (tmp_path / "duplication_chrT.json").write_text(json.dumps({"elements": {"ids": ["E1"]}}))
    (tmp_path / "motifs_chrT.json").write_text(
        json.dumps({"genes": {}, "elements": [{"id": "E1", "requires": [{"factor": "GATA4"}]}]})
    )
    assert _node("G1", None, "chrT", tmp_path)["id"] == "chrT:D1"
    assert _node("G1", {"start": 100}, "chrT", tmp_path)["id"] == "chrT:D1"
    assert _reader("G1", "chrT", tmp_path) == {"K562": "read", "HepG2": "silent"}
    els = _elements("G1", "chrT", tmp_path)
    assert [e["id"] for e in els] == ["E1"]
    e = els[0]
    assert (
        e["case"] == "syntax"
        and e["duplicated"]
        and e["requires"] == ["GATA4"]
        and e["human_fraction"] == 0.5
    )
    gene = {
        "symbol": "G1",
        "type": "protein_coding",
        "start": 100,
        "end": 900,
        "strand": "+",
        "transcripts": 1,
        "canonical": "T1",
        "exons": 3,
        "cds_segments": 3,
    }
    layers = {k: None for k in LAYERS}
    layers.update(
        {"gene": gene, "elements": els, "human_axis": els, "duplication": ["E1"], "reader": {"K562": "read"}}
    )
    d = {
        "symbol": "G1",
        "chrom": "chrT",
        "layers": layers,
        "unknown": ["protein", "origin", "paralogues", "requires", "node", "expression"],
    }
    text = render(d)
    assert "gene G1 {" in text and "element E1 {" in text
    assert "case syntax; inside a segmental duplication" in text
    assert "read in K562; silent in none" in text and "origin: not read for this locus" in text
    assert text.count("{") == text.count("}")


def test_a_poised_cell_type_is_named_and_not_dropped(tmp_path):
    """Before 2026-09-17 a poised cell appeared in neither group and simply vanished from the line.

    That is the worst of the three outcomes: silent would at least be a visible claim, and absent
    reads as "this cell type was never measured".
    """
    import json

    from genomeos.decompile import _reader

    (tmp_path / "reader_K562_chr21.json").write_text(
        json.dumps(
            {
                "cell_type": "K562",
                "silent_genes": ["TPTE", "BACH1"],  # poised genes are in here by design
                "poised_genes": ["BACH1"],
            }
        )
    )

    out = _reader("BACH1", "chr21", tmp_path)
    assert out == {"K562": "poised"}
    assert _reader("TPTE", "chr21", tmp_path) == {"K562": "silent"}
    assert _reader("APP", "chr21", tmp_path) == {"K562": "read"}


def test_a_result_without_poised_genes_reads_as_it_did(tmp_path):
    import json

    from genomeos.decompile import _reader

    (tmp_path / "reader_K562_chr21.json").write_text(
        json.dumps({"cell_type": "K562", "silent_genes": ["TPTE"]})
    )

    assert _reader("TPTE", "chr21", tmp_path) == {"K562": "silent"}
