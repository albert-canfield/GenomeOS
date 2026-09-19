# SPDX-License-Identifier: AGPL-3.0-or-later
"""Node comparison across species by symbol (offline)."""

import json

from genomeos.genome.mouse import compare_nodes, human_nodes_by_symbol


def test_compare_nodes_by_symbol(tmp_path):
    (tmp_path / "domains_chr21.json").write_text(
        json.dumps(
            {
                "domains": [
                    {"id": "chr21:D1", "genes": ["APP", "APP-DT"]},
                    {"id": "chr21:D2", "genes": ["SOD1", "SFRS15"]},
                ]
            }
        )
    )
    index, n = human_nodes_by_symbol(tmp_path)
    assert n == 2 and index["APP"] == "chr21:D1" and index["SOD1"] == "chr21:D2"
    index["FAR"] = "chr5:D9"
    mouse = [
        {"id": "chr16:D1", "start": 0, "end": 10, "coding_symbols": ["App", "Sod1"]},  # adjacent human nodes
        {"id": "chr16:D2", "start": 10, "end": 20, "coding_symbols": ["Sod1", "Sfrs15"]},  # one node
        {"id": "chr16:D3", "start": 20, "end": 30, "coding_symbols": ["App", "Gm1234"]},  # one match only
        {"id": "chr16:D4", "start": 30, "end": 40, "coding_symbols": []},
        {"id": "chr16:D5", "start": 40, "end": 50, "coding_symbols": ["App", "Far"]},  # another chromosome
    ]
    c = compare_nodes(mouse, index)
    assert (c["conserved"], c["split_adjacent"], c["split_scattered"], c["unmapped"]) == (1, 1, 1, 2)
    assert c["tested"] == 3 and c["fraction_conserved"] == 0.333 and c["fraction_same_neighbourhood"] == 0.667
    verdicts = [r["verdict"] for r in c["rows"]]
    assert verdicts == ["split_adjacent", "conserved", "unmapped", "unmapped", "split_scattered"]
    assert c["orthology"].startswith("inferred")
    # a curated orthology maps a renamed mouse gene onto its human orthologue
    c2 = compare_nodes(mouse, index, {"Gm1234": ["APP-DT"], "Far": ["SOD1"]})
    assert c2["orthology"].startswith("curated")
    assert [r["verdict"] for r in c2["rows"]] == [
        "split_adjacent",
        "conserved",
        "conserved",
        "unmapped",
        "split_adjacent",
    ]


def test_compara_orthology_from_a_small_dump(tmp_path, monkeypatch):
    import gzip

    from genomeos.genome import mouse

    # human gene table and a mouse GENCODE file the resolver reads
    genes = tmp_path / "genes.tsv"
    genes.write_text(
        "gene_id\tsymbol\tgene_type\tchrom\tstart\tend\tstrand\nENSG1.2\tAPP\tprotein_coding\tchr21\t1\t2\t+\n"
    )
    from genomeos.knowledge import homology

    monkeypatch.setattr(homology, "GENES_TSV", genes)
    ref = tmp_path / "ref"
    ref.mkdir()
    with gzip.open(ref / "gencode_vM25_chr16.gff3.gz", "wt") as fh:
        fh.write("##gff-version 3\n")
        fh.write("chr16\tHAVANA\tgene\t1\t2\t.\t+\t.\tID=ENSMUSG1.3;gene_id=ENSMUSG1.3;gene_name=App\n")
    monkeypatch.setattr(mouse, "REFERENCE", ref)
    header = ["gene_stable_id", "homology_type", "homology_gene_stable_id", "homology_species"]
    rows = iter(
        [
            header,
            ["ENSMUSG1", "ortholog_one2one", "ENSG1", "homo_sapiens"],
            ["ENSMUSG1", "ortholog_one2one", "ENSRNOG1", "rattus_norvegicus"],
            ["ENSMUSG1", "within_species_paralog", "ENSMUSG9", "mus_musculus"],
            ["ENSMUSG1", "ortholog_one2one", "ENSG7", "homo_sapiens"],  # human id unknown locally
        ]
    )
    p = mouse.fetch_compara_orthology(tmp_path, rows=rows)
    with gzip.open(p, "rt") as fh:
        text = fh.read()
    assert "2 human rows" in text and "1 human ids" in text
    assert mouse.load_orthology(tmp_path, source="compara") == {"App": ["APP"]}
    agreement = mouse.orthology_agreement(
        {"App": ["APP"], "Sod1": ["SOD1"]}, {"App": ["APP"]}, ["App", "Sod1", "X"]
    )
    assert agreement == {
        "genes": 3,
        "in_both": 1,
        "identical_when_in_both": 1.0,
        "mgi_only": 1,
        "compara_only": 0,
        "neither": 1,
    }
