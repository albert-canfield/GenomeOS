"""Blocks of a chromosome window, nested and to scale, with plausibility of moves."""

from pathlib import Path

import pytest

from genomeos.genome import Annotation, Genome, default_gencode
from genomeos.genome.blocks import blocks_for_window, check_move

CHRM = Path("data/reference/chrM.fa.gz")


@pytest.mark.skipif(not (CHRM.exists() and default_gencode({"chrM"})), reason="fetch chrM first")
def test_mitochondrial_blocks_nest_and_cover():
    ann = Annotation.from_gff3(default_gencode({"chrM"}), {"chrM"})
    seq = Genome.from_fasta(CHRM).chromosomes["chrM"].sequence
    r = blocks_for_window(ann, "chrM", 0, len(seq), seq, len(seq))
    c = r["counts"]
    assert c["gene"] == 37 and c["transcript"] == 37 and c["cds"] == 13
    by = {b["id"]: b for b in r["blocks"]}
    for b in r["blocks"]:
        if b["parent"]:
            p = by[b["parent"]]
            assert p["start"] <= b["start"] and b["end"] <= p["end"], (b["name"], p["name"])
    assert all(b["evidence"] in ("curated", "predicted", "none") for b in r["blocks"])
    unknown = [b for b in r["blocks"] if b["type"] == "unknown"]
    assert unknown and all(
        b["confidence"] == 0.0 for b in unknown
    )  # the D-loop / control region is unknown to us


@pytest.mark.skipif(not (CHRM.exists() and default_gencode({"chrM"})), reason="fetch chrM first")
def test_move_plausibility():
    ann = Annotation.from_gff3(default_gencode({"chrM"}), {"chrM"})
    seq = Genome.from_fasta(CHRM).chromosomes["chrM"].sequence
    r = blocks_for_window(ann, "chrM", 0, len(seq), seq, len(seq))
    genes = [b for b in r["blocks"] if b["type"] == "gene" and b["strand"] == "+"]
    co1 = next(b for b in genes if b["name"] == "MT-CO1")
    co2 = next(b for b in genes if b["name"] == "MT-CO2")
    # sliding CO1 onto CO2 (same strand) is not plausible
    v = check_move(r["blocks"], co1["id"], co2["start"] - 100, len(seq))
    assert not v["ok"] and any("overlaps" in x for x in v["reasons"])
    # moving off the chromosome is not plausible
    assert not check_move(r["blocks"], co1["id"], len(seq) - 10, len(seq))["ok"]
    # an exon may not leave its transcript; a shift by a non-multiple of 3 breaks the frame
    exon = next(
        b
        for b in r["blocks"]
        if b["type"] == "exon" and b["parent"].startswith("ENST") and b["name"] == "exon 1"
    )
    v = check_move(r["blocks"], exon["id"], exon["start"] + 1, len(seq))
    assert not v["ok"] and any("frame" in x or "leaves" in x for x in v["reasons"])


def test_reader_marks_nodes_and_genes(tmp_path):
    import json

    from genomeos.genome.blocks import Block, _apply_reader

    (tmp_path / "reader_K562_chr21.json").write_text(
        json.dumps(
            {
                "node_table": [
                    {"id": "chr21:D1", "open_fraction": 0.03, "peaks": 40},
                    {"id": "chr21:D2", "open_fraction": 0.0, "peaks": 0},
                ],
                "silent_node_ids": ["chr21:D2"],
                "silent_genes": ["TPTE"],
            }
        )
    )
    blocks = [
        Block("chr21:D1", "domain", 0, 1000, ".", "node D1", None, "inferred", 0.4, {}),
        Block("chr21:D2", "domain", 1000, 2000, ".", "node D2", None, "inferred", 0.4, {}),
        Block("g1", "gene", 10, 500, "+", "APP", None, "curated", 0.9, {"gene_type": "protein_coding"}),
        Block("g2", "gene", 1100, 1500, "-", "TPTE", None, "curated", 0.9, {"gene_type": "protein_coding"}),
        Block("g3", "gene", 1600, 1700, "-", "LINC1", None, "curated", 0.9, {"gene_type": "lncRNA"}),
    ]
    assert _apply_reader(blocks, "chr21", tmp_path) == ["K562"]
    assert blocks[0].attrs["K562_open_fraction"] == 0.03 and blocks[0].attrs["K562_node"] == "open"
    assert blocks[1].attrs["K562_node"] == "silent" and blocks[1].attrs["K562_peaks"] == 0
    assert blocks[2].attrs["K562_read"] == "read" and blocks[3].attrs["K562_read"] == "silent"
    assert "K562_read" not in blocks[4].attrs  # the reader speaks about coding genes only
    assert _apply_reader(blocks, "chr22", tmp_path) == []


def test_a_poised_gene_is_not_drawn_as_silent(tmp_path):
    """Poised arrives INSIDE silent_genes and must still be told apart from silent.

    The reader puts poised genes in silent_genes on purpose: every consumer reads "absent from this
    list" as read, so leaving them out would report a held gene as expressed. But poised is H3K27me3
    over a promoter that also carries H3K4me3 - a gene held ready, not one shut - and the Blocks lane
    drew all of them as silent until 2026-09-17. In one 4 Mb window of chr21 that was 23 genes.
    """
    import json

    from genomeos.genome.blocks import Block, _apply_reader

    (tmp_path / "reader_K562_chr21.json").write_text(
        json.dumps(
            {
                "node_table": [{"id": "chr21:D1", "open_fraction": 0.5, "peaks": 10}],
                "silent_node_ids": [],
                "silent_genes": ["TPTE", "BACH1"],  # the poised one is in here too, by design
                "poised_genes": ["BACH1"],
            }
        )
    )
    blocks = [
        Block("g1", "gene", 10, 500, "+", "APP", None, "curated", 0.9, {"gene_type": "protein_coding"}),
        Block("g2", "gene", 600, 900, "-", "TPTE", None, "curated", 0.9, {"gene_type": "protein_coding"}),
        Block("g3", "gene", 950, 990, "+", "BACH1", None, "curated", 0.9, {"gene_type": "protein_coding"}),
    ]

    _apply_reader(blocks, "chr21", tmp_path)

    assert blocks[0].attrs["K562_read"] == "read"
    assert blocks[1].attrs["K562_read"] == "silent"
    assert blocks[2].attrs["K562_read"] == "poised"


def test_a_reader_result_without_poised_genes_still_works(tmp_path):
    """Results written before poised_genes existed must not start reporting genes as read."""
    import json

    from genomeos.genome.blocks import Block, _apply_reader

    (tmp_path / "reader_K562_chr21.json").write_text(
        json.dumps({"node_table": [], "silent_node_ids": [], "silent_genes": ["TPTE"]})
    )
    blocks = [
        Block("g2", "gene", 600, 900, "-", "TPTE", None, "curated", 0.9, {"gene_type": "protein_coding"})
    ]

    _apply_reader(blocks, "chr21", tmp_path)

    assert blocks[0].attrs["K562_read"] == "silent"


def test_the_two_axes_reach_the_block_as_a_case(tmp_path):
    """Area J computed the case per block for weeks and it appeared nowhere until 2026-09-17.

    The result carries the pair and its confidence as an object, not a bare label, so a reader that
    took `hit["case"]` straight got a dict and showed nothing.
    """
    import json

    from genomeos.genome.blocks import Block, _apply_variation_case

    (tmp_path / "variation_chr21.json").write_text(
        json.dumps(
            {
                "blocks": [
                    {
                        "start": 0,
                        "end": 1000,
                        "mammal_fraction": 0.31,
                        "gnocchi": {"mean": 2.4, "bases": 1000},
                        "case": {
                            "case": "syntax",
                            "mammals": "held",
                            "humans": "constrained",
                            "confidence": 0.6,
                        },
                    }
                ]
            }
        )
    )
    blocks = [Block("u1", "unknown", 0, 1000, ".", "UNKNOWN", None, "inferred", 0.4, {})]

    _apply_variation_case(blocks, "chr21", tmp_path)

    assert blocks[0].attrs["case"] == "syntax"
    assert blocks[0].attrs["case_mammals"] == "held"
    assert blocks[0].attrs["case_humans"] == "constrained"
    assert blocks[0].attrs["mammal_fraction"] == 0.31
    assert blocks[0].attrs["gnocchi_mean"] == 2.4
    # the base COUNT does not come along: it is the field the 2026-09-17 bin-crediting fix corrected,
    # and a ratio that survived is not a licence to carry the absolute that did not
    assert "gnocchi_bases" not in blocks[0].attrs


def test_a_block_with_no_measurement_gets_no_case(tmp_path):
    import json

    from genomeos.genome.blocks import Block, _apply_variation_case

    (tmp_path / "variation_chr21.json").write_text(
        json.dumps({"blocks": [{"start": 0, "end": 10, "case": None}]})
    )
    blocks = [Block("u1", "unknown", 0, 10, ".", "UNKNOWN", None, "inferred", 0.4, {})]

    _apply_variation_case(blocks, "chr21", tmp_path)

    assert "case" not in blocks[0].attrs
