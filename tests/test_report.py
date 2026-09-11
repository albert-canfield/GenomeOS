from pathlib import Path

import pytest

from genomeos.report import to_markdown


def test_markdown_renders_every_section_present():
    rep = {
        "gene": "X",
        "chrom": "chr1",
        "sections": {
            "origin": {
                "locus": "chr1:1-100(+)",
                "strand": "+",
                "length_bp": 100,
                "type": "protein_coding",
                "transcripts": 2,
                "coding_transcripts": 1,
                "evidence": "curated: GENCODE",
            },
            "flow": {
                "canonical": "X-201",
                "exons": 3,
                "mrna_nt": 90,
                "utr5": 10,
                "cds_nt": 60,
                "utr3": 20,
                "protein_aa": 19,
                "codon_table": "standard",
                "identity_with_uniprot": 1.0,
                "uniprot_aa": 19,
                "evidence": "derived",
            },
            "individual": {
                "sample": "HG002",
                "variants_in_gene": 3,
                "coding_snvs": [
                    {
                        "pos": 5,
                        "change": "A>G",
                        "genotype": "0/1",
                        "consequence": "missense",
                        "hgvs_p": "p.Lys2Arg",
                    }
                ],
                "evidence": "measured",
            },
        },
    }
    md = to_markdown(rep)
    assert md.startswith("# X (chr1)") and "**Origin**" in md and "Identity with UniProt 100.0%" in md
    assert "p.Lys2Arg (0/1)" in md and "_measured_" in md
    assert "**Regulation**" not in md  # sections absent from the dict are not rendered


@pytest.mark.skipif(not Path("data/reference/chr21.fa.gz").exists(), reason="needs chr21")
def test_gene_report_app():
    from genomeos.genome import Annotation, IndexedGenome, default_gencode
    from genomeos.report import gene_report

    ann = Annotation.from_gff3(default_gencode({"chr21"}), {"chr21"})
    genome = IndexedGenome("data/reference/chr21.fa.gz")
    try:
        rep = gene_report("APP", "chr21", ann, genome)
    finally:
        genome.close()
    s = rep["sections"]
    assert s["origin"]["transcripts"] >= 10 and s["flow"]["protein_aa"] == 770
    assert "regulation" in s and s["regulation"]["enhancers_in_node"] > 0
    if "protein" in s:
        assert len(s["protein"]["partners_physical"]) == len(set(s["protein"]["partners_physical"]))
