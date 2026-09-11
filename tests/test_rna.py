from pathlib import Path

import pytest

from genomeos.molecules.rna import summarise_tissues, transcripts_of


def test_tissue_summary_patterns():
    flat = {f"t{i}": 10.0 for i in range(20)}
    s = summarise_tissues("X", "ENSG1.1", flat)
    assert (
        s["pattern"] == "expressed in all tissues"
        and s["median_tpm"] == 10.0
        and s["tissues_expressed"] == 20
    )
    enhanced = {**flat, "Brain": 200.0}
    assert summarise_tissues("X", "ENSG1.1", enhanced)["pattern"] == "tissue enhanced"
    assert summarise_tissues("X", "ENSG1.1", enhanced)["top"][0] == ("Brain", 200.0)
    restricted = {f"t{i}": 0.1 for i in range(20)} | {"Testis": 50.0}
    assert summarise_tissues("X", "ENSG1.1", restricted)["pattern"] == "tissue restricted"
    assert summarise_tissues("X", "ENSG1.1", {})["pattern"] == "not detected"


@pytest.mark.skipif(
    not Path("data/results/gencode_v50_chr21_chrM.gff3.gz").exists(), reason="needs chr21 models"
)
def test_transcripts_of_app():
    from genomeos.genome import Annotation, IndexedGenome, default_gencode

    ann = Annotation.from_gff3(default_gencode({"chr21"}), {"chr21"})
    genome = (
        IndexedGenome("data/reference/chr21.fa.gz") if Path("data/reference/chr21.fa.gz").exists() else None
    )
    r = transcripts_of("APP", ann, genome)
    assert r["count"] >= 10 and r["coding_isoforms"] >= 5
    canon = r["transcripts"][0]
    assert canon["canonical"] and canon["name"] == "APP-201" and canon["exons"] == 18
    if genome:
        assert canon["spliced_nt"] == 3583 and canon["protein_aa"] == 770
        genome.close()
    assert "protein_coding" in r["by_biotype"]


def test_summarise_isoforms_names_the_dominant_transcript_per_tissue():
    from genomeos.molecules.rna import summarise_isoforms

    rows = [
        {"transcriptId": "ENST1.1", "tissueSiteDetailId": "Brain", "median": 1.0},
        {"transcriptId": "ENST2.3", "tissueSiteDetailId": "Brain", "median": 9.0},
        {"transcriptId": "ENST1.1", "tissueSiteDetailId": "Liver", "median": 5.0},
        {"transcriptId": "ENST2.3", "tissueSiteDetailId": "Liver", "median": 0.0},
        {"transcriptId": "ENST1.1", "tissueSiteDetailId": "Empty", "median": 0.0},
    ]
    tx = [
        {
            "transcript": "ENST1",
            "name": "G-201",
            "canonical": True,
            "protein_length": 100,
            "biotype": "protein_coding",
        },
        {
            "transcript": "ENST2",
            "name": "G-202",
            "canonical": False,
            "protein_length": 80,
            "biotype": "protein_coding",
        },
    ]
    r = summarise_isoforms("G", "ENSG1.1", rows, tx)
    assert r["transcripts_measured"] == 2 and r["tissues_measured"] == 3 and r["canonical"] == "ENST1"
    assert r["canonical_dominant_in"] == 1 and list(r["tissues_where_another_isoform_dominates"]) == ["Brain"]
    assert r["tissues_where_another_isoform_dominates"]["Brain"]["share"] == 0.9
    assert list(r["isoforms"])[0] in ("ENST1", "ENST2") and r["isoforms"]["ENST2"]["tissues_dominant"] == 1
    assert "Empty" not in r["dominant_by_tissue"]  # nothing expressed there, nothing dominates
