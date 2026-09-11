import pytest

from genomeos.genome.lookup import features_at, normalise_record, parse_variant

RAW = {
    "most_severe_consequence": "missense_variant",
    "transcript_consequences": [
        {
            "gene_symbol": "TP53",
            "transcript_id": "ENST00000269305",
            "hgvsp": "ENSP00000269305.4:p.Arg175His",
            "hgvsc": "ENST00000269305.9:c.524G>A",
            "sift_prediction": "deleterious",
            "polyphen_prediction": "probably_damaging",
        }
    ],
    "colocated_variants": [
        {"id": "COSV52661038"},
        {
            "id": "rs28934578",
            "clin_sig": ["likely_pathogenic", "pathogenic"],
            "pubmed": [26900293, 34239995],
            "frequencies": {"T": {"gnomadg": 6.57e-06}},
        },
    ],
}


def test_parse_variant_forms():
    assert parse_variant("chr21:25897620 C>T") == ("chr21", 25897620, "C", "T")
    assert parse_variant("21 25,897,620 C T") == ("chr21", 25897620, "C", "T")
    assert parse_variant("MT:8993 T>G") == ("chrM", 8993, "T", "G")
    with pytest.raises(ValueError):
        parse_variant("APP A673T")


def test_normalise_record_keeps_what_is_known():
    n = normalise_record(RAW)
    assert n["gene"] == "TP53" and n["hgvsp"] == "p.Arg175His" and n["residue"] == 175
    assert n["aa_from"] == "R" and n["aa_to"] == "H"
    assert n["clinvar"] == ["likely_pathogenic", "pathogenic"] and n["dbsnp"] == ["rs28934578"]
    assert n["cosmic"] == ["COSV52661038"] and n["pubmed_count"] == 2 and n["gnomad_af"] == 6.57e-06


def test_features_at_residue():
    defn = {
        "sections": {
            "domains": {
                "items": {
                    "interpro": [],
                    "features": [
                        {"type": "Domain", "description": "DNA-binding", "start": 100, "end": 300},
                        {"type": "Region", "description": "Tetramer", "start": 320, "end": 360},
                    ],
                }
            },
            "modifications": {
                "items": [
                    {"type": "Modified residue", "description": "Phosphoserine", "start": 175, "end": 175}
                ]
            },
        }
    }
    hits = features_at(defn, 175)
    assert [h["type"] for h in hits] == ["Domain", "Modified residue"]
    assert features_at(defn, 500) == []
