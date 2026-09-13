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


def test_carriers_reads_a_local_individual(tmp_path, monkeypatch):
    from genomeos.genome import lookup as lk

    vcf = tmp_path / "HG002_chr21.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHG002\n"
        "chr21\t100\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\nchr21\t200\t.\tC\tT,G\t50\tPASS\t.\tGT\t1/1\n"
    )
    monkeypatch.setattr(lk, "INDIVIDUALS", {"HG002": (str(tmp_path / "HG002_{chrom}.vcf"), "test calls")})
    monkeypatch.setattr(lk, "INDIVIDUAL_FALLBACK", str(tmp_path / "none_{sample}_{chrom}.vcf"))
    from genomeos.genome import individuals as ind

    monkeypatch.setattr(ind, "ROOT", tmp_path / "nobody")  # imported genomes on this machine stay out
    assert lk.carriers("chr21", 100, "A", "G")[0]["genotype"] == "0/1"
    assert lk.carriers("chr21", 200, "C", "G")[0]["carries"] is True
    assert lk.carriers("chr21", 300, "A", "G")[0]["carries"] is False
    assert lk.carriers("chr22", 100, "A", "G") == []
