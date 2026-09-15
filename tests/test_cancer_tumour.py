from genomeos.cancer.tumour import (
    GERMLINE_AF,
    TumourVariant,
    grade,
    mutant_peptides,
    mutation_burden,
    normalise_vep,
)
from genomeos.genome.variants import Variant

VEP_TP53 = {
    "input": "17 7675088 . C T . . .",
    "most_severe_consequence": "missense_variant",
    "transcript_consequences": [
        {
            "gene_symbol": "TP53",
            "hgvsp": "ENSP00000269305.4:p.Arg175His",
            "hgvsc": "ENST00000269305.9:c.524G>A",
            "sift_prediction": "deleterious",
            "polyphen_prediction": "probably_damaging",
        }
    ],
    "colocated_variants": [
        {"id": "COSV52661038"},
        {"id": "rs28934578", "frequencies": {"T": {"gnomadg": 6.57e-06, "gnomade_afr": 0}}},
    ],
}
VEP_COMMON = {
    "input": "21 1 . A G . . .",
    "most_severe_consequence": "synonymous_variant",
    "transcript_consequences": [{"gene_symbol": "X", "hgvsp": "ENSP1:p.Leu10=", "hgvsc": "ENST1:c.30G>A"}],
    "colocated_variants": [{"id": "rs1", "frequencies": {"G": {"gnomadg": 0.42}}}],
}
KNOWLEDGE = {
    "study": "t",
    "samples": 100,
    "genes": {"TP53": {"frequency": 0.41, "hotspots": [["R175H", 30], ["R248Q", 20]]}},
}


def test_normalise_vep_keeps_the_fields_that_matter():
    n = normalise_vep(VEP_TP53)
    assert n["gene"] == "TP53" and n["consequence"] == "missense_variant"
    assert n["hgvsp"] == "p.Arg175His" and n["protein_change"] == "R175H" and n["residue"] == 175
    assert n["cosmic"] == ["COSV52661038"] and n["gnomad_af"] == 6.57e-06 and n["sift"] == "deleterious"
    c = normalise_vep(VEP_COMMON)
    assert c["protein_change"] == "L10=" and c["gnomad_af"] == 0.42


def test_grade_sets_aside_common_variants_and_ranks_hotspots_first():
    vs = [Variant("chr21", 0, "A", ("G",), gt=(0, 1)), Variant("chr17", 7675087, "C", ("T",), gt=(0, 1))]
    vep = {"chr21:0:A:G": normalise_vep(VEP_COMMON), "chr17:7675087:C:T": normalise_vep(VEP_TP53)}
    ranked = grade(vs, vep, KNOWLEDGE)
    assert ranked[0].gene == "TP53" and ranked[0].hotspot and not ranked[0].likely_germline
    assert ranked[0].pos == 7675088  # reported 1-based, as in the VCF
    assert ranked[0].score > 2.5
    assert ranked[1].likely_germline and ranked[1].gnomad_af >= GERMLINE_AF
    assert any("inherited" in e for e in ranked[1].evidence)
    mb = mutation_burden(ranked)
    assert mb["coding_somatic_variants"] == 1 and mb["evidence"] == "inferred"


def test_vep_line_is_one_based():
    from genomeos.cancer.tumour import _vep_line

    assert _vep_line(Variant("chr21", 34792137, "G", ("C",), gt=(0, 1))) == "21 34792138 . G C . . ."


def test_mutant_peptide_window():
    t = TumourVariant(
        "chr17",
        7675088,
        "C",
        "T",
        gene="TP53",
        consequence="missense_variant",
        protein_change="R175H",
        residue=175,
    )
    seq = "A" * 174 + "R" + "C" * 30
    p = mutant_peptides(t, seq, flank=3)
    assert p["wild_type"] == "AAARCCC" and p["mutant"] == "AAAHCCC" and p["window"] == [172, 178]
    assert "error" in mutant_peptides(t, "A" * 200)
    t.consequence = "stop_gained"
    assert mutant_peptides(t, seq) is None
