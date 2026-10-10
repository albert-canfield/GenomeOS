# SPDX-License-Identifier: AGPL-3.0-or-later
"""The trace on the tissue's dominant isoform, beside the canonical one (genomeos/flow/tissue.py)."""

import json

from genomeos.coords import Locus
from genomeos.flow import trace_gene_in_tissue
from genomeos.flow.tissue import (
    ESTABLISHED,
    UNKNOWN,
    agreement,
    cached_isoforms,
    combine,
    gene_census,
    match_tissue,
    tissue_isoform,
)
from genomeos.flow.trace import canonical_transcript, trace_gene
from genomeos.genome import Chromosome, Genome, Sequence
from genomeos.ir.model import Gene, Transcript

#   0..9 exon1 (GGG UTR + ATG AAA) | 9..19 intron | 19..25 exon2 (GAT CCC) | 25..35 intron |
#   35..44 exon3 (TGG TAA + CCC UTR)
SEQ = "GGGATGAAA" + "GTTTTTTTAG" + "GATCCC" + "GTTTTTTTAG" + "TGGTAACCC"
E1, E2, E3 = Locus("c", 0, 9), Locus("c", 19, 25), Locus("c", 35, 44)


def _genome():
    g = Genome("t")
    g.chromosomes["c"] = Chromosome("c", Sequence(SEQ))
    return g


def _transcripts():
    canonical = Transcript(
        id="ENST1.2", kind="transcript", exons=[E1, E2, E3],
        cds_segments=[Locus("c", 3, 9), E2, Locus("c", 35, 41)], tags=["Ensembl_canonical"],
    )  # fmt: skip
    skips_exon2 = Transcript(
        id="ENST2.5", kind="transcript", exons=[E1, E3], cds_segments=[Locus("c", 3, 9), Locus("c", 35, 41)]
    )
    non_coding = Transcript(id="ENST4.1", kind="transcript", exons=[E1, E2])
    return [canonical, skips_exon2, non_coding]


def _gene():
    return Gene(id="g1", kind="gene", symbol="G1", locus=Locus("c", 0, 44))


ISO = {
    "gene": "G1",
    "isoforms": {
        "ENST1": {"tissues": {"Liver": 10.0, "Brain_Cerebellum": 2.0, "Lung": 5.0, "Skin": 0.5}},
        "ENST2": {"tissues": {"Liver": 1.0, "Brain_Cerebellum": 8.0, "Lung": 4.8, "Skin": 0.2}},
        "ENST3": {"tissues": {"Heart": 9.0}},  # not in the local models
        "ENST4": {"tissues": {"Muscle": 7.0}},  # non-coding in the local models
    },
}


def test_canonical_choice_is_unchanged_by_the_refactor():
    txs = _transcripts()
    assert canonical_transcript(txs).id == "ENST1.2"
    tr = trace_gene(_genome(), _gene(), txs)
    assert tr.protein == "MKDPW"


def test_the_call_is_established_only_when_the_rule_holds():
    assert tissue_isoform(ISO, "Liver")["status"] == ESTABLISHED
    brain = tissue_isoform(ISO, "Brain_Cerebellum")
    assert brain["status"] == ESTABLISHED and brain["transcript"] == "ENST2" and brain["share"] == 0.8
    tie = tissue_isoform(ISO, "Lung")
    assert tie["status"] == UNKNOWN and tie["reason"] == "tie" and tie["transcript"] == "ENST1"
    low = tissue_isoform(ISO, "Skin")
    assert low["status"] == UNKNOWN and low["reason"] == "low_expression"
    assert tissue_isoform(ISO, "Kidney")["reason"] == "tissue_not_measured"
    assert tissue_isoform(None, "Liver")["reason"] == "no_isoform_data"
    assert tissue_isoform({"isoforms": {}}, "Liver")["reason"] == "no_isoform_data"
    # the rule's parameters move the call, and the call records them
    assert tissue_isoform(ISO, "Lung", tie_ratio=0.99)["status"] == ESTABLISHED
    assert tissue_isoform(ISO, "Skin", min_tpm=0.1)["rule"] == {"min_tpm": 0.1, "tie_ratio": 0.9}


def test_an_exact_tie_is_a_tie_whatever_the_order():
    iso = {"isoforms": {"ENST9": {"tissues": {"T": 3.0}}, "ENST8": {"tissues": {"T": 3.0}}}}
    call = tissue_isoform(iso, "T")
    assert call["status"] == UNKNOWN and call["reason"] == "tie" and call["transcript"] == "ENST8"


def test_tissue_names_match_exactly_or_folded_never_partially():
    names = ["Brain_Cerebellum", "Brain_Cortex", "Cells_EBV-transformed_lymphocytes"]
    assert match_tissue("Brain_Cerebellum", names) == "Brain_Cerebellum"
    assert match_tissue("brain cerebellum", names) == "Brain_Cerebellum"
    assert match_tissue("cells ebv transformed lymphocytes", names) == "Cells_EBV-transformed_lymphocytes"
    assert match_tissue("Brain", names) is None
    assert match_tissue("Cerebellum", names) is None


def test_the_trace_reports_the_tissue_isoform_beside_the_canonical():
    g, gene, txs = _genome(), _gene(), _transcripts()
    brain = trace_gene_in_tissue(g, gene, txs, "brain cerebellum", ISO)
    assert brain.status == ESTABLISHED and brain.isoform_id == "ENST2.5" and not brain.same_as_canonical
    assert brain.isoform.protein == "MKW" and brain.canonical.protein == "MKDPW"
    # exon 2 is skipped in the cerebellum: missense on the canonical, intronic on the tissue's isoform
    s = brain.substitute(20, "A", "T")
    assert s["canonical"]["consequence"] == "missense" and s["canonical"]["hgvs_p"] == "p.Asp3Val"
    assert s["consequence"] == "intron" and s["agreement"] == "class_changed"
    # exon 3: a stop on both, at residue 5 of the canonical and residue 3 of the tissue's protein
    s = brain.substitute(36, "G", "A")
    assert s["consequence"] == "nonsense" and s["hgvs_p"] == "p.Trp3Ter"
    assert s["canonical"]["hgvs_p"] == "p.Trp5Ter" and s["agreement"] == "residue_changed"
    # exon 1 reads the same on both
    assert brain.substitute(6, "A", "G")["agreement"] == "same"
    liver = trace_gene_in_tissue(g, gene, txs, "Liver", ISO)
    assert liver.same_as_canonical and liver.substitute(20, "A", "T")["agreement"] == "same"
    json.dumps(brain.to_dict())


def test_an_unknown_call_never_falls_back_to_the_canonical():
    g, gene, txs = _genome(), _gene(), _transcripts()
    for tissue, reason in (
        ("Lung", "tie"),
        ("Skin", "low_expression"),
        ("Kidney", "tissue_not_measured"),
        ("Heart", "dominant_not_in_local_models"),
        ("Muscle", "dominant_non_coding"),
    ):
        t = trace_gene_in_tissue(g, gene, txs, tissue, ISO)
        assert t.status == UNKNOWN and t.isoform is None and t.call["reason"] == reason, tissue
        s = t.substitute(20, "A", "T")
        assert s["consequence"] == UNKNOWN and s["hgvs_p"] is None and s["residue"] is None
        assert s["reason"] == reason and s["agreement"] is None
        assert s["canonical"]["consequence"] == "missense"  # the canonical stays beside it, labelled
    none = trace_gene_in_tissue(g, gene, txs, "Liver", None)
    assert none.status == UNKNOWN and none.call["reason"] == "no_isoform_data"
    assert trace_gene_in_tissue(g, gene, [], "Liver", ISO) is None


def test_agreement_classes():
    a = {"consequence": "missense", "residue": 5, "hgvs_p": "p.Ala5Thr", "hgvs_c": "c.13G>A"}
    assert agreement(a, dict(a)) == "same"
    assert agreement(a, {**a, "consequence": "synonymous"}) == "class_changed"
    assert agreement(a, {**a, "residue": 7}) == "residue_changed"
    assert agreement(a, {**a, "hgvs_p": "p.Gly5Arg"}) == "change_differs"


def test_gene_census_counts_per_tissue_with_denominators():
    g, gene, txs = _genome(), _gene(), _transcripts()
    variants = [(20, "A", "T"), (36, "G", "A"), (6, "A", "G"), (12, "T", "C"), (1, "G", "A")]
    tissues = ["Liver", "Brain_Cerebellum", "Lung", "Skin", "Heart", "Muscle", "Kidney"]
    r = gene_census(g, gene, txs, ISO, variants, tissues)
    # the intronic (12) and 5'UTR (1) changes are coding on no traced transcript: outside the universe
    assert r["variants"] == 3 and r["coding_on_canonical"] == 3 and r["transcripts_traced"] == 2
    brain = r["per_tissue"]["Brain_Cerebellum"]
    assert brain["variants"] == 3 and brain["established"] == 3 and brain["other_transcript"] == 3
    assert (brain["same"], brain["class_changed"], brain["residue_changed"], brain["disagree"]) == (
        1,
        1,
        1,
        2,
    )
    liver = r["per_tissue"]["Liver"]
    assert liver["same_transcript"] == 3 and liver["same"] == 3 and liver["disagree"] == 0
    for t, reason in (("Lung", "tie"), ("Skin", "low_expression"), ("Heart", "dominant_not_in_local_models"),
                      ("Muscle", "dominant_non_coding"), ("Kidney", "tissue_not_measured")):  # fmt: skip
        x = r["per_tissue"][t]
        assert x["established"] == 0 and x["unknown"][reason] == 3 and x["disagree"] == 0, t
    for x in r["per_tissue"].values():  # every universe variant is counted once per tissue
        assert x["established"] + sum(x["unknown"].values()) == x["variants"]
    json.dumps(r)
    both = combine([r, r], tissues)
    assert both["Brain_Cerebellum"]["disagree"] == 4 and both["Brain_Cerebellum"]["variants"] == 6
    assert both["Brain_Cerebellum"]["genes_disagreeing"] == 2
    assert both["Brain_Cerebellum"]["genes_established_other_transcript"] == 2
    assert both["Liver"]["genes_established_other_transcript"] == 0
    assert both["Lung"]["unknown"]["tie"] == 6


def test_a_variant_coding_only_on_the_tissue_isoform_enters_the_universe():
    # the tissue's isoform reads into what is the canonical's 3'UTR
    g, gene = _genome(), _gene()
    canonical = _transcripts()[0]
    longer = Transcript(
        id="ENST2.5", kind="transcript", exons=[E1, E3],
        cds_segments=[Locus("c", 3, 9), Locus("c", 35, 44)],
    )  # fmt: skip
    iso = {"isoforms": {"ENST1": {"tissues": {"T": 1.0}}, "ENST2": {"tissues": {"T": 9.0}}}}
    r = gene_census(g, gene, [canonical, longer], iso, [(42, "C", "A")], ["T"])
    assert r["variants"] == 1 and r["coding_on_canonical"] == 0
    assert r["per_tissue"]["T"]["class_changed"] == 1


def test_cached_isoforms_reads_the_cache_and_never_fetches(tmp_path):
    assert cached_isoforms("G1", tmp_path) is None
    (tmp_path / "gtex_isoforms_G1.json").write_text(json.dumps(ISO))
    assert cached_isoforms("g1", tmp_path)["gene"] == "G1"
    (tmp_path / "gtex_isoforms_BAD.json").write_text("{not json")
    assert cached_isoforms("BAD", tmp_path) is None
