# SPDX-License-Identifier: AGPL-3.0-or-later
"""The syntax candidates read one by one: open frames, copies of known proteins, the reading rules."""

import json
from types import SimpleNamespace

import pytest

from genomeos.attribution.candidates import (
    Intervals,
    binomial_tail,
    candidates,
    flank_genes,
    frame_score,
    merge,
    paralogues,
    reading,
    replicate_p,
    revcomp,
    set_tests,
    translate,
)

LO = {c: 1.0 for c in ("ATG", "GCC", "AAG", "GAG", "CTG")}  # a toy codon table: these read as coding


def test_translate_and_open_frames():
    orf = "ATG" + "GCC" * 20 + "AAG" * 5
    assert translate(orf).startswith("MA") and "*" not in translate(orf)
    fs = frame_score(orf, LO)
    assert fs["open_frame"] and fs["score"] == pytest.approx(1.0) and fs["peptide"].startswith("MAAA")
    # the same frame read from the other strand is found too
    assert frame_score(revcomp(orf), LO)["score"] == pytest.approx(1.0)
    # stops in every frame of both strands: no open frame, so no coding exon
    fs = frame_score("CTAG" * 30, LO)  # a stop every twelve bases in all six frames
    assert not fs["open_frame"] and fs["score"] is None


def test_intervals_and_merge():
    assert merge([(0, 10), (15, 20), (40, 50)], gap=5) == [(0, 20), (40, 50)]
    iv = Intervals([(0, 10, "a"), (5, 30, "b"), (100, 110, "c")])
    assert [v for _, _, v in iv.overlapping(8, 12)] == ["a", "b"]
    assert iv.covered(0, 40) == 30 and not iv.any(40, 100)


def test_binomial_tail():
    assert binomial_tail(0, 10, 0.1) == 1.0
    assert binomial_tail(1, 1, 0.05) == pytest.approx(0.05)
    assert binomial_tail(2, 3, 0.05) == pytest.approx(3 * 0.05**2 * 0.95 + 0.05**3)


def test_replicate_p():
    assert replicate_p(5, [1, 2, 3]) == pytest.approx(0.25)
    assert replicate_p(1, [1, 2, 3]) == pytest.approx(1.0)
    assert replicate_p(1, []) is None


def test_paralogues_finds_the_parent_protein(tmp_path):
    parent = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
    (tmp_path / "P1.json").write_text(
        json.dumps(
            {"gene": "P1", "sections": {"identity": {"items": {"accession": "Q1", "sequence": parent}}}}
        )
    )
    (tmp_path / "P2.json").write_text(
        json.dumps({"gene": "P2", "sections": {"identity": {"items": {"sequence": "W" * 40}}}})
    )
    out = paralogues({"copy": parent[5:25], "novel": "HHHHCCCCHHHHCCCC"}, proteins=tmp_path)
    assert out["copy"]["protein"] == "P1" and out["copy"]["share"] == pytest.approx(1.0)
    assert out["novel"]["protein"] is None and out["novel"]["share"] == 0.0


def test_candidates_filters_case_and_copies_and_refuses_a_short_list(tmp_path):
    rows = [
        {"start": 0, "end": 1000, "case": "syntax", "copy": False},
        {"start": 2000, "end": 3000, "case": "syntax", "copy": True},
        {"start": 4000, "end": 5000, "case": "relaxed", "copy": False},
    ]
    doc = {"candidates": rows, "real_unknown": {"by_case": {"syntax": {"blocks": 1}}}}
    (tmp_path / "organised_chr2.json").write_text(json.dumps(doc))
    got = candidates(tmp_path)
    assert [(r["chrom"], r["start"]) for r in got] == [("chr2", 0)]
    doc["real_unknown"]["by_case"]["syntax"]["blocks"] = 3
    (tmp_path / "organised_chr2.json").write_text(json.dumps(doc))
    with pytest.raises(ValueError):
        candidates(tmp_path)


def _gene(symbol, start, end, strand, kind="protein_coding"):
    return SimpleNamespace(
        symbol=symbol,
        type=kind,
        locus=SimpleNamespace(start=start, end=end, strand=SimpleNamespace(value=strand)),
    )


def test_flank_genes_says_which_end_faces_the_block():
    genes = [_gene("A", 0, 1000, "+"), _gene("B", 5000, 9000, "+"), _gene("C", 9500, 12000, "-")]
    ch = SimpleNamespace(genes=genes, gene_ends=sorted((g.locus.end, i) for i, g in enumerate(genes)))
    fl = flank_genes(ch, 1000, 5000)
    assert fl["left"]["gene"] == "A" and fl["left"]["facing"] == "3' end" and fl["left"]["distance"] == 0
    assert fl["right"]["gene"] == "B" and fl["right"]["facing"] == "5' end (TSS)"
    fl = flank_genes(ch, 9000, 9500)
    assert fl["right"]["gene"] == "C" and fl["right"]["facing"] == "3' end"


def _ev(**over):
    ev = {
        "chrom": "chr2",
        "start": 10_000,
        "end": 12_000,
        "flanks": {"left": None, "right": None},
        "segment_spans": [[10_500, 10_700]],
        "segments": 1,
        "ccre_classes": [],
        "ccre_classes_on_conserved": [],
        "reader": {"cells": 11, "open_on_element": []},
        "controls": {"mean_cells_open": 1.0, "coding_like_rate": 0.05},
        "ground_truth": {"vista": [], "mpra": [], "clinvar_pathogenic": [], "eqtl_genes": []},
        "repeats": {
            "conserved_bp": 200,
            "conserved_bp_by_repeat_age": {"unique": 200},
            "conserved_families": [],
        },
        "coding": {"segments": [{"open_frame": False, "score": None}]},
        "parser": {"exons_on_conserved": 0},
        "human_axis_footing": {"kilobases": 2, "kilobases_with_exon_bases": 0, "exon_bp_in_kilobases": 0},
        "paralogue": None,
    }
    ev.update(over)
    return ev


def test_reading_rules():
    r = reading(_ev())
    assert r["class"] == "unexplained" and r["confidence"] == 0.0
    assert any("open frame" in x for x in r["ruled_out"]) and any("DNase" in x for x in r["ruled_out"])

    primate = {
        "conserved_bp": 200,
        "conserved_bp_by_repeat_age": {"primate_specific": 150},
        "conserved_families": [],
    }
    assert reading(_ev(repeats=primate))["class"] == "alignment_artefact"
    rna = {
        "conserved_bp": 200,
        "conserved_bp_by_repeat_age": {"structured_rna": 180},
        "conserved_families": ["7SL"],
    }
    assert reading(_ev(repeats=rna))["group"] == "rna"

    open_cells = {"cells": 11, "open_on_element": ["a", "b", "c", "d"]}
    r = reading(_ev(reader=open_cells, ccre_classes=["dELS"], ccre_classes_on_conserved=["dELS"]))
    assert r["class"] == "regulatory" and r["confidence"] == pytest.approx(0.45)
    # open, but no more than the control windows are: not support
    r = reading(_ev(reader=open_cells, controls={"mean_cells_open": 3.5, "coding_like_rate": 0.05}))
    assert r["class"] == "unexplained"
    # a registry element off the conserved bases does not count
    r = reading(_ev(ccre_classes=["PLS"]))
    assert r["class"] == "unexplained" and any("off its conserved bases" in n for n in r["notes"])
    assert reading(_ev(ccre_classes=["PLS"], ccre_classes_on_conserved=["PLS"]))["class"] == "promoter"

    like = [{"open_frame": True, "score": 0.5, "coding_like": True}] * 2 + [
        {"open_frame": False, "score": None}
    ]
    r = reading(_ev(coding={"segments": like}, parser={"exons_on_conserved": 1}))
    assert r["class"] == "coding_exon" and r["group"] == "coding" and r["confidence"] == pytest.approx(0.45)
    # one coding-like segment among many is what the control rate gives
    many = like[:1] + [{"open_frame": False, "score": None}] * 20
    assert reading(_ev(coding={"segments": many}))["class"] == "unexplained"
    # coding-like segments inside transposon remnants are not counted
    old = [{**x, "in_repeat": "L1:L1ME2"} for x in like[:2]] + like[2:]
    r = reading(_ev(coding={"segments": old}))
    assert r["class"] == "unexplained" and any("transposon" in n for n in r["notes"])
    para = {
        "query": "conserved segment",
        "protein": "RPL7",
        "shared_kmers": 20,
        "query_kmers": 25,
        "share": 0.8,
    }
    assert reading(_ev(coding={"segments": like}, paralogue=para))["class"] == "coding_copy"
    short = {**para, "query_kmers": 10, "shared_kmers": 9}
    assert reading(_ev(coding={"segments": like}, paralogue=short))["class"] == "coding_exon"
    joined = {"exons_on_conserved": 1, "best_on_conserved": {"joins_annotated_exons_of": ["GNAT1"]}}
    r = reading(_ev(parser=joined))
    assert r["class"] == "coding_exon" and "GNAT1" in r["basis"][0]

    left = {"gene": "G", "type": "protein_coding", "strand": "+", "distance": 0, "facing": "3' end"}
    r = reading(_ev(flanks={"left": left, "right": None}))
    assert r["class"] == "transcript_extension"
    lnc = {**left, "type": "lncRNA"}
    assert reading(_ev(flanks={"left": lnc, "right": None}))["class"] == "unexplained"

    borrowed = {"kilobases": 2, "kilobases_with_exon_bases": 2, "exon_bp_in_kilobases": 300}
    base = reading(_ev(reader=open_cells, ccre_classes_on_conserved=["dELS"]))
    lower = reading(_ev(reader=open_cells, ccre_classes_on_conserved=["dELS"], human_axis_footing=borrowed))
    assert lower["human_axis_borrowed"] and lower["confidence"] == pytest.approx(base["confidence"] - 0.1)


def test_set_tests_reads_candidates_against_replicates():
    def row(cells, ccre, ctrl_cells):
        bp = {
            "window": 1000,
            "in_ccre": 100,
            "conserved": 100,
            "conserved_in_ccre": 10,
            "open": 50,
            "conserved_open": 5,
        }
        return {
            "reader": {"cells": 11, "open_on_element": ["x"] * cells},
            "ccres": ["E"] if ccre else [],
            "parser": {"exons_on_conserved": 0},
            "coding": {"segments_scored": 0, "segments": []},
            "bp": bp,
            "controls": {
                "segments_scored": 0,
                "segments_open_frame": 0,
                "segments_coding_like": 0,
                "bp": bp,
                "_per_window": [
                    {"cells_open": c, "ccre": False, "parser_exon_on_conserved": False} for c in ctrl_cells
                ],
            },
        }

    t = set_tests([row(6, True, [0, 1, 0]), row(4, True, [1, 0, 0])])
    assert t["replicates"] == 3
    assert t["reader_cells_open_mean"]["candidates"] == 5.0
    assert t["reader_cells_open_mean"]["p"] == pytest.approx(0.25)
    assert t["with_registry_element"]["controls"] == 0.0
