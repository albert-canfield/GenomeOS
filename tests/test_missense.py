# SPDX-License-Identifier: AGPL-3.0-or-later
"""AlphaMissense rows streamed for a person's variants, cached under the person, the inventory re-ranked."""

import json

from genomeos.genome import missense as am

TABLE = [
    "# Copyright\n",
    "#CHROM\tPOS\tREF\tALT\tgenome\tuniprot_id\ttranscript_id\tprotein_variant\tam_pathogenicity\tam_class\n",
    "chr21\t100\tA\tG\thg38\tP1\tENST1.4\tV2L\t0.9\tlikely_pathogenic\n",
    "chr21\t100\tA\tG\thg38\tP1\tENST2.1\tV9L\t0.2\tlikely_benign\n",
    "chr21\t100\tA\tT\thg38\tP1\tENST1.4\tV2F\t0.5\tambiguous\n",
    "chr22\t500\tC\tT\thg38\tP2\tENST3.2\tR5W\t0.7\tlikely_pathogenic\n",
]


def test_stream_keeps_only_the_wanted_variants():
    found, cost = am.stream_scores({("chr21", 100, "A", "G"), ("chr22", 999, "C", "T")}, rows=iter(TABLE))
    assert cost["rows_read"] == 4 and cost["variants_found"] == 1
    assert [e["transcript"] for e in found["chr21:100:A:G"]] == ["ENST1.4", "ENST2.1"]


def test_scores_for_caches_under_the_person(tmp_path):
    variants = [
        {"chrom": "chr21", "pos": 100, "ref": "A", "alt": "G"},
        {"chrom": "chr22", "pos": 500, "ref": "C", "alt": "T"},
    ]
    scores = am.scores_for("demo", variants, tmp_path, rows=iter(TABLE))
    assert set(scores) == {"chr21:100:A:G", "chr22:500:C:T"}
    cached = json.loads((tmp_path / "demo" / "alphamissense.json").read_text())
    assert cached["asked"] == ["chr21:100:A:G", "chr22:500:C:T"] and cached["last_pass"]["rows_read"] == 4
    # a second call with the same variants never streams (rows=None would fail offline)
    again = am.scores_for("demo", variants, tmp_path, rows=None)
    assert again == scores


def test_pick_prefers_the_persons_transcript_then_the_highest():
    entries = [
        {"transcript": "ENST1.4", "score": 0.9, "class": "likely_pathogenic", "protein_variant": "V2L"},
        {"transcript": "ENST2.1", "score": 0.2, "class": "likely_benign", "protein_variant": "V9L"},
    ]
    assert am.pick(entries, "ENST2")["on_transcript"] == "the person's"
    assert am.pick(entries, "ENST2")["score"] == 0.2
    assert am.pick(entries, None)["score"] == 0.9 and am.pick([], "ENST1") is None


def test_annotate_reranks_by_effect_and_counts_classes():
    missense = [
        {
            "gene": "A",
            "chrom": "chr21",
            "pos": 100,
            "ref": "A",
            "alt": "G",
            "transcript": "ENST1.4",
            "zygosity": "heterozygous",
            "site": "none",
        },
        {
            "gene": "B",
            "chrom": "chr22",
            "pos": 500,
            "ref": "C",
            "alt": "T",
            "transcript": None,
            "zygosity": "homozygous",
            "site": "site",
        },
        {
            "gene": "C",
            "chrom": "chr1",
            "pos": 7,
            "ref": "G",
            "alt": "A",
            "transcript": None,
            "zygosity": "homozygous",
            "site": "domain",
        },
    ]
    scores, _ = am.stream_scores({("chr21", 100, "A", "G"), ("chr22", 500, "C", "T")}, rows=iter(TABLE))
    out = am.annotate(missense, scores)
    assert [m["gene"] for m in out["missense_by_effect"]] == ["A", "B", "C"]
    assert out["missense_predicted"] == {
        "likely_pathogenic": 2,
        "ambiguous": 0,
        "likely_benign": 0,
        "unscored": 1,
    }
    assert out["missense_likely_pathogenic_homozygous"] == 1
    assert missense[0]["predicted"]["confidence"] == 0.7 and missense[2]["predicted"] is None
