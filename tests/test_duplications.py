"""Segmental duplications as curated copy-and-paste over the UNKNOWN blocks (area J step 3)."""

from genomeos.genome.duplications import (
    SegDup,
    band,
    covered,
    merge,
    pair_supported,
    parse_similar,
    summarise,
)


def test_helpers():
    assert merge([(5, 10), (0, 6), (20, 30)]) == [(0, 10), (20, 30)]
    assert covered([(0, 10), (20, 30)], 5, 25) == 10 and covered([(0, 10)], 50, 60) == 0
    assert band(0.99) == "young" and band(0.96) == "middle" and band(0.9) == "old"
    assert parse_similar("chr21:100-200 (jaccard 0.31)") == ("chr21", 100, 200)
    assert parse_similar("nonsense") is None


def test_summarise_blocks_pairs_and_elements():
    dups = [
        SegDup(1000, 3000, "chr21", 9000, 11000, 0.99, "+"),  # a young intra-chromosomal copy: blocks A and C
        SegDup(2000, 2500, "chr9", 100, 600, 0.93, "-"),  # an old copy from another chromosome, overlapping A
    ]
    blocks = [
        {"start": 0, "end": 4000, "length": 4000, "class": "unique_intergenic", "guess": {"tier": "neutral"}},
        {
            "start": 5000,
            "end": 6000,
            "length": 1000,
            "class": "unique_intergenic",
            "guess": {"tier": "neutral"},
        },
        {
            "start": 8000,
            "end": 12000,
            "length": 4000,
            "class": "mixed_intergenic",
            "guess": {"tier": "fossil"},
        },
    ]
    unknown = [
        {
            "start": 0,
            "end": 4000,
            "similar_to": ["chr21:8000-12000 (jaccard 0.4)", "chr21:5000-6000 (jaccard 0.2)"],
        },
        {"start": 8000, "end": 12000, "similar_to": ["chr21:0-4000 (jaccard 0.4)"]},
    ]
    elements = [{"id": "E1", "start": 1500, "end": 1700}, {"id": "E2", "start": 5100, "end": 5300}]
    s = summarise("chr21", dups, 20_000, blocks, unknown, elements)
    assert (
        s["pairs"] == 2 and s["intra_chromosomal"] == 1 and s["bands"] == {"young": 1, "middle": 0, "old": 1}
    )
    assert (
        s["duplicated_bp"] == 2000 and s["duplicated_fraction"] == 0.1
    )  # 1000-3000 merged; the other span elsewhere
    a, b, c = s["blocks"]
    assert a["duplicated_fraction"] == 0.5 and a["pairs"] == 2 and b["duplicated_fraction"] == 0.0
    assert c["duplicated_fraction"] == 0.0  # the partner span is not a track row of its own here
    assert s["by_tier"]["neutral"] == {
        "blocks": 2,
        "bp": 5000,
        "duplicated_bp": 2000,
        "blocks_mostly_duplicated": 1,
        "duplicated_fraction": 0.4,
    }
    # of the two heuristic pairs (A-C and A-B), the curated track supports A-C only
    assert s["similar_to_pairs"] == {"pairs": 2, "supported_by_curated_duplication": 1, "share": 0.5}
    assert pair_supported(dups, "chr21", (8000, 12000), (0, 4000)) and not pair_supported(
        dups, "chr21", (5000, 6000), (0, 4000)
    )
    assert s["elements"]["inside_a_duplication"] == 1 and s["elements"]["ids"] == ["E1"]


def test_distil_sums_chromosomes(tmp_path):
    import json

    from genomeos.genome.duplications import distil

    def result(chrom, dup_bp, frac, cu_frac):
        return {
            "result": f"duplication_{chrom}",
            "chrom": chrom,
            "pairs": 10,
            "intra_chromosomal": 4,
            "bands": {"young": 1, "middle": 2, "old": 7},
            "duplicated_bp": dup_bp,
            "duplicated_fraction": frac,
            "by_tier": {
                "constrained_unknown": {
                    "blocks": 2,
                    "bp": 1000,
                    "duplicated_bp": round(cu_frac * 1000),
                    "blocks_mostly_duplicated": 1,
                    "duplicated_fraction": cu_frac,
                }
            },
            "similar_to_pairs": {"pairs": 1, "supported_by_curated_duplication": 1, "share": 1.0},
            "elements": {"scored": 100, "inside_a_duplication": 2, "share": 0.02},
        }

    (tmp_path / "duplication_chrA.json").write_text(json.dumps(result("chrA", 1000, 0.1, 0.6)))
    (tmp_path / "duplication_chrB.json").write_text(json.dumps(result("chrB", 3000, 0.1, 0.2)))
    s = distil(tmp_path)
    assert s["chromosomes"] == 2 and s["pairs"] == 20 and s["bands"]["old"] == 14
    assert s["duplicated_bp"] == 4000 and s["duplicated_fraction_of_genome"] == 0.1
    cu = s["by_tier"]["constrained_unknown"]
    assert (
        cu["bp"] == 2000 and cu["duplicated_fraction"] == 0.4 and cu["share_blocks_mostly_duplicated"] == 0.5
    )
    assert s["similar_to_pairs"] == {"pairs": 2, "supported": 2, "share": 1.0}
    assert s["elements"] == {"scored": 200, "inside": 4, "share": 0.02}
    assert s["per_chromosome"]["chrA"]["constrained_unknown_duplicated"] == 0.6
