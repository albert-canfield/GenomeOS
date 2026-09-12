"""Motifs over promoters: indexed scanning of JASPAR profiles, enrichment and `requires:` (area J step 4)."""

import random

from genomeos.genome.motifs import (
    Index,
    dinucleotide_shuffle,
    enrichment,
    library_factors,
    parse_jaspar,
    requires,
    scan,
    shuffled,
)

GATA = """>MA0001.1\tGATA_TEST
A  [ 0 0 0 20 0 0 ]
C  [ 0 0 0 0 0 0 ]
G  [ 20 0 0 0 0 20 ]
T  [ 0 0 20 0 20 0 ]
>MA0002.1\tLONGONE
A  [ 20 0 0 0 0 0 0 0 0 20 ]
C  [ 0 20 0 0 0 0 0 0 20 0 ]
G  [ 0 0 20 0 0 0 0 20 0 0 ]
T  [ 0 0 0 20 20 20 20 0 0 0 ]
"""


def _code(kmer: str) -> int:
    code = 0
    for b in kmer:
        code = (code << 2) | "ACGT".index(b)
    return code


def test_parse_prepare_and_feasible_cores():
    motifs = parse_jaspar(GATA)
    assert [m.name for m in motifs] == ["GATA_TEST", "LONGONE"]
    assert motifs[0].width == 6 and motifs[1].width == 10
    m = motifs[0].prepare()
    assert m.core_k == 6 and m.maximum > m.threshold > m.minimum
    cores = m.feasible_cores()
    assert _code("GATATG") in cores and len(cores) < 4096  # the consensus passes, most k-mers do not
    big = motifs[1].prepare()
    assert big.core_k == 8 and 0 <= big.core_offset <= 2


def test_scan_finds_planted_sites_on_both_strands_and_enrichment_reads_them():
    rng = random.Random(3)
    motifs = [m.prepare() for m in parse_jaspar(GATA)]
    background = ["".join(rng.choice("ACGT") for _ in range(300)) for _ in range(6)]
    planted = list(background)
    planted[0] = planted[0][:100] + "GATATG" + planted[0][106:]  # a 6-mer: found, but chance hits tie it
    planted[1] = planted[1][:50] + "TGCAAAACGT" + planted[1][56:]  # reverse complement of LONGONE
    planted[2] = planted[2][:10] + "ACGTTTTGCA" + planted[2][20:]  # LONGONE forward
    hits = scan(motifs, planted)
    g = hits["GATA_TEST"]
    assert 0 in g and g[0][0] >= 0.85
    long = hits["LONGONE"]
    assert long[2] == (1.0, 10, 0) and long[1] == (1.0, 50, 1)
    assert set(long) <= {1, 2}  # a 10-mer does not appear by chance in 1,800 random bases
    assert scan(motifs, ["GATNTG" * 10]) == {}  # an N breaks the window and never scores
    ctrl = scan(motifs, shuffled(planted))
    e = enrichment(hits, ctrl, len(planted))
    assert e["LONGONE"]["sequences"] == 2 and e["LONGONE"]["enrichment"] > 1
    req = requires(
        {"GATA_TEST": (0.9, 100, 0), "LONGONE": (0.95, 10, 1)},
        {"GATA_TEST": {"enrichment": 1.0}, "LONGONE": {"enrichment": 3.0}},
    )
    assert [r["factor"] for r in req] == ["LONGONE"] and req[0]["strand"] == "-"


def test_index_packs_positions_and_library_factors_rank_enriched_ones():
    idx = Index(["ACGTACGTAC"], {4})
    entries = sorted(idx.tables[4][_code("ACGT")])
    # forward positions 0 and 4; the reverse complement GTACGTACGT holds ACGT at 2 and 6
    assert [(e >> 11) & 1 for e in entries] == [0, 0, 1, 1] and [e & 2047 for e in entries] == [0, 4, 2, 6]
    gene_hits = {f"G{i}": {"TFA": 0.9} for i in range(10)}
    for i in range(5):
        gene_hits[f"G{i}"]["TFB"] = 0.9
    libs = library_factors(
        gene_hits, {"lib.x": [f"G{i}" for i in range(5)], "lib.small": ["G0"]}, min_members=5
    )
    assert "lib.small" not in libs
    rows = libs["lib.x"]["factors"]
    assert [r["factor"] for r in rows] == ["TFB"] and rows[0]["ratio"] == 2.0  # TFA is everywhere: no ratio


def test_dinucleotide_shuffle_keeps_every_dinucleotide_count():
    from collections import Counter

    rng = random.Random(5)
    seq = "".join(rng.choice("ACGT") for _ in range(400)) + "CGCGCGCG" + "AAAATTTT"
    for _ in range(5):
        out = dinucleotide_shuffle(seq, rng)
        assert len(out) == len(seq) and out[0] == seq[0] and out[-1] == seq[-1]
        pairs = lambda s: Counter(zip(s, s[1:], strict=False))  # noqa: E731
        assert pairs(out) == pairs(seq)
    assert dinucleotide_shuffle("AC", rng) == "AC"


def test_distil_aggregates_requires_across_chromosomes_and_finds_operator_pairs(tmp_path):
    import json

    from genomeos.genome.motifs import distil

    def result(chrom, genes):
        return {
            "result": f"motifs_{chrom}",
            "chrom": chrom,
            "promoters": len(genes),
            "genes": {
                g: {"factors_hit": len(f), "requires": [{"factor": x, "score": 0.9} for x in f]}
                for g, f in genes.items()
            },
            "factor_enrichment": {},
            "cost": {"seconds": 1},
        }

    # library members carry A and B together; the rest of the genome carries A or B alone or neither
    lib = {f"L{i}": ["TFA", "TFB"] for i in range(8)}
    rest = {f"R{i}": (["TFA"] if i % 3 == 0 else ["TFB"] if i % 3 == 1 else ["TFC"]) for i in range(30)}
    (tmp_path / "motifs_chrA.json").write_text(
        json.dumps(result("chrA", {**lib, **dict(list(rest.items())[:15])}))
    )
    (tmp_path / "motifs_chrB.json").write_text(json.dumps(result("chrB", dict(list(rest.items())[15:]))))
    s = distil(tmp_path, members={"lib.ab": list(lib), "lib.tiny": ["L0"]})
    assert s["chromosomes"] == 2
    assert s["promoters"] == 38 and s["factor_share"]["TFA"] == round(18 / 38, 4)
    ab = s["libraries"]["lib.ab"]
    assert ab["members_placed"] == 8 and {f["factor"] for f in ab["factors"]} == {"TFA", "TFB"}
    pair = ab["operators"][0]
    assert set(pair["factors"]) == {"TFA", "TFB"} and pair["members_with_both"] == 8 and pair["ratio"] > 2
    assert "lib.tiny" not in s["libraries"]
