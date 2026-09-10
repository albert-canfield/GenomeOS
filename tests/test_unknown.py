"""Classifying the space between genes from sequence alone."""

import random

from genomeos.genome import Annotation
from genomeos.genome.unknown import (
    classify_block,
    investigate,
    kmer_table,
    load_patterns,
    tandem_fraction,
    unknown_blocks,
)

rng = random.Random(3)


def rand(n):
    # genome-like background: CpG-depleted and slightly AT-rich, so islands do not appear by chance
    return "".join(rng.choice("AACGTT") for _ in range(n)).replace("CG", "CA")


def test_pattern_file_and_named_classes():
    pats = load_patterns()
    names = {p.name for p in pats}
    assert {"telomeric_repeat", "cenp_b_box", "alu_core"} <= names
    tel = "TTAGGG" * 60 + rand(200)
    cls, ev, conf, f, hits = classify_block(tel, None, pats)
    assert cls == "telomere" and ev == "curated" and hits["telomeric_repeat"] >= 1


def test_tandem_low_complexity_orf_and_gap():
    pats = load_patterns()
    assert classify_block("N" * 500 + rand(50), None, pats)[0] == "gap"
    assert tandem_fraction("CACACACACACACACACACACA" * 10) > 0.9
    assert classify_block("CAG" * 400 + rand(300), None, pats)[0] == "tandem_repeat"
    codons = [
        a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")
    ]
    orf = "ATG" + "".join(rng.choice(codons) for _ in range(320)) + "TAA"
    cls = classify_block(rand(400) + orf + rand(400), None, pats)[0]
    assert cls == "long_orf"


def test_high_copy_kmers_reveal_an_interspersed_repeat():
    element = rand(300)
    genome = "".join(rand(rng.randint(200, 800)) + element for _ in range(80)) + rand(5000)
    counts = kmer_table(genome)
    pats = load_patterns()
    cls, ev, conf, f, _ = classify_block(element * 3 + rand(100), counts, pats)
    assert cls == "interspersed_repeat" and f["high_copy_fraction"] > 0.6
    cls2, _, _, f2, _ = classify_block(rand(2000), counts, pats)
    assert f2["high_copy_fraction"] < 0.05 and cls2 == "unclassified"


def test_investigate_orders_blocks_largest_first(tmp_path):
    ann = Annotation()
    seq = rand(6000)
    gff = tmp_path / "g.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c\tt\tgene\t1001\t1500\t.\t+\t.\tID=g1;gene_name=A;gene_type=protein_coding\n"
        "c\tt\tgene\t4001\t4200\t.\t-\t.\tID=g2;gene_name=B;gene_type=protein_coding\n"
    )
    ann = Annotation.from_gff3(gff, {"c"})
    blocks = unknown_blocks(ann, "c", len(seq), min_size=100)
    assert [b.length for b in blocks] == sorted((b.length for b in blocks), reverse=True)
    r = investigate(seq, ann, "c", min_size=100)
    assert r["unknown_blocks"] == 3 and r["unknown_bp"] == 1000 + 2500 + 1800
    assert set(r["by_class"]) >= {"unclassified"}
