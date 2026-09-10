from pathlib import Path

from genomeos.genome.regulatory import CCRE, ccre_index, count_in, load_ccres, summarise
from genomeos.genome.unknown import classify_block, load_patterns


def test_summary_and_density_queries():
    els = [
        CCRE("chr21", 100, 300, "a", "PLS", True),
        CCRE("chr21", 5000, 5300, "b", "dELS", False),
        CCRE("chr21", 9000, 9200, "c", "CTCF-only", True),
    ]
    s = summarise(els)
    assert s["elements"] == 3 and s["by_class"]["PLS"] == 1 and s["ctcf_bound"] == 2
    idx = ccre_index(els)
    assert count_in(idx, 0, 6000) == {"PLS": 1, "dELS": 1}
    assert count_in(idx, 8000, 20000) == {"CTCF-only": 1}


def test_chromatin_evidence_labels_a_block_regulatory():
    pats = load_patterns()
    seq = "".join("ACGTTGCA"[i % 8] for i in range(8000))
    cls, ev, conf, f, _ = classify_block(seq, None, pats, {"PLS": 1, "dELS": 2})
    assert cls == "regulatory" and ev == "curated" and f["ccre_per_10kb"] > 2
    assert classify_block(seq, None, pats, {"CTCF-only": 1})[0] != "regulatory"


def test_distilled_chr21_ccres_if_present():
    if not Path("data/results/ccres_chr21.bed.gz").exists():
        return
    els = load_ccres("chr21")
    assert len(els) > 10_000 and sum(1 for c in els if c.cls == "PLS") > 300
