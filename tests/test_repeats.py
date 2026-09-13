from genomeos.genome.repeats import Repeat, RepeatIndex, summarise
from genomeos.genome.unknown import classify_block, load_patterns


def test_repeat_index_coverage_by_class():
    reps = [
        Repeat(100, 400, "LINE", "L1", "L1PA2", 0.05),
        Repeat(350, 500, "SINE", "Alu", "AluY", 0.02),
        Repeat(900, 950, "Simple_repeat", "(CA)n", "(CA)n", 0.0),
    ]
    idx = RepeatIndex(reps)
    cov = idx.coverage(0, 1000)
    assert cov == {"LINE": 300, "SINE": 150, "Simple_repeat": 50}
    assert idx.coverage(380, 420) == {"LINE": 20, "SINE": 40}
    assert idx.coverage(600, 800) == {}
    s = summarise("c", reps, 1000)
    assert s["by_class"]["LINE"]["bp"] == 300 and s["repeat_fraction"] == 0.5
    assert s["top_families"]["LINE/L1"] == 300


def test_curated_repeat_coverage_beats_the_orf_rule():
    import random

    rng = random.Random(3)
    # a stop-free frame of 1,200 codons, as a LINE-1 ORF2 would give
    codons = [
        a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in ("TAA", "TAG", "TGA")
    ]
    seq = "ATG" + "".join(rng.choice(codons) for _ in range(1200)) + "TAA"
    pats = load_patterns()
    cls, ev, conf, f, _ = classify_block(seq, None, pats)
    assert cls == "long_orf"
    cls2, ev2, conf2, f2, _ = classify_block(seq, None, pats, None, {"LINE": int(len(seq) * 0.8)})
    assert cls2 == "interspersed_repeat_LINE" and ev2 == "curated" and conf2 == 0.9
    assert f2["interspersed_coverage"] == 0.8
    cls3, *_ = classify_block(seq, None, pats, None, {"Satellite": int(len(seq) * 0.9)})
    assert cls3 == "satellite_array"
    cls4, *_ = classify_block(seq, None, pats, None, {"Simple_repeat": int(len(seq) * 0.6)})
    assert cls4 == "tandem_repeat"
