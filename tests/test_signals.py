"""The genome's delimiters learned from data, and found again in raw sequence."""

from pathlib import Path

import pytest

from genomeos.genome import Annotation, IndexedGenome, Pwm, default_gencode, learn_signals, scan
from genomeos.genome.signals import SignalSet

CHR21 = Path("data/reference/chr21.fa.gz")
GFF = default_gencode({"chr21"}) or Path("missing")


def test_pwm_learns_consensus_and_scores():
    p = Pwm.learn("t", ["CAGGTAAGT", "AAGGTGAGT", "CAGGTAAGA", "GAGGTAAGT"], offset=3)
    assert p.consensus[3:5] == "GT" and p.length == 9
    assert p.score("CAGGTAAGT") > p.score("CCCCCCCCC")
    assert p.score("CAGGTAAGT") <= p.max_score() + 1e-9


@pytest.mark.skipif(not (CHR21.exists() and GFF.exists()), reason="fetch chr21 first")
def test_signals_learned_from_chr21_are_the_known_biology(tmp_path):
    ann = Annotation.from_gff3(GFF, {"chr21"})
    g = IndexedGenome(CHR21)
    sig = learn_signals(ann, g, "chr21")
    st = sig.stats
    assert st["donor_GT_fraction"] > 0.98 and st["acceptor_AG_fraction"] > 0.99
    assert sig.pwms["splice_donor"].consensus[3:5] == "GT"
    assert sig.pwms["splice_acceptor"].consensus[12:14] == "AG"
    assert sig.pwms["start_kozak"].consensus[6:9] == "ATG"
    stops = st["stop_codons"]
    canonical = sum(v for k, v in stops.items() if k in ("TAA", "TAG", "TGA"))
    assert canonical / sum(stops.values()) > 0.99  # one chr21 model ends off a stop codon
    assert 0.5 < st["polyA_signal_in_last_40nt"] < 0.9
    p = tmp_path / "sig.json"
    sig.save(p)
    assert SignalSet.load(p).pwms["splice_donor"].examples == sig.pwms["splice_donor"].examples
    # the learned donor PWM finds the real donors of APP in a raw window at the
    # measured 90%-recall threshold; APP is on the minus strand, so each exon
    # except the lowest one ends (in transcript order) at its forward start
    from genomeos.coords import Locus, Strand

    app = ann.gene("APP")
    m = ann.to_module("x")
    tx = next(t for t in m.entities[app.id].transcripts if "Ensembl_canonical" in t.tags)
    exons = sorted(tx.exons, key=lambda l: l.start)
    s, e = exons[0].start - 100, exons[3].end + 100
    window = g.fetch(Locus("chr21", s, e))
    sep = st["separability"]["splice_donor"]
    hits = [
        h
        for h in scan(window, sig, min_relative=sep["true_p10"])
        if h.signal == "splice_donor" and h.strand is Strand.MINUS
    ]
    true_donors = {ex.start - 1 - s for ex in exons[1:4]}
    found = sum(1 for d in true_donors if any(abs(h.pos - d) <= 1 for h in hits))
    assert found == 3, (found, sorted(true_donors), [h.pos for h in hits][:20])
    # the honest picture: a local signal alone is real but not sufficient
    assert 0.2 < sep["true_p10"] < 0.5 and 2 < sep["background_hits_per_kb_at_p10"] < 20
    g.close()
