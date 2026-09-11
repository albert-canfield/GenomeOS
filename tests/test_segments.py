from pathlib import Path

import pytest

from genomeos.genome.segments import SegmentParser, codon_log_odds
from genomeos.genome.signals import SignalSet

SIGNALS = Path("data/results/signals_chr21.json")


@pytest.mark.skipif(not SIGNALS.exists(), reason="needs the learned chr21 signals")
def test_parser_recovers_a_two_exon_gene_built_from_the_consensus_signals():
    sig = SignalSet.load(SIGNALS)
    d, a, k = sig.pwms["splice_donor"], sig.pwms["splice_acceptor"], sig.pwms["start_kozak"]
    # consensus windows: start (anchor = A of ATG at offset 6), donor (G of GT at offset 3),
    # acceptor (first exonic base at offset 14)
    kz, dn, ac = k.consensus.replace("n", "A"), d.consensus.replace("n", "A"), a.consensus.replace("n", "T")
    codons = ["GAG", "AAG", "CTG", "GCC", "GAC", "AAC", "ATC", "CTC"] * 14  # 112 codons, no stops
    exon1 = kz[: k.offset] + "ATG" + "".join(codons[:20]) + dn[: d.offset]  # ends with the 3 exonic bases
    intron = dn[d.offset :] + "T" * 40 + "A" * 30 + ac[: a.offset]  # GT…branch…polyY…AG
    exon2 = ac[a.offset :] + "GC" + "".join(codons[20:110]) + "TAA"  # 1 + 2 bases complete the codon
    background = "ACGT" * 400
    seq = background + exon1 + intron + exon2 + background
    lo = codon_log_odds(["".join(codons) * 3], background * 20)
    preds = SegmentParser(sig, lo, min_relative=0.5).parse("t", seq)
    assert preds, "no gene found"
    best = max(preds, key=lambda p: p.score)
    s0 = len(background) + k.offset  # position of the A of ATG
    e1 = s0 + 3 + 60 + d.offset  # donor: first intronic base (exon 1 ends with the 3 exonic consensus bases)
    a2 = e1 + len(intron)  # first base of exon 2
    stop_end = a2 + len(exon2)
    assert best.strand.value == "+"
    assert best.exons == [(s0, e1), (a2, stop_end)], best.exons
