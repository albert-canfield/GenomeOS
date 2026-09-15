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


def test_external_site_oracle_replaces_the_matrices_and_maps_both_strands():
    """With the true sites injected (feature c), the parser recovers the same gene from far weaker
    sequence signals; local() maps genomic marks to the parser's conventions on both strands."""
    from genomeos.genome.segments import Sequence
    from genomeos.predict.splice_sites import SpliceSites, score_of

    sig = SignalSet.load(SIGNALS)
    k = sig.pwms["start_kozak"]
    kz = k.consensus.replace("n", "A")
    codons = ["GAG", "AAG", "CTG", "GCC", "GAC", "AAC", "ATC", "CTC"] * 14
    # a donor and acceptor the matrices would reject: bare GT ... AG with no consensus context
    exon1 = kz[: k.offset] + "ATG" + "".join(codons[:20]) + "CCC"
    intron = "GT" + "TAATAGTGA" * 9 + "AG"  # stops in every frame: the intron cannot be read through
    exon2 = "GGC" + "".join(codons[20:110]) + "TAA"  # exon 1 is 66 bases, so exon 2 starts in frame
    background = "ACGT" * 400
    seq = background + exon1 + intron + exon2 + background
    lo = codon_log_odds(["".join(codons) * 3], background * 20)
    s0 = len(background) + k.offset
    e1 = len(background) + len(exon1)  # first intronic base (G of GT)
    a2 = e1 + len(intron)  # first exonic base of exon 2
    stop_end = a2 + len(exon2)
    assert not any(p.exons == [(s0, e1), (a2, stop_end)] for p in SegmentParser(sig, lo, 0.5).parse("t", seq))
    # genomic marks: donor on the last exonic base, acceptor on the first exonic base
    ss = SpliceSites("t", len(seq), cache=Path("/nonexistent"))
    ss.tracks["donor+"][e1 - 1] = 0.99
    ss.tracks["acceptor+"][a2] = 0.99
    preds = SegmentParser(sig, lo, 0.5, sites=ss.local).parse("t", seq)
    best = max(preds, key=lambda p: p.score)
    assert best.strand.value == "+" and best.exons == [(s0, e1), (a2, stop_end)], best.exons
    # the same gene on the minus strand: marks move to the reverse-complemented genome
    rc = str(Sequence(seq).reverse_complement())
    n = len(seq)
    ss2 = SpliceSites("t", n, cache=Path("/nonexistent"))
    ss2.tracks["donor-"][n - e1] = 0.99  # genomic-left base of the exon = its last exonic base on '-'
    ss2.tracks["acceptor-"][n - 1 - a2] = 0.99
    preds2 = SegmentParser(sig, lo, 0.5, sites=ss2.local).parse("t", rc)
    best2 = max(preds2, key=lambda p: p.score)
    assert best2.strand.value == "-" and best2.exons == [(n - stop_end, n - a2), (n - e1, n - s0)], (
        best2.exons
    )
    from genomeos.predict.splice_sites import BACKGROUND_LOG_ODDS

    assert abs(score_of(0.5) - BACKGROUND_LOG_ODDS) < 1e-9
    assert score_of(0.99) > score_of(0.5) > score_of(0.05) > 5


def test_markov_coding_prefers_the_frame_it_was_trained_on():
    from genomeos.genome.segments import MarkovCoding

    codons = ["GAG", "AAG", "CTG", "GCC", "GAC", "AAC", "ATC", "CTC"] * 40
    cds = "".join(codons)
    background = ("ACGT" * 300 + "TTGA" * 300 + "GCA" * 400) * 3
    m = MarkovCoding().learn([cds] * 3, background)
    in_frame = sum(m.codon(cds, i) for i in range(0, len(cds) - 8, 3))
    off_frame = sum(m.codon(cds, i) for i in range(1, len(cds) - 8, 3))
    assert in_frame > off_frame > -1e9
    assert m.codon(background, 100) < in_frame / (len(cds) / 3)  # background scores below the coding mean
