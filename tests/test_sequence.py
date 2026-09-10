from genomeos.genome import Locus, Sequence, Strand


def test_reverse_complement():
    assert Sequence("ATGC").reverse_complement() == "GCAT"
    assert Sequence("NNAC").reverse_complement() == "GTNN"


def test_gc_content_ignores_n():
    assert Sequence("GGCCNN").gc_content() == 1.0
    assert Sequence("ATATNN").gc_content() == 0.0


def test_locus_parse_and_str():
    l = Locus.parse("chr7:1,000,000-1,000,500(-)")
    assert l == Locus("chr7", 1_000_000, 1_000_500, Strand.MINUS)
    assert str(l) == "chr7:1000000-1000500(-)"
    assert l.length == 500


def test_telomeric_repeats():
    s = Sequence("TTAGGG" * 3 + "ACGT" + "CCCTAA" * 2)
    assert s.telomeric_repeats() == (3, 2)
