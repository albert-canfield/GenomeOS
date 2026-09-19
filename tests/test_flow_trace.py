from genomeos.coords import Locus, Strand
from genomeos.flow import trace
from genomeos.genome import Chromosome, Genome, Sequence
from genomeos.ir.model import Transcript


def _genome():
    g = Genome("t")
    #        5'UTR  exon1 CDS      intron       exon2 CDS+3'UTR
    chrom = "GGG" + "ATGAAA" + "GTTTTTTTAG" + "GATTAA" + "CCC"
    g.chromosomes["c"] = Chromosome("c", Sequence(chrom))
    return g


def test_plus_strand_trace_maps_bases_to_residues():
    g = _genome()
    t = Transcript(
        id="t1", kind="transcript", exons=[Locus("c", 0, 9), Locus("c", 19, 28)],
        cds_segments=[Locus("c", 3, 9), Locus("c", 19, 25)],
    )  # fmt: skip
    tr = trace(g, t, "G1")
    assert tr.mrna == "GGGAUGAAAGAUUAACCC"
    assert tr.utr5 == 3 and tr.utr3 == 3 and tr.cds == "AUGAAAGAUUAA"
    assert tr.protein == "MKD"
    # base 4 (the A of AAA) is codon 2 position 1 → Lys2
    r = tr.residue_of(6)
    assert r["region"] == "CDS" and r["codon"] == "AAA" and r["residue"] == 2 and r["aa3"] == "Lys"
    assert tr.residue_of(1)["region"] == "5'UTR"
    assert tr.residue_of(12)["region"] == "intron"
    assert tr.residue_of(26)["region"] == "3'UTR"
    assert tr.genomic_of_residue(3) == [19, 20, 21]
    # A>T at position 6 turns AAA into UAA: nonsense
    s = tr.substitute(6, "A", "T")
    assert s["consequence"] == "nonsense" and s["hgvs_p"] == "p.Lys2Ter" and s["hgvs_c"] == "c.4A>T"
    assert tr.substitute(8, "A", "G")["consequence"] == "synonymous"
    assert tr.substitute(19, "G", "A")["hgvs_p"] == "p.Asp3Asn"
    assert tr.substitute(12, "T", "C")["consequence"] == "intron"
    assert tr.substitute(6, "C", "T")["consequence"] == "reference mismatch"


def test_minus_strand_trace_reads_against_the_chromosome():
    g = Genome("t")
    # reverse complement of GGG ATGAAA GTTTTTTTAG GATTAA CCC
    fwd = "GGG" + "ATGAAA" + "GTTTTTTTAG" + "GATTAA" + "CCC"
    rc = Sequence(fwd).reverse_complement()
    g.chromosomes["c"] = Chromosome("c", rc)
    n = len(fwd)
    ex1 = Locus("c", n - 9, n, Strand.MINUS)
    ex2 = Locus("c", n - 28, n - 19, Strand.MINUS)
    t = Transcript(
        id="t1", kind="transcript", exons=[ex1, ex2],
        cds_segments=[Locus("c", n - 9, n - 3, Strand.MINUS), Locus("c", n - 25, n - 19, Strand.MINUS)],
    )  # fmt: skip
    tr = trace(g, t, "G1")
    assert tr.mrna == "GGGAUGAAAGAUUAACCC" and tr.protein == "MKD"
    # the genomic base that is the first A of AUG on the plus strand is a T
    pos = tr.genomic_of_residue(1)[0]
    assert str(rc)[pos] == "T"
    assert tr.residue_of(pos)["aa3"] == "Met"
    # plus-strand T>A here reads as A>U in the mRNA: AUG → UUG, start lost
    assert tr.substitute(pos, "T", "A")["consequence"] == "start_lost"


def test_trace_to_dict_is_json_ready():
    import json

    g = _genome()
    t = Transcript(id="t1", kind="transcript", exons=[Locus("c", 0, 9)])
    d = trace(g, t).to_dict()
    assert d["cds_start"] == -1 and d["protein_length"] == 0
    json.dumps(d)
