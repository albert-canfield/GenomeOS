from genomeos.genome import Genome, Locus, Sequence, Strand, read_fasta
from genomeos.ir import Transcript
from genomeos.runtime import (
    STANDARD_CODE,
    VERTEBRATE_MITOCHONDRIAL_CODE,
    find_orfs,
    splice,
    transcribe,
    translate,
)


def test_codon_table_is_complete_and_correct():
    assert len(STANDARD_CODE) == 64
    assert STANDARD_CODE["AUG"] == "M"
    assert STANDARD_CODE["UGG"] == "W"
    assert {STANDARD_CODE[c] for c in ("UAA", "UAG", "UGA")} == {"*"}
    assert VERTEBRATE_MITOCHONDRIAL_CODE["UGA"] == "W"
    assert VERTEBRATE_MITOCHONDRIAL_CODE["AGA"] == "*"


def test_transcribe_and_translate():
    dna = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"
    assert translate(transcribe(dna)) == "MAIVMGR"
    assert translate(transcribe(dna), to_stop=False) == "MAIVMGR*KGAR*"


def test_minus_strand_transcription():
    dna = Sequence("ATGAAATAA").reverse_complement()  # gene on the minus strand
    assert translate(transcribe(dna, Strand.MINUS)) == "MK"


def test_find_orfs_demo_genome():
    seqs = read_fasta("data/demo/demo.fa")
    seq = seqs["demo_synthetic"]
    orfs = list(find_orfs(seq, chrom="demo", min_aa=50))
    assert len(orfs) >= 1
    longest = max(orfs, key=lambda o: o.length_aa)
    assert longest.locus.start == 150 and longest.length_aa == 61
    assert longest.protein.startswith("M")


def test_splice_joins_exons():
    g = Genome("t")
    from genomeos.genome import Chromosome

    #           exon1      intron        exon2
    chrom = "ATGAAA" + "GTTTTTTTAG" + "TTTTAA"
    g.chromosomes["c"] = Chromosome("c", Sequence(chrom))
    t = Transcript(id="t1", kind="transcript", exons=[Locus("c", 0, 6), Locus("c", 16, 22)])
    mrna = splice(g, t)
    assert mrna == "AUGAAAUUUUAA"
    assert translate(mrna) == "MKF"


def test_selenocysteine_table_reads_uga_as_u_but_keeps_the_final_stop():
    from genomeos.runtime.central_dogma import SELENOCYSTEINE_CODE, translate_cds

    assert translate_cds("AUGUGAGGCUAA", SELENOCYSTEINE_CODE) == "MUG"
    assert translate_cds("AUGGGCUGA", SELENOCYSTEINE_CODE) == "MG"  # UGA as terminator at the end
    assert translate_cds("AUGUGAGGCUAA") == "M"  # standard table: UGA stops
