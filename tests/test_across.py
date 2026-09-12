"""One locus across species: alignment summary and conserved grammar (area J step 7)."""

from genomeos.genome.motifs import parse_jaspar
from genomeos.knowledge.across import LOCI, conserved_grammar, summarise_alignment


def test_summarise_alignment_counts_identity_coverage_and_locus():
    blocks = [
        {
            "alignments": [
                {
                    "species": "homo_sapiens",
                    "seq": "ACGT-ACGT",
                    "seq_region": "7",
                    "start": 1,
                    "end": 8,
                    "strand": 1,
                },
                {
                    "species": "mus_musculus",
                    "seq": "ACGTTACGA",
                    "seq_region": "5",
                    "start": 100,
                    "end": 108,
                    "strand": -1,
                },
            ]
        }
    ]
    s = summarise_alignment(blocks, "mus_musculus", 10)
    assert s["blocks"] == 1 and s["locus"] == "5:100-108(-)"
    assert s["human_bases_aligned"] == 8 and s["coverage"] == 0.8
    assert s["identity"] == round(7 / 8, 3) and s["other_sequence"] == "ACGTTACGA"
    assert s["human_aligned_sequence"] == "ACGT-ACGT"  # gaps kept: the same-site test needs the columns
    empty = summarise_alignment([], "gallus_gallus", 10)
    assert empty["coverage"] == 0.0 and empty["identity"] is None and empty["locus"] is None


def test_conserved_grammar_requires_the_same_aligned_site():
    motifs = [
        m.prepare()
        for m in parse_jaspar(
            ">M1\tETS_TEST\nA [ 0 0 0 0 20 20 ]\nC [ 20 0 0 0 0 0 ]\n"
            "G [ 0 20 20 0 0 0 ]\nT [ 0 0 0 20 0 0 ]\n"
        )
    ]
    site, filler = "CGGTAA", "ATATATATATATATATATATATAT"
    human = filler + site + filler
    pairs = {
        "mouse": (human, human),  # the same site at the same column
        "chicken": (human, filler + filler + site),  # a hit, but 24 bases away
        "fish": (human, filler + "TTTTTT" + filler),  # no hit
    }
    g = conserved_grammar(pairs, motifs)
    assert g["factors_in_human"] == 1 and g["species_compared"] == ["mouse", "chicken", "fish"]
    row = g["rows"][0]
    assert row["factor"] == "ETS_TEST" and row["same_site_in"] == ["mouse"] and row["everywhere"] is False
    assert g["everywhere"] == []
    # with the site aligned in every species it is the locus's grammar
    g2 = conserved_grammar({"mouse": (human, human), "chicken": (human, human)}, motifs)
    assert g2["everywhere"] == ["ETS_TEST"]
    # an insertion of 20 bases in the other species before the site: the columns absorb it
    ins = "G" * 20
    gapped_human = filler + "-" * 20 + site + filler
    g3 = conserved_grammar({"frog": (gapped_human, filler + ins + site + filler)}, motifs)
    assert g3["everywhere"] == ["ETS_TEST"]
    assert set(LOCI) == {"ZRS", "HERC2_OCA2"}


def test_column_maps_handle_gaps():
    from genomeos.knowledge.across import base_at_column, column_map

    assert column_map("AC-GT") == [0, 1, 3, 4]
    assert base_at_column("AC-GT") == {0: 0, 1: 1, 3: 2, 4: 3}


def test_cli_table_accepts_no_rows():
    from genomeos.cli import _table

    assert _table([], ["factor", "human"]).splitlines() == ["factor  human", "------  -----"]
    assert _table([{"factor": "ETS1", "human": 0.9}], ["factor", "human"]).splitlines()[2] == "ETS1    0.9  "
