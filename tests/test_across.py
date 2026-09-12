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
    g = conserved_grammar(pairs, motifs, families={})
    assert (
        g["factors_in_human"] == 1
        and g["families_in_human"] == 1
        and g["species_compared"] == ["mouse", "chicken", "fish"]
    )
    row = g["rows"][0]
    assert row["factor"] == "ETS_TEST" and row["same_site_in"] == ["mouse"] and row["everywhere"] is False
    assert g["everywhere"] == []
    # with the site aligned in every species it is the locus's grammar
    fam = {"ETS_TEST": {"family": "Ets-related", "class": "Tryptophan cluster factors"}}
    g2 = conserved_grammar({"mouse": (human, human), "chicken": (human, human)}, motifs, families=fam)
    assert g2["everywhere"] == ["ETS_TEST"] and g2["families_everywhere"] == {"Ets-related": ["ETS_TEST"]}
    # an insertion of 20 bases in the other species before the site: the columns absorb it
    ins = "G" * 20
    gapped_human = filler + "-" * 20 + site + filler
    g3 = conserved_grammar({"frog": (gapped_human, filler + ins + site + filler)}, motifs, families={})
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


def test_genome_coordinates_and_sites():
    from genomeos.knowledge.across import cluster_sites, to_genome

    blocks = [(0, 1000), (50, 2000)]  # 50 bases from genome 1000, then the next block from genome 2000
    assert to_genome(10, blocks) == 1010 and to_genome(60, blocks) == 2010 and to_genome(-1, blocks) is None
    held = [
        (100, 108, "HOXA9", "HOX", "mouse"),
        (100, 108, "HOXA9", "HOX", "chicken"),
        (104, 112, "CDX1", "HOX", "chicken"),
        (106, 114, "NKX6-1", "NK", "mouse"),
        (300, 310, "MEIS2", "TALE", "mouse"),  # held by mouse only
        (300, 310, "MEIS3", "TALE", "mouse"),  # a tie on start and width
        (500, 508, "ETS1", "Ets", "fish"),  # a poorly aligned species, not core
    ]
    sites = cluster_sites(held, core=["mouse", "chicken"])
    assert [(s["start"], s["end"]) for s in sites] == [(100, 114), (300, 310), (500, 508)]
    first, second, third = sites
    assert (
        first["everywhere"]
        and first["families"] == ["HOX", "NK"]
        and first["species"] == ["chicken", "mouse"]
    )
    assert first["factors_everywhere"] == ["HOXA9"]  # CDX1 and NKX6-1 hold it in one species each
    assert not second["everywhere"] and second["factors"] == ["MEIS2", "MEIS3"]
    assert not third["everywhere"] and third["species"] == ["fish"]


def test_summary_records_human_blocks():
    blocks = [
        {
            "alignments": [
                {
                    "species": "homo_sapiens",
                    "seq": "AC-GT",
                    "seq_region": "7",
                    "start": 101,
                    "end": 104,
                    "strand": 1,
                },
                {
                    "species": "danio_rerio",
                    "seq": "ACTGT",
                    "seq_region": "7",
                    "start": 5,
                    "end": 9,
                    "strand": 1,
                },
            ]
        },
        {
            "alignments": [
                {
                    "species": "homo_sapiens",
                    "seq": "TTT",
                    "seq_region": "7",
                    "start": 201,
                    "end": 203,
                    "strand": 1,
                },
                {
                    "species": "danio_rerio",
                    "seq": "TTA",
                    "seq_region": "7",
                    "start": 20,
                    "end": 22,
                    "strand": 1,
                },
            ]
        },
    ]
    s = summarise_alignment(blocks, "danio_rerio", 10)
    assert s["human_blocks"] == [(0, 100), (4, 200)]


def test_panel_statistics():
    from genomeos.knowledge.across import bh, fisher_greater, locus_group, mann_whitney_greater, panel_summary

    # Fisher: all 5 positives hold, none of 5 negatives: p = 1 / C(10, 5)
    assert abs(fisher_greater(5, 0, 0, 5) - 1 / 252) < 1e-12 and fisher_greater(0, 5, 5, 0) == 1.0
    assert (
        mann_whitney_greater([5, 6, 7, 8], [1, 2, 3, 4]) < 0.05 and mann_whitney_greater([1, 2], [5, 6]) > 0.5
    )
    assert mann_whitney_greater([], [1]) is None
    q = bh([0.01, 0.04, 0.03, 0.5])
    assert [round(x, 4) for x in q] == [0.04, 0.0533, 0.0533, 0.5]
    assert locus_group(["lb"], "positive") == "limb" and locus_group(["fb", "hb"], "positive") == "neural"
    assert locus_group(["lb", "fb"], "positive") is None and locus_group([], "negative") == "negative"

    def rec(i, group, held, sites, core=2):
        return {
            "id": f"{group}{i}",
            "group": group,
            "length": 1000,
            "core_species": ["mouse", "chicken"][:core],
            "strict_sites": sites,
            "families_held": held,
            "families_hit": held,
        }

    records = [rec(i, "limb", ["HOX", "SOX"] if i < 6 else ["HOX"], 12) for i in range(10)]  # SOX 6 of 10
    records += [rec(i, "negative", ["SOX"] if i % 2 else [], 6) for i in range(10)]
    records += [rec(99, "limb", ["HOX"], 30, core=1)]  # one core species only: not usable
    s = panel_summary(records)
    assert s["usable"] == {"limb": 10, "negative": 10}
    limb = s["groups"]["limb"]
    hox = next(r for r in limb["families"] if r["family"] == "HOX")
    assert hox["positives_held"] == 10 and hox["negatives_held"] == 0 and hox["q"] <= 0.05
    sox = next(r for r in limb["families"] if r["family"] == "SOX")
    assert sox["q"] > 0.05 and limb["density_p_greater"] < 0.05
    assert [r["family"] for r in limb["families_q05"]] == ["HOX"]
