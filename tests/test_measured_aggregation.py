# SPDX-License-Identifier: AGPL-3.0-or-later
"""R6 (external review of 2026-09-28): reporter tiles are kept, the label follows a declared rule.

The acceptance tests are the review's own words, pinned in `measured.R6_ACCEPTANCE`: one unusually
strong reporter tile cannot make an element active, and conflicting tiles stay inspectable. Synthetic
and local: no assay cache is read.
"""

from __future__ import annotations

from genomeos.attribution import measured
from genomeos.attribution.mpra import Element as MpraElement

CHROM = "chrT"


def tile(start, end, k562, strand="+", name=None):
    return MpraElement(
        CHROM, start, end, name or f"t{start}", activity={"K562": k562}, strand={"K562": strand}
    )


def row_for(layer, start, end, gene="AAA"):
    elements = [{"id": "E1", "start": start, "end": end, "predicted_coding": {"gene": gene}}]
    return measured.rows(CHROM, elements, layer)[0]


def test_acceptance_is_pinned():
    assert len(measured.R6_ACCEPTANCE) == 2
    assert "strong reporter tile" in measured.R6_ACCEPTANCE[0]


def test_one_strong_tile_among_inactive_ones_does_not_make_the_element_active():
    # five 200 bp tiles over one 250 bp element, each meeting the reciprocal-overlap rule
    tiles = [tile(1000 + 10 * i, 1200 + 10 * i, v) for i, v in enumerate([3.0, 0.2, 0.2, 0.2, 0.2])]
    layer = measured.Layer(chrom=CHROM, lentimpra=tiles)
    row = row_for(layer, 1000, 1250)
    m = row["measured"]["lentimpra"]
    assert m["label_by_cell"] == {"K562": "silent"}
    assert m["cells_active"] == [] and m["cells_silent"] == ["K562"]
    assert row["agreement"]["lentimpra"] == measured.DISAGREES
    assert len(m["tiles"]) == 5
    assert sorted(t["activity"]["K562"] for t in m["tiles"]) == [0.2, 0.2, 0.2, 0.2, 3.0]
    assert m[measured.REPORTER_MAX_FIELD] == {"K562": 3.0}
    assert m["max_activity_descriptive"] == 3.0
    assert m["activity"] == {"K562": 0.2}  # the median, descriptive
    assert m["tiles_by_cell"] == {"K562": 5} and m["tiles_active_by_cell"] == {"K562": 1}


def test_conflicting_tiles_are_neither_agreement_nor_disagreement_and_stay_listed():
    tiles = [tile(1000, 1200, 1.5, name="a"), tile(1040, 1240, 0.5, strand="-", name="b")]
    layer = measured.Layer(chrom=CHROM, lentimpra=tiles)
    row = row_for(layer, 1000, 1240)
    m = row["measured"]["lentimpra"]
    assert m["cells_conflicting"] == ["K562"] and not m["cells_active"] and not m["cells_silent"]
    a = row["agreement"]
    assert a["lentimpra"] == measured.TILES_CONFLICT
    assert a["assays_agreeing"] == 0 and a["assays_disagreeing"] == 0
    assert [(t["name"], t["strand"]["K562"], t["activity"]["K562"]) for t in m["tiles"]] == [
        ("a", "+", 1.5),
        ("b", "-", 0.5),
    ]
    text = measured.basis_text(row)
    assert "+1.50" in text and "+0.50" in text and "conflicting" in text
    counts = measured.reporter_counts([row])
    assert counts["lentimpra_tiles_conflict"] == 1
    assert counts["lentimpra_cells_where_tiles_disagree"] == 1
    census = measured.census(CHROM, [row], [{"id": "E1", "start": 1000, "end": 1240}], layer)
    assert census["agrees"]["lentimpra"] == 0 and census["disagrees"]["lentimpra"] == 0
    assert measured.pool([census])["lentimpra_tiles_conflict"] == 1


def test_a_majority_of_active_tiles_is_active():
    tiles = [tile(1000 + 10 * i, 1200 + 10 * i, v) for i, v in enumerate([1.2, 1.4, 0.1])]
    row = row_for(measured.Layer(chrom=CHROM, lentimpra=tiles), 1000, 1220)
    assert row["measured"]["lentimpra"]["cells_active"] == ["K562"]
    assert row["agreement"]["lentimpra"] == measured.AGREES


def test_a_one_tile_element_reads_as_before():
    row = row_for(measured.Layer(chrom=CHROM, lentimpra=[tile(1000, 1200, 1.3)]), 1000, 1250)
    m = row["measured"]["lentimpra"]
    assert m["activity"] == {"K562": 1.3} and m["cells_active"] == ["K562"]
    assert row["agreement"]["lentimpra"] == measured.AGREES
    text = measured.basis_text(row)
    assert "lentiMPRA log2(RNA/DNA) K562 +1.30, active in 1 of 1 cells at 1.0;" in text
    assert "median of" not in text


def test_every_assay_block_names_its_outcome_kind():
    assert set(measured.OUTCOME_KIND) == set(measured.ASSAYS)
    row = row_for(measured.Layer(chrom=CHROM, lentimpra=[tile(1000, 1200, 1.3)]), 1000, 1250)
    assert "integrated lentiviral" in row["measured"]["lentimpra"]["outcome_kind"]
    assert "episomal" not in row["measured"]["lentimpra"]["outcome_kind"]


def test_reporter_label_rule_edges():
    assert measured.reporter_label([1.0]) == "active"  # the threshold is inclusive
    assert measured.reporter_label([0.99]) == "silent"
    assert measured.reporter_label([5.0, 0.0]) == "conflicting"
    assert measured.reporter_label([5.0, 0.0, 0.0]) == "silent"


def test_satmut_repeats_are_read_and_their_disagreement_counted():
    primary = {
        1001: {"functional": True, "strong": True, "effect": 0.5},
        1002: {"functional": False, "strong": False, "effect": 0.0},
    }
    flip = {
        1001: {"functional": False, "strong": False, "effect": 0.0},
        1002: {"functional": False, "strong": False, "effect": 0.0},
    }
    e = measured.SatmutElement(CHROM, 1000, 1002, "S", 2, primary, repeat_bases=(("S-flip", flip),))
    got = measured.Layer._satmut_of(1000, 1002, [(e, 1.0)])
    assert [(r["experiment"], r["role"], r["bases_functional"]) for r in got["experiments_read"]] == [
        ("S", "primary", 1),
        ("S-flip", "repeat", 0),
    ]
    assert got["bases_measured_by_more_than_one_experiment"] == 2
    assert got["bases_where_experiments_disagree"] == 1
    assert got["bases_functional"] == 1  # the verdict stays on the primary, as declared
