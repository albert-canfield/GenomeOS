# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area I, the domain rule in a MAMMAL: the readers, the single-label rule, the census and the stop rule.

These are the lane's targeted tests. They build synthetic chromosomes and tiny files in `tmp_path`,
so none of them opens the git-ignored store: the 12 ATAC peak files, the 12 ENCODE RNA
quantifications, the 3 CTCF peak files, the cCRE registry and GENCODE vM21 are never read here.

What they pin is (1) that the node is `simpler_flysource`'s rule and nothing new, (2) that the
readers drop exactly what the registered exclusions say they drop and read only the columns the
manifest declares, (3) that `single_label_elements` implements the fly's "open in EXACTLY one
label" rule, (4) that the census is computable from positions alone -- no RNA value can reach it --
and (5) that the stop rule's verdict cannot be reached on a count that was not made, in either
direction.

The one test that reads a committed record uses `tracked_paths.must_be_present` and not a skip
(`tests/test_inert_skip_guards.py` docstring, `4f2652c`): the inventory is a file this lane commits,
so its absence is a broken checkout and must never read as "not run here". That test lands with
the inventory it reads, in the next commit, which is `1f0ee83`'s pattern: a guard keyed on a record
is written only once git tracks the record, so it is never an inert skip and never a red push.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos.attribution import simpler
from genomeos.attribution import simpler_flynode as fn
from genomeos.attribution import simpler_flysource as fs
from genomeos.attribution import simpler_mouse as sm
from genomeos.attribution.simpler import Chromosome, Gene

INVENTORY = Path("data/results/simpler_mouse_inventory.json")
REGISTRATION_FILE = Path("data/results/simpler_mouse_registration.json")


def _gene(gid: str, start: int, end: int, strand: str = "+", chrom: str = "chr1") -> Gene:
    return Gene(gid=gid, name=gid, chrom=chrom, start=start, end=end, strand=strand)


@pytest.fixture
def chroms() -> dict[str, Chromosome]:
    """chr1: a gene body around 10 kb, then two 5' ends far out to the right."""
    return {
        "chr1": Chromosome(
            [
                _gene("G1", 9_000, 11_000, "+"),
                _gene("G2", 200_000, 210_000, "+"),
                _gene("G3", 400_000, 410_000, "-"),
            ]
        )
    }


# ----------------------------------------------------------------- the node is the fly's, unchanged


def test_domain_of_is_flysources_domain_of(chroms):
    bounds = [50_000, 150_000, 350_000]
    for pos in (0, 60_000, 149_000, 300_000, 900_000):
        assert sm.domain_of(bounds, pos, sm.HALF_BIN) == fs.domain_of(bounds, pos, sm.HALF_BIN)


def test_boundary_rule_is_flysources_boundary_rule(chroms):
    bounds = [50_000, 150_000, 350_000]
    for pos in (30_000, 100_000, 250_000, 500_000):
        mine = sm.boundary_rule(chroms["chr1"], bounds, pos, sm.HALF_BIN)
        theirs = fs.boundary_rule(chroms["chr1"], bounds, pos, sm.HALF_BIN)
        assert (mine[0].gid if mine[0] else None) == (theirs[0].gid if theirs[0] else None)
        assert mine[1] == theirs[1]


def test_the_bin_is_flynodes_bin_and_is_not_retuned():
    assert (sm.BIN_BP, sm.HALF_BIN) == (fn.BIN_BP, fn.HALF_BIN)


def test_the_stop_rule_and_the_endpoint_are_the_worms():
    assert sm.STOP_RULE == simpler.MIN_DISCORDANT == 200
    assert sm.DETECTION_TPM == simpler.DETECTION_TPM == 5.0
    assert (sm.CLUSTER_BP, sm.BOOTSTRAP, sm.CONTROL_DRAWS, sm.ALPHA) == (
        simpler.CLUSTER_BP,
        simpler.BOOTSTRAP,
        simpler.CONTROL_DRAWS,
        simpler.ALPHA,
    )


def test_a_boundary_rule_gene_outside_the_domain_is_never_named(chroms):
    """The in-domain 5' end wins even when a nearer 5' end sits across a boundary."""
    gene, case = sm.boundary_rule(chroms["chr1"], [180_000], 150_000, sm.HALF_BIN)
    assert chroms["chr1"].nearest_tss(150_000).gid == "G2"
    assert case == fn.MOVES
    assert gene.gid == "G1"


def test_an_empty_domain_falls_back_to_nearest_tss(chroms):
    gene, case = sm.boundary_rule(chroms["chr1"], [250_000, 260_000], 255_000, sm.HALF_BIN)
    assert case == fn.EMPTY
    assert gene.gid == "G2"


# -------------------------------------------------------------- the readers read what is declared


def _bed(path: Path, rows: list[tuple[str, int, int]], extra: str = "") -> Path:
    body = "".join(f"{c}\t{s}\t{e}{extra}\n" for c, s, e in rows)
    if str(path).endswith(".gz"):
        path.write_bytes(gzip.compress(body.encode()))
    else:
        path.write_text(body)
    return path


def test_read_bed_spans_drops_contigs_that_are_not_in_the_chromosome_list(tmp_path):
    p = _bed(
        tmp_path / "p.bed.gz",
        [("chr1", 10, 20), ("chr1_GL456210_random", 10, 20), ("chrM", 10, 20), ("chrX", 5, 9)],
    )
    assert sm.read_bed_spans(p) == [("chr1", 10, 20), ("chrX", 5, 9)]


def test_read_bed_spans_reads_three_columns_and_ignores_the_narrowpeak_score(tmp_path):
    p = _bed(tmp_path / "p.bed.gz", [("chr1", 10, 20)], extra="\t.\t900\t.\t7.5\t3.1\t2.0\t5")
    assert sm.read_bed_spans(p) == [("chr1", 10, 20)]


def test_read_ccres_keeps_the_id_and_screens_own_class_only(tmp_path):
    p = tmp_path / "c.bed"
    p.write_text(
        "chr1\t100\t300\tEM10D1\tEM10E1\tdELS\n"
        "chr1\t400\t600\tEM10D2\tEM10E2\tPLS,CTCF-bound\n"
        "chrUn_GL456239\t1\t2\tEM10D3\tEM10E3\tdELS\n"
    )
    got = sm.read_ccres(p)
    assert [c.eid for c in got] == ["EM10E1", "EM10E2"]
    assert [c.kind for c in got] == [sm.DISTAL, sm.PROMOTER]
    assert got[0].mid == 200


@pytest.mark.parametrize(
    "klass,kind",
    [
        ("dELS", sm.DISTAL),
        ("dELS,CTCF-bound", sm.DISTAL),
        ("CTCF-only,CTCF-bound", sm.DISTAL),
        ("PLS", sm.PROMOTER),
        ("pELS,CTCF-bound", sm.PROMOTER),
        ("DNase-H3K4me3", sm.PROMOTER),
    ],
)
def test_the_promoter_distal_split_is_screens_class_and_not_a_basepair_distance(klass, kind):
    assert sm.CCRE("chr1", 0, 10, "e", klass).kind == kind


def test_boundaries_from_spans_takes_midpoints_sorts_and_deduplicates():
    bounds, summary = sm.boundaries_from_spans(
        [("chr1", 300, 400), ("chr1", 100, 200), ("chr1", 310, 390), ("chrX", 0, 100)]
    )
    assert bounds["chr1"] == [150, 350]
    assert bounds["chrX"] == [50]
    assert summary["boundaries"] == 3
    assert summary["chromosomes"] == 2


def test_boundaries_from_domains_takes_the_domain_edges_and_shares_an_abutting_one():
    bounds, summary = sm.boundaries_from_domains(
        [("chr1", 0, 100), ("chr1", 100, 300), ("chr1", 400, 500), ("chrUn_x", 0, 9)]
    )
    assert bounds["chr1"] == [0, 100, 300, 400, 500]
    assert summary["domains"] == 3
    assert summary["boundaries"] == 5
    assert summary["median_domain_bp"] == 100


def _book(path, sheets: dict[str, list[tuple]]):
    import openpyxl

    wb = openpyxl.Workbook()
    wb.remove(wb.active)
    for name, rows in sheets.items():
        ws = wb.create_sheet(name)
        ws.append(("chrom", "start", "end", "Compartment"))
        for r in rows:
            ws.append(r)
    wb.save(path)
    return path


def test_read_domain_sheet_reads_three_columns_and_never_the_compartment(tmp_path):
    p = _book(
        tmp_path / "b.xlsx",
        {"ncx NPC": [("chr1", 0, 100, "A"), ("chr2", 5, 9, "B"), ("chrUn_x", 1, 2, "A")]},
    )
    assert sm.read_domain_sheet(p, "ncx NPC", 3) == [("chr1", 0, 100), ("chr2", 5, 9)]


def test_read_domain_sheet_refuses_a_sheet_that_does_not_reproduce_its_recorded_row_count(tmp_path):
    p = _book(tmp_path / "b.xlsx", {"ncx NPC": [("chr1", 0, 100, "A")]})
    with pytest.raises(ValueError, match="the registration recorded 2398"):
        sm.read_domain_sheet(p, "ncx NPC", 2398)


def test_read_domain_sheet_refuses_a_workbook_without_the_registered_sheet(tmp_path):
    p = _book(tmp_path / "b.xlsx", {"something else": [("chr1", 0, 100, "A")]})
    with pytest.raises(ValueError, match="has no sheet"):
        sm.read_domain_sheet(p, "ncx NPC", 1)


def test_the_five_bonev_sheets_carry_the_counts_the_search_recorded():
    assert sm.BONEV_SHEETS["bonev_2017_ncx_npc"][:2] == ("ncx NPC", 2398)
    assert {v[1] for v in sm.BONEV_SHEETS.values()} == {2398, 2307, 2712, 2522, 2442}
    assert sm.PRIMARY_SET == "bonev_2017_ncx_npc"
    assert sm.PRIMARY_SET in sm.BONEV_SHEETS


def test_read_expression_finds_its_columns_by_header_name_not_position(tmp_path):
    a = tmp_path / "a.tsv"
    b = tmp_path / "b.tsv"
    a.write_text("gene_id\tTPM\tFPKM\nENSMUSG001.1\t9.0\t1\nERCC-00002\t999\t1\n")
    b.write_text("FPKM\tgene_id\tTPM\n1\tENSMUSG001.1\t0.5\n1\tENSMUSG002.3\t40\n")
    expr, summary = sm.read_expression({"forebrain": a, "liver": b}, ("forebrain", "liver"))
    assert expr["ENSMUSG001"] == (9.0, 0.5)
    assert expr["ENSMUSG002"] == (0.0, 40.0)
    assert "ERCC-00002" not in expr
    assert summary["rows_per_label"] == {"forebrain": 1, "liver": 2}


def test_read_expression_drops_spike_ins_and_anything_that_is_not_an_ensmusg_id(tmp_path):
    p = tmp_path / "a.tsv"
    p.write_text("gene_id\tTPM\ngSpikein_ERCC-1\t50\nENSMUSG9.1\t7\ntRNA-x\t3\n")
    expr, _ = sm.read_expression({"liver": p}, ("liver",))
    assert set(expr) == {"ENSMUSG9"}


def test_read_annotation_keeps_protein_coding_genes_and_strips_the_version(tmp_path):
    p = tmp_path / "g.gff3"
    p.write_text(
        "##gff-version 3\n"
        "chr1\tx\tgene\t100\t200\t.\t+\t.\tID=g1;gene_id=ENSMUSG001.5;gene_type=protein_coding;gene_name=Aaa\n"
        "chr1\tx\ttranscript\t100\t200\t.\t+\t.\tgene_id=ENSMUSG001.5;transcript_type=protein_coding\n"
        "chr1\tx\tgene\t300\t400\t.\t-\t.\tID=g2;gene_id=ENSMUSG002.1;gene_type=lncRNA;gene_name=Bbb\n"
        "chrUn_x\tx\tgene\t1\t2\t.\t+\t.\tID=g3;gene_id=ENSMUSG003.1;gene_type=protein_coding\n"
    )
    chroms, tss, names = sm.read_annotation(p)
    assert list(chroms) == ["chr1"]
    assert [g.gid for g in chroms["chr1"].genes] == ["ENSMUSG001"]
    assert tss == {"chr1": [100]}
    assert names == {"ENSMUSG001": "Aaa"}


def test_read_annotation_uses_the_3p_end_as_the_tss_of_a_minus_strand_transcript(tmp_path):
    p = tmp_path / "g.gff3"
    p.write_text(
        "chr1\tx\tgene\t100\t200\t.\t-\t.\tgene_id=ENSMUSG001.1;gene_type=protein_coding\n"
        "chr1\tx\ttranscript\t100\t200\t.\t-\t.\tgene_id=ENSMUSG001.1;transcript_type=protein_coding\n"
    )
    chroms, tss, _ = sm.read_annotation(p)
    assert chroms["chr1"].genes[0].tss == 200
    assert tss == {"chr1": [200]}


# ------------------------------------------------ the element rule: open in EXACTLY one of the labels


def _index(spans: list[tuple[str, int, int]]):
    by: dict[str, list[tuple[int, int]]] = {}
    for c, s, e in spans:
        by.setdefault(c, []).append((s, e))
    return sm.peak_index(by)


def test_overlaps_needs_one_base_and_is_half_open():
    idx = _index([("chr1", 100, 200)])
    assert sm.overlaps(idx, "chr1", 199, 250)
    assert sm.overlaps(idx, "chr1", 50, 101)
    assert not sm.overlaps(idx, "chr1", 200, 300)
    assert not sm.overlaps(idx, "chr1", 0, 100)
    assert not sm.overlaps(idx, "chr2", 100, 200)


def test_peak_index_merges_touching_spans():
    starts, ends = sm.peak_index({"chr1": [(10, 20), (15, 30), (100, 110)]})["chr1"]
    assert (starts, ends) == ([10, 100], [30, 110])


def test_single_label_elements_keeps_only_ccres_open_in_exactly_one_label():
    ccres = [
        sm.CCRE("chr1", 100, 200, "one", "dELS"),
        sm.CCRE("chr1", 300, 400, "two", "dELS"),
        sm.CCRE("chr1", 500, 600, "none", "dELS"),
    ]
    peaks = {
        "liver": _index([("chr1", 150, 160), ("chr1", 350, 360)]),
        "lung": _index([("chr1", 320, 330)]),
    }
    label_of, funnel = sm.single_label_elements(ccres, peaks, ("liver", "lung"))
    assert label_of == {"one": "liver"}
    assert funnel["open_in_exactly_one_label"] == 1
    assert funnel["open_in_more_than_one_label_excluded"] == 1
    assert funnel["open_in_no_label"] == 1
    assert funnel["ccres_on_the_chromosomes"] == 3


def test_elements_carries_the_single_label_and_splits_by_screens_class(chroms):
    ccres = [
        sm.CCRE("chr1", 300_000, 300_200, "d", "dELS"),
        sm.CCRE("chr1", 9_500, 9_700, "p", "PLS"),
    ]
    els, funnel = sm.elements(ccres, {"d": "liver", "p": "lung"}, chroms)
    assert sorted((e.kind, e.tissue) for e in els) == [
        (sm.DISTAL, "liver"),
        (sm.PROMOTER, "lung"),
    ]
    assert funnel["distal"] == 1
    assert funnel["promoter_class_positive_control"] == 1


def test_elements_drops_a_ccre_with_no_single_label(chroms):
    ccres = [sm.CCRE("chr1", 300_000, 300_200, "d", "dELS")]
    els, funnel = sm.elements(ccres, {}, chroms)
    assert els == []
    assert funnel["single_label_ccres"] == 0


# ------------------------------------- the census is computable from positions, and the stop rule


def _rows(chroms, pairs: list[tuple[int, str, str]]):
    els = [
        simpler.Element(chrom="chr1", start=mid - 100, end=mid + 100, kind=kind, tissue=lab)
        for mid, kind, lab in pairs
    ]
    return simpler.attribute(els, chroms)


def test_between_genes_is_the_distal_elements_in_no_gene_body(chroms):
    rows = _rows(chroms, [(10_000, sm.DISTAL, "liver"), (300_000, sm.DISTAL, "liver")])
    assert [r["mid"] for r in rows if r["domain_case"] == fn.HOST] == [10_000]
    main = sm.attribute(rows, chroms, {"chr1": [150_000]}, sm.HALF_BIN)
    assert [r["mid"] for r in sm.between_genes(main)] == [300_000]


def test_the_census_reaches_no_rna_value_at_all(chroms, monkeypatch):
    """`concordant` is the only door an RNA value can come through; the census must not open it."""
    monkeypatch.setattr(sm, "concordant", lambda *a, **k: pytest.fail("the census read an RNA value"))
    rows = _rows(chroms, [(300_000, sm.DISTAL, "liver"), (9_500, sm.PROMOTER, "liver")])
    main = sm.attribute(rows, chroms, {"chr1": [150_000]}, sm.HALF_BIN)
    cen = sm.census(main, {}, {"x": 1}, ("liver",))
    assert cen["between_genes"]["elements"] == 1
    assert cen["positive_control_promoters"] == 1


def test_the_census_verdict_is_no_source_qualifies_under_the_stop_rule(chroms):
    rows = _rows(chroms, [(300_000, sm.DISTAL, "liver")])
    main = sm.attribute(rows, chroms, {"chr1": [150_000]}, sm.HALF_BIN)
    cen = sm.census(main, {}, {}, ("liver",))
    assert cen["between_gene_elements_that_change_gene"] < sm.STOP_RULE
    assert cen["stop_rule_met"] is False
    assert cen["verdict"] == "NO_SOURCE_QUALIFIES"


def test_the_census_verdict_is_decidable_only_at_or_above_the_stop_rule():
    """Planted: a population built to move many elements must read DECIDABLE, and only then."""
    genes = [_gene("H0", 1_000, 2_000, "+")] + [
        _gene(f"H{i}", 500_000 + 10_000 * i, 500_000 + 10_000 * i + 500, "+") for i in range(1, 400)
    ]
    ch = {"chr1": Chromosome(genes)}
    mids = [400_000 + 1_000 * i for i in range(sm.STOP_RULE)]
    rows = _rows(ch, [(m, sm.DISTAL, "liver") for m in mids])
    bounds = {"chr1": [450_000]}
    main = sm.attribute(rows, ch, bounds, sm.HALF_BIN)
    cen = sm.census(main, {}, {}, ("liver",))
    moved = cen["between_gene_elements_that_change_gene"]
    assert cen["stop_rule_met"] is (moved >= sm.STOP_RULE)
    assert cen["verdict"] == ("DECIDABLE" if moved >= sm.STOP_RULE else "NO_SOURCE_QUALIFIES")
    assert cen["stop_rule"] == sm.STOP_RULE


def test_between_gene_moves_agrees_with_the_censuss_own_count(chroms):
    """The position-only cross-check in the inventory must not be a second, different count."""
    rows = _rows(
        chroms,
        [
            (300_000, sm.DISTAL, "liver"),
            (320_000, sm.DISTAL, "liver"),
            (10_000, sm.DISTAL, "liver"),
        ],
    )
    bounds = {"chr1": [150_000, 310_000]}
    main = sm.attribute(rows, chroms, bounds, sm.HALF_BIN)
    cen = sm.census(main, {}, {}, ("liver",))
    moves = fs.between_gene_moves(rows, chroms, bounds, sm.HALF_BIN)
    assert moves["between_gene_elements"] == cen["between_genes"]["elements"]
    assert moves["change_gene"] == cen["between_gene_elements_that_change_gene"]


# --------------------------------------------------------------- the endpoint and the three verdicts


def test_concordant_needs_the_threshold_and_the_highest_label():
    labels = ("a", "b", "c")
    assert sm.concordant((9.0, 1.0, 1.0), "a", labels)
    assert not sm.concordant((4.9, 1.0, 1.0), "a", labels)
    assert not sm.concordant((9.0, 9.0, 1.0), "a", labels)
    assert not sm.concordant(None, "a", labels)


def test_the_verdict_sentence_follows_the_interval_and_names_only_three_outcomes():
    assert "ahead of nearest TSS" in sm._verdict_sentence([0.5, 1.5])
    assert "Nearest TSS is ahead" in sm._verdict_sentence([-1.5, -0.5])
    assert "No difference detected" in sm._verdict_sentence([-0.5, 0.5])
    assert "No difference detected" in sm._verdict_sentence(None)


def test_compare_labels_a_zero_straddling_interval_as_no_difference(chroms):
    rows = _rows(chroms, [(300_000, sm.DISTAL, "liver"), (320_000, sm.DISTAL, "liver")])
    main = sm.attribute(rows, chroms, {"chr1": [150_000]}, sm.HALF_BIN)
    out = sm.compare(main, {}, "boundary", "nearest", ("liver",))
    assert out["delta_points"] == 0.0
    assert out["label"] == "NO_DIFFERENCE_DETECTED"


def test_a_gene_absent_from_the_rna_rows_is_not_concordant_and_is_counted(chroms):
    rows = _rows(chroms, [(300_000, sm.DISTAL, "liver")])
    main = sm.attribute(rows, chroms, {"chr1": [150_000]}, sm.HALF_BIN)
    out = sm.compare(main, {}, "boundary", "nearest", ("liver",))
    assert out["boundary_concordant"] == 0


# ----------------------------------------------------- the registration says what the brief required


def test_the_registration_names_the_stop_list_and_the_human_biosample_distinction():
    text = sm.REGISTRATION["stop_list_of_2026_09_27"]
    assert "READER_TERMS" in text
    assert "MOUSE TISSUE SOURCE IS NOT A HUMAN READER BIOSAMPLE" in text
    assert "No more biosamples" in text


def test_the_registration_fixes_the_population_the_endpoint_and_the_controls():
    reg = sm.REGISTRATION
    assert "no protein-coding gene body" in reg["population"]
    assert "transcripts per million" in reg["endpoint"]
    assert "rotated" in reg["beside_deciding_nothing"][0]
    assert "FALSIFIED" in reg["matched_control_and_falsifier"]
    assert set(reg["outcome_wording"]) == {
        "BOUNDARY_RULE_AHEAD",
        "NEAREST_AHEAD",
        "NO_DIFFERENCE_DETECTED",
    }
    assert reg["money"].startswith("none")


def test_the_registration_carries_the_two_invertebrates_registered_numbers():
    assert sm.WORM["between_gene_delta_points"] == 3.73
    assert sm.WORM["between_gene_ci95_points"] == [2.79, 4.74]
    assert sm.WORM["outron_excluded_delta_points"] == 4.47
    assert sm.FLY["brain_source_delta_points"] == 0.12
    assert sm.FLY["brain_source_ci95_points"] == [-0.25, 0.61]
    assert "+3.73" in sm.REGISTRATION["question"]
    assert "+0.12" in sm.REGISTRATION["question"]


def test_the_labels_are_twelve_distinct_tissues_and_brain_is_not_one_of_them():
    assert len(sm.LABELS) == len(set(sm.LABELS)) == 12
    assert "brain" not in sm.LABELS
    assert set(sm.ATAC_FILES) == set(sm.RNA_FILES) == set(sm.LABELS)
    assert sm.PRIMARY_SET not in sm.CTCF_FILES
    assert len(sm.CTCF_FILES) == 3


def test_every_declared_input_is_distinct_and_named_by_its_encode_accession():
    accs = [a for _e, a in sm.ATAC_FILES.values()]
    accs += [a for _e, a in sm.RNA_FILES.values()]
    accs += [a for _e, a, _w in sm.CTCF_FILES.values()]
    assert len(accs) == len(set(accs)) == 27
    assert all(a.startswith("ENCFF") for a in accs)


# ------------------------------------- the committed record, with it and not before it (4f2652c)


def test_the_committed_inventory_counts_the_between_gene_stratum_and_states_its_verdict():
    """`must_be_present`, not a skip: the inventory is tracked from this commit on."""
    tp.must_be_present(INVENTORY, was="the mouse inventory has not been written on this machine")
    inv = json.loads(INVENTORY.read_text())
    stop = inv["stop_rule"]
    assert stop["threshold"] == sm.STOP_RULE
    assert stop["met"] is (stop["between_gene_elements_that_change_gene"] >= sm.STOP_RULE)
    assert stop["verdict"] in ("DECIDABLE", "NO_SOURCE_QUALIFIES")
    assert inv["read_before_any_rna_value_was_read"] is True
    assert inv["alphagenome_requests"] == 0
    assert inv["census"]["between_genes"]["elements"] >= stop["between_gene_elements_that_change_gene"]
    assert "READER_TERMS" in inv["stop_list_of_2026_09_27"]


def test_the_committed_inventory_names_no_host_on_the_between_gene_population():
    """The whole point of this population: neither rule can name a gene body the element sits in."""
    tp.must_be_present(INVENTORY, was="the mouse inventory has not been written on this machine")
    bg = json.loads(INVENTORY.read_text())["census"]["between_genes"]
    assert bg["gene_structure_names_a_host"] == 0
    assert bg["differ_from_nearest_tss_where_the_nearest_gene_contains_the_element"] == 0


def test_every_committed_boundary_set_is_counted_on_the_same_population():
    tp.must_be_present(INVENTORY, was="the mouse inventory has not been written on this machine")
    cen = json.loads(INVENTORY.read_text())["census"]
    n = cen["between_genes"]["elements"]
    alts = cen["alternative_boundary_sets_between_genes"]
    assert len(alts) == len(sm.BONEV_SHEETS) + len(sm.CTCF_FILES) - 1
    assert all(v["elements"] == n for v in alts.values())
    assert set(cen["boundary_sets"]) == set(sm.BONEV_SHEETS) | set(sm.CTCF_FILES)


def test_the_committed_registration_was_written_before_any_rna_value_was_read():
    tp.must_be_present(REGISTRATION_FILE, was="the mouse registration has not been written on this machine")
    reg = json.loads(REGISTRATION_FILE.read_text())
    assert reg["written_before_any_rna_value_was_read"] is True
    assert reg["alphagenome_requests"] == 0
    r = reg["registration"]
    assert r["population"] == sm.REGISTRATION["population"]
    assert r["endpoint"] == sm.REGISTRATION["endpoint"]
    assert r["outcome"] == sm.REGISTRATION["outcome"]
    assert r["matched_control_and_falsifier"] == sm.REGISTRATION["matched_control_and_falsifier"]
    assert r["labels"] == list(sm.LABELS)
    assert "MOUSE TISSUE SOURCE IS NOT A HUMAN READER BIOSAMPLE" in r["stop_list_of_2026_09_27"]


def test_the_committed_registration_fixes_the_counts_every_set_reached():
    tp.must_be_present(REGISTRATION_FILE, was="the mouse registration has not been written on this machine")
    counts = json.loads(REGISTRATION_FILE.read_text())["registration"]["stop_rule_counts"]
    assert set(counts) == set(sm.BONEV_SHEETS) | set(sm.CTCF_FILES)
    assert min(counts.values()) >= sm.STOP_RULE
    assert counts[sm.PRIMARY_SET] == max(counts.values())


def test_the_registered_reading_names_the_PRIMARY_boundary_source_and_cannot_name_another():
    """The draft named the wrong source with the right count, and this is what catches that.

    `score`'s configuration filled the reading's boundary slot from a hand-written string left over
    from an earlier arrangement in which an ENCODE CTCF set was primary, so the reading would have
    called the 4,796 edges of Bonev's `ncx NPC` domains "CTCF peaks of embryonic-day-14.5 liver" --
    the primary set's count under a different set's name, the liver CTCF set carrying 49,312
    boundaries. The string is now derived from the primary set's own row, and this pins that.
    """
    short = sm.PRIMARY_SOURCE_SHORT
    sheet, domains, biosample = sm.BONEV_SHEETS[sm.PRIMARY_SET]
    assert f"`{sheet}`" in short
    assert f"{domains:,}" in short
    assert biosample in short
    assert "Bonev" in short
    # the primary set is an author-called Hi-C DOMAIN call, so no CTCF set's name may appear in it
    assert "CTCF" not in short
    for name, (_acc, _enc, ctcf_biosample) in sm.CTCF_FILES.items():
        assert name not in short
        assert ctcf_biosample not in short
    # and no OTHER Bonev sheet may be named either
    for name, (other_sheet, _n, other_biosample) in sm.BONEV_SHEETS.items():
        if name == sm.PRIMARY_SET:
            continue
        assert f"`{other_sheet}`" not in short
        assert other_biosample not in short


def test_the_committed_registration_agrees_with_the_primary_source_the_reading_names():
    tp.must_be_present(REGISTRATION_FILE, was="the mouse registration has not been written on this machine")
    reg = json.loads(REGISTRATION_FILE.read_text())
    sheet = sm.BONEV_SHEETS[sm.PRIMARY_SET][0]
    assert f"`{sheet}`" in reg["registration"]["boundary_set"]["primary"]
    assert "Bonev et al. 2017" in reg["registration"]["boundary_set"]["primary"]
    assert reg["census"]["boundary_sets"][sm.PRIMARY_SET]["sheet"] == sheet
    assert str(sm.BONEV_SHEETS[sm.PRIMARY_SET][1]) in reg["registration"]["boundary_set"]["primary"]
