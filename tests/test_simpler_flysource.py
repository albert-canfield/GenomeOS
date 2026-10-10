# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area I, the fly's between-gene question on a richer source: the rule at any bin, the counts, the verdict.

These are the lane's targeted tests. They build synthetic chromosomes, so none of them opens the
git-ignored store: the 519 MB RNA atlas, the DAR table and the boundary workbooks are never read
here. What they pin is (1) that the bin-parameterised rule IS simpler_flynode's rule when the bin is
simpler_flynode's, (2) that the readers drop what the exclusions say they drop, (3) that the density
sweep and the tailored ceiling behave as the inventory's argument needs them to, and (4) that the
verdict cannot pass a candidate on a count that was not made.
"""

from __future__ import annotations

import csv
import random

import pytest

from genomeos.attribution import simpler
from genomeos.attribution import simpler_fly as sf
from genomeos.attribution import simpler_flynode as fn
from genomeos.attribution import simpler_flysource as fs
from genomeos.attribution.simpler import Chromosome, Gene


def _gene(gid: str, start: int, end: int, strand: str = "+") -> Gene:
    return Gene(gid=gid, name=gid, chrom="chr2L", start=start, end=end, strand=strand)


@pytest.fixture
def chroms() -> dict[str, Chromosome]:
    """One arm: three genes with a gap wide enough to hold a boundary between each pair."""
    return {
        "chr2L": Chromosome(
            [
                _gene("a", 10_000, 12_000),
                _gene("b", 60_000, 62_000),
                _gene("c", 200_000, 202_000),
            ]
        )
    }


def _rows(mids: list[int], chroms: dict[str, Chromosome]) -> list[dict]:
    els = [simpler.Element("chr2L", m - 100, m + 100, sf.DISTAL, fs.LABELS[0]) for m in mids]
    return simpler.attribute(els, chroms)


# --------------------------------------------------------------- the rule is simpler_flynode's rule


def test_domain_of_at_flynodes_bin_is_flynodes_domain(chroms):
    bounds = [20_000, 100_000, 150_000]
    for pos in (5_000, 20_000, 30_000, 99_999, 100_250, 180_000, 400_000):
        assert fs.domain_of(bounds, pos, fn.HALF_BIN) == fn.domain_of(bounds, pos)


def test_boundary_rule_at_flynodes_bin_names_flynodes_gene(chroms):
    ch = chroms["chr2L"]
    bounds = [20_000, 100_000]
    for pos in range(1_000, 260_000, 937):
        mine = fs.boundary_rule(ch, bounds, pos, fn.HALF_BIN)
        theirs = fn.boundary_rule(ch, bounds, pos)
        assert (mine[0].gid if mine[0] else None) == (theirs[0].gid if theirs[0] else None)
        assert mine[1] == theirs[1]


def test_a_wider_bin_can_keep_the_nearest_5p_end_in_the_domain(chroms):
    """The only thing the bin changes: a 5' end just outside a narrow bin is inside a wide one."""
    ch = chroms["chr2L"]
    bounds = [12_600]
    pos = 30_000
    narrow, narrow_case = fs.boundary_rule(ch, bounds, pos, 100)
    wide, wide_case = fs.boundary_rule(ch, bounds, pos, 2_700)
    assert narrow.gid == "b" and narrow_case == fn.MOVES
    assert wide.gid == "a" and wide_case == fn.IN_DOMAIN


def test_empty_domain_falls_back_to_nearest_tss(chroms):
    ch = chroms["chr2L"]
    gene, case = fs.boundary_rule(ch, [30_000, 31_000], 30_500, 10)
    # The domain (29,990, 31,010) holds no 5' end, so the rule falls back to the arm's nearest.
    assert case == fn.EMPTY and gene.gid == "a"


# ------------------------------------------------------------------------------------- the readers


def test_read_boundary_table_prefixes_takes_midpoints_and_drops_other_arms():
    bounds, summary = fs.read_boundary_table(
        [("2L", 1_000, 1_800), ("chr3R", 5_000, 5_800), ("mitochondrion_genome", 1, 800)]
    )
    assert bounds == {"chr2L": [1_400], "chr3R": [5_400]}
    assert summary["boundaries"] == 2
    assert summary["span_widths_bp"] == {800: 2}
    assert summary["rows_off_the_used_arms"] == 1


def test_read_boundary_table_sorts_within_an_arm():
    bounds, _ = fs.read_boundary_table([("2L", 9_000, 9_800), ("2L", 1_000, 1_800)])
    assert bounds["chr2L"] == [1_400, 9_400]


def test_read_label_map_refuses_a_file_whose_declared_clusters_are_not_the_frozen_set(tmp_path):
    path = tmp_path / "cellTypeInfo.tsv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["CellType_lvl1", "RNAcluster"])
        w.writerow(["Astrocyte-like", "AST"])
    with pytest.raises(ValueError, match="frozen label set"):
        fs.read_label_map(path)


def test_read_label_map_ignores_cell_types_with_no_declared_counterpart(tmp_path):
    path = tmp_path / "cellTypeInfo.tsv"
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["CellType_lvl1", "RNAcluster"])
        for lab in fs.LABELS:
            w.writerow([lab, lab])
        w.writerow(["No counterpart", "-"])
        w.writerow(["Also none", ""])
    got = fs.read_label_map(path)
    assert "No counterpart" not in got and "Also none" not in got
    assert sorted(set(got.values())) == sorted(fs.LABELS)


def test_elements_keeps_only_regions_up_in_exactly_one_label(chroms):
    tss = {"chr2L": sorted(g.tss for g in chroms["chr2L"].genes)}
    dars = {
        "chr2L:100000-100500": {fs.LABELS[0]},
        "chr2L:120000-120500": {fs.LABELS[0], fs.LABELS[1]},
        "chrUextra:1-500": {fs.LABELS[0]},
        "chr2L:10000-10500": {fs.LABELS[2]},
    }
    els, funnel = fs.elements(dars, chroms, tss)
    kept = {e.mid: e for e in els}
    assert funnel["more_than_one_label"] == 1
    assert funnel["chromosome_not_used"] == 1
    assert 100_250 in kept and kept[100_250].kind == sf.DISTAL
    assert kept[10_250].kind == sf.PROMOTER
    assert funnel["regions_up_in_at_least_one_label"] == 4


def test_elements_carries_the_single_label_onto_the_element(chroms):
    tss = {"chr2L": sorted(g.tss for g in chroms["chr2L"].genes)}
    els, _ = fs.elements({"chr2L:100000-100500": {"AST"}}, chroms, tss)
    assert [e.tissue for e in els] == ["AST"]


# ------------------------------------------------------------------- the counts the inventory argues


def test_between_gene_moves_counts_only_elements_outside_a_gene_body(chroms):
    rows = _rows([11_000, 100_000], chroms)  # the first is inside gene 'a'
    got = fs.between_gene_moves(rows, chroms, {"chr2L": [50_000]}, fn.HALF_BIN)
    assert got["all_distal_elements"] == 2
    assert got["between_gene_elements"] == 1


def test_between_gene_moves_agrees_with_flynode_at_flynodes_bin(chroms):
    rows = _rows(list(range(5_000, 250_000, 1_111)), chroms)
    bounds = {"chr2L": [20_000, 70_000, 150_000]}
    mine = fs.between_gene_moves(rows, chroms, bounds, fn.HALF_BIN)
    theirs = fn.attribute(rows, chroms, bounds)
    between = fn.between_genes(theirs)
    assert mine["between_gene_elements"] == len(between)
    assert mine["change_gene"] == sum(1 for r in between if r["boundary_gene"] != r["nearest_gene"])


def test_place_uniform_spreads_by_arm_length_and_stays_inside_the_arm():
    arm = {"chr2L": 1_000_000, "chr3R": 3_000_000}
    bounds = fs.place_uniform(arm, 400, random.Random(1))
    assert sum(len(v) for v in bounds.values()) == 400
    assert len(bounds["chr3R"]) == 3 * len(bounds["chr2L"])
    for name, positions in bounds.items():
        assert positions == sorted(positions)
        assert all(fn.HALF_BIN <= p <= arm[name] - fn.HALF_BIN for p in positions)


def test_density_sweep_reports_min_mean_max_per_rung_and_a_stop_rule_flag(chroms):
    rows = _rows(list(range(20_000, 190_000, 2_003)), chroms)
    sweep = fs.density_sweep(rows, chroms, {"chr2L": 250_000}, counts=(4, 40), draws=3)
    assert [r["boundaries"] for r in sweep["rungs"]] == [4, 40]
    for rung in sweep["rungs"]:
        assert rung["change_gene_min"] <= rung["change_gene_mean"] <= rung["change_gene_max"]
        assert rung["mean_spacing_bp"] == 250_000 // rung["boundaries"]
    assert sweep["reaches_stop_rule_at_any_density"] is (
        sweep["best_single_draw_change_gene"] >= fs.STOP_RULE
    )
    assert sweep["draws_per_count"] == 3


def test_density_sweep_is_reproducible_from_its_seed(chroms):
    rows = _rows(list(range(20_000, 190_000, 2_003)), chroms)
    arm = {"chr2L": 250_000}
    first = fs.density_sweep(rows, chroms, arm, counts=(10,), draws=2)
    second = fs.density_sweep(rows, chroms, arm, counts=(10,), draws=2)
    assert first["rungs"] == second["rungs"]


def test_tailored_ceiling_places_one_boundary_per_gap_and_moves_at_least_as_many(chroms):
    rows = _rows([40_000, 150_000], chroms)
    ceiling = fs.tailored_ceiling(rows, chroms)
    none = fs.between_gene_moves(rows, chroms, {}, fn.HALF_BIN)
    assert ceiling["boundaries_placed"] <= ceiling["between_gene_elements"]
    assert ceiling["change_gene"] >= none["change_gene"]
    assert ceiling["min_gap_bp"] == fs.MIN_GAP_BP


def test_tailored_ceiling_skips_a_gap_narrower_than_the_bin(chroms):
    """An element 400 bp from its nearest 5' end: a boundary in that gap would sit on the 5' end."""
    rows = _rows([9_600], chroms)
    ceiling = fs.tailored_ceiling(rows, chroms)
    assert ceiling["between_gene_elements"] == 1
    assert ceiling["boundaries_placed"] == 0


# ------------------------------------------------------------------------------------- the verdict


def _cand(name: str, counted: bool, moves: int | None) -> dict:
    return {
        "name": name,
        "counted": counted,
        "between_gene_elements": 3_000,
        "between_gene_elements_that_change_gene": moves,
    }


def test_verdict_is_no_source_qualifies_when_nothing_counted_reaches_the_stop_rule():
    got = fs.verdict([_cand("few", True, fs.STOP_RULE - 1), _cand("fewer", True, 10)])
    assert got["outcome"] == "NO_SOURCE_QUALIFIES"
    assert got["qualifying"] == []
    assert got["best_counted"]["name"] == "few"
    assert got["candidates_counted"] == 2


def test_verdict_passes_a_candidate_that_reaches_the_stop_rule():
    got = fs.verdict([_cand("few", True, 10), _cand("enough", True, fs.STOP_RULE)])
    assert got["outcome"] == "SOURCE_QUALIFIES"
    assert got["qualifying"] == ["enough"]


def test_verdict_cannot_pass_a_candidate_whose_count_was_not_made():
    got = fs.verdict([_cand("big but unread", False, 10_000), _cand("few", True, 10)])
    assert got["outcome"] == "NO_SOURCE_QUALIFIES"
    assert got["candidates_not_counted"] == 1
    assert got["best_counted"]["name"] == "few"


def test_verdict_with_nothing_counted_names_no_best():
    got = fs.verdict([_cand("unread", False, None)])
    assert got["outcome"] == "NO_SOURCE_QUALIFIES"
    assert got["best_counted"] is None


def test_the_stop_rule_is_the_records_not_this_lanes():
    assert fs.STOP_RULE == fn.MIN_DISCORDANT == simpler.MIN_DISCORDANT == 200


# ---------------------------------------------------------------------------- the endpoint and wording


def test_concordant_needs_the_threshold_and_the_highest_label():
    n = len(fs.LABELS)
    high = [0.0] * n
    high[3] = fs.DETECTION_CPM + 1
    assert fs.concordant(tuple(high), fs.LABELS[3])
    assert not fs.concordant(tuple(high), fs.LABELS[4])
    low = [0.0] * n
    low[3] = fs.DETECTION_CPM - 0.01
    assert not fs.concordant(tuple(low), fs.LABELS[3])
    tied = [0.0] * n
    tied[3] = tied[4] = fs.DETECTION_CPM + 1
    assert not fs.concordant(tuple(tied), fs.LABELS[3])


def test_a_gene_absent_from_the_rna_rows_is_not_concordant():
    assert not fs.concordant(None, fs.LABELS[0])


def test_labels_are_unique_sorted_and_the_registration_carries_them():
    assert list(fs.LABELS) == sorted(set(fs.LABELS))
    assert fs.REGISTRATION["labels"] == fs.LABELS
    assert len(fs.LABELS) < fs.ACCESSIBILITY_CELL_TYPES_WITH_A_DECLARED_RNA_CLUSTER


def test_the_primary_boundary_set_is_the_only_one_marked_primary():
    primary = [k for k, v in fs.BOUNDARY_SETS.items() if v["primary"]]
    assert primary == [fs.PRIMARY_SET]
    assert fs.BOUNDARY_SETS[fs.PRIMARY_SET]["half_bin"] == fn.HALF_BIN


def test_every_registered_wording_template_formats_with_the_fill_keys():
    keys = set(_FILL)
    for name, template in fs.REGISTRATION["outcome_wording"].items():
        used = {part.split("}")[0].split(":")[0] for part in template.split("{")[1:] if "}" in part}
        assert used <= keys, f"{name} needs {sorted(used - keys)}"


def test_the_registration_names_the_route_the_source_and_both_licences():
    reg = fs.REGISTRATION
    assert "route (b)" in reg["route"]
    assert "Janssens" in reg["element_source"]
    assert "Davie" in reg["rna_source"]
    assert reg["element_licence"].startswith("not stated")
    assert "not inferred" in reg["element_licence"]


def test_the_registration_states_the_stage_gap_is_wider_not_narrower():
    assert "WIDER" in fs.REGISTRATION["boundary_set"]["stage"]


def test_the_tail_of_the_reading_refuses_the_causal_claim():
    tail = fs.REGISTRATION["outcome_wording"]["tail"]
    assert "Accessibility is not function" in tail
    assert "not that the element acts on it" in tail


#: The keys `_fill` produces. Kept here so a template that needs a key the run cannot supply fails
#: in the lane's own tests and not in the run that reads the endpoint.
_FILL = {
    "n",
    "kb",
    "kn",
    "cb",
    "cn",
    "k",
    "d",
    "lo",
    "hi",
    "mde",
    "kg",
    "cg",
    "k2",
    "d2",
    "lo2",
    "hi2",
    "N",
    "kA",
    "kbA",
    "cbA",
    "knA",
    "cnA",
    "dA",
    "loA",
    "hiA",
    "kgA",
    "cgA",
    "dgA",
    "logA",
    "higA",
    "cc",
    "p",
    "rm",
    "pr",
    "t",
    "nb",
    "miss",
    "fly",
}


def _cmp() -> dict:
    return {
        "n": 3265,
        "a_concordant": 100,
        "b_concordant": 90,
        "a_rate_pct": 3.1,
        "b_rate_pct": 2.8,
        "differ": 216,
        "delta_points": 0.31,
        "delta_ci95_points": [-0.2, 0.8],
    }


def _res() -> dict:
    return {
        "primary": _cmp(),
        "primary_vs_gene_structure": _cmp(),
        "all_distal": _cmp(),
        "all_distal_vs_gene_structure": _cmp(),
        "control": {"expected_rate": 0.05, "p_one_sided": 0.001},
        "rotation_beside": {"mean_rotated_rate_pct": 5.5, "p_one_sided": 0.6},
        "census": {"between_genes": {"mde_points_unclustered_vs_nearest_tss": 1.3}},
        "boundaries": 3121,
        "attributed_genes_missing_from_rna_rows": 7,
    }


def test_the_fill_keys_in_this_test_are_the_ones_the_run_produces():
    """If `_fill` grows or shrinks, this fails rather than letting a template silently break."""
    assert set(fs._fill(_res())) == _FILL


def test_reading_fills_every_part_and_leaves_no_placeholder():
    res = _res() | {
        "comparison": "NO_DIFFERENCE_DETECTED",
        "comparison_vs_gene_structure": "GENE_STRUCTURE_AHEAD",
        "falsifier": "NOT_FALSIFIED",
    }
    text = fs.reading(res)
    assert "{" not in text and "}" not in text
    assert "No difference detected between the boundary rule and nearest TSS." in text
    assert "the gene-structure rule is ahead of the boundary rule." in text
    assert "Accessibility is not function" in text
