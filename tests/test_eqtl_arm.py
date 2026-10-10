# SPDX-License-Identifier: AGPL-3.0-or-later
"""M5b: the fine-mapped eQTL arm, its exclusion, its matched windows, its test and its sizing
(genomeos/attribution/eqtl_arm.py)."""

from __future__ import annotations

import ast
import itertools
import random
import struct
from pathlib import Path

import numpy as np
import pytest

from genomeos.attribution import eqtl_arm as ea

DAPG_HOW = "fetched whole from UCSC (gbdb/hg38/gtex/eQtl/gtexDapg.bb) by scripts/eqtl_arm_calibrate.py's lane"


def blocks(*spans):
    return [{"chrom": c, "start": s, "end": e} for c, s, e in spans]


# ---- the exclusion -------------------------------------------------------------------------------


def test_the_exclusion_holds_every_base_of_a_block_and_its_margin_and_nothing_beyond():
    ex = ea.Exclusion(blocks(("chr1", 10_000, 20_000)), margin=1_000)
    assert ex.holds("chr1", 9_000) and ex.holds("chr1", 20_999) and ex.holds("chr1", 15_000)
    assert not ex.holds("chr1", 8_999) and not ex.holds("chr1", 21_000)
    assert not ex.holds("chr2", 15_000)
    assert ex.hits("chr1", 0, 9_001) and ex.hits("chr1", 20_999, 30_000)
    assert not ex.hits("chr1", 0, 9_000) and not ex.hits("chr1", 21_000, 30_000)
    assert ex.hits("chr1", 0, 100_000)  # an interval spanning the whole block


def test_the_exclusion_merges_blocks_whose_margins_meet():
    ex = ea.Exclusion(blocks(("chr1", 10_000, 20_000), ("chr1", 21_500, 30_000)), margin=1_000)
    assert ex.spans["chr1"] == [(9_000, 31_000)]
    assert ex.holds("chr1", 21_000) and ex.blocks == 2
    assert ex.digest() == ea.Exclusion(blocks(("chr1", 21_500, 30_000), ("chr1", 10_000, 20_000))).digest()


def test_the_default_margin_is_one_kilobase():
    assert ea.EXCLUSION_MARGIN == 1_000
    assert ea.Exclusion(blocks(("chr1", 5_000, 6_000))).spans["chr1"] == [(4_000, 7_000)]


# ---- the readers: position first -----------------------------------------------------------------


def record(chrom_id, start, end, fields):
    rest = "\t".join(fields).encode()
    return struct.pack("<III", chrom_id, start, end) + rest + b"\0"


def dapg_fields(pos_text, gene="ENSG1.1", name="GENE1", tissue="Whole_Blood", pip="0.9"):
    vals = {f: "x" for f in ea.DAPG_FIELDS[3:]}
    vals.update(eqtlPos=pos_text, geneId=gene, geneName=name, tissue=tissue, pip=pip)
    return [vals[f] for f in ea.DAPG_FIELDS[3:]]


def test_eqtl_pos_reads_the_tracks_one_based_position_as_a_zero_based_base():
    assert ea.eqtl_pos("chr21:5097921-5097921") == ("chr21", 5097920)
    assert ea.eqtl_pos("chr21:5097921") == ("chr21", 5097920)
    assert ea.eqtl_pos("garbage") is None


def test_a_record_whose_variant_lies_in_the_exclusion_is_dropped_before_its_fields_are_split():
    ex = ea.Exclusion(blocks(("chr1", 10_000, 20_000)), margin=1_000)
    pos_field = ea.DAPG_FIELDS.index("eqtlPos") - 3
    raw = (
        record(0, 5_000, 50_000, dapg_fields("chr1:15001-15001", name="INSIDE_SENTINEL"))
        + record(0, 5_000, 50_000, dapg_fields("chr1:30001-30001", name="OUTSIDE"))
        # the span crosses the block but the variant lies outside: kept, on the variant alone
        + record(0, 2_000, 30_000, dapg_fields("chr1:2001-2001", name="SPAN_CROSSES"))
        + record(0, 9_500, 9_600, dapg_fields("chr1:9501-9501", name="IN_THE_MARGIN"))
        + record(1, 5_000, 6_000, dapg_fields("chr1:5001-5001", name="OTHER_CHROM_ID"))
    )
    got = list(ea.decode_records(raw, 0, "chr1", ex, pos_field, ea.DAPG_FIELDS))
    assert [(p, r["geneName"]) for p, r in got] == [(30_000, "OUTSIDE"), (2_000, "SPAN_CROSSES")]
    assert not any("SENTINEL" in str(r) or "MARGIN" in str(r) for _, r in got)


def test_the_position_is_cut_out_before_the_remainder_is_decoded(monkeypatch):
    """Only records that passed the position check reach the full split: count the full decodes."""
    ex = ea.Exclusion(blocks(("chr1", 10_000, 20_000)), margin=0)
    pos_field = ea.DAPG_FIELDS.index("eqtlPos") - 3
    raw = b"".join(record(0, 1, 2, dapg_fields(f"chr1:{p + 1}-{p + 1}")) for p in (5, 15_000, 16_000, 25_000))
    seen = []
    real = ea.eqtl_pos

    def spy(text):
        seen.append(text)
        return real(text)

    monkeypatch.setattr(ea, "eqtl_pos", spy)
    got = list(ea.decode_records(raw, 0, "chr1", ex, pos_field, ea.DAPG_FIELDS))
    assert [p for p, _ in got] == [5, 25_000]
    assert len(seen) == 4 and all(t.startswith("chr1:") for t in seen)  # the position field and nothing else


def test_record_starts_reads_the_start_of_each_record_on_the_chromosome():
    raw = record(3, 10, 11, ["a"]) + record(4, 20, 21, ["b"]) + record(3, 30, 32, ["c\td"])
    assert ea.record_starts(raw, 3) == [10, 30]


def test_the_margin_is_the_one_every_earlier_dapg_read_used():
    tree = ast.parse(Path("scripts/eqtl_targets.py").read_text())
    margin = next(
        n.value.value
        for n in tree.body
        if isinstance(n, ast.Assign) and any(getattr(t, "id", None) == "MARGIN" for t in n.targets)
    )
    assert margin == ea.MARGIN


# ---- units and windows ---------------------------------------------------------------------------


def pair(start, end, gene="G", tss=100_000, regulated=True, chrom="chr1", cells=("K562",)):
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "gene": gene,
        "tss": tss,
        "regulated": regulated,
        "cells": cells,
    }


def test_scaled_gives_the_element_with_its_margin_or_a_centred_window():
    assert ea.scaled(1_000, 1_600, None) == (500, 2_100)
    assert ea.scaled(10_000, 10_600, 4_000) == (8_300, 12_300)
    assert ea.scaled(10_000, 10_600, 1_000) == (9_500, 11_100)  # never shorter than the element


def test_units_are_dropped_for_the_exclusion_the_tss_the_cis_window_and_a_negative_on_a_positive():
    ex = ea.Exclusion(blocks(("chr1", 50_000, 60_000)), margin=1_000)
    pairs = [
        pair(30_000, 30_500),  # kept
        pair(60_200, 60_700),  # touches the widened block
        pair(99_000, 99_500),  # holds the TSS with its guard
        pair(2_000_000, 2_000_500),  # beyond the cis window
        pair(30_100, 30_400, regulated=False),  # a negative on a positive of its gene
        pair(30_100, 30_400, gene="H", regulated=False),  # a negative of another gene: kept
    ]
    units, dropped = ea.build_units(pairs, ex, None, "element")
    assert [(u.start, u.gene) for u in units] == [(29_500, "G"), (29_600, "H")]
    assert dropped == {
        "touches_the_exclusion": 1,
        "holds_the_tss": 1,
        "beyond_the_cis_window": 1,
        "negative_overlapping_a_positive_of_its_gene": 1,
    }


def test_element_pairs_counts_a_pair_once_and_calls_it_positive_if_any_screen_did():
    class P:
        def __init__(self, cell, regulated, tss=500):
            self.chrom, self.start, self.end, self.gene = "chr1", 10, 20, "G"
            self.cell, self.regulated, self.tss = cell, regulated, tss

    got = ea.element_pairs([P("K562", False), P("GM12878", True), P("HCT116", False, tss=None)])
    assert got == [
        {
            "chrom": "chr1",
            "start": 10,
            "end": 20,
            "gene": "G",
            "tss": 500,
            "regulated": True,
            "cells": ("GM12878", "K562"),
        }
    ]


def density(step=100):
    """A synthetic common-variant count: one variant every `step` bases."""
    return lambda chrom, s, e: e // step - s // step


def test_every_drawn_window_meets_every_matching_rule():
    ex = ea.Exclusion(blocks(("chr1", 140_000, 150_000)), margin=1_000)
    u = ea.Unit("k", "chr1", 129_500, 131_000, "G", 100_000, True, ("K562",), (130_000, 130_500))
    count = density()
    u.common = count("chr1", u.start, u.end)
    avoid = [(68_000, 72_000)]
    got = ea.draw_windows(u, ea.WINDOWS + 1, count, ex, avoid, 10_000_000, random.Random(1))
    assert got is not None and len(got) == ea.WINDOWS + 1
    band = max(ea.DISTANCE_BAND * u.distance, ea.MIN_BAND_LENGTHS * u.length)
    for s, e in got:
        assert e - s == u.length
        assert abs(abs((s + e) // 2 - u.tss) - u.distance) <= band + 1
        assert abs(count("chr1", s, e) - u.common) <= ea.count_tolerance(u.common)
        assert not ex.hits("chr1", s, e) and not ea.holds_tss(s, e, u.tss)
        assert not ea.overlaps(s, e, u.start, u.end) and not ea.overlaps(s, e, *avoid[0])
    for (a, b), (c, d) in itertools.combinations(got, 2):
        assert not ea.overlaps(a, b, c, d)


def test_the_draw_gives_none_rather_than_an_unmatched_window():
    u = ea.Unit("k", "chr1", 129_500, 131_000, "G", 100_000, True, ("K562",), (130_000, 130_500))
    u.common = 1_000  # no window of this length can hold that many at the synthetic density
    assert ea.draw_windows(u, 3, density(), ea.Exclusion([]), [], 10_000_000, random.Random(1)) is None


def test_the_draw_is_reproducible_from_its_seed():
    u = ea.Unit("k", "chr1", 129_500, 131_000, "G", 100_000, True, ("K562",), (130_000, 130_500))
    u.common = 15
    a = ea.draw_windows(u, 5, density(), ea.Exclusion([]), [], 10_000_000, random.Random("s"))
    b = ea.draw_windows(u, 5, density(), ea.Exclusion([]), [], 10_000_000, random.Random("s"))
    assert a == b


def test_flanks_are_the_two_adjacent_windows_clear_of_the_exclusion_and_the_tss():
    ex = ea.Exclusion(blocks(("chr1", 5_000, 6_000)), margin=0)
    u = ea.Unit("k", "chr1", 7_000, 8_000, "G", 100_000, True, (), (7_000, 8_000))
    assert ea.flanks_of(u, ex, 1_000_000) == [(6_000, 7_000), (8_000, 9_000)]
    u2 = ea.Unit("k", "chr1", 5_500, 6_500, "G", 100_000, True, (), (5_500, 6_500))
    assert ea.flanks_of(u2, ex, 1_000_000) == [(6_500, 7_500)]


# ---- scoring -------------------------------------------------------------------------------------


def fine():
    return ea.Fine(
        pos=np.array([100, 150, 200, 900], dtype=np.int64),
        gene=np.array([0, 1, 0, 0], dtype=np.int32),
        tissue=np.array([0, 1, 1, 0], dtype=np.int16),
        pip=np.array([0.9, 0.8, 0.2, 0.6], dtype=np.float32),
        genes=[("ENSG1", "OLD1"), ("ENSG2", "B")],
        tissues=["Whole_Blood", "Liver"],
    )


def test_gene_codes_answer_to_the_tracks_name_and_the_current_symbol():
    codes = ea.gene_codes(fine(), {"ENSG1": "NEW1"})
    assert codes["OLD1"] == {0} and codes["NEW1"] == {0} and codes["B"] == {1}


def test_score_names_the_gene_only_at_the_threshold_and_keeps_its_tissues():
    f = fine()
    assert ea.score(f, {0}, 0, 160) == {
        "gene": True,
        "any": True,
        "gene_any_pip": True,
        "tissues": ["Whole_Blood"],
    }
    assert ea.score(f, {0}, 140, 300)["gene"] is False  # only a 0.2 row of gene 0 there
    assert ea.score(f, {0}, 140, 300)["gene_any_pip"] is True
    assert ea.score(f, set(), 0, 1_000)["gene"] is False
    assert ea.egene(f, {1}) and not ea.egene(f, set())


# ---- the test and the sizing ---------------------------------------------------------------------


def test_upper_tail_is_exact_against_enumeration():
    q = np.array([0.1, 0.5, 0.25, 1.0, 0.0, 0.3])
    for s in range(0, 8):
        brute = 0.0
        for bits in itertools.product((0, 1), repeat=len(q)):
            if sum(bits) >= s:
                brute += float(np.prod([p if b else 1 - p for p, b in zip(q, bits, strict=True)]))
        assert ea.upper_tail(q, s) == pytest.approx(brute, abs=1e-12)


def test_upper_tail_normal_branch_agrees_with_the_exact_one():
    rng = np.random.default_rng(0)
    q = rng.uniform(0.05, 0.3, 450)
    t = int(q.sum()) + 5
    approx = ea.upper_tail(q, t)  # 450 informative terms: the normal branch
    dist = np.zeros(len(q) + 1)
    dist[0] = 1.0
    for i, p in enumerate(q):
        dist[1 : i + 2] = dist[1 : i + 2] * (1 - p) + dist[0 : i + 1] * p
        dist[0] *= 1 - p
    assert approx == pytest.approx(float(dist[t:].sum()), abs=0.01)
    assert ea.upper_tail(q, 0) == 1.0 and ea.upper_tail(q, 451) == 0.0


def test_the_conditional_test_holds_its_level_when_the_block_is_another_window():
    rng = np.random.default_rng(1)
    rejections = 0
    trials = 400
    for _ in range(trials):
        rate = rng.uniform(0.0, 0.2, 150)
        hits = rng.random((150, ea.WINDOWS + 1)) < rate[:, None]
        if ea.conditional_p(hits[:, 0].astype(float), hits[:, 1:].sum(1).astype(float)) <= ea.ALPHA:
            rejections += 1
    assert rejections / trials <= ea.ALPHA + 0.03


def test_icc_is_one_for_identical_windows_and_near_zero_for_independent_ones():
    assert ea.icc(np.array([[1, 1, 1], [0, 0, 0], [1, 1, 1], [0, 0, 0]], dtype=float)) == pytest.approx(1.0)
    rng = np.random.default_rng(2)
    assert abs(ea.icc((rng.random((2000, 10)) < 0.3).astype(float))) < 0.02


def synthetic_rows(n, y_rate, w_rate, rng):
    return [
        {"y": bool(rng.random() < y_rate), "w": [bool(rng.random() < w_rate) for _ in range(ea.WINDOWS + 1)]}
        for _ in range(n)
    ]


def test_power_rises_with_the_share_right_and_holds_level_with_none_right():
    rng = np.random.default_rng(3)
    pos = synthetic_rows(400, 0.4, 0.05, rng)
    right, wrong = ea.arrays(pos), ea.arrays(pos, "w0")
    hi = ea.power(right, wrong, 100, 1.0, rng, 200)
    lo = ea.power(right, wrong, 100, 0.25, rng, 200)
    none = ea.power(right, wrong, 100, 0.0, rng, 400)
    assert hi > 0.95 and hi >= lo and none <= ea.ALPHA + 0.03


def test_cluster_interval_moves_units_of_a_gene_together():
    rng = np.random.default_rng(4)
    groups = [f"g{i // 5}" for i in range(200)]
    values = [float(i // 5 % 2) for i in range(200)]
    lo, hi = ea.cluster_interval(groups, values, rng, 500)
    assert lo < 0.5 < hi
    assert ea.cluster_interval(groups[:20], values[:20], rng, 50) is None  # fewer than ten genes


def sizes(power_half, fooled):
    half = {"power_at_the_available_blocks": power_half, "smallest_grid_blocks_reaching_80": 300}
    return {
        "window": {"by_share_right": {"0.5": half}},
        "tested_negative": {"rejection_rate_at_the_available_blocks_when_none_is_right": fooled},
    }


def test_the_fair_test_rule_needs_both_power_and_a_control_that_is_not_fooled():
    assert ea.fair_test(sizes(0.85, 0.06), 331)["fair"] is True
    low = ea.fair_test(sizes(0.5, 0.06), 331)
    assert low["fair"] is False and low["why_not"][0].startswith("(1)")
    fooled = ea.fair_test(sizes(0.9, 0.2), 331)
    assert fooled["fair"] is False and fooled["why_not"][0].startswith("(2)")
    assert ea.fair_test(sizes(0.9, 0.05), 0)["fair"] is False


def test_the_rule_text_names_both_conditions_and_its_bound():
    assert "f = 1/2" in ea.FAIR_TEST_RULE and "0.10" in ea.FAIR_TEST_RULE
    assert pytest.approx(2 * ea.ALPHA) == ea.FOOLED_BOUND


# ---- the blocks' attribution ---------------------------------------------------------------------


def test_block_attribution_takes_the_strongest_moving_element_inside_the_block():
    bl = [
        {"start": 1_000, "end": 5_000, "tier": "constrained_unknown"},
        {"start": 6_000, "end": 9_000, "tier": "constrained_unknown", "copy": True},
        {"start": 10_000, "end": 12_000, "tier": "regulatory"},
    ]

    def el(s, gene, lf, tissue="liver"):
        return {
            "start": s,
            "end": s + 200,
            "predicted": {"gene": gene, "log2_fold_change": lf, "tissue": tissue},
        }

    els = [el(1_100, "A", 0.15), el(2_000, "B", -0.4, "lung"), el(3_000, "C", 0.05), el(7_000, "D", 1.0)]
    got = ea.block_attribution("chr1", bl, els, 0.1)
    assert len(got) == 1
    assert (got[0]["gene"], got[0]["tissue"], got[0]["element"]) == ("B", "lung", [2_000, 2_200])


def test_available_counts_blocks_with_a_gene_near_enough_and_outside_the_block():
    rows_ = [
        {"chrom": "chr1", "start": 10_000, "end": 20_000, "gene": "A", "element": [12_000, 12_200]},
        {"chrom": "chr1", "start": 10_000, "end": 20_000, "gene": None, "element": None},
        {"chrom": "chr1", "start": 10_000, "end": 20_000, "gene": "FAR", "element": [12_000, 12_200]},
        {"chrom": "chr1", "start": 10_000, "end": 20_000, "gene": "IN", "element": [12_000, 12_200]},
        {"chrom": "chr1", "start": 10_000, "end": 20_000, "gene": "NONE", "element": [12_000, 12_200]},
    ]
    tss = {"A": [("chr1", 60_000)], "FAR": [("chr1", 5_000_000)], "IN": [("chr1", 18_000)]}
    got = ea.available(rows_, tss)
    assert got["naming_a_gene"] == 4
    assert got["usable_at_the_block_scale"] == 1 and got["usable_at_the_element_scale"] == 2
    assert got["not_usable"] == {
        "names_no_gene": 1,
        "beyond_the_cis_window": 1,
        "block_holds_the_tss": 1,
        "gene_without_a_tss_on_the_chromosome": 1,
    }


# ---- the track as held ---------------------------------------------------------------------------


@pytest.mark.needs_local_data("data/knowledge/eqtl_arm/gtexDapg.bb", how=DAPG_HOW)
def test_the_held_track_has_the_expected_fields_and_decodes_nothing_inside_an_exclusion():
    ex = ea.Exclusion(blocks(("chr21", 30_000_000, 35_000_000)), margin=0)
    f = ea.DapgTrack(ex).chromosome("chr21")
    assert len(f.pos) > 0 and np.all(np.diff(f.pos) >= 0)
    assert not np.any((f.pos >= 30_000_000) & (f.pos < 35_000_000))
    assert float(f.pip.min()) >= 0.0 and float(f.pip.max()) <= 1.0
    assert "Whole_Blood" in f.tissues
