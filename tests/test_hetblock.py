# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-hetblock: does scalehet's scale heterogeneity survive a locus-block null?

What is guarded here: that every statistic, threshold, floor, count and seed is scalehet's or magcal's
BY IDENTITY and not by copy, so the only thing this lane can change about scalehet's test is the
exchangeable unit; that the block null holds every chromosome's gene count EXACTLY and moves a gene as
one block; that with one gene per block the block null IS scalehet's null bit-identically, which is the
check that the unit is the only change; that the block unit is freshrna's own function and span; that
`decide_block` is planted in all four of its directions, including the MIXED case that claims neither
answer and the degeneracy cases that claim nothing; that a planted regional slope reads SURVIVES and a
planted single slope with correlated neighbours reads WITHDRAWN, so the test is not vacuous in either
direction; that every registered reading names the withdrawal plainly and none strengthens the gain,
proposes a regional model or calls a scale an explanation; and that the run refuses an uncommitted
registration, a moved input and a second run.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest
import tracked_paths as tp

from genomeos.attribution import freshrna as fr
from genomeos.attribution import hetblock as hb
from genomeos.attribution import magcal as mc
from genomeos.attribution import scalehet as sh

ROOT = Path(__file__).resolve().parents[1]
MODULE = ROOT / "genomeos/attribution/hetblock.py"
HARNESS = ROOT / "scripts/hetblock_harness.py"
REGISTER = ROOT / "scripts/hetblock_register.py"
RUN = ROOT / "scripts/hetblock_run.py"
SCALEHET_RESULT = ROOT / "data/results/scalehet_result.json"
sys.path.insert(0, str(ROOT / "scripts"))

# THE RECORDS THIS LANE WRITES ARE READ BY THEIR OWN TESTS, WHICH SHIP WITH THEM AND NEVER BEFORE THEM:
# data/results/hetblock_harness.json -> tests/test_hetblock_harness_record.py (committed with it)
# data/results/hetblock_registration.json -> tests/test_hetblock_registration_record.py (with it)
# data/results/hetblock_result.json -> tests/test_hetblock_result_record.py (with it)
# A test committed BEFORE its record makes every sha in between red, because tracked_paths.must_be_present
# correctly RAISES while the path is untracked; lane-scalehet's own test did that at 036eec6 and three
# tests failed on the sha the push captured. The mechanism is kept and only the ORDER changes.


def rows_on(chrom: str, n: int, slope: float, start: int = 0, span: int = 300_000) -> list[dict]:
    """n synthetic genes on one chromosome with D = slope * P, thirteen cells each, spread over blocks."""
    out = []
    for i in range(n):
        p = [((i + j + start) % 7) - 3.0 for j in range(13)]
        out.append(
            {
                "gene": f"g{chrom}{i}",
                "gene_id": f"ENSG{chrom}{i}",
                "symbol": f"S{chrom}{i}",
                "chrom": chrom,
                "start": i * span,
                "n_rules": 1,
                "n_inhibiting": 0,
                "P": p,
                "D": [slope * x for x in p],
            }
        )
    return out


def parts(rows: list[dict]) -> tuple[list[float], list[float], list[str], list[int]]:
    """The per-gene sums, the admitted chromosomes and their sizes, through scalehet's own functions."""
    counts = sh.gene_counts(rows)
    adm, _ = sh.admitted_chromosomes(counts)
    dps, pps = sh.gene_sums(rows)
    return dps, pps, adm, [counts[c] for c in adm]


def statistics_of(rows: list[dict]) -> dict:
    dps, pps, adm, _ = parts(rows)
    return sh.observed_statistics(dps, pps, [r["chrom"] for r in rows], adm)


# ---- what is carried by identity -----------------------------------------------------------------
def test_every_statistic_threshold_floor_and_seed_is_scalehets_or_magcals_object():
    assert hb.MIN_GENES is sh.MIN_GENES is mc.MIN_GENES
    assert hb.CONTROL_ALPHA is sh.CONTROL_ALPHA is mc.CONTROL_ALPHA
    assert hb.MIN_CHROMOSOMES is sh.MIN_CHROMOSOMES
    assert hb.FORM_INTERVAL is sh.FORM_INTERVAL is mc.FORM_INTERVAL
    assert hb.TOLERANCE is sh.TOLERANCE
    assert hb.CELLS is sh.CELLS is mc.CELLS
    assert hb.PERMUTATIONS == sh.PERMUTATIONS == {"count": 1000, "seed": 20261012}
    assert hb.CONTROL_ALPHA == 0.05


def test_the_alpha_is_the_laxer_0_05_and_the_stricter_number_only_reported():
    assert hb.CONTROL_ALPHA == 0.05
    assert pytest.approx(0.05 / 3, abs=1e-4) == hb.STRICTER_ALPHA
    assert "DECIDES NOTHING" in hb.ALPHA
    # no outcome may turn on the stricter number
    assert hb.decide_block(0.03, 0.03) == "SURVIVES_BLOCK_NULL"
    assert hb.decide_block(0.06, 0.06) == "WITHDRAWN_NOT_ESTABLISHED"


def test_scalehets_committed_figures_are_pinned_and_not_retyped():
    """scalehet's result is ANOTHER lane's already-committed record, so it is required, not skipped."""
    tp.must_be_present(SCALEHET_RESULT)
    d = json.loads(SCALEHET_RESULT.read_text())
    w = hb.SCALEHET_COMMITTED
    assert w["q_observed"] == d["q_observed"]
    assert w["r_observed"] == d["r_observed"]
    assert w["global_scale"] == d["global_scale"]
    assert w["p_q"] == d["p_q"] and w["p_r"] == d["p_r"]
    assert w["null_q_median"] == d["null_q_summary"]["median"]
    assert w["null_q_98_33"] == d["null_q_summary"]["percentile_98_33"]
    assert w["null_r_median"] == d["null_r_summary"]["median"]
    assert w["null_r_98_33"] == d["null_r_summary"]["percentile_98_33"]
    assert w["outcome"] == d["outcome"] == "HETEROGENEOUS"


def test_this_lane_edits_no_peers_file():
    tp.must_be_present(MODULE, HARNESS, REGISTER, RUN)
    for peer in (
        "genomeos/attribution/scalehet.py",
        "genomeos/attribution/magcal.py",
        "genomeos/attribution/scalesplit.py",
        "genomeos/attribution/freshrna.py",
        "genomeos/attribution/multrna.py",
    ):
        p = ROOT / peer
        if p.exists():
            assert hb.LANE not in p.read_text(), f"{peer} names this lane: it was edited"


def test_the_q_and_r_statistics_are_scalehets_function_objects():
    src = MODULE.read_text()
    assert "sh.range_stat" in src
    assert "def heterogeneity_q" not in src, "Q must not be reimplemented"
    assert "def range_stat" not in src, "R must not be reimplemented"
    assert "def permutation_p" not in src, "the p form must not be reimplemented"
    assert "def percentile_range" not in src, "the interval must not be reimplemented"


def test_the_null_draws_q_agrees_with_scalehets_function_within_tolerance():
    """The draw's Q is scalehet's numpy expression; it is the SAME SUM as heterogeneity_q."""
    import numpy as np

    rows = rows_on("chr1", 200, 0.2) + rows_on("chr2", 200, 0.25, start=3)
    dps, pps, adm, sizes = parts(rows)
    obs = statistics_of(rows)
    a, w, g = obs["scales"], obs["pp"], obs["global_scale"]
    numpy_form = float((np.array(w) * (np.array(a) - g) ** 2).sum())
    assert abs(numpy_form - sh.heterogeneity_q(a, w, g)) <= hb.TOLERANCE
    assert len(dps) == len(pps) == sum(sizes)


# ---- the block unit ------------------------------------------------------------------------------
def test_the_block_unit_is_freshrnas_own_function_and_span():
    assert hb.BLOCK_SPAN == 1_000_000
    assert hb.ROBUST_SPAN == 5_000_000
    rows = rows_on("chr1", 10, 0.2, span=300_000)
    blocks = hb.locus_blocks(rows, hb.BLOCK_SPAN)
    theirs = fr.locus_clusters(rows, hb.BLOCK_SPAN)
    assert sorted(blocks) == sorted(theirs.values())
    assert "freshrna.locus_clusters" in hb.BLOCK_UNIT


def test_a_row_without_a_start_is_refused_and_not_silently_read_as_zero():
    rows = rows_on("chr1", 4, 0.2)
    rows[2]["start"] = None
    with pytest.raises(ValueError, match="no start"):
        hb.locus_blocks(rows)


def test_a_block_never_straddles_two_chromosomes():
    rows = rows_on("chr1", 20, 0.2, span=100_000) + rows_on("chr2", 20, 0.2, span=100_000)
    for b in hb.locus_blocks(rows):
        assert len({rows[i]["chrom"] for i in b}) == 1


def test_the_blocks_are_a_partition_of_the_genes():
    rows = rows_on("chr1", 30, 0.2, span=400_000) + rows_on("chr7", 25, 0.2, span=400_000)
    flat = sorted(i for b in hb.locus_blocks(rows) for i in b)
    assert flat == list(range(len(rows)))


def test_a_coarser_span_gives_fewer_and_larger_blocks():
    rows = rows_on("chr1", 60, 0.2, span=200_000)
    fine = hb.block_sizes(hb.locus_blocks(rows, hb.BLOCK_SPAN))
    coarse = hb.block_sizes(hb.locus_blocks(rows, hb.ROBUST_SPAN))
    assert coarse["blocks"] < fine["blocks"]
    assert coarse["largest"] >= fine["largest"]


def test_element_components_with_no_shared_element_are_one_gene_each():
    rows = rows_on("chr1", 8, 0.2)
    blocks = hb.element_blocks(rows, {})
    assert sorted(blocks) == [[i] for i in range(len(rows))]
    assert hb.block_sizes(blocks)["largest"] == 1


def test_element_components_join_genes_that_share_an_element_and_chain_through_a_third():
    rows = rows_on("chr1", 4, 0.2)
    names = [r["gene"] for r in rows]
    # e1 drives genes 0 and OUTSIDER; e2 drives OUTSIDER and gene 2: the chain must join 0 and 2
    shared = {"chr1:e1": [names[0], "OUTSIDER"], "chr1:e2": ["OUTSIDER", names[2]]}
    blocks = hb.element_blocks(rows, shared)
    joined = [b for b in blocks if len(b) > 1]
    assert joined == [[0, 2]]


def test_the_selection_rule_prefers_the_element_components_whenever_the_graph_has_an_edge():
    assert hb.MIN_SHARED_ELEMENTS == 1
    assert "ELEMENT-SHARING COMPONENT if" in hb.BLOCK_SELECTION_RULE
    assert "440,558" in hb.BLOCK_UNIT_WHY_NOT_ELEMENTS
    assert "NOT ONE DRIVES TWO OF" in hb.BLOCK_UNIT_WHY_NOT_ELEMENTS


# ---- the null ------------------------------------------------------------------------------------
def test_the_block_null_holds_every_groups_gene_count_exactly():
    rows = rows_on("chr1", 150, 0.2, span=200_000) + rows_on("chr2", 120, 0.2, span=200_000, start=2)
    dps, pps, adm, sizes = parts(rows)
    blocks = hb.locus_blocks(rows)
    for drawn in hb.draw_group_sizes(blocks, sizes, 40, 7):
        assert drawn == sizes


def test_the_null_refuses_sizes_that_do_not_sum_to_the_genes_given():
    rows = rows_on("chr1", 120, 0.2)
    dps, pps, _, _ = parts(rows)
    with pytest.raises(ValueError, match="sum to"):
        hb.block_permutation_null(dps, pps, hb.locus_blocks(rows), [100, 19], 5, 1)


def test_the_null_refuses_blocks_that_are_not_a_partition():
    rows = rows_on("chr1", 10, 0.2)
    dps, pps, _, _ = parts(rows)
    with pytest.raises(ValueError, match="partition"):
        hb.block_permutation_null(dps, pps, [[0, 1], [1, 2]], [5, 5], 3, 1)


def test_ONE_GENE_PER_BLOCK_IS_SCALEHETS_NULL_BIT_IDENTICALLY():
    """The check that the exchangeable unit is the only thing this lane changes."""
    rows = rows_on("chr1", 140, 0.2) + rows_on("chr3", 110, 0.3, start=5)
    dps, pps, adm, sizes = parts(rows)
    mine = hb.block_permutation_null(
        dps, pps, [[i] for i in range(len(dps))], sizes, 50, sh.PERMUTATIONS["seed"]
    )
    theirs = sh.permutation_null(dps, pps, sizes, 50, sh.PERMUTATIONS["seed"])
    assert mine["q"] == theirs["q"]
    assert mine["r"] == theirs["r"]
    assert mine["global_scales"] == theirs["global_scales"]


def test_a_gene_travels_as_one_block_so_its_two_sums_are_never_separated():
    rows = rows_on("chr1", 110, 0.2) + rows_on("chr4", 110, 0.4, start=1)
    dps, pps, adm, sizes = parts(rows)
    blocks = hb.locus_blocks(rows)
    null = hb.block_permutation_null(dps, pps, blocks, sizes, 30, 3)
    # the global scale is the same in every draw -- the same genes are only redistributed -- up to
    # the order of the floating-point summation, which is what TOLERANCE is for
    g = sum(dps) / sum(pps)
    assert max(abs(x - g) for x in null["global_scales"]) <= hb.TOLERANCE
    assert max(null["global_scales"]) - min(null["global_scales"]) <= hb.TOLERANCE


def test_the_null_is_reproducible_at_its_seed_and_differs_at_another():
    rows = rows_on("chr1", 130, 0.2, span=150_000) + rows_on("chr5", 130, 0.35, start=2, span=150_000)
    dps, pps, adm, sizes = parts(rows)
    blocks = hb.locus_blocks(rows)
    a = hb.block_permutation_null(dps, pps, blocks, sizes, 25, 11)
    b = hb.block_permutation_null(dps, pps, blocks, sizes, 25, 11)
    c = hb.block_permutation_null(dps, pps, blocks, sizes, 25, 12)
    assert a["q"] == b["q"]
    assert a["q"] != c["q"]


def test_at_most_one_block_per_cut_is_split_and_the_count_is_reported():
    rows = rows_on("chr1", 120, 0.2, span=100_000) + rows_on("chr2", 120, 0.2, start=3, span=100_000)
    dps, pps, adm, sizes = parts(rows)
    null = hb.block_permutation_null(dps, pps, hb.locus_blocks(rows), sizes, 20, 5)
    rep = null["split_blocks_per_draw"]
    assert rep["cuts"] == len(sizes) - 1
    assert rep["max"] <= rep["cuts"]


def test_a_group_whose_sum_of_p_squared_is_zero_is_degenerate_not_a_scale_of_zero():
    rows = rows_on("chr1", 110, 0.2)
    flat = rows_on("chr2", 110, 0.2)
    for r in flat:
        r["P"] = [0.0] * 13
        r["D"] = [0.0] * 13
    allr = rows + flat
    dps, pps, adm, sizes = parts(allr)
    null = hb.block_permutation_null(dps, pps, hb.locus_blocks(allr), sizes, 5, 1)
    assert null["degenerate"] is None or "sum of P^2 is 0" in null["degenerate"]


# ---- the decision --------------------------------------------------------------------------------
def test_both_at_or_below_alpha_reads_survives():
    assert hb.decide_block(0.001, 0.01) == "SURVIVES_BLOCK_NULL"
    assert hb.decide_block(0.05, 0.05) == "SURVIVES_BLOCK_NULL"


def test_both_above_alpha_reads_withdrawn_not_established():
    assert hb.decide_block(0.2, 0.3) == "WITHDRAWN_NOT_ESTABLISHED"
    assert hb.decide_block(0.0501, 0.9) == "WITHDRAWN_NOT_ESTABLISHED"


def test_the_two_statistics_on_opposite_sides_read_mixed_and_claim_neither():
    assert hb.decide_block(0.01, 0.3) == "MIXED"
    assert hb.decide_block(0.3, 0.01) == "MIXED"


def test_r_alone_can_never_declare_either_answer():
    for p_r in (0.0001, 0.02, 0.2, 0.9):
        assert hb.decide_block(0.01, p_r) in ("SURVIVES_BLOCK_NULL", "MIXED")
        assert hb.decide_block(0.9, p_r) in ("WITHDRAWN_NOT_ESTABLISHED", "MIXED")


def test_a_degenerate_null_reads_inconclusive_and_cannot_be_overridden_by_a_small_p():
    assert hb.decide_block(0.0001, 0.0001, "a group's sum of P^2 is 0") == "INCONCLUSIVE"


def test_too_few_chromosomes_too_few_blocks_or_one_huge_block_are_degenerate():
    assert hb.degeneracy(9, 2443, 0.003) is not None
    assert hb.degeneracy(23, 99, 0.003) is not None
    assert hb.degeneracy(23, 2443, 0.11) is not None
    assert hb.degeneracy(23, 2443, 0.003) is None
    assert hb.MIN_BLOCKS == 100 == mc.MIN_GENES
    assert hb.MAX_BLOCK_SHARE == 0.10


def test_the_robust_clause_is_registered_in_advance_and_cannot_change_the_outcome():
    fields = {
        "q5": 1.0,
        "null_q5_median": 2.0,
        "null_q5_lo": 0.5,
        "null_q5_hi": 3.0,
        "n_blocks5": 591,
    }
    agrees = hb.robust_clause(0.01, 0.02, **fields)
    disagrees = hb.robust_clause(0.01, 0.2, **fields)
    assert "AGREES" in agrees
    assert "DOES NOT AGREE" in disagrees
    assert hb.decide_block(0.01, 0.01) == "SURVIVES_BLOCK_NULL"  # unchanged by either clause


# ---- not vacuous: planted in both directions -----------------------------------------------------
def planted(rows: list[dict], count: int = 200, seed: int = 5) -> str:
    """The registered procedure end to end on synthetic rows, through the registered functions."""
    dps, pps, adm, sizes = parts(rows)
    obs = sh.observed_statistics(dps, pps, [r["chrom"] for r in rows], adm)
    blocks = hb.locus_blocks(rows)
    null = hb.block_permutation_null(dps, pps, blocks, sizes, count, seed)
    sizes_of = hb.block_sizes(blocks)
    degen = hb.degeneracy(len(adm), sizes_of["blocks"], sizes_of["largest_share"]) or null["degenerate"]
    if degen:
        return hb.decide_block(1.0, 1.0, degen)
    return hb.decide_block(sh.permutation_p(obs["q"], null["q"]), sh.permutation_p(obs["r"], null["r"]))


def test_PLANTED_a_strong_regional_slope_survives_the_block_null():
    """A real per-chromosome difference must still read SURVIVES once blocks are the unit."""
    rows = (
        rows_on("chr1", 300, 0.10, span=250_000)
        + rows_on("chr2", 300, 0.30, start=2, span=250_000)
        + rows_on("chr3", 300, 0.50, start=4, span=250_000)
        + rows_on("chr4", 300, 0.70, start=6, span=250_000)
        + rows_on("chr5", 300, 0.90, start=8, span=250_000)
        + rows_on("chr6", 300, 1.10, start=10, span=250_000)
        + rows_on("chr7", 300, 1.30, start=12, span=250_000)
        + rows_on("chr8", 300, 1.50, start=1, span=250_000)
        + rows_on("chr9", 300, 1.70, start=3, span=250_000)
        + rows_on("chr10", 300, 1.90, start=5, span=250_000)
    )
    assert planted(rows) == "SURVIVES_BLOCK_NULL"


def test_PLANTED_one_slope_everywhere_does_not_read_survives():
    """With one slope generating every gene, the block null must not manufacture heterogeneity."""
    rows = []
    for k in range(10):
        rows += rows_on(f"chr{k + 1}", 200, 0.25, start=k, span=250_000)
    assert planted(rows) in ("WITHDRAWN_NOT_ESTABLISHED", "MIXED")


def test_PLANTED_a_slope_that_varies_BY_BLOCK_and_not_by_chromosome_is_not_read_as_regional():
    """Within-block correlation with no between-chromosome structure must not read SURVIVES.

    Each 1-Mb block carries its own slope, drawn in a repeating pattern that is identical on every
    chromosome, so neighbouring genes are strongly correlated while no chromosome differs from another in
    expectation. A gene-level null would see the block structure as between-chromosome spread; a block
    null must not.
    """
    slopes = [0.05, 0.45, 0.25, 0.65, 0.15, 0.55, 0.35, 0.75]
    rows = []
    for k in range(10):
        chrom = f"chr{k + 1}"
        for b, s in enumerate(slopes):
            for g in range(25):
                p = [((g + j + k) % 7) - 3.0 for j in range(13)]
                rows.append(
                    {
                        "gene": f"g{chrom}b{b}n{g}",
                        "gene_id": f"ENSG{chrom}b{b}n{g}",
                        "symbol": f"S{chrom}b{b}n{g}",
                        "chrom": chrom,
                        "start": b * 1_000_000 + g * 10_000,
                        "n_rules": 1,
                        "n_inhibiting": 0,
                        "P": p,
                        "D": [s * x for x in p],
                    }
                )
    assert planted(rows) != "SURVIVES_BLOCK_NULL"


# ---- the readings --------------------------------------------------------------------------------
NUMBERS = {
    "unit": "1-Mb locus block",
    "n_blocks": 2443,
    "genes": "16,496",
    "chroms": 23,
    "global_scale": "+0.2150",
    "q": "166.048",
    "null_q_median": "78.9561",
    "null_q_lo": "35.4193",
    "null_q_hi": "143.929",
    "p_q": "0.123877",
    "alpha": 0.05,
    "stricter": 0.0167,
    "r": "+0.1642",
    "null_r_median": "+0.0950",
    "null_r_lo": "+0.0521",
    "null_r_hi": "+0.1706",
    "p_r": "0.234765",
    "robust": "THE 5-Mb CLAUSE.",
}


@pytest.mark.parametrize("kind", ["SURVIVES_BLOCK_NULL", "WITHDRAWN_NOT_ESTABLISHED", "MIXED"])
def test_every_registered_reading_fills_from_the_run_numbers(kind):
    text = hb.reading(kind, **NUMBERS)
    assert "{" not in text and "}" not in text
    assert len(text) > 400


def test_the_inconclusive_reading_names_its_reason():
    text = hb.reading("INCONCLUSIVE", why="a group's sum of P^2 is 0")
    assert "a group's sum of P^2 is 0" in text
    assert "No heterogeneity verdict is given" in text


def test_the_full_reading_always_carries_the_prior_exposure_clause():
    for kind in ("SURVIVES_BLOCK_NULL", "WITHDRAWN_NOT_ESTABLISHED", "MIXED"):
        assert hb.PRIOR_EXPOSURE_CLAUSE in hb.full_reading(kind, **NUMBERS)
    assert hb.PRIOR_EXPOSURE_CLAUSE in hb.full_reading("INCONCLUSIVE", why="x")


def test_the_prior_exposure_is_stated_twice_and_names_the_nine_readings_and_the_rejection():
    assert len(hb.PRIOR_EXPOSURE) == 2
    both = " ".join(hb.PRIOR_EXPOSURE)
    for named in ("item 7", "7d", "7h", "capped-sum", "M2h", "M2i", "scalesplit", "recovered", "scalehet"):
        assert named in both
    assert "ALREADY SEEN POST HOC" in both
    assert "TESTED ONCE" in both and "REJECTED" in both
    clause = hb.PRIOR_EXPOSURE_CLAUSE
    assert "ALREADY SEEN POST HOC" in clause and "ACTED ON" in clause


def test_the_withdrawal_reading_says_withdrawn_plainly_and_is_not_softened():
    text = hb.OUTCOMES["WITHDRAWN_NOT_ESTABLISHED"]
    assert "IS WITHDRAWN TO NOT ESTABLISHED" in text
    assert "should no longer carry" in text
    assert "NOT SHOWN THAT ONE CONSTANT IS THE WHOLE CORRECTION EITHER" in text
    for soft in ("perhaps", "may still be", "arguably", "weak evidence that the magnitude varies"):
        assert soft not in text.lower()


def test_both_outcomes_are_written_in_advance_and_neither_is_a_failure():
    assert "REACHABLE OUTCOME" in hb.QUESTION
    assert "FULL ANSWER" in hb.QUESTION
    assert "not a failure to find something" in hb.QUESTION
    assert set(hb.OUTCOMES) == {
        "SURVIVES_BLOCK_NULL",
        "WITHDRAWN_NOT_ESTABLISHED",
        "MIXED",
        "INCONCLUSIVE",
    }


def test_no_reading_proposes_a_regional_model():
    for text in list(hb.OUTCOMES.values()) + [hb.WHAT_IS_NOT_ASKED, hb.STOP_RULE]:
        assert re.search(r"regional model", text, re.I)
        assert not re.search(r"we (should|will|now) fit (a )?per-region", text, re.I)
    assert "no regional model is proposed" in hb.OUTCOMES["SURVIVES_BLOCK_NULL"].lower()
    assert "NO REGIONAL MODEL IS PROPOSED AND NONE IS RULED OUT" in hb.OUTCOMES["WITHDRAWN_NOT_ESTABLISHED"]


def test_the_gain_is_never_strengthened_in_any_reading():
    for text in list(hb.OUTCOMES.values()) + [hb.QUESTION, hb.WHAT_IS_NOT_ASKED] + list(hb.LIMITATIONS):
        assert "about 2%" in text or "2%" not in text
        assert not re.search(r"explains? (the|most)", text, re.I)
    for text in hb.OUTCOMES.values():
        assert "+0.2146" in text
        assert "GLOBAL AVERAGE" in text


def test_no_reading_calls_a_scale_an_explanation_of_expression():
    for text in list(hb.OUTCOMES.values()) + list(hb.LIMITATIONS):
        assert not re.search(r"scale .{0,30}explains", text, re.I)
    assert "MAGNITUDE CORRECTION AND NEVER AN EXPLANATION" in " ".join(hb.LIMITATIONS)


def test_the_limitations_name_the_element_unit_pricing_nothing_and_the_split_blocks():
    joined = " ".join(hb.LIMITATIONS)
    assert "PRICES NOTHING HERE" in joined
    assert "THE CUTS SPLIT BLOCKS" in joined
    assert "NOT PRICED BY ANY PERMUTATION OF ADDRESSES" in joined
    assert "A FAILURE TO REJECT IS NOT A PROOF OF IDENTITY" in joined


# ---- the binding ---------------------------------------------------------------------------------
def test_every_bound_field_names_a_real_constant_and_registered_fields_returns_them_all():
    fields = hb.registered_fields()
    assert set(fields) == set(hb.BOUND)
    for field, name in hb.BOUND.items():
        assert hasattr(hb, name), f"BOUND[{field}] names {name}, which does not exist"
        assert fields[field] is not None


def test_the_binding_set_covers_every_decision_the_run_could_otherwise_move():
    for name in (
        "CONTROL_ALPHA",
        "MIN_GENES",
        "MIN_CHROMOSOMES",
        "MIN_BLOCKS",
        "MAX_BLOCK_SHARE",
        "BLOCK_SPAN",
        "ROBUST_SPAN",
        "MIN_SHARED_ELEMENTS",
        "PERMUTATIONS",
        "OUTCOMES",
        "ROBUST_CLAUSES",
        "DECISION",
        "FORM_INTERVAL",
    ):
        assert name in hb.BOUND.values(), f"{name} can move the reading and is not BOUND"


def test_the_run_refuses_a_second_run_an_unregistered_start_and_a_moved_input():
    src = RUN.read_text()
    tp.must_be_present(RUN)
    assert "already exists" in src
    assert "does not exist; register and commit it first" in src
    assert "a pinned input moved" in src
    assert "the registration differs from the module" in src


def test_the_harness_computes_no_block_statistic_before_the_registration():
    src = HARNESS.read_text()
    tp.must_be_present(HARNESS)
    assert "draw_group_sizes" in src, "the exact-count check must read sizes alone"
    assert "block_permutation_null" in src  # the singleton case only, against scalehet's own null
    assert src.count("hb.block_permutation_null") == 1
    assert "no Q, no R and no p of a 1-Mb or 5-Mb block null was computed" in src
