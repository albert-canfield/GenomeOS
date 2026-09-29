# SPDX-License-Identifier: AGPL-3.0-or-later
"""Milestone 1.3 clause 2: the reach and class descriptions (2026-09-28 registration, lane-tssreach).
Synthetic rows only: no archive, no cache, no request."""

from __future__ import annotations

import importlib.util
import random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("c2r", ROOT / "scripts" / "clause2_reach_control.py")
c2r = importlib.util.module_from_spec(spec)
spec.loader.exec_module(c2r)
mc = c2r.mc


def synthetic(seed: int = 1):
    rng = random.Random(seed)
    rows = []
    for i in range(400):
        s = i * 25_000 + rng.randrange(0, 5_000)
        rows.append(
            {
                "id": f"E{i}",
                "start": s,
                "end": s + 300,
                "predicted_coding": "G" if rng.random() < 0.4 else None,
                "predicted": {"gene": "G", "log2_fold_change": 0.9} if rng.random() < 0.5 else None,
            }
        )
    blocks = [
        {"start": 2_000_000 + k * 1_000_000, "end": 2_000_000 + k * 1_000_000 + 60_000} for k in range(5)
    ]
    for b in blocks:
        b["length"] = b["end"] - b["start"]
    return rows, blocks


def facts_for(rows, reach_of=lambda i: 3, class_of=lambda i: "dELS"):
    return [
        {
            "reach": reach_of(i),
            "stratum": c2r.stratum_of(reach_of(i)),
            "distance": 1000 * (i + 1),
            "genes": reach_of(i),
            "class": class_of(i),
            "yes": {q: bool(t(r)) for q, t in c2r.QUESTIONS.items()},
        }
        for i, r in enumerate(rows)
    ]


# --- negatives first -------------------------------------------------------------------------------


def test_a_failed_gate_reports_no_reach_or_class_figure():
    # the registered order: without d717b28's counts and faeb0da's per-element rates, nothing is read
    assert c2r.REPRODUCE["real_unknown"]["windows_drawn"] == 44_100
    assert c2r.UNMATCHED_TRIES == 4_000  # d717b28's cap, not the matched one
    assert "no reach or class figure is reported" in c2r.PRE_REGISTRATION["gate"]


def test_no_outcome_here_can_make_clause_2_pass():
    # registered before the run: both controls are descriptions and neither is a rescue
    assert "neither control can make clause 2 pass" in c2r.PRE_REGISTRATION["what_this_cannot_do"]
    for text in c2r.PRE_REGISTRATION["readings"].values():
        assert "passes" not in text or "does not pass" in text


def test_a_cell_too_small_in_either_arm_leaves_the_standardisation_and_is_reported():
    recs = [
        {
            "elements": 2,
            "carrying": 1,
            "block": {"n": 2, "by_class": {"dELS": {"n": 2, "names_a_coding_gene": 0}}},
            "window": {"n": 100, "by_class": {"dELS": {"n": 100, "names_a_coding_gene": 50}}},
        }
    ]
    out = c2r.standardise(recs, "by_class", "names_a_coding_gene", min_cell=30)
    assert out["cells_included"] == []  # the block arm holds 2, below the registered 30
    assert out["cells_excluded"]["dELS"]["block_elements"] == 2
    assert out["standardised_difference_points"] is None  # nothing standardised, nothing claimed


# --- the window the scorer used ---------------------------------------------------------------------


def test_the_window_is_the_scorer_s_own_and_not_an_assumption():
    from genomeos.predict import chromatin_tracks, splice_sites

    assert c2r.SCORER_WINDOW == splice_sites.WINDOW == chromatin_tracks.WINDOW == 1_048_576
    assert c2r.HALF_WINDOW == c2r.SCORER_WINDOW // 2 == mc.HALF_WINDOW == 524_288


def test_reach_is_counted_on_the_element_midpoint_not_the_block_midpoint():
    tss = [1_000_000, 1_200_000, 3_000_000]
    assert mc.tss_count(tss, 1_100_000, c2r.HALF_WINDOW) == 2
    assert mc.tss_count(tss, 2_000_000, c2r.HALF_WINDOW) == 0  # 524,288 short of both sides
    assert mc.tss_distance(tss, 2_000_000) == 800_000


def test_genes_in_window_counts_a_span_that_reaches_the_window_without_its_tss():
    starts, ends = sorted([100, 2_000_000]), sorted([1_500_000, 2_000_500])
    # the gene at 100..1,500,000 has its TSS far outside a window centred on 1,400,000 but is in it
    assert c2r.genes_in_window(starts, ends, 1_400_000) == 1
    assert c2r.genes_in_window(starts, ends, 1_000_000_000) == 0


def test_the_reach_strata_are_fixed_here_and_not_read_from_the_data():
    assert c2r.stratum_of(0) == "0"
    assert c2r.stratum_of(4) == "3-4"
    assert c2r.stratum_of(9) == "7-9"
    assert c2r.stratum_of(1000) == "25+"
    assert [s[0] for s in c2r.REACH_STRATA] == [0, 1, 2, 3, 5, 7, 10, 15, 25]


# --- the draw is d717b28's --------------------------------------------------------------------------


def test_the_draw_reproduces_the_matched_control_s_window_counts_exactly():
    rows, blocks = synthetic()
    preds = {"names_a_coding_gene": lambda r: bool(r.get("predicted_coding"))}
    want = mc.matched_windows(blocks, blocks, rows, preds, None, max_tries=c2r.UNMATCHED_TRIES)
    got = c2r.reach_windows(blocks, blocks, rows, facts_for(rows), max_tries=c2r.UNMATCHED_TRIES)
    assert [r["drawn"] for r in got] == [r["drawn"] for r in want]
    assert [r["carrying"] for r in got] == [r["carrying"] for r in want]
    assert [r["elements"] for r in got] == [r["elements"] for r in want]
    assert [r["windows_yes"]["names_a_coding_gene"] for r in got] == [
        r["windows_yes"]["names_a_coding_gene"] for r in want
    ]


def test_every_element_of_every_accepted_window_is_described_once():
    rows, blocks = synthetic()
    facts = facts_for(rows, reach_of=lambda i: 2)
    got = c2r.reach_windows(blocks, blocks, rows, facts, max_tries=c2r.UNMATCHED_TRIES)
    for r in got:
        assert r["window"]["reach"] == 2 * r["window"]["n"]
        assert sum(c["n"] for c in r["window"]["by_class"].values()) == r["window"]["n"]
        assert sum(c["n"] for c in r["window"]["by_stratum"].values()) == r["window"]["n"]


# --- the standardisation ----------------------------------------------------------------------------


def test_standardising_on_a_cell_with_no_imbalance_leaves_the_crude_difference_alone():
    # one stratum only: the weights can do nothing, so standardised == crude
    recs = [
        {
            "elements": 1,
            "carrying": 1,
            "block": {"n": 100, "by_stratum": {"0": {"n": 100, "names_a_coding_gene": 10}}},
            "window": {"n": 400, "by_stratum": {"0": {"n": 400, "names_a_coding_gene": 200}}},
        }
    ]
    out = c2r.standardise(recs, "by_stratum", "names_a_coding_gene")
    assert out["standardised_difference_points"] == -40.0
    assert out["window_rate_crude"] == 0.5


def test_standardising_removes_a_composition_difference_it_should_remove():
    # the block's elements sit in the low stratum, the windows' in the high one, and inside each
    # stratum the two arms agree: the crude gap is composition and the standardised gap is zero
    recs = [
        {
            "elements": 1,
            "carrying": 1,
            "block": {
                "n": 200,
                "by_stratum": {
                    "0": {"n": 100, "names_a_coding_gene": 5},
                    "25+": {"n": 100, "names_a_coding_gene": 80},
                },
            },
            "window": {
                "n": 1100,
                "by_stratum": {
                    "0": {"n": 100, "names_a_coding_gene": 5},
                    "25+": {"n": 1000, "names_a_coding_gene": 800},
                },
            },
        }
    ]
    out = c2r.standardise(recs, "by_stratum", "names_a_coding_gene")
    assert out["standardised_difference_points"] == 0.0
    assert out["window_rate_crude"] > out["window_rate_standardised_to_the_blocks"]


def test_a_gap_that_lives_inside_every_cell_survives_the_standardisation():
    recs = [
        {
            "elements": 1,
            "carrying": 1,
            "block": {
                "n": 200,
                "by_stratum": {
                    "0": {"n": 100, "names_a_coding_gene": 5},
                    "25+": {"n": 100, "names_a_coding_gene": 10},
                },
            },
            "window": {
                "n": 200,
                "by_stratum": {
                    "0": {"n": 100, "names_a_coding_gene": 45},
                    "25+": {"n": 100, "names_a_coding_gene": 50},
                },
            },
        }
    ]
    out = c2r.standardise(recs, "by_stratum", "names_a_coding_gene")
    assert out["standardised_difference_points"] == -40.0


# --- the registered readings ------------------------------------------------------------------------


def test_the_half_the_gap_threshold_is_the_one_the_earlier_controls_used():
    assert c2r.POOLED_GAP == -34.86 and c2r.PER_BLOCK_GAP == -28.77
    assert c2r.HALF_POOLED_GAP == -17.43


def test_far_lower_reach_that_explains_the_gap_reads_as_the_artefact_case():
    pooled = {"reach_ratio_block_over_window": 0.2, "zero_reach_share_difference_in_points": 30.0}
    r = c2r.reading(
        pooled,
        {"ci95_over_blocks": [-20.0, -10.0]},
        {"closes_more_than_half_the_gap": True, "standardised_difference_points": -5.0},
    )
    assert r["outcome"] == "reach_short_and_explains"
    assert r["clause_2"].startswith("stays not met")


def test_comparable_reach_reads_as_the_sequence_and_leaves_the_number_alone():
    pooled = {"reach_ratio_block_over_window": 0.99, "zero_reach_share_difference_in_points": 0.3}
    r = c2r.reading(
        pooled,
        {"ci95_over_blocks": [-0.4, 0.5]},
        {"closes_more_than_half_the_gap": False, "standardised_difference_points": -33.0},
    )
    assert r["outcome"] == "reach_comparable"
    assert "about the sequence, not the window" in r["text"]


def test_reach_that_differs_but_does_not_explain_is_its_own_reading():
    pooled = {"reach_ratio_block_over_window": 0.3, "zero_reach_share_difference_in_points": 12.0}
    r = c2r.reading(
        pooled,
        {"ci95_over_blocks": [-9.0, -6.0]},
        {"closes_more_than_half_the_gap": False, "standardised_difference_points": -28.0},
    )
    assert r["outcome"] == "reach_short_but_gap_survives"


def test_the_incoherent_corner_is_registered_as_an_error_not_a_finding():
    pooled = {"reach_ratio_block_over_window": 1.0, "zero_reach_share_difference_in_points": 0.1}
    r = c2r.reading(
        pooled,
        {"ci95_over_blocks": [-0.2, 0.2]},
        {"closes_more_than_half_the_gap": True, "standardised_difference_points": -2.0},
    )
    assert r["outcome"] == "reach_comparable_but_gap_closes"
    assert "implementation error" in r["text"]


def test_the_registration_names_both_controls_and_costs_no_request():
    for key in ("quantity_a_class", "quantity_b_reach", "estimators", "interval", "thresholds", "cost"):
        assert key in c2r.PRE_REGISTRATION
    assert c2r.PRE_REGISTRATION["cost"].startswith("0 AlphaGenome requests")
    assert c2r.PRE_REGISTRATION["registered"] == "2026-09-28"
