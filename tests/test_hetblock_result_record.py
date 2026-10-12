# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-hetblock's result record, guarded by the test that ships WITH it and never before it.

What is guarded here, and all of it WITHOUT naming which outcome came out, so the test would have
passed on either reading: that the result names the registration it was run under and the harness that
reproduced lane-scalehet's figures first; that its observed Q, R and global scale ARE scalehet's
committed values, because only the null changed; that the governing unit is the registered 1-Mb locus
block over the registered counts and is not degenerate; that the outcome name is exactly the function
of p_Q, p_R and ALPHA that `hetblock.decide_block` fixes, recomputed here from the record's own p
values; that the reading is the REGISTERED sentence for that outcome, filled and not rewritten, and
carries the prior-exposure clause; that the 5-Mb null entered only through a registered clause and did
not decide; that the element-sharing null is reported and priced nothing; and that NO REGIONAL MODEL
was fitted and the global +0.2146 was not strengthened.

This file reads data/results/hetblock_result.json and is committed in the same commit as it
(`22e12de`). No pytest.skip is keyed on that record (`4f2652c`).
"""

from __future__ import annotations

import json
from pathlib import Path

import tracked_paths as tp

from genomeos.attribution import hetblock as hb
from genomeos.attribution import scalehet as sh

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "data/results/hetblock_result.json"
REGISTRATION = ROOT / "data/results/hetblock_registration.json"
HARNESS = ROOT / "data/results/hetblock_harness.json"


def record() -> dict:
    tp.must_be_present(RECORD)
    return json.loads(RECORD.read_text())


def governing(d: dict) -> dict:
    return d["block_nulls"]["locus_block_1mb_GOVERNING"]


def test_the_result_is_tracked_and_is_step_four_run_under_the_committed_registration() -> None:
    d = record()
    tp.must_be_present(REGISTRATION, HARNESS)
    assert d["result"] == "hetblock_result"
    assert d["lane"] == hb.LANE
    assert d["step"].startswith("4")
    assert d["registration"] == "data/results/hetblock_registration.json"
    assert d["harness"] == "data/results/hetblock_harness.json"
    assert d["registration_date"] == hb.REGISTERED


def test_the_OBSERVED_statistics_ARE_scalehets_because_only_the_null_changed() -> None:
    d = record()
    committed = hb.SCALEHET_COMMITTED
    assert abs(d["observed"]["q"] - committed["q_observed"]) <= hb.TOLERANCE
    assert abs(d["observed"]["r"] - committed["r_observed"]) <= hb.TOLERANCE
    assert abs(d["observed"]["global_scale"] - committed["global_scale"]) <= hb.TOLERANCE
    assert d["admitted_count"] == committed["admitted_count"] == 23
    assert d["genes_admitted"] == committed["genes_admitted"]
    assert d["pairs_admitted"] == committed["pairs_admitted"]
    assert d["pairs_admitted"] == d["genes_admitted"] * len(hb.CELLS)


def test_the_governing_unit_is_the_registered_1mb_locus_block_and_is_not_degenerate() -> None:
    d = record()
    assert d["governing_unit"] == "locus_block_1mb"
    g = governing(d)
    sizes = g["block_sizes"]
    # the 1-Mb block counts the harness committed BEFORE the registration, on the admitted genes
    harness = json.loads(HARNESS.read_text())["exact_gene_counts"]["checked"]["one_mb"]
    assert harness["span"] == hb.BLOCK_SPAN
    assert g["blocks"] == sizes["blocks"] == harness["blocks"]
    assert sizes["largest"] == harness["largest"]
    assert sizes["genes"] == d["genes_admitted"]
    assert sizes["blocks"] >= hb.MIN_BLOCKS
    assert sizes["largest_share"] <= hb.MAX_BLOCK_SHARE
    assert d["admitted_count"] >= hb.MIN_CHROMOSOMES
    assert d["degenerate"] is None


def test_the_null_was_drawn_at_the_REGISTERED_count_and_seed_and_nothing_else() -> None:
    d = record()
    assert d["permutations"] == {"count": sh.PERMUTATIONS["count"], "seed": sh.PERMUTATIONS["seed"]}
    harness = json.loads(HARNESS.read_text())["exact_gene_counts"]["checked"]
    for key, span in (
        ("locus_block_1mb_GOVERNING", "one_mb"),
        ("locus_block_5mb_robustness_deciding_nothing", "five_mb"),
    ):
        null = d["block_nulls"][key]
        assert null["count"] == sh.PERMUTATIONS["count"]
        assert null["seed"] == sh.PERMUTATIONS["seed"]
        assert null["block_sizes"]["blocks"] == harness[span]["blocks"]
    # strictly coarser: the robustness null has fewer and larger blocks over the same genes
    fine = d["block_nulls"]["locus_block_1mb_GOVERNING"]["block_sizes"]
    coarse = d["block_nulls"]["locus_block_5mb_robustness_deciding_nothing"]["block_sizes"]
    assert coarse["blocks"] < fine["blocks"]
    assert coarse["largest"] >= fine["largest"]
    assert coarse["genes"] == fine["genes"] == d["genes_admitted"]


def test_the_OUTCOME_NAME_is_exactly_decide_blocks_function_of_the_recorded_p_VALUES() -> None:
    """Recomputed from the record's own numbers: nothing but p_Q, p_R, ALPHA and degeneracy moved it."""
    d = record()
    g = governing(d)
    assert hb.decide_block(g["p_q"], g["p_r"], d["degenerate"]) == d["outcome"]
    assert d["outcome"] in hb.OUTCOMES


def test_the_record_says_for_itself_which_statistic_clears_the_STRICTER_alpha_and_which_does_not() -> None:
    """Outcome-agnostic: the four recorded flags must be the arithmetic, so neither p can be read up."""
    for null in record()["block_nulls"].values():
        assert null["p_q_below_alpha"] == (null["p_q"] <= hb.CONTROL_ALPHA)
        assert null["p_r_below_alpha"] == (null["p_r"] <= hb.CONTROL_ALPHA)
        assert null["p_q_below_the_stricter_alpha"] == (null["p_q"] <= hb.STRICTER_ALPHA)
        assert null["p_r_below_the_stricter_alpha"] == (null["p_r"] <= hb.STRICTER_ALPHA)
        # and clearing the LAXER alpha never implies clearing the stricter one
        if null["p_q_below_the_stricter_alpha"]:
            assert null["p_q_below_alpha"]
        if null["p_r_below_the_stricter_alpha"]:
            assert null["p_r_below_alpha"]


def test_the_decision_used_the_LAXER_alpha_and_only_reported_the_stricter_one() -> None:
    d = record()
    g = governing(d)
    alpha = hb.CONTROL_ALPHA
    assert alpha == 0.05
    assert hb.STRICTER_ALPHA == 0.0167 < alpha
    if d["outcome"] == "SURVIVES_BLOCK_NULL":
        assert g["p_q"] <= alpha and g["p_r"] <= alpha
    elif d["outcome"] == "WITHDRAWN_NOT_ESTABLISHED":
        assert g["p_q"] > alpha and g["p_r"] > alpha
    assert str(hb.STRICTER_ALPHA) in d["reading"] or d["outcome"] == "INCONCLUSIVE"


def test_the_reading_is_the_REGISTERED_sentence_for_that_outcome_filled_and_not_rewritten() -> None:
    d = record()
    expected = hb.full_reading(d["outcome"], **d["reading_numbers"])
    assert d["reading"] == expected
    assert d["reading"].startswith(hb.OUTCOMES[d["outcome"]].split("{")[0])


def test_the_reading_always_carries_the_prior_exposure_clause_stated_twice() -> None:
    d = record()
    assert hb.PRIOR_EXPOSURE_CLAUSE in d["reading"]
    twice = " ".join(d["prior_exposure"]) + " " + d["prior_exposure_clause"]
    for named in ("7d", "7h", "capped-sum", "M2h", "M2i", "scalesplit", "recovered", "scalehet"):
        assert named in twice


def test_the_5mb_null_entered_only_through_a_registered_clause_and_decided_nothing() -> None:
    d = record()
    if d["outcome"] == "INCONCLUSIVE":
        return
    robust = d["reading_numbers"]["robust"]
    prefixes = [c.split("{")[0] for c in hb.ROBUST_CLAUSES.values()]
    assert any(robust.startswith(pre) for pre in prefixes), (
        "the 5-Mb clause is one of the two REGISTERED ones, filled and not written for this result"
    )
    # and it is the clause the registered RULE picks: the same side of ALPHA or not
    five = d["block_nulls"]["locus_block_5mb_robustness_deciding_nothing"]
    same_side = (five["p_q"] <= hb.CONTROL_ALPHA) == (governing(d)["p_q"] <= hb.CONTROL_ALPHA)
    wanted = "AGREES" if same_side else "DISAGREES"
    assert robust.startswith(hb.ROBUST_CLAUSES[wanted].split("{")[0])
    # and the name of the outcome is a function of the GOVERNING unit's p alone
    g = governing(d)
    assert hb.decide_block(g["p_q"], g["p_r"], d["degenerate"]) == d["outcome"]


def test_the_element_sharing_null_is_reported_and_priced_nothing_here() -> None:
    d = record()
    el = d["block_nulls"]["element_sharing_components_which_price_nothing_here"]
    assert el["block_sizes"]["blocks"] == d["genes_admitted"], "one gene per component: it IS the gene unit"
    assert el["block_sizes"]["largest"] == 1
    assert el["block_sizes"]["single_gene_blocks"] == d["genes_admitted"]
    assert "every component is a single gene" in el["note"]
    assert "IS scalehet's gene-level null" in el["note"]


def test_how_much_the_block_unit_widened_the_null_is_a_DESCRIPTIVE_and_not_a_verdict() -> None:
    w = record()["how_much_the_block_unit_widens_the_null"]
    assert w["gene_null_median_q"] == hb.SCALEHET_COMMITTED["null_q_median"]
    assert w["gene_level_p_q_scalehet_committed"] == hb.SCALEHET_COMMITTED["p_q"]
    assert "DESCRIPTIVE AND NOT A VERDICT" in w["note"]


def test_NO_REGIONAL_MODEL_WAS_FITTED_whatever_came_out() -> None:
    d = record()
    assert d["no_regional_model_was_fitted"]
    low = d["reading"].lower()
    assert "we propose" not in low
    assert "we fit" not in low
    for banned in ("per-region scale of", "regional model of", "we therefore model"):
        assert banned not in low


def test_the_global_scale_and_its_variance_share_are_carried_and_never_strengthened() -> None:
    d = record()
    assert "+0.2146" in d["reading"]
    assert "about 2% of the variance" in d["reading"]
    assert "GLOBAL AVERAGE" in d["reading"]
    low = d["reading"].lower()
    for overstatement in ("explains expression", "proves", "demonstrates that one", "settled"):
        assert overstatement not in low


def test_the_withdrawal_if_that_is_the_outcome_is_not_softened_and_says_so_plainly() -> None:
    d = record()
    if d["outcome"] != "WITHDRAWN_NOT_ESTABLISHED":
        return
    r = d["reading"]
    assert "SCALEHET'S REJECTION IS WITHDRAWN" in r
    assert "WITHDRAWN TO NOT ESTABLISHED" in r
    assert "should no longer carry" in r
    assert "NOT SHOWN THAT ONE CONSTANT IS THE WHOLE CORRECTION EITHER" in r
    for softener in ("may suggest", "weak evidence", "borderline", "nearly significant"):
        assert softener not in r.lower()


def test_the_limitations_and_the_stop_rule_are_carried_verbatim_from_the_registration() -> None:
    d = record()
    assert d["limitations"] == list(hb.LIMITATIONS)
    assert d["stop_rule"] == hb.STOP_RULE
    assert d["decision"] == hb.DECISION
    assert d["falsifier"] == hb.FALSIFIER


def test_the_run_declares_every_reader_its_inputs_relatively_and_cost_nothing() -> None:
    d = record()
    assert d["readers"]
    assert d["alphagenome_requests"] == 0
    assert "no key" in d["money"]
    for entry in d["result_manifest"]["inputs"]:
        assert not entry["path"].startswith("/"), f"input_entry paths are RELATIVE (d72f959): {entry['path']}"
        assert entry["sha256"]
