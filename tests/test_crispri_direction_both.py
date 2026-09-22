# SPDX-License-Identifier: AGPL-3.0-or-later
"""The upward arm's rules, above all the one that stops an unbalanced two-signed set reading as a pass."""

from __future__ import annotations

from genomeos.attribution import crispri_direction_both as cb


def _pair(measured_sign, predicted, cell="K562", gene="G", top=True):
    return {
        "cell": cell,
        "chrom": "chr1",
        "element": [1000, 1500],
        "gene": gene,
        "measured": -0.5 if measured_sign == "down" else 0.5,
        "measured_sign": measured_sign,
        "predicted": predicted,
        "predicted_sign": "down" if predicted < 0 else ("up" if predicted > 0 else "zero"),
        "was_top_target": top,
        "matches": 1,
        "matches_disagree": False,
    }


def _set(n_down, down_right, n_up, up_right, magnitude=0.5):
    out = []
    for i in range(n_down):
        out.append(_pair("down", -magnitude if i < down_right else magnitude, gene=f"D{i}"))
    for i in range(n_up):
        out.append(_pair("up", magnitude if i < up_right else -magnitude, gene=f"U{i}"))
    return out


def test_a_constant_sign_caller_scores_exactly_half_on_balanced_accuracy():
    """The whole reason balanced accuracy is the registered primary rather than raw agreement."""
    # 116 down / 36 up, the shape of the real set, answered "down" every time
    answered = _set(116, 116, 36, 0)
    assert cb.balanced(answered)["rate"] == 0.5
    # and it scores 0.763 on raw agreement, which is why raw agreement is never quoted against 0.5
    raw = cb.raw_agreement(answered)
    assert raw["rate"] == round(116 / 152, 4)
    assert raw["majority_class_rate"] == round(116 / 152, 4)
    # a constant-UP caller also scores exactly half, whichever sign it picks
    assert cb.balanced(_set(116, 0, 36, 36))["rate"] == 0.5


def test_balanced_accuracy_is_the_mean_of_the_two_sensitivities():
    answered = _set(10, 9, 10, 4)
    assert cb.sensitivity(answered, "down")["rate"] == 0.9
    assert cb.sensitivity(answered, "up")["rate"] == 0.4
    assert cb.balanced(answered)["rate"] == 0.65


def test_balanced_accuracy_refuses_a_one_signed_set():
    """The first half's set. It must not be given a number that implies a two-signed comparison."""
    ba = cb.balanced(_set(44, 41, 0, 0))
    assert ba["rate"] is None and ba["both_signs"] is False
    assert cb.judge_combined(ba)["verdict"] == "refused"


def test_a_predicted_zero_leaves_the_denominator_and_is_counted():
    answered = _set(4, 4, 0, 0) + [_pair("down", 0.0, gene="Z")]
    assert cb.sensitivity(answered, "down")["n"] == 4
    assert cb.raw_agreement(answered)["predicted_zero_excluded"] == 1


def test_the_upward_arm_fails_at_or_below_chance():
    """The registered outcome the lane exists to be able to return."""
    head = cb.sensitivity(_set(0, 0, 20, 10), "up")
    v = cb.judge_up(head, 0.46)
    assert v["verdict"] == "failed" and "downward effects only" in v["why"]
    # and below chance likewise
    assert cb.judge_up(cb.sensitivity(_set(0, 0, 20, 4), "up"), 0.46)["verdict"] == "failed"


def test_the_upward_arm_needs_both_margins_to_pass():
    """0.70 is not a pass on its own: it has to clear the model's own up-habit, interval and all."""
    seven_of_ten = cb.sensitivity(_set(0, 0, 10, 7), "up")
    assert seven_of_ten["rate"] == 0.7
    # the point clears 0.65 and the marginal by 0.24, but the Wilson lower bound does not clear 0.46
    assert cb.judge_up(seven_of_ten, 0.46)["verdict"] == "undecidable"
    # the same 0.70 against a model that says "up" more often fails the 0.10 margin instead
    assert cb.judge_up(cb.sensitivity(_set(0, 0, 200, 140), "up"), 0.62)["verdict"] == "undecidable"
    # and it passes only when both margins and the interval hold
    assert cb.judge_up(cb.sensitivity(_set(0, 0, 200, 190), "up"), 0.46)["verdict"] == "passed"


def test_the_combined_verdict_reads_the_registered_ladder():
    passed = cb.judge_combined(cb.balanced(_set(200, 190, 200, 180)))
    assert passed["verdict"] == "passed" and passed["band"] == "> 0.80"
    weak = cb.judge_combined(cb.balanced(_set(200, 190, 200, 40)))
    assert weak["band"] in ("0.55 to 0.65", "<= 0.55")
    assert weak["verdict"] in ("failed", "undecidable")
    assert cb.LADDER[weak["band"]] in weak["why"]


def test_a_failing_combined_set_says_the_headline_is_downward_only():
    """The band at or below 0.55 must say in words what it costs the passing result."""
    assert "downward effects only" in cb.LADDER["<= 0.55"]
    assert "survive" in cb.FALSIFIES or "does NOT survive" in cb.FALSIFIES


def test_value_for_gene_reads_a_gene_that_is_not_the_top_target():
    """The one fact the whole lane rests on: the cache keeps every gene, the compact table keeps one."""
    archive = {
        "E1": {
            "genes": [
                {"gene": "TOP", "by_cell": {"K562": -0.9}},
                {"gene": "OTHER", "by_cell": {"K562": 0.2, "GM12878": -0.1}},
            ]
        }
    }
    els = [{"id": "E1", "start": 1000, "end": 1500}]
    m = cb.value_for_gene(els, "chr1", "OTHER", "K562", archive)
    assert m["value"] == 0.2 and m["matches"] == 1
    # the compact-table route cannot answer the same pair
    compact = [
        {
            "id": "E1",
            "start": 1000,
            "end": 1500,
            "predicted": {"gene": "TOP"},
            "predicted_by_cell": {"K562": -0.9},
        }
    ]
    assert cb.top_target_of(compact, "OTHER", "K562") is False
    assert cb.top_target_of(compact, "TOP", "K562") is True
    # a cell with no stored value is not zero, it is unanswerable
    assert cb.value_for_gene(els, "chr1", "TOP", "GM12878", archive) is None


def test_value_for_gene_takes_the_largest_magnitude_and_flags_disagreement():
    archive = {
        "E1": {"genes": [{"gene": "G", "by_cell": {"K562": -0.2}}]},
        "E2": {"genes": [{"gene": "G", "by_cell": {"K562": 0.7}}]},
    }
    els = [{"id": "E1", "start": 1000, "end": 1500}, {"id": "E2", "start": 1100, "end": 1600}]
    m = cb.value_for_gene(els, "chr1", "G", "K562", archive)
    assert m["value"] == 0.7 and m["matches"] == 2 and m["matches_disagree"] is True


def test_the_magnitude_matched_sub_test_keeps_only_the_strong_band():
    answered = _set(4, 4, 4, 0, magnitude=0.5) + _set(6, 6, 6, 6, magnitude=0.02)
    mm = cb.magnitude_matched(answered)
    assert mm["threshold"] == 0.1
    assert mm["down"]["n"] == 4 and mm["up"]["n"] == 4
    assert mm["up"]["rate"] == 0.0  # the upward arm fails where the downward one does not: not magnitude


def test_the_magnitude_ladder_is_read_per_sign():
    lad = cb.by_magnitude(_set(4, 4, 4, 0, magnitude=0.5))
    assert lad["up"]["above_0.1"]["rate"] == 0.0
    assert lad["down"]["above_0.1"]["rate"] == 1.0
    assert lad["up"]["largest_error_magnitude"] == 0.5
    assert lad["down"]["largest_error_magnitude"] == 0.0


def test_the_registration_names_the_expected_direction_and_the_budget():
    assert "AT OR BELOW CHANCE" in cb.EXPECTED
    assert "0 of 3" in cb.EXPECTED
    assert "Requests: 0" in cb.PREREGISTERED
    assert "balanced accuracy" in cb.PREREGISTERED
    assert "no subset is searched" in cb.FALSIFIES
    assert set(cb.LADDER) == {"<= 0.55", "0.55 to 0.65", "0.65 to 0.80", "> 0.80"}


def test_the_costing_is_recomputed_from_the_tables_not_quoted_from_the_prior_lane():
    """The costing must come out of the benchmark rows, and it must say what was actually spent."""
    assert "One deletion request per distinct element" in cb.costing.__doc__
    assert "Reported whatever it comes to" in cb.costing.__doc__
    assert cb.PRIOR == {"k": 41, "n": 44, "rate": 0.9318}


def test_a_rederivation_mismatch_is_visible():
    prior = [
        {
            "chrom": "chr1",
            "element": [1000, 1500],
            "gene": "G",
            "cell": "K562",
            "predicted": -0.5,
            "predicted_sign": "down",
        }
    ]
    here = [{**_pair("down", 0.3, gene="G")}]
    import json
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "prior.json"
        p.write_text(json.dumps({"answered": prior}))
        out = cb.rederivation_check(here, p)
    assert out["checked"] and out["sign_changed"] == 1 and out["passes"] is False
    assert out["prior_headline"] == {"k": 41, "n": 44, "rate": 0.9318}
