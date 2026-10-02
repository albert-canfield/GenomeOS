# SPDX-License-Identifier: AGPL-3.0-or-later
"""The link-level direction test: its registration, its predicates and its estimators."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import cell2, fresh
from genomeos.attribution import direction_link as dl
from genomeos.attribution import increases as inc
from genomeos.attribution import repress2 as rp

REGISTRATION = Path("data/results/direction_link_registration.json")


# ---- the floors and the locus convention are imported, never chosen here -------------------------


def test_floors_are_imported_and_not_set_here():
    assert dl.POSITIVE_FLOOR is fresh.POSITIVE_FLOOR
    assert dl.LOCUS_FLOOR is fresh.LOCUS_FLOOR
    assert dl.LOCUS_FLOOR == cell2.POOLED_LOCUS_FLOOR
    assert dl.FLOORS["neither_chosen_here"] is True
    assert dl.FLOORS["applied"].endswith("never to the two arms pooled")


def test_locus_convention_is_cell2s_own_words_and_says_it_is_not_independence():
    assert dl.LOCUS_RULE is cell2.INDEPENDENT_LOCUS_RULE
    assert dl.LOCUS_SPAN == cell2.INDEPENDENT_LOCUS_SPAN == 1_000_000
    assert "not established biological independence" in dl.LOCUS_RULE
    assert dl.NOT_BIOLOGICAL_INDEPENDENCE is inc.NOT_BIOLOGICAL_INDEPENDENCE


def test_inherited_readings_are_imported_verbatim_and_not_restated():
    assert dl.INHERITED_REPRESS2_NO_GO is inc.INHERITED_NO_GO
    assert dl.INHERITED_REPRESS2_NO_GO == "the measured layer holds no repression call to test"
    assert dl.INHERITED_INCREASE_GO is inc.GO


def test_the_increase_population_is_recorded_as_unbalanced_and_partly_held_out():
    assert dl.INHERITED_FIGURES["increase_links_by_cell"] == {"K562": 39, "WTC11": 6, "HCT116": 3}
    assert sum(dl.INHERITED_FIGURES["increase_links_by_cell"].values()) == 48
    assert dl.INHERITED_FIGURES["increases_behind_the_links_by_split"]["heldout"] == 23
    assert "EVALUATION only" in dl.SPLIT_RULE
    assert "not balanced across cell types" in dl.NOT_BALANCED


def test_the_floors_are_on_counts_and_the_effect_sizes_are_recorded_beside_them():
    m = dl.INHERITED_FIGURES["largest_increase_per_link_over_the_40_links_that_file_shows"]
    assert (m["min"], m["median"], m["max"]) == (0.0102, 0.0987, 4.2101)
    assert (m["at_or_above_0.10"], m["of"]) == (20, 40)
    assert "on COUNTS" in dl.FLOORS_ARE_ON_COUNTS
    assert "not a statement that the effects behind it are large" in dl.FLOORS_ARE_ON_COUNTS


def test_the_decrease_arm_does_not_inherit_the_coverage_confusion():
    assert "EXTRACTOR'S OWN CONSTRUCTION and not model coverage" in dl.DECREASE_ARM_IS_NOT_MODEL_COVERAGE
    assert "and_on_an_attributed_element" in dl.DECREASE_ARM_IS_NOT_MODEL_COVERAGE
    # both arms are taken at one step, so neither arm's gate differs from the other's
    assert "same step for both arms" in dl.ELIGIBILITY


# ---- no rule-level test, and no cell substitution ------------------------------------------------


def test_the_registration_says_why_a_rule_level_test_is_impossible():
    r = dl.registration()
    assert r["level"] == "link level, not rule level"
    assert "exactly 1 has any compiled rule naming its own gene and 0" in r["why_not_rule_level"]
    assert dl.INHERITED_FIGURES["increases_with_a_rule_gated_on_their_own_cell"] == 0


def test_the_cached_cells_are_named_and_nothing_is_substituted():
    assert dl.CACHED_CELLS == ("HepG2", "IMR-90", "K562", "GM12878")
    assert dl.CELLS_NOT_CACHED == ("WTC11", "HCT116")
    assert dl.PRIMARY_CELL == "K562" and dl.PRIMARY_CELL in dl.CACHED_CELLS
    for cell in dl.CELLS_NOT_CACHED:
        assert cell not in dl.CACHED_CELLS
    assert "NO SUBSTITUTION IS MADE" in dl.CELL_AVAILABILITY
    assert "elements_hct116" in dl.CELL_AVAILABILITY and "NOT read here" in dl.CELL_AVAILABILITY


# ---- the predicates ------------------------------------------------------------------------------


def test_significant_is_repress2s_own_predicate():
    assert dl.significant({"outcome": rp.SIGNIFICANT_OUTCOMES[0]})
    assert dl.significant({"outcome": rp.SIGNIFICANT_OUTCOMES[1]})
    assert not dl.significant({"outcome": "not_significant_well_powered"})


@pytest.mark.parametrize(
    ("value", "expect"), [(-1.0, -1), (1e-9, 1), (-1e-9, -1), (0.0, None), (None, None), (4.21, 1)]
)
def test_sign_of_sends_a_zero_and_an_absent_value_to_no_sign(value, expect):
    assert dl.sign_of(value) == expect


@pytest.mark.parametrize(
    ("down", "up", "arm"),
    [(True, False, dl.DECREASES), (False, True, dl.INCREASES), (True, True, None), (False, False, None)],
)
def test_arm_of_puts_a_contested_link_in_neither_arm(down, up, arm):
    assert dl.arm_of(down, up) is arm


def test_the_measured_sign_follows_from_the_arm_only():
    assert dl.MEASURED_SIGN == {dl.DECREASES: -1, dl.INCREASES: 1}
    assert "No effect size, magnitude or threshold enters it" in dl.MEASURED_SIGN_CALL


# ---- the estimators -----------------------------------------------------------------------------


def _link(arm: str, agrees: bool, chrom: str, start: int, gene: str, cell: str = "K562"):
    sign = dl.MEASURED_SIGN[arm]
    return {
        "arm": arm,
        "agrees": agrees,
        "measured_sign": sign,
        "predicted_sign": sign if agrees else -sign,
        "chrom": chrom,
        "start": start,
        "end": start + 300,
        "gene": gene,
        "cell": cell,
    }


def _population(n_down: int, down_agree: int, n_up: int, up_agree: int):
    out = []
    for i in range(n_down):
        out.append(_link(dl.DECREASES, i < down_agree, f"chr{i % 20 + 1}", 5_000_000 * (i + 1), f"D{i}"))
    for i in range(n_up):
        start = 5_000_000 * (i + 1) + 2_000_000
        out.append(_link(dl.INCREASES, i < up_agree, f"chr{i % 20 + 1}", start, f"U{i}"))
    return out


def test_rate_and_difference_are_per_arm_and_none_when_an_arm_is_empty():
    pop = _population(4, 3, 2, 1)
    by_arm = {a: [r for r in pop if r["arm"] == a] for a in dl.ARMS}
    assert dl.rate(by_arm[dl.DECREASES]) == 0.75
    assert dl.rate(by_arm[dl.INCREASES]) == 0.5
    assert dl.difference(by_arm) == 0.25
    assert dl.rate([]) is None
    assert dl.difference({dl.DECREASES: by_arm[dl.DECREASES], dl.INCREASES: []}) is None


def test_the_gate_reports_both_margins_and_never_pools_the_arms():
    g = dl.gate(_population(10, 10, 0, 0)[:10], dl.DECREASES)
    assert g["answered_links"] == 10
    assert g["meets_both_floors"] is False
    assert g["short_by"] == [
        "answered links 10, short of 30 by 20",
        f"independent loci {g['independent_loci']}, short of 20 by {20 - g['independent_loci']}",
    ]
    assert "not established biological independence" in g["not_biological_independence"]


def test_the_gate_passes_only_at_or_above_both_imported_floors():
    links = [
        _link(dl.DECREASES, True, "chr1", 10_000_000 * (i + 1), f"G{i}") for i in range(dl.POSITIVE_FLOOR)
    ]
    g = dl.gate(links, dl.DECREASES)
    assert g["answered_links"] == dl.POSITIVE_FLOOR == 30
    assert g["independent_loci"] == 30 >= dl.LOCUS_FLOOR
    assert g["meets_both_floors"] is True and g["short_by"] == []


def test_the_bootstrap_unit_is_the_locus_and_links_at_one_locus_are_one_draw():
    # two links 1 kb apart are one locus under cell2's rule, so the bootstrap sees one unit
    close = [
        _link(dl.DECREASES, True, "chr1", 1_000_000, "A"),
        _link(dl.DECREASES, False, "chr1", 1_001_000, "B"),
    ]
    assert len(set(dl.loci_of(close))) == 1
    b = dl.bootstrap(_population(30, 25, 30, 10), n=50, seed=1)
    assert b["unit"] == "independent locus"
    assert b["draws"] == 50 and b["loci"] == len(set(dl.loci_of(_population(30, 25, 30, 10))))
    assert b["ci95"]["difference"][0] <= b["ci95"]["difference"][1]
    assert "over INDEPENDENT LOCI and not over links" in b["over_loci_not_links"]


def test_the_bootstrap_is_deterministic_at_a_seed():
    pop = _population(30, 25, 30, 10)
    assert dl.bootstrap(pop, n=100, seed=7) == dl.bootstrap(pop, n=100, seed=7)


def test_the_bootstrap_skips_a_draw_that_leaves_an_arm_empty_and_reports_how_many():
    pop = _population(30, 20, 1, 0)
    b = dl.bootstrap(pop, n=200, seed=3)
    assert b["skipped_draws"] > 0
    assert b["skipped_draws"] + len(dl.interval(b["ci95"]["difference"]) or []) > 0


def test_the_control_holds_the_arm_sizes_and_the_predicted_signs_fixed():
    pop = _population(40, 36, 40, 8)
    observed = dl.difference({a: [r for r in pop if r["arm"] == a] for a in dl.ARMS})
    c = dl.shuffle_control(pop, observed, n=300, seed=5)
    assert c["shuffles"] == 300
    assert c["observed_difference"] == round(observed, 4)
    assert c["p_one_sided"] is not None and 0.0 <= c["p_one_sided"] <= 1.0
    assert "HELD FIXED" in c["held_fixed"] and "WHAT VARIES" in c["held_fixed"]
    assert "narrower than a null that did" in c["limit"]
    assert "NOT centred on zero" in c["held_fixed"]


def test_the_shuffled_difference_is_centred_on_the_baseline_and_not_on_zero():
    """The trap the control exists to close: a lopsided caller earns a positive difference for free."""
    pop = _population(60, 57, 60, 6)
    d = dl.difference({a: [r for r in pop if r["arm"] == a] for a in dl.ARMS})
    c = dl.shuffle_control(pop, d, n=600, seed=11)
    # 111 of 120 predicted signs are -1 here, so the shuffle alone earns 2q - 1 = 0.85
    assert c["marginal_predicted_down_rate"] == round(111 / 120, 4)
    assert c["baseline"] == round(2 * 111 / 120 - 1, 4)
    assert abs(c["median_permuted_difference"] - c["baseline"]) < 0.1
    assert abs(d - c["baseline"]) < 0.11, "the raw difference here is the baseline, not skill"
    # and read against zero this population would have looked like a large effect
    assert d > 0.8


def test_the_excess_has_a_permutation_expectation_of_zero():
    pop = _population(60, 57, 60, 6)
    assert dl.marginal_down(pop) == 111 / 120
    assert dl.baseline(pop) == pytest.approx(2 * 111 / 120 - 1)
    assert dl.excess(pop) == pytest.approx(
        dl.difference({a: [r for r in pop if r["arm"] == a] for a in dl.ARMS}) - dl.baseline(pop)
    )
    assert dl.excess(pop) == pytest.approx(0.0), "equal arms: the excess is identically zero"


def test_the_registered_collapse_identity_holds_exactly():
    """COLLAPSE in code: excess == 2 * (balanced accuracy - 0.5) * (n_u - n_d) / (n_d + n_u).

    This is the fact that decides which statistic the reading is taken on, so it is pinned rather than
    described. It is checked over arm sizes and agreement counts that differ in every combination.
    """
    import random

    rng = random.Random(20261002)
    for _ in range(200):
        n_d, n_u = rng.randint(2, 40), rng.randint(2, 40)
        a, b = rng.randint(0, n_d), rng.randint(0, n_u)
        pop = _population(n_d, a, n_u, b)
        ba = dl.balanced_accuracy(pop)
        claim = 2 * (ba - dl.CHANCE) * (n_u - n_d) / (n_d + n_u)
        assert dl.excess(pop) == pytest.approx(claim)
        if n_d == n_u:
            assert dl.excess(pop) == pytest.approx(0.0)


def test_a_larger_decrease_arm_makes_a_positive_excess_mean_below_chance():
    """The consequence COLLAPSE registers: with n_u < n_d the difference's one-sided test inverts."""
    good = _population(80, 70, 20, 17)  # reads both signs well
    bad = _population(80, 24, 20, 6)  # reads both signs badly
    assert dl.balanced_accuracy(good) > dl.CHANCE and dl.excess(good) < 0
    assert dl.balanced_accuracy(bad) < dl.CHANCE and dl.excess(bad) > 0


def test_balanced_accuracy_is_the_statistic_a_constant_sign_caller_cannot_win():
    constant_down = [
        _link(dl.DECREASES, True, f"chr{i % 20 + 1}", 7_000_000 * (i + 1), f"D{i}") for i in range(50)
    ] + [
        _link(dl.INCREASES, False, f"chr{i % 20 + 1}", 7_000_000 * (i + 1) + 3_000_000, f"U{i}")
        for i in range(10)
    ]
    assert dl.balanced_accuracy(constant_down) == pytest.approx(dl.CHANCE)
    # and its raw pooled rate is 50/60, which is why no pooled rate is reported as a result
    assert sum(1 for r in constant_down if r["agrees"]) / len(constant_down) == pytest.approx(50 / 60)


def test_the_control_recovers_real_sign_reading_and_rejects_a_constant_sign_caller():
    strong = _population(60, 54, 30, 27)
    d = dl.difference({a: [r for r in strong if r["arm"] == a] for a in dl.ARMS})
    c = dl.shuffle_control(strong, d, n=600, seed=11)
    assert c["p_one_sided"] < dl.CONTROL_ALPHA, "the balanced-accuracy p is the one read"
    assert dl.bootstrap(strong, n=300, seed=4)["ci95"]["balanced_accuracy"][0] > dl.CHANCE
    # a caller at chance on both signs clears neither condition
    flat = _population(60, 30, 30, 15)
    d0 = dl.difference({a: [r for r in flat if r["arm"] == a] for a in dl.ARMS})
    assert dl.shuffle_control(flat, d0, n=600, seed=11)["p_one_sided"] >= dl.CONTROL_ALPHA
    assert dl.bootstrap(flat, n=300, seed=4)["ci95"]["balanced_accuracy"][0] <= dl.CHANCE


def test_the_control_reports_the_differences_own_p_beside_the_one_it_reads():
    pop = _population(60, 54, 30, 27)
    d = dl.difference({a: [r for r in pop if r["arm"] == a] for a in dl.ARMS})
    c = dl.shuffle_control(pop, d, n=300, seed=5)
    assert c["p_one_sided_difference"] is not None
    assert "not the reading" in c["p_one_sided_call"]
    assert c["ci95_permuted_balanced_accuracy"][0] < dl.CHANCE < c["ci95_permuted_balanced_accuracy"][1]


def test_the_control_is_deterministic_at_a_seed():
    pop = _population(30, 25, 30, 10)
    d = dl.difference({a: [r for r in pop if r["arm"] == a] for a in dl.ARMS})
    assert dl.shuffle_control(pop, d, n=200, seed=9) == dl.shuffle_control(pop, d, n=200, seed=9)


# ---- the readings -------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("ci", "p", "established"),
    [
        ([0.6, 0.8], 0.01, True),
        ([0.49, 0.8], 0.01, False),
        ([0.6, 0.8], 0.2, False),
        ([0.5, 0.8], 0.01, False),
        (None, 0.01, False),
        ([0.6, 0.8], None, False),
    ],
)
def test_both_conditions_are_required_for_the_established_reading(ci, p, established):
    """`ci` is the locus interval of balanced accuracy, read against CHANCE and not against zero."""
    r = dl.reading(ci, p)
    assert r["established"] is established
    assert r["reading"] == (dl.ESTABLISHED if established else dl.NOT_DETECTED)
    assert r["there_is_no_third"] is True


def test_there_are_exactly_two_readings_and_no_third():
    assert dl.READINGS == (dl.ESTABLISHED, dl.NOT_DETECTED)
    assert len(dl.READINGS) == 2 and dl.THERE_IS_NO_THIRD is True
    assert set(dl.WHAT_FOLLOWS) == {"established", "not_detected"}


def test_an_undetected_outcome_cannot_read_as_encouraging():
    for word in dl.FORBIDDEN_OF_A_SHORT_OR_UNDETECTED_OUTCOME:
        assert word in dl.NOT_DETECTED
    assert "NEVER as close, promising, nearly enough, a good start or enough for a pilot" in dl.NOT_DETECTED
    assert "not a measured absence" in dl.NOT_DETECTED
    assert "NEVER as close, promising, nearly enough, a good start or enough for a pilot" in dl.GATE_NO_GO
    assert dl.FORBIDDEN_OF_A_SHORT_OR_UNDETECTED_OUTCOME == (
        "close",
        "promising",
        "nearly enough",
        "a good start",
        "enough for a pilot",
    )


def test_the_registration_records_why_the_difference_is_not_the_reading():
    for phrase in (
        "2 * (balanced_accuracy - 0.5) * (n_u - n_d) / (n_d + n_u)",
        "212 against 48",
        "WORSE than chance",
        "exactly ONE skill number",
        "DESCRIPTION of the asymmetry",
    ):
        assert phrase in dl.COLLAPSE
    assert "a constant-sign caller - of either sign - scores exactly 0.5 on it by construction" in (
        dl.BALANCED_CALL
    )


def test_the_prior_direction_lanes_are_carried_verbatim_and_not_claimed_as_this_one():
    p = dl.PRIOR_DIRECTION
    assert p["downward_arm"]["rate"] == 0.8276 and p["upward_arm"]["rate"] == 0.4167
    assert p["balanced_accuracy"]["rate"] == 0.6221
    assert "UNDECIDABLE" in p["registered_verdict_carried_verbatim"]
    assert "neither confirms nor contradicts it" in p["not_this_lanes_population"]


def test_the_established_reading_claims_sign_agreement_and_no_mechanism():
    assert "SIGN AGREEMENT" in dl.ESTABLISHED
    for forbidden in ("silencer", "repressor", "repression mechanism"):
        assert f"not a {forbidden}" in dl.ESTABLISHED or forbidden in dl.ESTABLISHED
    assert "not a verdict on any compiled rule" in dl.ESTABLISHED


def test_what_the_result_cannot_establish_is_its_own_section():
    for phrase in (
        "a sign agreement is not a mechanism",
        "silencer",
        "one cell (K562)",
        "operational grouping",
        "undetected difference is not a measured absence",
        "WTC11 and HCT116 links are absent from the test rather than negative in it",
    ):
        assert phrase in dl.CANNOT_ESTABLISH


# ---- the committed registration -----------------------------------------------------------------


@pytest.mark.skipif(not REGISTRATION.exists(), reason="the registration has not been written")
def test_the_committed_registration_carries_what_the_run_may_not_move():
    r = json.loads(REGISTRATION.read_text())
    assert r["result"] == "direction_link_registration" and r["lane"] == "lane-direction"
    assert r["requests"] == 0
    assert r["floors"]["links"] == 30 and r["floors"]["independent_loci"] == 20
    assert r["floors"]["neither_chosen_here"] is True
    assert r["intervals"]["draws"] == dl.DRAWS and r["intervals"]["seed"] == dl.SEED
    assert r["baseline"] == dl.BASELINE and r["test_statistic"] == dl.TEST_STATISTIC
    assert "NOT on zero" in r["marginal_sign_rate_and_the_baseline"]
    assert r["chance"] == 0.5
    assert "balanced accuracy" in r["test_statistic"]
    assert r["why_the_difference_is_not_the_reading"] == dl.COLLAPSE
    assert r["prior_direction_lanes"]["registered_verdict_carried_verbatim"].startswith(
        "balanced accuracy 0.6221 is above chance but under 0.65"
    )
    assert r["control"]["shuffles"] == dl.SHUFFLES and r["control"]["seed"] == dl.SHUFFLE_SEED
    assert r["control"]["alpha"] == dl.CONTROL_ALPHA
    assert r["cells"]["cached_cell_tracks"] == list(dl.CACHED_CELLS)
    assert r["cells"]["cells_not_cached"] == list(dl.CELLS_NOT_CACHED)
    assert r["readings_before_the_outcome_was_seen"]["there_is_no_third"] is True
    assert r["readings_before_the_outcome_was_seen"]["established"] == dl.ESTABLISHED
    assert r["readings_before_the_outcome_was_seen"]["not_detected"] == dl.NOT_DETECTED
    assert r["follows_from"]["registered_go_carried_verbatim"] == inc.GO
    assert r["follows_from"]["and_behind_it"]["registered_no_go_carried_verbatim"] == inc.INHERITED_NO_GO
    assert r["what_the_run_will_write"] == "data/results/direction_link.json"


@pytest.mark.skipif(not REGISTRATION.exists(), reason="the registration has not been written")
def test_the_registration_declares_the_cache_group_with_every_member_named():
    r = json.loads(REGISTRATION.read_text())
    groups = [e for e in r["result_manifest"]["inputs"] if e.get("path") == "per_element_response_cache"]
    assert len(groups) == 1
    assert groups[0]["group"] is True
    assert groups[0]["files"] == len(groups[0]["members"]) >= 24
    assert all(m["sha256"] and m["bytes"] for m in groups[0]["members"])


# ---- the streaming reader of the per-element cache ------------------------------------------------


def _archive(path, n_elements=120, seed=7):
    """A synthetic archive shaped like the sweep's own, with braces and escapes inside strings."""
    import gzip
    import random

    rng = random.Random(seed)
    rec = {}
    for i in range(n_elements):
        eid = f"EH38E{1000000 + i}"
        genes = []
        for j in range(rng.randint(1, 40)):
            genes.append(
                {
                    "gene": rng.choice([f"G{j}", "odd{brace}name", 'quote\\"inside', f"LINC{j:05d}"]),
                    "n_tracks": 371,
                    "mean_log2fc": round(rng.uniform(-1, 1), 4),
                    "max_drop_tissue": "a } tissue { with braces",
                    "by_cell": {
                        "HepG2": round(rng.uniform(-1, 1), 4),
                        "IMR-90": 0.0,
                        "K562": round(rng.uniform(-1, 1), 4),
                        "GM12878": round(rng.uniform(-1, 1), 4),
                    },
                }
            )
        rec[eid] = {"id": eid, "chrom": "chrT", "start": i * 1000, "end": i * 1000 + 300, "genes": genes}
    with gzip.open(path, "wt") as fh:
        json.dump(rec, fh)
    return rec


@pytest.mark.parametrize("chunk", [1 << 10, 1 << 14, 1 << 23])
def test_the_cache_reader_returns_exactly_the_wanted_records_at_any_chunk_size(tmp_path, chunk):
    """The reader streams, so a record may straddle any number of chunk boundaries."""
    from scripts.direction_link import cached_records

    path = tmp_path / "chrT.json.gz"
    rec = _archive(path)
    import random

    wanted = set(random.Random(3).sample(sorted(rec), 25)) | {"EH38E9999999"}
    got = cached_records(path, wanted, chunk=chunk)
    assert set(got) == wanted - {"EH38E9999999"}
    assert all(got[k] == rec[k] for k in got)


def test_the_cache_reader_asks_for_nothing_when_nothing_is_wanted(tmp_path):
    from scripts.direction_link import cached_records

    path = tmp_path / "chrT.json.gz"
    _archive(path, n_elements=5)
    assert cached_records(path, set()) == {}
    assert cached_records(tmp_path / "absent.json.gz", {"EH38E1000000"}) == {}


# ---- the committed result -------------------------------------------------------------------------

RESULT = Path("data/results/direction_link.json")


@pytest.mark.skipif(not RESULT.exists(), reason="the test has not been run")
def test_the_committed_result_reconciles_with_the_population_it_follows_from():
    r = json.loads(RESULT.read_text())
    e = r["eligibility"]
    assert e["links"] == 260
    assert e["by_arm"] == {"decreases": 212, "increases": 48}
    assert e["contested"] == dl.INHERITED_FIGURES["contested_by_a_regulated_pair"] == 0
    # the two arms are lane-increase's own two counts at the same ladder step
    assert e["by_arm"]["increases"] == dl.INHERITED_FIGURES["significant_increases_on_an_attributed_element"]
    assert e["by_arm"]["decreases"] == dl.INHERITED_FIGURES["significant_decreases_on_an_attributed_element"]
    assert e["by_arm_and_cell"]["increases"] == dl.INHERITED_FIGURES["increase_links_by_cell"]


@pytest.mark.skipif(not RESULT.exists(), reason="the test has not been run")
def test_the_answerability_breakdown_is_exhaustive_at_its_own_denominator():
    a = json.loads(RESULT.read_text())["answerability"]
    assert a["cell_read"] == "K562"
    assert a["breakdown_reconciles"] is True
    assert a["answered"] + a["predicted_zero_excluded"] + a["absent"] == a["links_in_the_cell_read"]
    assert a["answered"] == 195 and a["absent"] == 35 and a["predicted_zero_excluded"] == 0
    for arm, n in (("decreases", 191), ("increases", 39)):
        b = a["by_arm"][arm]
        assert b["in_the_cell_read"] == n
        assert b["answered"] + b["predicted_zero_excluded"] + b["absent"] == n


@pytest.mark.skipif(not RESULT.exists(), reason="the test has not been run")
def test_no_cell_the_cache_does_not_carry_entered_any_denominator():
    a = json.loads(RESULT.read_text())["answerability"]
    out = a["unanswerable_because_the_cache_carries_no_such_cell"]
    assert out["links"] == 26
    # the run met a third uncached cell on the decrease arm that the increase population does not hold
    assert set(out["by_arm_and_cell"]["decreases"]) == {"HCT116", "Jurkat", "WTC11"}
    assert set(out["by_arm_and_cell"]["increases"]) == {"HCT116", "WTC11"}
    for arm in dl.ARMS:
        for cell in out["by_arm_and_cell"][arm]:
            assert cell not in dl.CACHED_CELLS
    assert "nothing_substituted" in out
    # and a cached cell other than the one registered as primary is counted, not pooled in
    assert a["links_in_another_cached_cell_not_read"]["by_arm_and_cell"]["decreases"] == {"GM12878": 4}


@pytest.mark.skipif(not RESULT.exists(), reason="the test has not been run")
def test_the_gate_refused_the_comparison_and_nothing_was_computed_past_it():
    r = json.loads(RESULT.read_text())
    g = r["gate"]
    assert g["both_arms_meet_both_floors"] is False
    assert g["reading"] == dl.GATE_NO_GO
    assert g["per_arm"]["decreases"]["meets_both_floors"] is True
    assert g["per_arm"]["increases"]["short_by"] == [
        "answered links 21, short of 30 by 9",
        "independent loci 13, short of 20 by 7",
    ]
    c = r["comparison"]
    assert c["taken"] is False
    assert c["why"] == dl.GATE_NO_GO
    for absent in ("balanced_accuracy", "control", "bootstrap", "difference_between_the_arms"):
        assert absent not in c, "the gate refused the comparison, so none of it may be in the result"
    note = c["the_rates_in_arms_are_on_the_record_not_read_as_a_result"]
    assert "they are NOT read as a result" in note
    assert "no sign shuffle was run" in note
    assert c["statistics_not_computed"]


@pytest.mark.skipif(not RESULT.exists(), reason="the test has not been run")
def test_the_result_carries_the_floors_the_locus_rule_and_what_it_cannot_establish():
    r = json.loads(RESULT.read_text())
    for arm in dl.ARMS:
        assert r["gate"]["per_arm"][arm]["floors"] == {"links": 30, "independent_loci": 20}
        assert (
            "not established biological independence"
            in (r["gate"]["per_arm"][arm]["not_biological_independence"])
        )
    assert r["requests"] == 0
    assert r["cannot_establish"] == dl.CANNOT_ESTABLISH
    assert r["locus_convention"]["span"] == 1_000_000
    assert r["prior_direction_lanes"]["registered_verdict_carried_verbatim"].startswith(
        "balanced accuracy 0.6221"
    )
