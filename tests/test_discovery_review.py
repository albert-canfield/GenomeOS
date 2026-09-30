# SPDX-License-Identifier: AGPL-3.0-or-later
"""A's adoption review, discovery only (lane-discover, 2026-09-29): union coverage, one row per changed
observation, unchanged denominators, the newly ambiguous rows, the stated-cell record and the shifted
screen, on synthetic intervals. The run reads the benchmark and the compiled programs; these tests do not."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

from genomeos.attribution import correctness as cr
from genomeos.attribution import measured as ms

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location("discovery_review", ROOT / "scripts/discovery_review.py")
dr = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = dr
_SPEC.loader.exec_module(dr)
pc = dr.pc
pa = dr.pa


def _pair(start, end, gene, cell="K562", effect=-0.3, significant=True, power=0.9, ds="D"):
    return ms.CrispriPair(
        chrom="chr1",
        start=start,
        end=end,
        gene=gene,
        cell=cell,
        dataset=ds,
        reference="R",
        regulated=significant and effect < 0,
        significant=significant,
        effect_size=effect,
        p_adjusted=0.01 if significant else 0.9,
        power_at_effect_size_20=power,
    )


# E0..E5: compiled predicted links on chr1; R0, R1: registry elements no link carries
ELEMENTS = [
    (100, 300, "G1", "activates"),  # E0: inside P0's 1000 bp interval, policy only
    (10100, 10300, "G2", "activates"),  # E1: inside P1, policy only
    (20100, 20300, "G3", "activates"),  # E2: inside P2, policy only
    (20500, 20700, "G3", "activates"),  # E3: inside P2, policy only
    (30000, 30300, "G4", "activates"),  # E4: P3 under both rules
    (40000, 40300, "G5", "activates"),  # E5: P4 under both rules
]
REGISTRY_ONLY = [
    (20200, 20400),  # R0: overlaps E2 inside P2, policy only
    (30300, 30400),  # R1: beside E4 inside P3, policy only
]
PAIRS = [
    _pair(0, 1000, "G1"),  # P0: newly attached, unique
    _pair(10000, 11000, "G2"),  # P1: newly attached, unique
    _pair(20000, 21000, "G3"),  # P2: newly attached, ambiguous over two links and R0
    _pair(30000, 30400, "G4"),  # P3: unique under the current rule, ambiguous under the policy
    _pair(40000, 40400, "G5", cell="GM12878"),  # P4: unique under both, unchanged
    _pair(50000, 50500, "G6"),  # P5: touches nothing
]
STATED = {
    ("E0", "G1"): frozenset({"K562"}),
    ("E1", "G2"): frozenset({"HepG2"}),
    ("E2", "G3"): frozenset({"HepG2"}),
    ("E3", "G3"): frozenset({"K-562"}),  # the same cell, written differently: norm_cell matches it
    ("E4", "G4"): frozenset({"K562"}),
    ("E5", "G5"): frozenset({"K562"}),
}


def _world(elements=ELEMENTS, registry_only=REGISTRY_ONLY):
    links = [pc.PredictedLink(i, "chr1", s, e, f"E{i}", g, a) for i, (s, e, g, a) in enumerate(elements)]
    reg = [(x.chrom, x.start, x.end, x.element) for x in links]
    reg += [("chr1", s, e, f"R{i}") for i, (s, e) in enumerate(registry_only)]
    return links, pa.ElementSet((x.chrom, x.start, x.end, x.index) for x in links), pa.ElementSet(reg)


def _both(pairs=PAIRS, world=None):
    links, link_set, registry = world or _world()
    cur = pc.observe(pairs, links, link_set, registry, pc.RULES[dr.CURRENT])
    pol = pc.observe(pairs, links, link_set, registry, pc.RULES[dr.POLICY])
    return cur, pol, links, registry


# --- union coverage: overlapping candidates are counted once ---------------------------------------
def test_union_coverage_counts_a_shared_base_once():
    # two candidates overlapping by 100 bp inside a 1000 bp interval: 500 bp covered, not 600
    assert dr.union_coverage(0, 1000, [(100, 400), (300, 600)]) == pytest.approx(0.5)
    # the same candidate twice is one candidate's worth
    assert dr.union_coverage(0, 1000, [(100, 300), (100, 300)]) == pytest.approx(0.2)
    # abutting intervals merge; a candidate past the interval's end is clipped to it
    assert dr.union_coverage(0, 1000, [(0, 100), (100, 200), (900, 1300)]) == pytest.approx(0.3)
    # nested candidates: the outer one's length
    assert dr.union_coverage(0, 1000, [(100, 700), (200, 300)]) == pytest.approx(0.6)
    assert dr.union_coverage(0, 1000, []) == 0.0
    assert dr.union_coverage(10, 10, [(0, 100)]) is None


def test_the_fraction_of_each_candidate_covered():
    assert dr.candidate_fraction(0, 1000, 900, 1100) == pytest.approx(0.5)
    assert dr.candidate_fraction(0, 1000, 100, 300) == pytest.approx(1.0)
    assert dr.candidate_fraction(0, 1000, 5, 5) is None


def test_a_changed_row_carries_every_admitted_candidate_and_the_union_coverage():
    cur, pol, links, registry = _both()
    rows = dr.changed_rows(cur, pol, registry, links, STATED)
    by = {r["id"]: r for r in rows}
    p2 = by[pc.observation_id(PAIRS[2])]
    # E2 [20100, 20300), R0 [20200, 20400) and E3 [20500, 20700) in [20000, 21000): 300 + 200 = 500 bp
    assert [x["id"] for x in p2["candidates"]] == ["E2", "E3", "R0"]
    assert p2["multiplicity"] == {dr.CURRENT: 0, dr.POLICY: 3}
    assert p2["union_coverage"] == {dr.CURRENT: 0.0, dr.POLICY: 0.5}
    assert all(x["fraction_covered"] == 1.0 for x in p2["candidates"])
    assert p2["change"] == dr.NEWLY_ATTACHED
    assert p2["status_under_new_rule"] == dr.UNRESOLVED
    assert p2["verdict_current_rule"] == pc.UNATTACHED
    assert p2["direction_description_policy"] == "agrees"
    # P3: E4 under both rules, R1 under the policy only; 300 + 100 = the whole 400 bp
    p3 = by[pc.observation_id(PAIRS[3])]
    assert p3["union_coverage"] == {dr.CURRENT: 0.75, dr.POLICY: 1.0}
    assert [(x["id"], x["current"]) for x in p3["candidates"]] == [("E4", True), ("R1", False)]


# --- observation-once and the unchanged denominators -----------------------------------------------
def test_each_changed_observation_appears_once_and_unchanged_ones_not_at_all():
    cur, pol, links, registry = _both()
    rows = dr.changed_rows(cur, pol, registry, links, STATED)
    ids = [r["id"] for r in rows]
    assert ids == [pc.observation_id(p) for p in PAIRS[:4]]
    assert len(set(ids)) == len(ids)
    # the same observation twice is refused, never counted twice
    with pytest.raises(ValueError):
        dr.changed_rows([*cur, cur[0]], [*pol, pol[0]], registry, links, STATED)
    with pytest.raises(ValueError):
        pc.observe([*PAIRS, PAIRS[0]], links, *_world()[1:], pc.RULES[dr.POLICY])


def test_the_denominators_are_unchanged_and_nothing_is_lost():
    cur, pol, links, registry = _both()
    rows = dr.changed_rows(cur, pol, registry, links, STATED)
    conds = dr.conditions(cur, pol, rows, len(links), len(links))
    one = conds["1_nothing_lost_denominators_unchanged"]
    assert one["current_attachments_lost"] == {"candidate_sets": 0, "attached_links": 0}
    assert one["denominators"]["observations"] == {
        dr.CURRENT: len(PAIRS),
        dr.POLICY: len(PAIRS),
        "same_observation_ids": True,
    }
    assert conds["2_each_observation_once_with_its_candidates_and_coverage"]["each_once"]
    assert conds["hold"]
    # a moved link denominator, or a rule that loses what the current rule attaches, fails condition 1
    assert not dr.conditions(cur, pol, rows, len(links), len(links) + 1)["hold"]
    swapped = dr.conditions(pol, cur, rows, len(links), len(links))
    assert swapped["1_nothing_lost_denominators_unchanged"]["current_attachments_lost"]["attached_links"] > 0
    assert not swapped["hold"]


def test_change_classes():
    cur, pol, links, registry = _both()
    got = [dr.change_class(c, p) for c, p in zip(cur, pol, strict=True)]
    assert got == [
        dr.NEWLY_ATTACHED,
        dr.NEWLY_ATTACHED,
        dr.NEWLY_ATTACHED,
        dr.CANDIDATES_ADDED_ATTACHED,
        None,
        None,
    ]
    assert dr.change_class(pol[3], cur[3]) == dr.LOST
    with pytest.raises(ValueError):
        dr.change_class(cur[0], pol[1])


# --- newly ambiguous: the historical verdict beside 'unresolved' -----------------------------------
def test_a_unique_assignment_that_becomes_ambiguous_keeps_its_historical_verdict():
    cur, pol, links, registry = _both()
    assert cur[3].resolution == pc.UNIQUE and cur[3].verdict == pc.SUPPORTED
    assert pol[3].resolution == pc.AMBIGUOUS_ONE
    rows = dr.changed_rows(cur, pol, registry, links, STATED)
    na = dr.newly_ambiguous(rows)
    assert na["n"] == 1
    (e,) = na["each"]
    assert e["id"] == pc.observation_id(PAIRS[3])
    assert e["historical_verdict_current_rule"] == pc.SUPPORTED
    assert e["status_under_new_rule"] == dr.UNRESOLVED
    assert na["historical_verdicts"] == {pc.SUPPORTED: 1}


# --- the stated-cell record: the 24/13 split, the 39 of 79 and the 39/54 check -----------------------
def test_stated_cells_come_from_the_target_claims_of_the_compiled_rules():
    claims = [
        cr.Claim("E0", "chr1", 100, 300, cr.TARGET, "G1", "G1", "K562"),
        cr.Claim("E0", "chr1", 100, 300, cr.TARGET, "G1", "G1", "HepG2"),
        cr.Claim("E0", "chr1", 100, 300, cr.ACTIVITY, "activates_target", "G1", "liver"),
        cr.Claim("E1", "chr1", 10100, 10300, cr.TARGET, "G2", "G2", "unknown"),
    ]
    stated, ruled = dr.stated_cells(claims)
    assert stated == {("E0", "G1"): frozenset({"K562", "HepG2"})}
    assert ruled == {("E0", "G1"), ("E1", "G2")}


def test_stated_cell_agreement():
    cur, pol, links, _ = _both()
    assert dr.agreement(pol[0], links, STATED) == dr.MATCH
    assert dr.agreement(pol[1], links, STATED) == dr.OTHER_ONLY
    assert dr.agreement(pol[2], links, STATED) == dr.MATCH  # E3's 'K-562' among the attached links
    assert dr.agreement(pol[4], links, STATED) == dr.OTHER_ONLY  # a GM12878 screen of a K562 claim
    assert dr.agreement(pol[0], links, {}) == dr.NO_STATED
    assert dr.agreement(pol[5], links, STATED) == dr.NOT_ATTACHED


def test_the_stated_cell_record_is_recomputed_and_checked():
    cur, pol, links, _ = _both()
    rec = dr.stated_cell_record(cur, pol, links, STATED, {"in_context": 1, "elsewhere_only": 1})
    got = rec["recomputed"]
    assert got["newly_supported_unique"]["n"] == 2
    assert got["newly_supported_unique"][dr.MATCH] == 1
    assert got["newly_supported_unique"][dr.OTHER_ONLY] == 1
    assert got["newly_supported_unique"]["cells_of_the_matches"] == {"K562": 1}
    assert got["newly_supported_ambiguous"] == {
        "n": 1,
        "with_a_stated_cell_match": 1,
        dr.MATCH: 1,
        dr.OTHER_ONLY: 0,
        dr.NO_STATED: 0,
    }
    assert got["current_rule_supported"]["n"] == 2
    assert got["current_rule_supported"][dr.MATCH] == 1
    assert got["current_rule_supported"]["not_a_match"] == 1
    assert rec["check_against_the_v2_v3_direction_split"]["equal"]
    assert not dr.stated_cell_record(cur, pol, links, STATED, {"in_context": 2, "elsewhere_only": 0})[
        "check_against_the_v2_v3_direction_split"
    ]["equal"]
    # the synthetic world is not the benchmark: the provisional 24/13 reading is not claimed here
    assert not rec["reproduces_the_provisional_reading"]
    assert rec["provisional_in_memory_reading"] == dr.PROVISIONAL


def test_the_provisional_reading_is_recorded_as_given():
    p = dr.PROVISIONAL
    assert p["newly_supported_unique"] == {"n": 37, dr.MATCH: 24, dr.OTHER_ONLY: 13, "matches_all_in": "K562"}
    assert p["newly_supported_ambiguous"] == {"n": 79, "with_a_stated_cell_match": 39}
    assert p["current_rule_supported"] == {"n": 93, dr.MATCH: 39, "not_a_match": 54}


# --- the shifted screen ----------------------------------------------------------------------------
@pytest.mark.parametrize(
    "real, shifted, want",
    [
        (3, [1, 1], dr.PASS),
        (3, [2, 2], dr.PASS),  # exactly 1.5
        (1, [1, 1], dr.FAIL),
        (2, [0, 0], dr.UNDEFINED),
        (0, [1, 0], dr.NOTHING),
        (0, [0, 0], dr.NOTHING),
    ],
)
def test_enrichment_reading(real, shifted, want):
    assert dr.enrichment(real, shifted, 100)["screen"] == want


def test_the_screen_uses_every_observation_of_its_stratum_and_the_mean_of_both_shifts():
    cur, pol, links, registry = _both()
    real = dr.flags(cur, pol)
    moved = []
    for shift in dr.SHIFTS:
        sp = dr.shifted(PAIRS, shift)
        assert all(q.end - q.start == p.end - p.start for p, q in zip(PAIRS, sp, strict=True))
        moved.append(
            dr.flags(*(pc.observe(sp, links, *_world()[1:], pc.RULES[n]) for n in (dr.CURRENT, dr.POLICY)))
        )
    scr = dr.screen(PAIRS, real, moved)
    rows = scr[dr.SCREEN_QUANTITY]
    assert rows["all_observations"]["observations"] == len(PAIRS)
    assert rows["all_observations"]["real"] == 3
    assert rows["D|K562"]["observations"] == 5 and rows["D|GM12878"]["observations"] == 1
    # moved 20 kb, P0 lands on E2 and E3, which link G3, not G1: nothing attaches to G1 there
    assert rows["all_observations"]["shifted"] == [0, 0]
    assert rows["all_observations"]["screen"] == dr.UNDEFINED
    for q in dr.ALSO_SCREENED:
        assert "screen" not in scr[q]["all_observations"]


def test_the_decision_reading():
    ok = {"hold": True}

    def scr(pooled, strata):
        rows = {"all_observations": {"screen": pooled}, **{k: {"screen": v} for k, v in strata.items()}}
        return {dr.SCREEN_QUANTITY: rows}

    assert dr.reading(ok, scr(dr.PASS, {"a": dr.PASS, "b": dr.UNDEFINED}))["reading"] == dr.MERITS
    r = dr.reading(ok, scr(dr.PASS, {"a": dr.PASS, "b": dr.FAIL}))
    assert r["reading"] == dr.RESTRICTED and r["strata_failing"] == ["b"]
    assert dr.reading(ok, scr(dr.FAIL, {"a": dr.PASS}))["reading"] == dr.DOES_NOT
    assert dr.reading({"hold": False}, scr(dr.PASS, {"a": dr.PASS}))["reading"] == dr.DOES_NOT


# --- summaries and the registration ----------------------------------------------------------------
def test_summary_counts_tenths_with_exactly_one_apart():
    s = dr.summary([0.0, 0.05, 0.5, 0.95, 1.0, 1.0])
    assert s["n"] == 6
    assert s["tenths"] == {"0.0-0.1": 2, "0.5-0.6": 1, "0.9-1.0": 1, "1.0": 2}
    assert dr.summary([]) == {"n": 0}
    with pytest.raises(ValueError):
        dr.summary([1.2])


def test_registration_uses_no_width_boundary_and_the_registered_wording():
    text = json.dumps(dr.REGISTRATION)
    assert "700" not in text
    assert "possible causal" not in text
    assert "all admitted registry candidates" in text
    assert dr.SCREEN_THRESHOLD == 1.5
    assert dr.SHIFTS == (20_000, -20_000)
    assert "not an estimate of false attribution" in dr.REGISTRATION["screen"]["what"]
    assert (
        "distance" in dr.REGISTRATION["screen"]["caveat"] and "density" in dr.REGISTRATION["screen"]["caveat"]
    )
    commits = {e["commit"] for e in dr.PRIOR_EXPOSURE}
    assert {"b7e4bf0", "dd8c49c", "6c0d39f", "0183ef5", "1253d69", "32af2cb"} <= commits
    assert any("in-memory" in e["read"] for e in dr.PRIOR_EXPOSURE)
