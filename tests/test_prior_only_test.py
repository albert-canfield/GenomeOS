# SPDX-License-Identifier: AGPL-3.0-or-later
"""The pilot's prior-only lead, tested on its own (lane-prior, 2026-09-29): the admissibility record.
Synthetic rows only: no benchmark file, no cache, no request."""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("pot", ROOT / "scripts" / "prior_only_test.py")
pot = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pot)

HARNESS_STUDIES = (
    "Gasperini2019",
    "HCT116",
    "K562_DC_TAP",
    "Klann",
    "Morris",
    "Nasser2021",
    "Reilly",
    "Schraivogel2020",
    "WTC11_DC_TAP",
    "Xie",
)


# --- negatives first -------------------------------------------------------------------------------


def test_every_held_out_file_study_is_inadmissible():
    # the benchmark kept a held-out positive only with H3K27ac at the element; the prior reads H3K27ac
    for src, v in pot.STUDIES.items():
        if v["file"] == pot.HELDOUT_FILE or pot.HELDOUT_FILE in v["file"]:
            assert v["class"] == pot.INADMISSIBLE, src


def test_the_other_half_of_the_lead_is_not_admissible():
    assert pot.STUDIES["crispri:Xie"]["class"] == pot.INADMISSIBLE
    assert pot.STUDIES["crispri:Morris"]["class"] == pot.INADMISSIBLE
    assert "crispri:Xie" not in pot.ADMISSIBLE_SOURCES


def test_no_admissible_study_is_unqualified():
    # both admissible screens chose their tested candidates on chromatin: the selection is stated
    assert set(pot.ADMISSIBLE_SOURCES) == {"crispri:Gasperini2019", "crispri:Schraivogel2020"}
    for src in pot.ADMISSIBLE_SOURCES:
        assert pot.STUDIES[src]["class"] == pot.ADMISSIBLE_STATED


# --- the record ------------------------------------------------------------------------------------


def test_every_harness_study_has_a_record_with_a_quoted_source():
    assert set(pot.STUDIES) == {f"crispri:{s}" for s in HARNESS_STUDIES}
    for src, v in pot.STUDIES.items():
        assert v["class"] in (pot.INADMISSIBLE, pot.ADMISSIBLE_STATED, pot.ADMISSIBLE), src
        assert v["quotes"], src
        for q in v["quotes"]:
            assert q["text"] and q["url"].startswith("https://"), src


def test_the_benchmark_evidence_names_the_filter_and_its_reproduction():
    held = pot.BENCHMARK_EVIDENCE["held_out_file"]
    assert any("no H3K27ac categories" in q["text"] for q in held["quotes"])
    assert "4,378 of 4,378" in held["reproduced_by_this_lane"]
    assert "no chromatin filter" in pot.BENCHMARK_EVIDENCE["training_file"]["filters"]


def test_local_counts_count_positives_outside_the_h3k27ac_categories():
    rows = {
        "crispri:A": [
            {"_file": "f", "Regulated": "TRUE", "elementChromatinCategory": "High H3K27ac"},
            {"_file": "f", "Regulated": "TRUE", "elementChromatinCategory": "CTCF element"},
            {"_file": "f", "Regulated": "FALSE", "elementChromatinCategory": "No H3K27ac"},
            {"_file": "f", "Regulated": "FALSE", "elementChromatinCategory": "H3K27ac"},
        ]
    }
    c = pot.local_counts(rows)["crispri:A"]["f"]
    assert (c["pairs"], c["positives"]) == (4, 2)
    assert c["positives_outside_h3k27ac_categories"] == 1
    assert c["pairs_outside_h3k27ac_categories"] == 2


def test_partition_counts_apply_the_harness_floors():
    ho = pot.ho
    ms = pot.ms
    units = []
    for i in range(30):
        outcome = ms.DECREASE if i < 21 else ms.NULL_INFORMATIVE
        units.append(ho.Unit("crispri:A", "chr2", i * 2_000_000, i * 2_000_000 + 500, outcome, gene=f"G{i}"))
    units.append(ho.Unit("crispri:A", "chr1", 0, 500, ms.DECREASE, gene="X"))
    got = pot.partition_counts({"crispri:A": tuple(units)}, ("chr2",))["crispri:A"]
    assert (got["units"], got["positives"], got["negatives"]) == (30, 21, 9)
    assert got["clears_floors"] is False  # 9 negatives, below the floor of 20


# --- the registration ------------------------------------------------------------------------------


def _cmp(src, a, b, interval, diff=0.1, av=0.5, bv=0.4):
    return {
        "source": src,
        "endpoint": "decrease",
        "a": a,
        "b": b,
        "status": "scored",
        "interval": interval,
        "difference": diff,
        "a_value": av,
        "b_value": bv,
    }


def _scores(src, share=0.97):
    return [
        {"source": src, "endpoint": "decrease", "labels": n, "status": "scored", "coverage": {"share": share}}
        for n in pot.LABELLINGS
    ]


def _comps(src, lead, activity, target, machinery, act_dist):
    return [
        _cmp(src, pot.PRIOR, pot.DISTANCE, lead),
        _cmp(src, pot.PRIOR_ACTIVITY, pot.PRIOR_DISTANCE, activity),
        _cmp(src, pot.PRIOR, pot.PRIOR_ACTIVITY, target),
        _cmp(src, pot.PRIOR_DISTANCE, pot.DISTANCE, machinery),
        _cmp(src, pot.PRIOR_ACTIVITY, pot.DISTANCE, act_dist),
        _cmp(src, pot.PRIOR, pot.UNCHANGED, [0.1, 0.2]),
    ]


SRC = "crispri:Gasperini2019"
UP, ZERO, DOWN = [0.05, 0.2], [-0.05, 0.2], [-0.2, -0.05]


def test_a_failed_reproduction_voids_the_run_whatever_the_scores():
    v = pot.assess(_comps(SRC, UP, UP, ZERO, ZERO, UP), _scores(SRC), reproduced=False)
    assert v["reading"] == "void_reproduction" and not v["replicated"]


def test_an_interval_across_zero_does_not_replicate():
    v = pot.assess(_comps(SRC, ZERO, UP, ZERO, ZERO, UP), _scores(SRC), reproduced=True)
    assert v["reading"] == "not_replicated"


def test_an_interval_below_zero_is_reversed():
    v = pot.assess(_comps(SRC, DOWN, UP, ZERO, ZERO, UP), _scores(SRC), reproduced=True)
    assert v["reading"] == "reversed"


def test_low_coverage_does_not_pass():
    v = pot.assess(_comps(SRC, UP, UP, ZERO, ZERO, UP), _scores(SRC, share=0.79), reproduced=True)
    assert v["reading"] == "not_replicated" and not v["endpoints"][0]["passes"]


def test_a_gain_activity_does_not_carry_is_not_an_activity_reading():
    # (i): activity on identical pairs includes zero; (iii): the machinery carries it
    v = pot.assess(_comps(SRC, UP, ZERO, ZERO, UP, ZERO), _scores(SRC), reproduced=True)
    assert v["reading"] == "replicated_not_activity"
    f = v["endpoints"][0]["falsifier"]
    assert f["i_activity_does_not_carry"] and f["iii_machinery_carries"]


def test_a_gain_the_compiled_target_carries_is_not_an_activity_reading():
    v = pot.assess(_comps(SRC, UP, UP, UP, ZERO, ZERO), _scores(SRC), reproduced=True)
    assert v["reading"] == "replicated_not_activity"
    assert v["endpoints"][0]["falsifier"]["ii_compiled_target_carries"]


def test_the_activity_reading_needs_the_lead_and_activity_above_zero():
    v = pot.assess(_comps(SRC, UP, UP, ZERO, ZERO, UP), _scores(SRC), reproduced=True)
    assert v["reading"] == "replicated_activity" and v["replicated"]


def test_no_scored_primary_endpoint_is_no_result():
    v = pot.assess([], [], reproduced=True)
    assert v["reading"] == "no_admissible_endpoint"


def test_reproduction_is_to_the_fourth_decimal():
    lead = pot.LEAD
    ok = _cmp(SRC, pot.PRIOR, pot.DISTANCE, [0.13354, 0.24909], 0.19271, 0.66012, 0.46738)
    assert pot.reproduces(ok, lead)
    off = _cmp(SRC, pot.PRIOR, pot.DISTANCE, [0.1336, 0.2491], 0.1927, 0.6601, 0.4674)
    assert not pot.reproduces(off, lead)


def test_the_registration_states_what_a_pass_would_not_mean():
    r = pot.registration()
    nots = " ".join(r["meaning"]["a_pass_would_not_mean"])
    assert "not a new finding" in r["readings"]["replicated_activity"].lower() or "new finding" in nots
    assert "ABC" in nots and "validation" in nots
    assert "never ABC" in r["abc_form"]
    assert r["primary_endpoints"] == [["crispri:Gasperini2019", "decrease"]]
    assert set(pot.FRESH).isdisjoint(pot.SEEN) and len(pot.FRESH) == 17


# --- the three prior labellings, on synthetic blocks, through the pilot's own code ------------------


def _fixed():
    pb, pl = pot.pb, pot.pl
    blocks = pb.Blocks("chrS", ["B1", "B2"], [20_000, 60_000], [20_500, 60_500], 500)
    peaks = pl.Peaks({"K562": [(60_000, 60_500)]})
    return pb.Fixed(
        {"chrS": blocks}, {"chrS": peaks}, {"chrS": {"G": 10_000}}, {"chrS": {"B1": ("G", 0.8, "K562")}}
    )


def _unit(s, e):
    return pot.ho.Unit("crispri:T", "chrS", s, e, pot.ms.DECREASE, gene="G", tss=10_000, cell="K562")


def _sig(x):
    import math

    return 1.0 / (1.0 + math.exp(-x))


def test_each_variant_removes_one_term_and_no_weight():
    import math

    fx = _fixed()
    u1, u2, out = _unit(20_100, 20_300), _unit(60_100, 60_300), _unit(90_000, 90_200)
    held = [u1, u2, out]
    got = {}
    for name in (pot.PRIOR_DISTANCE, pot.PRIOR_ACTIVITY, pot.PRIOR):
        lab, _ = pot.prior_labels(name, fx, held, "crispri:T", ("chrS",))
        got[name] = [lab.predict(u, "decrease") for u in held]
    # outside every cCRE all three abstain: the same coverage
    assert all(got[n][2] is None for n in got)
    d1, d2 = 10_250, 50_250
    t = lambda d: 3.0 - math.log1p(d / 5000.0)  # noqa: E731
    # (a) distance alone: activity constant at -1
    assert abs(got[pot.PRIOR_DISTANCE][0] - _sig(t(d1)) * _sig(-1.0)) < 1e-9
    assert got[pot.PRIOR_DISTANCE][0] > got[pot.PRIOR_DISTANCE][1]
    # (b) adds H3K27ac: B2 carries a K562 peak (+2), B1 not (-1)
    assert abs(got[pot.PRIOR_ACTIVITY][1] - _sig(t(d2)) * _sig(1.0)) < 1e-9
    assert abs(got[pot.PRIOR_ACTIVITY][0] - _sig(t(d1)) * _sig(-2.0)) < 1e-9
    # (c) adds the compiled target on B1: +1 + min(0.8, 1) to t, +0.5 to a (model cell K562)
    assert abs(got[pot.PRIOR][0] - _sig(t(d1) + 1.8) * _sig(-1.5)) < 1e-9
    assert got[pot.PRIOR][1] == got[pot.PRIOR_ACTIVITY][1]
    assert pot.pl.W["no_target"] == -3.0 and pot.pl.W["distance_d0"] == 5000.0
