# SPDX-License-Identifier: AGPL-3.0-or-later
"""Audit A (lane-place, 2026-09-29): the placement cascade, the registered policy and its statistics, on
synthetic intervals. The run itself reads the benchmark and the compiled programs; these tests do not."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

from genomeos.attribution import ablation as ab
from genomeos.attribution import measured as ms

ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location("placement_audit", ROOT / "scripts/placement_audit.py")
pa = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = pa
_SPEC.loader.exec_module(pa)


def _pair(start, end, gene="G", cell="K562", split="training", effect=-0.3, significant=True, power=0.9):
    return ms.CrispriPair(
        chrom="chr1",
        start=start,
        end=end,
        gene=gene,
        cell=cell,
        dataset="D",
        reference="R",
        regulated=significant and effect < 0,
        significant=significant,
        effect_size=effect,
        p_adjusted=0.01 if significant else 0.9,
        split=split,
        power_at_effect_size_20=power,
    )


def _raw(p: ms.CrispriPair, distance: int = 1000) -> dict[str, str]:
    return {"distanceToTSS": str(distance)}


# --- the two rules --------------------------------------------------------------------------------
def test_registered_threshold_is_the_measured_layers():
    assert pa.POLICY_THRESHOLD == ms.RECIPROCAL_OVERLAP
    reg = pa.registration()
    assert reg["policy"]["name"] == "min_side_half"
    assert set(reg["stop_rule"]) >= {"stop_a_nulls_as_fast", "stop_b_not_positional", "stop_c_wide_intervals"}
    assert reg["cause_order"] == list(pa.CAUSES)


def test_a_registry_element_inside_a_wider_tested_interval():
    # a 200 bp element wholly inside a 500 bp tested interval: 0.4 of the tested interval
    assert not pa.current_rule(100, 300, 0, 500)
    assert pa.policy_rule(100, 300, 0, 500)
    assert pa.min_side(100, 300, 0, 500) == 1.0


def test_equal_widths_agree_and_a_graze_is_refused_by_both():
    assert pa.current_rule(0, 300, 100, 400) and pa.policy_rule(0, 300, 100, 400)
    assert not pa.current_rule(0, 300, 290, 800) and not pa.policy_rule(0, 300, 290, 800)


def test_zero_width_meets_nothing():
    assert pa.min_side(10, 10, 0, 500) == 0.0
    assert not pa.policy_rule(0, 500, 10, 10)
    es = pa.ElementSet([("chr1", 0, 500, "a")])
    assert es.touching("chr1", 10, 10) == []


def test_element_set_touching():
    es = pa.ElementSet([("chr1", 0, 200, "a"), ("chr1", 300, 600, "b"), ("chr2", 0, 100, "c")])
    assert [i for _, _, i in es.touching("chr1", 150, 350)] == ["a", "b"]
    assert es.touching("chr1", 200, 300) == []  # half-open: 200 and 300 are not shared bases
    assert es.touching("chr3", 0, 100) == []
    assert (es.min_width, es.max_width) == (100, 300)


# --- the cascade ----------------------------------------------------------------------------------
H_NONE = pa.Hit(0, 0, ())
H_RULE = pa.Hit(1, 1, ("x",))
H_MIN = pa.Hit(1, 1, ())
H_PART = pa.Hit(1, 0, ())


@pytest.mark.parametrize(
    ("comp", "reg", "width", "want"),
    [
        (H_RULE, H_RULE, 500, pa.PLACED),
        (H_MIN, H_RULE, 500, pa.REGISTRY_WOULD),
        (H_NONE, H_NONE, 0, pa.ZERO_WIDTH),
        (H_NONE, H_MIN, 701, pa.UNREACHABLE),
        (H_NONE, H_MIN, 74, pa.UNREACHABLE),
        (H_NONE, H_MIN, 700, pa.WIDTH_MISMATCH),
        (H_NONE, H_MIN, 75, pa.WIDTH_MISMATCH),
        (H_NONE, H_PART, 500, pa.PARTIAL),
        (H_NONE, H_NONE, 500, pa.ABSENT),
    ],
)
def test_cause_cascade(comp, reg, width, want):
    assert pa.cause(comp, reg, width, 150, 350) == want


def test_exclusion_takes_the_state_furthest_along():
    swept = {"a": pa.EXCLUSION[3], "b": pa.EXCLUSION[2]}
    assert pa.exclusion(["a", "b"], set(), swept) == pa.EXCLUSION[2]
    assert pa.exclusion(["z"], set(), swept) == pa.EXCLUSION[4]  # never scored
    assert pa.exclusion(["a", "c"], {"c"}, swept) == pa.EXCLUSION[0]


# --- the links and C4's own numbers ---------------------------------------------------------------
def _links(pairs):
    return pa.links_of(pairs, [_raw(p) for p in pairs])


def test_a_link_is_positive_when_any_record_is_a_decrease():
    pairs = [
        _pair(0, 500),
        _pair(0, 500, effect=-0.01, significant=False),
        _pair(600, 1100, significant=False),
    ]
    links = _links(pairs)
    assert len(links) == 2
    by = {x.start: x for x in links}
    assert by[0].positive and len(by[0].outcomes) == 2
    assert not by[600].positive and by[600].group == ms.NULL_INFORMATIVE


def test_the_cascade_reproduces_ablation_placement():
    comp = [("chr1", 100, 400, "c1")]
    reg = comp + [("chr1", 1000, 1300, "r1"), ("chr1", 2100, 2300, "r2"), ("chr1", 5000, 5150, "r3")]
    pairs = [
        _pair(100, 400),  # placed
        _pair(1000, 1400, gene="H"),  # registry would
        _pair(2000, 2500, gene="I"),  # registry element inside, 0.4 of the tested interval
        _pair(8000, 8500, gene="J"),  # nothing there
        _pair(8000, 8500, gene="K", significant=False),  # a null, not a C4 link
    ]
    got = ab.placement(
        pairs, {"chr1": [(s, e) for _, s, e, _ in comp]}, {"chr1": [(s, e) for _, s, e, _ in reg]}
    )
    rows = pa.place(_links(pairs), pa.ElementSet(comp), pa.ElementSet(reg), pa.current_rule)
    pos = [r for r in rows if r.link.positive]
    assert got["measured_links"] == len(pos) == 4
    assert got["placed_training"] == sum(r.cause == pa.PLACED for r in pos) == 1
    assert got["not_placed_registry_would_training"] == sum(r.cause == pa.REGISTRY_WOULD for r in pos) == 1
    assert got["not_placed_no_registry_element_training"] == 2
    assert {r.link.gene: r.cause for r in rows} == {
        "G": pa.PLACED,
        "H": pa.REGISTRY_WOULD,
        "I": pa.WIDTH_MISMATCH,
        "J": pa.ABSENT,
        "K": pa.ABSENT,
    }
    policy = pa.place(_links(pairs), pa.ElementSet(comp), pa.ElementSet(reg), pa.policy_rule)
    assert {r.link.gene: r.cause for r in policy}["I"] == pa.REGISTRY_WOULD


def test_attrition_keeps_positives_and_nulls_apart():
    pairs = [_pair(0, 500), _pair(0, 500, gene="N", significant=False, power=0.1)]
    rows = pa.place(_links(pairs), pa.ElementSet([]), pa.ElementSet([]), pa.current_rule)
    a = pa.attrition(rows)
    assert a["positive"]["all"] == {pa.ABSENT: 1, "total": 1}
    assert a["null"]["training"] == {pa.ABSENT: 1, "total": 1}
    assert a[ms.NULL_INCONCLUSIVE]["all"]["total"] == 1
    assert pa.cause_families(rows) == {
        "positive": {f: int(f == "absent_from_both_sets") for f in pa.FAMILIES},
        "null": {f: int(f == "absent_from_both_sets") for f in pa.FAMILIES},
    }


# --- the statistics -------------------------------------------------------------------------------
def test_auroc_and_average_precision():
    s = np.array([0.9, 0.8, 0.3, 0.1])
    y = np.array([True, True, False, False])
    assert pa.auroc(s, y) == 1.0 and pa.average_precision(s, y) == 1.0
    assert pa.auroc(np.array([1.0, 1.0]), np.array([True, False])) == 0.5  # a tie counts half
    assert pa.auroc(s, np.array([False] * 4)) is None
    assert pa.average_precision(np.array([0.1, 0.9]), np.array([True, False])) == 0.5


def test_rescue_excess_and_bootstrap():
    pairs = [_pair(0, 500, gene=f"P{i}") for i in range(4)] + [
        _pair(0, 500, gene=f"N{i}", significant=False) for i in range(4)
    ]
    links = _links(pairs)
    empty = pa.place(links, pa.ElementSet([]), pa.ElementSet([]), pa.current_rule)
    full = pa.place(links, pa.ElementSet([("chr1", 0, 500, "e")]), pa.ElementSet([]), pa.current_rule)
    cols = pa.rescue_columns(empty, full)
    sums = {k: float(v.sum()) for k, v in cols.items()}
    assert sums == {"pos_un": 4, "null_un": 4, "pos_res": 4, "null_res": 4}
    assert pa.excess(sums) == 0.0
    lo, hi = pa.cluster_bootstrap([x.gene for x in links], cols, pa.excess)
    assert lo == hi == 0.0
    assert pa._stop(lo) and pa._stop(float("nan")) and not pa._stop(0.01)


def test_coverage_and_ambiguity():
    pairs = [_pair(0, 500), _pair(0, 500, gene="N", significant=False)]
    comp = pa.ElementSet([("chr1", 0, 300, "a"), ("chr1", 250, 500, "b")])
    rows = pa.place(_links(pairs), comp, pa.ElementSet([]), pa.policy_rule)
    assert pa.coverage(rows)["positive"]["coverage_compiled"] == 1.0
    assert pa.ambiguity(rows)["null"] == {"placed": 1, "carried_by_more_than_one_element": 1}


def test_compiled_conventions_find_moved_elements():
    reg = [("chr1", 0, 200, "a"), ("chr1", 300, 500, "b")]
    comp = [("chr1", 0, 200, "a"), ("chr1", 301, 500, "b"), ("chr1", 900, 1000, "z")]
    got = pa.compiled_conventions(comp, reg)
    assert (got["same_id_same_coordinates_in_registry"], got["same_id_other_coordinates"]) == (1, 1)
    assert got["id_not_in_registry"] == 1


def test_width_classes_and_the_descriptive_table():
    assert pa.width_class(74) == "narrower_than_half_the_narrowest"
    assert pa.width_class(75) == pa.width_class(700) == "reachable_width"
    assert pa.width_class(701) == "wider_than_twice_the_widest"
    pairs = [_pair(0, 1000), _pair(0, 1000, gene="N", significant=False)]
    reg = pa.ElementSet([("chr1", 100, 400, "a")])
    links = _links(pairs)
    cur = pa.place(links, reg, reg, pa.current_rule)
    pol = pa.place(links, reg, reg, pa.policy_rule)
    got = pa.by_width(cur, pol)
    assert got["positive"]["wider_than_twice_the_widest"]["rescued"] == 1
    assert got["null"]["wider_than_twice_the_widest"]["placed_current_rule"] == 0
