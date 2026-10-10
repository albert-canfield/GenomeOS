# SPDX-License-Identifier: AGPL-3.0-or-later
"""The panel gate for the polygenic scores lane-pgs left UNKNOWN (genome/pgs_panel.py)."""

from __future__ import annotations

from pathlib import Path

import pytest

from genomeos.genome import pgs_panel as pp
from genomeos.genome import pgs_reference as pr


def test_gate_is_frozen_at_the_registered_values():
    assert pp.GATE["min_people"] == 100
    assert pp.GATE["access"] == "open"
    assert pp.GATE["coverage_floor"] == pr.RULES["coverage_floor"] == 0.80
    assert set(pp.LEFT_UNKNOWN) == {"PGS000297", "PGS000018"}
    assert set(pp.CROSS_CHECK) == {"PGS000004", "PGS000115"}


def test_the_person_keeps_the_declared_group():
    assert pp.declared_group("HG002") == "asj"
    assert pp.PERSON == "HG002"


def test_labels_name_groups_only_through_the_panels_own_label():
    assert pp.group_of("aadr_ho", "Jew_Ashkenazi.HO") == "asj"
    assert pp.group_of("ebc_jew2010", "Ashkenazi_Jews") == "asj"
    assert pp.group_of("1kgp3", "CEU") == "nfe"
    assert pp.group_of("1kgp3", "FIN") == "fin"
    assert pp.group_of("1kgp3", "ACB") == "afr"
    # A label that is not a frozen code names no group, however close its people may be.
    assert pp.group_of("hgdp", "French") is None
    assert pp.group_of("aadr_ho", "Jew_Iraqi.HO") is None
    assert pp.group_of("sgdp", "CEU") is None


def test_the_person_and_relatives_and_related_people_are_never_counted():
    people = [
        {"id": "HG002", "label": "Ashkenazi"},
        {"id": "HG003", "label": "Ashkenazi"},
        {"id": "HG004", "label": "Ashkenazi"},
        {"id": "X1", "label": "Ashkenazi"},
        {"id": "X2", "label": "Ashkenazi", "related": True},
    ]
    assert [p["id"] for p in pp.counted(people)] == ["X1"]
    assert pp.group_counts("giab", people)["asj"] == 1


def test_evaluate_reads_access_before_size_and_needs_the_minimum():
    open_c = {"key": "a", "access": "open", "terms": "t", "data": "array genotypes"}
    shut = {"key": "b", "access": "controlled", "terms": "EGA", "data": "whole-genome sequence calls"}
    assert pp.evaluate(open_c, "asj", 100)["meets_size_access_label"]
    assert not pp.evaluate(open_c, "asj", 99)["meets_size_access_label"]
    e = pp.evaluate(shut, "asj", 128)
    assert not e["meets_size_access_label"] and e["fails"][0].startswith("access")
    assert "size: not counted" in pp.evaluate(open_c, "asj", None)["fails"]
    with pytest.raises(ValueError):
        pp.evaluate(open_c, "xyz", 5)
    assert pp.outcome([pp.evaluate(open_c, "asj", 7), pp.evaluate(shut, "asj", 128)], "asj") == (
        "NO_PANEL_QUALIFIES"
    )
    assert pp.outcome([pp.evaluate(open_c, "nfe", 404)], "nfe") == "PANEL_QUALIFIES"


def test_candidates_are_unique_and_carry_their_access():
    keys = [c["key"] for c in pp.CANDIDATES]
    assert len(keys) == len(set(keys))
    assert all(c["access"] in pp.ACCESS_LEVELS for c in pp.CANDIDATES)
    assert {"1kgp3", "hgdp", "sgdp", "hprc2", "giab", "aadr_ho", "tagc"} <= set(keys)


def test_committed_inventory_stops_the_lane_on_its_own_counts():
    import json

    import tracked_paths as tp

    path = Path("data/results/pgs_panel_inventory.json")
    tp.must_be_present(path)
    res = json.loads(path.read_text())
    assert res["gate"] == pp.GATE
    assert res["declared_group"] == "asj"
    assert res["result_manifest"]["complete"] is True
    assert res["scores_left_unknown"] == pp.LEFT_UNKNOWN

    # the outcome is the one the committed counts give under the committed gate
    evaluations = [
        pp.evaluate(
            {
                "key": c["candidate"],
                "access": c["access"],
                "terms": c["terms"],
                "data": c["data"],
            },
            "asj",
            c["declared_group_people"],
        )
        for c in res["candidates"]
    ]
    assert res["outcome"] == pp.outcome(evaluations, "asj") == "NO_PANEL_QUALIFIES"
    assert all(
        e["meets_size_access_label"] == c["meets_size_access_label"]
        for e, c in zip(evaluations, res["candidates"], strict=True)
    )

    # the largest open-access count is read from the candidates, and is below the minimum
    openly = [
        c["declared_group_people"]
        for c in res["candidates"]
        if c["access"] == "open" and c["declared_group_people"] is not None
    ]
    assert res["largest_count_open_access"]["people"] == max(openly) < pp.GATE["min_people"]

    # no panel person is named, and the reading says the person's UNKNOWN stands
    text = path.read_text()
    for name in ("ashkenazy", "GRC1211", "SAKD", "HGDP00001", "NA21142"):
        assert name not in text
    assert "UNKNOWN" in res["reading"] and "not a risk" in res["note"]
