"""The reference lineage's minute measured against two independent clocks."""

import json
from pathlib import Path

import pytest

from genomeos.organism.reference import ReferenceLineage, distil, parse_wormweb
from genomeos.organism.time_axis import _fit, compare
from genomeos.results import load_result


def _ref():
    node = lambda name, level, total, kids=(), tissue="": {  # noqa: E731
        "id": name.lower(),
        "name": name,
        "did": name,
        "data": {"levelDistance": level, "totalDistance": total, "deathDistance": 0, "type": tissue},
        "children": list(kids),
    }
    tree = node(
        "P0",
        0,
        0,
        [
            node(
                "AB",
                20,
                20,
                [
                    node("ABa", 40, 60, [node("ABal", 60, 120, tissue="neuron")]),
                    node("ABp", 40, 60, tissue="muscle"),
                ],
            ),
            node("P1", 30, 30, [node("EMS", 40, 70, tissue="intestine"), node("P2", 50, 80, tissue="repro")]),
        ],
    )
    return ReferenceLineage.from_dict(distil(parse_wormweb("var json = " + json.dumps(tree) + ";")))


def test_fit_recovers_a_known_scale():
    f = _fit([0, 100, 200, 300], [10, 80, 150, 220])
    assert f["slope"] == 0.7 and f["intercept_min"] == 10.0 and f["residual_sd_min"] == 0.0 and f["n"] == 4


def test_compare_fits_births_divisions_and_midlife_on_a_toy():
    ref = _ref()
    # a clock that runs at 0.5 reference minutes: births and divisions halved, 12 cells so the fits run
    lifetimes = {c.id: (c.born * 0.5, c.divides * 0.5) for c in ref.cells.values() if c.divides is not None}
    for k in range(10):  # pad with consistent synthetic cells so the fit has more than ten points
        lifetimes[f"synthetic{k}"] = (k * 5.0, k * 10.0)
    packer = {c.id: [(c.born + c.divides) / 4] * 6 for c in ref.cells.values() if c.divides is not None}
    for k in range(10):
        packer[f"synthetic{k}"] = [k * 7.5] * 6
    r = compare(ref, lifetimes, packer)
    assert r["atlas_births"] is None or 0.4 <= r["atlas_births"]["slope"] <= 0.6  # only reference cells count
    assert (
        r["packer_midlife"] is None or r["packer_midlife"]["n"] <= 5
    )  # synthetic cells are not in the reference
    assert "reference minute" in r["reading"] or r["reading"].startswith("not enough")


@pytest.mark.skipif(not Path("data/results/celegans_time_axis.json").exists(), reason="distil first")
def test_distilled_time_axis_says_the_reference_runs_slow():
    r = load_result("celegans_time_axis")
    assert r["packer_midlife"]["n"] > 50 and 0.6 <= r["packer_midlife"]["slope"] <= 0.8
    assert r["atlas_births"]["n"] > 300 and r["atlas_divisions"]["n"] > 300
    assert 0.55 <= r["atlas_divisions"]["slope"] <= 0.85
    assert 0.55 <= r["minutes_per_reference_minute"] <= 0.8 and "runs slow" in r["reading"]
