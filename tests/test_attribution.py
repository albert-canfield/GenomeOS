# SPDX-License-Identifier: AGPL-3.0-or-later
"""The 98%: bigWig sections decoded, constraint distributed over blocks, guesses tiered."""

import json
import os
import struct
from array import array

import pytest

from genomeos.attribution.bigwig import (
    BEDGRAPH,
    FIXED_STEP,
    VARIABLE_STEP,
    BigWig,
    IntervalStats,
    LeafItem,
    coalesce,
    decode_section,
)
from genomeos.attribution.budget import TIERS, build, distil, guess, tallies
from genomeos.attribution.constraint import Element, elements_over_blocks


def _section(kind: int, start: int, items, step: int = 1, span: int = 1) -> bytes:
    head = struct.pack("<IIIIIBBH", 13, start, start + 10 * len(items), step, span, kind, 0, len(items))
    if kind == VARIABLE_STEP:
        body = b"".join(struct.pack("<If", s, v) for s, v in items)
    elif kind == BEDGRAPH:
        body = b"".join(struct.pack("<IIf", s, e, v) for s, e, v in items)
    else:  # fixedStep, and any unknown type for the error path
        body = array("f", items).tobytes()
    return head + body


def test_decode_three_section_types():
    runs = decode_section(_section(FIXED_STEP, 100, [1.0, 2.0, 3.0]))
    assert len(runs) == 1 and runs[0][0] == 100 and runs[0][1] == 1 and list(runs[0][2]) == [1.0, 2.0, 3.0]
    runs = decode_section(_section(VARIABLE_STEP, 0, [(5, 0.5), (9, 2.5)], span=2))
    assert [(s, n, list(v)) for s, n, v in runs] == [(5, 2, [0.5]), (9, 2, [2.5])]
    runs = decode_section(_section(BEDGRAPH, 0, [(10, 20, 4.0)]))
    assert [(s, n, list(v)) for s, n, v in runs] == [(10, 10, [4.0])]
    with pytest.raises(ValueError):
        decode_section(_section(7, 0, [1.0]))


def test_accumulate_counts_only_inside_intervals_and_above_threshold():
    ivs = [(100, 105), (110, 112)]
    stats = [IntervalStats(*iv) for iv in ivs]
    vals = array("f", [0.0, 3.0, 3.0, -1.0, 2.27, 9.0, 9.0, 9.0, 9.0, 9.0, 1.0, 5.0, 9.0])  # bases 100..112
    BigWig._accumulate(100, 1, vals, ivs, [e for _, e in ivs], stats, 2.27)
    a, b = stats
    assert (a.bases, a.above, a.maximum) == (
        5,
        3,
        3.0,
    )  # 3.0, 3.0, 2.27 reach the threshold; 9.0 lies outside
    assert a.mean == pytest.approx((0 + 3 + 3 - 1 + 2.27) / 5)
    assert (b.bases, b.above, b.maximum) == (2, 1, 5.0)
    assert stats[0].as_dict()["fraction_above"] == 0.6


def test_coalesce_joins_neighbours_and_splits_far_blocks():
    items = [LeafItem(0, 0, 0, 0, off, 100) for off in (0, 100, 250, 10_000_000)]
    groups = coalesce(items, max_gap=100, max_size=1000)
    assert [[i.offset for i in g] for g in groups] == [[0, 100, 250], [10_000_000]]


def test_elements_over_blocks_clips_and_counts():
    els = [Element(5, 15, 10), Element(20, 30, 50), Element(95, 105, 7)]
    res = elements_over_blocks(els, [(0, 25), (100, 200)])
    assert res[0] == {"n": 2, "bp": 15, "max_lod": 50, "fraction": 0.6}
    assert res[1] == {"n": 1, "bp": 5, "max_lod": 7, "fraction": 0.05}


def test_guess_tiers():
    assert guess("gap", 1.0, None, None)["tier"] == "structural"
    assert guess("centromere", 0.5, {"fraction_above": 0.0}, None)["tier"] == "structural"
    assert guess("interspersed_repeat_LINE", 0.9, {"fraction_above": 0.01}, None)["tier"] == "fossil"
    g = guess("interspersed_repeat_LINE", 0.9, {"fraction_above": 0.12}, None)
    assert g["tier"] == "constrained_unknown" and "exapted" in g["label"]
    assert guess("regulatory", 0.7, {"fraction_above": 0.08}, None)["tier"] == "regulatory"
    assert guess("unique_intergenic", 0.3, {"fraction_above": 0.01}, {"n": 0})["tier"] == "neutral"
    assert (
        guess("unique_intergenic", 0.3, {"fraction_above": 0.07}, {"n": 0})["tier"] == "constrained_unknown"
    )
    assert guess("unclassified", 0.0, {"fraction_above": 0.035}, {"n": 4})["tier"] == "constrained_unknown"
    assert guess("long_orf", 0.4, {"fraction_above": 0.0}, None)["tier"] == "neutral"
    unmeasured = guess("unique_intergenic", 0.3, None, None)
    assert unmeasured["confidence"] <= 0.3 and "not measured" in unmeasured["label"]
    assert all(
        guess(c, 0.5, {"fraction_above": 0.5}, None)["tier"] in TIERS
        for c in ("regulatory", "mixed_intergenic", "gap")
    )


def _unknown():
    return {
        "unknown_bp": 300,
        "blocks": [
            {"start": 0, "end": 100, "length": 100, "class": "gap", "evidence": "curated", "confidence": 1.0},
            {
                "start": 100,
                "end": 200,
                "length": 100,
                "class": "regulatory",
                "evidence": "curated",
                "confidence": 0.7,
            },
            {
                "start": 200,
                "end": 300,
                "length": 100,
                "class": "unique_intergenic",
                "evidence": "inferred",
                "confidence": 0.3,
            },
        ],
    }


def test_build_without_network_and_tallies():
    out = build("chrTest", phylop=False, elements=False, unknown=_unknown(), length=1000)
    assert [r["guess"]["tier"] for r in out["blocks"]] == ["structural", "regulatory", "constrained_unknown"]
    assert out["by_tier"]["structural"]["fraction_of_unknown"] == pytest.approx(1 / 3, abs=1e-4)
    assert out["by_tier"]["structural"]["fraction_of_chromosome"] == 0.1
    assert out["constrained_fraction"] is None  # nothing measured
    assert out["sources"] == {"phylop": None, "elements": None}


def test_tallies_and_distil(tmp_path):
    rows = [
        {"length": 100, "class": "gap", "phylop": None, "guess": {"tier": "structural", "confidence": 1.0}},
        {
            "length": 200,
            "class": "regulatory",
            "phylop": {"above": 20, "bases": 200},
            "guess": {"tier": "regulatory", "confidence": 0.7},
        },
        {
            "length": 700,
            "class": "unique_intergenic",
            "phylop": {"above": 7, "bases": 700},
            "guess": {"tier": "neutral", "confidence": 0.6},
        },
    ]
    t = tallies(rows, 2000)
    assert t["constrained_bp"] == 27 and t["measured_bp"] == 900 and t["constrained_fraction"] == 0.03
    assert t["by_class"]["regulatory"]["constrained_fraction"] == 0.1
    assert t["guessed_fraction"] == 1.0
    res = {
        "result": "budget_chr21",
        "chromosome_length": 2000,
        "unknown_bp": 1000,
        "cost": {"phylop": {"mb_fetched": 3.5}},
        **t,
    }
    (tmp_path / "budget_chr21.json").write_text(json.dumps(res))
    (tmp_path / "budget_chr22.json").write_text(json.dumps({**res, "result": "budget_chr22"}))
    d = distil(tmp_path)
    assert d["chromosomes"] == 2 and d["genome_bp"] == 4000 and d["unknown_bp"] == 2000
    assert d["constrained_bp"] == 54 and d["by_tier"]["neutral"]["bp"] == 1400
    assert d["by_tier"]["neutral"]["fraction_of_genome"] == 0.35
    assert d["phylop_mb_fetched"] == 7.0


@pytest.mark.skipif(
    not os.environ.get("GENOMEOS_LIVE"), reason="set GENOMEOS_LIVE=1 to read the Zoonomia track over HTTP"
)
def test_live_zoonomia_range():  # pragma: no cover - network
    from genomeos.attribution.constraint import PHYLOP_241_URL

    bw = BigWig(PHYLOP_241_URL)
    (st,) = bw.summarise("chr21", [(25_880_000, 25_881_000)], 2.27)
    assert st.bases == 1000 and 0.1 < st.fraction_above < 0.3  # APP exon neighbourhood
