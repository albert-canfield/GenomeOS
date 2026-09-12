"""The second axis of constraint: within-human variation next to the mammalian one (area J step 1)."""

import json
import os

import pytest

from genomeos.attribution.bigwig import BigWig, IntervalStats
from genomeos.attribution.variation import (
    CASE_ORDER,
    ELEMENT_MAMMAL_MIN,
    GNOCCHI_THRESHOLD,
    HUMAN_MIN_FRACTION,
    aggregate,
    case_of,
    classify_blocks,
    classify_elements,
    elements_of,
    merge,
    non_overlapping,
    stats_dict,
    tally_blocks,
    tally_elements,
    vista_by_status,
)


def _st(start, end, bases, above, maximum, total=None):
    total = total if total is not None else maximum * bases
    return IntervalStats(start=start, end=end, bases=bases, total=total, maximum=maximum, above=above)


def test_merge_and_non_overlapping():
    assert merge([(10, 20), (5, 12), (30, 40), (40, 45)]) == [(5, 20), (30, 45)]
    items = [{"start": 0, "end": 10}, {"start": 5, "end": 12}, {"start": 12, "end": 20}]
    kept, dropped = non_overlapping(items)
    assert [k["start"] for k in kept] == [0, 12] and dropped == 1


def test_stats_dict_and_cases():
    assert stats_dict(None) is None
    assert stats_dict(_st(0, 1000, 0, 0, float("-inf"))) is None
    d = stats_dict(_st(0, 4000, 4000, 3000, 4.4))
    assert d["fraction_above"] == 0.75 and d["strong"] is True
    # the four cases, at the block bars
    assert case_of(0.10, 0.60)["case"] == "syntax"
    assert case_of(0.01, 0.00)["case"] == "tolerant"
    assert case_of(0.01, 0.60)["case"] == "recent"
    assert case_of(0.10, 0.00)["case"] == "relaxed"
    # an element uses its own mammalian bar
    assert case_of(0.10, 0.60, ELEMENT_MAMMAL_MIN)["case"] == "recent"
    # unmeasured on either axis gives no case; confidence is a reading, capped at 0.6
    assert case_of(None, 0.5) is None and case_of(0.1, None) is None
    c = case_of(1.0, 1.0)
    assert c["confidence"] == 0.6 and c["mammals"] == "constrained" and c["humans"] == "constrained"
    assert case_of(0.05, HUMAN_MIN_FRACTION)["confidence"] == 0.3


def test_classify_and_tally_blocks():
    blocks = [
        {
            "start": 0,
            "end": 1000,
            "length": 1000,
            "class": "gap",
            "phylop": None,
            "guess": {"tier": "structural"},
        },
        {
            "start": 1000,
            "end": 5000,
            "length": 4000,
            "class": "unique_intergenic",
            "phylop": {"fraction_above": 0.12},
            "guess": {"tier": "constrained_unknown"},
        },
        {
            "start": 5000,
            "end": 9000,
            "length": 4000,
            "class": "unique_intergenic",
            "phylop": {"fraction_above": 0.0},
            "guess": {"tier": "neutral"},
        },
    ]
    stats = [None, _st(1000, 5000, 4000, 0, 0.5), _st(5000, 9000, 4000, 4000, 3.0)]
    rows = classify_blocks(blocks, stats)
    assert rows[0]["case"] is None and rows[0]["gnocchi"] is None
    assert rows[1]["case"]["case"] == "relaxed" and rows[2]["case"]["case"] == "recent"
    t = tally_blocks(rows)
    assert t["structural"]["human_constrained_fraction"] is None
    assert t["constrained_unknown"]["human_constrained_fraction"] == 0.0
    assert t["neutral"]["human_constrained_fraction"] == 1.0
    assert t["neutral"]["cases"]["recent"] == {"blocks": 1, "bp": 4000}
    assert set(t["neutral"]["cases"]) == set(CASE_ORDER)


def test_elements_loaded_and_classified(tmp_path):
    (tmp_path / "constrained_targets_chrT.json").write_text(
        json.dumps(
            {
                "elements": [
                    {
                        "id": "E1",
                        "start": 100,
                        "end": 400,
                        "constrained_fraction": 0.5,
                        "predicted_coding": {"gene": "G1", "log2_fold_change": -0.4},
                    }
                ]
            }
        )
    )
    (tmp_path / "enhancer_targets_chrT.json").write_text(
        json.dumps(
            {
                "elements": [
                    {"id": "E1", "start": 100, "end": 400, "predicted_coding": {"gene": "G1"}},
                    {"id": "E2", "start": 900, "end": 1200, "predicted_coding": {}},
                ]
            }
        )
    )
    els = elements_of("chrT", tmp_path)
    assert [e["id"] for e in els] == ["E1", "E2"]
    assert els[0]["origin"] == "constrained" and els[0]["target"] == "G1"
    assert els[1]["mammal_fraction"] is None
    els[1]["mammal_fraction"] = 0.0
    rows = classify_elements(els, [_st(100, 400, 1000, 1000, 5.0), _st(900, 1200, 1000, 0, -1.0)])
    assert rows[0]["case"]["case"] == "syntax" and rows[0]["gnocchi"]["strong"]
    assert rows[1]["case"]["case"] == "tolerant"
    t = tally_elements(rows)
    assert t["by_case"]["syntax"] == {"elements": 1, "name_a_gene": 1, "strong": 1, "name_a_gene_share": 1.0}
    assert t["by_case"]["tolerant"]["name_a_gene"] == 0 and t["unmeasured"] == 0


def test_aggregate_and_vista():
    sts = [_st(0, 1000, 1000, 1000, 4.5), _st(2000, 3000, 1000, 0, 1.0), _st(4000, 5000, 0, 0, float("-inf"))]
    a = aggregate(sts)
    assert a["track_bases"] == 2000 and a["fraction_above"] == 0.5
    assert a["share_intervals_constrained"] == 0.5 and a["share_intervals_strong"] == 0.5
    rows = [{"status": "positive"}, {"status": "negative"}, {"status": "negative"}]
    v = vista_by_status(rows, sts)
    assert v["positive"]["fraction_above"] == 1.0 and v["negative"]["fraction_above"] == 0.0


@pytest.mark.skipif(
    not os.environ.get("GENOMEOS_LIVE"), reason="set GENOMEOS_LIVE=1 to read the Gnocchi track over HTTP"
)
def test_live_gnocchi_range():  # pragma: no cover - network
    from genomeos.attribution.variation import GNOCCHI_URL

    bw = BigWig(GNOCCHI_URL)
    (st,) = bw.summarise("chr21", [(40_000_000, 40_020_000)], GNOCCHI_THRESHOLD)
    bw.close()
    assert st.bases > 0 and st.maximum > 4.0  # a kilobase in the top percentile sits here
