# SPDX-License-Identifier: AGPL-3.0-or-later
"""The not-open profile types the rules the context-evidence census counted, and nothing else.

Every check here exists to stop this profile describing a different population, or a different band,
from the reading it claims to describe: the rule enumeration is held to `context_evidence.rule_loci`
locus for locus, the chr21 state counts to the registered chr21 census, the effect band to the field
the cached prediction already carries, the distance bands to `executor._band`, and the
informativeness of a cell to the tolerance the base-rate result fixed before any share was computed.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import context_evidence as ce
from genomeos.attribution import executor as ex
from genomeos.attribution import not_open_profile as prof
from genomeos.predict.enhancer_target import STRONG_EFFECT
from tests.committed_data import committed, must_be_committed

CHR21_RELATIVE = "data/results/budget_chr21.json"
CHR21 = Path(CHR21_RELATIVE)
CENSUS21_RELATIVE = "data/results/context_evidence_chr21.json"
CENSUS21 = Path(CENSUS21_RELATIVE)
BASERATE_RELATIVE = "data/results/context_evidence_baserate.json"
BASERATE = Path(BASERATE_RELATIVE)
# A MIXED guard, SPLIT rather than converted whole. The condition named two paths of opposite
# kinds: `data/results/budget_chr21.json`, which git TRACKS, and `ce.TRACK_METADATA`
# (`data/cache/entex/alphagenome_track_metadata_copy.csv`), which git IGNORES. One skip reason
# covered both, so a deleted committed budget was indistinguishable from a machine without the
# track-metadata cache. The ignored half keeps its skip -- by name, as `tests/local_data.py`
# requires -- and the tracked half FAILS.
needs_chr21 = pytest.mark.skipif(
    not ce.TRACK_METADATA.exists(),
    reason=(
        f"NOT RUN HERE, not passed: {ce.TRACK_METADATA} is git-ignored machine-local data and is "
        f"absent from this checkout"
    ),
)


@pytest.fixture(autouse=True)
def _chr21s_budget_is_committed() -> None:
    """The tracked half of the old mixed guard: git tracks chr21's budget, so its absence FAILS."""
    must_be_committed(CHR21_RELATIVE)


def test_the_bands_are_imported_and_not_restated():
    assert prof.STRONG_EFFECT is STRONG_EFFECT
    assert prof.INNER_BAND_OF is ex._band
    assert prof.OUTER_LIMIT == ex.WIDE_MAX_DISTANCE
    # the three bands the project already uses, in its own words, plus one residual row above its top
    assert prof.DISTANCE_BANDS[:3] == (ex._band(0), ex._band(10_000), ex._band(100_000))
    assert len(set(prof.DISTANCE_BANDS)) == 4


def test_distance_band_delegates_below_the_project_s_top_and_is_residual_above():
    for d in (0, 1, 4_999, 5_000, 49_999, 50_000, prof.OUTER_LIMIT - 1):
        assert prof.distance_band(d) == ex._band(d)
    assert prof.distance_band(prof.OUTER_LIMIT) == prof.OUTER_BAND
    assert prof.distance_band(10_000_000) == prof.OUTER_BAND
    assert prof.distance_band(None) == prof.NO_DISTANCE


def test_no_name_this_profile_introduces_carries_a_verdict():
    names = [
        *prof.SOURCES,
        *prof.EFFECT_BANDS,
        *prof.DISTANCE_BANDS,
        prof.OUTER_BAND,
        *prof.Rule(
            source=prof.SOURCE_PREDICTED,
            element="E",
            chrom="chr21",
            start=1,
            end=2,
            cell="K562",
            gene="G",
            element_class="enhancer",
            activity_axis="activates_target",
            origin_axis="unique",
            effect=-1.0,
            effect_band=prof.STRONG,
            distance=1000,
        )
        .row()
        .keys(),
    ]
    for name in names:
        for stem in ce.FORBIDDEN_IN_A_STATE_NAME:
            assert stem not in name.lower(), f"{name!r} carries the verdict stem {stem!r}"


def test_the_limitations_are_the_reading_s_own_words():
    lim = prof.limitations()
    assert lim["not_closed"] == ce.NOT_CLOSED
    assert lim["not_validation"] == ce.NOT_VALIDATION
    # the base-rate figure may never be quoted without the label that it is not independent
    assert "not independent" in lim["base_rate_is_descriptive_and_not_independent"]
    assert "no verdict is moved" in lim["descriptive_only"]


@committed(BASERATE_RELATIVE)
def test_informativeness_is_read_from_the_base_rate_result():
    r = json.loads(BASERATE.read_text())
    info = prof.informativeness()
    assert info.tolerance == r["tolerance"]
    assert info.difference == {c: v["difference"] for c, v in r["per_assigned_cell"].items()}
    for cell, diff in info.difference.items():
        assert info.informative(cell) == (diff > info.tolerance)
    # the cell whose assignment carries no information about openness is not called informative
    assert not info.informative("SK-N-SH")
    assert info.informative(None) is False
    # the three cells lane-context2 singled out are all cells the base-rate result carries
    assert set(prof.CELLS_LANE_CONTEXT2_SINGLED_OUT) <= set(info.difference)


def test_an_axis_value_is_the_compiler_s_own_join():
    from genomeos.attribution import compile as cp

    axes = {
        "activity": [["activates_target"]],
        "molecular_role": [["enhancer_like"], list(cp.REPRESSION_ALTERNATIVES)],
    }
    # the value is exactly what the compiler writes after the axis name on the element
    written = {line.split(": ", 1)[0]: line.split(": ", 1)[1] for line in cp.axis_lines(axes)}
    for name, value in written.items():
        assert prof.axis_value(axes, name) == value
    # `,` where every group holds, `|` where a group is unresolved alternatives
    assert prof.axis_value(axes, "activity") == "activates_target"
    assert prof.axis_value(axes, "molecular_role") == "enhancer_like, " + "|".join(cp.REPRESSION_ALTERNATIVES)
    assert prof.axis_value(axes, "origin") == "unknown"


@needs_chr21
def test_every_rule_carries_the_compiled_class_and_axes():
    for r in prof.rules("chr21"):
        assert r.element_class
        assert r.activity_axis and r.origin_axis
        assert r.row()["activity_axis"] == r.activity_axis


@needs_chr21
def test_the_rules_are_the_ones_the_census_counted():
    mine = prof.rules("chr21")
    loci = ce.rule_loci("chr21")
    assert [(r.cell, r.start, r.end) for r in mine] == loci
    assert all(r.source in prof.SOURCES for r in mine)
    assert all(r.gene for r in mine)


@needs_chr21
def test_the_effect_band_is_the_prediction_s_own_field():
    for r in prof.rules("chr21"):
        if r.source != prof.SOURCE_PREDICTED:
            assert r.effect_band is None and r.effect is None
            continue
        assert r.effect is not None
        # the band is read off the cached prediction; it must still agree with the cut-off it names
        assert r.effect_band == (prof.STRONG if abs(r.effect) >= prof.STRONG_EFFECT else prof.WEAK)


@needs_chr21
@committed(CENSUS21_RELATIVE)
def test_chr21_states_reproduce_the_registered_census():
    census = json.loads(CENSUS21.read_text())
    readers = ce.Readers()
    table = ce.mapping()
    counts: dict[str, int] = dict.fromkeys(ce.STATES, 0)
    for r in prof.rules("chr21"):
        state, _rest = ce.parse_value(ce.state_for(r.cell, "chr21", r.start, r.end, readers, table))
        counts[state] += 1
    assert sum(counts.values()) == census["rules"]
    assert counts == census["per_state"]
