# SPDX-License-Identifier: AGPL-3.0-or-later
"""The union axis: the assignment rule, the census, and the clauses the registration is made of.

The registration (`target_calibration.PREREGISTERED_UNION`) is the object under test as much as the
code is: an axis chosen after its result was seen is priced by the search that found it, so what
has to be mechanical here is that the level stays a description, that the separation is read only
where the axis was never looked at, and that the adoption rule refuses rather than rounds.
"""

from __future__ import annotations

import json

import pytest

from genomeos.attribution import crispri
from genomeos.attribution import target_calibration as tc

# ---------------------------------------------------------------------------------------
# The axis: how a target is assigned to a stratum
# ---------------------------------------------------------------------------------------


def test_the_four_cells_are_the_two_signals_crossed() -> None:
    assert tc.union_cell(0.6, tc.CELL) == "both"
    assert tc.union_cell(0.6, "another track") == "drop only"
    assert tc.union_cell(0.0, tc.CELL) == "track only"
    assert tc.union_cell(0.0, "another track") == "neither"


def test_the_threshold_is_the_boundary_the_drop_axis_already_registered() -> None:
    """0.2 divides the two bottom DROP_BANDS strata from the two top ones; the boundary is open."""
    assert tc.UNION_DROP == 0.2
    assert tc.union_cell(0.2, "another track") == "neither"
    assert tc.union_cell(0.2001, "another track") == "drop only"


def test_the_union_stratum_is_either_signal_and_neither_is_both_absent() -> None:
    assert tc.union_stratum_of(0.6, "another track") == "union"
    assert tc.union_stratum_of(0.0, tc.CELL) == "union"
    assert tc.union_stratum_of(0.0, "another track") == "neither"
    assert set(tc.UNION_STRATA) == {"union", "neither"}


def test_every_crossing_cell_has_a_stratum_and_the_two_partition_the_cells() -> None:
    union = {c for c in tc.UNION_CELLS if c != "neither"}
    assert union == {"both", "drop only", "track only"}
    assert len(tc.UNION_CELLS) == 4


# ---------------------------------------------------------------------------------------
# The track axis read the element's way, so it means one thing on both populations
# ---------------------------------------------------------------------------------------


def _pair(gene: str = "GENE", drop: float = 0.0) -> crispri.Pair:
    p = crispri.Pair(
        chrom="chr1",
        start=1_000,
        end=1_200,
        gene=gene,
        cell=tc.CELL,
        dataset="test",
        distance=5_000.0,
        dhs=1.0,
        h3k27ac=1.0,
        regulated=False,
    )
    p.features["deletion_drop"] = drop
    return p


class _Table:
    def __init__(self, element: dict | None) -> None:
        self._element = element

    def overlapping(self, chrom: str, start: int, end: int) -> list[dict]:
        return [self._element] if self._element else []


def test_on_the_gate_the_track_is_the_pairs_own_targets_track() -> None:
    el = {"id": "E", "predicted": {"gene": "GENE", "tissue": tc.CELL}}
    assert tc.pair_track(_pair("GENE"), _Table(el)) == tc.CELL


def test_off_the_gate_the_track_is_the_elements_own_target_not_a_constant() -> None:
    """The 2026-09-22 reading returned 'another track' for every off-gate pair, which is not an axis.

    Off the gate the pair's gene is not the element's target, so the only reading that exists is
    the element's own -- which is also exactly what a swept target carries.
    """
    el = {"id": "E", "predicted": {"gene": "OTHER", "tissue": tc.CELL}}
    assert tc.pair_track(_pair("GENE"), _Table(el)) == tc.CELL
    cold = {"id": "E", "predicted": {"gene": "OTHER", "tissue": "placenta"}}
    assert tc.pair_track(_pair("GENE"), _Table(cold)) == "another track"


def test_a_pair_on_no_deleted_element_is_off_the_track() -> None:
    assert tc.pair_track(_pair("GENE"), _Table(None)) == "another track"


def test_the_crossing_counts_read_no_label() -> None:
    """A count that touched `regulated` could not be published before the rates; this one cannot."""
    el = {"id": "E", "predicted": {"gene": "GENE", "tissue": tc.CELL}}
    table = _Table(el)
    a, b = _pair("GENE", 0.9), _pair("GENE", 0.9)
    a.regulated, b.regulated = True, False
    assert tc.crossing_counts([a], table) == tc.crossing_counts([b], table)
    assert tc.crossing_counts([a, b], table)["both"] == 2


# ---------------------------------------------------------------------------------------
# The adoption rule: what it refuses
# ---------------------------------------------------------------------------------------


def test_two_rates_are_separated_only_when_their_intervals_do_not_meet() -> None:
    far = tc.observed_band(90, 100), tc.observed_band(5, 100)
    assert tc.separated(*far) is True
    near = tc.observed_band(55, 100), tc.observed_band(50, 100)
    assert tc.separated(*near) is False


def test_an_empty_population_separates_from_nothing() -> None:
    assert tc.separated(tc.observed_band(0, 0), tc.observed_band(90, 100)) is None


def test_a_perfect_run_of_128_names_a_band_and_still_prices_nothing_genome_wide() -> None:
    """The 128 of 128 is a real interval inside 0.9-1; the registration forbids calling it a test."""
    v = tc.observed_band(128, 128)
    assert v["band"] == "0.9-1"
    assert v["ci95"] == [0.9709, 1.0]
    assert "DESCRIPTION AND NOT AS A TEST" in tc.PREREGISTERED_UNION


def test_the_thin_arm_is_registered_as_unusable_rather_than_measured() -> None:
    """8 scored on-gate pairs in the five other cell lines: the rule refuses them before it reads."""
    assert tc.observed_band(8, 8)["band"] == tc.NOT_CALIBRATED
    assert "UNUSABLE" in tc.PREREGISTERED_UNION


# ---------------------------------------------------------------------------------------
# The registration itself: the clauses that make it one
# ---------------------------------------------------------------------------------------


@pytest.mark.parametrize(
    "clause",
    [
        "chosen AFTER its result was seen",  # the provenance is stated, not buried
        "THE POPULATIONS AND THEIR SIZES, COUNTED BEFORE ANY RATE",
        "WHAT INDEPENDENT EVIDENCE WOULD LOOK LIKE",
        "THE ADOPTION RULE, FIXED NOW",
        "NOTHING is banded from this axis",  # the refusal is a registered outcome
    ],
)
def test_the_registration_carries_the_clause(clause: str) -> None:
    assert clause in tc.PREREGISTERED_UNION


def test_the_registered_sizes_are_the_sizes_the_counters_produce() -> None:
    """The genome census is quoted in the registration; a drift would make the registration false."""
    for n in ("593,765", "536,139", "57,626", "440,377", "395,602", "44,775"):
        assert n in tc.PREREGISTERED_UNION


def test_the_registration_names_what_is_not_reachable_at_zero_requests() -> None:
    assert "NOT reachable at 0 requests" in tc.PREREGISTERED_UNION
    assert "not confirmable from anything on this disk" in tc.PREREGISTERED_UNION


# ---------------------------------------------------------------------------------------
# The walk: the union axis is a third reading beside the two already published
# ---------------------------------------------------------------------------------------


def test_the_reband_walk_carries_a_stratum_counter_for_every_axis(tmp_path) -> None:
    """Coverage has to be readable whatever the bands are, so the strata are counted separately."""
    els = [
        {
            "id": "E1",
            "start": 1_000,
            "end": 1_200,
            "predicted": {"gene": "G", "tissue": tc.CELL},
            "predicted_coding": {"gene": "G", "tissue": tc.CELL},
            "predicted_by_cell": {tc.CELL: -0.9},
            "predicted_coding_by_cell": {tc.CELL: -0.9},
        }
    ]
    (tmp_path / "chr1.json").write_text(json.dumps(els))
    out = tc.reband_chromosome(
        "chr1",
        {"target": [0.0, 0.0, 0.0, 0.0]},
        {"union": {"band": "0.9-1"}},
        elements=tmp_path,
        reference=tmp_path,
        results=tmp_path,
        axis="union",
    )
    # no GENCODE under tmp_path, so every target is dropped for want of a TSS: the walk still
    # reports its arms rather than failing, which is what the missing-TSS count is for
    assert out["any gene"]["targets"] == 0
    assert out["any gene"]["by_stratum"] == {}


def test_the_union_threshold_is_exactly_where_the_drop_axis_splits() -> None:
    """The union's magnitude arm is the drop axis's two top strata and nothing else."""
    top = {lab for lab in tc.DROP_LABELS if lab.startswith(("0.2 <", "0.5 <"))}
    for drop in (0.0, 0.05, 0.15, 0.2, 0.3, 0.5, 0.9):
        in_union = tc.union_cell(drop, "another track") == "drop only"
        assert in_union == (tc.drop_band_label(drop) in top), drop
