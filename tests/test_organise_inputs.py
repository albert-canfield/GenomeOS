# SPDX-License-Identifier: AGPL-3.0-or-later
"""`organise.inputs()` declares what a block-organiser result depends on, and a declaration that
shrinks when an input is missing is the defect these tests plant.

Every absence here is PLANTED in a `tmp_path` results directory. A test that asserted the raise by
pointing at a path that happens not to exist on the machine running it would pass for a reason that
has nothing to do with the code: it would go green on this machine today and silently stop testing
anything the moment the file arrived. So each test writes the result, writes or withholds the file
the result points at, and asserts on the message.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import organise
from genomeos.attribution.targets import RUNS

RUN = RUNS[0]


def _budget(d: Path, chrom: str = "chr21") -> None:
    """The one input `inputs()` requires: enough of a budget record for a declaration to exist."""
    doc = {"blocks": [{"start": 0, "end": 1000, "length": 1000, "class": "unknown"}]}
    (d / f"budget_{chrom}.json").write_text(json.dumps(doc))


def _summary_with_pointer(d: Path, where: str, chrom: str = "chr21") -> Path:
    p = d / f"{RUN}_{chrom}.json"
    p.write_text(json.dumps({"elements_where": where, "elements_total": 1}))
    return p


# --- site 1: the pointed-at table ------------------------------------------------------------


def test_a_pointer_whose_table_is_absent_raises_and_names_the_path(tmp_path: Path) -> None:
    """PLANTED: the summary is there and declares a table; the table is not written."""
    _budget(tmp_path)
    summary = _summary_with_pointer(tmp_path, "data/results/elements/absent_chr21.json.gz")
    with pytest.raises(FileNotFoundError) as e:
        organise.inputs("chr21", tmp_path)
    msg = str(e.value)
    assert "data/results/elements/absent_chr21.json.gz" in msg, msg
    assert str(summary) in msg, msg


def test_the_same_pointer_is_declared_when_the_table_is_planted(tmp_path: Path) -> None:
    """The positive control: identical to the test above except the table exists, so the raise is
    the absence and not the shape of the record. Without this, a function that raised on every
    pointer would pass the test above."""
    _budget(tmp_path)
    table = tmp_path / "elements_chr21.json.gz"
    table.write_bytes(b"")
    _summary_with_pointer(tmp_path, str(table))
    got = organise.inputs("chr21", tmp_path)
    assert table in got
    assert got.count(table) == 1


def test_a_pointer_resolved_against_the_repository_root_is_still_declared(tmp_path: Path) -> None:
    """A relative pointer is resolved against `results_dir.parent.parent` before it is judged, so
    the raise must not fire on a table that is merely named relatively."""
    results = tmp_path / "data" / "results"
    results.mkdir(parents=True)
    _budget(results)
    table = tmp_path / "data" / "elements" / "chr21.json.gz"
    table.parent.mkdir(parents=True)
    table.write_bytes(b"")
    _summary_with_pointer(results, "data/elements/chr21.json.gz")
    assert table in organise.inputs("chr21", results)


def test_a_summary_with_no_pointer_at_all_declares_only_itself(tmp_path: Path) -> None:
    """Absence of an `elements_where` field is not a missing input: it means the run keeps its
    elements inline. Only a pointer that was MADE and does not resolve raises."""
    _budget(tmp_path)
    (tmp_path / f"{RUN}_chr21.json").write_text(json.dumps({"elements_total": 0}))
    got = organise.inputs("chr21", tmp_path)
    assert got == [tmp_path / "budget_chr21.json", tmp_path / f"{RUN}_chr21.json"]


# --- site 2: the budget pair, which `blocks()` cannot do without --------------------------------


def test_neither_budget_file_raises_and_names_both_paths(tmp_path: Path) -> None:
    """PLANTED: a results directory with a variation reading and no budget at all. `read_axes`
    returns None here and `blocks()` raises, so a declaration would describe a join that cannot be
    built."""
    (tmp_path / "variation_chr21.json").write_text(json.dumps({"blocks": []}))
    with pytest.raises(FileNotFoundError) as e:
        organise.inputs("chr21", tmp_path)
    msg = str(e.value)
    assert str(tmp_path / "budget_axes_chr21.json") in msg, msg
    assert str(tmp_path / "budget_chr21.json") in msg, msg


@pytest.mark.parametrize("name", ["budget_chr21", "budget_axes_chr21"])
def test_either_budget_file_alone_is_enough(tmp_path: Path, name: str) -> None:
    """The pair is an either/or in `read_axes`, not two required inputs, so neither one alone may
    raise. This is the test that would fail if the fix had been applied per path."""
    (tmp_path / f"{name}.json").write_text(json.dumps({"blocks": []}))
    assert organise.inputs("chr21", tmp_path) == [tmp_path / f"{name}.json"]


def test_an_absent_per_chromosome_reading_is_a_smaller_join_and_not_a_raise(tmp_path: Path) -> None:
    """`blocks()` reads variation, duplication and the attribution runs tolerantly (`or {}` and
    `run_elements`), so a chromosome that has not had them run is declared short on purpose. If this
    ever starts raising, a partial genome has become unwritable and the change went too far."""
    _budget(tmp_path)
    got = organise.inputs("chr21", tmp_path)
    assert got == [tmp_path / "budget_chr21.json"]
    for absent in ("variation_chr21", "duplication_chr21", *(f"{r}_chr21" for r in RUNS)):
        assert not (tmp_path / f"{absent}.json").exists()


# --- site 3: the result's own `inputs` block, which was a constant --------------------------------


def _runs(d: Path, chrom: str = "chr21") -> None:
    for r in RUNS:
        (d / f"{r}_{chrom}.json").write_text(json.dumps({"elements": []}))


def test_the_declaration_names_the_budget_file_that_actually_answered(tmp_path: Path) -> None:
    """PLANTED both ways. `read_axes` prefers `budget_axes_<chrom>`; the old block named
    `budget_<chrom>` as a constant, so it was wrong on every chromosome where both are present --
    which is all 24. A constant cannot pass both halves of this test."""
    from genomeos.attribution.budget import read_axes
    from genomeos.attribution.organise import _declared

    blocks = {"blocks": [], "by_tier": {}}
    (tmp_path / "budget_chr21.json").write_text(json.dumps(blocks))
    stored_only = _declared("chr21", tmp_path, read_axes("chr21", tmp_path))
    assert stored_only["budget"] == "budget_chr21"

    (tmp_path / "budget_axes_chr21.json").write_text(json.dumps(blocks))
    both = _declared("chr21", tmp_path, read_axes("chr21", tmp_path))
    assert both["budget"] == "budget_axes_chr21", both
    assert both["budget_branch"] == "budget_axes_chr21"
    assert both["budget_pair_considered"] == ["budget_axes_chr21", "budget_chr21"]


def test_every_attribution_run_in_RUNS_is_declared_not_a_hardcoded_two(tmp_path: Path) -> None:
    """PLANTED: all three runs written. The old block listed `constrained_targets` and
    `enhancer_targets` literally and omitted `enhancer_targets_all`, which carries most of the
    elements and the only pointer at a local table. This test fails against any literal list."""
    from genomeos.attribution.organise import _declared

    _budget(tmp_path)
    _runs(tmp_path)
    got = _declared("chr21", tmp_path, {"axes_source": "budget_chr21"})
    assert got["elements"] == [f"{r}_chr21" for r in RUNS]
    assert "enhancer_targets_all_chr21" in got["elements"]


def test_an_absent_conditional_input_is_recorded_and_not_claimed(tmp_path: Path) -> None:
    """PLANTED: variation and duplication withheld, the runs present. The old block named them
    whether or not they were there, so the record claimed a human axis and a copy flag a reader
    could not have had. Neither dropping them nor claiming them is right: the branch is recorded."""
    from genomeos.attribution.organise import _declared

    _budget(tmp_path)
    _runs(tmp_path)
    got = _declared("chr21", tmp_path, {"axes_source": "budget_chr21"})
    assert got["variation"] is None and got["duplication"] is None
    assert set(got["inputs_optional_absent"]) == {"variation_chr21", "duplication_chr21"}
    for stem, note in got["inputs_optional_absent"].items():
        assert "not on this machine" in note, (stem, note)


def test_an_absent_run_is_recorded_and_the_present_ones_still_declared(tmp_path: Path) -> None:
    """PLANTED: one run withheld. A chromosome that has not had a run is a smaller join, so it must
    not raise -- but which run is missing has to be in the record."""
    from genomeos.attribution.organise import _declared

    _budget(tmp_path)
    for r in RUNS[1:]:
        (tmp_path / f"{r}_chr21.json").write_text(json.dumps({"elements": []}))
    got = _declared("chr21", tmp_path, {"axes_source": "budget_chr21"})
    assert got["elements"] == [f"{r}_chr21" for r in RUNS[1:]]
    assert f"{RUNS[0]}_chr21" in got["inputs_optional_absent"]


def test_the_pointed_at_table_reaches_the_declaration(tmp_path: Path) -> None:
    """PLANTED: a run summary pointing at a table that is written. The table is git-ignored and
    lives only on the machine that made the result, so a declaration that omits it cannot be
    rebuilt against. It must appear by path in `declared`."""
    from genomeos.attribution.organise import _declared

    _budget(tmp_path)
    table = tmp_path / "elements_chr21.json.gz"
    table.write_bytes(b"")
    _summary_with_pointer(tmp_path, str(table))
    got = _declared("chr21", tmp_path, {"axes_source": "budget_chr21"})
    assert str(table) in got["declared"]


def test_the_declaration_raises_when_a_pointed_at_table_is_absent(tmp_path: Path) -> None:
    """PLANTED: the raise from site 1 must propagate to the writer rather than be swallowed while
    the result's own block is assembled. `targets.run_elements` already raises on this same
    condition, so a declaration that tolerated it would be laxer than the read it describes."""
    from genomeos.attribution.organise import _declared

    _budget(tmp_path)
    _summary_with_pointer(tmp_path, "data/results/elements/absent_chr21.json.gz")
    with pytest.raises(FileNotFoundError, match="absent_chr21"):
        _declared("chr21", tmp_path, {"axes_source": "budget_chr21"})
