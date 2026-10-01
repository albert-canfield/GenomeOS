# SPDX-License-Identifier: AGPL-3.0-or-later
"""The reporting logic of scripts/re2g_paired.py, on synthetic pairs only.

No prediction file and no benchmark table is read. What is checked is that the three kinds of zero stay
apart, that the coverage reconciliation does not confuse two different counts, and that the import
closure keeps the case the directories keep.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from genomeos.attribution import crispri

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def paired():
    spec = importlib.util.spec_from_file_location("re2g_paired", ROOT / "scripts" / "re2g_paired.py")
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def make(answered: bool, drop: float, covered: bool, regulated: bool) -> crispri.Pair:
    p = crispri.Pair(
        chrom="chr1",
        start=100,
        end=600,
        gene="G",
        cell="K562",
        dataset="synthetic",
        distance=10_000.0,
        dhs=1.0,
        h3k27ac=1.0,
        regulated=regulated,
        weight=1.0,
    )
    p.covered = covered
    p.features = {"deletion_answered": 1.0 if answered else 0.0, "deletion_drop": drop, "top_target": 0.0}
    return p


def test_the_three_kinds_of_zero_are_counted_apart(paired):
    rows = [
        make(True, 0.5, True, True),  # answered, non-zero
        make(True, 0.0, True, False),  # answered, zero
        make(False, 0.0, True, True),  # never answered, element in the registry
        make(False, 0.0, False, False),  # never answered, no overlapping element
    ]
    out = paired.zero_fill_split(rows)
    assert out["answered_nonzero"]["pairs"] == 1
    assert out["answered_zero"]["pairs"] == 1
    assert out["unanswered_filled"]["pairs"] == 2
    split = out["unanswered_filled_split"]
    assert split["no_overlapping_registry_element"]["pairs"] == 1
    assert split["element_in_the_registry_but_the_sweep_was_silent_about_this_gene"]["pairs"] == 1


def test_an_answered_zero_is_never_folded_into_the_unanswered_bucket(paired):
    """The project's rule: a value never looked at is not a measured value of zero."""
    out = paired.zero_fill_split([make(True, 0.0, True, False)])
    assert out["answered_zero"]["pairs"] == 1
    assert out["unanswered_filled"]["pairs"] == 0


def test_every_bucket_carries_its_own_positive_count(paired):
    rows = [make(True, 0.5, True, True), make(True, 0.5, True, False)]
    out = paired.zero_fill_split(rows)
    assert out["answered_nonzero"] == {"pairs": 2, "positives": 1, "weighted_positives": 1.0}


def test_the_split_accounts_for_every_pair(paired):
    rows = [
        make(True, 0.5, True, True),
        make(True, 0.0, True, False),
        make(False, 0.0, True, False),
        make(False, 0.0, False, True),
    ]
    out = paired.zero_fill_split(rows)
    assert (
        out["answered_nonzero"]["pairs"] + out["answered_zero"]["pairs"] + out["unanswered_filled"]["pairs"]
    ) == len(rows)


def test_the_reconciliation_compares_coverage_with_coverage_not_with_answered(paired):
    """174 is the headline's COVERAGE count; comparing it against the answered count would be the bug."""
    rows = [make(False, 0.0, False, False)] + [make(False, 0.0, True, False) for _ in range(3)]
    # `coverage_arms_k562_heldout` is a TOP-LEVEL key of crispri_published.json, not nested under
    # `heldout_published_pairs`. Reading it in the wrong place raised KeyError and wrote no result.
    headline = {
        "coverage_arms_k562_heldout": {
            "arm1_all_pairs": {"uncovered_pairs": 1, "gain": 0.1361},
            "arm2_coverage_matched": {"median_gain": 0.1449},
            "arm3_coverage_indicator_control": {"gain": 0.0008},
        }
    }
    out = paired.coverage_reconciliation(rows, headline)
    assert out["ours_uncovered_pairs"] == 1
    assert out["uncovered_agrees"] is True
    assert out["ours_unanswered_pairs"] == 4  # a different count, and not an error
    assert "not a join difference" in out["these_are_two_different_counts"]


def test_the_import_closure_keeps_the_case_the_directory_keeps(paired, tmp_path):
    """On a case-insensitive filesystem a bare is_file would accept pkg/Thing.py for pkg/thing.py."""
    (tmp_path / "pkg").mkdir()
    (tmp_path / "pkg" / "thing.py").write_text("")
    assert paired._is_file_exactly(tmp_path, ["pkg", "thing.py"]) is True
    assert paired._is_file_exactly(tmp_path, ["pkg", "Thing.py"]) is False
    assert paired._is_file_exactly(tmp_path, ["pkg", "absent.py"]) is False


def test_the_computed_closure_lists_only_files_that_are_really_there(paired):
    for rel in paired.counting_path():
        assert (ROOT / rel).is_file(), rel
        parts = rel.split("/")
        assert paired._is_file_exactly(ROOT, parts), rel


def test_the_closure_reaches_this_lanes_module_and_the_crispri_module(paired):
    path = paired.counting_path()
    assert "genomeos/attribution/re2g.py" in path
    assert "genomeos/attribution/crispri.py" in path


def test_the_shared_input_block_names_the_shared_and_the_selection_correlated_input(paired):
    text = " ".join(paired.SHARED_INPUTS.values())
    assert "DNase" in text and "H3K27ac" in text
    assert "shared input" in paired.SHARED_INPUTS["the_honest_statement"]
    assert "not refitted" in text or "Nothing is refitted" in text


def test_every_file_this_lane_wrote_is_claimed_as_its_own(paired):
    """A file of this lane's listed as another session's would make the cleanliness record misleading."""
    for rel in (
        "genomeos/attribution/re2g.py",
        "scripts/re2g_register.py",
        "scripts/re2g_paired.py",
        "scripts/re2g_k562_only_range.py",
        "scripts/re2g_identity.py",
        "scripts/re2g_likeforlike.py",
        "tests/test_re2g.py",
        "tests/test_re2g_paired.py",
    ):
        assert rel in paired.OWN_CODE, rel
        assert (ROOT / rel).is_file(), rel


def test_the_coverage_arms_key_is_read_where_the_committed_file_really_keeps_it():
    """Guards the KeyError this lane hit: the block is top-level, not under heldout_published_pairs."""
    import json

    payload = json.loads((ROOT / "data/results/crispri_published.json").read_text())
    assert "coverage_arms_k562_heldout" in payload
    assert "coverage_arms_k562_heldout" not in payload["heldout_published_pairs"]
    assert payload["coverage_arms_k562_heldout"]["arm1_all_pairs"]["uncovered_pairs"] == 174
