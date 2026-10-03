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
import tracked_paths as tp

from genomeos import manifest as mf
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
    assert mf._is_file_exactly(tmp_path, ["pkg", "thing.py"]) is True
    assert mf._is_file_exactly(tmp_path, ["pkg", "Thing.py"]) is False
    assert mf._is_file_exactly(tmp_path, ["pkg", "absent.py"]) is False


def test_the_computed_closure_lists_only_files_that_are_really_there(paired):
    for rel in mf.counting_path(paired.ENTRY, paired.ROOT):
        assert (ROOT / rel).is_file(), rel
        parts = rel.split("/")
        assert mf._is_file_exactly(ROOT, parts), rel


def test_the_closure_reaches_this_lanes_module_and_the_crispri_module(paired):
    path = mf.counting_path(paired.ENTRY, paired.ROOT)
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


def test_the_paired_result_leads_with_the_notice_that_withdraws_its_own_positive(paired):
    """The first key a reader meets must be the one pointing at the binding, negative comparison."""
    import json

    path = ROOT / "data/results/re2g_paired.json"
    # TRACKED, so it HAS been generated in every checkout of this commit; the skip could never fire
    tp.must_be_present(path, was="re2g_paired.json has not been generated in this checkout")
    payload = json.loads(path.read_text())
    assert next(iter(payload)) in ("result", "date", "read_this_first")
    notice = payload["read_this_first"]
    assert "MAY NOT SAY THE DELETION MODEL RANKS BETTER" in notice
    assert "MAY NOT BE CITED ALONE" in notice
    assert "re2g_likeforlike.json" in notice
    # the figure itself is not copied in, so the two results cannot drift apart
    assert "0.0593" not in notice


def test_the_like_for_like_result_leads_with_its_own_binding_outcome():
    import json

    path = ROOT / "data/results/re2g_likeforlike.json"
    # TRACKED: see the note above; a missing committed result is a broken tree, not an ungenerated one
    tp.must_be_present(path, was="re2g_likeforlike.json has not been generated in this checkout")
    payload = json.loads(path.read_text())
    notice = payload["read_this_first"]
    assert "THE BINDING RESULT OF THIS LANE" in notice
    assert payload["paired_delta"]["primary_k562"]["reading"] in notice


# --- no registered reading word on an unregistered population ---------------------------------------

READING_WORDS = ("ranks better than ENCODE-rE2G on these pairs", "no difference detected", "ranks worse")


def _likeforlike_module():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "re2g_likeforlike", ROOT / "scripts" / "re2g_likeforlike.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _results_with_registered_sets():
    """Each committed result of this lane with the registered population set taken FROM CODE.

    The set is not read out of the result file: a result that simply omits the field would otherwise
    excuse itself, which is exactly how the breach this guards against survived one review.
    """
    import json

    from genomeos.attribution import re2g

    registered = {
        # re2g.POPULATIONS registers all three, each with its own `why`.
        "re2g_paired": ("primary_k562", "secondary_pooled", "gm12878"),
        # re2g.SECOND_REGISTRATION registers the K562 primary and nothing else.
        "re2g_likeforlike": tuple(_likeforlike_module().REGISTERED_POPULATIONS),
    }
    assert set(re2g.POPULATIONS) == {"primary", "secondary", "gm12878"}
    out = []
    for name, names in registered.items():
        path = ROOT / "data/results" / f"{name}.json"
        if not path.is_file():
            continue
        payload = json.loads(path.read_text())
        if "paired_delta" in payload:
            out.append((name, payload, names))
    return out


def test_only_a_registered_population_carries_a_reading_word():
    """General over whatever populations a result has, so a new descriptive arm cannot reintroduce this.

    A descriptive arm carrying "ranks better than ENCODE-rE2G on these pairs" gets lifted as a finding.
    It happened once in re2g_likeforlike.json, three lines below a primary reading "no difference
    detected", and this is what fails if it happens again. The registered set comes from code, so
    omitting the field in the result is not a way to pass.
    """
    results = _results_with_registered_sets()
    if not results:
        pytest.skip("no result of this lane has been generated in this checkout")
    for name, payload, registered in results:
        for population, block in payload["paired_delta"].items():
            if not isinstance(block, dict) or "delta_auprc" not in block:
                continue
            if population in registered:
                continue
            assert block.get("reading") is None, (
                f"{name}: descriptive population {population!r} carries a reading word "
                f"{block.get('reading')!r}"
            )
            assert block.get("descriptive"), f"{name}: {population!r} is not labelled descriptive"
            arm = block.get("dnase_plus_distance_alone") or {}
            assert arm.get("reading") is None, (
                f"{name}: descriptive population {population!r} has an arm carrying a reading word"
            )


def test_the_top_level_reading_block_holds_the_registered_population_only():
    for name, payload, registered in _results_with_registered_sets():
        assert set(payload["reading"]) <= set(registered), name
        for population in payload["paired_delta"]:
            if population not in registered:
                assert population not in payload["reading"], (
                    f"{name}: descriptive {population!r} is in the reading block"
                )


def test_a_descriptive_arm_keeps_its_estimate_and_interval_and_loses_only_the_word():
    """The fix is a labelling fix: the figures stay so a reader can still see them."""
    for name, payload, registered in _results_with_registered_sets():
        for population, block in payload["paired_delta"].items():
            if population in registered or not isinstance(block, dict):
                continue
            if "delta_auprc" not in block:
                continue
            assert block["delta_auprc"] is not None, f"{name}: {population!r} lost its point estimate"
            assert "ci95" in block, f"{name}: {population!r} lost its interval field"


def test_the_baseline_arm_on_the_registered_population_is_labelled_pre_specified_not_post_hoc():
    """It is in the registration code with its own reading before any figure of it existed."""
    import json

    path = ROOT / "data/results/re2g_likeforlike.json"
    # TRACKED: see the note above; a missing committed result is a broken tree, not an ungenerated one
    tp.must_be_present(path, was="re2g_likeforlike.json has not been generated in this checkout")
    payload = json.loads(path.read_text())
    arm = payload["paired_delta"]["primary_k562"]["dnase_plus_distance_alone"]
    assert "pre-specified secondary arm" in arm["label"]
    assert "post hoc" not in arm["label"].lower()
    assert arm["reading"] in READING_WORDS
