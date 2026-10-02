# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fourth increment's refusal, and the re-distillation's registered frame.

These tests run without the git-ignored local data: they exercise the rules that make the refusal
mean something -- the reconciliation that refuses a payload, the anchor check that refuses a payload
whose figures are not comparable with the committed count, the frame that cannot be streamed
against unregistered, and the distillation that cannot be written over the pinned one.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name: str):
    sys.path.insert(0, str(ROOT))
    sys.path.insert(0, str(ROOT / "scripts"))
    spec = importlib.util.spec_from_file_location(f"_{name}", ROOT / f"scripts/{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def i4():
    return _load("response_map_increment4_count")


@pytest.fixture(scope="module")
def frame_mod():
    return _load("eqtl_crispri_frame")


# --- the rule that refuses three of the four candidates ----------------------------------------


def test_the_native_locus_rule_is_increment_3s_own_words(i4):
    """Carried, not paraphrased: paraphrasing it is how it would be relaxed."""
    c3 = _load("response_map_increment3_count")
    assert i4.NATIVE_LOCUS_RULE == c3.NATIVE_LOCUS_RULE
    assert "reporter assay reads a copy of the sequence outside the locus" in i4.NATIVE_LOCUS_RULE
    assert i4.REPORTER_WORDING == (
        "they read a sequence, or the bases inside it, not the native locus after a perturbation of it"
    )


def test_the_committed_count_still_says_what_this_count_starts_from(i4):
    """If increment 3's own figures move, this refusal's starting point has gone with them."""
    p = ROOT / i4.INCREMENT_3_COUNT
    if not p.exists():  # the results are git-ignored in a clean worktree
        pytest.skip("increment 3's count is not in this checkout")
    c = json.loads(p.read_text())
    k = {x["kind"]: x for x in c["kinds_of_evidence"]}
    assert c["denominator"]["addable_loci_re_derived_here"] == 110
    assert c["denominator"]["independent_loci_of_the_measured_layer"] == 473
    assert k["eqtl"]["loci_on_the_existing_110"] == 6
    assert k["eqtl"]["further"]["further_independent_loci"] == 264
    for kind in ("lentimpra", "vista", "satmut"):
        assert k[kind]["native_locus_observations_on_the_existing_110"] == 0
        assert k[kind]["an_observation_of_the_native_locus"] is False


def test_the_withdrawn_number_is_carried_with_both_figures(i4):
    """The point of the field is that a reader can see the overstatement, not be told of it."""
    w = i4.WITHDRAWN
    assert w["the_number"] == 92
    assert w["the_honest_figure"] == 6
    assert w["the_number"] // w["the_honest_figure"] == 15
    assert "proximity arm" in w["why_it_is_wrong"]
    assert "198" in w["the_shape_it_is"]


def test_the_pointer_field_trap_is_in_the_result_not_only_the_docstring(i4):
    assert "elements_where" in i4.POINTER_FIELD_NOTE
    assert "cache_target_effect" in i4.POINTER_FIELD_NOTE
    assert "17 published results" in i4.POINTER_FIELD_NOTE


# --- the reconciliation, shown refusing ---------------------------------------------------------


def _arms_stub(deep_cause="adds_no_locus_under_the_rules_as_registered"):
    """The smallest shape `reconcile` reads, with two loci populations that balance."""
    return {
        "mine": [0, 0, 1],
        "best": {0: "all_three_present", 1: "no_reader_state"},
        "eqtl_groups": {7, 8, 9},
        "eqtl_reached": {7},
        "eqtl_further": {8, 9},
        "deepening": {"refused_under_cause": deep_cause},
        "widening": {"refused_under_cause": "one_native_locus_kind_and_no_possible_chain"},
    }


def _blocks_stub():
    return [
        {"candidate": "lentimpra", "refused_under_cause": "zero_native_locus_observations"},
        {"candidate": "vista", "refused_under_cause": "zero_native_locus_observations"},
        {"candidate": "satmut", "refused_under_cause": "zero_native_locus_observations"},
        {"candidate": "eqtl"},
    ]


def test_reconciliation_holds_when_every_arm_has_a_cause(i4):
    rec = i4.reconcile(_arms_stub(), _blocks_stub())
    assert rec["arms_they_divide_into"] == 5
    assert rec["arms_accounted_for"] == 5
    assert rec["arms_built"] == 0
    assert rec["holds"] is True


def test_reconciliation_refuses_when_an_arm_loses_its_cause(i4):
    """An arm may never leave the count unnamed, and the payload is not served when one does."""
    blocks = _blocks_stub()
    blocks.pop()  # the eQTL candidate, and with it both its arms
    rec = i4.reconcile(_arms_stub(), blocks)
    assert rec["arms_accounted_for"] == 3
    assert rec["holds"] is False


def test_reconciliation_refuses_when_the_eqtl_loci_do_not_add_up(i4):
    a = _arms_stub()
    a["eqtl_further"] = {8}  # one locus has gone missing between the two arms
    rec = i4.reconcile(a, _blocks_stub())
    assert rec["eqtl_loci_accounted_for"] == 2
    assert rec["eqtl_loci_of_this_kind"] == 3
    assert rec["holds"] is False


def test_main_raises_rather_than_saving_a_payload_that_does_not_reconcile(i4, monkeypatch):
    """The refusal is in the writer, not in a reviewer's attention."""
    monkeypatch.setattr(i4, "arms", lambda *a, **k: _arms_stub())
    monkeypatch.setattr(i4, "candidate_blocks", lambda a: _blocks_stub()[:2])

    def _never(*a, **k):
        raise AssertionError("save_result must not be reached")

    monkeypatch.setattr(i4, "save_result", _never)
    with pytest.raises(i4.RefusedError, match="reconciliation does not hold"):
        i4.main()


def test_main_raises_when_the_anchors_do_not_reproduce(i4, monkeypatch):
    """A figure not comparable with the committed count may not be published beside it."""
    monkeypatch.setattr(i4, "arms", lambda *a, **k: _arms_stub())
    monkeypatch.setattr(i4, "candidate_blocks", lambda a: _blocks_stub())
    monkeypatch.setattr(
        i4,
        "anchors",
        lambda a, root=ROOT: {"all_reproduce": False, "anchors_differing": {"pairs": [1, 2]}},
    )

    def _never(*a, **k):
        raise AssertionError("save_result must not be reached")

    monkeypatch.setattr(i4, "save_result", _never)
    with pytest.raises(i4.RefusedError, match="anchors do not reproduce"):
        i4.main()


# --- the re-distillation's frame ----------------------------------------------------------------


def test_the_margin_is_imported_and_not_chosen(frame_mod):
    eqtl_targets = _load("eqtl_targets")
    assert frame_mod.MARGIN == eqtl_targets.MARGIN == 500


def test_the_floor_is_a_number_and_above_what_increment_3_built(frame_mod):
    assert frame_mod.FLOOR_ADDABLE_LOCI_ON_A_SHOWN_ELEMENT == 20
    assert frame_mod.FLOOR_ADDABLE_LOCI_ON_A_SHOWN_ELEMENT > 6


def test_the_hits_do_not_go_where_the_pinned_ones_are(frame_mod):
    from genomeos.attribution import eqtl

    assert frame_mod.KNOWLEDGE != eqtl.KNOWLEDGE
    assert frame_mod.KNOWLEDGE.as_posix() == "data/knowledge/gtex_crispri"


def test_intervals_carry_the_margin_on_each_side(frame_mod):
    iv = frame_mod.intervals({"E1": ("chr1", 10_000, 10_200)})
    assert iv.n == 1
    # `Intervals.at` tests `start < pos <= end` on the widened interval (9_500, 10_700], which is
    # a 1-based position against a half-open 0-based interval. The boundary is asserted both ways
    # rather than described, because a margin off by one would still have produced a plausible count.
    assert iv.at("chr1", 9_600) == ["chr1:10000-10200"]  # inside the left margin
    assert iv.at("chr1", 9_501) == ["chr1:10000-10200"]  # its first included position
    assert iv.at("chr1", 9_500) == []  # one before it
    assert iv.at("chr1", 10_700) == ["chr1:10000-10200"]  # its last included position
    assert iv.at("chr1", 10_701) == []  # one past it
    assert iv.at("chr1", 9_400) == []


def test_two_pairs_on_one_element_index_one_interval(frame_mod):
    """The frame is of elements, so an element many pairs attach to is indexed once."""
    iv = frame_mod.intervals(
        {"E1": ("chr1", 10_000, 10_200), "E2": ("chr1", 10_000, 10_200), "E3": ("chr2", 5, 50)}
    )
    assert iv.n == 2


def test_distil_refuses_without_a_registration(frame_mod, monkeypatch, tmp_path, capsys):
    """Nothing is streamed against an unregistered frame."""
    monkeypatch.setattr(frame_mod, "ROOT", tmp_path)
    monkeypatch.setattr(frame_mod, "frame", lambda *a, **k: {"att_iv": {"E1": ("chr1", 10_000, 10_200)}})
    monkeypatch.setattr(sys, "argv", ["eqtl_crispri_frame.py", "--distil"])

    def _never(*a, **k):
        raise AssertionError("eqtl.distil must not be reached")

    monkeypatch.setattr(frame_mod.eqtl, "distil", _never)
    with pytest.raises(frame_mod.RefusedError, match="registered before"):
        frame_mod.main()


def test_distil_refuses_when_the_frame_has_moved_since_registration(frame_mod, monkeypatch, tmp_path):
    """A stream against a frame that is not the registered one proves nothing about it."""
    (tmp_path / "data/results").mkdir(parents=True)
    (tmp_path / "data/results/eqtl_crispri_frame_registration.json").write_text(
        json.dumps({"frame": {"intervals_indexed": 99}})
    )
    monkeypatch.setattr(frame_mod, "ROOT", tmp_path)
    monkeypatch.setattr(frame_mod, "frame", lambda *a, **k: {"att_iv": {"E1": ("chr1", 10_000, 10_200)}})
    monkeypatch.setattr(sys, "argv", ["eqtl_crispri_frame.py", "--distil"])

    def _never(*a, **k):
        raise AssertionError("eqtl.distil must not be reached")

    monkeypatch.setattr(frame_mod.eqtl, "distil", _never)
    with pytest.raises(frame_mod.RefusedError, match="frame has moved"):
        frame_mod.main()


def test_distil_refuses_to_overwrite_the_pinned_distillation(frame_mod, monkeypatch, tmp_path):
    """data/knowledge/gtex is a sha256-pinned input of four committed results."""
    (tmp_path / "data/results").mkdir(parents=True)
    (tmp_path / "data/results/eqtl_crispri_frame_registration.json").write_text(
        json.dumps({"frame": {"intervals_indexed": 1}})
    )
    monkeypatch.setattr(frame_mod, "ROOT", tmp_path)
    monkeypatch.setattr(frame_mod, "frame", lambda *a, **k: {"att_iv": {"E1": ("chr1", 10_000, 10_200)}})
    monkeypatch.setattr(frame_mod, "KNOWLEDGE", frame_mod.eqtl.KNOWLEDGE)
    monkeypatch.setattr(sys, "argv", ["eqtl_crispri_frame.py", "--distil"])

    def _never(*a, **k):
        raise AssertionError("eqtl.distil must not be reached")

    monkeypatch.setattr(frame_mod.eqtl, "distil", _never)
    with pytest.raises(frame_mod.RefusedError, match="sha256-pinned"):
        frame_mod.main()


def test_count_refuses_before_anything_has_been_distilled(frame_mod, monkeypatch, tmp_path):
    monkeypatch.setattr(frame_mod, "KNOWLEDGE", tmp_path / "empty")
    (tmp_path / "empty").mkdir()
    monkeypatch.setattr(frame_mod, "frame", lambda *a, **k: {"att_iv": {"E1": ("chr1", 10_000, 10_200)}})
    monkeypatch.setattr(sys, "argv", ["eqtl_crispri_frame.py", "--count"])
    with pytest.raises(frame_mod.RefusedError, match="run --distil first"):
        frame_mod.main()


def test_assessed_and_found_are_different_fields(frame_mod):
    """Conflating them is how an unassessed element becomes an element with no eQTL."""
    assert "BY CONSTRUCTION" in frame_mod.ASSESSED_BY_CONSTRUCTION
    assert "may not be reported as one" in frame_mod.ASSESSED_BY_CONSTRUCTION


def test_the_retained_hit_rule_sets_no_cutoff_of_its_own(frame_mod):
    assert "GTEx's significance call is the one used" in frame_mod.RETAINED_HIT_RULE
    assert "sets no cut-off of its own" in frame_mod.RETAINED_HIT_RULE


def test_retained_rows_and_attachments_are_separate_units(frame_mod):
    """A row inside two overlapping elements attaches twice, so rows and attachments differ.

    The first version of this block subtracted a row count from an attachment count and published a
    margin-only figure of -1,220. Only rows may be subtracted from rows.
    """
    src = (ROOT / "scripts/eqtl_crispri_frame.py").read_text()
    assert '"attachments_element_by_record"' in src
    assert '"distinct_rows_attached_inside_an_element"' in src
    assert '"distinct_rows_retained_by_the_distillation"' in src
    assert '"distinct_rows_retained_in_the_margin_only"' in src
    # the window is one window, so the margin-only figure is 0 by construction and says so
    assert "BY " in src and "CONSTRUCTION" in src
    assert "not a finding" in src
    # the defect is named in the file, so the next reader does not have to rediscover why
    assert "-1,220" in src


def test_a_row_in_two_overlapping_elements_counts_once_as_a_row_and_twice_as_an_attachment():
    """The arithmetic the -1,220 came from, in isolation, with the two units kept apart."""
    row = {"tissue": "Lung", "pos": 1_000, "gene_id": "ENSG1", "slope": -0.4}
    hits_by_chrom = {"chr1": [row]}
    # two overlapping compiled elements both contain position 1,000
    with_hit = {"E1": [row], "E2": [row]}
    att_iv = {"E1": ("chr1", 900, 1_100), "E2": ("chr1", 950, 1_050)}

    def key(chrom, h):
        return (h["tissue"], chrom, h["pos"], h["gene_id"])

    attachments = sum(len(v) for v in with_hit.values())
    retained_rows = {key(c, h) for c, v in hits_by_chrom.items() for h in v}
    attached_rows = {key(att_iv[e][0], h) for e, v in with_hit.items() for h in v}
    assert attachments == 2
    assert len(retained_rows) == len(attached_rows) == 1
    assert len(retained_rows) - len(attached_rows) == 0  # and never negative
    assert retained_rows - attached_rows == set()


def test_the_attachment_window_is_the_registered_one_and_not_the_bare_element(i4, frame_mod):
    """`hits_in` widens by the distil margin itself, so the window is the element plus 500 each side.

    An earlier draft of both writers said the attachment was "the element ITSELF". It is not, and
    the margin-only figure of 0 in the re-distillation is that identity rather than a finding.
    """
    c3 = _load("response_map_increment3_count")
    rows = [{"pos": 1_000, "tissue": "Lung", "gene_id": "G", "slope": 0.1, "pval": 1e-9}]
    # a variant 400 bases outside the element still attaches, because the rule carries the margin
    assert c3.hits_in(rows, 1_400, 1_500) == rows
    assert c3.hits_in(rows, 1_501, 1_600) == []  # 501 out: past the margin
    assert c3.EQTL_MARGIN == frame_mod.MARGIN == 500
    assert "widened by 500 bases each side" in i4.WITHDRAWN["the_rule_that_gives_6"]
    assert "not about a tolerance" in i4.WITHDRAWN["the_window_is_not_what_separates_92_from_6"]
