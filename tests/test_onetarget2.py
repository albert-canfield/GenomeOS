# SPDX-License-Identifier: AGPL-3.0-or-later
"""The remaining one-target consumers' registration: its guards, each planted to REFUSE.

Every guard in `genomeos/attribution/onetarget2.py` is exercised twice here: once on the real census
(it passes), and once on a planted census or a planted tree that it must REFUSE. The plants are
built in memory or in tmp_path, so no consumer module is touched by this file.

Mutation results are recorded in each test's docstring: removing the guard's own raise makes the
named test fail. Twelve guards were found in this project that could never fire; none of these five
is one of them, and the drift guard caught this registration's own file on its first run.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from genomeos.attribution import onetarget2 as ot

# ---------------------------------------------------------------------------
# the real census
# ---------------------------------------------------------------------------


def test_the_real_census_passes_every_guard():
    ot.check_all()


def test_every_identity_module_names_its_own_falsifier():
    """A falsifier per module is the row's own wording, and it is what makes the wave readable."""
    for c in ot.CENSUS:
        assert c["falsifier"].strip(), c["module"]
        assert c["class"] in ot.CLASSES, c["module"]
        assert c["status"] in ot.STATUSES, c["module"]
        assert c["claims_today"].strip(), c["module"]
        assert c["reads"].strip(), c["module"]
    assert len({c["falsifier"] for c in ot.CENSUS}) == len(ot.CENSUS)


def test_no_module_is_in_two_lists():
    """A module classified twice would be counted twice and reasoned about once."""
    lists = {
        "census": {c["module"] for c in ot.CENSUS},
        "record_only": {m["module"] for m in ot.RECORD_ONLY},
        "already_moved": {m["module"] for m in ot.ALREADY_MOVED},
        "excluded_by_name": {m["module"] for m in ot.EXCLUDED_BY_NAME},
    }
    names = list(lists)
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            assert not (lists[a] & lists[b]), f"{a} and {b} share {sorted(lists[a] & lists[b])}"


def test_the_invariant_is_carried_with_its_committed_figures_and_its_own_limit():
    """The 2026-09-27 reading is quoted, not strengthened: its own exception is carried too."""
    assert "849,469" in ot.THE_INVARIANT
    assert "0 head disagreements" in ot.THE_INVARIANT
    assert "onetarget_consumers.json" in ot.THE_INVARIANT
    assert "not a theorem" in ot.THE_INVARIANT
    assert "the coding head is not" in ot.THE_INVARIANT


def test_the_two_invariant_modules_are_registered_as_not_to_be_moved():
    """Moving a module that cannot change would move a number with no reason to move it."""
    inv = ot.by_class("invariant")
    assert inv, "the finding that two consumers are invariant is part of this registration"
    for c in ot.CENSUS:
        if c["class"] == "invariant":
            assert "NOT to be moved" in c["falsifier"], c["module"]


def test_both_closure_consumers_keep_their_falsifier_registered_and_unfired():
    for c in ot.CENSUS:
        if c["status"] == "frozen_closure":
            assert c["module"] in ot.CLOSURE_FROZEN
            assert "NOT RUN THIS WAVE" in c["falsifier"], c["module"]


# ---------------------------------------------------------------------------
# guard 1: the census may not drift from the tree
# ---------------------------------------------------------------------------


def _tree(tmp: Path, rel: str, body: str) -> Path:
    p = tmp / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(body)
    return p


HEAD_LINE = '    pc = e.get("predicted' + '_coding") or {}\n'


def test_the_scan_finds_a_head_read_and_ignores_a_mere_mention(tmp_path):
    _tree(tmp_path, "genomeos/a.py", "x = 1\n" + HEAD_LINE)
    _tree(tmp_path, "genomeos/b.py", 'KEYS = ("predicted", "predicted' + '_coding")\n')
    found = ot.head_read_sites(tmp_path)
    assert found == {"genomeos/a.py": [2]}


def test_PLANTED_an_undeclared_head_reading_module_is_REFUSED(tmp_path):
    """Mutation: delete the `unknown` raise in check_census_matches_tree and this test fails."""
    _tree(tmp_path, "genomeos/brand_new_consumer.py", HEAD_LINE)
    with pytest.raises(ot.CensusDriftError) as e:
        ot.check_census_matches_tree(tmp_path)
    assert "brand_new_consumer.py" in str(e.value)
    assert "in no list" in str(e.value)


def _faithful_tree(tmp: Path, shift: dict[str, int] | None = None) -> None:
    """A synthetic tree whose head reads sit on exactly the lines the census declares."""
    shift = shift or {}
    for c in ot.CENSUS:
        want = sorted(set(c["lines"]) | set(c.get("not_head_lines") or {}))
        want = [n + shift.get(c["module"], 0) for n in want]
        lines = ["\n"] * (max(want) + 1)
        for n in want:
            lines[n - 1] = HEAD_LINE
        _tree(tmp, c["module"], "".join(lines))


@pytest.mark.xfail(
    strict=True,
    reason="SUPERSEDED by the moves of 2026-10-03, kept rather than edited (precedent dbb5d4a). "
    "`_faithful_tree` builds its synthetic tree from CENSUS's REGISTERED lines, and three modules "
    "have since moved, so the tree it builds no longer describes them and the drift guard refuses "
    "it -- correctly. The live form of this check is "
    "`test_the_tree_built_from_the_EFFECTIVE_lines_is_accepted` below, which builds from "
    "`effective_lines()`. This one is left in place, failing strictly, so that the supersession is "
    "visible and so that a future change which made the registered lines describe the tree again "
    "would turn this red.",
)
def test_the_faithful_tree_is_accepted(tmp_path):
    """The plant below differs from this by one line, so the refusal is about the line and nothing else."""
    _faithful_tree(tmp_path)
    ot.check_census_matches_tree(tmp_path)


def test_PLANTED_a_declared_module_whose_line_moved_is_REFUSED(tmp_path):
    """A line that moved may be a peer's edit; a line that vanished may be a dropped claim.

    Mutation: delete the line-set raise in check_census_matches_tree and this test fails.
    """
    mod = ot.CENSUS[0]["module"]
    _faithful_tree(tmp_path, {mod: 1})
    with pytest.raises(ot.CensusDriftError) as e:
        ot.check_census_matches_tree(tmp_path)
    assert mod in str(e.value)
    assert "accounts for lines" in str(e.value)


def test_PLANTED_a_declared_module_that_reads_no_head_at_all_is_REFUSED(tmp_path):
    """Mutation: delete the `mod not in found` raise and this test fails."""
    for m in ot.CENSUS:
        _tree(tmp_path, m["module"], "x = 1\n")
    with pytest.raises(ot.CensusDriftError) as e:
        ot.check_census_matches_tree(tmp_path)
    assert "reads no compact head" in str(e.value)


@pytest.mark.xfail(
    strict=True,
    reason="SUPERSEDED by the moves of 2026-10-03, kept rather than edited (precedent dbb5d4a). "
    "This plant claims a REGISTERED line twice, and for a module that has since moved the "
    "registered line is no longer one the tree matches, so the overlap clause is not the clause "
    "that fires. A plant against a stale rule is worse than no plant: it would pass while testing "
    "nothing. The live form is "
    "`test_PLANTED_an_EFFECTIVE_line_declared_both_as_a_head_and_as_not_one_is_REFUSED` below.",
)
def test_PLANTED_a_line_declared_both_as_a_head_and_as_not_one_is_REFUSED(monkeypatch):
    """Mutation: delete the overlap raise in check_census_matches_tree and this test fails."""
    bad = tuple(
        {**c, "not_head_lines": {c["lines"][0]: "claimed twice"}} if i == 0 else c
        for i, c in enumerate(ot.CENSUS)
    )
    monkeypatch.setattr(ot, "CENSUS", bad)
    with pytest.raises(ot.CensusDriftError) as e:
        ot.check_census_matches_tree()
    assert "both as a head read and as not one" in str(e.value)


# ---------------------------------------------------------------------------
# guard 2: one falsifier per module, never one for the set
# ---------------------------------------------------------------------------


def test_PLANTED_two_modules_sharing_one_falsifier_is_REFUSED(monkeypatch):
    """This is the row's own point: a shared falsifier lets one module pass on another's evidence.

    Mutation: delete the `f in seen` raise in check_one_falsifier_per_module and this test fails.
    """
    shared = "the consumers' figures move"
    bad = tuple({**c, "falsifier": shared} for c in ot.CENSUS)
    monkeypatch.setattr(ot, "CENSUS", bad)
    with pytest.raises(ot.SharedFalsifierError) as e:
        ot.check_one_falsifier_per_module()
    assert "one falsifier verbatim" in str(e.value)


def test_PLANTED_a_module_with_no_falsifier_is_REFUSED(monkeypatch):
    """Mutation: delete the empty-falsifier raise and this test fails."""
    bad = tuple({**c, "falsifier": "   "} if i == 3 else c for i, c in enumerate(ot.CENSUS))
    monkeypatch.setattr(ot, "CENSUS", bad)
    with pytest.raises(ot.SharedFalsifierError) as e:
        ot.check_one_falsifier_per_module()
    assert "has no falsifier" in str(e.value)


# ---------------------------------------------------------------------------
# guard 3: no closure file may be marked movable
# ---------------------------------------------------------------------------


def test_PLANTED_a_closure_file_marked_free_is_REFUSED(monkeypatch):
    """sender_closure() hashes the working tree, so even an uncommitted edit refuses the send.

    Mutation: delete the raise in check_no_closure_file_is_moving and this test fails.
    """
    bad = tuple({**c, "status": "free"} if c["module"] in ot.CLOSURE_FROZEN else c for c in ot.CENSUS)
    monkeypatch.setattr(ot, "CENSUS", bad)
    with pytest.raises(ot.ClosureFileMovingError) as e:
        ot.check_no_closure_file_is_moving()
    assert "import closure" in str(e.value)


def test_the_declared_closure_files_really_are_in_the_senders_closure():
    """A frozen list that named the wrong files would protect nothing."""
    from genomeos.attribution import astrorun

    closure = set(astrorun.sender_closure())
    for m in ot.CLOSURE_FROZEN:
        assert m in closure, m
    for c in ot.CENSUS:
        if c["status"] != "frozen_closure":
            assert c["module"] not in closure, c["module"]


# ---------------------------------------------------------------------------
# guard 4: no preview
# ---------------------------------------------------------------------------


def test_PLANTED_a_field_a_run_would_fill_is_REFUSED():
    """Mutation: delete the FORBIDDEN_FIELDS raise in check_no_preview and this test fails."""
    with pytest.raises(ot.PreviewInRegistrationError) as e:
        ot.check_no_preview({"result": "x", "census": [{"module": "m", "verdict": "moved"}]})
    assert "a field a run would fill" in str(e.value)


def test_PLANTED_a_measurement_in_a_registration_is_REFUSED():
    """Mutation: delete the float raise in check_no_preview and this test fails."""
    with pytest.raises(ot.PreviewInRegistrationError) as e:
        ot.check_no_preview({"result": "x", "scope_counts": {"identity": 18, "share": 0.62}})
    assert "is a float" in str(e.value)


def test_a_top_level_result_stamp_is_not_a_preview():
    """The registry writes `result` holding this registration's own NAME, which is not an outcome."""
    ot.check_no_preview({"result": "onetarget2_registration", "date": "2026-10-03"})


def test_the_registered_payload_holds_no_preview_and_no_measurement():
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "onetarget2_register", Path(__file__).resolve().parents[1] / "scripts/onetarget2_register.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    p = mod.payload()
    ot.check_no_preview(p)
    assert p["what_the_run_will_write"] == ot.WILL_WRITE
    for forbidden in ("results", "figures", "findings", "measured_figures"):
        assert forbidden not in p


# ---------------------------------------------------------------------------
# guard 5: the registration may not misstate its own scope
# ---------------------------------------------------------------------------


def test_PLANTED_a_scope_figure_that_has_drifted_is_REFUSED(monkeypatch):
    """The prose quotes these figures, so a drift is a registration that misstates its own scope.

    Mutation: delete the raise in check_scope_counts and this test fails.
    """
    monkeypatch.setattr(ot, "CENSUS", ot.CENSUS[:-1])
    with pytest.raises(ot.CensusDriftError) as e:
        ot.check_scope_counts()
    assert "scope counts" in str(e.value)


def test_the_prose_figure_and_the_declared_figure_are_the_same_number():
    assert ot.SCOPE_COUNTS["identity"] == 18
    assert "18 modules whose own" in ot.ROW_SAYS_FIVE_THE_TREE_SAYS_OTHERWISE
    assert f"leaving {ot.SCOPE_COUNTS['identity_free_today']} " in (ot.ROW_SAYS_FIVE_THE_TREE_SAYS_OTHERWISE)
    assert f"holds {ot.SCOPE_COUNTS['scripts_out_of_scope']} further" in (
        ot.ROW_SAYS_FIVE_THE_TREE_SAYS_OTHERWISE
    )


# ---------------------------------------------------------------------------
# guard 6: the one piece of prior evidence is checked against the file
# ---------------------------------------------------------------------------


def test_the_quoted_invariant_is_what_the_committed_result_holds():
    hc = ot.check_the_quoted_invariant()
    assert hc["compared"] == ot.PRIOR_EVIDENCE_COMPARED
    assert hc["head_disagreements"] == ot.PRIOR_EVIDENCE_DISAGREEMENTS


def test_PLANTED_a_misquoted_prior_figure_is_REFUSED(tmp_path):
    """Mutation: delete the `got != want` raise in check_the_quoted_invariant and this test fails."""
    import json

    p = tmp_path / "onetarget_consumers.json"
    p.write_text(
        json.dumps(
            {
                ot.PRIOR_EVIDENCE_ARM: {
                    "head_control": {
                        "compared": ot.PRIOR_EVIDENCE_COMPARED,
                        "head_disagreements": 7,
                        "coding_head_disagreements": 0,
                    }
                }
            }
        )
    )
    with pytest.raises(ot.MisquotedPriorEvidenceError) as e:
        ot.check_the_quoted_invariant(p)
    assert "Re-read the committed result" in str(e.value)


def test_PLANTED_a_missing_prior_evidence_file_is_REFUSED(tmp_path):
    """Mutation: delete the `not p.exists()` raise and this test fails."""
    with pytest.raises(ot.MisquotedPriorEvidenceError) as e:
        ot.check_the_quoted_invariant(tmp_path / "absent.json")
    assert "not on this machine" in str(e.value)


def test_PLANTED_an_invariant_that_states_no_figure_is_REFUSED(monkeypatch):
    """Mutation: delete the THE_INVARIANT spelling raise and this test fails."""
    monkeypatch.setattr(ot, "THE_INVARIANT", "the heads agree")
    with pytest.raises(ot.MisquotedPriorEvidenceError) as e:
        ot.check_the_quoted_invariant()
    assert "does not state the" in str(e.value)


def test_the_registration_says_what_it_cannot_establish_and_authorises_nothing():
    assert len(ot.CANNOT_ESTABLISH) >= 5
    assert "authorises no conclusion" in ot.AUTHORISES_NO_CONCLUSION
    assert "DO NOT EXIST" in ot.AUTHORISES_NO_CONCLUSION
    assert len(ot.REFUSALS) == 6


# ---------------------------------------------------------------------------
# added 2026-10-03 with the first three moves. Nothing above is edited: the two tests the moves
# superseded are left in place as strict expected failures with their reason beside them, and these
# are their live forms, built from `effective_lines()` instead of from the registered lines.
# ---------------------------------------------------------------------------


def _effective_tree(tmp: Path, shift: dict[str, int] | None = None) -> None:
    """A synthetic tree whose head reads sit on the lines the census EXPECTS TODAY."""
    shift = shift or {}
    eff = ot.effective_lines()
    for c in ot.CENSUS:
        mod = c["module"]
        want = sorted(set(eff[mod]) | set(c.get("not_head_lines") or {}))
        want = [n + shift.get(mod, 0) for n in want]
        lines = ["\n"] * (max(want) + 1)
        for n in want:
            lines[n - 1] = HEAD_LINE
        _tree(tmp, mod, "".join(lines))


def test_the_tree_built_from_the_EFFECTIVE_lines_is_accepted(tmp_path):
    """The live form of the superseded faithful-tree check; the plant below differs by one line."""
    _effective_tree(tmp_path)
    ot.check_census_matches_tree(tmp_path)


def test_PLANTED_an_EFFECTIVE_line_that_moved_is_REFUSED(tmp_path):
    """Mutation: delete the line-set raise in check_census_matches_tree and this test fails."""
    mod = ot.CENSUS[0]["module"]
    _effective_tree(tmp_path, {mod: 1})
    with pytest.raises(ot.CensusDriftError) as e:
        ot.check_census_matches_tree(tmp_path)
    assert mod in str(e.value) and "accounts for lines" in str(e.value)


def test_PLANTED_an_EFFECTIVE_line_declared_both_as_a_head_and_as_not_one_is_REFUSED(monkeypatch):
    """Mutation: delete the overlap raise in check_census_matches_tree and this test fails."""
    eff = ot.effective_lines()
    bad = tuple(
        {**c, "not_head_lines": {eff[c["module"]][0]: "claimed twice"}} if i == 0 else c
        for i, c in enumerate(ot.CENSUS)
    )
    monkeypatch.setattr(ot, "CENSUS", bad)
    with pytest.raises(ot.CensusDriftError) as e:
        ot.check_census_matches_tree()
    assert "both as a head read and as not one" in str(e.value)


def test_effective_lines_is_the_registered_lines_until_a_module_moves():
    eff = ot.effective_lines()
    for c in ot.CENSUS:
        mod = c["module"]
        if mod in ot.MOVED_LINES:
            assert eff[mod] == tuple(ot.MOVED_LINES[mod]["now"])
        else:
            assert eff[mod] == ot.REGISTERED_LINES[mod]


# ---- the moves' own refusals, kept apart from the registered ones --------------------------------


def test_the_moves_refusals_are_kept_apart_from_the_registered_ones():
    """The registration's own REFUSALS tuple is not rewritten because the moves needed more."""
    assert len(ot.REFUSALS) == 6
    assert len(ot.REFUSALS_ADDED_BY_THE_MOVES) == 2
    assert not set(ot.REFUSALS) & set(ot.REFUSALS_ADDED_BY_THE_MOVES)


def test_the_registered_lines_are_what_the_registration_committed():
    assert ot.declared() == ot.REGISTERED_LINES
    ot.check_registered_lines_are_not_edited()


def test_PLANTED_editing_a_registered_line_to_match_a_moved_tree_is_REFUSED(monkeypatch):
    """The temptation after a move is to update the census; that replaces the claim.

    Mutation: delete the `declared() != REGISTERED_LINES` raise and this test fails.
    """
    bad = tuple({**c, "lines": (c["lines"][0] + 19,)} if i == 0 else c for i, c in enumerate(ot.CENSUS))
    monkeypatch.setattr(ot, "CENSUS", bad)
    with pytest.raises(ot.CensusDriftError) as e:
        ot.check_registered_lines_are_not_edited()
    assert "no longer match what 24adf33 registered" in str(e.value)


def test_PLANTED_a_move_record_that_restates_the_claim_is_REFUSED(monkeypatch):
    """Mutation: delete the MOVED_LINES['was'] raise and this test fails."""
    mod = next(iter(ot.MOVED_LINES))
    monkeypatch.setattr(ot, "MOVED_LINES", {mod: {**ot.MOVED_LINES[mod], "was": (999,)}})
    with pytest.raises(ot.CensusDriftError) as e:
        ot.check_registered_lines_are_not_edited()
    assert "may not restate the claim" in str(e.value)


def test_every_moved_module_is_an_identity_module_and_says_what_moved():
    identity = set(ot.by_class("identity"))
    for mod, rec in ot.MOVED_LINES.items():
        assert mod in identity, mod
        assert rec["what"].strip(), mod
        assert rec["was"] == ot.REGISTERED_LINES[mod], mod


# ---- one refusal for the one invariant the whole wave stands on ----------------------------------


def test_a_clean_census_passes_the_invariant_check():
    ot.check_head_invariant({"elements": 10, "head_disagrees": 0, "coding_head_disagrees": 3}, "x")


def test_PLANTED_an_any_gene_head_disagreement_STOPS_THE_WAVE():
    """Mutation: delete the raise in check_head_invariant and this test fails."""
    with pytest.raises(ot.InvariantRefutedError) as e:
        ot.check_head_invariant({"elements": 10, "head_disagrees": 1}, "closure.window_elements")
    assert "THE WAVE STOPS" in str(e.value) and "closure.window_elements" in str(e.value)


def test_a_coding_head_disagreement_alone_is_counted_and_is_NOT_a_refutation():
    """The registered limit doing work: 10 of 4,794 had no predicted_coding where the cache named one."""
    ot.check_head_invariant({"elements": 4794, "coding_head_disagrees": 10}, "x")
