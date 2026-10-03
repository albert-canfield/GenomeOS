# SPDX-License-Identifier: AGPL-3.0-or-later
"""The guard on the predicted `silencer_like` key, and on prose about it.

Before this file every guard enforcing "a direction is never a function" sat on the MEASURED path
(`increase_links`, `hct116_count`, `direction_link`, `measured`) or on the COMPILED label (the
programme headers, the grammar axis). The PREDICTED key had none and prose had none, which is the
layer that generates the 230,839 calls and where all eight overstatements corrected in 7072f30
lived -- one of them in the predicted module's own docstring.

Each guard here is written to FIRE. The planted cases are not decoration: a guard whose subject
cannot be broken in a test has not been shown to reach its subject.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.attribution import silencerdoctrine as sd

ROOT = Path(__file__).resolve().parent.parent

#: The seven sentences 7072f30 corrected, verbatim as they stood at a5a43b3. This is the negative
#: set: the guard exists if and only if it catches all seven. They live in this file rather than in
#: the corpus, and `tests/` is not part of the corpus, so quoting them here offends nothing.
THE_SEVEN_CORRECTED_SENTENCES = (
    # docs/ALPHAGENOME.md -- no direction gloss at all
    "43 elements that the registry calls\nenhancer-like behave as silencers for the gene they move.",
    # docs/NODES-READER-WRITER.md -- the sentence runs over a line break
    'the registry\'s\n"enhancer-like" class, defined by chromatin marks, includes elements that\n'
    "behave as silencers for the gene they move.",
    # docs/PROGRESS.md:67 -- inside a parenthesis of a long bullet, no gloss
    "job `enhancer_targets_chr21` scored 200 distal enhancers (63.5% move a gene by 0.1 log2; "
    "43 behave as silencers);",
    # docs/PROGRESS.md:421
    "And 43 of the 127 elements that the ENCODE registry calls enhancer-like behave as silencers: "
    "expression rises when they are deleted.",
    # docs/PROGRESS.md:469 -- the inference made explicit with a "so"
    "1,157 elements that the ENCODE registry calls enhancer-like rise when deleted, so they act as "
    "silencers: the registry names a class of element, not a direction of effect.",
    # docs/ROADMAP.md:358
    "Two findings worth as much as the headline: the nearest-gene heuristic names the wrong gene "
    "almost a third of the time, and 1,157\nelements the registry calls enhancer-like behave as "
    "silencers, rising when\ndeleted.",
    # docs/ROADMAP.md:3558 -- a table cell whose own sentence carries no subject at all
    "90.2% of those sit inside the element's CTCF node and 71.2% are exactly the nearest TSS; "
    "1,157 behave as silencers. Needs a key;",
)


def _scan_string(tmp_path: Path, text: str, name: str = "planted.md") -> list[dict]:
    p = tmp_path / name
    p.write_text(text)
    return sd.scan_text(name, p)


# ---- the doctrine this guard enforces ------------------------------------------------------------


def test_the_grammar_still_forbids_the_inference_the_predicted_key_is():
    """The source sentence. Nothing asserted it before this file, in tests or in code."""
    grammar = (ROOT / "genomeos/lang/grammar.py").read_text()
    assert '"silencer": "a repressive element; never inferred from a direction of effect alone"' in grammar
    # and the axis the key actually belongs to is glossed as a direction and nothing else
    assert (
        '"represses_target": "removing it raises its target\'s expression (predicted or measured)"'
    ) in grammar


def test_the_predicted_module_still_emits_the_key_this_guard_watches():
    """If the emitter renames the key or stops deriving it from the direction field, say so here.

    `genomeos/predict/enhancer_target.py` is inside the frozen closure and is read, never written.
    """
    src = (ROOT / "genomeos/predict/enhancer_target.py").read_text()
    assert f'"{sd.KEY}": repress,' in src, "the key this guard watches is no longer emitted"
    assert 'repress += p["action"] == "represses"' in src, (
        "the count is no longer taken from the direction field alone, which is the whole reason it "
        "is a direction and not a function"
    )
    assert sd.KEY == "silencer_like"


# ---- part 1: the predicted key over the committed results ----------------------------------------


def test_the_population_carrying_the_key_is_the_three_predicted_families_plus_silenceragree():
    """77 tracked results carry the key; 75 of them are the predicted summary.

    Established by walking every tracked `data/results/*.json` for the key at any depth, not from
    any published figure. The 75 is three families of 24 chromosomes plus a genome-wide file each.
    The other two are this check's own outputs, which carry the key as a measured share.
    """
    report = sd.key_report(ROOT)
    assert report["by_family"] == {
        "constrained_targets": 25,
        "enhancer_targets": 25,
        "enhancer_targets_all": 25,
        "other": 2,
    }
    assert report["results_carrying_the_key"] == 77
    assert report["occurrences"] == 341
    others = [r for r in sd.results_carrying_the_key(ROOT) if sd.family(r) is None]
    assert sorted(others) == [
        "data/results/silenceragree.json",
        "data/results/silenceragree_registration.json",
    ]


def test_every_committed_result_carrying_the_key_keeps_it_a_count_of_a_direction():
    assert sd.key_report(ROOT)["violations"] == []


def test_the_count_is_recomputed_from_the_element_rows_wherever_they_are_committed():
    """48 of the 77 carry their own per-element rows, and there the count is exact, not bounded.

    Both equalities hold on all 48: the summary count equals the number of rows the model calls
    `represses`, AND it equals the number whose recorded log2 fold change is POSITIVE. The second is
    what makes the key a direction rather than a function in the committed bytes: it counts elements
    whose deletion RAISES the target.
    """
    report = sd.key_report(ROOT)
    assert report["recounted_from_their_own_element_rows"] == 48
    assert report["violations"] == []


@pytest.mark.parametrize(
    "value",
    [0.377, "silencer", True, None, -1],
    ids=["a_rate", "a_label", "a_bool", "nothing", "negative"],
)
def test_the_key_may_not_hold_anything_but_a_non_negative_int(value):
    bad = sd.key_violations("planted.json", {"summary": {sd.KEY: value, "n": 10}})
    assert bad, f"{value!r} passed as a count of elements"


def test_a_count_larger_than_the_population_beside_it_is_refused():
    bad = sd.key_violations("planted.json", {"summary": {sd.KEY: 11, "with_predicted_target": 10}})
    assert any("exceeds its population" in b for b in bad)


def test_the_named_target_population_must_be_strong_plus_weak():
    bad = sd.key_violations(
        "planted.json",
        {"summary": {sd.KEY: 3, "strong": 2, "weak": 5, "with_predicted_target": 9}},
    )
    assert any("not the named-target population" in b for b in bad)


def test_a_sibling_naming_a_role_is_refused():
    bad = sd.key_violations("planted.json", {"summary": {sd.KEY: 3, "n": 9, "silencer": True}})
    assert any("names a role" in b for b in bad)
    bad = sd.key_violations("planted.json", {"summary": {sd.KEY: 3, "n": 9, "molecular_role": "silencer"}})
    assert any("names a role" in b for b in bad)


def _planted_rows(actions_and_fcs):
    return {
        "summary": {sd.KEY: sum(1 for a, _ in actions_and_fcs if a == "represses")},
        "elements": [{"predicted": {"action": a, "log2_fold_change": fc}} for a, fc in actions_and_fcs],
    }


def test_the_recount_fires_when_the_summary_disagrees_with_its_own_rows():
    doc = _planted_rows([("represses", 0.4), ("activates", -0.4)])
    doc["summary"][sd.KEY] = 2
    bad = sd.recount_violations("planted.json", doc)
    assert any("rows say represses" in b for b in bad)


def test_the_recount_fires_when_a_represses_row_records_a_fall_not_a_rise():
    """The direction half. A row called `represses` whose target FELL breaks the key's meaning."""
    doc = _planted_rows([("represses", -0.4), ("activates", -0.9)])
    bad = sd.recount_violations("planted.json", doc)
    assert any("POSITIVE log2_fold_change" in b for b in bad)
    assert any("raises the target" in b for b in bad)


def test_the_recount_is_not_vacuous_on_a_well_formed_plant():
    assert (
        sd.recount_violations("planted.json", _planted_rows([("represses", 0.4), ("activates", -0.4)])) == []
    )


# ---- part 2: prose -------------------------------------------------------------------------------


def test_the_tracked_corpus_asserts_no_silencer_function_anywhere():
    """The positive set: the tree as 7072f30 left it, with the one frozen line excused."""
    assert sd.scan_tree(ROOT) == []


@pytest.mark.parametrize(
    "sentence", THE_SEVEN_CORRECTED_SENTENCES, ids=[f"corrected_{i}" for i in range(1, 8)]
)
def test_each_of_the_seven_corrected_sentences_is_caught(tmp_path, sentence):
    hits = _scan_string(tmp_path, sentence)
    assert hits, "a sentence 7072f30 had to correct would pass this guard"
    assert hits[0]["tier"] == "FUNCTION_ASSERTION"


def test_a_sentence_that_crosses_a_line_break_is_still_caught(tmp_path):
    """One of the seven ran over a line break, so a line-at-a-time scan would have missed it."""
    hits = _scan_string(tmp_path, "the predicted elements\nbehave as silencers for their gene.\n")
    assert len(hits) == 1 and hits[0]["phrase"] == "behave as silencers"


@pytest.mark.parametrize(
    "sentence",
    [
        "The element behaves as a silencer for that gene.",
        "These elements act as silencers.",
        "The element functions as a repressor of MYC.",
        "It serves as a silencer in K562.",
        "The element operates as a silencer.",
        "It works as a repressor.",
    ],
)
def test_the_function_assertion_tier_catches_every_as_a_silencer_form(tmp_path, sentence):
    assert _scan_string(tmp_path, sentence)


def test_the_function_assertion_tier_needs_no_subject_test(tmp_path):
    """Deliberate. All eight offences took this form and none of it is ever true of this layer, so
    it is refused with nothing about the predicted layer in the block at all."""
    assert _scan_string(tmp_path, "Something behaves as a silencer.")


def test_the_role_assertion_tier_fires_only_where_the_subject_is_the_predicted_key(tmp_path):
    predicted = "Of the predicted elements, 1,157 are silencers."
    assert _scan_string(tmp_path, predicted)[0]["tier"] == "ROLE_ASSERTION"
    # the same shape about an assayed element is a different subject and stays allowed
    assert _scan_string(tmp_path, "The ReSE fragment is a silencer by assay.", "a.md") == []


def test_the_role_assertion_tier_does_not_fire_on_the_doctrine_being_stated(tmp_path):
    """The project states the doctrine in sentences of exactly this shape, and they must survive."""
    for refusal in (
        "A predicted increase is never a silencer and never a repressor.",
        "A deletion direction is not a silencer label.",
        "The predicted element is silencer-like, not a silencer.",
    ):
        assert _scan_string(tmp_path, refusal, "doctrine.md") == []


def test_the_hedged_word_is_never_flagged(tmp_path):
    for ok in (
        "These predicted elements are silencer-like.",
        "The predicted count is silencer_like.",
        "1,157 predicted elements are silencer-like, meaning expression rises on deletion.",
    ):
        assert _scan_string(tmp_path, ok, "ok.md") == []


def test_python_code_is_not_scanned_only_its_comments_and_strings(tmp_path):
    src = "behaves_as_a_silencer = 1  # a name, not a sentence\nx = behaves_as_a_silencer\n"
    assert sd.scan_text("planted.py", _write(tmp_path, "planted.py", src)) == []
    commented = "# The element behaves as a silencer for that gene.\nx = 1\n"
    assert sd.scan_text("planted.py", _write(tmp_path, "c.py", commented))


def _write(tmp_path: Path, name: str, text: str) -> Path:
    p = tmp_path / name
    p.write_text(text)
    return p


def test_a_python_docstring_is_scanned(tmp_path):
    src = '"""A rise means\nthe element behaves as a silencer for that gene.\n"""\n'
    hits = sd.scan_text("planted.py", _write(tmp_path, "d.py", src))
    assert len(hits) == 1 and hits[0]["line"] == 2


# ---- the one allowlist entry ---------------------------------------------------------------------


def test_the_allowlist_holds_exactly_the_one_frozen_line():
    assert sd.ALLOWED == (("genomeos/predict/enhancer_target.py", 11),)


def test_the_frozen_line_is_excused_by_path_and_line_and_not_by_the_pattern():
    """With the allowlist off, the pattern DOES catch it, and it is the only offence in the tree."""
    unfiltered = sd.scan_tree(ROOT, allow=())
    assert [(h["path"], h["line"], h["tier"]) for h in unfiltered] == [
        ("genomeos/predict/enhancer_target.py", 11, "FUNCTION_ASSERTION")
    ]


def test_the_allowlist_entry_still_excuses_a_live_offence():
    """An allowlist entry covering nothing is a hole. This is the check that it still covers one."""
    assert sd.allowlist_still_earns_its_place(ROOT) == []


def test_the_allowlisted_file_is_on_the_paid_studys_own_computed_closure():
    """The allowlist's REASON, machine-checked rather than taken as prose.

    The reason the line is not corrected is that the file is inside a frozen paid-study closure. The
    closure is the computed import closure recorded in that study's registration, so the claim is
    checkable: if the file ever leaves it, the allowlist loses its justification and this fails.
    """
    reg = json.loads((ROOT / sd.CLOSURE_REGISTRATION).read_text())
    closure = reg["result_manifest"]["code_cleanliness"]["counting_path"]
    assert sd.ALLOWED[0][0] in closure
    assert "import closure" in reg["result_manifest"]["code_cleanliness"]["counting_path_is_computed"]


def test_the_allowlist_excuses_one_line_and_not_the_file(tmp_path):
    """A second offence in the frozen file, on any other line, is still refused."""
    src = (ROOT / sd.ALLOWED[0][0]).read_text()
    planted = src + "\n# The element behaves as a silencer for that gene.\n"
    hits = sd.scan_text(sd.ALLOWED[0][0], _write(tmp_path, "frozen.py", planted))
    assert [h["line"] for h in hits] == [len(planted.splitlines())]


# ---- the corpus itself ---------------------------------------------------------------------------


def test_the_corpus_reaches_every_kind_of_file_the_brief_names():
    """A guard over a corpus that silently shrank would pass forever."""
    tracked = sd._tracked(ROOT, sd.CORPUS)
    for required in (
        "docs/ROADMAP.md",
        "README.md",
        "genomeos/predict/enhancer_target.py",
        "genomeos/lang/grammar.py",
    ):
        assert required in tracked
    assert any(p.startswith("genomeos/web/static/") for p in tracked)
    assert any(p.startswith("data/demo/") and p.endswith(".bio") for p in tracked)
    assert len(tracked) > 250


def test_the_guards_own_module_is_in_the_corpus_and_passes_it():
    """No exclusion for the file that defines the patterns: it is scanned like everything else."""
    tracked = sd._tracked(ROOT, sd.CORPUS)
    assert "genomeos/attribution/silencerdoctrine.py" in tracked
    own = ROOT / "genomeos/attribution/silencerdoctrine.py"
    assert sd.scan_text("genomeos/attribution/silencerdoctrine.py", own) == []
