# SPDX-License-Identifier: AGPL-3.0-or-later
"""The compiled header's `when: cell_type` sentence, and that the false one cannot return.

The superseded sentence said every rule is gated on a cell the element had been assayed or scored
within. For a predicted rule that is false - the compiled cell is the argmax over the scorer's
RNA-seq track axis - and on chr21 those are 5,174 of 5,176 rules. These tests assert:

* the exact superseded sentence is absent from the compiler and from the new program, and present
  only in the quarantine module that exists to migrate v1;
* the agreed replacement is present, as one sentence, with both counts;
* the new program's rule lines are byte-identical to v1's, so no figure moved;
* v1's bytes are untouched and still match the sha256 `data/results/label_gene.json` pins;
* the generator REFUSES a header whose counts disagree with the rules underneath it.

Each one is planted against before it is trusted: `test_counterfactual_*` restore the defect in a
copy and show the assertion fails, so a test that passes either way is caught here.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import pytest

from genomeos.attribution import compile as ac
from genomeos.attribution import reheader as rh

REPO = Path(__file__).resolve().parent.parent
COMPILER = REPO / "genomeos" / "attribution" / "compile.py"
V1 = REPO / "data" / "organisms" / "human" / "noncoding_chr21.bio"
V2 = REPO / "data" / "organisms" / "human" / "noncoding_chr21_v2.bio"

#: the sha256 `data/results/label_gene.json` and `label_gene_registration.json` name for v1
V1_PINNED_SHA256 = "778615512beaa46ddbb02f3cbee8502ceb1ff8d07dd43ce553e086b29e035192"

#: the agreed wording, as one sentence with chr21's counts
AGREED = (
    "The cell named is the tissue whose predicted expression of the target gene changed most when "
    "the element is deleted (over 371 RNA-seq tracks; 5,174 rules), or the cell a CRISPRi screen "
    "measured (2 rules)."
)


def _unwrap(text: str) -> str:
    """The header's comment lines as running prose, so a wrapped sentence can be matched whole."""
    out = []
    for ln in text.splitlines():
        if ln.startswith("rule ") or ln.startswith("element "):
            break
        if ln.startswith("# "):
            out.append(ln[2:])
        elif ln != "#":
            out.append("\n")
    return re.sub(r"[ \t]+", " ", " ".join(out))


# ---- the superseded sentence is gone from the code ----


def test_the_superseded_sentence_is_absent_from_the_compiler() -> None:
    assert rh.CELL_SENTENCE_SUPERSEDED not in COMPILER.read_text()


def test_the_superseded_sentence_lives_only_in_the_quarantine_module() -> None:
    carriers = sorted(
        p.relative_to(REPO).as_posix()
        for p in (REPO / "genomeos").rglob("*.py")
        if rh.CELL_SENTENCE_SUPERSEDED in p.read_text()
    )
    assert carriers == ["genomeos/attribution/reheader.py"]


def test_no_compiler_header_can_emit_the_superseded_sentence() -> None:
    for pred, exp in ((5174, 2), (0, 0), (1, 0), (999999, 7)):
        paragraph = " ".join(ln[2:] for ln in ac.cell_sentence_lines(pred, exp))
        assert rh.CELL_SENTENCE_SUPERSEDED not in paragraph


# ---- the new sentence is there, and says what was agreed ----


def test_the_compiler_emits_the_agreed_sentence_with_computed_counts() -> None:
    paragraph = " ".join(ln[2:] for ln in ac.cell_sentence_lines(5174, 2))
    assert AGREED in paragraph


def test_the_counts_are_computed_and_not_typed() -> None:
    paragraph = " ".join(ln[2:] for ln in ac.cell_sentence_lines(12, 3))
    assert "12 rules), or the cell a CRISPRi screen measured (3 rules)." in paragraph
    assert "5,174" not in paragraph and "(2 rules)" not in paragraph


def test_the_compiler_comment_lines_stay_within_the_program_width() -> None:
    for pred, exp in ((5174, 2), (123456789, 987654321)):
        assert all(len(ln) <= 100 for ln in ac.cell_sentence_lines(pred, exp))


# ---- the new program ----


def test_v1_bytes_are_untouched_and_still_match_the_pinned_sha256() -> None:
    import hashlib

    assert hashlib.sha256(V1.read_bytes()).hexdigest() == V1_PINNED_SHA256
    pin = json.loads((REPO / "data" / "results" / "label_gene.json").read_text())
    named = [
        i
        for i in pin["result_manifest"]["inputs"]
        if i.get("path", "").endswith("data/organisms/human/noncoding_chr21.bio")
    ]
    assert [i["sha256"] for i in named] == [V1_PINNED_SHA256]


def test_the_new_program_carries_the_agreed_sentence_and_not_the_old_one() -> None:
    text = V2.read_text()
    assert rh.CELL_SENTENCE_SUPERSEDED not in text
    assert "measured or predicted in" not in text
    assert AGREED in _unwrap(text)


def test_no_figure_moved_the_rule_lines_are_byte_identical() -> None:
    v1 = [ln for ln in V1.read_text().splitlines() if ln.startswith("rule ")]
    v2 = [ln for ln in V2.read_text().splitlines() if ln.startswith("rule ")]
    assert len(v1) == 5176
    assert [a for a, b in zip(v1, v2, strict=True) if a != b] == []
    assert v1 == v2


def test_the_new_program_differs_from_v1_in_the_cell_paragraph_and_nowhere_else() -> None:
    import difflib

    v1, v2 = V1.read_text().splitlines(), V2.read_text().splitlines()
    removed = [ln for ln in difflib.unified_diff(v1, v2, n=0, lineterm="") if ln.startswith("-#")]
    added = [ln for ln in difflib.unified_diff(v1, v2, n=0, lineterm="") if ln.startswith("+#")]
    assert [ln[1:] for ln in removed] == list(rh.CELL_PARAGRAPH_SUPERSEDED)
    assert [ln[1:] for ln in added] == ac.cell_sentence_lines(5174, 2)
    changed = [
        ln
        for ln in difflib.unified_diff(v1, v2, n=0, lineterm="")
        if ln[:1] in "+-" and ln[:3] not in ("---", "+++") and not ln.startswith(("-#", "+#"))
    ]
    assert changed == []


def test_the_new_program_is_the_generator_applied_to_v1() -> None:
    assert rh.reheader_cell_sentence(V1.read_text()) == V2.read_text()


def test_the_counts_written_into_the_new_program_are_its_own() -> None:
    assert ac.count_rules(V2.read_text()) == (5174, 2)


# ---- the generator refuses rather than ship a header that disagrees with its body ----


def test_the_generator_refuses_a_program_with_an_unaccounted_rule_tier() -> None:
    text = rh.CELL_PARAGRAPH_SUPERSEDED[0] + "\n"
    text += "\n".join(rh.CELL_PARAGRAPH_SUPERSEDED[1:]) + "\n"
    text += 'rule E1 activates G { evidence: curated "a tier the cell sentence does not speak for" }\n'
    with pytest.raises(ValueError, match="only some of them"):
        rh.reheader_cell_sentence(text)


PREDICTED_RULE = (
    'rule E1 activates G { when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" }'
)


def _synthetic(n: int = 1) -> str:
    """A minimum program: the superseded paragraph and `n` predicted rule lines, nothing else."""
    return "\n".join([*rh.CELL_PARAGRAPH_SUPERSEDED, *([PREDICTED_RULE] * n)]) + "\n"


def test_the_generator_refuses_when_the_replacement_adds_a_rule(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        rh, "cell_sentence_lines", lambda p, e: ["# a paragraph that smuggles a rule", PREDICTED_RULE]
    )
    with pytest.raises(ValueError, match="rule counts differ from the ones written into it"):
        rh.reheader_cell_sentence(_synthetic())


def test_the_generator_refuses_when_a_rule_line_moves(monkeypatch: pytest.MonkeyPatch) -> None:
    """The last guard, reached when the counts agree but the rule lines do not."""
    monkeypatch.setattr(rh, "count_rules", lambda text: (1, 0))
    monkeypatch.setattr(
        rh, "cell_sentence_lines", lambda p, e: ["# a paragraph that smuggles a rule", PREDICTED_RULE]
    )
    with pytest.raises(ValueError, match=r"\+1 were added or dropped; a header fix changes no rule"):
        rh.reheader_cell_sentence(_synthetic())


def test_the_generator_refuses_a_program_without_the_superseded_paragraph() -> None:
    with pytest.raises(ValueError, match="not in this program, verbatim"):
        rh.reheader_cell_sentence(V2.read_text())


def test_the_compiler_refuses_a_header_whose_counts_disagree_with_its_rules() -> None:
    src = COMPILER.read_text()
    assert "refusing to write a header that contradicts its own rules" in src
    assert "n_predicted_rules = len(elements)" in src


# ---- the planted counterfactuals: each assertion above is shown to fail on the defect ----


def test_counterfactual_restoring_the_old_sentence_fails_the_program_test(tmp_path: Path) -> None:
    """Put v1's paragraph back into a copy of v2; the program assertions must then fail."""
    restored = V2.read_text().replace(
        "\n".join(ac.cell_sentence_lines(5174, 2)), "\n".join(rh.CELL_PARAGRAPH_SUPERSEDED)
    )
    assert restored != V2.read_text()
    copy = tmp_path / "noncoding_chr21_counterfactual.bio"
    copy.write_text(restored)
    text = copy.read_text()
    # the two assertions of test_the_new_program_carries_the_agreed_sentence_and_not_the_old_one
    assert rh.CELL_SENTENCE_SUPERSEDED in text, "the old sentence is back, as planted"
    assert AGREED not in _unwrap(text), "the agreed sentence is gone, as planted"
    # and the rule lines are untouched, so the counterfactual isolates the header alone
    v1 = [ln for ln in V1.read_text().splitlines() if ln.startswith("rule ")]
    assert [ln for ln in text.splitlines() if ln.startswith("rule ")] == v1


def test_counterfactual_restoring_the_old_sentence_fails_the_compiler_test(tmp_path: Path) -> None:
    """Put the sentence back into a copy of the compiler; the compiler assertion must then fail."""
    restored = COMPILER.read_text().replace(
        "CELL_GATE_SENTENCE = (",
        f'_RESTORED_DEFECT = "{rh.CELL_SENTENCE_SUPERSEDED}"\nCELL_GATE_SENTENCE = (',
        1,
    )
    assert restored != COMPILER.read_text()
    copy = tmp_path / "compile_counterfactual.py"
    copy.write_text(restored)
    assert rh.CELL_SENTENCE_SUPERSEDED in copy.read_text(), "the compiler assertion would now fail"
    assert rh.CELL_SENTENCE_SUPERSEDED not in COMPILER.read_text(), "and passes on the real file"


# ---- the compiler's own refusal, and a drift planted in it ---------------------------------------

PRED = 'rule E1 activates G { when: cell_type = K562; evidence: predicted "AlphaGenome deletion, K562" }'
EXPT = 'rule E1_measured activates G { when: cell_type = K562; evidence: experimental "a CRISPRi screen" }'


def test_check_cell_counts_accepts_a_header_that_matches_its_rules() -> None:
    ac.check_cell_counts(["# header", PRED, PRED, EXPT], 2, 1, "chrT")


@pytest.mark.parametrize(
    ("declared", "where"),
    [((3, 1), "declares 3 predicted"), ((2, 0), "0 experimental"), ((0, 0), "declares 0 predicted")],
)
def test_check_cell_counts_refuses_a_declared_count_that_drifted(
    declared: tuple[int, int], where: str
) -> None:
    with pytest.raises(ValueError, match="contradicts its own rules"):
        ac.check_cell_counts(["# header", PRED, PRED, EXPT], *declared, "chrT")


def test_check_cell_counts_refuses_a_rule_tier_the_sentence_does_not_speak_for() -> None:
    other = 'rule E2 activates G { evidence: curated "GENCODE, a tier the sentence omits" }'
    with pytest.raises(ValueError, match="only some of the rules it speaks for"):
        ac.check_cell_counts(["# header", PRED, other], 1, 0, "chrT")


def test_compile_chromosome_calls_the_check_with_counts_it_derived() -> None:
    src = COMPILER.read_text()
    assert "n_predicted_rules = len(elements)" in src
    assert "check_cell_counts(lines, n_predicted_rules, n_measured_rules, chrom)" in src


def _module_from(source: str, name: str):
    """Load an edited copy of the compiler as a module, so a drift can be planted in it."""
    import importlib.util

    spec = importlib.util.spec_from_loader(name, loader=None)
    mod = importlib.util.module_from_spec(spec)
    mod.__file__ = str(COMPILER)
    exec(compile(source, name, "exec"), mod.__dict__)
    return mod


def test_counterfactual_a_drifted_predicted_count_makes_the_compiler_refuse_chr21() -> None:
    """Plant +1 on the compiler's own predicted count: compiling chr21 must then refuse."""
    src = COMPILER.read_text()
    old = "n_predicted_rules = len(elements)"
    assert src.count(old) == 1
    drifted = _module_from(src.replace(old, f"{old} + 1", 1), "compile_with_a_drifted_count")
    with pytest.raises(ValueError, match="declares 5175 predicted and 2 experimental"):
        drifted.compile_chromosome("chr21")
    # and the real compiler, with the count derived, does not refuse
    assert ac.compile_chromosome("chr21").count("\nrule ") == 5176
