# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The 2026-10-03 licence-header sweep added comments and nothing else, proved blob against blob.

WHY A TEST AND NOT A NOTE. Fifteen engine files gained an `SPDX-License-Identifier: Apache-2.0`
header on 2026-10-03 (Albert's item 3b; `tests/test_engine_licence_headers.py` holds the record of
the defect). EIGHT OF THEM ARE IN THE 52-FILE IMPORT CLOSURE of the AstroREG-2 sender, and that
closure is the subject of a supervisor sign-off gating a paid study of 1,232 AlphaGenome requests:
the sign-off's condition is that the closure which runs is the closure that was dry-run reviewed,
checked as `genomeos.attribution.astrorun.closure_sha256(sender_closure())`. A digest over file
CONTENTS cannot tell a comment from a behaviour change -- which is exactly why it is a good check,
and also why a header sweep moves it. The digest therefore has to be re-issued, and what makes that
safe is not anyone's word that only comments moved. It is this.

WHAT IS ASSERTED, and both halves are needed:

1. **AST identity.** `ast.dump(ast.parse(old)) == ast.dump(ast.parse(new))`. The default `ast.dump`
   excludes line numbers, so inserting a comment line passes while any change to a statement, an
   expression or an import fails. The point worth stating: comments are not in the AST but DOCSTRINGS
   ARE, so a header written into a module docstring would FAIL here. Every one of these is a `#`
   comment above the docstring, which leaves the docstring as the module's first statement.
2. **Pure insertion.** No pre-existing line removed, none modified, every added line a `#` comment or
   blank. AST identity alone would permit DELETING a comment, and some comments are operative:
   `# noqa`, `# type: ignore` and `# pragma: no cover` all change tooling behaviour while being
   invisible to the AST. Pure insertion alone would permit inserting a comment that... is only a
   comment, but it says nothing about the rest of the file, which is what (1) is for.

BLOB TO BLOB, NOT BLOB TO WORKING TREE. Both sides are named by their git blob sha and read with
`git cat-file`, so this test is a claim about one historical change and stays true for ever. Reading
the "after" side from the working tree instead would have turned a one-off proof into a standing ban
on ever editing `genomeos/runtime/grn.py` again, which is not this lane's to impose: a later, real
change to one of these files is a legitimate change that moves the digest and needs the sign-off
re-issued on its own terms, and that is the coordinator's business and not a red here.

WHAT THIS DOES NOT CLAIM. It does not say the digest is right, that the sign-off was re-issued, or
that the other 44 files of the closure are untouched -- the sender's own machinery owns all three. It
says only that these 15 diffs carry no code.
"""

from __future__ import annotations

import ast
import difflib
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

#: The two lines inserted, in the house form `genomeos/coords.py` has carried since it was written.
HEADER = (
    "# SPDX-License-Identifier: Apache-2.0",
    "# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.",
)


@dataclass(frozen=True)
class Change:
    """One file, the blob before the sweep, the blob after it, and whether the sender imports it."""

    path: str
    before: str
    after: str
    in_sender_closure: bool


#: Measured, not transcribed: `git rev-parse HEAD:<path>` for the before side and
#: `git hash-object <path>` for the after side, at the commit that made the change. The closure
#: column is `genomeos.attribution.astrorun.sender_closure()` membership, read the same day.
#: `genomeos/coords.py` is the odd one: it ALREADY declared Apache-2.0, on line 6 under its module
#: docstring, where the five-line detector in `tests/test_engine_licence_headers.py` could not see it.
#: The window was deliberately not widened -- the plant there exists to keep a tag below the fifth
#: line unread -- so a second, identical tag was inserted at the top and line 6 left alone. Two
#: identical tags in one file is already the practice: `genomeos/bio.py` has one at line 1 and another
#: at line 23. A MOVE would have been a deletion plus an insertion and would have failed check (2).
CHANGES = (
    Change(
        "genomeos/coords.py",
        "506829eac37e115d02067d9dbff6f794886748a8",
        "b50069620ee9293dd7c6929c2bfe2d348dbcf1bb",
        True,
    ),
    Change(
        "genomeos/ir/model.py",
        "738ab2d3df886da1ca0fbe086728c183e2337af4",
        "34a0aca62e56b923c11cd0697a6422d0213fb1d0",
        True,
    ),
    Change(
        "genomeos/lang/parser.py",
        "dbe39faf6c2e031040c1acf5487c4aab0c45c8ac",
        "169b935b658500758a66735de0683a224a9aca98",
        True,
    ),
    Change(
        "genomeos/runtime/boolean.py",
        "5e93b6dc3f434e1f1728cf5a9906bd23bcfc0c9e",
        "97e2e6c8aad88520c9961dc7ee98723e73c1a9c6",
        False,
    ),
    Change(
        "genomeos/runtime/cell.py",
        "dd6e7ef9c90562168495c97118824105ef125d37",
        "7bc094a9b397a73b3604c40b8dd24d2540544fa6",
        True,
    ),
    Change(
        "genomeos/runtime/central_dogma.py",
        "75387adc2f5a7926b4941b0a70927bad62bf0e15",
        "f5efd3d20c19120080301e54eb96eb265cfc4e46",
        True,
    ),
    Change(
        "genomeos/runtime/compose.py",
        "86872c2fac896154c69561942aff138ae6a6177b",
        "e07ecb7792f89f38b9f21ae6906ae53d6536ca9d",
        False,
    ),
    Change(
        "genomeos/runtime/debugger.py",
        "e4944b1d130188f30b4de73e922ecbc79ae82ac7",
        "9b570830fae7be7da5bfd3d68d69ca28e33dab16",
        False,
    ),
    Change(
        "genomeos/runtime/gastrulation.py",
        "3116b0c4ff2c73d4222b9aa2223fe4ccf10f4faa",
        "994a444984b4750672511436e75d9dace3e39e62",
        False,
    ),
    Change(
        "genomeos/runtime/grn.py",
        "afa38086e822d7fb9302e5b393ed6e43ec330ce8",
        "53db9bbff056279e8edbc11d2b63564b5876020d",
        True,
    ),
    Change(
        "genomeos/runtime/sbml.py",
        "4977f158f506526b806fa8bfbf57334b5d7cc0da",
        "dce46be56a9a92367cd089364bd20dfb6fe964a5",
        False,
    ),
    Change(
        "genomeos/runtime/segmentation.py",
        "df0f2ebcf3dd631ef915ed293609c0787a8698a2",
        "1cdf76434490800970889f397cdd18535cec1da8",
        False,
    ),
    Change(
        "genomeos/runtime/spatial.py",
        "3490ad8cd23f141244ff01e440fe7b387fc1431d",
        "e8d551f3b058eb684dd3a88bc695f8e6bc50119a",
        False,
    ),
    Change(
        "genomeos/runtime/uncertainty.py",
        "926dff5ed68fd3cd038286e8cbb2ed5b33f8481e",
        "2c1a4f698ea02a0ea01d38ee39a7c6938046c26c",
        True,
    ),
    Change(
        "genomeos/runtime/variant_effect.py",
        "945c7ed3f468062460a1c83fc7f984942572ba3c",
        "76f15cb80651a806140eec9845826b4348df7fe3",
        True,
    ),
)


def _blob(sha: str) -> str:
    """The bytes git holds under that sha, decoded. Immutable by construction."""
    out = subprocess.run(["git", "cat-file", "blob", sha], capture_output=True, cwd=ROOT, check=False)
    assert out.returncode == 0, f"git has no blob {sha}: {out.stderr.decode(errors='replace')}"
    return out.stdout.decode()


def _requires_git() -> None:
    """Skips by name where the objects cannot be there at all, rather than passing quietly.

    The condition is `.git`, following `tests/test_number_provenance.py`: a tree without it holds no
    object store, so there is nothing to read either side from. This is not an engine test -- it
    imports no `genomeos` module, so `package_engine.engine_tests()` does not select it -- but a
    source export or a shallow copy is a tree it could still be run in.
    """
    if not (ROOT / ".git").exists():
        pytest.skip("no .git in this tree, so the blobs this test names cannot be read from it")


def ast_is_identical(before: str, after: str) -> bool:
    """Check (1). `ast.dump` without `include_attributes` carries no line numbers."""
    return ast.dump(ast.parse(before)) == ast.dump(ast.parse(after))


def insertions_and_removals(before: str, after: str) -> tuple[list[str], list[str]]:
    """Every line the diff adds and every line it removes, in that order."""
    diff = list(difflib.ndiff(before.splitlines(), after.splitlines()))
    return (
        [line[2:] for line in diff if line.startswith("+ ")],
        [line[2:] for line in diff if line.startswith("- ")],
    )


def is_pure_comment_insertion(before: str, after: str) -> bool:
    """Check (2). Nothing removed or modified, and every added line a comment or blank."""
    added, removed = insertions_and_removals(before, after)
    if removed or not added:
        return False
    return all(not line.strip() or line.lstrip().startswith("#") for line in added)


@pytest.mark.parametrize("change", CHANGES, ids=[c.path for c in CHANGES])
def test_the_diff_is_ast_identical_and_a_pure_comment_insertion(change: Change) -> None:
    """Both proofs, per file, over the two blobs the sweep names."""
    _requires_git()
    before, after = _blob(change.before), _blob(change.after)
    assert before != after, f"{change.path}: the two blobs are the same, so this proves nothing"
    assert ast_is_identical(before, after), (
        f"{change.path}: the syntax tree CHANGED between {change.before[:7]} and {change.after[:7]}. "
        "A licence header may not change code. If this file needed a real edit, it is no longer "
        "licensing hygiene: it is a change to reviewed code, and for the eight files in the "
        "AstroREG-2 sender's closure that is the coordinator's call."
    )
    added, removed = insertions_and_removals(before, after)
    assert not removed, (
        f"{change.path}: {len(removed)} line(s) were removed or modified: {removed[:3]}. Not a pure "
        "insertion. AST identity would not have caught this, because `# noqa`, `# type: ignore` and "
        "`# pragma: no cover` are operative and invisible to the AST."
    )
    assert is_pure_comment_insertion(before, after), f"{change.path}: added non-comment lines: {added}"
    assert list(HEADER) == added, f"{change.path}: added lines are not the two house header lines: {added}"


def test_the_population_is_what_it_says_it_is() -> None:
    """A table of fifteen rows proves nothing if it is secretly fourteen, or names one file twice."""
    assert len(CHANGES) == 15, len(CHANGES)
    assert len({c.path for c in CHANGES}) == 15, "a path appears twice"
    assert len({c.before for c in CHANGES}) == 15 and len({c.after for c in CHANGES}) == 15
    assert sum(c.in_sender_closure for c in CHANGES) == 8, (
        "eight of the fifteen are in the sender's closure; that count is why this file exists"
    )
    for c in CHANGES:
        assert (ROOT / c.path).is_file(), f"{c.path} is not in the tree"


def test_PLANTED_both_checks_refuse_what_they_are_for(tmp_path: Path) -> None:
    """Each check shown refusing, because a check that cannot be made to fail is a formality.

    Hermetic: every case is a string in this function. Nothing is written into `genomeos/`, and no
    blob is created, so a peer running this suite against the shared checkout sees no file of mine.
    """
    # An operative comment the AST cannot see. Spelled from parts so this file's own linter does
    # not read the test DATA as a directive aimed at this line.
    operative = "# " + "noqa" + ": E501"
    original = f'"""A module."""\n\nVALUE = 1  {operative}\n'
    header = "\n".join(HEADER) + "\n"

    good = header + original
    assert ast_is_identical(original, good) and is_pure_comment_insertion(original, good)

    # (1) refuses a changed statement, although the diff is otherwise an insertion of comments.
    code_changed = header + original.replace("VALUE = 1", "VALUE = 2")
    assert not ast_is_identical(original, code_changed), "a changed constant must fail check (1)"

    # (1) refuses a header written into the docstring, which is the trap worth planting: a docstring
    # IS in the AST, so this is a code change however much it looks like a comment.
    tagged_docstring = '"""A module.\n\nSPDX-License-Identifier: Apache-2.0\n"""'
    in_docstring = original.replace('"""A module."""', tagged_docstring)
    assert not ast_is_identical(original, in_docstring), "a tag inside the docstring must fail check (1)"
    assert is_pure_comment_insertion(original, in_docstring) is False, "nor is it a comment insertion"

    # (2) refuses a DELETED comment, which check (1) cannot see at all. A lint suppression is the real
    # case: dropping one turns the suppression off while the syntax tree is untouched.
    noqa_dropped = header + original.replace(f"  {operative}", "")
    assert ast_is_identical(original, noqa_dropped), "the AST really is blind to this, which is the point"
    assert not is_pure_comment_insertion(original, noqa_dropped), "a dropped comment must fail check (2)"

    # (2) refuses a MOVED tag, which is why `coords.py` gained a second one instead of moving its own.
    moved = f'"""A module."""\n{header}\nVALUE = 1  {operative}\n'
    assert not is_pure_comment_insertion(header + original, moved), "a move is a deletion plus an insertion"

    # and a file that changed nothing at all is not a proof either
    assert not is_pure_comment_insertion(original, original), "an empty diff must not pass as an insertion"
