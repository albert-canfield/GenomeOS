# SPDX-License-Identifier: AGPL-3.0-or-later
"""Every file the engine package ships must say Apache-2.0, and the ones that do not are PINNED.

Why this exists (2026-10-02). `tests/test_engine_boundary.py` enforces decision D40 by walking
IMPORTS: no Apache-2.0 engine file may import the AGPL-3.0 application. That is half of the boundary.
The other half was checked nowhere: `scripts/package_engine.py` copies the engine into a package
whose pyproject declares **Apache-2.0**, with no licence filter of any kind. So a file carrying an
`AGPL-3.0-or-later` header is copied, unchanged, into an Apache-2.0 wheel. An engine that imports
nothing from the application and ships AGPL files inside an Apache package is still a licensing
failure, and it is the kind discovered on release rather than in a suite.

THE POPULATION IS THE PACKAGER'S OWN, IMPORTED AND NOT COPIED. `package_engine.engine_files()` is
called directly, so this test cannot drift from what is actually shipped -- the same reason the
cell-provenance module imports a constant rather than restating it. Re-deriving the list here would
have missed two things that only reading the function shows: it also ships NON-.py files (the seven
`std/*.bio` modules), and it DOES exclude `__pycache__`, so the 62 `.pyc` files under the engine
packages are not shipped. Both of those were guessed wrong before the function was read.

MEASURED WHEN THIS WAS WRITTEN, over that population: **5 files declare `AGPL-3.0-or-later`**, 15
`.py` files declare nothing, and 7 `.bio` files declare nothing. The three sets are pinned APART
because they are different problems. An AGPL file inside an Apache package is a CONFLICT. A file
declaring nothing takes the package's licence by default, which is probably the intent for the
runtime modules but is nowhere stated. The `.bio` modules have no SPDX convention at all -- their
headers are prose comments -- so whether they should carry one is a question and not a defect.

RESOLVED 2026-10-03, by Albert, and the pins are edited in the same commit as the resolution because
that is what they were built to force. Three of the four numbers above changed and one did not:

  - **The 5 AGPL files are gone from the engine.** They moved to `genomeos/provenance/`, a package
    `PACKAGES` does not name, so the packager no longer copies them anywhere. `DECLARE_AGPL` is now
    EMPTY, which is the strongest form this pin has ever had: the equality now says the engine ships
    no AGPL file at all, and the first one to appear fails. The two AGPL *tests* went with them, by
    consequence rather than by edit -- `engine_tests()` selects by the imports it reads, and a test
    importing `genomeos.provenance` is no longer a test importing only the engine. That set is pinned
    here too now, at empty, because it was the same conflict going through a door this file did not
    watch.
  - **14 of the 15 header-less `.py` files now declare `Apache-2.0`**, in the two-line house form, at
    the top of the file where every other engine file carries it. `DECLARE_NOTHING_PY` is EMPTY.
  - **The 15th, `genomeos/coords.py`, was never header-less.** It has declared `Apache-2.0` since it
    was written, on line 6, under a five-line module docstring, so `_declared()` read the first five
    lines and did not see it: this pin held a mislabel, and the file was not the thing that was wrong.
    WIDENING THE WINDOW WAS CONSIDERED AND REJECTED. A window that skipped the module docstring would
    have read line 6, and it would have cost the property
    `test_PLANTED_the_detector_reads_each_of_the_three_states` exists for -- that a tag below the fifth
    line is NOT read -- for one file's sake. So the detector is untouched and the plant keeps its full
    force. `coords.py` instead has a second, identical `Apache-2.0` comment INSERTED at the top, its
    line 6 left exactly where it was. Two identical tags in one file is already the practice here:
    `genomeos/bio.py` carries one at line 1 and another at line 23. The insertion also keeps that
    file's diff a PURE INSERTION of comment lines, which all 15 share and which matters for the 8 of
    them that are in the AstroREG-2 sender's reviewed 52-file import closure -- a digest over file
    contents cannot tell a comment from a behaviour change, so the proof that nothing but comments
    moved is in `tests/test_licence_headers_are_insertions.py`.
  - **The 7 `.bio` files are untouched.** There is still no SPDX convention for BioLang modules, so
    the open question stays open and nothing was invented for it.

WHAT THIS TEST DOES AND DOES NOT DECIDE. It does not move a file, change a header, or alter the
packager: each of those is a licensing decision, and licensing decisions are Albert's. It PINS the
known sets EXACTLY, so that:

  - an AGPL file entering the engine fails immediately, which is the regression that matters. It
    was "a SIXTH" while the pin held five; since 2026-10-03 the pin is empty and it is the FIRST;
  - RESOLVING one forces this pin to be edited in the same commit, so it cannot rot into a permanent
    allowlist that quietly grows. An allowlist nobody has to update is how a violation becomes the
    status quo, and these are equalities rather than subset checks for exactly that reason.

The pins are a record of a defect, not permission for it.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

APACHE = "SPDX-License-Identifier: Apache-2.0"

#: Shipped files that DECLARE the application's licence. EMPTY since 2026-10-03: the five that were
#: here moved to `genomeos/provenance/` and are shipped by nothing. Pinned exactly, as an equality on
#: the empty set, so the FIRST one to appear fails rather than the sixth.
DECLARE_AGPL: tuple[str, ...] = ()

#: Shipped engine TESTS that declare the application's licence. The same conflict through a door this
#: file did not watch until the move: `package_engine.engine_tests()` selects test files by their
#: imports and copies them into the Apache-2.0 package as its own suite (LICENSING.md names them a
#: separate Apache entry), and nothing checked their headers. Two were AGPL --
#: `tests/test_rule_number_sources.py` and `tests/test_rule_number_sources_findings.py` -- and both
#: left the selection when the modules they import moved. EMPTY, by equality, for the same reason.
DECLARE_AGPL_TESTS: tuple[str, ...] = ()

#: Shipped `.py` files declaring NO licence. EMPTY since 2026-10-03: 14 of the 15 were given the
#: two-line Apache-2.0 house header, and the 15th (`genomeos/coords.py`) had one all along on line 6 --
#: a mislabel this pin inherited from the five-line window, corrected by inserting a second identical
#: tag at the top of that file rather than by widening the window. Pinned separately from the AGPL
#: set, as it always was, so that a conflict can never be mistaken for the milder state of taking the
#: package's licence by default.
DECLARE_NOTHING_PY: tuple[str, ...] = ()

#: Shipped engine TESTS declaring NO licence. The milder state on the tests side, pinned apart from
#: the conflict for the same reason as on the files side. One today: LICENSING.md states that the
#: selected tests are Apache-2.0, "decided by the copyright holder on 2026-09-28", and 30 of the 31
#: say so in their own first lines while this one says nothing and takes it by default. Recorded, not
#: resolved -- putting a header on a file is a licensing act and licensing acts are Albert's.
DECLARE_NOTHING_TESTS = ("tests/test_unstated_confidence.py",)

#: Shipped non-`.py` files: the standard library's BioLang modules. They carry prose comment headers
#: and there is no SPDX convention for `.bio` at all, so this is an open question, not a defect.
DECLARE_NOTHING_BIO = (
    "genomeos/std/ageing.bio",
    "genomeos/std/cell_types.bio",
    "genomeos/std/development.bio",
    "genomeos/std/human_stages.bio",
    "genomeos/std/human_turnover.bio",
    "genomeos/std/methylation.bio",
    "genomeos/std/signalling.bio",
)


def _packager():
    """`scripts/package_engine.py`, loaded by path because `scripts/` is not an importable package."""
    path = ROOT / "scripts" / "package_engine.py"
    spec = importlib.util.spec_from_file_location("package_engine_for_licence_test", path)
    assert spec and spec.loader, path
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def shipped() -> list[Path]:
    """Exactly what the packager copies, from the packager's own function."""
    return list(_packager().engine_files())


def _declared(path: Path) -> str | None:
    head = "".join(path.read_text(errors="replace").splitlines(keepends=True)[:5])
    for line in head.splitlines():
        if "SPDX-License-Identifier:" in line:
            return line.split("SPDX-License-Identifier:", 1)[1].strip()
    return None


def _rel(p: Path) -> str:
    return str(p.relative_to(ROOT))


def test_the_engine_ships_no_file_declaring_the_application_s_licence_beyond_the_pinned_five():
    """The conflict, pinned by equality so resolving one forces this list to be edited."""
    found = sorted(_rel(p) for p in shipped() if (_declared(p) or "").startswith("AGPL"))
    assert found == sorted(DECLARE_AGPL), (
        "the set of shipped engine files declaring AGPL has CHANGED.\n"
        f"  found:  {found}\n"
        f"  pinned: {sorted(DECLARE_AGPL)}\n"
        "An addition is a new licensing conflict: scripts/package_engine.py copies these into a "
        "package declaring Apache-2.0 with no licence filter, so the file ships under a licence it "
        "does not claim. A removal means one was resolved -- update this pin in the same commit, "
        "because a pin nobody has to edit becomes a permanent allowlist."
    )


def test_PLANTED_an_agpl_file_entering_the_engine_fails_this_check(tmp_path):
    """The counterfactual. A pin that cannot be shown to refuse is a list, not a check.

    It was a SIXTH file when the pin held five. Since the five moved out on 2026-10-03 the pin is the
    empty set and the planted file is the FIRST, which is the same plant doing more work: there is no
    longer any member whose presence could make a planted one look ordinary.

    The plant is hermetic: a file is written under `tmp_path` and `shipped()` is monkeypatched to
    return the real population plus that file. Nothing is written into `genomeos/`, because a peer
    running a suite against this shared checkout would see such a file and go red for a reason that
    is about no tree anyone intends to ship -- which is the defect class this project spent today
    removing.
    """
    import pytest

    planted = tmp_path / "genomeos" / "lang" / "planted_sixth.py"
    planted.parent.mkdir(parents=True)
    planted.write_text("# SPDX-License-Identifier: AGPL-3.0-or-later\nVALUE = 1\n")
    assert (_declared(planted) or "").startswith("AGPL"), "the detector must see the planted header"

    real = shipped()

    # The comparison the pinned check makes, run over the planted population. `_rel` is not used
    # here because the planted file is outside ROOT; the point is the SET equality, not the spelling.
    def found_agpl(files: list[Path]) -> list[str]:
        return sorted(
            p.name if not str(p).startswith(str(ROOT)) else _rel(p)
            for p in files
            if (_declared(p) or "").startswith("AGPL")
        )

    before = found_agpl(real)
    assert before == sorted(DECLARE_AGPL), "the real population must match the pin before planting"
    assert before == [], "the pin is empty today, so the real population must be empty too"

    after = found_agpl([*real, planted])
    assert after != sorted(DECLARE_AGPL), "planting an AGPL file must change the set"
    assert "planted_sixth.py" in after, "the planted file must be the thing that changed it"
    with pytest.raises(AssertionError):
        assert after == sorted(DECLARE_AGPL)


def test_PLANTED_the_detector_reads_each_of_the_three_states(tmp_path):
    """A detector that answered the same for every header would make all the pins vacuous."""
    cases = {
        "agpl.py": ("# SPDX-License-Identifier: AGPL-3.0-or-later\n", "AGPL-3.0-or-later"),
        "apache.py": ("# SPDX-License-Identifier: Apache-2.0\n", "Apache-2.0"),
        "bare.py": ('"""No licence header at all."""\n', None),
        "late.py": ("\n" * 9 + "# SPDX-License-Identifier: Apache-2.0\n", None),
    }
    for name, (text, expected) in cases.items():
        p = tmp_path / name
        p.write_text(text)
        assert _declared(p) == expected, name
    # `late.py` is the one worth stating: a tag below the fifth line is NOT read, so a file cannot
    # satisfy this check by burying its licence where a reader would not look for it.
    # UNCHANGED ON 2026-10-03, deliberately. `genomeos/coords.py` declares Apache-2.0 on line 6, under
    # a five-line docstring, and this window was NOT widened to reach it: a tag below the fifth line
    # staying unread is the property this plant exists for, and it is worth more than the one file.
    # `coords.py` got a second, identical tag INSERTED at the top instead, which is what
    # `genomeos/bio.py` already does -- one at line 1 and another at line 23, same licence, both
    # comments. Insertion also keeps that file's diff a pure insertion, which matters because it is
    # one of the 52 files in the AstroREG-2 sender's reviewed import closure.


def shipped_tests() -> list[Path]:
    """Exactly the test files the packager copies, from the packager's own function.

    Imported and not re-derived, for the reason the module docstring gives about `engine_files()`:
    `engine_tests()` selects by READING each test's imports, so a list restated here would drift the
    moment a test gained or lost one.
    """
    return list(_packager().engine_tests())


def test_the_engine_ships_no_TEST_declaring_the_application_s_licence():
    """The door this file did not watch until 2026-10-03, pinned the same way and apart.

    The package ships a test suite of its own: `engine_tests()` picks the tests that import only the
    engine and nothing from `scripts/`, and they are copied into the Apache-2.0 package exactly as the
    modules are. LICENSING.md makes them a separate Apache-2.0 entry, "taking precedence over the
    `tests/` entry above", decided by the copyright holder on 2026-09-28 -- so an AGPL header on one of
    them is the same conflict as an AGPL header on a module, and the pins above could not see it,
    because their population is `engine_files()`.

    Two were AGPL when this was written: `tests/test_rule_number_sources.py` and
    `tests/test_rule_number_sources_findings.py`. Neither was edited. They left the selection because
    the modules they import moved to `genomeos.provenance`, and `engine_tests()` reads imports -- which
    is why the fix for the modules fixed the tests without anyone touching them.
    """
    found = sorted(_rel(p) for p in shipped_tests() if (_declared(p) or "").startswith("AGPL"))
    assert found == sorted(DECLARE_AGPL_TESTS), (
        "the set of shipped engine TESTS declaring AGPL has CHANGED.\n"
        f"  found:  {found}\n"
        f"  pinned: {sorted(DECLARE_AGPL_TESTS)}\n"
        "An addition ships an AGPL test inside a package declaring Apache-2.0. A removal means one "
        "was resolved -- update this pin in the same commit, for the same reason as the one above."
    )
    assert shipped_tests(), "no tests were selected at all, so this check would pass on nothing"


def test_the_shipped_tests_are_either_apache_or_pinned_as_declaring_nothing():
    """The positive half on the tests side, so the two pins above are not vacuous."""
    pinned = set(DECLARE_AGPL_TESTS) | set(DECLARE_NOTHING_TESTS)
    assert len(pinned) == len(DECLARE_AGPL_TESTS) + len(DECLARE_NOTHING_TESTS), (
        "a test appears in both pins, so one of them is not saying what it claims"
    )
    silent = sorted(_rel(p) for p in shipped_tests() if _declared(p) is None)
    assert silent == sorted(DECLARE_NOTHING_TESTS), (
        "the set of shipped engine tests carrying NO SPDX header has changed.\n"
        f"  found:  {silent}\n"
        f"  pinned: {sorted(DECLARE_NOTHING_TESTS)}\n"
        "LICENSING.md says the selected tests are Apache-2.0; add the header rather than extending "
        "this pin, and edit the pin in the same commit."
    )
    offenders = [
        _rel(p)
        for p in shipped_tests()
        if _rel(p) not in pinned
        and APACHE not in "".join(p.read_text(errors="replace").splitlines(keepends=True)[:5])
    ]
    assert not offenders, f"shipped engine tests outside both pins must carry '{APACHE}': {offenders}"


def test_PLANTED_an_agpl_test_entering_the_engine_suite_fails_this_check(tmp_path):
    """Hermetic, for the reason the other plant gives: nothing is written into `tests/`."""
    import pytest

    planted = tmp_path / "test_planted_agpl.py"
    planted.write_text("# SPDX-License-Identifier: AGPL-3.0-or-later\ndef test_x():\n    assert True\n")
    assert (_declared(planted) or "").startswith("AGPL"), "the detector must see the planted header"

    real = shipped_tests()

    def found_agpl(files: list[Path]) -> list[str]:
        return sorted(
            p.name if not str(p).startswith(str(ROOT)) else _rel(p)
            for p in files
            if (_declared(p) or "").startswith("AGPL")
        )

    assert found_agpl(real) == sorted(DECLARE_AGPL_TESTS) == [], "the pin must hold before planting"
    after = found_agpl([*real, planted])
    assert "test_planted_agpl.py" in after, "the planted file must be the thing that changed the set"
    with pytest.raises(AssertionError):
        assert after == sorted(DECLARE_AGPL_TESTS)


def test_the_shipped_files_declaring_nothing_are_pinned_too():
    """Milder than the conflict and counted apart, so the two can never be confused."""
    found = sorted(_rel(p) for p in shipped() if _declared(p) is None)
    expected = sorted((*DECLARE_NOTHING_PY, *DECLARE_NOTHING_BIO))
    assert found == expected, (
        "the set of shipped engine files carrying NO SPDX header has changed.\n"
        f"  found:  {found}\n"
        f"  pinned: {expected}\n"
        "A `.py` file here takes the package's licence by default, which is probably the intent for "
        "the runtime modules but is stated nowhere per file: add a header rather than extending this "
        "pin. A `.bio` file here is an open question -- there is no SPDX convention for BioLang "
        "modules and their headers are prose."
    )


def test_every_other_shipped_file_says_apache_and_the_pins_are_disjoint():
    """The positive half: outside the pins, the rule holds with no exceptions."""
    pinned = set(DECLARE_AGPL) | set(DECLARE_NOTHING_PY) | set(DECLARE_NOTHING_BIO)
    assert len(pinned) == len(DECLARE_AGPL) + len(DECLARE_NOTHING_PY) + len(DECLARE_NOTHING_BIO), (
        "a file appears in more than one pin, so one of them is not saying what it claims"
    )
    offenders = [
        _rel(p)
        for p in shipped()
        if _rel(p) not in pinned
        and APACHE not in "".join(p.read_text(errors="replace").splitlines(keepends=True)[:5])
    ]
    assert not offenders, (
        f"shipped engine files outside every pin must carry '{APACHE}' in their first five lines: {offenders}"
    )


def test_the_pins_name_files_the_packager_still_ships():
    """A pin naming a file that is gone, or no longer shipped, would hide a real violation."""
    shipped_rel = {_rel(p) for p in shipped()}
    stale = [
        rel for rel in (*DECLARE_AGPL, *DECLARE_NOTHING_PY, *DECLARE_NOTHING_BIO) if rel not in shipped_rel
    ]
    shipped_test_rel = {_rel(p) for p in shipped_tests()}
    stale += [rel for rel in (*DECLARE_AGPL_TESTS, *DECLARE_NOTHING_TESTS) if rel not in shipped_test_rel]
    assert not stale, (
        f"pinned files the packager no longer ships, so the pin is stale and hides nothing: {stale}"
    )


def test_the_packager_applies_no_licence_filter_which_is_why_this_file_exists():
    """The premise, read from the source rather than assumed.

    If someone gives the packager a licence filter, this fails and this whole file should be
    reconsidered -- the right outcome, because the pins would then guard something that no longer
    happens.

    The question is asked of `engine_files()` ITSELF and not of the module, because the module has
    an SPDX header of its own and writes three more into the files it generates. A first draft here
    asserted `"SPDX" not in src` over the whole file and failed immediately on line 1 -- the
    module's own `# SPDX-License-Identifier: AGPL-3.0-or-later`. Scope the question to the thing
    whose behaviour is the premise.
    """
    import inspect

    src = inspect.getsource(_packager().engine_files)
    for word in ("SPDX", "licence", "license", "Apache", "AGPL"):
        assert word not in src, (
            f"package_engine.engine_files() now mentions {word!r}, so it may select files by "
            "licence. This file's premise is that it copies every engine file unfiltered, and the "
            "pins below exist only because of that. Re-read the function."
        )
    assert "__pycache__" in src, (
        "engine_files() no longer excludes __pycache__, so the 62 .pyc files under the engine "
        "packages would be shipped. That is a packaging defect in its own right and this file's "
        "population would silently include them."
    )
