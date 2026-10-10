# SPDX-License-Identifier: AGPL-3.0-or-later
"""Every Python file under `genomeos/` and `scripts/` declares one licence, and it is its side's.

WHY THIS EXISTS (2026-10-10, ROADMAP section 5 item 15). The codebase-map refresh (`fc21c29`) found
that 94 of the application modules under `genomeos/` declared no licence and that no test checked
PRESENCE. `tests/test_engine_licence_headers.py` already covers the engine: it reads what
`scripts/package_engine.py` ships and pins every shipped file that does not say Apache-2.0. Nothing
covered the application, so a new module could arrive without a header and nothing would notice.
LICENSING.md already decided the split: the engine set is Apache-2.0, and everything else under
`genomeos/`, together with `scripts/`, is AGPL-3.0-or-later. This file checks that each file says so.

THE SIDE COMES FROM THE PACKAGER, NOT FROM A COPY OF IT. The engine set is
`package_engine.engine_files()` filtered to `.py`. It is called here rather than restated, for the
reason the engine-licence test gives: a list copied here would drift from what is actually shipped.
Every `.py` file outside that set is application.

WHAT COUNTS AS A DECLARATION. A real `#` comment token reading `SPDX-License-Identifier: <id>`, found
with `tokenize`. A docstring or string that mentions the tag is not a declaration:
`genomeos/provenance/__init__.py` quotes it in prose on line 12, and `scripts/package_engine.py`
writes Apache tags into the files it generates. The first declaration has to come before the
module's first statement, not counting the module docstring. That allows the line-1 form and the
form below a module docstring (`genomeos/__init__.py`, `genomeos/cli.py`), and it refuses a tag
buried in code. This rule is deliberately looser than the five-line window in the engine test. That
window guards the shipped package; this file guards presence across the whole application.

THE FILES WITHOUT A HEADER ARE PINNED, BY EQUALITY AND WITH A REASON EACH. They are pinned in three
separate sets because they wait for three different things:

  - `PENDING_SIDE`: files whose side, engine or application, is still open in Albert's
    licence-placement question. EMPTY for `.py`. The 2026-10-03 resolution moved the five AGPL files
    out of the engine and gave the 14 header-less engine files Apache-2.0. What remains open is the
    SPDX convention for the 7 `genomeos/std/*.bio` files. They are not Python, so they are outside
    this population, and they stay pinned in the engine test's `DECLARE_NOTHING_BIO`.
  - `PENDING_CLOSURE`: application files in the AstroREG-2 sender's 52-file import closure.
    `astrorun.sender_closure()` hashes the WORKING TREE, and the sign-off for the paid run names
    that digest. So even a comment moves it, and run 2 cannot resume until someone signs again.
    Whether these headers go in before or after run 2 is on Albert's queue with the run and its
    signer, so this lane did not edit them.
  - `DEFERRED_HELD`: files a live lane holds on the work board. A later commit heads them when that
    lane finishes.

Each pin is an equality, not a subset. An addition (a new header-less file) fails here. So does a
removal, which means someone resolved a pin, and that commit has to edit the pin too. A pin that
nobody has to update becomes a permanent allowlist.
"""

from __future__ import annotations

import ast
import importlib.util
import io
import re
import tokenize
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

ENGINE_LICENCE = "Apache-2.0"
APPLICATION_LICENCE = "AGPL-3.0-or-later"

#: Files whose side is open in Albert's licence-placement question (night-shift ledger, "decide:
#: licence placement for the Apache engine"). Empty for `.py` since 2026-10-03; see the docstring.
PENDING_SIDE: dict[str, str] = {}

_CLOSURE = (
    "in the AstroREG-2 sender's import closure: sender_closure() hashes the working tree and the "
    "run's sign-off names that digest, so a header needs a fresh sign-off; Albert's call with run 2"
)

#: Header-less application files in the AstroREG-2 sender's import closure, measured 2026-10-10 as
#: the intersection of `astrorun.sender_closure()` (52 files) with the files declaring nothing.
PENDING_CLOSURE: dict[str, str] = {
    p: _CLOSURE
    for p in (
        "genomeos/attribution/astrorun.py",
        "genomeos/genome/__init__.py",
        "genomeos/genome/anatomy.py",
        "genomeos/genome/annotation.py",
        "genomeos/genome/bgzf.py",
        "genomeos/genome/domains.py",
        "genomeos/genome/fasta.py",
        "genomeos/genome/fetch.py",
        "genomeos/genome/genome.py",
        "genomeos/genome/index.py",
        "genomeos/genome/regulation.py",
        "genomeos/genome/regulatory.py",
        "genomeos/genome/repeats.py",
        "genomeos/genome/sequence.py",
        "genomeos/genome/signals.py",
        "genomeos/genome/unknown.py",
        "genomeos/genome/variants.py",
        "genomeos/knowledge/__init__.py",
        "genomeos/knowledge/cells.py",
        "genomeos/knowledge/go.py",
        "genomeos/knowledge/ligand_receptor.py",
        "genomeos/knowledge/obo.py",
        "genomeos/knowledge/reactome.py",
        "genomeos/lib/__init__.py",
        "genomeos/lib/catalog.py",
        "genomeos/lib/membership.py",
        "genomeos/predict/__init__.py",
        "genomeos/predict/alphagenome_adapter.py",
        "genomeos/results.py",
        "scripts/astroreg2_send.py",
    )
}

#: Header-less application files a live lane held on the work board when the sweep ran.
DEFERRED_HELD: dict[str, str] = {
    "genomeos/forge/__init__.py": "genomeos/forge/ held by lane-forgecal on 2026-10-10",
    "genomeos/forge/design.py": "genomeos/forge/ held by lane-forgecal on 2026-10-10",
}

#: Files that carry the same tag twice by design: one at the top and one older one below the module
#: docstring that the engine test's five-line window could not see. Both are comments and both name
#: the same licence. `tests/test_licence_headers_are_insertions.py` records why a second tag was
#: inserted instead of moving the first.
TWO_IDENTICAL_TAGS: dict[str, str] = {
    "genomeos/bio.py": "lines 1 and 26, both Apache-2.0",
    "genomeos/coords.py": "lines 1 and 8, both Apache-2.0; the second predates the sweep of 2026-10-03",
}

_TAG = re.compile(r"^#\s*SPDX-License-Identifier:\s*(\S+)\s*$")


def _packager():
    """`scripts/package_engine.py`, loaded by path because `scripts/` is not an importable package."""
    path = ROOT / "scripts" / "package_engine.py"
    spec = importlib.util.spec_from_file_location("package_engine_for_presence_test", path)
    assert spec and spec.loader, path
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _rel(p: Path) -> str:
    return p.relative_to(ROOT).as_posix()


def engine_py() -> set[str]:
    """The Apache side, as the packager ships it, `.py` only."""
    return {_rel(p) for p in _packager().engine_files() if p.suffix == ".py"}


def population() -> list[Path]:
    """Every `.py` under `genomeos/` and `scripts/`, without bytecode caches."""
    files = [*(ROOT / "genomeos").rglob("*.py"), *(ROOT / "scripts").rglob("*.py")]
    return sorted(p for p in files if "__pycache__" not in p.parts)


def declarations(text: str) -> list[tuple[int, str]]:
    """(line, licence id) for every comment token that is an SPDX tag; strings never count."""
    found = []
    for tok in tokenize.generate_tokens(io.StringIO(text).readline):
        if tok.type == tokenize.COMMENT:
            m = _TAG.match(tok.string)
            if m:
                found.append((tok.start[0], m.group(1)))
    return found


def first_statement_line(text: str) -> int | None:
    """The line of the module's first statement that is not its docstring, or None if it has none."""
    body = ast.parse(text).body
    first = body[0] if body else None
    if (
        isinstance(first, ast.Expr)
        and isinstance(first.value, ast.Constant)
        and isinstance(first.value.value, str)
    ):
        body = body[1:]
    if not body:
        return None
    node = body[0]
    decorators = getattr(node, "decorator_list", None) or []
    return min([node.lineno, *(d.lineno for d in decorators)])


def problems(rel: str, text: str, engine: set[str]) -> list[str]:
    """What is wrong with one file's declaration; empty when it is right. Pins are not applied here."""
    tags = declarations(text)
    if not tags:
        return ["declares no licence"]
    out = []
    ids = {lic for _, lic in tags}
    if len(ids) > 1:
        out.append(f"declares {len(ids)} different licences: {sorted(ids)}")
    if len(tags) > 1 and rel not in TWO_IDENTICAL_TAGS:
        out.append(f"carries {len(tags)} SPDX lines (lines {[n for n, _ in tags]}); one is the rule")
    want = ENGINE_LICENCE if rel in engine else APPLICATION_LICENCE
    side = "engine" if rel in engine else "application"
    if tags[0][1] != want:
        out.append(f"declares {tags[0][1]} but is on the {side} side, which is {want}")
    first = first_statement_line(text)
    if first is not None and tags[0][0] > first:
        out.append(f"its first SPDX line ({tags[0][0]}) comes after its first statement (line {first})")
    return out


@pytest.fixture(scope="module")
def census() -> dict[str, list[str]]:
    engine = engine_py()
    return {_rel(p): problems(_rel(p), p.read_text(), engine) for p in population()}


def test_the_header_less_files_are_exactly_the_pinned_ones(census):
    """Presence, as an equality on the files declaring nothing."""
    silent = sorted(rel for rel, probs in census.items() if probs == ["declares no licence"])
    pinned = sorted({**PENDING_SIDE, **PENDING_CLOSURE, **DEFERRED_HELD})
    assert silent == pinned, (
        "the set of Python files declaring NO licence has changed.\n"
        f"  new, not pinned:   {sorted(set(silent) - set(pinned))}\n"
        f"  pinned, now headed: {sorted(set(pinned) - set(silent))}\n"
        f"A new file needs '# SPDX-License-Identifier: {APPLICATION_LICENCE}' as its first line "
        f"(or '{ENGINE_LICENCE}' inside the engine set); do not extend a pin for it. A pinned file "
        "that is now headed was resolved: remove it from its pin in the same commit."
    )


def test_every_declared_licence_is_single_placed_and_matches_its_side(census):
    """Exactly one licence, in the header, and the right one for the file's side."""
    bad = {rel: probs for rel, probs in census.items() if probs and probs != ["declares no licence"]}
    assert not bad, "files whose licence declaration is wrong:\n" + "\n".join(
        f"  {rel}: {'; '.join(probs)}" for rel, probs in sorted(bad.items())
    )


def test_the_pins_are_disjoint_named_and_point_at_real_files():
    """A pin naming a missing file, or one file in two pins, would hide a real violation."""
    sets = (PENDING_SIDE, PENDING_CLOSURE, DEFERRED_HELD)
    assert sum(len(s) for s in sets) == len(set().union(*sets)), "a file appears in more than one pin"
    for pin in (*sets, TWO_IDENTICAL_TAGS):
        for rel, why in pin.items():
            assert (ROOT / rel).is_file(), f"{rel} is pinned but is not in the tree"
            assert why.strip(), f"{rel} is pinned without a reason"
    engine = engine_py()
    assert not (set(PENDING_CLOSURE) | set(DEFERRED_HELD)) & engine, "an engine file is pinned as application"
    assert set(TWO_IDENTICAL_TAGS) <= engine, "the two-tag files are both engine files"


def test_the_closure_pin_names_files_still_in_the_sender_closure():
    """The closure reason holds only while the file is in the closure; one that leaves gets its header."""
    from genomeos.attribution import astrorun

    closure = set(astrorun.sender_closure())
    left = sorted(set(PENDING_CLOSURE) - closure)
    assert not left, f"no longer in the AstroREG-2 sender closure, so head them and unpin: {left}"


def test_the_population_and_the_engine_set_are_what_they_say():
    """A check over an empty or mis-sided population passes on nothing."""
    engine = engine_py()
    rels = {_rel(p) for p in population()}
    assert len(engine) >= 30, f"the engine set looks wrong: {len(engine)} .py files"
    assert engine <= rels, f"engine files outside the population: {sorted(engine - rels)}"
    assert {"genomeos/bio.py", "genomeos/coords.py", "genomeos/version.py"} <= engine
    assert "genomeos/cli.py" in rels and "genomeos/cli.py" not in engine
    assert "scripts/package_engine.py" in rels and len(rels) - len(engine) > 400


def test_PLANTED_each_rule_refuses_what_it_is_for():
    """Each rule shown refusing on a string, so none is a formality. Nothing is written to the tree."""
    engine = {"genomeos/runtime/x.py"}
    agpl = f"# SPDX-License-Identifier: {APPLICATION_LICENCE}\n"
    apache = f"# SPDX-License-Identifier: {ENGINE_LICENCE}\n"
    doc = '"""A module."""\n'
    code = "VALUE = 1\n"

    assert problems("genomeos/x.py", agpl + doc + code, engine) == []
    assert problems("genomeos/x.py", doc + agpl + code, engine) == [], "below the docstring is allowed"
    assert problems("scripts/x.py", "#!/usr/bin/env python3\n" + agpl + doc, engine) == []
    assert problems("genomeos/runtime/x.py", apache + doc, engine) == []

    assert problems("genomeos/x.py", doc + code, engine) == ["declares no licence"]
    tagged_doc = '"""A module.\n\nSPDX-License-Identifier: AGPL-3.0-or-later\n"""\n'
    assert problems("genomeos/x.py", tagged_doc + code, engine) == ["declares no licence"], (
        "a tag inside a docstring is prose, not a declaration"
    )
    in_string = 'HEADER = "# SPDX-License-Identifier: Apache-2.0"\n'
    assert problems("genomeos/x.py", agpl + in_string, engine) == [], "a tag inside a string is data"

    assert any("application side" in s for s in problems("genomeos/x.py", apache + doc, engine))
    assert any("engine side" in s for s in problems("genomeos/runtime/x.py", agpl + doc, engine))
    assert any("different licences" in s for s in problems("genomeos/x.py", agpl + apache + doc, engine))
    assert any("one is the rule" in s for s in problems("genomeos/x.py", agpl + doc + agpl, engine))
    assert any("after its first statement" in s for s in problems("genomeos/x.py", doc + code + agpl, engine))

    # The pinned two-tag files pass because they are pinned, and only for the count.
    assert problems("genomeos/bio.py", apache + doc + apache, {"genomeos/bio.py"}) == []
    assert problems("genomeos/bio.py", apache + doc + agpl, {"genomeos/bio.py"}) != []
