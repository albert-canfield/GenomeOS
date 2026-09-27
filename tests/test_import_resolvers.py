# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""`import protein:X` in both worlds: with GenomeOS installed, and with the engine alone.

The engine ships no proteome. GenomeOS declares its resolver as an entry point, so `bio` finds it
wherever GenomeOS is installed without anyone calling `register()`; with no resolver anywhere the
import stops with an error naming what is missing. The packaged engine's side of this is run in an
interpreter with no GenomeOS by `scripts/package_engine.py` (`tests/test_engine_package.py`).
"""

from __future__ import annotations

import importlib.metadata

import pytest

from genomeos.bio import main as bio
from genomeos.lang import parse
from genomeos.lang import parser as engine
from genomeos.lang.parser import BioLangError


@pytest.fixture
def no_registered():
    """IMPORT_RESOLVERS as a standalone `bio` starts: empty, nothing registered by the CLI. Emptied
    and refilled in place, since modules hold the dict itself by name."""
    saved = dict(engine.IMPORT_RESOLVERS)
    engine.IMPORT_RESOLVERS.clear()
    yield engine.IMPORT_RESOLVERS
    engine.IMPORT_RESOLVERS.clear()
    engine.IMPORT_RESOLVERS.update(saved)


def test_genomeos_declares_the_protein_resolver_as_an_entry_point():
    found = [e for e in importlib.metadata.entry_points(group=engine.RESOLVER_GROUP) if e.name == "protein"]
    assert [e.value for e in found] == ["genomeos.lib.biolang:resolve_protein"], (
        "reinstall the project (uv sync) so its entry points are current"
    )


def test_with_genomeos_installed_the_import_resolves_without_register(no_registered):
    from genomeos.lib.biolang import resolve_protein

    m = parse("module t\nimport protein:TP53\n")
    assert m.entities["TP53"].accession == "P04637"
    assert no_registered["protein"] is resolve_protein  # loaded once, then registered


def test_with_no_resolver_anywhere_the_import_names_what_is_missing(no_registered, monkeypatch):
    monkeypatch.setattr(importlib.metadata, "entry_points", lambda **_: [])
    with pytest.raises(BioLangError) as e:
        parse("module t\nimport protein:TP53\n")
    msg = str(e.value)
    assert "no resolver registered for import scheme 'protein'" in msg
    assert "the engine ships none" in msg and engine.RESOLVER_GROUP in msg
    assert no_registered == {}  # nothing was invented to fill the gap


def test_a_registered_resolver_comes_before_an_installed_one(no_registered):
    no_registered["protein"] = lambda symbol: f'protein {symbol} {{ evidence: curated "t"; confidence: 0.5 }}'
    m = parse("module t\nimport protein:ZZZ\n")
    assert m.entities["ZZZ"].confidence == 0.5


def test_bio_reports_an_unresolvable_import_as_an_error_not_a_crash(tmp_path, capsys):
    prog = tmp_path / "p.bio"
    prog.write_text("module p\nimport nothing:here\n")
    assert bio(["check", str(prog)]) == 2
    err = capsys.readouterr().err
    assert err.startswith("bio: no resolver registered for import scheme 'nothing'")
    assert "the engine ships none" not in err  # only schemes the engine knows get the longer reason
