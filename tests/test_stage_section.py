# SPDX-License-Identifier: AGPL-3.0-or-later
"""Staging one section of a shared file: nothing outside it may reach the index, and an unplaceable
section is refused rather than appended."""

import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("stage_section", Path("scripts/stage_section.py"))
ss = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ss)

BASE = [
    "# The file\n",
    "## Theirs\n",
    "their committed line\n",
    "## Mine\n",
    "my committed line\n",
    "## What comes next\n",
    "tail\n",
]
WORK = [
    "# The file\n",
    "## Theirs\n",
    "their committed line\n",
    "THEIR HALF-WRITTEN DRAFT\n",
    "## Mine\n",
    "my new line\n",
    "## What comes next\n",
    "tail\n",
]


def test_another_lanes_uncommitted_text_cannot_reach_the_index():
    out = ss.splice(BASE, WORK, "## Mine")
    assert "THEIR HALF-WRITTEN DRAFT\n" not in out
    assert "my new line\n" in out and "my committed line\n" not in out
    assert out == BASE[:3] + ["## Mine\n", "my new line\n"] + BASE[5:]


def test_a_section_already_committed_stages_as_no_change():
    assert ss.splice(BASE, BASE, "## Mine") == BASE


def test_the_section_ends_at_the_next_heading_not_the_end_of_file():
    out = ss.splice(BASE, WORK, "## Theirs")
    # taking their section carries their draft, which is why each lane names only its own heading
    assert "THEIR HALF-WRITTEN DRAFT\n" in out
    assert out[out.index("## Mine\n") + 1] == "my committed line\n"  # mine is untouched


def test_a_heading_absent_from_the_working_copy_is_refused():
    with pytest.raises(ss.RefusedError, match="nothing of yours"):
        ss.splice(BASE, WORK, "## Not a section of mine")


def test_a_section_new_to_the_committed_file_is_refused_without_an_anchor():
    """It must refuse, not append: the end of the file is where another lane is also writing."""
    work = BASE[:5] + ["## Brand new\n", "first draft\n"] + BASE[5:]
    with pytest.raises(ss.RefusedError, match="name the heading it goes before"):
        ss.splice(BASE, work, "## Brand new")


def test_a_new_section_with_an_anchor_goes_before_it():
    work = BASE[:5] + ["## Brand new\n", "first draft\n"] + BASE[5:]
    out = ss.splice(BASE, work, "## Brand new", "## What comes next")
    assert out.index("## Brand new\n") == out.index("## What comes next\n") - 2
    assert "THEIR HALF-WRITTEN DRAFT\n" not in out


def test_an_anchor_that_is_not_in_the_committed_file_is_refused():
    work = BASE[:5] + ["## Brand new\n", "first draft\n"] + BASE[5:]
    with pytest.raises(ss.RefusedError, match="anchor"):
        ss.splice(BASE, work, "## Brand new", "## No such heading")
