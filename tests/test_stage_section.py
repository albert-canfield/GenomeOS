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


# A roadmap-shaped file: areas are `###` under one `##`, which is where the tool was wrong.
ROADMAP = [
    "# GenomeOS roadmap\n",
    "## The areas\n",
    "### A. BioLang and BioVM\n",
    "area A, owned by another lane\n",
    "### E. From one cell to an organism\n",
    "area E, committed text\n",
    "### J. Grammar by comparison\n",
    "area J, owned by another lane\n",
    "## Data jobs\n",
    "the jobs table\n",
]


def test_a_subsection_stops_at_the_next_subsection_not_the_next_top_level_heading():
    """The bug: ending at the next `## ` swept areas E through J into one splice, carrying two other
    lanes' text out under the staging lane's commit."""
    work = list(ROADMAP)
    work[5] = "area E, edited by me\n"
    work[3] = "AREA A, ANOTHER LANE'S HALF-WRITTEN DRAFT\n"
    work[7] = "AREA J, ANOTHER LANE'S HALF-WRITTEN DRAFT\n"
    out = ss.splice(ROADMAP, work, "### E. From one cell to an organism")
    assert "area E, edited by me\n" in out
    assert "AREA A, ANOTHER LANE'S HALF-WRITTEN DRAFT\n" not in out
    assert "AREA J, ANOTHER LANE'S HALF-WRITTEN DRAFT\n" not in out
    assert out == ROADMAP[:5] + ["area E, edited by me\n"] + ROADMAP[6:]


def test_a_subsection_that_is_last_in_its_parent_stops_at_the_next_top_level_heading():
    work = list(ROADMAP)
    work[7] = "area J, edited by me\n"
    out = ss.splice(ROADMAP, work, "### J. Grammar by comparison")
    assert out == ROADMAP[:7] + ["area J, edited by me\n"] + ROADMAP[8:]
    assert "the jobs table\n" in out  # the `## Data jobs` section after it is untouched


def test_a_top_level_section_still_swallows_its_own_subsections():
    work = list(ROADMAP)
    work[3] = "area A, edited by me\n"
    out = ss.splice(ROADMAP, work, "## The areas")
    assert "area A, edited by me\n" in out  # `## The areas` legitimately contains every `###` area
    assert out[-2:] == ROADMAP[-2:]  # and stops at `## Data jobs`


def test_a_hash_at_the_start_of_a_line_in_a_code_block_is_not_a_heading():
    base = [
        "## Mine\n",
        "```\n",
        "# explicit paths only\n",
        "git add x\n",
        "```\n",
        "the line after the block\n",
        "## Theirs\n",
        "their text\n",
    ]
    work = list(base)
    work[5] = "the line after the block, edited\n"
    out = ss.splice(base, work, "## Mine")
    assert "the line after the block, edited\n" in out
    assert out[:5] == base[:5] and out[6:] == base[6:]


def test_a_heading_argument_that_is_not_a_heading_is_refused():
    with pytest.raises(ss.RefusedError, match="not a markdown heading"):
        ss.splice(ROADMAP, ROADMAP, "area E, committed text")
