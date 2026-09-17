# SPDX-License-Identifier: AGPL-3.0-or-later
"""Section 5, the one ordered list across every area, as the Progress tab reads it."""

from __future__ import annotations

from pathlib import Path

from genomeos import roadmap

TEXT = """# GenomeOS

Reviewed and reorganised on 2026-09-17 against the code.

## 5. Next steps, consolidated

**Where things stand.** A paragraph that is not a numbered item.

1. ~~The endpoint that ran.~~ **Done 2026-09-17 (`2a07c80`): it returned +0.075.**
2. **The panel's chr2.** Decides the depletion at full size. No key.
3a. ~~The 69 syntax blocks.~~ **Cancelled unrun.**
3b. **Measurement over the real unknown.** 721 of 882 blocks hold nothing measured,
   which is the binding constraint.
4. **Disk.** The project is 7.3 GB of a 460 GB disk.

The list below is the consolidated order as it stood on 2026-09-11, kept as history.

1. Whole-genome data jobs: done on 2026-09-11.
2. An older item nobody should see in the tab.

## 6. Milestones
"""


def test_the_live_list_is_parsed_and_the_history_below_it_is_not() -> None:
    steps = roadmap.parse_next_steps(TEXT)

    assert [s["number"] for s in steps] == ["1", "2", "3a", "3b", "4"]


def test_a_struck_through_item_is_done_and_the_rest_are_next() -> None:
    by_number = {s["number"]: s for s in roadmap.parse_next_steps(TEXT)}

    assert by_number["1"]["state"] == "done"
    assert by_number["3a"]["state"] == "done"
    assert by_number["2"]["state"] == "next"
    assert by_number["3b"]["state"] == "next"


def test_the_whole_entry_is_kept_because_the_reasons_are_the_point() -> None:
    """The tab shows a title and hides the rest behind a disclosure; the rest must exist."""
    step = {s["number"]: s for s in roadmap.parse_next_steps(TEXT)}["3b"]

    assert "721 of 882" in step["text"]
    assert "binding constraint" in step["text"]
    assert step["title"].startswith("Measurement over the real unknown")


def test_a_title_is_the_first_sentence_without_its_markup() -> None:
    step = {s["number"]: s for s in roadmap.parse_next_steps(TEXT)}["4"]

    assert step["title"] == "Disk."
    assert "**" not in step["title"]


def test_a_section_without_a_history_marker_keeps_every_item() -> None:
    text = TEXT.replace(
        "The list below is the consolidated order as it stood on 2026-09-11, kept as history.", ""
    )

    assert len(roadmap.parse_next_steps(text)) == 7


def test_the_real_roadmap_parses_and_every_item_carries_a_state() -> None:
    steps = roadmap.load(Path("."))["next_steps"]

    assert steps, "docs/ROADMAP.md section 5 should have items"
    assert all(s["state"] in ("done", "next") for s in steps)
    assert all(s["number"] and s["text"] for s in steps)


def test_a_step_that_says_done_in_bold_is_done_wherever_it_says_it() -> None:
    """The lanes write "PAR polarity rules: **done 2026-09-14**" as often as they start with "Done".

    Reading only the opening word left area E showing nine planned steps of which five were finished,
    and the Progress tab overstating what is left is worse than understating it: that tab is what gets
    read to decide what to do next.
    """
    assert roadmap.step_state("PAR polarity rules: **done 2026-09-14** (above), and 11 of 11.") == "done"
    assert roadmap.step_state("Glia from factors: **done 2026-09-14**; the rest waits.") == "done"


def test_a_date_is_required_so_prose_about_being_done_does_not_count() -> None:
    assert roadmap.step_state("Something that is **not done** yet, 2026-09-14") == "planned"
    assert roadmap.step_state("A measured contact table for the embryo.") == "planned"
    assert roadmap.step_state("What **would be done** on 2026-09-14 if the data existed") == "planned"


def test_the_open_count_excludes_the_steps_it_knows_are_finished() -> None:
    """It used to count every step in the Next section, finished or not."""
    text = """# R

## 3. Areas

### E. From one cell to an organism

- **Goal.** A worm.
- **Next.** 1. PAR polarity rules: **done 2026-09-14** (above). 2. A measured contact table.
  3. Glia from factors: **done 2026-09-14**.

---
"""
    area = next(a for a in roadmap.parse_areas(text) if a["letter"] == "E")

    assert len(area["next"]) == 3
    assert area["counts"]["next"] == 1
    assert area["counts"]["steps_done"] == 2
