# SPDX-License-Identifier: AGPL-3.0-or-later
"""The knockdown-field documentation census: its instrument, and its counterfactuals.

A census whose whole output is the word NOT ESTABLISHED rests on counts of zero, and a zero from a
broken instrument reads exactly like a zero from an absent field. So the tests here prove BOTH
directions: that the counter finds a field name when one is there, that the quote finder fails when
the sentence is not there, and that the verdict follows the counts rather than the prose around them.

Nothing here opens a study file or reads a knockdown value, and nothing here needs the network: the
one test that wants a real producer document skips when the cache has not been filled by a run.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _module():
    spec = importlib.util.spec_from_file_location(
        "n1_kd_semantics_census", ROOT / "scripts" / "n1_kd_semantics_census.py"
    )
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def test_the_counter_finds_a_field_name_that_is_present() -> None:
    """The positive control. `occurrences` is the instrument every NOT ESTABLISHED verdict rests on,
    so it must be shown to return a nonzero count for a name that IS in the text."""
    m = _module()
    text = "mean_pop.cells['fold_expr'] = x  # and again: fold_expr"
    assert m.occurrences(text, ("fold_expr", "pct_expr")) == {"fold_expr": 2, "pct_expr": 0}


def test_the_counter_finds_the_2019_names_in_the_real_producer_archive() -> None:
    """The non-vacuity proof against the actual document. The census reports zero for the four 2022
    spellings in `Perturbseq_GI`; the same instrument over the same archive must find the 2019
    family, or the zeros would be evidence of nothing."""
    m = _module()
    archive = m.CACHE / m.DOCS["code_perturbseq_gi"]["file"]
    if not archive.exists():
        pytest.skip("the producer code archive is not cached; run the census script to fetch it")
    text = m.archive_text(archive)
    found = m.occurrences(text, ("pct_first_expr", "control_first_expr", "fold_first_expr"))
    assert all(v > 0 for v in found.values()), found
    # and the 2022 spellings really are absent from the same bytes
    assert m.occurrences(text, m.FIELDS_2022) == dict.fromkeys(m.FIELDS_2022, 0)


def test_the_quote_finder_fails_on_a_sentence_that_is_not_there() -> None:
    """The negative control for the other instrument: a producer sentence the census requires must
    not be reported present when the document does not carry it."""
    m = _module()
    assert m.find_quote("Knockdown was computed as the ratio of mean expression", "ratio of mean")
    assert not m.find_quote("nothing about knockdown here", m.QUOTES["knockdown_definition"])


def test_the_quote_finder_tolerates_line_breaks_but_not_missing_words() -> None:
    m = _module()
    assert m.find_quote("an on-target\n  knockdown of at\tleast 60%", m.QUOTES["producer_threshold_60"])
    assert not m.find_quote("an on-target knockdown of at least", m.QUOTES["producer_threshold_60"])


def _docs(**over: str) -> dict[str, str]:
    base = dict.fromkeys(
        (
            "pmc9380471",
            "figshare_20029387",
            "figshare_21632564",
            "code_perturbseq_gi",
            "code_guide_calling",
        ),
        "",
    )
    base.update(over)
    return base


def test_the_verdict_follows_the_counts_not_the_prose() -> None:
    """The counterfactual on the census itself: plant a producer document that DOES name a column and
    the field must stop being reported as absent everywhere. The verdict is driven by the evidence."""
    m = _module()
    clean = m.census(_docs())
    assert clean["fields_absent_from_every_producer_document"] == list(m.FIELDS_2022)
    assert [f["verdict"] for f in clean["fields"]] == ["NOT ESTABLISHED"] * len(m.FIELDS_2022)

    planted = m.census(_docs(figshare_20029387="fold_expr: the ratio of perturbed to control mean"))
    assert "fold_expr" not in planted["fields_absent_from_every_producer_document"]
    by_name = {f["field"]: f for f in planted["fields"]}
    assert by_name["fold_expr"]["counts"]["deposit_descriptions"] == 1
    assert by_name["pct_expr"]["counts"]["deposit_descriptions"] == 0


def test_the_census_never_chooses_a_threshold() -> None:
    """The lane's hard constraint, as an assertion: the census reports what the producers chose and
    chooses nothing itself."""
    m = _module()
    body = m.census(_docs())
    assert body["answer"]["not_chosen_here"].startswith("no threshold is chosen")
    levels = [t["level"] for t in body["quantity"]["producers_own_thresholds"]]
    assert levels == [
        "at least 30% on-target knockdown, or the target not detected",
        "at least 60% on-target knockdown",
    ]
    assert body["answer"]["verdict"].startswith("NO")


def test_the_quantity_is_established_and_the_columns_are_not() -> None:
    """The distinction the census exists to keep, stated as a test so a later edit cannot blur it."""
    m = _module()
    body = m.census(_docs())
    assert body["quantity"]["verdict"] == "ESTABLISHED"
    assert all(f["verdict"] == "NOT ESTABLISHED" for f in body["fields"])
    assert body["correction_to_the_committed_record"]["verdict"] == "NOT ESTABLISHED for the 2022 file"
