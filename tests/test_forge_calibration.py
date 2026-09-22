# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registered facts behind genomeos/forge/calibration.py's refusal.

The refusal rests on two claims that a later edit could quietly falsify: that the emitted
confidence is single-valued, and that the repository holds no design whose answer is unknown to
the module it runs over. Each is pinned here, so that making BioForge's confidence calibratable
breaks these tests rather than leaving docs/BIOFORGE-CONFIDENCE.md silently wrong.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from genomeos.forge.calibration import (
    CENSUS,
    CONFIDENCE_SITES,
    PREREGISTRATION,
    THE_CONSTANT,
    VERDICT,
)
from genomeos.lang import parse_file
from genomeos.organism.forge import run_designs

DESIGN_FILES = [Path(p) for p in CENSUS["design_block_files"]]


@lru_cache(maxsize=1)
def _run() -> dict[str, tuple[int, int, float | None, float | None, float, bool]]:
    """Every design in the repository, run once: the search is the expensive part of these tests."""
    out: dict[str, tuple[int, int, float | None, float | None, float, bool]] = {}
    for path in DESIGN_FILES:
        for r in run_designs(parse_file(str(path)), None, seed=0):
            ex = r.to_experiment()
            out[r.design.name] = (
                r.evaluations,
                sum(1 for c in r.candidates if c.feasible),
                r.best.loss if r.best else None,
                ex.confidence if ex else None,
                r.design.confidence,
                r.solved,
            )
    return out


def test_the_emitted_confidence_is_single_valued() -> None:
    """One distinct value means one reliability bin, which is why no curve exists at any n."""
    per_design = _run()
    assert len(per_design) == CENSUS["design_blocks_in_data"] == 4
    confidences = {v[3] for v in per_design.values()}
    assert confidences == {THE_CONSTANT}, (
        f"the calibration refusal assumes a degenerate predictor; it now takes {confidences}."
        " If the confidence has been made to vary, re-open area H's third item and re-read"
        " docs/BIOFORGE-CONFIDENCE.md, whose verdict no longer follows"
    )


def test_the_confidence_ignores_everything_the_search_computes() -> None:
    """Designs that differ by 4x in evaluations and 11x in feasible candidates get the same number."""
    per_design = _run()
    evaluations = {v[0] for v in per_design.values()}
    feasible = {v[1] for v in per_design.values()}
    assert len(evaluations) > 1 and len(feasible) > 1, "the inputs no longer vary; the pin is vacuous"
    assert {v[3] for v in per_design.values()} == {THE_CONSTANT}
    # the design's own stated confidence is parsed and then discarded by to_experiment()
    stated = {v[4] for v in per_design.values()}
    assert stated != {THE_CONSTANT} and len(stated) > 1, (
        "the shipped designs state 0.5 and 0.6; if they no longer differ from the emitted 0.3 the"
        " 'stated confidence is discarded' claim can no longer be demonstrated from the data"
    )


def test_every_outcome_is_positive_so_there_is_no_negative_case() -> None:
    """A constant outcome column is the second fatal degeneracy, independent of the first."""
    solved = {v[5] for v in _run().values()}
    assert solved == {True}, (
        "a design now fails, which is the first negative case in the repository. That is what"
        " obstruction 1 in WHAT_WOULD_CHANGE_IT asked for: area H's third item can be re-opened"
    )


def test_the_constant_is_still_written_where_the_docs_say_it_is() -> None:
    """Seven sites, one literal. The doc quotes these line numbers; they must stay true."""
    root = Path(__file__).resolve().parent.parent
    for site in CONFIDENCE_SITES:
        rel, _, lineno = site.rpartition(":")
        line = (root / rel).read_text().splitlines()[int(lineno) - 1]
        assert str(THE_CONSTANT) in line, f"{site} no longer writes {THE_CONSTANT}: {line.strip()!r}"


def test_the_registration_and_the_verdict_are_present_and_agree() -> None:
    assert VERDICT == "UNCHECKABLE_BY_CONSTRUCTION"
    assert PREREGISTRATION["registered_expectation"].startswith("refusal")
    for key in ("what_would_count_as_calibrated", "what_would_count_as_uncheckable", "baselines"):
        assert PREREGISTRATION[key]
    assert PREREGISTRATION["the_population_clause"]["quoted_for"]
    assert CENSUS["fixtures_pairing_a_design_with_a_published_outcome"] == 0
