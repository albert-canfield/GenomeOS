# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registered facts behind genomeos/forge/calibration.py's refusal.

The refusal rests on two claims that a later edit could quietly falsify: that the emitted
confidence is single-valued, and that the repository holds no design whose answer is unknown to
the module it runs over. Each is pinned here, so that making BioForge's confidence calibratable
breaks these tests rather than leaving docs/BIOFORGE-CONFIDENCE.md silently wrong.
"""

from __future__ import annotations

import re
from functools import lru_cache
from pathlib import Path

from genomeos.forge.calibration import (
    CEILING_SITES,
    CENSUS,
    PREREGISTRATION,
    RETIRED_SITES,
    THE_CONSTANT,
    VERDICT,
)
from genomeos.lang import parse_file
from genomeos.organism.forge import UNSTATED_CONFIDENCE, run_designs

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


def test_the_emitted_confidence_is_no_longer_the_constant_and_states_no_probability() -> None:
    """Review R4 retired the 0.3 on the emitted experiment. The old pin asserted it was still there;
    this one asserts the replacement: no number is stated, and no probability is claimed."""
    per_design = _run()
    assert len(per_design) == CENSUS["design_blocks_in_data"] == 4
    confidences = {v[3] for v in per_design.values()}
    assert confidences == {UNSTATED_CONFIDENCE} and THE_CONSTANT not in confidences


def test_the_confidence_ignores_everything_the_search_computes() -> None:
    """Designs that differ by 4x in evaluations and 11x in feasible candidates get the same number."""
    per_design = _run()
    evaluations = {v[0] for v in per_design.values()}
    feasible = {v[1] for v in per_design.values()}
    assert len(evaluations) > 1 and len(feasible) > 1, "the inputs no longer vary; the pin is vacuous"
    assert {v[3] for v in per_design.values()} == {UNSTATED_CONFIDENCE}
    # the design's own stated confidence is parsed and still not carried onto the answer
    stated = {v[4] for v in per_design.values()}
    assert stated != {UNSTATED_CONFIDENCE} and len(stated) > 1


def test_every_outcome_is_positive_so_there_is_no_negative_case() -> None:
    """A constant outcome column is the second fatal degeneracy, independent of the first."""
    solved = {v[5] for v in _run().values()}
    assert solved == {True}, (
        "a design now fails, which is the first negative case in the repository. That is what"
        " obstruction 1 in WHAT_WOULD_CHANGE_IT asked for: area H's third item can be re-opened"
    )


def test_the_retired_literal_is_gone_and_the_ceilings_are_where_the_docs_say() -> None:
    """Pinned by pattern, not by line: the literal that created the number is gone, and the ceilings
    that only lower an evidence-quality score are counted per file."""
    root = Path(__file__).resolve().parent.parent
    for site in RETIRED_SITES:
        rel = site.rpartition(":")[0]
        assert "confidence=0.3" not in (root / rel).read_text().replace(" ", "")
    ceiling = re.compile(r"min\([\w.]*confidence, 0\.3\)")
    for rel, n in CEILING_SITES.items():
        assert len(ceiling.findall((root / rel).read_text())) == n, rel


def test_the_registration_and_the_verdict_are_present_and_agree() -> None:
    assert VERDICT == "UNCHECKABLE_BY_CONSTRUCTION"
    assert PREREGISTRATION["registered_expectation"].startswith("refusal")
    for key in ("what_would_count_as_calibrated", "what_would_count_as_uncheckable", "baselines"):
        assert PREREGISTRATION[key]
    assert PREREGISTRATION["the_population_clause"]["quoted_for"]
    assert CENSUS["fixtures_pairing_a_design_with_a_published_outcome"] == 0
