# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area H, 2026-10-10: the registration that decides whether BioForge's design answers have a calibration
record, pinned before the census that applies it (genomeos/forge/calibration_inventory.py)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from genomeos.forge.calibration_inventory import (
    MAX_INTERVAL_WIDTH,
    MINIMUM_INDEPENDENT_TRIPLES,
    REGISTRATION,
    clopper_pearson,
    minimum_triples,
    worst_width,
)


def test_the_exact_interval_reproduces_the_published_numbers() -> None:
    """The same [0.3976, 1] and [0.1581, 1] that docs/BIOFORGE-CONFIDENCE.md §4(c) printed at 4/4 and 2/2."""
    assert clopper_pearson(4, 4) == pytest.approx((0.3976, 1.0), abs=5e-5)
    assert clopper_pearson(2, 2) == pytest.approx((0.1581, 1.0), abs=5e-5)
    lo, hi = clopper_pearson(19, 38)
    assert lo == pytest.approx(1 - hi, abs=1e-9) and 1 / 3 < lo < hi < 2 / 3
    with pytest.raises(ValueError):
        clopper_pearson(3, 2)


def test_the_minimum_is_the_smallest_n_whose_worst_interval_is_narrower_than_a_third() -> None:
    assert MAX_INTERVAL_WIDTH == 1 / 3
    assert minimum_triples() == MINIMUM_INDEPENDENT_TRIPLES == 38
    assert worst_width(38) < 1 / 3 <= worst_width(37)
    assert REGISTRATION["the_minimum"]["n"] == MINIMUM_INDEPENDENT_TRIPLES


def test_the_registration_names_every_circularity_rule_and_refuses_a_rate_on_circular_triples() -> None:
    assert set(REGISTRATION["circular_if_any"]) == {
        "same_source",
        "acceptance_check",
        "authored_frame",
        "evidence_is_the_model",
    }
    assert "compute no rate" in REGISTRATION["decision"]
    assert "never counted as an independent outcome" in REGISTRATION["never"]


# The census rules on toys, before the census is run on the repository (area H, 2026-10-10).
_TOY = """
module toy.blind
import bio.std.development
organism T {{ root: R; cell_type: Zygote; resolution: populations; factors: A, B }}
timer cycle {{ duration: 60 min; when: cell_type = Zygote }}
decision grow {{
  action: divide; when: cell_type = Zygote, A = present; fraction: 1.0
  evidence: experimental "{rule}"; confidence: 0.5
}}
design no_growth {{
  knockout_any_of: {frame}; target: count at 300 min = 1; until: 300 min
  evidence: {design_evidence}; confidence: 0.4
}}
"""


def _toy_record(rule: str, frame: str, design_evidence: str) -> dict:
    from genomeos.forge.calibration_inventory import design_record, model_evidence
    from genomeos.lang import parse
    from genomeos.organism.forge import run_design

    m = parse(_TOY.format(rule=rule, frame=frame, design_evidence=design_evidence))
    return design_record("toy.bio", m, run_design(m, m.designs[0], seed=0), model_evidence(m), {}, [])


def test_citation_keys_read_the_first_author_and_year_of_each_citation() -> None:
    from genomeos.forge.calibration_inventory import citation_keys, cites

    src = "Kemphues et al. 1988, Cell 52:311; van der Flier et al. 2009, Cell 136:903; a note, no year"
    assert citation_keys(src) == [("Kemphues", "1988"), ("van der Flier", "2009")]
    assert cites('"source": "Jenny et al. 2002, EMBO J 21:6338"', ("Jenny", "2002"))
    assert not cites("Jenny et al. 2003", ("Jenny", "2002"))


def test_a_blind_frame_and_an_outside_citation_make_a_toy_triple_independent() -> None:
    """Independence is reachable: the rules are not written so that nothing can pass them."""
    rec = _toy_record(
        "Smith et al. 2001", "A, B", 'experimental "Jones et al. 2010: without A nothing divides"'
    )
    assert rec["answer"] == "-A" and rec["circular"] == [] and rec["independent"]


@pytest.mark.parametrize(
    ("rule", "frame", "design_evidence", "fires"),
    [
        (
            "Jones et al. 2010 (A drives division)",
            "A, B",
            'experimental "Jones et al. 2010: without A nothing divides"',
            "same_source",
        ),
        (
            "Smith et al. 2001",
            "A",
            'experimental "Jones et al. 2010: without A nothing divides"',
            "authored_frame",
        ),
        ("Smith et al. 2001", "A, B", 'curated "the module requires A to divide"', "evidence_is_the_model"),
    ],
)
def test_each_rule_fires_on_its_own(rule: str, frame: str, design_evidence: str, fires: str) -> None:
    rec = _toy_record(rule, frame, design_evidence)
    assert fires in rec["circular"] and not rec["independent"]


def test_a_shipped_experiment_in_an_organism_program_is_always_an_acceptance_check() -> None:
    from genomeos.forge.calibration_inventory import experiment_record, model_evidence
    from genomeos.lang import parse

    src = _TOY.format(rule="Smith et al. 2001", frame="A, B", design_evidence='curated "x"').replace(
        "design no_growth",
        "stage S { from: 0 min }\nexperiment ko_a { knockout: A; until: 300 min;"
        ' assert: count at 300 min = 1; evidence: experimental "Jones et al. 2010" }\ndesign no_growth',
    )
    m = parse(src)
    rec = experiment_record("toy.bio", m, m.experiments[0], model_evidence(m), {})
    assert rec["published"] and rec["expressible_as_a_design_question"]
    assert rec["circular"] == ["acceptance_check"] and not rec["independent"]


# The census result, data/results/forge_calibration_inventory.json, written 2026-10-10 by
# scripts/forge_calibration_inventory.py on a clean worktree of e06b1ca (lane-forgecal).
RESULT = Path(__file__).resolve().parent.parent / "data" / "results" / "forge_calibration_inventory.json"


def test_the_committed_inventory_is_what_the_code_computes() -> None:
    """Shares the design run with tests/test_forge_design_answers.py through design_answers.results_for."""
    from genomeos.forge.calibration_inventory import build_inventory, census_of

    written = json.loads(RESULT.read_text(encoding="utf-8"))
    assert census_of(written) == json.loads(json.dumps(build_inventory()))
    assert written["result_manifest"]["complete"] and written["money"].startswith("none")


def test_no_triple_is_independent_so_no_calibration_is_registered() -> None:
    """Exact counts of 2026-10-10 (e06b1ca). Every published outcome in the repository is also an acceptance
    check of the model it would be scored on, so the shortfall is the whole registered minimum."""
    inv = json.loads(RESULT.read_text(encoding="utf-8"))
    c = inv["counts"]
    assert (c["design_blocks"], c["design_triples_independent"], c["design_answers_tied"]) == (4, 0, 1)
    assert c["design_triples_circular_by_rule"] == {
        "same_source": 4,
        "acceptance_check": 3,
        "authored_frame": 4,
        "evidence_is_the_model": 2,
    }
    assert (c["experiment_blocks"], c["experiments_with_a_published_citation"]) == (39, 37)
    assert (c["expressible_as_inverse_design_questions"], c["set_aside_published_experiments"]) == (27, 10)
    assert c["inverse_triples_circular_by_rule"]["same_source"] == 27
    assert c["inverse_triples_circular_by_rule"]["acceptance_check"] == 27
    assert (c["independent_triples"], c["minimum_independent_triples"], c["shortfall"]) == (0, 38, 38)
    assert inv["decision"].startswith("shortfall of 38: no calibration is registered")
    # the three design answers that match a recorded outcome match the model's own acceptance experiment
    matched = [d for d in inv["design_triples"] if d["answer_matches_the_recorded_outcome"]]
    assert len(matched) == 3 and all("acceptance_check" in d["circular"] for d in matched)
