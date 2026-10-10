# SPDX-License-Identifier: AGPL-3.0-or-later
"""Area H, 2026-10-10: the CI path for the shipped `design` programs.

`bio test` runs experiment blocks and asserts, never design blocks, so until this file a change to a
module could move a design's answer unseen (haematopoiesis_designs.bio's two answers were pinned by no
test; only `solved` was). This file runs every design block under the CI `bio test` roots through
BioForge (genomeos/forge/design_answers.py) and holds each answer, its ties and its search effort to the
record below. It checks stability only: the record says what BioForge answers, not that the answer is
right, and an answer that reproduces its own model's acceptance phenotype is not a prediction.

To change a pin, change the module or the design on purpose, rerun, and write the new record here with
the date and the commit that moved it.
"""

from __future__ import annotations

from genomeos.forge import design_answers
from genomeos.forge.design_answers import (
    NOT_A_CLAIM,
    answer_records,
    design_programs,
    results_for,
)
from genomeos.lang import parse
from genomeos.organism.forge import run_design

CELEGANS = "data/organisms/celegans/designs.bio"
BLOOD = "data/organisms/human/haematopoiesis_designs.bio"

# Pinned 2026-10-10 by lane-forgecal at seed 0 on the tree of the commit that adds this file.
PINNED = [
    {
        "program": CELEGANS,
        "design": "two_intestinal_founders",
        "answer": "-POP-1",
        "solved": True,
        "loss": 0.0,
        "tied_with_the_answer": [],
        "tie_broken_by": None,
        "evaluations": 7,
        "feasible_candidates": 2,
        "varies_knobs": False,
    },
    {
        "program": CELEGANS,
        "design": "two_ems_founders_no_germline",
        "answer": "-PIE-1",
        "solved": True,
        "loss": 0.0,
        "tied_with_the_answer": [],
        "tie_broken_by": None,
        "evaluations": 6,
        "feasible_candidates": 6,
        "varies_knobs": False,
    },
    {
        "program": BLOOD,
        "design": "lose_neutrophils_keep_the_rest",
        "answer": "-GFI1",
        "solved": True,
        "loss": 0.0,
        "tied_with_the_answer": [],
        "tie_broken_by": None,
        "evaluations": 11,
        "feasible_candidates": 2,
        "varies_knobs": False,
    },
    {
        # -IKZF1 and -TCF3 both reach the targets at loss 0 with one knockout; the answer is IKZF1 only
        # because the final sort falls back to the label. The pin holds the tie so a rename cannot flip it
        # unseen.
        "program": BLOOD,
        "design": "no_lymphocytes_at_all",
        "answer": "-IKZF1",
        "solved": True,
        "loss": 0.0,
        "tied_with_the_answer": ["-TCF3"],
        "tie_broken_by": "label order (forge.run_design's final sort)",
        "evaluations": 29,
        "feasible_candidates": 22,
        "varies_knobs": False,
    },
]


def test_every_design_block_under_the_ci_roots_is_pinned() -> None:
    """A new design block anywhere bio test looks fails here until its answer is pinned."""
    assert design_programs() == [CELEGANS, BLOOD]
    found = [(p, r.design.name) for p in design_programs() for r in results_for(p)]
    assert found == [(r["program"], r["design"]) for r in PINNED]


def test_every_shipped_design_gives_its_pinned_answer() -> None:
    assert answer_records() == PINNED


def test_the_record_claims_no_accuracy() -> None:
    forbidden = ("accuracy", "probability", "correct", "match", "published", "calibrat")
    for rec in answer_records():
        assert not [k for k in rec if any(word in k for word in forbidden)], rec
    assert "not a check that the answer is right" in NOT_A_CLAIM


def test_a_knob_free_answer_does_not_depend_on_the_seed() -> None:
    """The seed reaches run_design in two places: the knob search's draws and Body(seed, means=True),
    which makes no random generator. A design with no `vary` knob therefore cannot depend on the seed;
    shown here end to end on the cheap C. elegans designs, and on a toy whose timers would spread with a
    seed if `means` did not hold them at their means. A shipped design that gains a knob is re-run."""
    for p in design_programs():
        for r in results_for(p):
            if r.design.vary or p == CELEGANS:
                again = design_answers.answer_record(
                    p, next(x for x in results_for(p, 1) if x.design == r.design)
                )
                assert again == design_answers.answer_record(p, r)
    src = """
module toy.design
import bio.std.development
organism T { root: R; cell_type: Zygote; resolution: populations }
timer cycle { duration: 60 min; sd: 20 min; when: cell_type = Zygote }
decision grow { action: divide; when: cell_type = Zygote; fraction: 1.0 }
design pace { target: count at 300 min in 30..34; until: 300 min }
"""
    m = parse(src)

    def losses(means: bool) -> set:
        runs = [run_design(m, m.designs[0], seed=seed, means=means) for seed in (0, 1, 7)]
        return {tuple((c.label(), c.loss) for c in r.candidates) for r in runs}

    assert len(losses(means=True)) == 1  # what BioForge runs: the seed changes nothing
    assert len(losses(means=False)) > 1  # and the toy is not vacuous: with spread timers it would
