# SPDX-License-Identifier: AGPL-3.0-or-later
"""Can BioForge's predicted confidence be calibrated at all? The registration and the refusal.

Area H's third open item asks "how often does a design's answer match the published outcome",
which presumes the number beside the answer is a probability. This module establishes what the
number is, counts what could be checked against it, registers the comparison before running it,
and records the verdict. The verdict is a refusal, and the refusal is the result.

**What the number is.** `organism/forge.py:143` writes `confidence=0.3` as a literal onto the
`Experiment` a design emits. It is not derived from `loss`, `feasible`, `evaluations`, the number
of tied feasible candidates or the margin to the runner-up, all of which the search computes and
then discards. Six further sites cap an existing confidence at the same constant
(`organism/forge.py:66,70`, `forge/design.py:59,63,66`, and `attribution/budget.py:60` for an
unrelated quantity). A cap can only lower a number and is defensible as a ceiling; the literal at
line 143 is different in kind, because it creates a number out of nothing and attaches it to a
freshly predicted entity, where a reader is most likely to read it as a posterior. The design
block's own stated confidence (0.5 and 0.6 in the shipped programs) is parsed, stored on the
`Design`, and then ignored by `to_experiment()`.

So the emitted confidence is SINGLE-VALUED by construction. That is not a small-sample problem: a
reliability diagram needs at least two populated bins, and a predictor with one distinct value has
one bin at every n. **No calibration curve for this quantity can exist at any sample size**, and
that conclusion is read off the source, not off the data.

**What was left to check.** With the curve ruled out, one weaker statement survives: a single base
rate, "among the answers BioForge emits, what fraction match the published outcome, and does 0.3
lie inside its interval". `CENSUS` counts what that could be run on. The answer is 2 designs whose
published outcome is recorded in the program file with a citation, or 4 if two designs whose
answers appear only in prose are admitted. Both populations are named below and neither is pooled.

**Why the base rate was computed and then refused.** At 4 of 4 the exact 95% interval is
[0.3976, 1.0000] and excludes 0.3, so the arithmetic says the stated confidence is too low and
invites raising it. That number is not published here, because the population it is computed on is
the population in which the answer was placed before the search ran. `THE_GATE` records the
proof: the design block for `two_intestinal_founders` states its published outcome as *"Lin et al.
1995: without POP-1, MS takes the E fate"*, and the rule the search must traverse to reach that
outcome, `fate_E` in `data/organisms/celegans/founders.bio:47`, carries the evidence string *"Lin
et al. 1995 (without POP-1, MS takes the E fate)"* - the same clause of the same paper, hand
written at confidence 0.9. The search is not predicting POP-1; it is evaluating a rule authored
from the answer, over a candidate list into which the author put the answer. 4 of 4 measures that
the forward simulator reads its own rules correctly. It carries no information about a design
question BioForge has not been given the answer to.

This is the fifth instance this month of the class named in docs/LESSONS.md on 2026-09-22: a curve
is honest on the population it was fitted on and nowhere else, and the gate that chose that
population is usually invisible. Here the gate is a hand-written candidate list and a rule copied
from the citation, and the quantity under test was never a curve to begin with.
"""

from __future__ import annotations

from typing import Any

# The constant, and every site that writes it. Read off the source on 2026-09-22.
THE_CONSTANT = 0.3
CONFIDENCE_SITES: dict[str, str] = {
    "genomeos/organism/forge.py:143": "confidence=0.3 - CREATES the number on the emitted Experiment;"
    " the only one a user of `genomeos grow --design` reads",
    "genomeos/organism/forge.py:66": "min(timer.confidence, 0.3) - ceiling on a knob the search moved",
    "genomeos/organism/forge.py:70": "min(decision.confidence, 0.3) - ceiling",
    "genomeos/forge/design.py:59": "min(gene.confidence, 0.3) - ceiling",
    "genomeos/forge/design.py:63": "min(rule.confidence, 0.3) - ceiling",
    "genomeos/forge/design.py:66": "min(parameter.confidence, 0.3) - ceiling",
    "genomeos/attribution/budget.py:60": "min(confidence, 0.3) - the same literal for an unrelated"
    " quantity in another area, which is what a provenance marker looks like and what a probability"
    " does not",
}
WHAT_IT_IS = (
    "a provenance marker, not a probability: one hand-set constant meaning 'this came from a search,"
    " not from a bench', used at the same value by an unrelated area for an unrelated quantity"
)
DISCARDED_BY_THE_CONFIDENCE = (
    "loss",
    "feasible",
    "evaluations",
    "number of feasible candidates found",
    "margin between the best and the second-best candidate",
    "the design block's own stated confidence",
)

# Counted before anything was designed or measured. 2026-09-22.
CENSUS: dict[str, Any] = {
    "design_blocks_in_data": 4,
    "design_block_files": (
        "data/organisms/celegans/designs.bio",
        "data/organisms/human/haematopoiesis_designs.bio",
    ),
    "design_blocks_in_tests_as_toys": 2,
    "with_published_outcome_cited_in_the_file": 2,
    "with_answer_recorded_only_in_prose": 2,
    "fixtures_pairing_a_design_with_a_published_outcome": 0,
    "tests_that_read_a_published_outcome_and_compare_it_to_a_design_answer": 0,
    "adjacent_but_not_design_shaped": (
        "31 `experiment` blocks in data/, 29 of them with a published citation, and 58 `# test:`"
        " lines across 21 of 65 .bio files. These are perturbation/outcome pairs, not design/outcome"
        " pairs: an experiment states the perturbation and checks the simulator, while a design must"
        " FIND the perturbation. They cannot stand in for the missing population, because in an"
        " experiment there is no answer for BioForge to get wrong"
    ),
}

PREREGISTRATION: dict[str, Any] = {
    "written": (
        "2026-09-22, after the four design blocks were counted and run and after the emitted"
        " confidence was read off the source, and BEFORE any match rate was compared with 0.3 and"
        " before any interval was computed. The census and the degeneracy are inputs to the design of"
        " the check, not outcomes of it; what is registered here is what the match rate would have to"
        " look like to be publishable"
    ),
    "the_quantity_under_test": (
        "the confidence `organism/forge.py:143` attaches to the predicted Experiment a design emits,"
        " as read by a user of `genomeos grow --design`"
    ),
    "what_is_compared": (
        "for each design block with a published outcome, the perturbation BioForge returns as its"
        " best feasible candidate, against the perturbation the cited publication reports. A match is"
        " exact set equality of the knockout set; anything else is a miss"
    ),
    "baselines": {
        "trivial_modal": (
            "always predict the modal outcome. Every design in the repository is solved, so the modal"
            " outcome is 'the answer is right' and a predictor that answers 'right' unconditionally"
            " scores the same as BioForge. Registered as the baseline BioForge must BEAT, not tie"
        ),
        "uniform_from_the_candidate_list": (
            "pick uniformly from the design's own knockout_any_of at its own at_most. Per design that"
            " is 1/6, 1/5, 1/10 and 1/28; all four by chance is 1/8400"
        ),
        "constant_confidence": (
            "quote 0.3 for everything, which is what the code already does. The calibration question"
            " is whether ANY predictor here separates from this one, and a single-valued predictor"
            " cannot"
        ),
    },
    "what_would_count_as_calibrated": (
        "(a) the emitted confidence takes at least two distinct values over the designs, so that at"
        " least two reliability bins are populated, AND (b) within each populated bin the observed"
        " match rate's 95% interval contains the bin's stated confidence, AND (c) the ordering of the"
        " bins agrees with the ordering of their observed rates. All three, or the word calibrated is"
        " not used"
    ),
    "what_would_count_as_uncheckable": (
        "any one of these, each of which is fatal on its own: (a) the predictor is single-valued, so"
        " there is one bin and no curve exists at ANY n; (b) the outcome is single-valued, so the"
        " label set has no negative case and no rate can be estimated from it; (c) the 95%"
        " Clopper-Pearson interval on the observed rate is wider than 0.5 at the available n; (d) the"
        " candidate list or the module rules were authored from the same publication the answer is"
        " scored against, in which case the search retrieves rather than predicts and the population"
        " is not the one the confidence is quoted for"
    ),
    "the_population_clause": {
        "quoted_for": (
            "EVERY feasible answer BioForge emits, for any design block, over any organism program,"
            " for any target. The constant is unconditional: nothing in the code narrows it"
        ),
        "checkable_against": (
            "the 2 design blocks whose published outcome is cited in the program file, or the 4 if"
            " the two whose answers live only in prose are admitted. Reported separately and never"
            " pooled, per the 2026-09-17 entry in docs/LESSONS.md"
        ),
        "and_these_are_not_the_same_population": (
            "the second is a hand-built subset of the first, selected by an author who held the"
            " published answer while writing both the candidate list and the module rule that leads"
            " to it. A rate measured on it is not a rate for the population the number is quoted for"
        ),
    },
    "registered_expectation": (
        "refusal. The predictor was already known to be a literal before this was written, so (a)"
        " under what_would_count_as_uncheckable was expected to fire. What was NOT known in advance"
        " is whether the surviving one-bin base rate would look publishable, and the registration"
        " exists so that a tempting number could not be promoted after the fact"
    ),
}

THE_GATE: dict[str, Any] = {
    "claim": (
        "the published outcome a design is scored against is the same sentence, from the same paper,"
        " as the evidence on the module rule the search traverses to produce it"
    ),
    "worked_example": {
        "design": "two_intestinal_founders, data/organisms/celegans/designs.bio:14",
        "its_stated_outcome": 'curated "Lin et al. 1995: without POP-1, MS takes the E fate"',
        "the_rule_it_must_traverse": "fate_E, data/organisms/celegans/founders.bio:47",
        "that_rule_s_evidence": 'experimental "Lin et al. 1995 (without POP-1, MS takes the E fate);'
        ' Maduro et al. 2005 (END-1/END-3)", confidence 0.9',
        "so": "the answer is hand written into the model at confidence 0.9 and the search is scored"
        " for returning it. The same holds for GFI1, whose requirement in GMP_to_Myeloblast cites"
        " Hock et al. 2003 (GFI1: neutropenia), which is the phenotype the design asks for",
    },
    "and_the_candidate_list": (
        "every knockout_any_of contains the published factor, because the author put it there. A"
        " design whose answer is absent from its own candidate list does not exist in this repository"
    ),
}

VERDICT = "UNCHECKABLE_BY_CONSTRUCTION"
VERDICT_TEXT = (
    "BioForge's predicted confidence cannot be calibrated, and the obstruction is not the sample"
    " size. It is a single hand-set constant, so it has one reliability bin at every n and no curve"
    " over it can exist. The weaker one-bin base rate is available at n=2 or n=4 and is withheld,"
    " because both populations were authored from the publications they are scored against. The"
    " honest statement is that the number has never been a probability: it is a provenance marker"
    " reading 'predicted, not validated', and it should be read and documented as one"
)

WHAT_WOULD_CHANGE_IT: tuple[str, ...] = (
    "1. A NEGATIVE CASE. Every design in the repository is solved and every answer is right, so the"
    " outcome column is constant and no rate can be estimated from it. At least one design whose"
    " published answer is 'no single perturbation achieves this', or whose candidate list does not"
    " contain the answer, is needed before the question has two sides.",
    "2. A BLIND CANDIDATE LIST. Fix knockout_any_of by a rule written before the design is chosen -"
    " for instance every factor the module declares - rather than by curation around the answer."
    " This is the loci-benchmark move: the frame is a function of the data, not of the result.",
    "3. DESIGNS WHOSE ANSWER IS NOT IN THE MODULE. So long as the rule that yields the outcome cites"
    " the paper the outcome comes from, the search retrieves. A held-out publication, whose finding"
    " was not used to author any rule, is the only kind of design that tests prediction.",
    "4. A CONFIDENCE THAT VARIES. The search already computes quantities that differ sharply across"
    " the four designs and throws all of them away: feasible candidates found (2, 6, 2, 22) and"
    " evaluations (7, 6, 11, 29). Deriving the confidence from one of those would give it more than"
    " one bin - which is necessary for calibration and nowhere near sufficient, since the derived"
    " number would then need calibrating on a population that does not yet exist. Item 1 comes first.",
    "5. n. With a varying predictor and both outcomes present, distinguishing a stated 0.3 from a"
    " true 0.9 by exact binomial test at 80% power needs n >= 6; a 95% interval narrower than the"
    " 0.6024 that n=4 gives at its best needs more. The repository holds 2 to 4. n is the smallest of"
    " the five obstructions and the only one that more of the same work would fix.",
)
