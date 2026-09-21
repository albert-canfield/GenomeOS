# SPDX-License-Identifier: Apache-2.0
# Part of the BioLang engine (language, IR, VM, standard library); see LICENSING.md.
"""Pre-registration for stage 2's gate (a), written before the instrument and before any run.

§5.3 states the gate: *"Expressing a costly protein must reduce the output of unrelated genes sharing
the pool, in the documented direction and order of magnitude: Ceroni et al. 2015 (Nat Methods 12:415)
in E. coli and Frei et al. 2020 (Nat Commun 11:4641) in mammalian cells. Falsified if no reduction
appears, or if it appears with the pool switched off."*

Two of those three clauses can be tested with what this project holds. The third cannot, and saying
so now — before a number exists — is the difference between a registered limitation and an excuse
found afterwards.

**What is askable.** Direction, monotonicity and the falsifier are properties of the model's own
arithmetic against its own declared capacities, so they need no external data. They are also the
clauses that can fail: an implementation that scaled nothing, or scaled everything regardless of
capacity, fails them.

**What is NOT askable, registered as unaskable rather than attempted.** "The documented order of
magnitude" requires Ceroni's and Frei's measured curves, and this project holds neither: there is no
burden dataset under `data/knowledge`, and nothing in the repository carries a capacity-monitor
reading or a mammalian resource-competition coefficient. A number invented for the comparison would
be a number invented to pass. So the magnitude clause is recorded as **not yet askable**, with what
it would take: the per-construct output of an unrelated reporter at two or more expression levels
from either paper, digitised with its axis stated.

This is the same distinction §18 drew for the model's 1 Mb reach — *a wrong answer and an unasked
question arrive in the same shape*, and only the second one is free to detect.
"""

from __future__ import annotations

from typing import Any

#: expression levels of the costly protein, spanning unconstrained to heavily constrained. Fixed
#: here so the curve is not chosen after seeing which shape it makes.
LEVELS: tuple[float, ...] = (1e3, 1e4, 1e5, 1e6, 5e6, 1e7, 5e7)

#: the unrelated gene's demand, held constant across every level: it is the thing that must fall
#: without its own demand changing, or the fall says nothing about sharing.
UNRELATED_DEMAND = 100.0

#: the falsifier arm's capacity multiplier. Large enough that no level above can exhaust the pool.
UNCONSTRAINED_MULTIPLIER = 1e6

PRE_REGISTRATION: dict[str, Any] = {
    "written": (
        "2026-09-21, before the instrument existed and before any level was run; committed in this "
        "file before the script that runs it"
    ),
    "claim": (
        "expressing a costly protein reduces the output of an unrelated gene that shares its pools, "
        "and does so only because the pool is the constraint"
    ),
    "design": (
        "one program, two genes and two proteins sharing a ribosome pool and an ATP pool "
        "(data/demo/stage2_pools.bio). The unrelated gene's own demand is held at "
        f"{UNRELATED_DEMAND} at every level, so any change in its output comes from sharing rather "
        f"than from itself. The costly protein is expressed at {len(LEVELS)} fixed levels, "
        f"{LEVELS[0]:.0e} to {LEVELS[-1]:.0e}, chosen before the run"
    ),
    "outcomes": {
        "direction": (
            "the unrelated gene's factor at the highest level is strictly below its factor at the "
            "lowest. A reduction is the documented direction; an increase or no change fails"
        ),
        "monotonicity": (
            "the factor never rises as the costly protein is expressed harder. One reversal fails "
            "it: a resource that is competed for does not become more available under more demand"
        ),
        "falsifier": (
            "with the pool multiplied by "
            f"{UNCONSTRAINED_MULTIPLIER:.0e} and the same demands, every factor stays 1.0. This is "
            "the clause §5.3 names explicitly, and it is the one a bug passes: a model that cuts an "
            "unrelated gene whatever the capacity is not reproducing a burden"
        ),
        "magnitude": (
            "NOT ASKABLE with what this project holds, and registered so before the run. It needs "
            "Ceroni's or Frei's measured output of an unrelated reporter at two or more expression "
            "levels, digitised with its axis stated; no such dataset is in the repository. The run "
            "reports the model's own curve so a later comparison has something to be held against, "
            "and claims nothing about agreement"
        ),
    },
    "verdict_rule": (
        "pass only if direction, monotonicity and the falsifier all hold. Any one failing fails the "
        "gate. The magnitude clause cannot contribute to a pass, because an unasked question is not "
        "evidence"
    ),
    "cannot_do": (
        "this is the model against its own declared capacities, not against a measurement. It shows "
        "that the implementation produces a burden of the documented SHAPE where a shared pool is "
        "short, and it cannot show that the shape matches any organism. The capacities themselves "
        "are cited but they are HeLa-scale round numbers, so the levels at which the curve bends are "
        "a property of the demo program rather than of a cell"
    ),
}
