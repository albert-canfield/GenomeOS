# SPDX-License-Identifier: AGPL-3.0-or-later
"""The cause classes for direction rule v2's disagreeing rules, fixed before any of them is read.

    uv run python scripts/contradiction_classes.py

Direction rule v2 (`genomeos/attribution/direction_v2.py`, registered `5491177`, implemented
`4e1cc41`) reads the sign in a compiled rule's own assigned cell and withholds a direction
wherever that cell's retained values do not carry one. Over the 440,377 predicted-layer rules of
24 chromosomes it left 408,127 unresolved, and the overwhelming majority of those are the
*absence* of a per-cell value rather than a disagreement in one. Two of its tallies are not
absence:

* `signs_disagree_in_cell`, an `UNRESOLVED_REASONS` entry: the rule's assigned cell retains
  values that point both ways;
* `by_cell_value_opposite_to_the_selected_extreme`, the `SELECTION_EXCLUDED_DIAGNOSTIC`:
  `by_cell[cell]` has the opposite sign to the extreme v1 selected.

This module fixes the classes that can explain such a disagreement, and the rule that picks the
primary one, **before** any of those rules is read. It does not take the classes from this lane's
own head: it imports them from `scripts/repression_trace.py`, where they were registered on
2026-09-29 (`5d14d37`, pinned by `tests/test_repression_trace.py`, `0762073`). A second vocabulary
for the same question would be worse than none.

One class is added for this population and is marked as added. Nothing in the imported five, the
imported sixth or the imported primary rule is restated, widened or softened here.

Reads nothing. Computes nothing. Classifies nothing. 0 model requests, no network, no money, no
cache read: `main()` prints the registration and says that no rule has been classified.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from genomeos.attribution import direction_v2 as dv  # noqa: E402

_spec = importlib.util.spec_from_file_location("repression_trace", ROOT / "scripts/repression_trace.py")
assert _spec is not None and _spec.loader is not None
rt = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(rt)

NAME = "contradiction_classes"

# ==================================================================================================
# The registration (2026-10-02, lane-contradict). Committed before any of the disagreeing rules was
# read, and before any classification code exists.
# ==================================================================================================
REGISTERED = "2026-10-02"
STATUS = (
    "a registration of cause classes for an internal development reading of existing compiled "
    "rules, not a validation. It says nothing about whether any direction is right"
)
QUESTION = (
    "for each predicted-layer rule where direction rule v2's retained evidence points against the "
    "direction v1 asserts, which registered cause class explains the disagreement, and does any "
    "rule survive every class as a genuine contradiction"
)

# ---- the population, taken from the committed code's own identifiers ------------------------------

#: The two tallies that are disagreement rather than absence, named by the committed identifiers
#: that produce them rather than by any prose about them. The counts are not restated here: this
#: module takes no count on trust and computes none.
POPULATION = (
    {
        "tally": "signs_disagree_in_cell",
        "where": "genomeos/attribution/direction_v2.UNRESOLVED_REASONS, decided in direction_v2._v2",
        "means": "the rule's assigned cell retains at least two numerically distinct non-zero "
        "values and predict.enhancer_target._cell_summary reports signs_disagree of them",
    },
    {
        "tally": "by_cell_value_opposite_to_the_selected_extreme",
        "where": "scripts/direction_v2.chromosome_counts, registered as "
        "direction_v2.SELECTION_EXCLUDED_DIAGNOSTIC",
        "means": "by_cell[cell] is present, is not zero, and has the opposite sign to the extreme "
        "v1 selected (predicted_coding['log2_fold_change'])",
    },
)

#: The diagnostic's registered reading, carried word for word from the registration that owns it,
#: because this lane may not strengthen it: it is a property of the cache, never a count of wrong
#: rules, and the variant that would reverse those calls was not adopted.
DIAGNOSTIC_READING = dv.SELECTION_EXCLUDED_DIAGNOSTIC

#: What v2 can and cannot do to a direction, carried from its own registration. It bounds what any
#: class assignment here can ever licence.
V2_CANNOT_REVERSE = dv.FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE

POPULATION_NOT_ASSUMED_DISJOINT = (
    "the two tallies are counted separately by scripts/direction_v2.py and their sum is NOT "
    "registered here as a number of distinct rules. Whether a rule can be counted under both is a "
    "question about the committed code and the cache, settled by reading them and never by adding "
    "the two counts together. This clause is written before either is read so that a sum cannot be "
    "adopted later for convenience"
)

NO_ENUMERATION_ASSUMED = (
    "no file is assumed to name these rules. scripts/direction_v2.py prints Counter totals and "
    "records no element identifier for any unresolved call or diagnostic hit, so whether the "
    "population exists on disk as names is itself a finding and not a premise"
)

# ---- the classes ---------------------------------------------------------------------------------

#: The registered vocabulary, imported and not copied: the same object repression_trace registered.
#: `tests/test_contradiction_classes.py` asserts it is that object by `is`.
INHERITED_CLASSES = rt.CLASSES

SIGN, CELL, CONTRADICTION, INDIRECT, JUDGE, SPLIT = (
    rt.SIGN,
    rt.CELL,
    rt.CONTRADICTION,
    rt.INDIRECT,
    rt.JUDGE,
    rt.SPLIT,
)
PRESENT, ABSENT, UNKNOWN = rt.PRESENT, rt.ABSENT, rt.UNKNOWN

#: The one class added for this population, and the only entry of this module's own making.
WITHIN_CELL_SPLIT = "g"

ADDED_CLASSES = {
    WITHIN_CELL_SPLIT: {
        "name": "direction split within the assigned cell",
        "definition": "class (f) restricted to the rule's own assigned cell: the values the cache "
        "retains for that one cell name point both ways, so the direction the rule asserts is the "
        "sign of whichever of them the selection took. (f) as registered is the split across all "
        "tracks of all cells (max_drop_log2fc < 0 < max_rise_log2fc); this is the narrower case "
        "where the opposing values are both labelled with the cell the rule names",
        "decided_by": "present when direction_v2.retained_values for the rule's assigned cell "
        "yields at least two numerically distinct non-zero values whose signs disagree under "
        "predict.enhancer_target._cell_summary - which is direction_v2._v2's own "
        "`signs_disagree_in_cell` condition, imported and not restated",
        "added_after_first_reading": False,
        "added_for_this_population": True,
        "why_added": "neither (b) nor (f) covers it. (b) is a mismatch between the claim's cell and "
        "a deciding observation's cell, and here both values carry the rule's own cell. (f) is "
        "registered over all tracks and its decided_by reads max_drop and max_rise regardless of "
        "cell, so a within-cell split satisfies (f)'s test only incidentally and (f)'s definition "
        "not at all. It is registered as a narrowing of (f) and inherits (f)'s standing: recorded "
        "beside, never primary, because it says where a direction came from and not that anything "
        "disagrees with a measurement",
    },
}

CLASSES = {**INHERITED_CLASSES, **ADDED_CLASSES}

# ---- the assignment rule -------------------------------------------------------------------------

#: The imported order, which this module does not change.
INHERITED_PRIMARY_RULE = rt.PRIMARY_RULE

PRIMARY_RULE = (
    "repression_trace.PRIMARY_RULE unchanged - (a), then (e), then (c), then (b), with (d) and (f) "
    "recorded beside and never primary - and (g) added beside (f) and likewise never primary. "
    "repression_trace.primary() is called for the decision; this module adds no ordering of its "
    "own and overrides none"
)

CONTRADICTION_IS_RESIDUAL = (
    "a rule may be called a genuine contradiction, class (c), only when no other class covers it. "
    "(c) is reached by elimination and is never the first reading: a sign that differs because of a "
    "convention, a cell mismatch or a split across tracks is not a contradiction"
)

#: The clause that decides most of this population in advance, and the reason it is registered here
#: rather than discovered later.
CONTRADICTION_NEEDS_A_MEASUREMENT = (
    "class (c) as registered requires a measurement: `the model's own value for the named gene on "
    "the measured cell's track has the opposite sign to the measured effect`, and its decided_by "
    "reads `opposite in sign to the measurement`. A predicted-layer rule carries no measured "
    "effect at all - scripts/direction_v2.chromosome_counts excludes the measured layer in those "
    "words, because a measured rule carries a screen's effect and no AlphaGenome direction - so for "
    "every rule of this population (c) is UNKNOWN by construction: there is no measurement for a "
    "model value to be opposite to. A disagreement between two model values is therefore NOT class "
    "(c) under the registered definition, however large it is. Calling it one would widen a "
    "registered class, which this lane may not do, and the class to widen is not (c) but a new one "
    "a lane with a measurement in hand may register"
)

ASSIGNMENT_RULE = (
    "for each rule of POPULATION: read each class's decided_by against the committed code and the "
    "rule's own retained values, record present, absent or unknown per class, and take the primary "
    "class from repression_trace.primary(). Where the decided_by of a class cannot be evaluated "
    "without a file this lane may not open, the class is recorded UNKNOWN and the reason named; an "
    "unevaluable class is never recorded absent, because absent is a reading and unknown is not"
)

NO_CHANGE = (
    "no rule, label, direction or verdict is changed, and none is proposed for change. What to do "
    "with a genuine contradiction is a separate decision and is not taken here: a surviving rule is "
    "named and nothing further"
)

CELL_GROUPING_CAVEAT = (
    "where loci are grouped, cell2.group is an operational grouping, NOT established biological independence"
)

MAY_BE_CALLED = (
    "the cause classes, fixed in advance, for the compiled rules where direction rule v2's own "
    "retained evidence points against the direction v1 asserts"
)
MAY_NOT_BE_CALLED = (
    "a validation of the directions, the model or the compiler",
    "a measurement of accuracy, in either direction",
    "a rate: the population is selected for disagreeing and no share of anything is estimated from it",
    "evidence that any direction is wrong: v2 can withhold a direction and never reverse one "
    "(V2_CANNOT_REVERSE), so no reading here licences a reversal",
)


def registration() -> dict[str, Any]:
    """The registered constants, as a result file would state them."""
    return {
        "name": NAME,
        "registered": REGISTERED,
        "status": STATUS,
        "question": QUESTION,
        "population": list(POPULATION),
        "population_not_assumed_disjoint": POPULATION_NOT_ASSUMED_DISJOINT,
        "no_enumeration_assumed": NO_ENUMERATION_ASSUMED,
        "diagnostic_reading": DIAGNOSTIC_READING,
        "v2_cannot_reverse": V2_CANNOT_REVERSE,
        "classes": CLASSES,
        "classes_inherited_from": "scripts/repression_trace.py, registered 2026-09-29 (5d14d37)",
        "classes_added_here": list(ADDED_CLASSES),
        "inherited_primary_rule": INHERITED_PRIMARY_RULE,
        "primary_rule": PRIMARY_RULE,
        "contradiction_is_residual": CONTRADICTION_IS_RESIDUAL,
        "contradiction_needs_a_measurement": CONTRADICTION_NEEDS_A_MEASUREMENT,
        "assignment_rule": ASSIGNMENT_RULE,
        "no_change": NO_CHANGE,
        "cell_grouping_caveat": CELL_GROUPING_CAVEAT,
        "may_be_called": MAY_BE_CALLED,
        "may_not_be_called": list(MAY_NOT_BE_CALLED),
    }


def primary(classes: dict[str, str]) -> str:
    """The primary class, from repression_trace's own order. (d), (f) and (g) are never primary."""
    if classes.get(WITHIN_CELL_SPLIT) not in (None, PRESENT, ABSENT, UNKNOWN):
        raise ValueError(f"not a reading: {classes[WITHIN_CELL_SPLIT]!r}")
    return rt.primary(classes)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="the registered cause classes, before any rule is read")
    ap.add_argument("--json", action="store_true", help="print the registration as JSON")
    args = ap.parse_args(argv)
    if args.json:
        print(json.dumps(registration(), indent=1, sort_keys=True))
        return 0
    print(f"{NAME}: registered {REGISTERED}")
    print(STATUS)
    print()
    for key, c in CLASSES.items():
        mark = " [added for this population]" if c.get("added_for_this_population") else ""
        print(f"  ({key}) {c['name']}{mark}")
    print()
    print(PRIMARY_RULE)
    print()
    print(CONTRADICTION_IS_RESIDUAL)
    print(CONTRADICTION_NEEDS_A_MEASUREMENT)
    print()
    print("no rule has been classified: this commit carries the registration and nothing else.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
