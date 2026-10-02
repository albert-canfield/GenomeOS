# SPDX-License-Identifier: AGPL-3.0-or-later
"""The reading of direction v2's disagreeing rules, under the classes registered first at c6183c2.

    uv run python scripts/contradiction_reading.py

`scripts/contradiction_classes.py` fixed the cause classes and the assignment rule before any of
these rules was read, and is left exactly as it was committed: this file carries the reading and
nothing of the registration.

The reading is a derivation from committed code rather than from data, because the population is
not enumerated on disk and the cache that would enumerate it may not be opened here. Every
statement below is decided by reading `genomeos/attribution/direction_v2.py`,
`scripts/direction_v2.py` and `genomeos/predict/enhancer_target.py`, and each branch of the
derivation is pinned by a planted row in `tests/test_contradiction_reading.py`.

Reads nothing: every value here is a constant. 0 model requests, no network, no money, no cache
read. No rule, label, direction or verdict is changed or proposed for change.
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


def _load(name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, ROOT / f"scripts/{name}.py")
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


#: the registration, imported and not restated. This file adds no class and changes no ordering.
cc = _load("contradiction_classes")

NAME = "contradiction_reading"
READ_ON = "2026-10-02"
REGISTRATION_COMMITTED_AT = "c6183c2"
READ_AFTER = (
    "the registration was committed first, at c6183c2, with no rule read, no class assigned and "
    "no classification code written"
)

# ---- what could not be read, and why --------------------------------------------------------------

ENUMERATION_ABSENT = (
    "no committed artefact names a single rule of this population, so the names cannot be "
    "reported. scripts/direction_v2.py accumulates `Counter` totals keyed by axis and by reason "
    "and writes no file at all: it has no json.dump, no output path, and it records no element "
    "identifier, locus, gene or cell for any unresolved call or diagnostic hit. No file under "
    "data/results carries any of the eight tokens of direction_v2.UNRESOLVED_REASONS or the "
    "diagnostic's key. The counts therefore exist only as the stdout of one run and as prose in "
    "the ROADMAP row, and recovering the names would require compile._attributed over 24 "
    "chromosomes together with enhancer_target.load_cached per element - the cache read this lane "
    "is forbidden. A count without its names is what was available, and it is reported as that"
)

# ---- the derivation -------------------------------------------------------------------------------

SUBSET_DERIVATION = (
    "the diagnostic's condition implies v2's `signs_disagree_in_cell` in the ordinary case, so the "
    "two tallies overlap and their sum is not a count of distinct rules. The diagnostic fires when "
    "by_cell[cell] is present, is not zero and has the opposite sign to "
    "predicted_coding['log2_fold_change']. That value is predict_target's own winner, and "
    "predict_target takes it with its tissue from one cached gene row: `log2_fold_change` is "
    "max_drop_log2fc with tissue max_drop_tissue, or max_rise_log2fc with tissue max_rise_tissue. "
    "direction_v2._v1 assigns that same tissue as the cell, and retained_values reads the same "
    "row. So whenever the cached row's winning extreme still carries the assigned cell, the cell "
    "retains that extreme AND by_cell[cell] of the opposite sign: two values, both non-zero, "
    "necessarily numerically distinct, and signs_disagree under enhancer_target._cell_summary. "
    "_v2's checks in their committed order then reach `signs_disagree_in_cell`. Two branches "
    "escape it, and both are named rather than assumed away: a third retained value of exactly "
    "0.0 reaches `zero_value_in_cell` first, which the ROADMAP row reports as 1 genome-wide; and a "
    "cached row whose max_drop_tissue and max_rise_tissue both differ from the tissue the "
    "prediction result recorded leaves by_cell[cell] as the only retained value, which reaches "
    "`one_value_only`. Which branch each diagnostic hit took cannot be decided without the cache"
)

NEITHER_TALLY_CONTAINS_THE_OTHER = (
    "`signs_disagree_in_cell` also occurs with no diagnostic hit at all, when the cell retains "
    "both max_drop and max_rise and the row keeps no by_cell value for it. So the overlap is not "
    "an identity either, and the only safe statement about the two counts is that they are "
    "counted separately and may not be added"
)

ONE_VALUE_ONLY_GLOSS_IS_TOO_STRONG = (
    "a correction to the gloss carried in the ROADMAP row, not to any code: it reads "
    "`one_value_only - the rule's assigned cell is not one of the four lines that kept a per-cell "
    "value, so the only value for that cell IS v1's own selected extreme`. The second clause does "
    "not follow from the first. On the second escape branch above the single retained value is "
    "by_cell[cell] and not the selected extreme, and it is then a value of the opposite sign - a "
    "disagreement on disk sitting inside the bucket the row reads as absence. How many rules are "
    "on that branch is unknown here: it is bounded above by the diagnostic's own count and cannot "
    "be read off any committed artefact"
)

# ---- the class assignment -------------------------------------------------------------------------

#: The per-class reading over the whole population, keyed by the registered class letters. Not one
#: of the four classes the primary rule can rank is decidable here, and the reason is the same for
#: each: the imported vocabulary was registered for a claim disagreeing with a MEASUREMENT, and
#: this population is a model disagreeing with itself.
POPULATION_CLASSES = {
    cc.SIGN: (
        cc.UNKNOWN,
        "(a)'s decided_by recomputes the label side from the cached answer with predict_target and "
        "the compiler's verb rule, and the measurement side from the raw benchmark row. Both need "
        "files this lane may not open, and the measurement side has no row to read at all. "
        "Unevaluable, so unknown and not absent",
    ),
    cc.CELL: (
        cc.UNKNOWN,
        "(b) is present when no deciding observation's cell equals the claim's cell. A "
        "predicted-layer rule has no deciding observation, so the class has nothing to compare "
        "and is unevaluable. The artefact-to-artefact cell mismatch in SUBSET_DERIVATION's second "
        "escape branch is NOT (b): (b) is registered over a claim and an observation, and this "
        "lane may not widen it to cover a prediction result disagreeing with a cache row",
    ),
    cc.CONTRADICTION: (
        cc.UNKNOWN,
        "(c) by the registration's CONTRADICTION_NEEDS_A_MEASUREMENT: unknown by construction for "
        "every rule of this population, because the class is defined against a measured effect and "
        "the predicted layer carries none",
    ),
    cc.INDIRECT: (
        cc.UNKNOWN,
        "(d) needs the model's any-gene prediction, the run's compact verdict or another gene's "
        "significant change in the same screen. None is readable here. Never primary in any case",
    ),
    cc.JUDGE: (
        cc.UNKNOWN,
        "(e) needs a verdict of correctness.verdict_of to re-run and compare. These rules carry no "
        "verdict unless S4 happened to judge them, which is not readable here",
    ),
    cc.SPLIT: (
        cc.UNKNOWN,
        "(f) is present when the cached answer records both a drop and a rise for the named gene "
        "regardless of cell. Readable only from the cache. Never primary",
    ),
    cc.WITHIN_CELL_SPLIT: (
        cc.PRESENT,
        "(g) is present for every rule of the `signs_disagree_in_cell` tally by definition: its "
        "decided_by IS _v2's condition for that reason, so membership in the tally and presence of "
        "the class are the same fact. For a diagnostic hit it is present on the ordinary branch of "
        "SUBSET_DERIVATION and unknown on the two escape branches",
    ),
}

PRIMARY_CLASS_UNREACHABLE = (
    "the imported primary rule names no primary class for this population, and that is the "
    "reading rather than a gap in it. The rule ranks (a), (e), (c) and (b); all four are unknown "
    "above, and repression_trace.primary() raises on a reading with none of them present. The one "
    "class that is present, (g), is registered as never primary - it says where a direction came "
    "from, not that anything disagrees with a measurement. The cause is structural: four of the "
    "six inherited classes presuppose a deciding observation or a judge's verdict, and this "
    "population has neither. It is a model disagreeing with itself inside one cell, which is the "
    "case the inherited vocabulary was not registered for"
)

GENUINE_CONTRADICTIONS = (
    "none, and the reason must be read exactly. Not one rule of this population survives "
    "elimination as a class (c) genuine contradiction, because (c) was never reachable: its "
    "registered definition requires a measured effect for a model value to be opposite to, and a "
    "predicted-layer rule has no measurement. This is NOT a finding that the disagreements are "
    "explained away, and it is NOT a finding that no contradiction exists in the compiled genome. "
    "It is that the registered class cannot be satisfied by evidence of this kind, so the question "
    "the brief asks cannot be answered in the vocabulary it asked for. What IS established is that "
    "every rule of the `signs_disagree_in_cell` tally is class (g): the cache retains values of "
    "both signs for the rule's own cell, and the direction the compiled program asserts is the "
    "sign of whichever one the selection took"
)

NO_WIDENING = (
    "two widenings were available and both are refused. (c) was not widened from a model value "
    "against a measurement to a model value against another model value, and (b) was not widened "
    "from a claim-against-observation cell mismatch to a prediction-result-against-cache-row one. "
    "Either would have produced an answer in the registered vocabulary by changing what the "
    "vocabulary says, which is the one thing a lane classifying against a registration may not do"
)

WHAT_WOULD_SETTLE_IT = (
    "two things, neither of them this lane's to do. First, an enumeration: scripts/direction_v2.py "
    "recording the element id, locus, gene and assigned cell of each unresolved call and each "
    "diagnostic hit into a result file, which makes the population readable without re-deriving it "
    "and settles the overlap by listing it. Second, for class (c) to be reachable at all, a "
    "measured effect on the same element, gene and cell - the measured layer, which this "
    "population excludes by construction. Neither is a change to a rule, and no rule is edited, "
    "proposed for edit or named for reversal here"
)

MAY_NOT_BE_CALLED = (
    "a validation of the directions, the model or the compiler",
    "a measurement of accuracy, in either direction",
    "a clean bill of health: no genuine contradiction was found because the class could not be "
    "reached, not because candidates were eliminated by evidence",
    "a count of rules: no rule is named and no count is taken here",
)


def reading() -> dict[str, Any]:
    """The derivation, as a result file would state it. Reads nothing: every value is a constant."""
    return {
        "name": NAME,
        "read_on": READ_ON,
        "read_after": READ_AFTER,
        "registration": cc.NAME,
        "registration_committed_at": REGISTRATION_COMMITTED_AT,
        "enumeration_absent": ENUMERATION_ABSENT,
        "subset_derivation": SUBSET_DERIVATION,
        "neither_tally_contains_the_other": NEITHER_TALLY_CONTAINS_THE_OTHER,
        "one_value_only_gloss_is_too_strong": ONE_VALUE_ONLY_GLOSS_IS_TOO_STRONG,
        "population_classes": {k: {"reading": r, "why": w} for k, (r, w) in POPULATION_CLASSES.items()},
        "primary_class_unreachable": PRIMARY_CLASS_UNREACHABLE,
        "genuine_contradictions": GENUINE_CONTRADICTIONS,
        "no_widening": NO_WIDENING,
        "what_would_settle_it": WHAT_WOULD_SETTLE_IT,
        "rules_named": [],
        "rules_named_why_empty": ENUMERATION_ABSENT,
        "may_not_be_called": list(MAY_NOT_BE_CALLED),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="the reading of direction v2's disagreeing rules")
    ap.add_argument("--json", action="store_true", help="print the registration and the reading as JSON")
    args = ap.parse_args(argv)
    if args.json:
        print(json.dumps({"registration": cc.registration(), "reading": reading()}, indent=1, sort_keys=True))
        return 0
    print(f"{NAME}: read {READ_ON}, under the classes registered at {REGISTRATION_COMMITTED_AT}")
    print()
    for key, c in cc.CLASSES.items():
        print(f"  ({key}) {c['name']}: {POPULATION_CLASSES[key][0]}")
    print()
    print(ENUMERATION_ABSENT)
    print()
    print(SUBSET_DERIVATION)
    print()
    print(NEITHER_TALLY_CONTAINS_THE_OTHER)
    print()
    print(ONE_VALUE_ONLY_GLOSS_IS_TOO_STRONG)
    print()
    print(PRIMARY_CLASS_UNREACHABLE)
    print()
    print(GENUINE_CONTRADICTIONS)
    print()
    print(NO_WIDENING)
    print()
    print("no rule is named, because no committed artefact names one; no rule is changed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
