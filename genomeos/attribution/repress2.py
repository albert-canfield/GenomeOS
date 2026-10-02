# SPDX-License-Identifier: AGPL-3.0-or-later
"""The compiler's repression calls: a blinded eligibility count, then one registered mechanism.

Two independent routes point at the same place. `data/results/not_open_profile.json` reports that of
the assessable rules whose element axis is `represses_target`, **25,943 of 30,480 (0.8511)** put their
element in a cell where reader v1 does not detect it open, against **29,141 of 51,103 (0.5702)** for
`activates_target`. Separately, the S4 trace found **all four** judged repressions wrong in direction -
a population of four, so a count and never a rate. Nothing has yet held a repression call against a
measured population.

This module does three things and no more.

**Gate 1, before any direction is read.** How many measured pairs with a *significant signed effect*
are attached to a `represses_target` rule, and to an `activates_target` rule, and over how many
independent loci. The eligibility predicate reads whether a pair's outcome is one of the two
significant labels and never which of the two, so the count is taken without reading a sign. Both
floors are imported, not chosen: `fresh.LOCUS_FLOOR` (which is `cell2.POOLED_LOCUS_FLOOR`) and
`fresh.POSITIVE_FLOOR`. Below either floor this records a no-go with the count and no direction is
read. The project's record says only 212 measured rule links exist in all, so a population too small
to decide is an expected and legitimate outcome here, not a failure of the lane.

**One mechanism, registered as a hypothesis before anything was computed.** A compiled direction is
the sign of the single most extreme of 371 tracks. Where the true effect is small that sign is close
to a coin flip, so a rise in some unrelated tissue becomes `represses`. `MECHANISM` states it,
`PREDICTION` states what it predicts, and `REFUTATION` states, in advance, the finding that would
refute it. The prediction is a conjunction of four strict comparisons and any one of them failing to
rise refutes it; there is no wording in which a flat or falling share reads as support.

**A 0-request internal check.** Each attributed element already carries
`predicted_coding_by_cell`, the same deletion sweep's signed value on each of the four retained cell
lines (`crispri.MODEL_CELLS`). For repress and activate rules separately, how often does the compiled
call's sign agree with the element's own sign in those four cells. This is a consistency check between
two readings of the same model, so a call that disagrees in every retained cell is an outlier of this
reading and not by itself a mechanism.

**What this module does not do.** It changes no rule, no compiled label and no verdict, it deletes
nothing, it refits nothing and it makes no model request. The openness axis remains a consistency
check between two readings of the same ENCODE chromatin: `not_open_in_reader` means *not detected open
at the reader's registered call*, never *closed* (`context_evidence.NOT_CLOSED`), and reader v1's
peaks and AlphaGenome's training set are the same chromatin, so nothing here is independent evidence
and nothing here is validation (`context_evidence.NOT_VALIDATION`). Every band, floor and grouping is
imported from where the project already fixed it; this module introduces no threshold of its own.
"""

from __future__ import annotations

from typing import Any

from genomeos.attribution import cell2, fresh
from genomeos.attribution import context_evidence as ce
from genomeos.attribution import measured as ms
from genomeos.attribution import not_open_profile as nop
from genomeos.attribution.crispri import MODEL_CELLS
from genomeos.attribution.crispri_direction import wilson

# ---- the two activity axis values this lane separates -------------------------------------------

#: The compiler's own axis values, as `compile.element_axes` writes them. A rule is counted under one
#: of these only when its element's activity axis is *exactly* that value: the axis can hold more than
#: one term at once (`compile.measured_axes` joins them with a comma) and such a rule is neither a
#: repression rule nor an activation rule for this lane's purposes. It is counted in its own row.
REPRESSES = "represses_target"
ACTIVATES = "activates_target"
ACTIVITIES = (REPRESSES, ACTIVATES)
ACTIVITY_CALL = (
    "the element's own `activity:` axis, read off the compiled element by "
    "not_open_profile.axis_value, which reads compile.element_axes and compile.measured_axes "
    "unchanged. A rule counts as a repression rule only when that axis is exactly "
    f"{REPRESSES!r} and as an activation rule only when it is exactly {ACTIVATES!r}; a rule whose "
    "axis holds more than one term is counted in its own row and is in neither population"
)

#: Why the measured layer's axis is not that layer's rule direction, stated here because a count that
#: mixed the two would be wrong and the asymmetry is easy to miss. `measured.rule_links` emits a link
#: only from a pair the benchmark calls `Regulated`, which is `Significant AND EffectSize < 0`, so
#: every CRISPRi link's action is `activates`; a significant *increase* raises no rule at all and
#: reaches the program only as `represses_target` on the element axis (`compile.measured_axes`).
MEASURED_AXIS_IS_NOT_THE_LINK_DIRECTION = (
    "on a measured element the activity axis and the rule direction are different things: "
    "measured.rule_links emits a link only from a pair the benchmark calls Regulated, which is "
    "Significant AND EffectSize < 0, so every compiled CRISPRi link's action is `activates`, while a "
    "significant increase raises no rule and reaches the program only as `represses_target` on the "
    "element's activity axis. Counts are therefore reported per rule source as well as per axis, and "
    "no count here reads a measured element's axis as the direction of a measured link"
)

# ---- gate 1: the blinded eligibility count ------------------------------------------------------

#: The two outcome labels `measured.CrispriPair.outcome` gives a significant pair. Membership in this
#: pair of labels is the whole eligibility predicate: the gate asks whether a pair's outcome is one of
#: the two and never which of the two, so no direction and no sign enters the count.
SIGNIFICANT_OUTCOMES = (ms.DECREASE, ms.INCREASE)

#: Everything the gate is permitted to read about a measured pair. The blinding condition of the gate:
#: an effect size, a p-value, a power column, a model score or a predicted value may not enter it.
GATE_FIELDS = ("gene", "cell", "chrom", "start", "end", "outcome_is_one_of_the_two_significant_labels")

GATE_CALL = (
    "a measured pair is eligible when measured.CrispriPair.outcome is one of the two labels a "
    f"significant pair takes, {SIGNIFICANT_OUTCOMES}. The gate reads membership in that pair of "
    "labels and never which member, so it is taken without reading a direction or a sign. This is "
    "the same two-sided significance crispri_direction.signed uses, read off the cached row's own "
    "outcome field rather than recomputed"
)

#: Both floors are imported. `fresh.LOCUS_FLOOR` is `cell2.POOLED_LOCUS_FLOOR`, so a locus count here
#: means exactly what it means in the second-cell-type lane and the two can be compared.
LOCUS_FLOOR = fresh.LOCUS_FLOOR
POSITIVE_FLOOR = fresh.POSITIVE_FLOOR
FLOORS_CALL = (
    f"independent loci at least {LOCUS_FLOOR} (fresh.LOCUS_FLOOR, which is "
    f"cell2.POOLED_LOCUS_FLOOR) and eligible measured links at least {POSITIVE_FLOOR} "
    "(fresh.POSITIVE_FLOOR). Both are imported from where the project fixed them before this lane "
    "existed; neither is set here and neither may be moved after a count"
)

#: The locus convention, imported from `cell2` rather than restated. It is an operational grouping for
#: deciding whether a population has enough distinct places in the genome to measure, and it is **not**
#: established biological independence; nothing here may be cited as an independence guarantee.
LOCUS_RULE = cell2.INDEPENDENT_LOCUS_RULE
LOCUS_SPAN = cell2.INDEPENDENT_LOCUS_SPAN
NOT_BIOLOGICAL_INDEPENDENCE = (
    "the grouping is cell2's, operational and not established biological independence: two links it "
    "calls one locus may be two regulatory events and two it calls separate may share a domain, a TAD "
    "or a trans factor"
)

#: The two readings of gate 1, fixed before the count. There is no third.
GATE_PASS = "at or above both registered floors: this population can carry a direction reading"
GATE_NO_GO = (
    "below a registered floor: a no-go. The count goes on the record, no direction is read on this "
    "population, no rate is reported and the floor is not moved"
)
GATE_READINGS = (GATE_PASS, GATE_NO_GO)

# ---- the mechanism, registered as a hypothesis before anything was computed ---------------------

MECHANISM = (
    "a compiled direction is the sign of the single most extreme of 371 AlphaGenome tracks (the S4 "
    "trace, 641909e): the element's `predicted_coding` block keeps the largest |log2 fold change| "
    "over the tracks and the compiler reads `activates` or `represses` off that one value's sign. "
    "Where the element's true effect is small, which track is the most extreme is close to arbitrary "
    "and its sign is close to a coin flip, so a rise in some unrelated tissue is read as "
    "`represses_target`. On this account a large part of the repression population is not a repression "
    "mechanism at all but the sign of a small number chosen by a maximum"
)

PREDICTION = (
    "if that is what is happening, repression calls concentrate where the effect is weak and where "
    "the reader does not detect the element open. Over the assessable rules of the predicted layer, "
    "crossed by effect band (not_open_profile's own `strength`: strong at or above "
    "enhancer_target.STRONG_EFFECT, weak below) and by reader state (context_evidence.state_for: "
    "open_in_reader against not_open_in_reader), the repression share rises along both axes: "
    "weak above strong at each reader state, and not_open above open in each band. That is four "
    "strict comparisons, and the prediction is their conjunction"
)

REFUTATION = (
    "a repression share that does not rise as the effect weakens and as openness falls refutes the "
    "prediction. Concretely: the prediction is refuted if any one of its four comparisons - "
    "weak vs strong within open_in_reader, weak vs strong within not_open_in_reader, not_open vs open "
    "within the strong band, not_open vs open within the weak band - is not a strict rise in the "
    "repression share. A flat or falling share is a refutation and is reported as one; there is no "
    "wording in which it reads as support. Each share is reported with its Wilson 95% interval and "
    "the population it was counted over, so a rise inside the noise can be seen for what it is, but "
    "the reading is decided on the four directions as registered above and on nothing else"
)

#: The two readings of the prediction, fixed before the table. There is no third and no "partly".
PREDICTION_SUPPORTED = (
    "the repression share rises along both axes, as the mechanism predicts: all four registered "
    "comparisons are strict rises. This is consistent with the mechanism and establishes it about no "
    "rule; a descriptive concentration is not a demonstration that the sign came from a maximum"
)
PREDICTION_REFUTED = (
    "the repression share does not rise along at least one axis, so the prediction as registered does "
    "not hold and the mechanism is refuted on this population. A refuted mechanism is a result"
)
PREDICTION_READINGS = (PREDICTION_SUPPORTED, PREDICTION_REFUTED)

#: The four comparisons, named so the result can report each one by name rather than by position.
COMPARISONS = (
    ("weak_above_strong_within_open_in_reader", "weak", "strong", ce.STATE_OPEN, ce.STATE_OPEN),
    ("weak_above_strong_within_not_open_in_reader", "weak", "strong", ce.STATE_NOT_OPEN, ce.STATE_NOT_OPEN),
    ("not_open_above_open_within_strong", "strong", "strong", ce.STATE_NOT_OPEN, ce.STATE_OPEN),
    ("not_open_above_open_within_weak", "weak", "weak", ce.STATE_NOT_OPEN, ce.STATE_OPEN),
)

# ---- the 0-request internal check ---------------------------------------------------------------

#: The four cell lines the deletion sweep was scored in, imported from `crispri` rather than restated.
RETAINED_CELLS = MODEL_CELLS
#: The field on the attributed element that carries their signed values. The same precedence
#: `crispri.deletion_values` uses, coding first; this check reads the coding field only, because the
#: compiled rule's own direction is read off `predicted_coding`.
BY_CELL_FIELD = "predicted_coding_by_cell"

CONCORDANCE_CALL = (
    "the compiled call's sign is the sign of the element's predicted_coding log2 fold change, which "
    "is what the compiler reads the direction off. Beside it, the element's own "
    f"{BY_CELL_FIELD} holds the same sweep's signed value on each of {list(RETAINED_CELLS)}. A cell "
    "agrees when its value has the same sign; a value of exactly zero has no sign and agrees with "
    "neither, and is counted under its own name. Only rules of the predicted layer are counted: a "
    "measured rule carries the screen's effect size and no AlphaGenome effect of its own"
)

#: What a disagreement in every retained cell is, and what it is not.
OUTLIER_CALL = (
    "a call that disagrees with the element's own sign in every one of the four retained cells is an "
    "outlier of this reading: the one track the maximum selected points one way and all four retained "
    "lines point the other. It is counted as an outlier and not as a mechanism"
)
CONCORDANCE_IS_NOT_VALIDATION = (
    "this is a consistency check between two readings of the same model on the same element - the "
    "most extreme track against four named cell lines of the same sweep - so it is never independent "
    "evidence and never validation. The compiled call's own tissue is in general not one of the four, "
    "so a disagreement may be the cell-specific answer being different rather than the call being "
    "wrong, and no rule is marked wrong by this count"
)

# ---- gate 1 -------------------------------------------------------------------------------------


def significant(pair: dict[str, Any]) -> bool:
    """Whether one cached CRISPRi pair carries a significant signed effect, read without its sign.

    The only predicate the gate applies. It asks whether the pair's outcome is one of the two labels a
    significant pair takes and never which of the two, so no direction enters the count.
    """
    return pair.get("outcome") in SIGNIFICANT_OUTCOMES


def eligible_links(
    rules: list[nop.Rule], measured_rows: list[dict[str, Any]], activity: str, match_cell: bool
) -> list[cell2.LocusKey]:
    """The blinded locus keys of every (rule, measured pair) link this activity axis carries.

    One key per (rule, eligible pair) whose element, gene and - when `match_cell` - cell agree. The
    key is `cell2.LocusKey`, which is the only path into the grouping and carries nothing but the
    cell, the chromosome, the element interval and the gene.

    `match_cell` is reported both ways and never pooled. A rule is gated by `when: cell_type`, so a
    measurement in another cell is not a reading of that rule in its own cell; a count that ignored
    the cell would be larger and would be answering a different question.
    """
    by_element: dict[str, list[dict[str, Any]]] = {}
    for row in measured_rows:
        pairs = (row["measured"].get("crispri") or {}).get("pairs") or []
        if pairs:
            by_element.setdefault(row["id"], []).extend(pairs)
    out: list[cell2.LocusKey] = []
    for r in rules:
        if r.activity_axis != activity:
            continue
        for pair in by_element.get(_element_id(r), ()):
            if pair["gene"] != r.gene:
                continue
            if match_cell and pair["cell"] != r.cell:
                continue
            if significant(pair):
                out.append(cell2.LocusKey(r.cell, r.chrom, r.start, r.end, r.gene))
    return out


def _element_id(rule: nop.Rule) -> str:
    """The attributed element a rule sits on. A measured rule's element carries the `_measured`
    suffix `not_open_profile.rules` gives it; the measured row's own id is the bare one."""
    return rule.element[: -len("_measured")] if rule.element.endswith("_measured") else rule.element


def gate(keys: list[cell2.LocusKey], population: str) -> dict[str, Any]:
    """One population against both registered floors. No direction has been read at this point."""
    links = len(keys)
    loci = cell2.count_loci(keys, LOCUS_SPAN)
    met = links >= POSITIVE_FLOOR and loci >= LOCUS_FLOOR
    short = []
    if links < POSITIVE_FLOOR:
        short.append(
            f"eligible measured links {links}, short of {POSITIVE_FLOOR} by {POSITIVE_FLOOR - links}"
        )
    if loci < LOCUS_FLOOR:
        short.append(f"independent loci {loci}, short of {LOCUS_FLOOR} by {LOCUS_FLOOR - loci}")
    return {
        "population": population,
        "eligible_measured_links": links,
        "distinct_measured_pairs_behind_them": len(set(keys)),
        "independent_loci": loci,
        "floors": {"independent_loci": LOCUS_FLOOR, "eligible_measured_links": POSITIVE_FLOOR},
        "meets_both_floors": met,
        "reading": GATE_PASS if met else GATE_NO_GO,
        "short_by": short,
        "independent_locus_rule": LOCUS_RULE,
        "not_biological_independence": NOT_BIOLOGICAL_INDEPENDENCE,
    }


# ---- the descriptive table ----------------------------------------------------------------------


def share_row(repress: int, activate: int) -> dict[str, Any]:
    """The repression share of one cell of the table, with its population and Wilson interval."""
    assessable = repress + activate
    w = wilson(repress, assessable)
    return {
        "represses_target": repress,
        "activates_target": activate,
        "assessable_rules": assessable,
        "repression_share": None if not assessable else round(repress / assessable, 4),
        "repression_share_ci95": w["ci95"],
        "population_of_that_share": (
            f"the {assessable} rules of this cell of the table whose activity axis is exactly "
            f"{REPRESSES!r} or exactly {ACTIVATES!r}"
        ),
    }


def table(counts: dict[tuple[str, str, str], int]) -> dict[str, Any]:
    """The repression share by effect band and reader state, from counts keyed (band, state, axis)."""
    out: dict[str, Any] = {}
    for band in nop.EFFECT_BANDS:
        for state in (ce.STATE_OPEN, ce.STATE_NOT_OPEN):
            out[f"{band}|{state}"] = share_row(
                counts.get((band, state, REPRESSES), 0), counts.get((band, state, ACTIVATES), 0)
            )
    return out


def comparisons(cells: dict[str, Any]) -> dict[str, Any]:
    """Each registered comparison, by name, as a strict rise or not. Decided on direction only."""
    out: dict[str, Any] = {}
    for name, hi_band, lo_band, hi_state, lo_state in COMPARISONS:
        hi = cells[f"{hi_band}|{hi_state}"]
        lo = cells[f"{lo_band}|{lo_state}"]
        a, b = hi["repression_share"], lo["repression_share"]
        out[name] = {
            "higher_side": f"{hi_band}|{hi_state}",
            "lower_side": f"{lo_band}|{lo_state}",
            "higher_share": a,
            "lower_share": b,
            "rises": None if a is None or b is None else a > b,
            "difference": None if a is None or b is None else round(a - b, 4),
        }
    return out


def prediction_reading(compared: dict[str, Any]) -> dict[str, Any]:
    """The reading of the prediction, in the wording registered before the table was computed."""
    rises = {k: v["rises"] for k, v in compared.items()}
    empty = [k for k, v in rises.items() if v is None]
    supported = not empty and all(rises.values())
    return {
        "per_comparison": rises,
        "comparisons_on_an_empty_cell": empty,
        "supported": supported,
        "reading": PREDICTION_SUPPORTED if supported else PREDICTION_REFUTED,
        "which_failed": [k for k, v in rises.items() if v is False],
        "refutation_condition_as_registered": REFUTATION,
    }


# ---- the 0-request internal check ---------------------------------------------------------------

AGREES, DISAGREES, NO_SIGN, NOT_SCORED = "agrees", "disagrees", "no_sign_value_is_zero", "not_scored"


def cell_signs(call: float, by_cell: dict[str, Any] | None) -> dict[str, str]:
    """How each retained cell's own signed value stands to the compiled call's sign."""
    out = {}
    for cell in RETAINED_CELLS:
        v = (by_cell or {}).get(cell)
        if v is None:
            out[cell] = NOT_SCORED
        elif float(v) == 0.0:
            out[cell] = NO_SIGN
        else:
            out[cell] = AGREES if (float(v) > 0) == (call > 0) else DISAGREES
    return out


def concordance(rows: list[dict[str, str]]) -> dict[str, Any]:
    """Sign concordance across the four retained cells over one activity axis's rules."""
    scored = [r for r in rows if all(v != NOT_SCORED for v in r.values())]
    agree_counts: dict[int, int] = {}
    for r in scored:
        n = sum(v == AGREES for v in r.values())
        agree_counts[n] = agree_counts.get(n, 0) + 1
    outliers = sum(1 for r in scored if all(v == DISAGREES for v in r.values()))
    unanimous = sum(1 for r in scored if all(v == AGREES for v in r.values()))
    return {
        "rules": len(rows),
        "rules_scored_in_all_four_retained_cells": len(scored),
        "rules_not_scored_in_all_four": len(rows) - len(scored),
        "cells_agreeing_with_the_compiled_sign": {str(k): agree_counts.get(k, 0) for k in range(5)},
        "agrees_in_all_four": unanimous,
        "disagrees_in_all_four": outliers,
        "disagrees_in_all_four_share": None if not scored else round(outliers / len(scored), 4),
        "population_of_that_share": (
            f"the {len(scored)} rules of this axis whose element carries a value in all four of "
            f"{list(RETAINED_CELLS)}"
        ),
        "per_cell_agreement": {cell: sum(1 for r in scored if r[cell] == AGREES) for cell in RETAINED_CELLS},
        "per_cell_value_exactly_zero": {
            cell: sum(1 for r in scored if r[cell] == NO_SIGN) for cell in RETAINED_CELLS
        },
    }


# ---- the registration ---------------------------------------------------------------------------


def limitations() -> dict[str, str]:
    """Imported from the readings this lane describes, not restated."""
    return {
        "not_closed": ce.NOT_CLOSED,
        "not_validation": ce.NOT_VALIDATION,
        "descriptive_only": (
            "the effect-band by reader-state table is descriptive. It locates where repression calls "
            "sit; it establishes about no rule that the rule is wrong, and it moves no verdict."
        ),
        "concordance_is_not_validation": CONCORDANCE_IS_NOT_VALIDATION,
        "measured_axis_is_not_the_link_direction": MEASURED_AXIS_IS_NOT_THE_LINK_DIRECTION,
        "not_biological_independence": NOT_BIOLOGICAL_INDEPENDENCE,
        "four_is_a_count": (
            "the S4 trace's four judged repressions, all wrong in direction, are a count of four and "
            "never a rate. Nothing here turns them into one."
        ),
    }


def registration() -> dict[str, Any]:
    """Every call this lane makes, from the code itself, committed before any result is computed."""
    return {
        "lane": "lane-repress2",
        "status": (
            "registered before any count, any table and any direction was read; nothing in this "
            "payload is a result"
        ),
        "question": (
            "do the compiler's repression calls hold against a measured population, and is the "
            "repression population concentrated where the most-extreme-track mechanism says it "
            "should be"
        ),
        "gate_1": {
            "what_is_counted": (
                "measured pairs with a significant signed effect attached to a rule of each activity "
                "axis, and the independent loci they fall in"
            ),
            "call": GATE_CALL,
            "blinded_fields": list(GATE_FIELDS),
            "significant_outcome_labels": list(SIGNIFICANT_OUTCOMES),
            "floors": FLOORS_CALL,
            "locus_floor": LOCUS_FLOOR,
            "link_floor": POSITIVE_FLOOR,
            "independent_locus_rule": LOCUS_RULE,
            "independent_locus_span": LOCUS_SPAN,
            "readings": list(GATE_READINGS),
            "below_the_floor": (
                "the count goes on the record as a no-go and the lane continues only as the "
                "descriptive part. The floor is not lowered and no rate is reported"
            ),
            "cell_matching": (
                "reported both ways and never pooled: a rule is gated by `when: cell_type`, so a "
                "measurement in another cell is not a reading of that rule in its own cell"
            ),
        },
        "mechanism": MECHANISM,
        "prediction": PREDICTION,
        "refutation": REFUTATION,
        "prediction_readings": list(PREDICTION_READINGS),
        "comparisons": [c[0] for c in COMPARISONS],
        "activity_call": ACTIVITY_CALL,
        "effect_band_call": nop.EFFECT_BAND_CALL,
        "strong_effect": nop.STRONG_EFFECT,
        "effect_bands": list(nop.EFFECT_BANDS),
        "state_call": (
            "genomeos.attribution.context_evidence.state_for, imported unchanged, with its own "
            "mapping table and reader v1's own openness call. No threshold, mapping or registered "
            "definition is changed here"
        ),
        "internal_check": {
            "call": CONCORDANCE_CALL,
            "retained_cells": list(RETAINED_CELLS),
            "field": BY_CELL_FIELD,
            "outlier": OUTLIER_CALL,
            "not_validation": CONCORDANCE_IS_NOT_VALIDATION,
        },
        "limitations": limitations(),
        "introduces_no_threshold": (
            "every band, floor, grouping and state in this lane is imported from where the project "
            "already fixed it: not_open_profile's effect band and enhancer_target.STRONG_EFFECT, "
            "context_evidence.state_for, cell2's locus rule and span, fresh.LOCUS_FLOOR and "
            "fresh.POSITIVE_FLOOR, crispri.MODEL_CELLS. This lane sets no value of its own"
        ),
        "changes_nothing": (
            "no rule, no compiled label and no verdict is changed, nothing is deleted and nothing is refitted"
        ),
        "alphagenome_requests": 0,
        "money": "none: no request is made",
        "code": {
            "module": "genomeos/attribution/repress2.py",
            "register": "scripts/repress_register.py",
            "run": "scripts/repress_population.py",
            "tests": "tests/test_repress2.py",
        },
    }
