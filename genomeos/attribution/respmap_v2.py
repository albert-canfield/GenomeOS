# SPDX-License-Identifier: AGPL-3.0-or-later
"""Direction rule v2, recorded per published assertion of the response map.

`scripts/direction_v2.py` accumulates two `collections.Counter`s and writes no file: it names no
element, locus, gene or cell, so no committed artefact carries the v2 class of any one assertion.
A class that cannot be read cannot be reported, and a per-assertion class may not be re-derived
from a Counter total. This module records one row per assertion - element, gene, cell, the v1
action the map published, the v2 class and, when v2 did not resolve, its reason - so the count is
read off committed rows and not off a tally.

It applies `attribution.direction_v2` unchanged. It scores nothing, re-scores nothing, changes no
rule, no threshold and no published figure, makes 0 model requests and reads only files already on
disk: the committed response-map payload and the finished deletion sweep's own per-element cache.
"""

from __future__ import annotations

from typing import Any

from genomeos.attribution import direction_v2 as dv
from genomeos.predict import enhancer_target as et

#: The result this module's run writes, and the published payload it reads.
RESULT = "respmap_direction_v2"
SOURCE_RESULT = "response_map_increment2"

QUESTION = (
    "of the response map's published `predicted` assertions, how many does direction rule v2 leave "
    "unresolved, and how many does it resolve to the opposite action from the one the map "
    "published? Counted per assertion from committed rows, never from a Counter total"
)

POPULATION = (
    "every assertion of the published increment-2 response-map payload "
    "(`data/results/response_map_increment2.json`) whose `status` is `predicted`. That payload's "
    "own `counts.by_status` gives 166 of its 527 assertions as `predicted` and 361 as `observed`, "
    "and `observed` is reserved for what an experiment measured, so the 361 are outside this "
    "question by their status and not by any choice made here. The denominator is therefore 166, "
    "fixed by the published payload before this lane read a single direction"
)

INCREMENT_3_IS_NOT_THIS_POPULATION = (
    "increment 3 is a different population and is not counted here: the committed "
    "`data/results/response_map_increment3.json` carries 574 assertions over 337 entities in 93 "
    "chains, of which 73 are `predicted`. The figures 527 assertions, 567 entities, 195 chains and "
    "166 `predicted` are increment 2's, as its own `counts` block records. Recorded because the "
    "brief this lane was given attached increment 2's counts to increment 3's name"
)

#: The v2 class of one assertion. Exhaustive and mutually exclusive: every row carries exactly one.
CLASSES = (
    "resolved_agrees_with_published",  # v2 resolved and named the action the map published
    "resolved_opposite_to_published",  # v2 resolved and named the other action: an outright flip
    "unresolved",  # v2 withheld a direction; the row also carries the reason
)

CLASS_RULE = (
    "the published action is the assertion's own `measurement.direction`, which is "
    "`quoted.predicted_coding.action` verbatim, and v1 reads that same field "
    "(`direction_v2._v1`), so no v1 call is recomputed from a value here. The v2 call is "
    "`direction_v2.direction_call(predicted_coding, V2, row=...)` on the cached gene row for the "
    "element's named target, applied unchanged. `activates` and `represses` are the same axis as "
    "`activates` and `inhibits`: the token is compared on the axis, never on its spelling"
)

FLIP_IS_REPORTED_AS_MEASURED = (
    "the flip count is reported as measured and not as expected. "
    f"{dv.FLIPS_ARE_STRUCTURALLY_IMPOSSIBLE} - that is registered in advance, so a measured 0 here "
    "is the arithmetic of the rule and carries no agreement between v1 and v2; and a measured "
    "non-zero would refute the registered reading rather than this count, which is why each row "
    "also records whether the published value is numerically among the cell's retained values"
)

ABSENCE_GLOSS_CORRECTED = (
    "`one_value_only` may NOT be glossed as `v1's own selected extreme was the one value retained`. "
    "On one escape branch the single retained value is `by_cell[cell]` and not the published value, "
    "and it can carry the opposite sign - a disagreement sitting inside the bucket that reads as "
    "absence. Genome-wide the count of that branch is unknown and bounded above by the "
    "selection-excluded diagnostic. Over these 166 assertions it is not unknown: every row records "
    "the retained values with the field each came from, so the branch is counted exactly, as "
    "`unresolved_one_value_only_by_cell_opposite_to_published`"
)

EARLIER_TALLIES_NOT_TOUCHED = (
    "the two earlier genome-wide tallies - 14 `signs_disagree_in_cell` and 10 "
    "selection-excluded diagnostic hits - overlap and may not be added. Neither is read, reused, "
    "adjusted or re-derived here, and neither is a denominator or a numerator of anything in this "
    "result. The counts in this result are taken over this result's own 166 rows and nowhere else, "
    "so they may not be added to those two either"
)

DIAGNOSTIC_IS_NOT_A_V2_OUTPUT = dv.SELECTION_EXCLUDED_DIAGNOSTIC

VALIDATES_NOTHING = (
    f"{dv.COUNTS_VALIDATE_NOTHING}. {dv.NOT_VALIDATION}. This result says what the model's own rule "
    "says about the model's own published assertions. It does not say whether either is right, it "
    "measures no accuracy, it holds no v2 call against any measured outcome, and nothing in it is "
    "evidence that an unresolved element has no action"
)

UNRESOLVED_MEANS = dv.UNRESOLVED_MEANS

NO_RECOMMENDATION = (
    "this result is an input to the deferred decision on whether BioLang gains an action token for "
    "an unresolved direction, and it makes no recommendation about that decision. The decision "
    "changes every compiled program and is not this lane's to argue; only the number it would need "
    "is produced here"
)

NOTHING_RESCORED = (
    "no element is re-scored, no rule is edited, no threshold is moved and no published figure is "
    "changed. `MAGNITUDE_FLOOR` is still the imported `predict.enhancer_target.MIN_EFFECT` and "
    "consistency is still unanimity; this module defines no threshold of its own. The response-map "
    "payload is read and never written"
)

CANNOT_ESTABLISH = (
    "a v2 class is not a measurement of the element's action, and an unresolved class is not "
    "absence of regulation",
    "no accuracy, agreement or improvement is measured: v2 reads the same AlphaGenome sweep v1 "
    "reads, restricted to one cell name",
    "no per-locus independence is claimed: these 166 assertions are not grouped here, and "
    "`attribution.cell2.group` would be an operational grouping and NOT established biological "
    "independence if they were",
    "the cell on a row is the cell the compiled rule assigned, not a cell any experiment in the "
    "response map read: every one of these assertions carries `context.cell` as `not recorded`",
    "the counts do not generalise beyond this population: they are increment 2's published "
    "predicted assertions and not the genome-wide predicted layer",
)


def predicted_assertions(payload: dict[str, Any]) -> list[dict[str, Any]]:
    """The payload's own `predicted` assertions, in the order it published them."""
    return [a for a in payload.get("assertions") or [] if a.get("status") == "predicted"]


def locus_of(assertion: dict[str, Any]) -> tuple[str, str]:
    """The chromosome and the compiled element id, read off the assertion's own id.

    The id is built by `response_map2._compiled_rule` as `r2|<chrom>:<start>-<end>|<element>|<gene>`
    and the element field is checked against `quoted.id`, so a reshaped id fails loudly here rather
    than silently classifying the wrong element.
    """
    parts = assertion["id"].split("|")
    if len(parts) != 4:
        raise ValueError(f"not a response-map assertion id: {assertion['id']!r}")
    chrom = parts[1].split(":")[0]
    element = parts[2]
    quoted = (assertion.get("quoted") or {}).get("id")
    if quoted is not None and quoted != element:
        raise ValueError(f"id names element {element!r} but quotes {quoted!r}")
    return chrom, element


def _axis(token: str | None) -> str:
    """The activity axis of an action token, so `represses` and `inhibits` are one axis."""
    return "activates_target" if token == dv.ACTIVATES else "represses_target"


def classify(assertion: dict[str, Any], row: dict[str, Any] | None = None) -> dict[str, Any]:
    """One assertion's row: what it published, what v2 says, and why v2 said it.

    `row` is the cached gene row; when it is not given it is read from the deletion cache on disk.
    """
    chrom, element = locus_of(assertion)
    pc = (assertion.get("quoted") or {}).get("predicted_coding") or {}
    gene = pc.get("gene")
    cell = pc.get("tissue") or ""
    published = assertion.get("measurement", {}).get("direction")
    v1 = dv.direction_call(pc, dv.V1)
    if row is None:
        row = dv.cached_row(chrom, element, gene)
    v2 = dv.direction_call(pc, dv.V2, row=row)
    if not v2.resolved:
        klass = "unresolved"
    elif _axis(v2.action) == _axis(published):
        klass = "resolved_agrees_with_published"
    else:
        klass = "resolved_opposite_to_published"
    retained = list(dv.retained_values(row, cell)) if row is not None else []
    selected = pc.get("log2_fold_change")
    values = [v for _, v in retained]
    by_cell = [v for src, v in retained if src == "by_cell"]
    opposite = bool(
        by_cell
        and selected is not None
        and by_cell[0] != 0.0
        and float(selected) != 0.0
        and (by_cell[0] > 0) != (float(selected) > 0)
    )
    return {
        "assertion_id": assertion["id"],
        "chrom": chrom,
        "element": element,
        "gene": gene,
        "cell": cell,
        "published_action": published,
        "published_axis": _axis(published),
        "v1_action": v1.action,
        "published_value": selected,
        "v2_class": klass,
        "v2_action": v2.action if v2.resolved else None,
        "v2_unresolved_reason": v2.reason,
        "cached_row_found": row is not None,
        "retained_values": [{"field": src, "value": v} for src, v in retained],
        "published_value_is_among_retained": selected is not None and float(selected) in set(values),
        # ABSENCE_GLOSS_CORRECTED: the escape branch, counted and not glossed.
        "only_retained_is_by_cell_opposite_to_published": bool(
            len(retained) == 1 and retained[0][0] == "by_cell" and opposite
        ),
        # DIAGNOSTIC_IS_NOT_A_V2_OUTPUT: a property of the cache, reported beside the class.
        "by_cell_value_opposite_to_published": opposite,
        "cell_is_a_cached_track": cell in et.CELLS,
    }


def counts(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The counts this result reports, each one over these rows and never pooled with another."""
    by_class = {k: sum(1 for r in rows if r["v2_class"] == k) for k in CLASSES}
    by_reason = {
        reason: sum(1 for r in rows if r["v2_unresolved_reason"] == reason)
        for reason in dv.UNRESOLVED_REASONS
    }
    by_axis: dict[str, dict[str, int]] = {}
    for axis in ("activates_target", "represses_target"):
        on_axis = [r for r in rows if r["published_axis"] == axis]
        by_axis[axis] = {"published": len(on_axis)} | {
            k: sum(1 for r in on_axis if r["v2_class"] == k) for k in CLASSES
        }
    return {
        "denominator": len(rows),
        "by_class": by_class,
        "unresolved_by_reason": {k: v for k, v in by_reason.items() if v},
        "unresolved_reasons_with_no_assertion": sorted(k for k, v in by_reason.items() if not v),
        "by_published_axis": by_axis,
        "no_cached_row": sum(1 for r in rows if not r["cached_row_found"]),
        "published_value_among_retained": sum(1 for r in rows if r["published_value_is_among_retained"]),
        "unresolved_one_value_only_by_cell_opposite_to_published": sum(
            1 for r in rows if r["only_retained_is_by_cell_opposite_to_published"]
        ),
        "selection_excluded_diagnostic_hits": sum(
            1 for r in rows if r["by_cell_value_opposite_to_published"]
        ),
        "cells_that_are_cached_tracks": sum(1 for r in rows if r["cell_is_a_cached_track"]),
    }


def emit(payload: dict[str, Any]) -> dict[str, Any]:
    """Every published predicted assertion of `payload`, with its v2 class, and the counts."""
    rows = [classify(a) for a in predicted_assertions(payload)]
    out = counts(rows)
    total = sum(out["by_class"].values())
    if total != out["denominator"]:
        raise ValueError(f"{total} classified rows against {out['denominator']} assertions")
    return {"rows": rows, "counts": out}
