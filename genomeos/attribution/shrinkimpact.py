# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Leaf-by-leaf comparison of a committed result against a recomputation on the complete stores.

Why this exists. A declared input group built as a filesystem glob filtered by `.exists()` records
whatever happened to be on the machine, so the declaration shrinks silently on a partial store
instead of refusing (`d72f959`; `lane-globcensus` counted 18 such call sites). `lane-globcensus`
stage 2 then named six committed results whose declared group count is below the chromosome
complement. A count below the complement is NOT by itself a shrink: it can equally be the run's
legitimate scope, or -- as this lane measured for `closure_cell_inventory` -- a path-prefix
GROUPING that files one chromosome's input under another group's label. Separating the three needs
the counts recomputed with the complete stores in place and compared leaf by leaf with what was
committed, which is what this module's caller does.

What it deliberately does not do. It never rewrites a committed result and it never re-runs a
registered test. A leaf that moves is reported with its old and its new value, and a correction is
written as its own additive result. A recomputed count is what the complete stores give; it does
not by itself overturn a registered reading.
"""

from __future__ import annotations

from typing import Any

#: A leaf is any value that is not a container. `None` is a leaf, and so is an EMPTY container:
#: an empty dict that becomes a dict with one key is a leaf that moved, and flattening an empty
#: container to nothing at all would hide exactly that.
_CONTAINERS = (dict, list)


def leaves(obj: Any, prefix: str = "") -> dict[str, Any]:
    """Flatten `obj` to {dotted path: leaf value}, list members addressed by index.

    The path of a top-level scalar is its key; nested keys are joined with `.`; a list member is
    `path[i]`. The separator is not escaped, so a key that itself contains a dot can in principle
    collide with a nesting -- in these results no key does, and `compare` reports the key sets it
    saw, so a collision shows up as a path present on one side only rather than as a silent match.
    """
    if isinstance(obj, dict) and obj:
        out: dict[str, Any] = {}
        for key, value in obj.items():
            out.update(leaves(value, f"{prefix}.{key}" if prefix else str(key)))
        return out
    if isinstance(obj, list) and obj:
        out = {}
        for i, value in enumerate(obj):
            out.update(leaves(value, f"{prefix}[{i}]"))
        return out
    return {prefix: obj}


def compare(committed: Any, recomputed: Any, *, label: str) -> dict[str, Any]:
    """Compare two payloads leaf by leaf and say exactly which leaves moved.

    `moved` carries the old and the new value of every leaf present on both sides with different
    values. `only_in_committed` and `only_in_recomputed` carry the paths present on one side only:
    those are not "moved" -- a leaf that one side does not have has no old or no new value -- and
    they are counted separately so that a recomputation with a different SHAPE cannot be read as
    agreement. `leaves_compared` is the denominator: the number of paths present on both sides.
    """
    left = leaves(committed)
    right = leaves(recomputed)
    shared = sorted(set(left) & set(right))
    moved = [
        {"path": path, "committed": left[path], "recomputed": right[path]}
        for path in shared
        if left[path] != right[path]
    ]
    report = {
        "label": label,
        "leaves_compared": len(shared),
        "moved": moved,
        "moved_count": len(moved),
        "only_in_committed": sorted(set(left) - set(right)),
        "only_in_recomputed": sorted(set(right) - set(left)),
    }
    report["reading"] = reading(report)
    return report


def reading(report: dict[str, Any]) -> str:
    """One sentence, and it never calls a shape difference agreement."""
    label = report["label"]
    compared = report["leaves_compared"]
    moved = report["moved_count"]
    extra = len(report["only_in_committed"]) + len(report["only_in_recomputed"])
    if moved == 0 and extra == 0:
        return (
            f"{label}: no number moves -- all {compared} recomputed leaves equal the committed "
            f"ones, and neither side carries a leaf the other does not"
        )
    if moved == 0:
        return (
            f"{label}: none of the {compared} leaves present on both sides moved, but {extra} "
            f"leaves are present on one side only, so this is NOT 'no number moves'"
        )
    return f"{label}: {moved} of {compared} leaves MOVED" + (
        f", and {extra} leaves are present on one side only" if extra else ""
    )


def every_reading_is_no_movement(reports: list[dict[str, Any]]) -> bool:
    """True only when every report is a clean no-movement: used to gate the summary wording."""
    return bool(reports) and all(
        r["moved_count"] == 0 and not r["only_in_committed"] and not r["only_in_recomputed"] for r in reports
    )
