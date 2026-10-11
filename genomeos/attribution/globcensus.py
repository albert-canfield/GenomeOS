# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""globcensus (2026-10-11): declared input groups built by a glob filtered on `.exists()`.

THE DEFECT, in the supervisor's framing rather than either of the two we tried first. A manifest's
declared input GROUP is sometimes built like this:

    mf.files_entry("reader v1 DNase narrowPeak files", [p for p in peak_files(CHROMS) if p.exists()])

The `.exists()` filter means the declaration is a function of the FILESYSTEM at write time rather than
of the commit. On a partial store the group does not refuse and does not warn: it SILENTLY SHRINKS, and
the result still looks complete, with nothing anywhere saying a smaller set was declared than intended.

ABSENCE IS SILENT. That is the whole of it, and it is narrower and more interesting than the two
descriptions this lane began with -- one of us called the artefacts faulty, the other called the tooling
incapable, and both were wrong. The declaration is honest about what it saw; what is missing is any
record that it saw less than the writer meant.

WHAT THIS CENSUS IS: static, a list and a count. No repairs, no rewrites, and no result's bytes touched.
It is INTERIM: the counts below are measured, and the classification of each VARYING group into
"legitimate scope difference" or "silent shrink" is NOT done. That distinction needs each result's own
declared scope read against its group, which is the next step and is not claimed here.
"""

from __future__ import annotations

import ast
import collections
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent.parent

#: The two trees that can write a manifest.
SEARCHED = ("genomeos", "scripts")


def files_entry_calls(root: Path | None = None) -> list[dict[str, Any]]:
    """Every `files_entry(...)` call, by AST, with whether its arguments mention `exists`.

    BY AST AND NOT BY GREP, and that is not fussiness: a one-line `grep` for `files_entry.*exists()`
    misses a call whose filter sits on a later line, and exactly one of the call sites is of that shape
    (`scripts/trio_q100_phase.py:72`, where the call spans four lines and the filter is on line 74).
    The brief's own count of 17 came from the one-line form; the true figure is 18. This is the same
    instrument class the project has been finding all week, here inside the measurement of the defect.
    """
    base = root if root is not None else ROOT
    out: list[dict[str, Any]] = []
    for tree in SEARCHED:
        for p in sorted((base / tree).rglob("*.py")):
            try:
                parsed = ast.parse(p.read_text())
            except SyntaxError:
                continue
            for node in ast.walk(parsed):
                if not isinstance(node, ast.Call):
                    continue
                name = getattr(node.func, "attr", None) or getattr(node.func, "id", None)
                if name != "files_entry":
                    continue
                src = ast.unparse(node)
                out.append(
                    {
                        "path": str(p.relative_to(base)),
                        "line": node.lineno,
                        "exists_filtered": "exists" in src,
                        "spans_lines": (node.end_lineno or node.lineno) - node.lineno + 1,
                    }
                )
    return out


def declared_groups(results_dir: Path | None = None) -> dict[str, list[tuple[str, Any]]]:
    """Per group name, the `(result, member count)` of every committed declaration of it."""
    base = results_dir if results_dir is not None else ROOT / "data/results"
    groups: dict[str, list[tuple[str, Any]]] = collections.defaultdict(list)
    for f in sorted(base.glob("*.json")):
        try:
            d = json.loads(f.read_text())
        except Exception:
            continue
        for item in (d.get("result_manifest") or {}).get("inputs") or []:
            if isinstance(item, dict) and item.get("group") is True:
                groups[str(item.get("path"))].append((f.stem, item.get("files")))
    return dict(groups)


def varying(groups: dict[str, list[tuple[str, Any]]]) -> dict[str, list[Any]]:
    """The group names whose member count is NOT the same in every result that declares them.

    A varying count is the census's real question and NOT by itself a defect: a result that genuinely
    covered 23 chromosomes should declare 23. What the count cannot distinguish is that case from a
    result that MEANT 24 and silently got 23 because a store was partial. Naming them is this census's
    deliverable; separating them is the next step.
    """
    return {
        k: sorted({n for _, n in v}, key=lambda x: (x is None, x))
        for k, v in groups.items()
        if len({n for _, n in v}) > 1
    }


#: What the census does not claim, written here so a later reader cannot widen it.
DOES_NOT_CLAIM = (
    "not that any committed result IS wrong; not that a varying member count is a defect (a result"
    " that covered fewer chromosomes should declare fewer); not that the constant-count groups are"
    " safe, since a group constant at a shrunken value would look identical to one constant at its"
    " full value; and nothing about groups built WITHOUT an exists() filter, which are the majority"
)
