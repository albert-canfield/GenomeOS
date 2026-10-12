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


# ---------------------------------------------------------------------------------------- stage 2
# Stage 1 counted and deliberately stopped: it named 6 groups whose member count varies between the
# results declaring them, and refused to call any of them a shrink from a count alone. Stage 2 does the
# distinguishing check stage 1 said was needed -- each group read against ITS OWN result's declared
# scope, and then against what the store holds now.
#
# AND THE COUNT ALONE WOULD HAVE BEEN WRONG, which is why stage 1 refused it: six of the seven short
# declarations are missing chrY, and a group of 22 against a scope of 24 is exactly what a deliberate
# chr1-22 analysis looks like. The verdict turns on WHICH chromosome is absent, not on how many.

#: A short declaration is only a SHRINK when the scope it falls short of is the result's own, AND the
#: missing member is present in the store now -- which together say the store was partial at write time
#: rather than the analysis being narrower on purpose.
SHRINK_TEST = (
    "a declaration is a SILENT SHRINK when (1) its member count is below the result's own"
    " parameters.chromosomes, (2) the missing chromosomes are present in the store today, and (3) no"
    " exclusion in that result names them. All three are required: (1) alone is satisfied by any"
    " deliberately narrower analysis, and (2) alone by a store that has merely grown since"
)

#: Measured 2026-10-12. Six declarations in six committed results meet all three conditions.
SHRINKS_FOUND = (
    {
        "result": "direction_link",
        "group": "per_element_response_cache",
        "scope": 24,
        "declared": 22,
        "missing": ("chr13", "chrY"),
        "written": "2026-10-02",
    },
    {
        "result": "hct116_count",
        "group": "per_element_response_cache",
        "scope": 24,
        "declared": 22,
        "missing": ("chr13", "chrY"),
        "written": "2026-10-02",
    },
    {
        "result": "reptest_answerable",
        "group": "per_element_response_cache",
        "scope": 24,
        "declared": 22,
        "missing": ("chr13", "chrY"),
        "written": "2026-10-02",
    },
    {
        "result": "rescore_cost",
        "group": "per_element_response_cache",
        "scope": 24,
        "declared": 22,
        "missing": ("chr13", "chrY"),
        "written": "2026-10-02",
    },
    {
        "result": "closure_cell2_inventory",
        "group": "GENCODE 50 per-chromosome annotation",
        "scope": 24,
        "declared": 23,
        "missing": ("chr21",),
        "written": "2026-10-10",
    },
    {
        "result": "closure_cell_inventory",
        "group": "GENCODE 50 per-chromosome annotation",
        "scope": 24,
        "declared": 22,
        "missing": ("chr21", "chrY"),
        "named_in_exclusions": ("chrY",),
        "written": "2026-10-10",
    },
)

#: The reading, and it is the one the framing predicted rather than the one the lane went looking for.
#: NARROWED by its own guard before it was committed: the first draft said "NOT ONE of them says so",
#: and the silence test refused it, because closure_cell_inventory DOES name one of its two.
STAGE_2_READING = (
    "THE DEFECT IS REALISED. Six committed results each declare parameters.chromosomes = 24 and an"
    " input group of 22 or 23, and in every case the missing chromosomes are on disk TODAY -- so the"
    " store was partial when they were written and has since been completed. Each result therefore"
    " STATES a 24-chromosome scope and was computed on fewer. FIVE OF THE SIX are wholly silent: their"
    " exclusions name guide counts, signs and thresholds and never an absent chromosome. THE SIXTH is"
    " partly explained -- closure_cell_inventory states 'chrY: the library's design carries no guide"
    " there', which is a declared exclusion and not a shrink, but it says nothing about chr21, which"
    " is. So 9 of the 10 missing chromosomes across the six are unexplained, and chrY in that one"
    " result is the only one any result accounts for. This is not one incident either: two stores"
    " (data/knowledge/alphagenome/elements and data/reference) on two dates (2026-10-02 and"
    " 2026-10-10), eight days apart"
)

#: What stage 2 still does not establish, kept separate from what it does.
STAGE_2_DOES_NOT_CLAIM = (
    "not that any of the six results' NUMBERS are wrong -- a two-chromosome shortfall may move a"
    " genome-wide figure or may not, and this census computes no figure; not that chrY's absence is"
    " the same kind of thing as chr13's or chr21's, since chrY is routinely empty; not that the other"
    " declarations are safe, because a group constant at a shrunken value remains invisible here; and"
    " nothing about results whose parameters declare no chromosome scope at all, of which there are"
    " six and which cannot be classified by this test"
)


def is_shrink(declared: int, scope: int, missing_now_on_disk: bool, named_in_exclusions: bool) -> bool:
    """SHRINK_TEST as a function: all three conditions, none sufficient alone."""
    return declared < scope and missing_now_on_disk and not named_in_exclusions
