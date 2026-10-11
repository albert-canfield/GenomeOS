#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""globcensus: the count, written once. Static, a list and a count, no repairs.

    uv run python scripts/globcensus_count.py

INTERIM BY CONSTRUCTION: it counts call sites and declared group sizes. It does NOT classify a varying
count into "legitimate scope difference" or "silent shrink" -- that needs each result's own declared
scope read against its group, and claiming it from a count alone is the error this census exists to
describe.

Input paths are handed to input_entry REPO-RELATIVE: input_entry records str(path) verbatim while
files_entry relativises, so an absolute path would record this machine's checkout. manifest.py is inside
the paid study's frozen closure and cannot be fixed, so a relative path at the call site is the only way
not to add another unrebuildable result.
"""

from __future__ import annotations

from pathlib import Path

from genomeos import manifest as mf
from genomeos.attribution import globcensus as gc
from genomeos.results import save_result

RESULT = "globcensus"
OWN_CODE = (
    "genomeos/attribution/globcensus.py",
    "scripts/globcensus_count.py",
    "tests/test_globcensus.py",
)


def manifest(sources: list[Path]) -> dict:
    inputs = [mf.input_entry(p, partition=None) for p in sources]
    return {
        "sources": [
            {"accession": "the committed result registry, data/results/*.json", "version": "as committed"},
            {"accession": "the genomeos/ and scripts/ trees, read by AST", "version": "as committed"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "searched_trees": list(gc.SEARCHED),
            "method": "ast.parse over every .py, then every files_entry call; group declarations read"
            " from result_manifest.inputs where group is true",
            "why_not_grep": "a one-line pattern cannot see a call whose filter sits on a later line",
        },
        "exclusions": [
            "static only: nothing was run, no store was read, no result's bytes were touched, and no"
            " varying count is classified"
        ],
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        "partitions": {
            "none": "one population, the committed registry; no partition is taken and the absence of"
            " one is not a claim of independence"
        },
    }


def main() -> int:
    calls = gc.files_entry_calls()
    filtered = [c for c in calls if c["exists_filtered"]]
    groups = gc.declared_groups()
    var = gc.varying(groups)
    declarations = sum(len(v) for v in groups.values())
    reader = groups.get("reader v1 DNase narrowPeak files", [])

    payload = {
        "question": (
            "which committed results declare an input group built by a glob filtered on .exists(), so"
            " the declaration silently SHRINKS on a partial store instead of refusing"
        ),
        "the_defect": (
            "ABSENCE IS SILENT. files_entry(name, [p for p in paths if p.exists()]) makes the declared"
            " group a function of the filesystem at write time rather than of the commit. On a partial"
            " store it neither refuses nor warns, and the result still reads as complete"
        ),
        "interim": (
            "counts are measured; the classification of each varying group into legitimate scope or"
            " silent shrink is NOT done and is not claimed here"
        ),
        "call_sites": {
            "files_entry_calls_total": len(calls),
            "exists_filtered": len(filtered),
            "the_brief_said": 17,
            "why_the_brief_undercounted": (
                "17 is what a one-line grep for files_entry.*exists() sees. The AST sees 18:"
                " scripts/trio_q100_phase.py:72 spans four lines with the filter on line 74"
            ),
            "sites": filtered,
        },
        "committed_declarations": {
            "group_declarations": declarations,
            "distinct_group_names": len(groups),
            "constant_member_count": len(groups) - len(var),
            "varying_member_count": len(var),
        },
        "the_negative_the_lane_was_opened_on": {
            "group": "reader v1 DNase narrowPeak files",
            "results_declaring_it": len(reader),
            "member_counts_seen": sorted({n for _, n in reader}, key=lambda x: (x is None, x)),
            "reading": (
                "312 in every one of the 13: NO SHRINK ANYWHERE. The defect is real as a mechanism and"
                " so far UNREALISED for this group. The 13-against-312 reported on 2026-10-10 was a"
                " worktree of this lane's own making, not any artefact"
            ),
        },
        "groups_whose_count_varies": var,
        "a_varying_count_is_not_itself_a_defect": (
            "a result that covered 23 chromosomes should declare 23. What a count cannot do is tell"
            " that from a result that meant 24 and silently got 23"
        ),
        "does_not_claim": gc.DOES_NOT_CLAIM,
        "the_limit_that_would_flatter_us": (
            "a group CONSTANT at a shrunken value is indistinguishable here from one constant at its"
            " full value, so the constant-count groups are NOT shown safe."
            " test_PLANTED_a_group_constant_at_a_shrunken_value_is_INVISIBLE_here holds that as a test"
        ),
    }
    sources = [Path("data/results")]
    print("saved", save_result(RESULT, payload, manifest=manifest(sources)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
