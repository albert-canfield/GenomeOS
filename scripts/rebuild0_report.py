#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-rebuild0: the committed record of the 0-difference rebuild of M2h, and two tool findings.

    uv run python scripts/rebuild0_report.py

Written because a verdict that exists only in a session is not a record. The coordinator ruled this
belongs in data/results with provenance rather than in a commit message, since it is a measurement with
inputs and a verdict and should be queryable and rebuildable.

INPUT PATHS ARE REPO-RELATIVE ON PURPOSE: `genomeos.manifest.input_entry` records `str(path)` verbatim
while its sibling `files_entry` relativises, so an absolute path would record this machine's checkout
and make this result unrebuildable elsewhere -- the defect that leaves committed results unrebuildable
(lane-silencerebuild, 2026-10-03). manifest.py is inside the paid study's frozen closure, so a relative
path at the call site is the only way not to add another.
"""

from __future__ import annotations

from pathlib import Path

from genomeos import manifest as mf
from genomeos.results import save_result

RESULT = "rebuild0_report"
MAGCAL_RESULT = Path("data/results/magcal_result.json")
MAXFORM_REG = Path("data/results/maxform_registration.json")

OWN_CODE = ("scripts/rebuild0_report.py",)


def manifest() -> dict:
    inputs = [mf.input_entry(MAGCAL_RESULT, partition=None), mf.input_entry(MAXFORM_REG, partition=None)]
    return {
        "sources": [
            {
                "accession": "data/results/magcal_result.json (M2h, written at fbce775)",
                "version": "pinned by sha256",
            },
            {
                "accession": "data/results/maxform_registration.json (written at c7eda3d)",
                "version": "pinned by sha256",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "tool": "scripts/manifest_rebuild.py",
            "method": "a clean worktree at each artefact's own recorded code.git_sha, stores and declared"
            " machine-local inputs linked read-only, fresh venv from the committed uv.lock, offline",
        },
        "exclusions": [
            "no RNA and no model was run here: this records the verdicts of two rebuilds and two defects"
            " found in the rebuild tool itself"
        ],
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        "partitions": {
            "none": "two artefacts reported side by side, not pooled; no partition is taken and the"
            " absence of one is not a claim of independence"
        },
    }


def payload() -> dict:
    return {
        "question": "does M2h's result rebuild from committed bytes, and does maxform's registration?",
        "the_sole_positive_reproduces": {
            "artefact": "data/results/magcal_result.json",
            "added_in_commit": "ac5580b",
            "its_own_declared_code_git_sha": "fbce775",
            "rebuilt": True,
            "rebuilt_sha256": "10c5c3250f0484a93e364557e7c5c2bc1c9be4ff6ddf4121ff6875f989e1ed35",
            "differences": 0,
            "must_hold_failures": 0,
            "claim_as_it_may_be_stated": (
                "ZERO differences over 484 COMPARED leaves, with 4 exempted leaves named below and 16"
                " never compared. NOT 'reproduces exactly'"
            ),
            "why_it_matters": (
                "every form in this series so far reads WORSE than no rule (7, 7d, 7h, and the capped"
                " sum at 0 genes). ac5580b is the ONLY result that reverses the pattern, so this is the"
                " check on the sole positive rather than bookkeeping"
            ),
            "registration_before_run_HOLDS_BY_COMMIT_ORDER": (
                "the result's own declared code.git_sha IS fbce775, the commit that added"
                " magcal_registration.json, so the registration was committed before the code that"
                " wrote the result ran. Contrast maxform below, whose census, registration and result"
                " landed in one commit"
            ),
        },
        "leaf_accounting": {
            "total": 500,
            "compared": 484,
            "not_compared": 16,
            "the_sum_closes": "484 + 16 = 500",
            "leaves_set_aside": 4,
            "set_aside_is_NOT_an_addend": (
                "the 4 are a SUBSET of the 484: leaves that were compared, differed, and were"
                " exempted. 484 + 4 = 488 adds a subset back onto its own superset, which is why that"
                " sum appears not to close"
            ),
            "the_16_by_class": {
                "result_manifest_code_block": 14,
                "date": 1,
                "timing": 1,
                "model_dependencies": 0,
                "resource_fields": 0,
            },
        },
        "the_4_set_aside_by_exact_path": [
            {
                "path": "/result_manifest/code_cleanliness/dirty",
                "class": "environment",
                "committed": True,
                "rebuilt": False,
            },
            {
                "path": "/result_manifest/code_cleanliness/foreign_uncommitted_code",
                "class": "environment",
                "committed": "7 items",
                "rebuilt": 0,
            },
            {
                "path": "/result_manifest/traced_inputs/opens",
                "class": "run, KNOWN-LOSSY",
                "committed": 395,
                "rebuilt": 816,
            },
            {"path": "/seconds", "class": "timing", "value": 8.5},
        ],
        "written_with_a_dirty_tree": {
            "code_cleanliness_dirty": True,
            "foreign_uncommitted_files_at_write_time": 7,
            "note": (
                "correctly exempted as an environment field -- it describes the tree the run happened"
                " in and not anything the run computed -- and recorded here where a later reader sees it"
            ),
        },
        "maxform_registration_COULD_NOT_REACH_A_VERDICT": {
            "artefact": "data/results/maxform_registration.json",
            "rebuilt": False,
            "inputs_declared": 62,
            "inputs_checked": 61,
            "inputs_unchecked": 1,
            "unavailable": ["input data/results/maxform_census.json is absent"],
            "cause_from_the_record": (
                "the registration declares code.git_sha c7eda3d (21:57:20) and dirty_result_paths"
                " ['data/results/maxform_census.json']; that census first entered git at b3c560b"
                " (22:07:45). link_machine_local_inputs cannot supply it either, because it links only"
                " git-IGNORED declared inputs and the census is a TRACKED file. So the declared sha"
                " PREDATES its input's existence in git and no checkout of that sha can rebuild it"
            ),
            "this_is_the_finished_form": (
                "'could not reach a verdict, and here is the dependency it named' -- not a comparison."
                " An earlier attempt of mine reported 53 differences; those measured the METHOD (a raw"
                " worktree plus link_stores_read_only.py, which links the three stores and not the"
                " manifest's declared machine-local inputs) and are WITHDRAWN"
            ),
        },
        "finding_1_a_known_lossy_exemption_is_hiding_a_2x_difference": {
            "leaf": "/result_manifest/traced_inputs/opens",
            "committed": 395,
            "rebuilt": 816,
            "the_tools_own_words": (
                "KNOWN-LOSSY EXEMPTION, not noise: setting it aside WILL HIDE a genuine change in the"
                " set of files this writer opens outside data/"
            ),
            "reading": (
                "the exemption is doing exactly what its own text warns of. Under the rule of"
                " 2026-10-02 a known-lossy leaf must be COMPARED WITHIN A STATED TOLERANCE, not"
                " skipped; only ENVIRONMENT_FIELDS and run-resource keys are exempt, by exact key. It"
                " names no path so it cannot say WHICH files, and a 2x difference in opens is not noise"
            ),
            "not_repaired_here": (
                "manifest_rebuild.py is another lane's verification tool: reported, not patched"
            ),
        },
        "finding_2_the_reconciles_field_cannot_be_false": {
            "field": "leaves.reconciles",
            "code": "total == compared + (total - compared)",
            "why_vacuous": "the right-hand side IS total by arithmetic, so the field is True for any input",
            "measured": (
                "returns True for total=10 compared=99, for total=0 compared=7, and for total=3"
                " compared=-5 -- including cases a real check must refuse"
            ),
            "consequence": (
                "the accounting does close, but BY DEFINITION rather than by verification, and"
                " 'reconciles: true' is not evidence of anything. I had cited it as evidence"
            ),
            "the_non_vacuous_check_it_reaches_for": (
                "that the ENUMERATED dropped-leaf lists -- ignored_keys, the manifest code block,"
                " model_dependencies, the timing fields, the resource fields -- SUM to not_compared"
            ),
            "found_by": "the coordinator, genomeos-99; verified here first-hand rather than accepted",
            "not_repaired_here": "another lane's verification tool: reported, not patched",
        },
    }


def main() -> int:
    print("saved", save_result(RESULT, payload(), manifest=manifest()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
