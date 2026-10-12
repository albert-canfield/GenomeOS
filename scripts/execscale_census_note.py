#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""execscale census note: which committed programmes differ from the census, and by whose commit.

    uv run python scripts/execscale_census_note.py        # write data/results/execscale_census_note.json

WHAT THIS IS, AND WHAT IT IS NOT. `data/results/execscale_census.json` (`81d3949`, 16:52 on
2026-10-12) wrote down the committed state before the calibrated magnitude existed. It is NOT
rewritten here and not one byte of it is touched: a census is a record of one tree, and a record that
is edited when the tree moves is not a record. This is a SIBLING note beside it, written additively,
saying which of its 58 programmes the current HEAD no longer matches and which commit moved each one.

WHY IT IS NEEDED. `71f9d47` (19:20, the coordinator's restore of the `share:` generator, honouring a
registered falsifier) changed data/organisms/human/haematopoiesis.bio after the census was taken. A
reader who compares the census to HEAD therefore finds a difference, and nothing in the census says
whether that difference is a legitimate commit or a lane's uncommitted slip. This note answers that
question for every programme, from git, and nothing more.

WHAT IT CLAIMS: a comparison of bytes, and for each difference the commits that touched that path
after the census. NOTHING ABOUT THE MAGNITUDE, nothing about RNA, nothing fitted, no model request and
no network. It does not judge whether a change was right -- it names the commit so a reader can.

THE GROUP DIGEST, stated because it is the number a reader will trip on. The census's `group_sha256`
is EXACTLY reproducible from the 58 blobs at `81d3949` in the record's own order, so the record is
sound. Two things make it look otherwise, and both are recorded here rather than argued:

  1. the working tree's group digest differs, by the one programme `71f9d47` moved;
  2. the ORDER matters and the census's order is not a global sort. `programs()` concatenates
     per-tree sorted lists over three trees, while `manifest.files_entry` sorts its names before
     hashing, so the same 58 files give one digest in the census's order and another through
     `files_entry`. Anything that reaches a group digest by the sorted route computes a different
     number from a correct record.

`manifest.group_digest` and the duplicate at `scripts/execscale_census.py:125` AGREE, given the same
list in the same order; that is checked in tests/test_execscale_rescoped.py and recorded here as a
reading. Neither is repaired by this script: the duplicate's removal is on the timed list.

Input paths go to `input_entry` REPO-RELATIVE, so the record does not pin this machine's checkout.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.results import RESULTS_DIR, save_result

RESULT = "execscale_census_note"
CENSUS = "execscale_census"
#: the commit the census record was taken on, and the commit that moved a programme after it
CENSUS_COMMIT = "81d3949"
OWN_CODE = (
    "scripts/execscale_census_note.py",
    "scripts/execscale_census.py",
    "tests/test_execscale_rescoped.py",
)
REPO = Path(__file__).resolve().parent.parent


def _run(*args: str) -> str:
    return subprocess.run(list(args), cwd=REPO, capture_output=True, text=True, check=True).stdout.strip()


def group_digest(names: list[str]) -> str:
    """The group digest of the names given, IN THE ORDER GIVEN, by the engine's implementation.

    `manifest.group_digest` and not a third copy: this note is about which of two numbers a reader is
    holding, so it must not add another implementation to be wrong in.
    """
    digest, _total, _members = mf.group_digest(names)
    return digest


def comparison(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """Every programme of the census against HEAD and against the working tree, with its commits."""
    census = json.loads((results_dir / f"{CENSUS}.json").read_text())
    rows = census["programs"]
    names = [row["path"] for row in rows]
    moved: list[dict[str, Any]] = []
    uncommitted: list[str] = []
    for row in rows:
        path = row["path"]
        on_disk = (REPO / path).read_bytes()
        disk_sha = hashlib.sha256(on_disk).hexdigest()
        head_sha = hashlib.sha256(
            subprocess.run(["git", "show", f"HEAD:{path}"], cwd=REPO, capture_output=True, check=True).stdout
        ).hexdigest()
        if disk_sha != head_sha:
            uncommitted.append(path)
        if disk_sha == row["sha256"]:
            continue
        log = _run("git", "log", "--oneline", f"{CENSUS_COMMIT}..HEAD", "--", path)
        moved.append(
            {
                "path": path,
                "census_sha256": row["sha256"],
                "census_bytes": row["bytes"],
                "now_sha256": disk_sha,
                "now_bytes": len(on_disk),
                "commits_after_the_census": log.splitlines(),
                "is_committed": disk_sha == head_sha,
            }
        )
    return {
        "programs_in_the_census": len(rows),
        "programs_identical_to_the_census": len(rows) - len(moved),
        "programs_that_moved": moved,
        "programs_with_an_uncommitted_difference": uncommitted,
        "group_sha256_in_the_census": census["group_sha256"],
        "group_sha256_of_the_working_tree_in_the_census_order": group_digest(names),
        "group_sha256_of_the_working_tree_sorted_by_name": group_digest(sorted(names)),
        "head": _run("git", "rev-parse", "HEAD"),
    }


def readings() -> dict[str, Any]:
    """The three readings this note is written to carry, each a sentence a reader can check."""
    return {
        "the_census_record_is_sound": (
            "its group_sha256 is exactly reproducible from the 58 blobs at 81d3949 in the record's own"
            " order, and every per-file digest matches those blobs. The record is a record of one tree"
            " and that tree still answers for it. tests/test_execscale_rescoped.py::test_the_recorded_"
            "group_digest_is_REPRODUCIBLE_from_the_tree_it_was_taken_on is the check."
        ),
        "the_two_group_digest_implementations_agree": (
            "manifest.group_digest and the duplicate at scripts/execscale_census.py:125 return the same"
            " digest for the same list in the same order, on the committed group of 58 programmes."
            " The duplicate is NOT removed here: removing it is a removal inside a check_staged window,"
            " so it is reported for the timed list instead."
        ),
        "a_group_digest_carries_the_LIST_as_well_as_the_CONTENT": (
            "both implementations hash name + b'\\0' + bytes per member IN THE ORDER GIVEN and neither"
            " sorts internally, so membership and order move the digest while no per-file digest moves."
            " The 58 per-file digests alone therefore do not carry the census's claim: a programme"
            " added, removed or reordered is invisible to them. And because manifest.files_entry SORTS"
            " its names before hashing while execscale_census.programs() concatenates per-tree sorted"
            " lists over three trees, the two routes to a group's digest disagree by construction --"
            " which is how a checker can report that only the group digest differs on a tree whose"
            " every per-file digest matches."
        ),
    }


def manifest() -> dict[str, Any]:
    inputs = [mf.input_entry(Path("data/results") / f"{CENSUS}.json", partition=None)]
    inputs.append(
        mf.files_entry(
            "the committed .bio programmes the census covers, as this note reads them from disk",
            [row["path"] for row in json.loads((RESULTS_DIR / f"{CENSUS}.json").read_text())["programs"]],
        )
    )
    return {
        "sources": [
            {"accession": "data/results/execscale_census.json", "version": "81d3949, unmodified"},
            {"accession": "this checkout's git history", "version": "HEAD at the time of writing"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "census_commit": CENSUS_COMMIT,
            "method": (
                "sha256 per programme on disk against the census row and against the HEAD blob; for"
                " each difference, git log CENSUS_COMMIT..HEAD on that path; group digests by"
                " manifest.group_digest in two orders"
            ),
        },
        "exclusions": [
            "the census record is not rewritten and not one byte of it is touched",
            "no model request, no key, no network, no money; nothing is fitted and no RNA is read",
            "no judgement of whether a change was right: the commit is named so a reader can judge",
        ],
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        "partitions": {
            "none": "one population, the 58 programmes of the census; no partition is taken or implied"
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args(argv)
    payload = {
        "question": (
            "which of the census's 58 committed programmes does HEAD no longer match, and which"
            " commit moved each one?"
        ),
        "claim": (
            "a byte comparison and a git attribution, nothing else: every programme of"
            " data/results/execscale_census.json against its bytes on disk now, and for each"
            " difference the commits after 81d3949 that touched that path. The census is not"
            " rewritten. NOTHING about the magnitude, nothing about RNA, nothing fitted."
        ),
        "does_not_claim": (
            "it says nothing about whether a change was correct, nothing about RNA, and nothing about"
            " the magnitude the census was taken for; it fits nothing and reads no measurement"
        ),
        "readings": readings(),
        **comparison(),
    }
    print("saved", save_result(RESULT, payload, manifest=manifest()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
