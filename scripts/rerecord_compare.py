# SPDX-License-Identifier: AGPL-3.0-or-later
"""Compare a result re-recorded under a new name against the committed file it is written beside.

    uv run python scripts/rerecord_compare.py NEW.json OLD.json [--json]

Four results were re-recorded on 2026-10-02 under new names because their manifests declared groups of
input files under a label with no member list, which no second environment can open or hash, while the
committed bytes could not be replaced: each one's sha256 is a declared input of other committed results.
The acceptance test of such a re-record is that **no figure moved**, and the way to show that is leaf by
leaf rather than by quoting a chosen list of headline numbers: a comparison that checks the numbers
somebody thought to name says nothing about the ones they did not.

So this reads both files with `scripts/manifest_rebuild.py`'s own `comparable`, `diff`, `leaves` and
`leaf_reconciliation` -- the same functions that decide a rebuild's verdict, not a second opinion about
them -- and sorts every difference into one of two piles:

* the fields a re-record under a new name **must** change, each by exact path or by one documented
  prefix (RENAME_FIELDS below): the name itself, the `supersedes` block that says which committed file
  this one is written beside, the `code_cleanliness` block that describes the tree the run happened in,
  and the `members`/`group` keys whose absence is the whole reason for the re-record;
* everything else, which is a figure that moved. One is enough to fail: the exit status is 1 and the
  caller is expected to stop and report it rather than adjust anything to make the two agree.

`MUST_HOLD` (own code committed, no foreign uncommitted code on the counting path) is checked on both
files and reported, because `code_cleanliness` is exempted above and that exemption must not be able to
hide a run whose own code was uncommitted.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import manifest_rebuild as mr  # noqa: E402

from genomeos import manifest as mf  # noqa: E402

#: The only differences a re-record under a new name may carry, by exact path or by prefix. A prefix
#: entry ends in `/` and covers the keys below that path; `/result` is the name itself. Anything not
#: matched here is reported as a figure that moved, so this list is the whole of what is forgiven and
#: adding to it is a visible change rather than a judgement call at reading time.
RENAME_FIELDS: tuple[str, ...] = (
    "/result",
    "/result_manifest/supersedes",
    "/result_manifest/supersedes/",
    "/result_manifest/code_cleanliness",
    "/result_manifest/code_cleanliness/",
    "/code_cleanliness",
    "/code_cleanliness/",
)

#: Keys that may appear on a declared input in the new file and not in the old one: `files_entry` did
#: not name a group's members until 2026-10-02, and naming them is the point of the re-record.
MEMBER_KEYS = ("members", "group")


def _input_member_key(path: str) -> bool:
    """Whether a difference path is a declared input gaining its member list (or the marker that says
    it is a group), at any depth below that key: `/result_manifest/inputs[7]/members[3]/sha256`."""
    for container in ("/result_manifest/inputs[", "/inputs["):
        if not path.startswith(container):
            continue
        rest = path[path.index("]", len(container)) + 1 :] if "]" in path else ""
        for key in MEMBER_KEYS:
            if rest == f"/{key}" or rest.startswith(f"/{key}/") or rest.startswith(f"/{key}["):
                return True
    return False


def is_rename_field(path: str) -> bool:
    """Whether this difference is one the rename itself accounts for."""
    if path in RENAME_FIELDS:
        return True
    if any(p.endswith("/") and path.startswith(p) for p in RENAME_FIELDS):
        return True
    return _input_member_key(path)


def classify(differences: list[str]) -> tuple[list[str], list[str]]:
    """(figures that moved, differences the rename accounts for). The path is everything before the
    first `:` of a `diff()` line, as `environment_differences` reads it too."""
    moved, accounted = [], []
    for d in differences:
        (accounted if is_rename_field(d.split(":")[0]) else moved).append(d)
    return moved, accounted


def group_accounting(payload: dict[str, Any]) -> dict[str, Any]:
    """What this manifest's inputs declare: how many entries, how many are groups, how many of those
    name their members, and how many member files are named with a sha256 of their own."""
    inputs = ((payload.get(mf.KEY) or {}).get("inputs")) or []
    groups = [i for i in inputs if isinstance(i, dict) and mr.is_group(i)]
    labelled = [
        i
        for i in inputs
        if isinstance(i, dict)
        and not mr.is_group(i)
        and "files" in i
        and not Path(str(i.get("path"))).exists()
    ]
    members = [m for g in groups for m in (g.get("members") or [])]
    return {
        "inputs_declared": len(inputs),
        "groups": len(groups),
        "groups_naming_their_members": sum(1 for g in groups if g.get("members")),
        "member_files_named": len(members),
        "member_files_with_their_own_sha256": sum(
            1 for m in members if isinstance(m, dict) and m.get("sha256") and m.get("bytes") is not None
        ),
        "entries_that_are_a_label_no_path_resolves": [str(i.get("path")) for i in labelled],
        "declared_paths": len(mr.declared_paths(inputs)),
    }


def compare(new: Path, old: Path) -> dict[str, Any]:
    """The whole report. `old` is the committed file, `new` the one written beside it."""
    a, b = json.loads(old.read_text()), json.loads(new.read_text())
    differences = mr.diff(mr.comparable(a), mr.comparable(b))
    moved, accounted = classify(differences)
    return {
        "committed": {
            "file": old.as_posix(),
            "result": a.get("result"),
            "date": a.get("date"),
            "leaves": mr.leaf_reconciliation(a),
            "inputs": group_accounting(a),
        },
        "re_recorded": {
            "file": new.as_posix(),
            "result": b.get("result"),
            "date": b.get("date"),
            "leaves": mr.leaf_reconciliation(b),
            "inputs": group_accounting(b),
        },
        "leaves_compared": mr.leaves(mr.comparable(a)),
        "differences": len(differences),
        "figures_that_moved": moved,
        "accounted_for_by_the_rename": accounted,
        "timing_fields_ignored": mr.timing_paths(a),
        "resource_fields_ignored": mr.resource_paths(a),
        "must_hold_failures": mr.must_hold_failures(a, b),
        "no_figure_moved": not moved,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("new", type=Path, help="the result re-recorded under a new name")
    ap.add_argument("old", type=Path, help="the committed result it is written beside")
    ap.add_argument("--json", action="store_true", help="the whole report as JSON")
    args = ap.parse_args(argv)

    r = compare(args.new, args.old)
    if args.json:
        print(json.dumps(r, indent=2))
    else:
        for side in ("committed", "re_recorded"):
            s = r[side]
            lv, gi = s["leaves"], s["inputs"]
            print(f"{side:12s} {s['result']} ({s['date']}) {s['file']}")
            print(f"             leaves {lv['compared']} compared of {lv['total']}, {lv['not_compared']} not")
            print(
                f"             inputs {gi['inputs_declared']} declared, {gi['groups']} groups, "
                f"{gi['groups_naming_their_members']} naming members, {gi['member_files_named']} member "
                f"files ({gi['member_files_with_their_own_sha256']} with sha256 and bytes), "
                f"{gi['declared_paths']} resolvable paths"
            )
            if gi["entries_that_are_a_label_no_path_resolves"]:
                print(
                    f"             {len(gi['entries_that_are_a_label_no_path_resolves'])} entries are a "
                    f"label no path resolves: {gi['entries_that_are_a_label_no_path_resolves']}"
                )
        print(f"leaves compared      {r['leaves_compared']}")
        print(f"differences          {r['differences']}")
        print(f"  accounted for      {len(r['accounted_for_by_the_rename'])} (the rename itself)")
        print(f"  figures that moved {len(r['figures_that_moved'])}")
        for d in r["figures_that_moved"][:50]:
            print(f"    {d}")
        print(f"must_hold_failures   {len(r['must_hold_failures'])}")
        for d in r["must_hold_failures"]:
            print(f"    {d}")
        print(f"no_figure_moved      {r['no_figure_moved']}")
    return 0 if r["no_figure_moved"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
