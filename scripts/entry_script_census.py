# SPDX-License-Identifier: AGPL-3.0-or-later
"""Which committed results were written by a script outside this repository. READ-ONLY.

    uv run python scripts/entry_script_census.py              # the README-cited population
    uv run python scripts/entry_script_census.py --all        # every data/results/*.json as well
    uv run python scripts/entry_script_census.py --json
    uv run python scripts/entry_script_census.py --write      # save it as a result (--write implies --all)

Why (2026-10-02, lane-entrypoint). `genomeos.manifest.code_cleanliness` compared git's dirty list
against the paths the CALLER DECLARED and never looked at the script that was running, so a result
could assert `own_code_is_committed: True` while the bytes came from a file in no commit. Measured on
`data/results/organised_chr21.json`
(working tree at d4d0449, 19:41; superseded by e1def2e, 20:41): `code.argv[0]` was
`/private/tmp/.../scratchpad/write.py`, the
counting path holds 39 `genomeos/` modules and no script, and the manifest says the code was
committed. The check is fixed for new writes; this census says how far the old ones reach.

UNTIL 2026-10-03 THIS SCRIPT ONLY PRINTED, so its numbers were in no result and no row could quote
them; `--write` saves the whole census through `genomeos.results.save_result`, with a cleanliness block,
this committed script counted as the entry, and every examined result hashed as a declared input.

THE 918 RESULTS THAT RECORD NO `code.argv` ARE NOT ONE FINDING, and splitting them is why this lane
exists rather than a tidy-up. A result written BEFORE `manifest._argv` existed records nothing because
there was nothing to record: that is history, and calling it a breach would be false. A result written
AFTER it is a live gap. So `split_by_boundary` finds the commit that introduced the field by `git log`
on the FUNCTION rather than by guessing a date, places each unrecorded result by ITS OWN COMMITTED
HISTORY and never by its `date` field (which the writer supplies and could be anything), and names a
result it cannot place in its own class instead of folding it into whichever side is smaller.

NOTHING IS REWRITTEN AND NO RESULT'S BYTES ARE TOUCHED. Several results' sha256 ARE registrations in
`data/results/manifest_headlines.json`, and rewriting one would break a pin. The census reads and
reports; what to do about what it finds is not its decision.

HOW THE VERDICT IS READ, and it is read from the record rather than recomputed. `manifest.code.argv`
is written by `manifest._argv`, which makes argv[0] RELATIVE to the repository root when it is under
it and leaves it absolute when it is not. So an absolute argv[0] is itself the evidence that the
script was outside the repository at write time -- no present-tree lookup is involved, and a script
that has since been committed, moved or deleted cannot change the answer. The one case where that
inference does not hold is named rather than folded in: an absolute path that lies UNDER this
checkout's root was inside the repository and recorded absolutely (a writer that passed a resolved
`__file__` before `_argv` relativised, or a differently-rooted checkout), and it is counted apart
under `inside_recorded_absolute` instead of being read as outside.

"README-cited" is the definition `tests/test_headline_registry_matches_readme.py` already uses, not a
new one: the entries of `manifest_headlines.json` (`rebuilt` + `pending`), which is the record of what
the README's claims rest on, together with every `*.json` the README names that is a file in
data/results. Both halves are reported, and a registered name with no file on this disk is reported
as `absent` rather than skipped.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULTS = Path("data/results")
README = Path("README.md")
HEADLINES = RESULTS / "manifest_headlines.json"
ROOT = Path(__file__).resolve().parent.parent

RESULT = "entry_script_census"
#: This script, as the entry whose transitive import closure is the counting path of its own result
#: (genomeos.manifest.counting_path). Named rather than derived from __file__ so the closure is the
#: same however the script is invoked.
ENTRY = "scripts/entry_script_census.py"
OWN_CODE = (ENTRY, "tests/test_entry_script_census.py")

#: The field whose arrival is the boundary, and the one function that has ever written it.
ARGV_FIELD = "result_manifest.code.argv, written by genomeos.manifest._argv"

#: The verdict for each recorded argv[0]. `outside_repository` is the one the census was asked for.
VERDICTS = {
    "outside_repository": "argv[0] is an absolute path that is not under this checkout: the script "
    "that wrote the bytes was outside the repository and no commit holds it",
    "inside_recorded_absolute": "argv[0] is absolute but under this checkout's root, so it was "
    "inside the repository and merely recorded absolutely; NOT counted as outside",
    "inside_repository": "argv[0] is a repository-relative path, which is what manifest._argv writes "
    "for a script under the root",
    "inline_source": "argv[0] is '-c': the source was never in a file",
    "stdin": "argv[0] is '-': the source was piped in",
    "test_runner": "argv[0] is a test runner's entry point",
    "unrecorded": "the manifest records no code.argv, so nothing says which script ran",
    "absent": "the result is registered as README-cited but is not a file on this disk",
}

_RUNNERS = ("pytest", "py.test", "unittest")


def readme_cited() -> tuple[list[str], dict[str, list[str]]]:
    """The README-cited result names, and where each citation comes from."""
    where: dict[str, list[str]] = {}
    try:
        record = json.loads(HEADLINES.read_text())
    except OSError:
        record = {}
    for entry in list(record.get("rebuilt") or []) + list(record.get("pending") or []):
        name = entry.get("result")
        if name:
            where.setdefault(name, []).append(f"manifest_headlines ({entry.get('quoted_in')})")
    try:
        text = README.read_text()
    except OSError:
        text = ""
    for filename in sorted(set(re.findall(r"[a-z0-9_]+\.json", text))):
        if (RESULTS / filename).is_file():
            where.setdefault(filename.removesuffix(".json"), []).append("README.md names the file")
    return sorted(where), where


def verdict_for(argv0: Any) -> str:
    if argv0 is None:
        return "unrecorded"
    argv0 = str(argv0)
    if argv0 == "-c":
        return "inline_source"
    if argv0 == "-":
        return "stdin"
    if Path(argv0).name in _RUNNERS or argv0.endswith(("pytest/__main__.py", "unittest/__main__.py")):
        return "test_runner"
    if not Path(argv0).is_absolute():
        return "inside_repository"
    try:
        Path(argv0).relative_to(ROOT)
    except ValueError:
        return "outside_repository"
    return "inside_recorded_absolute"


def bytes_state(p: Path) -> str:
    """Whether the result FILE's own bytes are committed: `committed`, `modified`, `untracked`, `unknown`.

    A finding needs this column. The 24 `organised_chr*.json` files this census flags are WORKING-TREE
    modifications whose committed versions are older, legacy-shaped results with no `result_manifest`
    at all, and another lane is regenerating them from a committed entry point. A file about to be
    discarded is not a published result that certified itself falsely, and a census that printed the
    two the same way would invite its reader to count 24 published breaches where there are none.
    """
    try:
        listed = subprocess.run(
            ["git", "ls-files", "-z", "--", str(p)], capture_output=True, text=True, timeout=30, check=True
        ).stdout
        if str(p) not in [x for x in listed.split("\0") if x]:
            return "untracked"
        status = subprocess.run(
            ["git", "status", "--porcelain", "-z", "--", str(p)],
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    return "modified" if [x for x in status.split("\0") if x] else "committed"


def read_one(name: str) -> dict[str, Any]:
    """One result's recorded entry script. Reads the file and nothing else; never writes."""
    p = RESULTS / f"{name}.json"
    if not p.is_file():
        return {"result": name, "verdict": "absent", "argv0": None}
    try:
        d = json.loads(p.read_text())
    except (OSError, ValueError) as e:
        return {"result": name, "verdict": "unrecorded", "argv0": None, "unreadable": str(e)}
    m = d.get("result_manifest") or {}
    code = m.get("code") or {}
    argv = code.get("argv")
    argv0 = argv[0] if isinstance(argv, list) and argv else None
    clean = m.get("code_cleanliness") or {}
    return {
        "result": name,
        "verdict": verdict_for(argv0),
        "argv0": argv0,
        "claims_own_code_is_committed": clean.get("own_code_is_committed"),
        "carries_entry_script": isinstance(clean.get("entry_script"), dict),
        "counting_path_count": clean.get("counting_path_count"),
        "bytes_state": bytes_state(p),
    }


def _git(*args: str) -> str | None:
    """git, read-only, from the repository root. None when git cannot answer, never a guess."""
    try:
        out = subprocess.run(
            ["git", *args], capture_output=True, text=True, timeout=120, check=True, cwd=ROOT
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout


def _is_ancestor(older: str, newer: str) -> bool | None:
    """Whether `older` is an ancestor of `newer`; None when git cannot answer."""
    try:
        done = subprocess.run(
            ["git", "merge-base", "--is-ancestor", older, newer],
            capture_output=True,
            text=True,
            timeout=120,
            cwd=ROOT,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return True if done.returncode == 0 else (False if done.returncode == 1 else None)


def _oldest_commit_line(text: str | None) -> tuple[str, str, str] | None:
    """The oldest `%H\t%ad\t%s` line in a git log, ignoring any diff the log also printed."""
    if not text:
        return None
    rows = [ln.split("\t", 2) for ln in text.splitlines() if re.match(r"^[0-9a-f]{40}\t", ln)]
    rows = [r for r in rows if len(r) == 3]
    return (rows[-1][0], rows[-1][1], rows[-1][2]) if rows else None


def argv_boundary() -> dict[str, Any]:
    """The commit that made a result able to record `code.argv`, read from git and not from a date.

    Asked two ways that have to agree, because the whole before/after split rests on this one answer.
    `git log -L :_argv:genomeos/manifest.py` follows the FUNCTION through the file's history and its
    oldest entry is the commit that wrote it; the pickaxe on the call site answers the same question
    from the other end. When the two disagree there is no boundary and the split says so rather than
    picking the one that reads better -- a wrong boundary would move results between "history" and
    "a live gap", which is exactly the distinction this census exists to draw.
    """
    by_function = _oldest_commit_line(
        _git("log", "-L", ":_argv:genomeos/manifest.py", "--format=%H\t%ad\t%s", "--date=iso")
    )
    by_call_site = _oldest_commit_line(
        _git(
            "log",
            "-S",
            '"argv": _argv(',
            "--format=%H\t%ad\t%s",
            "--date=iso",
            "--",
            "genomeos/manifest.py",
        )
    )
    agree = bool(by_function and by_call_site and by_function[0] == by_call_site[0])
    found = {
        "git log -L :_argv:genomeos/manifest.py (oldest)": by_function[0] if by_function else None,
        "git log -S '\"argv\": _argv(' -- genomeos/manifest.py (oldest)": (
            by_call_site[0] if by_call_site else None
        ),
    }
    return {
        "field": ARGV_FIELD,
        "sha": by_function[0] if agree else None,
        "date": by_function[1] if agree else None,
        "subject": by_function[2] if agree else None,
        "found_by": found,
        "the_two_queries_agree": agree,
        "how": (
            "read from git log on the function and on its call site, never from a date; when the two "
            "disagree there is no boundary and every result is left unplaceable with that reason"
        ),
    }


#: What each side of the boundary means. `after_the_field_existed` is the only class that is a defect.
SIDES = {
    "before_the_field_existed": "the commit that added this result is an ancestor of the boundary "
    "commit, so the result was written when nothing could have recorded code.argv. HISTORY, not a "
    "breach: there was no field to fill.",
    "after_the_field_existed": "the commit that added this result is a descendant of the boundary "
    "commit, so the field existed when the result was written and is empty anyway. THE FINDING: a "
    "live gap.",
    "added_in_the_boundary_commit": "the boundary commit itself added this result, so the field and "
    "the result entered the record together and neither side can be shown. UNPLACEABLE.",
    "never_committed": "no commit holds this file, so it has no committed history to place it by and "
    "its `date` field is the writer's word. UNPLACEABLE, and it is also not a published result.",
    "unplaceable_no_add_commit": "git tracks the file but reports no commit that added it, so the "
    "earliest commit containing it could not be read. UNPLACEABLE.",
    "unplaceable_unrelated_history": "the commit that added this result is neither an ancestor nor a "
    "descendant of the boundary commit. UNPLACEABLE.",
    "unplaceable_no_boundary": "the boundary itself could not be read, so no result can be placed.",
}
#: Which classes are unplaceable. Named here rather than inferred, so none is ever folded into
#: whichever side happens to be smaller.
UNPLACEABLE = (
    "added_in_the_boundary_commit",
    "never_committed",
    "unplaceable_no_add_commit",
    "unplaceable_unrelated_history",
    "unplaceable_no_boundary",
)


def _tree_paths(ref: str) -> set[str] | None:
    out = _git("ls-tree", "-r", "--name-only", ref, "--", str(RESULTS))
    return None if out is None else {ln for ln in out.splitlines() if ln}


def _tracked_paths() -> set[str]:
    out = _git("ls-files", "--", str(RESULTS))
    return {ln for ln in (out or "").splitlines() if ln}


def _first_add_commit(rel: str) -> str | None:
    out = _git("log", "--reverse", "--diff-filter=A", "--format=%H", "--", rel)
    lines = [ln for ln in (out or "").splitlines() if ln]
    return lines[0] if lines else None


def split_by_boundary(names: list[str]) -> dict[str, Any]:
    """Each unrecorded result placed before or after the commit that made `code.argv` recordable.

    The side is read from the RESULT'S OWN COMMITTED HISTORY and never from its `date` field, which
    the writer supplies and could be anything. A file present in the boundary commit's PARENT tree was
    committed before the boundary, which is the same answer `git log --diff-filter=A` gives for it and
    is one git call for the whole set instead of one per file; anything else is asked per file.
    """
    boundary = argv_boundary()
    sha = boundary["sha"]
    tracked = _tracked_paths()
    parent_tree = _tree_paths(f"{sha}^") if sha else None
    boundary_tree = _tree_paths(sha) if sha else None
    if parent_tree is None and sha:
        parent_tree = set()  # a root commit has no parent: nothing was committed before it

    sides: dict[str, list[str]] = {}
    for name in names:
        rel = f"{RESULTS.as_posix()}/{name}.json"
        if not sha or boundary_tree is None:
            side = "unplaceable_no_boundary"
        elif rel not in tracked:
            side = "never_committed"
        elif rel in (parent_tree or set()):
            side = "before_the_field_existed"
        elif rel in boundary_tree:
            side = "added_in_the_boundary_commit"
        else:
            add = _first_add_commit(rel)
            if add is None:
                side = "unplaceable_no_add_commit"
            elif _is_ancestor(sha, add):
                side = "after_the_field_existed"
            elif _is_ancestor(add, sha):
                side = "before_the_field_existed"
            else:
                side = "unplaceable_unrelated_history"
        sides.setdefault(side, []).append(name)

    counts = {k: len(v) for k, v in sorted(sides.items())}
    unplaceable = {k: v for k, v in sorted(sides.items()) if k in UNPLACEABLE}
    return {
        "boundary": boundary,
        "examined": len(names),
        "counts": counts,
        "the_finding": {
            "class": "after_the_field_existed",
            "count": counts.get("after_the_field_existed", 0),
            "results": sorted(sides.get("after_the_field_existed", [])),
            "reading": "results written after the field existed that record no entry script at all: "
            "the field was there to fill and is empty, so nothing names the code that made the bytes",
        },
        "the_history": {
            "class": "before_the_field_existed",
            "count": counts.get("before_the_field_existed", 0),
            "reading": "results written before the field existed. They record nothing because there "
            "was nothing to record; calling this a breach would be false.",
        },
        "the_unplaceable": {
            "count": sum(len(v) for v in unplaceable.values()),
            "by_class": {k: sorted(v) for k, v in unplaceable.items()},
            "reading": "neither side can be shown from committed history, so each is named in its own "
            "class and none is folded into whichever side is smaller",
        },
        "sides_are": SIDES,
        "the_side_is_read_from": (
            "the result's own committed history (the boundary commit's parent tree, or git log "
            "--diff-filter=A for the file), never its `date` field, which the writer supplies"
        ),
    }


def census(names: list[str]) -> dict[str, Any]:
    rows = [read_one(n) for n in names]
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["verdict"]] = counts.get(row["verdict"], 0) + 1
    outside = [r for r in rows if r["verdict"] == "outside_repository"]
    return {
        "examined": len(rows),
        "counts": counts,
        "outside_repository": outside,
        "outside_repository_count": len(outside),
        "falsely_certified": [r for r in outside if r.get("claims_own_code_is_committed") is True],
        "falsely_certified_and_committed": [
            r
            for r in outside
            if r.get("claims_own_code_is_committed") is True and r.get("bytes_state") == "committed"
        ],
        "rows": rows,
        "verdicts_are": VERDICTS,
        "nothing_was_written": "this script reads result files and writes none of them",
    }


#: The registered readings this artefact carries, word for word. They are quoted and never restated:
#: a census that recomputed them would be free to report whichever number its own population gave.
READINGS = {
    "falsely_certified": "0 of 1,109 committed results falsely certified",
    "falsely_certified_means": (
        "no committed result asserts `own_code_is_committed: True` while naming a script outside the "
        "repository. It does NOT mean every result is rebuildable."
    ),
    "organised_chr": (
        "the 24 `organised_chr*` results that prompted the entry-script work were NEVER PUBLISHED in "
        "that state: every committed version was legacy-shaped with no `result_manifest` at all, and "
        "the census's own `bytes_state` column recorded all 24 as `modified`, none `committed`."
    ),
    "organised_chr_do_not_say": (
        'do not repeat the stronger claim that they "asserted own_code_is_committed: True", which '
        "was a working-tree state superseded at e1def2e"
    ),
}


def population(names: list[str]) -> dict[str, Any]:
    """What the examined set is, counted rather than described.

    The printed table called this "every committed result" and the set is a glob of the WORKING TREE,
    so a result no commit holds is in it. The denominator is named here because the before/after split
    turns on exactly that difference: a file git does not track has no committed history to place it by.
    """
    tracked = _tracked_paths()
    rels = [f"{RESULTS.as_posix()}/{n}.json" for n in names]
    inside = [r for r in rels if r in tracked]
    return {
        "glob": "data/results/*.json in the working tree, listed once at the start of this run",
        "files": len(rels),
        "tracked_by_git": len(inside),
        "untracked": len(rels) - len(inside),
        "it_is_not": (
            "'every committed result': the glob is the working tree, and an untracked result is in it"
        ),
        "the_count_moves": (
            "other lanes write results into this directory while the census runs, so two runs minutes "
            "apart give different totals; the number is only ever the number at this run's git_sha"
        ),
    }


def stdin_written(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The results whose argv[0] is '-': the source was piped in and exists in no file."""
    found = [r for r in rows if r["verdict"] == "stdin"]
    return {
        "count": len(found),
        "results": [
            {
                "result": r["result"],
                "argv0": r["argv0"],
                "claims_own_code_is_committed": r["claims_own_code_is_committed"],
                "bytes_state": r["bytes_state"],
            }
            for r in found
        ],
        "reading": (
            "NEITHER claims `own_code_is_committed`, so neither is falsely certified; neither is "
            "rebuildable, because the source was piped in and exists in no file"
        ),
        "not_this_lane_s": (
            "lane-headlinewriter is working on these two; this census counts them and touches nothing"
        ),
    }


def artefact(every: list[str] | None = None) -> dict[str, Any]:
    """The whole census as the result records it. Reads result files and git; writes none of them.

    `every` is the examined set, listed once by the caller so that the files the manifest HASHES are
    the same files these counts were taken from: the directory is written into by other lanes while
    this runs, and a second glob would hash a different set from the one it reports on.
    """
    names, where = readme_cited()
    every = sorted(p.stem for p in RESULTS.glob("*.json")) if every is None else list(every)
    cited = census(names)
    everything = census(every)
    unrecorded = [r["result"] for r in everything["rows"] if r["verdict"] == "unrecorded"]
    return {
        "status": (
            "a read-only census of which script wrote each result in data/results, with the results "
            "that record no entry script split by whether the field existed when they were written"
        ),
        "population": population(every),
        "readme_cited": cited,
        "all_results": everything,
        "unrecorded_split": split_by_boundary(unrecorded),
        "stdin_written": stdin_written(everything["rows"]),
        "citations": where,
        "readings_carried_word_for_word": READINGS,
        "nothing_was_written": (
            "this census reads result files and writes none of them: several results' sha256 are "
            "registrations in data/results/manifest_headlines.json and rewriting one would break a "
            "pin. The census reads and reports; what to do about what it finds is not its decision."
        ),
        "itself": (
            f"data/results/{RESULT}.json is a result in data/results like any other, so this run "
            "counts whatever version of it was already on disk -- none at all on the first run -- and "
            "never the bytes it is about to write"
        ),
        "examined_files": [f"{RESULTS.as_posix()}/{n}.json" for n in every],
    }


def manifest_for(every: list[str]) -> dict[str, Any]:
    """The manifest for the census's own result: the examined results hashed one by one."""
    return {
        "sources": [
            {
                "accession": f"data/results/*.json in this checkout ({len(every)} files)",
                "version": "the working tree at this run's git_sha; each file pinned here by sha256",
            },
            {
                "accession": "git history of this repository (genomeos/manifest.py and data/results)",
                "version": "read at this run; the boundary commit is named in unrecorded_split",
            },
        ],
        "inputs": [
            mf.files_entry(
                "data/results/*.json as examined by this census",
                [RESULTS / f"{n}.json" for n in every],
                partition=None,
                role="read for result_manifest.code.argv and result_manifest.code_cleanliness only; "
                "never written",
            )
        ],
        "assembly": "n/a: this census reads result manifests and git history, no genomic sequence",
        "coordinates": "n/a: no genomic interval is read or reported",
        "parameters": {
            "argv_field": ARGV_FIELD,
            "verdicts": sorted(VERDICTS),
            "sides": sorted(SIDES),
            "unplaceable_classes": list(UNPLACEABLE),
            "alphagenome_requests": 0,
            "results_rewritten": 0,
        },
        "exclusions": [
            "README.md is read to find the names it cites and is not an input to any count",
            "a registered name with no file on this disk is reported as `absent`, never skipped",
            "no result's bytes are read for anything but its own manifest, and none is written",
        ],
        "partitions": "n/a: a census of a record, not an evaluation",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--all", action="store_true", help="every data/results/*.json as a second table")
    ap.add_argument("--json", action="store_true")
    ap.add_argument(
        "--write",
        action="store_true",
        help=f"save the whole census as data/results/{RESULT}.json (implies --all)",
    )
    args = ap.parse_args(argv)

    if args.write:
        every = sorted(p.stem for p in RESULTS.glob("*.json"))
        payload = artefact(every)
        path = save_result(RESULT, payload, manifest=manifest_for(every))
        split = payload["unrecorded_split"]
        b = split["boundary"]
        print(f"written: {path}")
        print(
            f"population: {payload['population']['files']} files, "
            f"{payload['population']['tracked_by_git']} tracked, "
            f"{payload['population']['untracked']} untracked"
        )
        print(f"boundary for {b['field']}: {b['sha']} {b['date']}")
        print(f"unrecorded: {split['examined']}")
        for side, n in split["counts"].items():
            print(f"  {side:34s} {n}")
        print(f"  THE FINDING  after the field existed: {split['the_finding']['count']}")
        print(f"  HISTORY      before the field existed: {split['the_history']['count']}")
        print(f"  UNPLACEABLE  {split['the_unplaceable']['count']}")
        return 0

    names, where = readme_cited()
    cited = census(names)
    everything = census(sorted(p.stem for p in RESULTS.glob("*.json"))) if args.all else None
    if args.json:
        print(json.dumps({"readme_cited": cited, "all_results": everything, "citations": where}, indent=1))
        return 0

    print(f"README-cited results examined: {cited['examined']}")
    for verdict, n in sorted(cited["counts"].items()):
        print(f"  {verdict:28s} {n}")
    print(f"\nargv[0] OUTSIDE the repository: {cited['outside_repository_count']}")
    for row in cited["outside_repository"]:
        print(f"  {row['result']}")
        print(f"    argv[0]               {row['argv0']}")
        print(f"    own_code_is_committed {row['claims_own_code_is_committed']}")
        print(f"    the result's bytes   {row['bytes_state']}")
        print(f"    cited in              {'; '.join(where.get(row['result'], []))}")
    if cited["falsely_certified"]:
        print(
            f"\nof those, asserting own_code_is_committed True: "
            f"{len(cited['falsely_certified'])} "
            f"({', '.join(r['result'] for r in cited['falsely_certified'])})"
        )
    if everything is not None:
        print(f"\nevery committed result ({everything['examined']} files), for scale:")
        # The line above is lane-entrypoint's and is left exactly as committed; this one corrects its
        # denominator rather than replacing it, because the glob is the WORKING TREE and a result no
        # commit holds is in it. `population()` counts the two apart in the written result.
        pop = population(sorted(p.stem for p in RESULTS.glob("*.json")))
        print(
            f"  (that glob is the working tree: {pop['tracked_by_git']} of those "
            f"{pop['files']} files are tracked by git and {pop['untracked']} are untracked, "
            "so the denominator is not 'every committed result')"
        )
        for verdict, n in sorted(everything["counts"].items()):
            print(f"  {verdict:28s} {n}")
        print(f"  outside and asserting own_code_is_committed True: {len(everything['falsely_certified'])}")
        print(
            f"  of those, with the result's own bytes COMMITTED: "
            f"{len(everything['falsely_certified_and_committed'])}"
        )
        for row in everything["falsely_certified"]:
            print(f"    {row['result']:28s} bytes {row['bytes_state']:10s} {row['argv0']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
