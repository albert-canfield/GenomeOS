# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The verdict of a check, written down as a file that names the tree it judged.

WHY THIS EXISTS. Four times in one night a check's exit code lied, each time in a different
hand: a guard piped into `tail` exited 0 over a REFUSED and the `&&` after it ran; a
`scripts/check.sh ... | tail` reported 0 over four real failures; a work-board command was
reported as exiting 0 when it exits 2, because the reading was taken through a pipe; and a
backgrounded check returned 0 while its log said `19 failed`. A pipeline's `$?` is the last
command's, a background launch's `$?` is the launcher's, and neither is the checker's.

A fifth fault is worse and no pipe is involved. The check takes nine to twelve minutes and
several sessions commit to this one checkout, so the tree moves underneath a run. One run
reported four failures that existed in NEITHER the tree it started on nor the tree it ended
on: pytest had collected the old test file and executed it against a script a peer committed
nine minutes later. An exit code cannot say WHICH TREE it judged, so even an honest one is
not a verdict about anything in particular.

So the verdict is a file, and the file names its subject.

WHICH TREE HASH IS RECORDED, AND WHY THAT ONE. The hash recorded is the hash of the WORKING
TREE -- every tracked and every untracked-but-not-ignored file as it sits on disk -- and not
the hash of the index and not the hash of HEAD. That is the right hash because it is the only
one that describes what the check actually read. `ruff check .` reads files on disk; `pytest`
COLLECTS FILES ON DISK, including a test file nobody has staged; `bio test` reads programs on
disk. The index describes what someone intends to commit, which may be a subset, a superset or
neither, and HEAD describes what was committed before any of the editing. Recording either of
those would let a run be credited to content it never opened, which is the exact mistake this
file exists to stop. The one happy consequence is that in a clean checkout of a commit -- which
is how scripts/pre-push.sh runs the check -- the working-tree hash IS that commit's tree hash,
so a consumer can compare it directly against the thing it is about to push.

The hash is taken TWICE, before the legs and after them, and both are recorded. Equal hashes
mean the run had one subject. Unequal hashes mean the tree moved mid-run, which is the
nine-minute hazard, and the file then says so instead of offering a verdict about a tree that
no longer existed by the time the counts were printed.

WHAT KIND OF RED, recorded in the file. A verdict that only says "red" cannot say whether the tree is
at fault, and on 2026-10-02 five red files were caused by nothing in the tree: two by a `.json` handed
to ruff as a lint argument, two by `git check-ignore` being unable to answer past a symlinked data
store, and a fifth that refused a push for the same reason. Each read exactly like a real failure
until somebody opened the log. So every file now carries `error_class`, and the list of classes with
it, so a reader a month later needs neither the log nor this docstring: `test` is a failed assertion,
`test_setup` is pytest red with errors and no failures -- tests that could not start -- `tooling` is a
non-test leg, `tree_moved` is the hazard above, `unknown` is a red with nothing to attribute it to.
It is a reading aid and nothing more: `require` still refuses anything that is not green, whatever
class it is, because a red verdict of any class cannot answer for a push.
The four false reds of 2026-10-02 are listed in docs/FALSE-RED-VERDICTS.md, with what the record
holds and what it does not.

HOW THE FILE IS WRITTEN. To a temporary path in the same directory, flushed and fsynced, and
then `os.replace`d into position. A reader therefore sees either the previous complete file or
the new complete file, never half of one, and a run that is killed before `write` leaves the
previous file in place rather than a truncated thing that reads as a pass. A killed run leaves
no verdict for its own tree, and a missing verdict is not a pass.

HOW IT IS READ. `require --tree T` answers for one tree and one tree only. It looks for the
file named after T; a verdict for any other tree is not an answer to the question, however
green it is. On success it prints one token, GENOMEOS_VERDICT_OK, and exits 0; shell callers
are asked to require BOTH, because a shell that takes only the exit code is back in the trap
this file is about.
"""

from __future__ import annotations

import argparse
import contextlib
import json
import os
import platform
import re
import socket
import subprocess
import sys
import tempfile
import time
from pathlib import Path

SCHEMA = "genomeos-check-status/1"
OK_TOKEN = "GENOMEOS_VERDICT_OK"  # printed by `require` on success, and required by shell callers

# `require` exits with this when the answer is no. Distinct from 1 and 2 so a caller can tell a
# refusal from the tool having failed to run at all.
EXIT_REFUSED = 3


def _git(*args: str, cwd: Path, env: dict[str, str] | None = None) -> str:
    full = dict(os.environ)
    if env:
        full.update(env)
    out = subprocess.run(
        ["git", *args],
        cwd=str(cwd),
        env=full,
        check=True,
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


def repo_root(start: Path | None = None) -> Path:
    return Path(_git("rev-parse", "--show-toplevel", cwd=start or Path.cwd()))


def status_dir(root: Path) -> Path:
    """Where verdicts live: one directory shared by the checkout and all of its worktrees.

    `--git-common-dir` and not `--git-dir`, because scripts/pre-push.sh runs the check inside a
    temporary worktree whose per-worktree git dir the pushing shell will never look in, while the
    common dir is the same path for both. Under .git, because a verdict is per-machine runtime
    state like data/jobs: it is about one worktree at one moment and means nothing elsewhere, so
    it must not be committable and must not show up in `git status`.
    """
    override = os.environ.get("GENOMEOS_CHECK_STATUS_DIR")
    if override:
        return Path(override)
    common = _git("rev-parse", "--path-format=absolute", "--git-common-dir", cwd=root)
    return Path(common) / "genomeos-check"


def working_tree_hash(root: Path) -> str:
    """The hash of the working tree: tracked and untracked-not-ignored files, as they are on disk.

    Built in a THROWAWAY index file in a temporary directory, never the repository's own index.
    That matters in a checkout several sessions share: the shared-checkout rule is that nothing
    stages a peer's half-written file, and nothing here does -- GIT_INDEX_FILE points at a path
    under a mkdtemp that is deleted on the way out, so `git add -A` here cannot affect what any
    session is about to commit. The assertion below is there so a later edit cannot quietly drop
    the override and turn this reader into a writer.
    """
    with tempfile.TemporaryDirectory(prefix="genomeos-treehash-") as td:
        index = Path(td) / "index"
        # raised and not asserted: `python -O` strips an assert, and a check that can be compiled
        # away is not a check. These two cost nothing and they are the whole safety of the `add`.
        if index.exists():
            raise RuntimeError("the throwaway index must start absent")
        env = {"GIT_INDEX_FILE": str(index)}
        _git("add", "-A", ".", cwd=root, env=env)
        if not index.exists():
            raise RuntimeError(f"git wrote no index at {index}; GIT_INDEX_FILE did not take effect")
        return _git("write-tree", cwd=root, env=env)


# pytest's own summary line, e.g. "19 failed, 3239 passed, 4 skipped in 412.33s" or "3258 passed in
# 400.11s". Counting words rather than positions, because the order and the set both vary with the
# run and a positional parse would read a green run as unparsable the first time nothing is skipped.
_COUNT_WORDS = (
    "passed",
    "failed",
    "error",
    "errors",
    "skipped",
    "xfailed",
    "xpassed",
    "deselected",
    "warnings",
    "warning",
)
_COUNT_RE = re.compile(r"(\d+)\s+(" + "|".join(_COUNT_WORDS) + r")\b")
# the summary is the last line that both reports a count and reports a duration
_SUMMARY_RE = re.compile(r"\bin\s+[\d.]+\s*s(econds)?\b")


def parse_pytest_counts(text: str) -> dict[str, int] | None:
    """Counts from pytest's summary line, or None when there is no summary line to read.

    None is not zero and must never be reported as a pass: it is what a killed run, a crashed
    interpreter or a changed pytest looks like, and the caller's job is to say so.
    """
    summary = None
    for line in text.splitlines():
        stripped = line.strip().strip("=").strip()
        if _COUNT_RE.search(stripped) and _SUMMARY_RE.search(stripped):
            summary = stripped
    if summary is None:
        return None
    counts = {"passed": 0, "failed": 0, "errors": 0, "skipped": 0}
    for number, word in _COUNT_RE.findall(summary):
        key = "errors" if word in ("error", "errors") else word
        if key in ("warnings", "warning"):
            continue
        counts[key] = counts.get(key, 0) + int(number)
    return counts


def decide(exit_code: int, tree_begin: str, tree_end: str, counts: dict[str, int] | None) -> tuple[str, str]:
    """The verdict and the one-line reason for it. Green requires every reason to hold at once."""
    if tree_begin != tree_end:
        return "tree_moved", (
            f"the working tree changed during the run: began {tree_begin}, ended {tree_end}; "
            "the counts belong to neither tree in full"
        )
    if exit_code != 0:
        return "red", f"the check exited {exit_code}"
    if counts is None:
        return "unreadable", (
            "the check exited 0 but no pytest summary line could be read, so the counts are unknown"
        )
    if counts.get("failed", 0) or counts.get("errors", 0):
        return "red", (
            f"the check exited 0 but the counts say {counts.get('failed', 0)} failed and "
            f"{counts.get('errors', 0)} errors"
        )
    if counts.get("passed", 0) <= 0:
        return "unreadable", "the check exited 0 but no test passed, so there is nothing to call green"
    return "green", f"{counts['passed']} passed, 0 failed, 0 errors, {counts.get('skipped', 0)} skipped"


#: WHAT KIND OF THING MADE A RUN RED, recorded so a tooling-caused red is distinguishable from a test
#: red IN THE RECORD ITSELF. Until 2026-10-02 it was not: five red status files that day were caused by
#: nothing in the tree, and each read exactly like a real failure until somebody opened the log.
#:
#: Four of them were argument and environment faults -- two where a `.json` was handed to ruff as a
#: lint argument (96 "errors" in a result file, the run dead at the ruff-check leg), and two where
#: `git check-ignore` could not answer past a symlinked data store and 20 tests ERRORED at setup. The
#: fifth refused a push. The distinction below is the one that separates them at a glance:
#:
#:   test        a test FAILED. Something the suite asserts about the code did not hold.
#:   test_setup  pytest was red with ERRORS and NO failures: tests could not START. Fixtures,
#:               collection, markers, missing data, a tool that could not answer -- nothing was
#:               asserted wrongly, so this is where an environment or tooling fault lands.
#:   tooling     a leg that is not a test refused or failed: lint, format, shellcheck, startup.
#:   tree_moved  the working tree changed mid-run, so the counts belong to no one tree.
#:   unknown     red or unreadable with nothing in the file to attribute it to. Never a pass.
#:
#: A run with BOTH failures and errors is `test`, deliberately: a real failure exists and must be
#: fixed whatever else is wrong, and a class that hid it behind the errors would be the softer reading.
ERROR_CLASSES = ("", "test", "test_setup", "tooling", "tree_moved", "unknown")

#: Legs of scripts/check.sh that are not tests. `biolang` is absent on purpose: `bio test` is the
#: BioLang programs testing themselves, so its red is a test red and is classified as one.
NON_TEST_LEGS = ("startup", "ruff-check", "ruff-format", "shellcheck", "report")


def classify(verdict: str, failed_leg: str, counts: dict[str, int] | None) -> str:
    """Which of ERROR_CLASSES this verdict is. Empty for green; never empty for anything else."""
    if verdict == "green":
        return ""
    if verdict == "tree_moved":
        return "tree_moved"
    if failed_leg in NON_TEST_LEGS and failed_leg:
        return "tooling"
    if counts:
        if counts.get("failed", 0):
            return "test"
        if counts.get("errors", 0):
            return "test_setup"
    return "unknown"


def _atomic_write_json(path: Path, payload: dict) -> None:
    """Either the old file or the whole new one: write beside the target, fsync, then rename.

    Same directory so the rename is within one filesystem and therefore atomic; fsync before the
    rename so a crash cannot leave a present-but-empty file where a verdict is expected.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp-status-", suffix=".json", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        os.unlink(tmp)
        raise


def status_path(root: Path, tree: str) -> Path:
    return status_dir(root) / f"status-{tree}.json"


def latest_path(root: Path) -> Path:
    return status_dir(root) / "latest.json"


def write_status(
    root: Path,
    *,
    exit_code: int,
    tree_begin: str,
    scope: str,
    scope_files: list[str],
    started_at: float,
    pytest_log: Path | None = None,
    failed_leg: str = "",
    note: str = "",
) -> dict:
    tree_end = working_tree_hash(root)
    log_text = ""
    if pytest_log and pytest_log.exists():
        log_text = pytest_log.read_text(errors="replace")
    counts = parse_pytest_counts(log_text) if log_text else None
    verdict, reason = decide(exit_code, tree_begin, tree_end, counts)
    payload = {
        "schema": SCHEMA,
        "verdict": verdict,
        "reason": reason,
        "error_class": classify(verdict, failed_leg, counts),
        "error_classes": list(ERROR_CLASSES),
        "exit_code": exit_code,
        "tree_begin": tree_begin,
        "tree_end": tree_end,
        "tree_changed": tree_begin != tree_end,
        "counts": counts,
        "counts_read": counts is not None,
        "scope": scope,
        "scope_files": scope_files,
        "failed_leg": failed_leg,
        "note": note,
        "pytest_log": str(pytest_log) if pytest_log else "",
        "started_at": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(started_at)),
        "finished_at": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime()),
        "seconds": round(time.time() - started_at, 1),
        "pid": os.getpid(),
        "ppid": os.getppid(),
        "host": socket.gethostname(),
        "platform": platform.platform(),
        "repo": str(root),
    }
    # the per-tree file first, so that the moment `latest.json` advertises a tree the file it
    # names is already complete and a reader following the pointer cannot miss
    _atomic_write_json(status_path(root, tree_begin), payload)
    _atomic_write_json(latest_path(root), payload)
    _prune(status_dir(root), keep={status_path(root, tree_begin), latest_path(root)})
    return payload


PRUNE_AFTER_SECONDS = 14 * 24 * 3600


def _prune(directory: Path, *, keep: set[Path], now: float | None = None) -> list[Path]:
    """Drop verdicts and test logs older than a fortnight. Never the two just written.

    A verdict for a tree nobody will push again is dead weight, but deleting on a schedule is also
    how a consumer could be made to find nothing; so this only ever removes what is already far too
    old to answer for anything being pushed, and a failure to delete is not allowed to fail a run.
    """
    now = time.time() if now is None else now
    removed = []
    for path in list(directory.glob("status-*.json")) + list(directory.glob("pytest-*.log")):
        if path in keep:
            continue
        try:
            if now - path.stat().st_mtime > PRUNE_AFTER_SECONDS:
                path.unlink()
                removed.append(path)
        except OSError:
            continue
    return removed


def tree_diff_summary(root: Path, a: str, b: str, limit: int = 12) -> list[str]:
    """The paths by which two trees differ, so a hash mismatch is a message and not a mystery.

    A refusal that only prints two forty-character hashes tells the reader nothing they can act on,
    and a guard whose reason is thrown away has already cost this project a commit. Best effort: if
    either tree is not in this object store the list is empty and the refusal still stands.
    """
    try:
        out = subprocess.run(
            ["git", "diff", "--name-status", a, b],
            cwd=str(root),
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (subprocess.CalledProcessError, OSError):
        return []
    lines = [line for line in out.splitlines() if line]
    if len(lines) > limit:
        return lines[:limit] + [f"... and {len(lines) - limit} more paths"]
    return lines


def read_status(root: Path, tree: str) -> dict | None:
    path = status_path(root, tree)
    try:
        text = path.read_text()
    except FileNotFoundError:
        return None
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return None


def require(root: Path, tree: str, *, scope: str | None, out=None, err=None) -> int:
    """Answer, for ONE tree, whether a complete green verdict exists. Anything else refuses."""
    # resolved on the call and not in the signature: a default bound at import time keeps a
    # reference to the original stream, so a caller that has redirected sys.stdout -- a test
    # harness, a log capture -- would not see what this printed
    out = sys.stdout if out is None else out
    err = sys.stderr if err is None else err
    tree = tree.strip()
    if not re.fullmatch(r"[0-9a-f]{40}", tree):
        print(f"verdict: REFUSED: {tree!r} is not a tree hash", file=err)
        return EXIT_REFUSED
    status = read_status(root, tree)
    if status is None:
        print(f"verdict: REFUSED: no check verdict exists for tree {tree}.", file=err)
        print(f"  looked in {status_path(root, tree)}", file=err)
        other = None
        with contextlib.suppress(OSError, json.JSONDecodeError):
            other = json.loads(latest_path(root).read_text())
        if other:
            judged = other.get("tree_begin") or ""
            print(
                f"  the most recent verdict in this checkout judged tree {judged} "
                f"({other.get('verdict')}) -- a DIFFERENT tree, so it does not answer for this one.",
                file=err,
            )
            differences = tree_diff_summary(root, judged, tree)
            if differences:
                print(f"  the two trees differ in {len(differences)} path(s):", file=err)
                for line in differences:
                    print(f"    {line}", file=err)
        print("  A missing verdict is not a pass. Run scripts/check.sh on this tree.", file=err)
        return EXIT_REFUSED
    if status.get("schema") != SCHEMA:
        print(f"verdict: REFUSED: {status_path(root, tree)} is not a {SCHEMA} file", file=err)
        return EXIT_REFUSED
    for field in ("tree_begin", "tree_end"):
        if status.get(field) != tree:
            print(
                f"verdict: REFUSED: the verdict file for {tree} says {field}="
                f"{status.get(field)}; it did not judge this tree throughout.",
                file=err,
            )
            return EXIT_REFUSED
    if scope and status.get("scope") != scope:
        print(
            f"verdict: REFUSED: the verdict for {tree} has scope {status.get('scope')!r}, "
            f"not {scope!r}. A run that linted only some files is not a verdict on the project.",
            file=err,
        )
        return EXIT_REFUSED
    if status.get("verdict") != "green":
        print(
            f"verdict: REFUSED: the verdict for {tree} is {status.get('verdict')!r}: {status.get('reason')}",
            file=err,
        )
        # said here as well as in the file, because the reader of a refused push is the person who
        # most needs to know whether the red is about the code at all
        if status.get("error_class"):
            print(
                f"  error_class {status['error_class']!r}, failed_leg "
                f"{status.get('failed_leg') or '(none)'!r}"
                + (" -- NOT a test failure" if status["error_class"] != "test" else ""),
                file=err,
            )
        if status.get("note"):
            print(f"  note: {status['note']}", file=err)
        if status.get("pytest_log"):
            print(f"  the check's test output is in {status['pytest_log']}", file=err)
        return EXIT_REFUSED
    counts = status.get("counts") or {}
    print(
        f"{OK_TOKEN} {tree} passed={counts.get('passed')} failed={counts.get('failed')} "
        f"errors={counts.get('errors')} skipped={counts.get('skipped')} "
        f"scope={status.get('scope')} finished={status.get('finished_at')}",
        file=out,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="genomeos.verdict", description=__doc__.splitlines()[0])
    parser.add_argument("--repo", default=None, help="repository root (default: the one containing cwd)")
    sub = parser.add_subparsers(dest="cmd", required=True)

    sub.add_parser("tree", help="print the hash of the working tree as it is on disk")
    sub.add_parser("status-dir", help="print the directory verdict files are written to")

    p_write = sub.add_parser("write", help="write the status file for a finished run")
    p_write.add_argument("--exit-code", type=int, required=True)
    p_write.add_argument("--tree-begin", required=True, help="the hash printed by `tree` before the legs ran")
    p_write.add_argument("--scope", default="project", choices=["project", "files"])
    p_write.add_argument("--scope-file", action="append", default=[])
    p_write.add_argument("--started-at", type=float, required=True, help="unix seconds the run began")
    p_write.add_argument("--pytest-log", default="")
    p_write.add_argument("--failed-leg", default="")
    p_write.add_argument("--note", default="")

    p_req = sub.add_parser("require", help="refuse unless a green verdict exists for exactly this tree")
    p_req.add_argument("--tree", required=True)
    p_req.add_argument("--scope", default=None, choices=["project", "files"])

    p_show = sub.add_parser("show", help="print the verdict file for a tree, or the most recent one")
    p_show.add_argument("--tree", default=None)

    args = parser.parse_args(argv)
    root = Path(args.repo) if args.repo else repo_root()

    if args.cmd == "tree":
        print(working_tree_hash(root))
        return 0
    if args.cmd == "status-dir":
        print(status_dir(root))
        return 0
    if args.cmd == "write":
        payload = write_status(
            root,
            exit_code=args.exit_code,
            tree_begin=args.tree_begin.strip(),
            scope=args.scope,
            scope_files=args.scope_file,
            started_at=args.started_at,
            pytest_log=Path(args.pytest_log) if args.pytest_log else None,
            failed_leg=args.failed_leg,
            note=args.note,
        )
        print(
            f"check: verdict {payload['verdict']} for tree {payload['tree_begin']} "
            f"({payload['reason']}) -> {status_path(root, payload['tree_begin'])}"
        )
        return 0
    if args.cmd == "require":
        return require(root, args.tree, scope=args.scope)
    if args.cmd == "show":
        path = latest_path(root) if args.tree is None else status_path(root, args.tree)
        try:
            sys.stdout.write(path.read_text())
        except FileNotFoundError:
            print(f"verdict: no file at {path}", file=sys.stderr)
            return EXIT_REFUSED
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
