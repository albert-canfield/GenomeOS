# SPDX-License-Identifier: AGPL-3.0-or-later
"""Refuse a commit that would undo work someone else just pushed.

    GIT_INDEX_FILE=/tmp/idx-$$ uv run python scripts/check_staged.py
    GIT_INDEX_FILE=/tmp/idx-$$ uv run python scripts/check_staged.py --since 3.days

Several sessions write this checkout at once. A lane that stages a whole file stages the copy it read,
and if HEAD moved between the read and the commit, the commit silently reverts whatever arrived in
between. On 2026-09-15 that cost a complete feature: a panel commit staged `docs/BIOLANG-GRAMMAR.md`,
`genomeos/lang/`, `genomeos/ir/`, `genomeos/runtime/` and `tests/test_share.py` from a stale read and
took out the whole `share:` construct, its grammar row, its documentation and its 155 lines of tests.
Nobody noticed until the lane that wrote it came back. It was the third such accident that day.

`scripts/stage_section.py` prevents this for a markdown section. This catches the rest, structurally,
by asking one question of the staged tree rather than trusting anyone's care:

    does this commit delete lines that a recent commit added?

For every path that differs between HEAD and the index, it takes the lines the commit would remove and
compares them against the lines each recent commit ADDED to that same path. An honest edit of your own
text hits nothing. A stale-base staging hits the peer commit exactly, and is named with its sha, its
subject and the lines themselves.

It is deliberately a stale-base detector and not a merge policeman. Deleting your own lines, or lines
older than the window, is ordinary work and passes. Rewriting a line that a peer added minutes ago is
sometimes right too -- so a finding is a refusal you can override with `--force` once you have looked,
not a lock. What it removes is the case where nobody looked at all.

A staged deletion of a file that exists in HEAD is always reported: an accidental one is unrecoverable
by the author who lost it, and a deliberate one costs a flag.

**Do not pipe it.** `check_staged.py | tail -2` reports the refusal and then exits 0, because a shell
pipeline's status is the last command's, so the `&&` that was meant to stop the commit runs it instead.
This session did exactly that on 2026-09-16 and committed through a refusal it had printed and read.
Run it bare, or redirect to a file and grep that.

Exit 0 clean, 2 on a finding, 1 on a usage error.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys

# Lines that carry no authorship: blank, a lone brace, a bare fence. Counting them would make every
# reformatting look like a revert.
NOISE = {"", "{", "}", "(", ")", "[", "]", "```", '"""', "'''", "*/", "#", "//", "--", "|"}
MIN_LEN = 3
DEFAULT_SINCE = "2.days"
MAX_COMMITS = 60
SAMPLES = 4


def _git(*args: str) -> str:
    """Run git and return stdout; a failure is empty, since every caller treats it as 'nothing there'."""
    done = subprocess.run(["git", *args], capture_output=True, text=True)
    return done.stdout if done.returncode == 0 else ""


def meaningful(lines: list[str]) -> set[str]:
    """The lines worth attributing: stripped, not noise, long enough to be somebody's writing."""
    out = set()
    for line in lines:
        s = line.strip()
        if s in NOISE or len(s) < MIN_LEN:
            continue
        out.add(s)
    return out


def _diff_lines(*args: str) -> tuple[set[str], set[str]]:
    """The meaningful lines added and removed by a diff."""
    text = _git("diff", "--no-color", "--unified=0", *args)
    added, removed = [], []
    for line in text.splitlines():
        if line.startswith("+++") or line.startswith("---"):
            continue
        if line.startswith("+"):
            added.append(line[1:])
        elif line.startswith("-"):
            removed.append(line[1:])
    return meaningful(added), meaningful(removed)


def staged_paths(index: str | None) -> list[tuple[str, str]]:
    """(status, path) for everything the index changes against HEAD."""
    env = dict(os.environ)
    if index:
        env["GIT_INDEX_FILE"] = index
    done = subprocess.run(
        ["git", "diff-index", "--cached", "--name-status", "HEAD"],
        capture_output=True,
        text=True,
        env=env,
    )
    if done.returncode != 0:
        raise SystemExit(f"cannot read the index: {done.stderr.strip()}")
    rows = []
    for line in done.stdout.splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            rows.append((parts[0][0], parts[-1]))
    return rows


def recent_commits(path: str, since: str) -> list[tuple[str, str, str]]:
    """(sha, author, subject) of the commits that touched this path inside the window, newest first."""
    text = _git(
        "log",
        f"--since={since}",
        f"--max-count={MAX_COMMITS}",
        "--format=%H%x00%an%x00%s",
        "HEAD",
        "--",
        path,
    )
    out = []
    for line in text.splitlines():
        sha, _, rest = line.partition("\0")
        author, _, subject = rest.partition("\0")
        if sha:
            out.append((sha, author, subject))
    return out


def reverted_by(path: str, removed: set[str], since: str) -> list[dict[str, object]]:
    """Recent commits whose added lines this staging would take back out."""
    findings = []
    for sha, author, subject in recent_commits(path, since):
        added, _ = _diff_lines(f"{sha}^!", "--", path)
        overlap = removed & added
        if overlap:
            findings.append(
                {
                    "sha": sha[:7],
                    "author": author,
                    "subject": subject,
                    "lines": sorted(overlap, key=len, reverse=True),
                }
            )
    return findings


def check(index: str | None, since: str) -> list[str]:
    """Every reason to refuse this commit, in the words its author needs to act on."""
    problems = []
    if index:
        # every read below is against this index; the per-commit diffs ignore it, so one setting serves
        os.environ["GIT_INDEX_FILE"] = index
    for status, path in staged_paths(index):
        if status == "D":
            problems.append(f"{path}: staged as DELETED, and HEAD has it. Intended? --force says so.")
            continue
        if status == "A":
            continue  # a file HEAD does not have cannot revert anything
        _, removed = _diff_lines("HEAD", "--cached", "--", path)
        if not removed:
            continue
        for found in reverted_by(path, removed, since):
            lines = found["lines"]
            assert isinstance(lines, list)
            shown = "\n".join(f"        - {line[:100]}" for line in lines[:SAMPLES])
            more = f"\n        ... and {len(lines) - SAMPLES} more" if len(lines) > SAMPLES else ""
            problems.append(
                f"{path}: removes {len(lines)} lines that {found['sha']} added "
                f"({found['author']}: {found['subject'][:60]})\n{shown}{more}"
            )
    return problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--index", default=os.environ.get("GIT_INDEX_FILE"), help="the index to check")
    ap.add_argument("--since", default=DEFAULT_SINCE, help="how far back a peer commit counts")
    ap.add_argument("--force", action="store_true", help="report but do not refuse")
    args = ap.parse_args(argv)

    problems = check(args.index, args.since)
    if not problems:
        print("staged tree takes nothing back out")
        return 0
    print("REFUSED: this commit would undo work that is already in HEAD\n")
    for problem in problems:
        print(f"    {problem}\n")
    print(
        "    Stage the HEAD copy plus your own hunk instead (scripts/stage_section.py for a markdown\n"
        "    section), reading the base and writing the index in one call so HEAD cannot move between.\n"
        "    If you have looked and the removal is right, re-run with --force."
    )
    return 0 if args.force else 2


if __name__ == "__main__":
    sys.exit(main())
