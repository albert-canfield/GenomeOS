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

For every path that differs between HEAD and the index, it takes the lines the commit would remove,
SUBTRACTS the lines the same commit adds back, and compares what is left against the lines each recent
commit ADDED to that same path. An honest edit of your own text hits nothing. A stale-base staging hits
the peer commit exactly, and is named with its sha, its subject and the lines themselves.

The subtraction is what makes the question "is this content gone", not "did a line move". Comparison is
on the stripped text, so git reports a line put inside an `if`/`else` as a removal plus an addition of
the same content; before 2026-10-02 only the removal was read and that honest shape was refused as a
revert of the line it preserves. Content the diff puts back, at any indentation and anywhere in the
file, is not undone. The subtraction is exact, so a line that merely resembles a held one -- a variable
renamed, a comment appended -- is still a removal and is still refused.

It is deliberately a stale-base detector and not a merge policeman. Deleting your own lines, or lines
older than the window, is ordinary work and passes. Rewriting a line that a peer added minutes ago is
sometimes right too -- so a finding is a refusal you can override with `--force` once you have looked,
not a lock. What it removes is the case where nobody looked at all.

Since 2026-10-02 it asks a second, unrelated question of the same index, of every staged
`data/results/*.json` (`unreproducible_results`):

    could any commit ever reproduce this number?

A result stamps the revision it was written at, so when the code that wrote it was uncommitted that
sha names a tree WITHOUT that code and nothing reproduces the file. That is the one fault here which
makes a number permanently unverifiable rather than merely awkward. Three shapes are refused: the
writer's own code uncommitted, a peer's uncommitted file on the import closure that computed the
number, and no `result_manifest` at all unless the name is on the committed legacy allowlist that
`genomeos.results.save_result` reads. Writing such a result locally stays possible; committing it
does not.

A staged deletion of a file that exists in HEAD is always reported: an accidental one is unrecoverable
by the author who lost it, and a deliberate one costs a flag.

**Do not pipe it, and do not silence it.** Two ways to defeat this tool, both used by the session that
wrote it, on consecutive days:

- `check_staged.py | tail -2` reports the refusal and exits 0, because a pipeline's status is the last
  command's, so the `&&` that was meant to stop the commit runs it instead (2026-09-16: a commit went
  through a refusal that had been printed and read);
- `check_staged.py >/dev/null` keeps the exit code and throws away the reason, so an honest `&&` chain
  stops with nothing on screen and the refusal looks like a silent failure of something else
  (2026-09-17: a commit was written twice because its author could not see why the first died).

The first of those happened a third time on 2026-09-17, to the same session, in a compound command
that staged, checked and committed in one line to save a round trip. The warning above was read that
morning and did not survive being in a hurry. **If you must pipe it, run the shell with
`set -o pipefail`** — a refusal then exits 2 through the pipe instead of 0, which is one word between
a guard that works and a guard that watches. Better still, give it its own line and read it.

Run it bare. If output must be captured, send it to a file and read the file.

Exit 0 clean, 2 on a finding, 1 on a usage error.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

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
    if not staged_paths(index):
        # an empty commit is what happens when a push that looked like it failed had actually landed
        # and the work was committed again on top: 49edf1e on 2026-09-17, an exact duplicate of
        # 4fc7b71 with the same message and no diff. Harmless, untidy, and invisible until someone
        # reads the log.
        return ["the staged tree is identical to HEAD: this commit would be empty"]
    if index:
        # every read below is against this index; the per-commit diffs ignore it, so one setting serves
        os.environ["GIT_INDEX_FILE"] = index
    for status, path in staged_paths(index):
        if status == "D":
            problems.append(f"{path}: staged as DELETED, and HEAD has it. Intended? --force says so.")
            continue
        if status == "A":
            continue  # a file HEAD does not have cannot revert anything
        put_back, removed = _diff_lines("HEAD", "--cached", "--", path)
        # Git reports a moved or re-indented line as a removal plus an addition of the same text, and
        # `meaningful` already compares stripped, so content this same diff puts back -- at any
        # indentation, anywhere in the file -- is not undone and has nothing to protect. Without this
        # subtraction, putting a peer's line inside an `if`/`else` was refused as a revert of it, and a
        # false refusal is pressure to reach for `--force`, which is the flag that can really lose work.
        # The subtraction is exact on the stripped text: a line that merely resembles a held one --
        # a renamed variable, a comment appended -- is a different string and is still a removal.
        removed -= put_back
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


# `from genomeos.a.b import c` and `import genomeos.a.b`: the dotted module a staged file needs.
_IMPORT = re.compile(r"^\s*(?:from\s+(genomeos(?:\.\w+)*)\s+import\b|import\s+(genomeos(?:\.\w+)*))", re.M)


def missing_imports(index: str | None) -> list[str]:
    """Staged Python that imports a genomeos module the committed tree will not contain.

    38087a1 on 2026-09-28 committed a test importing `genomeos.certainty` while the module was still
    untracked. In the shared tree everything passed, because the file was there; in the clean worktree
    the pre-push hook builds, ruff sorted the unknown module as third-party (I001) and the test could
    not import, and every lane's push stopped behind it. The committed tree is the index, so the
    question is whether the index holds the module.
    """
    env = dict(os.environ)
    if index:
        env["GIT_INDEX_FILE"] = index
    listed = subprocess.run(["git", "ls-files", "--cached"], capture_output=True, text=True, env=env)
    tracked = set(listed.stdout.splitlines())
    problems = []
    for status, path in staged_paths(index):
        if status == "D" or not path.endswith(".py"):
            continue
        shown = subprocess.run(["git", "show", f":{path}"], capture_output=True, text=True, env=env)
        for match in _IMPORT.finditer(shown.stdout):
            module = match.group(1) or match.group(2)
            base = module.replace(".", "/")
            if f"{base}.py" in tracked or f"{base}/__init__.py" in tracked:
                continue
            problems.append(
                f"{path}: imports {module}, but neither {base}.py nor {base}/__init__.py is in the "
                "commit. Stage the module with it: a clean checkout cannot import it."
            )
    return problems


#: A result lives directly under this prefix -- `data/results/<name>.json`, and `<name>` is what the
#: legacy allowlist lists.
RESULTS_PREFIX = "data/results/"
#: The two fields of `result_manifest.code_cleanliness` that decide whether a commit can reproduce the
#: numbers. Named here only so a refusal can print the field a reader must go and look at; the values
#: are written by the one shared `genomeos.manifest.code_cleanliness`.
OWN_COMMITTED = "own_code_is_committed"
FOREIGN_ON_PATH = "foreign_uncommitted_code_on_the_counting_path"
#: The marker lines tests/test_check_staged.py cuts between to show, by removal, what gets through
#: without these conditions. Keep the loop inside them.
CONDITIONS_BEGIN = "    # --- the three result conditions ---"
CONDITIONS_END = "    # --- end of the three result conditions ---"


def _repo_root() -> Path:
    """The tree being judged."""
    top = _git("rev-parse", "--show-toplevel").strip()
    return Path(top) if top else Path.cwd()


def contract() -> tuple[str, Path, frozenset[str]]:
    """The manifest key, where the legacy allowlist lives, and the names on it -- all three read from
    the modules `save_result` itself uses.

    `genomeos.results.legacy_names` is the function that decides who is exempt when a result is
    written and `genomeos.results.LEGACY_ALLOWLIST` is where it looks, so the gate and the writer
    cannot disagree: a second copy of 955 names in this file would drift from the committed list the
    first time anybody touched either. The list is read from the tree under test, the code to read it
    from this script's own repository, because a fixture repository holds results and an allowlist but
    no package. A list that cannot be read is empty, so enforcement fails closed and a name that is
    not on it must carry the contract.

    The import is deliberately inside the function: `scripts/` is this script's sys.path entry, not
    the repository root, and a module-level import of `genomeos` would make a usage error of every
    run from somewhere else -- including the runs that have nothing to do with results.
    """
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from genomeos.manifest import KEY
    from genomeos.results import LEGACY_ALLOWLIST, legacy_names

    return KEY, LEGACY_ALLOWLIST, legacy_names(_repo_root() / LEGACY_ALLOWLIST)


def unreproducible_results(index: str | None) -> list[str]:
    """Staged results that no commit reproduces, in the three shapes seen on 2026-10-02.

    A result stamps the revision it was written at. When the code that wrote it was uncommitted, that
    sha names a tree WITHOUT that code and no commit anywhere reproduces the file. Of every fault here
    this is the only one that makes a number permanently unverifiable rather than merely awkward,
    because everything else can at least be re-derived. It happened twice that night:
    data/results/response_map_increment3_count.json at 3360c49 stamped 727c936 while naming its own
    script under `own_uncommitted_code`, and an earlier result did the same at 42bbc4e.

    So, of every staged `data/results/*.json`:

    (a) `code_cleanliness.own_code_is_committed` false -- the writer's own code is not in the tree the
        result points at;
    (b) `code_cleanliness.foreign_uncommitted_code_on_the_counting_path` non-empty -- another lane's
        uncommitted file sat on the import closure that computed the number. Seven results carried
        this that night, while one lane held two modules that sit on every result's counting path;
    (c) no `result_manifest` at all, unless the name is on the legacy allowlist.

    (c) is not a check that was missing. `tests/test_manifest_enforced.py` caught
    data/results/astroreg_request_plan.json at 1f880b2 and did not fail to refuse it -- it refused it
    fifteen minutes into a push, in a clean worktree, by which time the repair had landed at 5d14d65
    and the verdict was about a tree that no longer existed. The question is the same one; asking it of
    the index is what makes the answer act on anything.

    Writing such a result locally stays possible, deliberately: trial runs need it, and that night a
    lane wrote a trial result and then deleted it so nothing produced from uncommitted code could be
    taken for the real one. The refusal belongs at the commit that publishes the file, not at the write.

    `own_code_is_committed` is compared against `False` and not read for falsity, because the shared
    function writes a bool, and an absent or null field is a manifest that never answered the question
    -- which `save_result` refuses on the key set for every non-legacy name before anything reaches
    here. For the same reason a manifest with no `code_cleanliness` block passes (a) and (b) silently:
    that gap is the writer's, it is already closed there, and guessing at it here would refuse the
    legacy results whose manifests predate the block.
    """
    problems: list[str] = []
    rows = [
        path
        for status, path in staged_paths(index)
        if status != "D"
        and path.startswith(RESULTS_PREFIX)
        and path.endswith(".json")
        and "/" not in path[len(RESULTS_PREFIX) :]
    ]
    if not rows:
        return problems
    env = dict(os.environ)
    if index:
        env["GIT_INDEX_FILE"] = index
    # --- the three result conditions ---
    # inside the markers with the conditions it serves, so the copy the counterfactual test makes of
    # this file needs no package at all to run
    key, listed_at, legacy = contract()
    for path in rows:
        # the staged blob, not the working copy: what this commit would publish
        shown = subprocess.run(["git", "show", f":{path}"], capture_output=True, text=True, env=env)
        try:
            body = json.loads(shown.stdout)
        except ValueError as exc:
            problems.append(f"{path}: not readable as JSON ({exc}), so its `{key}` cannot be read")
            continue
        name = path[len(RESULTS_PREFIX) : -len(".json")]
        given = body.get(key) if isinstance(body, dict) else None
        if not isinstance(given, dict):
            if name in legacy:
                continue  # written before the contract, and the committed list says so
            problems.append(
                f"{path}: no `{key}`, and `{name}` is not on the legacy allowlist ({listed_at}), so "
                "nothing in this commit says which code and which inputs produced these numbers. "
                "astroreg_request_plan.json went in this way at 1f880b2."
            )
            continue
        clean = given.get("code_cleanliness")
        if not isinstance(clean, dict):
            continue
        if clean.get(OWN_COMMITTED) is False:
            own = clean.get("own_uncommitted_code") or []
            named = ", ".join(str(p) for p in own) or "not named"
            problems.append(
                f"{path}: `{key}.code_cleanliness.{OWN_COMMITTED}` is false, so the revision it stamps "
                f"({str(clean.get('git_sha'))[:7]}) names a tree WITHOUT the code that wrote it and no "
                f"commit reproduces this file. Uncommitted own code: {named}."
            )
        foreign = clean.get(FOREIGN_ON_PATH)
        if isinstance(foreign, list) and foreign:
            problems.append(
                f"{path}: `{key}.code_cleanliness.{FOREIGN_ON_PATH}` is not empty, so another lane's "
                "uncommitted code sat on the import closure that computed these numbers and the "
                f"stamped tree does not contain it: {', '.join(str(p) for p in foreign)}."
            )
    # --- end of the three result conditions ---
    return problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--index", default=os.environ.get("GIT_INDEX_FILE"), help="the index to check")
    ap.add_argument("--since", default=DEFAULT_SINCE, help="how far back a peer commit counts")
    ap.add_argument("--force", action="store_true", help="report but do not refuse")
    args = ap.parse_args(argv)

    problems = check(args.index, args.since)
    unimportable = missing_imports(args.index)
    unreproducible = unreproducible_results(args.index)
    if unimportable:
        # not a revert, so not the message below; and not something --force should wave through,
        # because the fix is always to stage the module
        print("REFUSED: this commit imports a module it does not contain\n")
        for problem in unimportable:
            print(f"    {problem}\n")
        return 2
    if unreproducible:
        # a separate heading because it is a separate question: not "is a peer's work going out" but
        # "could anyone ever get this number back from a commit"
        print("REFUSED: this commit adds a result no commit reproduces\n")
        for problem in unreproducible:
            print(f"    {problem}\n")
        print(
            "    Commit the code, re-run the result, and stage the two together, so the revision the\n"
            "    result stamps names a tree that contains the code that wrote it. Writing such a result\n"
            "    locally stays possible on purpose -- trial runs need it -- and committing it does not.\n"
            "    If you have looked and the result is right as it stands, re-run with --force.\n"
        )
    if not problems:
        if not unreproducible:
            print("staged tree takes nothing back out")
            return 0
        return 0 if args.force else 2
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
