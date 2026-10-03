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

Since 2026-10-03 it asks a third question, of the work board this time rather than of git
(`staged_held_by_a_live_lane`):

    is another lane working on this file right now?

The two questions above are both about lines that are already COMMITTED. Nothing here used to say
anything about including a peer's UNCOMMITTED additions, and on 2026-10-03 that cost a red HEAD: a
lane staged a whole test file carrying another lane's in-flight tests and committed it under its own
message -- guard clean, exit 0, no `--force` -- and the tests went in without the imports they needed,
so HEAD failed on `NameError`. The information was already in the checkout at the moment it happened:
the holder had declared the file on the work board, and `.claude/hooks/shared_checkout_guard.py`
already says so when you EDIT such a file. It was the commit that never asked. So a staged path that a
live lane other than you declares in its `files` is refused, named with its holder, its task and its
age. There is deliberately no `--force` for it: the holder is a person to ask, not a flag.

Alongside the refusal it prints a NOTICE, which refuses nothing, for a SHARED file whose whole
working copy this commit stages (`whole_copy_notices`). A file counts as shared two ways: it is on
the list below, or a live lane holds it. The list alone was the hook's whole test on 2026-10-03 and
it has eight entries, none of them a test file, which is the second reason nothing fired that day.
Being held is the other way, so the list stops being the only way a file can be shared -- and the
list stays, because a file can be shared without anyone holding it today. The notice goes quiet by
itself when you do the right thing: it fires only while the index for that path is byte-identical to
the working copy, so staging the HEAD copy plus your own hunk silences it.

The subtraction above has one exception, also 2026-10-03 (`live_put_back`). Because it compares the
stripped text of whole lines and takes content back from ANYWHERE in the file, a removed line pasted
verbatim into a comment or into a docstring stopped counting as removed: the guard read the quotation
as the line. Quoting superseded PROSE as history is this project's practice and must keep working, so
the distinction is made on the added text's structure and not on how it looks. In a Python file the
staged blob is tokenized, and a line that no token other than a comment or a string literal touches
is not code the program runs; a removal that was live code in HEAD is cancelled only by live code in
the index. A moved comment, a reworded docstring and every non-Python file keep the old behaviour
exactly. What this cannot tell apart is written at `live_code_text` and `live_put_back`.

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
import io
import json
import os
import re
import subprocess
import sys
import tokenize
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


#: Token types that carry no live code. A line that none of the OTHER token types touches is a
#: comment line, or a line inside a string literal -- a docstring's interior, a triple-quoted block.
#: Read by name with a default, because the f-string token types only exist from python 3.12.
_QUOTED_TOKENS = frozenset(
    getattr(tokenize, name, -1)
    for name in (
        "COMMENT",
        "STRING",
        "NL",
        "NEWLINE",
        "INDENT",
        "DEDENT",
        "ENCODING",
        "ENDMARKER",
        "FSTRING_START",
        "FSTRING_MIDDLE",
        "FSTRING_END",
    )
)
#: The marker lines tests/test_check_staged.py cuts between to show, by removal, that a removed line
#: quoted into a docstring goes through without this block.
QUOTED_BEGIN = "        # --- the quoted-cancellation condition ---"
QUOTED_END = "        # --- end of the quoted-cancellation condition ---"


def live_code_text(source: str) -> set[str] | None:
    """The stripped text of every line of this Python source that holds live code.

    How the two are told apart, stated plainly: the source is tokenized, and a line is live code if
    any token that is NOT a comment, a string literal or layout begins or ends on it. A line nothing
    else touches is a comment, or the inside of a string literal, and the program does not run it.
    This is the structure of the text and not a guess at its shape -- no regex for a leading `#`,
    which would miss a docstring's interior, and no indentation rule, which would miss everything.

    `None` means the source could not be tokenized. The caller then keeps the old behaviour, because
    a refusal resting on a file this script could not read would be a guess, and this script's own
    history says a false refusal is pressure towards `--force`.

    What it cannot tell apart:

    - a line inside a string that IS the program's output -- generated code in a template, a
      `textwrap.dedent` block, a SQL or shell literal. Moving a line into such a literal is a real
      move, and this reads it as a quotation and refuses it. That refusal goes to the coordinator;
      it does not go into a comment;
    - a line of an f-string that carries an interpolation, which brings NAME tokens onto the line
      and so reads as live code even though its literal text is not;
    - which quotation was deliberate. It does not try: a removed line of code is reported whether it
      was quoted on purpose or pasted in by accident.
    """
    lines = source.splitlines()
    live: set[int] = set()
    try:
        for token in tokenize.generate_tokens(io.StringIO(source).readline):
            if token.type in _QUOTED_TOKENS:
                continue
            for number in range(token.start[0], token.end[0] + 1):
                live.add(number)
    except (tokenize.TokenError, IndentationError, SyntaxError):
        return None
    return {lines[n - 1].strip() for n in live if 1 <= n <= len(lines)}


def live_put_back(path: str, index: str | None, put_back: set[str], removed: set[str]) -> set[str]:
    """`put_back`, with a quotation no longer able to cancel the removal of a line of code.

    A text that was live code in HEAD has to be live code in the index to count as put back. Every
    other text -- a moved comment, a reworded docstring, a prose sentence quoted as history -- is
    cancelled from anywhere in the file, exactly as before.

    What it cannot tell apart, beyond `live_code_text`'s list:

    - anything that is not a `.py` file. A removed line of shell, markdown, JSON or JavaScript
      pasted into a comment still cancels its own removal. Python is the language this script can
      read without guessing, and the three findings of 2026-10-03 were in Python and in markdown;
    - a text that exists in HEAD both as live code and as prose. It is judged as code, because
      being code somewhere in HEAD is what the question is about;
    - a quotation pasted at the line's OWN original indentation. git then matches it against the
      line it replaces and reports no removal at all, so nothing reaches this subtraction and there
      is nothing to refuse. What counts as removed is the diff's judgement, not this script's, and
      that case is outside both. A quotation of code inside prose is almost always re-indented,
      which is the case this closes.
    """
    if not path.endswith(".py"):
        return put_back
    env = dict(os.environ)
    if index:
        env["GIT_INDEX_FILE"] = index
    staged = subprocess.run(["git", "show", f":{path}"], capture_output=True, text=True, env=env)
    head = subprocess.run(["git", "show", f"HEAD:{path}"], capture_output=True, text=True, env=env)
    if staged.returncode != 0 or head.returncode != 0:
        return put_back
    was_code = live_code_text(head.stdout)
    is_code = live_code_text(staged.stdout)
    if was_code is None or is_code is None:
        return put_back
    code_removals = removed & was_code
    return (put_back - code_removals) | (put_back & code_removals & is_code)


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
        # --- the quoted-cancellation condition ---
        # a line that comes back only inside a comment or a string literal has not come back
        put_back = live_put_back(path, index, put_back, removed)
        # --- end of the quoted-cancellation condition ---
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


#: An entry in one of these states is somebody's current work. The states are written by
#: `genomeos.work.board`, the reader the board's own `genomeos work` command uses, which demotes an
#: entry to `stale` after six hours and to `abandoned` after two days. The freshness window is
#: therefore the board's own and not a second copy of it here, exactly as the legacy allowlist below
#: is read from the module that writes it rather than copied.
LIVE_STATES = frozenset({"working", "waiting"})
#: The marker lines tests/test_check_staged.py cuts between to show, by removal, that a peer's
#: in-flight file goes in silently without this block.
HELD_BEGIN = "    # --- the live-hold condition ---"
HELD_END = "    # --- end of the live-hold condition ---"


#: Files that several sessions write at once, so staging one whole is how a peer's hunk travels.
#: `.claude/hooks/shared_checkout_guard.py` has its own copy of this list on its edit path; the two
#: are deliberately independent, because the hook is untracked, needs Albert's approval to change
#: and cannot be imported from a tracked tool. Each fails open on its own, and neither is the only
#: way a file counts as shared here -- the work board is the second way.
SHARED = frozenset(
    {
        "genomeos/cli.py",
        "genomeos/jobs.py",
        "genomeos/web/server.py",
        "genomeos/web/static/index.html",
        "README.md",
        "CONTRIBUTING.md",
        "docs/ROADMAP.md",
        "docs/PROGRESS.md",
        "docs/ATTRIBUTION.md",
        "docs/LESSONS.md",
    }
)
#: The marker lines tests/test_check_staged.py cuts between to show, by removal, that a held file
#: gets no notice at all without this block.
WHOLE_COPY_BEGIN = "    # --- the whole-copy condition ---"
WHOLE_COPY_END = "    # --- end of the whole-copy condition ---"


def live_board(root: Path) -> list[dict[str, object]]:
    """Every work-board entry that is somebody's current work.

    `genomeos.work.board` is imported rather than re-implemented: it is the function the `genomeos
    work` command and the Progress tab read with, so the gate cannot disagree with the board about
    what counts as live. The import is inside the function, and only reached when the board
    directory exists, for the same reason `contract()` imports inside itself -- `scripts/` is this
    script's sys.path entry, not the repository root, and most runs of this script have nothing to
    do with the board.
    """
    if not (root / "data" / "work").is_dir():
        return []
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from genomeos.work import board

    return [e for e in board(root) if e.get("state") in LIVE_STATES]


def declared(entry: dict[str, object]) -> tuple[set[str], set[str]]:
    """The (files, directories) an entry claims.

    `genomeos work --files` splits on commas, so a space-separated list given as one shell argument
    is stored as one element with spaces in it -- there are such entries on the board today. Both
    separators are read here, because an entry that names four files in one string is claiming four
    files whichever way it was typed.
    """
    files: set[str] = set()
    dirs: set[str] = set()
    for raw in entry.get("files") or []:  # type: ignore[union-attr]
        for item in str(raw).replace(",", " ").split():
            (dirs if item.endswith("/") else files).add(item)
    return files, dirs


def is_me(who: str, lane: str, message_name: str) -> bool:
    """Whether a board entry's owner is the lane that is committing.

    Two ways to say so, and the second needs no new habit: `--lane` (which scripts/commit_own.sh
    passes from `-L` or `GENOMEOS_LANE`), or the committing message file's name, which commit_own.sh
    already REQUIRES to carry the lane as a literal substring. So the same test that decides the
    message file is yours decides the board entry is yours, and a lane that follows the naming rule
    is recognised with nothing to configure.
    """
    if not who:
        return False
    return who == lane or (bool(message_name) and who in message_name)


def staged_held_by_a_live_lane(index: str | None, lane: str, message_name: str) -> list[str]:
    """Staged paths another live lane is working on, which this commit would take whole.

    Only an EXACT file claim refuses. A directory claim (`tests/`, `scripts/`) raises the whole-copy
    notice instead and not a refusal, deliberately: several entries on this board hold `tests/` or
    `scripts/` for a day at a time, and refusing on those would stop every other lane from
    committing any test or any script. A named file is a claim on that file; a directory is a claim
    on a region, and the region is where people work in parallel. That is the hole this leaves: a
    holder who declared only the directory is not protected by the refusal, only by the notice.
    """
    problems: list[str] = []
    root = _repo_root()
    board = live_board(root)
    if not board:
        return problems
    # --- the live-hold condition ---
    for status, path in staged_paths(index):
        if status == "D":
            continue
        for entry in board:
            who = str(entry.get("who") or "")
            if is_me(who, lane, message_name):
                continue
            files, _dirs = declared(entry)
            if path not in files:
                continue
            age = entry.get("age")
            minutes = round(float(age) / 60) if isinstance(age, (int, float)) else "?"
            task = str(entry.get("task") or "no task recorded")
            problems.append(f"{path}: held on the work board by {who}, {minutes} min ago: {task[:140]}")
    # --- end of the live-hold condition ---
    return problems


def whole_copy_notices(index: str | None) -> list[str]:
    """Shared files whose whole working copy this commit stages, with every uncommitted hunk in it.

    Not a refusal. Most whole-copy stagings of a shared file are honest -- all the hunks are yours --
    and refusing them would be refusing ordinary work, so this says what is going in and leaves the
    decision where it belongs. The condition is three things at once: the path is shared (listed, or
    held by a live lane), the index for it is byte-identical to the working copy, and HEAD differs
    from it. The middle one is what makes the notice stop by itself: stage the HEAD copy plus your own
    hunk and the index is no longer the working copy, so nothing is said.
    """
    notices: list[str] = []
    root = _repo_root()
    board = live_board(root)
    if index:
        os.environ["GIT_INDEX_FILE"] = index
    # --- the whole-copy condition ---
    for status, path in staged_paths(index):
        if status not in {"M", "A"}:
            continue
        because = []
        if path in SHARED:
            because.append("a shared file")
        for entry in board:
            who = str(entry.get("who") or "")
            age = entry.get("age")
            minutes = round(float(age) / 60) if isinstance(age, (int, float)) else "?"
            task = str(entry.get("task") or "no task recorded")[:90]
            files, dirs = declared(entry)
            inside = next((d for d in sorted(dirs) if path.startswith(d)), None)
            if path in files:
                because.append(f"held on the work board by {who} ({minutes} min ago: {task})")
            elif inside:
                because.append(
                    f"inside {inside}, held on the work board by {who} ({minutes} min ago: {task})"
                )
        if not because:
            continue
        # the working copy, not HEAD: an index that differs from the worktree is a partial staging,
        # which is the thing this notice asks for, so it is not asked for again
        if _git("diff", "--name-only", "--", path).strip():
            continue
        hunks = _git("diff", "--no-color", "--unified=0", "HEAD", "--cached", "--", path).count("\n@@")
        if not hunks:
            continue
        notices.append(
            f"{path} is {'; '.join(because)}, and this stages its whole working copy: all "
            f"{hunks} uncommitted hunk(s), a peer's included. Read `git diff HEAD -- {path}`. If any "
            "hunk is not yours, stage the HEAD copy plus your own hunk instead "
            "(scripts/stage_section.py for a markdown section)."
        )
    # --- end of the whole-copy condition ---
    return notices


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
    ap.add_argument("--lane", default=os.environ.get("GENOMEOS_LANE", ""), help="the lane committing")
    ap.add_argument(
        "--message-name", default="", help="the commit message file's name, which carries the lane"
    )
    args = ap.parse_args(argv)

    notices = whole_copy_notices(args.index)
    held = staged_held_by_a_live_lane(args.index, args.lane, args.message_name)
    problems = check(args.index, args.since)
    for notice in notices:
        # first, and above any refusal, because it is the sentence that explains the refusal when
        # there is one and the only sentence there is when there is not
        print(f"NOTICE: {notice}\n")
    unimportable = missing_imports(args.index)
    unreproducible = unreproducible_results(args.index)
    if held:
        # Not a revert and not a result, so neither heading below; and no --force, because the thing
        # to do about a file somebody is working on is to talk to them. The lanes never force in this
        # project -- the coordinator does, after asking the holder -- so a flag here would only move
        # the decision back to the lane that is already not looking.
        print("REFUSED: this commit takes a file another lane is working on right now\n")
        for problem in held:
            print(f"    {problem}\n")
        print(
            "    Their lines are not committed yet, so staging the file whole puts their work in\n"
            "    progress into your commit under your message. On 2026-10-03 that landed a test file\n"
            "    without the imports it needed and HEAD went red on NameError.\n"
            "    Two routes, and no flag: hand the file to the coordinator, who checks with the holder\n"
            "    and commits it as one change (the holder releasing it with `genomeos work update\n"
            "    --files` is the same route from the other end); or stage only your own hunk -- the HEAD\n"
            "    copy plus your lines, scripts/stage_section.py for a markdown section -- and leave\n"
            "    theirs in the worktree.\n"
        )
        return 2
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
