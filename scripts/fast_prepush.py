#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fast pre-push: the TARGETED tests, run TWICE -- once in this tree, once with no data stores.

    uv run python scripts/fast_prepush.py             # both legs, on the window described below
    uv run python scripts/fast_prepush.py --plan      # print the window and the selection, run nothing
    uv run python scripts/fast_prepush.py --only store-free

    # what scripts/pre-push.sh calls, from a clean worktree of the commit being pushed
    uv run python scripts/fast_prepush.py --decide --remote-ref refs/heads/dev
    uv run python scripts/fast_prepush.py --gate --sha SHA --tree TREE --remote-ref refs/heads/dev
    uv run python scripts/fast_prepush.py --require-fast --tree TREE

WHAT IT IS FOR, and it is one defect class and not a second full gate. On 2026-10-03 three CI runs went
red -- 3 failures, then 6 -- on failures the LOCAL FULL SUITE COULD NOT SEE: every one was a test that
assumes this machine's git-ignored data (`data/knowledge`, `data/cache`, `data/reference`). This
checkout has those stores, CI is a fresh checkout and does not, so the suite is green here and red
there by construction. A lesson was written about the class and the very next commit reproduced it. The
full gate costs 17 to 26 minutes and still misses it, because the environment it runs in is the one
that hides it. The leg below removes the environment instead of adding tests, and it catches exactly
this class, in about a minute, at the hand that wrote the code.

WIRED INTO `scripts/pre-push.sh` ON 2026-10-10 (ROADMAP section 5, E1), as the supervisor ruled on
2026-10-03 (ROADMAP-HISTORY, the row beside `1b2e768`; CONTRIBUTING step 4). Until then this paragraph
said the script did not replace the gate, and it did not: nothing called it. Now a push to `dev` takes
THE FAST GATE -- `ruff check` and `ruff format --check` over the project, then ALWAYS_RUN and the window
below in BOTH legs -- and its verdict is a status file of scope "fast" (`fast-<tree>.json`, beside and
never in place of the project verdicts). THE FULL SUITE STAYS WHERE IT DECIDES SOMETHING, and `--decide`
below sends the push to `scripts/check.sh` instead whenever one of these holds:

  * the daily CI is RED or UNREAD: `<common>/genomeos-check/daily-ci`, written by
    `scripts/last_green_ci.py --daily`, is absent, says anything but `success`, or was read more than
    DAILY_READ_MAX_SECONDS ago -- the full gate until a reading is green;
  * the push is not to `refs/heads/dev`: `main`, and the `promote-<sha>` branch that
    `scripts/promote_main.sh --prepare` pushes, so the sha promoted to main has its local full suite;
  * `GENOMEOS_FULL_GATE=1` is in the environment: the switch for a sha a PAID run will be sent from,
    and for anyone who wants the full suite on a dev push.

THE TRADE, stated as the ruling states it: a red in a test file that did not change, and is not in
ALWAYS_RUN, does not block a dev push. The daily CI and the promotion's full suite are where it is found.

THE WINDOW IS SINCE THE LAST GREEN CI, NOT SINCE ORIGIN -- and that is the whole design, not a detail.
A file changed in a commit that already reached origin while its CI run was red, cancelled or unread is
INVISIBLE to a since-origin window. That is precisely how the six got out: they were pushed, CI went
red, the pushes continued, and from then on `origin/dev..HEAD` no longer contained the files that were
breaking. Measured on this checkout, 2026-10-03: the since-origin window names 4 test files and the
since-last-green window names 113. On a normal day with green CI behind you the two are the same set
and the wider rule costs nothing.

NO NETWORK CALL HERE. The last green sha is read from a local file, `<git common dir>/genomeos-check/
last-green-ci`, which `scripts/last_green_ci.py` writes; that script's docstring holds the reason and
the linked-worktree measurement behind the path. A hook that asks GitHub has a rate limit and an auth
failure in the way of every push.

THE THREE FALLBACK CASES, AND TWO OF THEM ARE NOT THE SAME THING:

  MISSING  -- a fresh clone, or the first push after this lands. Nothing is known and nothing is wrong.
              Fall back to the since-origin window plus ALWAYS_RUN, print a notice, DO NOT BLOCK.
  STALE    -- older than 48 h. USE ITS SHA ANYWAY. A stale record still names a sha CI once went green
              on, so the window from it is WIDER than the since-origin window, strictly safer, and
              available with no network. Narrowing on staleness would be "we know less, so we check
              less", which is backwards; and a notice that recurs on every push is read once and
              invisible by the fifth time, so the staleness warning must not be the only consequence.
              `tests/test_fast_prepush.py` asserts the FILE COUNT for this case and not the message,
              because a hook that printed the warning and narrowed anyway would pass a message test.
  FRESH    -- the normal path: the window starts at the recorded sha.

Two further cases are real and are handled with the missing one, each under its own name so a reader is
never told "missing" about a file that is there: UNKNOWN-SHA (a record naming a commit this checkout
does not have, e.g. after a reclone) and NO-BASE (no remote-tracking ref either, so there is nothing to
diff against at all -- then the selection is ALWAYS_RUN alone and it still does not block).

THE STORE-FREE LEG HAS TWO ASSERTIONS AND BOTH ARE REQUIRED BEFORE A PASS IS TRUSTED:

  1. THE THREE STORES ARE GENUINELY ABSENT. A leg that silently found them present would pass
     everything and report nothing wrong -- a guard keyed to a condition it does not verify is not a
     guard, it is a green light wired to nothing. This is the assertion that makes the leg's green mean
     something, so its failure is a REFUSAL and not a pass.
  2. THE WORKTREE'S SOURCE IS WHAT GOT IMPORTED. `PYTHONPATH` precedence is not a promise: an editable
     install of this project pointing at the main checkout, a stale `.pth`, or a `__pycache__` would
     each let the leg test the tree it was trying to exclude. So `genomeos.__file__` is read inside the
     leg's interpreter and required to be under the worktree.

The worktree is removed through `scripts/remove_worktree.py`, NEVER `rm -rf`: that script exists
because `rm -rf` on a path that may hold a clone of a store is the one operation that can destroy the
machine's only copy of 11.2 GB that is not in git. `--force` is passed because a pytest run leaves
`__pycache__` and `.pytest_cache` behind and git refuses a worktree with untracked files; the stores are
never linked into this worktree, so there is nothing of the machine's inside it to force away.

TEMPORARY DIRECTORIES COME FROM `tempfile.mkdtemp`. Not `mktemp -t`: GNU mktemp refuses a `-t` template
with fewer than three trailing X's, which once made `scripts/push_own.sh` unrunnable on Linux while
green on every Mac. Python's mkdtemp has no template to get wrong.

CHANGED HELPER MODULES. Seven of the files under `tests/` are imported rather than collected. They are
handled by what they are, and the reasoning is printed in the plan rather than left implicit:
`tests/always_run.py` carries four tests of its own and pytest collects an explicitly named file
whatever `python_files` says, so it is named directly; a helper with no tests of its own is replaced by
the test modules that IMPORT it, because that is where a change to it can fail; and `tests/conftest.py`
is loaded by every pytest run there is, so a change to it is exercised by both legs already and gets a
notice instead of an expansion. Nothing here is an exemption: a changed file is either run or its
importers are.
"""

from __future__ import annotations

import argparse
import ast
import json
import os
import re
import socket
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from datetime import UTC, datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import last_green_ci  # noqa: E402  (a sibling script, imported for the ONE path both of them use)

ROOT = Path(__file__).resolve().parent.parent

#: the three git-ignored machine-local stores whose absence the store-free leg both creates and asserts
STORES: tuple[str, ...] = ("data/knowledge", "data/cache", "data/reference")
#: older than this and the record is stale. It is still USED; see the docstring.
STALE_SECONDS = 48 * 3600
#: the remote-tracking ref the fallback window is measured from
REMOTE_REF = "origin/dev"
#: the pathspec the window is taken over
TESTS_PATHSPEC = "tests/*.py"
#: the one ref whose push MAY take the fast gate; every other ref the hook checks takes the full one
DEV_REF = "refs/heads/dev"
#: set (to anything but 0) and a push takes the full gate whatever else holds: the paid-path switch
FULL_GATE_ENV = "GENOMEOS_FULL_GATE"
#: a daily reading of CI older than this is UNREAD. "Read first thing each morning": this morning's
#: reading is current, and yesterday morning's lapses a day later, with two hours' slack for a morning
#: that starts later than the one before. One number, here, so the hook and its tests cannot disagree.
DAILY_READ_MAX_SECONDS = 26 * 3600
#: the fast verdict's own schema and scope. NOT genomeos/verdict.py's schema, so that even a fast file
#: copied to a project verdict's name is refused there by its schema check before its scope is read.
FAST_SCHEMA = "genomeos-fast-gate/1"
FAST_SCOPE = "fast"
#: printed by `--require-fast` on success, and required by the shell caller as well as the exit code
FAST_OK_TOKEN = "GENOMEOS_FAST_OK"
#: `--require-fast` exits with this when the answer is no, as genomeos/verdict.py's `require` does
EXIT_REFUSED = 3
#: fast verdicts older than this are pruned when a new one is written, as the project verdicts are
FAST_PRUNE_AFTER_SECONDS = 14 * 24 * 3600


def git(root: Path, *args: str, check: bool = True) -> str:
    r = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True)
    if check and r.returncode != 0:
        raise RuntimeError(f"git {' '.join(args)} exited {r.returncode}: {r.stderr.strip()}")
    return r.stdout.strip()


def always_run(root: Path = ROOT) -> tuple[str, ...]:
    """ALWAYS_RUN, read from `tests/always_run.py` itself so this file holds no second copy of the list.

    Read from `MEMBERS` and not from `ALWAYS_RUN`, because `ALWAYS_RUN` is built by a `tuple(...)`
    comprehension over `MEMBERS` and `ast.literal_eval` cannot evaluate a call -- the first version of
    this function tried and raised. `MEMBERS` is a literal tuple of dict literals, so it evaluates.

    Read by AST and not by importing: this runs before any leg, in whatever interpreter the hook has,
    and the list must not depend on `tests/` being importable from here. A file that has changed shape
    enough that neither name reads falls back to RUNNING the module, which prints the list -- so a
    rename there costs a subprocess and never costs the list.
    """
    source = (root / "tests" / "always_run.py").read_text()
    tree = ast.parse(source)
    for node in tree.body:
        targets = (
            [node.target]
            if isinstance(node, ast.AnnAssign)
            else (node.targets if isinstance(node, ast.Assign) else [])
        )
        if not any(getattr(t, "id", "") == "MEMBERS" for t in targets):
            continue
        try:
            members = ast.literal_eval(node.value)
        except ValueError:
            break
        paths = tuple(m["path"] for m in members if isinstance(m, dict) and m.get("path"))
        if paths:
            return paths
    printed = subprocess.run(
        [sys.executable, str(root / "tests" / "always_run.py")], cwd=root, capture_output=True, text=True
    )
    paths = tuple(printed.stdout.split())
    if not paths:
        raise RuntimeError(
            f"tests/always_run.py yielded no ALWAYS_RUN members, by MEMBERS or by running it "
            f"(exit {printed.returncode}: {printed.stderr.strip()[-200:]})"
        )
    return paths


# ----------------------------------------------------------------------------- the window


@dataclass
class Window:
    case: str
    base: str | None
    notices: list[str] = field(default_factory=list)
    record: dict | None = None
    block: bool = False
    age_seconds: float | None = None

    def as_dict(self) -> dict:
        return {
            "case": self.case,
            "base": self.base,
            "age_seconds": self.age_seconds,
            "notices": self.notices,
            "record": self.record,
            "block": self.block,
        }


def read_last_green(root: Path) -> tuple[dict | None, Path, str | None]:
    """(the record, where it was looked for, why there is none). Never raises on a bad file."""
    path = last_green_ci.last_green_path(root)
    if not path.exists():
        return None, path, "no such file"
    try:
        record = json.loads(path.read_text())
    except (OSError, ValueError) as e:
        return None, path, f"unreadable ({e.__class__.__name__}: {e})"
    if not isinstance(record, dict) or not record.get("sha"):
        return None, path, "readable but names no sha"
    return record, path, None


def since_origin_base(root: Path, remote_ref: str = REMOTE_REF) -> str | None:
    """The fallback base: the merge base of the remote-tracking ref and HEAD, or None if there is none.

    The merge base and not the ref itself, so a diverged local branch names the commits that are
    actually new here rather than every commit either side has.
    """
    if subprocess.run(
        ["git", "rev-parse", "--verify", "-q", remote_ref], cwd=root, capture_output=True
    ).returncode:
        return None
    base = subprocess.run(["git", "merge-base", remote_ref, "HEAD"], cwd=root, capture_output=True, text=True)
    return base.stdout.strip() or None


def resolve_window(
    root: Path = ROOT,
    now: datetime | None = None,
    remote_ref: str = REMOTE_REF,
    stale_seconds: int = STALE_SECONDS,
) -> Window:
    """Which commit the window starts at, and under which of the cases. Never blocks on its own."""
    record, path, why = read_last_green(root)

    def fallback(case: str, notice: str) -> Window:
        base = since_origin_base(root, remote_ref)
        if base is None:
            return Window(
                case="no-base",
                base=None,
                notices=[
                    notice,
                    f"fast-prepush: and {remote_ref} does not exist here either, so there is no commit "
                    f"to measure a window from. The selection is ALWAYS_RUN alone. Not blocking: this "
                    f"is a fact about the checkout, not about the code being pushed.",
                ],
                record=record,
            )
        return Window(
            case=case,
            base=base,
            notices=[
                notice,
                f"fast-prepush: falling back to the since-{remote_ref} window "
                f"({_short(root, base)}..HEAD) plus ALWAYS_RUN. Not blocking.",
            ],
            record=record,
        )

    if record is None:
        return fallback(
            "missing",
            f"fast-prepush: NOTICE: no last-green-CI record ({why}) at {path}. Nothing is known about "
            f"which sha CI last passed, and nothing is wrong -- a fresh clone looks exactly like this. "
            f"Write one with `uv run python scripts/last_green_ci.py --from-cache`.",
        )

    sha = str(record["sha"])
    if subprocess.run(
        ["git", "rev-parse", "--verify", "-q", f"{sha}^{{commit}}"], cwd=root, capture_output=True
    ).returncode:
        return fallback(
            "unknown-sha",
            f"fast-prepush: NOTICE: the last-green-CI record at {path} names {sha[:12]}, which is not a "
            f"commit in this checkout, so no window can be measured from it. Fetch it, or rewrite the "
            f"record.",
        )

    age = _age_seconds(record, now or datetime.now(UTC))
    if age is not None and age > stale_seconds:
        return Window(
            case="stale",
            base=sha,
            age_seconds=age,
            record=record,
            notices=[
                f"fast-prepush: WARNING: the last-green-CI record is {age / 3600:.1f} h old "
                f"(over {stale_seconds / 3600:.0f} h) -- CI's green sha has probably moved. USING ITS "
                f"SHA ANYWAY and WIDENING the window: a stale record still names a sha CI once went "
                f"green on, so the window from it contains the since-{remote_ref} window and then "
                f"some. Narrowing here would check LESS for knowing less. Refresh it with "
                f"`uv run python scripts/last_green_ci.py --from-cache`.",
            ],
        )

    return Window(
        case="fresh",
        base=sha,
        age_seconds=age,
        record=record,
        notices=[
            f"fast-prepush: window from the last green CI sha {record.get('short') or sha[:7]}"
            + (f" (run {record['run_id']})" if record.get("run_id") else "")
            + (f", recorded {age / 3600:.1f} h ago" if age is not None else ""),
        ],
    )


def _age_seconds(record: dict, now: datetime) -> float | None:
    """How old the RECORD is, by when it was written, falling back to the run's own date.

    When it was written, because that is what "the file is stale" means: the question is how long ago
    somebody last read CI, not how long ago the run happened. A record written now about a run from
    last week is fresh information about a stale repository state, and the window it gives is correct.
    """
    for key in ("written", "date"):
        raw = record.get(key)
        if not raw:
            continue
        try:
            when = datetime.fromisoformat(str(raw).replace("Z", "+00:00"))
        except ValueError:
            continue
        if when.tzinfo is None:
            when = when.replace(tzinfo=UTC)
        return (now - when).total_seconds()
    return None


def _short(root: Path, sha: str) -> str:
    return git(root, "rev-parse", "--short", sha, check=False) or sha[:7]


def changed_test_files(root: Path, base: str | None) -> list[str]:
    """The test files the window covers. An empty window is an empty list, never everything."""
    if base is None:
        return []
    out = git(root, "diff", "--name-only", base, "HEAD", "--", TESTS_PATHSPEC)
    return sorted(p for p in out.splitlines() if p)


# ----------------------------------------------------------------------------- the selection


def has_tests(path: Path) -> bool:
    """Whether pytest would collect at least one test from this file when NAMED EXPLICITLY.

    Named explicitly matters: `python_files` does not apply to a path given on the command line, which
    is why `tests/always_run.py` -- four tests, a name pytest would never match -- is collectible here
    and is not collected by a bare `pytest tests/`. Measured with `--collect-only`: 4 tests collected.
    """
    try:
        tree = ast.parse(path.read_text())
    except (OSError, SyntaxError):
        return False
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)) and node.name.startswith("test_"):
            return True
    return False


def importers_of(root: Path, stem: str) -> list[str]:
    """The `tests/test_*.py` files that import the helper module `stem`, by reading their import lines."""
    pattern = re.compile(
        rf"^\s*(?:from\s+(?:tests\.)?{re.escape(stem)}\s+import|import\s+(?:tests\.)?{re.escape(stem)}\b)",
        re.M,
    )
    out = []
    for p in sorted((root / "tests").glob("test_*.py")):
        try:
            if pattern.search(p.read_text()):
                out.append(str(p.relative_to(root)))
        except OSError:
            continue
    return out


@dataclass
class Selection:
    files: list[str]
    notices: list[str] = field(default_factory=list)
    changed: list[str] = field(default_factory=list)
    helpers: dict[str, list[str]] = field(default_factory=dict)

    def as_dict(self) -> dict:
        return {
            "files": self.files,
            "count": len(self.files),
            "changed_in_window": self.changed,
            "changed_count": len(self.changed),
            "helpers": self.helpers,
            "notices": self.notices,
        }


def select(root: Path, changed: list[str]) -> Selection:
    """The targeted set: the changed test files, the importers of changed helpers, and ALWAYS_RUN."""
    files: set[str] = set()
    notices: list[str] = []
    helpers: dict[str, list[str]] = {}
    for rel in changed:
        path = root / rel
        if not path.exists():
            notices.append(f"fast-prepush: {rel} changed in the window and is not in the tree now; skipped")
            continue
        if has_tests(path):
            files.add(rel)
            continue
        stem = Path(rel).stem
        if stem == "conftest":
            helpers[rel] = []
            notices.append(
                f"fast-prepush: {rel} changed. It carries no tests and cannot be collected, and it is "
                f"LOADED BY EVERY pytest run, so both legs already exercise it. Nothing added for it."
            )
            continue
        pulled = importers_of(root, stem)
        helpers[rel] = pulled
        files.update(pulled)
        notices.append(
            f"fast-prepush: {rel} changed. It carries no tests of its own, so the {len(pulled)} test "
            f"module(s) that import it are in the set instead: that is where a change to it fails."
        )
    for member in always_run(root):
        files.add(member)
    return Selection(files=sorted(files), notices=notices, changed=list(changed), helpers=helpers)


# ----------------------------------------------------------------------------- the legs


#: every child process is THIS interpreter, so no leg ever resolves an environment of its own
PYTHON = sys.executable


def leg_env(root: Path, tree: Path | None) -> dict[str, str]:
    """The environment a leg runs pytest in: `tree` first on the import path, and never a sync.

    `uv` IS NOT IN THE INNER LOOP AT ALL, and that is deliberate. `scripts/pre-push.sh` runs its
    clean-worktree check through `uv run` with `UV_PROJECT_ENVIRONMENT` pointed back at this
    checkout's `.venv`, because the thing it launches is `scripts/check.sh`, a shell script that
    calls `uv` itself. Here the launcher is already inside the right interpreter, so the legs are run
    with `sys.executable` and there is no environment to point anything at, nothing to resolve and
    nothing to get wrong on another machine. That matters for this script above all others: a hook
    whose job is to catch CI-only failures must not itself behave differently in CI, and
    "UV_PROJECT_ENVIRONMENT names a venv that is not where this machine put it" is exactly such a
    difference. `UV_NO_SYNC` is still set for any grandchild that does reach uv, because `uv sync`
    mutates the one environment every session on this machine shares and has removed an extra
    mid-suite.
    """
    env = dict(os.environ)
    env["UV_NO_SYNC"] = "1"
    env.pop("UV_PROJECT_ENVIRONMENT", None)
    if tree is not None:
        env["PYTHONPATH"] = str(tree)
    else:
        env.pop("PYTHONPATH", None)
    return env


def assert_stores_absent(tree: Path) -> dict:
    """REQUIRED before a store-free pass is trusted: the three stores are genuinely not there.

    A leg that found them present would pass everything and report nothing wrong, so this is a
    refusal and never a pass. `exists()` follows symlinks on purpose: a store reached through a link
    is present as far as every test is concerned, and a link is exactly how the stores get into the
    other worktree in this repository.
    """
    present = [s for s in STORES if (tree / s).exists()]
    return {
        "checked": list(STORES),
        "present": present,
        "ok": not present,
        "why": (
            f"REFUSED: {', '.join(present)} exist under {tree}, so this leg is NOT store-free and its "
            f"green would mean nothing. A guard keyed to a condition it does not verify is not a guard."
        )
        if present
        else f"none of {', '.join(STORES)} exists under {tree}",
    }


def assert_import_path(root: Path, tree: Path) -> dict:
    """REQUIRED: the worktree's own source is what an import of `genomeos` resolves to in this leg.

    MEASURED RATHER THAN REASONED, and the measurement corrected the reason. The claim this was
    written for was "PYTHONPATH precedence puts the worktree first". Probed on 2026-10-03 against a
    real store-free worktree, DROPPING PYTHONPATH ENTIRELY CHANGES NOTHING: the worktree still wins,
    because the leg runs with its cwd inside the worktree and that is what goes first on `sys.path`.
    So PYTHONPATH is belt and braces here, not the mechanism, and saying otherwise would have
    credited this guard with work the cwd was doing.

    What the assertion does catch is real and was also measured: against a tree that carries no
    source of its own, `import genomeos` resolves to `/.../GenomeOS/genomeos/__init__.py` -- THE MAIN
    CHECKOUT, which is reachable from this venv -- and this refuses. That is the live hazard: the leg
    silently testing the very tree it exists to exclude. `tests/test_fast_prepush.py` keeps that
    refusal running, so the assertion is known to fire and not merely known to pass.
    """
    probe = "import json,genomeos,sys;print(json.dumps({'file': genomeos.__file__, 'prefix': sys.path[:3]}))"
    r = subprocess.run(
        [PYTHON, "-c", probe],
        cwd=tree,
        env=leg_env(root, tree),
        capture_output=True,
        text=True,
    )
    if r.returncode != 0:
        return {
            "ok": False,
            "why": f"REFUSED: the import probe exited {r.returncode}: {r.stderr.strip()[-400:]}",
        }
    try:
        got = json.loads(r.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return {"ok": False, "why": f"REFUSED: the import probe printed nothing readable: {r.stdout!r}"}
    imported = Path(got["file"]).resolve()
    inside = str(imported).startswith(str(tree.resolve()) + os.sep)
    return {
        "ok": inside,
        "imported": str(imported),
        "worktree": str(tree.resolve()),
        "sys_path_head": got.get("prefix"),
        "why": (
            f"genomeos imported from {imported}"
            if inside
            else f"REFUSED: `import genomeos` resolved to {imported}, which is NOT under {tree.resolve()}. "
            f"This leg would have tested the tree it is meant to exclude."
        ),
    }


def run_pytest(root: Path, tree: Path | None, files: list[str], label: str) -> dict:
    """One pytest run. Exit 5 (nothing collected) is RED: files were named and nothing ran."""
    if not files:
        return {"leg": label, "ok": True, "skipped": "nothing to run", "exit": None, "seconds": 0.0}
    cwd = tree or root
    started = time.monotonic()
    r = subprocess.run(
        [PYTHON, "-m", "pytest", "-q", "-p", "no:cacheprovider", *files],
        cwd=cwd,
        env=leg_env(root, tree),
        capture_output=True,
        text=True,
    )
    seconds = time.monotonic() - started
    tail = (r.stdout + r.stderr).strip().splitlines()
    summary = next(
        (ln for ln in reversed(tail) if re.search(r"\d+ (passed|failed|error)", ln)), tail[-1] if tail else ""
    )
    return {
        "leg": label,
        "ok": r.returncode == 0,
        "exit": r.returncode,
        "seconds": round(seconds, 1),
        "summary": re.sub(r"\x1b\[[0-9;]*m", "", summary)[:300],
        "cwd": str(cwd),
        "nothing_collected": r.returncode == 5,
        "output_tail": "\n".join(tail[-40:]) if r.returncode != 0 else "",
    }


def store_free_leg(
    root: Path, files: list[str], sha: str = "HEAD", plant: dict[str, str] | None = None
) -> dict:
    """A detached worktree of `sha` with no stores in it, the two assertions, then the targeted tests.

    `plant` writes the named files into the worktree before the leg runs, and PRODUCTION PASSES
    NOTHING. It exists for one reason: a leg whose blocking has never been demonstrated is a leg
    nobody has tested, and the only place the demonstration is honest is inside a real store-free
    worktree -- a test file carries the project's `needs_local_data` wiring only when it sits under
    the worktree's own `tests/`, because that is how pytest finds a conftest. So
    `tests/test_fast_prepush.py` plants two files there, one reading a git-ignored store WITHOUT the
    marker and one WITH it, and requires the first to make this leg red and the second not to. The
    parameter adds files to a throwaway worktree; it cannot remove, skip or weaken anything, and the
    two assertions below run against the planted tree exactly as they do against an unplanted one.
    """
    parent = tempfile.mkdtemp(prefix="genomeos-storefree-")
    tree = Path(parent) / "tree"
    out: dict = {"worktree": str(tree)}
    added = subprocess.run(
        ["git", "worktree", "add", "-q", "--detach", str(tree), sha],
        cwd=root,
        capture_output=True,
        text=True,
    )
    out["worktree_added"] = added.returncode == 0
    if added.returncode != 0:
        out["ok"] = False
        out["why"] = f"REFUSED: could not make a store-free worktree of {sha}: {added.stderr.strip()}"
        return out
    try:
        for rel, text in (plant or {}).items():
            target = tree / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(text)
        out["planted"] = sorted(plant or {})
        out["stores"] = assert_stores_absent(tree)
        if not out["stores"]["ok"]:
            out["ok"] = False
            out["why"] = out["stores"]["why"]
            return out
        out["import_path"] = assert_import_path(root, tree)
        if not out["import_path"]["ok"]:
            out["ok"] = False
            out["why"] = out["import_path"]["why"]
            return out
        out["pytest"] = run_pytest(root, tree, files, "store-free")
        out["ok"] = out["pytest"]["ok"]
        return out
    finally:
        out["removal"] = remove_worktree(root, tree)
        if not out["removal"]["ok"]:
            out.setdefault("why", "")
            out["why"] += (
                f" (and the worktree {tree} was NOT removed: {out['removal'].get('detail', '')[:300]}. "
                f"It is reported and left alone -- never `rm -rf` a worktree; see "
                f"scripts/remove_worktree.py)"
            )


def remove_worktree(root: Path, tree: Path) -> dict:
    """Through `scripts/remove_worktree.py` and nothing else. `--force`: pytest leaves untracked files."""
    r = subprocess.run(
        [PYTHON, "scripts/remove_worktree.py", str(tree), "--force"],
        cwd=root,
        env=leg_env(root, None),
        capture_output=True,
        text=True,
    )
    return {"ok": r.returncode == 0, "exit": r.returncode, "detail": (r.stdout + r.stderr).strip()[-1200:]}


# ----------------------------------------------------------------------------- the gate: which one


def read_daily(root: Path) -> tuple[dict | None, Path, str | None]:
    """(the daily reading of CI, where it was looked for, why there is none). Never raises on a bad file."""
    path = last_green_ci.daily_read_path(root)
    if not path.exists():
        return None, path, "no such file"
    try:
        record = json.loads(path.read_text())
    except (OSError, ValueError) as e:
        return None, path, f"unreadable ({e.__class__.__name__}: {e})"
    if not isinstance(record, dict) or not record.get("conclusion") or not record.get("read_at"):
        return None, path, "readable but names no conclusion or no read_at"
    return record, path, None


def gate_decision(
    root: Path,
    remote_ref: str,
    env: dict[str, str] | None = None,
    now: datetime | None = None,
    max_age: int = DAILY_READ_MAX_SECONDS,
) -> dict:
    """Which gate a push takes. "fast" only when NOTHING below asks for the full suite.

    Every reason found is listed and not only the first, so the push says all of what it knows. It does
    not raise on anything it reads: a reading it cannot make is itself a reason for the full gate, and
    `scripts/pre-push.sh` starts from the full gate and leaves it only on this function's explicit
    "fast", so a crash here costs nine minutes and never a check.
    """
    env = dict(os.environ) if env is None else env
    now = now or datetime.now(UTC)
    full: list[str] = []
    asked = str(env.get(FULL_GATE_ENV, "")).strip()
    if asked not in ("", "0"):
        full.append(
            f"{FULL_GATE_ENV}={asked} asks for the full suite (the switch for a promotion or a paid-path sha)"
        )
    if remote_ref != DEV_REF:
        full.append(
            f"the push is to {remote_ref}, not {DEV_REF}: main and the promote-<sha> branches keep the full "
            f"suite, so main's bar does not move"
        )
    record, path, why = read_daily(root)
    fast_reason = ""
    if record is None:
        full.append(
            f"the daily CI is UNREAD: no reading ({why}) at {path}. Record one after reading CI: "
            f"`uv run python scripts/last_green_ci.py --daily CONCLUSION --sha SHA --run-id ID`"
        )
    else:
        conclusion = str(record.get("conclusion")).strip().lower()
        short = record.get("short") or str(record.get("sha") or "")[:7]
        if conclusion != last_green_ci.GREEN:
            full.append(
                f"the daily CI is RED: the reading taken {record.get('read_at')} says {conclusion!r} on "
                f"{short}. The full gate holds until a reading is green."
            )
        age = _age_seconds({"written": record.get("read_at")}, now)
        if age is None:
            full.append(f"the daily reading's read_at {record.get('read_at')!r} is not a time: UNREAD")
        elif age > max_age:
            full.append(
                f"the daily CI is UNREAD today: the last reading was taken {age / 3600:.1f} h ago, over the "
                f"{max_age / 3600:.0f} h limit"
            )
        elif age < -600:
            full.append(f"the daily reading's read_at {record.get('read_at')} is in the future: UNREAD")
        else:
            fast_reason = (
                f"a dev push, and the daily CI was read {conclusion} on {short} {age / 3600:.1f} h ago"
            )
    return {
        "mode": "full" if full else "fast",
        "reasons": full or [fast_reason],
        "remote_ref": remote_ref,
        "daily": record,
        "daily_path": str(path),
    }


# ----------------------------------------------------------------------------- the gate: the fast one


def fast_status_dir(root: Path) -> Path:
    """The directory the project verdicts live in, honouring the same override genomeos/verdict.py does."""
    override = os.environ.get("GENOMEOS_CHECK_STATUS_DIR")
    if override:
        return Path(override)
    return last_green_ci.git_common_dir(root) / last_green_ci.DIRNAME


def fast_status_path(root: Path, tree: str) -> Path:
    """`fast-<tree>.json`, and NEVER `status-<tree>.json` or `latest.json`.

    genomeos/verdict.py's `require` reads `status-<tree>.json` and nothing else, and the project's
    status files are keyed by TREE ALONE, so a fast verdict written under that name would do two wrong
    things at once: it could be read as the project verdict for its tree, and it would DESTROY the real
    project verdict for the same tree. Under its own name it can do neither; its own schema and its
    scope "fast" refuse it there even if somebody copies it across.
    """
    return fast_status_dir(root) / f"fast-{tree}.json"


def tree_hash(root: Path) -> str:
    """The working-tree hash genomeos/verdict.py records, so the fast verdict names its tree the same way."""
    r = subprocess.run(
        [PYTHON, "-m", "genomeos.verdict", "tree"],
        cwd=root,
        env=leg_env(root, root),
        capture_output=True,
        text=True,
    )
    return r.stdout.strip() if r.returncode == 0 else ""


def ruff_leg(root: Path, args: list[str], label: str) -> dict:
    """One ruff run over the project, as `scripts/check.sh` runs it, with this interpreter's ruff."""
    started = time.monotonic()
    r = subprocess.run(
        [PYTHON, "-m", "ruff", *args], cwd=root, env=leg_env(root, root), capture_output=True, text=True
    )
    tail = (r.stdout + r.stderr).strip().splitlines()
    return {
        "leg": label,
        "ok": r.returncode == 0,
        "exit": r.returncode,
        "seconds": round(time.monotonic() - started, 1),
        "summary": tail[-1][:300] if tail else "",
        "output_tail": "\n".join(tail[-40:]) if r.returncode != 0 else "",
    }


def write_fast_status(root: Path, payload: dict) -> Path:
    """Atomically, beside the project verdicts and under a name they never read; prune the old ones."""
    path = fast_status_path(root, payload["tree"])
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp-fast-", suffix=".json", dir=str(path.parent))
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True, default=str)
            fh.write("\n")
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    except BaseException:
        os.unlink(tmp)
        raise
    now = time.time()
    for old in path.parent.glob("fast-*.json"):
        try:
            if old != path and now - old.stat().st_mtime > FAST_PRUNE_AFTER_SECONDS:
                old.unlink()
        except OSError:
            continue
    return path


def gate(root: Path, sha: str, tree: str, remote_ref: str = DEV_REF, now: datetime | None = None) -> dict:
    """THE FAST GATE, as ruled on 2026-10-03: lint, format, then ALWAYS_RUN and the window in both legs.

    Run by `scripts/pre-push.sh` inside its clean worktree of the commit being pushed, so `root` is that
    worktree, its stores are the read-only clones, and the normal leg runs THERE (with the import-path
    assertion, so it cannot silently test the shared checkout). Stops at the first red leg, as
    `scripts/check.sh` does. Writes its verdict to `fast-<tree>.json` with scope "fast", green only if
    every leg is green AND the tree measured before and after is the tree being pushed.
    """
    started = time.time()
    result: dict = {
        "schema": FAST_SCHEMA,
        "scope": FAST_SCOPE,
        "project_verdict": False,
        "what_this_is": (
            "the FAST gate's verdict for one dev push (ruff check, ruff format --check, ALWAYS_RUN and the "
            "test files changed since the last green CI, normal and store-free). It is NOT a project "
            "verdict: genomeos/verdict.py never reads this file, and the full suite runs in the daily CI, "
            "before the promotion to main and on any paid path."
        ),
        "sha": sha,
        "tree": tree,
        "remote_ref": remote_ref,
        "legs": [],
        "failed_leg": "",
        "root": str(root),
        "host": socket.gethostname(),
        "pid": os.getpid(),
        "started_at": datetime.fromtimestamp(started, UTC).isoformat(),
    }
    result["tree_begin"] = tree_hash(root)
    legs: list[dict] = result["legs"]
    files: list[str] = []
    for args, label in ((["check", "."], "ruff-check"), (["format", "--check", "."], "ruff-format")):
        leg = ruff_leg(root, args, label)
        legs.append(leg)
        print(f"fast-gate: {label}: exit {leg['exit']} in {leg['seconds']}s -- {leg['summary']}")
        if not leg["ok"]:
            result["failed_leg"] = label
            print(leg["output_tail"], file=sys.stderr)
            break
    if not result["failed_leg"]:
        p = plan(root, now=now)
        result["plan"] = p
        for notice in p["window"]["notices"] + p["selection"]["notices"]:
            print(notice)
        files = p["selection"]["files"]
        print(
            f"fast-gate: window {p['window']['case']} from {(p['window']['base'] or 'none')[:12]}: "
            f"{p['selection']['changed_count']} changed test file(s), {len(files)} to run with ALWAYS_RUN"
        )
        imported = assert_import_path(root, root)
        result["import_path"] = imported
        if not imported["ok"]:
            result["failed_leg"] = "normal"
            print(imported["why"], file=sys.stderr)
        else:
            normal = run_pytest(root, root, files, "normal")
            legs.append(normal)
            print(
                f"fast-gate: normal leg: exit {normal['exit']} in {normal['seconds']}s "
                f"-- {normal.get('summary')}"
            )
            if not normal["ok"]:
                result["failed_leg"] = "normal"
                print(normal.get("output_tail", ""), file=sys.stderr)
    if not result["failed_leg"]:
        sf = store_free_leg(root, files, sha=sha)
        legs.append({"leg": "store-free", **sf})
        if sf.get("stores"):
            print(f"fast-gate: store-free assertion: {sf['stores']['why']}")
        if sf.get("import_path"):
            print(f"fast-gate: import-path assertion: {sf['import_path']['why']}")
        if sf.get("pytest"):
            print(
                f"fast-gate: store-free leg: exit {sf['pytest']['exit']} in {sf['pytest']['seconds']}s "
                f"-- {sf['pytest'].get('summary', '')}"
            )
        if not sf["ok"]:
            result["failed_leg"] = "store-free"
            print(sf.get("why", ""), file=sys.stderr)
            if sf.get("pytest"):
                print(sf["pytest"].get("output_tail", ""), file=sys.stderr)
    result["tree_end"] = tree_hash(root)
    same_tree = result["tree_begin"] == result["tree_end"] == tree
    if result["failed_leg"]:
        result["verdict"], result["reason"] = "red", f"the {result['failed_leg']} leg is red"
    elif not same_tree:
        result["verdict"] = "red"
        result["reason"] = (
            f"the tree moved or was not the one pushed: begin {result['tree_begin'] or 'unreadable'}, "
            f"end {result['tree_end'] or 'unreadable'}, pushed {tree}"
        )
    else:
        result["verdict"], result["reason"] = "green", "every fast leg green on the tree being pushed"
    result["finished_at"] = datetime.now(UTC).isoformat()
    result["seconds"] = round(time.time() - started, 1)
    result["path"] = str(write_fast_status(root, result))
    print(f"fast-gate: verdict {result['verdict']} ({result['reason']}) -> {result['path']}")
    return result


def require_fast(root: Path, tree: str, out=None, err=None) -> int:
    """Answer, for ONE tree, whether a complete green FAST verdict exists. Anything else refuses."""
    out = sys.stdout if out is None else out
    err = sys.stderr if err is None else err
    tree = tree.strip()
    if not re.fullmatch(r"[0-9a-f]{40}", tree):
        print(f"fast-gate: REFUSED: {tree!r} is not a tree hash", file=err)
        return EXIT_REFUSED
    path = fast_status_path(root, tree)
    try:
        status = json.loads(path.read_text())
    except (OSError, ValueError) as e:
        print(f"fast-gate: REFUSED: no readable fast verdict for tree {tree} at {path} ({e})", file=err)
        return EXIT_REFUSED
    problems = []
    if status.get("schema") != FAST_SCHEMA:
        problems.append(f"schema {status.get('schema')!r}, not {FAST_SCHEMA!r}")
    if status.get("scope") != FAST_SCOPE:
        problems.append(f"scope {status.get('scope')!r}, not {FAST_SCOPE!r}")
    for key in ("tree", "tree_begin", "tree_end"):
        if status.get(key) != tree:
            problems.append(f"{key}={status.get(key)!r}, not this tree")
    if status.get("verdict") != "green":
        problems.append(f"verdict {status.get('verdict')!r}: {status.get('reason')}")
    if problems:
        print(f"fast-gate: REFUSED for tree {tree}: " + "; ".join(problems), file=err)
        return EXIT_REFUSED
    selection = (status.get("plan") or {}).get("selection") or {}
    window = (status.get("plan") or {}).get("window") or {}
    print(
        f"{FAST_OK_TOKEN} {tree} scope={FAST_SCOPE} files={selection.get('count')} "
        f"window={window.get('case')} finished={status.get('finished_at')}",
        file=out,
    )
    return 0


# ----------------------------------------------------------------------------- the whole path


def plan(root: Path = ROOT, now: datetime | None = None) -> dict:
    window = resolve_window(root, now=now)
    changed = changed_test_files(root, window.base)
    selection = select(root, changed)
    return {"window": window.as_dict(), "selection": selection.as_dict()}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", action="store_true", help="print the window and the selection, run nothing")
    ap.add_argument("--only", choices=("normal", "store-free"), help="run one leg only")
    ap.add_argument("--sha", default="HEAD", help="the commit the store-free worktree is made of")
    ap.add_argument("--root", type=Path, default=ROOT)
    ap.add_argument("--json", action="store_true", help="print the whole result as JSON at the end")
    ap.add_argument("--decide", action="store_true", help="print which gate a push takes, GATE fast|full")
    ap.add_argument("--gate", action="store_true", help="run the fast gate and write its fast-<tree> verdict")
    ap.add_argument("--require-fast", action="store_true", help="refuse unless a green fast verdict exists")
    ap.add_argument("--tree", help="the tree of the commit being pushed (--gate, --require-fast)")
    ap.add_argument("--remote-ref", default=DEV_REF, help="the ref being pushed to (--decide, --gate)")
    args = ap.parse_args(argv)
    root = args.root

    if args.decide:
        d = gate_decision(root, args.remote_ref)
        # ONE first line the hook matches on; the reasons after it are for the person pushing
        print(f"GATE {d['mode']}: {len(d['reasons'])} reason(s)")
        for reason in d["reasons"]:
            print(f"  - {reason}")
        return 0
    if args.require_fast:
        if not args.tree:
            ap.error("--require-fast needs --tree")
        return require_fast(root, args.tree)
    if args.gate:
        sha = git(root, "rev-parse", "--verify", f"{args.sha}^{{commit}}")
        tree = args.tree or git(root, "rev-parse", f"{sha}^{{tree}}")
        result = gate(root, sha, tree, args.remote_ref)
        return 0 if result["verdict"] == "green" else 1

    p = plan(root)
    for notice in p["window"]["notices"] + p["selection"]["notices"]:
        print(notice)
    print(
        f"fast-prepush: window {p['window']['case']}, base "
        f"{(p['window']['base'] or 'none')[:12]}: {p['selection']['changed_count']} changed test file(s), "
        f"{p['selection']['count']} to run (ALWAYS_RUN included)"
    )
    if args.plan:
        print(json.dumps(p, indent=1))
        return 0

    files = p["selection"]["files"]
    result: dict = {"plan": p, "legs": []}
    status = 0
    if args.only != "store-free":
        normal = run_pytest(root, None, files, "normal")
        result["legs"].append(normal)
        print(
            f"fast-prepush: normal leg: exit {normal['exit']} in {normal['seconds']}s "
            f"-- {normal.get('summary', '')}"
        )
        if not normal["ok"]:
            status = 1
            print(normal.get("output_tail", ""), file=sys.stderr)
    if args.only != "normal":
        sf = store_free_leg(root, files, sha=args.sha)
        result["legs"].append(sf)
        if sf.get("stores"):
            print(f"fast-prepush: store-free assertion: {sf['stores']['why']}")
        if sf.get("import_path"):
            print(f"fast-prepush: import-path assertion: {sf['import_path']['why']}")
        if sf.get("pytest"):
            print(
                f"fast-prepush: store-free leg: exit {sf['pytest']['exit']} in "
                f"{sf['pytest']['seconds']}s -- {sf['pytest'].get('summary', '')}"
            )
        print(f"fast-prepush: worktree removed through scripts/remove_worktree.py: {sf['removal']['ok']}")
        if not sf["ok"]:
            status = 1
            print(sf.get("why", ""), file=sys.stderr)
            if sf.get("pytest"):
                print(sf["pytest"].get("output_tail", ""), file=sys.stderr)
    result["ok"] = status == 0
    if args.json:
        print(json.dumps(result, indent=1, default=str))
    print(
        "fast-prepush: " + ("both legs green" if status == 0 else "RED -- see above; the push is not ready")
    )
    return status


if __name__ == "__main__":
    sys.exit(main())
