#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""The fast pre-push: the TARGETED tests, run TWICE -- once in this tree, once with no data stores.

    uv run python scripts/fast_prepush.py             # both legs, on the window described below
    uv run python scripts/fast_prepush.py --plan      # print the window and the selection, run nothing
    uv run python scripts/fast_prepush.py --only store-free

WHAT IT IS FOR, and it is one defect class and not a second full gate. On 2026-10-03 three CI runs went
red -- 3 failures, then 6 -- on failures the LOCAL FULL SUITE COULD NOT SEE: every one was a test that
assumes this machine's git-ignored data (`data/knowledge`, `data/cache`, `data/reference`). This
checkout has those stores, CI is a fresh checkout and does not, so the suite is green here and red
there by construction. A lesson was written about the class and the very next commit reproduced it. The
full gate costs 17 to 26 minutes and still misses it, because the environment it runs in is the one
that hides it. The leg below removes the environment instead of adding tests, and it catches exactly
this class, in about a minute, at the hand that wrote the code.

It does NOT replace `scripts/pre-push.sh`. That gate still runs, unchanged, and still refuses a push
without a status-file verdict for the exact tree. This is a fast path in front of it.

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
    args = ap.parse_args(argv)
    root = args.root

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
