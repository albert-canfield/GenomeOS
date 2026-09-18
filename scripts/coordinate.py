# SPDX-License-Identifier: AGPL-3.0-or-later
"""One screen of who holds what, so a session does not have to ask another session.

    uv run python scripts/coordinate.py            # the whole board
    uv run python scripts/coordinate.py --mine genomeos-0e
    uv run python scripts/coordinate.py --json

**Why this exists.** On 2026-09-17 two sessions produced 96 commits and roughly twenty long messages
to each other. The verification in those messages was worth it — it caught seven wrong numbers — but
almost none of the *coordination* was: who owns which file, what has landed, what is unclaimed, and
whether anything is about to collide are all answerable from artefacts already in the repository. A
session that asks another session those questions is paying for a round trip to read a file it can
read itself, and the answer arrives as a summary, which this week established is the least reliable
form of a fact.

So this joins what already exists and invents nothing: `roadmap.load` for areas, owners and §5's
ordered steps, `work.board` for who is working and how stale they are, `work.uncommitted` for the
working tree, and git for worktrees and unmerged branches. Every number it prints is read at the
moment it prints, and nothing here is a place to *write* status — the board and the roadmap stay the
only writable records, because a third one would drift from both.

**What it will not do.** It does not assign work, because the roadmap does that; it does not judge
whether a lane is behind, because the board's own note is the lane's to write; and it does not delete
anything. It reports, and the collisions it reports are the ones a human or a session has to resolve.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path
from typing import Any

from genomeos import roadmap, work

STALE_HOURS = 3.0  # a board entry older than this is reported as stale rather than as current work


def git(root: Path, *args: str) -> str:
    done = subprocess.run(["git", *args], cwd=root, capture_output=True, text=True, check=False, env=_env())
    return done.stdout.strip()


def _env() -> dict[str, str]:
    import os

    e = dict(os.environ)
    e.pop("GIT_INDEX_FILE", None)  # never read another session's private index
    return e


def worktrees(root: Path) -> list[dict[str, Any]]:
    """Agent worktrees, whether each is locked, and whether its branch is merged into dev."""
    out: list[dict[str, Any]] = []
    block: dict[str, Any] = {}
    for line in git(root, "worktree", "list", "--porcelain").splitlines() + [""]:
        if not line.strip():
            if block.get("path") and block.get("branch"):
                out.append(block)
            block = {}
            continue
        key, _, value = line.partition(" ")
        if key == "worktree":
            block = {"path": value, "locked": False, "branch": ""}
        elif key == "branch":
            block["branch"] = value.replace("refs/heads/", "")
        elif key == "locked":
            block["locked"] = True
    for w in out:
        merged = subprocess.run(
            ["git", "merge-base", "--is-ancestor", w["branch"], "dev"],
            cwd=root,
            capture_output=True,
            env=_env(),
        )
        w["merged_into_dev"] = merged.returncode == 0
        # content-identical is the check that matters: today every "unmerged" agent branch had been
        # rebased and reworded, so ancestry and patch-id both said unmerged while nothing differed
        w["files_differing_from_dev"] = len(
            [x for x in git(root, "diff", "--name-only", "dev", w["branch"]).splitlines() if x]
        )
    return out


def collisions(board: list[dict[str, Any]], dirty: list[dict[str, str]]) -> dict[str, Any]:
    """Two kinds: two live sessions declaring one file, and a dirty file nobody live has declared."""
    live = [e for e in board if e["state"] != "done" and e["age"] / 3600.0 <= STALE_HOURS]
    claimed: dict[str, list[str]] = {}
    for e in board:
        for f in e.get("files") or []:
            claimed.setdefault(f, []).append(e["who"])
    shared = {f: who for f, who in claimed.items() if len({*who}) > 1}
    live_files = {f for e in live for f in (e.get("files") or [])}
    unclaimed_dirty = [d["path"] for d in dirty if d["path"] not in live_files]
    return {
        "files_claimed_by_more_than_one_session": shared,
        "dirty_files_no_live_session_has_declared": unclaimed_dirty,
        "reading": (
            "a file claimed twice is a merge waiting to happen; a dirty file nobody live has declared"
            " belongs to a session that has stopped, and committing it guesses on their behalf"
        ),
    }


def report(root: Path) -> dict[str, Any]:
    rm = roadmap.load(root)
    board = work.board(root)
    dirty = work.uncommitted(root)
    wts = worktrees(root)
    live = [e for e in board if e["state"] != "done" and e["age"] / 3600.0 <= STALE_HOURS]
    stale = [e for e in board if e["state"] != "done" and e["age"] / 3600.0 > STALE_HOURS]
    steps = rm["next_steps"]
    owners_live = {e["who"] for e in live}
    # an area whose declared owner has no live entry is work nobody is holding right now
    unheld = [
        {"area": a["letter"], "title": a["title"], "owners": a.get("owners") or []}
        for a in rm["areas"]
        if a["counts"].get("planned") and not (set(a.get("owners") or []) & owners_live)
    ]
    return {
        "generated": time.strftime("%Y-%m-%d %H:%M"),
        "dev": git(root, "log", "--oneline", "-1", "dev"),
        "behind_origin": git(root, "rev-list", "--count", "dev..origin/dev") or "0",
        "ahead_of_origin": git(root, "rev-list", "--count", "origin/dev..dev") or "0",
        "sessions_live": [
            {"who": e["who"], "area": e["area"], "task": e["task"][:90], "files": e.get("files") or []}
            for e in live
        ],
        "sessions_stale": [
            {"who": e["who"], "hours": round(e["age"] / 3600.0, 1), "task": e["task"][:70]} for e in stale
        ],
        "steps": [
            {"number": s["number"], "progress": s["progress"], "title": s["title"][:80]} for s in steps
        ],
        "steps_by_progress": {
            k: sum(1 for s in steps if s["progress"] == k) for k in ("done", "partial", "blocked", "planned")
        },
        "areas_with_planned_work_and_nobody_live": unheld,
        "worktrees": wts,
        "worktree_branches_with_content_not_on_dev": [
            w["branch"] for w in wts if w["files_differing_from_dev"]
        ],
        "uncommitted": [d["path"] for d in dirty],
        "collisions": collisions(board, dirty),
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--json", action="store_true")
    ap.add_argument("--mine", default="", help="a session name, to print what it holds first")
    args = ap.parse_args(argv)

    r = report(Path("."))
    if args.json:
        print(json.dumps(r, indent=1))
        return 0

    print(f"dev  {r['dev']}")
    if r["ahead_of_origin"] != "0" or r["behind_origin"] != "0":
        print(f"     ahead of origin {r['ahead_of_origin']}, behind {r['behind_origin']}")
    print()
    if args.mine:
        held = [s for s in r["sessions_live"] if s["who"] == args.mine]
        print(f"YOU ({args.mine}): {held[0]['task'] if held else 'no live board entry'}")
        if held and held[0]["files"]:
            print(f"     files: {', '.join(held[0]['files'])}")
        print()

    print(f"LIVE ({len(r['sessions_live'])} session(s) within {STALE_HOURS}h)")
    for s in r["sessions_live"]:
        print(f"  {s['who']:<13} [{s['area'] or '-'}] {s['task']}")
    if not r["sessions_live"]:
        print("  nobody")
    print(f"\nSTALE ({len(r['sessions_stale'])}) — their files are not free, but nobody is watching them")
    for s in r["sessions_stale"][:8]:
        print(f"  {s['who']:<13} {s['hours']:>6.1f}h  {s['task']}")

    print("\nROADMAP §5")
    for s in r["steps"]:
        print(f"  {s['number']:<4} {s['progress']:<8} {s['title']}")
    print(f"  totals: {r['steps_by_progress']}")

    if r["areas_with_planned_work_and_nobody_live"]:
        print("\nAREAS WITH PLANNED WORK AND NOBODY LIVE")
        for a in r["areas_with_planned_work_and_nobody_live"]:
            print(f"  {a['area']}  {a['title'][:66]}  (owners: {', '.join(a['owners']) or 'none named'})")

    print(f"\nWORKTREES ({len(r['worktrees'])})")
    for w in r["worktrees"]:
        flag = "locked" if w["locked"] else "      "
        diff = w["files_differing_from_dev"]
        note = f"{diff} file(s) differ from dev" if diff else "content identical to dev"
        print(f"  {flag}  {w['branch'][:46]:<46} {note}")

    c = r["collisions"]
    print("\nCOLLISIONS")
    if c["files_claimed_by_more_than_one_session"]:
        for f, who in c["files_claimed_by_more_than_one_session"].items():
            print(f"  claimed twice: {f} — {', '.join(sorted(set(who)))}")
    if c["dirty_files_no_live_session_has_declared"]:
        for f in c["dirty_files_no_live_session_has_declared"]:
            print(f"  dirty, undeclared: {f}")
    if not any(
        c[k] for k in ("files_claimed_by_more_than_one_session", "dirty_files_no_live_session_has_declared")
    ):
        print("  none")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
