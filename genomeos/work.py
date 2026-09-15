"""The work board: who is working on what right now, so the Progress tab says it instead of anyone asking.

Several contributors work on one checkout at once. Each keeps one entry, written by
`genomeos work start | update | done`, as data/work/<who>.json (git-ignored runtime state, like
data/jobs): the task, its roadmap area, the files held and what comes next. The board adds what
the checkout itself shows: files changed but not committed, and commits of the last day.
"""

from __future__ import annotations

import json
import re
import subprocess
import time
from pathlib import Path
from typing import Any

STALE_AFTER = 6 * 3600  # an entry not touched for six hours is shown as stale, not hidden
KEEP_DONE = 24 * 3600  # finished entries stay visible for a day
_WHO = re.compile(r"^[A-Za-z0-9_.-]{1,64}$")


def _dir(root: Path) -> Path:
    return root / "data" / "work"


def _path(root: Path, who: str) -> Path:
    if not _WHO.match(who):
        raise ValueError(f"not a contributor name: {who!r}")
    return _dir(root) / f"{who}.json"


def read(root: Path, who: str) -> dict[str, Any] | None:
    p = _path(root, who)
    try:
        return json.loads(p.read_text())
    except (OSError, json.JSONDecodeError):
        return None


def _write(root: Path, entry: dict[str, Any]) -> dict[str, Any]:
    p = _path(root, entry["who"])
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = p.with_suffix(".tmp")
    tmp.write_text(json.dumps(entry, indent=2) + "\n")
    tmp.replace(p)
    return entry


def _files(files: list[str] | None) -> list[str]:
    return sorted({x.strip() for f in (files or []) for x in f.split(",") if x.strip()})


def start(
    root: Path,
    who: str,
    task: str,
    area: str | None = None,
    files: list[str] | None = None,
    next_step: str | None = None,
    now: float | None = None,
) -> dict[str, Any]:
    """Begin a task: replaces the contributor's entry."""
    now = time.time() if now is None else now
    return _write(
        root,
        {
            "who": who,
            "task": task,
            "area": (area or "").upper() or None,
            "files": _files(files),
            "next": next_step,
            "state": "working",
            "note": None,
            "started": now,
            "updated": now,
            "finished": None,
        },
    )


def update(
    root: Path,
    who: str,
    note: str | None = None,
    task: str | None = None,
    area: str | None = None,
    files: list[str] | None = None,
    next_step: str | None = None,
    state: str | None = None,
    now: float | None = None,
) -> dict[str, Any]:
    """Change part of an entry (progress note, files, next step, waiting or working)."""
    entry = read(root, who)
    if entry is None:
        if not task:
            raise ValueError(f"{who} has no entry; start one with a task")
        return start(root, who, task, area, files, next_step, now)
    if state is not None and state not in ("working", "waiting"):
        raise ValueError("state is working or waiting; use done to finish")
    for key, val in (("note", note), ("task", task), ("next", next_step), ("state", state)):
        if val is not None:
            entry[key] = val
    if area is not None:
        entry["area"] = area.upper() or None
    if files is not None:
        entry["files"] = _files(files)
    entry["updated"] = time.time() if now is None else now
    return _write(root, entry)


def done(root: Path, who: str, note: str | None = None, now: float | None = None) -> dict[str, Any]:
    """Finish the task: the entry stays on the board as finished for a day."""
    entry = read(root, who)
    if entry is None:
        raise ValueError(f"{who} has no entry")
    now = time.time() if now is None else now
    entry.update(state="done", finished=now, updated=now)
    if note:
        entry["note"] = note
    return _write(root, entry)


def board(root: Path, now: float | None = None) -> list[dict[str, Any]]:
    """Every entry, working first, with its age; finished entries older than a day are dropped."""
    now = time.time() if now is None else now
    out = []
    d = _dir(root)
    for p in sorted(d.glob("*.json")) if d.exists() else []:
        try:
            e = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if e.get("state") == "done" and now - (e.get("finished") or 0) > KEEP_DONE:
            continue
        updated = e.get("updated")
        e["age"] = round(now - updated) if isinstance(updated, (int, float)) else 0
        if e.get("state") != "done" and e["age"] > STALE_AFTER:
            e["state"] = "stale"
        out.append(e)
    order = {"working": 0, "waiting": 1, "stale": 2, "done": 3}
    return sorted(out, key=lambda e: (order.get(e.get("state"), 9), e["age"]))


def _git(root: Path, *args: str) -> str:
    try:
        r = subprocess.run(  # noqa: S603
            ["git", *args], cwd=root, capture_output=True, text=True, timeout=10, check=False
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    return r.stdout if r.returncode == 0 else ""


def uncommitted(root: Path) -> list[dict[str, str]]:
    """Files changed in the checkout and not committed yet: work in flight, whoever holds it."""
    out = []
    for line in _git(root, "status", "--porcelain", "--untracked-files=all").splitlines():
        if len(line) < 4:
            continue
        code, path = line[:2], line[3:].strip().strip('"')
        if " -> " in path:
            path = path.split(" -> ", 1)[1]
        out.append(
            {"path": path, "change": "new" if code == "??" else "deleted" if "D" in code else "modified"}
        )
    return out


def recent_commits(root: Path, hours: int = 24, limit: int = 30) -> list[dict[str, Any]]:
    """Commits of the last hours on the current branch, the subject cut at its first clause."""
    raw = _git(root, "log", f"--since={hours}.hours", f"-{limit}", "--format=%h%x09%ct%x09%s")
    out = []
    for line in raw.splitlines():
        parts = line.split("\t", 2)
        if len(parts) != 3:
            continue
        sha, ts, subject = parts
        head = re.split(r"(?<=[a-z0-9)])[:;]\s", subject, maxsplit=1)[0]
        out.append({"sha": sha, "time": int(ts), "title": head[:160], "subject": subject})
    return out
