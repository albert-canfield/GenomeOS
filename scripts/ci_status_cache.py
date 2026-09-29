#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Write the CI status the Progress tab shows, so that no web request waits on GitHub.

Reading CI takes several seconds of `gh` calls, too slow for a page load, so this script reads it
once and writes a small git-ignored file, data/cache/ci_status.json, stamped with the time it was
written. /api/state reads that file and shows its age beside every figure it takes from it, so a
stale reading says it is stale.

What it records:
- the last CI runs of .github/workflows/ci.yml, and for the latest completed one its jobs, the
  pytest summary line of the `test` job and the tests that failed;
- the latest green run, and how many completed runs in a row have not been green;
- the tips of origin/main and origin/dev, and the latest completed `test` check on main's sha;
- the promotion gate's verdict: scripts/promote_main.sh in its dry run, which fetches dev and main
  into the remote-tracking refs, reads GitHub, prints what it would do and changes nothing. This
  script never asks the gate for anything but a dry run.

Usage: uv run python scripts/ci_status_cache.py [--out PATH] [--no-promote] [--limit N]
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from collections.abc import Callable
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "cache" / "ci_status.json"
WORKFLOW = "ci.yml"
RUN_FIELDS = "databaseId,headSha,headBranch,event,status,conclusion,createdAt,updatedAt,url"
#: the only way this script calls the promotion gate: a dry run for the newest green revision
PROMOTE = ("scripts/promote_main.sh", "--dry-run", "--latest-green")

Runner = Callable[[list[str]], tuple[int, str, str]]


def run(cmd: list[str], timeout: int = 180) -> tuple[int, str, str]:
    """A command's exit code, stdout and stderr; a missing program or a timeout is an exit code too."""
    try:
        p = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True, timeout=timeout)
    except (OSError, subprocess.SubprocessError) as e:
        return 127, "", str(e)
    return p.returncode, p.stdout, p.stderr


_STAMP = re.compile(r"\d{4}-\d\d-\d\dT[\d:.]+Z\s(.*)$")
_SUMMARY = re.compile(r"^=*\s*(\d+ (?:failed|passed)\b.*? in [\d.]+s\b.*?)\s*=*$")
_FAILURE = re.compile(r"^_{3,} (\S.*?) _{3,}$")  # a failure's section header: ____ test_name ____


def pytest_summary(log: str) -> tuple[str | None, list[str]]:
    """The last pytest summary line of a job log and the tests that failed, with the timestamps cut.

    A failed test is named by its FAILED line where pytest prints one, and by its section header
    in the failures section otherwise (CI runs pytest with -rs, which leaves the FAILED lines out).
    """
    summary, failed = None, []
    for line in log.splitlines():
        m = _STAMP.search(line)
        body = (m.group(1) if m else line).strip()
        name = body[len("FAILED ") :] if body.startswith("FAILED ") else None
        header = _FAILURE.match(body)
        name = name or (header.group(1) if header else None)
        if name and name not in failed:
            failed.append(name)
        s = _SUMMARY.match(body)
        if s:
            summary = s.group(1)
    return summary, failed


def _run_row(r: dict) -> dict:
    return {
        "id": r.get("databaseId"),
        "sha": r.get("headSha"),
        "branch": r.get("headBranch"),
        "event": r.get("event"),
        "status": r.get("status"),
        "conclusion": r.get("conclusion"),
        "created": r.get("createdAt"),
        "url": r.get("url"),
    }


def collect(sh: Runner = run, promote: bool = True, limit: int = 20) -> dict:
    """Everything the status file holds, gathered through `sh` so a test can stand in for GitHub."""
    assert "--push" not in PROMOTE and "--dry-run" in PROMOTE
    errors: list[str] = []
    out: dict = {"workflow": WORKFLOW}

    code, repo, err = sh(["gh", "repo", "view", "--json", "nameWithOwner", "-q", ".nameWithOwner"])
    out["repo"] = repo.strip() if code == 0 else None
    if code != 0:
        errors.append(f"gh repo view: {err.strip()[:300]}")

    if promote:
        code, so, se = sh(list(PROMOTE))
        out["promote_dry_run"] = {
            "args": list(PROMOTE[1:]),
            "exit": code,
            "output": (so + se).strip()[-1500:],
        }

    code, raw, err = sh(
        ["gh", "run", "list", "--workflow", WORKFLOW, "--limit", str(limit), "--json", RUN_FIELDS]
    )
    runs = [_run_row(r) for r in json.loads(raw)] if code == 0 and raw.strip() else []
    if code != 0:
        errors.append(f"gh run list: {err.strip()[:300]}")
    runs.sort(key=lambda r: r["created"] or "", reverse=True)
    out["runs"] = runs
    done = [r for r in runs if r["status"] == "completed"]
    latest = dict(done[0]) if done else None
    if latest:
        code, raw, err = sh(["gh", "run", "view", str(latest["id"]), "--json", "jobs"])
        jobs = json.loads(raw).get("jobs", []) if code == 0 and raw.strip() else []
        if code != 0:
            errors.append(f"gh run view: {err.strip()[:300]}")
        latest["jobs"] = [{"name": j.get("name"), "conclusion": j.get("conclusion")} for j in jobs]
        test = next((j for j in jobs if j.get("name") == "test"), None)
        latest["test_summary"], latest["failed_tests"] = None, []
        if test:
            code, log, err = sh(["gh", "run", "view", "--job", str(test.get("databaseId")), "--log"])
            if code == 0:
                latest["test_summary"], latest["failed_tests"] = pytest_summary(log)
            else:
                errors.append(f"gh run view --log: {err.strip()[:300]}")
    out["latest_completed"] = latest
    out["latest_success"] = next((r for r in done if r["conclusion"] == "success"), None)
    streak = 0
    for r in done:
        if r["conclusion"] == "success":
            break
        streak += 1
    out["consecutive_failures"] = streak

    refs = {}
    for name in ("main", "dev"):
        code, sha, _ = sh(["git", "rev-parse", "--verify", "--quiet", f"refs/remotes/origin/{name}"])
        refs[name] = sha.strip() if code == 0 else None
    out["refs"] = refs
    main_check = None
    if out["repo"] and refs.get("main"):
        api = f"repos/{out['repo']}/commits/{refs['main']}/check-runs?check_name=test&filter=all&per_page=100"
        code, raw, err = sh(["gh", "api", api])
        if code == 0 and raw.strip():
            checks = [
                c
                for c in json.loads(raw).get("check_runs", [])
                if (c.get("app") or {}).get("slug") == "github-actions" and c.get("status") == "completed"
            ]
            checks.sort(key=lambda c: c.get("completed_at") or "")
            c = checks[-1] if checks else None
            main_check = {
                "sha": refs["main"],
                "conclusion": c.get("conclusion") if c else None,
                "completed_at": c.get("completed_at") if c else None,
                "url": c.get("html_url") if c else None,
                "completed_runs": len(checks),
            }
        else:
            errors.append(f"gh api check-runs: {err.strip()[:300]}")
    out["main_check"] = main_check
    out["errors"] = errors
    return out


def write(status: dict, path: Path = OUT, now: datetime | None = None) -> Path:
    """The status with the time it was written, replacing the old file in one rename."""
    stamp = (now or datetime.now(UTC)).strftime("%Y-%m-%dT%H:%M:%SZ")
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps({"written": stamp, **status}, indent=1) + "\n")
    tmp.replace(path)
    return path


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    ap.add_argument("--out", type=Path, default=OUT)
    ap.add_argument("--no-promote", action="store_true", help="skip the promotion gate's dry run")
    ap.add_argument("--limit", type=int, default=20, help="CI runs to list")
    a = ap.parse_args(argv)
    status = collect(promote=not a.no_promote, limit=a.limit)
    path = write(status, a.out)
    latest = status.get("latest_completed") or {}
    print(
        f"{path}: latest {latest.get('conclusion')} on {str(latest.get('sha') or '')[:10]}"
        f" ({latest.get('test_summary')}); {status['consecutive_failures']} not green in a row;"
        f" main test {(status.get('main_check') or {}).get('conclusion')};"
        f" promotion dry run exit {(status.get('promote_dry_run') or {}).get('exit')}"
    )
    for e in status["errors"]:
        print(f"  {e}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
