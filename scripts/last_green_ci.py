#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Record the sha CI last went GREEN on, so the fast pre-push can widen its window with no network call.

    uv run python scripts/last_green_ci.py --sha 747b7ff --run-id 18123456789 --date 2026-10-01T09:14:00Z
    uv run python scripts/last_green_ci.py --from-cache      # from data/cache/ci_status.json's latest_success
    uv run python scripts/last_green_ci.py --show            # print what is recorded, change nothing

    # the DAILY READING, which decides whether a dev push may take the fast gate at all
    uv run python scripts/last_green_ci.py --daily success --sha 747b7ff --run-id 1812345678
    uv run python scripts/last_green_ci.py --daily failure --sha 747b7ff --run-id 18123456789
    uv run python scripts/last_green_ci.py --daily-from-cache   # the latest completed dev run in the cache

THE DAILY READING, added 2026-10-10 when the fast gate was wired (ROADMAP section 5, E1). The ruling of
2026-10-03 keeps the full suite where it decides something, and the first place it names is "the daily
CI, read first thing each morning -- if it is red or unread, pushes revert to the full gate until it is
green". A push hook cannot ask GitHub (the reason above), so the reading is written down here, beside
the last-green record, as `<common>/genomeos-check/daily-ci`: the conclusion of the latest completed CI
run on dev, its sha and run id, and `read_at`, the moment CI was actually read. `scripts/fast_prepush.py
--decide` reads it and takes the FULL gate when it is absent (unread), when its conclusion is anything
but `success` (red, cancelled, skipped: not green), or when `read_at` is older than its limit (unread
today). A green reading also moves the last-green record to the same sha, because it is one. A red one
leaves the last-green record where it is: the window still starts at the last sha CI passed on.

`--daily-from-cache` takes `read_at` from the cache's own `written` stamp and NOT from the clock: the
cache is a reading taken when `scripts/ci_status_cache.py` ran, and copying a ten-day-old cache now must
not be recorded as a reading taken this morning.

WHY A FILE AND NOT A `gh` CALL. The fast pre-push's window is "test files changed since the last sha
with a GREEN CI run", and the obvious way to get that sha is to ask GitHub. That would put a network
round trip, a rate limit and an authentication failure mode on EVERY push, in a hook whose whole value
is that it costs about a minute. So the hook never asks: it reads this file, and this script -- run by
whoever is already reading CI, which in practice is the coordinator -- is what puts the answer there.

WHERE IT GOES, and this is the part that is easy to get wrong. Under the GIT COMMON DIRECTORY, beside
the verdict status files, at `<common>/genomeos-check/last-green-ci`. Not `<worktree>/.git/...`:
INSIDE A LINKED WORKTREE `.git` IS A REGULAR FILE, not a directory. Measured 2026-10-03 on this
checkout, git 2.50.0, against the live worktree a peer had open:

    $ file .../scratchpad/sweep_wt/.git
    ASCII text                                     # 73 bytes, a `gitdir:` pointer
    $ git -C .../sweep_wt rev-parse --git-dir
    /Users/.../GenomeOS/.git/worktrees/sweep_wt    # per-worktree, nobody else looks here
    $ git -C .../sweep_wt rev-parse --git-common-dir
    /Users/.../GenomeOS/.git                       # the same path the main checkout gets

A path built as `<root>/.git/genomeos-check/last-green-ci` is therefore a path under a FILE when the
reader is in a worktree: it can never exist, so the reader would see "missing" forever and silently
narrow its window every time -- which is the same shape of mistake, a platform assumption nobody
probed, that put three CI runs red on 2026-10-03. `genomeos/verdict.py:status_dir` had already written
this rule down for the status files; this file obeys it and `tests/test_fast_prepush.py` asserts the
two directories are the same one.

IT RECORDS, IT DOES NOT JUDGE. Nothing here decides what "green" means. The caller has read CI and
names the sha, the run id and the run's date; this writes those three plus the time of writing, and
the reader in scripts/fast_prepush.py decides what to do with the age. `--from-cache` is the one
convenience: it takes them out of `data/cache/ci_status.json`, which `scripts/ci_status_cache.py`
already writes from `gh`, so the coordinator's existing CI read can feed this with no second request.
It reads that file from disk and makes no request of its own.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

#: the one name of the file, so the reader and the writer cannot drift apart
FILENAME = "last-green-ci"
#: the daily reading of CI, which `scripts/fast_prepush.py --decide` reads; see the docstring
DAILY_FILENAME = "daily-ci"
#: the one conclusion that is green; every other GitHub conclusion keeps the full gate
GREEN = "success"
#: the directory under the git COMMON dir, shared by the checkout and every worktree of it
DIRNAME = "genomeos-check"
#: where `scripts/ci_status_cache.py` leaves its reading of CI, for `--from-cache`
CI_STATUS_CACHE = ROOT / "data" / "cache" / "ci_status.json"


def git_common_dir(root: Path) -> Path:
    """The git directory shared by this checkout and all of its worktrees, as an absolute path.

    `--git-common-dir`, never `--git-dir` and never `<root>/.git`; the docstring above says why and
    shows the measurement. `--path-format=absolute` because the plain form prints a RELATIVE path
    (".git") from the main checkout, which resolves against the caller's cwd rather than the root.
    """
    out = subprocess.run(
        ["git", "rev-parse", "--path-format=absolute", "--git-common-dir"],
        cwd=root,
        capture_output=True,
        text=True,
        check=True,
    )
    return Path(out.stdout.strip())


def last_green_path(root: Path = ROOT) -> Path:
    return git_common_dir(root) / DIRNAME / FILENAME


def write(root: Path, sha: str, run_id: str | None, date: str | None, now: datetime | None = None) -> dict:
    """Write the record and return it. The sha is resolved against this repository first.

    A sha this checkout does not have is refused rather than written: a record naming an object
    nobody here can diff against would read as "fresh" and then fail at the diff, which is a worse
    failure than no record at all -- the reader's missing-file path is designed and harmless.
    """
    full = subprocess.run(
        ["git", "rev-parse", "--verify", f"{sha}^{{commit}}"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    if full.returncode != 0:
        raise SystemExit(
            f"last_green_ci: {sha!r} is not a commit in this checkout, so nothing was written "
            f"({full.stderr.strip()}). Fetch it first; a record naming an absent object is worse "
            f"than no record."
        )
    record = {
        "sha": full.stdout.strip(),
        "short": subprocess.run(
            ["git", "rev-parse", "--short", full.stdout.strip()], cwd=root, capture_output=True, text=True
        ).stdout.strip(),
        "run_id": str(run_id) if run_id is not None else None,
        "date": date,
        "written": (now or datetime.now(UTC)).isoformat(),
        "written_by": "scripts/last_green_ci.py",
    }
    path = last_green_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(record, indent=1) + "\n")
    tmp.replace(path)  # one atomic rename, so a reader never sees half a record
    return record


def _run_fields(run: dict) -> tuple[str | None, str | None, str | None]:
    """(sha, run id, date) of one run row, in EITHER spelling of the keys.

    `gh run list --json` spells them headSha/databaseId/updatedAt, and `scripts/ci_status_cache.py`
    writes its rows through `_run_row`, which spells them sha/id/created. Until 2026-10-10 this read only
    the first spelling, so against the real cache file it found no sha and refused every time (checked
    on this checkout's cache of 2026-09-30: `latest_success` carries `sha` and no `headSha`).
    """
    sha = run.get("headSha") or run.get("sha")
    rid = run.get("databaseId", run.get("id"))
    date = run.get("updatedAt") or run.get("createdAt") or run.get("created")
    return (str(sha) if sha else None), (str(rid) if rid is not None else None), date


def from_cache(cache: Path = CI_STATUS_CACHE) -> tuple[str, str | None, str | None]:
    """(sha, run id, date) out of `scripts/ci_status_cache.py`'s file. Reads disk; no request."""
    d = json.loads(cache.read_text())
    run = d.get("latest_success")
    sha, rid, date = _run_fields(run) if isinstance(run, dict) else (None, None, None)
    if not sha:
        raise SystemExit(
            f"last_green_ci: {cache} records no latest_success with a sha, so there is no green "
            f"sha to write. Run scripts/ci_status_cache.py, or pass --sha yourself."
        )
    return sha, rid, date


def daily_read_path(root: Path = ROOT) -> Path:
    return git_common_dir(root) / DIRNAME / DAILY_FILENAME


def write_daily(
    root: Path,
    conclusion: str,
    sha: str,
    run_id: str | None,
    date: str | None,
    read_at: datetime | None = None,
    now: datetime | None = None,
) -> dict:
    """Write the daily reading, and on a green one move the last-green record to the same sha.

    The sha is resolved first and a sha this checkout lacks is refused before ANYTHING is written, for
    the reason `write` gives: a reading the gate cannot connect to a commit is worse than none, and
    none is the designed, harmless case (the full gate).
    """
    conclusion = conclusion.strip().lower()
    if not conclusion:
        raise SystemExit("last_green_ci: --daily needs a conclusion (success, failure, cancelled, ...)")
    full = subprocess.run(
        ["git", "rev-parse", "--verify", f"{sha}^{{commit}}"], cwd=root, capture_output=True, text=True
    )
    if full.returncode != 0:
        raise SystemExit(
            f"last_green_ci: {sha!r} is not a commit in this checkout, so no reading was written "
            f"({full.stderr.strip()}). Fetch it first."
        )
    stamp = now or datetime.now(UTC)
    record = {
        "conclusion": conclusion,
        "green": conclusion == GREEN,
        "sha": full.stdout.strip(),
        "short": full.stdout.strip()[:7],
        "run_id": str(run_id) if run_id is not None else None,
        "date": date,
        "read_at": (read_at or stamp).isoformat(),
        "written": stamp.isoformat(),
        "written_by": "scripts/last_green_ci.py --daily",
    }
    if record["green"]:
        write(root, record["sha"], run_id, date, now=read_at or stamp)
    path = daily_read_path(root)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(record, indent=1) + "\n")
    tmp.replace(path)
    return record


def daily_from_cache(cache: Path = CI_STATUS_CACHE) -> tuple[str, str, str | None, str | None, datetime]:
    """(conclusion, sha, run id, date, read_at) of the latest COMPLETED run on dev in the cache.

    `read_at` is the cache's own `written` stamp: when CI was read, not when this copied it.
    """
    d = json.loads(cache.read_text())
    runs = [r for r in d.get("runs") or [] if isinstance(r, dict)]
    on_dev = [
        r
        for r in runs
        if r.get("status") == "completed" and (r.get("branch") or r.get("headBranch")) == "dev"
    ]
    on_dev.sort(key=lambda r: str(r.get("created") or r.get("createdAt") or ""), reverse=True)
    if not on_dev:
        raise SystemExit(
            f"last_green_ci: {cache} lists no completed CI run on dev, so there is no daily reading to "
            f"record. Run scripts/ci_status_cache.py, or pass --daily CONCLUSION --sha SHA yourself."
        )
    run = on_dev[0]
    sha, rid, date = _run_fields(run)
    written = d.get("written")
    if not sha or not run.get("conclusion") or not written:
        raise SystemExit(
            f"last_green_ci: the latest dev run in {cache} lacks a sha, a conclusion or the cache's "
            f"written stamp ({sha!r}, {run.get('conclusion')!r}, {written!r}); nothing recorded."
        )
    read_at = datetime.fromisoformat(str(written).replace("Z", "+00:00"))
    if read_at.tzinfo is None:
        read_at = read_at.replace(tzinfo=UTC)
    return str(run["conclusion"]), sha, rid, date, read_at


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sha", help="the sha CI last went green on")
    ap.add_argument("--run-id", help="the CI run's id, so the record can be traced back to a run")
    ap.add_argument("--date", help="the run's timestamp, ISO 8601")
    ap.add_argument("--from-cache", action="store_true", help="take all three from data/cache/ci_status.json")
    ap.add_argument("--show", action="store_true", help="print the record and change nothing")
    ap.add_argument(
        "--daily",
        metavar="CONCLUSION",
        help="record the daily reading of CI: success, failure, cancelled, ...",
    )
    ap.add_argument(
        "--daily-from-cache", action="store_true", help="the daily reading from data/cache/ci_status.json"
    )
    ap.add_argument("--root", type=Path, default=ROOT)
    args = ap.parse_args(argv)

    if args.show:
        out = {}
        for name, path in (("last_green", last_green_path(args.root)), ("daily", daily_read_path(args.root))):
            out[name] = {"path": str(path), "exists": path.exists()}
            if path.exists():
                out[name]["record"] = json.loads(path.read_text())
        print(json.dumps(out, indent=1))
        return 0 if out["last_green"]["exists"] else 1

    if args.daily_from_cache:
        conclusion, sha, run_id, date, read_at = daily_from_cache()
        record = write_daily(args.root, conclusion, sha, run_id, date, read_at=read_at)
        print(json.dumps({"path": str(daily_read_path(args.root)), "record": record}, indent=1))
        return 0
    if args.daily:
        if not args.sha:
            ap.error("--daily needs --sha (and --run-id and --date when known)")
        record = write_daily(args.root, args.daily, args.sha, args.run_id, args.date)
        print(json.dumps({"path": str(daily_read_path(args.root)), "record": record}, indent=1))
        return 0

    if args.from_cache:
        sha, run_id, date = from_cache()
    elif args.sha:
        sha, run_id, date = args.sha, args.run_id, args.date
    else:
        ap.error("pass --sha (with --run-id and --date), or --from-cache, or --show")
        return 64

    record = write(args.root, sha, run_id, date)
    print(json.dumps({"path": str(last_green_path(args.root)), "record": record}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
