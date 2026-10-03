#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Record the sha CI last went GREEN on, so the fast pre-push can widen its window with no network call.

    uv run python scripts/last_green_ci.py --sha 747b7ff --run-id 18123456789 --date 2026-10-01T09:14:00Z
    uv run python scripts/last_green_ci.py --from-cache      # from data/cache/ci_status.json's latest_success
    uv run python scripts/last_green_ci.py --show            # print what is recorded, change nothing

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


def from_cache(cache: Path = CI_STATUS_CACHE) -> tuple[str, str | None, str | None]:
    """(sha, run id, date) out of `scripts/ci_status_cache.py`'s file. Reads disk; no request."""
    d = json.loads(cache.read_text())
    run = d.get("latest_success")
    if not isinstance(run, dict) or not run.get("headSha"):
        raise SystemExit(
            f"last_green_ci: {cache} records no latest_success with a headSha, so there is no green "
            f"sha to write. Run scripts/ci_status_cache.py, or pass --sha yourself."
        )
    rid = run.get("databaseId")
    return (
        str(run["headSha"]),
        (str(rid) if rid is not None else None),
        run.get("updatedAt") or run.get("createdAt"),
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sha", help="the sha CI last went green on")
    ap.add_argument("--run-id", help="the CI run's id, so the record can be traced back to a run")
    ap.add_argument("--date", help="the run's timestamp, ISO 8601")
    ap.add_argument("--from-cache", action="store_true", help="take all three from data/cache/ci_status.json")
    ap.add_argument("--show", action="store_true", help="print the record and change nothing")
    ap.add_argument("--root", type=Path, default=ROOT)
    args = ap.parse_args(argv)

    if args.show:
        path = last_green_path(args.root)
        if not path.exists():
            print(json.dumps({"path": str(path), "exists": False}, indent=1))
            return 1
        print(
            json.dumps({"path": str(path), "exists": True, "record": json.loads(path.read_text())}, indent=1)
        )
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
