# SPDX-License-Identifier: AGPL-3.0-or-later
"""Find the commits a branch lost: reflog transitions where the old tip is not an ancestor of the new.

    uv run python scripts/check_history.py              # dev
    uv run python scripts/check_history.py --ref main --json

Several sessions move the same branch. The documented way is `git update-ref refs/heads/dev NEW OLD`,
which fails if someone else moved it first. Drop the third argument — as one session did on every
commit of 2026-09-17 without noticing — and a peer commit landing in the seconds between `rev-parse
HEAD` and `update-ref` is overwritten. There is no conflict and no error. The commit is simply not in
the history, and nobody finds out until they go looking for work that should be there.

This is the check that works afterwards, and it is the only one: for each consecutive pair in the
reflog, ask whether the old tip is an ancestor of the new one. If it is not, the branch moved
sideways and whatever was only on the old tip left the history.

**Not every hit is a loss.** An amend or a rebase is also a sideways move, and the honest signal is
weak: a pair whose subjects match is almost always the same work under a new hash, and a pair whose
subjects differ is the shape that loses somebody's commit. Both are reported, labelled, and neither is
called a verdict — `relanded` says whether a commit with the old subject is on the branch today, which
is what turns "this happened" into "this was recovered".

On this repository when it was written: 598 transitions on dev, 6 sideways, 3 of them with different
subjects (2026-09-11 and two on 2026-09-12), and all 6 re-landed. None on 2026-09-17.
"""

from __future__ import annotations

import argparse
import json
import subprocess


def git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, check=False).stdout.strip()


def is_ancestor(old: str, new: str) -> bool:
    return (
        subprocess.run(["git", "merge-base", "--is-ancestor", old, new], capture_output=True).returncode == 0
    )


def subject(commit: str) -> str:
    return git("log", "--format=%s", "-1", commit)


def when(commit: str) -> str:
    return git("log", "--format=%ad", "--date=short", "-1", commit)


def relanded(ref: str, subj: str, old: str) -> str | None:
    """A commit on the branch today carrying the lost tip's subject, if there is one."""
    if not subj:
        return None
    for line in git("log", "--format=%H %s", ref).splitlines():
        h, _, s = line.partition(" ")
        if s == subj and h != old:
            return h[:7]
    return None


def sideways(ref: str) -> list[dict]:
    tips = [t for t in git("reflog", "show", ref, "--format=%H").splitlines() if t]
    out = []
    for new, old in zip(tips, tips[1:], strict=False):  # reflog is newest first
        if new == old or is_ancestor(old, new):
            continue
        s_old, s_new = subject(old), subject(new)
        out.append(
            {
                "old": old[:7],
                "new": new[:7],
                "date": when(old),
                "old_subject": s_old,
                "new_subject": s_new,
                # the weak signal, stated as weak: same subject is an amend, different is the shape
                # that overwrites somebody's commit
                "looks_like": "amend or rewrite" if s_old == s_new else "a ref moved over a commit",
                "relanded_as": relanded(ref, s_old, old),
            }
        )
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ref", default="dev")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    rows = sideways(args.ref)
    total = max(0, len(git("reflog", "show", args.ref, "--format=%H").splitlines()) - 1)
    if args.json:
        print(json.dumps({"ref": args.ref, "transitions": total, "sideways": rows}, indent=1))
        return 0

    print(f"{args.ref}: {total} reflog transitions, {len(rows)} not fast-forward")
    lost = 0
    for r in rows:
        flag = "!" if r["looks_like"] != "amend or rewrite" else " "
        print(f" {flag} {r['date']}  {r['old']} -> {r['new']}  ({r['looks_like']})")
        print(f"     lost: {r['old_subject'][:88]}")
        if r["looks_like"] != "amend or rewrite":
            lost += 1
            print(f"     kept: {r['new_subject'][:88]}")
        print(f"     re-landed as: {r['relanded_as'] or 'NOT FOUND ON THIS BRANCH'}")
    if rows:
        missing = [r for r in rows if not r["relanded_as"]]
        print(
            f"\n{lost} moved a ref over a different commit; "
            f"{len(missing)} have no commit with that subject on {args.ref} today."
        )
        if missing:
            print("Those are the ones to look at: the work may only exist in the reflog.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
