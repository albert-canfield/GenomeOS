# SPDX-License-Identifier: AGPL-3.0-or-later
"""Remove a git worktree that holds read-only store clones, without `rm` and without a recursive chmod.

    uv run python scripts/remove_worktree.py WORKTREE            # check, unlock the recorded clones, remove
    uv run python scripts/remove_worktree.py WORKTREE --check    # run every refusal check, change nothing
    uv run python scripts/remove_worktree.py WORKTREE --force    # pass --force to git, for untracked scratch

THE INCIDENT. `scripts/link_stores_read_only.py` clones `data/reference`, `data/knowledge` and
`data/cache` into a verdict worktree at mode 0444/0555, so a verification run cannot modify the inputs
it verifies against. `git worktree remove` then cannot delete the worktree, and on 2026-10-02 a lane
therefore fell back to `chmod -R u+w` over the whole worktree followed by `rm -rf`. It did that
carefully -- it checked the clones were distinct inodes from the real stores first and re-checked the
stores afterwards, and it disclosed all of it -- but both halves are the hazard. A recursive chmod
restores write on paths nobody recorded, and `rm -rf` on a path that may contain a link or a clone into
a store is the one operation that can destroy the machine's only copy of 11.2 GB that is not in git.
The rule in force since is: NO WORKTREE IS REMOVED BY `rm`. This is the tool behind that rule.

MEASURED, 2026-10-03, git 2.50.0, a planted repo with genuine `cp -c` clones locked to 0444/0555 --
because a tool whose reason for existing is asserted rather than measured is a tool nobody can check:

  git worktree remove <wt>            exit 255   error: failed to delete '<wt>': Permission denied
  git worktree remove --force <wt>    exit 128   fatal: '<wt>' is not a working tree

So the premise holds, and the second line is the part that was not expected and is worse than the
premise. The FAILED removal is NOT atomic: before it reached the read-only `data/` it had already
deleted the worktree's tracked files and its `.git` file and DEREGISTERED the worktree, so the retry
cannot find a working tree at all. What is left on disk is an orphan directory that `git worktree
remove` can no longer touch by any flag -- which is exactly the dead end that leaves a lane with no git
route and reaching for `rm -rf`. Hence the order here is fixed and is the whole design: UNLOCK FIRST,
then hand a worktree git can actually delete to git. `tests/test_remove_worktree.py::test_d_*` keeps
that measurement running rather than quoting it.

WHAT IT REFUSES, and each refusal exists because its absence is how a store gets destroyed:

1. A path git does not list as a worktree. Taken from `git worktree list --porcelain` and nowhere else.
   An orphan left by a half-done removal lands here too: it is REPORTED to its owner, never removed,
   because this tool has no way to tell an orphan from a directory that merely looks like one.
2. The main worktree. `git worktree remove` refuses it anyway; naming it here means the refusal says so.
3. A recorded path whose inode is a REAL STORE's. If the worktree holds the store itself rather than a
   clone of it, restoring write reaches the machine's only copy and removing the worktree removes the
   store. Every recorded path is compared against EVERY real store, not just its own, and the check is
   `(st_dev, st_ino)` from `os.stat`, which follows, so an indirection cannot hide the identity.
4. A recorded path that is a symlink. A symlink carries its TARGET's mode, so a chmod through one
   reaches whatever it points at; this is the defect `link_stores_read_only` was written to remove and
   it must not come back through the remover.
5. A recorded path that contains, or sits inside, a real store. Inode equality does not catch a
   recorded SUBDIRECTORY of a store, and containment does. SAID PLAINLY: no test reaches this check.
   A mutation run on 2026-10-03 took it out and `tests/test_remove_worktree.py` stayed green, and that
   is not a gap in the tests but a fact about directories -- APFS has no directory hard links, so for
   two directory paths "same inode" and "same path once resolved" are one fact, and check 3 fires
   first on every case this one would catch. The reachable remainder is a recorded path INSIDE a store
   while also inside the worktree, which needs a store to live inside the worktree, which needs the
   worktree to contain the repository root, which check 2 refuses. It is kept because it costs one
   `realpath` and fails closed, and it is named here so nobody reads it as a measured guard.
6. A recorded path that is not inside the worktree, or whose live inode is not the one recorded. Either
   means the record no longer describes what is there, and a record that is not trusted is not acted on.

It touches nothing but the paths `link_stores_read_only` RECORDED it cloned, and it prints each real
store's inode and mode before and after, so "the stores are untouched" is a reading and not a promise.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import stat
import subprocess
import sys
from pathlib import Path
from typing import Any

_spec = importlib.util.spec_from_file_location(
    "lsro", Path(__file__).resolve().parent / "link_stores_read_only.py"
)
lsro = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(lsro)


class RefusedError(Exception):
    """A check said no. Every message names the path, because a refusal nobody can act on is a stop."""


def registered_worktrees(root: Path) -> list[Path]:
    """The worktrees git lists, resolved. The ONLY source of what may be removed."""
    out = subprocess.run(["git", "worktree", "list", "--porcelain"], cwd=root, capture_output=True, text=True)
    if out.returncode != 0:
        raise RefusedError(f"`git worktree list` failed in {root}: {out.stderr.strip()}")
    return [
        Path(line[len("worktree ") :]).resolve()
        for line in out.stdout.splitlines()
        if line.startswith("worktree ")
    ]


def real_stores(root: Path) -> dict[str, Path]:
    """The machine's one copy of each store that is actually on this machine."""
    return {n: (root / "data" / n) for n in lsro.STORES if (root / "data" / n).is_dir()}


def store_state(root: Path) -> dict[str, dict[str, Any]]:
    """Each real store's inode and mode, for the before/after reading this tool prints itself."""
    state = {}
    for name, p in real_stores(root).items():
        st = p.stat()
        state[name] = {
            "path": str(p),
            "inode": st.st_ino,
            "device": st.st_dev,
            "mode": oct(stat.S_IMODE(st.st_mode)),
        }
    return state


def _inside(child: Path, parent: Path) -> bool:
    """Whether `child` is `parent` or below it, with every symlink on both paths resolved first."""
    c, p = os.path.realpath(child), os.path.realpath(parent)
    return c == p or c.startswith(p.rstrip(os.sep) + os.sep)


def check(root: Path, worktree: Path) -> dict[str, Any]:
    """Every refusal check, touching nothing. Raises `RefusedError` naming the offending path."""
    worktree = worktree.resolve()
    registered = registered_worktrees(root)
    if worktree == root.resolve():
        raise RefusedError(f"{worktree} is the MAIN worktree, not one of git's removable worktrees")
    if worktree not in registered:
        listed = "\n  ".join(str(p) for p in registered) or "(none)"
        raise RefusedError(
            f"{worktree} is not in `git worktree list`, so this tool will not remove it. If a failed "
            f"removal left it behind, it is an orphan: report it to its owner. git lists:\n  {listed}"
        )

    stores = real_stores(root)
    record = lsro.read_manifest(worktree, root)
    plan: list[dict[str, Any]] = []
    for entry in record["clones"]:
        path = Path(entry["path"])
        # The symlink check comes FIRST, before containment, so that a link out of the worktree is
        # diagnosed by its mechanism rather than by where it happens to land. `_inside` resolves
        # symlinks, so it would catch this one too -- but it would say "not inside the worktree", which
        # does not tell the reader that the hazard is a mode borrowed from the target.
        if path.is_symlink():
            raise RefusedError(
                f"recorded clone {path} is a SYMLINK: its mode is its target's, so restoring write "
                f"through it would reach {os.path.realpath(path)}"
            )
        if not path.exists():
            if not _inside(path, worktree):
                raise RefusedError(
                    f"recorded clone {path} is NOT inside {worktree}; the record is not trusted"
                )
            plan.append({"store": entry.get("store"), "path": str(path), "state": "recorded but absent"})
            continue
        # The INODE check runs before the weaker complaints about where the path sits, so that the most
        # dangerous condition is always the one reported. A record pointing at a real store is both
        # "the same inode as the store" and "outside the worktree"; the first is what a reader needs to
        # be told, and until this was reordered the second arrived first and buried it.
        live = os.stat(path)
        for _name, real in stores.items():
            rst = os.stat(real)
            if (live.st_dev, live.st_ino) == (rst.st_dev, rst.st_ino):
                raise RefusedError(
                    f"REFUSING {path}: it is the SAME INODE as the real store {real} "
                    f"(device {rst.st_dev}, inode {rst.st_ino}). This worktree holds the store itself "
                    f"rather than a clone of it, so removing it would remove the store."
                )
            if _inside(path, real) or _inside(real, path):
                raise RefusedError(f"REFUSING {path}: it contains, or sits inside, the real store {real}")
        if not _inside(path, worktree):
            raise RefusedError(f"recorded clone {path} is NOT inside {worktree}; the record is not trusted")
        if "inode" in entry and (live.st_dev, live.st_ino) != (entry["device"], entry["inode"]):
            raise RefusedError(
                f"REFUSING {path}: the linker recorded device {entry['device']} inode {entry['inode']} "
                f"and what is there now is device {live.st_dev} inode {live.st_ino}, so the record no "
                f"longer describes this path"
            )
        plan.append(
            {
                "store": entry.get("store"),
                "path": str(path),
                "state": "a clone, distinct inode from the real store",
                "inode": live.st_ino,
                "real_store_inode": os.stat(stores[entry["store"]]).st_ino
                if entry.get("store") in stores
                else None,
            }
        )
    return {
        "worktree": str(worktree),
        "manifest": str(lsro.manifest_path(worktree, root)),
        "recorded_clones": len(record["clones"]),
        "plan": plan,
        "stores_before": store_state(root),
    }


def remove(root: Path, worktree: Path, force: bool = False, check_only: bool = False) -> dict[str, Any]:
    """Check, unlock only the recorded clones, then let git remove the worktree. Never `rm`."""
    out = check(root, worktree)
    out["checked_only"] = check_only
    if check_only:
        out["stores_after"] = out["stores_before"]
        out["ok"] = True
        return out

    unlocked = []
    for item in out["plan"]:
        if item["state"].startswith("a clone"):
            unlocked.append({"path": item["path"], "entries": lsro.unlock(Path(item["path"]))})
    out["unlocked"] = unlocked

    argv = ["git", "worktree", "remove", *(["--force"] if force else []), str(Path(worktree).resolve())]
    r = subprocess.run(argv, cwd=root, capture_output=True, text=True)
    out["git_command"] = " ".join(argv)
    out["git_exit"] = r.returncode
    out["git_stderr"] = r.stderr.strip()
    out["git_stdout"] = r.stdout.strip()
    out["removed"] = r.returncode == 0 and not Path(worktree).resolve().exists()
    out["stores_after"] = store_state(root)
    out["stores_unchanged"] = out["stores_after"] == out["stores_before"]
    if out["removed"]:
        # the record describes a worktree that no longer exists; one file, unlinked, never a tree walk
        Path(out["manifest"]).unlink(missing_ok=True)
        out["manifest_removed"] = True
    else:
        # git refused. Report its exact words and STOP: escalating from here is what this tool exists
        # to prevent, and the escalation available -- `rm -rf` -- is the one that can destroy a store.
        out["refused_by_git"] = r.stderr.strip() or f"git exited {r.returncode} without a message"
    out["ok"] = bool(out["removed"]) and bool(out["stores_unchanged"])
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("worktree", type=Path)
    ap.add_argument("--force", action="store_true", help="pass --force to git, for untracked scratch files")
    ap.add_argument("--check", action="store_true", help="run every check and change nothing")
    args = ap.parse_args(argv)
    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    try:
        r = remove(root, args.worktree, force=args.force, check_only=args.check)
    except RefusedError as e:
        json.dump({"ok": False, "refused": str(e), "stores": store_state(root)}, sys.stdout, indent=1)
        print()
        return 2
    json.dump(r, sys.stdout, indent=1)
    print()
    return 0 if r["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
