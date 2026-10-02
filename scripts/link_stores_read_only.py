# SPDX-License-Identifier: AGPL-3.0-or-later
"""Put the git-ignored data stores into a verdict worktree READ-ONLY, and take the lock off again.

    uv run python scripts/link_stores_read_only.py WORKTREE            # link
    uv run python scripts/link_stores_read_only.py WORKTREE --unlock   # before removing the worktree

The amended acceptance rule (2026-10-02, the supervisor's ruling). The project's rule was "a status-file
verdict from a worktree of the committed tree", and it was unachievable as stated: `data/cache`,
`data/reference` and `data/knowledge` are git-ignored machine-local data, absent from a worktree by
definition, and 21 tests therefore FAILED in one -- 16 in `tests/test_context_evidence.py`, every one
because the cell-to-biosample mapping is read from a table on this machine and is never guessed from
names. A verdict decided by the environment rather than by the code is the defect class the project
spent the day removing, reappearing inside the rule adopted to escape it.

The amended rule: the verdict is a worktree of the committed tree WITH the git-ignored stores linked
READ-ONLY, so a test that needs local data runs there; and where the stores are genuinely absent -- CI --
such a test skips by name under `needs_local_data` and never fails (`tests/local_data.py`).

WHY A CLONE AND NOT A SYMLINK, which is what `scripts/pre-push.sh` did until today. A symlink carries
its TARGET's mode. The stores' directories are mode 0755, so a write through the link reached the only
copy of the store on this machine -- measured on 2026-10-02:
`open("<worktree>/data/reference/HG002_chr1.vcf.gz", "r+b")` through the link SUCCEEDED, and a probe file
created through it appeared in the real store. A verification run that can modify the inputs it verifies
against is the strongest form of the failure class being removed, because it can corrupt data rather
than merely misreport. `cp -c` asks APFS for a clone: no bytes are copied, both sides share them until
one is written, and the clone is a SEPARATE INODE -- so taking the write bits off afterwards takes them
off the worktree's side only and never off the machine's one copy. It is the mechanism
`scripts/manifest_rebuild.py: link_read_only` already uses per file; this is the same thing per store.

MEASURED COST, because it runs on every push: `data/knowledge`, 6.5 GB and 45,326 files, clones in 7.3 s
and locks in 1.2 s, and free disk did not move. A clone that silently became a real copy would show up
as both -- minutes, and gigabytes -- which is why the figures are printed rather than assumed.

WHAT IS OUT OF SCOPE AND WHY, stated rather than left to be noticed: `data/results` is NOT linked.
It is ignored by PATTERN with `!` re-includes and 848 of its files are committed, so it cannot be cloned
over a worktree's own tracked copies without destroying the tree identity the verdict rests on. Linking
only its ignored members is a per-path job (`manifest_rebuild.link_machine_local_inputs` does exactly
that for one result's declared inputs) and no test among the 21 needs one, so it is left undone and said.

TREE IDENTITY. A verdict must name the committed tree, `tree_begin == tree_end ==` it. A clone is a real
directory where a symlink was, so this is not self-evident and is not assumed: `--check-tree` prints the
worktree's hash before and after linking, and `tests/test_link_stores_read_only.py` asserts they are
equal. `.gitignore` carries the three stores WITHOUT a trailing slash, which matches a directory as well
as a file, and the comment there records that the spelling was chosen for exactly this reason.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any

#: The git-ignored stores linked whole. `data/results` is deliberately not here; see the docstring.
STORES = ("reference", "knowledge", "cache")

#: Directories keep traverse and list, and lose create. Files lose write. Both are the kernel's to
#: enforce, which is what makes this stronger than the audit hook in `scripts/rebuild_write_guard.py`:
#: a C library that opens a file itself is invisible to a hook and is refused by a mode.
DIR_MODE = 0o555
FILE_MODE = 0o444


def tree_hash(where: Path) -> str:
    """The working-tree hash the verdict names, as `genomeos.verdict` computes it."""
    r = subprocess.run(
        [sys.executable, "-m", "genomeos.verdict", "tree"], cwd=where, capture_output=True, text=True
    )
    return r.stdout.strip() if r.returncode == 0 else f"unavailable: {r.stderr.strip()[-200:]}"


def lock(path: Path) -> int:
    """Take the write bits off `path` and everything under it. Returns how many entries were touched."""
    n = 0
    for dirpath, dirnames, filenames in os.walk(path, topdown=False, followlinks=False):
        for name in filenames:
            p = Path(dirpath) / name
            if not p.is_symlink():
                os.chmod(p, FILE_MODE)
                n += 1
        for name in dirnames:
            os.chmod(Path(dirpath) / name, DIR_MODE)
            n += 1
    os.chmod(path, DIR_MODE)
    return n + 1


def unlock(path: Path) -> int:
    """Put the write bits back, so removing the worktree cannot be refused by a read-only directory.

    Every entry under here is a clone of its own, so this touches nothing the machine keeps.
    """
    n = 0
    for dirpath, dirnames, filenames in os.walk(path, topdown=True, followlinks=False):
        os.chmod(dirpath, 0o755)
        for name in filenames:
            p = Path(dirpath) / name
            if not p.is_symlink():
                os.chmod(p, 0o644)
                n += 1
        n += len(dirnames)
    return n


def clone(src: Path, dst: Path) -> tuple[bool, float]:
    """(whether APFS cloned it, seconds). A fallback to a real copy is REPORTED, never silent: a clone
    costs nothing and a copy costs the store's size, and the two read the same in a report otherwise."""
    t = time.monotonic()
    cloned = subprocess.run(["cp", "-c", "-R", str(src), str(dst)], capture_output=True).returncode == 0
    if not cloned:
        shutil.copytree(src, dst, symlinks=True)
    return cloned, round(time.monotonic() - t, 2)


def free_bytes(path: Path) -> int:
    """Free bytes on the filesystem `path` is on, or -1 when it cannot be asked. Never raises: a disk
    reading is a figure in a report and must not be able to stop the linking it is reporting on."""
    for candidate in (path, *path.parents):
        try:
            st = os.statvfs(candidate)
        except OSError:
            continue
        return st.f_bavail * st.f_frsize
    return -1


def link(root: Path, worktree: Path, check_tree: bool = False) -> dict[str, Any]:
    """Clone each store into `worktree` read-only. Idempotent: a store already there is left alone."""
    out: dict[str, Any] = {"worktree": str(worktree), "stores": [], "failures": []}
    if check_tree:
        out["tree_before"] = tree_hash(worktree)
    # the directory is made BEFORE the disk is measured: a worktree that does not exist yet cannot be
    # asked how much room it has, and that reading is a figure in a report, not a precondition
    (worktree / "data").mkdir(parents=True, exist_ok=True)
    before = free_bytes(worktree)
    for name in STORES:
        src, dst = root / "data" / name, worktree / "data" / name
        if not src.is_dir():
            out["stores"].append({"store": name, "linked": False, "why": "not on this machine"})
            continue
        if dst.exists() or dst.is_symlink():
            out["stores"].append({"store": name, "linked": False, "why": "already in the worktree"})
            continue
        try:
            cloned, seconds = clone(src, dst)
            entries = lock(dst)
        except OSError as e:
            out["failures"].append(f"{name}: {e}")
            continue
        probe = _write_is_refused(dst)
        if not probe:
            out["failures"].append(
                f"{name}: a write into the linked store was NOT refused after locking it, so this "
                f"worktree could modify the inputs it verifies against"
            )
        out["stores"].append(
            {
                "store": name,
                "linked": True,
                "cloned": cloned,
                "copied_not_cloned": not cloned,
                "entries_locked": entries,
                "seconds": seconds,
                "a_write_into_it_is_refused": probe,
            }
        )
    out["free_bytes_before"] = before
    out["free_bytes_after"] = free_bytes(worktree)
    out["free_bytes_spent"] = before - free_bytes(worktree)
    if check_tree:
        out["tree_after"] = tree_hash(worktree)
        out["tree_unchanged"] = out["tree_before"] == out["tree_after"]
    out["ok"] = not out["failures"]
    return out


def _write_is_refused(store: Path) -> bool:
    """Whether creating a file in the locked store really fails. Asked of the kernel, not assumed."""
    probe = store / ".link_stores_read_only_probe"
    try:
        with open(probe, "w"):
            pass
    except OSError:
        return True
    probe.unlink(missing_ok=True)
    return False


def unlink(worktree: Path) -> dict[str, Any]:
    """Put the write bits back on every linked store, so `git worktree remove` cannot be refused."""
    out: dict[str, Any] = {"worktree": str(worktree), "unlocked": []}
    for name in STORES:
        dst = worktree / "data" / name
        if dst.is_dir() and not dst.is_symlink():
            out["unlocked"].append({"store": name, "entries": unlock(dst)})
    out["ok"] = True
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("worktree", type=Path)
    ap.add_argument("--unlock", action="store_true", help="put the write bits back before removal")
    ap.add_argument("--check-tree", action="store_true", help="print the worktree hash before and after")
    args = ap.parse_args(argv)
    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    r = (
        unlink(args.worktree.resolve())
        if args.unlock
        else link(root, args.worktree.resolve(), args.check_tree)
    )
    json.dump(r, sys.stdout, indent=1)
    print()
    return 0 if r["ok"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
