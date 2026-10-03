# SPDX-License-Identifier: AGPL-3.0-or-later
"""Put the git-ignored data stores into a verdict worktree READ-ONLY, and take the lock off again.

    uv run python scripts/link_stores_read_only.py WORKTREE            # link
    uv run python scripts/link_stores_read_only.py WORKTREE --unlock   # before removing the worktree

TO REMOVE SUCH A WORKTREE, use `scripts/remove_worktree.py` and not `rm`. `git worktree remove` cannot
delete a worktree whose clones are at 0444 -- MEASURED 2026-10-03, git 2.50.0: exit 255,
`error: failed to delete '<wt>': Permission denied` -- and the failed attempt deregisters the worktree
first, so `--force` then answers `fatal: is not a working tree` and what is on disk is an orphan git
can no longer touch. That dead end is what made a lane fall back to `chmod -R u+w` and `rm -rf` on
2026-10-02. The remover unlocks only the paths THIS SCRIPT RECORDED it cloned (see MANIFEST_DIR) and
refuses any whose inode is a real store's.

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
import hashlib
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

#: Where the linker RECORDS what it cloned, and the ONLY record `scripts/remove_worktree.py` trusts.
#: Added 2026-10-03, after a lane met the hole this leaves: `git worktree remove` cannot delete a
#: worktree whose clones are at 0444, so the lane restored write with `chmod -R u+w` over the whole
#: worktree and then `rm -rf`. Both halves are the hazard. A recursive chmod restores write on paths
#: nobody recorded, and `rm -rf` on a path that may contain a link or a clone into a store is the one
#: operation that can destroy the machine's only copy. A remover that trusts a RECORD instead touches
#: exactly the paths the linker says it made, and can check each one's inode against the real store's
#: before it touches it. The record lives in the repository's COMMON git directory and NOT inside the
#: worktree: a file inside the worktree would be untracked content, and `tree_before == tree_after` is
#: precisely what this script exists to preserve. It is keyed by the resolved worktree path, so two
#: worktrees never share a record and a stale record can never be read as another worktree's.
MANIFEST_DIR = "genomeos_linked_stores"

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


def manifest_path(worktree: Path, root: Path) -> Path:
    """The record file for `worktree`, under the common git directory. Pure: it creates nothing.

    The common directory is asked of the WORKTREE first and of `root` second. Either answers the same
    thing for a real worktree, and asking the worktree first is what makes this work when `root` is not
    a git repository at all -- which `link` is called with (a root with no stores on this machine is
    named rather than failed), so a record that could only be placed from `root` would turn that into a
    crash. Raises when neither is a repository: there is then no worktree for anyone to remove.
    """
    last: Exception | None = None
    for cwd in (worktree, root):
        try:
            common = subprocess.check_output(
                ["git", "rev-parse", "--path-format=absolute", "--git-common-dir"],
                cwd=cwd,
                text=True,
                stderr=subprocess.DEVNULL,
            ).strip()
        except (OSError, subprocess.CalledProcessError) as e:
            last = e
            continue
        key = hashlib.sha256(str(worktree.resolve()).encode()).hexdigest()[:16]
        return Path(common) / MANIFEST_DIR / f"{key}.json"
    raise RuntimeError(f"no git common directory from {worktree} or {root}: {last}")


def read_manifest(worktree: Path, root: Path) -> dict[str, Any]:
    """What the linker recorded for `worktree`, or an empty record. A malformed file reads as EMPTY
    rather than raising: the remover's refusals must come from its checks, and an unreadable record
    means nothing is known to be a clone, which is the safe reading, not a reason to guess."""
    try:
        r = json.loads(manifest_path(worktree, root).read_text())
    except (OSError, ValueError, RuntimeError):
        return {"worktree": str(worktree), "clones": []}
    return (
        r
        if isinstance(r, dict) and isinstance(r.get("clones"), list)
        else {
            "worktree": str(worktree),
            "clones": [],
        }
    )


def write_manifest(worktree: Path, root: Path, clones: list[dict[str, Any]]) -> Path:
    """Record `clones` for `worktree`, MERGED over whatever an earlier run recorded for the same
    worktree: `link` is idempotent and leaves a store it finds already there alone, so a second run
    must not drop the first run's record of it. One entry per store, this run's winning."""
    p = manifest_path(worktree, root)
    p.parent.mkdir(parents=True, exist_ok=True)
    merged = {c["store"]: c for c in read_manifest(worktree, root)["clones"] if "store" in c}
    merged.update({c["store"]: c for c in clones})
    record = {
        "worktree": str(worktree),
        "root": str(root),
        "recorded_at": time.time(),
        "clones": [merged[k] for k in sorted(merged)],
    }
    p.write_text(json.dumps(record, indent=1) + "\n")
    return p


def _recorded(store: str, src: Path, dst: Path, cloned: bool, entries: int) -> dict[str, Any]:
    """One clone's record. The inodes of BOTH sides go in, so a remover can refuse a record whose
    worktree side has since been replaced by the real store rather than trusting the path string."""
    s, d = src.stat(), dst.stat()
    return {
        "store": store,
        "path": str(dst),
        "source": str(src),
        "inode": d.st_ino,
        "device": d.st_dev,
        "source_inode": s.st_ino,
        "source_device": s.st_dev,
        "cloned": cloned,
        "entries_locked": entries,
    }


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
    recorded: list[dict[str, Any]] = []
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
        recorded.append(_recorded(name, src, dst, cloned, entries))
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
    # The record is written LAST, so a run that dies midway leaves a locked clone that NOTHING records.
    # That direction is the safe one and is the reason the order is this way round rather than the
    # other: `scripts/remove_worktree.py` unlocks only recorded paths, so an unrecorded clone makes it
    # refuse or makes git refuse, and either way a human is told. The opposite order would leave a
    # record of a clone that does not exist, which is a record inviting a chmod of whatever is there.
    try:
        out["manifest"] = str(write_manifest(worktree, root, recorded))
    except (OSError, RuntimeError) as e:
        out["manifest"] = None
        out["manifest_why"] = str(e)
        # A clone that was locked and NOT recorded is a trap: `scripts/remove_worktree.py` unlocks only
        # recorded paths, so nobody can unlock it and nobody is told why. That is a failure. Nothing
        # cloned means nothing to record, and no worktree anyone needs to remove, so that is not.
        if recorded:
            out["failures"].append(f"the record of {len(recorded)} locked clone(s) could not be written: {e}")
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
