# SPDX-License-Identifier: AGPL-3.0-or-later
"""Refuse, by mechanism, any write into the data stores a rebuild verifies a result against.

The defect this closes, measured on 2026-10-02 and not reasoned: `scripts/manifest_rebuild.py` links the
machine-local inputs under `data/results` as APFS clones at mode 0444, so a writer gets PermissionError
through them -- but it SYMLINKS `data/reference`, `data/knowledge` and `data/cache` into the worktree, and
a symlink carries the target's mode. `open(<worktree>/data/reference/HG002_chr1.vcf.gz, "r+b")` through
such a link succeeds, the real store's directories are mode 0755, and a file created through the link is
created in the one and only copy on this machine. That was true of every rebuild this tool has run. A
verification tool that can modify the inputs it checks a result against is the worst form of the failure
class the project spent the day removing, because it can corrupt data rather than merely misreport.

The mechanism is an audit hook, the same instrument `genomeos.manifest`'s tracer uses and installed the
same way: once per process, left in place because an audit hook cannot be removed, and inert until it is
armed. Armed, it raises PermissionError before a byte moves.

THE AUDIT BEHIND THIS FILE, AND THE WINDOW IT DOES NOT COVER. Because the tool COULD have written to
those stores, it may already have, and a guard that only prevents future writes leaves that unknown. So
the three stores were audited on 2026-10-02: 45,734 files walked, 1,028 modified at or after 2026-09-28
(830 on 09-28, 103 on 09-29, 1 on 09-30, 13 on 10-01, 81 on 10-02) across 20 top-level directories, and
EVERY one of those directories has a named committed writer. The two with no writer under `scripts/` or
`genomeos/` are `data/cache/review_screenshots` (five PNGs of the web UI, named in docs/ROADMAP.md) and
`data/cache/paired_enhancer` (two downloaded publication supplements); two `.DS_Store` files are Finder's.
The sharpest figure: `data/reference`, 0 of 179 files modified in the entire window -- and that is the
store the probe demonstrated writable. The thing that could have been written was not written.

The standing caveat, which is not a footnote. "Every modification is explained by a known writer" is
ESTABLISHED. "No rebuild ever wrote to an input store" is NOT, and cannot be established from mtimes at
all. The reason is exact: several of those directories are caches that the results' OWN writers populate,
and a rebuild runs the result's own command -- so a cache write by the original run and a cache write by
a rebuild of that run are THE SAME WRITER AND THE SAME PATH, and an mtime cannot tell them apart.
`data/cache/holdout`, 142 files, is the concrete instance. Three further limits belong with it: an mtime
can be rewritten; a write that restored identical bytes leaves no trace in mtime or in content; and
attribution by searching for a path in a committed writer shows that such a writer EXISTS, not that it
rather than a rebuild made that file. Everything outside those three directories is out of scope,
`data/results` and `data/organisms` included.

So the honest window: a rebuild verdict taken from the commit that added this file onward is covered by
the mechanism; an earlier one rests on an assumption that is reasonable and unevidenced. This guard makes
it evidenced from now on, NOT retrospectively, and `manifest_rebuild` prints that window on every report
rather than leaving it here.

WHAT IS REFUSED, and this is NARROWER than "any write outside the worktree" on purpose. A rebuild
legitimately writes outside its worktree all the time -- a temporary file, the uv cache, a library's
config directory -- so refusing every outside write would break the rebuilds rather than protect the
stores. What is refused is a write whose path RESOLVES INTO one of the protected stores, by real path, so
the symlink route and a direct absolute path to the real store are refused alike.

WHICH EVENTS, by name (`WATCHED`): a named set, not "every write". `open` in a writing mode or with a
writing flag, and the filesystem calls that change a path without opening it -- remove, rename, mkdir,
rmdir, truncate, chmod, utime, link, symlink, and shutil's copyfile, move and rmtree. An event not on
that list is not watched, and saying so is the point of writing the list down.

WHAT THIS IS NOT. An audit hook sees `io.open`, `os.open` and `pathlib` alike, which is why it is the
instrument and not a monkeypatch of `builtins.open` -- a peer's patch of `builtins.open` missed 3,234 of
3,259 reads. But it only sees what the PYTHON runtime does. A C library that opens a file itself --
htslib, and so pysam, which this project streams BAMs and VCFs through -- is INVISIBLE to it, and so is
any other extension module, any subprocess the guard was not put into, and `mmap` writes to a mapping
taken before the hook armed. So this is a guard against Python-level writes and must never be described
as absolute.

THE TWO MECHANISMS ARE COMPLEMENTARY, NOT REDUNDANT, and this is the part a reader should take away.
The hook NAMES the write -- it says which event, which path and which store, before a byte moves, which
is what makes a refusal diagnosable instead of a mysterious EACCES. The MODE PREVENTS the write, and
only the mode reaches C: `cp -c` plus 0444 is enforced by the kernel against every writer in any
language, htslib included, and it covers exactly the case the hook cannot see. So the stores a verdict
worktree reads are linked at 0444 (`scripts/link_stores_read_only.py`), and the hook is what tells the
reader which code tried to write and where.

THE CHILD PROCESS. `manifest_rebuild` runs the result's own command in a SUBPROCESS, and an audit hook
installed in the parent does not reach it -- which is where the risk actually lives, since the parent only
reads and hashes. So the parent also writes a `sitecustomize.py` beside the worktree that imports THIS
file and calls `install_from_env`, and puts that directory on the child's PYTHONPATH; CPython's `site`
imports `sitecustomize` at startup, so the hook is in before the first line of the result's writer runs.
There is one implementation and both processes use it. Whether the child really armed is not assumed: the
child writes `MARKER` into the guard directory and the parent reports its absence as a reason no verdict
may be given.
"""

from __future__ import annotations

import contextlib
import os
import sys
from pathlib import Path
from typing import Any

#: The protected store roots, `os.pathsep`-separated real paths, as the parent passes them to the child.
ENV_ROOTS = "GENOMEOS_REBUILD_PROTECTED_STORES"

#: Where this file is, so a generated `sitecustomize.py` can import it without `genomeos` on the path:
#: the worktree is at an older revision and its `genomeos` package may predate this guard entirely.
ENV_GUARD = "GENOMEOS_REBUILD_WRITE_GUARD"

#: The file the child writes once the hook is armed. Evidence, so the parent never assumes it.
MARKER = "write-guard-armed"

#: What a per-process read record is called. The pid is in the name, which is what keeps two processes
#: from writing one file.
READS_PREFIX = "reads-"

#: Where that evidence goes. A separate variable from ENV_GUARD because ENV_GUARD names the guard FILE,
#: and deriving the directory from it put the marker in `scripts/` -- a stray file in the repository and
#: no evidence where the parent looked for it. Caught by the test that asks for the evidence by path.
ENV_MARKER_DIR = "GENOMEOS_REBUILD_WRITE_GUARD_DIR"

#: Open flags that mean the file is being written. Kept identical to `genomeos.manifest._WRITE_FLAGS`;
#: `tests/test_rebuild_write_guard.py` pins the two against each other over a table of modes and flags,
#: because this file cannot import `genomeos` (see ENV_GUARD) and two copies that nothing compares are
#: how ten copies of the cleanliness block drifted apart.
_WRITE_FLAGS = os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_APPEND | os.O_TRUNC

#: Audit events that change a path, and which of each event's arguments name the path being changed.
#: `open` is not here: a mode and a flag decide whether it writes, which `is_write` does.
#: `os.rename` and `shutil.move` name both, because the source is removed as well as the target written.
WATCHED: dict[str, tuple[int, ...]] = {
    "os.remove": (0,),
    "os.rename": (0, 1),
    "os.mkdir": (0,),
    "os.rmdir": (0,),
    "os.truncate": (0,),
    "os.chmod": (0,),
    "os.utime": (0,),
    "os.link": (1,),
    "os.symlink": (1,),
    "shutil.copyfile": (1,),
    "shutil.move": (0, 1),
    "shutil.rmtree": (0,),
}

# --- read recording (2026-10-02) ------------------------------------------------------------------
#
# The guard above refuses writes and records nothing, so reads go in HERE rather than into a second
# instrument: the project would otherwise have two tracers to keep in step, and the hook is already
# installed and already sees every open.
#
# Why it is needed, measured end to end on a live result rather than reasoned: `manifest_rebuild` is
# BLIND BY CONSTRUCTION to any file a manifest does not name. `constrained_unknown_targets` rebuilt at
# 193 of 193 inputs declared and checked, 193 opened and hashed, 3,858 of 3,896 leaves compared, no
# must_hold failures, "0 differences AFTER SETTING ASIDE 1 leaves" -- while a tracer on the same run saw
# a 194TH read, `data/results/unknown_chr21.json`, pinned by no declared sha256. Two instruments
# disagreed about one run and the flattering one was the clean bill of health.
#
# THE TWO WAYS A RECORDER GIVES A PASS IT HAS NOT EARNED, both of which happened, and what is done here
# instead:
#
# (a) PATCHING `builtins.open`. It leaves `io.open` -- and therefore all of `pathlib`, and pandas --
#     invisible: 25 reads recorded where there were 3,259, with one declared input reported "never read"
#     when it had been read 3,258 times. An audit hook on the `open` EVENT sees `io.open`, `os.open` and
#     `pathlib` alike, which is the whole reason the instrument is a hook and not a patch.
#
# (b) ONE RECORD DUMPED AT `atexit`. Under `--workers 2` a pool child that read nothing reached exit
#     last and OVERWROTE the record, and the result came back with 0 paths opened -- satisfying "every
#     read is declared" VACUOUSLY. So the record is APPEND-ONLY, written with `os.write` at the moment
#     of the read, to a file named for the process that made it. A silent child cannot overwrite a
#     reader's record because it never opens the reader's file, and a record that was never written is
#     an ABSENT file, which is a refusal rather than an empty pass.
#
# THE LIMIT, in the code and not only in a report: an audit hook sees what the PYTHON runtime does. A C
# library that opens a file itself -- htslib, and so pysam, which this project streams BAMs and VCFs
# through -- is INVISIBLE to it. So this records Python-level reads, and a claim that every read was
# declared is a claim about Python-level reads only. It is not a claim about every byte the process read.
#
# WHAT IS RECORDED: one line per DISTINCT path, the first time it is read, as `pid\tread\tpath`. Not a
# count: how many times a path was read is NOT recorded, because a run that opens the same table a
# million times would otherwise write a million lines. The set is what the comparison needs -- declared
# against read -- and the count's absence is said rather than left to be assumed.

#: Where the per-process read records go, as the parent passes it to the child.
ENV_RECORD_DIR = "GENOMEOS_REBUILD_READ_RECORD_DIR"

#: Distinct paths remembered before deduping stops. A bound, so a process opening millions of distinct
#: paths cannot grow the set without limit; past it every read is written, so nothing is lost but volume.
_DEDUPE_CAP = 1 << 16

#: The append-only file descriptor for this process's reads, and the paths already written.
_RECORD_FD: int | None = None
_SEEN: set[str] = set()

#: The armed roots. Empty means the hook is installed and inert, which is the state a test session and
#: every other process is left in.
_ROOTS: tuple[str, ...] = ()
_INSTALLED = False


def is_write(mode: Any, flags: Any) -> bool:
    """Whether an `open` audit event opened its file to write it. `mode` is the string `open` was given,
    or None for `os.open`, where only the flags say. The same test as `genomeos.manifest._is_write`."""
    if isinstance(mode, str):
        return any(c in mode for c in "wax+")
    return bool(isinstance(flags, int) and flags & _WRITE_FLAGS)


def _as_path(value: Any) -> str | None:
    """The string a path argument names, or None when it names no path (an integer file descriptor)."""
    if isinstance(value, int) or value is None:
        return None
    if isinstance(value, bytes):
        try:
            return os.fsdecode(value)
        except (UnicodeDecodeError, ValueError):
            return None
    try:
        return os.fspath(value) if not isinstance(value, str) else value
    except TypeError:
        return None


def inside(path: str, roots: tuple[str, ...]) -> str | None:
    """The protected root `path` resolves into, or None. Real paths on both sides, so a symlink into a
    store is refused exactly as an absolute path to it is, and a path that does not exist yet still
    resolves through the directories above it."""
    try:
        real = os.path.realpath(path)
    except (OSError, ValueError):
        return None
    for root in roots:
        if real == root or real.startswith(root + os.sep):
            return root
    return None


def offence(event: str, args: tuple[Any, ...], roots: tuple[str, ...]) -> tuple[str, str] | None:
    """(the path, the store it is inside) when this event would change a protected store, else None.

    Pure, and taking its roots as an argument, so both cases of the counterfactual can be put to it
    without arming anything.
    """
    if not roots:
        return None
    if event == "open":
        path, mode, flags = (list(args) + [None, None, None])[:3]
        if not is_write(mode, flags):
            return None
        indices: tuple[int, ...] = (0,)
        candidates = [path]
    elif event in WATCHED:
        indices = WATCHED[event]
        candidates = [args[i] if i < len(args) else None for i in indices]
    else:
        return None
    for value in candidates:
        p = _as_path(value)
        if p is None:
            continue
        root = inside(p, roots)
        if root is not None:
            return p, root
    return None


def record_read(path: str) -> None:
    """Append one line for `path`, now, with `os.write`. No buffer, no atexit, no shared file.

    Failures are swallowed: a recorder that raised inside an audit hook would break every open in the
    process. An unwritten record shows up as a record with nothing in it, and `manifest_rebuild`
    withholds the verdict by name on that rather than reading it as "no undeclared reads".
    """
    if _RECORD_FD is None:
        return
    if len(_SEEN) < _DEDUPE_CAP:
        if path in _SEEN:
            return
        _SEEN.add(path)
    with contextlib.suppress(OSError):
        os.write(_RECORD_FD, f"{os.getpid()}\tread\t{path}\n".encode())


def _hook(event: str, args: tuple[Any, ...]) -> None:
    if event == "open" and _RECORD_FD is not None:
        path, mode, flags = (list(args) + [None, None, None])[:3]
        if not is_write(mode, flags):
            p = _as_path(path)
            if p is not None:
                record_read(p)
    found = offence(event, args, _ROOTS)
    if found is None:
        return
    path, root = found
    raise PermissionError(
        f"a rebuild may not write to the store it verifies a result against: {event} on {path} resolves "
        f"into {root}, which is linked into the worktree read-only. THIS IS NOT A VERIFICATION FAILURE "
        f"and it says nothing about whether the result reproduces: the rebuilt command tried to WRITE to "
        f"an input store, which is most often a writer caching into data/cache beside its real work. "
        f"Before 2026-10-02 that write went through to the only copy of the store on this machine, "
        f"silently, which is why it is refused rather than allowed. Rerun the writer in the checkout if "
        f"the cache needs filling; the rebuild must not be the thing that fills it"
    )


def record_to(directory: str | os.PathLike[str] | None) -> str | None:
    """Open this process's append-only read record in `directory`, or stop recording when None.

    `O_APPEND` and a name carrying the pid are the whole mechanism: two processes never write the same
    file, and every write lands at the end, so nothing a process recorded can be overwritten by a
    process that recorded nothing. There is no `atexit`, no buffer and no single shared file, because a
    record assembled at exit is what a silent pool child overwrote.
    """
    global _RECORD_FD, _SEEN
    if _RECORD_FD is not None:
        with contextlib.suppress(OSError):
            os.close(_RECORD_FD)
        _RECORD_FD = None
    _SEEN = set()
    if directory is None:
        return None
    d = os.fspath(directory)
    os.makedirs(d, exist_ok=True)
    path = os.path.join(d, f"{READS_PREFIX}{os.getpid()}.tsv")
    _RECORD_FD = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    return path


def reads_recorded(directory: str | os.PathLike[str]) -> dict[str, Any]:
    """Every distinct path recorded in `directory`, with the processes that recorded it.

    An ABSENT or empty record is reported as `recorded: False`, never as "no reads": the difference is
    the whole of case (b), where a record nobody wrote satisfied "every read is declared" vacuously.
    """
    d = Path(os.fspath(directory))
    files = sorted(d.glob(f"{READS_PREFIX}*.tsv")) if d.is_dir() else []
    paths: set[str] = set()
    pids: set[int] = set()
    for f in files:
        with contextlib.suppress(OSError):
            for line in f.read_text(errors="replace").splitlines():
                parts = line.split("\t", 2)
                if len(parts) == 3 and parts[1] == "read":
                    if parts[0].isdigit():
                        pids.add(int(parts[0]))
                    paths.add(parts[2])
    return {
        "recorded": bool(paths),
        "record_files": [f.name for f in files],
        "processes": sorted(pids),
        "paths": sorted(paths),
        "distinct_paths": len(paths),
        "counts_not_recorded_because": (
            "one line per DISTINCT path, so how many times a path was read is not recorded: a run that "
            "opens one table a million times would otherwise write a million lines"
        ),
        "limit": (
            "Python-level reads only. An audit hook sees io.open, os.open and pathlib alike, which is "
            "why it is a hook and not a patch of builtins.open -- that patch recorded 25 reads where "
            "there were 3,259. It does NOT see a C library opening a file itself, so htslib and pysam "
            "reads are invisible and 'every read was declared' is a claim about Python-level reads"
        ),
    }


def arm(roots: object) -> tuple[str, ...]:
    """Install the hook if it is not in yet and arm it over `roots`. Returns the roots now armed.

    Installed once and never removed, because `sys.addaudithook` cannot be undone; `disarm` empties the
    roots instead, which leaves the hook in place and inert. The pattern is `genomeos.manifest`'s
    `trace_begin`, deliberately, so there is one way audit hooks are handled in this repository.
    """
    global _ROOTS, _INSTALLED
    if isinstance(roots, str):
        roots = roots.split(os.pathsep)
    real = []
    for r in roots or ():
        p = _as_path(r)
        if p and p.strip():
            real.append(os.path.realpath(p))
    _ROOTS = tuple(sorted(set(real), key=lambda s: -len(s)))
    if not _INSTALLED:
        sys.addaudithook(_hook)
        _INSTALLED = True
    return _ROOTS


def disarm() -> None:
    """Empty the roots and stop recording. The hook stays installed and passes everything, which is its
    resting state: `sys.addaudithook` cannot be undone, so this is the only state that can be restored."""
    global _ROOTS
    _ROOTS = ()
    record_to(None)


def armed() -> tuple[str, ...]:
    return _ROOTS


def install_from_env() -> tuple[str, ...]:
    """Arm from ENV_ROOTS and leave MARKER beside the guard, so the parent has evidence rather than a
    hope. Called by the generated `sitecustomize.py` in the rebuild's child process."""
    roots = arm(os.environ.get(ENV_ROOTS, ""))
    record_to(os.environ.get(ENV_RECORD_DIR) or None)
    where = os.environ.get(ENV_MARKER_DIR)
    if roots and where:
        marker = os.path.join(where, MARKER)
        try:
            with open(marker, "a") as f:
                f.write(f"{os.getpid()}\t{os.pathsep.join(roots)}\n")
        except OSError:
            pass
    return roots


#: What the parent writes beside the worktree for the child to import at startup. It imports THIS file by
#: its absolute path, so the worktree's own (older) `genomeos` package is never involved.
SITECUSTOMIZE = """# SPDX-License-Identifier: AGPL-3.0-or-later
# Written by scripts/manifest_rebuild.py for the duration of one rebuild. It arms the write guard in the
# rebuilt command's own process, which an audit hook in the parent cannot reach.
import importlib.util
import os

_guard = os.environ.get("GENOMEOS_REBUILD_WRITE_GUARD")
if _guard:
    _spec = importlib.util.spec_from_file_location("genomeos_rebuild_write_guard", _guard)
    _mod = importlib.util.module_from_spec(_spec)
    _spec.loader.exec_module(_mod)
    _mod.install_from_env()
"""


def write_sitecustomize(directory: str | os.PathLike[str], guard: str | os.PathLike[str]) -> dict[str, Any]:
    """Put `sitecustomize.py` in `directory` and return the environment that arms the child.

    `directory` is OUTSIDE the worktree on purpose: a file inside it would make the rebuilt result's own
    cleanliness block record an untracked path, which is a difference this tool would then have to
    explain away.
    """
    directory = os.fspath(directory)
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "sitecustomize.py"), "w") as f:
        f.write(SITECUSTOMIZE)
    return {
        ENV_GUARD: os.path.abspath(os.fspath(guard)),
        ENV_MARKER_DIR: os.path.abspath(directory),
        ENV_RECORD_DIR: os.path.abspath(directory),
    }


def marker_pids(directory: str | os.PathLike[str]) -> list[int]:
    """The process ids that recorded an armed guard in `directory`; empty when none did."""
    try:
        with open(os.path.join(os.fspath(directory), MARKER)) as f:
            text = f.read()
    except OSError:
        return []
    out = []
    for line in text.splitlines():
        head = line.split("\t", 1)[0].strip()
        if head.isdigit():
            out.append(int(head))
    return out
