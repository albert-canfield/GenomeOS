# SPDX-License-Identifier: AGPL-3.0-or-later
"""Rebuild a result from its manifest in a clean second checkout, and say what differs or is missing.

    uv run python scripts/manifest_rebuild.py RESULT.json --where DIR [--venv fresh|shared] [--keep]

Review item R9's acceptance: a second environment reconstructs a representative result from its
manifest, or names exactly which dependency is unavailable. The steps, each of which can stop the
rebuild with the reason:

1. read `result_manifest` from RESULT.json (a result without one cannot be rebuilt: said so);
2. `git worktree add --detach` the recorded `code.git_sha` under DIR; uncommitted code at write time
   (`code.dirty_code_paths`) is named, since the commit alone did not run;
3. link the local data stores (data/reference, data/knowledge, data/cache) read-only, as the
   pre-push hook does, link the machine-local inputs under data/results read-only as well (see
   `link_machine_local_inputs`), and check every input's sha256 against the manifest, naming any input
   that is absent or has different bytes. Every declared input is accounted for: a group of files read
   together is opened and hashed member by member, and an input the tool cannot resolve stops the
   rebuild by name. `inputs_declared`, `inputs_checked` and `inputs_unchecked` carry the denominator,
   so an input that was not opened cannot leave the list and read as one that matched;
   `inputs_satisfied_from_machine_local_paths` says how much of the rebuild rested on this machine;
4. run `code.argv` in the worktree, in a fresh environment from the committed uv.lock (`--venv
   fresh`, offline) or the checkout's own (`--venv shared`);
5. compare the rebuilt result with RESULT.json field by field, ignoring only `date` and the
   `code` block of the manifest (which records the run, not the result).
   Since 2026-09-28 wall-clock timing keys are ignored too, at any depth (`seconds`, `*_seconds`,
   `seconds_*`, `*_per_second`); the report lists their paths under `timing_fields_ignored`.
   Since 2026-10-02 the run-resource readings named in `RESOURCE_KEYS` are ignored too, by exact key,
   and listed under `resource_fields_ignored`.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

from genomeos import manifest as mf

STORES = ("reference", "knowledge", "cache")
IGNORED = ("date",)

#: Where a result's machine-local inputs live. data/results is git-ignored as a rule (docs/DATA.md), so a
#: result computed from other results declares inputs that no clean worktree can hold: reader v1's cached
#: DNase peak sets, 13 biosamples by 23 chromosomes, are on this machine and in no commit. Until 2026-10-02
#: the rebuild linked only the three data stores, so those inputs were reported absent one by one and the
#: rebuild stopped before any comparison: response_map_increment2 reached 111 of 410 inputs opened, and
#: context_evidence hit the same wall. They are linked here instead, read-only, and hashed like any other.
MACHINE_LOCAL_DIR = "data/results"

#: What an input's `sha256` says when nothing was hashed to produce it. Such an input cannot be
#: checked, which is a reason to stop and say so, not a reason to pass.
NOT_HASHED = "n/a"

#: The fields a rebuild in a clean worktree is allowed to differ on, by **exact path**, because they
#: describe the tree the run happened in and not anything the run computed. A result that honestly records
#: the other lanes' outstanding files in a shared checkout cannot reproduce those lists in a clean worktree,
#: where there are none; under the unamended rule such a result could only ever fail, while one that omitted
#: the record would pass. Adopted 2026-10-02 after the coordinator raised it and the supervisor ruled.
ENVIRONMENT_FIELDS = tuple(
    f"/{container}{field}"
    for container in ("code_cleanliness/", "result_manifest/code_cleanliness/", "result_manifest/code/")
    for field in (
        "dirty",
        "foreign_uncommitted_code",
        "dirty_code_paths",
        "untracked_code_paths",
        "dirty_result_paths",
    )
)

#: These must hold on **both** sides, so the exemption above can never excuse a result whose own code was
#: uncommitted, or one with a foreign uncommitted file on the path that produced its numbers.
MUST_HOLD: tuple[tuple[str, Any], ...] = (
    ("own_code_is_committed", True),
    ("foreign_uncommitted_code_on_the_counting_path", []),
)


def no_verdict_reason(declared: int, checked: int) -> str | None:
    """Why no verdict may be reported over these inputs, or None when one may.

    The one rule the tool must never be able to break: a clean verdict over inputs it never opened. The
    project's lesson of 2026-10-01 is that a verification tool can fail in the direction that flatters, and
    this is where that would happen -- a report that compared the fields of two files and said nothing about
    the bytes behind them reads as a reproduction. So the decision lives in one function, is asked for both
    before the run and again before the verdict is written, and says yes only when at least one input was
    declared and every declared input was opened and hashed.
    """
    if declared <= 0:
        return (
            "the manifest declares no inputs, so no bytes were opened: a comparison of fields would say "
            "nothing about what the run read, and no verdict is reported"
        )
    if checked != declared:
        return f"{checked} of {declared} declared inputs were opened and hashed: no verdict is reported"
    return None


def leaves(x: Any, where: str = "") -> int:
    """How many leaf values a comparison covers. A count of top-level keys can hide a nested difference:
    `cell2_eligibility.json` has 14 top-level keys and several hundred leaves."""
    if isinstance(x, dict):
        return sum(leaves(v, f"{where}/{k}") for k, v in x.items())
    if isinstance(x, list):
        return sum(leaves(v, f"{where}[{i}]") for i, v in enumerate(x))
    return 1


def leaf_reconciliation(original: dict[str, Any]) -> dict[str, Any]:
    """Where every leaf of the original went, so a count carries its denominator.

    `compared` plus `not_compared` equals `total`: a comparison reported without the leaves it dropped
    invites the reader to assume it covered the file. The dropped leaves are the ones `comparable()`
    removes before any diff -- the date, the manifest's own `code` block, `model_dependencies` and every
    timing key -- and they are named here by cause rather than merely counted.
    """
    total = leaves(original)
    compared = leaves(comparable(original))
    return {
        "total": total,
        "compared": compared,
        "not_compared": total - compared,
        "not_compared_because": {
            "ignored_keys": list(IGNORED),
            "manifest_code_block": "result_manifest.code: the stamp of the run, not a value it computed",
            "model_dependencies": "recorded after the fact for results made before 2026-09-28",
            "timing_keys": "every key is_timing() matches, listed under timing_fields_ignored",
            "resource_keys": (
                "the exact keys in RESOURCE_KEYS, a run-resource reading rather than a value the run "
                "computed, listed under resource_fields_ignored"
            ),
        },
        "reconciles": total == compared + (total - compared),
    }


def environment_differences(found: list[str]) -> tuple[list[str], list[str]]:
    """(real differences, environment differences). A path counts as environmental only if it is exactly
    one of ENVIRONMENT_FIELDS, or an element of one of those lists; no pattern matching."""
    real, env = [], []
    for d in found:
        path = d.split(":")[0]
        base = path.split("[")[0]
        (env if base in ENVIRONMENT_FIELDS else real).append(d)
    return real, env


def must_hold_failures(original: dict[str, Any], rebuilt: dict[str, Any]) -> list[str]:
    """Every place either side breaks MUST_HOLD. Checked at any depth, on both sides."""
    out = []

    def walk(x: Any, side: str, where: str = "") -> None:
        if isinstance(x, dict):
            for k, v in x.items():
                for field, wanted in MUST_HOLD:
                    if k == field and v != wanted:
                        out.append(f"{side}{where}/{k}: {v!r}, must be {wanted!r}")
                walk(v, side, f"{where}/{k}")
        elif isinstance(x, list):
            for i, v in enumerate(x):
                walk(v, side, f"{where}[{i}]")

    walk(original, "original")
    walk(rebuilt, "rebuilt")
    return out


def is_timing(key: Any) -> bool:
    """A wall-clock key: `seconds`, or a snake_case name with a `seconds` token or ending `per_second`.
    `second` alone (an ordinal, as in second_endpoint) and `duration` (often biological) are compared."""
    if not isinstance(key, str):
        return False
    tokens = key.lower().split("_")
    return "seconds" in tokens or tokens[-2:] == ["per", "second"]


#: Run-resource readings, exempt from comparison by **exact key**. Added 2026-10-02 on the supervisor's
#: ruling, and deliberately not through ENVIRONMENT_FIELDS, which is the code-cleanliness list and stays
#: that way. A peak memory figure is of the same class as a wall-clock second: what the machine spent to
#: produce the result, not anything the result asserts. Two otherwise identical rebuilds differed in this
#: one leaf and nothing else -- `/cost/peak_memory_mb` 1142.5 against 1148.8 in placement_audit, and
#: `/compute/peak_rss_mb` 1300.5 against 1298.2 in context_contrast_feasibility. Exact keys only: a pattern
#: such as `*_mb` would one day swallow a leaf a run did compute, and the point of a short list is that
#: adding to it is a decision somebody makes on the record.
RESOURCE_KEYS = ("peak_rss_mb", "peak_memory_mb")


def is_resource(key: Any) -> bool:
    """A run-resource reading, matched by exact key and never by pattern, so a name that merely looks like
    one -- `peak_signal_mb`, `peak_rss_mb_per_cell`, `memory_mb_budget` -- is compared like any other leaf."""
    return isinstance(key, str) and key in RESOURCE_KEYS


def _ignored_key(key: Any) -> bool:
    """A key comparable() drops: a wall-clock reading or a run-resource one. Both are reported by path,
    under timing_fields_ignored and resource_fields_ignored, so a reader sees what was set aside."""
    return is_timing(key) or is_resource(key)


def _strip_timing(x: Any) -> Any:
    if isinstance(x, dict):
        return {k: _strip_timing(v) for k, v in x.items() if not _ignored_key(k)}
    if isinstance(x, list):
        return [_strip_timing(v) for v in x]
    return x


def _paths_where(x: Any, pred: Any, where: str = "") -> list[str]:
    if isinstance(x, dict):
        out = []
        for k, v in x.items():
            out.extend([f"{where}/{k}"] if pred(k) else _paths_where(v, pred, f"{where}/{k}"))
        return out
    if isinstance(x, list):
        return [p for i, v in enumerate(x) for p in _paths_where(v, pred, f"{where}[{i}]")]
    return []


def timing_paths(x: Any, where: str = "") -> list[str]:
    """Every path at which comparable() drops a timing key."""
    return _paths_where(x, is_timing, where)


def resource_paths(x: Any, where: str = "") -> list[str]:
    """Every path at which comparable() drops a run-resource key."""
    return _paths_where(x, is_resource, where)


def diff(a: Any, b: Any, where: str = "") -> list[str]:
    """Every path at which two JSON values differ."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                out.append(f"{where}/{k}: only in {'rebuilt' if k in b else 'original'}")
            else:
                out.extend(diff(a[k], b[k], f"{where}/{k}"))
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{where}: {len(a)} items vs {len(b)}"]
        return [d for i, (x, y) in enumerate(zip(a, b, strict=True)) for d in diff(x, y, f"{where}[{i}]")]
    return [] if a == b else [f"{where}: {a!r} vs {b!r}"[:300]]


def comparable(payload: dict[str, Any]) -> dict[str, Any]:
    out = {k: v for k, v in payload.items() if k not in IGNORED}
    m = dict(out.get(mf.KEY) or {})
    m.pop("code", None)
    # recorded after the fact for results made before 2026-09-28 (R9 follow-up), so a rerun at their
    # revision cannot write it: it states the inputs' provenance, not a value the run computes
    m.pop("model_dependencies", None)
    if m:
        out[mf.KEY] = m
    out = _strip_timing(out)
    return out


def is_group(entry: dict[str, Any]) -> bool:
    """Whether a declared input stands for several files read together rather than for one path.

    `files_entry` marks a group since 2026-10-02. Before that a group carried `files` and so was
    indistinguishable from `input_entry` on a directory, which carries `files` too; such an entry is
    not treated as a group here, and falls to the unresolvable branch, which says so by name.
    """
    return entry.get("group") is True or isinstance(entry.get("members"), list)


def declared_paths(inputs: Any) -> list[str]:
    """Every repository-relative path a manifest's inputs name, a group's members included.

    A label is not a path, so a group contributes its members and not its own `path`; a group recorded
    before its members were named contributes nothing, which is why such a group remains a named failure
    in `_check_group` rather than something this could quietly satisfy.
    """
    out: list[str] = []
    for i in inputs if isinstance(inputs, list) else []:
        if not isinstance(i, dict):
            continue
        if is_group(i):
            out += [
                m["path"]
                for m in i.get("members") or []
                if isinstance(m, dict) and isinstance(m.get("path"), str)
            ]
        elif isinstance(i.get("path"), str):
            out.append(i["path"])
    return out


def git_ignored(root: Path, paths: list[str]) -> set[str]:
    """Which of `paths` git ignores in `root`: the machine-local ones, answered for all of them at once.

    A path git does not ignore is never linked from this machine, whatever is sitting at it: for a tracked
    path the committed tree at the rebuild's sha is what the rebuild must read, and a path that is neither
    tracked nor ignored is not a machine-local store but a stray file. Both are left to be named absent.
    """
    if not paths:
        return set()
    r = subprocess.run(
        ["git", "-C", str(root), "check-ignore", "--stdin"],
        input="\n".join(paths),
        capture_output=True,
        text=True,
    )
    return {line for line in r.stdout.splitlines() if line}


def link_read_only(src: Path, dst: Path) -> None:
    """Link one machine-local file into the worktree so the rebuild can read it and cannot write it.

    `cp -c` asks APFS for a clone: no bytes are copied, both sides share them until one is written, and the
    clone is a separate inode -- so the 0o444 that follows takes the write bits off the worktree's side only
    and never off the machine's one and only copy. That is the reason for a clone rather than a symlink or a
    hard link: both of those carry the machine's own mode, so `open(path, "w")` through either would truncate
    the machine's file, and the only way to refuse it would be to chmod the original, which other sessions
    share. Through this link such an open fails with PermissionError before a byte is written, and the
    rebuild cannot write a result into the directory it is reading its inputs from. On a filesystem with no
    clone support the fallback is a plain copy, read-only in the same way; a worktree copy that is still
    writable after the chmod is an error, not a warning.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    if subprocess.run(["cp", "-c", str(src), str(dst)], capture_output=True).returncode:
        shutil.copyfile(src, dst)
    dst.chmod(0o444)
    if dst.stat().st_mode & 0o222:
        raise OSError(f"{dst} is writable after chmod 0444: the rebuild could overwrite the input it reads")


def link_machine_local_inputs(
    worktree: Path, root: Path, inputs: Any, output: str | None = None
) -> dict[str, Any]:
    """Link the git-ignored data/results inputs this manifest declares into the worktree, read-only.

    Only a path the manifest declares is linked, only under MACHINE_LOCAL_DIR, only when git ignores it in
    this checkout, and only when the worktree does not already hold that path -- a file the commit carries
    is read from the commit. Nothing here decides whether an input matches: the bytes linked are hashed
    afterwards against the sha256 the manifest declares, so a machine-local input whose bytes have since
    changed is a reported difference and not a silent pass. `output` is the result the rebuild will write,
    which is this run's product rather than its input and is never linked over.
    """
    out: dict[str, Any] = {
        "linked": [],
        "absent_on_this_machine": [],
        "not_linked_because_git_does_not_ignore_it": [],
        "not_linked_because_it_is_the_output": [],
        "failures": [],
    }
    want = sorted(
        {
            p
            for p in declared_paths(inputs)
            if not os.path.isabs(p)
            and ".." not in Path(p).parts
            and (p == MACHINE_LOCAL_DIR or p.startswith(MACHINE_LOCAL_DIR + "/"))
            and not (worktree / p).exists()
        }
    )
    ignored = git_ignored(root, want)
    for p in want:
        if output is not None and p == output:
            out["not_linked_because_it_is_the_output"].append(p)
            continue
        if p not in ignored:
            out["not_linked_because_git_does_not_ignore_it"].append(p)
            continue
        src = root / p
        if not src.exists():
            out["absent_on_this_machine"].append(p)
            continue
        try:
            if src.is_dir():
                for f in sorted(x for x in src.rglob("*") if x.is_file()):
                    link_read_only(f, worktree / p / f.relative_to(src))
            else:
                link_read_only(src, worktree / p)
        except OSError as e:
            out["failures"].append(f"{p}: {e}")
            continue
        out["linked"].append(p)
    return out


def unlock_links(worktree: Path, paths: list[str]) -> None:
    """Put the write bits back on the linked copies, so removing the worktree cannot be refused by them.
    Each is a clone of its own, so this touches nothing the machine keeps."""
    for p in paths:
        target = worktree / p
        for f in sorted(x for x in target.rglob("*") if x.is_file()) if target.is_dir() else [target]:
            if f.is_file():
                os.chmod(f, 0o644)


def check_inputs(
    worktree: Path, inputs: Any, root: Path, machine_local: set[str] | None = None
) -> dict[str, Any]:
    """Open and hash every declared input, a group member by member, and account for all of them.

    The contract the caller rests on: `unavailable` is empty only when at least one input was declared
    and every declared input was opened and hashed, and `entries` holds one record per declared input,
    so nothing can leave the list silently. An input this cannot resolve is a named failure.

    Written 2026-10-02 to replace a loop that resolved every input by `entry["path"]`. Because
    `files_entry` records a label there and not a path, each group was reported "absent" with no file
    opened, and its entry never reached `inputs`: a list with no denominator, which read as though
    every input had matched. The direction of that error flattered the result -- a result whose inputs
    were grouped was excused as resting on an unavailable dependency, when its files were on disk and
    could have been compared.
    """
    out: dict[str, Any] = {
        "declared": 0,
        "checked": 0,
        "files_opened": 0,
        # how many declared inputs were opened and hashed from a git-ignored, machine-local path rather
        # than from the committed tree: the share of the rebuild that rests on this machine
        "machine_local_opened": 0,
        "entries": [],
        "unavailable": [],
        "checked_outside_the_worktree": [],
    }
    local = machine_local or set()
    if not isinstance(inputs, list) or not inputs:
        out["unavailable"].append(
            "the manifest declares no inputs, so nothing was opened or hashed: a comparison of this "
            "result's fields would say nothing about the bytes the run read"
        )
        return out
    out["declared"] = len(inputs)
    linked = [str((root / "data" / d).resolve()) for d in STORES]
    for i in inputs:
        if not isinstance(i, dict) or not isinstance(i.get("path"), str) or not i["path"]:
            out["entries"].append({"path": None, "sha256_matches": False, "problem": "no path"})
            out["unavailable"].append(f"a declared input carries no path: {i!r}"[:200])
            continue
        path = i["path"]
        if is_group(i):
            _check_group(worktree, i, out, local)
        elif i.get("sha256") == NOT_HASHED:
            out["entries"].append(
                {"path": path, "sha256_matches": False, "problem": f"sha256 is {NOT_HASHED!r}"}
            )
            out["unavailable"].append(
                f"input {path} records sha256 {NOT_HASHED!r}: nothing was hashed when the result was "
                "written, so there is nothing for a rebuild to check these bytes against"
            )
        else:
            _check_one_path(worktree, i, linked, out, local)
    return out


def _check_group(
    worktree: Path, i: dict[str, Any], out: dict[str, Any], local: set[str] | None = None
) -> None:
    """One grouped input: every member opened and hashed on its own, then the group's own digest."""
    path, members = i["path"], i.get("members")
    named = [m for m in members or [] if isinstance(m, dict) and isinstance(m.get("path"), str)]
    if not members or len(named) != len(members):
        out["entries"].append(
            {
                "path": path,
                "group": True,
                "members_declared": i.get("files"),
                "members_opened_and_hashed": 0,
                "sha256_matches": False,
                "problem": "the group names no members",
            }
        )
        out["unavailable"].append(
            f"input {path} is a group of {i.get('files', 'an unstated number of')} files recorded under "
            "a label whose members are not named, so not one of them could be opened: it was written "
            "before files_entry named its members (2026-10-02) and cannot be checked from this manifest"
        )
        return
    names = sorted(m["path"] for m in named)
    absent = [n for n in names if not (worktree / n).exists()]
    if absent:
        out["entries"].append(
            {
                "path": path,
                "group": True,
                "members_declared": len(names),
                "members_opened_and_hashed": 0,
                "absent_members": absent[:20],
                "sha256_matches": False,
            }
        )
        out["unavailable"].append(
            f"input {path}: {len(absent)} of {len(names)} files in the group are absent, first {absent[:3]}"
        )
        return
    try:
        digest, size, seen = mf.group_digest(names, root=worktree)
    except OSError as e:
        out["entries"].append(
            {
                "path": path,
                "group": True,
                "members_declared": len(names),
                "members_opened_and_hashed": 0,
                "sha256_matches": False,
                "problem": f"a member cannot be read: {e}",
            }
        )
        out["unavailable"].append(f"input {path}: a file in the group cannot be read: {e}")
        return
    out["files_opened"] += len(seen)
    recorded = {m["path"]: m.get("sha256") for m in named}
    differing = [s["path"] for s in seen if recorded.get(s["path"]) not in (None, s["sha256"])]
    ok = digest == i.get("sha256")
    from_machine = sorted(set(names) & (local or set()))
    out["entries"].append(
        {
            "path": path,
            "group": True,
            "members_declared": len(names),
            "members_opened_and_hashed": len(seen),
            "members_with_different_bytes": differing,
            "members_from_machine_local_paths": len(from_machine),
            "sha256_matches": ok and not differing,
            "bytes": size,
        }
    )
    if from_machine:
        # a group counts once, as one declared input, whichever of its members came from this machine
        out["machine_local_opened"] += 1
    if differing:
        out["unavailable"].append(
            f"input {path}: {len(differing)} of {len(names)} files in the group have different bytes, "
            f"first {differing[:3]}"
        )
        return
    if not ok:
        out["unavailable"].append(
            f"input {path}: the group's own digest differs from the manifest's (sha256 {digest[:12]})"
        )
        return
    out["checked"] += 1


def _check_one_path(
    worktree: Path, i: dict[str, Any], linked: list[str], out: dict[str, Any], local: set[str] | None = None
) -> None:
    """One input recorded as a single path, which may be a file or a directory."""
    path = i["path"]
    outside: dict[str, Any] | None = None
    if os.path.isabs(path):
        # `worktree / path` leaves the worktree when path is absolute, so the bytes hashed would not be
        # the bytes the rebuild reads -- except under data/reference, data/knowledge and data/cache,
        # which are linked into the worktree and so are the same file.
        real = str(Path(path).resolve())
        store = next((s for s in linked if real == s or real.startswith(s + os.sep)), None)
        if store is None:
            out["entries"].append(
                {"path": path, "sha256_matches": False, "problem": "absolute, outside the linked stores"}
            )
            out["unavailable"].append(
                f"input {path} is an absolute path outside the linked data stores: the bytes there are "
                "not the bytes this worktree reads, so it cannot be checked against this manifest"
            )
            return
        outside = {"path": path, "linked_store": store}
    p = worktree / path
    if not p.exists():
        # `files` on an entry with no member list reads as a group recorded under a label before
        # files_entry named its members; saying so is the difference between an input that is missing
        # and an input this tool cannot resolve.
        reads_as_group = isinstance(i.get("files"), int) and i["files"] > 1
        why = (
            f"is absent, and reads as a group of {i['files']} files recorded under a label before "
            "files_entry named its members (2026-10-02), so no file could be opened"
            if reads_as_group
            else "is absent"
        )
        out["entries"].append(
            {
                "path": path,
                "sha256_matches": False,
                "problem": "unresolved group" if reads_as_group else "absent",
            }
        )
        out["unavailable"].append(f"input {path} {why}")
        return
    try:
        digest, size, count = mf.sha256_of(p)
    except OSError as e:
        out["entries"].append({"path": path, "sha256_matches": False, "problem": f"cannot be read: {e}"})
        out["unavailable"].append(f"input {path} cannot be read: {e}")
        return
    if outside is not None:
        out["checked_outside_the_worktree"].append(outside)
    ok = digest == i.get("sha256")
    out["files_opened"] += count
    entry = {"path": path, "sha256_matches": ok, "bytes": size, "files_opened": count}
    if path in (local or set()):
        # opened and hashed, and counted here, whether or not it matched: a machine-local input whose
        # bytes differ is still an input this rebuild took from this machine, and the difference is
        # reported below. The count says what the rebuild rested on, not what passed.
        entry["machine_local"] = True
        out["machine_local_opened"] += 1
    out["entries"].append(entry)
    if not ok:
        out["unavailable"].append(f"input {path} has different bytes (sha256 {digest[:12]})")
        return
    out["checked"] += 1


def rebuild(
    result: Path, where: Path, venv: str = "fresh", keep: bool = False, inputs_only: bool = False
) -> dict[str, Any]:
    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    original = json.loads(result.read_text())
    m = original.get(mf.KEY)
    if not isinstance(m, dict):
        return {"rebuilt": False, "unavailable": [f"{result} carries no {mf.KEY}: nothing to rebuild from"]}
    code = m.get("code") or {}
    sha, argv = code.get("git_sha"), code.get("argv")
    report: dict[str, Any] = {
        "result": str(result),
        "git_sha": sha,
        "argv": argv,
        "dirty_code_at_write": code.get("dirty_code_paths", []),
        "unavailable": [],
    }
    if not sha or not argv:
        report["unavailable"].append("manifest has no git_sha or argv")
        return {**report, "rebuilt": False}
    wt = where / f"rebuild-{sha[:10]}"
    subprocess.run(["git", "-C", str(root), "worktree", "add", "-q", "--detach", str(wt), sha], check=True)
    links: dict[str, Any] = {"linked": []}
    try:
        for d in STORES:
            if (root / "data" / d).is_dir() and not (wt / "data" / d).exists():
                (wt / "data" / d).symlink_to(root / "data" / d)
        links = link_machine_local_inputs(wt, root, m.get("inputs"), output=f"data/results/{result.name}")
        found = check_inputs(wt, m.get("inputs"), root, machine_local=set(links["linked"]))
        report["inputs"] = found["entries"]
        report["inputs_declared"] = found["declared"]
        report["inputs_checked"] = found["checked"]
        report["inputs_unchecked"] = found["declared"] - found["checked"]
        report["files_opened_and_hashed"] = found["files_opened"]
        report["inputs_satisfied_from_machine_local_paths"] = found["machine_local_opened"]
        report["machine_local_inputs"] = {
            "count": found["machine_local_opened"],
            "of_declared": found["declared"],
            "linked_read_only_from": str(root / MACHINE_LOCAL_DIR),
            "how": "a clone at mode 0o444: the rebuild reads these bytes and a writer gets PermissionError",
            "paths": links["linked"],
            **{k: v for k, v in links.items() if k != "linked" and v},
        }
        report["unavailable"] += [
            f"a machine-local input could not be linked read-only: {f}" for f in links["failures"]
        ]
        if found["checked_outside_the_worktree"]:
            report["checked_outside_the_worktree"] = found["checked_outside_the_worktree"]
        report["unavailable"] += found["unavailable"]
        if report["unavailable"]:
            return {**report, "rebuilt": False}
        if inputs_only:
            # Steps 1 to 3 only. Item 10's capacity rule allows two lanes at once to read the per-element
            # cache, whose reader peaks near 2.9 GB, so a rebuild whose command reads it cannot be run by a
            # third lane. The inputs can still be linked and hashed, which costs a megabyte at a time, and
            # what that buys is the input figure and nothing else: no command ran, so nothing was compared.
            report["unavailable"].append(
                "--inputs-only: the command was not run, so no field was compared and no verdict is "
                "reported; the input figures above are all this says"
            )
            return {**report, "rebuilt": False}
        env = {k: v for k, v in os.environ.items() if k not in ("VIRTUAL_ENV", "UV_PROJECT_ENVIRONMENT")}
        if venv == "shared":
            env.update(UV_NO_SYNC="1", UV_PROJECT_ENVIRONMENT=str(root / ".venv"))
        else:
            sync = subprocess.run(
                ["uv", "sync", "--frozen", "--offline", "--quiet"],
                cwd=wt,
                env=env,
                capture_output=True,
                text=True,
            )
            report["venv"] = "fresh, from the committed uv.lock, offline"
            if sync.returncode:
                report["unavailable"].append(
                    f"fresh environment: uv sync --offline failed: {sync.stderr[-400:]}"
                )
                return {**report, "rebuilt": False}
        env["PYTHONPATH"] = str(wt)
        # The committed uv.lock is used as it is. Without UV_FROZEN `uv run` re-locks (the lock at
        # fe0880a still names genomeos 0.9.0 against pyproject's 1.0.0) and the rewritten lock makes
        # the second checkout dirty, which the rebuilt manifest then records as a one-byte difference.
        env.update(UV_FROZEN="1", UV_OFFLINE="1")
        run = subprocess.run(["uv", "run", "python", *argv], cwd=wt, env=env, capture_output=True, text=True)
        report["exit"] = run.returncode
        report["stdout_tail"] = run.stdout[-1500:]
        if run.returncode:
            report["unavailable"].append(f"the command failed: {run.stderr[-800:]}")
            return {**report, "rebuilt": False}
        name = original.get("result") or result.stem
        written = wt / "data" / "results" / f"{name}.json"
        if written.name != result.name:  # the manifest was read from a copy under another name
            shutil.copyfile(written, wt / "data" / "results" / result.name)
        rebuilt = json.loads((wt / "data" / "results" / result.name).read_text())
        report["rebuilt_sha256"] = mf.sha256_of(written)[0]
        found = diff(comparable(original), comparable(rebuilt))
        real, env = environment_differences(found)
        report["differences"] = real
        report["environment_fields_ignored"] = env
        report["must_hold_failures"] = must_hold_failures(original, rebuilt)
        report["differences"] += report["must_hold_failures"]
        report["timing_fields_ignored"] = sorted(set(timing_paths(original)) | set(timing_paths(rebuilt)))
        report["resource_fields_ignored"] = sorted(
            set(resource_paths(original)) | set(resource_paths(rebuilt))
        )
        report["fields_compared"] = len(comparable(original))
        report["leaves"] = leaf_reconciliation(original)
        report["leaves_compared"] = report["leaves"]["compared"]
        # despite its name this key compares fields (date and the code block ignored); the next is bytes
        # (timing keys are ignored as well, and listed in timing_fields_ignored)
        why = no_verdict_reason(report["inputs_declared"], report["inputs_checked"])
        if why:
            # unreachable through the branch above; asked again so that no later edit can write a verdict
            # for a result whose inputs were not all opened, the failure this tool had on 2026-10-02
            report["unavailable"].append(why)
            return {**report, "rebuilt": False}
        report["identical_bytes_except_date_and_run"] = not report["differences"]
        report["identical_bytes"] = written.read_bytes() == result.read_bytes()
        return {**report, "rebuilt": True}
    finally:
        if not keep:
            unlock_links(wt, links["linked"])
            subprocess.run(["git", "-C", str(root), "worktree", "remove", "--force", str(wt)], check=False)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("result", type=Path)
    ap.add_argument("--where", type=Path, required=True, help="directory for the second checkout")
    ap.add_argument("--venv", choices=("fresh", "shared"), default="fresh")
    ap.add_argument("--keep", action="store_true", help="leave the worktree in place")
    ap.add_argument(
        "--inputs-only",
        action="store_true",
        help="link and hash every declared input, then stop: no command is run and no verdict is reported",
    )
    args = ap.parse_args(argv)
    r = rebuild(args.result.resolve(), args.where.resolve(), args.venv, args.keep, args.inputs_only)
    json.dump({k: v for k, v in r.items() if k != "stdout_tail"}, sys.stdout, indent=1)
    print()
    return 0 if r.get("rebuilt") and not r.get("differences") else 1


if __name__ == "__main__":
    raise SystemExit(main())
