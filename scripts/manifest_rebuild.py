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
   pre-push hook does, and check every input's sha256 against the manifest, naming any input that
   is absent or has different bytes. Every declared input is accounted for: a group of files read
   together is opened and hashed member by member, and an input the tool cannot resolve stops the
   rebuild by name. `inputs_declared`, `inputs_checked` and `inputs_unchecked` carry the denominator,
   so an input that was not opened cannot leave the list and read as one that matched;
4. run `code.argv` in the worktree, in a fresh environment from the committed uv.lock (`--venv
   fresh`, offline) or the checkout's own (`--venv shared`);
5. compare the rebuilt result with RESULT.json field by field, ignoring only `date` and the
   `code` block of the manifest (which records the run, not the result).
   Since 2026-09-28 wall-clock timing keys are ignored too, at any depth (`seconds`, `*_seconds`,
   `seconds_*`, `*_per_second`); the report lists their paths under `timing_fields_ignored`.
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


def _strip_timing(x: Any) -> Any:
    if isinstance(x, dict):
        return {k: _strip_timing(v) for k, v in x.items() if not is_timing(k)}
    if isinstance(x, list):
        return [_strip_timing(v) for v in x]
    return x


def timing_paths(x: Any, where: str = "") -> list[str]:
    """Every path at which comparable() drops a timing key."""
    if isinstance(x, dict):
        out = []
        for k, v in x.items():
            out.extend([f"{where}/{k}"] if is_timing(k) else timing_paths(v, f"{where}/{k}"))
        return out
    if isinstance(x, list):
        return [p for i, v in enumerate(x) for p in timing_paths(v, f"{where}[{i}]")]
    return []


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


def check_inputs(worktree: Path, inputs: Any, root: Path) -> dict[str, Any]:
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
        "entries": [],
        "unavailable": [],
        "checked_outside_the_worktree": [],
    }
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
            _check_group(worktree, i, out)
        elif i.get("sha256") == NOT_HASHED:
            out["entries"].append(
                {"path": path, "sha256_matches": False, "problem": f"sha256 is {NOT_HASHED!r}"}
            )
            out["unavailable"].append(
                f"input {path} records sha256 {NOT_HASHED!r}: nothing was hashed when the result was "
                "written, so there is nothing for a rebuild to check these bytes against"
            )
        else:
            _check_one_path(worktree, i, linked, out)
    return out


def _check_group(worktree: Path, i: dict[str, Any], out: dict[str, Any]) -> None:
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
    out["entries"].append(
        {
            "path": path,
            "group": True,
            "members_declared": len(names),
            "members_opened_and_hashed": len(seen),
            "members_with_different_bytes": differing,
            "sha256_matches": ok and not differing,
            "bytes": size,
        }
    )
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


def _check_one_path(worktree: Path, i: dict[str, Any], linked: list[str], out: dict[str, Any]) -> None:
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
    out["entries"].append({"path": path, "sha256_matches": ok, "bytes": size, "files_opened": count})
    if not ok:
        out["unavailable"].append(f"input {path} has different bytes (sha256 {digest[:12]})")
        return
    out["checked"] += 1


def rebuild(result: Path, where: Path, venv: str = "fresh", keep: bool = False) -> dict[str, Any]:
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
    try:
        for d in STORES:
            if (root / "data" / d).is_dir() and not (wt / "data" / d).exists():
                (wt / "data" / d).symlink_to(root / "data" / d)
        found = check_inputs(wt, m.get("inputs"), root)
        report["inputs"] = found["entries"]
        report["inputs_declared"] = found["declared"]
        report["inputs_checked"] = found["checked"]
        report["inputs_unchecked"] = found["declared"] - found["checked"]
        report["files_opened_and_hashed"] = found["files_opened"]
        if found["checked_outside_the_worktree"]:
            report["checked_outside_the_worktree"] = found["checked_outside_the_worktree"]
        report["unavailable"] += found["unavailable"]
        if report["unavailable"]:
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
        report["fields_compared"] = len(comparable(original))
        report["leaves"] = leaf_reconciliation(original)
        report["leaves_compared"] = report["leaves"]["compared"]
        # despite its name this key compares fields (date and the code block ignored); the next is bytes
        # (timing keys are ignored as well, and listed in timing_fields_ignored)
        if report["inputs_checked"] != report["inputs_declared"] or not report["inputs_declared"]:
            # unreachable through the branch above; kept so that no later edit can write a verdict for a
            # result whose inputs were not all opened, which is the failure this tool had on 2026-10-02
            report["unavailable"].append(
                f"{report['inputs_checked']} of {report['inputs_declared']} declared inputs were opened "
                "and hashed: no verdict is reported"
            )
            return {**report, "rebuilt": False}
        report["identical_bytes_except_date_and_run"] = not report["differences"]
        report["identical_bytes"] = written.read_bytes() == result.read_bytes()
        return {**report, "rebuilt": True}
    finally:
        if not keep:
            subprocess.run(["git", "-C", str(root), "worktree", "remove", "--force", str(wt)], check=False)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("result", type=Path)
    ap.add_argument("--where", type=Path, required=True, help="directory for the second checkout")
    ap.add_argument("--venv", choices=("fresh", "shared"), default="fresh")
    ap.add_argument("--keep", action="store_true", help="leave the worktree in place")
    args = ap.parse_args(argv)
    r = rebuild(args.result.resolve(), args.where.resolve(), args.venv, args.keep)
    json.dump({k: v for k, v in r.items() if k != "stdout_tail"}, sys.stdout, indent=1)
    print()
    return 0 if r.get("rebuilt") and not r.get("differences") else 1


if __name__ == "__main__":
    raise SystemExit(main())
