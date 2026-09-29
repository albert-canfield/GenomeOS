# SPDX-License-Identifier: AGPL-3.0-or-later
"""Which results already carry which part of the provenance contract (review item R9).

    uv run python scripts/manifest_census.py            # the table
    uv run python scripts/manifest_census.py --json     # per file, for tooling
    uv run python scripts/manifest_census.py --missing assembly   # the files lacking one field

Reads every data/results/*.json through the tolerant reader (genomeos/manifest.py) and counts, per
required field, how many results declare it in a manifest, how many carry an older top-level key
that covers part of it ("legacy"), and how many say nothing. Nothing is rewritten: the historical
results stay as they were written. A legacy key is weak evidence (a `genome` key may name a FASTA
path, not a build), so the declared column is the only one that meets the contract.

    uv run python scripts/manifest_census.py --enforcement   # review item 12 S6, from git history

Item 12 S6 (second external review, 2026-09-28) found that a new result failing the contract was
written before it was refused, and that a retry then found the file, took it for historical and only
warned. `--enforcement` counts what that could have let in: every caller of `save_result` and of the
manifest stamp (read from the source, not assumed), which result names were historical when the
contract arrived (`CONTRACT_COMMIT`: the tree before it, plus the git-ignored results on this disk
written before it), which names came after, and whether any result, on disk or in any commit since,
is a new name without a complete manifest, which is what the retry path leaves behind.
It also lists writers that write into data/results without `save_result`, found by a pattern (a
lower bound: a writer that builds its path another way is not seen).
"""

from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
from pathlib import Path

from genomeos import manifest as mf
from genomeos.results import RESULTS_DIR, legacy_names, quarantine_dir


def census(results_dir: Path = RESULTS_DIR) -> dict:
    per_file: dict[str, dict] = {}
    unreadable: list[str] = []
    for p in sorted(results_dir.glob("*.json")):
        try:
            payload = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError, UnicodeDecodeError):
            unreadable.append(p.name)
            continue
        r = mf.read(payload)
        per_file[p.stem] = {"declared": r["declared"], "complete": r["complete"], "fields": r["fields"]}
    fields = {
        f: {
            "declared": sum(1 for v in per_file.values() if v["fields"][f] == "declared"),
            "legacy": sum(1 for v in per_file.values() if v["fields"][f] == "legacy"),
            "none": sum(1 for v in per_file.values() if v["fields"][f] is None),
        }
        for f in mf.REQUIRED
    }
    legacy = legacy_names()
    unlisted = [n for n in per_file if n not in legacy]
    qdir = quarantine_dir(results_dir)
    return {
        "results": len(per_file),
        "legacy_allowlist": len(legacy),
        "on_allowlist": len(per_file) - len(unlisted),
        "not_on_allowlist": len(unlisted),
        "not_on_allowlist_incomplete": sorted(n for n in unlisted if not per_file[n]["complete"]),
        "quarantined": sorted(p.stem for p in qdir.glob("*.json")) if qdir.is_dir() else [],
        "unreadable": unreadable,
        "with_manifest": sum(1 for v in per_file.values() if v["declared"]),
        "complete": sum(1 for v in per_file.values() if v["complete"]),
        "no_field_at_all": sum(1 for v in per_file.values() if not any(v["fields"].values())),
        "fields": fields,
        "files": per_file,
    }


# --- item 12 S6: enforcement census, read from git history -------------------------------------------

#: fe0880a, the commit that introduced the contract (review R9, lane-manifest, 2026-09-28 00:12 +0100).
#: A result name is historical when it was in the registry before this commit.
CONTRACT_COMMIT = mf.CONTRACT_COMMIT
ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOTS = ("genomeos", "scripts", "tests")


def _git(*args: str, root: Path = ROOT) -> str:
    return subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True, check=True).stdout


def historical_names(root: Path = ROOT, commit: str = CONTRACT_COMMIT) -> dict[str, str]:
    """Result names in the registry before the contract, each with where it was found: "tracked" (a
    data/results/*.json in the tree of the contract commit's parent) or "ignored-local" (a git-ignored
    data/results/*.json on this disk, last written before the contract commit)."""
    before = f"{commit}^"
    out: dict[str, str] = {}
    for line in _git("ls-tree", "--name-only", before, "data/results/", root=root).splitlines():
        if line.endswith(".json"):
            out[Path(line).stem] = "tracked"
    t0 = int(_git("show", "-s", "--format=%ct", commit, root=root).strip())
    ignored = _git("ls-files", "--others", "--ignored", "--exclude-standard", "data/results/", root=root)
    for line in ignored.splitlines():
        p = root / line
        if line.endswith(".json") and "/" not in line[len("data/results/") :] and p.stat().st_mtime < t0:
            out.setdefault(p.stem, "ignored-local")
    return out


def _calls(tree: ast.AST, names: set[str]) -> list[ast.Call]:
    out = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            f = node.func
            name = f.id if isinstance(f, ast.Name) else f.attr if isinstance(f, ast.Attribute) else None
            if name in names:
                out.append(node)
    return out


def _results_targets(tree: ast.AST, text: str) -> set[str]:
    """Names assigned a path that points into data/results (RESULTS_DIR / ..., "data/results/...",
    ... / "data" / "results" / ...)."""
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign):
            src = ast.get_source_segment(text, node.value) or ""
            if "RESULTS_DIR /" in src or "data/results/" in src or '"results" /' in src:
                out |= {t.id for t in node.targets if isinstance(t, ast.Name)}
    return out


def callers(root: Path = ROOT) -> dict:
    """Every call of save_result, of the manifest stamp and code_revision in tracked Python source, and
    the writers that write into data/results directly (pattern: a name assigned a path into
    data/results, then .write_text/.write_bytes on it, or open(name, "w"))."""
    files = [f for f in _git("ls-files", "--", *(f"{d}/*.py" for d in SOURCE_ROOTS), root=root).splitlines()]
    save_sites: list[dict] = []
    stamp_sites: list[str] = []
    direct: list[str] = []
    for f in files:
        if f == "genomeos/results.py":
            continue
        text = (root / f).read_text()
        if (
            "save_result" not in text
            and "stamp" not in text
            and "code_revision" not in text
            and "results" not in text
        ):
            continue
        try:
            tree = ast.parse(text)
        except SyntaxError:
            continue
        mentions_key = "result_manifest" in text or "mf.KEY" in text or "manifest.KEY" in text
        for c in _calls(tree, {"save_result"}):
            kw = {k.arg for k in c.keywords}
            save_sites.append(
                {
                    "file": f,
                    "line": c.lineno,
                    "manifest_argument": "manifest" in kw or len(c.args) >= 4,
                    "file_mentions_result_manifest": mentions_key,
                    "strict": next((ast.unparse(k.value) for k in c.keywords if k.arg == "strict"), None),
                    "results_dir_given": "results_dir" in kw or len(c.args) >= 3,
                }
            )
        if f != "genomeos/manifest.py":
            for c in _calls(tree, {"stamp", "code_revision"}):
                fn = c.func
                owner = ast.unparse(fn.value) if isinstance(fn, ast.Attribute) else ""
                if owner in ("mf", "manifest", "genomeos.manifest") or isinstance(fn, ast.Name):
                    stamp_sites.append(f"{f}:{c.lineno}")
        targets = _results_targets(tree, text)
        if targets and f.startswith(("scripts/", "genomeos/")):
            for c in _calls(tree, {"write_text", "write_bytes", "open"}):
                fn = c.func
                hit = (
                    isinstance(fn, ast.Attribute)
                    and fn.attr in ("write_text", "write_bytes")
                    and isinstance(fn.value, ast.Name)
                    and fn.value.id in targets
                ) or (
                    isinstance(fn, ast.Name)
                    and fn.id == "open"
                    and c.args
                    and isinstance(c.args[0], ast.Name)
                    and c.args[0].id in targets
                    and any(isinstance(a, ast.Constant) and "w" in str(a.value) for a in c.args[1:])
                )
                if hit:
                    direct.append(f"{f}:{c.lineno}")
    by_dir: dict[str, int] = {}
    for s in save_sites:
        by_dir[s["file"].split("/")[0]] = by_dir.get(s["file"].split("/")[0], 0) + 1
    return {
        "save_result_sites": len(save_sites),
        "save_result_files": len({s["file"] for s in save_sites}),
        "save_result_sites_by_dir": by_dir,
        "sites_passing_manifest_argument": sum(s["manifest_argument"] for s in save_sites),
        "sites_in_files_naming_result_manifest": sum(
            s["file_mentions_result_manifest"] and not s["manifest_argument"] for s in save_sites
        ),
        "sites_passing_strict": [
            f"{s['file']}:{s['line']} strict={s['strict']}" for s in save_sites if s["strict"]
        ],
        "sites_with_results_dir": sum(s["results_dir_given"] for s in save_sites),
        "stamp_or_code_revision_sites": stamp_sites,
        "direct_writers_into_results": direct,
        "sites": save_sites,
    }


def enforcement(root: Path = ROOT, commit: str = CONTRACT_COMMIT) -> dict:
    """What the write-then-refuse and retry-warns path could have let into the registry."""
    hist = historical_names(root, commit)
    on_disk = {p.stem: p for p in sorted((root / "data" / "results").glob("*.json"))}
    new = sorted(n for n in on_disk if n not in hist)
    status: dict[str, str] = {}
    for n, p in on_disk.items():
        r = mf.read(json.loads(p.read_text()))
        status[n] = "complete" if r["complete"] else "declared_incomplete" if r["declared"] else "none"
    since = f"{commit}^..HEAD"
    added = sorted(
        {
            Path(x).stem
            for x in _git(
                "log",
                "--format=",
                "--name-only",
                "--diff-filter=A",
                since,
                "--",
                "data/results/*.json",
                root=root,
            ).splitlines()
            if x.strip()
        }
    )
    versions = incomplete_versions = 0
    incomplete_committed: list[str] = []
    for n in added:
        if n in hist:
            continue
        path = f"data/results/{n}.json"
        for c in _git("log", "--format=%H", since, "--", path, root=root).split():
            try:
                blob = _git("show", f"{c}:{path}", root=root)
            except subprocess.CalledProcessError:
                continue  # the commit deleted it
            versions += 1
            if not mf.read(json.loads(blob))["complete"]:
                incomplete_versions += 1
                incomplete_committed.append(f"{c[:7]} {n}")
    return {
        "contract_commit": commit,
        "historical": len(hist),
        "historical_tracked": sum(v == "tracked" for v in hist.values()),
        "historical_ignored_local": sum(v == "ignored-local" for v in hist.values()),
        "on_disk": len(on_disk),
        "on_disk_historical": sum(n in hist for n in on_disk),
        "on_disk_new": len(new),
        "on_disk_new_complete": sum(status[n] == "complete" for n in new),
        "on_disk_new_incomplete": [n for n in new if status[n] != "complete"],
        "on_disk_historical_with_manifest": sum(status[n] != "none" for n in on_disk if n in hist),
        "on_disk_historical_incomplete_manifest": sorted(
            n for n in on_disk if n in hist and status[n] == "declared_incomplete"
        ),
        "new_names_added_since_contract": sum(n not in hist for n in added),
        "committed_versions_of_new_names": versions,
        "committed_versions_incomplete": incomplete_versions,
        "committed_incomplete": incomplete_committed,
        "untracked_not_ignored_results": [
            x
            for x in _git(
                "ls-files", "--others", "--exclude-standard", "data/results/", root=root
            ).splitlines()
            if x.endswith(".json")
        ],
    }


def print_enforcement(e: dict, c: dict) -> None:
    say = print
    say(f"callers: save_result called at {c['save_result_sites']} sites in {c['save_result_files']} files")
    say(f"  by directory {c['save_result_sites_by_dir']}")
    say(f"  {c['sites_passing_manifest_argument']} pass a manifest argument")
    say(f"  {c['sites_in_files_naming_result_manifest']} more sit in a file naming result_manifest")
    say(f"  {c['sites_with_results_dir']} pass a results directory")
    say(f"  strict passed at: {c['sites_passing_strict'] or 'none'}")
    say(f"  stamp or code_revision called elsewhere: {c['stamp_or_code_revision_sites'] or 'none'}")
    direct = c["direct_writers_into_results"]
    say(f"  writes into data/results without save_result (pattern candidates): {len(direct)}")
    for d in direct:
        say(f"    {d}")
    say(f"registry: {e['on_disk']} results on disk")
    say(f"  {e['historical']} names historical at {e['contract_commit'][:7]}:")
    say(f"    {e['historical_tracked']} tracked in its parent")
    say(f"    {e['historical_ignored_local']} git-ignored on this disk, last written before it")
    say(f"  on disk: {e['on_disk_historical']} historical, {e['on_disk_new']} new")
    say(f"  new with a complete manifest: {e['on_disk_new_complete']}")
    say(f"  new without one: {len(e['on_disk_new_incomplete'])} {e['on_disk_new_incomplete'] or ''}")
    say(f"  historical rewritten with a manifest: {e['on_disk_historical_with_manifest']}")
    say(f"    incomplete (the warn path): {e['on_disk_historical_incomplete_manifest']}")
    say(f"  new names committed since the contract: {e['new_names_added_since_contract']}")
    say(f"    their committed versions: {e['committed_versions_of_new_names']}")
    say(f"    incomplete: {e['committed_versions_incomplete']} {e['committed_incomplete'] or ''}")
    say(f"  untracked, not ignored results: {e['untracked_not_ignored_results'] or 'none'}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, default=RESULTS_DIR)
    ap.add_argument("--json", action="store_true", help="print the whole census as JSON")
    ap.add_argument("--missing", choices=mf.REQUIRED, help="list the results that say nothing on one field")
    ap.add_argument(
        "--enforcement", action="store_true", help="item 12 S6: callers and the registry since the contract"
    )
    args = ap.parse_args(argv)
    if args.enforcement:
        e, c = enforcement(), callers()
        if args.json:
            json.dump(
                {"enforcement": e, "callers": {k: v for k, v in c.items() if k != "sites"}},
                sys.stdout,
                indent=1,
            )
            print()
        else:
            print_enforcement(e, c)
        return 0
    c = census(args.dir)
    if args.json:
        json.dump(c, sys.stdout, indent=1)
        print()
        return 0
    if args.missing:
        for name, v in c["files"].items():
            if v["fields"][args.missing] is None:
                print(name)
        return 0
    print(
        f"{c['results']} results read ({len(c['unreadable'])} unreadable); "
        f"{c['with_manifest']} carry a manifest, {c['complete']} meet the contract, "
        f"{c['no_field_at_all']} say nothing on any field"
    )
    bad, quarantined = c["not_on_allowlist_incomplete"], c["quarantined"]
    print(
        f"legacy allowlist {c['legacy_allowlist']} names: {c['on_allowlist']} results on it, "
        f"{c['not_on_allowlist']} not, of which {len(bad)} incomplete{': ' + ', '.join(bad) if bad else ''}; "
        f"{len(quarantined)} quarantined{': ' + ', '.join(quarantined) if quarantined else ''}"
    )
    print(f"{'field':<12} {'declared':>9} {'legacy':>7} {'none':>6}")
    for f, n in c["fields"].items():
        print(f"{f:<12} {n['declared']:>9} {n['legacy']:>7} {n['none']:>6}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
