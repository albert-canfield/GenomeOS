# SPDX-License-Identifier: AGPL-3.0-or-later
"""Which committed results were written by a script outside this repository. READ-ONLY.

    uv run python scripts/entry_script_census.py              # the README-cited population
    uv run python scripts/entry_script_census.py --all        # every data/results/*.json as well
    uv run python scripts/entry_script_census.py --json

Why (2026-10-02, lane-entrypoint). `genomeos.manifest.code_cleanliness` compared git's dirty list
against the paths the CALLER DECLARED and never looked at the script that was running, so a result
could assert `own_code_is_committed: True` while the bytes came from a file in no commit. Measured on
`data/results/organised_chr21.json`
(working tree at d4d0449, 19:41; superseded by e1def2e, 20:41): `code.argv[0]` was
`/private/tmp/.../scratchpad/write.py`, the
counting path holds 39 `genomeos/` modules and no script, and the manifest says the code was
committed. The check is fixed for new writes; this census says how far the old ones reach.

NOTHING IS REWRITTEN AND NO RESULT'S BYTES ARE TOUCHED. Several results' sha256 ARE registrations in
`data/results/manifest_headlines.json`, and rewriting one would break a pin. The census reads and
reports; what to do about what it finds is not its decision.

HOW THE VERDICT IS READ, and it is read from the record rather than recomputed. `manifest.code.argv`
is written by `manifest._argv`, which makes argv[0] RELATIVE to the repository root when it is under
it and leaves it absolute when it is not. So an absolute argv[0] is itself the evidence that the
script was outside the repository at write time -- no present-tree lookup is involved, and a script
that has since been committed, moved or deleted cannot change the answer. The one case where that
inference does not hold is named rather than folded in: an absolute path that lies UNDER this
checkout's root was inside the repository and recorded absolutely (a writer that passed a resolved
`__file__` before `_argv` relativised, or a differently-rooted checkout), and it is counted apart
under `inside_recorded_absolute` instead of being read as outside.

"README-cited" is the definition `tests/test_headline_registry_matches_readme.py` already uses, not a
new one: the entries of `manifest_headlines.json` (`rebuilt` + `pending`), which is the record of what
the README's claims rest on, together with every `*.json` the README names that is a file in
data/results. Both halves are reported, and a registered name with no file on this disk is reported
as `absent` rather than skipped.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
from pathlib import Path
from typing import Any

RESULTS = Path("data/results")
README = Path("README.md")
HEADLINES = RESULTS / "manifest_headlines.json"
ROOT = Path(__file__).resolve().parent.parent

#: The verdict for each recorded argv[0]. `outside_repository` is the one the census was asked for.
VERDICTS = {
    "outside_repository": "argv[0] is an absolute path that is not under this checkout: the script "
    "that wrote the bytes was outside the repository and no commit holds it",
    "inside_recorded_absolute": "argv[0] is absolute but under this checkout's root, so it was "
    "inside the repository and merely recorded absolutely; NOT counted as outside",
    "inside_repository": "argv[0] is a repository-relative path, which is what manifest._argv writes "
    "for a script under the root",
    "inline_source": "argv[0] is '-c': the source was never in a file",
    "stdin": "argv[0] is '-': the source was piped in",
    "test_runner": "argv[0] is a test runner's entry point",
    "unrecorded": "the manifest records no code.argv, so nothing says which script ran",
    "absent": "the result is registered as README-cited but is not a file on this disk",
}

_RUNNERS = ("pytest", "py.test", "unittest")


def readme_cited() -> tuple[list[str], dict[str, list[str]]]:
    """The README-cited result names, and where each citation comes from."""
    where: dict[str, list[str]] = {}
    try:
        record = json.loads(HEADLINES.read_text())
    except OSError:
        record = {}
    for entry in list(record.get("rebuilt") or []) + list(record.get("pending") or []):
        name = entry.get("result")
        if name:
            where.setdefault(name, []).append(f"manifest_headlines ({entry.get('quoted_in')})")
    try:
        text = README.read_text()
    except OSError:
        text = ""
    for filename in sorted(set(re.findall(r"[a-z0-9_]+\.json", text))):
        if (RESULTS / filename).is_file():
            where.setdefault(filename.removesuffix(".json"), []).append("README.md names the file")
    return sorted(where), where


def verdict_for(argv0: Any) -> str:
    if argv0 is None:
        return "unrecorded"
    argv0 = str(argv0)
    if argv0 == "-c":
        return "inline_source"
    if argv0 == "-":
        return "stdin"
    if Path(argv0).name in _RUNNERS or argv0.endswith(("pytest/__main__.py", "unittest/__main__.py")):
        return "test_runner"
    if not Path(argv0).is_absolute():
        return "inside_repository"
    try:
        Path(argv0).relative_to(ROOT)
    except ValueError:
        return "outside_repository"
    return "inside_recorded_absolute"


def bytes_state(p: Path) -> str:
    """Whether the result FILE's own bytes are committed: `committed`, `modified`, `untracked`, `unknown`.

    A finding needs this column. The 24 `organised_chr*.json` files this census flags are WORKING-TREE
    modifications whose committed versions are older, legacy-shaped results with no `result_manifest`
    at all, and another lane is regenerating them from a committed entry point. A file about to be
    discarded is not a published result that certified itself falsely, and a census that printed the
    two the same way would invite its reader to count 24 published breaches where there are none.
    """
    try:
        listed = subprocess.run(
            ["git", "ls-files", "-z", "--", str(p)], capture_output=True, text=True, timeout=30, check=True
        ).stdout
        if str(p) not in [x for x in listed.split("\0") if x]:
            return "untracked"
        status = subprocess.run(
            ["git", "status", "--porcelain", "-z", "--", str(p)],
            capture_output=True,
            text=True,
            timeout=30,
            check=True,
        ).stdout
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    return "modified" if [x for x in status.split("\0") if x] else "committed"


def read_one(name: str) -> dict[str, Any]:
    """One result's recorded entry script. Reads the file and nothing else; never writes."""
    p = RESULTS / f"{name}.json"
    if not p.is_file():
        return {"result": name, "verdict": "absent", "argv0": None}
    try:
        d = json.loads(p.read_text())
    except (OSError, ValueError) as e:
        return {"result": name, "verdict": "unrecorded", "argv0": None, "unreadable": str(e)}
    m = d.get("result_manifest") or {}
    code = m.get("code") or {}
    argv = code.get("argv")
    argv0 = argv[0] if isinstance(argv, list) and argv else None
    clean = m.get("code_cleanliness") or {}
    return {
        "result": name,
        "verdict": verdict_for(argv0),
        "argv0": argv0,
        "claims_own_code_is_committed": clean.get("own_code_is_committed"),
        "carries_entry_script": isinstance(clean.get("entry_script"), dict),
        "counting_path_count": clean.get("counting_path_count"),
        "bytes_state": bytes_state(p),
    }


def census(names: list[str]) -> dict[str, Any]:
    rows = [read_one(n) for n in names]
    counts: dict[str, int] = {}
    for row in rows:
        counts[row["verdict"]] = counts.get(row["verdict"], 0) + 1
    outside = [r for r in rows if r["verdict"] == "outside_repository"]
    return {
        "examined": len(rows),
        "counts": counts,
        "outside_repository": outside,
        "outside_repository_count": len(outside),
        "falsely_certified": [r for r in outside if r.get("claims_own_code_is_committed") is True],
        "falsely_certified_and_committed": [
            r
            for r in outside
            if r.get("claims_own_code_is_committed") is True and r.get("bytes_state") == "committed"
        ],
        "rows": rows,
        "verdicts_are": VERDICTS,
        "nothing_was_written": "this script reads result files and writes none of them",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--all", action="store_true", help="every data/results/*.json as a second table")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    names, where = readme_cited()
    cited = census(names)
    everything = census(sorted(p.stem for p in RESULTS.glob("*.json"))) if args.all else None
    if args.json:
        print(json.dumps({"readme_cited": cited, "all_results": everything, "citations": where}, indent=1))
        return 0

    print(f"README-cited results examined: {cited['examined']}")
    for verdict, n in sorted(cited["counts"].items()):
        print(f"  {verdict:28s} {n}")
    print(f"\nargv[0] OUTSIDE the repository: {cited['outside_repository_count']}")
    for row in cited["outside_repository"]:
        print(f"  {row['result']}")
        print(f"    argv[0]               {row['argv0']}")
        print(f"    own_code_is_committed {row['claims_own_code_is_committed']}")
        print(f"    the result's bytes   {row['bytes_state']}")
        print(f"    cited in              {'; '.join(where.get(row['result'], []))}")
    if cited["falsely_certified"]:
        print(
            f"\nof those, asserting own_code_is_committed True: "
            f"{len(cited['falsely_certified'])} "
            f"({', '.join(r['result'] for r in cited['falsely_certified'])})"
        )
    if everything is not None:
        print(f"\nevery committed result ({everything['examined']} files), for scale:")
        for verdict, n in sorted(everything["counts"].items()):
            print(f"  {verdict:28s} {n}")
        print(f"  outside and asserting own_code_is_committed True: {len(everything['falsely_certified'])}")
        print(
            f"  of those, with the result's own bytes COMMITTED: "
            f"{len(everything['falsely_certified_and_committed'])}"
        )
        for row in everything["falsely_certified"]:
            print(f"    {row['result']:28s} bytes {row['bytes_state']:10s} {row['argv0']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
