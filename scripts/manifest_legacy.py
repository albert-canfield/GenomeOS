# SPDX-License-Identifier: AGPL-3.0-or-later
"""The legacy allowlist of the result registry, generated once from git history (item 12 S6).

    uv run python scripts/manifest_legacy.py            # write data/results_legacy.txt; refuses if it exists
    uv run python scripts/manifest_legacy.py --check    # compare the committed list with git history

`genomeos/results.py` lets a result without a complete manifest into data/results/ only when its
name is on this list: the names that were in the registry before the contract arrived (fe0880a,
`manifest.CONTRACT_COMMIT`). Until item 12 S6 "historical" meant "the file already exists", so a
retry of a failed new result found its own file and passed. The list is the census's historical
set (`manifest_census.historical_names`): every data/results/*.json in the tree of the contract
commit's parent ("tracked"), and every git-ignored data/results/*.json on the generating machine
last written before the contract commit ("ignored-local").

It is generated once and never extended by code: this script refuses to overwrite it, and a name
is added only by a reviewed commit that says why. `--check` reports any difference from history;
an ignored-local name absent from this machine is reported as not checkable here, not as wrong.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from manifest_census import ROOT, historical_names  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos.results import LEGACY_ALLOWLIST, LEGACY_COUNT  # noqa: E402

SOURCES = ("tracked", "ignored-local")


def render(names: dict[str, str], commit: str = mf.CONTRACT_COMMIT) -> str:
    counts = {s: sum(v == s for v in names.values()) for s in SOURCES}
    head = [
        "# Legacy allowlist of the result registry (item 12 S6), generated 2026-09-28 by",
        "# scripts/manifest_legacy.py: the result names in data/results/ before the manifest contract",
        f"# arrived at {commit[:7]}. {counts['tracked']} tracked in the tree of {commit[:7]}^,",
        f"# {counts['ignored-local']} git-ignored on the generating machine and last written before it.",
        "# genomeos/results.py lets a result without a complete manifest into the registry only when",
        "# its name is here (it warns); anything else is quarantined. Generated once from git history",
        "# and never extended by code: add a name only in a reviewed commit that says why.",
        "# docs/DATA.md, 'The result registry: what enters it and what is quarantined'.",
        "# name<TAB>source",
    ]
    return "\n".join(head + [f"{n}\t{names[n]}" for n in sorted(names)]) + "\n"


def read(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            name, _, source = line.partition("\t")
            out[name.strip()] = source.strip()
    return out


def check(path: Path, root: Path = ROOT) -> list[str]:
    """Differences between the committed list and git history, as sentences; empty when they agree."""
    listed = read(path)
    hist = historical_names(root)
    problems = []
    if len(listed) != LEGACY_COUNT:
        problems.append(f"{len(listed)} names listed, {LEGACY_COUNT} registered")
    for n, s in sorted(listed.items()):
        if s not in SOURCES:
            problems.append(f"{n}: unknown source {s!r}")
        elif s == "tracked" and hist.get(n) != "tracked":
            problems.append(f"{n}: listed as tracked, not in the tree before the contract")
        elif s == "ignored-local" and n in hist and hist[n] != "ignored-local":
            problems.append(f"{n}: listed as ignored-local, history says {hist[n]}")
    for n, s in sorted(hist.items()):
        if n not in listed:
            problems.append(f"{n}: historical ({s}) but not listed")
    return problems


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=ROOT / LEGACY_ALLOWLIST)
    ap.add_argument("--check", action="store_true", help="compare the committed list with git history")
    args = ap.parse_args(argv)
    if args.check:
        listed = read(args.out)
        problems = check(args.out)
        here = historical_names()
        absent = sum(1 for n, s in listed.items() if s == "ignored-local" and n not in here)
        print(f"{len(listed)} names listed; {len(problems)} differences from git history", end="")
        print(f"; {absent} ignored-local names not on this machine (not checkable here)")
        for p in problems:
            print(f"  {p}")
        return 1 if problems else 0
    if args.out.exists():
        print(f"{args.out} exists: the allowlist is generated once and never rewritten", file=sys.stderr)
        return 2
    names = historical_names()
    if len(names) != LEGACY_COUNT:
        print(f"history gives {len(names)} names, {LEGACY_COUNT} registered: not written", file=sys.stderr)
        return 3
    args.out.write_text(render(names))
    print(f"{args.out}: {len(names)} names")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
