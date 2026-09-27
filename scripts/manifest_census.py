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
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from genomeos import manifest as mf
from genomeos.results import RESULTS_DIR


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
    return {
        "results": len(per_file),
        "unreadable": unreadable,
        "with_manifest": sum(1 for v in per_file.values() if v["declared"]),
        "complete": sum(1 for v in per_file.values() if v["complete"]),
        "no_field_at_all": sum(1 for v in per_file.values() if not any(v["fields"].values())),
        "fields": fields,
        "files": per_file,
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", type=Path, default=RESULTS_DIR)
    ap.add_argument("--json", action="store_true", help="print the whole census as JSON")
    ap.add_argument("--missing", choices=mf.REQUIRED, help="list the results that say nothing on one field")
    args = ap.parse_args(argv)
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
    print(f"{'field':<12} {'declared':>9} {'legacy':>7} {'none':>6}")
    for f, n in c["fields"].items():
        print(f"{f:<12} {n['declared']:>9} {n['legacy']:>7} {n['none']:>6}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
