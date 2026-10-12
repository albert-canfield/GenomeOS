# SPDX-License-Identifier: AGPL-3.0-or-later
"""Write data/results/forge_calibration_inventory.json: BioForge's calibration inventory (area H, 2026-10-10).

Applies the registration in genomeos/forge/calibration_inventory.py (committed in 5df76f2 before this
census was coded) to every design and experiment block under the CI `bio test` roots. Runs the four
shipped designs once (about a minute); poses no new design question to BioForge; fetches nothing.

    uv run python scripts/forge_calibration_inventory.py [--check]

`--check` writes nothing and exits 1 when the committed file's census (everything but its
`result_manifest`) differs from what the code computes now over the programs the file records; programs
added under the roots since it was written are printed by name, outside that comparison.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.forge import calibration_inventory as ci  # noqa: E402
from genomeos.forge.design_answers import PROGRAM_ROOTS, SEED  # noqa: E402
from genomeos.results import save_result  # noqa: E402

ENTRY = "scripts/forge_calibration_inventory.py"
OWN_CODE = (
    "genomeos/forge/calibration_inventory.py",
    "genomeos/forge/design_answers.py",
    "genomeos/organism/forge.py",
    ENTRY,
    "tests/test_forge_calibration_inventory.py",
)
NAME = "forge_calibration_inventory"
OUT = ROOT / "data" / "results" / f"{NAME}.json"


def result_manifest(reads: list[str], scripts: list[str]) -> dict[str, Any]:
    inputs = [
        mf.files_entry("the organism programs the census parsed, with their imports", reads, partition=None)
    ]
    inputs += [
        {
            **mf.input_entry(ROOT / s, partition=None, read_as="gate registration named by a program"),
            "path": s,
        }
        for s in scripts
    ]
    return {
        "sources": [
            {
                "accession": "the BioLang programs under " + ", ".join(PROGRAM_ROOTS) + " in this repository",
                "version": "the committed files, by sha256 in inputs",
            },
            {
                "accession": "the citations recorded in those programs' evidence strings and in the gate"
                " registrations they name",
                "version": "as committed; no publication was fetched or read",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "n/a: organism programs; no genome sequence or coordinate is read",
        "coordinates": "n/a: no genomic coordinates are read or written",
        "parameters": {
            "seed": SEED,
            "program_roots": list(PROGRAM_ROOTS),
            "confidence_level": ci.CONFIDENCE_LEVEL,
            "max_interval_width": ci.MAX_INTERVAL_WIDTH,
            "minimum_independent_triples": ci.MINIMUM_INDEPENDENT_TRIPLES,
            "registration_commit": "5df76f2",
        },
        "exclusions": [
            "design and experiment blocks outside the CI bio test roots (the compiled programs under"
            " data/knowledge are git-ignored and are not CI programs)",
            "published experiments the design grammar cannot express (additions at a later time,"
            " ablations given as an environment, programs with no organism): counted and set aside",
        ],
        "partitions": "n/a: an inventory, nothing is trained, fitted or held out",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    ap.add_argument("--check", action="store_true", help="compare with the committed file; write nothing")
    args = ap.parse_args()
    mf.trace_begin()
    t0 = time.time()
    if args.check:
        if not OUT.is_file():
            print("forge calibration inventory: no committed file")
            return 1
        written = json.loads(OUT.read_text(encoding="utf-8"))
        inv = json.loads(json.dumps(ci.build_inventory(ROOT, programs=written["programs_scanned"])))
        same = ci.census_of(written) == inv
        print("forge calibration inventory: " + ("matches the committed file" if same else "DIFFERS"))
        added = ci.programs_added_since(written, ROOT)
        if added:
            print(f"programs added since it was written, outside its census: {added}")
        return 0 if same else 1
    inv = json.loads(json.dumps(ci.build_inventory(ROOT)))
    scripts = sorted({s for p in inv["programs_scanned"] for s in ci.named_scripts(ROOT, p)})
    inv["seconds"] = round(time.time() - t0, 1)
    inv["money"] = "none: nothing is fetched, no key is read"
    path = save_result(NAME, inv, manifest=result_manifest(mf.traced_reads(), scripts))
    c = inv["counts"]
    print(
        f"wrote {path}: {c['design_blocks']} design triples,"
        f" {c['experiment_blocks']} experiment blocks"
        f" ({c['expressible_as_inverse_design_questions']} expressible as design questions);"
        f" independent {c['independent_triples']} of a registered minimum {c['minimum_independent_triples']}."
        f" {inv['decision']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
