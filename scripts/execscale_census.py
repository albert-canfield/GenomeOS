#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""execscale census: the committed state, registered BEFORE the calibrated magnitude is added.

    uv run python scripts/execscale_census.py            # write data/results/execscale_census.json
    uv run python scripts/execscale_census.py --verify    # recompute and compare to that record

WHY A CENSUS COMES FIRST. The next commit adds the magnitude M2h measured (`ac5580b`, scale +0.2146
on the multiplicative composition's centred prediction) to the runtime as an OPT-IN parameter whose
default is 1.0. The whole safety property of that change is that at the default NOTHING moves, and a
property nobody measured before the change cannot be checked after it. So the committed state is
digested and written down here, with the claim stated in advance and the method of checking it named.

THE CLAIM (`CLAIM` below), held to byte level: with the parameter at its default, every committed
.bio programme is byte-identical and every committed census count is unchanged.

HOW IT IS VERIFIED (`VERIFICATION` below): sha256 over the bytes of each programme and over the group,
and the `bio test` counts recomputed by running the same checks in process. A digest comparison, not
a reading of diffs or a glance at a count: `--verify` recomputes every field and prints the first
field that differs, or PASS.

Input paths go to `input_entry` REPO-RELATIVE, so the record does not pin this machine's checkout.
"""

from __future__ import annotations

import argparse
import contextlib
import hashlib
import io
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

from genomeos import bio
from genomeos import manifest as mf
from genomeos.lang import parse_file
from genomeos.results import RESULTS_DIR, save_result

RESULT = "execscale_census"
OWN_CODE = (
    "scripts/execscale_census.py",
    "scripts/execscale_demo.py",
    "tests/test_execscale.py",
    "genomeos/runtime/grn.py",
    "genomeos/attribution/compile.py",
)
#: the paths `bio test` is run over, in this order; the lane's acceptance runs the same three
TEST_PATHS = ("data/demo", "genomeos/std", "data/organisms")
#: the programmes whose bytes the claim is about: every committed .bio under these trees
PROGRAM_TREES = ("data/organisms", "genomeos/std", "data/demo")

CLAIM = (
    "With the calibrated magnitude at its default of 1.0, NOTHING changes: every committed .bio"
    " programme under data/organisms, genomeos/std and data/demo is byte-identical to the bytes"
    " digested here, and every count in this census -- the `bio test` file count and checks-passed"
    " total, each programme's entity and rule count, and each programme's declared `# test:` lines,"
    " chr21's `# test: rules == 5176` among them -- is unchanged. The parameter is opt-in: a run gets"
    " the measured magnitude only by asking for it by name, and no default run and no committed"
    " programme asks"
)
VERIFICATION = (
    "sha256 over the bytes of each programme, plus one group digest over all of them in a fixed order,"
    " plus the counts recomputed by running `bio test`'s own checks in process. `--verify` recomputes"
    " every field of this record and compares field by field, printing the first mismatch with both"
    " values. Byte digests and integer counts, never a reading of a diff and never an eyeballed count:"
    " the error this guards against is a change that looks inert in a diff and moves a number"
)
DOES_NOT_CLAIM = (
    "this census says nothing about RNA, about whether the calibrated composition predicts expression,"
    " or about whether the default SHOULD be 1.0. It is a record of the committed state and of one"
    " safety property of an addition to it. The RNA evidence for the magnitude is M2h's (`ac5580b`)"
    " and M2i's (`f6e71c6`) and is not re-read here"
)
#: the magnitude the next commit makes available, quoted here so the census names what it guards
#: against: the number is FIXED, read from `ac5580b`, and nothing in this lane fits anything
MAGNITUDE = {
    "value": 0.2146,
    "source_commit": "ac5580b",
    "reading": (
        "M2h: on genes held out of the fit, the multiplicative composition scaled by the one factor"
        " fitted on the training genes came closer to measured RNA than the no-rule baseline (held-out"
        " Delta MSE +0.04339 [98.33% CI +0.0379, +0.04894] over 108,264 gene-cell pairs, 8,328 genes;"
        " control p = 0.000999; fitted scale +0.2146 [98.33% CI +0.2007, +0.2283]), about 2% of the"
        " measured variance"
    ),
    "it_is_a_global_average": (
        "lane-scalehet/hetblock2 (`d9e3ddd`) measured that the magnitude is NOT one global number: it"
        " varies by chromosome and that heterogeneity survives a locus-block null. +0.2146 is a global"
        " AVERAGE over a varying magnitude and must be described that way"
    ),
    "nothing_is_fitted_here": "the constant is read from ac5580b; this lane fits nothing",
}

_TOTAL = re.compile(r"^(\d+) files, (\d+)/(\d+) checks passed$")
_PER_FILE = re.compile(r"^([✓✗]) (\S+): module (\S+), (\d+) entities, (\d+) rules, (\d+)/(\d+) checks$")


def programs() -> list[Path]:
    """Every committed .bio programme the claim covers, in a fixed order."""
    out: list[Path] = []
    for tree in PROGRAM_TREES:
        out += sorted(Path(tree).rglob("*.bio"))
    return out


#: this checkout, found from this script's own location rather than from the working directory
REPO = Path(__file__).resolve().parent.parent


def as_repo_path(path: Path) -> Path:
    """`path` as this repository names it: relative to REPO when it is inside REPO, else unchanged.

    The group digest hashes a programme's NAME beside its bytes, and the name must be the repository's
    so that the digest is the same wherever the checkout sits. The census writer already hands relative
    names; this makes an absolute path under REPO hash the same, which is what the check comparing the
    record to the tree needs. A path outside REPO (a temporary file in a test) is left as it is.
    """
    return Path(os.path.relpath(path, REPO)) if path.is_absolute() and path.is_relative_to(REPO) else path


def group_digest(paths: list[Path]) -> str:
    """One digest over the group: each programme's repository path and bytes, in the order given."""
    h = hashlib.sha256()
    for p in paths:
        p = as_repo_path(p)
        h.update(str(p).encode())
        h.update(b"\0")
        h.update(p.read_bytes())
    return h.hexdigest()


def bio_test_counts() -> dict[str, Any]:
    """`bio test TEST_PATHS` run in process: the totals and the per-file line for every programme."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        code = bio.main(["test", *TEST_PATHS])
    totals: dict[str, Any] = {"exit_code": code}
    per_file: dict[str, dict[str, Any]] = {}
    for line in buf.getvalue().splitlines():
        m = _TOTAL.match(line)
        if m:
            totals["files"] = int(m.group(1))
            totals["checks_passed"] = int(m.group(2))
            totals["checks_total"] = int(m.group(3))
            continue
        m = _PER_FILE.match(line)
        if m:
            per_file[m.group(2)] = {
                "mark": m.group(1),
                "module": m.group(3),
                "entities": int(m.group(4)),
                "rules": int(m.group(5)),
                "checks_ok": int(m.group(6)),
                "checks_total": int(m.group(7)),
            }
    totals["paths"] = list(TEST_PATHS)
    return {"totals": totals, "per_file": per_file}


def declared_tests(path: Path) -> list[str]:
    """The programme's own `# test:` lines, verbatim and in file order."""
    return [line.strip() for line in path.read_text().splitlines() if line.strip().startswith("# test:")]


def census() -> dict[str, Any]:
    """Every field the claim is held to, recomputed from the working tree."""
    paths = programs()
    rows = []
    for p in paths:
        digest, size, _ = mf.sha256_of(p)
        module = parse_file(p)
        rows.append(
            {
                "path": str(p),
                "sha256": digest,
                "bytes": size,
                "lines": len(p.read_text().splitlines()),
                "module": module.name,
                "entities": len(module.entities),
                "rules": len(module.rules),
                "declared_test_lines": declared_tests(p),
            }
        )
    bt = bio_test_counts()
    return {
        "programs": rows,
        "program_count": len(rows),
        "group_sha256": group_digest(paths),
        "compiled_rule_totals": {r["path"]: r["rules"] for r in rows},
        "bio_test": bt,
        "chr21": {str(p): declared_tests(p) for p in paths if "chr21" in p.name},
    }


#: data files a programme itself reads while it parses (the C. elegans contact tables), declared
#: because the counting path opens them even though the claim is about the .bio bytes
READ_WHILE_PARSING = (
    "data/organisms/celegans/contacts_acvu.tsv",
    "data/organisms/celegans/contacts_embryo.tsv",
    "data/organisms/celegans/contacts_embryo_cshaper.tsv",
)


def manifest() -> dict:
    inputs = [mf.input_entry(p, partition=None) for p in programs()]
    inputs.append(
        mf.files_entry("the contact tables the C. elegans programmes read as they parse", READ_WHILE_PARSING)
    )
    return {
        "sources": [
            {"accession": "the committed .bio programmes of this checkout", "version": "as committed"},
            {"accession": "genomeos.bio's own `# test:` checks, run in process", "version": "as committed"},
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "test_paths": list(TEST_PATHS),
            "program_trees": list(PROGRAM_TREES),
            "method": "sha256 per programme and over the group; counts from bio.main(['test', ...])",
            "magnitude_guarded_against": MAGNITUDE["value"],
            "magnitude_source": MAGNITUDE["source_commit"],
        },
        "exclusions": [
            "no model request, no key, no network, no money; nothing is fitted and no RNA is read"
        ],
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
        "partitions": {
            "none": "one population, the committed programmes; no partition is taken and none is implied"
        },
    }


def _flat(d: Any, prefix: str = "") -> dict[str, Any]:
    out: dict[str, Any] = {}
    if isinstance(d, dict):
        for k, v in d.items():
            out.update(_flat(v, f"{prefix}.{k}" if prefix else str(k)))
    elif isinstance(d, list):
        for i, v in enumerate(d):
            out.update(_flat(v, f"{prefix}[{i}]"))
    else:
        out[prefix] = d
    return out


def verify(results_dir: Path = RESULTS_DIR) -> int:
    """Recompute the census and compare it field by field to the committed record."""
    record = json.loads((results_dir / f"{RESULT}.json").read_text())
    want = _flat({k: record[k] for k in CENSUS_KEYS})
    got = _flat(census())
    bad = [k for k in sorted(set(want) | set(got)) if want.get(k, "<absent>") != got.get(k, "<absent>")]
    if bad:
        for k in bad[:20]:
            print(f"MOVED {k}: committed {want.get(k, '<absent>')!r} -> now {got.get(k, '<absent>')!r}")
        print(f"FAIL {len(bad)} field(s) differ from {RESULT}.json")
        return 1
    print(f"PASS {len(got)} fields identical to {RESULT}.json; group sha256 {got['group_sha256']}")
    return 0


#: the keys of the record the claim is held to; the prose keys carry no number and are not compared
CENSUS_KEYS = (
    "programs",
    "program_count",
    "group_sha256",
    "compiled_rule_totals",
    "bio_test",
    "chr21",
)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify", action="store_true", help="compare the working tree to the record")
    args = ap.parse_args(argv)
    if args.verify:
        return verify()
    payload = {
        "question": (
            "what is the committed state that the calibrated magnitude, added as an opt-in parameter"
            " with default 1.0, must leave untouched?"
        ),
        "claim": CLAIM,
        "verification": VERIFICATION,
        "magnitude": MAGNITUDE,
        "does_not_claim": DOES_NOT_CLAIM,
        **census(),
    }
    print("saved", save_result(RESULT, payload, manifest=manifest()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
