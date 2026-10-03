#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Measure the registered check. Run ONLY after data/results/silenceragree_registration.json is
committed.

The script refuses to run if the registration is not committed, because the whole value of this
lane is that the rule was fixed before the number existed.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import silenceragree as sa  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "silenceragree"
ENTRY = "scripts/silenceragree_run.py"
REGISTRATION = "data/results/silenceragree_registration.json"
OWN_CODE = (
    "genomeos/attribution/silenceragree.py",
    "scripts/silenceragree_fetch.py",
    "scripts/silenceragree_register.py",
    "scripts/silenceragree_run.py",
    "tests/test_silenceragree.py",
)


def registration_commit() -> str:
    out = subprocess.run(  # noqa: S603
        ["git", "log", "-1", "--format=%H", "--", REGISTRATION],  # noqa: S607
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    sha = out.stdout.strip()
    if not sha:
        raise SystemExit(
            f"{REGISTRATION} is not committed. The rule must be fixed before the number exists; "
            "run scripts/silenceragree_register.py and commit it ALONE first."
        )
    return sha


#: The contract's shape for `coordinates`, which `save_result` requires as {base, interval} and
#: which `result_manifest()` first stated as prose. ADDED rather than edited in place: the guard
#: refused removing the prose line -- it is work already in HEAD -- and this project supersedes
#: additively. The prose sentence is kept verbatim below as the note, so nothing it said is lost and
#: the declaration a reader of the written file sees is this one.
COORDINATES = {
    "base": 0,
    "interval": "half-open",
    "note": "0-based half-open for the project's elements as stored; the NCBI genomicinfo span is "
    "used as given, which costs at most one base at each end of a fragment whose median length is "
    "192 bp, and overlap is tested as element_start < rese_end and element_end > rese_start. The "
    "overlap in bases is recorded on every joined row, so a stricter rule can be applied afterwards "
    "without re-joining.",
}


#: The 23 committed cCRE subsets the control matches its STRATA on. Declared as a group because the
#: first write was quarantined for reading all 23 without declaring one of them: the file that
#: decides every stratum boundary was invisible to a rebuild.
def ccre_group():
    return mf.files_entry(
        "ENCODE cCRE v3 per-chromosome subsets, read for the element CLASS the control matches on",
        sorted((ROOT / "data/results").glob("ccres_chr*.bed.gz")),
        partition="the registry's own committed cCRE subsets",
    )


def result_manifest(reg_sha: str) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "ReSE screen-validated silencers in RefSeq Functional Elements, "
                "NCBI Gene esummary, query 'ReSE[All Fields] AND human[orgn]'",
                "version": "fetched 2026-10-03; 5,478 records; NCBI places no use restrictions",
            },
            {
                "accession": "Pang B, Snyder MP. Systematic identification of silencers in human "
                "cells. Nat Genet 2020;52(3):254-263",
                "version": "10.1038/s41588-020-0578-5",
            },
            {
                "accession": "ENCODE cCRE v3 element classes, as committed in data/results/ccres_chr*.bed.gz",
                "version": "as committed",
            },
        ],
        "inputs": [
            mf.input_entry(sa.RESE_CACHE, partition="the external source; data/cache is never committed"),
            mf.input_entry(
                sa.ELEMENTS_DIR, partition="the project's scored elements; data/knowledge is never committed"
            ),
            mf.input_entry(
                sa.CACHE_DIR,
                partition="the AlphaGenome element cache, read for the per-gene drop/rise pair that "
                "makes the tie bias countable; data/knowledge is never committed",
            ),
            mf.input_entry(ROOT / REGISTRATION, partition="the registration this run is bound by"),
        ],
        "assembly": "GRCh38; the ReSE intervals are joined on NCBI genomicinfo and never on the "
        "GRCh37 coordinates in each record's description",
        "coordinates": "overlap tested as element_start < rese_end and element_end > rese_start",
        "parameters": {
            "registration_commit": reg_sha,
            "control_draws": sa.CONTROL["draws"],
            "control_seed": sa.CONTROL["seed"],
            "band_threshold_points": sa.BAND_THRESHOLD_POINTS,
            "near_zero": sa.NEAR_ZERO,
        },
        "exclusions": [
            "480 records naming ReSE only outside the description (genomicinfo is a parent span)",
            "4 with no genomicinfo, 2 unplaced on GRCh38",
            "3,408 fragments overlapping no scored element",
            "403 fragments whose overlapping element has no named predicted target",
            "SilencerDB entirely: no licence is stated anywhere in it",
        ],
        "partitions": "none: an external check against a source the model never saw",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    reg_sha = registration_commit()
    payload = sa.readings()
    payload["registration"] = {"file": REGISTRATION, "commit": reg_sha}
    payload["base_rate"] = sa.base_rate()
    payload["bands_as_registered"] = sa.BANDS
    payload["what_this_licenses"] = sa.WHAT_AGREEMENT_LICENSES
    payload["what_this_does_not_license"] = sa.WHAT_AGREEMENT_DOES_NOT_LICENSE
    payload["tie_bias"] = sa.TIE_BIAS
    payload["four_cell_line_negative"] = sa.FOUR_LINE_NEGATIVE
    payload["result_manifest"] = result_manifest(reg_sha)
    # Both lines SUPERSEDE what result_manifest() stated, additively: the guard refuses removing a
    # line that is already in HEAD, so the prose `coordinates` and the short input list stay where
    # they are and the correct values are written over them here. What lands in the file is this.
    payload["result_manifest"]["coordinates"] = COORDINATES
    payload["result_manifest"]["inputs"].append(ccre_group())
    # ADDED, not an edit: the contract wants `partitions` as a dict or a string opening "n/a:", and
    # result_manifest() states it as prose that the validator refused. That line is already in HEAD
    # and the guard refuses removing it, so the correct value is written over it here. The substance
    # is unchanged -- there is no evaluation partition, because the model never saw this source as a
    # label at all; what is NOT claimed is that it never saw the source's REGIONS, which it did.
    payload["result_manifest"]["partitions"] = (
        "n/a: no evaluation partition. This is an external check against a curated assay, not a "
        "held-out split of anything the model was fitted on. Whether the ReSE screen's data entered "
        "the training corpus is NOT ESTABLISHED and is registered as unknown, so this field records "
        "the absence of a partition and must not be read as a claim of independence."
    )
    path = save_result(NAME, payload)
    print(f"wrote {path}, sha256 {hashlib.sha256((ROOT / path).read_bytes()).hexdigest()}")
    print(json.dumps({k: v for k, v in payload.items() if k.startswith("reading")}, indent=1))
    print(json.dumps(payload["falsifier_control_disagrees_with_itself"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
