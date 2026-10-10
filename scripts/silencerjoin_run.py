#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Run the ReSE join under the committed decision, and write data/results/silencerjoin.json.

Refuses to run unless data/results/silencerjoin_decision.json is committed AND its rule and cell-line
mapping are, field for field, the ones the module would apply now: the rule was fixed before the join
and may not move under it.

Reads the derived direction call and never writes it: no scored element file is opened for writing.
Free: 0 model requests, no key, no money, no network.
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
from genomeos.attribution import silencerjoin as sj  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "silencerjoin"
ENTRY = "scripts/silencerjoin_run.py"
DECISION = sj.DECISION
OWN_CODE = (
    "genomeos/attribution/silencerjoin.py",
    "scripts/silencerjoin_decide.py",
    "scripts/silencerjoin_run.py",
    "tests/test_silencerjoin.py",
)


def decision_commit() -> str:
    out = subprocess.run(  # noqa: S603
        ["git", "log", "-1", "--format=%H", "--", DECISION],  # noqa: S607
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    sha = out.stdout.strip()
    if not sha:
        raise SystemExit(
            f"{DECISION} is not committed. The join rule must be fixed before the join runs; run "
            "scripts/silencerjoin_decide.py and commit it first."
        )
    return sha


def decision_still_applies() -> dict[str, Any]:
    """The committed decision, refused if its rule or mapping differs from what would be applied."""
    committed = json.loads((ROOT / DECISION).read_text())
    if committed["rule"] != sj.rule():
        raise SystemExit("the committed decision's rule differs from the module's: refusing to join")
    if committed["cell_line_mapping"]["rese_to_scored_line"] != sj.CELL_LINE_MAP:
        raise SystemExit("the committed decision's cell-line mapping differs from the module's")
    return committed


def result_manifest(decision_sha: str) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "ReSE screen-validated silencers in RefSeq Functional Elements, NCBI Gene "
                "esummary, query 'ReSE[All Fields] AND human[orgn]'",
                "version": "fetched 2026-10-03 by scripts/silenceragree_fetch.py; 5,478 records; public "
                "domain: NCBI places no restrictions on the use or distribution of the data it holds",
            },
            {
                "accession": "Pang B, Snyder MP. Systematic identification of silencers in human cells. "
                "Nat Genet 2020;52(3):254-263",
                "version": "10.1038/s41588-020-0578-5",
            },
            {
                "accession": "the project's scored elements and their derived direction call, "
                "data/knowledge/alphagenome/all_elements (the AlphaGenome deletion sweep)",
                "version": "as on disk; pinned by sha256 under inputs",
            },
        ],
        "inputs": [
            mf.input_entry(sj.RESE_CACHE, partition="the external source; data/cache is never committed"),
            mf.input_entry(
                sj.ELEMENTS_DIR,
                partition="the project's scored elements, read for coordinates and, after the join, for "
                "the derived direction; data/knowledge is never committed",
            ),
            mf.input_entry(ROOT / DECISION, partition="the decision this join is bound by"),
        ],
        "assembly": "GRCh38; the ReSE intervals are NCBI genomicinfo and never the GRCh37 coordinates "
        "in each record's description",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "a fragment is [chrstart, chrstop + 1) from genomicinfo; elements are [start, end) "
            "as stored; shared bases are min(ends) - max(starts)",
        },
        "parameters": {
            "decision_commit": decision_sha,
            "role_rule": sj.ROLE_RULE,
            "role_overlap": sj.ROLE_OVERLAP,
            "reciprocal_sensitivity": list(sj.RECIPROCAL_SENSITIVITY),
            "min_side_sensitivity": sj.MIN_SIDE_SENSITIVITY,
            "cell_line_map": sj.CELL_LINE_MAP,
        },
        "exclusions": [
            "480 records naming ReSE only outside the description (genomicinfo is a parent span)",
            "4 with no genomicinfo, 2 unplaced on GRCh38",
            "SilencerDB entirely: no licence is stated anywhere in it",
        ],
        "partitions": "n/a: an annotation join against a curated assay, not an evaluation; nothing is "
        "scored against a held-out split",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    sha = decision_commit()
    decision_still_applies()
    payload = sj.build()
    payload["decision"] = {"file": DECISION, "commit": sha}
    payload["readings"] = {
        "role": {"molecular_role": sj.ROLE, "evidence_status": sj.EVIDENCE_STATUS, "evidence": sj.EVIDENCE},
        "unknown_reasons": sj.REASONS,
    }
    payload["what_the_role_licenses"] = sj.WHAT_THE_ROLE_LICENSES
    payload["what_the_role_does_not_license"] = sj.WHAT_THE_ROLE_DOES_NOT_LICENSE
    payload["derived_call_unchanged"] = (
        "the derived direction call is read from the scored element files and never written; the "
        "committed predicted results are untouched, and this result stands beside them"
    )
    payload["cost"] = "none: 0 model requests, no key, no money, no network"
    payload["result_manifest"] = result_manifest(sha)
    path = save_result(NAME, payload)
    print(f"wrote {path}, sha256 {hashlib.sha256((ROOT / path).read_bytes()).hexdigest()}")
    show = {
        k: payload[k]
        for k in ("population", "by_cell_line", "derived_call_on_the_role_elements", "sensitivity")
    }
    print(json.dumps(show, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
