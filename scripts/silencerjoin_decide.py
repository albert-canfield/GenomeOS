#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Write the ReSE join decision -- overlap rule, cell-line mapping, what is reported -- BEFORE the join.

Writes data/results/silencerjoin_decision.json, committed alone and before scripts/silencerjoin_run.py
may run. It reads the ReSE records and nothing else: no scored element is opened, so the decision
cannot carry a join outcome.

Free: the ReSE records are the probe's NCBI cache, read in place. No model request, no key, no money.
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import silencerjoin as sj  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "silencerjoin_decision"
ENTRY = "scripts/silencerjoin_decide.py"
OWN_CODE = (
    "genomeos/attribution/silencerjoin.py",
    "scripts/silencerjoin_decide.py",
    "scripts/silencerjoin_run.py",
    "tests/test_silencerjoin.py",
)

#: Keys a join would fill. The decision holding any of them, at any depth, is refused before writing.
FORBIDDEN_IN_DECISION = (
    "elements",
    "by_cell_line",
    "elements_given_the_role",
    "derived_call_on_the_role_elements",
    "sensitivity",
    "pairs_sharing_a_base",
    "population",
)


def forbidden(doc: Any, path: str = "") -> list[str]:
    """Every key, at any depth, that only a join could fill."""
    out: list[str] = []
    if isinstance(doc, dict):
        for k, v in doc.items():
            if k in FORBIDDEN_IN_DECISION:
                out.append(f"{path}/{k}")
            out.extend(forbidden(v, f"{path}/{k}"))
    elif isinstance(doc, list):
        for i, v in enumerate(doc):
            out.extend(forbidden(v, f"{path}[{i}]"))
    return out


def result_manifest() -> dict[str, Any]:
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
        ],
        "inputs": [
            mf.input_entry(sj.RESE_CACHE, partition="the external source; data/cache is never committed"),
        ],
        "assembly": "GRCh38; the ReSE intervals are NCBI genomicinfo and never the GRCh37 coordinates "
        "in each record's description",
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": "a fragment is [chrstart, chrstop + 1) from genomicinfo; elements are [start, end) "
            "as stored",
        },
        "parameters": {
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
        "partitions": "n/a: a decision written before a join; no evaluation partition exists",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    payload = sj.decision()
    bad = forbidden(payload)
    if bad:
        raise SystemExit(f"the decision holds fields only a join could fill: {bad}")
    payload["stop_list"] = (
        "checked against ROADMAP section 5, the stop list of 2026-09-27 reinstated: this adds no "
        "biosample, no area J work and no AlphaGenome request"
    )
    payload["cost"] = "none: 0 model requests, no key, no money, no network (the cache is read in place)"
    payload["result_manifest"] = result_manifest()
    path = save_result(NAME, payload)
    print(f"wrote {path}, sha256 {hashlib.sha256((ROOT / path).read_bytes()).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
