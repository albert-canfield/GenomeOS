# SPDX-License-Identifier: AGPL-3.0-or-later
"""Record, after the result, the three facts about this lane that git shows and the files do not.

    uv run --frozen python scripts/label_gene_disclose.py

A DISCLOSURE, not an amendment: it is written after the figures are known, so it is not a
pre-registration and nothing in it is presented as one. It changes no band, no floor, no reading,
no control and no label rule, and the result's reading stays (c), the data cannot tell.

What it records: which of this lane's three files were written while a peer's fix to
genomeos/results.py sat uncommitted, with each file's sha and time and why the two that carry no
figure are disclosed rather than rewritten; that this lane's label universe is enumerated from its
own declared input and resolves no CURIE and reads no ontology file; and that the bands were
committed fifteen minutes before any value of control 1 existed, with an explicit list of what had
and had not been seen when the registration was written.

0 model requests, no network, no money.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import label_gene as lg  # noqa: E402

RESULT = "label_gene_disclosure"
ABOUT = (
    Path("data/results/label_gene_registration.json"),
    Path("data/results/label_gene_registration_amendment_1.json"),
    Path("data/results/label_gene.json"),
)

OWN_CODE = (
    "genomeos/attribution/label_gene.py",
    "scripts/label_gene_register.py",
    "scripts/label_gene_amend_1.py",
    "scripts/label_gene_count.py",
    "scripts/label_gene_disclose.py",
    "tests/test_label_gene.py",
)


def payload() -> dict[str, Any]:
    entries = [
        mf.input_entry(path, partition="a file of this lane that this disclosure is about")
        for path in ABOUT
        if path.exists()
    ]
    body = lg.disclosure_payload()
    body["result_manifest"] = {
        "sources": [
            {
                "accession": "this lane's own registration, amendment and result",
                "version": "as committed at 47835d1, f579926 and 6ab0b1f, digested in inputs",
            }
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "n/a: a disclosure about this lane's own files; it reads no sequence",
        "coordinates": "n/a: a disclosure about this lane's own files; it reads no coordinate",
        "parameters": {
            "band_at_control": list(lg.BAND_AT_CONTROL),
            "band_far_above": lg.BAND_FAR_ABOVE,
            "minimum_clusters_genes": lg.MIN_CLUSTERS,
            "nothing_moved": "every band and floor above is identical to 47835d1's, and the "
            "reading of the result stays (c)",
        },
        "exclusions": [
            "no design change: this file amends nothing and may not be read as a registration",
            "no recount: the result at 6ab0b1f was already computed in a worktree of a tree "
            "containing d33146a",
        ],
        "partitions": {"disclosure": "carries no observation and no denominator"},
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    return body


def main() -> int:
    from genomeos.results import save_result

    path = save_result(RESULT, payload())
    print(f"disclosed: {path}")
    print("  (1) the result at 6ab0b1f ran the results.py its stamped tree holds: no recount")
    print("  (2) the registration's and amendment's stamps are stale or pre-fix and are NOT")
    print("      rewritten, because repairing them after the figures are known is the worse fault")
    print("  (3) the label universe is enumerated from the declared input: no CURIE, no ontology")
    print("  (4) the bands were committed before any value of control 1 existed")
    print("no design changed; the reading stays (c). 0 model requests, no network, no money")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
