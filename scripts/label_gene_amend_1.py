# SPDX-License-Identifier: AGPL-3.0-or-later
"""Write AMENDMENT_1 of the label-gene registration as its own file, before any count.

    uv run --frozen python scripts/label_gene_amend_1.py

The landed registration at 47835d1 is not rewritten: this repository records an amendment as a
separate file naming what it amends, the way gene_identity recorded its three. The amendment is
additive - the bands, the floors, the four readings and (c)'s words are repeated here unchanged so
a reader can check that none moved - and it is pre-count: when it is committed, no concordance
exists in the tree.

0 model requests, no network, no money.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import label_gene as lg  # noqa: E402

RESULT = "label_gene_registration_amendment_1"
REGISTRATION = Path("data/results/label_gene_registration.json")

OWN_CODE = (
    "genomeos/attribution/label_gene.py",
    "scripts/label_gene_register.py",
    "scripts/label_gene_amend_1.py",
    "scripts/label_gene_count.py",
    "tests/test_label_gene.py",
)


def payload() -> dict[str, Any]:
    entries: list[dict[str, Any]] = []
    if REGISTRATION.exists():
        entries.append(
            mf.input_entry(
                REGISTRATION,
                partition="the registration this file amends, committed at 47835d1 and not "
                "rewritten by this amendment",
            )
        )
    body = lg.amendment_1_payload()
    body["result_manifest"] = {
        "sources": [
            {
                "accession": "this lane's own registration, data/results/label_gene_registration.json",
                "version": "as committed at 47835d1, digested in inputs",
            }
        ],
        "inputs": entries,
        "input_count": len(entries),
        "assembly": "n/a: an amendment to a design; it reads no sequence and no interval",
        "coordinates": "n/a: an amendment to a design; it reads no coordinate",
        "parameters": {
            "top_genes_for_concentration": lg.TOP_GENES_FOR_CONCENTRATION,
            "band_at_control": list(lg.BAND_AT_CONTROL),
            "band_far_above": lg.BAND_FAR_ABOVE,
            "minimum_clusters_genes": lg.MIN_CLUSTERS,
            "nothing_moved": "every band and floor above is identical to 47835d1's",
        },
        "exclusions": [
            "nothing of the registration is withdrawn, weakened or exempted by this amendment",
            "the gene-equal-weight figure is excluded from selecting the reading: it is a "
            "SENSITIVITY and decides nothing",
        ],
        "partitions": {"amendment": "a design document; it carries no observation and no denominator"},
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    return body


def main() -> int:
    from genomeos.results import save_result

    path = save_result(RESULT, payload())
    print(f"amended: {path}")
    print("  (1) the concentration of within-gene pairs in the top 5 genes is PRINTED")
    print("  (2) a gene-equal-weight concordance is a SENSITIVITY and decides nothing")
    print("  (3) if the two disagree, the result must say so in those words")
    print("no concordance computed; 0 model requests, no network, no money")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
