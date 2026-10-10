# SPDX-License-Identifier: AGPL-3.0-or-later
"""Write the covered-CURIE set and the committed selection rule to
`data/results/oodinventory_rule.json`, BEFORE any public source is searched.

The point of running this first is that the rule cannot then be bent to whatever the search
turns up. The file carries the whole sorted covered set and its sha256, so the candidate list
written later names the set it was judged against.

    uv run python scripts/oodinventory_curies.py

Reads one local CSV: the copy of AlphaGenome's published track metadata. No network, no model
request, no biosample, nothing fitted.
"""

from __future__ import annotations

import collections
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import oodinventory as ood  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "oodinventory_rule"
ENTRY = "scripts/oodinventory_curies.py"
OWN_CODE = (
    "genomeos/attribution/oodinventory.py",
    "scripts/oodinventory_curies.py",
    "tests/test_oodinventory.py",
)


def result_manifest(read_now: list[str]) -> dict[str, Any]:
    """The contract's block. One input, read through the module, recorded by the audit hook."""
    declared = [p for p in read_now if p.endswith("alphagenome_track_metadata_copy.csv")]
    if not declared:  # the hook can miss a read it did not see; declare the path we opened
        declared = [str(ood.METADATA_PATH)]
    inputs = [
        mf.input_entry(p, partition="not an evaluation partition: the catalogue being counted")
        for p in declared
    ]
    return {
        "sources": [
            {
                "accession": "AlphaGenome published output-track metadata, local copy at "
                "data/cache/entex/alphagenome_track_metadata_copy.csv",
                "version": "as held in this checkout; sha256 in metadata_sha256 of this result. "
                "It is a copy of a published track list, NOT a copy of the training set.",
            }
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "n/a: no coordinate is read. This result counts ontology terms in a metadata "
        "table and states a selection rule.",
        "coordinates": "n/a: no interval is read or written",
        "parameters": {
            "rna_output_value": ood.RNA_OUTPUT,
            "match_field": ood.MATCH_FIELD,
            "selection_rule_version": ood.SELECTION_RULE["version"],
        },
        "exclusions": [
            "every row whose output is not rna_seq (chip_histone, cage, dnase, atac, "
            "contact_maps, procap): the question is about the RNA training domain",
            "nonzero_mean and every other numeric column: nothing here reads a signal value",
            "no biosample is added and no out-of-domain test is run (2026-09-27 stop list)",
        ],
        "partitions": "n/a: nothing is fitted and nothing is held out; this is a count and a rule",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    rows = ood.rna_rows(ROOT / ood.METADATA_PATH)
    curies = ood.covered_curies(ROOT / ood.METADATA_PATH)
    record: dict[str, Any] = {
        "lane": "lane-oodinventory",
        "stage": "selection rule, committed before any public source is searched",
        "metadata_file": str(ood.METADATA_PATH),
        "metadata_sha256": mf.sha256_of(ROOT / ood.METADATA_PATH)[0],
        "rna_row_count": len(rows),
        "rna_row_strands": dict(collections.Counter(r["strand"] for r in rows)),
        "rna_distinct_biosample_names": len({r["biosample_name"] for r in rows}),
        "rna_rows_with_empty_curie": sum(1 for r in rows if not r["ontology_curie"].strip()),
        "rna_rows_with_accession": sum(1 for r in rows if r["name"].startswith("ENC")),
        "accessions_anywhere_in_file": sorted(
            {r["name"] for r in _all_rows() if r["name"].startswith("ENC")}
        ),
        "covered_curie_count": len(curies),
        "covered_curie_digest": ood.curie_digest(curies),
        "covered_curie_prefixes": dict(collections.Counter(c.split(":")[0] for c in sorted(curies))),
        "covered_curies": sorted(curies),
        "selection_rule": ood.SELECTION_RULE,
        "wording_limit": ood.SELECTION_RULE["wording_limit"],
        "alphagenome_requests": 0,
        "money": "none: one local metadata file is read",
        "seconds": round(time.time() - started, 2),
    }
    record["result_manifest"] = result_manifest(mf.traced_reads())
    path = save_result(NAME, record)
    print(
        f"wrote {path}: {len(rows)} rna_seq rows, {len(curies)} distinct {ood.MATCH_FIELD}, "
        f"{record['rna_rows_with_accession']} with an accession"
    )
    return 0


def _all_rows() -> list[dict[str, str]]:
    """Every row, not only the RNA ones, so the result can name the file's only accessions."""
    import csv

    with (ROOT / ood.METADATA_PATH).open(newline="") as handle:
        return list(csv.DictReader(handle))


if __name__ == "__main__":
    raise SystemExit(main())
