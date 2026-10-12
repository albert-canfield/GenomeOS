#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Write `data/results/pseudobulk_rule.json`: the inventory rule AND THE MINIMUM, before counting.

THIS RUNS AND COMMITS BEFORE ANYTHING IS COUNTED. It reads the two committed records that fix the
covered set -- `data/results/freshrna_result.json` and `data/results/freshrna_sources.json` -- names
the thirteen covered cell types with the bulk experiment already on disk for each, and states the
MINIMUM the comparison needs before it is meaningful: how many cell types and how many genes
(`genomeos.attribution.pseudobulk.MINIMUM`). Nothing is searched here and no measured value is read.

The ONLY thing learned before this rule was written is ENCODE's own `assay_title` vocabulary over
released transcription experiments -- which titles exist for single-nucleus and single-cell RNA --
recorded in `vocabulary_learned_first` below. No count of qualifying cell types was taken first.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import pseudobulk as pb  # noqa: E402
from genomeos.results import save_result  # noqa: E402

ENTRY = "scripts/pseudobulk_rule.py"
OWN_CODE = (
    "genomeos/attribution/pseudobulk.py",
    "scripts/pseudobulk_rule.py",
    "scripts/pseudobulk_inventory.py",
    "tests/test_pseudobulk.py",
)

VOCABULARY_LEARNED_FIRST = {
    "what": (
        "ENCODE's own assay_title facet over type=Experiment&status=released"
        "&assay_slims=Transcription, read once before this rule was written, to learn WHICH assay "
        "titles exist for single-nucleus and single-cell RNA. A vocabulary, not a count: the query "
        "was made with limit=0 and no biosample breakdown, so no number of qualifying cell types "
        "was read before the minimum below was fixed."
    ),
    "titles_found": {"snRNA-seq": 395, "scRNA-seq": 30, "long read scRNA-seq": 64},
    "titles_taken": list(pb.SINGLE_CELL_ASSAYS),
}


def result_manifest(covered: list[dict[str, Any]], read_now: list[str]) -> dict[str, Any]:
    inputs = [
        mf.input_entry(
            p,
            partition="not an evaluation partition: a committed record that fixes the covered set",
        )
        for p in sorted(set(read_now))
        if "/data/" in p.replace("\\", "/") or p.startswith("data/")
    ]
    return {
        "sources": [
            {
                "accession": "data/results/freshrna_result.json and "
                "data/results/freshrna_sources.json, M2i's committed records (`f6e71c6`)",
                "version": "as committed in this checkout; pinned by sha256 in pinned_sha256",
            },
            {
                "accession": "ENCODE portal assay_title facet over released human transcription "
                "experiments; vocabulary only, no biosample breakdown and no file fetched",
                "version": f"queried {time.strftime('%Y-%m-%d')}; counts in "
                f"vocabulary_learned_first; {pb.ENCODE_LICENCE}",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "n/a: no coordinate is read. This record states a rule and a minimum.",
        "coordinates": "n/a: no interval is read or written",
        "parameters": {
            "selection_rule_version": pb.SELECTION_RULE["version"],
            "match_field": pb.MATCH_FIELD,
            "single_cell_assays": list(pb.SINGLE_CELL_ASSAYS),
            "minimum_cell_types": pb.MINIMUM["cell_types"],
            "minimum_genes": pb.MINIMUM["genes"],
            "covered_cell_count": len(covered),
            "scale_if_used": pb.SCALE,
        },
        "exclusions": [
            "every non-human experiment, and every experiment not status=released",
            "genetically modified and treated biosamples: a knockdown is not the cell type's own RNA",
            "every controlled-access source: GTEx individual-level data, dbGaP, HCA managed access",
            "no biosample is added, no AlphaGenome request is made, no deletion is requested "
            "(2026-09-27 stop list, reinstated)",
            "no experimental data file is fetched by this step",
        ],
        "partitions": "n/a: nothing is fitted and nothing is held out; this record states a rule",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    covered = pb.covered_cells()
    record: dict[str, Any] = {
        "result": "pseudobulk_rule",
        "date": time.strftime("%Y-%m-%d"),
        "lane": "lane-pseudobulk",
        "stage": "the rule AND THE MINIMUM, committed before anything is counted",
        "question": pb.SELECTION_RULE["question"],
        "why_this_lane_exists": (
            "M2h (`ac5580b`) found the compiled rules fit measured RNA better when scaled by one "
            "fitted constant, held-out Delta MSE +0.04339 [98.33% CI +0.0379, +0.04894], scale "
            "+0.2146, about 2% of variance; M2i (`f6e71c6`) replicated it with nothing fitted on "
            "RNA fresh to this project, +0.03727 [+0.03202, +0.04299]. lane-trainexp measured that "
            "independence from AlphaGenome's training data is UNKNOWN FOR 13 OF 13 (`18de375`), and "
            "the only way past that is a cell type OUTSIDE AlphaGenome's training domain. "
            "lane-oodinventory (`c441b90` rule, `171f5bf` result) found exactly ONE such candidate "
            "with open paired data, bile duct UBERON:0002394, and ITS RNA IS SINGLE-NUCLEUS ONLY. "
            "The supervisor ruled a test there NOT INTERPRETABLE AS IT STANDS: it would change the "
            "cell type AND the assay modality at once, so a result would be ambiguous whichever way "
            "it fell. This lane is the prerequisite that would make it interpretable, and it needs "
            "no biosample lift."
        ),
        "symmetry": (
            "If pseudobulk reproduces bulk on this endpoint, the bile-duct test becomes "
            "interpretable and the supervisor's recommendation to Albert can change. If it does "
            "not, the out-of-domain test cannot be done with that candidate and the biosample "
            "question closes without a lift. Both outcomes are useful and this rule prefers "
            "neither."
        ),
        "selection_rule": pb.SELECTION_RULE,
        "minimum": pb.MINIMUM,
        "minimum_stated_before_counting": (
            "This record is written and committed before the search runs. The two numbers in "
            "`minimum` -- 4 cell types and 1,000 genes -- are fixed here, with their reasons, and "
            "the search that follows cannot move them."
        ),
        "vocabulary_learned_first": VOCABULARY_LEARNED_FIRST,
        "covered_cells": covered,
        "covered_cell_count": len(covered),
        "covered_cell_keys": [c["cell"] for c in covered],
        "bulk_side_is_already_on_disk": (
            "All thirteen covered cell types satisfy the bulk clause by construction: M2i read one "
            "open ENCODE bulk RNA-seq experiment for each, and the files are on disk under "
            "data/cache/freshrna/rna/ with their sha256 in data/results/freshrna_sources.json. The "
            "open question is the single-nucleus clause alone."
        ),
        "what_the_search_will_write": "data/results/pseudobulk_inventory.json",
        "scale_if_used": pb.SCALE,
        "scale_origin": (
            "M2h's fitted scale from `ac5580b`, used in this lane ONLY as the fixed constant "
            "+0.2146. Nothing is refitted here, and this lane fits nothing beyond what its "
            "registration declares."
        ),
        "wording_limit": pb.WORDING_LIMIT,
        "no_biosample_added": pb.SELECTION_RULE["no_biosample"],
        "alphagenome_requests": 0,
        "money": "$0.00: no model request, no key, no data needing an application",
        "pinned_sha256": {
            "data/results/freshrna_result.json": pb.sha256_of(pb.FRESHRNA_RESULT),
            "data/results/freshrna_sources.json": pb.sha256_of(pb.FRESHRNA_SOURCES),
            "data/results/magcal_result.json": pb.sha256_of(pb.MAGCAL_RESULT),
            "data/results/oodinventory_candidates.json": pb.sha256_of(pb.OODINVENTORY_RESULT),
        },
    }
    record["seconds"] = round(time.time() - started, 2)
    record["result_manifest"] = result_manifest(covered, mf.traced_reads())
    path = save_result("pseudobulk_rule", record)
    print(json.dumps({"written": str(path), "minimum": pb.MINIMUM["cell_types"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
