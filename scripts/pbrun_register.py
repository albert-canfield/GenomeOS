#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registration, written and committed BEFORE any comparison value is read.

The 1,000-gene clause cleared first (`3ef7e82`, 8,051 genes in all four cell types on both sides),
so this lane may run -- and this script is what fixes, in advance, every choice a reader could
otherwise suspect of having been made after a value was seen: the pseudobulk aggregation stated
exactly, the four cell types, the gene population, the endpoint and the statistic carried verbatim
from 7d, 7h, M2h and M2i where they apply, the one deviation named AS a deviation, what agreement
means, what disagreement means, the equivalence margin and where its number comes from, the
falsifier, and what "inconclusive" is so that an ambiguous result cannot be read later as a weak
agreement.

It reads committed records and no matrix. No measured value of either arm is opened here.
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
from genomeos.attribution import pbrun as pr  # noqa: E402
from genomeos.attribution import pseudobulk as pb  # noqa: E402
from genomeos.results import save_result  # noqa: E402

ENTRY = "scripts/pbrun_register.py"
OWN_CODE = (
    "genomeos/attribution/pbrun.py",
    "scripts/pbrun_count.py",
    "scripts/pbrun_register.py",
    "tests/test_pbrun.py",
)


def result_manifest(read_now: list[str], parameters: dict[str, Any]) -> dict[str, Any]:
    inputs = [
        mf.input_entry(p, partition="M2h's HELD-OUT side, as the count record pins it")
        for p in sorted(set(read_now))
        if "/data/" in p.replace("\\", "/") or p.startswith("data/")
    ]
    return {
        "sources": [
            {
                "accession": (
                    "data/results/pbrun_gene_count.json (this lane's committed count), "
                    "data/results/magcal_result.json (M2h), "
                    "data/results/multrna_result.json (7h's genes table)"
                ),
                "version": "as committed in this checkout; pinned by sha256 in pinned_sha256",
            }
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "n/a: a registration. No coordinate and no value is read.",
        "coordinates": "n/a: no interval is read or written",
        "parameters": parameters,
        "exclusions": [
            "every measured value: this is the registration and it reads none",
            "every cell type outside lane-scsweep's four",
            "every gene outside the 8,051 the committed count record pins",
        ],
        "partitions": {
            "read": (
                "M2h's HELD-OUT side of the multiplicative form, intersected to the 8,051 genes that "
                "resolve in all four cell types on both sides"
            ),
            "not_read": "M2h's TRAINING genes: the scale is the constant already fitted on them",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE),
        "entry": ENTRY,
        "own_code": list(OWN_CODE),
    }


def main() -> int:
    started = time.time()
    count = json.loads(pr.GENE_COUNT_PATH.read_text())
    if count["answer"] != "QUALIFIES":
        raise SystemExit(
            f"refused: the committed count answers {count['answer']!r}, so there is nothing to "
            "register and the lane ends there"
        )
    cells = count["qualifying_cells"]
    population = count["qualifying_genes"]
    margin = pr.margin()

    record: dict[str, Any] = {
        "result": "pbrun_registration",
        "date": time.strftime("%Y-%m-%d"),
        "lane": "lane-pbrun",
        "roadmap_item": (
            "area I: does pseudobulked single-cell or single-nucleus RNA reproduce bulk RNA-seq on "
            "M2h's endpoint, in cell types the reader already covers that have both?"
        ),
        "stage": (
            "step 2 of 3: registered BEFORE any comparison value is read, after the 1,000-gene "
            "clause cleared at 8,051 of 8,328"
        ),
        "question": (
            "Read M2h's endpoint twice on the same genes and the same four cell types -- once with "
            "the measured side taken from bulk RNA-seq and once with it taken from pseudobulked "
            "single-cell or single-nucleus RNA -- at M2h's FIXED scale +0.2146. Do the two readings "
            "agree?"
        ),
        "what_it_gates": (
            "AGREEMENT licenses the modality swap FOR THIS ENDPOINT ONLY and makes the bile-duct "
            "out-of-domain test (lane-oodinventory, `171f5bf`) interpretable, so the supervisor's "
            "recommendation to Albert can change. DISAGREEMENT closes Albert's biosample question "
            "without a lift. INCONCLUSIVE licenses nothing."
        ),
        "cells": cells,
        "cells_source": pr.CELLS_SOURCE,
        "cell_count": len(cells),
        "modality_per_cell": {e["cell"]: e["modality"] for e in count["per_cell"]},
        "population": (
            "exactly the 8,051 gene symbols of `qualifying_genes` in "
            "data/results/pbrun_gene_count.json, pinned by its sha256: M2h's HELD-OUT side of the "
            "multiplicative form (8,328 genes of 7h's committed population under M2h's split on the "
            "gene symbol) intersected with the genes that resolve, by lane-pseudobulk's clause word "
            "for word, in every one of the four cell types on BOTH sides. Each gene contributes its "
            "FOUR (gene, cell type) pairs. No filter on expression and no filter on any measured "
            "value; no gene is admitted that 7h did not admit."
        ),
        "population_genes": len(population),
        "population_pairs": len(population) * len(cells),
        "gene_clause_verbatim": pr.GENE_CLAUSE,
        "minimum": pr.MINIMUM,
        "sc_selection": pr.SC_SELECTION,
        "aggregation": pr.AGGREGATION,
        "carried_verbatim": pr.CARRIED,
        "deviation_from_m2h": pr.DEVIATION_FROM_M2H,
        "scale": pr.SCALE,
        "scale_origin": (
            "M2h's fitted constant +0.2145911381525326 from `ac5580b`, per_form.multiplicative.scale "
            "of data/results/magcal_result.json, used as a FIXED CONSTANT. Nothing is fitted in this "
            "lane: a test parses this lane's module and fails on any fitter in it or on any "
            "right-hand side for SCALE other than lane-pseudobulk's own constant."
        ),
        "nothing_is_fitted": (
            "no scale, no intercept, no per-gene or per-cell parameter, no normalisation chosen from "
            "a value, no threshold. The two arms differ in ONE thing, the measured side."
        ),
        "metric": pr.CARRIED["metric"],
        "endpoint": (
            "per arm, Delta MSE over the population's pairs at SCALE. The REGISTERED READING is: the "
            "single-cell arm's Delta MSE with its gene bootstrap interval, the bulk arm's beside it, "
            "and the PAIRED difference (single-cell minus bulk) with its own interval, read against "
            "the margin."
        ),
        "bootstrap": pr.BOOTSTRAP,
        "form_interval": [0.8333, 99.1667],
        "form_interval_words": "98.33% CI",
        "control": pr.CONTROL,
        "permutations": pr.PERMUTATIONS,
        "control_alpha": 0.05,
        "min_genes": pr.CARRIED["min_genes"],
        "margin": margin,
        "margin_rule": pr.MARGIN_RULE,
        "which_governs": pr.CARRIED["which_governs_m2h"],
        "decision": pr.DECISION,
        "agreement_means": pr.DECISION["agreement"],
        "disagreement_means": pr.DECISION["disagreement"],
        "inconclusive_means": pr.DECISION["inconclusive"],
        "falsifier": pr.DECISION["falsifier"],
        "stop_rule": pr.DECISION["stop_rule"],
        "harness": pr.HARNESS,
        "tolerance": pr.TOLERANCE,
        "prior_exposure": pr.EXPOSURE,
        "wording_limit": pr.WORDING_LIMIT,
        "limitations": pr.LIMITATIONS,
        "descriptives_not_verdicts": [
            "each arm's Delta MSE at scale 1, the uncalibrated size 7d and 7h read",
            "each arm's Delta MSE restricted to the pairs of ONE cell type, so the single-NUCLEUS "
            "ovary arm can be read on its own",
            "the Pearson correlation between the two arms' D over all pairs, and per cell type",
            "each cell type's nuclei or cells counted, and its pseudobulk library size",
        ],
        "what_the_run_will_write": "data/results/pbrun_result.json",
        "the_code_that_runs_it": (
            "scripts/pbrun_run.py, which must be committed before it may write into the registry; "
            "it reads this registration and refuses on any disagreement with it (HARNESS)"
        ),
        "no_biosample_added": pb.SELECTION_RULE["no_biosample"],
        "stop_list_check": (
            "docs/ROADMAP.md section 5, 'The stop list of 2026-09-27, reinstated': NO biosample "
            "added, no AlphaGenome request, no deletion request, nothing in area J, no node-model "
            "work, no BioForge calibration."
        ),
        "alphagenome_requests": 0,
        "money": "$0.00: no model request, no key, no data needing an application",
        "pinned_sha256": {
            "data/results/pbrun_gene_count.json": pr.sha256_of(pr.GENE_COUNT_PATH),
            "data/results/magcal_result.json": pr.sha256_of(pr.MAGCAL_RESULT),
            "data/results/multrna_result.json": pr.sha256_of(pr.MULTRNA_RESULT),
            "data/results/freshrna_sources.json": pr.sha256_of(pr.FRESHRNA_SOURCES),
        },
    }
    record["seconds"] = round(time.time() - started, 2)
    record["result_manifest"] = result_manifest(
        mf.traced_reads(),
        {
            "scale": pr.SCALE,
            "margin": margin,
            "resamples": pr.BOOTSTRAP["resamples"],
            "bootstrap_seed": pr.BOOTSTRAP["seed"],
            "permutations": pr.PERMUTATIONS,
            "population_genes": len(population),
            "cell_count": len(cells),
            "sc_selection_version": pr.SC_SELECTION["version"],
        },
    )
    path = save_result("pbrun_registration", record)
    print(
        json.dumps(
            {
                "written": str(path),
                "cells": cells,
                "population_genes": len(population),
                "margin": margin,
                "scale": pr.SCALE,
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
