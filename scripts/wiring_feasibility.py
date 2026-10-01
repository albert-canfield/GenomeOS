# SPDX-License-Identifier: AGPL-3.0-or-later
"""The wiring diagnostic's feasibility: can matched rewirings be drawn at all?

Writes data/results/wiring_feasibility.json.

    uv run --frozen python scripts/wiring_feasibility.py

Measures the matching only, against genomeos/attribution/wiring.THRESHOLDS, which were committed before
this script was first run. For the K562 CRISPRi training pairs it assesses the three nulls under each
scheme of wiring.SCHEMES, and the primary null on held-out K562; the first scheme in order whose primary
assessment meets every threshold is the one a registration may use, and none is the no-go.

No model is fitted, no AUPRC or other performance figure of any assignment is computed, and no model
request is made. Reads the two benchmark tables, the sweep's deletion table and per-element cache
(as crispri.score_published does), and, for the expression covariate, the non-targeting rows of
Replogle et al.'s K562 raw pseudobulk file already cached for N1 (md5 checked against Figshare first).
"""

from __future__ import annotations

import csv
import gzip
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import crispri, wiring  # noqa: E402
from genomeos.attribution import n1_perturb_response as n1  # noqa: E402
from genomeos.attribution.measured import CRISPRI_SPLIT_OF  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "wiring_feasibility"
RAW_H5AD = Path("data/cache/n1/K562_gwps_raw_bulk_01.h5ad")


def ensembl_ids() -> tuple[dict[str, str], list[str]]:
    """Each benchmark gene symbol's Ensembl ID (version stripped), from both tables; a symbol given two IDs
    is left out and named."""
    seen: dict[str, set[str]] = defaultdict(set)
    for name in (crispri.TRAINING, crispri.HELDOUT):
        with gzip.open(crispri.KNOWLEDGE / name, "rt") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                gid = (r.get("measuredGeneEnsemblId") or "").split(".")[0]
                if gid.startswith("ENSG"):
                    seen[r["measuredGeneSymbol"]].add(gid)
    one = {s: next(iter(ids)) for s, ids in seen.items() if len(ids) == 1}
    return one, sorted(s for s, ids in seen.items() if len(ids) > 1)


@mf.depends_on_models("alphagenome")  # the sweep's model, as far as the disk says
def manifest(training: list[crispri.Pair], heldout: list[crispri.Pair], raw_sha: str) -> dict:
    inputs = [
        mf.input_entry(crispri.KNOWLEDGE / name, partition=CRISPRI_SPLIT_OF[name], pairs=len(pairs))
        for name, pairs in ((crispri.TRAINING, training), (crispri.HELDOUT, heldout))
    ]
    inputs.append(mf.input_entry(crispri.ELEMENTS, partition=None, role="DeletionTable: elements by overlap"))
    inputs.append(
        mf.input_entry(
            crispri.ELEMENT_CACHE, partition=None, role="ElementCache: every gene's deletion value"
        )
    )
    inputs.append(
        {
            "path": str(RAW_H5AD),
            "sha256": raw_sha,
            "md5": n1.SOURCE["files"]["raw"]["md5"],
            "partition": None,
            "role": "K562 expression covariate: mean raw pseudobulk over the non-targeting rows",
        }
    )
    return {
        "sources": [
            {
                "accession": "EngreitzLab/CRISPR_comparison resources/crispr_data, EPCrisprBenchmark "
                "training_K562 and heldout_5_cell_types (Gschwind et al.)",
                "version": "main branch, unpinned upstream; fetched 2026-09-16; pinned here by sha256",
                "url": crispri.BASE_URL,
            },
            {
                "accession": "ENCODE SCREEN cCREs scored by AlphaGenome deletion (all-element table and "
                "per-element response cache)",
                "version": "AlphaGenome as served during the 2026-09 all-element sweep (unpinned); pinned "
                "here by sha256",
            },
            {
                "accession": "Replogle et al. 2022 processed Perturb-seq, Figshare+ 20029387, "
                "K562_gwps_raw_bulk_01.h5ad (non-targeting rows only)",
                "version": "v1",
                "doi": n1.SOURCE["doi"],
                "license": n1.SOURCE["license"],
            },
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "thresholds": wiring.THRESHOLDS,
            "regulated_floor": wiring.REGULATED_FLOOR,
            "schemes": list(wiring.SCHEMES),
            "nulls": list(wiring.NULLS),
            "draws": wiring.FEASIBILITY_DRAWS,
            "seed": wiring.SEED,
            "derange_tries": wiring.DERANGE_TRIES,
            "activity_bins": wiring.ACTIVITY_BINS,
            "expression_bins": wiring.EXPRESSION_BINS,
            "min_distance": crispri.MIN_DISTANCE,
            "model_cells": list(crispri.MODEL_CELLS),
            "alphagenome_requests": 0,
            "performance_metrics_computed": 0,
        },
        "exclusions": [
            "pairs outside the population (uncovered, or covered with no cached value for the gene) are "
            "never rewired; they are counted under pairs_outside_population and keep their real values",
            "held-out pairs in cell types other than K562 are not assessed",
            "a gene symbol the benchmark gives two Ensembl IDs has no expression value (named in "
            "expression_coverage)",
        ],
        "partitions": {
            "training": "the 10,356 K562 training pairs: the primary pair set",
            "heldout_k562": "the 1,918 held-out K562 pairs: a secondary pair set, primary null only",
        },
    }


def main() -> int:
    t0 = time.time()
    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    crispri.annotate(training + heldout, crispri.DeletionTable(), crispri.ElementCache())
    md5, raw_sha = n1.file_digests(RAW_H5AD)
    if md5 != n1.SOURCE["files"]["raw"]["md5"]:
        print(f"refused: {RAW_H5AD} md5 {md5} is not Figshare's {n1.SOURCE['files']['raw']['md5']}")
        return 2
    by_gene_id = wiring.control_expression(RAW_H5AD)
    ensembl_of, conflicts = ensembl_ids()
    expression = wiring.expression_by_symbol(ensembl_of, by_gene_id)
    heldout_k562 = [p for p in heldout if p.cell == "K562"]
    sets = {"training": training, "heldout_k562": heldout_k562}

    coverage = {}
    for name, pairs in sets.items():
        genes = {p.gene for p in pairs if wiring.in_population(p)}
        coverage[name] = {
            "population_genes": len(genes),
            "with_expression": sum(g in expression for g in genes),
        }
    coverage["symbols_with_two_ensembl_ids"] = conflicts
    coverage["measured_genes_in_file"] = len(by_gene_id)

    assessments = []
    for scheme in wiring.SCHEMES:
        for null in wiring.NULLS:
            assessments.append(wiring.assess(training, expression, "training", null, scheme))
        assessments.append(wiring.assess(heldout_k562, expression, "heldout_k562", "element_kept", scheme))
    primary = [a for a in assessments if a["pair_set"] == "training" and a["null"] == "element_kept"]
    chosen = wiring.choose(primary)
    decision = {
        "rule": "the first scheme in SCHEMES order whose primary assessment (element_kept, training) meets "
        "every threshold; none is the no-go",
        "chosen_scheme": chosen,
        "no_go": chosen is None,
        "primary_by_scheme": {a["scheme"]["name"]: a["verdict"]["adequate"] for a in primary},
        "secondaries_at_chosen_scheme": None
        if chosen is None
        else {
            f"{a['pair_set']}/{a['null']}": a["verdict"]["adequate"]
            for a in assessments
            if a["scheme"]["name"] == chosen["name"]
            and not (a["pair_set"] == "training" and a["null"] == "element_kept")
        },
    }
    payload = {
        "status": "feasibility of the matching only; no performance metric of any assignment was computed",
        "lane": "lane-wiring",
        "question": "do the existing element-to-gene deletion assignments outperform matched rewired ones on "
        "the already-examined K562 CRISPRi pairs",
        "nulls": wiring.NULLS,
        "schemes": list(wiring.SCHEMES),
        "thresholds": wiring.THRESHOLDS,
        "regulated_floor": wiring.REGULATED_FLOOR,
        "thresholds_note": "THRESHOLDS, SCHEMES and REGULATED_FLOOR were committed before this script was "
        "first run; see the code stamp below",
        "expression_coverage": coverage,
        "decision": decision,
        "assessments": assessments,
        "alphagenome_requests": 0,
    }
    payload[mf.KEY] = manifest(training, heldout, raw_sha)
    path = save_result(RESULT, payload)
    for a in assessments:
        v = a["verdict"]
        failed = ", ".join(c["check"] for c in v["failed"][:6])
        print(
            f"{a['scheme']['name']} {a['pair_set']:13s} {a['null']:13s} adequate={v['adequate']}"
            + (f"  failed: {failed}{' ...' if len(v['failed']) > 6 else ''}" if failed else "")
        )
    print(f"decision: {decision['chosen_scheme']} no_go={decision['no_go']}")
    print(f"({time.time() - t0:.0f} s) -> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
