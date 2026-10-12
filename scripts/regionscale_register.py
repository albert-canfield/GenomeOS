# SPDX-License-Identifier: AGPL-3.0-or-later
"""lane-regionscale, step 2: the registration, written and committed BEFORE any held-out value is read.

    uv run python scripts/regionscale_register.py

It writes `data/results/regionscale_registration.json`: every BOUND field of
`genomeos/attribution/regionscale.py` verbatim, the sha256 of every input the run will read, and the
declared readers. It refuses to run unless the committed harness exists and reports
HARNESS_REPRODUCES, because a registration written on a harness that did not reproduce M2h's numbers
would be a registration of an unknown quantity.

It fits nothing, scores nothing and reads no held-out value. The one thing in it that comes from the
table is a set of GENE COUNTS, which depend on no measured value
(`regionscale.GENE_COUNTS_ARE_NOT_MEASUREMENTS`) and which the harness has already computed and
committed; they are repeated here because the HELD-OUT DESIGN's feasibility is part of what is being
registered.

No model request, no network, no money, 0 AlphaGenome requests.
"""

from __future__ import annotations

import json
import os
import sys
import time
from datetime import date
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import regionscale as rs  # noqa: E402
from genomeos.results import save_result  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from regionscale_harness import (  # noqa: E402
    OWN_CODE,
    READERS,
    SOURCE_FILES,
    load_rows,
)

RESULT = "regionscale_registration"
OUT = "data/results/regionscale_registration.json"
HARNESS_RESULT = "data/results/regionscale_harness.json"


def refuse_unless_harness_reproduced() -> dict[str, Any]:
    if Path(OUT).exists():
        raise SystemExit(f"refused: {OUT} already exists; a registration is written once (STOP_RULE)")
    if not Path(HARNESS_RESULT).exists():
        raise SystemExit(f"refused: {HARNESS_RESULT} does not exist; run and commit the harness first")
    h = json.loads(Path(HARNESS_RESULT).read_text())
    if h["outcome"] != "HARNESS_REPRODUCES":
        raise SystemExit("refused: the committed harness did not reproduce M2h's numbers")
    return h


def manifest(entries: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {"accession": "data/results/multrna_result.json (7h, lane-multrna)", "version": "pinned"},
            {"accession": "data/results/magcal_result.json (M2h, lane-magcal)", "version": "pinned"},
            {
                "accession": "data/results/scalehet_result.json (lane-scalehet)",
                "version": "pinned; its rejection of a single global scale is why this lane exists",
            },
            {
                "accession": "data/results/hetblock_result.json (lane-hetblock2)",
                "version": "pinned; its locus-block re-test of that rejection, carried word for word",
            },
            {
                "accession": "data/results/regionscale_harness.json (this lane, step 1)",
                "version": "pinned; the registration refuses to exist without HARNESS_REPRODUCES",
            },
            {
                "accession": "ENCODE polyA plus RNA-seq gene quantifications, V29, 13 reader biosamples",
                "version": "read by 7d, 7h and M2h; this lane reads their committed rows, not the files",
            },
        ],
        "inputs": entries,
        "input_count": len(entries),
        "readers": READERS,
        "reader_count": len(READERS),
        "assembly": "GRCh38; genes as 7h matched them, by gene_id",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {
            "form": rs.FORM,
            "regional_unit": "chromosome",
            "min_genes_each_half": rs.MIN_GENES,
            "min_regions": rs.MIN_REGIONS,
            "form_interval": list(rs.FORM_INTERVAL),
            "nominal_interval": list(rs.NOMINAL_INTERVAL),
            "bootstrap": {"resamples": rs.BOOTSTRAP["resamples"], "seed": rs.BOOTSTRAP["seed"]},
            "permutations": rs.PERMUTATIONS,
            "control_alpha": rs.CONTROL_ALPHA,
            "tolerance": rs.TOLERANCE,
            "m2h_tolerance": rs.M2H_TOLERANCE,
        },
        "exclusions": [
            "none: the population is 7h's, unchanged. A chromosome below the ADMISSION floors is NOT"
            " excluded -- its held-out genes stay in the scored set and are scored by the global scale"
            " under both arms"
        ],
        "partitions": {
            "fit": "per admitted chromosome, THAT chromosome's training genes alone; globally, all"
            " training genes",
            "scored": "M2h's own held-out genes, all of them, identical for both arms",
            "nothing_is_fitted_here": "this script fits nothing and reads no held-out value",
        },
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }


def main() -> int:
    t0 = time.time()
    harness = refuse_unless_harness_reproduced()
    pinned_files = (*SOURCE_FILES, HARNESS_RESULT)
    entries = [mf.input_entry(p, partition="pinned by sha256") for p in pinned_files]
    rows = load_rows()
    counts = rs.split_counts(rows)
    adm, fall = rs.admitted_regions(counts)
    reg: dict[str, Any] = {
        "registration": RESULT,
        "lane": rs.LANE,
        "step": "2 of 4: the registration, before any held-out value of this design is read",
        "date": date.today().isoformat(),
        "registered": rs.REGISTERED,
        "harness_record": {
            "result": HARNESS_RESULT,
            "outcome": harness["outcome"],
            "m2h_reproduces": harness["m2h_reproduction_FIRST"]["reproduces"],
            "m2h_reading": harness["m2h_reproduction_FIRST"]["reading"],
        },
        # NOTE the ORDER and the NAME. `registered_fields()` carries a BOUND field called "harness"
        # -- the registered HARNESS text -- so the provenance block above is called "harness_record".
        # Written as "harness" it was silently OVERWRITTEN by this spread and the registration lost the
        # proof that M2h reproduced; a test reading this file caught it.
        **rs.registered_fields(),
        "the_held_out_design_is_feasible": {
            "admitted": adm,
            "admitted_count": len(adm),
            "fallback_to_global_both_arms": fall,
            "min_genes_each_half": rs.MIN_GENES,
            "min_regions": rs.MIN_REGIONS,
            "testable": rs.is_testable(adm),
            "these_are_counts": rs.GENE_COUNTS_ARE_NOT_MEASUREMENTS,
        },
        "pinned_sha256": {e["path"]: e["sha256"] for e in entries},
        "readers": READERS,
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "what_is_not_in_this_file": (
            "no per-chromosome scale, no Delta MSE of either arm, no DELTA, no bootstrap interval, no"
            " permutation and no p value. Those are step 3, scripts/regionscale_run.py, which refuses"
            " to start unless this file exists, binds every field above word for word, and finds every"
            " pinned sha256 unmoved"
        ),
        "seconds": round(time.time() - t0, 1),
    }
    print("saved", save_result(RESULT, reg, manifest=manifest(entries)))
    print(
        f"registered: unit=chromosome, admitted={len(adm)} of {len(counts)},"
        f" fallback={fall}, testable={rs.is_testable(adm)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
