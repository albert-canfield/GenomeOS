# SPDX-License-Identifier: AGPL-3.0-or-later
"""Gate 2's metric recounted as registered: a merge judged against each merged block's own starting
label (the registered run compared it with the first block's only).

    uv run python scripts/pilot_metric_recount.py

A correction of a counting defect found after gate 2's registered run (docs/ATTRIBUTION.md, "The
result: the pilot beats the unchanged labels ..."). It re-solves every held-out source exactly as
`scripts/pilot_biological_gate.py` did, first reproduces the run's own counts with the run's counter
(`pilot_bio.MERGE_COUNT_AS_RUN`) and stops if they differ, then counts as registered. Nothing is
scored again: no labelling's score, pass or reading can change. Writes
data/results/pilot_metric_recount.json.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import holdout as ho  # noqa: E402
from genomeos.attribution import pilot as pl  # noqa: E402
from genomeos.attribution import pilot_bio as pb  # noqa: E402
from genomeos.results import RESULTS_DIR, save_result  # noqa: E402

NAME = "pilot_metric_recount"
RUN = "pilot_biological_gate"
KEYS = ("committed", "validated", "errors", "neutral", "untested")


@mf.depends_on_models("alphagenome")  # the pilot's target prior reads the compiled layer
def manifest(parameters: dict[str, Any]) -> dict[str, Any]:
    run = RESULTS_DIR / f"{RUN}.json"
    return {
        "sources": [
            {
                "accession": f"this repository: data/results/{RUN}.json (the registered run)",
                "version": "pinned by sha256",
            },
            {
                "accession": "the inputs of that run, as its manifest lists them",
                "version": "pinned by sha256",
            },
        ],
        "inputs": [mf.input_entry(run, partition=None)],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": parameters,
        "exclusions": ["nothing is scored: only the committed corrections' validation is counted again"],
        "partitions": {
            "validation": "chromosomes " + ", ".join(pb.VALIDATION),
            "held_out_source": "each source in turn, as in the registered run",
        },
    }


def main() -> None:
    t0 = time.process_time()
    run = json.loads((RESULTS_DIR / f"{RUN}.json").read_text())
    fixed = pb.Fixed.load(pb.VALIDATION)
    units = ho.all_units()
    refs = ho.crispri_references()
    cache: dict[str, pl.Outcome] = {}
    per: dict[str, Any] = {}
    totals = {"as_run": dict.fromkeys(KEYS, 0), "as_registered": dict.fromkeys(KEYS, 0)}
    for src in pb.HELD_OUT:
        view = ho.evidence(src, units, refs)
        resp = pb.responsiveness(view)
        held = [u for u in units[src] if u.chrom in pb.VALIDATION]
        solved: dict[str, pb.Solved] = {}
        for chrom in pb.VALIDATION:
            obs, _ = pb.view_obs(view, fixed, chrom)
            rel = pb.relevant_blocks([u for u in held if u.chrom == chrom], fixed.blocks[chrom])
            solved[chrom] = pb.solve("joint", chrom, rel, obs, fixed, resp, cache)
        vep = pb.VALIDATION_ENDPOINT[ho.assay_of(src)]
        pb.MERGE_COUNT_AS_RUN[0] = True
        try:
            as_run = pb.validate(solved, held, vep)
        finally:
            pb.MERGE_COUNT_AS_RUN[0] = False
        registered_run = run["sources"][src]["corrections"]["pilot"]
        if any(as_run[k] != registered_run[k] for k in KEYS):
            raise SystemExit(
                f"{src}: the re-solve does not reproduce the registered run's counts; nothing written"
            )
        fixed_count = pb.validate(solved, held, vep)
        per[src] = {
            "as_run": {k: as_run[k] for k in KEYS},
            "as_registered": {k: fixed_count[k] for k in KEYS},
            "by_kind_as_registered": fixed_count["by_kind"],
            "examples_as_registered": [x for x in fixed_count["examples"] if x["change"]["kind"] == "merge"][
                :10
            ],
        }
        for k in KEYS:
            totals["as_run"][k] += as_run[k]
            totals["as_registered"][k] += fixed_count[k]
        print(src, per[src]["as_run"], per[src]["as_registered"], flush=True)
    cpu_h = run["metric"]["pilot_cpu_hours"]
    payload = {
        "question": "gate 2's validated corrections per compute-hour, with a merge judged against each "
        "merged block's own starting label as registered",
        "status": pb.STATUS,
        "correction_of": RUN,
        "defect": "the registered run's counter compared a merged part with its first block's starting label "
        "only, so the second block's change was never counted",
        "reproduced_the_run": True,
        "per_source": per,
        "totals": totals,
        "metric_as_registered": {
            **totals["as_registered"],
            "pilot_cpu_hours": cpu_h,
            "validated_per_cpu_hour": round(totals["as_registered"]["validated"] / cpu_h, 2),
            "errors_per_cpu_hour": round(totals["as_registered"]["errors"] / cpu_h, 2),
            "cpu_hours_from": "the registered run's measured pilot CPU time; this recount's time is beside",
        },
        "recount_cpu_seconds": round(time.process_time() - t0, 1),
        "alphagenome_requests": 0,
    }
    p = save_result(
        NAME, payload, manifest=manifest({"merge_count": "each block against its own starting label"})
    )
    print("totals", totals, "saved", p, flush=True)


if __name__ == "__main__":
    main()
