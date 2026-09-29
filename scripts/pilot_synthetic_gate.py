# SPDX-License-Identifier: AGPL-3.0-or-later
"""Item 12 S3 / item 13, the coherence pilot's first gate: a constructed neighbourhood whose true labels
need two simultaneous corrections, where no single correction fits or even helps. The pilot must find
the pair; if it cannot, it stops and nothing biological is scored.

    uv run python scripts/pilot_synthetic_gate.py

Registered 2026-09-29 before it was run (genomeos/attribution/pilot.py SYNTH_*; docs/ATTRIBUTION.md
"Item 13 pilot registered"). No data is read and no model is called. Writes
data/results/pilot_synthetic_gate.json.
"""

from __future__ import annotations

import resource
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import pilot as pl  # noqa: E402
from genomeos.results import save_result  # noqa: E402

NAME = "pilot_synthetic_gate"
BUDGET_CPU_SECONDS = 600  # registered: past ten CPU minutes the gate is reported as over budget


def manifest(parameters: dict[str, Any]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "this repository: genomeos/attribution/pilot.py, the registered generator "
                "of synthetic neighbourhoods (no external data)",
                "version": "pinned by sha256",
            }
        ],
        "inputs": [mf.input_entry(ROOT / "genomeos/attribution/pilot.py", partition=None)],
        "assembly": "n/a: synthetic neighbourhoods on a synthetic chromosome, no genome",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": parameters,
        "exclusions": [
            "generated draws that fail the registered validity check are counted by reason and not scored"
        ],
        "partitions": "n/a: a synthetic solver gate, no evidence is held out",
    }


def main() -> None:
    t0, c0 = time.time(), time.process_time()
    res = pl.run_synthetic_gate()
    cpu = time.process_time() - c0
    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)
    params = {
        "registration": pl.registration(),
        "seed": pl.SYNTH_SEED,
        "instances_per_motif": pl.SYNTH_INSTANCES,
        "max_draws_per_motif": pl.SYNTH_MAX_DRAWS,
        "pass_share": pl.SYNTH_PASS_SHARE,
        "motifs": list(pl.SYNTH_MOTIFS),
        "family_motif": pl.SYNTH_FAMILY_MOTIF,
        "compute_budget_cpu_seconds": BUDGET_CPU_SECONDS,
    }
    payload = {
        "question": "does the pilot find a repair that needs two simultaneous corrections, where no single "
        "correction fits the evidence or lowers the energy, and does it return indistinguishable "
        "alternatives as one family with a measurement that separates them",
        "status": "a solver gate on constructed neighbourhoods: it tests the search, not biology",
        "pass_rule": pl.SYNTH_PASS_RULE,
        "validity": pl.SYNTH_VALIDITY,
        "controls": pl.SYNTH_CONTROLS,
        "falsifier": pl.SYNTH_FALSIFIER,
        **res,
        "within_budget": cpu <= BUDGET_CPU_SECONDS,
        "compute": {
            "cpu_seconds": round(cpu, 2),
            "wall_seconds": round(time.time() - t0, 2),
            "peak_rss_mb": round(peak_mb, 1),
            "measured_with": "time.process_time and getrusage, this process",
        },
        "alphagenome_requests": 0,
    }
    p = save_result(NAME, payload, manifest=manifest(params))
    summary = {
        m: (v["generated_passed"], v["generated"], v["hand_built"]["passed"])
        for m, v in res["motifs"].items()
    }
    print(res["verdict"], summary)
    print("saved", p, round(cpu, 1), "cpu s", flush=True)


if __name__ == "__main__":
    main()
