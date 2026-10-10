# SPDX-License-Identifier: AGPL-3.0-or-later
"""AstroREG-2's registered post-run analyses. REFUSES unless the run's ledger shows the run complete.

    uv run python scripts/astroreg2_analyse.py --run 2

The gate comes first and reads the LEDGER ONLY: until its last line is the `run_complete` the sender
writes, no answer is opened and the refusal says why. Run 2 today stopped at request 364 of 1,232 with
no such line, so today this script refuses.

Once the run is complete, and after the supervisor's read over every answer, it computes item (h) as
`astroargmax_registration.json` and its amendment 1 fix it, and reports the primary claim's assembly as
REFUSED with each unregistered term named (genomeos/attribution/astropost.py). Nothing is sent, no key is
read, and no ledger is written.
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

import astroreg2_send  # noqa: E402

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import astropost, astrorun  # noqa: E402

RESULT = "astroreg2_analysis"
ENTRY = "scripts/astroreg2_analyse.py"
OWN_CODE = (
    "genomeos/attribution/astropost.py",
    "scripts/astroreg2_analyse.py",
    "tests/test_astropost.py",
)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--run",
        type=int,
        required=True,
        help="which run's ledger and answers to analyse; the ledger is derived from it, never chosen",
    )
    args = ap.parse_args()
    t0 = time.time()
    ledger = astrorun.ledger_path_for_run(args.run)
    try:
        gate = astropost.check_the_run_is_complete(ledger)
    except astropost.AnalysisRefusedError as e:
        print(f"REFUSED: {e}")
        return 2
    print(f"the run is complete on its ledger: {gate['charges']} charged, {gate['answered']} answered")

    plan = astroreg2_send.reviewed_plan()
    cache_root = astrorun.ROOT_FOR_BLOBS / astrorun.ASTROREG2_CACHE_ROOT
    try:
        answers = astropost.load_the_answers(plan, ledger, cache_root)
    except astropost.AnalysisRefusedError as e:
        print(f"REFUSED: {e}")
        return 2
    chrom_of = {r["element"]: r["chrom"] for r in plan}
    out: dict[str, Any] = {
        "status": "AstroREG-2's registered post-run analyses, on a complete run",
        "run_id": args.run,
        "gate": gate,
        "item_h": astropost.item_h(answers, chrom_of),
        "primary_claim": astropost.primary_assembly(),
    }
    out[mf.KEY] = manifest(ledger, plan)
    from genomeos.results import save_result

    print(f"item (h), {astropost.DECIDING_SET}: {out['item_h']['the_reading']['words']}")
    print(f"-> {save_result(RESULT, out)} ({time.time() - t0:.0f} s)")
    return 0


def manifest(ledger: Path, plan: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "sources": [
            {
                "accession": "AlphaGenome answers of the AstroREG-2 run, as written by its sender",
                "version": "the request cache root, read only for elements the run's ledger answered",
            }
        ],
        "inputs": [
            mf.input_entry(ledger, partition=None, role="the run's ledger: complete, and what it answered"),
            mf.input_entry(astroreg2_send.PLAN, partition=None, role="the reviewed request list"),
            mf.input_entry(
                astropost.ASTROARGMAX_REGISTRATION, partition=None, role="item (h)'s registration"
            ),
            mf.input_entry(astropost.ASTROARGMAX_AMENDMENT_1, partition=None, role="its amendment 1"),
            mf.input_entry(
                astropost.ASTROREG2_REGISTRATION, partition=None, role="AstroREG-2's registered claim"
            ),
        ],
        "assembly": "GRCh38",
        "coordinates": {"base": 0, "interval": "half-open"},
        "parameters": {"requests_in_the_reviewed_list": len(plan)},
        "cache_root": str(astrorun.ASTROREG2_CACHE_ROOT),
        "exclusions": [
            "no answer is read unless the ledger's last line is run_complete and not a pilot stop",
            "the primary claim's assembly is refused with its unregistered terms named",
        ],
        "partitions": {"astrocyte": "the answered elements of the reviewed list"},
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


if __name__ == "__main__":
    raise SystemExit(main())
