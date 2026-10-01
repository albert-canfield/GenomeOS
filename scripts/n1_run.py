# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run N1 once, after the owner authorises measurement access: the calibration gate, then the analysis.

    uv run --frozen python scripts/n1_run.py --normalized PATH --raw PATH --authorisation "who, when"

Not to be run before that authorisation. It applies data/results/n1_registration.json exactly as frozen
(genomeos/attribution/n1_perturb_response.py `run`), checks both files' md5 against Figshare before any
read, records their sha256, and writes data/results/n1_result.json. It refuses when that result already
exists: the experiment is applied once, and a failed gate or an insufficient coverage is its result.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import n1_perturb_response as n1  # noqa: E402
from genomeos.results import RESULTS_DIR, load_result, save_result  # noqa: E402

RESULT = "n1_result"
REGISTRATION = "n1_registration"


def execute(
    normalized: Path,
    raw: Path,
    authorisation: str,
    results_dir: Path = RESULTS_DIR,
    expected: dict | None = None,
) -> dict:
    """The run itself. `expected` replaces Figshare's md5 only for synthetic files in tests."""
    raise SystemExit(
        "refused: superseded by amendment 1 (data/results/n1_registration_amendment_1.json); "
        "run scripts/n1_run_v2.py"
    )
    if (results_dir / f"{RESULT}.json").exists():
        raise SystemExit(f"refused: {results_dir / RESULT}.json exists; N1 is applied once")
    if not authorisation.strip():
        raise SystemExit("refused: --authorisation must name who authorised measurement access, and when")
    reg = load_result(REGISTRATION, results_dir)
    if not reg or "plan" not in reg:
        raise SystemExit(f"refused: no {REGISTRATION} in {results_dir}")
    reader = n1.H5adPseudobulk(normalized, raw, expected=expected)
    out = n1.run(reg["plan"], reader)
    files = n1.SOURCE["files"]
    payload = {
        "authorisation": authorisation,
        "registration_commit": (reg.get("result_manifest") or {}).get("code", {}).get("git_sha"),
        "files_sha256": reader.sha256,
        "datasets_read": reader.datasets_read,
        **out,
        "result_manifest": {
            "sources": [
                {
                    "accession": "Replogle et al. 2022 processed Perturb-seq, Figshare+ 20029387",
                    "version": "v1",
                    "doi": n1.SOURCE["doi"],
                    "license": "CC BY 4.0",
                }
            ],
            "inputs": [
                {
                    "path": str(normalized),
                    "sha256": reader.sha256["normalized"],
                    "md5": files["normalized"]["md5"],
                    "partition": "calibration and evaluation",
                },
                {
                    "path": str(raw),
                    "sha256": reader.sha256["raw"],
                    "md5": files["raw"]["md5"],
                    "partition": "calibration (strata only)",
                },
                mf.input_entry(results_dir / f"{REGISTRATION}.json", partition=None),
            ],
            "assembly": "GRCh38",
            "coordinates": "n/a: genes are matched by Ensembl ID; no interval is read at run time",
            "parameters": dict(n1.CONSTANTS),
            "exclusions": [
                {
                    "what": "genes uncalibrated at the gate",
                    "count": len(out.get("gate", {}).get("genes_uncalibrated", [])),
                }
            ],
            "partitions": {
                "calibration": "the 585 non-targeting rows",
                "evaluation": "the knockdown-eligible candidate rows",
            },
        },
    }
    save_result(RESULT, payload, results_dir)
    return payload


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--normalized", type=Path, required=True, help="K562_gwps_normalized_bulk_01.h5ad")
    ap.add_argument("--raw", type=Path, required=True, help="K562_gwps_raw_bulk_01.h5ad")
    ap.add_argument("--authorisation", required=True, help="who authorised measurement access, and when")
    args = ap.parse_args(argv)
    out = execute(args.normalized, args.raw, args.authorisation)
    print(f"N1: {out['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
