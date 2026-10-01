# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run N1 under amendment 1, once, after the owner authorises measurement access.

    uv run --frozen python scripts/n1_run_v2.py --normalized PATH --raw PATH --authorisation "who, when"

Before any measurement is opened it checks the freeze: the sha256 of the original registration, of
genomeos/attribution/n1_perturb_response.py and of this runner, and the constants, against
data/results/n1_registration_amendment_1.json; any difference refuses. Then it checks both files' md5
against Figshare and records their sha256, and applies `run_v2`: identities, the calibration gate on the
non-targeting rows, and, while the knockdown rule is unsupported (amendment 1, item 2), a stop before any
candidate row is read. It writes data/results/n1_result_amendment_1.json and refuses when that exists.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import n1_perturb_response as n1  # noqa: E402
from genomeos.results import RESULTS_DIR, load_result, save_result  # noqa: E402

RESULT = "n1_result_amendment_1"
REGISTRATION = "n1_registration"
RUNNER = Path(__file__).resolve()
MODULE = Path(n1.__file__).resolve()


def execute(
    normalized: Path,
    raw: Path,
    authorisation: str,
    results_dir: Path = RESULTS_DIR,
    expected: dict | None = None,
) -> dict:
    """The run itself. `expected` replaces Figshare's md5 only for synthetic files in tests."""
    if (results_dir / f"{RESULT}.json").exists():
        raise SystemExit(f"refused: {results_dir / RESULT}.json exists; N1 is applied once")
    if not authorisation.strip():
        raise SystemExit("refused: --authorisation must name who authorised measurement access, and when")
    amendment = load_result(n1.AMENDMENT, results_dir)
    original = load_result(REGISTRATION, results_dir)
    if not amendment or not original or "plan" not in original:
        raise SystemExit(f"refused: {n1.AMENDMENT} and {REGISTRATION} must both be in {results_dir}")
    problems = n1.freeze_problems(amendment, results_dir / f"{REGISTRATION}.json", MODULE, RUNNER)
    if problems:
        raise SystemExit("refused, the code or constants differ from the freeze: " + "; ".join(problems))
    reader = n1.H5adPseudobulkV2(normalized, raw, expected=expected)
    out = n1.run_v2(original["plan"], reader)
    files = n1.SOURCE["files"]
    payload = {
        "authorisation": authorisation,
        "amendment_commit": (amendment.get("result_manifest") or {}).get("code", {}).get("git_sha"),
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
                    "partition": "calibration",
                },
                {
                    "path": str(raw),
                    "sha256": reader.sha256["raw"],
                    "md5": files["raw"]["md5"],
                    "partition": "calibration (strata only)",
                },
                mf.input_entry(results_dir / f"{REGISTRATION}.json", partition=None),
                mf.input_entry(results_dir / f"{n1.AMENDMENT}.json", partition=None),
            ],
            "assembly": "GRCh38",
            "coordinates": "n/a: genes are matched by Ensembl ID; no interval is read at run time",
            "parameters": dict(n1.CONSTANTS_V2),
            "exclusions": [
                {
                    "what": "genes uncalibrated at the gate",
                    "count": len(out.get("gate", {}).get("genes_uncalibrated", [])),
                }
            ],
            "partitions": {
                "calibration": "the 585 non-targeting rows; no candidate row is read under amendment 1"
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
    print(f"N1 (amendment 1): {out['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
