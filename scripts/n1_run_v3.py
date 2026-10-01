# SPDX-License-Identifier: AGPL-3.0-or-later
"""Run N1 under amendment 2, once, after the owner authorises measurement access.

    uv run --frozen python scripts/n1_run_v3.py --ad PATH --raw PATH --authorisation "who, when"

--ad is the producers' "anderson-darling p-values, BH-corrected.csv.gz" (Figshare+ 21632564) and --raw is
K562_gwps_raw_bulk_01.h5ad (Figshare+ 20029387). Before any data file is opened it checks the freeze: the
sha256 of the original registration, of amendments 1 and 2, of genomeos/attribution/n1_perturb_response.py
and of this runner, and the constants, against data/results/n1_registration_amendment_2.json; any
difference refuses. Then it checks both files' md5 against Figshare, records their sha256, and applies
`run_v3`. It writes data/results/n1_result_amendment_2.json and refuses when that exists.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import n1_perturb_response as n1  # noqa: E402
from genomeos.results import RESULTS_DIR, load_result, save_result  # noqa: E402

RESULT = "n1_result_amendment_2"
REGISTRATION = "n1_registration"
RUNNER = Path(__file__).resolve()
MODULE = Path(n1.__file__).resolve()


def execute(
    ad: Path,
    raw: Path,
    authorisation: str,
    results_dir: Path = RESULTS_DIR,
    expected: dict | None = None,
) -> dict:
    """The run itself. `expected` replaces Figshare's md5 values only for synthetic files in tests."""
    if (results_dir / f"{RESULT}.json").exists():
        raise SystemExit(f"refused: {results_dir / RESULT}.json exists; N1 is applied once")
    if not authorisation.strip():
        raise SystemExit("refused: --authorisation must name who authorised measurement access, and when")
    names = (REGISTRATION, n1.AMENDMENT, n1.AMENDMENT_2)
    regs = {name: load_result(name, results_dir) for name in names}
    if not all(regs.values()) or "plan" not in regs[REGISTRATION]:
        raise SystemExit(f"refused: {', '.join(names)} must all be in {results_dir}")
    amendment = regs[n1.AMENDMENT_2]
    files = {
        "original_registration": results_dir / f"{REGISTRATION}.json",
        "amendment_1": results_dir / f"{n1.AMENDMENT}.json",
        "module": MODULE,
        "runner": RUNNER,
    }
    problems = n1.freeze_problems_v3(amendment, files)
    if problems:
        raise SystemExit("refused, the code or constants differ from the freeze: " + "; ".join(problems))
    expected = expected or {}
    ad_file = n1.PublishedAdFile(ad, expected_md5=expected.get("ad", n1.AD_SOURCE["md5"]))
    raw_file = n1.H5adRawControls(raw, expected=expected.get("raw"))
    out = n1.run_v3(regs[REGISTRATION]["plan"], amendment["ad_header_digest"], raw_file, ad_file)
    payload = {
        "authorisation": authorisation,
        "amendment_commit": (amendment.get("result_manifest") or {}).get("code", {}).get("git_sha"),
        "files_sha256": {"ad": ad_file.sha256, "raw": raw_file.sha256["raw"]},
        "datasets_read": raw_file.datasets_read + [f"ad: {p}" for p in ad_file.passes],
        **out,
        "result_manifest": {
            "sources": [
                {
                    "accession": "Replogle et al. 2022 supplemental files, Figshare+ 21632564",
                    "version": "v1",
                    "license": "CC0",
                },
                {
                    "accession": "Replogle et al. 2022 processed Perturb-seq, Figshare+ 20029387",
                    "version": "v1",
                    "license": "CC BY 4.0",
                },
            ],
            "inputs": [
                {
                    "path": str(ad),
                    "sha256": ad_file.sha256,
                    "md5": n1.AD_SOURCE["md5"],
                    "partition": "evaluation",
                },
                {
                    "path": str(raw),
                    "sha256": raw_file.sha256["raw"],
                    "md5": n1.SOURCE["files"]["raw"]["md5"],
                    "partition": "strata",
                },
                *[mf.input_entry(results_dir / f"{name}.json", partition=None) for name in names],
            ],
            "assembly": "GRCh38",
            "coordinates": "n/a: genes are matched by Ensembl ID; no interval is read at run time",
            "parameters": dict(n1.CONSTANTS_V3),
            "exclusions": [{"what": k, "count": v} for k, v in out.get("coverage", {}).items()],
            "partitions": {
                "strata": "the 585 non-targeting rows of the raw pseudobulk",
                "evaluation": "the 56 principal columns",
            },
        },
    }
    save_result(RESULT, payload, results_dir)
    return payload


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--ad", type=Path, required=True, help="anderson-darling p-values, BH-corrected.csv.gz")
    ap.add_argument("--raw", type=Path, required=True, help="K562_gwps_raw_bulk_01.h5ad")
    ap.add_argument("--authorisation", required=True, help="who authorised measurement access, and when")
    args = ap.parse_args(argv)
    out = execute(args.ad, args.raw, args.authorisation)
    print(f"N1 (amendment 2): {out['status']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
