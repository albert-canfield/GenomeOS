#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""Apply N1b's registration once, after it is committed, and write the result.

The registration is checked before a byte of `X` is asked for: it must exist, it must be committed,
and its sha256 is recorded in the result. If it is uncommitted or absent the run stops and reads
nothing. The frozen code's shas are checked the same way.
"""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import n1b_calibration as n1b  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from scripts import n1b_fetch, n1b_register  # noqa: E402

NAME = "n1b_calibration"
OUT = ROOT / f"data/results/{NAME}.json"
REG = n1b_register.OUT
ENTRY = "scripts/n1b_run.py"
OWN_CODE = n1b_register.OWN_CODE


def git(*args: str) -> str:
    return subprocess.run(  # noqa: S603
        ["git", *args],  # noqa: S607
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    ).stdout.strip()


def freeze_problems() -> list[str]:
    problems = []
    if not REG.exists():
        return ["the registration does not exist; nothing is read"]
    rel = str(REG.relative_to(ROOT))
    if git("diff", "HEAD", "--name-only", "--", rel):
        problems.append(f"{rel} is uncommitted; the registration must be committed before the run")
    if not git("log", "-1", "--format=%H", "--", rel):
        problems.append(f"{rel} has no commit")
    reg = json.loads(REG.read_text())
    for path, rec in reg["frozen_code"].items():
        if n1b_register.sha256_text(ROOT / path) != rec["sha256"]:
            problems.append(f"{path} has changed since the registration recorded its sha256")
    if reg["frozen_code_uncommitted"]:
        problems.append(
            f"the registration was written over uncommitted code: {reg['frozen_code_uncommitted']}"
        )
    return problems


def live_rss_bytes() -> int:
    import resource

    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(usage if sys.platform == "darwin" else usage * 1024)


def _inconclusive(intervals: dict) -> list[dict]:
    """Scopes whose interval SPANS a band edge: the point estimate sits on one side, the interval on
    both, so the data cannot tell whether the band holds. Reported apart from the verdict, because a
    point estimate inside a band with an interval crossing it is not the band holding."""
    bands = {"variance": (n1b.VAR_LO, n1b.VAR_HI), "tail": (n1b.TAIL_LO, n1b.TAIL_HI)}
    out = []
    for name, iv in intervals.items():
        for leg, (lo, hi) in bands.items():
            low, high = iv[leg]["low"], iv[leg]["high"]
            if low is None or high is None:
                continue
            inside = low >= lo and high <= hi
            outside = high < lo or low > hi
            if not inside and not outside:
                out.append(
                    {
                        "scope": name,
                        "leg": leg,
                        "interval": [low, high],
                        "band": [lo, hi],
                        "reading": n1b.INCONCLUSIVE_WORDING,
                        "identical_resample_share": iv[leg]["identical_resample_share"],
                    }
                )
    return out


def main() -> int:
    import numpy as np

    mf.trace_begin()
    problems = freeze_problems()
    if problems:
        for p in problems:
            print(f"REFUSED: {p}")
        return 1
    reg = json.loads(REG.read_text())
    probe = json.loads(n1b_fetch.PROBE.read_text())
    layout = probe["x_layout"]
    nt = probe["non_targeting_rows"]
    if nt != reg["access_plan"]["rows_to_read"]:
        print("REFUSED: the rows on disk are not the rows the registration named")
        return 1

    t0 = time.time()
    x, fetch_log = n1b_fetch.fetch_rows(nt, layout)
    fetch_s = time.time() - t0
    counts = [probe["obs"]["num_cells_filtered"][p] for p in range(len(nt))]
    usable = n1b.usable_rows(counts)
    keep = usable["rows_used"]
    n_used = [counts[p] for p in keep]
    expr_all = n1b_fetch.raw_expression([nt[p] for p in keep], n_used)
    t = n1b.t_from_rows(x[keep], n_used)
    del x
    bad = n1b.uncalibratable_genes(t, expr_all)
    good = [j for j in range(t.shape[1]) if j not in set(bad)]
    t = np.ascontiguousarray(t[:, good])
    expr = [float(expr_all[j]) for j in good]
    strata = n1b.expression_quintiles(expr)
    fig = n1b.gene_figures(t)
    scopes: dict[str, list[int]] = {"overall": list(range(len(good)))}
    stratum_cols = {s: [j for j, st in enumerate(strata) if st == s] for s in range(n1b.STRATA)}
    for s, cols in stratum_cols.items():
        scopes[f"quintile {s}"] = cols
    overall = n1b.summarise(fig, scopes["overall"], t)
    recs = []
    for s, cols in stratum_cols.items():
        rec = n1b.summarise(fig, cols, t)
        ex = [expr[j] for j in cols]
        recs.append(
            {
                "stratum": s,
                "control_expression_min": min(ex) if ex else None,
                "control_expression_max": max(ex) if ex else None,
                "control_expression_median": float(np.median(ex)) if ex else None,
                **rec,
            }
        )
    points = {"overall": overall, **{f"quintile {r['stratum']}": r for r in recs}}
    intervals = n1b.cluster_bootstrap(t, scopes, points)
    verdict = n1b.decide(overall, recs)
    inconclusive = _inconclusive(intervals)
    out = {
        "study": n1b.STUDY,
        "registration": {
            "path": str(REG.relative_to(ROOT)),
            "sha256": n1b_register.sha256_text(REG),
            "commit": git("log", "-1", "--format=%H", "--", str(REG.relative_to(ROOT))),
            "its_own_reading_of_a_fail": reg["what_a_fail_means"]["verdict"],
        },
        "not_a_resumption_of_n1": reg["not_a_resumption_of_n1"],
        "outcome_exposure": reg["outcome_exposure"],
        "rows": {
            "non_targeting_read": len(nt),
            "usable": len(keep),
            "unusable": len(usable["rows_dropped"]),
            "unusable_reasons": usable["rows_dropped"][:3],
            "cells_pooled_over_usable_rows": int(sum(n_used)),
            "cells_per_row_min_median_max": [
                min(n_used),
                float(np.median(n_used)),
                max(n_used),
            ],
        },
        "genes": {
            "in_file": layout["cols"],
            "uncalibratable": len(bad),
            "calibrated": len(good),
            "why_a_gene_leaves": "any non-finite T, or a non-finite control expression",
        },
        "stratifier": reg["stratification"],
        "assumptions": list(n1b.ASSUMPTIONS),
        "bands_are_relative": reg["bands_are_relative"],
        "overall": overall,
        "strata": recs,
        "intervals": intervals,
        "verdict": verdict,
        "inconclusive_scopes": inconclusive,
        "inconclusive_wording": reg["inconclusive_reading"],
        "power_stated_in_advance": reg["power"],
        "fetch": {
            "requests": len(fetch_log),
            "bytes": sum(e["bytes"] for e in fetch_log),
            "seconds": round(fetch_s, 1),
            "log": fetch_log,
            "headers_read_before_the_registration": {
                "requests": len(probe["fetch_log"]),
                "bytes": probe["bytes_fetched"],
                "log": probe["fetch_log"],
            },
            "md5_whole_file": n1b.MD5_NOT_VERIFIABLE,
        },
        "local_data_declared_not_committed": {
            "rule": "Albert has ruled 'data/cache stays local'; nothing below is committed",
            "data/cache/n1/K562_gwps_raw_bulk_01.h5ad": {
                "sha256": n1b_fetch.sha256_file(n1b_fetch.RAW_LOCAL),
                "md5_verified_against_figshare": True,
                "role": "the stratifier only",
            },
            "data/cache/n1b/probe.json": {"sha256": n1b_fetch.sha256_file(n1b_fetch.PROBE)},
        },
        "cost": reg["cost"],
        "memory": {
            "live_peak_rss_bytes": live_rss_bytes(),
            "caveat": "this live reading is valid only in a process that does nothing substantial "
            "first; it measures THIS process and is reported, not asserted on. The tests inject their "
            "ceiling instead of reading a process-wide figure.",
            "matrix_bytes": int(t.nbytes),
        },
    }
    out["result_manifest"] = {
        "sources": [
            {
                "accession": "Replogle et al. 2022 processed Perturb-seq, Figshare+ 20029387 "
                "(K562_gwps_normalized_bulk_01.h5ad: the 585 non-targeting rows of X by exact HTTP "
                "range request, plus the headers and three obs columns the probe read)",
                "version": "10.25452/figshare.plus.20029387.v1, CC BY 4.0",
            },
            {
                "accession": "Replogle et al. 2022 processed Perturb-seq, Figshare+ 20029387 "
                "(K562_gwps_raw_bulk_01.h5ad: the usable non-targeting rows, for the per-gene control "
                "expression that forms the strata and for nothing else)",
                "version": "10.25452/figshare.plus.20029387.v1, CC BY 4.0",
            },
            {
                "accession": "GenomeOS n1b_calibration_registration, committed before any byte of X",
                "version": out["registration"]["commit"] or "uncommitted",
            },
        ],
        "inputs": [
            mf.input_entry(REG, partition="the registration this run applies"),
            mf.input_entry(n1b_fetch.PROBE, partition="the file's own index and obs columns"),
            mf.input_entry(n1b_fetch.RAW_LOCAL, partition="the stratifier"),
        ],
        "assembly": "n/a: pseudobulk rows and gene columns; no coordinate is read",
        "coordinates": "n/a: no genomic interval is read",
        "parameters": {
            **n1b.CONSTANTS,
            "primary_var_summary": n1b.PRIMARY_VAR_SUMMARY,
            "primary_tail_summary": n1b.PRIMARY_TAIL_SUMMARY,
            "pass_rule": n1b.PASS_RULE,
            "x_rows_read": len(nt),
            "x_bytes_read": sum(e["bytes"] for e in fetch_log),
            "x_requests": len(fetch_log),
        },
        "exclusions": [
            "every row that is not non-targeting. NO FACTOR ROW WAS READ, AT ALL",
            f"{len(usable['rows_dropped'])} non-targeting rows carry no finite num_cells_filtered, so "
            "T cannot be formed on them; their bytes were fetched so the exclusion is auditable",
            f"{len(bad)} genes left before any figure was formed: a non-finite T or a non-finite "
            "control expression",
            "obs/control_expr, NaN on all 585 non-targeting rows, so quintiles of it cannot be formed",
        ],
        "partitions": {
            f"quintile {r['stratum']}": f"{r['genes']} genes, control expression "
            f"{r['control_expression_min']:.4g} to {r['control_expression_max']:.4g}"
            for r in recs
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }
    path = save_result(NAME, out)
    print(f"wrote {path}")
    print(f"PASSED: {verdict['passed']}")
    print(f"overall variance {overall['var_primary']:.4f}  tail {overall['tail_primary']:.5f}")
    for r in recs:
        print(
            f"  quintile {r['stratum']}: var {r['var_primary']:.4f}  tail {r['tail_primary']:.5f}  "
            f"genes {r['genes']}  expr median {r['control_expression_median']:.4g}"
        )
    print(f"failed scopes: {verdict['failed_scopes']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
