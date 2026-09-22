# SPDX-License-Identifier: AGPL-3.0-or-later
"""Calibrate the predicted enhancer-to-gene targets against the CRISPRi screens, then band the sweep.

    uv run python scripts/target_calibration.py              # the calibration and the whole sweep
    uv run python scripts/target_calibration.py --training    # the fit only, no held-out pair scored
    uv run python scripts/target_calibration.py --chroms chr21,chr22

Reads the cached benchmark tables (data/knowledge/crispri), the all-enhancer deletion table
(data/knowledge/alphagenome/all_elements), the registry classes (data/results/ccres_*.bed.gz) and
GENCODE v50 (data/reference), and saves data/results/target_calibration.json. No model request, no
network, and the committed sweep summaries are not touched.
"""

from __future__ import annotations

import argparse
import time
from typing import Any

from genomeos.attribution import crispri
from genomeos.attribution import target_calibration as tc
from genomeos.results import save_result


def show_reliability(rows: list[dict[str, Any]]) -> None:
    print(
        f"    {'bin':>4} {'pairs':>6} {'cov reg':>8} {'cov not':>8} {'predicted':>10} "
        f"{'observed':>9} {'95% interval':>17}"
    )
    for i, r in enumerate(rows, 1):
        creg, cnot = r.get("coverage_regulated", [0, 0, None]), r.get("coverage_not_regulated", [0, 0, None])
        lo, hi = r["ci95"]
        print(
            f"    {i:>4} {r['pairs']:>6} {str(creg[2]):>8} {str(cnot[2]):>8} "
            f"{r['mean_predicted']:>10.4f} {r['observed']:>9.4f} "
            f"{f'{lo:.4f}-{hi:.4f}':>17} {'' if r['consistent'] else 'OUTSIDE'}"
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--training", action="store_true", help="fit and report on training pairs only")
    ap.add_argument("--chroms", default="", help="comma-separated chromosomes for the sweep pass")
    args = ap.parse_args()
    t0 = time.time()

    training, heldout = crispri.load(crispri.TRAINING), crispri.load(crispri.HELDOUT)
    table = crispri.DeletionTable()
    if args.training:
        crispri.annotate(training, table)
        found = tc.add_features(training, table)
        train = tc.scored(training)
        labels = [p.regulated for p in train]
        top = [p for p in training if p.covered and p.features.get("top_target")]
        print(f"training pairs scored: {len(train)} ({sum(labels)} regulated); missing: {found}")
        print(f"coverage: {tc.coverage(training, found)['arms']}")
        w = tc.fit(train, tc.SWEEP_FEATURES)
        cal = tc.calibration(
            tc.predict(w, train, tc.SWEEP_FEATURES),
            labels,
            tc.BINS,
            [tc.stratum(p) for p in train],
            tc.stratum_coverage(training),
        )
        print(
            f"in-sample ECE {cal['reliability']['ece']}, AUPRC {cal['auprc']}, "
            f"bins {cal['reliability']['bins_consistent']}/{cal['reliability']['bins']}, "
            f"HL {cal['reliability']['hosmer_lemeshow']}"
        )
        show_reliability(cal["reliability"]["rows"])
        loco = [0.0] * len(train)
        for c in sorted({p.chrom for p in train}):
            rows = [p for p in train if p.chrom != c]
            wc = tc.fit(rows, tc.SWEEP_FEATURES)
            idx = [i for i, p in enumerate(train) if p.chrom == c]
            for i, v in zip(idx, tc.predict(wc, [train[i] for i in idx], tc.SWEEP_FEATURES), strict=True):
                loco[i] = v
        lc = tc.calibration(
            loco, labels, tc.BINS, [tc.stratum(p) for p in train], tc.stratum_coverage(training)
        )
        print(
            f"leave-chromosome-out ECE {lc['reliability']['ece']}, Brier {lc['brier']}, "
            f"AUPRC {lc['auprc']}, bins {lc['reliability']['bins_consistent']}/{lc['reliability']['bins']}, "
            f"HL {lc['reliability']['hosmer_lemeshow']}"
        )
        show_reliability(lc["reliability"]["rows"])
        raw = tc.calibration(tc.raw_confidence(train), labels, tc.BINS)
        print(
            f"effect size as a confidence: ECE {raw['reliability']['ece']}, Brier {raw['brier']}, "
            f"bins {raw['reliability']['bins_consistent']}/{raw['reliability']['bins']}, "
            f"HL {raw['reliability']['hosmer_lemeshow']}"
        )
        for row in tc.by_drop_band(top):
            print(f"    drop {row['drop']:>20}: {row['regulated']}/{row['pairs']} = {row['rate']}")
        print(f"training only, nothing held out scored ({time.time() - t0:.0f} s)")
        return 0

    result = tc.score(training, heldout, table)
    print(f"pre-registered: {result['preregistered']}")
    print(f"populations: {result['populations']}")
    for name, cal in result["heldout_k562"].items():
        r = cal["reliability"]
        print(
            f"  held-out K562 {name:30s} ECE {r['ece']:.4f}  MCE {r['mce']:.4f}  "
            f"Brier {cal['brier']:.5f}  bins inside {r['bins_consistent']}/{r['bins']}  "
            f"AUPRC {cal['auprc']}"
        )
    show_reliability(result["heldout_k562"]["calibrated, sweep features"]["reliability"]["rows"])
    print(f"checks: {result['checks']}")
    print(f"verdict: {result['verdict']}")
    for band in result["measured_by_drop_band"]["K562 pooled (both, after the held-out test)"]:
        print(
            f"  predicted target, drop {band['drop']:>20}: {band['regulated']}/{band['pairs']} "
            f"= {band['rate']} {band['ci95']}"
        )

    chroms = tuple(c for c in args.chroms.split(",") if c) or tc.CHROMS
    weights = {
        "sweep": tc.fit(tc.scored(training), tc.SWEEP_FEATURES),
        "target": tc.fit([p for p in tc.scored(training) if p.features["top_target"]], tc.TARGET_FEATURES),
    }
    swept = tc.sweep(weights, chroms, progress=lambda m: print(f"  {m}"))
    result["genome_wide"] = swept["genome_wide"]
    result["per_chromosome"] = swept["per_chromosome"]
    for label, v in swept["genome_wide"].items():
        print(f"  sweep, {label}: {v['targets']} targets, bands {v['predicted_target_bands']}")
        print(f"    missing: {v['missing']}")
    path = save_result("target_calibration", result)
    print(f"saved {path} ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
