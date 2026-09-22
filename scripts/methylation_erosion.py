#!/usr/bin/env python3
"""Does the methylation engine's solo-WCGW prediction hold in the WGBS this project already has?

The engine (`genomeos/runtime/methylation.py`) runs the published mechanism: replication makes a
methylated CpG hemimethylated, UHRF1 + DNMT1 restore it with a fidelity below 1, and both
maintenance and de novo methylation get *worse* as the local CpG density falls. Run it and an
isolated CpG erodes while a clustered one holds. That is a prediction about real methylomes, and
the project holds eight ENCODE WGBS biosamples to check it against.

This is the application side of the test: it reads sequence and the binned WGBS caches, which the
engine is not allowed to know about, and calls the engine only through its parameters and its
trajectory. The claim below is registered in code before any methylation value is read
(the class thresholds are fixed from *sequence* alone, with `calibrate`).

    python scripts/methylation_erosion.py calibrate [--chrom chr21 ...]   sequence classes only
    python scripts/methylation_erosion.py test [--chrom chr21 chr22 ...]  the claim, and the fit
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from genomeos.coords import Locus  # noqa: E402
from genomeos.genome import epigenome as ep  # noqa: E402
from genomeos.genome.index import IndexedGenome  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from genomeos.runtime.methylation import (  # noqa: E402
    CpGContext,
    MethylationParams,
    divisions_for_level,
    expected_trajectory,
    fidelity_for_level,
)

# ------------------------------------------------------------------------------------------
# The claim, registered before the data is read
# ------------------------------------------------------------------------------------------

LONG_CULTURED = ("K562", "HepG2", "GM12878")
INTACT = ("H1", "hepatocyte", "CD14-positive monocyte")
ALSO_REPORTED = ("SK-N-SH", "IMR-90")

CHROMS = ("chr19", "chr20", "chr21", "chr22")

CLAIM = {
    "name": "solo_wcgw_erosion_ordering",
    "chromosomes": list(CHROMS),
    "registered": "2026-09-16",
    "engine": "genomeos/runtime/methylation.py (CpG methylation state machine)",
    "mechanism": (
        "Maintenance (UHRF1 + DNMT1) and de novo methylation (DNMT3A/B) are less efficient at an "
        "isolated CpG than in a cluster, so each division loses more methylation at a solo CpG "
        "than at a dense one and the two run apart as divisions accumulate."
    ),
    "prediction": (
        "In WGBS, 200 bp windows whose CpGs are all isolated and A/T-flanked (solo-WCGW, Zhou et "
        "al. 2018) are less methylated than windows whose CpGs all have a close neighbour, and the "
        "difference is LARGER in the long-cultured lines (K562, HepG2, GM12878) than in the three "
        "biosamples with an intact methylome (H1, hepatocyte, CD14-positive monocyte)."
    ),
    "decision_rule": (
        "The ordering holds if min(gap over K562, HepG2, GM12878) > max(gap over H1, hepatocyte, "
        "monocyte), where gap = dense methylation minus solo-WCGW methylation in points. Reported "
        "for the whole chromosome set and, as the region control, within 50 kb windows that hold "
        "both classes. The gap relative to the dense level (gap / dense) is a secondary metric, "
        "registered here because a line whose whole methylome sits near zero cannot show a large "
        "gap in points."
    ),
    "falsified_if": (
        "any long-cultured line shows a smaller gap than any intact biosample, or the gap is not "
        "positive where the engine says it must be."
    ),
    "long_cultured": list(LONG_CULTURED),
    "intact": list(INTACT),
    "also_reported_outside_the_claim": list(ALSO_REPORTED),
}

# Window classes, fixed from sequence alone (`calibrate`) before any methylation value is read.
SOLO_WINDOW_BP = 35  # Zhou et al. 2018: an isolated CpG has no other CpG within 35 bp
MIN_SOLO_CPGS = 2  # a solo-WCGW window: at least this many CpGs, every one of them solo-WCGW
MIN_DENSE_CPGS = 6  # a dense window: at least this many CpGs, none of them solo
MAX_DENSE_OBS_EXP = 0.6  # ... and not CpG-island-like (Gardiner-Garden & Frommer 1987), because
#                          islands are unmethylated in every cell type and would not be a control
MIN_COVERED_CALLS = 2  # a window counts only where WGBS covered 2+ CpGs at 5+ reads
REGION_BP = 50_000  # the region control: gaps computed inside one 50 kb window

# The assumptions the fit needs and the literature does not supply, stated twice so the answer can
# be read against both. Nothing here is measured; everything here is labelled inferred.
DE_NOVO_SCENARIOS = {
    "low": {"de_novo_dense": 0.01, "de_novo_solo": 0.002},
    "high": {"de_novo_dense": 0.05, "de_novo_solo": 0.01},
}
LITERATURE_FIDELITY = 0.97  # the middle of the published 0.95-0.99 (Laird 2004, Genereux 2005)
FIT_DIVISIONS = (20, 50, 100, 200, 500, 1000)
NEIGHBOUR_SCALE = 2.0  # CpGs; UNKNOWN in the literature, assumed
SOLO_CTX = CpGContext("solo_wcgw", 0, "W")
DENSE_CTX = CpGContext("dense", MIN_DENSE_CPGS + 2, "S")


# ------------------------------------------------------------------------------------------
# Sequence: which 200 bp windows are solo-WCGW and which are dense
# ------------------------------------------------------------------------------------------


@dataclass(slots=True)
class Classes:
    chrom: str
    solo: list[int]  # bin indices
    dense: list[int]
    counted: dict[str, int]


def classify(chrom: str, reference: Path) -> Classes:
    """Split a chromosome's 200 bp bins into solo-WCGW windows and dense windows, from sequence."""
    with IndexedGenome(reference / f"{chrom}.fa.gz") as g:
        length = g.lengths[chrom]
        seq = str(g.fetch(Locus(chrom, 0, length))).upper()
    positions: list[int] = []
    flank_w: list[bool] = []
    i = seq.find("CG")
    while i != -1:
        if i >= 1 and i + 2 < length:
            positions.append(i)
            flank_w.append(seq[i - 1] in "AT" and seq[i + 2] in "AT")
        i = seq.find("CG", i + 1)
    # neighbours within +/- 35 bp, from the sorted position list
    n = len(positions)
    solo_flag = [False] * n
    lo = 0
    for k in range(n):
        while positions[k] - positions[lo] > SOLO_WINDOW_BP:
            lo += 1
        hi = k
        while hi + 1 < n and positions[hi + 1] - positions[k] <= SOLO_WINDOW_BP:
            hi += 1
        solo_flag[k] = (hi - lo) == 0  # only itself inside the window
    nb = length // ep.BIN + 1
    cpgs = [0] * nb
    solo_wcgw = [0] * nb
    solo_any = [0] * nb
    for k, pos in enumerate(positions):
        b = pos // ep.BIN
        cpgs[b] += 1
        if solo_flag[k]:
            solo_any[b] += 1
            solo_wcgw[b] += flank_w[k]
    solo: list[int] = []
    dense: list[int] = []
    for b in range(nb):
        if cpgs[b] < MIN_SOLO_CPGS:
            continue
        start, end = b * ep.BIN, min((b + 1) * ep.BIN, length)
        window = seq[start:end]
        if "N" in window:
            continue
        if cpgs[b] >= MIN_SOLO_CPGS and solo_wcgw[b] == cpgs[b]:
            solo.append(b)
        elif cpgs[b] >= MIN_DENSE_CPGS and solo_any[b] == 0:
            c, gc = window.count("C"), window.count("G")
            obs_exp = cpgs[b] * len(window) / (c * gc) if c and gc else math.inf
            if obs_exp < MAX_DENSE_OBS_EXP:
                dense.append(b)
    return Classes(
        chrom,
        solo,
        dense,
        {
            "bins": nb,
            "cpgs": n,
            "solo_cpgs": sum(solo_flag),
            "solo_wcgw_cpgs": sum(1 for k in range(n) if solo_flag[k] and flank_w[k]),
            "solo_windows": len(solo),
            "dense_windows": len(dense),
        },
    )


# ------------------------------------------------------------------------------------------
# The measurement
# ------------------------------------------------------------------------------------------


def _weighted(cols: list, bins: list[int]) -> tuple[float | None, int, int]:
    """Call-weighted methylation over these bins: (fraction, windows measured, covered calls)."""
    frac = covered = windows = 0
    for b in bins:
        if b >= len(cols[1]):
            continue
        cov = cols[1][b]
        if cov < MIN_COVERED_CALLS:
            continue
        covered += cov
        frac += cols[2][b]
        windows += 1
    return (frac / covered / 1000 if covered else None), windows, covered


def _region_gap(cols: list, solo: list[int], dense: list[int]) -> tuple[float | None, int]:
    """Mean of (dense - solo) inside 50 kb windows holding both classes, weighted by solo windows."""
    per = ep.BIN
    step = REGION_BP // per
    groups: dict[int, dict[str, list[int]]] = {}
    for name, bins in (("solo", solo), ("dense", dense)):
        for b in bins:
            groups.setdefault(b // step, {"solo": [], "dense": []})[name].append(b)
    num = den = 0.0
    regions = 0
    for g in groups.values():
        s, sn, _ = _weighted(cols, g["solo"])
        d, dn, _ = _weighted(cols, g["dense"])
        if s is None or d is None or sn == 0 or dn == 0:
            continue
        num += (d - s) * sn
        den += sn
        regions += 1
    return (num / den if den else None), regions


def measure(chroms: list[str], reference: Path) -> dict:
    manifest = ep.load_manifest()
    if not manifest:
        raise SystemExit("no epigenome manifest: run `genomeos epigenome manifest` first")
    cells = [c for c, r in manifest["cell_types"].items() if "cpg" in r.get("methylation", {})]
    classes = {c: classify(c, reference) for c in chroms}
    out: dict[str, dict] = {}
    for cell in cells:
        acc = {
            "solo_frac": 0.0,
            "solo_cov": 0,
            "solo_windows": 0,
            "dense_frac": 0.0,
            "dense_cov": 0,
            "dense_windows": 0,
            "region_num": 0.0,
            "region_den": 0.0,
            "regions": 0,
            "chroms": [],
        }
        for chrom in chroms:
            cols = ep.load_methylation_profile(cell, chrom)
            if cols is None:
                continue
            cl = classes[chrom]
            for name, bins in (("solo", cl.solo), ("dense", cl.dense)):
                frac = covered = windows = 0
                for b in bins:
                    if b < len(cols[1]) and cols[1][b] >= MIN_COVERED_CALLS:
                        covered += cols[1][b]
                        frac += cols[2][b]
                        windows += 1
                acc[f"{name}_frac"] += frac
                acc[f"{name}_cov"] += covered
                acc[f"{name}_windows"] += windows
            gap, regions = _region_gap(cols, cl.solo, cl.dense)
            if gap is not None:
                n_solo = sum(1 for b in cl.solo if b < len(cols[1]) and cols[1][b] >= MIN_COVERED_CALLS)
                acc["region_num"] += gap * n_solo
                acc["region_den"] += n_solo
                acc["regions"] += regions
            acc["chroms"].append(chrom)
        if not acc["solo_cov"] or not acc["dense_cov"]:
            continue
        solo = acc["solo_frac"] / acc["solo_cov"] / 1000
        dense = acc["dense_frac"] / acc["dense_cov"] / 1000
        out[cell] = {
            "solo_wcgw": round(solo, 4),
            "dense": round(dense, 4),
            "gap": round(dense - solo, 4),
            "gap_over_dense": round((dense - solo) / dense, 4) if dense else None,
            "gap_within_50kb": (
                round(acc["region_num"] / acc["region_den"], 4) if acc["region_den"] else None
            ),
            "regions_with_both": acc["regions"],
            "solo_windows_measured": acc["solo_windows"],
            "dense_windows_measured": acc["dense_windows"],
            "solo_calls": acc["solo_cov"],
            "dense_calls": acc["dense_cov"],
            "chromosomes": acc["chroms"],
        }
    return {"cells": out, "classes": {c: classes[c].counted for c in chroms}}


# ------------------------------------------------------------------------------------------
# What the engine has to assume to reproduce those numbers
# ------------------------------------------------------------------------------------------


def params(scenario: str, solo_f: float = 0.9, dense_f: float = LITERATURE_FIDELITY) -> MethylationParams:
    return MethylationParams(
        maintenance_fidelity_dense=dense_f,
        maintenance_fidelity_solo_wcgw=solo_f,
        maintenance_fidelity_solo_scgs=solo_f,
        de_novo_dense=DE_NOVO_SCENARIOS[scenario]["de_novo_dense"],
        de_novo_solo=DE_NOVO_SCENARIOS[scenario]["de_novo_solo"],
        tet_erasure=0.0,
        neighbour_scale=NEIGHBOUR_SCALE,
        assumptions=(f"de novo scenario {scenario}", "no active TET turnover", "tet_erasure = 0"),
    )


def fit(cells: dict, start: dict[str, float], scenario: str) -> dict:
    """For each biosample: the divisions its dense level implies at literature fidelity, and the
    solo-versus-dense fidelity the two levels imply at a range of division counts. Inferred."""
    p = params(scenario)
    out: dict[str, dict] = {}
    for cell, row in cells.items():
        implied = divisions_for_level(row["dense"], p, DENSE_CTX, initial_level=start["dense"])
        per_n: dict[str, dict] = {}
        for n in FIT_DIVISIONS:
            f_dense = fidelity_for_level(
                row["dense"], n, p, DENSE_CTX, knob="dense", initial_level=start["dense"]
            )
            f_solo = None
            if f_dense is not None:
                q = params(scenario, dense_f=f_dense)
                f_solo = fidelity_for_level(
                    row["solo_wcgw"], n, q, SOLO_CTX, knob="solo", initial_level=start["solo_wcgw"]
                )
            per_n[str(n)] = {
                "maintenance_fidelity_dense": round(f_dense, 5) if f_dense is not None else None,
                "maintenance_fidelity_solo_wcgw": round(f_solo, 5) if f_solo is not None else None,
                "deficit": (
                    round(f_dense - f_solo, 5) if f_dense is not None and f_solo is not None else None
                ),
            }
        out[cell] = {
            "divisions_implied_at_literature_fidelity": implied,
            "per_division_count": per_n,
        }
    return {
        "inferred": True,
        "start_level": start,
        "de_novo_scenario": {"name": scenario, **DE_NOVO_SCENARIOS[scenario]},
        "literature_fidelity_dense": LITERATURE_FIDELITY,
        "neighbour_scale": NEIGHBOUR_SCALE,
        "note": (
            "Inferred, not measured. The start level is the intact reference's own measured level, "
            "the de novo rates and the neighbour scale are assumptions (UNKNOWN in the literature), "
            "and TET turnover is set to zero, so every number here is conditional on those four."
        ),
        "per_cell": out,
    }


def verdict(cells: dict) -> dict:
    def gaps(names: tuple[str, ...], key: str) -> dict[str, float]:
        return {n: cells[n][key] for n in names if n in cells and cells[n].get(key) is not None}

    out = {}
    for key in ("gap", "gap_over_dense", "gap_within_50kb"):
        cultured, intact = gaps(LONG_CULTURED, key), gaps(INTACT, key)
        if len(cultured) < len(LONG_CULTURED) or len(intact) < len(INTACT):
            out[key] = {"holds": None, "reason": "a biosample of the claim is missing"}
            continue
        lo, hi = min(cultured.values()), max(intact.values())
        out[key] = {
            "holds": bool(lo > hi),
            "min_long_cultured": round(lo, 4),
            "max_intact": round(hi, 4),
            "margin": round(lo - hi, 4),
            "long_cultured": cultured,
            "intact": intact,
            "every_gap_positive": all(v > 0 for v in {**cultured, **intact}.values()),
        }
    return out


# ------------------------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("command", choices=("calibrate", "test"))
    ap.add_argument("--chrom", nargs="*", default=list(CHROMS))
    ap.add_argument("--reference", default="data/reference")
    ap.add_argument("--reference-start", default="H1", help="the biosample whose levels start the fit")
    args = ap.parse_args(argv)
    reference = Path(args.reference)

    if args.command == "calibrate":
        for chrom in args.chrom:
            print(chrom, json.dumps(classify(chrom, reference).counted))
        return 0

    m = measure(list(args.chrom), reference)
    cells = m["cells"]
    if args.reference_start not in cells:
        raise SystemExit(f"{args.reference_start} has no WGBS in the cache; cannot start the fit")
    start = {
        "dense": cells[args.reference_start]["dense"],
        "solo_wcgw": cells[args.reference_start]["solo_wcgw"],
    }
    payload = {
        "claim": CLAIM,
        "chromosomes": list(args.chrom),
        "window_classes": {
            "bin_bp": ep.BIN,
            "solo_window_bp": SOLO_WINDOW_BP,
            "min_solo_cpgs": MIN_SOLO_CPGS,
            "min_dense_cpgs": MIN_DENSE_CPGS,
            "max_dense_obs_exp": MAX_DENSE_OBS_EXP,
            "min_covered_calls": MIN_COVERED_CALLS,
            "region_bp": REGION_BP,
            "counts": m["classes"],
        },
        "measured": cells,
        "verdict": verdict(cells),
        "inferred_fit": {
            name: fit(cells, start, name)
            for name in DE_NOVO_SCENARIOS  # both de novo scenarios
        },
        "engine_prediction": {
            "note": "the engine's own trajectory under the demo program's assumptions",
            "divisions": 100,
            "solo_wcgw": round(
                expected_trajectory(SOLO_CTX, params("high", solo_f=0.90), 100, initial_level=0.95)[-1], 4
            ),
            "dense": round(
                expected_trajectory(DENSE_CTX, params("high", solo_f=0.90), 100, initial_level=0.95)[-1], 4
            ),
        },
    }
    path = save_result("methylation_erosion", payload)
    print(json.dumps({"measured": cells, "verdict": payload["verdict"]}, indent=2))
    print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
