# SPDX-License-Identifier: AGPL-3.0-or-later
"""How much a 1 kb track's whole-bin counting overstates an absolute figure, on chr21.

    uv run python scripts/bigwig_bin_shift.py [--chrom chr21]

The bigWig reader counts every base of every bin an interval touches. For Zoonomia phyloP, one
value per base, that is exact. For gnomAD Gnocchi, one value per kilobase, it is the share of
touched kilobases -- the honest reading of a measurement with no finer resolution, and the reading
`variation.py` documents -- but it overstates any absolute count about the interval itself, worst
for intervals shorter than a bin, where a single kilobase is credited in full.

`IntervalStats` now carries both readings from one pass: `bases`/`above` count touched bins,
`overlap_bases`/`overlap_above` count only the bases inside the interval. This script reports the
two columns side by side over one chromosome's UNKNOWN blocks and its registry elements, so the
shift is published before any result is recomputed against it. The chromosome was named before the
numbers were read: chr21, the chromosome every reading in this project is scored on first.

Saves data/results/bigwig_bin_shift.json. One read of the track, both columns.
"""

from __future__ import annotations

import argparse
import gzip
import time
from pathlib import Path
from typing import Any

from genomeos.attribution import organise, variation
from genomeos.results import save_result

RESULTS = Path("data/results")
BIN = 1000  # the Gnocchi track's resolution, for the "shorter than a bin" split


def registry_elements(chrom: str, results: Path = RESULTS) -> list[tuple[int, int]]:
    """Non-overlapping registry element spans of one chromosome, from the committed cCRE track."""
    p = results / f"ccres_{chrom}.bed.gz"
    if not p.exists():
        return []
    spans = []
    with gzip.open(p, "rt") as f:
        for line in f:
            parts = line.split("\t")
            if len(parts) < 3 or parts[0] != chrom:
                continue
            spans.append((int(parts[1]), int(parts[2])))
    spans.sort()
    out: list[tuple[int, int]] = []
    for s, e in spans:  # the summariser wants non-overlapping intervals
        if out and s < out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
            continue
        out.append((s, e))
    return out


def columns(intervals: list[tuple[int, int]], stats: list[Any]) -> dict[str, Any]:
    """Both readings over one set of intervals, and the same split for intervals under one bin."""
    rows = [(e - s, st) for (s, e), st in zip(intervals, stats, strict=True) if st and st.bases]
    if not rows:
        return {"intervals": 0}

    def fold(sel: list[tuple[int, Any]]) -> dict[str, Any]:
        if not sel:
            return {"intervals": 0}
        length = sum(n for n, _ in sel)
        touched = sum(st.bases for _, st in sel)
        overlap = sum(st.overlap_bases for _, st in sel)
        above = sum(st.above for _, st in sel)
        over_above = sum(st.overlap_above for _, st in sel)
        return {
            "intervals": len(sel),
            "interval_bp": length,
            "touched_bin_bases": touched,
            "overlap_bases": overlap,
            "overstatement": round(touched / overlap, 3) if overlap else None,
            "constrained_bp_touched_bins": above,
            "constrained_bp_overlap": over_above,
            "constrained_bp_overstatement": round(above / over_above, 3) if over_above else None,
            "fraction_above_touched_bins": round(above / touched, 4) if touched else None,
            "fraction_above_overlap": round(over_above / overlap, 4) if overlap else None,
            "overlap_bases_never_exceed_length": overlap <= length,
        }

    return {
        "all": fold(rows),
        "shorter_than_one_bin": fold([r for r in rows if r[0] < BIN]),
        "at_least_one_bin": fold([r for r in rows if r[0] >= BIN]),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr21")
    args = ap.parse_args()
    chrom = args.chrom
    t0 = time.time()
    blocks = [(b["start"], b["end"]) for b in organise.blocks(chrom)]
    elements = registry_elements(chrom)
    out: dict[str, Any] = {
        "chromosome_named_in_advance": chrom,
        "track": variation.GNOCCHI_URL,
        "track_resolution_bp": BIN,
        "threshold": variation.GNOCCHI_THRESHOLD,
        "reading": (
            "touched-bin counting is the share of kilobases an interval touches and is what "
            "variation.py documents; overlap counting is the interval's own bases and is what an "
            "absolute figure needs. Ratios are almost unchanged, absolute counts are not."
        ),
    }
    for name, intervals in (("unknown_blocks", blocks), ("registry_elements", elements)):
        if not intervals:
            out[name] = {"refused": "no intervals of this kind cached for this chromosome"}
            continue
        stats, cost = variation.gnocchi_over(chrom, intervals)
        out[name] = {"cost": cost, **columns(intervals, stats)}
        print(f"{chrom} {name}: {len(intervals)} intervals, {cost}", flush=True)
        for key in ("all", "shorter_than_one_bin", "at_least_one_bin"):
            c = out[name].get(key, {})
            if c.get("intervals"):
                print(
                    f"  {key:22s} n {c['intervals']:6d}  length {c['interval_bp']:>10,}  "
                    f"touched {c['touched_bin_bases']:>10,}  overlap {c['overlap_bases']:>10,}  "
                    f"x{c['overstatement']}  constrained {c['constrained_bp_touched_bins']:>9,} "
                    f"-> {c['constrained_bp_overlap']:>9,}  fraction {c['fraction_above_touched_bins']} "
                    f"-> {c['fraction_above_overlap']}",
                    flush=True,
                )
    out["seconds"] = round(time.time() - t0, 1)
    path = save_result("bigwig_bin_shift", out)
    print(f"-> {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
