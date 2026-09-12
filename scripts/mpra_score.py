# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured enhancer activity (ENCODE4 lentiMPRA, K562 / HepG2 / WTC11) against GenomeOS, one chromosome.

    uv run python scripts/mpra_score.py --chrom chr21             # registry, reader, constraint, model
    uv run python scripts/mpra_score.py --chrom chr21 --no-model  # without AlphaGenome

Saves data/results/mpra_<chrom>.json.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

from genomeos.attribution import mpra
from genomeos.attribution.constraint import phylop_over_blocks
from genomeos.genome.reader import PeakIndex, load_peaks
from genomeos.genome.regulatory import load_ccres
from genomeos.jobs import heartbeat
from genomeos.predict import status
from genomeos.results import save_result


def beat(name: str, msg: str | None = None) -> None:
    """Proof of life for this chromosome's job and for the genome job that may be running it."""
    heartbeat(name)
    if os.environ.get("GENOMEOS_JOB"):
        heartbeat(os.environ["GENOMEOS_JOB"])
    if msg:
        print(msg, flush=True)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chrom", default="chr21")
    ap.add_argument("--no-model", action="store_true")
    args = ap.parse_args()
    chrom, name = args.chrom, f"mpra_{args.chrom}"
    elements = mpra.load(chrom)
    if not elements:
        print(f"{chrom}: no lentiMPRA elements")
        return 1
    ccres = load_ccres(chrom)
    peaks = {cell: PeakIndex(load_peaks(cell, chrom)) for cell in mpra.MODELLED if load_peaks(cell, chrom)}
    intervals, keep = [], []
    last_end = -1
    for e in elements:  # the bigWig summariser wants non-overlapping intervals
        if e.start >= last_end:
            intervals.append((e.start, e.end))
            keep.append(e)
            last_end = e.end
    beat(name)
    t0 = time.time()
    stats, cost = phylop_over_blocks(
        chrom,
        intervals,
        progress=lambda d, n, b: beat(name, f"{chrom}: phyloP {d}/{n} blocks, {b / 1e6:.0f} MB"),
    )
    print(f"{chrom}: phyloP over {len(keep):,} elements in {time.time() - t0:.0f} s ({cost})", flush=True)
    rows = mpra.annotate(keep, ccres, peaks, stats)
    model_cost = None
    if not args.no_model:
        st = status()
        if not st["enabled"]:
            print(f"AlphaGenome is disabled: {st['reason']}; scoring without the model")
        else:
            from genomeos.genome import IndexedGenome
            from genomeos.predict.chromatin_tracks import PredictedDnase
            from genomeos.predict.splice_sites import client_factory

            g = IndexedGenome(f"data/reference/{chrom}.fa")
            pd_ = PredictedDnase(chrom, g.lengths[chrom])
            g.close()

            def progress(msg: str) -> None:
                beat(name)
                print(msg, flush=True)

            means = pd_.means([(e.start, e.end, e.key) for e in keep], client_factory, progress)
            mpra.attach_predicted(rows, means)
            model_cost = {"requests": pd_.requests, "seconds": round(pd_.seconds, 1)}
    summary = mpra.summarise(rows)
    out = {
        "chrom": chrom,
        "elements": len(elements),
        "annotated": len(rows),
        "reader_cells": sorted(peaks),
        "summary": summary,
        "rows": rows,
        "phylop_cost": cost,
        "model_cost": model_cost,
        "evidence": f"{mpra.EVIDENCE}; experimental: ENCODE DNase peaks (reader); measured: Zoonomia phyloP; "
        "curated: ENCODE cCREs; predicted: AlphaGenome DNase per cell line",
    }
    save_result(name, out)
    for cell, s in summary["cells"].items():
        same = (s.get("predicted_dnase_vs_activity") or {}).get("same_cell") or {}
        rd = (s.get("reader_open_predicts_active") or {}).get("same_cell") or {}
        print(
            f"done {cell}: {s['elements']} elements, active {s['active']}; "
            f"reader precision {rd.get('precision')} recall {rd.get('recall')}; "
            f"predicted DNase rho {same.get('spearman')} "
            f"(top quartile active {same.get('active_in_top_quartile')}, "
            f"bottom {same.get('active_in_bottom_quartile')})"
        )
    if "specific_elements" in summary:
        print(f"specific: {summary['specific_elements']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
