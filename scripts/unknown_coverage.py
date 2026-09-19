# SPDX-License-Identifier: AGPL-3.0-or-later
"""What is actually MEASURED over the unknown space, block by block. No model, no requests.

    uv run python scripts/unknown_coverage.py [--chroms chr21,chr22] [--no-save]

Three lanes reached the same wall from different directions on 2026-09-17: the deletion sweep finished
and its per-tier readings did not survive standardisation, a calibration of those predictions failed on
the prevalence of the screens that measured it, and a motif model that transfers to VISTA cannot beat
its own shuffle swept blind. What the unknown space is short of is not another model over it. The
transfer lane put a number on it — 161 of the 882 blocks of the real unknown hold a measured element
and 721 hold none — and this is that number for every tier, with the shape an assay would have to have.

It reads the measured layers this project already holds, each through its own loader:

- **lentiMPRA** (ENCODE4 K562/HepG2/WTC11): the elements an episomal reporter measured;
- **VISTA**: in-vivo enhancer tests, positive and negative alike, since a negative is a measurement;
- **CRISPRi** (the ENCODE enhancer-gene benchmark, training and held-out): the elements a screen
  perturbed in their own chromosomes.

For every UNKNOWN block it reports which layers touch it and how much of its sequence they cover, and
for every tier it reports the blocks and megabases that nothing has measured. That last number is the
one an assay is designed against.
"""

from __future__ import annotations

import argparse
import bisect
import time
from collections import defaultdict
from statistics import median
from typing import Any

from genomeos.attribution import crispri, mpra, organise, vista
from genomeos.results import save_result

CRISPRI_FILES = (
    "EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz",
    "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz",
)


def label_of(block: dict[str, Any]) -> str:
    """The same labelling the matched reading uses: tier, and for the real unknown its case."""
    tier = block.get("tier")
    if tier == "constrained_unknown":
        return (
            "constrained_unknown_copy"
            if block.get("copy")
            else f"real_unknown_{block.get('case') or 'uncased'}"
        )
    return str(tier)


def measured_spans(chrom: str) -> dict[str, list[tuple[int, int]]]:
    """Every measured interval on one chromosome, by the assay that measured it."""
    spans: dict[str, list[tuple[int, int]]] = {}
    try:
        spans["lentimpra"] = sorted((e.start, e.end) for e in mpra.load(chrom))
    except Exception:  # noqa: BLE001 - a layer that is not fetched is reported as absent, not fatal
        spans["lentimpra"] = []
    try:
        spans["vista"] = sorted((e.start, e.end) for e in vista.load_loci(chrom))
    except Exception:  # noqa: BLE001
        spans["vista"] = []
    crispri_spans: list[tuple[int, int]] = []
    for name in CRISPRI_FILES:
        try:
            crispri_spans.extend((p.start, p.end) for p in crispri.load(name) if p.chrom == chrom)
        except Exception:  # noqa: BLE001
            continue
    spans["crispri"] = sorted(set(crispri_spans))
    return spans


def overlap_bp(spans: list[tuple[int, int]], lo: int, hi: int) -> int:
    """Bases of a block covered by these intervals, counting any base once."""
    if not spans:
        return 0
    starts = [s for s, _ in spans]
    i = max(0, bisect.bisect_left(starts, lo) - 1)
    covered, last = 0, lo
    while i < len(spans) and spans[i][0] < hi:
        s, e = max(spans[i][0], lo), min(spans[i][1], hi)
        if e > s:
            s = max(s, last)
            if e > s:
                covered += e - s
                last = e
        i += 1
    return covered


def read_chromosome(chrom: str) -> list[dict[str, Any]]:
    """Every UNKNOWN block of one chromosome with the measured layers that touch it."""
    spans = measured_spans(chrom)
    rows = []
    for b in organise.blocks(chrom):
        lo, hi = b["start"], b["end"]
        by_assay = {name: overlap_bp(v, lo, hi) for name, v in spans.items()}
        rows.append(
            {
                "chrom": chrom,
                "block": f"{chrom}:{lo}-{hi}",
                "label": label_of(b),
                "length": b["length"],
                "measured_bp": by_assay,
                "assays": sorted(k for k, v in by_assay.items() if v),
            }
        )
    return rows


def summarise(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Per label: how many blocks anything has measured, and the megabases nothing has."""
    by_label: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_label[r["label"]].append(r)
    out = {}
    for label, v in sorted(by_label.items()):
        touched = [r for r in v if r["assays"]]
        untouched = [r for r in v if not r["assays"]]
        out[label] = {
            "blocks": len(v),
            "blocks_measured": len(touched),
            "blocks_untouched": len(untouched),
            "share_measured": round(len(touched) / len(v), 4) if v else None,
            "bp": sum(r["length"] for r in v),
            "measured_bp": sum(sum(r["measured_bp"].values()) for r in v),
            "untouched_mb": round(sum(r["length"] for r in untouched) / 1e6, 2),
            "by_assay_blocks": {
                assay: sum(1 for r in v if r["measured_bp"].get(assay))
                for assay in ("lentimpra", "vista", "crispri")
            },
            "median_untouched_block_kb": (
                round(median([r["length"] for r in untouched]) / 1000, 2) if untouched else None
            ),
        }
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--no-save", action="store_true", help="print without overwriting the genome's result")
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]

    t0 = time.time()
    rows: list[dict[str, Any]] = []
    for chrom in chroms:
        got = read_chromosome(chrom)
        rows.extend(got)
        print(
            f"{chrom}: {len(got)} blocks, {sum(1 for r in got if r['assays'])} with a measurement", flush=True
        )

    per_label = summarise(rows)
    real = [r for r in rows if r["label"].startswith("real_unknown_")]
    out = {
        "result": "unknown_coverage",
        "chromosomes": chroms,
        "assays": ["lentimpra", "vista", "crispri"],
        "by_label": per_label,
        "real_unknown": {
            "blocks": len(real),
            "blocks_measured": sum(1 for r in real if r["assays"]),
            "untouched_mb": round(sum(r["length"] for r in real if not r["assays"]) / 1e6, 2),
            "measured_bp": sum(sum(r["measured_bp"].values()) for r in real),
            "bp": sum(r["length"] for r in real),
        },
        "reading": (
            "a block is 'measured' if any of the three assays overlaps it by a single base, which is the "
            "most generous definition available and still leaves most of the space untouched. Coverage "
            "is not evidence of function; it is the precondition for asking about it"
        ),
        "seconds": round(time.time() - t0, 1),
    }
    if args.no_save:
        print("not saved (--no-save)")
    else:
        print(f"saved {save_result(out['result'], out)}")
    print(f"\n{'label':26} {'blocks':>7} {'measured':>9} {'share':>7} {'untouched Mb':>13}")
    for label, v in per_label.items():
        print(
            f"{label:26} {v['blocks']:7} {v['blocks_measured']:9} "
            f"{v['share_measured']:7} {v['untouched_mb']:13}"
        )
    r = out["real_unknown"]
    print(
        f"\nreal unknown: {r['blocks_measured']} of {r['blocks']} blocks measured, "
        f"{r['measured_bp'] / 1e6:.2f} Mb of {r['bp'] / 1e6:.2f} Mb covered, {r['untouched_mb']} Mb untouched"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
