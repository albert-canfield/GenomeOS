# SPDX-License-Identifier: AGPL-3.0-or-later
"""How much of the untouched 13.77 Mb a reporter assay could not measure in principle.

    uv run python scripts/measurability.py [--chroms chr21,chr22] [--no-save]

`scripts/unknown_coverage.py` found that 680 of the 882 real-unknown blocks -- 13.77 Mb -- have never
been touched by lentiMPRA, VISTA or the CRISPRi benchmark, and that at 300 bp this tiles into 45,900
oligos. 45,900 is arithmetic: 13,770,000 divided by 300. This script measures the part of it a library
could actually order and read back, per block and pooled, and charges every excluded base to a named
reason. It is the opposite of a prediction over the unknown space: it sharpens an artefact of the
instrument so that the design excludes or flags sequence before oligos are paid for.

`genomeos.attribution.measurability` holds the definitions -- what counts as unmeasurable, in two
families, with every threshold a named constant and every arbitrary one reported at three values. This
script only joins the local layers: the block organiser, RepeatMasker, the curated segmental
duplications, the bgzip-indexed reference and GENCODE for the distance to a coding TSS.

Nothing is fetched. Every block that could not be assessed is counted with the layer that was missing.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path
from typing import Any

from genomeos.attribution import measurability as meas
from genomeos.results import save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
REFERENCE = Path("data/reference")
PEER_READING = {"blocks": 680, "mb": 13.77, "oligos": 45_900}  # unknown_coverage, 2026-09-17


def coding_tss(chrom: str) -> list[int]:
    """Every protein-coding transcription start on the chromosome, sorted; [] when unannotated."""
    from genomeos.genome import Annotation, default_gencode

    gff = default_gencode({chrom})
    if gff is None:
        return []
    ann = Annotation.from_gff3(gff, {chrom})
    return sorted(
        (g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start)
        for g in ann.genes.values()
        if g.type == "protein_coding" and g.locus.chrom == chrom
    )


def mpra_loci(chrom: str, cap: int) -> list[tuple[int, int]]:
    """`OLIGO` bp windows centred on lentiMPRA elements: sequence a reporter already measured."""
    from genomeos.attribution import mpra

    try:
        elements = sorted({(e.start, e.end) for e in mpra.load(chrom)})
    except Exception:  # noqa: BLE001 - an unfetched layer contributes nothing to the check
        return []
    out = []
    for s, e in sorted(meas.sample(elements, cap)):
        lo = max(0, (s + e) // 2 - meas.OLIGO // 2)
        out.append((lo, lo + meas.OLIGO))
    return out


def _skipped(chrom: str, blocks: list[dict[str, Any]], missing: str) -> dict[str, Any]:
    return {
        "rows": [],
        "mpra": [],
        "blocks_seen": len(blocks),
        "unassessed": [
            {"block": f"{chrom}:{b['start']}-{b['end']}", "chrom": chrom, "missing": missing} for b in blocks
        ],
    }


def run_chromosome(chrom: str, mpra_cap: int) -> dict[str, Any]:
    """Every untouched real-unknown block of one chromosome, assessed, plus the calibration windows."""
    from genomeos.coords import Locus
    from genomeos.genome import IndexedGenome, reference_fasta

    unassessed: list[dict[str, str]] = []
    blocks = meas.untouched_blocks(chrom)
    fasta = reference_fasta(chrom, REFERENCE)
    if not fasta.exists():
        return _skipped(chrom, blocks, "reference")
    calibration_loci = mpra_loci(chrom, mpra_cap) if mpra_cap else []
    regions = [(b["start"], b["end"]) for b in blocks] + calibration_loci
    by_class = meas.repeat_spans_by_class(chrom, regions)
    dups = meas.segdup_spans(chrom)
    if by_class is None or dups is None:
        return _skipped(chrom, blocks, "rmsk" if by_class is None else "superdups")
    groups = meas.class_groups(by_class)
    segdups, young = dups
    tss = coding_tss(chrom)
    genome = IndexedGenome(str(fasta))
    rows, mp = [], []
    try:
        for b in blocks:
            seq = str(genome.fetch(Locus(chrom, b["start"], b["end"])))
            if len(seq) != b["end"] - b["start"]:
                unassessed.append(
                    {"block": f"{chrom}:{b['start']}-{b['end']}", "chrom": chrom, "missing": "sequence"}
                )
                continue
            mid = (b["start"] + b["end"]) // 2
            row = meas.assess_block(b, seq, groups, segdups, young, meas.nearest(tss, mid))
            row["label"] = b["label"]
            rows.append(row)
        for lo, hi in calibration_loci:
            seq = str(genome.fetch(Locus(chrom, lo, hi))).upper()
            if len(seq) < meas.OLIGO:
                continue
            mp.append(
                {
                    "chrom": chrom,
                    "start": lo,
                    "tss_distance": meas.nearest(tss, (lo + hi) // 2),
                    # `non_unique_in_block` has no block to be non-unique in, so its share is 0.0 by
                    # construction here and that one reason is not calibrated; the other six are.
                    "window": meas.window_metrics(seq, lo, hi, groups, segdups, young, 0.0),
                }
            )
    finally:
        genome.close()
    if not tss:
        for r in rows:
            r["tss_distance"] = None
    return {"rows": rows, "mpra": mp, "unassessed": unassessed, "blocks_seen": len(blocks)}


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--mpra-cap", type=int, default=meas.MPRA_CALIBRATION_MAX // len(CHROMS))
    ap.add_argument("--min-usable", type=int, default=3, help="oligos below which a block cannot be tiled")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else CHROMS

    t0 = time.time()
    rows: list[dict[str, Any]] = []
    mpra_rows: list[dict[str, Any]] = []
    unassessed: list[dict[str, str]] = []
    seen = 0
    for chrom in chroms:
        try:
            got = run_chromosome(chrom, args.mpra_cap)
        except FileNotFoundError as exc:
            unassessed.append({"block": f"{chrom}:*", "chrom": chrom, "missing": str(exc)})
            print(f"{chrom}: not assessed ({exc})", flush=True)
            continue
        rows.extend(got["rows"])
        mpra_rows.extend(got["mpra"])
        unassessed.extend(got["unassessed"])
        seen += got["blocks_seen"]
        print(
            f"{chrom}: {got['blocks_seen']} untouched blocks, {len(got['rows'])} assessed, "
            f"{sum(r['oligos'] for r in got['rows'])} oligos",
            flush=True,
        )

    p = meas.pooled(rows, min_usable=args.min_usable)
    by_case: dict[str, Any] = {}
    for label in sorted({r["label"] for r in rows}):
        sel = [r for r in rows if r["label"] == label]
        by_case[label] = {**meas.pooled(sel, min_usable=args.min_usable), "covariates": meas.covariates(sel)}

    usable_blocks = [r for r in rows if meas.block_counts(r)["attributable"] >= args.min_usable]
    thin_blocks = [r for r in rows if meas.block_counts(r)["attributable"] < args.min_usable]
    out = {
        "result": "measurability_real_unknown",
        "chromosomes": chroms,
        "definition": meas.__doc__.split("\n\n")[1] if meas.__doc__ else "",
        "family_a": list(meas.FAMILY_A),
        "family_b": list(meas.FAMILY_B),
        "arithmetic_reasons": list(meas.ARITHMETIC),
        "precedence": list(meas.PRECEDENCE),
        "peer_reading": PEER_READING,
        "reproduction_check": {
            "untouched_blocks_found": seen,
            "untouched_mb_found": round(sum(r["length"] for r in rows) / 1e6, 2),
            "expected": PEER_READING,
        },
        "pooled": p,
        "by_case": by_case,
        "window_covariates": meas.window_covariates(rows),
        "block_covariates": {
            "tileable": meas.covariates(usable_blocks),
            f"below_{args.min_usable}_usable_oligos": meas.covariates(thin_blocks),
        },
        "sensitivity": meas.sensitivity(rows),
        "calibration_lentimpra": meas.calibration(mpra_rows),
        "coverage": {
            "blocks_seen": seen,
            "blocks_assessed": len(rows),
            "blocks_unassessed": len(unassessed),
            "unassessed_by_missing_layer": {
                k: sum(1 for u in unassessed if u["missing"] == k)
                for k in sorted({u["missing"] for u in unassessed})
            },
            "blocks_without_tss_distance": sum(1 for r in rows if r["tss_distance"] is None),
        },
        "blocks": [{**{k: v for k, v in r.items() if k != "windows"}, **meas.block_counts(r)} for r in rows],
        "reading": (
            "a count of what the instrument cannot reach, not a claim about function. Family A cannot "
            "be synthesised or resolved; Family B can be ordered but its answer belongs to a sequence "
            "rather than to a locus, so a library may keep it flagged rather than drop it"
        ),
        "seconds": round(time.time() - t0, 1),
    }
    if args.no_save:
        print("not saved (--no-save)")
    else:
        print(f"saved {save_result(out['result'], out)}")

    print(
        f"\n{p['blocks']} blocks, {p['bp'] / 1e6:.2f} Mb, {p['oligos_tiled']} oligos at "
        f"{meas.OLIGO} bp step {meas.STEP} (peer's arithmetic: {PEER_READING['oligos']})"
    )
    print(f"{'reason':26} {'windows':>9} {'Mb':>8} {'share':>7}")
    for name in (*meas.FAMILY_A, *meas.FAMILY_B):
        n = p["charged_windows"][name]
        print(f"{name:26} {n:9} {n * meas.OLIGO / 1e6:8.3f} {n / max(1, p['oligos_tiled']):7.2%}")
    for name in meas.ARITHMETIC:
        bp = p["excluded_bp"][name]
        print(f"{name:26} {'-':>9} {bp / 1e6:8.3f} {bp / max(1, p['bp']):7.2%}")
    print(
        f"{'usable (both families)':26} {p['oligos_attributable']:9} "
        f"{p['attributable_bp'] / 1e6:8.3f} {p['oligos_attributable'] / max(1, p['oligos_tiled']):7.2%}"
    )
    print(
        f"\nsynthesisable oligos {p['oligos_synthesisable']}, attributable to one locus "
        f"{p['oligos_attributable']}, blocks under {args.min_usable} usable oligos "
        f"{p['blocks_below_min_usable']} of {p['blocks']}"
    )
    print(f"base accounting: {p['bp_accounted']} of {p['bp']} bp, unaccounted {p['bp_unaccounted']}")
    c = out["calibration_lentimpra"]
    print(
        f"lentiMPRA check: {c['excluded_family_a']} of {c['windows']} already-measured windows would be "
        f"excluded by Family A ({c['false_exclusion_rate']}), {c['flagged_family_b']} flagged by "
        f"Family B ({c['family_b_rate']}); GC {c['covariates']['median_gc']}, TSS "
        f"{(c['covariates']['median_tss_distance'] or 0) / 1000:.1f} kb"
    )
    print(f"\n{'group':30} {'n':>6} {'len':>8} {'GC':>6} {'TSS kb':>8}")
    for name, cov in (
        ("tileable blocks", out["block_covariates"]["tileable"]),
        (
            f"below {args.min_usable} usable",
            out["block_covariates"][f"below_{args.min_usable}_usable_oligos"],
        ),
    ):
        tss = cov["median_tss_distance"]
        print(
            f"{name:30} {cov['n']:6} {cov['median_length'] or 0:8.0f} {cov['median_gc'] or 0:6.3f} "
            f"{(tss / 1000 if tss else 0):8.1f}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
