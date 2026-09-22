# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured mappability for every oligo of the library's test arm, against the repeat proxy.

    uv run python scripts/mappability.py [--chroms chr21,chr22] [--ks 24,36] [--no-save]

The library's test arm (`data/knowledge/library/unknown_library.tsv`, 91,919 oligos of 300 bp over the
real unknown's blocks) is already the set that survived the synthesis rules. What it does not have
measured is whether each oligo's answer could be attributed to one locus: that figure is still a
proxy over interspersed-repeat and segmental-duplication content, and it swings 19,084 to 28,563
oligos as the repeat cut-off moves.

This script reads Umap multi-read mappability for GRCh38 at k = 24, 36, 50 and 100, one k per column
and never pooled, through the project's bigWig range reader -- no file is downloaded, and every run
prints the bytes it moved. The two arms are reported separately and cross-tabulated; the cell the
design cares about is the oligo the proxy passes and the track calls unmappable.

`genomeos.attribution.mappability` holds the tracks, the rule and its three values. This script only
joins the local layers: the manifest, RepeatMasker, the curated segmental duplications and the block
list of the 680 untouched blocks that `measurability_real_unknown` assessed.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path
from typing import Any

from genomeos.attribution import mappability as mp
from genomeos.attribution import measurability as meas
from genomeos.results import RESULTS_DIR, save_result

MANIFEST = Path("data/knowledge/library/unknown_library.tsv")
ARM = "test"
CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]


def block_span(block: str) -> tuple[str, int, int]:
    chrom, _, rest = block.partition(":")
    lo, _, hi = rest.partition("-")
    return chrom, int(lo), int(hi)


def untouched_blocks(results_dir: Path = RESULTS_DIR) -> set[str]:
    """The 680 blocks `measurability_real_unknown` assessed: the published proxy's own universe."""
    p = results_dir / "measurability_real_unknown.json"
    if not p.exists():
        return set()
    return {b["block"] for b in json.loads(p.read_text()).get("blocks", [])}


def load_manifest(path: Path = MANIFEST, arm: str = ARM) -> tuple[dict[str, list[dict[str, Any]]], dict]:
    """Every oligo of one arm by chromosome, plus the checks the reading of missing values rests on.

    A missing Umap value is read as a measured zero rather than as an unassessed base, which is only
    sound if no oligo can contain an assembly gap. The N count and the length check are therefore
    part of the output and not an assumption.
    """
    by_chrom: dict[str, list[dict[str, Any]]] = {}
    n_bases = wrong_length = 0
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = {name: i for i, name in enumerate(header)}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if parts[cols["arm"]] != arm:
                continue
            chrom = parts[cols["chrom"]]
            start, end = int(parts[cols["start"]]), int(parts[cols["end"]])
            seq = parts[cols["sequence"]].upper()
            n_bases += len(seq) - sum(seq.count(b) for b in "ACGT")
            if end - start != mp.OLIGO or len(seq) != mp.OLIGO:
                wrong_length += 1
            block = parts[cols["block"]]
            tss = parts[cols["tss_distance"]]
            by_chrom.setdefault(chrom, []).append(
                {
                    "id": parts[cols["id"]],
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "block": block,
                    "block_length": block_span(block)[2] - block_span(block)[1],
                    "gc": float(parts[cols["gc"]]) if parts[cols["gc"]] else None,
                    "tss_distance": int(tss) if tss else None,
                }
            )
    checks = {
        "arm": arm,
        "oligos": sum(len(v) for v in by_chrom.values()),
        "n_bases_in_arm": n_bases,
        "oligos_not_one_oligo_long": wrong_length,
        "reading": (
            "with no N base in the arm, a base the track carries no value for is unmappable and not "
            "an assembly gap; that is what licenses counting missing values as zeros"
        ),
    }
    return by_chrom, checks


def add_proxy(chrom: str, oligos: list[dict[str, Any]]) -> str | None:
    """Family B's two fractions per oligo, from the local tracks; the missing layer's name on failure."""
    blocks = sorted({block_span(o["block"])[1:] for o in oligos})
    by_class = meas.repeat_spans_by_class(chrom, list(blocks))
    if by_class is None:
        return "rmsk"
    dups = meas.segdup_spans(chrom)
    if dups is None:
        return "superdups"
    interspersed = meas.class_groups(by_class).get("interspersed", [])
    segdups = dups[0]
    for o in oligos:
        o["interspersed_fraction"] = round(meas.fraction_covered(interspersed, o["start"], o["end"]), 4)
        o["segdup_fraction"] = round(meas.fraction_covered(segdups, o["start"], o["end"]), 4)
    return None


def per_track(
    k: int,
    by_chrom: dict[str, list[dict[str, Any]]],
    args: argparse.Namespace,
) -> dict[str, Any]:
    """One k: read the track chromosome by chromosome, then count. Nothing crosses k."""
    track = mp.TRACKS[k]
    rows: list[dict[str, Any]] = []
    unassessed: dict[str, int] = {}
    stats = {"bytes_fetched": 0, "range_requests": 0, "seconds": 0.0, "oligos_from_cache": 0}
    per_chrom: dict[str, Any] = {}
    for chrom in sorted(by_chrom, key=lambda c: CHROMS.index(c) if c in CHROMS else 99):
        oligos = by_chrom[chrom]
        res = mp.read_chromosome(
            k,
            chrom,
            [(o["start"], o["end"]) for o in oligos],
            byte_cap=args.byte_cap,
            refresh=args.refresh,
        )
        stats["bytes_fetched"] += res.bytes_fetched
        stats["range_requests"] += res.requests
        stats["seconds"] += res.seconds
        stats["oligos_from_cache"] += res.from_cache
        per_chrom[chrom] = res.as_dict()
        if not res.in_track:
            unassessed["chromosome_not_in_track"] = unassessed.get("chromosome_not_in_track", 0) + len(oligos)
            continue
        if res.cap_reached:
            unassessed["byte_cap_reached"] = unassessed.get("byte_cap_reached", 0) + res.unread
        for o in oligos:
            if "interspersed_fraction" not in o:
                continue  # already counted under the missing repeat layer
            m = res.rows.get((o["start"], o["end"]))
            if m is None:
                unassessed["no_track_row"] = unassessed.get("no_track_row", 0) + 1
                continue
            rows.append({**o, "bases": m["bases"], "above": m["above"], "total": m["total"]})
        print(
            f"  k{k} {chrom}: {len(oligos)} oligos, {res.bytes_fetched / 1e6:.1f} MB, "
            f"{res.seconds:.0f} s, {res.from_cache} cached",
            flush=True,
        )
    untouched = untouched_blocks()
    sub = [r for r in rows if r["block"] in untouched]
    seen_blocks = {r["block"] for r in sub}
    zero = sum(1 for r in rows if r["bases"] == 0)
    return {
        **track.as_dict(),
        "read": {**stats, "seconds": round(stats["seconds"], 1), "per_chromosome": per_chrom},
        "base_cut": mp.BASE_CUT,
        "min_fraction": args.min_fraction,
        "arms": mp.arms(rows, args.min_fraction),
        "crosstab": mp.crosstab(rows, args.min_fraction),
        "sensitivity": mp.sensitivity(rows),
        "proxy_sensitivity": mp.proxy_sensitivity(rows, min_fraction=args.min_fraction),
        "groups": mp.groups(rows, args.min_fraction),
        "matched_proxy_against_track": mp.proxy_against_track(rows, args.min_fraction),
        "untouched_680_blocks": {
            "oligos": len(sub),
            "blocks": len(seen_blocks),
            "blocks_expected": len(untouched),
            "blocks_without_a_test_oligo": sorted(untouched - seen_blocks),
            "blocks_without_a_test_oligo_reason": (
                "the library's test arm stops at chrX, so a chrY block of the 680 has no oligo to "
                "assess; it is named here rather than counted as mappable or unmappable"
            ),
            "arms": mp.arms(sub, args.min_fraction),
            "crosstab": mp.crosstab(sub, args.min_fraction),
            "sensitivity": mp.sensitivity(sub),
            "groups": mp.groups(sub, args.min_fraction),
        },
        "coverage": {
            "oligos_assessed": len(rows),
            "oligos_unassessed": sum(unassessed.values()),
            "unassessed_by_category": unassessed,
            "oligos_in_the_680_untouched_blocks": len(sub),
            "oligos_in_previously_assayed_blocks": len(rows) - len(sub),
            "oligos_with_no_value_anywhere": zero,
            "oligos_with_no_value_reading": (
                "counted as unmappable, not unassessed: Umap writes nothing where mappability is "
                "zero and the arm has no N base"
            ),
        },
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--ks", default=",".join(str(k) for k in mp.DEFAULT_KS))
    ap.add_argument("--arm", default=ARM)
    ap.add_argument("--min-fraction", type=float, default=mp.MAPPABLE_FRACTION)
    ap.add_argument("--byte-cap", type=int, default=mp.BYTE_CAP)
    ap.add_argument("--refresh", action="store_true", help="ignore the cached per-oligo summaries")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    ks = [int(x) for x in args.ks.split(",") if x]

    t0 = time.time()
    if not MANIFEST.exists():
        print(f"no manifest at {MANIFEST}: nothing to assess")
        return 1
    by_chrom, checks = load_manifest(MANIFEST, args.arm)
    if args.chroms:
        wanted = set(args.chroms.split(","))
        by_chrom = {c: v for c, v in by_chrom.items() if c in wanted}
    print(f"{checks['oligos']} {args.arm} oligos, {len(by_chrom)} chromosomes; {checks}", flush=True)

    missing_layer: dict[str, int] = {}
    for chrom, oligos in by_chrom.items():
        missing = add_proxy(chrom, oligos)
        if missing:
            missing_layer[missing] = missing_layer.get(missing, 0) + len(oligos)
            print(f"{chrom}: proxy not computed ({missing} not cached)", flush=True)

    tracks = {}
    for k in ks:
        print(f"k{k}: {mp.TRACKS[k].url}", flush=True)
        tracks[f"k{k}"] = per_track(k, by_chrom, args)

    out: dict[str, Any] = {
        "result": "mappability_real_unknown",
        "question": (
            "the library's 'attributable to one locus' figure rested on a repeat proxy; this measures "
            "mappability itself and reports the two arms separately"
        ),
        "assembly": mp.ASSEMBLY,
        "track_kind": mp.TRACK_KIND,
        "citation": mp.CITATION,
        "manifest": str(MANIFEST),
        "manifest_checks": checks,
        "chromosomes": sorted(by_chrom),
        "chromosomes_absent_from_the_arm": [c for c in CHROMS if c not in by_chrom],
        "oligo_bp": mp.OLIGO,
        "rule": (
            f"a base is uniquely readable at multi-read mappability >= BASE_CUT {mp.BASE_CUT}; an oligo "
            f"is mappable when at least MAPPABLE_FRACTION of its {mp.OLIGO} bases are, reported at "
            f"{list(mp.SENSITIVITY_FRACTIONS)}"
        ),
        "proxy_definition": (
            f"measurability Family B: interspersed_fraction >= {meas.INTERSPERSED_MAX} or "
            f"segdup_fraction >= {meas.SEGDUP_MAX}"
        ),
        "byte_cap_per_track": args.byte_cap,
        "cache_dir": str(mp.CACHE_DIR),
        "cache_bytes": mp.cache_bytes(),
        "cache_cap": mp.CACHE_CAP,
        "tracks": tracks,
        "proxy_layer_missing": missing_layer,
        "not_merged": (
            "the mappability arm is never added to or averaged with the repeat proxy, and no number "
            "here pools two values of k"
        ),
        "seconds": round(time.time() - t0, 1),
    }
    if args.no_save:
        print("not saved (--no-save)")
    else:
        print(f"saved {save_result(out['result'], out)}")

    total_bytes = sum(t["read"]["bytes_fetched"] for t in tracks.values())
    print(
        f"\nbytes read {total_bytes / 1e6:.1f} MB over {len(ks)} tracks, "
        f"cache {mp.cache_bytes() / 1e6:.1f} MB"
    )
    print(f"\n{'k':>4} {'oligos':>8} {'proxy':>8} {'track':>8} {'both':>8} {'agree':>7}")
    for k in ks:
        t = tracks[f"k{k}"]
        a, x = t["arms"], t["crosstab"]
        print(
            f"{k:>4} {a['oligos_assessed']:8} {a['attributable_by_repeat_proxy']:8} "
            f"{a['attributable_by_mappability_track']:8} {a['attributable_by_both_arms']:8} "
            f"{x['agreement'] or 0:7.3f}"
        )
    print(f"\ncross-tabulation at min_fraction {args.min_fraction}")
    print(f"{'k':>4} " + " ".join(f"{c:>32}" for c in mp.CELLS))
    for k in ks:
        x = tracks[f"k{k}"]["crosstab"]
        print(f"{k:>4} " + " ".join(f"{x[c]:>32}" for c in mp.CELLS))
    print(
        f"\nthe 680 untouched blocks alone, the published proxy's own universe "
        f"(min_fraction {args.min_fraction})"
    )
    print(f"{'k':>4} {'oligos':>8} {'proxy':>8} {'track':>8} {'both':>8} {'dangerous':>10}")
    for k in ks:
        u = tracks[f"k{k}"]["untouched_680_blocks"]
        a = u["arms"]
        print(
            f"{k:>4} {u['oligos']:8} {a['attributable_by_repeat_proxy']:8} "
            f"{a['attributable_by_mappability_track']:8} {a['attributable_by_both_arms']:8} "
            f"{u['crosstab']['proxy_passes_track_unmappable']:10}"
        )
    print("\nsensitivity: the rule at three values")
    for k in ks:
        for s in tracks[f"k{k}"]["sensitivity"]:
            print(
                f"k{k:<4} min_fraction {s['min_fraction']:<5} mappable "
                f"{s['attributable_by_mappability_track']:8} dangerous cell "
                f"{s['crosstab']['proxy_passes_track_unmappable']:8}"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
