# SPDX-License-Identifier: AGPL-3.0-or-later
"""Segmental duplications: the genome's copy-and-paste, read from UCSC's curated track.

`genomicSuperDups` (Bailey and Eichler's whole-genome assembly comparison, kept current on
hg38) lists every pair of regions over 1 kb of non-repeat-masked sequence that align at
90% identity or more, with the partner locus and the fraction of matching bases. Fetched
per chromosome through the UCSC track API, as RepeatMasker is, and kept as a compact result.

What it is for (ROADMAP area J step 3, ATTRIBUTION step 2's "organising evidence"):

- every UNKNOWN block of the budget gets its duplicated fraction and its partners, so a
  block that is a copy reads as a copy before anything is attributed to it;
- the UNKNOWN classifier's `similar_to` pairs (shared 20-mers between the largest blocks,
  a heuristic) are held against the curated pairs, which says how far to trust it;
- the scored regulatory elements that lie inside a duplication are counted: a duplicated
  enhancer is regulatory code pasted twice, and the deletion model scored one copy.

Identity dates a copy: 0.98 and above is a young duplication (primate), below 0.95 older;
the summary keeps the three bands. Evidence is `curated` for the pairs and `inferred` for
what a block's duplication implies.
"""

from __future__ import annotations

import json
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, load_result, save_result

UCSC = "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=genomicSuperDups;chrom={chrom};maxItemsOutput=1000000"
YOUNG = 0.98  # fraction of matching bases from which a copy reads as young (primate-era)
OLD = 0.95
EVIDENCE = {
    "pairs": "curated: UCSC hg38 genomicSuperDups (segmental duplications >= 1 kb, >= 90% identity)",
    "blocks": (
        "inferred: a block's duplicated fraction from the curated pairs; what the copy does is not said"
    ),
}


@dataclass(frozen=True)
class SegDup:
    start: int
    end: int
    other_chrom: str
    other_start: int
    other_end: int
    frac_match: float
    strand: str

    def as_row(self) -> list:
        return [
            self.start,
            self.end,
            self.other_chrom,
            self.other_start,
            self.other_end,
            round(self.frac_match, 4),
            self.strand,
        ]

    @classmethod
    def from_row(cls, r: list) -> SegDup:
        return cls(r[0], r[1], r[2], r[3], r[4], r[5], r[6])


def rows_path(chrom: str, results_dir: Path = RESULTS_DIR) -> Path:
    return results_dir / f"superdups_{chrom}.bed.gz"


def save_rows(chrom: str, dups: list[SegDup], results_dir: Path = RESULTS_DIR) -> Path:
    """The curated pairs as a local compressed BED (raw track rows stay out of the repository)."""
    import gzip

    p = rows_path(chrom, results_dir)
    p.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(p, "wt") as fh:
        for d in dups:
            fh.write("\t".join(str(x) for x in d.as_row()) + "\n")
    return p


def load_rows(chrom: str, results_dir: Path = RESULTS_DIR) -> list[SegDup]:
    import gzip

    p = rows_path(chrom, results_dir)
    if not p.exists():
        return []
    out = []
    with gzip.open(p, "rt") as fh:
        for line in fh:
            s, e, oc, os_, oe, fm, st = line.rstrip("\n").split("\t")
            out.append(SegDup(int(s), int(e), oc, int(os_), int(oe), float(fm), st))
    return out


def fetch(chrom: str, timeout: int = 300) -> list[SegDup]:
    req = urllib.request.Request(
        UCSC.format(chrom=chrom), headers={"User-Agent": "GenomeOS/0.9 (duplications)"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        d = json.load(r)
    rows = d.get("genomicSuperDups", [])
    if isinstance(rows, dict):
        rows = rows.get(chrom, [])
    out = [
        SegDup(
            x["chromStart"],
            x["chromEnd"],
            x["otherChrom"],
            x["otherStart"],
            x["otherEnd"],
            x["fracMatch"],
            x["strand"],
        )
        for x in rows
    ]
    out.sort(key=lambda s: (s.start, s.end))
    return out


def merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], e))
        else:
            out.append((s, e))
    return out


def covered(merged: list[tuple[int, int]], start: int, end: int) -> int:
    """Bases of [start, end) inside the merged intervals."""
    total = 0
    for s, e in merged:
        if e <= start:
            continue
        if s >= end:
            break
        total += min(e, end) - max(s, start)
    return total


def band(frac: float) -> str:
    return "young" if frac >= YOUNG else "middle" if frac >= OLD else "old"


def parse_similar(s: str) -> tuple[str, int, int] | None:
    """'chr21:100-200 (jaccard 0.31)' → ('chr21', 100, 200)."""
    try:
        loc = s.split()[0]
        chrom, span = loc.split(":")
        a, b = span.split("-")
        return chrom, int(a), int(b)
    except (ValueError, IndexError):
        return None


def pair_supported(dups: list[SegDup], chrom: str, a: tuple[int, int], b: tuple[int, int]) -> bool:
    """Is there a curated pair whose one span overlaps block a and whose other span overlaps block b?"""
    for d in dups:
        here = d.start < a[1] and a[0] < d.end
        there = d.other_chrom == chrom and d.other_start < b[1] and b[0] < d.other_end
        if here and there:
            return True
        here = d.start < b[1] and b[0] < d.end
        there = d.other_chrom == chrom and d.other_start < a[1] and a[0] < d.other_end
        if here and there:
            return True
    return False


def summarise(
    chrom: str,
    dups: list[SegDup],
    chrom_length: int | None,
    budget_blocks: list[dict] | None,
    unknown_blocks: list[dict] | None,
    elements: list[dict] | None,
) -> dict[str, Any]:
    merged = merge([(d.start, d.end) for d in dups])
    dup_bp = sum(e - s for s, e in merged)
    bands = {"young": 0, "middle": 0, "old": 0}
    intra = 0
    for d in dups:
        bands[band(d.frac_match)] += 1
        intra += d.other_chrom == chrom
    out: dict[str, Any] = {
        "pairs": len(dups),
        "intra_chromosomal": intra,
        "inter_chromosomal": len(dups) - intra,
        "bands": bands,
        "duplicated_bp": dup_bp,
        "duplicated_fraction": round(dup_bp / chrom_length, 4) if chrom_length else None,
    }
    # per UNKNOWN block, and per tier
    rows = []
    by_tier: dict[str, dict] = {}
    if budget_blocks:
        for b in budget_blocks:
            c = covered(merged, b["start"], b["end"])
            partners = sorted(
                {
                    f"{d.other_chrom}:{d.other_start}-{d.other_end}"
                    for d in dups
                    if d.start < b["end"] and b["start"] < d.end
                }
            )
            frac = round(c / b["length"], 4) if b["length"] else 0.0
            tier = b["guess"]["tier"]
            t = by_tier.setdefault(
                tier, {"blocks": 0, "bp": 0, "duplicated_bp": 0, "blocks_mostly_duplicated": 0}
            )
            t["blocks"] += 1
            t["bp"] += b["length"]
            t["duplicated_bp"] += c
            t["blocks_mostly_duplicated"] += frac >= 0.5
            rows.append(
                {
                    "start": b["start"],
                    "end": b["end"],
                    "class": b["class"],
                    "tier": tier,
                    "duplicated_fraction": frac,
                    "pairs": len(partners),
                    "partners": partners[:8],
                }
            )
        for t in by_tier.values():
            t["duplicated_fraction"] = round(t["duplicated_bp"] / t["bp"], 4) if t["bp"] else None
    out["blocks"] = rows
    out["by_tier"] = by_tier
    # the heuristic against the curated pairs
    if unknown_blocks:
        pairs = set()
        for b in unknown_blocks:
            for s in b.get("similar_to") or []:
                p = parse_similar(s)
                if p and p[0] == chrom:
                    a, c = (b["start"], b["end"]), (p[1], p[2])
                    pairs.add((min(a, c), max(a, c)))
        supported = sum(pair_supported(dups, chrom, a, c) for a, c in pairs)
        out["similar_to_pairs"] = {
            "pairs": len(pairs),
            "supported_by_curated_duplication": supported,
            "share": round(supported / len(pairs), 3) if pairs else None,
        }
    if elements:
        inside = [e for e in elements if covered(merged, e["start"], e["end"]) > 0]
        out["elements"] = {
            "scored": len(elements),
            "inside_a_duplication": len(inside),
            "share": round(len(inside) / len(elements), 3) if elements else None,
            "ids": [e["id"] for e in inside[:50]],
        }
    return out


def elements_of(chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict]:
    seen: set[str] = set()
    out = []
    for name in ("constrained_targets", "enhancer_targets"):
        r = load_result(f"{name}_{chrom}", results_dir) or {}
        for e in r.get("elements", []):
            if e["id"] not in seen:
                seen.add(e["id"])
                out.append({"id": e["id"], "start": e["start"], "end": e["end"]})
    return out


def run_and_save(chrom: str, results_dir: Path = RESULTS_DIR, dups: list[SegDup] | None = None) -> dict:
    t0 = time.time()
    if dups is None:
        dups = load_rows(chrom, results_dir) or fetch(chrom)
    save_rows(chrom, dups, results_dir)
    budget = load_result(f"budget_{chrom}", results_dir) or {}
    unknown = load_result(f"unknown_{chrom}", results_dir) or {}
    out = summarise(
        chrom,
        dups,
        budget.get("chromosome_length"),
        budget.get("blocks"),
        unknown.get("blocks"),
        elements_of(chrom, results_dir),
    )
    out.update(
        {
            "chrom": chrom,
            "rows_file": str(rows_path(chrom, results_dir)),
            "row_fields": ["start", "end", "other_chrom", "other_start", "other_end", "frac_match", "strand"],
            "bands_definition": {"young": f">= {YOUNG}", "middle": f"{OLD} to {YOUNG}", "old": f"< {OLD}"},
            "evidence": EVIDENCE,
            "cost": {"seconds": round(time.time() - t0, 1)},
        }
    )
    save_result(f"duplication_{chrom}", out, results_dir)
    return out


def distil(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    """The genome-wide summary of every saved duplication_<chrom> result."""
    import glob

    per_chrom: dict[str, dict] = {}
    by_tier: dict[str, dict] = {}
    bands = {"young": 0, "middle": 0, "old": 0}
    pairs = intra = dup_bp = 0
    genome_bp = 0
    similar = {"pairs": 0, "supported": 0}
    elements = {"scored": 0, "inside": 0}
    for f in sorted(glob.glob(str(results_dir / "duplication_chr*.json"))):
        r = load_result(Path(f).stem, results_dir)
        if not r or "by_tier" not in r:
            continue
        chrom = r["chrom"]
        pairs += r["pairs"]
        intra += r["intra_chromosomal"]
        dup_bp += r["duplicated_bp"]
        if r.get("duplicated_fraction"):
            genome_bp += round(r["duplicated_bp"] / r["duplicated_fraction"])
        for k in bands:
            bands[k] += r["bands"][k]
        for t, v in r["by_tier"].items():
            d = by_tier.setdefault(
                t, {"blocks": 0, "bp": 0, "duplicated_bp": 0, "blocks_mostly_duplicated": 0}
            )
            for k in d:
                d[k] += v[k]
        sp = r.get("similar_to_pairs") or {}
        similar["pairs"] += sp.get("pairs", 0)
        similar["supported"] += sp.get("supported_by_curated_duplication", 0)
        el = r.get("elements") or {}
        elements["scored"] += el.get("scored", 0)
        elements["inside"] += el.get("inside_a_duplication", 0)
        per_chrom[chrom] = {
            "pairs": r["pairs"],
            "duplicated_fraction": r.get("duplicated_fraction"),
            "constrained_unknown_duplicated": (r["by_tier"].get("constrained_unknown") or {}).get(
                "duplicated_fraction"
            ),
        }
    for d in by_tier.values():
        d["duplicated_fraction"] = round(d["duplicated_bp"] / d["bp"], 4) if d["bp"] else None
        d["share_blocks_mostly_duplicated"] = (
            round(d["blocks_mostly_duplicated"] / d["blocks"], 3) if d["blocks"] else None
        )
    similar["share"] = round(similar["supported"] / similar["pairs"], 3) if similar["pairs"] else None
    elements["share"] = round(elements["inside"] / elements["scored"], 3) if elements["scored"] else None
    return {
        "chromosomes": len(per_chrom),
        "pairs": pairs,
        "intra_chromosomal": intra,
        "bands": bands,
        "duplicated_bp": dup_bp,
        "duplicated_fraction_of_genome": round(dup_bp / genome_bp, 4) if genome_bp else None,
        "by_tier": by_tier,
        "similar_to_pairs": similar,
        "elements": elements,
        "per_chromosome": per_chrom,
        "evidence": EVIDENCE,
    }
