"""Curated repeat annotation (RepeatMasker via the UCSC track API).

Sequence patterns find young repeats; RepeatMasker's library finds all of
them, including copies 30% diverged from their consensus, and names their
class (LINE, SINE, LTR, DNA, Satellite, Simple_repeat, Low_complexity) and
family. GenomeOS fetches one chromosome's rows from the public UCSC API,
keeps a compact BED (start, end, class, family, name, divergence) and a
class summary, and uses the intervals to give UNKNOWN blocks a curated
repeat class. Stream, distil, discard: chr21 is 20 MB of JSON once, 1 MB
kept.
"""

from __future__ import annotations

import bisect
import gzip
import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

UCSC = "https://api.genome.ucsc.edu/getData/track?genome=hg38;track=rmsk;chrom={chrom};maxItemsOutput=1000000"
RESULTS = Path("data/results")
EVIDENCE = "RepeatMasker annotation of hg38 (UCSC rmsk track, Dfam/Repbase library)"
INTERSPERSED = ("LINE", "SINE", "LTR", "DNA", "Retroposon", "RC")


@dataclass(frozen=True, slots=True)
class Repeat:
    start: int
    end: int
    cls: str
    family: str
    name: str
    divergence: float  # fraction of bases diverged from the consensus

    @property
    def length(self) -> int:
        return self.end - self.start


def fetch_repeats(chrom: str, timeout: int = 300) -> list[Repeat]:
    req = urllib.request.Request(UCSC.format(chrom=chrom), headers={"User-Agent": "GenomeOS/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        d = json.load(r)
    rows = d.get("rmsk", [])
    if isinstance(rows, dict):  # some tracks key rows by chromosome
        rows = rows.get(chrom, [])
    out = [
        Repeat(
            x["genoStart"],
            x["genoEnd"],
            x["repClass"].rstrip("?"),
            x["repFamily"].rstrip("?"),
            x["repName"],
            x.get("milliDiv", 0) / 1000,
        )
        for x in rows
    ]
    out.sort(key=lambda r: r.start)
    return out


def save_repeats(chrom: str, repeats: list[Repeat]) -> Path:
    RESULTS.mkdir(parents=True, exist_ok=True)
    p = RESULTS / f"rmsk_{chrom}.bed.gz"
    with gzip.open(p, "wt") as fh:
        for r in repeats:
            fh.write(f"{r.start}\t{r.end}\t{r.cls}\t{r.family}\t{r.name}\t{r.divergence:.3f}\n")
    return p


def load_repeats(chrom: str) -> list[Repeat]:
    p = RESULTS / f"rmsk_{chrom}.bed.gz"
    if not p.exists():
        return []
    out = []
    with gzip.open(p, "rt") as fh:
        for line in fh:
            s, e, c, f, n, d = line.rstrip("\n").split("\t")
            out.append(Repeat(int(s), int(e), c, f, n, float(d)))
    return out


def summarise(chrom: str, repeats: list[Repeat], chrom_length: int | None = None) -> dict[str, Any]:
    by_class: dict[str, dict[str, Any]] = {}
    by_family: dict[str, int] = {}
    for r in repeats:
        d = by_class.setdefault(r.cls, {"copies": 0, "bp": 0, "divergence_sum": 0.0})
        d["copies"] += 1
        d["bp"] += r.length
        d["divergence_sum"] += r.divergence
        by_family[f"{r.cls}/{r.family}"] = by_family.get(f"{r.cls}/{r.family}", 0) + r.length
    for d in by_class.values():
        d["mean_divergence"] = round(d["divergence_sum"] / d["copies"], 3)
        del d["divergence_sum"]
    total = sum(r.length for r in repeats)
    return {
        "chrom": chrom,
        "copies": len(repeats),
        "repeat_bp": total,
        "repeat_fraction": round(total / chrom_length, 4) if chrom_length else None,
        "by_class": dict(sorted(by_class.items(), key=lambda kv: -kv[1]["bp"])),
        "top_families": dict(sorted(by_family.items(), key=lambda kv: -kv[1])[:15]),
        "evidence": EVIDENCE,
        "confidence": 0.9,
    }


class RepeatIndex:
    """Per-class base coverage of any interval, from sorted repeat intervals."""

    def __init__(self, repeats: list[Repeat]) -> None:
        self.repeats = sorted(repeats, key=lambda r: r.start)
        self.starts = [r.start for r in self.repeats]
        self.max_len = max((r.length for r in self.repeats), default=0)

    def __bool__(self) -> bool:
        return bool(self.repeats)

    def coverage(self, start: int, end: int) -> dict[str, int]:
        """Bases of [start, end) covered by each repeat class."""
        out: dict[str, int] = {}
        i = bisect.bisect_left(self.starts, start - self.max_len)
        for r in self.repeats[i:]:
            if r.start >= end:
                break
            ov = min(r.end, end) - max(r.start, start)
            if ov > 0:
                out[r.cls] = out.get(r.cls, 0) + ov
        return out


def repeat_index(chrom: str) -> RepeatIndex:
    return RepeatIndex(load_repeats(chrom))
