"""Candidate regulatory elements from ENCODE (SCREEN registry of cCREs).

ENCODE's registry classifies ~1 million human elements from chromatin data
(DNase, H3K4me3, H3K27ac, CTCF ChIP): promoter-like (PLS), proximal and
distal enhancer-like (pELS, dELS), CTCF-only, and DNase-H3K4me3. This is
experimental evidence about the UNKNOWN space that no sequence pattern can
give. Following stream-distil-discard, the 64 MB genome-wide file is
streamed once and only the rows of the chromosomes asked for are kept.
"""

from __future__ import annotations

import gzip
import io
import urllib.request
from dataclasses import dataclass
from pathlib import Path

CCRE_URL = "https://downloads.wenglab.org/V3/GRCh38-cCREs.bed"
RESULTS = Path("data/results")
EVIDENCE = "ENCODE SCREEN registry of cCREs v3 (chromatin: DNase, H3K4me3, H3K27ac, CTCF)"

CLASS_LABEL = {
    "PLS": "promoter-like",
    "pELS": "proximal enhancer-like",
    "dELS": "distal enhancer-like",
    "CTCF-only": "CTCF site (insulator / loop anchor)",
    "DNase-H3K4me3": "open chromatin with promoter mark",
}


@dataclass(frozen=True, slots=True)
class CCRE:
    chrom: str
    start: int
    end: int
    id: str
    cls: str  # PLS, pELS, dELS, CTCF-only, DNase-H3K4me3
    ctcf_bound: bool

    @property
    def label(self) -> str:
        return CLASS_LABEL.get(self.cls, self.cls)


def _parse(line: str) -> CCRE | None:
    f = line.rstrip("\n").split("\t")
    if len(f) < 6:
        return None
    parts = f[5].split(",")
    return CCRE(f[0], int(f[1]), int(f[2]), f[4], parts[0], "CTCF-bound" in parts)


def stream_ccres(chroms: set[str], url: str = CCRE_URL, progress=None) -> list[CCRE]:
    """Stream the registry over HTTP and keep the rows of `chroms`; nothing is written to disk."""
    req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
    out: list[CCRE] = []
    n = 0
    with urllib.request.urlopen(req, timeout=300) as resp:  # noqa: S310
        for line in io.TextIOWrapper(io.BufferedReader(resp, 1 << 20), encoding="ascii", errors="replace"):
            n += 1
            if progress and n % 200_000 == 0:
                progress(n)
            if line.split("\t", 1)[0] in chroms:
                c = _parse(line)
                if c:
                    out.append(c)
    return out


def save_ccres(chrom: str, elements: list[CCRE]) -> Path:
    RESULTS.mkdir(parents=True, exist_ok=True)
    p = RESULTS / f"ccres_{chrom}.bed.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(f"# ENCODE cCREs v3, {chrom} subset distilled by GenomeOS from {CCRE_URL}\n")
        for c in elements:
            fh.write(f"{c.chrom}\t{c.start}\t{c.end}\t{c.id}\t{c.cls}\t{int(c.ctcf_bound)}\n")
    return p


def load_ccres(chrom: str) -> list[CCRE]:
    p = RESULTS / f"ccres_{chrom}.bed.gz"
    if not p.exists():
        return []
    out = []
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            out.append(CCRE(f[0], int(f[1]), int(f[2]), f[3], f[4], f[5] == "1"))
    return out


def summarise(elements: list[CCRE]) -> dict:
    by: dict[str, int] = {}
    bp: dict[str, int] = {}
    for c in elements:
        by[c.cls] = by.get(c.cls, 0) + 1
        bp[c.cls] = bp.get(c.cls, 0) + (c.end - c.start)
    return {
        "elements": len(elements),
        "by_class": by,
        "bp_by_class": bp,
        "ctcf_bound": sum(1 for c in elements if c.ctcf_bound),
        "evidence": EVIDENCE,
    }


def ccre_index(elements: list[CCRE]) -> list[tuple[int, int, str]]:
    """Sorted (start, end, cls) for fast density queries."""
    return sorted((c.start, c.end, c.cls) for c in elements)


def count_in(index: list[tuple[int, int, str]], start: int, end: int) -> dict[str, int]:
    """cCREs overlapping [start, end) by class (binary search on start)."""
    import bisect

    starts = [x[0] for x in index]
    i = bisect.bisect_left(starts, start - 5000)
    out: dict[str, int] = {}
    while i < len(index) and index[i][0] < end:
        s, e, cls = index[i]
        if e > start:
            out[cls] = out.get(cls, 0) + 1
        i += 1
    return out
