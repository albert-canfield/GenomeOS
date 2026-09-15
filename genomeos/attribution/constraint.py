# SPDX-License-Identifier: AGPL-3.0-or-later
"""Evolutionary constraint over UNKNOWN blocks, from two public sources.

1. Zoonomia phyloP over 241 placental mammals (Christmas et al. 2023, Science), read
   per base from UCSC's bigWig with range requests, never downloaded. A base with
   phyloP >= 2.27 is "constrained" at 5% FDR, the paper's own threshold; 10.7% of
   the genome passes it. Per block we keep the bases seen, the mean, the maximum and
   the constrained fraction.
2. The 100-vertebrate phastCons conserved elements (UCSC `phastConsElements100way`),
   fetched per chromosome through the UCSC REST API and cached as a small BED under
   `data/knowledge/constraint` (a cache, git-ignored). Per block: how many elements,
   how many bases they cover, and the highest lod score.

Both are evidence of *purifying selection*, which is the nearest thing to a
measurement of "this matters for the organism" that exists for every base. Neither
says what the base does; that is what the rest of the package is for.
"""

from __future__ import annotations

import gzip
import json
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from genomeos.attribution.bigwig import BigWig, IntervalStats

PHYLOP_241_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/cactus241way/cactus241way.phyloP.bw"
PHYLOP_THRESHOLD = 2.27  # Zoonomia: FDR 5% for constraint
UCSC_API = "https://api.genome.ucsc.edu/getData/track"
ELEMENTS_TRACK = "phastConsElements100way"
CACHE = Path("data/knowledge/constraint")
WINDOW = 5_000_000  # one API request per window of conserved elements
PAUSE = 1.0  # seconds between API requests; UCSC asks for restraint


@dataclass(frozen=True)
class Element:
    start: int
    end: int
    lod: int


def phylop_over_blocks(
    chrom: str,
    intervals: list[tuple[int, int]],
    threshold: float = PHYLOP_THRESHOLD,
    url: str = PHYLOP_241_URL,
    progress=None,
) -> tuple[list[IntervalStats], dict]:
    """phyloP summaries for the intervals of one chromosome; returns the stats and what it cost."""
    bw = BigWig(url)
    try:
        t0 = time.time()
        stats = bw.summarise(chrom, intervals, threshold, progress=progress)
        cost = {
            "requests": bw.src.requests,
            "mb_fetched": round(bw.src.bytes_fetched / 1e6, 1),
            "seconds": round(time.time() - t0, 1),
        }
    finally:
        bw.close()
    return stats, cost


def _api(params: dict, timeout: int = 180) -> dict:
    q = ";".join(f"{k}={v}" for k, v in params.items())
    req = urllib.request.Request(
        f"{UCSC_API}?{q}", headers={"User-Agent": "GenomeOS/0.9 (conserved elements)"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        return json.load(r)


def fetch_elements(
    chrom: str,
    chrom_length: int,
    track: str = ELEMENTS_TRACK,
    cache: Path = CACHE,
    progress=None,
) -> list[Element]:
    """Conserved elements of one chromosome, from the cache or the UCSC API in 5 Mb windows."""
    cache.mkdir(parents=True, exist_ok=True)
    p = cache / f"{track}_{chrom}.bed.gz"
    if p.exists():
        return load_elements(p)
    out: list[Element] = []
    seen: set[tuple[int, int]] = set()
    for start in range(0, chrom_length, WINDOW):
        end = min(start + WINDOW, chrom_length)
        d = _api(
            {
                "genome": "hg38",
                "track": track,
                "chrom": chrom,
                "start": start,
                "end": end,
                "maxItemsOutput": 1_000_000,
            }
        )
        for it in d.get(track, []):
            key = (it["chromStart"], it["chromEnd"])
            if key in seen:
                continue
            seen.add(key)
            lod = int(str(it.get("name", "lod=0")).split("=")[-1] or 0)
            out.append(Element(it["chromStart"], it["chromEnd"], lod))
        if progress:
            progress(end, chrom_length, len(out))
        time.sleep(PAUSE)
    out.sort(key=lambda e: e.start)
    with gzip.open(p, "wt") as fh:
        fh.write(f"# UCSC hg38 {track}, {chrom}, fetched by GenomeOS through {UCSC_API}\n")
        for e in out:
            fh.write(f"{chrom}\t{e.start}\t{e.end}\t{e.lod}\n")
    return out


def load_elements(p: Path) -> list[Element]:
    out = []
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            out.append(Element(int(f[1]), int(f[2]), int(f[3])))
    return out


def elements_over_blocks(elements: list[Element], intervals: list[tuple[int, int]]) -> list[dict]:
    """Per interval: conserved elements overlapping it, bases covered (clipped), highest lod.

    Elements and intervals are both sorted by start; intervals do not overlap each other.
    """
    order = sorted(range(len(intervals)), key=lambda i: intervals[i])
    res: list[dict] = [{"n": 0, "bp": 0, "max_lod": 0} for _ in intervals]
    j = 0
    for i in order:
        s, e = intervals[i]
        # elements are short (median under 100 bp), so a forward pointer with a small look-back is enough
        while j < len(elements) and elements[j].end <= s:
            j += 1
        k = j
        while k < len(elements) and elements[k].start < e:
            el = elements[k]
            lo, hi = max(s, el.start), min(e, el.end)
            if lo < hi:
                r = res[i]
                r["n"] += 1
                r["bp"] += hi - lo
                r["max_lod"] = max(r["max_lod"], el.lod)
            k += 1
    for i, (s, e) in enumerate(intervals):
        res[i]["fraction"] = round(res[i]["bp"] / (e - s), 4) if e > s else 0.0
    return res
