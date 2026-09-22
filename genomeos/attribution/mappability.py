# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured mappability for the library's oligos, in place of the repeat proxy.

`genomeos/attribution/measurability.py` splits the real unknown's oligos into a family that cannot
be synthesised and a family that can be ordered but cannot be **attributed to one locus**. The
second family was measured with a proxy: how much of the window RepeatMasker calls an interspersed
repeat, and how much of it lies in a curated segmental duplication. That proxy is not the quantity
the design needs. It swings the attributable count from 19,084 to 28,563 oligos as the repeat share
cut-off moves from 0.25 to 0.75 -- a 50% swing on an arbitrary number -- and the one rule that did
ask about uniqueness, `non_unique_in_block`, compared a window only against its own block and so
fired on 23 windows: a lower bound, not a measurement.

This module reads the quantity itself. **Umap multi-read mappability** (Karimzadeh, Ernst, Kundaje
and Hoffman 2018, *Nucleic Acids Research* 46:e120, doi:10.1093/nar/gky677) gives every base of
GRCh38 the share of the length-`k` reads overlapping it that align to exactly one place in the
genome. A base at 1.0 can be read back from anywhere; a base at 0.0 is covered by no uniquely
placeable read at all. That is attribution, measured against the whole assembly rather than against
one 10 kb block.

**k is part of every number here and is never pooled across tracks.** A 24-mer and a 100-mer answer
different questions about a 300 bp oligo, and the four tracks are read and reported in four separate
columns. The tracks are UCSC's copies of the Hoffman lab files, `k{k}.Umap.MultiTrackMappability.bw`
under `/gbdb/hg38/hoffmanMappability/`, GRCh38 (UCSC hg38) throughout.

**Nothing is downloaded.** The four files are 2.0, 1.6, 1.3 and 0.9 GB; they are read through
`genomeos.attribution.bigwig` by HTTP range request, header and R-tree first and then the data
sections that overlap the oligos, and every run reports the bytes it moved against `BYTE_CAP`. Only
the per-oligo summaries are kept, in `data/knowledge/mappability/`, which `cache_bytes` holds under
`CACHE_CAP`.

**What "mappable" means for a 300 bp oligo** is one named rule, reported at three values so the
reader sees the sensitivity rather than a chosen number:

    a base is uniquely readable when its multi-read mappability is at least BASE_CUT = 1.0, and an
    oligo is mappable when at least MAPPABLE_FRACTION of its 300 bases are uniquely readable.

`BASE_CUT` is 1.0 and not a fraction because anything below 1.0 already means some read covering the
base has a second home in the genome, which is precisely the failure the design is trying to avoid;
the tolerance the design actually has is over *how much of the oligo* may fail, and that is
`MAPPABLE_FRACTION`, reported at 0.50, 0.90 and 0.99.

The track is **not merged with the proxy**. Both arms are counted, and the cross-tabulation between
them is the result: the cell that matters is the oligo the proxy passes and the track calls
unmappable, because the design would order that oligo believing its answer could be placed.

One property of the track has to be stated before any count is read from it: Umap stores no value
where mappability is zero, so a window with no value at all is measured as unmappable and not as
unassessed. That reading is only safe because the library's test arm contains **not one N base**
(`scripts/mappability.py` re-checks it), so a missing value cannot be an assembly gap.
"""

from __future__ import annotations

import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from genomeos.attribution import measurability as meas
from genomeos.attribution.bigwig import BigWig
from genomeos.compare import Strata, standardised

# --- the tracks -------------------------------------------------------------------------------
ASSEMBLY = "GRCh38 (UCSC hg38)"
TRACK_KIND = "Umap multi-read mappability (single-read mappability averaged over the k-mers covering a base)"
CITATION = (
    "Karimzadeh M, Ernst C, Kundaje A, Hoffman MM 2018, Umap and Bismap: quantifying genome and "
    "methylome mappability, Nucleic Acids Research 46:e120, doi:10.1093/nar/gky677"
)
HOST = "https://hgdownload.soe.ucsc.edu/gbdb/hg38/hoffmanMappability"


@dataclass(frozen=True)
class Track:
    """One mappability track. `size_bytes` is what the file weighs, not what is ever fetched."""

    k: int
    url: str
    size_bytes: int
    kind: str = TRACK_KIND
    assembly: str = ASSEMBLY
    citation: str = CITATION

    @property
    def name(self) -> str:
        return f"Umap k{self.k}"

    def as_dict(self) -> dict[str, Any]:
        return {
            "track": self.name,
            "k": self.k,
            "kind": self.kind,
            "assembly": self.assembly,
            "url": self.url,
            "file_bytes": self.size_bytes,
        }


# sizes checked against the Content-Range of a 64-byte request on 2026-09-17; a host that disagrees
# with one of these is a different file and the reader's own size check will refuse it
TRACKS: dict[int, Track] = {
    24: Track(24, f"{HOST}/k24.Umap.MultiTrackMappability.bw", 2_028_863_031),
    36: Track(36, f"{HOST}/k36.Umap.MultiTrackMappability.bw", 1_647_208_796),
    50: Track(50, f"{HOST}/k50.Umap.MultiTrackMappability.bw", 1_331_314_604),
    100: Track(100, f"{HOST}/k100.Umap.MultiTrackMappability.bw", 864_604_710),
}
DEFAULT_KS = (24, 36, 50, 100)

# --- the rule ---------------------------------------------------------------------------------
BASE_CUT = 1.0  # a base is uniquely readable only when every k-mer covering it maps to one place
MAPPABLE_FRACTION = 0.90  # share of an oligo's bases that must be uniquely readable
SENSITIVITY_FRACTIONS = (0.50, 0.90, 0.99)  # the rule reported at three values, not tuned to one

# --- the proxy this replaces, at its own three values -----------------------------------------
PROXY_INTERSPERSED = (0.25, 0.50, 0.75)  # measurability.SENSITIVITY["interspersed_max"]

# --- budgets ----------------------------------------------------------------------------------
BYTE_CAP = 96 << 20  # per track and run; a run that reaches it stops and names the unread oligos
BATCH = 2_000  # oligos per summarise call, so the byte cap is checked often enough to bind
CACHE_DIR = Path("data/knowledge/mappability")
CACHE_CAP = 100 << 20  # the whole cache; only per-oligo summaries are kept, never track bytes

OLIGO = meas.OLIGO

# the four cells of the cross-tabulation, named so neither arm can be read as the other
CELLS = (
    "both_flag",  # proxy flags it, track calls it unmappable
    "proxy_flags_track_mappable",  # the proxy's false alarms against the measurement
    "proxy_passes_track_unmappable",  # the dangerous cell: would be ordered as attributable
    "neither_flags",
)

# Bins for the matched comparison. Block length because an oligo's own length is 300 by
# construction; GC and distance to a coding TSS because those are the covariates that confounded
# every comparison this project corrected in the week before this was written.
STRATA = Strata(
    length=(1_000, 3_000, 10_000, 30_000),
    gc=(0.30, 0.35, 0.40, 0.45, 0.55),
    tss_distance=(1_000, 10_000, 50_000, 200_000, 1_000_000),
)


# --- reading the track ------------------------------------------------------------------------


@dataclass
class ReadResult:
    """What one chromosome's read cost and what it could not reach."""

    rows: dict[tuple[int, int], dict[str, Any]] = field(default_factory=dict)
    bytes_fetched: int = 0
    requests: int = 0
    seconds: float = 0.0
    in_track: bool = True
    cap_reached: bool = False
    unread: int = 0
    from_cache: int = 0

    def as_dict(self) -> dict[str, Any]:
        return {
            "oligos_read": len(self.rows),
            "oligos_from_cache": self.from_cache,
            "oligos_unread": self.unread,
            "bytes_fetched": self.bytes_fetched,
            "range_requests": self.requests,
            "seconds": round(self.seconds, 1),
            "chromosome_in_track": self.in_track,
            "byte_cap_reached": self.cap_reached,
        }


def cache_path(k: int, chrom: str, cache_dir: Path = CACHE_DIR) -> Path:
    return cache_dir / f"umap_k{k}_{chrom}.tsv"


def cache_bytes(cache_dir: Path = CACHE_DIR) -> int:
    if not cache_dir.is_dir():
        return 0
    return sum(p.stat().st_size for p in cache_dir.glob("*.tsv"))


def load_cache(k: int, chrom: str, cache_dir: Path = CACHE_DIR) -> dict[tuple[int, int], dict[str, Any]]:
    """Per-oligo summaries kept from an earlier run: start, end, bases with a value, bases at the
    cut, and the sum of the values. Nothing here is track bytes."""
    p = cache_path(k, chrom, cache_dir)
    if not p.exists():
        return {}
    out: dict[tuple[int, int], dict[str, Any]] = {}
    with p.open() as fh:
        header = fh.readline()
        if not header.startswith("start\t"):
            return {}
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 5:
                continue
            s, e = int(parts[0]), int(parts[1])
            out[(s, e)] = {
                "start": s,
                "end": e,
                "bases": int(parts[2]),
                "above": int(parts[3]),
                "total": float(parts[4]),
            }
    return out


def save_cache(
    k: int,
    chrom: str,
    rows: dict[tuple[int, int], dict[str, Any]],
    cache_dir: Path = CACHE_DIR,
    cap: int = CACHE_CAP,
) -> bool:
    """Write the summaries; refuse and say so when the cache would pass `cap`."""
    cache_dir.mkdir(parents=True, exist_ok=True)
    p = cache_path(k, chrom, cache_dir)
    body = ["start\tend\tbases\tabove\ttotal"]
    for (s, e), r in sorted(rows.items()):
        body.append(f"{s}\t{e}\t{r['bases']}\t{r['above']}\t{r['total']:.6g}")
    text = "\n".join(body) + "\n"
    already = p.stat().st_size if p.exists() else 0
    if cache_bytes(cache_dir) - already + len(text.encode()) > cap:
        return False
    p.write_text(text)
    return True


def read_chromosome(
    k: int,
    chrom: str,
    intervals: Sequence[tuple[int, int]],
    source: str | Path | None = None,
    byte_cap: int = BYTE_CAP,
    batch: int = BATCH,
    cache_dir: Path | None = CACHE_DIR,
    refresh: bool = False,
) -> ReadResult:
    """Summarise `intervals` of one chromosome against the k-mer track, by range request.

    Intervals are sorted and de-duplicated here because the bigWig reader requires them not to
    overlap. Oligos already in the cache cost no bytes. The read stops when `byte_cap` is reached
    and the oligos it did not get to are counted, never silently zeroed.
    """
    track = TRACKS[k]
    want = sorted(set(intervals))
    res = ReadResult()
    cached = {} if refresh or cache_dir is None else load_cache(k, chrom, cache_dir)
    todo = [iv for iv in want if iv not in cached]
    for iv in want:
        if iv in cached:
            res.rows[iv] = cached[iv]
    res.from_cache = len(res.rows)
    if not todo:
        return res
    t0 = time.time()
    bw = BigWig(source if source is not None else track.url)
    try:
        if chrom not in bw.chroms:
            res.in_track = False
            res.unread = len(todo)
            res.bytes_fetched = bw.src.bytes_fetched
            res.requests = bw.src.requests
            res.seconds = time.time() - t0
            return res
        done = 0
        for i in range(0, len(todo), batch):
            if bw.src.bytes_fetched >= byte_cap:
                res.cap_reached = True
                break
            chunk = todo[i : i + batch]
            for iv, st in zip(chunk, bw.summarise(chrom, list(chunk), BASE_CUT), strict=True):
                res.rows[iv] = {
                    "start": iv[0],
                    "end": iv[1],
                    "bases": st.bases,
                    "above": st.above,
                    "total": round(st.total, 4),
                }
            done += len(chunk)
        res.unread = len(todo) - done
        res.bytes_fetched = bw.src.bytes_fetched
        res.requests = bw.src.requests
    finally:
        bw.close()
    res.seconds = time.time() - t0
    if cache_dir is not None and res.rows:
        save_cache(k, chrom, {**cached, **res.rows}, cache_dir)
    return res


# --- the rule, applied ------------------------------------------------------------------------


def window_bp(row: dict[str, Any]) -> int:
    return row["end"] - row["start"]


def unique_fraction(row: dict[str, Any]) -> float:
    """Share of the oligo's bases at or above `BASE_CUT`.

    The denominator is the oligo's length and not the bases the track carries a value for: Umap
    writes nothing where mappability is zero, so a base with no value is a measured zero.
    """
    return row["above"] / window_bp(row)


def mean_mappability(row: dict[str, Any]) -> float:
    """Mean multi-read mappability over the oligo, missing values counted as the zeros they are."""
    return row["total"] / window_bp(row)


def is_mappable(row: dict[str, Any], min_fraction: float = MAPPABLE_FRACTION) -> bool:
    return unique_fraction(row) >= min_fraction


def proxy_flagged(
    row: dict[str, Any],
    interspersed_max: float = meas.INTERSPERSED_MAX,
    segdup_max: float = meas.SEGDUP_MAX,
) -> bool:
    """Family B of `measurability`: the repeat and duplication proxy for "not one locus"."""
    return row["interspersed_fraction"] >= interspersed_max or row["segdup_fraction"] >= segdup_max


# --- counting ---------------------------------------------------------------------------------


def crosstab(
    rows: Sequence[dict[str, Any]],
    min_fraction: float = MAPPABLE_FRACTION,
    interspersed_max: float = meas.INTERSPERSED_MAX,
    segdup_max: float = meas.SEGDUP_MAX,
) -> dict[str, Any]:
    """The 2x2 between the repeat proxy and the measured track, with both disagreements named."""
    cells = dict.fromkeys(CELLS, 0)
    for r in rows:
        p = proxy_flagged(r, interspersed_max, segdup_max)
        m = is_mappable(r, min_fraction)
        if p and not m:
            cells["both_flag"] += 1
        elif p and m:
            cells["proxy_flags_track_mappable"] += 1
        elif m:
            cells["neither_flags"] += 1
        else:
            cells["proxy_passes_track_unmappable"] += 1
    n = len(rows)
    agree = cells["both_flag"] + cells["neither_flags"]
    return {
        "oligos": n,
        "min_fraction": min_fraction,
        "interspersed_max": interspersed_max,
        "segdup_max": segdup_max,
        **cells,
        "agreement": round(agree / n, 4) if n else None,
        "proxy_flags_track_mappable_share": round(cells["proxy_flags_track_mappable"] / n, 4) if n else None,
        "proxy_passes_track_unmappable_share": (
            round(cells["proxy_passes_track_unmappable"] / n, 4) if n else None
        ),
    }


def arms(
    rows: Sequence[dict[str, Any]],
    min_fraction: float = MAPPABLE_FRACTION,
    interspersed_max: float = meas.INTERSPERSED_MAX,
    segdup_max: float = meas.SEGDUP_MAX,
) -> dict[str, Any]:
    """The attributable count once per arm, and their conjunction, each named.

    Three numbers and not one: the repeat proxy's arm, the mappability track's arm, and the oligos
    both arms keep. They are never added or averaged, because they are answers to two different
    questions and the whole point of measuring the track was that the proxy was standing in for it.
    """
    n = len(rows)
    by_proxy = sum(1 for r in rows if not proxy_flagged(r, interspersed_max, segdup_max))
    by_track = sum(1 for r in rows if is_mappable(r, min_fraction))
    both = sum(
        1 for r in rows if not proxy_flagged(r, interspersed_max, segdup_max) and is_mappable(r, min_fraction)
    )
    return {
        "oligos_assessed": n,
        "attributable_by_repeat_proxy": by_proxy,
        "attributable_by_mappability_track": by_track,
        "attributable_by_both_arms": both,
        "share_by_repeat_proxy": round(by_proxy / n, 4) if n else None,
        "share_by_mappability_track": round(by_track / n, 4) if n else None,
        "share_by_both_arms": round(both / n, 4) if n else None,
    }


def sensitivity(
    rows: Sequence[dict[str, Any]], fractions: Sequence[float] = SENSITIVITY_FRACTIONS
) -> list[dict[str, Any]]:
    """The rule at each of its three values: the arm, and the cross-tabulation, per value."""
    out = []
    for f in fractions:
        a = arms(rows, f)
        x = crosstab(rows, f)
        cells = {k: x[k] for k in CELLS}
        out.append({"min_fraction": f, **a, "crosstab": cells, "agreement": x["agreement"]})
    return out


def proxy_sensitivity(
    rows: Sequence[dict[str, Any]],
    values: Sequence[float] = PROXY_INTERSPERSED,
    min_fraction: float = MAPPABLE_FRACTION,
) -> list[dict[str, Any]]:
    """The proxy's own swing, held against a track that does not move with it."""
    out = []
    for v in values:
        a = arms(rows, min_fraction, interspersed_max=v)
        x = crosstab(rows, min_fraction, interspersed_max=v)
        out.append(
            {
                "interspersed_max": v,
                "attributable_by_repeat_proxy": a["attributable_by_repeat_proxy"],
                "attributable_by_mappability_track": a["attributable_by_mappability_track"],
                "crosstab": {k: x[k] for k in CELLS},
            }
        )
    return out


# --- groups -----------------------------------------------------------------------------------


def group_covariates(rows: Iterable[dict[str, Any]]) -> dict[str, Any]:
    """The three covariates every group in this project prints, via `measurability.covariates`.

    An oligo's own length is `OLIGO` by construction, so the length distribution reported is the
    distribution of the blocks the oligos were tiled from, which is the length that varies.
    """
    shaped = [{"length": r["block_length"], "gc": r["gc"], "tss_distance": r["tss_distance"]} for r in rows]
    cov = meas.covariates(shaped)
    return {
        "oligos": cov.pop("n"),
        "oligo_bp": OLIGO,
        "median_block_length": cov.pop("median_length"),
        "block_length_quartiles": cov.pop("length_quartiles"),
        **cov,
    }


def groups(
    rows: Sequence[dict[str, Any]], min_fraction: float = MAPPABLE_FRACTION
) -> dict[str, dict[str, Any]]:
    """Covariates for every group this reports, including all four cells of the cross-tabulation.

    No tier is used as a control anywhere: the groups here are cut by the two rules themselves.
    """
    out: dict[str, list[dict[str, Any]]] = {name: [] for name in CELLS}
    out["track_mappable"] = []
    out["track_unmappable"] = []
    out["proxy_flagged"] = []
    out["proxy_passed"] = []
    for r in rows:
        p = proxy_flagged(r)
        m = is_mappable(r, min_fraction)
        out["track_mappable" if m else "track_unmappable"].append(r)
        out["proxy_flagged" if p else "proxy_passed"].append(r)
        if p and not m:
            out["both_flag"].append(r)
        elif p and m:
            out["proxy_flags_track_mappable"].append(r)
        elif m:
            out["neither_flags"].append(r)
        else:
            out["proxy_passes_track_unmappable"].append(r)
    return {k: group_covariates(v) for k, v in out.items()}


def proxy_against_track(
    rows: Sequence[dict[str, Any]], min_fraction: float = MAPPABLE_FRACTION
) -> dict[str, Any]:
    """Does the proxy predict the measurement once length, GC and TSS distance are held fixed?

    `genomeos.compare.standardised`, not a rate against a rate: targets are the oligos the proxy
    flags, controls are the oligos it passes, and the outcome is the track calling the oligo
    unmappable. The matched difference, the matched n, the targets dropped for want of a control and
    both arms' covariate medians all come back in the same dict.
    """
    shaped = [
        {
            "length": r["block_length"],
            "gc": r["gc"],
            "tss_distance": r["tss_distance"],
            "track_unmappable": not is_mappable(r, min_fraction),
            "proxy": proxy_flagged(r),
        }
        for r in rows
        if r["tss_distance"] is not None and r["gc"] is not None
    ]
    out = standardised(
        [r for r in shaped if r["proxy"]],
        [r for r in shaped if not r["proxy"]],
        STRATA,
        hit="track_unmappable",
    )
    out["excluded_without_a_covariate"] = len(rows) - len(shaped)
    out["reading"] = (
        "a positive matched difference means the proxy's flag does carry information about measured "
        "mappability at equal block length, GC and TSS distance; it does not make the two the same "
        "quantity, and the cross-tabulation is what says how far apart they are"
    )
    return out
