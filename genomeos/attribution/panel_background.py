# SPDX-License-Identifier: AGPL-3.0-or-later
"""What the human panel's matched background is made of, and how much of its offset that explains.

On 2026-09-16 the panel's GC- and replication-timing-matched background turned out not to be a
neutral baseline: every non-coding tier reads above it (regulatory 23 of 23 chromosomes, fossil 22,
neutral 21), a systematic +6 to +11 points that belongs to the control rather than to the tier. That
reading said what the offset *is not*; it did not say what the background *holds*. This module
measures it.

The background is rebuilt exactly as `human_panel.build_background` builds it -- the kilobases of the
hg38 grid that lie inside the panel's alignment blocks, with canonical CDS bases and the events
inside them removed, binned by GC and by replication timing -- and every kilobase is then annotated
from local data only:

- coding bases: GENCODE canonical CDS (what the panel already removes), CDS of *any* transcript, and
  exons of any gene including non-coding ones and pseudogenes;
- conserved elements: phastCons 100-way, from the cache `attribution/constraint.py` writes;
- constrained bases: Zoonomia phyloP over 241 mammals, sampled by range request when asked for;
- segmental duplication and the RepeatMasker classes, from the per-chromosome BED tracks;
- registry coverage: ENCODE cCREs;
- distance to the nearest coding TSS, GC, and the replication timing the panel matched on.

The offset is decomposed the one way that is not a story: the background is rebuilt with each
candidate taken out of it, through the panel's own `build_background`, and the tier ratios are
recomputed. What a covariate explains is how far it moves the ratio towards 1; what is left over is
reported as left over rather than chased with a fifth covariate.

Two-group comparisons go through `genomeos.compare.standardised`, never through a ratio of means.
"""

from __future__ import annotations

import gzip
import json
import random
from array import array
from bisect import bisect_left, bisect_right
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from statistics import median
from typing import Any

from genomeos.attribution.human_panel import (
    BACKGROUND_BIN,
    CACHE,
    DELETED,
    GC_STRATA,
    MIN_RECURRING,
    PRESENT,
    REPLACED,
    TIERS,
    EventIndex,
    Panel,
    build_background,
    coding_genes,
    events,
    measure,
    merge_intervals,
    rt_edges,
    rt_over,
    rt_stratum,
)
from genomeos.compare import Strata, standardised
from genomeos.results import RESULTS_DIR, load_result

CELL = 100  # feature coverage is kept per 100 bases, the resolution the panel already stores GC at
PROMOTER = 2_000  # a base this far from a coding TSS or nearer is promoter-proximal
CONSTRAINT_CACHE = Path("data/knowledge/constraint")
TRACKS = Path("data/results")
REPEAT_CLASSES = ("SINE", "LINE", "LTR", "DNA", "Satellite", "Simple_repeat", "Low_complexity")
GROUPS = ("cds", *TIERS)
CANDIDATES = ("cds_any", "exon_any", "conserved_elements", "promoter", "segmental_duplication")
TOGETHER = ("exon_any", "conserved_elements", "promoter")  # the three that are not repeats
WINDOW_CAP = 120_000  # windows per chromosome in the standardised comparison; longer ones subsample
TSS_EDGES = (1_000, 5_000, 20_000, 100_000)
SPLIT_AT = 200  # a kilobase holding this many bases of a feature is counted as holding it
FAR = 10**9  # the distance reported where a chromosome has no coding TSS at all
SEED = 20260917

# -- coverage, which is not a covariate ----------------------------------------------------------
# A covariate describes the sequence; coverage describes what was done to it, and standardising on
# the first while the arms differ in the second manufactures a difference out of effort. In this lane
# there is one effort axis and this is it: an event is *recurring* only when two assemblies carry the
# minority state, so a window that only 40 of the panel's 89 assemblies reached cannot show what a
# window all 89 reached can. COVERAGE_MEASURE names it; COVERAGE_EDGES stratifies it; the kilobases
# under COVERAGE_MIN are the ones the coverage-conditioned background leaves out.
COVERAGE_MEASURE = "assemblies informing the window, as a fraction of the panel's assemblies"
COVERAGE_EDGES = (0.5, 0.8, 0.95, 0.99)
COVERAGE_MIN = 0.95
COVERAGE_KEY = "without_low_assembly_coverage"
# On a sex chromosome most assemblies inform nothing, so the coverage cut would leave a background of
# scraps and a rate built from them means nothing. Below this share of kilobases kept, the arm reports
# itself unassessed with the reason rather than returning a number.
COVERAGE_KEEP_MIN = 0.5
COVERAGE_TOO_THIN = "conditioning on coverage leaves under half the background's kilobases"

# The two sentences this measurement is allowed to end on, as data rather than as prose, so that a
# later reader cannot upgrade the weaker one to the stronger one by paraphrase.
CLAIM_WITHOUT_COVERAGE = "the covariates I measured do not explain it"
CLAIM_WITH_COVERAGE = "nothing I measured explains it, including coverage"

SOURCES = {
    "coding": "GENCODE 50 (GRCh38.p14): canonical CDS, all-transcript CDS, all-gene exons, coding TSS",
    "conserved_elements": "UCSC phastConsElements100way, the cache attribution/constraint.py writes",
    "constrained_bases": "Zoonomia cactus241way phyloP >= 2.27 (5% FDR), read by range, sampled",
    "segmental_duplication": "UCSC genomicSuperDups, per-chromosome BED (genomeos duplications)",
    "repeats": "UCSC RepeatMasker, per-chromosome BED, by class",
    "registry": "ENCODE cCREs v3, per-chromosome BED",
    "replication_timing": "ENCODE UW Repli-seq lifted to hg38, the timing the panel matched on",
}


# ================================================================================================
# coverage arrays: bases of a feature per 100-base cell
# ================================================================================================
def cover_cells(intervals: list[tuple[int, int]], cells: int, cell: int = CELL) -> array:
    """Bases covered per cell, for one feature. Intervals are merged first, so no base counts twice."""
    a = array("B", bytes(cells))
    clipped = [(max(0, s), min(e, cells * cell)) for s, e in intervals if e > s]
    for s, e in merge_intervals([iv for iv in clipped if iv[1] > iv[0]]):
        i, j = s // cell, (e - 1) // cell
        if i == j:
            a[i] = min(cell, a[i] + (e - s))
            continue
        a[i] = min(cell, a[i] + ((i + 1) * cell - s))
        for k in range(i + 1, j):
            a[k] = cell
        a[j] = min(cell, a[j] + (e - j * cell))
    return a


def cells_bp(a: array, start: int, end: int, cell: int = CELL) -> int:
    """Bases of the feature inside [start, end); the range is rounded out to whole cells."""
    return sum(a[start // cell : (end - 1) // cell + 1])


def span_bp(a: array, s: int, e: int) -> int:
    """Feature bases in [s, e), at cell resolution and never more than the span itself."""
    return min(cells_bp(a, s, e), e - s)


def informing(panel: Panel) -> list[int]:
    """Per alignment block: how many assemblies have a line over it, whatever that line says.

    Present, deleted and replaced all inform; a padded `.` does not. This is the panel's measurement
    effort over a block, and it is what `coverage_over` turns into a per-window fraction.
    """
    return [
        st.count(PRESENT) + st.count(DELETED) + sum(st.count(c) for c in REPLACED)
        for st in panel.block_status
    ]


def coverage_over(panel: Panel, per_block: list[int], s: int, e: int) -> tuple[float | None, int]:
    """The mean fraction of the panel's assemblies informing [s, e), and the aligned bases behind it."""
    seen = bp = 0
    for i in panel.blocks_over(s, e):
        ov = min(e, panel.block_end[i]) - max(s, panel.block_start[i])
        if ov > 0:
            bp += ov
            seen += ov * per_block[i]
    return (seen / (bp * panel.n) if bp and panel.n else None), bp


def nearest(sorted_points: list[int], pos: int) -> int:
    """Distance from `pos` to the nearest point, or FAR when the chromosome has none."""
    if not sorted_points:
        return FAR
    i = bisect_left(sorted_points, pos)
    best = FAR
    if i < len(sorted_points):
        best = sorted_points[i] - pos
    if i:
        best = min(best, pos - sorted_points[i - 1])
    return best


# ================================================================================================
# local tracks and gene models
# ================================================================================================
def _bed(path: Path, start_col: int, end_col: int) -> list[tuple[int, int]]:
    out: list[tuple[int, int]] = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            try:
                out.append((int(f[start_col]), int(f[end_col])))
            except (IndexError, ValueError):
                continue
    return out


def gene_features(chrom: str) -> dict[str, Any]:
    """All-transcript CDS, all-gene exons and the coding TSS positions, from GENCODE."""
    from genomeos.genome.annotation import Annotation, default_gencode

    gff = default_gencode({chrom})
    if gff is None:
        return {}
    ann = Annotation.from_gff3(gff, {chrom})
    cds_any: list[tuple[int, int]] = []
    exon_any: list[tuple[int, int]] = []
    tss: list[int] = []
    for g in ann.genes.values():
        if g.locus.chrom != chrom:
            continue
        for t in g.transcripts.values():
            exon_any += [(x.start, x.end) for x in t.exons]
            cds_any += [(c.start, c.end) for c, _ in t.cds]
        if g.type == "protein_coding":
            tss.append(g.locus.start if "MINUS" not in str(g.locus.strand) else g.locus.end)
    return {"cds_any": merge_intervals(cds_any), "exon_any": merge_intervals(exon_any), "tss": sorted(tss)}


def local_tracks(chrom: str, tracks: Path = TRACKS, cache: Path = CONSTRAINT_CACHE) -> dict[str, Any]:
    """Every per-chromosome track this measurement reads; the ones absent are named, not zeroed."""
    out: dict[str, Any] = {"missing": []}
    for name, (p, s, e) in {
        "segmental_duplication": (tracks / f"superdups_{chrom}.bed.gz", 0, 1),
        "registry": (tracks / f"ccres_{chrom}.bed.gz", 1, 2),
        "conserved_elements": (cache / f"phastConsElements100way_{chrom}.bed.gz", 1, 2),
    }.items():
        if not p.exists():
            out["missing"].append(name)
            continue
        out[name] = merge_intervals(_bed(p, s, e))
    rmsk = tracks / f"rmsk_{chrom}.bed.gz"
    if not rmsk.exists():
        out["missing"].append("repeats")
        return out
    rows: list[tuple[int, int, str]] = []
    with gzip.open(rmsk, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 2:
                rows.append((int(f[0]), int(f[1]), f[2]))
    out["repeats"] = merge_intervals([(s, e) for s, e, _ in rows])
    for k in REPEAT_CLASSES:
        out[f"repeat_{k}"] = merge_intervals([(s, e) for s, e, c in rows if c == k])
    return out


# ================================================================================================
# the rebuilt panel, exactly as build() sets it up
# ================================================================================================
@dataclass
class Rebuilt:
    """One chromosome's store, its events, its timing and the blocks the committed result holds."""

    chrom: str
    panel: Panel
    evs: list[Any]
    index: EventIndex
    rt: dict[int, float]
    edges: list[float]
    cds_all: list[tuple[int, int]]
    blocks: list[dict[str, Any]]
    length: int
    committed: dict[str, Any] = field(default_factory=dict)
    _measured: dict[str, dict[str, Any]] = field(default_factory=dict)

    def tier_intervals(self, tier: str) -> list[tuple[int, int]]:
        if tier == "cds":
            return self.cds_all
        return [(b["start"], b["end"]) for b in self.blocks if b["tier"] == tier]

    def measured(self, tier: str) -> dict[str, Any]:
        """The panel's numbers over a tier, computed once: only the background changes below."""
        if tier not in self._measured:
            ivs = self.tier_intervals(tier)
            self._measured[tier] = measure(self.panel, ivs, self.index.over(ivs)) if ivs else {}
        return self._measured[tier]


def rebuild(chrom: str, cache: Path = CACHE, results_dir: Path | None = None) -> Rebuilt:
    """Load the store and reproduce the panel's inputs: no network, no model, nothing re-derived."""
    results_dir = results_dir or RESULTS_DIR
    panel = Panel(cache / chrom)
    committed = load_result(f"human_panel_{chrom}", results_dir) or {}
    budget = load_result(f"budget_{chrom}", results_dir) or {}
    length = budget.get("chromosome_length") or panel.block_end[-1]
    evs = events(panel, [(0, length)])
    rt_path = cache / chrom / "replication_timing.json"
    rt = {int(k): v for k, v in json.loads(rt_path.read_text()).items()} if rt_path.exists() else {}
    genes = coding_genes(chrom)
    return Rebuilt(
        chrom=chrom,
        panel=panel,
        evs=evs,
        index=EventIndex(evs),
        rt=rt,
        edges=rt_edges(rt),
        cds_all=merge_intervals([iv for _, ivs in genes for iv in ivs]),
        blocks=committed.get("blocks", []),
        length=length,
        committed=committed,
    )


def reproduces_the_panel(rb: Rebuilt, built: dict[str, Any]) -> dict[str, Any]:
    """Does the rebuilt background give back the committed result's own ratios?

    Nothing below means anything if the rebuild is not the panel's. The committed
    `pooled[tier].ratio_gc_rt` is the number the tier baseline was read from, and this says how far
    each rebuilt ratio sits from it. Tiers the committed result does not carry are named, not zeroed.
    """
    pooled = rb.committed.get("pooled") or {}
    diffs, absent = {}, []
    for name in GROUPS:
        if name not in pooled or built.get(name) is None:
            absent.append(name)
            continue
        diffs[name] = round(built[name] - pooled[name]["ratio_gc_rt"], 4)
    return {
        "tiers_checked": len(diffs),
        "tiers_not_in_the_committed_result": absent,
        "largest_absolute_difference": max((abs(v) for v in diffs.values()), default=None),
        "per_tier": diffs,
    }


def tier_ratios(rb: Rebuilt, exclude: list[tuple[int, int]]) -> dict[str, float | None]:
    """Every tier's ratio against a background built with `exclude` taken out of it."""
    bg = build_background(rb.panel, exclude, rb.evs, rb.rt, rb.edges)
    out: dict[str, float | None] = {}
    for name in GROUPS:
        ivs = rb.tier_intervals(name)
        if not ivs:
            continue
        exp = bg.expected(rb.panel, ivs)
        out[name] = round(rb.measured(name)["recurring_events"] / exp, 4) if exp else None
    out["background_per_kb"] = round(1000 * bg.overall, 4)
    out["excluded_bases"] = sum(e - s for s, e in merge_intervals(exclude))
    return out


# ================================================================================================
# the background's kilobases, and what each of them holds
# ================================================================================================
def background_bins(
    rb: Rebuilt, feats: dict[str, array], tss: list[int]
) -> tuple[list[dict], dict[str, int]]:
    """One row per kilobase the panel's background is built from, and a named count for the rest.

    The keep rule is `build_background`'s own: a kilobase counts towards a stratum's rate when at
    least half of it survives the canonical-CDS cut and the store has a GC value for it.
    """
    panel, w = rb.panel, BACKGROUND_BIN
    per_block = informing(panel)
    aligned: Counter = Counter()
    for i in range(len(panel.block_start)):
        s, e = panel.block_start[i], panel.block_end[i]
        for b in range(s // w, (e - 1) // w + 1):
            lo, hi = max(s, b * w), min(e, (b + 1) * w)
            if lo < hi:
                aligned[b] += hi - lo
    ev_count: Counter = Counter()
    for ev in rb.evs:
        if ev.minor >= MIN_RECURRING:
            ev_count[ev.start // w] += 1
    rows: list[dict] = []
    dropped: Counter = Counter()
    for b in sorted(aligned):
        s, e = b * w, (b + 1) * w
        kept = max(0, aligned[b] - span_bp(feats["cds_canonical"], s, e))
        if kept < w / 2:
            dropped["below_half_a_kilobase_after_the_coding_cut"] += 1
            continue
        gc = panel.gc_fraction([(s, e)])
        if gc is None:
            dropped["no_gc_in_the_store"] += 1
            continue
        st = rt_stratum(rt_over(rb.rt, [(s, e)]), rb.edges)
        row = {
            "bin": b,
            "start": s,
            "end": e,
            "aligned": aligned[b],
            "kept": kept,
            "gc": round(gc, 4),
            "gc_stratum": bisect_right(GC_STRATA, gc),
            "rt_stratum": st,
            "tss": nearest(tss, s + w // 2),
            "coverage": coverage_over(panel, per_block, s, e)[0],
            "recurring": ev_count.get(b, 0),
        }
        for name, arr in feats.items():
            row[name] = span_bp(arr, s, e)
        rows.append(row)
    return rows, dict(dropped)


NO_TRACK = "no track for this chromosome"
NO_BASES = "the group has no bases on this chromosome"


def _features(shares: dict[str, tuple[int, int]], absent: list[str], why: str = NO_TRACK) -> dict[str, Any]:
    """Every feature with its denominator beside it: a share of nothing is not the same as no share.

    A feature whose track this chromosome does not carry gets `share: null` and `bases_assessed: 0`,
    never 0.0, so "measured and absent" and "never looked" cannot be read as the same number. `why`
    says which of the two an unassessed feature is.
    """
    out: dict[str, Any] = {}
    for f, (covered, assessed) in shares.items():
        share = covered / assessed if assessed else None
        out[f] = {
            "share": round(share, 4) if share is not None else None,
            "bases_assessed": assessed,
            "reads": "not assessed" if share is None else _reads_share(share),
        }
        if not assessed:
            out[f]["why"] = why
    for f in absent:
        out[f] = {"share": None, "bases_assessed": 0, "reads": "not assessed", "why": NO_TRACK}
    return out


def _reads_share(share: float) -> str:
    """The cell in words, so that measured-and-empty cannot be read as never-looked-at."""
    if share == 0:
        return "measured, and none of it"
    return f"measured, {round(100 * share, 2)}% of the bases"


def composition(
    rows: list[dict], names: list[str], absent: list[str] | None = None, weight: str = "kept"
) -> dict[str, Any]:
    """A group of kilobases: its covariates and its feature shares, base-weighted, with denominators."""
    n = sum(r[weight] for r in rows)
    if not rows or not n:
        return {
            "kilobases": len(rows),
            "bases": n,
            "features": _features({f: (0, 0) for f in names}, absent or [], NO_BASES),
        }
    return {
        "kilobases": len(rows),
        "bases": n,
        "gc_median": round(median(r["gc"] for r in rows), 4),
        "tss_median": int(median(r["tss"] for r in rows)),
        "tss_assessed": sum(1 for r in rows if r["tss"] < FAR),
        "coverage_median": round(median(r["coverage"] for r in rows if r["coverage"] is not None), 4),
        "coverage_assessed": sum(1 for r in rows if r["coverage"] is not None),
        "recurring_per_kb": round(1000 * sum(r["recurring"] for r in rows) / n, 3),
        "features": _features({f: (sum(min(r[f], r[weight]) for r in rows), n) for f in names}, absent or []),
    }


def block_composition(
    rb: Rebuilt,
    tier: str,
    feats: dict[str, array],
    tss: list[int],
    names: list[str],
    absent: list[str] | None = None,
) -> dict[str, Any]:
    """The same table for a tier, over its own blocks rather than over the background's kilobases."""
    ivs = rb.tier_intervals(tier)
    if not ivs:
        return {
            "blocks": 0,
            "bases": 0,
            "features": _features({f: (0, 0) for f in names}, absent or [], NO_BASES),
        }
    m = rb.measured(tier)
    aligned = m["aligned_bases"] or 1
    spans = [e - s for s, e in ivs]
    gcs = [g for g in (rb.panel.gc_fraction([iv]) for iv in ivs) if g is not None]
    dists = [nearest(tss, (s + e) // 2) for s, e in ivs]
    total = sum(spans)
    return {
        "blocks": len(ivs),
        "bases": m["bases"],
        "aligned_bases": m["aligned_bases"],
        "length_median": int(median(spans)),
        "gc_median": round(median(gcs), 4) if gcs else None,
        "gc_assessed": len(gcs),
        "tss_median": int(median(dists)),
        "tss_assessed": sum(1 for d in dists if d < FAR),
        # the tier's own measurement effort: the share of assembly-bases that informed it at all
        "coverage_median": round(1 - m["missing_share"], 4) if m.get("missing_share") is not None else None,
        "coverage_assessed": len(ivs) if m.get("missing_share") is not None else 0,
        "recurring_per_kb": round(1000 * m["recurring_events"] / aligned, 3),
        "features": _features(
            {f: (sum(span_bp(feats[f], s, e) for s, e in ivs), total) for f in names}, absent or []
        ),
    }


def feature_splits(rows: list[dict], names: list[str], absent: list[str] | None = None) -> dict[str, Any]:
    """For each feature: the rate of the background's kilobases that hold it against those that do not.

    This is the descriptive half of the decomposition. A feature that makes the background quiet shows
    up here as a lower rate on the holding side, and as a share of the background large enough to move
    the average; a feature that is either rare or no quieter cannot explain the offset whatever else
    it correlates with. `kilobases_assessed` is the denominator: where a track is absent it is 0 and
    every rate is null, which is not the same as a feature that is present and never co-occurs.
    """
    out: dict[str, Any] = {}
    total = sum(r["kept"] for r in rows) or 1
    for f in names:
        with_, without = [r for r in rows if r[f] >= SPLIT_AT], [r for r in rows if r[f] < SPLIT_AT]
        nw = sum(r["kept"] for r in with_)
        nn = sum(r["kept"] for r in without)
        out[f] = {
            "kilobases_assessed": len(rows),
            "kilobases_with": len(with_),
            "share_of_background_bases": round(nw / total, 4),
            "with_per_kb": round(1000 * sum(r["recurring"] for r in with_) / nw, 3) if nw else None,
            "without_per_kb": round(1000 * sum(r["recurring"] for r in without) / nn, 3) if nn else None,
        }
    for f in absent or []:
        out[f] = {
            "kilobases_assessed": 0,
            "kilobases_with": None,
            "share_of_background_bases": None,
            "with_per_kb": None,
            "without_per_kb": None,
            "why": "no track for this chromosome",
        }
    return out


# ================================================================================================
# the standardised comparison, at the resolution where a rate becomes a share
# ================================================================================================
def window_rows(rb: Rebuilt, tss: list[int], cap: int = WINDOW_CAP) -> tuple[list[dict], dict[str, Any]]:
    """100-base windows with a boolean: does this window carry a recurring event?

    A rate cannot be standardised as a share, so the unit is made small enough that "carries one" is
    informative -- about eight recurring events per aligned kilobase, so roughly half the windows.
    Long chromosomes are subsampled by a fixed step, which is reported rather than hidden.

    Each window carries its coverage as well as its covariates: how much of the panel reached it is
    not a property of the sequence and must be held fixed separately from GC, timing and distance.
    """
    panel = rb.panel
    per_block = informing(panel)
    ev_cell: Counter = Counter()
    for ev in rb.evs:
        if ev.minor >= MIN_RECURRING:
            ev_cell[ev.start // CELL] += 1
    cells = sorted(
        {
            c
            for i in range(len(panel.block_start))
            for c in range(panel.block_start[i] // CELL, (panel.block_end[i] - 1) // CELL + 1)
        }
    )
    step = max(1, -(-len(cells) // cap))
    rows: list[dict] = []
    thin = no_coverage = 0
    for c in cells[::step]:
        s, e = c * CELL, (c + 1) * CELL
        cov, aligned = coverage_over(panel, per_block, s, e)
        gc = panel.gc_fraction([(s, e)])
        if aligned < CELL / 2 or gc is None:
            thin += 1
            continue
        if cov is None:
            no_coverage += 1
            continue
        rows.append(
            {
                "start": s,
                "length": aligned,
                "gc": round(gc, 4),
                "rt": rt_stratum(rt_over(rb.rt, [(s, e)]), rb.edges) or 0,
                "tss": nearest(tss, s + CELL // 2),
                "coverage": round(cov, 4),
                "recurring": bool(ev_cell.get(c, 0)),
            }
        )
    return rows, {
        "step": step,
        "cells_in_panel_blocks": len(cells),
        "windows": len(rows),
        "windows_too_thin_or_without_gc": thin,
        "windows_without_a_coverage_value": no_coverage,
        "coverage_measure": COVERAGE_MEASURE,
    }


def _windows_in(rows: list[dict], ivs: list[tuple[int, int]]) -> list[dict]:
    starts = [r["start"] for r in rows]
    out = []
    for s, e in merge_intervals(ivs):
        out += rows[bisect_left(starts, s) : bisect_left(starts, e)]
    return out


def standardised_against_background(
    rb: Rebuilt, rows: list[dict], cds: list[tuple[int, int]]
) -> dict[str, Any]:
    """Each tier against the panel's own control -- every non-coding window -- three ways.

    Holding GC and replication timing, which is the match the panel makes; then also holding distance
    to the nearest coding TSS, which is a covariate; then also holding coverage, which is not. The
    third is the one that matters: standardising on covariates is not conditioning on effort, and a
    difference that survives the first two and dies on the third was always an artefact of effort.
    `compare.standardised` carries the medians and the dropped count in all three.
    """
    coding = {r["start"] for r in _windows_in(rows, [(s - CELL + 1, e) for s, e in cds])}
    controls = [r for r in rows if r["start"] not in coding]
    ladder = {
        "gc_and_timing": Strata(length=(CELL - 1,), gc=GC_STRATA, rt=(0.5, 1.5)),
        "gc_timing_and_tss": Strata(length=(CELL - 1,), gc=GC_STRATA, rt=(0.5, 1.5), tss=TSS_EDGES),
        "gc_timing_tss_and_coverage": Strata(
            length=(CELL - 1,),
            gc=GC_STRATA,
            rt=(0.5, 1.5),
            tss=TSS_EDGES,
            coverage=COVERAGE_EDGES,
        ),
    }
    out: dict[str, Any] = {
        "controls": len(controls),
        "coding_windows_excluded": len(coding),
        "coverage_measure": COVERAGE_MEASURE,
    }
    for tier in TIERS:
        targets = [r for r in _windows_in(rows, rb.tier_intervals(tier)) if r["start"] not in coding]
        if not targets:
            continue
        out[tier] = {k: standardised(targets, controls, s, hit="recurring") for k, s in ladder.items()}
    return out


# ================================================================================================
# constrained bases, sampled
# ================================================================================================
def phylop_sample(chrom: str, rows: list[dict], n: int = 0, seed: int = SEED) -> dict[str, Any]:
    """Constrained-base fraction over a random sample of the background's kilobases.

    phyloP is read by range request and a whole chromosome would be a download in all but name, so it
    is sampled: the number of kilobases is the caller's, the sampling error is reported, and the
    kilobases not sampled are a counted category rather than a zero.
    """
    if n <= 0 or not rows:
        return {"kilobases_sampled": 0, "kilobases_not_sampled": len(rows), "why": "not requested"}
    from statistics import pstdev

    from genomeos.attribution.constraint import phylop_over_blocks

    rng = random.Random(seed)
    take = sorted(rng.sample(rows, min(n, len(rows))), key=lambda r: r["start"])
    stats, cost = phylop_over_blocks(chrom, [(r["start"], r["end"]) for r in take])
    seen = sum(s.bases for s in stats)
    above = sum(s.above for s in stats)
    per_kb = [s.above / s.bases for s in stats if s.bases]
    frac = above / seen if seen else None
    return {
        "kilobases_sampled": len(take),
        "kilobases_not_sampled": len(rows) - len(take),
        "kilobases_with_no_phylop": len(take) - len(per_kb),
        "bases_read": seen,
        "constrained_fraction": round(frac, 4) if frac is not None else None,
        # the kilobase is the sampling unit, so the error is taken over kilobases, not over bases
        "standard_error": round(pstdev(per_kb) / len(per_kb) ** 0.5, 4) if len(per_kb) > 1 else None,
        "cost": cost,
    }


def tier_constraint(rb: Rebuilt) -> dict[str, Any]:
    """The phyloP constrained fraction each tier's blocks already carry, from the committed result.

    The budget measured it block by block; nothing is re-read here. Blocks where it was never measured
    are counted, not treated as zero, so a tier with no constrained bases and a tier never looked at
    cannot be confused.
    """
    out: dict[str, Any] = {}
    for tier in TIERS:
        rows = [b for b in rb.blocks if b["tier"] == tier]
        vals = [(b["mammal_fraction"], b["end"] - b["start"]) for b in rows if b.get("mammal_fraction")]
        bp = sum(w for _, w in vals)
        out[tier] = {
            "blocks": len(rows),
            "blocks_assessed": len(vals),
            "constrained_fraction": round(sum(v * w for v, w in vals) / bp, 4) if bp else None,
        }
    return out


# ================================================================================================
# one chromosome
# ================================================================================================
def _ivs(tracks: dict, genes: dict, name: str) -> list[tuple[int, int]]:
    if name in ("cds_any", "exon_any"):
        return genes.get(name, [])
    if name == "promoter":
        return merge_intervals([(max(0, p - PROMOTER), p + PROMOTER) for p in genes.get("tss", [])])
    return tracks.get(name, [])


def feature_arrays(rb: Rebuilt, tracks: dict, genes: dict) -> tuple[dict[str, array], list[int]]:
    """The coverage arrays every table reads, including one per tier so a kilobase knows where it is."""
    cells = rb.length // CELL + 2
    sets: dict[str, list[tuple[int, int]]] = {"cds_canonical": rb.cds_all}
    for name in ("cds_any", "exon_any", "promoter"):
        sets[name] = _ivs(tracks, genes, name)
    for name in ("conserved_elements", "segmental_duplication", "registry", "repeats"):
        if name in tracks:
            sets[name] = tracks[name]
    for c in REPEAT_CLASSES:
        if f"repeat_{c}" in tracks:
            sets[f"repeat_{c}"] = tracks[f"repeat_{c}"]
    for t in TIERS:
        sets[f"in_{t}"] = rb.tier_intervals(t)
    return {k: cover_cells(v, cells) for k, v in sets.items()}, genes.get("tss", [])


def analyse(
    chrom: str,
    cache: Path = CACHE,
    results_dir: Path | None = None,
    phylop: int = 0,
    tracks_dir: Path = TRACKS,
    constraint_dir: Path = CONSTRAINT_CACHE,
) -> dict[str, Any]:
    """Measure one chromosome's background composition and decompose its offset."""
    rb = rebuild(chrom, cache, results_dir)
    tracks = local_tracks(chrom, tracks_dir, constraint_dir)
    genes = gene_features(chrom)
    feats, tss = feature_arrays(rb, tracks, genes)
    names = [k for k in feats if k != "cds_canonical"]
    absent = _unassessed(tracks, genes)
    rows, dropped = background_bins(rb, feats, tss)
    outside = [r for r in rows if not any(r[f"in_{t}"] for t in TIERS)]
    comp: dict[str, Any] = {
        "background": composition(rows, names, absent),
        "background_outside_every_block": composition(outside, names, absent),
    }
    for g in GROUPS:
        comp[g] = block_composition(rb, g, feats, tss, names, absent)
    offset: dict[str, Any] = {"as_the_panel_builds_it": tier_ratios(rb, rb.cds_all)}
    for name in CANDIDATES:
        ivs = _ivs(tracks, genes, name)
        if not ivs:
            offset[f"without_{name}"] = {"assessed": False, "why": _why(tracks, genes, name)}
            continue
        offset[f"without_{name}"] = {"assessed": True, **tier_ratios(rb, merge_intervals(rb.cds_all + ivs))}
    together = [n for n in TOGETHER if _ivs(tracks, genes, n)]
    every = merge_intervals(rb.cds_all + [iv for n in together for iv in _ivs(tracks, genes, n)])
    offset["without_any_of_them"] = {
        "assessed": bool(together),
        "covariates_included": together,
        "covariates_left_out": [n for n in TOGETHER if n not in together],
        **tier_ratios(rb, every),
    }
    thin = [(r["start"], r["end"]) for r in rows if (r["coverage"] or 0) < COVERAGE_MIN]
    enough = rows and len(rows) - len(thin) >= COVERAGE_KEEP_MIN * len(rows)
    offset[COVERAGE_KEY] = {
        "assessed": bool(enough),
        "is_coverage_not_a_covariate": True,
        "measure": COVERAGE_MEASURE,
        "threshold": COVERAGE_MIN,
        "kilobases": len(rows),
        "kilobases_left_out": len(thin),
        **(tier_ratios(rb, merge_intervals(rb.cds_all + thin)) if enough else {"why": COVERAGE_TOO_THIN}),
    }
    win, win_cost = window_rows(rb, tss)
    stand = standardised_against_background(rb, win, rb.cds_all)
    held = [stand[t]["gc_timing_tss_and_coverage"]["targets_matched"] for t in TIERS if t in stand]
    return {
        "chrom": chrom,
        "assemblies": rb.panel.n,
        "evidence": SOURCES,
        # the window comparison holds coverage on every chromosome; the exclusion arm is the coarser
        # second way of doing it and is allowed to be unassessed without weakening the sentence
        "claim_available": CLAIM_WITH_COVERAGE if held and all(held) else CLAIM_WITHOUT_COVERAGE,
        "coverage_conditioned": {
            "in_the_window_comparison": bool(held and all(held)),
            "by_rebuilding_the_background": bool(enough),
        },
        "reproduces_the_panel": reproduces_the_panel(rb, offset["as_the_panel_builds_it"]),
        "composition": comp,
        "feature_splits": feature_splits(rows, names, absent),
        "offset": offset,
        "explained": {t: explained(offset, t) for t in TIERS},
        "constrained_bases": {
            "background_sampled": phylop_sample(chrom, rows, phylop),
            "tiers_from_the_budget": tier_constraint(rb),
        },
        "standardised": stand,
        "coverage": {
            "kilobases_in_panel_blocks": len(rows) + sum(dropped.values()),
            "kilobases_measured": len(rows),
            "kilobases_not_measured": dropped,
            # measured, and kept by the panel, but matched on GC alone because no timing reached them
            "kilobases_measured_without_replication_timing": sum(1 for r in rows if r["rt_stratum"] is None),
            "covariates_assessed": [n for n in CANDIDATES if _ivs(tracks, genes, n)],
            "covariates_not_assessed": {
                n: _why(tracks, genes, n) for n in CANDIDATES if not _ivs(tracks, genes, n)
            },
            "tracks_absent": tracks["missing"],
            "gencode_read": bool(genes),
            "coding_tss_read": len(tss),
            "windows": win_cost,
        },
    }


def _unassessed(tracks: dict, genes: dict) -> list[str]:
    """Feature names with no data on this chromosome: they get a null share, never a zero one."""
    out = list(tracks["missing"])
    if "repeats" in tracks["missing"]:
        out += [f"repeat_{c}" for c in REPEAT_CLASSES]
    if not genes:
        out += ["cds_any", "exon_any", "promoter"]
    return out


def _why(tracks: dict, genes: dict, name: str) -> str:
    if name in ("cds_any", "exon_any", "promoter") and not genes:
        return "GENCODE not read for this chromosome"
    if name in tracks["missing"]:
        return "no track for this chromosome"
    return "read, and empty on this chromosome"


def explained(offset: dict[str, Any], tier: str) -> dict[str, Any]:
    """How far towards 1 each exclusion moves a tier's ratio, as a share of its distance from 1.

    Every entry says whether the covariate was assessed. A covariate that was measured and moves
    nothing reads `{"assessed": true, "explains": 0.0}`; one that could not be measured reads
    `{"assessed": false, "explains": null}` with the reason. The two are not interchangeable.
    """
    base = offset["as_the_panel_builds_it"].get(tier)
    out: dict[str, Any] = {"ratio_as_built": base}
    if base is None or base == 1:
        out["why"] = "the tier is absent, or already exactly at the background"
        return out
    for key, row in offset.items():
        if key == "as_the_panel_builds_it" or not isinstance(row, dict):
            continue
        if not row.get("assessed"):
            out[key] = {"assessed": False, "explains": None, "reads": "not assessed", "why": row.get("why")}
            continue
        r = row.get(tier)
        share = round((base - r) / (base - 1), 4) if r is not None else None
        out[key] = {
            "assessed": True,
            "ratio": r,
            "explains": share,
            "reads": reads_explained(share),
        }
    return out


def reads_explained(share: float | None, floor: float = 0.02) -> str:
    """The decomposition cell in words. Accounting for none is a sentence, not a zero."""
    if share is None:
        return "not assessed"
    if abs(share) < floor:
        return "measured, accounts for none"
    if share < 0:
        return f"measured, widens the offset by {round(-100 * share)}%"
    return f"measured, accounts for {round(100 * share)}% of the offset"
