# SPDX-License-Identifier: AGPL-3.0-or-later
"""What a reporter assay could not measure in the 13.77 Mb nothing has measured yet.

`scripts/unknown_coverage.py` established on 2026-09-17 that 0.45% of the unknown space has ever
been touched by an assay this project holds, and that the untouched part of the real unknown --
680 blocks, 13.77 Mb -- tiles at 300 bp into about 45,900 oligos. That number is an upper bound
taken from the arithmetic 13,770,000 / 300. This module measures how much of it a library could
actually order, and refuses to predict anything.

**What counts as "cannot be measured in principle", fixed here before any block was read.** Two
families, because they fail for different reasons and a library treats them differently:

*Family A, not synthesisable or not resolvable* -- the oligo cannot exist, or cannot be told apart
from another oligo of the same library:

1. `assembly_gap`: more than `N_MAX` of the window is N. Those bases are not in the reference; there
   is no sequence to order.
2. `gc_extreme`: GC outside [`GC_LOW`, `GC_HIGH`]. Array synthesis and the emulsion PCR that
   amplifies the pool both drop out at the extremes, so the oligo is absent or under-represented
   rather than measured.
3. `homopolymer`: a run of one base at least `HOMOPOLYMER` long. Phosphoramidite coupling slips on
   long runs and the read-back cannot count the run's length, so the member's identity is uncertain.
4. `tandem_low_complexity`: at least `TANDEM_MAX` of the window is RepeatMasker Simple_repeat,
   Low_complexity or Satellite. Tandem arrays are the classic synthesis and assembly failure, and a
   satellite window is also unplaceable.
5. `non_unique_in_block`: at least `DUP_KMER_MAX` of the window's `KMER`-mers occur more than once
   inside its own block, so two windows of the same tiling are the same molecule and the reads
   cannot be assigned to one of them.

*Family B, synthesisable but not attributable to one locus* -- the assay would return a number, and
the number would be a property of a sequence rather than of a place in the genome:

6. `segmental_duplication`: at least `SEGDUP_MAX` of the window lies in a curated segmental
   duplication (UCSC genomicSuperDups, >= 1 kb at >= 90% identity). The same sequence exists
   elsewhere, so a positive cannot be assigned to this locus and no CRISPRi or deletion follow-up
   can target it.
7. `interspersed_repeat`: at least `INTERSPERSED_MAX` of the window is LINE, SINE, LTR, DNA,
   Retroposon or RC. The window is one of thousands of near-identical genomic copies.

Two more categories are arithmetic, not sequence, and are named so that the base accounting closes:
`shorter_than_one_oligo` (a block under `OLIGO` bp yields no window at all) and `tiling_remainder`
(the tail of a block that no window covers at `STEP`).

Every threshold above is a named constant with its reason, and each arbitrary one is reported at two
or three values (`SENSITIVITY`) instead of being chosen. Family B is reported beside Family A, never
folded into it, because excluding a repeat-derived window is a decision about attribution and a
library may legitimately take it either way.

Nothing here is a model of function and nothing is extrapolated: the module reads sequence, the
RepeatMasker track and the duplication track, and counts.
"""

from __future__ import annotations

import bisect
import gzip
import random
from collections.abc import Iterable, Sequence
from pathlib import Path
from statistics import median
from typing import Any

from genomeos.genome.duplications import YOUNG, covered, merge
from genomeos.results import RESULTS_DIR

# --- the library's geometry, taken from the design this measures -------------------------------
OLIGO = 300  # the peer design's oligo length (scripts/unknown_library.py); 13.77 Mb / 300 = 45,900
STEP = 300  # non-overlapping, the tiling whose arithmetic produced 45,900; overlap is stated, not assumed

# --- Family A thresholds ----------------------------------------------------------------------
N_MAX = 0.10  # share of a window that may be N; above it there is no sequence to order
GC_LOW = 0.25  # below this an oligo drops out of array synthesis and pool PCR
GC_HIGH = 0.75  # above this the same, plus secondary structure the polymerase stalls on
HOMOPOLYMER = 10  # run of one base from which coupling slips and read-back miscounts the run
TANDEM_MAX = 0.50  # share of a window in a tandem/low-complexity/satellite repeat
KMER = 40  # k for within-block uniqueness: longer than any motif, shorter than an oligo
DUP_KMER_MAX = 0.50  # share of a window's k-mers that recur in its own block

# --- Family B thresholds ----------------------------------------------------------------------
SEGDUP_MAX = 0.50  # share of a window inside a curated segmental duplication
INTERSPERSED_MAX = 0.50  # share of a window inside an interspersed repeat family

TANDEM_CLASSES = ("Simple_repeat", "Low_complexity", "Satellite")
INTERSPERSED_CLASSES = ("LINE", "SINE", "LTR", "DNA", "Retroposon", "RC")

FAMILY_A = ("assembly_gap", "gc_extreme", "homopolymer", "tandem_low_complexity", "non_unique_in_block")
FAMILY_B = ("segmental_duplication", "interspersed_repeat")
ARITHMETIC = ("shorter_than_one_oligo", "tiling_remainder")
REASONS = (*FAMILY_A, *FAMILY_B)
# The order a window is charged to exactly one reason. Physical impossibility first, then
# synthesis, then resolution, then attribution: a window that fails several is charged to the
# first, and the raw per-reason counts are reported beside the exclusive ones.
PRECEDENCE = REASONS

SENSITIVITY: dict[str, tuple] = {
    # each arbitrary threshold at the values the result reports instead of one chosen value
    "n_max": (0.0, 0.10, 0.50),
    "gc": ((0.25, 0.75), (0.30, 0.70), (0.20, 0.80)),
    "homopolymer": (8, 10, 12),
    "tandem_max": (0.25, 0.50, 0.75),
    "dup_kmer_max": (0.25, 0.50, 0.90),
    "segdup_max": (0.01, 0.50, 1.00),
    "interspersed_max": (0.25, 0.50, 0.75),
}

MPRA_CALIBRATION_MAX = 20_000  # lentiMPRA elements sampled for the false-exclusion check
SEED = 17


class Thresholds:
    """One complete set of cut-offs, so the same windows can be read at several of them."""

    __slots__ = (
        "dup_kmer_max",
        "gc_high",
        "gc_low",
        "homopolymer",
        "interspersed_max",
        "n_max",
        "segdup_max",
        "tandem_max",
    )

    def __init__(
        self,
        n_max: float = N_MAX,
        gc_low: float = GC_LOW,
        gc_high: float = GC_HIGH,
        homopolymer: int = HOMOPOLYMER,
        tandem_max: float = TANDEM_MAX,
        dup_kmer_max: float = DUP_KMER_MAX,
        segdup_max: float = SEGDUP_MAX,
        interspersed_max: float = INTERSPERSED_MAX,
    ) -> None:
        self.n_max = n_max
        self.gc_low = gc_low
        self.gc_high = gc_high
        self.homopolymer = homopolymer
        self.tandem_max = tandem_max
        self.dup_kmer_max = dup_kmer_max
        self.segdup_max = segdup_max
        self.interspersed_max = interspersed_max

    def as_dict(self) -> dict[str, Any]:
        return {k: getattr(self, k) for k in self.__slots__}


DEFAULT = Thresholds()


# --- sequence measurements --------------------------------------------------------------------


def n_fraction(seq: str) -> float:
    """Share of the window that is not A, C, G or T: assembly gap or ambiguity."""
    if not seq:
        return 1.0
    acgt = sum(seq.count(b) for b in "ACGT")
    return 1.0 - acgt / len(seq)


def gc_fraction(seq: str) -> float | None:
    """GC over the called bases only; None when nothing is called."""
    acgt = sum(seq.count(b) for b in "ACGT")
    if not acgt:
        return None
    return (seq.count("G") + seq.count("C")) / acgt


def longest_homopolymer(seq: str) -> int:
    """Length of the longest run of one base, N runs included (an N run is still a run)."""
    best = run = 0
    prev = ""
    for b in seq:
        run = run + 1 if b == prev else 1
        prev = b
        if run > best:
            best = run
    return best


def duplicated_kmer_fraction(block_seq: str, offset: int, length: int, k: int = KMER) -> float:
    """Share of a window's k-mers that occur more than once in the block it is tiled from.

    Counted over the whole block once by the caller would be cheaper; this signature keeps the
    function testable on its own, and `kmer_counts` is the shared path the chromosome pass uses.
    """
    counts = kmer_counts(block_seq, k)
    return duplicated_share(block_seq, counts, offset, length, k)


def kmer_counts(block_seq: str, k: int = KMER) -> dict[str, int]:
    counts: dict[str, int] = {}
    for i in range(len(block_seq) - k + 1):
        km = block_seq[i : i + k]
        counts[km] = counts.get(km, 0) + 1
    return counts


def duplicated_share(
    block_seq: str, counts: dict[str, int], offset: int, length: int, k: int = KMER
) -> float:
    """Share of the window's k-mers whose count in the block exceeds one."""
    window = block_seq[offset : offset + length]
    n = len(window) - k + 1
    if n <= 0:
        return 0.0
    dup = sum(1 for i in range(n) if counts.get(window[i : i + k], 0) > 1)
    return dup / n


# --- interval measurements --------------------------------------------------------------------


def fraction_covered(merged: list[tuple[int, int]], start: int, end: int) -> float:
    """Share of [start, end) inside already-merged intervals."""
    if end <= start:
        return 0.0
    return covered(merged, start, end) / (end - start)


def class_groups(by_class: dict[str, list[tuple[int, int]]]) -> dict[str, list[tuple[int, int]]]:
    """The two repeat groups the exclusions use, merged, plus everything else under `other_repeat`."""
    groups: dict[str, list[tuple[int, int]]] = {"tandem": [], "interspersed": [], "other_repeat": []}
    for cls, spans in by_class.items():
        key = (
            "tandem"
            if cls in TANDEM_CLASSES
            else "interspersed"
            if cls in INTERSPERSED_CLASSES
            else "other_repeat"
        )
        groups[key].extend(spans)
    return {k: merge(v) for k, v in groups.items()}


# --- one window -------------------------------------------------------------------------------


def window_metrics(
    seq: str,
    start: int,
    end: int,
    groups: dict[str, list[tuple[int, int]]],
    segdups: list[tuple[int, int]],
    young_segdups: list[tuple[int, int]],
    dup_kmer: float,
) -> dict[str, Any]:
    """Everything measured about one oligo window, with no cut-off applied yet."""
    return {
        "start": start,
        "end": end,
        "n_fraction": round(n_fraction(seq), 4),
        "gc": None if gc_fraction(seq) is None else round(gc_fraction(seq) or 0.0, 4),
        "homopolymer": longest_homopolymer(seq),
        "tandem_fraction": round(fraction_covered(groups.get("tandem", []), start, end), 4),
        "interspersed_fraction": round(fraction_covered(groups.get("interspersed", []), start, end), 4),
        "segdup_fraction": round(fraction_covered(segdups, start, end), 4),
        "young_segdup_fraction": round(fraction_covered(young_segdups, start, end), 4),
        "dup_kmer_fraction": round(dup_kmer, 4),
    }


def failed_reasons(m: dict[str, Any], t: Thresholds = DEFAULT) -> list[str]:
    """Every reason a window fails, in precedence order; empty when it is usable."""
    out = []
    if m["n_fraction"] > t.n_max:
        out.append("assembly_gap")
    gc = m["gc"]
    if gc is None or gc < t.gc_low or gc > t.gc_high:
        out.append("gc_extreme")
    if m["homopolymer"] >= t.homopolymer:
        out.append("homopolymer")
    if m["tandem_fraction"] >= t.tandem_max:
        out.append("tandem_low_complexity")
    if m["dup_kmer_fraction"] >= t.dup_kmer_max:
        out.append("non_unique_in_block")
    if m["segdup_fraction"] >= t.segdup_max:
        out.append("segmental_duplication")
    if m["interspersed_fraction"] >= t.interspersed_max:
        out.append("interspersed_repeat")
    return out


def charged_reason(m: dict[str, Any], t: Thresholds = DEFAULT) -> str:
    """The single reason a window is charged to, or `usable`."""
    reasons = failed_reasons(m, t)
    return reasons[0] if reasons else "usable"


def verdict(m: dict[str, Any], t: Thresholds = DEFAULT) -> str:
    """`usable`, `unattributable` (Family B only) or `unsynthesisable` (any Family A reason)."""
    reasons = set(failed_reasons(m, t))
    if reasons & set(FAMILY_A):
        return "unsynthesisable"
    if reasons & set(FAMILY_B):
        return "unattributable"
    return "usable"


# --- one block --------------------------------------------------------------------------------


def tiling(start: int, end: int, oligo: int = OLIGO, step: int = STEP) -> list[tuple[int, int]]:
    """Every window of the tiling that fits wholly inside [start, end)."""
    if end - start < oligo:
        return []
    return [(s, s + oligo) for s in range(start, end - oligo + 1, step)]


def assess_block(
    block: dict[str, Any],
    block_seq: str,
    groups: dict[str, list[tuple[int, int]]],
    segdups: list[tuple[int, int]],
    young_segdups: list[tuple[int, int]],
    tss_distance: int | None,
    oligo: int = OLIGO,
    step: int = STEP,
) -> dict[str, Any]:
    """One untouched block: its own composition, its windows' metrics and its tiling arithmetic."""
    lo, hi = block["start"], block["end"]
    seq = block_seq.upper()
    windows = tiling(lo, hi, oligo, step)
    counts = kmer_counts(seq, KMER) if windows else {}
    metrics = [
        window_metrics(
            seq[s - lo : e - lo],
            s,
            e,
            groups,
            segdups,
            young_segdups,
            duplicated_share(seq, counts, s - lo, e - s, KMER),
        )
        for s, e in windows
    ]
    tiled_bp = sum(min(e, hi) - s for s, e in windows) if windows else 0
    return {
        "chrom": block["chrom"],
        "block": f"{block['chrom']}:{lo}-{hi}",
        "case": block.get("case"),
        "length": hi - lo,
        "gc": None if gc_fraction(seq) is None else round(gc_fraction(seq) or 0.0, 4),
        "n_fraction": round(n_fraction(seq), 4),
        "tandem_fraction": round(fraction_covered(groups.get("tandem", []), lo, hi), 4),
        "interspersed_fraction": round(fraction_covered(groups.get("interspersed", []), lo, hi), 4),
        "other_repeat_fraction": round(fraction_covered(groups.get("other_repeat", []), lo, hi), 4),
        "repeat_fraction": round(
            fraction_covered(
                merge(
                    groups.get("tandem", []) + groups.get("interspersed", []) + groups.get("other_repeat", [])
                ),
                lo,
                hi,
            ),
            4,
        ),
        "segdup_fraction": round(fraction_covered(segdups, lo, hi), 4),
        "young_segdup_fraction": round(fraction_covered(young_segdups, lo, hi), 4),
        "tss_distance": tss_distance,
        "oligos": len(windows),
        "tiled_bp": tiled_bp,
        "remainder_bp": (hi - lo) - tiled_bp,
        "windows": metrics,
    }


def block_counts(row: dict[str, Any], t: Thresholds = DEFAULT) -> dict[str, Any]:
    """Per-block tallies at one threshold set: usable, attributable and the charged reasons."""
    charged: dict[str, int] = dict.fromkeys(REASONS, 0)
    usable = attributable = 0
    for m in row["windows"]:
        r = charged_reason(m, t)
        if r == "usable":
            usable += 1
            attributable += 1
        else:
            charged[r] += 1
            if r in FAMILY_B:
                usable += 1  # Family B is synthesisable; it fails attribution, not synthesis
    return {
        "oligos": row["oligos"],
        "synthesisable": usable,
        "attributable": attributable,
        "charged": charged,
    }


# --- track and sequence access ----------------------------------------------------------------


def repeat_spans_by_class(
    chrom: str, regions: list[tuple[int, int]], results_dir: Path = RESULTS_DIR
) -> dict[str, list[tuple[int, int]]] | None:
    """RepeatMasker intervals overlapping `regions`, by class, streamed from the local BED.

    Returns None when the chromosome's track is not cached, so the caller can count the block as
    unassessed rather than silently repeat-free.
    """
    p = results_dir / f"rmsk_{chrom}.bed.gz"
    if not p.exists():
        return None
    wanted = merge(regions)
    if not wanted:
        return {}
    lo, hi = wanted[0][0], wanted[-1][1]
    starts = [x[0] for x in wanted]
    out: dict[str, list[tuple[int, int]]] = {}
    with gzip.open(p, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            s, e, cls = int(parts[0]), int(parts[1]), parts[2]
            if e <= lo or s >= hi:
                continue
            if not _overlaps(wanted, starts, s, e):
                continue
            out.setdefault(cls, []).append((s, e))
    return {k: merge(v) for k, v in out.items()}


def _overlaps(merged: list[tuple[int, int]], starts: list[int], s: int, e: int) -> bool:
    """Whether [s, e) touches any of the merged intervals whose starts are `starts`."""
    i = bisect.bisect_right(starts, s) - 1
    return any(0 <= j < len(merged) and merged[j][0] < e and s < merged[j][1] for j in (i, i + 1))


def segdup_spans(chrom: str, results_dir: Path = RESULTS_DIR) -> tuple[list, list] | None:
    """(all curated duplications, those at >= YOUNG identity), merged; None when not cached."""
    from genomeos.genome import duplications

    if not duplications.rows_path(chrom, results_dir).exists():
        return None
    rows = duplications.load_rows(chrom, results_dir)
    return merge([(d.start, d.end) for d in rows]), merge(
        [(d.start, d.end) for d in rows if d.frac_match >= YOUNG]
    )


CRISPRI_FILES = (
    "EPCrisprBenchmark_combined_data.training_K562.GRCh38.tsv.gz",
    "EPCrisprBenchmark_combined_data.heldout_5_cell_types.GRCh38.tsv.gz",
)


def assay_spans(chrom: str) -> dict[str, list[tuple[int, int]]]:
    """The three measured layers on one chromosome, merged, by the assay that measured them.

    The same three loaders and the same one-base rule as `scripts/unknown_coverage.py`, restated
    here so this module does not depend on a script: the count of untouched blocks it produces is
    printed against that reading's 680 blocks / 13.77 Mb as a reproduction check.
    """
    from genomeos.attribution import crispri, mpra, vista

    spans: dict[str, list[tuple[int, int]]] = {}
    try:
        spans["lentimpra"] = merge([(e.start, e.end) for e in mpra.load(chrom)])
    except Exception:  # noqa: BLE001 - an unfetched layer is reported absent, not fatal
        spans["lentimpra"] = []
    try:
        spans["vista"] = merge([(e.start, e.end) for e in vista.load_loci(chrom)])
    except Exception:  # noqa: BLE001
        spans["vista"] = []
    rows: list[tuple[int, int]] = []
    for name in CRISPRI_FILES:
        try:
            rows.extend((p.start, p.end) for p in crispri.load(name) if p.chrom == chrom)
        except Exception:  # noqa: BLE001
            continue
    spans["crispri"] = merge(rows)
    return spans


def real_unknown_label(block: dict[str, Any]) -> str | None:
    """`real_unknown_<case>` for a constrained_unknown block that is not a copy, else None."""
    if block.get("tier") != "constrained_unknown" or block.get("copy"):
        return None
    return f"real_unknown_{block.get('case') or 'uncased'}"


def untouched_blocks(chrom: str) -> list[dict[str, Any]]:
    """The real-unknown blocks of one chromosome that no assay overlaps by a single base."""
    from genomeos.attribution import organise

    spans = assay_spans(chrom)
    out = []
    for b in organise.blocks(chrom):
        label = real_unknown_label(b)
        if label is None:
            continue
        if any(covered(v, b["start"], b["end"]) for v in spans.values()):
            continue
        out.append({**b, "chrom": chrom, "label": label})
    return out


def nearest(sorted_points: Sequence[int], pos: int) -> int | None:
    """Distance from `pos` to the closest of the sorted points; None when there are none."""
    if not sorted_points:
        return None
    i = bisect.bisect_left(sorted_points, pos)
    return min((abs(sorted_points[j] - pos) for j in (i - 1, i) if 0 <= j < len(sorted_points)), default=None)


# --- pooled summaries -------------------------------------------------------------------------


def _median(values: Iterable[float]) -> float | None:
    vals = [v for v in values if v is not None]
    return round(median(vals), 4) if vals else None


def covariates(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The three covariates every group in this project must print: length, GC, TSS distance."""
    lengths = [r["length"] for r in rows]
    return {
        "n": len(rows),
        "median_length": _median(lengths),
        "length_quartiles": _quartiles(lengths),
        "median_gc": _median(r["gc"] for r in rows),
        "median_tss_distance": _median(r["tss_distance"] for r in rows),
        "tss_distance_quartiles": _quartiles([r["tss_distance"] for r in rows if r["tss_distance"]]),
        "without_tss_distance": sum(1 for r in rows if r["tss_distance"] is None),
    }


def _quartiles(values: list[float]) -> list[float] | None:
    """[q25, median, q75] by nearest rank, truncating: `v[int(q * (n - 1))]`."""
    v = sorted(x for x in values if x is not None)
    if not v:
        return None
    return [round(v[int(q * (len(v) - 1))], 2) for q in (0.25, 0.5, 0.75)]


def window_covariates(rows: list[dict[str, Any]], t: Thresholds = DEFAULT) -> dict[str, Any]:
    """Per charged reason: how many windows, their GC, and the blocks and TSS distances behind them.

    Window length is `OLIGO` by construction, so the length a group is reported at is the length of
    the blocks its windows came from.
    """
    groups: dict[str, list[tuple[dict, dict]]] = {}
    for r in rows:
        for m in r["windows"]:
            groups.setdefault(charged_reason(m, t), []).append((r, m))
    out = {}
    for name in ("usable", *REASONS):
        pairs = groups.get(name, [])
        out[name] = {
            "windows": len(pairs),
            "bp": len(pairs) * OLIGO,
            "median_gc": _median(m["gc"] for _, m in pairs),
            "median_block_length": _median(r["length"] for r, _ in pairs),
            "median_tss_distance": _median(r["tss_distance"] for r, _ in pairs),
            "blocks": len({r["block"] for r, _ in pairs}),
        }
    return out


def pooled(rows: list[dict[str, Any]], t: Thresholds = DEFAULT, min_usable: int = 3) -> dict[str, Any]:
    """The whole set at one threshold set, with the base accounting closed to the last base."""
    charged: dict[str, int] = dict.fromkeys(REASONS, 0)
    raw: dict[str, int] = dict.fromkeys(REASONS, 0)
    oligos = synthesisable = attributable = 0
    short_bp = remainder_bp = 0
    below: list[str] = []
    for r in rows:
        c = block_counts(r, t)
        oligos += c["oligos"]
        synthesisable += c["synthesisable"]
        attributable += c["attributable"]
        for k, v in c["charged"].items():
            charged[k] += v
        for m in r["windows"]:
            for name in failed_reasons(m, t):
                raw[name] += 1
        if r["oligos"] == 0:
            short_bp += r["length"]
        else:
            remainder_bp += r["remainder_bp"]
        if c["attributable"] < min_usable:
            below.append(r["block"])
    total_bp = sum(r["length"] for r in rows)
    excluded_bp = {name: charged[name] * OLIGO for name in REASONS}
    excluded_bp["shorter_than_one_oligo"] = short_bp
    excluded_bp["tiling_remainder"] = remainder_bp
    accounted = attributable * OLIGO + sum(excluded_bp.values())
    return {
        "blocks": len(rows),
        "bp": total_bp,
        "oligo": OLIGO,
        "step": STEP,
        "overlap": OLIGO - STEP,
        "thresholds": t.as_dict(),
        "oligos_tiled": oligos,
        "oligos_synthesisable": synthesisable,
        "oligos_attributable": attributable,
        "charged_windows": charged,
        "reasons_raw_windows": raw,
        "excluded_bp": excluded_bp,
        "attributable_bp": attributable * OLIGO,
        "bp_accounted": accounted,
        "bp_unaccounted": total_bp - accounted,
        "blocks_below_min_usable": len(below),
        "min_usable_oligos": min_usable,
        "blocks_below_examples": below[:20],
    }


def sensitivity(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Every arbitrary threshold moved one at a time, so no number rests on a chosen cut-off."""
    out = []
    for key, values in SENSITIVITY.items():
        for v in values:
            kw: dict[str, Any] = {}
            if key == "gc":
                kw = {"gc_low": v[0], "gc_high": v[1]}
                label = f"gc {v[0]}-{v[1]}"
            else:
                kw = {key: v}
                label = f"{key} {v}"
            p = pooled(rows, Thresholds(**kw))
            out.append(
                {
                    "threshold": label,
                    "changed": key,
                    "value": list(v) if isinstance(v, tuple) else v,
                    "oligos_synthesisable": p["oligos_synthesisable"],
                    "oligos_attributable": p["oligos_attributable"],
                    "blocks_below_min_usable": p["blocks_below_min_usable"],
                }
            )
    return out


def calibration(mpra_rows: list[dict[str, Any]], t: Thresholds = DEFAULT) -> dict[str, Any]:
    """The same Family A rules over sequence a reporter assay has already measured successfully.

    lentiMPRA ordered these elements and got a number out of them, so a rule that excludes many of
    them is too strict. This is a plausibility check on the cut-offs, NOT a matched comparison: the
    lentiMPRA elements are ENCODE cCREs and differ from the untouched blocks in GC and in distance
    to a coding TSS, which is why both are printed.
    """
    charged: dict[str, int] = dict.fromkeys(REASONS, 0)
    charged["usable"] = 0
    for r in mpra_rows:
        charged[charged_reason(r["window"], t)] += 1
    n = len(mpra_rows)
    a = sum(charged[k] for k in FAMILY_A)
    b = sum(charged[k] for k in FAMILY_B)
    return {
        "windows": n,
        "excluded_family_a": a,
        "false_exclusion_rate": round(a / n, 4) if n else None,
        "flagged_family_b": b,
        "family_b_rate": round(b / n, 4) if n else None,
        "not_calibrated": ["non_unique_in_block"],  # no block to be non-unique in; 0 by construction
        "charged": charged,
        "covariates": {
            "median_gc": _median(r["window"]["gc"] for r in mpra_rows),
            "median_tss_distance": _median(r["tss_distance"] for r in mpra_rows),
            "median_length": OLIGO,
        },
        "reading": (
            "measured elements the same rules would have thrown away; a high rate means the cut-offs "
            "are too strict, not that the elements are unmeasurable"
        ),
    }


def sample(items: list, cap: int, seed: int = SEED) -> list:
    if len(items) <= cap:
        return items
    return random.Random(seed).sample(items, cap)
