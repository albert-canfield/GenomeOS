# SPDX-License-Identifier: AGPL-3.0-or-later
"""Milestone 1.3, clause 2, item 9b (M5b): a fine-mapped eQTL arm, calibrated outside the 882 blocks.

M5 (`measured_unknown`, 7906142) found that public element perturbation reaches 13 of the 882
real-unknown blocks against a coding gene, 62 short of the smallest calibrated design, and that the one
source whose tested universe reaches the sizing is GTEx v8 fine-mapped cis-eQTL (DAP-G): 879 blocks hold
a common variant GTEx could have tested. Nobody had measured that arm's window rate, its clustering or
its sensitivity, and the CRISPRi sizing (lane-s2, `cab70d9`) describes a gene-level perturbation readout
and nothing else. This module measures them, on sequence OUTSIDE the blocks only, and sizes the arm.

THE EXCLUSION COMES FIRST. `Exclusion` is the 882 blocks widened by `EXCLUSION_MARGIN` on each side.
`DapgTrack` decodes the GTEx DAP-G bigBed record by record and cuts the variant's position out of a
record BEFORE any other field is split: a record whose variant lies inside the exclusion is dropped
there, its gene, tissue and PIP never parsed, never counted and never stored. No calibration unit or
window may touch the exclusion either (`Exclusion.hits`), so nothing in a calibration can be a value
from inside a block. A count of the dropped records is deliberately not kept: how many fine-mapped rows
fall in the blocks is itself a fact about the blocks' eQTLs.

WHAT IS CALIBRATED. The CRISPRi benchmark (`crispri.TRAINING`, `crispri.HELDOUT`) says which element
regulates which gene by measurement. On its positive pairs (`Regulated` TRUE) and on its tested
negatives, at two scales (the element with the 500-base margin the earlier DAP-G lanes used, and a
window of the blocks' median length centred on the element), this module measures:

- sensitivity: how often a fine-mapped variant (DAP-G PIP >= `PIP`) inside a unit names the unit's own
  measured gene, in any GTEx tissue and in the tissue matched to the screen's cell;
- the window rate: the same at `WINDOWS` (four) windows per unit, plus a fifth drawn as a pseudo-block,
  matched on length, on the count of dbSNP 155 common variants inside (the proxy M5 used for the
  variants GTEx tested), and on distance to the gene's TSS, drawn outside the exclusion, the TSS and
  every positive element of the same gene;
- the naming rate: of the units holding a fine-mapped variant for any gene, the share naming the unit's
  own gene;
- clustering: the within-unit intraclass correlation of the windows' hits, the correlation of a unit's
  hit with its windows' mean, and the hit rate of the two adjacent same-length flanks given the unit's
  own hit (a fine-mapped signal spills across neighbouring sequence through linkage);
- the sizing: by simulation, resampling the calibrated units, how many blocks each compared against its
  own `WINDOWS` windows a test needs for 80% probability of detection at one-sided alpha 0.05, when a
  share f of the blocks' attributions is right, stated as the CRISPRi sizing was.

An eQTL is an association between a variant and a gene's expression across people. It is not a
perturbation, and nothing here is mechanism. Fine-mapping localises an association to a variant with a
posterior; it does not say which element acts.

FAIR_TEST_RULE is fixed here, in the committed code, before the calibration is run: it decides from the
sizing alone whether a test inside the blocks may be registered.
"""

from __future__ import annotations

import bisect
import hashlib
import json
import math
import random
import re
import struct
import zlib
from collections import defaultdict
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from genomeos.attribution.bigwig import coalesce
from genomeos.attribution.executor import DAPG_PIP
from genomeos.attribution.human_panel import GBDB, BigBed

STORE = Path("data/knowledge/eqtl_arm")
DAPG_URL = f"{GBDB}/gtex/eQtl/gtexDapg.bb"
COMMON_URL = f"{GBDB}/snp/dbSnp155Common.bb"
DAPG_FILE = STORE / "gtexDapg.bb"
COMMON_FILE = STORE / "dbSnp155Common.bb"
COMMON_CACHE = STORE / "common_positions"

#: fine-mapped: `executor.DAPG_PIP`, the threshold every DAP-G lane has used
PIP = DAPG_PIP
#: the element margin of every earlier DAP-G read (scripts/eqtl_targets.MARGIN, through
#: scripts/response_map_increment3_count.EQTL_MARGIN); a test checks the two agree
MARGIN = 500
#: the blocks are widened by this on each side before anything is decoded or drawn near them
EXCLUSION_MARGIN = 1_000
#: GTEx v8 tested variants within 1 Mb of a gene's TSS
CIS_WINDOW = 1_000_000
#: a window may not hold the gene's TSS, nor come this close to it (the blocks lie between gene bodies)
TSS_GUARD = 1_000
#: matched windows per unit used as its control, plus one more drawn as a pseudo-block for the null.
#: Four, not ten: a design probe that read coordinates and common-variant counts only (no eQTL value)
#: found that ten disjoint windows within a distance band fit around under half of the positives at the
#: element scale and about a tenth at the block scale, the losses concentrated near the TSS
WINDOWS = 4
#: a window's distance to the TSS lies within this share of the unit's, the band being at least
#: `MIN_BAND_LENGTHS` unit lengths wide on each side so that the windows fit
DISTANCE_BAND = 0.25
MIN_BAND_LENGTHS = 3
#: a window's common-variant count lies within this share of the unit's (at least one variant)
COUNT_BAND = 0.25
#: draws allowed for each window before the unit is called unmatched
DRAW_TRIES = 1_000
SEED = 20261010
#: the screen's cell to the GTEx tissue whose fine-mapping is matched to it; a judgement, stated
TISSUE_OF_CELL: dict[str, str | None] = {
    "K562": "Whole_Blood",
    "Jurkat": "Whole_Blood",
    "GM12878": "Cells_EBV-transformed_lymphocytes",
    "HCT116": "Colon_Transverse",
    "WTC11": None,
}
#: the sizing
ALPHA = 0.05
POWER = 0.80
SHARES = (1.0, 0.75, 0.5, 0.25, 0.1)
N_GRID = (25, 50, 75, 100, 150, 200, 300, 400, 531, 750, 882, 1000, 1500, 2000, 3000, 5000)
EXPERIMENTS = 2_000
BOOTSTRAPS = 2_000

#: the most a test may reject when no attribution is right but the blocks behave as tested
#: non-regulating elements: twice the nominal level
FOOLED_BOUND = 2 * ALPHA
FAIR_TEST_RULE = (
    "Fixed in the committed code before the calibration is run. A test inside the blocks is registered "
    "only if, at the calibrated rates and at the scale the test would use, BOTH hold for the blocks "
    "available to it (those whose attribution names a gene whose TSS lies on the same chromosome within "
    "the 1 Mb cis window and outside the unit): (1) the conditional test reaches 80% probability of "
    "detection at one-sided alpha 0.05 when half of the attributions are right and the rest behave as "
    "matched windows (f = 1/2); and (2) when none is right but the blocks behave as CRISPRi-tested "
    "elements that do not regulate the named gene, it rejects at most twice the nominal level (0.10), "
    "so the matched-window control is not fooled by sequence of the kind a screen tests. If neither "
    "scale meets both, no test is registered and the reason is recorded with the numbers. The rule "
    "reads the sizing table and the count of available blocks, nothing else."
)


# ---- the exclusion -------------------------------------------------------------------------------


class Exclusion:
    """The blocks, widened by `margin` on each side and merged: nothing inside may be decoded or drawn."""

    def __init__(self, rows: Iterable[Mapping[str, Any]], margin: int = EXCLUSION_MARGIN) -> None:
        by: dict[str, list[tuple[int, int]]] = defaultdict(list)
        n = 0
        for r in rows:
            by[r["chrom"]].append((max(0, int(r["start"]) - margin), int(r["end"]) + margin))
            n += 1
        self.margin = margin
        self.blocks = n
        self.spans: dict[str, list[tuple[int, int]]] = {}
        for c, ivs in by.items():
            merged: list[list[int]] = []
            for s, e in sorted(ivs):
                if merged and s <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], e)
                else:
                    merged.append([s, e])
            self.spans[c] = [(s, e) for s, e in merged]
        self._starts = {c: [s for s, _ in v] for c, v in self.spans.items()}

    def holds(self, chrom: str, pos: int) -> bool:
        """Whether the base `pos` lies inside the widened blocks."""
        starts = self._starts.get(chrom)
        if not starts:
            return False
        i = bisect.bisect_right(starts, pos) - 1
        return i >= 0 and pos < self.spans[chrom][i][1]

    def hits(self, chrom: str, start: int, end: int) -> bool:
        """Whether [start, end) shares a base with the widened blocks."""
        starts = self._starts.get(chrom)
        if not starts or end <= start:
            return False
        i = bisect.bisect_left(starts, end) - 1
        return i >= 0 and self.spans[chrom][i][1] > start

    def bases(self) -> int:
        return sum(e - s for v in self.spans.values() for s, e in v)

    def digest(self) -> str:
        text = json.dumps({"margin": self.margin, "spans": self.spans}, sort_keys=True)
        return hashlib.sha256(text.encode()).hexdigest()


# ---- the readers: position first, everything else after the exclusion ----------------------------


def eqtl_pos(text: str) -> tuple[str, int] | None:
    """A DAP-G record's `eqtlPos` ('chr1:12345-12345' or 'chr1:12345', 1-based) as (chrom, 0-based
    base), the reading `human_panel._eqtl_pos` makes; None if malformed."""
    m = re.fullmatch(r"(\w+):(\d+)(?:-\d+)?", text.strip())
    return (m.group(1), int(m.group(2)) - 1) if m else None


@dataclass
class Fine:
    """One chromosome's DAP-G records outside the exclusion, as sorted arrays."""

    pos: np.ndarray  # int64, the variant's 0-based base
    gene: np.ndarray  # int32 code into `genes`
    tissue: np.ndarray  # int16 code into `tissues`
    pip: np.ndarray  # float32
    genes: list[tuple[str, str]]  # (unversioned Ensembl id, name as the track gives it)
    tissues: list[str]

    def slice(self, start: int, end: int) -> slice:
        return slice(
            int(np.searchsorted(self.pos, start, "left")), int(np.searchsorted(self.pos, end, "left"))
        )


DAPG_FIELDS = (
    "chrom",
    "chromStart",
    "chromEnd",
    "name",
    "score",
    "strand",
    "thickStart",
    "thickEnd",
    "reserved",
    "blockCount",
    "blockSizes",
    "chromStarts",
    "eqtlPos",
    "eqtlName",
    "genePos",
    "geneName",
    "geneId",
    "geneStrand",
    "tissue",
    "pip",
    "clusterPip",
)


def decode_records(
    raw: bytes, chrom_id: int, chrom: str, exclusion: Exclusion, pos_field: int, fields: Sequence[str]
) -> Iterable[tuple[int, dict[str, str]]]:
    """The records of one uncompressed bigBed data block whose variant lies OUTSIDE the exclusion.

    The order is the guarantee: the record's text is split only as far as the position field, the
    position is checked, and only a record that passes has the rest of its fields split. `pos_field` is
    the position's index among the fields after the first three."""
    p = 0
    n = len(raw)
    while p + 12 <= n:
        c, _s, _e = struct.unpack_from("<III", raw, p)
        z = raw.index(b"\0", p + 12)
        rest = raw[p + 12 : z]
        p = z + 1
        if c != chrom_id:
            continue
        head = rest.split(b"\t", pos_field + 1)
        if len(head) <= pos_field:
            continue
        got = eqtl_pos(head[pos_field].decode("ascii", "replace"))
        if got is None or got[0] != chrom or exclusion.holds(chrom, got[1]):
            continue  # dropped on its position alone; nothing else of it is split or kept
        vals = rest.decode("utf-8", "replace").split("\t")
        yield got[1], dict(zip(fields[3:], vals, strict=False))


class DapgTrack:
    """The GTEx v8 DAP-G bigBed, held locally, read one chromosome at a time behind the exclusion."""

    def __init__(self, exclusion: Exclusion, path: Path = DAPG_FILE) -> None:
        self.bb = BigBed(path)
        if tuple(self.bb.fields) != DAPG_FIELDS:
            raise ValueError(f"unexpected DAP-G fields: {self.bb.fields}")
        self.exclusion = exclusion
        self.pos_field = DAPG_FIELDS.index("eqtlPos") - 3

    def chromosome(self, chrom: str) -> Fine:
        bb = self.bb
        empty = Fine(
            np.zeros(0, np.int64),
            np.zeros(0, np.int32),
            np.zeros(0, np.int16),
            np.zeros(0, np.float32),
            [],
            [],
        )
        if chrom not in bb.chroms:
            return empty
        cid, size = bb.chroms[chrom]
        gene_code: dict[tuple[str, str], int] = {}
        tissue_code: dict[str, int] = {}
        rows: list[tuple[int, int, int, float]] = []
        for group in coalesce(bb.leaf_items(chrom, [(0, size)]), max_size=16 << 20):
            base = group[0].offset
            blob = bb.src.read(base, group[-1].offset + group[-1].size - base)
            for it in group:
                raw = blob[it.offset - base : it.offset - base + it.size]
                if bb.header.uncompress_buf_size:
                    raw = zlib.decompress(raw)
                for pos, rec in decode_records(raw, cid, chrom, self.exclusion, self.pos_field, DAPG_FIELDS):
                    try:
                        pip = float(rec["pip"])
                    except (KeyError, ValueError):
                        continue
                    key = ((rec.get("geneId") or "").split(".")[0], rec.get("geneName") or "")
                    g = gene_code.setdefault(key, len(gene_code))
                    t = tissue_code.setdefault(rec.get("tissue") or "", len(tissue_code))
                    rows.append((pos, g, t, pip))
        if not rows:
            return empty
        rows.sort()
        arr = np.array(rows, dtype=np.float64)
        return Fine(
            pos=arr[:, 0].astype(np.int64),
            gene=arr[:, 1].astype(np.int32),
            tissue=arr[:, 2].astype(np.int16),
            pip=arr[:, 3].astype(np.float32),
            genes=list(gene_code),
            tissues=list(tissue_code),
        )


class CommonVariants:
    """dbSNP 155 Common item positions per chromosome (no other field), cached as arrays.

    A common variant is not an eQTL value; it is the proxy M5 used for a variant GTEx could have
    tested (minor allele 1% or more in a major population), counted here to match windows."""

    def __init__(self, path: Path = COMMON_FILE, cache: Path | None = COMMON_CACHE) -> None:
        self.path = path
        self.cache = cache
        self._bb: BigBed | None = None
        self._pos: dict[str, np.ndarray] = {}

    def positions(self, chrom: str) -> np.ndarray:
        if chrom in self._pos:
            return self._pos[chrom]
        cached = self.cache / f"{chrom}.npy" if self.cache else None
        if cached is not None and cached.exists():
            arr = np.load(cached)
        else:
            arr = self._decode(chrom)
            if cached is not None:
                cached.parent.mkdir(parents=True, exist_ok=True)
                np.save(cached, arr)
        self._pos[chrom] = arr
        return arr

    def _decode(self, chrom: str) -> np.ndarray:
        if self._bb is None:
            self._bb = BigBed(self.path)
        bb = self._bb
        if chrom not in bb.chroms:
            return np.zeros(0, np.int64)
        cid, size = bb.chroms[chrom]
        out: list[int] = []
        for group in coalesce(bb.leaf_items(chrom, [(0, size)]), max_size=64 << 20):
            base = group[0].offset
            blob = bb.src.read(base, group[-1].offset + group[-1].size - base)
            for it in group:
                raw = blob[it.offset - base : it.offset - base + it.size]
                if bb.header.uncompress_buf_size:
                    raw = zlib.decompress(raw)
                out.extend(record_starts(raw, cid))
        arr = np.array(out, dtype=np.int64)
        arr.sort()
        return arr

    def count(self, chrom: str, start: int, end: int) -> int:
        pos = self.positions(chrom)
        return int(np.searchsorted(pos, end, "left") - np.searchsorted(pos, start, "left"))


def record_starts(raw: bytes, chrom_id: int) -> list[int]:
    """The start of every record of one uncompressed bigBed data block on `chrom_id`; nothing else."""
    out = []
    p = 0
    n = len(raw)
    while p + 12 <= n:
        c, s, _e = struct.unpack_from("<III", raw, p)
        p = raw.index(b"\0", p + 12) + 1
        if c == chrom_id:
            out.append(s)
    return out


# ---- the calibration units -----------------------------------------------------------------------


@dataclass
class Unit:
    """One measured element-gene pair at one scale, with its matched windows once drawn."""

    key: str
    chrom: str
    start: int
    end: int
    gene: str
    tss: int
    regulated: bool
    cells: tuple[str, ...]
    element: tuple[int, int]
    common: int = 0
    windows: list[tuple[int, int]] = field(default_factory=list)
    flanks: list[tuple[int, int]] = field(default_factory=list)

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def mid(self) -> int:
        return (self.start + self.end) // 2

    @property
    def distance(self) -> int:
        return abs(self.mid - self.tss)


def element_pairs(pairs: Iterable[Any]) -> list[dict[str, Any]]:
    """The benchmark's valid pairs, one per (element, measured gene): positive if any screen called it.

    `pairs` are `crispri.Pair`s. A pair tested in several screens is one unit; its cells are kept."""
    by: dict[tuple[str, int, int, str], dict[str, Any]] = {}
    for p in pairs:
        if p.tss is None:
            continue
        k = (p.chrom, p.start, p.end, p.gene)
        rec = by.setdefault(
            k,
            {
                "chrom": p.chrom,
                "start": p.start,
                "end": p.end,
                "gene": p.gene,
                "tss": p.tss,
                "regulated": False,
                "cells": set(),
            },
        )
        rec["regulated"] = rec["regulated"] or bool(p.regulated)
        rec["cells"].add(p.cell)
    return [dict(v, cells=tuple(sorted(v["cells"]))) for _, v in sorted(by.items())]


def scaled(start: int, end: int, length: int | None) -> tuple[int, int]:
    """The element with `MARGIN` on each side (length None), or a window of `length` centred on it."""
    s, e = start - MARGIN, end + MARGIN
    if length is None or length <= e - s:
        return max(0, s), e
    mid = (start + end) // 2
    s = mid - length // 2
    return max(0, s), max(0, s) + length


def holds_tss(start: int, end: int, tss: int) -> bool:
    return start - TSS_GUARD <= tss < end + TSS_GUARD


def overlaps(a0: int, a1: int, b0: int, b1: int) -> bool:
    return a0 < b1 and b0 < a1


def build_units(
    pairs: list[dict[str, Any]], exclusion: Exclusion, length: int | None, scale: str
) -> tuple[list[Unit], dict[str, int]]:
    """The units of one scale, and how many pairs each rule dropped (counted on coordinates only)."""
    positives = positive_spans(pairs)
    dropped: dict[str, int] = defaultdict(int)
    out = []
    for p in pairs:
        s, e = scaled(p["start"], p["end"], length)
        mid = (s + e) // 2
        if abs(mid - p["tss"]) > CIS_WINDOW:
            dropped["beyond_the_cis_window"] += 1
            continue
        if exclusion.hits(p["chrom"], s, e):
            dropped["touches_the_exclusion"] += 1
            continue
        if holds_tss(s, e, p["tss"]):
            dropped["holds_the_tss"] += 1
            continue
        if not p["regulated"] and any(
            overlaps(p["start"] - MARGIN, p["end"] + MARGIN, a, b)
            for a, b in positives[(p["chrom"], p["gene"])]
        ):
            dropped["negative_overlapping_a_positive_of_its_gene"] += 1
            continue
        key = f"{scale}:{p['chrom']}:{p['start']}-{p['end']}:{p['gene']}"
        out.append(
            Unit(
                key, p["chrom"], s, e, p["gene"], p["tss"], p["regulated"], p["cells"], (p["start"], p["end"])
            )
        )
    return out, dict(dropped)


def positive_spans(pairs: list[dict[str, Any]]) -> dict[tuple[str, str], list[tuple[int, int]]]:
    """Every measured positive element of each gene, with the margin: no window may touch one."""
    out: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    for p in pairs:
        if p["regulated"]:
            out[(p["chrom"], p["gene"])].append((p["start"] - MARGIN, p["end"] + MARGIN))
    return out


def count_tolerance(common: int) -> int:
    return max(1, round(COUNT_BAND * common))


def draw_windows(
    unit: Unit,
    n: int,
    count: Any,
    exclusion: Exclusion,
    avoid: Sequence[tuple[int, int]],
    chrom_size: int,
    rng: random.Random,
) -> list[tuple[int, int]] | None:
    """`n` disjoint windows matched to the unit, or None when the draw cannot find them.

    Matched on length (equal), distance from the window's midpoint to the gene's TSS (within
    `band` of the unit's, either side of the TSS), and the count of common variants inside
    (within `count_tolerance`). A window never touches the exclusion, the TSS with its guard, the unit,
    another window, or any measured positive element of the same gene. `count(chrom, s, e)` counts."""
    length, d = unit.length, unit.distance
    band = max(DISTANCE_BAND * d, MIN_BAND_LENGTHS * length)
    tol = count_tolerance(unit.common)
    taken = [(unit.start, unit.end)]
    out: list[tuple[int, int]] = []
    misses = 0
    while len(out) < n:
        if misses >= DRAW_TRIES:
            return None  # one window could not be found in its budget: the unit is unmatched
        misses += 1
        dd = rng.uniform(max(0.0, d - band), d + band)
        mid = int(round(unit.tss + rng.choice((-1, 1)) * dd))
        s = mid - length // 2
        e = s + length
        if s < 0 or e > chrom_size or holds_tss(s, e, unit.tss) or exclusion.hits(unit.chrom, s, e):
            continue
        if any(overlaps(s, e, a, b) for a, b in taken) or any(overlaps(s, e, a, b) for a, b in avoid):
            continue
        if abs(count(unit.chrom, s, e) - unit.common) > tol:
            continue
        out.append((s, e))
        taken.append((s, e))
        misses = 0
    return out


def flanks_of(unit: Unit, exclusion: Exclusion, chrom_size: int) -> list[tuple[int, int]]:
    """The two same-length windows immediately beside the unit, where they are clear of the exclusion
    and the TSS guard."""
    out = []
    for s, e in ((unit.start - unit.length, unit.start), (unit.end, unit.end + unit.length)):
        if (
            s >= 0
            and e <= chrom_size
            and not exclusion.hits(unit.chrom, s, e)
            and not holds_tss(s, e, unit.tss)
        ):
            out.append((s, e))
    return out


# ---- scoring -------------------------------------------------------------------------------------


def gene_codes(fine: Fine, symbol_of_id: Mapping[str, str]) -> dict[str, set[int]]:
    """Each name a gene goes by -> its codes in the track: the track's own name (GENCODE v26) and the
    symbol the project's GENCODE gives its Ensembl id."""
    out: dict[str, set[int]] = defaultdict(set)
    for code, (gid, name) in enumerate(fine.genes):
        if name:
            out[name].add(code)
        sym = symbol_of_id.get(gid)
        if sym:
            out[sym].add(code)
    return out


def score(fine: Fine, codes: set[int], start: int, end: int) -> dict[str, Any]:
    """What the fine-mapping says inside [start, end) about one gene (given by its codes)."""
    sl = fine.slice(start, end)
    pip = fine.pip[sl]
    is_fine = pip >= PIP
    if codes:
        is_gene = np.isin(fine.gene[sl], np.fromiter(codes, dtype=np.int32))
    else:
        is_gene = np.zeros(len(pip), dtype=bool)
    both = is_fine & is_gene
    return {
        "gene": bool(both.any()),
        "any": bool(is_fine.any()),
        "gene_any_pip": bool(is_gene.any()),
        "tissues": sorted({fine.tissues[t] for t in fine.tissue[sl][both]}),
    }


def egene(fine: Fine, codes: set[int]) -> bool:
    """Whether the gene has a fine-mapped variant anywhere on the chromosome outside the exclusion."""
    if not codes:
        return False
    mask = np.isin(fine.gene, np.fromiter(codes, dtype=np.int32))
    return bool((fine.pip[mask] >= PIP).any())


def scored_unit(unit: Unit, fine: Fine, codes: set[int]) -> dict[str, Any]:
    """A unit, its windows and its flanks scored against its own gene; coordinates kept, rows not."""
    own = score(fine, codes, unit.start, unit.end)
    wins = [score(fine, codes, s, e) for s, e in unit.windows]
    fl = [score(fine, codes, s, e) for s, e in unit.flanks]
    mapped = sorted({t for c in unit.cells if (t := TISSUE_OF_CELL.get(c))})
    return {
        "key": unit.key,
        "chrom": unit.chrom,
        "gene": unit.gene,
        "regulated": unit.regulated,
        "cells": list(unit.cells),
        "length": unit.length,
        "distance": unit.distance,
        "common": unit.common,
        "y": own["gene"],
        "y_any": own["any"],
        "y_gene_any_pip": own["gene_any_pip"],
        "tissues": own["tissues"],
        "matched_tissues": mapped,
        "y_matched_tissue": bool(set(own["tissues"]) & set(mapped)) if mapped else None,
        "w": [x["gene"] for x in wins],
        "w_any": [x["any"] for x in wins],
        "w_gene_any_pip": [x["gene_any_pip"] for x in wins],
        "w_matched_tissue": [bool(set(x["tissues"]) & set(mapped)) for x in wins] if mapped else None,
        "flanks": [x["gene"] for x in fl],
        "egene": egene(fine, codes),
        "gene_in_track": bool(codes),
    }


# ---- summaries -----------------------------------------------------------------------------------


def _mean(xs: Sequence[float]) -> float | None:
    return round(float(np.mean(xs)), 6) if len(xs) else None


def window_rate(rows: Sequence[Mapping[str, Any]], field_: str = "w") -> float | None:
    """The mean, over units, of the share of a unit's matched windows that hit (all WINDOWS + 1 of them:
    every one was drawn under the same rules, and which of them a test would hold out is chosen at
    random in the sizing)."""
    return _mean([float(np.mean(r[field_])) for r in rows if r[field_]])


def icc(matrix: np.ndarray) -> float | None:
    """One-way ANOVA intraclass correlation of a units x windows 0/1 matrix (equal group sizes)."""
    if matrix.ndim != 2 or matrix.shape[0] < 2 or matrix.shape[1] < 2:
        return None
    g, k = matrix.shape
    means = matrix.mean(1)
    grand = matrix.mean()
    msb = k * ((means - grand) ** 2).sum() / (g - 1)
    msw = ((matrix - means[:, None]) ** 2).sum() / (g * (k - 1))
    den = msb + (k - 1) * msw
    return round(float((msb - msw) / den), 6) if den > 0 else None


def cluster_interval(
    groups: Sequence[str], values: Sequence[float], rng: np.random.Generator, n: int = BOOTSTRAPS
) -> list[float] | None:
    """A 95% percentile interval for the mean of `values`, resampling their GROUPS (the genes): units
    sharing a gene move together, because they share its eQTLs."""
    if not len(values):
        return None
    ug, inv = np.unique(np.asarray(groups), return_inverse=True)
    g = len(ug)
    if g < 10:
        return None
    vals = np.asarray(values, dtype=float)
    sums = np.bincount(inv, weights=vals, minlength=g)
    counts = np.bincount(inv, minlength=g).astype(float)
    est = np.empty(n)
    for b in range(n):
        w = np.bincount(rng.integers(0, g, g), minlength=g)
        est[b] = (w @ sums) / (w @ counts)
    return [round(float(np.quantile(est, 0.025)), 6), round(float(np.quantile(est, 0.975)), 6)]


DISTANCE_BANDS = ((0, 10_000), (10_000, 50_000), (50_000, 200_000), (200_000, CIS_WINDOW + 1))


def summarise(
    rows: Sequence[Mapping[str, Any]], rng: np.random.Generator, bootstraps: int = BOOTSTRAPS
) -> dict[str, Any]:
    """Every calibrated quantity of one group of units (positives or negatives) at one scale."""
    n = len(rows)
    if not n:
        return {"units": 0}

    def sens(rs: Sequence[Mapping[str, Any]]) -> float | None:
        return _mean([float(r["y"]) for r in rs])

    def diff(rs: Sequence[Mapping[str, Any]]) -> float | None:
        a, b = sens(rs), window_rate(rs)
        return None if a is None or b is None else a - b

    w = np.array([r["w"] for r in rows], dtype=float)
    y = np.array([r["y"] for r in rows], dtype=float)
    groups = [f"{r['chrom']}:{r['gene']}" for r in rows]
    holding_any = [r for r in rows if r["y_any"]]
    win_any = sum(sum(r["w_any"]) for r in rows)
    win_named = sum(sum(a and b for a, b in zip(r["w_any"], r["w"], strict=True)) for r in rows)
    mapped = [r for r in rows if r["y_matched_tissue"] is not None]
    tissues: dict[str, int] = defaultdict(int)
    for r in rows:
        for t in r["tissues"]:
            tissues[t] += 1
    flank_hit = [r for r in rows if r["y"] and r["flanks"]]
    flank_miss = [r for r in rows if not r["y"] and r["flanks"]]
    out: dict[str, Any] = {
        "units": n,
        "genes": len({(r["chrom"], r["gene"]) for r in rows}),
        "units_with_a_hit": int(y.sum()),
        "rate_at_the_units": sens(rows),
        "rate_at_the_units_interval_over_genes": cluster_interval(groups, y, rng, bootstraps),
        "rate_at_the_windows": window_rate(rows),
        "rate_at_the_windows_interval_over_genes": cluster_interval(groups, w.mean(1), rng, bootstraps),
        "rate_at_the_first_drawn_window": _mean([float(r["w"][0]) for r in rows]),
        "difference_units_minus_windows": None if (d := diff(rows)) is None else round(d, 6),
        "difference_interval_over_genes": cluster_interval(groups, y - w.mean(1), rng, bootstraps),
        "naming": {
            "units_holding_a_fine_mapped_variant_for_any_gene": len(holding_any),
            "of_them_naming_the_units_own_gene": sum(1 for r in holding_any if r["y"]),
            "windows_holding_one_for_any_gene": int(win_any),
            "of_them_naming_the_units_gene": int(win_named),
        },
        "matched_tissue": {
            "units_whose_cell_maps_to_a_tissue": len(mapped),
            "with_a_hit_in_that_tissue": sum(1 for r in mapped if r["y_matched_tissue"]),
            "rate_at_the_units": _mean([float(r["y_matched_tissue"]) for r in mapped]),
            "rate_at_the_windows": window_rate(mapped, "w_matched_tissue"),
        },
        "units_with_a_hit_by_tissue": dict(sorted(tissues.items(), key=lambda kv: (-kv[1], kv[0]))),
        "any_pip_descriptive": {
            "meaning": "a variant of the track for the gene at any PIP (a member of a DAP-G signal "
            "cluster), which is LD-tagged association, not localised",
            "rate_at_the_units": _mean([float(r["y_gene_any_pip"]) for r in rows]),
            "rate_at_the_windows": window_rate(rows, "w_gene_any_pip"),
        },
        "fine_mapped_egene": {
            "units_whose_gene_has_a_fine_mapped_variant_outside_the_exclusion": sum(
                1 for r in rows if r["egene"]
            ),
            "units_whose_gene_is_in_the_track_at_all": sum(1 for r in rows if r["gene_in_track"]),
            "rate_at_the_units_among_egene_units": sens([r for r in rows if r["egene"]]),
            "rate_at_the_windows_among_egene_units": window_rate([r for r in rows if r["egene"]]),
        },
        "by_distance": {
            f"{lo}-{hi}": {
                "units": len(b),
                "rate_at_the_units": sens(b),
                "rate_at_the_windows": window_rate(b),
            }
            for lo, hi in DISTANCE_BANDS
            for b in [[r for r in rows if lo <= r["distance"] < hi]]
        },
        "clustering": {
            "icc_of_the_windows_within_a_unit": icc(w),
            "design_effect_of_a_units_window_mean": None
            if (rho := icc(w)) is None
            else round(1 + (w.shape[1] - 1) * rho, 4),
            "correlation_of_the_units_hit_with_its_windows_mean": None
            if y.std() == 0 or w.mean(1).std() == 0
            else round(float(np.corrcoef(y, w.mean(1))[0, 1]), 6),
            "flank_rate_when_the_unit_hits": _mean([float(np.mean(r["flanks"])) for r in flank_hit]),
            "flank_rate_when_it_does_not": _mean([float(np.mean(r["flanks"])) for r in flank_miss]),
            "units_with_flanks": {"hit": len(flank_hit), "no_hit": len(flank_miss)},
        },
    }
    terc = np.quantile([r["common"] for r in rows], [1 / 3, 2 / 3]) if n >= 3 else None
    if terc is not None:
        cuts = [(-1, terc[0]), (terc[0], terc[1]), (terc[1], 1e12)]
        out["by_common_variant_tertile"] = {
            f"{lo:g}<c<={hi:g}": {
                "units": len(b),
                "rate_at_the_units": sens(b),
                "rate_at_the_windows": window_rate(b),
            }
            for lo, hi in cuts
            for b in [[r for r in rows if lo < r["common"] <= hi]]
        }
    return out


# ---- the test the sizing is for, and the sizing --------------------------------------------------


def upper_tail(probs: np.ndarray, observed: int) -> float:
    """P(S >= observed) for S a sum of independent Bernoulli(probs): exact below 400 informative
    terms, normal with a continuity correction above."""
    ones = int((probs >= 1).sum())
    q = probs[(probs > 0) & (probs < 1)]
    t = observed - ones
    if t <= 0:
        return 1.0
    if t > len(q):
        return 0.0
    if len(q) <= 400:
        dist = np.zeros(len(q) + 1)
        dist[0] = 1.0
        for i, p in enumerate(q):
            dist[1 : i + 2] = dist[1 : i + 2] * (1 - p) + dist[0 : i + 1] * p
            dist[0] *= 1 - p
        return float(min(1.0, dist[t:].sum()))
    mu = q.sum()
    sd = math.sqrt(float((q * (1 - q)).sum()))
    z = (t - 0.5 - mu) / sd
    return 0.5 * math.erfc(z / math.sqrt(2))


def conditional_p(y: np.ndarray, wsum: np.ndarray, k: int = WINDOWS) -> float:
    """The registered statistic's one-sided p-value: under the null a block is exchangeable with its
    `k` windows, so given h hits among its k + 1 intervals the block is the hit with probability
    h / (k + 1), independently across blocks; S is the number of blocks that hit."""
    h = y + wsum
    return upper_tail(h / (k + 1), int(y.sum()))


Source = tuple[np.ndarray | None, np.ndarray]


def draw(source: Source, idx: np.ndarray, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    """Blocks drawn from a calibrated source: (hit at the block, hits among its `WINDOWS` controls).

    A source is (hits at the units, or None, and a units x (WINDOWS + 1) matrix of window hits). With
    units, the block is the unit and one of its windows, chosen at random, is left out, so a block has
    `WINDOWS` controls. Without (the window null), one window chosen at random is the block and the
    other `WINDOWS` its controls, which makes that null exchangeable by construction: which window
    was drawn first carries no weight."""
    y, w = source
    w = w[idx]
    j = rng.integers(0, w.shape[1], len(idx))
    left = w[np.arange(len(idx)), j]
    block = left if y is None else y[idx]
    return block, w.sum(1) - left


def power(
    right: Source,
    wrong: Source,
    n: int,
    share: float,
    rng: np.random.Generator,
    experiments: int = EXPERIMENTS,
    weights: np.ndarray | None = None,
) -> float:
    """Probability that the conditional test rejects at `ALPHA` with `n` blocks of which a share
    `share` behave as a calibrated positive (`right`, resampled with `weights` if given) and the rest
    as `wrong`; each block brings its own `WINDOWS` windows."""
    n_right, n_wrong = len(right[1]), len(wrong[1])
    rejected = 0
    for _ in range(experiments):
        r = int(rng.binomial(n, share))
        ir = rng.choice(n_right, r, p=weights) if weights is not None else rng.integers(0, n_right, r)
        iw = rng.integers(0, n_wrong, n - r)
        yr, sr = draw(right, ir, rng)
        yw, sw = draw(wrong, iw, rng)
        if conditional_p(np.concatenate([yr, yw]), np.concatenate([sr, sw])) <= ALPHA:
            rejected += 1
    return rejected / experiments


def source(rows: Sequence[Mapping[str, Any]], units: bool = True) -> Source:
    """A calibrated group as a sizing source: its units' hits (or None for the window null) and its
    window hits."""
    w = np.array([r["w"] for r in rows], dtype=float).reshape(len(rows), WINDOWS + 1)
    return (np.array([r["y"] for r in rows], dtype=float) if units else None), w


def band_of(distance: int) -> int:
    return next(i for i, (lo, hi) in enumerate(DISTANCE_BANDS) if lo <= distance < hi)


def distance_weights(rows: Sequence[Mapping[str, Any]], target: Mapping[int, float]) -> dict[str, Any]:
    """Resampling weights that give the positives the blocks' distribution over `DISTANCE_BANDS`
    (`target`: band index -> share of the blocks). A band the blocks occupy and the positives do not
    cannot be represented; its share is reported and the rest renormalised."""
    bands = np.array([band_of(int(r["distance"])) for r in rows])
    have = {b: int((bands == b).sum()) for b in range(len(DISTANCE_BANDS))}
    missing = sum(s for b, s in target.items() if s > 0 and have.get(b, 0) == 0)
    w = np.array([target.get(int(b), 0.0) / have[int(b)] for b in bands], dtype=float)
    total = w.sum()
    return {
        "weights": (w / total) if total > 0 else None,
        "target_share_not_represented": round(missing, 4),
        "positives_by_band": {f"{lo}-{hi}": have[i] for i, (lo, hi) in enumerate(DISTANCE_BANDS)},
    }


def sizing(
    positives: Sequence[Mapping[str, Any]],
    negatives: Sequence[Mapping[str, Any]],
    available: int,
    rng: np.random.Generator,
    grid: Sequence[int] = N_GRID,
    shares: Sequence[float] = SHARES,
    experiments: int = EXPERIMENTS,
    target_bands: Mapping[int, float] | None = None,
) -> dict[str, Any]:
    """Power over the grid and at the available count, per share right and per null model; and, under
    the window null, with the positives reweighted to the blocks' distances to their genes."""
    right = source(positives)
    nulls: dict[str, Source | None] = {
        "window": source(positives, units=False),
        "tested_negative": source(negatives) if negatives else None,
    }
    ns = sorted(set(grid) | ({available} if available > 0 else set()))
    out: dict[str, Any] = {}
    for name, wrong in nulls.items():
        if wrong is None or not len(wrong[1]) or not len(right[1]):
            out[name] = None
            continue
        p0 = float(np.mean(wrong[1])) if wrong[0] is None else float(np.mean(wrong[0]))
        sens = float(np.mean(right[0]))
        win = float(np.mean(np.concatenate([right[1], wrong[1]])))
        table = {}
        for f in shares:
            pw = {n: round(power(right, wrong, n, f, rng, experiments), 4) for n in ns}
            reach = [n for n in sorted(pw) if pw[n] >= POWER and n in grid]
            block_rate = f * sens + (1 - f) * p0
            table[str(f)] = {
                "block_rate_expected": round(block_rate, 6),
                "block_to_window_ratio_expected": round(block_rate / win, 4) if win > 0 else None,
                "power_by_blocks": {str(n): pw[n] for n in ns},
                "smallest_grid_blocks_reaching_80": reach[0] if reach else None,
                "power_at_the_available_blocks": pw.get(available),
            }
        fpr = {n: round(power(right, wrong, n, 0.0, rng, experiments), 4) for n in ns}
        out[name] = {
            "wrong_block_rate": round(p0, 6),
            "right_block_rate": round(sens, 6),
            "window_rate": round(win, 6),
            "by_share_right": table,
            "rejection_rate_when_no_attribution_is_right": {str(k): v for k, v in fpr.items()},
            "rejection_rate_at_the_available_blocks_when_none_is_right": fpr.get(available),
        }
    if target_bands is not None and out.get("window") is not None:
        dw = distance_weights(positives, target_bands)
        wts = dw.pop("weights")
        rows_ = {}
        if wts is not None:
            sens_w = float(wts @ right[0])
            for f in (1.0, 0.5):
                pw = {n: round(power(right, nulls["window"], n, f, rng, experiments, wts), 4) for n in ns}
                reach = [n for n in sorted(pw) if pw[n] >= POWER and n in grid]
                rows_[str(f)] = {
                    "power_by_blocks": {str(n): pw[n] for n in ns},
                    "smallest_grid_blocks_reaching_80": reach[0] if reach else None,
                    "power_at_the_available_blocks": pw.get(available),
                }
            dw["right_block_rate_reweighted"] = round(sens_w, 6)
        out["window_distance_reweighted"] = {
            "meaning": "descriptive, not part of the rule: the window null with the positives resampled "
            "to the blocks' distribution of distance to the named gene's TSS",
            "target_share_by_band": {
                f"{lo}-{hi}": round(target_bands.get(i, 0.0), 4) for i, (lo, hi) in enumerate(DISTANCE_BANDS)
            },
            **dw,
            "by_share_right": rows_,
        }
    return out


def fair_test(sizes: Mapping[str, Any], available: int) -> dict[str, Any]:
    """FAIR_TEST_RULE applied to one scale's sizing and its count of available blocks."""
    why = []
    win, neg = sizes.get("window"), sizes.get("tested_negative")
    if available <= 0:
        why.append("no block is available to the test at this scale")
    if win is None:
        why.append("the window null could not be calibrated")
    else:
        half = win["by_share_right"]["0.5"]
        got, need = half["power_at_the_available_blocks"], half["smallest_grid_blocks_reaching_80"]
        if got is None or got < POWER:
            why.append(
                f"(1) {available} available blocks give {got} probability of detection at f = 1/2; 80% "
                f"needs {need if need is not None else f'more than {max(N_GRID)}'} blocks"
            )
    fooled = None
    if neg is None:
        why.append("the tested-negative null could not be calibrated")
    else:
        fooled = neg["rejection_rate_at_the_available_blocks_when_none_is_right"]
        if fooled is None or fooled > FOOLED_BOUND:
            why.append(
                f"(2) with none right and the blocks behaving as tested non-regulating elements, the "
                f"test rejects {fooled} of the time at {available} blocks, above {FOOLED_BOUND}"
            )
    return {"fair": not why, "available_blocks": available, "fooled_rate": fooled, "why_not": why}


# ---- the blocks' attribution: what the test would use, read without any eQTL value ---------------


def block_attribution(
    chrom: str, blocks: list[dict[str, Any]], elements: list[dict[str, Any]], min_log2: float
) -> list[dict[str, Any]]:
    """Each real-unknown block of one chromosome with the gene its attribution names most strongly.

    Membership is `measured_unknown.population_of`'s (an element belongs to the block holding its
    midpoint); the gene is the one whose predicted deletion change is largest among the block's
    elements at |log2| >= `min_log2`, with that element's tissue and coordinates."""
    els = sorted(elements, key=lambda e: (e["start"] + e["end"]) // 2)
    mids = [(e["start"] + e["end"]) // 2 for e in els]
    out = []
    for b in blocks:
        if b.get("tier") != "constrained_unknown" or b.get("copy"):
            continue
        lo, hi = bisect.bisect_left(mids, b["start"]), bisect.bisect_left(mids, b["end"])
        best = None
        for e in els[lo:hi]:
            pred = e.get("predicted") or {}
            lf = pred.get("log2_fold_change")
            if not pred.get("gene") or lf is None or abs(lf) < min_log2:
                continue
            if best is None or abs(lf) > abs(best[0]):
                best = (lf, pred["gene"], pred.get("tissue"), (e["start"], e["end"]))
        out.append(
            {
                "chrom": chrom,
                "start": b["start"],
                "end": b["end"],
                "gene": best[1] if best else None,
                "tissue": best[2] if best else None,
                "log2": best[0] if best else None,
                "element": list(best[3]) if best else None,
            }
        )
    return out


def available(
    rows: Sequence[Mapping[str, Any]], tss_of: Mapping[str, Sequence[tuple[str, int]]]
) -> dict[str, Any]:
    """How many blocks a test could use, per scale, and why the others cannot (coordinates only)."""
    reasons: dict[str, int] = defaultdict(int)
    usable: dict[str, list[int]] = {"block": [], "element": []}
    for r in rows:
        if not r["gene"]:
            reasons["names_no_gene"] += 1
            continue
        mid = (r["start"] + r["end"]) // 2
        cands = [t for c, t in tss_of.get(r["gene"], ()) if c == r["chrom"]]
        if not cands:
            reasons["gene_without_a_tss_on_the_chromosome"] += 1
            continue
        tss = min(cands, key=lambda t: abs(t - mid))
        if abs(tss - mid) > CIS_WINDOW:
            reasons["beyond_the_cis_window"] += 1
            continue
        if holds_tss(r["start"], r["end"], tss):
            reasons["block_holds_the_tss"] += 1
        else:
            usable["block"].append(abs(tss - mid))
        es, ee = scaled(r["element"][0], r["element"][1], None)
        if not holds_tss(es, ee, tss):
            usable["element"].append(abs(tss - (es + ee) // 2))
    return {
        "blocks": len(rows),
        "naming_a_gene": sum(1 for r in rows if r["gene"]),
        "usable_at_the_block_scale": len(usable["block"]),
        "usable_at_the_element_scale": len(usable["element"]),
        "not_usable": dict(reasons),
        "distance_quartiles": {
            k: [int(x) for x in np.quantile(v, [0.25, 0.5, 0.75])] if v else None for k, v in usable.items()
        },
        "distance_band_shares": {
            k: {i: sum(1 for d in v if band_of(int(d)) == i) / len(v) for i in range(len(DISTANCE_BANDS))}
            if v
            else None
            for k, v in usable.items()
        },
    }
