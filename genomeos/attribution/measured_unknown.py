# SPDX-License-Identifier: AGPL-3.0-or-later
"""Milestone 1.3, clause 2, from public measurement: which public assays reach the real unknown.

The clause asks for "the constrained-unknown blocks attributed to a gene and a tissue". On 2026-09-27
the blocks named a coding gene 16.2 points below length-matched random windows, and every control since
has left the model arm below its windows; the measured arm the project could run (`c17eedc`/`a7f207f`)
found 3 blocks with an element CRISPRi ever tested against a coding gene. Permanent non-goal D8 rules out
synthesis, so the clause can only move on measurement somebody else has already made and published.

This module is the census of that measurement, and nothing more. It holds:

- `population`, the 882 real-unknown blocks (organiser tier `constrained_unknown`, copies out) with
  their scored elements and the gene and tissue the attribution names, checked against the counts every
  clause-2 lane has used (882 blocks, 531 carrying a scored element, 3,280 elements);
- `Chain`, a minimal UCSC chain lift (hg19 to GRCh38) at liftOver's default of 95% of bases mapped,
  for the public tables that were never re-released on GRCh38;
- `census`, which asks of one source's tested intervals, by coordinates and tested gene only, how many
  blocks they reach and how.

**No readout value enters any count.** A source is reduced to `Tested` rows (chromosome, start, end,
and, where the assay names one, the gene it was measured against and the cell it was measured in) by
parsers that select those columns by name and never touch an effect size, a p-value, an activity or a
call. Whether an element *did* anything is not asked here; only whether it was *asked*.

**What "enough" means.** The calibrated power simulation (lane-s2, `019a9ab` registered, `cab70d9`
result) gives, for 80% probability of detection at equal cost and with each block compared against its
own windows, these designs in compared blocks and tested elements: a block-to-window ratio of a half
needs all 531 blocks and 1,062 elements; a quarter 200 and 400; a tenth 75 and 300; none 100 and 200; a
three-quarters ratio is infeasible. Those figures were calibrated on CRISPRi element-gene readouts (the
benchmark's own sensitivity, intraclass correlation and window anchor), so they apply to a gene-level
perturbation endpoint and to nothing else. `FLOOR_BLOCKS` is the smallest of them: below 75 compared
blocks no calibrated design reaches 80% at any ratio.
"""

from __future__ import annotations

import bisect
import gzip
import json
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from genomeos.attribution import organise
from genomeos.attribution.measured import RECIPROCAL_OVERLAP

CHROMS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
ELEMENTS = Path("data/knowledge/alphagenome/all_elements")
STORE = Path("data/knowledge/measured_unknown")
#: the sweep's own threshold for "this deletion moved a gene" (constrained_unknown_targets.MIN_LOG2)
MIN_LOG2 = 0.1
#: the counts every clause-2 lane has used; a population that differs is a different population
GATE = {"blocks": 882, "carrying": 531, "elements": 3280}
#: lane-s2's calibrated designs (`cab70d9`): true block-to-window ratio -> (compared blocks, tested
#: elements) for 80% probability of detection at equal cost; None where no design reaches it
SIZING: dict[str, tuple[int, int] | None] = {
    "three_quarters": None,
    "half": (531, 1062),
    "quarter": (200, 400),
    "tenth": (75, 300),
    "none": (100, 200),
}
FLOOR_BLOCKS = 75
FLOOR_ELEMENTS = 200
#: liftOver's default minMatch: an interval lifts only when this share of its bases map on one chain
MIN_MATCH = 0.95

#: what a source's readout is about, which decides whether it can speak to clause 2 at all
ENDPOINTS = {
    "gene": "an element perturbed in its own chromosome and a named gene's expression measured",
    "cell": "an element perturbed in its own chromosome and the cell's growth or survival measured",
    "activity": "a sequence's reporter activity measured in a cell, outside or inside its locus, "
    "with no gene named",
    "variant_gene": "natural variants associated with a named gene's expression in a named tissue",
    "trait": "natural variants associated with an organism-level trait, with no gene or tissue named",
}
#: the endpoints the clause's "attributed to a gene" can be scored against
GENE_ENDPOINTS = ("gene", "variant_gene")


@dataclass(frozen=True, slots=True)
class Tested:
    """One interval a public assay tested, reduced to coordinates and the gene and cell it names."""

    chrom: str
    start: int  # 0-based, half-open
    end: int
    gene: str | None = None
    cell: str | None = None


# ---- the population ---------------------------------------------------------------------------------


def _elements(chrom: str, root: Path = ELEMENTS) -> list[dict[str, Any]]:
    p = root / f"{chrom}.json"
    return json.loads(p.read_text()) if p.exists() else []


def _moves(e: dict[str, Any]) -> bool:
    pred = e.get("predicted") or {}
    log2 = pred.get("log2_fold_change")
    return bool(pred.get("gene")) and log2 is not None and abs(log2) >= MIN_LOG2


def population_of(
    chrom: str, blocks: list[dict[str, Any]], elements: list[dict[str, Any]]
) -> list[dict[str, Any]]:
    """The real-unknown blocks of one chromosome, each with its scored elements and what the attribution
    names. An element belongs to the block that holds its midpoint, as `constrained_unknown_targets`."""
    els = sorted(elements, key=lambda e: (e["start"] + e["end"]) // 2)
    mids = [(e["start"] + e["end"]) // 2 for e in els]
    out = []
    for b in blocks:
        if b.get("tier") != "constrained_unknown" or b.get("copy"):
            continue
        lo, hi = bisect.bisect_left(mids, b["start"]), bisect.bisect_left(mids, b["end"])
        got = els[lo:hi]
        coding = [e["predicted_coding"] for e in got if e.get("predicted_coding")]
        out.append(
            {
                "chrom": chrom,
                "start": b["start"],
                "end": b["end"],
                "length": b["end"] - b["start"],
                "case": b.get("case"),
                "elements": [(e["start"], e["end"]) for e in got],
                "named_genes": sorted(
                    {(e["predicted"] or {}).get("gene") for e in got if _moves(e)} - {None}
                ),
                "named_cells": sorted(
                    {(e["predicted"] or {}).get("tissue") for e in got if _moves(e)} - {None}
                ),
                "named_coding_genes": sorted({c.get("gene") for c in coding} - {None}),
                "named_coding_cells": sorted({c.get("tissue") for c in coding} - {None}),
            }
        )
    return out


def population(chroms: Iterable[str] = CHROMS) -> list[dict[str, Any]]:
    """The 882 real-unknown blocks genome-wide, read through the organiser and the element archive."""
    rows: list[dict[str, Any]] = []
    for c in chroms:
        rows += population_of(c, organise.blocks(c), _elements(c))
    return rows


def gate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """The population against the counts every clause-2 lane used; `ok` False means a different one."""
    got = {
        "blocks": len(rows),
        "carrying": sum(1 for r in rows if r["elements"]),
        "elements": sum(len(r["elements"]) for r in rows),
    }
    return {"expected": dict(GATE), "got": got, "ok": got == GATE}


# ---- the lift -----------------------------------------------------------------------------------------


class Chain:
    """A UCSC chain file held as aligned blocks, enough to lift an interval the way liftOver does at its
    default minMatch: the interval lifts onto the one chain that maps the most of its bases, only if that
    chain maps at least `MIN_MATCH` of them, and becomes the span of the bases it maps."""

    def __init__(self) -> None:
        # target chrom -> sorted list of (t_start, t_end, chain_id, q_chrom, q_start, q_strand, q_size)
        self.blocks: dict[str, list[tuple[int, int, int, str, int, str, int]]] = defaultdict(list)
        self._starts: dict[str, list[int]] = {}
        self._longest: dict[str, int] = {}

    @classmethod
    def from_lines(cls, lines: Iterable[str]) -> Chain:
        self = cls()
        head: list[str] | None = None
        t_pos = q_pos = 0
        for raw in lines:
            line = raw.strip()
            if not line:
                continue
            f = line.split()
            if f[0] == "chain":
                head = f
                t_pos, q_pos = int(f[5]), int(f[10])
                continue
            if head is None:
                continue
            size = int(f[0])
            t_chrom, q_chrom, q_strand, q_size, cid = head[2], head[7], head[9], int(head[8]), int(head[12])
            self.blocks[t_chrom].append((t_pos, t_pos + size, cid, q_chrom, q_pos, q_strand, q_size))
            if len(f) == 3:
                t_pos += size + int(f[1])
                q_pos += size + int(f[2])
            else:
                head = None
        for c, bl in self.blocks.items():
            bl.sort()
            self._starts[c] = [b[0] for b in bl]
            self._longest[c] = max(b[1] - b[0] for b in bl)
        return self

    @classmethod
    def from_file(cls, path: Path) -> Chain:
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt") as fh:
            return cls.from_lines(fh)

    def lift(self, chrom: str, start: int, end: int) -> tuple[str, int, int] | None:
        """[start, end) on the target build onto the query build, or None when it does not lift."""
        bl = self.blocks.get(chrom)
        if not bl or end <= start:
            return None
        starts = self._starts[chrom]
        i = max(0, bisect.bisect_left(starts, start - self._longest[chrom]))
        mapped: dict[int, list[tuple[str, int, int]]] = defaultdict(list)
        while i < len(bl) and bl[i][0] < end:
            t0, t1, cid, qc, q0, qs, qsize = bl[i]
            s, e = max(t0, start), min(t1, end)
            if e > s:
                a, b = q0 + (s - t0), q0 + (e - t0)
                if qs == "-":
                    a, b = qsize - b, qsize - a
                mapped[cid].append((qc, a, b))
            i += 1
        if not mapped:
            return None
        best = max(mapped.values(), key=lambda segs: sum(b - a for _, a, b in segs))
        if sum(b - a for _, a, b in best) < MIN_MATCH * (end - start):
            return None
        return best[0][0], min(a for _, a, _ in best), max(b for _, _, b in best)


def lift_all(chain: Chain, rows: Iterable[Tested]) -> tuple[list[Tested], int]:
    """Every row lifted, and how many did not lift."""
    out, lost = [], 0
    for r in rows:
        got = chain.lift(r.chrom, r.start, r.end)
        if got is None:
            lost += 1
            continue
        out.append(Tested(got[0], got[1], got[2], r.gene, r.cell))
    return out, lost


# ---- parsers: coordinates and the named gene, selected by name, never a value ---------------------------


def region(text: str) -> tuple[str, int, int] | None:
    """'chr1:100-200' (optionally with ':strand') as (chrom, start, end); None for anything else."""
    try:
        chrom, rest = text.split(":", 1)
        lo, hi = rest.split(":", 1)[0].split("-", 1)
        return chrom, int(lo), int(hi)
    except ValueError:
        return None


def bed_rows(lines: Iterable[str], cell: str | None = None) -> Iterator[Tested]:
    """The first three columns of a BED file; headers, track lines and malformed rows are skipped."""
    for line in lines:
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 3:
            continue
        try:
            s, e = int(f[1]), int(f[2])
        except ValueError:
            continue
        if f[0].startswith("chr") and e > s:
            yield Tested(f[0], s, e, None, cell)


def table_rows(
    lines: Iterable[str],
    chrom: str,
    start: str,
    end: str,
    gene: str | None = None,
    cell: str | None = None,
    sep: str = "\t",
) -> Iterator[Tested]:
    """A headed table reduced to the named coordinate columns and, if given, the gene column. Only these
    columns are looked up; every other field of the row is split and dropped unread."""
    it = iter(lines)
    header = next(it, "").rstrip("\n").split(sep)
    col = {name: i for i, name in enumerate(header)}
    missing = [n for n in (chrom, start, end, gene) if n is not None and n not in col]
    if missing:
        raise KeyError(f"columns not in header: {missing}")
    ic, i0, i1 = col[chrom], col[start], col[end]
    ig = col[gene] if gene is not None else None
    for line in it:
        f = line.rstrip("\n").split(sep)
        if len(f) <= max(ic, i0, i1, ig or 0):
            continue
        try:
            s, e = int(float(f[i0])), int(float(f[i1]))
        except ValueError:
            continue
        g = f[ig].strip() if ig is not None else None
        if f[ic].startswith("chr") and e > s:
            yield Tested(f[ic], s, e, g if g and g not in ("NA", ".") else None, cell)


def table_or_bed_rows(
    lines: list[str], chrom: str, start: str, end: str, gene: str | None = None, cell: str | None = None
) -> list[Tested]:
    """`table_rows` on a headed table; the first three columns of a file published without its header
    (ENCODE ships some guide tables that way). Anything else raises, rather than guess."""
    try:
        return list(table_rows(lines, chrom, start, end, gene=gene, cell=cell))
    except KeyError:
        first = (lines[0] if lines else "").split("\t")
        if len(first) >= 3 and first[0].startswith("chr") and first[1].isdigit() and first[2].isdigit():
            return list(bed_rows(lines, cell=cell))
        raise


def region_column_rows(lines: Iterable[str], column: str, sep: str = "\t") -> Iterator[Tested]:
    """A headed table whose one coordinate column holds 'chr:start-end' strings."""
    it = iter(lines)
    header = next(it, "").rstrip("\n").split(sep)
    if column not in header:
        raise KeyError(f"column not in header: {column}")
    i = header.index(column)
    for line in it:
        f = line.rstrip("\n").split(sep)
        r = region(f[i]) if len(f) > i else None
        if r and r[0].startswith("chr") and r[2] > r[1]:
            yield Tested(*r)


def fasta_name_rows(lines: Iterable[str]) -> Iterator[tuple[Tested, str, str]]:
    """Guide FASTA whose names carry 'chrom_start_end_strand_...': (interval, strand, spacer)."""
    name = None
    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:]
            continue
        if not name or not line:
            continue
        f = name.split("_")
        if len(f) >= 4 and f[0].startswith("chr"):
            try:
                s, e = int(f[1]), int(f[2])
            except ValueError:
                name = None
                continue
            yield Tested(f[0], min(s, e), max(s, e) + 1), f[3], line.upper()
        name = None


# ---- the census -----------------------------------------------------------------------------------------


def reciprocal(a0: int, a1: int, b0: int, b1: int, fraction: float = RECIPROCAL_OVERLAP) -> bool:
    """The measured layer's rule: the overlap covers `fraction` of each interval."""
    ov = min(a1, b1) - max(a0, b0)
    return ov > 0 and ov >= fraction * (a1 - a0) and ov >= fraction * (b1 - b0)


def _block_index(rows: list[dict[str, Any]]) -> dict[str, tuple[list[int], list[dict[str, Any]]]]:
    by: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in rows:
        by[r["chrom"]].append(r)
    out = {}
    for c, v in by.items():
        v.sort(key=lambda r: r["start"])
        out[c] = ([r["start"] for r in v], v)
    return out


def census(
    rows: list[dict[str, Any]], tested: Iterable[Tested], coding: set[str] | None = None
) -> dict[str, Any]:
    """What one source's tested intervals reach among the blocks, by coordinates and tested gene only.

    `blocks_touched`: a tested interval shares at least one base with the block (the most generous
    rule, `unknown_coverage`'s). `blocks_element_matched`: a tested interval matches one of the block's
    scored elements at the measured layer's 0.5 reciprocal overlap, the rule the measured arm used.
    `tested_inside`: tested intervals whose midpoint lies in a block, the elements a design could count.
    Gene fields are filled only where the source names the gene each interval was measured against."""
    index = _block_index(rows)
    touched: dict[tuple[str, int], set[str]] = defaultdict(set)
    matched: dict[tuple[str, int], int] = defaultdict(int)
    inside: dict[tuple[str, int], int] = defaultdict(int)
    genes: dict[tuple[str, int], set[str]] = defaultdict(set)
    n = 0
    for t in tested:
        n += 1
        got = index.get(t.chrom)
        if not got:
            continue
        starts, blocks = got
        i = bisect.bisect_right(starts, t.end - 1) - 1
        while i >= 0:
            b = blocks[i]
            if b["end"] <= t.start:
                break  # blocks do not overlap one another, so nothing earlier can reach t
            key = (b["chrom"], b["start"])
            touched[key].add(f"{t.start}-{t.end}")
            if t.gene:
                genes[key].add(t.gene)
            if b["start"] <= (t.start + t.end) // 2 < b["end"]:
                inside[key] += 1
            if any(reciprocal(t.start, t.end, e0, e1) for e0, e1 in b["elements"]):
                matched[key] += 1
            i -= 1
    by_key = {(r["chrom"], r["start"]): r for r in rows}
    carrying = {k for k, r in by_key.items() if r["elements"]}
    named = {k: set(r["named_coding_genes"]) | set(r["named_genes"]) for k, r in by_key.items()}
    coding_blocks = {k for k, g in genes.items() if coding is not None and g & coding}
    attributed = {k for k, g in genes.items() if g & named.get(k, set())}
    out: dict[str, Any] = {
        "tested_intervals": n,
        "blocks_touched": len(touched),
        "blocks_touched_carrying": len(set(touched) & carrying),
        "blocks_element_matched": len(matched),
        "elements_matched": sum(matched.values()),
        "tested_inside": sum(inside.values()),
        "blocks_with_tested_inside": {
            f"at_least_{k}": sum(1 for v in inside.values() if v >= k) for k in (1, 2, 4)
        },
        "blocks_touched_by_case": _by_case(by_key, set(touched)),
    }
    if any(genes.values()):
        out["blocks_tested_against_a_coding_gene"] = len(coding_blocks) if coding is not None else None
        out["blocks_tested_against_an_attributed_gene"] = len(attributed)
    out["touched_blocks"] = sorted(f"{c}:{s}" for c, s in touched)
    out["coding_tested_blocks"] = sorted(f"{c}:{s}" for c, s in coding_blocks)
    return out


def _by_case(by_key: dict[tuple[str, int], dict[str, Any]], keys: set[tuple[str, int]]) -> dict[str, int]:
    out: dict[str, int] = defaultdict(int)
    for k in keys:
        out[str(by_key[k].get("case"))] += 1
    return dict(sorted(out.items()))


def shortfall(blocks: int, endpoint: str) -> dict[str, Any]:
    """Where a source stands against the calibrated designs. Only a gene-level endpoint is sized there;
    for any other the shortfall is the whole requirement, because its readout cannot score a gene."""
    gene_level = endpoint in GENE_ENDPOINTS
    usable = blocks if gene_level else 0
    half = SIZING["half"]
    return {
        "gene_level_endpoint": gene_level,
        "blocks_usable_for_the_clause": usable,
        "short_of_the_floor": max(0, FLOOR_BLOCKS - usable),
        "short_of_a_half": max(0, (half[0] if half else 0) - usable),
        "designs_within_reach": [k for k, v in SIZING.items() if v is not None and usable >= v[0]],
    }
