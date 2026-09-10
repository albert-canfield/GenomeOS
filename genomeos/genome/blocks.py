"""Blocks of a chromosome window for the 2-D block viewer/editor.

A block is any delimited stretch we can name: gene, transcript, exon, CDS,
UTR, CpG island, repeat, telomeric block, assembly gap, and UNKNOWN (the
intergenic space between annotated genes, which is not "nothing": it is
where regulation lives and where our knowledge is thinnest). Blocks nest
through `parent`; every block carries evidence and a confidence so the
viewer can show what is curated, what is inferred, and what is unknown.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

from genomeos.genome.sequence import Sequence

COARSE_WINDOW = 3_000_000  # above this, only genes / unknown / gaps / islands are returned
ELEMENT_WINDOW = 5_000_000  # sequence elements are scanned only for windows up to this size


@dataclass(slots=True)
class Block:
    id: str
    type: str
    start: int
    end: int
    strand: str = "."
    name: str = ""
    parent: str | None = None
    evidence: str = "none"
    confidence: float = 0.0
    attrs: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "id": self.id,
            "type": self.type,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "name": self.name,
            "parent": self.parent,
            "evidence": self.evidence,
            "confidence": self.confidence,
            "attrs": self.attrs,
        }


def blocks_for_window(
    annotation,
    chrom: str,
    start: int,
    end: int,
    seq: Sequence | str | None = None,
    chrom_length: int | None = None,
) -> dict[str, Any]:
    blocks: list[Block] = []
    coarse = (end - start) > COARSE_WINDOW
    genes = [
        g
        for g in annotation.genes.values()
        if g.locus.chrom == chrom and g.locus.end > start and g.locus.start < end
    ]
    genes.sort(key=lambda g: g.locus.start)
    for g in genes:
        gid = g.id
        blocks.append(
            Block(
                gid,
                "gene",
                g.locus.start,
                g.locus.end,
                g.locus.strand.value,
                g.symbol,
                None,
                "curated",
                0.95,
                {"gene_type": g.type, "transcripts": len(g.transcripts)},
            )
        )
        if coarse:
            continue
        for t in g.transcripts.values():
            canonical = "Ensembl_canonical" in t.tags or len(g.transcripts) == 1
            blocks.append(
                Block(
                    t.id,
                    "transcript",
                    t.locus.start,
                    t.locus.end,
                    t.locus.strand.value,
                    t.name,
                    gid,
                    "curated",
                    0.9 if canonical else 0.7,
                    {"transcript_type": t.type, "canonical": canonical, "exons": len(t.exons)},
                )
            )
            cds = sorted((c for c, _ in t.cds), key=lambda l: l.start)
            for i, e in enumerate(sorted(t.exons, key=lambda l: l.start)):
                blocks.append(
                    Block(
                        f"{t.id}:exon{i + 1}",
                        "exon",
                        e.start,
                        e.end,
                        e.strand.value,
                        f"exon {i + 1}",
                        t.id,
                        "curated",
                        0.9,
                    )
                )
            for i, c in enumerate(cds):
                blocks.append(
                    Block(
                        f"{t.id}:cds{i + 1}",
                        "cds",
                        c.start,
                        c.end,
                        c.strand.value,
                        f"CDS {i + 1}",
                        t.id,
                        "curated",
                        0.9,
                    )
                )
    # UNKNOWN: the space between genes, evidence none, role unknown. The minimum
    # size scales with the window so small genomes (mtDNA's control region) show theirs.
    min_unknown = min(1000, max(50, (end - start) // 50))
    cursor = start
    for g in genes:
        if g.locus.start - cursor >= min_unknown:
            blocks.append(
                Block(
                    f"unknown:{cursor}-{g.locus.start}",
                    "unknown",
                    cursor,
                    g.locus.start,
                    ".",
                    "UNKNOWN",
                    None,
                    "none",
                    0.0,
                    {"role": "unknown", "note": "intergenic: regulation lives here"},
                )
            )
        cursor = max(cursor, g.locus.end)
    if end - cursor >= 1000:
        blocks.append(
            Block(
                f"unknown:{cursor}-{end}",
                "unknown",
                cursor,
                end,
                ".",
                "UNKNOWN",
                None,
                "none",
                0.0,
                {"role": "unknown"},
            )
        )
    # sequence elements
    if seq is not None and (end - start) <= ELEMENT_WINDOW:
        s = str(Sequence(str(seq)))
        for m in re.finditer(r"N{100,}", s):
            blocks.append(
                Block(
                    f"gap:{start + m.start()}",
                    "gap",
                    start + m.start(),
                    start + m.end(),
                    ".",
                    "gap",
                    None,
                    "curated",
                    1.0,
                )
            )
        window = 200
        isl: list[list[int]] = []
        for i in range(0, len(s) - window + 1, window):
            w = s[i : i + window]
            c, g = w.count("C"), w.count("G")
            if c and g and (c + g) / window > 0.5 and w.count("CG") * window / (c * g) > 0.6:
                if isl and isl[-1][1] == i:
                    isl[-1][1] = i + window
                else:
                    isl.append([i, i + window])
        for a, b in isl:
            if b - a >= 400:
                blocks.append(
                    Block(
                        f"cpg:{start + a}",
                        "cpg_island",
                        start + a,
                        start + b,
                        ".",
                        "CpG island",
                        None,
                        "predicted",
                        0.6,
                        {"note": "GC > 50%, CpG obs/exp > 0.6; promoter marker in mammals"},
                    )
                )
        for m in re.finditer(
            r"A{12,}|C{12,}|G{12,}|T{12,}|(?:AC){10,}|(?:AG){10,}|(?:AT){10,}|(?:CA){10,}|(?:CT){10,}|(?:GA){10,}|(?:GT){10,}|(?:TA){10,}|(?:TC){10,}|(?:TG){10,}",
            s,
        ):
            blocks.append(
                Block(
                    f"rep:{start + m.start()}",
                    "repeat",
                    start + m.start(),
                    start + m.end(),
                    ".",
                    m.group()[:6] + "…",
                    None,
                    "predicted",
                    0.7,
                )
            )
        for m in re.finditer(r"(?:TTAGGG){3,}|(?:CCCTAA){3,}", s):
            blocks.append(
                Block(
                    f"tel:{start + m.start()}",
                    "telomere",
                    start + m.start(),
                    start + m.end(),
                    ".",
                    "telomeric repeat",
                    None,
                    "curated",
                    0.9,
                )
            )
    _apply_unknown_classes(blocks, chrom)
    return {
        "chrom": chrom,
        "start": start,
        "end": end,
        "length": chrom_length,
        "coarse": coarse,
        "blocks": [b.to_dict() for b in blocks],
        "counts": _counts(blocks),
    }


def _apply_unknown_classes(blocks: list[Block], chrom: str) -> None:
    """If `genomeos unknown` has classified this chromosome, label the UNKNOWN blocks with its findings."""
    from genomeos.results import load_result

    r = load_result(f"unknown_{chrom}")
    if not r:
        return
    by_start = {(b["start"], b["end"]): b for b in r["blocks"]}
    for b in blocks:
        if b.type != "unknown":
            continue
        hit = by_start.get((b.start, b.end))
        if not hit:
            # window-clipped unknown spans: find the classified block containing this one
            hit = next((c for c in r["blocks"] if c["start"] <= b.start and b.end <= c["end"]), None)
        if hit and hit["class"] != "unclassified":
            b.name = f"UNKNOWN · {hit['class']}"
            b.evidence = hit["evidence"]
            b.confidence = hit["confidence"]
            b.attrs.update(
                {
                    "class": hit["class"],
                    **{f"f_{k}": v for k, v in hit["features"].items()},
                    "patterns": ", ".join(f"{k}×{v}" for k, v in hit["patterns"].items()) or "",
                }
            )
            if hit.get("similar_to"):
                b.attrs["similar_to"] = "; ".join(hit["similar_to"][:3])


def _counts(blocks: list[Block]) -> dict[str, int]:
    out: dict[str, int] = {}
    for b in blocks:
        out[b.type] = out.get(b.type, 0) + 1
    return out


def check_move(
    blocks: list[dict[str, Any]], block_id: str, new_start: int, chrom_length: int | None
) -> dict[str, Any]:
    """Plausibility of moving a block (with its children) to a new start. Returns verdict + reasons."""
    by_id = {b["id"]: b for b in blocks}
    b = by_id.get(block_id)
    if not b:
        return {"ok": False, "reasons": ["unknown block"]}
    delta = new_start - b["start"]
    reasons: list[str] = []
    length = b["end"] - b["start"]
    if new_start < 0 or (chrom_length and new_start + length > chrom_length):
        reasons.append("outside the chromosome")
    parent = by_id.get(b["parent"]) if b.get("parent") else None
    if parent and (new_start < parent["start"] or new_start + length > parent["end"]):
        reasons.append(f"leaves its parent {parent['type']} {parent['name']}")
    siblings = [
        o
        for o in blocks
        if o.get("parent") == b.get("parent") and o["id"] != block_id and o["type"] == b["type"]
    ]
    for o in siblings:
        if b["type"] in ("gene",) and o["strand"] != b["strand"]:
            continue
        if new_start < o["end"] and new_start + length > o["start"]:
            reasons.append(
                f"overlaps {o['type']} {o['name']}" + (" on the same strand" if b["type"] == "gene" else "")
            )
            break
    if b["type"] == "exon" and delta % 3 != 0:
        reasons.append("shifts the reading frame (move by a multiple of 3 to keep the frame, if coding)")
    return {
        "ok": not reasons,
        "reasons": reasons,
        "delta": delta,
        "note": "a plausible move is a design proposal: evidence predicted, confidence 0.3",
    }
