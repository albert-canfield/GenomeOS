"""Regulatory elements as BioIR entities, with the genes they reach.

ENCODE tells us where the elements are and what kind they are; it does not
say which gene each one regulates. GenomeOS assigns targets with explicit
bases, each with its own evidence:

  promoter-like (PLS)        → the gene whose transcription start lies within
                               PROMOTER_REACH: "promoter of", curated 0.8
  enhancer-like (pELS, dELS) → the coding genes whose TSS lies in the same
                               domain (between CTCF boundaries), nearest
                               first: "nearest TSS in domain" 0.4, others
                               "same domain" 0.25; all inferred
  CTCF-only                  → insulator, no targets
  DNase-H3K4me3              → open chromatin, nearest TSS within 2 kb, 0.5

Reach is bounded by the domain because enhancer–promoter contacts rarely
cross a TAD boundary; that is the whole point of the node model.
"""

from __future__ import annotations

import bisect
from typing import Any

from genomeos.coords import Locus, Strand
from genomeos.genome.domains import Domain, infer_domains
from genomeos.genome.regulatory import CCRE
from genomeos.genome.regulatory import EVIDENCE as ENCODE_EVIDENCE
from genomeos.ir.model import Evidence, EvidenceKind, RegulatoryElement

PROMOTER_REACH = 1_000
OPEN_CHROMATIN_REACH = 2_000
CLS = {
    "PLS": "promoter",
    "pELS": "enhancer",
    "dELS": "enhancer",
    "CTCF-only": "insulator",
    "DNase-H3K4me3": "open_chromatin",
}


def _tss(g) -> int:
    return g.locus.end - 1 if g.locus.strand is Strand.MINUS else g.locus.start


def assign_targets(
    chrom: str, ccres: list[CCRE], annotation, domains: list[Domain]
) -> list[RegulatoryElement]:
    """Turn ENCODE elements into RegulatoryElement entities with targets."""
    genes = [g for g in annotation.genes.values() if g.locus.chrom == chrom]
    coding = sorted((g for g in genes if g.type == "protein_coding"), key=_tss)
    tss_all = sorted(((_tss(g), g) for g in genes), key=lambda x: x[0])
    tss_pos = [t for t, _ in tss_all]
    dom_starts = [d.start for d in domains]

    def domain_at(pos: int) -> Domain | None:
        i = bisect.bisect_right(dom_starts, pos) - 1
        return domains[i] if 0 <= i < len(domains) and domains[i].start <= pos < domains[i].end else None

    def nearest(pos: int, reach: int) -> tuple[Any, int] | None:
        i = bisect.bisect_left(tss_pos, pos)
        best = None
        for j in (i - 1, i):
            if 0 <= j < len(tss_all):
                d = abs(tss_all[j][0] - pos)
                if d <= reach and (best is None or d < best[1]):
                    best = (tss_all[j][1], d)
        return best

    out: list[RegulatoryElement] = []
    for c in ccres:
        mid = (c.start + c.end) // 2
        cls = CLS.get(c.cls, "unknown")
        dom = domain_at(mid)
        el = RegulatoryElement(
            id=c.id,
            kind="regulatory_element",
            locus=Locus(chrom, c.start, c.end),
            cls=cls,
            domain=dom.id if dom else "",
            source=ENCODE_EVIDENCE,
            evidence=Evidence(EvidenceKind.CURATED, "ENCODE cCRE class"),
            confidence=0.9,
        )
        if cls == "promoter":
            hit = nearest(mid, PROMOTER_REACH)
            if hit:
                el.targets = [
                    {"gene": hit[0].symbol, "distance": hit[1], "basis": "promoter of", "confidence": 0.8}
                ]
        elif cls == "open_chromatin":
            hit = nearest(mid, OPEN_CHROMATIN_REACH)
            if hit:
                el.targets = [
                    {
                        "gene": hit[0].symbol,
                        "distance": hit[1],
                        "basis": "open chromatin at TSS",
                        "confidence": 0.5,
                    }
                ]
        elif cls == "enhancer" and dom:
            lo = bisect.bisect_left([_tss(g) for g in coding], dom.start)
            cands = []
            for g in coding[lo:]:
                t = _tss(g)
                if t >= dom.end:
                    break
                cands.append((abs(t - mid), g))
            cands.sort(key=lambda x: x[0])
            el.targets = [
                {
                    "gene": g.symbol,
                    "distance": d,
                    "basis": "nearest TSS in domain" if k == 0 else "same domain",
                    "confidence": 0.4 if k == 0 else 0.25,
                }
                for k, (d, g) in enumerate(cands[:8])
            ]
        out.append(el)
    return out


def regulation_of(
    symbol: str, chrom: str, ccres: list[CCRE], annotation, chrom_length: int
) -> dict[str, Any]:
    """The regulatory input of one gene: its promoter elements, the enhancers that can
    reach it (same domain), and the insulators that bound the domain."""
    domains = infer_domains(chrom, chrom_length, ccres, annotation)
    elements = assign_targets(chrom, ccres, annotation, domains)
    g = annotation.gene(symbol)
    tss = _tss(g)
    dom = next((d for d in domains if d.start <= tss < d.end), None)
    mine = [e for e in elements if any(t["gene"] == g.symbol for t in e.targets)]
    promoters = [e for e in mine if e.cls in ("promoter", "open_chromatin")]
    enhancers = [e for e in mine if e.cls == "enhancer"]
    nearest = [e for e in enhancers if e.targets[0]["gene"] == g.symbol]
    insulators = [
        e
        for e in elements
        if e.cls == "insulator"
        and dom
        and dom.start - 5000 <= e.locus.start <= dom.end + 5000
        and (abs(e.locus.start - dom.start) < 5000 or abs(e.locus.end - dom.end) < 5000)
    ]
    inside = [e for e in enhancers if g.locus.start <= e.locus.start < g.locus.end]

    def row(e: RegulatoryElement) -> dict[str, Any]:
        t = next(t for t in e.targets if t["gene"] == g.symbol)
        return {
            "id": e.id,
            "class": e.cls,
            "start": e.locus.start,
            "end": e.locus.end,
            "distance": t["distance"],
            "basis": t["basis"],
            "confidence": t["confidence"],
            "intragenic": g.locus.start <= e.locus.start < g.locus.end,
        }

    return {
        "gene": g.symbol,
        "chrom": chrom,
        "tss": tss,
        "strand": g.locus.strand.value,
        "domain": dom.to_dict() if dom else None,
        "promoters": [row(e) for e in promoters],
        "enhancers_in_domain": len(enhancers),
        "enhancers_nearest_to_this_gene": len(nearest),
        "enhancers_inside_gene": len(inside),
        "enhancers": sorted((row(e) for e in enhancers), key=lambda r: r["distance"])[:200],
        "insulators_bounding": [{"id": e.id, "start": e.locus.start, "end": e.locus.end} for e in insulators],
        "competing_genes": dom.genes if dom else [],
        "evidence": {
            "elements": "curated: " + ENCODE_EVIDENCE,
            "targets": "inferred: reach bounded by the CTCF domain, nearest TSS first",
            "domain": "inferred: CTCF-only boundaries, no Hi-C",
        },
    }
