"""Nodes above genes: domains inferred from CTCF boundaries.

The genome folds into topologically associating domains (TADs) whose
boundaries are marked by CTCF sites held by cohesin. Without Hi-C data we
approximate boundaries with ENCODE's CTCF-only elements: sites where CTCF
binds and nothing else does, which is what insulators look like in chromatin
data. Domains are the intervals between boundaries, merged below a minimum
size. Each domain lists the genes (by transcription start) and the
regulatory elements it contains. Evidence: inferred, confidence 0.4; real
domain calls need Hi-C and would come in as curated.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from statistics import median
from typing import Any

from genomeos.genome.regulatory import CCRE

EVIDENCE = (
    "inferred: CTCF-only ENCODE elements used as boundary proxies "
    "(TAD boundaries are CTCF/cohesin sites); no Hi-C"
)
MIN_DOMAIN = 50_000
MERGE_BOUNDARIES_WITHIN = 5_000


@dataclass(slots=True)
class Domain:
    id: str
    chrom: str
    start: int
    end: int
    genes: list[str] = field(default_factory=list)
    coding_genes: int = 0
    promoters: int = 0
    enhancers: int = 0
    ctcf_inside: int = 0
    confidence: float = 0.4

    @property
    def length(self) -> int:
        return self.end - self.start

    def to_dict(self, compact: bool = False) -> dict[str, Any]:
        """`compact` keeps the first 12 gene symbols and a count (the saved result stays small)."""
        genes = self.genes[:12] if compact else self.genes
        return {
            "id": self.id,
            "chrom": self.chrom,
            "start": self.start,
            "end": self.end,
            "length": self.length,
            "genes": genes,
            "genes_count": len(self.genes),
            "coding_genes": self.coding_genes,
            "promoters": self.promoters,
            "enhancers": self.enhancers,
            "ctcf_inside": self.ctcf_inside,
            "evidence": EVIDENCE,
            "confidence": self.confidence,
        }


def boundaries_from_ccres(ccres: list[CCRE], merge_within: int = MERGE_BOUNDARIES_WITHIN) -> list[int]:
    sites = sorted((c.start + c.end) // 2 for c in ccres if c.cls == "CTCF-only")
    out: list[int] = []
    for s in sites:
        if out and s - out[-1] <= merge_within:
            out[-1] = (out[-1] + s) // 2
        else:
            out.append(s)
    return out


def infer_domains(
    chrom: str, length: int, ccres: list[CCRE], annotation=None, min_size: int = MIN_DOMAIN
) -> list[Domain]:
    bounds = [0, *boundaries_from_ccres(ccres), length]
    # merge intervals shorter than min_size into their neighbour
    edges = [bounds[0]]
    for b in bounds[1:]:
        if b - edges[-1] < min_size and len(edges) > 1:
            continue
        edges.append(b)
    if edges[-1] != length:
        edges[-1] = length
    domains = [
        Domain(f"{chrom}:D{i + 1}", chrom, a, b)
        for i, (a, b) in enumerate(zip(edges, edges[1:], strict=False))
    ]
    # assign elements
    import bisect

    starts = [d.start for d in domains]

    def dom_at(pos: int) -> Domain | None:
        i = bisect.bisect_right(starts, pos) - 1
        return domains[i] if 0 <= i < len(domains) and domains[i].start <= pos < domains[i].end else None

    for c in ccres:
        d = dom_at((c.start + c.end) // 2)
        if not d:
            continue
        if c.cls == "PLS":
            d.promoters += 1
        elif c.cls in ("pELS", "dELS"):
            d.enhancers += 1
        elif c.cls == "CTCF-only":
            d.ctcf_inside += 1
    if annotation is not None:
        for g in annotation.genes.values():
            if g.locus.chrom != chrom:
                continue
            tss = g.locus.end - 1 if g.locus.strand.value == "-" else g.locus.start
            d = dom_at(tss)
            if d:
                d.genes.append(g.symbol)
                if g.type == "protein_coding":
                    d.coding_genes += 1
    return domains


def summarise(domains: list[Domain]) -> dict[str, Any]:
    sizes = [d.length for d in domains]
    genes = [d.coding_genes for d in domains]
    return {
        "domains": len(domains),
        "size_median": int(median(sizes)) if sizes else None,
        "size_max": max(sizes) if sizes else None,
        "coding_genes_per_domain_median": median(genes) if genes else None,
        "domains_without_coding_genes": sum(1 for g in genes if g == 0),
        "largest_gene_count": max(genes) if genes else None,
        "enhancers_per_domain_median": median(d.enhancers for d in domains) if domains else None,
        "evidence": EVIDENCE,
        "confidence": 0.4,
    }
