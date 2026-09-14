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

import gzip
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path
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
    evidence: str = EVIDENCE

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
            "evidence": self.evidence,
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
    chrom: str,
    length: int,
    ccres: list[CCRE],
    annotation=None,
    min_size: int = MIN_DOMAIN,
    orientation: dict[str, set[str]] | None = None,
) -> list[Domain]:
    """Nodes between boundaries.

    By default the boundaries are CTCF-only elements, the caller every committed node rests on,
    and the output is exactly what it has always been. With `orientation` (element id -> motif
    strands, from `ctcf_motif_strands`) the boundaries are the orientation-aware ones of
    `oriented_boundaries` instead, and each domain carries ORIENTED_EVIDENCE."""
    if orientation is not None:
        doms = _domains_from(
            chrom, length, ccres, annotation, min_size, oriented_boundaries(ccres, orientation)
        )
        return _with_evidence(doms, ORIENTED_EVIDENCE)
    return _domains_from(chrom, length, ccres, annotation, min_size, boundaries_from_ccres(ccres))


def _domains_from(
    chrom: str, length: int, ccres: list[CCRE], annotation, min_size: int, boundaries: list[int]
) -> list[Domain]:
    bounds = [0, *boundaries, length]
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


# ------------------------------------------------------------------------------------------
# CTCF motif orientation at the boundaries
# ------------------------------------------------------------------------------------------

#: JASPAR's CTCF core profile. Its forward strand reads GCCACCAGGGGGCGC, the CCACNAGGTGGCAG core
#: that Rao et al. 2014 (Cell 159:1665) call forward: a loop is anchored by a forward (+) site at its
#: upstream end and a reverse (-) site at its downstream end, the convergent rule of loop extrusion.
CTCF_PROFILE = "MA0139"
ORIENTATION_EVIDENCE = (
    "inferred: JASPAR MA0139 best hits inside CTCF-only ENCODE elements; convergent = the domain's "
    "upstream boundary carries a + site and its downstream boundary a - site (Rao et al. 2014)"
)


def boundary_clusters(
    ccres: list[CCRE], merge_within: int = MERGE_BOUNDARIES_WITHIN
) -> list[tuple[int, list[str]]]:
    """`boundaries_from_ccres` with the CTCF-only element ids that were merged into each boundary."""
    sites = sorted(((c.start + c.end) // 2, c.id) for c in ccres if c.cls == "CTCF-only")
    out: list[tuple[int, list[str]]] = []
    for s, cid in sites:
        if out and s - out[-1][0] <= merge_within:
            out[-1] = ((out[-1][0] + s) // 2, [*out[-1][1], cid])
        else:
            out.append((s, [cid]))
    return out


def domain_edges(length: int, boundaries: list[int], min_size: int = MIN_DOMAIN) -> list[int]:
    """The edges `infer_domains` keeps, chromosome ends included."""
    bounds = [0, *boundaries, length]
    edges = [bounds[0]]
    for b in bounds[1:]:
        if b - edges[-1] < min_size and len(edges) > 1:
            continue
        edges.append(b)
    if edges[-1] != length:
        edges[-1] = length
    return edges


def orient_boundaries(
    length: int,
    ccres: list[CCRE],
    strands: dict[str, set[str]],
    min_size: int = MIN_DOMAIN,
) -> list[dict[str, Any]]:
    """Each interior domain edge with the motif strands of its CTCF-only elements and a class.

    `strands` maps a CTCF-only element id to the strands its motif hits lie on ({"+"}, {"-"},
    both, or empty). A boundary between domains k and k+1 anchors k+1 with a + site and k with a
    - site. A domain is convergent when its upstream edge has a + site and its downstream edge a -
    site, divergent when its upstream edge has only - and its downstream edge only +. Classes:
    convergent (an edge of at least one convergent domain), divergent (an edge of a divergent
    domain and of no convergent one), motif (a site, neither), none (no site)."""
    clusters = boundary_clusters(ccres)
    edges = domain_edges(length, [p for p, _ in clusters], min_size)
    by_pos: dict[int, list[str]] = {}
    for p, ids in clusters:
        by_pos.setdefault(p, []).extend(ids)
    rows = []
    for e in edges:
        ids = by_pos.get(e, []) if 0 < e < length else []
        st: set[str] = set()
        for i in ids:
            st |= strands.get(i, set())
        rows.append({"position": e, "elements": ids, "plus": "+" in st, "minus": "-" in st})
    conv = [False] * len(rows)
    div = [False] * len(rows)
    for k in range(len(rows) - 1):  # domain k spans edges k .. k+1
        up, down = rows[k], rows[k + 1]
        if up["plus"] and down["minus"]:
            conv[k] = conv[k + 1] = True
        elif up["minus"] and not up["plus"] and down["plus"] and not down["minus"]:
            div[k] = div[k + 1] = True
    out = []
    for k, r in enumerate(rows):
        if not 0 < r["position"] < length:
            continue
        has = r["plus"] or r["minus"]
        cls = "convergent" if conv[k] else "divergent" if div[k] else "motif" if has else "none"
        sites = "both" if r["plus"] and r["minus"] else "+" if r["plus"] else "-" if r["minus"] else "none"
        out.append({**r, "class": cls, "sites": sites})
    return out


# ------------------------------------------------------------------------------------------
# The orientation-aware caller: an alternative node set beside the CTCF-only one
# ------------------------------------------------------------------------------------------

ORIENTED_EVIDENCE = (
    "inferred: boundaries where consecutive CTCF motifs (JASPAR MA0139, single-strand hits inside ENCODE "
    "elements with CTCF ChIP support) flip from reverse to forward, the loop-extrusion rule; no Hi-C"
)
STRAND_CACHE = Path("data/knowledge/ctcf_strands")


def ctcf_site_elements(ccres: list[CCRE]) -> list[CCRE]:
    """Elements with CTCF ChIP support: the CTCF-only class and any class flagged CTCF-bound.

    The CTCF-only class alone leaves out every CTCF site that overlaps a promoter or an enhancer;
    the registry flags those as CTCF-bound, and an extrusion barrier does not care what else binds."""
    return [c for c in ccres if c.cls == "CTCF-only" or c.ctcf_bound]


def ctcf_motif_strands(
    chrom: str,
    ccres: list[CCRE],
    fetch: Callable[[int, int], str],
    relative: float | None = None,
    cache: Path | None = STRAND_CACHE,
    tag: str = "",
) -> dict[str, set[str]]:
    """Element id -> the strands of its MA0139 hits, for every element with CTCF ChIP support.

    `fetch(start, end)` returns the forward sequence. Cached per chromosome under
    data/knowledge/ctcf_strands (a cache, git-ignored) because a genome's worth of scanning takes
    minutes."""
    from genomeos.genome import motifs as mo

    rel = mo.RELATIVE_SCORE if relative is None else relative
    p = cache / f"{tag}{chrom}_{rel}.tsv.gz" if cache is not None else None
    if p is not None and p.exists():
        out: dict[str, set[str]] = {}
        with gzip.open(p, "rt") as fh:
            for line in fh:
                cid, st = line.rstrip("\n").split("\t")
                out[cid] = set(st) if st != "." else set()
        return out
    profiles = [
        m for m in mo.load_motifs(relative=min(rel, mo.RELATIVE_SCORE)) if m.id.startswith(CTCF_PROFILE)
    ]
    out = {}
    for c in ctcf_site_elements(ccres):
        hits = mo.all_hits(profiles, fetch(c.start, c.end), min_relative=rel)
        out[c.id] = {"+" if h[3] == 0 else "-" for h in hits}
    if p is not None:
        p.parent.mkdir(parents=True, exist_ok=True)
        part = p.with_name(p.name + ".part")
        with gzip.open(part, "wt") as fh:
            for cid, st in out.items():
                fh.write(f"{cid}\t{''.join(sorted(st)) or '.'}\n")
        part.replace(p)
    return out


def oriented_boundaries(
    ccres: list[CCRE], strands: dict[str, set[str]], merge_within: int = MERGE_BOUNDARIES_WITHIN
) -> list[int]:
    """Boundary positions where consecutive single-strand CTCF sites flip from reverse to forward.

    Sites are the elements of `ctcf_site_elements` whose motif hits lie on one strand; an element
    with hits on both strands cannot say which way it points and is skipped. A boundary sits midway
    between a reverse site and the forward site that follows it (the pair points away from each
    other, as sites flanking measured boundaries do); boundaries closer than `merge_within` are
    merged as `boundaries_from_ccres` merges CTCF-only sites."""
    sites = sorted(
        ((c.start + c.end) // 2, next(iter(strands[c.id])))
        for c in ctcf_site_elements(ccres)
        if len(strands.get(c.id, ())) == 1
    )
    flips = [
        (a + b) // 2 for (a, sa), (b, sb) in zip(sites, sites[1:], strict=False) if sa == "-" and sb == "+"
    ]
    out: list[int] = []
    for f in flips:
        if out and f - out[-1] <= merge_within:
            out[-1] = (out[-1] + f) // 2
        else:
            out.append(f)
    return out


def _with_evidence(domains: list[Domain], evidence: str) -> list[Domain]:
    for d in domains:
        d.evidence = evidence
    return domains
