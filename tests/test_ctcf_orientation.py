"""CTCF motif orientation at inferred node boundaries (genome/domains.py)."""

from __future__ import annotations

from genomeos.genome.domains import (
    boundaries_from_ccres,
    boundary_clusters,
    domain_edges,
    infer_domains,
    orient_boundaries,
)
from genomeos.genome.regulatory import CCRE


def _ctcf(pos: int, name: str) -> CCRE:
    return CCRE("chrT", pos - 100, pos + 100, name, "CTCF-only", True)


def test_clusters_and_edges_reproduce_infer_domains():
    ccres = [_ctcf(p, f"E{i}") for i, p in enumerate([100_000, 103_000, 200_000, 220_000, 400_000, 600_000])]
    clusters = boundary_clusters(ccres)
    assert [p for p, _ in clusters] == boundaries_from_ccres(ccres)
    assert clusters[0][1] == ["E0", "E1"]  # merged within 5 kb
    edges = domain_edges(800_000, [p for p, _ in clusters])
    doms = infer_domains("chrT", 800_000, ccres)
    assert edges == [doms[0].start, *[d.end for d in doms]]
    assert 220_000 not in edges  # closer than the minimum domain to the previous edge


def test_convergent_divergent_and_unmarked_edges():
    # edges at 200k, 400k, 600k, 800k on a 1 Mb chromosome
    ccres = [_ctcf(p, f"B{p // 1000}") for p in (200_000, 400_000, 600_000, 800_000)]
    strands = {"B200": {"+"}, "B400": {"-"}, "B600": {"-"}, "B800": {"+"}}
    rows = {r["position"]: r for r in orient_boundaries(1_000_000, ccres, strands)}
    # domain 200k-400k: + then - -> convergent; domain 600k-800k: - then + -> divergent
    assert rows[200_000]["class"] == "convergent" and rows[400_000]["class"] == "convergent"
    assert rows[600_000]["class"] == "divergent" and rows[800_000]["class"] == "divergent"
    strands["B600"] = set()
    rows = {r["position"]: r for r in orient_boundaries(1_000_000, ccres, strands)}
    assert rows[600_000]["class"] == "none" and rows[800_000]["class"] == "motif"
    assert rows[200_000]["sites"] == "+" and rows[600_000]["sites"] == "none"
    both = {**strands, "B400": {"+", "-"}}
    rows = {r["position"]: r for r in orient_boundaries(1_000_000, ccres, both)}
    assert rows[400_000]["sites"] == "both" and rows[400_000]["class"] == "convergent"
