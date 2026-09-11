"""The protein knowledge graph, built from compiled definitions (no network).

Nodes are proteins (by gene symbol), pathways, tissues and domains; edges
carry the relation and the evidence of the section they came from:

    protein —[associates, STRING score, physical or not]— protein
    protein —[member of]— pathway            (Reactome, curated)
    protein —[expressed in, nTPM]— tissue    (HPA, experimental)
    protein —[has domain]— InterPro entry    (curated)

`build()` reads every cached definition under data/knowledge/proteins and
returns the graph; `neighbourhood()` cuts the part around one protein for
the Molecules tab; `summarise()` is the distilled result kept per run.
"""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from genomeos.molecules.compiler import CACHE


class Graph:
    def __init__(self) -> None:
        self.nodes: dict[str, dict[str, Any]] = {}
        self.edges: list[dict[str, Any]] = []
        self._adj: dict[str, list[int]] = defaultdict(list)

    def add_node(self, nid: str, kind: str, **attrs: Any) -> None:
        if nid not in self.nodes:
            self.nodes[nid] = {"id": nid, "kind": kind, **attrs}

    def add_edge(self, a: str, b: str, rel: str, evidence: str, confidence: float, **attrs: Any) -> None:
        self.edges.append(
            {"a": a, "b": b, "rel": rel, "evidence": evidence, "confidence": confidence, **attrs}
        )
        i = len(self.edges) - 1
        self._adj[a].append(i)
        self._adj[b].append(i)

    def degree(self, nid: str, rel: str | None = None) -> int:
        return sum(1 for i in self._adj.get(nid, []) if rel is None or self.edges[i]["rel"] == rel)

    def neighbourhood(self, nid: str, max_nodes: int = 60) -> dict[str, Any]:
        """The node, its direct neighbours, and edges among them (ranked by confidence)."""
        if nid not in self.nodes:
            return {"nodes": [], "edges": [], "centre": nid}
        idx = sorted(self._adj.get(nid, []), key=lambda i: -self.edges[i]["confidence"])
        keep = {nid}
        for i in idx:
            e = self.edges[i]
            keep.add(e["b"] if e["a"] == nid else e["a"])
            if len(keep) >= max_nodes:
                break
        edges = [e for e in self.edges if e["a"] in keep and e["b"] in keep]
        return {"centre": nid, "nodes": [self.nodes[n] for n in keep], "edges": edges}


def _load_definitions(cache_dir: Path = CACHE) -> list[dict[str, Any]]:
    out = []
    for p in sorted(cache_dir.glob("*.json")):
        try:
            d = json.loads(p.read_text())
        except json.JSONDecodeError:
            continue
        if (d.get("sections", {}).get("identity") or {}).get("items"):
            out.append(d)
    return out


def build(cache_dir: Path = CACHE, min_score: float = 0.7) -> Graph:
    g = Graph()
    defs = _load_definitions(cache_dir)
    compiled = {d["gene"] for d in defs}
    for d in defs:
        s = d["sections"]
        ident = s["identity"]["items"]
        gene = d["gene"]
        g.add_node(
            gene,
            "protein",
            accession=ident.get("accession"),
            name=ident.get("name"),
            length=ident.get("length"),
            compiled=True,
        )
        for pw in (s.get("pathways") or {}).get("items") or []:
            g.add_node(pw["id"], "pathway", name=pw.get("name"))
            g.add_edge(gene, pw["id"], "member_of", "curated: Reactome", 0.85)
        for dom in ((s.get("domains") or {}).get("items") or {}).get("interpro", []):
            g.add_node(dom["id"], "domain", name=dom.get("name"))
            g.add_edge(gene, dom["id"], "has_domain", "curated: InterPro", 0.9)
        ex = (s.get("expression") or {}).get("items") or {}
        for tissue, v in (ex.get("tissue_ntpm") or {}).items():
            if isinstance(v, int | float):
                g.add_node(f"tissue:{tissue}", "tissue", name=tissue)
                g.add_edge(gene, f"tissue:{tissue}", "expressed_in", "experimental: HPA nTPM", 0.8, ntpm=v)
        for it in (s.get("interactions") or {}).get("items") or []:
            if it["score"] < min_score:
                continue
            partner = it["partner"]
            g.add_node(partner, "protein", compiled=partner in compiled)
            # one edge per unordered pair
            if partner in compiled and partner < gene:
                continue
            g.add_edge(
                gene,
                partner,
                "associates",
                "predicted: STRING" + (" (experimental channel)" if it["physical_evidence"] else ""),
                0.7 if it["physical_evidence"] else 0.5,
                score=it["score"],
                physical=it["physical_evidence"],
            )
    return g


def summarise(g: Graph) -> dict[str, Any]:
    kinds = Counter(n["kind"] for n in g.nodes.values())
    rels = Counter(e["rel"] for e in g.edges)
    proteins = [n for n in g.nodes.values() if n["kind"] == "protein" and n.get("compiled")]
    deg = sorted(((g.degree(n["id"], "associates"), n["id"]) for n in proteins), reverse=True)
    pw_size = Counter(e["b"] for e in g.edges if e["rel"] == "member_of")
    biggest = [(g.nodes[p].get("name"), c) for p, c in pw_size.most_common(8)]
    # connected components over association edges among compiled proteins
    parent = {n["id"]: n["id"] for n in proteins}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for e in g.edges:
        if e["rel"] == "associates" and e["a"] in parent and e["b"] in parent:
            parent[find(e["a"])] = find(e["b"])
    comps = Counter(find(p["id"]) for p in proteins)
    return {
        "nodes": len(g.nodes),
        "edges": len(g.edges),
        "node_kinds": dict(kinds),
        "edge_kinds": dict(rels),
        "compiled_proteins": len(proteins),
        "physical_associations": sum(1 for e in g.edges if e["rel"] == "associates" and e.get("physical")),
        "hubs": [{"gene": n, "associations": d} for d, n in deg[:10]],
        "isolated_proteins": sum(1 for d, _ in deg if d == 0),
        "largest_component": max(comps.values()) if comps else 0,
        "components": len(comps),
        "biggest_pathways": biggest,
        "evidence": "edges carry the evidence of their source section: Reactome, InterPro (curated), "
        "HPA (experimental), STRING (predicted)",
    }


_CACHE: dict[str, Any] = {}


def cached_build(cache_dir: Path = CACHE, min_score: float = 0.7) -> Graph:
    """The graph over the whole local proteome, rebuilt only when the definition cache changed
    (19,000 definitions take about 2.5 s to load; a request should not pay that every time)."""
    files = list(cache_dir.glob("*.json"))
    sig = (len(files), max((f.stat().st_mtime for f in files), default=0.0), min_score)
    if _CACHE.get("sig") != sig:
        _CACHE["graph"] = build(cache_dir, min_score)
        _CACHE["sig"] = sig
    return _CACHE["graph"]
