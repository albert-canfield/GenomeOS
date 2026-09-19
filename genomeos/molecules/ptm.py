# SPDX-License-Identifier: AGPL-3.0-or-later
"""Post-translational state: the modified residues UniProt records, and who writes them.

The protein model separates what a protein is (ProteinDefinition) from one protein in one place at
one time (ProteinState), and until now nothing populated the modifications of a state. UniProt's
"Modified residue", "Glycosylation", "Lipidation", "Disulfide bond" and "Cross-link" features carry
the site, the chemistry and, often, the writer ("Phosphoserine; by CDK5, PRPK, AMPK, NUAK1 and ATM").
This module reads them out of the compiled definitions: per protein, the sites by class with their
writers; genome-wide, the writer → substrate edges (a kinase's substrates, an acetyltransferase's) as
a curated layer of the knowledge graph, and the counts as a committed summary.

Curated evidence (UniProt, largely from the literature it cites); a modification is a site that can be
modified, not a measurement that it is modified in a given cell.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path
from typing import Any

CACHE = Path("data/knowledge/proteins")
INDEX = CACHE / "_ptm_index.json"
EVIDENCE = (
    "curated: UniProt sequence features (modified residue, glycosylation, lipidation, cross-link, disulfide)"
)
CLASSES = (
    ("phospho", ("phospho",)),
    ("acetyl", ("acetyl",)),
    ("methyl", ("methyl",)),
    ("ubiquitin", ("ubiquitin", "glycyl lysine isopeptide")),
    ("sumo", ("sumo",)),
    ("glyco", ("glycos", "glycan", "glcnac", "fucos", "mannos", "galactos", "xylos")),
    ("lipid", ("myristoyl", "palmitoyl", "farnesyl", "geranyl", "gpi-anchor", "lipoyl", "octanoyl")),
    ("disulfide", ("disulfide",)),
    ("hydroxy", ("hydroxy",)),
    ("nitro", ("nitro",)),
    ("adp-ribosyl", ("adp-ribosyl",)),
    ("citrulline", ("citrulline",)),
    ("sulfo", ("sulfo",)),
)
_BY = re.compile(r";\s*by\s+([^;.]+)")
_AUTO = re.compile(r"\bautocatal", re.I)


def classify(description: str, feature_type: str) -> str:
    d = description.lower()
    if feature_type == "Disulfide bond":
        return "disulfide"
    if feature_type == "Glycosylation":
        return "glyco"
    if feature_type == "Lipidation":
        return "lipid"
    for cls, keys in CLASSES:
        if any(k in d for k in keys):
            return cls
    return "other"


def writers_of(description: str) -> list[str]:
    """The enzymes named after "by": 'Phosphoserine; by CDK5, PRPK, AMPK, NUAK1 and ATM' → 5 symbols."""
    m = _BY.search(description)
    if not m:
        return []
    text = m.group(1)
    if _AUTO.search(text) or "autocatalysis" in text.lower():
        return ["(self)"]
    parts = re.split(r",|\band\b|/|\bor\b", text)
    out = []
    for p in parts:
        p = p.strip().strip(".")
        p = re.sub(r"\s*\(.*?\)", "", p)
        p = re.sub(r"^(isoform|in vitro|in )\b.*", "", p).strip()
        if p and len(p) <= 15 and re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9\-]*", p):
            out.append(p.upper())
    return out


def sites(defn: dict[str, Any]) -> list[dict[str, Any]]:
    """The modifiable sites of one compiled definition, with class and writers."""
    out = []
    for f in (defn["sections"].get("modifications") or {}).get("items") or []:
        desc = f.get("description") or ""
        out.append(
            {
                "type": f.get("type"),
                "start": f.get("start"),
                "end": f.get("end"),
                "class": classify(desc, f.get("type") or ""),
                "description": desc.split(";")[0].strip(),
                "writers": writers_of(desc),
            }
        )
    return out


def states(defn: dict[str, Any]) -> list[dict[str, Any]]:
    """ProteinState-shaped records: one per modifiable site, the writer as the condition."""
    pid = defn.get("id") or f"gene:{defn['gene']}"
    out = []
    for s in sites(defn):
        out.append(
            {
                "protein": pid,
                "modifications": [f"{s['description']} at {s['start']}"],
                "written_by": s["writers"],
                "class": s["class"],
                "evidence": EVIDENCE,
                "confidence": 0.8,
            }
        )
    return out


def build_index(cache_dir: Path = CACHE, out: Path | None = None) -> dict[str, Any]:
    """Genome-wide: sites per class, writers and their substrates, from every cached definition."""
    per_class: Counter = Counter()
    writers: dict[str, Counter] = {}
    proteins_with = 0
    proteins = 0
    site_total = 0
    per_protein: dict[str, dict[str, int]] = {}
    for p in sorted(cache_dir.glob("*.json")):
        if p.name.startswith("_"):
            continue
        try:
            d = json.loads(p.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if "sections" not in d:
            continue
        proteins += 1
        ss = sites(d)
        if not ss:
            continue
        proteins_with += 1
        site_total += len(ss)
        counts = Counter(s["class"] for s in ss)
        per_class.update(counts)
        per_protein[d["gene"]] = dict(counts)
        for s in ss:
            for w in s["writers"]:
                if w != "(self)":
                    writers.setdefault(w, Counter())[d["gene"]] += 1
    index = {
        "proteins": proteins,
        "proteins_with_sites": proteins_with,
        "sites": site_total,
        "by_class": dict(per_class.most_common()),
        "writers": {w: dict(c) for w, c in writers.items()},
        "per_protein": per_protein,
        "evidence": EVIDENCE,
    }
    out = out or INDEX
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(index))
    return index


def load_index(path: Path = INDEX) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text()) if path.exists() else None
    except (OSError, json.JSONDecodeError):
        return None


def summary(index: dict[str, Any], top: int = 20) -> dict[str, Any]:
    """The committed summary: counts, and the writers with the most substrates."""
    writers = index.get("writers", {})
    ranked = sorted(writers.items(), key=lambda kv: -len(kv[1]))
    return {
        "proteins": index["proteins"],
        "proteins_with_sites": index["proteins_with_sites"],
        "sites": index["sites"],
        "by_class": index["by_class"],
        "writers": len(writers),
        "writer_edges": sum(len(c) for c in writers.values()),
        "top_writers": [
            {"writer": w, "substrates": len(c), "sites": sum(c.values())} for w, c in ranked[:top]
        ],
        "evidence": index["evidence"],
    }
