# SPDX-License-Identifier: AGPL-3.0-or-later
"""The proteome as a library: the compiled human proteome distilled into one packaged table.

The full compiled definitions (about 290 MB, one JSON per gene, rebuildable with the proteome
job) stay local. This module keeps what a program needs about every protein in a few
megabytes shipped inside GenomeOS: accession, length, one-line function, location, domains,
pathways, physically supported partners, tissue pattern, diseases, structures, and the
confidence per field. `distil()` builds it from the local cache; `load()` reads the packaged
copy; `block()` writes one protein as a BioLang block on demand.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

PACKAGED = Path(__file__).parent / "data" / "proteome.json.gz"
FIELDS = (
    "accession",
    "name",
    "length",
    "existence",
    "function",
    "location",
    "domains",
    "pathways",
    "partners",
    "writers",
    "tissue_pattern",
    "cell_type_pattern",
    "diseases",
    "structures_experimental",
    "alphafold",
    "symbol_match",
    "chrom",
)
EVIDENCE = {
    "accession": "curated: UniProtKB/Swiss-Prot",
    "function": "curated: UniProt",
    "domains": "curated: InterPro",
    "pathways": "curated: Reactome",
    "partners": "predicted: STRING, experimental channel ≥ 0.4",
    "writers": "curated: UniProt modified residues, 'by' the named enzyme",
    "tissue_pattern": "experimental: Human Protein Atlas",
    "diseases": "curated: UniProt",
    "structures_experimental": "experimental: PDB",
    "alphafold": "predicted: AlphaFold DB",
}


def _row(d: dict[str, Any], chrom_of: dict[str, str]) -> dict[str, Any] | None:
    from genomeos.molecules.compiler import writer_counts

    s = d.get("sections", {})
    ident = (s.get("identity") or {}).get("items") or {}
    if not ident.get("accession"):
        return None
    fn = (s.get("function") or {}).get("items") or {}
    dom = (s.get("domains") or {}).get("items") or {}
    ex = (s.get("expression") or {}).get("items") or {}
    inter = (s.get("interactions") or {}).get("items") or []
    summary = (fn.get("summary") or [""])[0]
    return {
        "accession": ident["accession"],
        "name": ident.get("name"),
        "length": ident.get("length"),
        "existence": (ident.get("existence") or "").split(":")[0].strip() or None,
        "function": summary[:240] + ("…" if len(summary) > 240 else ""),
        "location": (fn.get("location") or [])[:3],
        "domains": [x["name"] for x in dom.get("interpro", []) if x.get("name")][:10],
        "pathways": [x["id"] for x in (s.get("pathways") or {}).get("items") or []][:12],
        "partners": list(dict.fromkeys(x["partner"] for x in inter if x.get("physical_evidence")))[:10],
        "writers": dict(sorted(writer_counts(d).items(), key=lambda kv: (-kv[1], kv[0]))[:12]),
        "tissue_pattern": ex.get("tissue_specificity"),
        "cell_type_pattern": ex.get("cell_type_specificity"),
        "diseases": [x["name"] for x in (s.get("diseases") or {}).get("items") or [] if x.get("name")][:5],
        "structures_experimental": len((s.get("structures_experimental") or {}).get("items") or []),
        "alphafold": bool((s.get("structures_predicted") or {}).get("items")),
        "symbol_match": ident.get("symbol_match", True),
        "chrom": chrom_of.get(d["gene"]),
    }


def distil(cache_dir: Path = Path("data/knowledge/proteins"), out: Path = PACKAGED) -> dict[str, Any]:
    """Build the packaged table from the local cache; returns a summary."""
    chrom_of: dict[str, str] = {}
    for p in Path("data/results").glob("proteome_chr*.json"):
        d = json.loads(p.read_text())
        for sym in d.get("per_gene", {}):
            chrom_of[sym] = d["chrom"]
    table: dict[str, dict[str, Any]] = {}
    for p in sorted(cache_dir.glob("*.json")):
        try:
            d = json.loads(p.read_text())
        except json.JSONDecodeError:
            continue
        row = _row(d, chrom_of)
        if row:
            table[d["gene"]] = row
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "proteins": table,
        "fields": FIELDS,
        "evidence": EVIDENCE,
        "note": "distilled from the compiled definitions; the full record is the local knowledge cache",
    }
    with gzip.open(out, "wt") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    return {"proteins": len(table), "bytes": out.stat().st_size, "path": str(out)}


_LOADED: dict[str, Any] = {}


def load(path: Path = PACKAGED) -> dict[str, dict[str, Any]]:
    if "table" not in _LOADED:
        if not path.exists():
            return {}
        with gzip.open(path, "rt") as fh:
            _LOADED["table"] = json.load(fh)["proteins"]
    return _LOADED["table"]


def get(symbol: str) -> dict[str, Any] | None:
    return load().get(symbol.upper())


def block(symbol: str) -> str:
    """One BioLang `protein` block from the packaged table (no network, no cache needed)."""
    r = get(symbol)
    if not r:
        return f"# {symbol.upper()}: not in the packaged proteome\n"
    props = [f"accession: {r['accession']}"]
    if r["domains"]:
        props.append(f"domains: {', '.join(r['domains'])}")
    if r["pathways"]:
        props.append(f"pathways: {', '.join(r['pathways'])}")
    if r["partners"]:
        props.append(f"interactions: {', '.join(r['partners'])}")
    props.append(
        'evidence: curated "UniProt; InterPro; Reactome; STRING physical channel (packaged proteome)"'
    )
    props.append(f"confidence: {0.9 if r['symbol_match'] else 0.5}")
    head = f"# {r['name']} ({r['length']} aa) — {r['function'][:110]}" if r["function"] else f"# {r['name']}"
    from genomeos.molecules.compiler import writer_rules

    rules = writer_rules(symbol.upper(), r.get("writers") or {})
    body = head + "\n" + f"protein {symbol.upper()} {{\n" + "".join(f"  {x};\n" for x in props) + "}\n"
    return body + "".join(x + "\n" for x in rules)


def summary() -> dict[str, Any]:
    t = load()
    n = len(t)
    if not n:
        return {"proteins": 0}
    return {
        "proteins": n,
        "with_function": sum(1 for r in t.values() if r["function"]),
        "with_pathways": sum(1 for r in t.values() if r["pathways"]),
        "with_partners": sum(1 for r in t.values() if r["partners"]),
        "with_experimental_structure": sum(1 for r in t.values() if r["structures_experimental"]),
        "with_disease": sum(1 for r in t.values() if r["diseases"]),
        "symbol_mismatch": sum(1 for r in t.values() if not r["symbol_match"]),
        "evidence": EVIDENCE,
    }
