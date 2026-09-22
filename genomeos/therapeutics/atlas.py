"""Healthy-tissue RNA for the whole surface universe, shipped with the package.

Every candidate the pipeline scores needs the same question answered: how much
of this protein does a healthy person already carry, and where. Until now that
was one Human Protein Atlas request per gene, which meant the offline path
returned "no cached HPA record", normal-tissue safety was unknown, and the
score was capped for missing safety data on every candidate. It also meant
that scanning for a target rather than starting from an altered gene was not
possible at all: five thousand requests is not a scan.

The Atlas answers a whole protein class in one request, so the class that
matters here — the 5,573 genes HPA predicts to be membrane proteins, which is
the universe a circulating binder could ever reach — is fetched once, distilled
to the tissues GenomeOS weighs, and shipped gzipped inside the package. The raw
download is never kept. What is committed is the small summary under
data/results; what is shipped is the table.

The rows keep the Atlas's own column names, so `normal_profile` reads a local
row and a freshly fetched one through the same code and cannot tell them apart.

It is a population-level healthy-tissue measurement and never a patient's
tumour, it covers 20 tissues and not the body, and a gene outside the membrane
class is absent from the table rather than absent from the body: all three are
carried in the summary and in the evidence of every answer built from it.
"""

from __future__ import annotations

import gzip
import json
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

HPA_API = "https://www.proteinatlas.org/api/search_download.php"

#: The protein classes fetched. The membrane class is the universe a
#: circulating binder could reach; the CD markers are added because that is
#: where the antigens of the approved cell therapies live (CD19, TNFRSF17),
#: and a few of them are not in the membrane class.
QUERIES: tuple[str, ...] = (
    "protein_class:Predicted membrane proteins",
    "protein_class:CD markers",
)

PACKAGED = Path(__file__).resolve().parent / "data" / "normal_tissue.json.gz"
RESULT = "normal_tissue_atlas"

#: Columns kept from each Atlas row, beyond the per-tissue nTPM values.
KEEP = (
    "Gene",
    "Ensembl",
    "RNA tissue specificity",
    "RNA tissue distribution",
    "Subcellular main location",
    "Subcellular location",
)

#: Protein-class flags worth a boolean each; the full class list is dropped.
CLASS_FLAGS = {
    "Predicted membrane proteins": "membrane",
    "CD markers": "cd_marker",
    "FDA approved drug targets": "approved_drug_target",
    "Cancer-related genes": "cancer_related",
}


def _fetch(query: str, columns: str, timeout: int = 300) -> list[dict[str, Any]]:
    url = f"{HPA_API}?search={urllib.parse.quote(query)}&format=json&columns={columns}&compress=no"
    req = urllib.request.Request(url, headers={"User-Agent": "GenomeOS/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        return json.load(r)


def _distil_row(row: dict[str, Any], tissue_keys: tuple[str, ...]) -> dict[str, Any]:
    """One Atlas row, kept to what the pipeline reads, with its own column names."""
    out: dict[str, Any] = {k: row.get(k) for k in KEEP if row.get(k) is not None}
    classes = set(row.get("Protein class") or [])
    for name, flag in CLASS_FLAGS.items():
        if name in classes:
            out[flag] = True
    for key in tissue_keys:
        raw = row.get(key)
        if raw in (None, ""):
            continue
        try:
            out[key] = round(float(raw), 1)
        except (TypeError, ValueError):
            continue
    return out


def distil(tissues: dict[str, tuple[str, float]], queries: tuple[str, ...] = QUERIES, log=None) -> dict:
    """Fetch the protein classes, keep the tissues GenomeOS weighs, discard the rest."""
    from .providers import hpa_column

    columns = ",".join(["g", "eg", "pc", "scl", "scml", "rnats", "rnatd", *tissues])
    tissue_keys = tuple(hpa_column(field) for field in tissues)
    genes: dict[str, dict[str, Any]] = {}
    per_query: dict[str, int] = {}
    for q in queries:
        rows = _fetch(q, columns)
        per_query[q] = len(rows)
        if log:
            print(f"  {q}: {len(rows)} genes", file=log, flush=True)
        for row in rows:
            symbol = (row.get("Gene") or "").upper()
            if symbol and symbol not in genes:
                genes[symbol] = _distil_row(row, tissue_keys)
        time.sleep(0.5)
    measured = sum(1 for g in genes.values() if any(k in g for k in tissue_keys))
    return {
        "genes": genes,
        "meta": {
            "source": "Human Protein Atlas consensus tissue RNA (proteinatlas.org), CC BY-SA 4.0",
            "queries": list(queries),
            "genes_per_query": per_query,
            "tissues": [name for name, _w in tissues.values()],
            "genes": len(genes),
            "genes_with_tissue_values": measured,
            "level": "population-level healthy tissue; never a patient's tumour",
            "limits": [
                "20 tissues, not the body: a tissue not queried is not a tissue where the gene is absent",
                "a gene outside the membrane and CD-marker classes is absent from this table, which "
                "is not evidence that it is absent from healthy tissue",
                "RNA is not protein, and nTPM across tissues is not receptors per cell",
            ],
        },
    }


def write(table: dict, path: Path = PACKAGED) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        json.dump(table, fh, separators=(",", ":"))
    return path


_CACHE: dict[str, Any] | None = None


def load(path: Path = PACKAGED) -> dict[str, Any] | None:
    """The shipped table, read once per process."""
    global _CACHE  # noqa: PLW0603 - one process-wide read of a packaged file
    if _CACHE is None:
        if not path.exists():
            return None
        with gzip.open(path, "rt") as fh:
            _CACHE = json.load(fh)
    return _CACHE


def row(gene: str, path: Path = PACKAGED) -> dict[str, Any] | None:
    """One gene's healthy-tissue row, in the Atlas's own column names, or None."""
    table = load(path)
    if not table:
        return None
    return (table.get("genes") or {}).get(gene.upper())


def universe(path: Path = PACKAGED) -> list[str]:
    """Every gene in the table: the surface universe a binder could reach."""
    table = load(path)
    return sorted((table or {}).get("genes") or {})


def summary(table: dict) -> dict[str, Any]:
    """The committed record of a distillation run: counts, not the table."""
    meta = dict(table["meta"])
    genes = table["genes"]
    flagged = {
        flag: sum(1 for g in genes.values() if g.get(flag)) for flag in sorted(set(CLASS_FLAGS.values()))
    }
    with_specificity = sum(1 for g in genes.values() if g.get("RNA tissue specificity"))
    return {
        "result": RESULT,
        **meta,
        "protein_class_flags": flagged,
        "genes_with_specificity_class": with_specificity,
        "examples": {
            g: {
                "specificity": genes[g].get("RNA tissue specificity"),
                "highest_queried_tissue": _highest(genes[g]),
            }
            for g in ("ERBB2", "EGFR", "CD19", "TNFRSF17", "MSLN", "FOLH1")
            if g in genes
        },
        "note": (
            "The table itself ships gzipped inside the package "
            "(genomeos/therapeutics/data/normal_tissue.json.gz); the raw download is discarded."
        ),
    }


def _highest(gene_row: dict[str, Any]) -> dict[str, Any] | None:
    values = {
        k[len("Tissue RNA - ") : -len(" [nTPM]")]: v
        for k, v in gene_row.items()
        if k.startswith("Tissue RNA - ") and isinstance(v, int | float)
    }
    if not values:
        return None
    tissue = max(values, key=lambda t: values[t])
    return {"tissue": tissue, "nTPM": values[tissue]}
