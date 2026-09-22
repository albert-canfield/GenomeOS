# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""DepMap, streamed and distilled: never a matrix on disk, only the columns asked for.

The Cancer Dependency Map is the largest public statement of what a cancer cell
cannot live without: a genome-wide CRISPR knockout screen in over a thousand
cell lines, with the same lines genotyped for mutations, copy number and
expression. It is also large enough that downloading it is a mistake — the
release is tens of gigabytes and a single release directory would not fit next
to the rest of this project.

So nothing is kept. A file is streamed over HTTP, the header is read once, the
column indices of the genes actually asked for are computed, and every row is
discarded after those few numbers are taken out of it. What lands on disk is a
small per-dataset cache under `data/knowledge/therapeutics/depmap/`, gzipped,
holding the requested genes and nothing else, with the release, the file, the
figshare id and the licence recorded beside the numbers.

The portal's own API is behind a bot check, so the files are addressed by their
figshare ids, which is where the release is published and is the citable form.

Four facts are read, each in its own file and each kept apart:

| Dataset | File | Says |
|---|---|---|
| gene effect | `CRISPRGeneEffect.csv` | Chronos gene effect, 0 = no effect, -1 = median common essential |
| damaging mutation | `OmicsSomaticMutationsMatrixDamaging.csv` | 0/1/2 per gene per line |
| absolute copy number | `OmicsAbsoluteCNGene.csv` | PureCN absolute copies, so 0 is a homozygous deletion |
| signatures | `OmicsSignatures.csv` | MSIsensor2 MSI score, ploidy, LoH fraction, WGD, CIN, aneuploidy |

A fifth, `OmicsExpressionProteinCodingGenesTPMLogp1.csv`, carries log2(TPM+1)
per gene per line and is the tumour side of the marker-selectivity work.

Cell lines are not patients. They carry the alteration and the dependency in
the same cells, which is what makes a dependency test possible at all, and they
have no stroma, no immune system and no tissue architecture. Every result built
on this module says so.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import urllib.request
from pathlib import Path
from typing import Any

FIGSHARE_FILE = "https://ndownloader.figshare.com/files/{}"
CACHE = Path("data/knowledge/therapeutics/depmap")
UA = {"User-Agent": "GenomeOS/0.1 (therapeutics)"}

#: The release every number in this module comes from. Pinned, not "latest":
#: a dependency test that moves when the release moves is not reproducible.
RELEASE = "DepMap 24Q4 Public"
DOI = "10.25452/figshare.plus.27993248.v1"
ARTICLE = "27993248"
LICENCE = "DepMap Public 24Q4, Broad Institute, CC BY 4.0"
CITATION = "DepMap, Broad (2024). DepMap 24Q4 Public. Figshare+. doi:10.25452/figshare.plus.27993248.v1"

#: figshare file ids within that release, with what each one is.
FILES: dict[str, tuple[str, str, str]] = {
    "gene_effect": (
        "51064667",
        "CRISPRGeneEffect.csv",
        "Chronos gene effect, copy-number corrected and screen-quality corrected",
    ),
    "damaging": (
        "51065747",
        "OmicsSomaticMutationsMatrixDamaging.csv",
        "1 or 2 where a gene carries at least one likely loss-of-function mutation",
    ),
    "absolute_cn": (
        "51065303",
        "OmicsAbsoluteCNGene.csv",
        "PureCN absolute copy number per gene; segmental-duplication genes are masked",
    ),
    "signatures": (
        "51065726",
        "OmicsSignatures.csv",
        "model-level MSIsensor2 MSI score, ploidy, LoH fraction, WGD, CIN, aneuploidy",
    ),
    "model": (
        "51065297",
        "Model.csv",
        "one row per cell line: lineage, Oncotree disease, origin",
    ),
    "expression": (
        "51065489",
        "OmicsExpressionProteinCodingGenesTPMLogp1.csv",
        "log2(TPM+1) per protein-coding gene per line, RSEM unstranded",
    ),
}

LIMITS: tuple[str, ...] = (
    "cell lines, not patient tumours: no stroma, no immune system, no tissue architecture, "
    "and a growth-selected subset of the disease",
    "a CRISPR gene effect is a measurement of proliferation in culture over days to weeks, "
    "not of what a drug would do in a patient",
    "the damaging-mutation matrix is a call about a variant's predicted consequence, not a "
    "measured loss of protein",
    "genes overlapping segmental duplications are masked in the copy-number matrices, so an "
    "absent copy-number value is not a copy-number-neutral gene",
)


class DepMapUnavailableError(Exception):
    """Raised when a DepMap file cannot be streamed and no cache covers the request."""


def symbol(column: str) -> str:
    """`"PRMT5 (10419)"` -> `"PRMT5"`; DepMap matrices name columns that way."""
    return column.split(" (")[0].strip().strip('"').upper()


def _path(dataset: str, cache_dir: Path) -> Path:
    return cache_dir / f"{dataset}.json.gz"


def _read(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        with gzip.open(path, "rt") as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError):
        return None


def _write(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt") as fh:
        json.dump(payload, fh, separators=(",", ":"))


def _stream(file_id: str, timeout: int = 180):
    url = FIGSHARE_FILE.format(file_id)
    req = urllib.request.Request(url, headers=UA)
    return urllib.request.urlopen(req, timeout=timeout)  # noqa: S310


def stream_wide(
    dataset: str,
    genes: set[str],
    decimals: int = 4,
    log=None,
) -> dict[str, Any]:
    """Stream a models-by-genes matrix and keep only the requested gene columns.

    The header gives the column order once; after that every row is split, the
    handful of wanted indices are read, and the row is thrown away. Nothing but
    the kept numbers is ever held in memory or written to disk.
    """
    file_id, name, description = FILES[dataset]
    wanted = {g.upper() for g in genes}
    values: dict[str, dict[str, float]] = {}
    index: dict[int, str] = {}
    duplicates: list[str] = []
    try:
        with _stream(file_id) as response:
            text = io.TextIOWrapper(response, encoding="utf-8", newline="")
            header = text.readline().rstrip("\r\n").split(",")
            seen: set[str] = set()
            for i, column in enumerate(header[1:], start=1):
                gene = symbol(column)
                if gene not in wanted:
                    continue
                if gene in seen:
                    duplicates.append(gene)
                    continue
                seen.add(gene)
                index[i] = gene
            columns = sorted(index)
            for line in text:
                if not line.strip():
                    continue
                fields = line.rstrip("\r\n").split(",")
                model = fields[0].strip('"')
                row: dict[str, float] = {}
                for i in columns:
                    raw = fields[i] if i < len(fields) else ""
                    if raw == "" or raw == "NA":
                        continue
                    try:
                        row[index[i]] = round(float(raw), decimals)
                    except ValueError:
                        continue
                if row:
                    values[model] = row
                if log and len(values) % 250 == 0 and row:
                    print(f"    {name}: {len(values)} models", file=log, flush=True)
    except OSError as e:
        raise DepMapUnavailableError(f"{name} ({file_id}): {type(e).__name__}: {str(e)[:160]}") from e
    found = sorted(set(index.values()))
    return {
        "dataset": dataset,
        "file": name,
        "figshare_file": file_id,
        "what": description,
        "release": RELEASE,
        "doi": DOI,
        "licence": LICENCE,
        "requested": sorted(wanted),
        "genes": found,
        "missing": sorted(wanted - set(found)),
        "duplicate_columns": sorted(set(duplicates)),
        "models": len(values),
        "values": values,
    }


def stream_narrow(dataset: str, key: str, columns: tuple[str, ...]) -> dict[str, Any]:
    """Stream a small one-row-per-model table (Model.csv, OmicsSignatures.csv).

    `OmicsSignatures.csv` leaves its index column unnamed, so a table asked for
    by `ModelID` silently produced nothing and the MSI-high control came back
    "untestable" rather than wrong. The key now falls back to the first column
    and the resolved name is recorded, and an empty table raises instead of
    being cached.
    """
    file_id, name, description = FILES[dataset]
    rows: dict[str, dict[str, Any]] = {}
    resolved = key
    try:
        with _stream(file_id) as response:
            reader = csv.DictReader(io.TextIOWrapper(response, encoding="utf-8", newline=""))
            fields = reader.fieldnames or []
            if key not in fields and fields:
                resolved = fields[0]
            for row in reader:
                ident = (row.get(resolved) or "").strip()
                if not ident:
                    continue
                rows[ident] = {c: (row.get(c) or "").strip() for c in columns}
    except OSError as e:
        raise DepMapUnavailableError(f"{name} ({file_id}): {type(e).__name__}: {str(e)[:160]}") from e
    if not rows:
        raise DepMapUnavailableError(
            f"{name} ({file_id}) returned no rows keyed by {resolved!r}; the table was not cached"
        )
    return {
        "dataset": dataset,
        "file": name,
        "figshare_file": file_id,
        "what": description,
        "release": RELEASE,
        "doi": DOI,
        "licence": LICENCE,
        "key": resolved,
        "key_requested": key,
        "columns": list(columns),
        "models": len(rows),
        "values": rows,
    }


def matrix(
    dataset: str,
    genes: set[str],
    net: bool = True,
    cache_dir: Path = CACHE,
    decimals: int = 4,
    log=None,
) -> dict[str, Any]:
    """The requested genes of a DepMap matrix, from the cache when it covers them.

    The cache records which genes were *asked for*, not only which were found,
    so a gene that DepMap does not carry is not re-fetched on every run.
    """
    path = _path(dataset, cache_dir)
    cached = _read(path)
    wanted = {g.upper() for g in genes}
    if cached is not None and cached.get("release") == RELEASE:
        if wanted <= set(cached.get("requested") or ()):
            cached["from_cache"] = True
            return cached
        wanted |= set(cached.get("requested") or ())
    if not net:
        raise DepMapUnavailableError(
            f"no cached {dataset} covering {len(wanted)} genes and network disabled "
            f"(cache: {path}{'' if cached else ', absent'})"
        )
    if log:
        print(f"  streaming {FILES[dataset][1]} for {len(wanted)} genes", file=log, flush=True)
    table = stream_wide(dataset, wanted, decimals=decimals, log=log)
    _write(path, table)
    table["from_cache"] = False
    return table


def table(
    dataset: str,
    key: str,
    columns: tuple[str, ...],
    net: bool = True,
    cache_dir: Path = CACHE,
    log=None,
) -> dict[str, Any]:
    """A small per-model table, from the cache when it has the columns."""
    path = _path(dataset, cache_dir)
    cached = _read(path)
    if (
        cached is not None
        and cached.get("release") == RELEASE
        and set(columns) <= set(cached["columns"])
        and cached.get("values")
    ):
        cached["from_cache"] = True
        return cached
    if not net:
        raise DepMapUnavailableError(f"no cached {dataset} and network disabled (cache: {path})")
    if log:
        print(f"  streaming {FILES[dataset][1]}", file=log, flush=True)
    got = stream_narrow(dataset, key, columns)
    _write(path, got)
    got["from_cache"] = False
    return got


def models(net: bool = True, cache_dir: Path = CACHE, log=None) -> dict[str, dict[str, Any]]:
    """Cell line -> lineage, Oncotree disease and subtype."""
    got = table(
        "model",
        "ModelID",
        ("OncotreeLineage", "OncotreePrimaryDisease", "OncotreeSubtype", "StrippedCellLineName"),
        net=net,
        cache_dir=cache_dir,
        log=log,
    )
    return got["values"]


def signatures(net: bool = True, cache_dir: Path = CACHE, log=None) -> dict[str, dict[str, Any]]:
    """Cell line -> MSI score and the other model-level signatures."""
    got = table(
        "signatures",
        "ModelID",
        ("MSIScore", "Ploidy", "LoHFraction", "WGD", "CIN", "Aneuploidy"),
        net=net,
        cache_dir=cache_dir,
        log=log,
    )
    return got["values"]


def provenance() -> dict[str, Any]:
    """What every result built on this module has to carry with it."""
    return {
        "release": RELEASE,
        "doi": DOI,
        "figshare_article": ARTICLE,
        "licence": LICENCE,
        "citation": CITATION,
        "files": {k: {"figshare_file": v[0], "name": v[1], "what": v[2]} for k, v in FILES.items()},
        "storage": "streamed over HTTP and discarded; only the requested gene columns are cached",
        "limits": list(LIMITS),
    }
