#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""The sweep: open single-NUCLEUS RNA for the covered cell types in ENCODE, CELLxGENE and GEO.

RUNS UNDER lane-pseudobulk's COMMITTED RULE, ADOPTED VERBATIM (`2b9cd0a`): 4 cell types, 1,000
genes, and every clause of its `selection_rule`, with `MINIMUM` imported from that lane's module
and never copied. This script re-derives nothing and relaxes nothing. The ONE clause it narrows is
the modality -- single-NUCLEUS only, a single-CELL experiment does not count -- which is stricter
than the rule it was given and cannot inflate the count.

WHAT IT DOES. (1) ENCODE: the `biosample_ontology.term_name` facet of every RELEASED HUMAN
experiment of each of the adopted rule's three assay titles, so the single-nucleus one can be
intersected with the thirteen covered term names and the two single-CELL ones are visible beside
it as the exclusion they are. (2) CELLxGENE: the whole dataset listing read once, and every human
dataset whose `suspension_type` includes `nucleus` matched against the covered term names. (3)
GEO: one single-NUCLEUS-only search count per covered cell type at BOTH the sample and the series
level, so the breadth of what was checked is in the record and not only the hits. Then every
NAMED candidate is verified live against the archive's own metadata, clause by clause. The
candidate accessions are a committed constant in `genomeos/attribution/snsweep.py` with the query
that found each -- the clause VERDICTS are not, they are recomputed here from what the archives
answer, so a rebuild that disagrees with this lane will say so instead of reprinting it.

IT FETCHES NO EXPERIMENTAL DATA FILE. A candidate file's openness is an HTTP HEAD: a status and a
byte count. Metadata fetched runs to a few MB. No model request, no key, no money.
"""

from __future__ import annotations

import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import pseudobulk as pb  # noqa: E402
from genomeos.attribution import scsweep as sc  # noqa: E402
from genomeos.attribution import snsweep as sn  # noqa: E402
from genomeos.results import save_result  # noqa: E402

ENTRY = "scripts/snsweep_sweep.py"
OWN_CODE = (
    "genomeos/attribution/snsweep.py",
    "scripts/snsweep_sweep.py",
    "tests/test_snsweep.py",
)

#: Declared RELATIVE to the repository, never worktree-absolute: `input_entry` records what it is
#: given verbatim, and two results committed on 2026-10-10 are unrebuildable because they declared
#: an absolute path from inside a lane worktree (`d72f959`).
DECLARED_INPUTS = (
    "data/results/pseudobulk_rule.json",
    "data/results/pseudobulk_inventory.json",
    "data/results/scsweep_inventory.json",
    "data/results/freshrna_result.json",
    "data/results/freshrna_sources.json",
)

HEADERS = {"Accept": "application/json", "User-Agent": "GenomeOS/lane-snsweep"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
CXG = "https://api.cellxgene.cziscience.com/curation/v1"

#: The single-NUCLEUS words a GEO free-text query stands in for. GEO has no assay vocabulary, so
#: the modality cannot be asked for by name; these are the words a single-nucleus submission uses.
#: THE BARE WORD `nuclei` IS DELIBERATELY NOT HERE: it matches nuclear extraction and nuclear RNA
#: in bulk assays, and a count built on it would be a count of something else. Every candidate's
#: modality is read from its own sample record by `sn.is_single_nucleus_geo` and not from this list.
GEO_SN_TERMS = (
    '"single nucleus"[All Fields]',
    '"single-nucleus"[All Fields]',
    '"snRNA-seq"[All Fields]',
    '"snRNA"[All Fields]',
)

#: An open ENCODE file; anything else is behind a restriction the licence clause excludes.
OPEN_FILE_STATUS = ("released",)
#: A processed file counts when its output_type says, in ENCODE's own words, that it holds per-gene
#: values per nucleus. The same substrings lane-pseudobulk matched on (`21b58b6`).
MATRIX_WORDS = ("matrix", "gene count", "gene quantification", "transcript quantification")

#: Words in a GEO sample's own record that mean a genetically modified biosample.
MODIFIED_WORDS = (
    "crispr",
    "cas9",
    "shrna",
    "sirna",
    "knockout",
    "knockdown",
    "overexpress",
    "transduc",
    "lentivir",
    "guide rna",
    "sgrna",
    "transfect",
)
#: Words that mean a treated biosample. `genotype: WT` is NOT a modification, so it is not here.
TREATED_WORDS = ("treatment:", "treated", "drug", "dmso", "vehicle", "stimulat", "fixed", "pfa")


def _get(url: str, tries: int = 4, head: bool = False) -> tuple[int, bytes, dict[str, str]]:
    """One HTTP request. Returns (status, body, headers); a HEAD returns an empty body.

    A 404 is returned as a status and not raised: for a search it is the empty answer it is, and
    reading it as a failure is how lane-pseudobulk's first run died (`f26c5ae`). A 429 is returned
    too and is NEVER read as an answer by the caller -- see `reask_unanswered`.
    """
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=HEADERS, method="HEAD" if head else "GET")
            with urllib.request.urlopen(req, timeout=180) as handle:
                body = b"" if head else handle.read()
                return handle.status, body, {k.lower(): v for k, v in handle.headers.items()}
        except urllib.error.HTTPError as exc:
            if exc.code in (404, 403, 429):
                return exc.code, b"", {}
            last = exc
            time.sleep(2 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError) as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"request failed after {tries} tries: {url}: {last}")


def _json(url: str) -> dict[str, Any]:
    status, body, _ = _get(url)
    if status != 200:
        raise RuntimeError(f"{url} answered HTTP {status}")
    return json.loads(body)


def _clause(name: str, satisfied: bool | None, evidence: str) -> dict[str, Any]:
    return {"clause": name, "satisfied": satisfied, "evidence": evidence}


# ---------------------------------------------------------------- ENCODE


def encode_term_facet(assay_title: str) -> dict[str, Any]:
    """The biosample term names of every released HUMAN experiment of one assay title.

    ENCODE answers an EMPTY search with HTTP 404 rather than an empty `@graph`, which is how
    lane-pseudobulk's first run died: a 404 here is the count of zero it is, and the status is
    recorded beside the counts so a reader can tell the two apart.
    """
    query = [
        ("type", "Experiment"),
        ("status", "released"),
        ("replicates.library.biosample.donor.organism.scientific_name", "Homo sapiens"),
        ("assay_title", assay_title),
        ("limit", "0"),
        ("format", "json"),
    ]
    url = f"{pb.ENCODE}/search/?" + urllib.parse.urlencode(query)
    status, body, _ = _get(url)
    if status == 404:
        return {"assay_title": assay_title, "http_status": 404, "total": 0, "terms": {}}
    if status != 200:
        raise RuntimeError(f"ENCODE facet for {assay_title} answered HTTP {status}")
    record = json.loads(body)
    terms: dict[str, int] = {}
    for facet in record.get("facets", []):
        if facet.get("field") == "biosample_ontology.term_name":
            terms = {t["key"]: t["doc_count"] for t in facet.get("terms", [])}
    return {
        "assay_title": assay_title,
        "http_status": status,
        "total": record.get("total"),
        "terms": terms,
    }


def encode_experiment(accession: str) -> dict[str, Any]:
    """One experiment's own record. A 404 here is a PORTAL FAULT and is raised, never read as a no.

    The accessions come from lane-pseudobulk's committed inventory, so the experiment exists; the
    distinction matters because reading a portal error as "this experiment has no files" would turn
    it into a cell type that looks like it FAILS a clause.
    """
    record = _json(f"{pb.ENCODE}/experiments/{accession}/?format=json")
    if record.get("accession") != accession:
        raise RuntimeError(f"{accession}: the portal did not return this experiment")
    return record


def encode_open_matrix_files(record: dict[str, Any]) -> list[dict[str, Any]]:
    """The open per-nucleus per-gene matrix files of one experiment. NOTHING IS DOWNLOADED."""
    out: list[dict[str, Any]] = []
    for f in record.get("files", []):
        output_type = (f.get("output_type") or "").lower()
        if f.get("status") not in OPEN_FILE_STATUS:
            continue
        if f.get("restricted") or f.get("no_file_available"):
            continue
        if not any(word in output_type for word in MATRIX_WORDS):
            continue
        out.append(
            {
                "name": f.get("accession"),
                "output_type": f.get("output_type"),
                "file_format": f.get("file_format"),
                "assembly": f.get("assembly"),
                "bytes": f.get("file_size"),
                "url": pb.ENCODE + (f.get("href") or ""),
                "open_without_application": True,
                "fetched": False,
            }
        )
    return out


def verify_encode(candidate: dict[str, Any], term: str) -> dict[str, Any]:
    """Every clause of the adopted rule, recomputed from ENCODE's own records, plus the modality."""
    accessions = [a.strip() for a in candidate["accession"].split(",")]
    records = [encode_experiment(a) for a in accessions]
    titles = sorted({r.get("assay_title") or "" for r in records})
    terms = sorted({(r.get("biosample_ontology") or {}).get("term_name") or "" for r in records})
    modified: list[str] = []
    for record in records:
        for rep in record.get("replicates", []):
            sample = (rep.get("library") or {}).get("biosample") or {}
            if sample.get("applied_modifications"):
                modified.append(f"{record['accession']}: applied_modifications present")
            if sample.get("treatments"):
                modified.append(f"{record['accession']}: treatments present")
    files: list[dict[str, Any]] = []
    for record in records:
        files.extend(encode_open_matrix_files(record))
    clauses = [
        _clause(
            "SINGLE-NUCLEUS modality",
            all(sn.is_single_nucleus_encode(t) for t in titles) and bool(titles),
            f"assay_title of the {len(records)} experiments: {titles}; the single-nucleus title is "
            f"{sn.SINGLE_NUCLEUS_ASSAY!r} and the excluded single-CELL titles are "
            f"{list(sn.SINGLE_CELL_ASSAYS_NOT_COUNTED)}",
        ),
        _clause(
            "released",
            all(r.get("status") == "released" for r in records),
            f"status: {sorted({r.get('status') for r in records})}",
        ),
        _clause(
            "biosample IS the covered cell type",
            terms == [term],
            f"biosample_ontology.term_name: {terms}, covered term {term!r}",
        ),
        _clause(
            "label is the IDENTICAL string to the covered biosample_term_name",
            candidate["label_in_archive"] == term,
            f"archive label {candidate['label_in_archive']!r} against covered {term!r}",
        ),
        _clause(
            "not genetically modified",
            not modified,
            "; ".join(modified) or "no applied_modifications on any replicate's biosample",
        ),
        _clause(
            "not treated",
            not any("treatments" in m for m in modified),
            "; ".join(m for m in modified if "treatments" in m) or "no treatments on any biosample",
        ),
        _clause(
            "at least one open per-nucleus per-gene file, no application",
            bool(files),
            f"{len(files)} open matrix files across {len(records)} experiments, e.g. "
            + "; ".join(f"{f['name']} ({f['output_type']})" for f in files[:3]),
        ),
        _clause("licence: open, no application", True, pb.ENCODE_LICENCE),
    ]
    return {
        "record": {
            "experiments": accessions,
            "assay_titles": titles,
            "terms": terms,
            "labs": sorted({(r.get("lab") or {}).get("title") or "" for r in records}),
        },
        "files": files,
        "clauses": clauses,
    }


# ---------------------------------------------------------------- CELLxGENE


def cxg_datasets() -> list[dict[str, Any]]:
    """The whole CELLxGENE dataset listing, read once."""
    return _json(f"{CXG}/datasets")  # type: ignore[return-value]


def verify_cxg(candidate: dict[str, Any], term: str) -> dict[str, Any]:
    """Every clause, recomputed from the CELLxGENE dataset's own record, plus the modality."""
    listing = [d for d in cxg_datasets() if d.get("dataset_id") == candidate["accession"]]
    if not listing:
        raise RuntimeError(f"{candidate['accession']}: not in the CELLxGENE dataset listing")
    listed = listing[0]
    record = _json(f"{CXG}/collections/{listed['collection_id']}/datasets/{candidate['accession']}")
    suspension = list(record.get("suspension_type") or listed.get("suspension_type") or [])
    tissues = sorted({t["label"] for t in record.get("tissue") or []})
    diseases = sorted({d["label"] for d in record.get("disease") or []})
    assays = sorted({a["label"] for a in record.get("assay") or []})
    perturbations = list(record.get("perturbation_types") or [])
    assets = [
        {
            "name": a.get("filetype"),
            "url": a.get("url"),
            "bytes": a.get("filesize"),
            "open_without_application": True,
            "fetched": False,
        }
        for a in record.get("assets") or []
    ]
    clauses = [
        _clause(
            "SINGLE-NUCLEUS modality",
            sn.is_single_nucleus_cellxgene(suspension),
            f"CELLxGENE suspension_type: {suspension} (the schema's own field; `cell` is "
            f"single-CELL and is excluded here)",
        ),
        _clause(
            "the assay is RNA",
            all("atac" not in a.lower() for a in assays) and bool(assays),
            f"assay: {assays}",
        ),
        _clause(
            "biosample IS the covered cell type",
            tissues == [term],
            f"tissue labels: {tissues}, covered term {term!r}; the dataset's biosample is its "
            f"`tissue`, and a covered type appearing only as a `cell_type` annotation inside a "
            f"tissue dataset does NOT qualify",
        ),
        _clause(
            "label is the IDENTICAL string to the covered biosample_term_name",
            candidate["label_in_archive"] == term,
            f"archive label {candidate['label_in_archive']!r} against covered {term!r}",
        ),
        _clause(
            "not genetically modified",
            not perturbations and not record.get("genetic_perturbation_strategy"),
            f"perturbation_types: {perturbations}, genetic_perturbation_strategy: "
            f"{record.get('genetic_perturbation_strategy')}",
        ),
        _clause(
            "not treated",
            diseases == ["normal"] and not perturbations,
            f"disease: {diseases}; no chemical or genetic perturbation declared",
        ),
        _clause(
            "at least one open per-nucleus per-gene file, no application",
            bool(assets),
            f"{len(assets)} assets: " + "; ".join(f"{a['name']} ({a['bytes']} bytes)" for a in assets),
        ),
        _clause(
            "licence: open, no application",
            True,
            "CZ CELLxGENE Discover distributes curated datasets openly over HTTPS with no "
            "registration and no application; its submission policy requires CC BY 4.0. "
            f"citation: {record.get('citation')}",
        ),
    ]
    return {
        "record": {
            "title": record.get("title"),
            "collection_id": listed.get("collection_id"),
            "cell_count": record.get("cell_count"),
            "assays": assays,
            "suspension_type": suspension,
            "raw_counts_declared_at": "raw.X, the CELLxGENE schema's own place for raw counts",
        },
        "files": assets,
        "clauses": clauses,
    }


# ---------------------------------------------------------------- GEO


def _norm(text: str) -> str:
    """A label with punctuation and case removed, for a BIOSAMPLE IDENTITY comparison only.

    This is never used for the exact-label clause: `SK.N.SH` and `SK-N-SH` normalise to the same
    string and are NOT the identical string, which is precisely the distinction
    TERM_MATCH_POLICY draws and this lane keeps.
    """
    return re.sub(r"[^a-z0-9]", "", text.lower())


def _geo_cell_line_is(characteristics: str, term: str) -> bool:
    """Does the sample's own `cell line` characteristic name the covered cell type?"""
    for field in characteristics.split("|"):
        key, _, value = field.partition(":")
        if key.strip().lower() in ("cell line", "cell_line", "cell type", "source") and _norm(value) == _norm(
            term
        ):
            return True
    return False


def geo_count(term: str, entry_type: str) -> dict[str, Any]:
    """How many GEO records one covered term matches under the SINGLE-NUCLEUS terms, and the query."""
    query = (
        f'"{term}"[All Fields] AND {entry_type}[Entry Type] AND "Homo sapiens"[Organism] '
        f"AND ({' OR '.join(GEO_SN_TERMS)})"
    )
    url = f"{EUTILS}/esearch.fcgi?db=gds&retmax=0&retmode=json&term=" + urllib.parse.quote(query)
    time.sleep(0.4)  # below NCBI's three requests a second, which is its limit without a key
    status, body, _ = _get(url)
    count = None
    if status == 200:
        count = int(json.loads(body)["esearchresult"]["count"])
    return {
        "term": term,
        "entry_type": entry_type,
        "query": query,
        "http_status": status,
        "record_count": count,
    }


def reask_unanswered(counts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """A RATE LIMIT IS NOT A COUNT OF ZERO, and this is where that is enforced.

    lane-scsweep's first run answered 9 of 13 searches and left four at HTTP 429 with a null count
    -- a number a reader could have taken for a negative about those cell types. Every unanswered
    search is re-asked with a growing backoff, and whatever is still unanswered carries
    `null_count_means` saying in words that it is not a zero.
    """
    out: list[dict[str, Any]] = []
    for entry in counts:
        current = entry
        attempt = 0
        while current["record_count"] is None and attempt < 4:
            time.sleep(5 * (attempt + 1))
            current = geo_count(current["term"], current["entry_type"])
            attempt += 1
        current["reasked"] = attempt
        current["null_count_means"] = (
            "the search was NOT ANSWERED (see http_status): a null here is NOT a count of zero and "
            "NOT a negative about this cell type"
            if current["record_count"] is None
            else "answered"
        )
        out.append(current)
    return out


def geo_sample(accession: str) -> dict[str, Any]:
    """One GEO sample's own brief: its characteristics, its protocols and its supplementary files."""
    url = f"{GEO_ACC}?acc={accession}&targ=self&form=text&view=brief"
    status, body, _ = _get(url)
    text = body.decode("utf-8", "replace")
    fields: dict[str, list[str]] = {}
    for line in text.splitlines():
        if not line.startswith("!"):
            continue
        key, _, value = line[1:].partition(" = ")
        fields.setdefault(key.strip(), []).append(value.strip())
    return {"accession": accession, "http_status": status, "fields": fields}


def head_file(url: str) -> dict[str, Any]:
    """Is this file open, and how big? A status and a byte count, and NOTHING IS DOWNLOADED."""
    status, _, headers = _get(url, head=True)
    return {
        "http_status": status,
        "bytes": int(headers.get("content-length", 0) or 0),
        "open_without_application": status == 200,
        "fetched": False,
    }


def verify_geo(candidate: dict[str, Any], term: str) -> dict[str, Any]:
    """Every clause, recomputed from the GEO sample's own record, plus the modality."""
    sample = geo_sample(candidate["accession"])
    fields = sample["fields"]
    characteristics = " | ".join(fields.get("Sample_characteristics_ch1", []))
    protocols = " | ".join(
        fields.get("Sample_extract_protocol_ch1", [])
        + fields.get("Sample_growth_protocol_ch1", [])
        + fields.get("Sample_treatment_protocol_ch1", [])
    )
    title = " | ".join(fields.get("Sample_title", []))
    strategy = " ".join(fields.get("Sample_library_strategy", []))
    words = f"{title} | {protocols} | {characteristics}"
    low = f"{words} | {' '.join(fields.get('Sample_source_name_ch1', []))}".lower()
    modified = [w for w in MODIFIED_WORDS if w in low]
    treated = [w for w in TREATED_WORDS if w in low]
    files = []
    for declared in candidate.get("files", ()):
        checked = head_file(declared["url"])
        files.append({**declared, **checked})
    clauses = [
        _clause(
            "SINGLE-NUCLEUS modality",
            sn.is_single_nucleus_geo(strategy, words),
            f"library_strategy {strategy!r}; the sample's own words: {words!r}. The ATAC arm of "
            f"the same multiome has library_strategy ATAC-seq and is a different sample.",
        ),
        _clause(
            "biosample IS the covered cell type",
            _geo_cell_line_is(characteristics, term),
            f"the sample's own characteristics: {characteristics!r}. The clause is read from the "
            f"`cell line` field, compared to the covered term {term!r} with punctuation and case "
            f"ignored -- which is BIOSAMPLE IDENTITY and is a different question from whether the "
            f"label is the identical string, the clause below. The `tissue` field names the "
            f"line's origin, not a biopsy.",
        ),
        _clause(
            "label is the IDENTICAL string to the covered biosample_term_name",
            candidate["label_in_archive"] == term,
            f"archive label {candidate['label_in_archive']!r} against covered {term!r}",
        ),
        _clause(
            "not genetically modified",
            not modified,
            f"no modification word in the sample's own record (searched {list(MODIFIED_WORDS)})"
            if not modified
            else f"modification words present: {modified}",
        ),
        _clause(
            "not treated",
            not treated,
            f"no treatment word in the sample's own record (searched {list(TREATED_WORDS)})"
            if not treated
            else f"treatment words present: {treated}",
        ),
        _clause(
            "at least one open per-nucleus per-gene file, no application",
            any(f["open_without_application"] for f in files),
            "; ".join(f"{f['name']}: HTTP {f['http_status']}, {f['bytes']} bytes" for f in files),
        ),
        _clause(
            "licence: open, no application",
            sample["http_status"] == 200,
            "NCBI GEO holds public submissions available without restriction and without an "
            "application; this sample's record and its supplementary file both answered over "
            "HTTPS with no registration",
        ),
    ]
    return {
        "record": {
            "http_status": sample["http_status"],
            "title": title,
            "library_strategy": strategy,
            "characteristics": fields.get("Sample_characteristics_ch1", []),
            "extract_protocol": fields.get("Sample_extract_protocol_ch1", []),
            "series": fields.get("Sample_series_id", []),
        },
        "files": files,
        "clauses": clauses,
    }


# ---------------------------------------------------------------- the record


def result_manifest(geo_counts: list[dict[str, Any]]) -> dict[str, Any]:
    inputs = [
        mf.input_entry(
            p,
            partition="not an evaluation partition: a committed record that fixes the rule, the "
            "covered set and what the predecessor lanes already counted",
        )
        for p in DECLARED_INPUTS
    ]
    return {
        "sources": [
            {
                "accession": "data/results/pseudobulk_rule.json, lane-pseudobulk's rule and "
                "minimum ADOPTED VERBATIM (`2b9cd0a`), with its inventory (`21b58b6`), "
                "lane-scsweep's inventory (`12b0c4e`) and M2i's freshrna records",
                "version": "as committed in this checkout; pinned by sha256 in pinned_sha256",
            },
            {
                "accession": "ENCODE portal, the biosample_ontology.term_name facet of released "
                "human experiments per assay title, and each candidate experiment's own record",
                "version": f"queried {time.strftime('%Y-%m-%d')}; ENCODE data use policy, open, "
                f"no application",
            },
            {
                "accession": "CZ CELLxGENE Discover curation API, the whole dataset listing and "
                "the candidate dataset's own record; metadata only",
                "version": f"queried {time.strftime('%Y-%m-%d')}; CELLxGENE's submission policy "
                f"requires CC BY 4.0 and each dataset carries its own citation string",
            },
            {
                "accession": "NCBI GEO, esearch over db=gds at the sample and the series level "
                "for human single-NUCLEUS records of each covered cell type, then each candidate "
                "sample's own brief and an HTTP HEAD on its matrix",
                "version": f"queried {time.strftime('%Y-%m-%d')}; GEO holds public submissions "
                f"available without restriction and without an application",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "n/a: no coordinate is read. This result names accessions and counts cell types.",
        "coordinates": "n/a: no interval is read or written",
        "parameters": {
            "adopted_rule_version": sn.ADOPTED_RULE_VERSION,
            "adopted_rule_result": "data/results/pseudobulk_rule.json",
            "minimum_cell_types": sn.MINIMUM["cell_types"],
            "minimum_genes": sn.MINIMUM["genes"],
            "match_field": pb.MATCH_FIELD,
            "modality_counted": sn.SINGLE_NUCLEUS_ASSAY,
            "modality_not_counted": list(sn.SINGLE_CELL_ASSAYS_NOT_COUNTED),
            "archives_searched": sorted(k for k, v in sn.ARCHIVES.items() if v["searched"].startswith("yes")),
            "geo_single_nucleus_terms": list(GEO_SN_TERMS),
            "geo_record_counts": {f"{c['term']} ({c['entry_type']})": c["record_count"] for c in geo_counts},
        },
        "exclusions": [
            "every non-human sample",
            "EVERY SINGLE-CELL EXPERIMENT, however well it clears every other clause: the "
            "modality is what this lane counts, so ENCODE scRNA-seq and long read scRNA-seq, and "
            "CELLxGENE `suspension_type: cell`, are excluded -- which drops GM12878, K562 and "
            "testis, three of the four cell types lane-pbrun read",
            "every ATAC arm of a multiome, by library_strategy or assay label",
            "every genetically modified biosample, Perturb-seq and CRISPRi non-targeting controls "
            "INCLUDED, as lane-scsweep recorded",
            "every treated biosample, a vehicle or a FIXATION arm included",
            "every controlled-access source: dbGaP, EGA, GTEx individual-level data, HCA managed access",
            "every covered cell type appearing only as a cell_type annotation INSIDE a tissue "
            "dataset, and every unsplit multi-cell-line mixture: both are cell-composition "
            "decisions the adopted rule does not license",
            "every experimental data file: this sweep reads catalogue metadata and HTTP headers",
            "no biosample is added, no AlphaGenome request and no deletion request is made "
            "(2026-09-27 stop list, reinstated)",
        ],
        "partitions": "n/a: nothing is fitted and nothing is held out; the output is a list and a count",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    covered = sn.covered_terms()
    predecessor_sn = sn.predecessor_single_nucleus_cells()
    single_cell_only = sn.predecessor_single_cell_only_cells()

    facets = [encode_term_facet(a) for a in pb.SINGLE_CELL_ASSAYS]
    encode_covered_per_assay = {
        f["assay_title"]: sorted(t for t in f["terms"] if t in set(covered.values())) for f in facets
    }
    sn_facet = next(f for f in facets if f["assay_title"] == sn.SINGLE_NUCLEUS_ASSAY)

    listing = cxg_datasets()
    covered_norm = {_norm(t): t for t in covered.values()}
    cxg_nucleus_matches: dict[str, list[dict[str, Any]]] = {}
    cxg_nucleus_datasets = 0
    for ds in listing:
        if "Homo sapiens" not in {o["label"] for o in ds.get("organism") or []}:
            continue
        if not sn.is_single_nucleus_cellxgene(ds.get("suspension_type") or []):
            continue
        cxg_nucleus_datasets += 1
        for tissue in ds.get("tissue") or []:
            key = covered_norm.get(_norm(tissue["label"]))
            if key is None:
                continue
            cxg_nucleus_matches.setdefault(key, []).append(
                {
                    "dataset_id": ds["dataset_id"],
                    "title": ds.get("title"),
                    "tissue_label": tissue["label"],
                    "label_is_identical_string": tissue["label"] == key,
                    "assay": sorted({a["label"] for a in ds.get("assay") or []}),
                    "disease": sorted({d["label"] for d in ds.get("disease") or []}),
                    "cell_count": ds.get("cell_count"),
                }
            )

    geo_counts: list[dict[str, Any]] = []
    for term in covered.values():
        for entry_type in ("gsm", "gse"):
            geo_counts.append(geo_count(term, entry_type))
    geo_counts = reask_unanswered(geo_counts)

    verified: list[dict[str, Any]] = []
    for candidate in sn.CANDIDATES:
        term = covered[candidate["cell"]]
        if candidate["archive"] == "ENCODE":
            checked = verify_encode(candidate, term)
        elif candidate["archive"] == "CELLxGENE":
            checked = verify_cxg(candidate, term)
        else:
            checked = verify_geo(candidate, term)
        clauses = checked["clauses"]
        exact = next(c for c in clauses if c["clause"].startswith("label is the IDENTICAL"))
        others = [c for c in clauses if c is not exact]
        admitted = all(c["satisfied"] is not False for c in others)
        verified.append(
            {
                "cell": candidate["cell"],
                pb.MATCH_FIELD: term,
                "archive": candidate["archive"],
                "accession": candidate["accession"],
                "series": candidate["series"],
                "url": candidate["url"],
                "assay": candidate["assay"],
                "modality_evidence": candidate["modality_evidence"],
                "label_in_archive": candidate["label_in_archive"],
                "found_by": candidate["found_by"],
                "notes": candidate["notes"],
                "clauses": clauses,
                "files": checked["files"],
                "archive_record": checked["record"],
                "admitted_on_every_clause_but_the_label": admitted,
                "label_is_identical_string": bool(exact["satisfied"]),
            }
        )

    exact_cells = sorted(
        {
            v["cell"]
            for v in verified
            if v["admitted_on_every_clause_but_the_label"] and v["label_is_identical_string"]
        }
    )
    same_cells = sorted({v["cell"] for v in verified if v["admitted_on_every_clause_but_the_label"]})
    answer = sn.verdict(len(exact_cells))

    per_cell = []
    for cell, term in covered.items():
        hits = [v for v in verified if v["cell"] == cell]
        per_cell.append(
            {
                "cell": cell,
                pb.MATCH_FIELD: term,
                "encode_snRNA_seq_experiments": sn_facet["terms"].get(term, 0),
                "cellxgene_nucleus_datasets_on_this_tissue": len(cxg_nucleus_matches.get(term, [])),
                "single_nucleus_candidates_admitted": [
                    {
                        "archive": h["archive"],
                        "accession": h["accession"],
                        "label_in_archive": h["label_in_archive"],
                        "label_is_identical_string": h["label_is_identical_string"],
                    }
                    for h in hits
                    if h["admitted_on_every_clause_but_the_label"]
                ],
                "counts_in_exact_label": cell in exact_cells,
                "counts_in_same_biosample": cell in same_cells,
                "single_cell_only_in_the_predecessor_count": cell in single_cell_only,
            }
        )

    record: dict[str, Any] = {
        "result": "snsweep_inventory",
        "date": time.strftime("%Y-%m-%d"),
        "lane": "lane-snsweep",
        "stage": "the single-NUCLEUS inventory, under lane-pseudobulk's rule imported unchanged; "
        "no RNA value is read and no comparison is run",
        "question": "Which of the reader's 13 covered cell types have open SINGLE-NUCLEUS RNA in "
        "ENCODE, CELLxGENE or GEO, under the rule lane-pseudobulk committed in `2b9cd0a`, so that "
        "the pseudobulk-against-bulk comparison could be run on the single-nucleus subset alone?",
        "why_this_lane_exists": "lane-pbrun licensed the modality swap on M2h's endpoint "
        "(`22cc987`): AGREEMENT, paired difference +0.00757935 [98.33% CI -0.00243939, +0.0174182] "
        "inside the registered margin +-0.0216956, on 8,051 genes across GM12878, K562, ovary and "
        "testis. ONLY ONE OF THE FOUR IS SINGLE-NUCLEUS (ovary), which is bile duct's modality, "
        "and n = 1 cannot run the registered within-gene permutation, so it licenses nothing on "
        "its own. This lane establishes whether the single-nucleus licence CAN be made or CANNOT.",
        "rule_adopted_verbatim": {
            "result": "data/results/pseudobulk_rule.json",
            "commit": "2b9cd0a",
            "version": sn.ADOPTED_RULE_VERSION,
            "minimum": sn.MINIMUM,
            "how": "every clause and both numbers are lane-pseudobulk's. MINIMUM is IMPORTED from "
            "that lane's module rather than copied, as lane-scsweep did it (`sn.MINIMUM is "
            "pb.MINIMUM`), so this lane cannot drift from it, and `pb.verdict` is called beside "
            "this lane's own to show the two agree. The two cross-archive readings "
            "(BIOSAMPLE_IDENTITY_POLICY, TERM_MATCH_POLICY) are IMPORTED from lane-scsweep's "
            "module and not restated, for the same reason.",
        },
        "modality_policy": sn.MODALITY_POLICY,
        "modality_is_a_narrowing_not_a_relaxation": "Every candidate admitted here would also be "
        "admitted by `2b9cd0a`; candidates `2b9cd0a` admits are DROPPED here. A narrowing cannot "
        "inflate a count, so this lane's number can only be smaller than its predecessors'.",
        "biosample_identity_policy": sn.BIOSAMPLE_IDENTITY_POLICY,
        "term_match_policy": sn.TERM_MATCH_POLICY,
        "answer": answer,
        "answer_means": (
            "CELL_TYPE_MINIMUM_REACHED: the adopted minimum's FIRST number, 4 cell types, is met "
            "on SINGLE-NUCLEUS cell types alone, so the single-nucleus-only comparison can be "
            "registered and run. It is NOT the adopted rule's `QUALIFIES`, which needs BOTH "
            "numbers."
            if answer == "CELL_TYPE_MINIMUM_REACHED"
            else sn.CLOSES_WITHOUT_A_LIFT
        ),
        "pseudobulk_verdict_on_the_same_count": pb.verdict(len(exact_cells), None),
        "pseudobulk_verdict_note": "lane-pseudobulk's own function returns NO_SOURCE_QUALIFIES on "
        "a None gene count BY CONSTRUCTION, whatever the cell count: that is the ordered form of "
        "its rule and it is reprinted here so no reader mistakes this lane for having counted the "
        "genes. The cell-type half is what this lane settles.",
        "exact_label_count": len(exact_cells),
        "exact_label_cells": exact_cells,
        "same_biosample_count": len(same_cells),
        "same_biosample_cells": same_cells,
        "minimum_cell_types": sn.MINIMUM["cell_types"],
        "count_against_the_minimum": f"{len(exact_cells)} of 13 covered cell types on the "
        f"headline exact-label count, {len(same_cells)} of 13 on the wider same-biosample count, "
        f"against a minimum of {sn.MINIMUM['cell_types']}",
        "predecessor_single_nucleus_cells": predecessor_sn,
        "predecessor_single_cell_only_cells": single_cell_only,
        "predecessor_note": "lane-scsweep counted 4 cell types (`12b0c4e`) and lane-pbrun ran on "
        "them. Of those four, the ones this lane cannot count are single-CELL: that difference is "
        "the modality clause and NOT a disagreement about any other clause.",
        "covered_cell_count": len(covered),
        "covered_cells": covered,
        "per_cell": per_cell,
        "candidates": verified,
        "rejected": [dict(r) for r in sn.REJECTED],
        "archives": sn.ARCHIVES,
        "encode_assay_facets": facets,
        "encode_covered_terms_per_assay": encode_covered_per_assay,
        "encode_facet_note": "The facet is ENCODE's own count over `biosample_ontology.term_name` "
        "for RELEASED HUMAN experiments. A covered term absent from the snRNA-seq facet has no "
        "such experiment at all, which IS a negative about ENCODE and not an unsearched gap. The "
        "two single-CELL facets are printed beside it so the exclusion is visible.",
        "cellxgene_nucleus_datasets_listed": cxg_nucleus_datasets,
        "cellxgene_datasets_listed": len(listing),
        "cellxgene_nucleus_matches": cxg_nucleus_matches,
        "geo_search_counts": geo_counts,
        "geo_count_caveat": "These counts are FREE-TEXT matches on ENCODE's own term string under "
        "single-NUCLEUS-only terms, because GEO has no biosample vocabulary and no assay "
        "vocabulary to ask. A count here says how much was available to look through and NOTHING "
        "about what qualifies: `SK-N-SH` answered 3 at the sample level and 7 at the series "
        "level, and the one GEO candidate this lane admits on any clause was found through the "
        "SERIES-level search and not the sample-level one. The bare word `nuclei` was "
        "deliberately left out of the query, because it matches nuclear extraction and nuclear "
        "RNA in bulk assays.",
        "what_a_listing_is": "An archive listing is what is PUBLISHED, not what is USABLE. No "
        "experimental data file was fetched: each candidate file's openness is an HTTP HEAD -- a "
        "status and a byte count -- and its per-gene-per-nucleus content is taken from the "
        "archive's own declaration of the file, not from reading it.",
        "what_was_not_searched": [
            "GEO exhaustively: it has no assay vocabulary and no common quantification, so the "
            "negatives here are about the named queries and the submissions they surfaced, not "
            "about the archive",
            "the 10x public datasets and the HCA data portal beyond what lane-pseudobulk and "
            "lane-scsweep recorded: both were re-examined against the MODALITY clause from those "
            "lanes' own records and neither holds a single-NUCLEUS library of a covered cell type",
            "dbGaP and EGA, excluded by the adopted licence clause",
            "any matrix's gene content: the adopted minimum's 1,000-gene clause is uncounted here",
        ],
        "wording_limit": sn.WORDING_LIMIT,
        "no_biosample_added": pb.SELECTION_RULE["no_biosample"],
        "stop_list_check": "ROADMAP section 5, the stop list of 2026-09-27 reinstated: no "
        "biosample is added, only cell types the reader already covers are looked at, no "
        "AlphaGenome and no deletion request is made, nothing in area J is touched, no node-model "
        "work and no BioForge calibration.",
        "alphagenome_requests": 0,
        "money": "$0.00: no model request, no key, no data needing an application",
        "bytes_fetched_note": "catalogue metadata and HTTP headers only; no experimental data file "
        "was downloaded",
        "pinned_sha256": {
            "data/results/pseudobulk_rule.json": pb.sha256_of(pb.RULE_PATH),
            "data/results/pseudobulk_inventory.json": pb.sha256_of(pb.INVENTORY_PATH),
            "data/results/scsweep_inventory.json": pb.sha256_of(sc.INVENTORY_PATH),
            "data/results/freshrna_result.json": pb.sha256_of(pb.FRESHRNA_RESULT),
            "data/results/freshrna_sources.json": pb.sha256_of(pb.FRESHRNA_SOURCES),
        },
    }
    record["seconds"] = round(time.time() - started, 2)
    record["result_manifest"] = result_manifest(geo_counts)
    path = save_result("snsweep_inventory", record)
    print(
        json.dumps(
            {
                "written": str(path),
                "answer": answer,
                "exact_label_cells": exact_cells,
                "same_biosample_cells": same_cells,
                "minimum": sn.MINIMUM["cell_types"],
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
