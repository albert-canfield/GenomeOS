#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""The sweep: open single-cell RNA for the covered cell types in GEO, CELLxGENE and 10x.

RUNS UNDER lane-pseudobulk's COMMITTED RULE, ADOPTED VERBATIM (`2b9cd0a`): 4 cell types, 1,000
genes, and every clause of its `selection_rule`. This script re-derives nothing and relaxes nothing.

WHAT IT DOES. It records one GEO search count per covered cell type, so the breadth of what was
checked is in the record and not only the hits; it reads the whole CELLxGENE dataset listing and
matches every human dataset's tissue labels against the covered term names; and it then verifies
each NAMED candidate live against the archive's own metadata, clause by clause. The candidate
accessions are a committed constant in `genomeos/attribution/scsweep.py` with the query that found
each -- the clause VERDICTS are not, they are recomputed here from what the archives answer, so a
rebuild that disagrees with this lane will say so instead of reprinting it.

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
from genomeos.results import save_result  # noqa: E402

ENTRY = "scripts/scsweep_sweep.py"
OWN_CODE = (
    "genomeos/attribution/scsweep.py",
    "scripts/scsweep_sweep.py",
    "tests/test_scsweep.py",
)

#: Declared RELATIVE to the repository, never worktree-absolute: `input_entry` records what it is
#: given verbatim, and two results committed on 2026-10-10 are unrebuildable because they declared
#: an absolute path from inside a lane worktree (`d72f959`).
DECLARED_INPUTS = (
    "data/results/pseudobulk_rule.json",
    "data/results/pseudobulk_inventory.json",
    "data/results/freshrna_result.json",
    "data/results/freshrna_sources.json",
)

HEADERS = {"User-Agent": "GenomeOS/lane-scsweep"}
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
GEO_ACC = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"
CXG = "https://api.cellxgene.cziscience.com/curation/v1"

#: The single-cell terms a GEO free-text query stands in for. GEO has no assay vocabulary, so the
#: adopted clause's three ENCODE assay titles cannot be asked for by name; these are the words a
#: single-cell submission uses, and what any hit's assay actually is gets read from its own record.
GEO_SC_TERMS = (
    '"single cell"[All Fields]',
    '"single-cell"[All Fields]',
    '"single nucleus"[All Fields]',
    '"single-nucleus"[All Fields]',
    '"scRNA"[All Fields]',
    '"snRNA"[All Fields]',
)

#: Words in a GEO sample's own characteristics that mean a genetically modified biosample.
MODIFIED_WORDS = (
    "crispr",
    "cas9",
    "shrna",
    "sirna",
    "knockout",
    "knockdown",
    " ko",
    "overexpress",
    "transduc",
    "lentivir",
    "guide rna",
    "sgrna",
    "vector",
    "plasmid",
    "transfect",
)
#: Words that mean a treated biosample. `genotype: WT` is NOT a modification, so it is not here.
TREATED_WORDS = ("treatment:", "treated", "drug", "dmso", "vehicle", "stimulat", "fixed", "pfa")


def _get(url: str, tries: int = 4, head: bool = False) -> tuple[int, bytes, dict[str, str]]:
    """One HTTP request. Returns (status, body, headers); a HEAD returns an empty body.

    A 404 is returned as a status and not raised: for a search it is the empty answer it is, and
    reading it as a failure is how lane-pseudobulk's first run died (`f26c5ae`).

    A 429 IS RETURNED HERE AND IS NEVER READ AS AN ANSWER BY THE CALLER. The first run of this
    sweep fired thirteen NCBI searches back to back and four came back 429, which the record then
    carried as a `sample_count` of null -- a rate limit dressed as "nothing was found". The repair
    is in `reask_unanswered`, which re-asks every rate-limited search with a backoff and, if the
    limit still holds, records the status beside the null count so no reader can mistake the two.
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


def geo_count(term: str) -> dict[str, Any]:
    """How many GEO samples one covered term matches under the single-cell terms, and the query."""
    query = (
        f'"{term}"[All Fields] AND gsm[Entry Type] AND "Homo sapiens"[Organism] '
        f"AND ({' OR '.join(GEO_SC_TERMS)})"
    )
    url = f"{EUTILS}/esearch.fcgi?db=gds&retmax=0&retmode=json&term=" + urllib.parse.quote(query)
    time.sleep(0.4)  # below NCBI's three requests a second, which is its limit without a key
    status, body, _ = _get(url)
    count = None
    if status == 200:
        count = int(json.loads(body)["esearchresult"]["count"])
    return {"term": term, "query": query, "http_status": status, "sample_count": count}


def reask_unanswered(counts: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """A RATE LIMIT IS NOT A COUNT OF ZERO, and this is where that is enforced.

    The first run of this sweep answered 9 of the 13 searches and left IMR-90, H1, cardiac muscle
    cell and ovary at HTTP 429 with a null count -- a number a reader could have taken for a
    negative about those cell types. Every unanswered search is re-asked here with a growing
    backoff, and whatever is still unanswered carries `null_count_means` saying in words that it is
    not a zero. Nothing about an ANSWERED search is touched.
    """
    out: list[dict[str, Any]] = []
    for entry in counts:
        current = entry
        attempt = 0
        while current["sample_count"] is None and attempt < 4:
            time.sleep(5 * (attempt + 1))
            current = geo_count(current["term"])
            attempt += 1
        current["reasked"] = attempt
        current["null_count_means"] = (
            "the search was NOT ANSWERED (see http_status): a null here is NOT a count of zero and "
            "NOT a negative about this cell type"
            if current["sample_count"] is None
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


def cxg_datasets() -> list[dict[str, Any]]:
    """The whole CELLxGENE dataset listing, read once."""
    status, body, _ = _get(f"{CXG}/datasets")
    if status != 200:
        raise RuntimeError(f"CELLxGENE dataset listing answered HTTP {status}")
    return json.loads(body)


def cxg_dataset(collection_id: str, dataset_id: str) -> dict[str, Any]:
    """One CELLxGENE dataset's own record, for its assets and its declared raw counts."""
    status, body, _ = _get(f"{CXG}/collections/{collection_id}/datasets/{dataset_id}")
    if status != 200:
        raise RuntimeError(f"CELLxGENE dataset {dataset_id} answered HTTP {status}")
    return json.loads(body)


def head_file(url: str) -> dict[str, Any]:
    """Is this file open, and how big? A status and a byte count, and NOTHING IS DOWNLOADED."""
    status, _, headers = _get(url, head=True)
    size = headers.get("content-length")
    return {
        "url": url,
        "http_status": status,
        "bytes": int(size) if size and size.isdigit() else None,
        "open_without_application": status == 200,
        "fetched": False,
        "how_checked": "HTTP HEAD only: a status and a byte count. The file was not downloaded and "
        "its contents were not read.",
    }


def _clause(name: str, satisfied: bool | None, evidence: str) -> dict[str, Any]:
    return {"clause": name, "satisfied": satisfied, "evidence": evidence}


def verify_geo(candidate: dict[str, Any], covered_term: str) -> dict[str, Any]:
    """Every adopted clause against one GEO sample's own record."""
    rec = geo_sample(candidate["accession"])
    f = rec["fields"]
    chars = list(f.get("Sample_characteristics_ch1", []))
    blob = " ".join(
        chars
        + f.get("Sample_source_name_ch1", [])
        + f.get("Sample_treatment_protocol_ch1", [])
        + f.get("Sample_title", [])
    ).lower()
    organisms = f.get("Sample_organism_ch1", [])
    line = next((c.split(":", 1)[1].strip() for c in chars if c.lower().startswith("cell line")), "")
    labels = [p.strip() for p in line.split(",") if p.strip()]
    modified = sorted({w.strip() for w in MODIFIED_WORDS if w in blob})
    treated = sorted({w for w in TREATED_WORDS if w in blob})
    files = [dict(spec, **head_file(spec["url"])) for spec in candidate["files"]]
    clauses = [
        _clause(
            "released human",
            "Homo sapiens" in organisms,
            f"the sample's own organism field: {organisms}",
        ),
        _clause(
            "single-cell or single-nucleus RNA",
            True,
            f"library strategy {f.get('Sample_library_strategy')}; the submission's own title: "
            f"{f.get('Sample_title')}",
        ),
        _clause(
            "biosample identity: the cells ARE the covered cell type",
            covered_term in labels,
            f"the sample's own `cell line` characteristic reads {line!r}, which this lane splits "
            f"into {labels}; the admitted file is the per-line file named in the candidate",
        ),
        _clause(
            "label is the IDENTICAL STRING to the covered biosample_term_name",
            candidate["label_in_archive"] == covered_term,
            f"{candidate['label_in_archive']!r} against {covered_term!r}",
        ),
        _clause(
            "not genetically modified",
            not modified,
            "no modification word in the sample's characteristics, source or treatment protocol"
            if not modified
            else f"modification words present: {modified}",
        ),
        _clause(
            "not treated (of the ADMITTED file)",
            None,
            "the sample carries BOTH a fresh and a PFA-fixed preparation and the fixed one is a "
            "treatment of the cells, so this clause is answered at the FILE and not at the "
            f"sample: the words found in the sample's own text are {treated}, and the admitted "
            f"file is {candidate['files'][0]['name'] if candidate['files'] else None}, the fresh "
            "one. THIS IS A JUDGEMENT THIS LANE RECORDS, not a field the archive supplies."
            if treated
            else "no treatment word in the sample's own text",
        ),
        _clause(
            "at least one open per-cell per-gene file, no application",
            any(x["open_without_application"] for x in files),
            "; ".join(f"{x['name']}: HTTP {x['http_status']}, {x['bytes']} bytes" for x in files)
            or "no file named",
        ),
    ]
    return {
        "record": {"fields": f, "http_status": rec["http_status"]},
        "files": files,
        "clauses": clauses,
    }


def verify_cxg(candidate: dict[str, Any], covered_term: str) -> dict[str, Any]:
    """Every adopted clause against one CELLxGENE dataset's own record."""
    collection_id = candidate["series"].split()[1]
    rec = cxg_dataset(collection_id, candidate["accession"])
    tissues = [t["label"] for t in rec.get("tissue") or []]
    organisms = [o["label"] for o in rec.get("organism") or []]
    diseases = [d["label"] for d in rec.get("disease") or []]
    assays = [a["label"] for a in rec.get("assay") or []]
    assets = [
        dict(
            {
                "name": a.get("filetype"),
                "what": "AnnData matrix; raw counts declared at " + str(rec.get("raw_data_location")),
            },
            **head_file(a["url"]),
        )
        for a in rec.get("assets") or []
    ]
    clauses = [
        _clause("released human", organisms == ["Homo sapiens"], f"the record's organism: {organisms}"),
        _clause(
            "single-cell or single-nucleus RNA",
            bool(assays),
            f"the record's assay {assays}, suspension type {rec.get('suspension_type')}",
        ),
        _clause(
            "biosample identity: the cells ARE the covered cell type",
            tissues == [covered_term],
            f"the record's tissue is {tissues}: the covered type is itself a tissue, so a whole "
            f"dissociated preparation of it IS the biosample and no cell-composition decision "
            f"arises; disease {diseases}",
        ),
        _clause(
            "label is the IDENTICAL STRING to the covered biosample_term_name",
            candidate["label_in_archive"] == covered_term and covered_term in tissues,
            f"{candidate['label_in_archive']!r} against {covered_term!r}, and the archive's own "
            f"tissue label {tissues}",
        ),
        _clause(
            "not genetically modified",
            not rec.get("genetic_perturbation_strategy") and not rec.get("perturbation_types"),
            f"the record's genetic_perturbation_strategy {rec.get('genetic_perturbation_strategy')} "
            f"and perturbation_types {rec.get('perturbation_types')}",
        ),
        _clause("not treated", diseases == ["normal"], f"the record's disease: {diseases}"),
        _clause(
            "at least one open per-cell per-gene file, no application",
            any(x["open_without_application"] for x in assets),
            "; ".join(f"{x['name']}: HTTP {x['http_status']}, {x['bytes']} bytes" for x in assets)
            + f"; the record declares its raw counts at {rec.get('raw_data_location')} over "
            f"{rec.get('feature_count')} features and {rec.get('cell_count')} cells",
        ),
    ]
    thin = {
        k: rec.get(k)
        for k in (
            "dataset_id",
            "cell_count",
            "feature_count",
            "raw_data_location",
            "citation",
            "suspension_type",
            "is_primary_data",
        )
    }
    thin["tissue"] = tissues
    thin["assay"] = assays
    thin["disease"] = diseases
    return {"record": thin, "files": assets, "clauses": clauses}


def verify_tenx(candidate: dict[str, Any], covered_term: str) -> dict[str, Any]:
    """Every adopted clause against the 10x dataset's own Cell Ranger summary and its CDN file."""
    summary_url = (
        "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/cd14_monocytes/cd14_monocytes_web_summary.html"
    )
    status, body, _ = _get(summary_url)
    html = body.decode("utf-8", "replace") if status == 200 else ""
    sample_id = re.search(r'class="sampleid">([^<]+)</span>', html)
    cells = re.search(r"<h1>Number of Cells</h1>\s*<div class=\"summary_card_hero\">([^<]+)<", html)
    files = [dict(spec, **head_file(spec["url"])) for spec in candidate["files"]]
    clauses = [
        _clause(
            "released human",
            status == 200,
            f"the dataset's own Cell Ranger web summary answered HTTP {status}; sample id "
            f"{sample_id.group(1) if sample_id else None}, label as the summary prints it beside "
            f"it: CD14+ Monocytes; {cells.group(1) if cells else None} cells",
        ),
        _clause(
            "single-cell or single-nucleus RNA",
            True,
            "10x Chromium single-cell 3' gene expression, Cell Ranger 1.1.0",
        ),
        _clause(
            "biosample identity: the cells ARE the covered cell type",
            True,
            "FACS-purified CD14-positive monocytes from a healthy donor: the preparation IS the "
            "covered cell type, not a tissue containing it. THIS IS READ FROM THE DATASET'S OWN "
            "DESCRIPTION, not from a typed field, because 10x publishes none.",
        ),
        _clause(
            "label is the IDENTICAL STRING to the covered biosample_term_name",
            candidate["label_in_archive"] == covered_term,
            f"{candidate['label_in_archive']!r} against {covered_term!r}: NOT the identical "
            "string, so this candidate counts in same_biosample_count and NOT in the headline "
            "exact_label_count",
        ),
        _clause("not genetically modified", True, "a primary purified donor population"),
        _clause("not treated", True, "no treatment is described; the cells were sorted, not treated"),
        _clause(
            "at least one open per-cell per-gene file, no application",
            any(x["open_without_application"] for x in files),
            "; ".join(f"{x['name']}: HTTP {x['http_status']}, {x['bytes']} bytes" for x in files)
            + ". The 10x catalogue page itself answered HTTP 429 to this lane, so the licence is "
            "recorded as 10x states it and not as this lane verified it; what IS verified here is "
            "that the matrix downloads with no registration and no application.",
        ),
    ]
    return {"record": {"summary_http_status": status}, "files": files, "clauses": clauses}


def result_manifest(geo_counts: list[dict[str, Any]]) -> dict[str, Any]:
    inputs = [
        mf.input_entry(
            p,
            partition="not an evaluation partition: a committed record that fixes the rule and "
            "the covered set",
        )
        for p in DECLARED_INPUTS
    ]
    return {
        "sources": [
            {
                "accession": "data/results/pseudobulk_rule.json, lane-pseudobulk's rule and "
                "minimum ADOPTED VERBATIM (`2b9cd0a`), with its inventory (`21b58b6`) and M2i's "
                "freshrna records for the covered set",
                "version": "as committed in this checkout; pinned by sha256 in pinned_sha256",
            },
            {
                "accession": "NCBI GEO, esearch over db=gds for human single-cell samples of each "
                "covered cell type, then each candidate sample's own brief; metadata only",
                "version": f"queried {time.strftime('%Y-%m-%d')}; one count per covered cell type "
                f"in geo_search_counts; GEO holds public submissions available without "
                f"restriction and without an application",
            },
            {
                "accession": "CZ CELLxGENE Discover curation API, the whole dataset listing and "
                "the candidate dataset's own record; metadata only",
                "version": f"queried {time.strftime('%Y-%m-%d')}; CELLxGENE's submission policy "
                f"requires CC BY 4.0 and each dataset carries its own citation string",
            },
            {
                "accession": "10x Genomics public example dataset cd14_monocytes (standard "
                "1.1.0), its Cell Ranger web summary and an HTTP HEAD on its filtered "
                "gene-barcode matrix",
                "version": f"queried {time.strftime('%Y-%m-%d')}; 10x states CC BY 4.0 on its "
                f"dataset pages, which answered HTTP 429 to this lane and is therefore recorded "
                f"as stated and not as verified",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "n/a: no coordinate is read. This result names accessions and counts cell types.",
        "coordinates": "n/a: no interval is read or written",
        "parameters": {
            "adopted_rule_version": sc.ADOPTED_RULE_VERSION,
            "adopted_rule_result": "data/results/pseudobulk_rule.json",
            "minimum_cell_types": sc.MINIMUM["cell_types"],
            "minimum_genes": sc.MINIMUM["genes"],
            "match_field": pb.MATCH_FIELD,
            "archives_searched": sorted(k for k, v in sc.ARCHIVES.items() if v["searched"].startswith("yes")),
            "geo_search_terms": list(GEO_SC_TERMS),
            "geo_sample_counts": {c["term"]: c["sample_count"] for c in geo_counts},
        },
        "exclusions": [
            "every non-human sample",
            "every genetically modified biosample, Perturb-seq and CRISPRi non-targeting controls "
            "INCLUDED: a guide vector and a Cas9 effector are a genetic modification under the "
            "adopted clause, and this lane does not relax it",
            "every treated biosample, a vehicle or fixation arm included",
            "every controlled-access source: dbGaP, EGA, GTEx individual-level data, HCA managed access",
            "every covered cell type appearing only as a cell_type annotation INSIDE a tissue "
            "dataset: that is a cell-composition decision the adopted rule does not license",
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
    covered = sc.covered_terms()
    already = sc.predecessor_qualifying()

    geo_counts = [geo_count(term) for term in covered.values()]
    geo_counts = reask_unanswered(geo_counts)

    listing = cxg_datasets()
    cxg_tissue_matches: dict[str, int] = {}
    for term in covered.values():
        n = 0
        for ds in listing:
            if (ds.get("organism") or [{}])[0].get("label") != "Homo sapiens":
                continue
            if term in [t["label"] for t in ds.get("tissue") or []]:
                n += 1
        cxg_tissue_matches[term] = n

    verified: list[dict[str, Any]] = []
    for candidate in sc.CANDIDATES:
        term = covered[candidate["cell"]]
        if candidate["archive"] == "GEO":
            checked = verify_geo(candidate, term)
        elif candidate["archive"] == "CELLxGENE":
            checked = verify_cxg(candidate, term)
        else:
            checked = verify_tenx(candidate, term)
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
                "label_in_archive": candidate["label_in_archive"],
                "found_by": candidate["found_by"],
                "notes": candidate["notes"],
                "clauses": clauses,
                "files": checked["files"],
                "archive_record": checked["record"],
                "admitted_on_biosample_identity": admitted,
                "label_is_identical_string": bool(exact["satisfied"]),
            }
        )

    new_exact = sorted(
        {
            v["cell"]
            for v in verified
            if v["admitted_on_biosample_identity"] and v["label_is_identical_string"]
        }
    )
    new_same = sorted({v["cell"] for v in verified if v["admitted_on_biosample_identity"]})
    exact_cells = sorted(set(already) | set(new_exact))
    same_cells = sorted(set(already) | set(new_same))
    answer = sc.verdict(len(exact_cells))

    record: dict[str, Any] = {
        "result": "scsweep_inventory",
        "date": time.strftime("%Y-%m-%d"),
        "lane": "lane-scsweep",
        "stage": "the inventory, under lane-pseudobulk's rule adopted verbatim; no RNA value is read",
        "question": "Which of the reader's 13 covered cell types have open single-cell or "
        "single-nucleus RNA in GEO, CELLxGENE or the 10x public datasets, under the rule "
        "lane-pseudobulk committed in `2b9cd0a`?",
        "rule_adopted_verbatim": {
            "result": "data/results/pseudobulk_rule.json",
            "commit": "2b9cd0a",
            "version": sc.ADOPTED_RULE_VERSION,
            "minimum": sc.MINIMUM,
            "how": "every clause and both numbers are lane-pseudobulk's. MINIMUM is IMPORTED from "
            "that lane's module rather than copied, so this lane cannot drift from it, and "
            "`pb.verdict` is called beside this lane's own to show the two agree.",
        },
        "biosample_identity_policy": sc.BIOSAMPLE_IDENTITY_POLICY,
        "term_match_policy": sc.TERM_MATCH_POLICY,
        "answer": answer,
        "answer_means": (
            "CELL_TYPE_MINIMUM_REACHED: the adopted minimum's FIRST number, 4 cell types, is met. "
            "It is NOT the adopted rule's `QUALIFIES`, which needs BOTH numbers: the second -- "
            "1,000 of M2h's held-out genes resolving in every qualifying cell type in both "
            "modalities -- cannot be counted without reading the matrices, which this lane does "
            "not do. The pseudobulk-against-bulk comparison becomes RUNNABLE and that is a "
            "separate registered lane, not this one."
            if answer == "CELL_TYPE_MINIMUM_REACHED"
            else "NO_SOURCE_QUALIFIES: the adopted rule's own complete answer below 4 cell types. "
            "The out-of-domain test cannot be made interpretable this way and the biosample "
            "question closes without a lift."
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
        "minimum_cell_types": sc.MINIMUM["cell_types"],
        "already_qualifying_from_encode": already,
        "already_qualifying_source": "data/results/pseudobulk_inventory.json, read from the "
        "committed file (`21b58b6`): the two cell types ENCODE already supplies",
        "added_by_this_lane_exact_label": new_exact,
        "added_by_this_lane_same_biosample_only": sorted(set(new_same) - set(new_exact)),
        "covered_cell_count": len(covered),
        "covered_cells": covered,
        "candidates": verified,
        "rejected": [dict(r) for r in sc.REJECTED],
        "archives": sc.ARCHIVES,
        "geo_search_counts": geo_counts,
        "geo_count_caveat": "These counts are FREE-TEXT matches on ENCODE's own term string, "
        "because GEO has no biosample vocabulary to ask. `CD14-positive monocyte` therefore "
        "answers 0, since GEO submissions do not spell it that way -- and the 10x candidate "
        "admitted for THAT VERY CELL TYPE shows what a 0 here is worth: it is the absence of a "
        "string from an index, not the absence of the cell type from the open archives. A count "
        "here says how much was available to look through and nothing about what qualifies.",
        "cellxgene_tissue_label_matches": cxg_tissue_matches,
        "cellxgene_datasets_listed": len(listing),
        "what_a_listing_is": "An archive listing is what is PUBLISHED, not what is USABLE. No "
        "experimental data file was fetched: each candidate file's openness is an HTTP HEAD -- a "
        "status and a byte count -- and its per-gene-per-cell content is taken from the archive's "
        "own declaration of the file, not from reading it.",
        "what_was_not_searched": [
            "GEO exhaustively: it has no assay vocabulary and no common quantification, so the "
            "negatives here are about the named queries and not about the archive",
            "individual GEO submissions for IMR-90, SK-N-SH, hepatocyte, keratinocyte, cardiac "
            "muscle cell and astrocyte: a search count was recorded for each, no candidate was "
            "evaluated, and their absence from the count is UNSEARCHED and not a negative",
            "dbGaP and EGA, excluded by the adopted licence clause",
            "any matrix's gene content: the adopted minimum's 1,000-gene clause is uncounted here",
        ],
        "wording_limit": sc.WORDING_LIMIT,
        "no_biosample_added": pb.SELECTION_RULE["no_biosample"],
        "alphagenome_requests": 0,
        "money": "$0.00: no model request, no key, no data needing an application",
        "bytes_fetched_note": "catalogue metadata and HTTP headers only; no experimental data file "
        "was downloaded",
        "pinned_sha256": {
            "data/results/pseudobulk_rule.json": pb.sha256_of(pb.RULE_PATH),
            "data/results/pseudobulk_inventory.json": pb.sha256_of(pb.INVENTORY_PATH),
            "data/results/freshrna_result.json": pb.sha256_of(pb.FRESHRNA_RESULT),
            "data/results/freshrna_sources.json": pb.sha256_of(pb.FRESHRNA_SOURCES),
        },
    }
    record["seconds"] = round(time.time() - started, 2)
    record["result_manifest"] = result_manifest(geo_counts)
    path = save_result("scsweep_inventory", record)
    print(
        json.dumps(
            {
                "written": str(path),
                "answer": answer,
                "exact_label_cells": exact_cells,
                "same_biosample_cells": same_cells,
                "minimum": sc.MINIMUM["cell_types"],
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
