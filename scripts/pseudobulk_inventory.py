#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""The inventory: which covered cell types have BOTH open bulk RNA-seq and open single-nucleus RNA?

RUNS AFTER `data/results/pseudobulk_rule.json` IS COMMITTED and before any measured value is read.
This script reads catalogue metadata only: it fetches no experimental data file, makes no model
request and spends nothing.

HOW IT COUNTS. The bulk clause is satisfied by construction for all thirteen covered cell types --
M2i already holds one open ENCODE bulk RNA-seq experiment for each, on disk with its sha256 -- so
the only open question is the single-nucleus clause. For each of the three ENCODE single-cell assay
titles the portal's released human experiments are listed exhaustively, their
`biosample_ontology.term_name` is matched EXACTLY against the covered term names, and a matching
experiment is admitted only if the biosample carries no applied modification and no treatment and
the experiment has at least one OPEN processed file from which a per-gene value per nucleus or per
cell can be read. The output types actually found are recorded rather than assumed.

THE ANSWER MAY BE `NO_SOURCE_QUALIFIES`, which the committed rule declares acceptable and COMPLETE.
"""

from __future__ import annotations

import json
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
from genomeos.results import save_result  # noqa: E402

ENTRY = "scripts/pseudobulk_inventory.py"
OWN_CODE = (
    "genomeos/attribution/pseudobulk.py",
    "scripts/pseudobulk_rule.py",
    "scripts/pseudobulk_inventory.py",
    "tests/test_pseudobulk.py",
)

HEADERS = {"Accept": "application/json", "User-Agent": "GenomeOS/lane-pseudobulk"}

#: An open ENCODE file: anything else is behind a restriction the licence clause excludes.
OPEN_FILE_STATUS = ("released",)

#: A processed file counts for the single-cell clause when its output_type says, in ENCODE's own
#: words, that it holds per-gene values per cell or nucleus. Matched on substrings so the portal's
#: exact phrasing is DISCOVERED and recorded, not assumed.
MATRIX_WORDS = ("matrix", "gene count", "gene quantification", "transcript quantification")


# ENCODE ANSWERS AN EMPTY SEARCH WITH HTTP 404, NOT WITH AN EMPTY `@graph`. Measured here: the first
# run of this script died after four retries on `assay_title=long read scRNA-seq` restricted to
# released human experiments, which is not a portal failure but a count of zero. `fetch` therefore
# reads a 404 as the empty result it is, and records `encode_http_status` so the record says that the
# zero came from a 404 rather than from an answered search with no rows.
#
# THAT MAKES ONE OTHER CHECK NECESSARY, and `open_matrix_files` carries it: a per-experiment fetch
# whose accession came from a search that had just returned it must NOT read a 404 as "this
# experiment has no files", because that would turn a portal error into a cell type that looks like
# it fails the single-cell clause. It asserts the record came back with its own accession in it.
def fetch(url: str, tries: int = 4) -> dict[str, Any]:
    """One ENCODE portal JSON request, retried on a transient failure. No key, no cost."""
    last: Exception | None = None
    for attempt in range(tries):
        try:
            req = urllib.request.Request(url, headers=HEADERS)
            with urllib.request.urlopen(req, timeout=180) as handle:
                return json.load(handle)
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return {"@graph": [], "total": 0, "encode_http_status": 404}
            last = exc
            time.sleep(2 * (attempt + 1))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"ENCODE request failed after {tries} tries: {url}: {last}")


def search(assay_title: str) -> list[dict[str, Any]]:
    """Every released human experiment of one single-cell assay title, with the fields the rule needs."""
    fields = [
        "accession",
        "assay_title",
        "status",
        "biosample_ontology.term_name",
        "biosample_ontology.term_id",
        "biosample_ontology.classification",
        "biosample_summary",
        "lab.title",
        "date_released",
        "replicates.library.biosample.applied_modifications",
        "replicates.library.biosample.treatments",
    ]
    query = [
        ("type", "Experiment"),
        ("status", "released"),
        ("replicates.library.biosample.donor.organism.scientific_name", "Homo sapiens"),
        ("assay_title", assay_title),
        ("limit", "all"),
        ("format", "json"),
    ] + [("field", f) for f in fields]
    url = f"{pb.ENCODE}/search/?" + urllib.parse.urlencode(query)
    return fetch(url).get("@graph", [])


def modified_or_treated(experiment: dict[str, Any]) -> list[str]:
    """Why a biosample is rejected, in words, or an empty list when it is the cell type's own RNA."""
    reasons: list[str] = []
    for rep in experiment.get("replicates", []):
        sample = (rep.get("library") or {}).get("biosample") or {}
        if sample.get("applied_modifications"):
            reasons.append("applied_modifications present: a genetically modified biosample")
        if sample.get("treatments"):
            reasons.append("treatments present: not the untreated cell type's own RNA")
    return sorted(set(reasons))


def open_matrix_files(accession: str) -> tuple[list[dict[str, Any]], list[str]]:
    """The open per-cell matrix files of one experiment, and every output_type seen in it."""
    record = fetch(f"{pb.ENCODE}/experiments/{accession}/?format=json")
    if record.get("accession") != accession:
        raise RuntimeError(
            f"{accession}: the portal did not return this experiment "
            f"(encode_http_status {record.get('encode_http_status')}). The accession came from a "
            f"search that had just returned it, so this is a portal fault and it must NOT be read "
            f"as 'the experiment has no files', which would make the cell type look like it fails "
            f"the single-cell clause"
        )
    seen: list[str] = []
    out: list[dict[str, Any]] = []
    for f in record.get("files", []):
        output_type = f.get("output_type") or ""
        seen.append(output_type)
        if f.get("status") not in OPEN_FILE_STATUS:
            continue
        if f.get("restricted") or f.get("no_file_available"):
            continue
        if not any(word in output_type.lower() for word in MATRIX_WORDS):
            continue
        out.append(
            {
                "accession": f.get("accession"),
                "output_type": output_type,
                "file_format": f.get("file_format"),
                "assembly": f.get("assembly"),
                "genome_annotation": f.get("genome_annotation"),
                "size_bytes": f.get("file_size"),
                "url": pb.ENCODE + (f.get("href") or ""),
            }
        )
    return out, sorted(set(seen))


OTHER_SOURCES = [
    {
        "source": "Tabula Sapiens (CZ Biohub)",
        "licence": "CC BY 4.0, open, no application",
        "clause": "single_cell clause: PASSES on licence and modality, FAILS on the covered set",
        "detail": (
            "Tabula Sapiens publishes single-cell RNA across human tissues under CC BY 4.0. Its "
            "annotations are CL cell-type labels inside tissues, and none of them is an ENCODE "
            "immortalised cell line, so the six cell lines of the covered set (K562, HepG2, "
            "GM12878, IMR-90, H1, SK-N-SH) cannot be matched at all. For the tissue members of the "
            "covered set the labels are cell types WITHIN the tissue rather than the tissue, so "
            "matching one to ENCODE's `ovary` or `testis` bulk would require a cell-composition "
            "decision this lane has not registered and which is itself the confound under test."
        ),
    },
    {
        "source": "Human Cell Atlas Data Portal",
        "licence": "open releases are open; matched tiers are managed access",
        "clause": "licence clause for the managed tiers; covered set for the open ones",
        "detail": (
            "The same finding lane-oodinventory recorded (`171f5bf`), re-examined against THIS "
            "lane's clauses: this lane needs no matched chromatin accessibility, so HCA's open "
            "single-cell RNA releases are not excluded by that. They are primary human tissues and "
            "organoids annotated by CL cell type, so they meet no ENCODE cell line in the covered "
            "set, and a tissue-level match carries the same composition decision as Tabula Sapiens."
        ),
    },
    {
        "source": "GTEx single-nucleus (snRNA-seq pilot, eight tissues)",
        "licence": "summary matrices open; individual-level data dbGaP-controlled",
        "clause": "licence clause, and the covered set",
        "detail": (
            "GTEx's single-nucleus pilot covers lung, breast, oesophagus, heart, skeletal muscle, "
            "prostate, skin and oesophagus muscularis. None of those is a covered cell type under "
            "the exact-term match, and the individual-level data is dbGaP-controlled, which the "
            "licence clause excludes."
        ),
    },
    {
        "source": "GEO / SRA",
        "licence": "varies by submission",
        "clause": "recorded, not searched exhaustively",
        "detail": (
            "GEO is an archive of independent submissions with no common pipeline or common "
            "quantification. A single-nucleus submission on a covered cell line would be "
            "admissible under this lane's clauses, and the limitation is recorded plainly: this "
            "inventory did not search GEO exhaustively, so its answer is about ENCODE and the four "
            "named atlases and is not a claim that no such submission exists anywhere."
        ),
    },
]


def result_manifest(read_now: list[str], counts: dict[str, int]) -> dict[str, Any]:
    inputs = [
        mf.input_entry(
            p, partition="not an evaluation partition: a committed record that fixes the covered set"
        )
        for p in sorted(set(read_now))
        if "/data/" in p.replace("\\", "/") or p.startswith("data/")
    ]
    return {
        "sources": [
            {
                "accession": "data/results/pseudobulk_rule.json, this lane's committed rule and "
                "minimum, with data/results/freshrna_result.json and freshrna_sources.json",
                "version": "as committed in this checkout; pinned by sha256 in pinned_sha256",
            },
            {
                "accession": "ENCODE portal search API, released human experiments of "
                + ", ".join(pb.SINGLE_CELL_ASSAYS)
                + "; metadata only, no experimental file fetched",
                "version": f"queried {time.strftime('%Y-%m-%d')}; experiment counts per assay in "
                f"encode_experiment_counts; {pb.ENCODE_LICENCE}",
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": "n/a: no coordinate is read. This result counts experiments and names accessions.",
        "coordinates": "n/a: no interval is read or written",
        "parameters": {
            "selection_rule_version": pb.SELECTION_RULE["version"],
            "match_field": pb.MATCH_FIELD,
            "single_cell_assays": list(pb.SINGLE_CELL_ASSAYS),
            "minimum_cell_types": pb.MINIMUM["cell_types"],
            "minimum_genes": pb.MINIMUM["genes"],
            "encode_experiment_counts": counts,
        },
        "exclusions": [
            "every non-human experiment, and every experiment not status=released",
            "every biosample with an applied modification or a treatment",
            "every file that is not released, or is restricted, or has no file available",
            "every controlled-access source: GTEx individual-level data, dbGaP, HCA managed access",
            "every experimental data file: this search reads catalogue metadata only",
            "no biosample is added and no out-of-domain test is run (2026-09-27 stop list)",
        ],
        "partitions": "n/a: nothing is fitted and nothing is held out; the output is a list",
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


def main() -> int:
    mf.trace_begin()
    started = time.time()
    covered = pb.covered_cells()
    by_term = {c[pb.MATCH_FIELD]: c for c in covered}

    counts: dict[str, int] = {}
    hits: list[dict[str, Any]] = []
    all_terms: dict[str, set[str]] = {}
    for assay in pb.SINGLE_CELL_ASSAYS:
        rows = search(assay)
        counts[assay] = len(rows)
        terms = set()
        for row in rows:
            term = (row.get("biosample_ontology") or {}).get("term_name") or ""
            terms.add(term)
            if term in by_term:
                hits.append({"assay_title": assay, "term": term, "experiment": row})
        all_terms[assay] = terms

    per_cell: list[dict[str, Any]] = []
    output_types_seen: set[str] = set()
    for entry in covered:
        term = entry[pb.MATCH_FIELD]
        mine = [h for h in hits if h["term"] == term]
        admitted: list[dict[str, Any]] = []
        rejected: list[dict[str, Any]] = []
        for hit in mine:
            exp = hit["experiment"]
            accession = exp["accession"]
            reasons = modified_or_treated(exp)
            if reasons:
                rejected.append({"experiment": accession, "assay_title": hit["assay_title"], "why": reasons})
                continue
            files, seen = open_matrix_files(accession)
            output_types_seen.update(seen)
            if not files:
                rejected.append(
                    {
                        "experiment": accession,
                        "assay_title": hit["assay_title"],
                        "why": [
                            "no open processed file holds a per-gene value per nucleus or per cell; "
                            "output types present: " + ", ".join(seen)
                        ],
                    }
                )
                continue
            admitted.append(
                {
                    "experiment": accession,
                    "experiment_url": f"{pb.ENCODE}/experiments/{accession}/",
                    "assay_title": hit["assay_title"],
                    "biosample_summary": exp.get("biosample_summary"),
                    "lab": (exp.get("lab") or {}).get("title"),
                    "date_released": exp.get("date_released"),
                    "files": files,
                }
            )
        per_cell.append(
            {
                "cell": entry["cell"],
                pb.MATCH_FIELD: term,
                "bulk": {
                    "assay": entry["bulk_assay"],
                    "experiment": entry["bulk_experiment"],
                    "files": entry["bulk_files"],
                    "licence": pb.ENCODE_LICENCE,
                    "satisfied": True,
                },
                "single_cell": {
                    "admitted": admitted,
                    "rejected": rejected,
                    "satisfied": bool(admitted),
                },
                "has_both": bool(admitted),
            }
        )

    qualifying = [c for c in per_cell if c["has_both"]]
    answer = pb.verdict(len(qualifying), None if len(qualifying) < pb.MINIMUM["cell_types"] else 0)

    record: dict[str, Any] = {
        "result": "pseudobulk_inventory",
        "date": time.strftime("%Y-%m-%d"),
        "lane": "lane-pseudobulk",
        "stage": "the inventory the committed rule was written for, before any value is read",
        "rule_result": "data/results/pseudobulk_rule.json",
        "question": pb.SELECTION_RULE["question"],
        "minimum": pb.MINIMUM,
        "answer": answer,
        "qualifying_cell_count": len(qualifying),
        "qualifying_cells": [c["cell"] for c in qualifying],
        "covered_cell_count": len(covered),
        "per_cell": per_cell,
        "encode_experiment_counts": counts,
        "encode_terms_per_assay": {k: len(v) for k, v in all_terms.items()},
        "encode_output_types_seen": sorted(output_types_seen),
        "other_sources_searched": OTHER_SOURCES,
        "bulk_clause_satisfied_for_all_thirteen": True,
        "licence": pb.ENCODE_LICENCE,
        "wording_limit": pb.WORDING_LIMIT,
        "no_biosample_added": pb.SELECTION_RULE["no_biosample"],
        "alphagenome_requests": 0,
        "money": "$0.00: no model request, no key, no data needing an application",
        "bytes_fetched": 0,
        "bytes_note": "catalogue metadata only; no experimental data file was downloaded",
        "pinned_sha256": {
            "data/results/pseudobulk_rule.json": pb.sha256_of(pb.RULE_PATH),
            "data/results/freshrna_result.json": pb.sha256_of(pb.FRESHRNA_RESULT),
            "data/results/freshrna_sources.json": pb.sha256_of(pb.FRESHRNA_SOURCES),
        },
    }
    record["seconds"] = round(time.time() - started, 2)
    record["result_manifest"] = result_manifest(mf.traced_reads(), counts)
    path = save_result("pseudobulk_inventory", record)
    print(
        json.dumps(
            {
                "written": str(path),
                "answer": answer,
                "qualifying": record["qualifying_cells"],
                "counts": counts,
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
