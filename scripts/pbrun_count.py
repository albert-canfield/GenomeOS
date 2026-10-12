#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""Count the 1,000-gene clause FIRST, before any registration and before any comparison value.

lane-scsweep (`12b0c4e`) reached lane-pseudobulk's 4-cell-type minimum and returned
`CELL_TYPE_MINIMUM_REACHED` rather than `QUALIFIES`, because the second number -- 1,000 of M2h's
held-out genes resolving in every qualifying cell type in BOTH modalities -- was uncounted. This
script counts it, and the lane's whole next step hangs on the answer: below 1,000 the answer is
`NO_SOURCE_QUALIFIES`, the modality swap stays unlicensed, the bile-duct out-of-domain test stays
uninterpretable, and Albert's biosample question closes without a lift.

WHAT IS READ AND WHAT IS NOT. On the bulk side, the gene_id column of the exact ENCODE RSEM tables
M2i already read (`data/cache/freshrna/rna/`), and nothing else from them: no TPM, no count, no
value. On the single-cell side, the GENE AXIS ONLY of each chosen matrix -- a STARsolo features.tsv,
an RSEM gene_id column, an h5ad var index -- and no matrix entry at all. So this script reads no
measured value on either side, which is why it is allowed to run before the registration.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import pbrun as pr  # noqa: E402
from genomeos.attribution import pseudobulk as pb  # noqa: E402
from genomeos.results import save_result  # noqa: E402

ENTRY = "scripts/pbrun_count.py"
OWN_CODE = (
    "genomeos/attribution/pbrun.py",
    "scripts/pbrun_count.py",
    "tests/test_pbrun.py",
)

#: The single-cell sources, per cell type: the archive, the file on disk relative to the repository,
#: how its gene axis is read and which identifier it carries. GM12878 and ovary are derived from
#: lane-pseudobulk's committed inventory by SC_SELECTION['encode']; K562 and testis are the one file
#: lane-scsweep admitted for each.
GEO_K562 = pr.CACHE_DIR / "GSM9224207_Fresh_K562.h5ad"
CXG_TESTIS = pr.CACHE_DIR / "09e3ae14-f355-4ed7-ad81-787b587c8a16.h5ad"

NON_ENCODE_SOURCES: dict[str, dict[str, Any]] = {
    "K562": {
        "archive": "GEO",
        "accession": "GSM9224207",
        "series": "GSE306708",
        "assay": "inCAGE single-cell RNA-seq (droplet)",
        "modality": "single-cell",
        "path": str(GEO_K562),
        "url": (
            "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9224nnn/GSM9224207/suppl/GSM9224207_Fresh_K562.h5ad"
        ),
        "reader": "h5ad",
        "gene_axis": "var/_index",
        "counts_at": "X",
        "identifier": pr.SYMBOL,
        "licence": (
            "NCBI GEO holds public submissions available without restriction and without an "
            "application; no controlled-access tier is involved"
        ),
    },
    "testis": {
        "archive": "CELLxGENE Discover",
        "accession": "ec6ef004-5f30-4a40-b282-980cca20a561",
        "series": "Tabula Sapiens - Testis",
        "assay": "10x 3' v3 single-cell RNA-seq",
        "modality": "single-cell",
        "path": str(CXG_TESTIS),
        "url": "https://datasets.cellxgene.cziscience.com/09e3ae14-f355-4ed7-ad81-787b587c8a16.h5ad",
        "reader": "h5ad",
        "gene_axis": "raw/var/_index",
        "counts_at": "raw/X",
        "identifier": pr.ENSEMBL,
        "licence": (
            "CZ CELLxGENE Discover distributes curated datasets openly over HTTPS with no "
            "registration and no application; its submission policy requires CC BY 4.0"
        ),
    },
}

#: The output types SC_SELECTION['encode'] prefers, best first.
ENCODE_OUTPUT_PRIORITY = ("sparse gene count matrix of unique reads", "gene quantifications")
#: An ENCODE scRNA-seq 'gene quantifications' file is ONE CELL, so all of an experiment's qualifying
#: files contribute; a STARsolo matrix is a whole experiment, so exactly one is taken.
ONE_FILE_PER_EXPERIMENT = ("sparse gene count matrix of unique reads",)


def annotation_key(value: str | None) -> tuple[int, str]:
    """GENCODE's version as a number, so V29 sorts above V24 and above a missing annotation."""
    text = (value or "").lstrip("Vv")
    return (int(text) if text.isdigit() else -1, value or "")


def encode_sources(cell: str) -> list[dict[str, Any]]:
    """SC_SELECTION['encode'] applied to lane-pseudobulk's committed inventory, for one cell type."""
    inventory = json.loads(pr.PSEUDOBULK_INVENTORY.read_text())
    entry = next(e for e in inventory["per_cell"] if e["cell"] == cell)
    chosen: list[dict[str, Any]] = []
    for experiment in entry["single_cell"]["admitted"]:
        files = [f for f in experiment["files"] if f.get("assembly") == "GRCh38"]
        if not files:
            continue
        best = max(annotation_key(f.get("genome_annotation")) for f in files)
        files = [f for f in files if annotation_key(f.get("genome_annotation")) == best]
        picked: list[dict[str, Any]] = []
        for output_type in ENCODE_OUTPUT_PRIORITY:
            same = sorted((f for f in files if f["output_type"] == output_type), key=lambda f: f["accession"])
            if same:
                picked = same[:1] if output_type in ONE_FILE_PER_EXPERIMENT else same
                break
        for f in picked:
            is_tar = f["file_format"] == "tar"
            chosen.append(
                {
                    "archive": "ENCODE",
                    "accession": f["accession"],
                    "experiment": experiment["experiment"],
                    "assay": experiment["assay_title"],
                    "modality": "single-nucleus"
                    if experiment["assay_title"] == "snRNA-seq"
                    else "single-cell",
                    "output_type": f["output_type"],
                    "assembly": f["assembly"],
                    "genome_annotation": f["genome_annotation"],
                    "path": str(pr.CACHE_DIR / f"{f['accession']}{'.tar.gz' if is_tar else '.tsv'}"),
                    "url": f["url"],
                    "bytes_declared": f["size_bytes"],
                    "reader": "starsolo_features_tsv" if is_tar else "rsem_gene_id_column",
                    "gene_axis": "GeneFull_Ex50pAS/filtered/features.tsv" if is_tar else "gene_id column",
                    "counts_at": "matrix.mtx" if is_tar else "expected_count column",
                    "identifier": pr.ENSEMBL,
                    "licence": pb.ENCODE_LICENCE,
                }
            )
    return chosen


def single_cell_sources(cells: list[str]) -> dict[str, list[dict[str, Any]]]:
    """Every contributing single-cell file per cell type, by SC_SELECTION and nothing else."""
    out: dict[str, list[dict[str, Any]]] = {}
    for cell in cells:
        if cell in NON_ENCODE_SOURCES:
            out[cell] = [dict(NON_ENCODE_SOURCES[cell])]
        else:
            out[cell] = encode_sources(cell)
    return out


def gene_axis_counts(source: dict[str, Any]) -> dict[str, int]:
    """The rows per identifier of one contributing file's GENE AXIS. No matrix entry is read."""
    path = ROOT / source["path"]
    reader = source["reader"]
    if reader == "rsem_gene_id_column":
        return pr.row_counts_rsem_tsv(path)
    if reader == "starsolo_features_tsv":
        return pr.row_counts_starsolo_tar(path)
    if reader == "h5ad":
        return pr.row_counts_h5ad(path, source["gene_axis"])
    raise ValueError(f"no reader named {reader!r}")


def bulk_sources(cells: list[str]) -> dict[str, list[dict[str, Any]]]:
    """The exact bulk files M2i read, per cell type, from its committed sources record."""
    biosamples = json.loads(pr.FRESHRNA_SOURCES.read_text())["biosamples"]
    out: dict[str, list[dict[str, Any]]] = {}
    for cell in cells:
        entry = biosamples[cell]
        out[cell] = [
            {
                "accession": f["accession"],
                "path": f["path"],
                "sha256": f["sha256"],
                "experiment": entry["experiment"],
                "assay": entry["assay"],
            }
            for f in entry["files"]
        ]
    return out


def result_manifest(read_now: list[str], parameters: dict[str, Any]) -> dict[str, Any]:
    inputs = [
        mf.input_entry(
            p,
            partition=(
                "M2h's HELD-OUT side where the input is 7h's genes table; not an evaluation "
                "partition for the rest"
            ),
        )
        for p in sorted(set(read_now))
        if "/data/" in p.replace("\\", "/") or p.startswith("data/")
    ]
    return {
        "sources": [
            {
                "accession": (
                    "data/results/multrna_result.json (7h's genes table), "
                    "data/results/magcal_result.json (M2h), "
                    "data/results/freshrna_result.json and freshrna_sources.json (M2i), "
                    "data/results/pseudobulk_inventory.json, data/results/scsweep_inventory.json"
                ),
                "version": "as committed in this checkout; pinned by sha256 in pinned_sha256",
            },
            {
                "accession": (
                    "ENCODE scRNA-seq and snRNA-seq processed files of GM12878 and ovary, as named "
                    "in single_cell_sources"
                ),
                "version": f"downloaded {time.strftime('%Y-%m-%d')}; {pb.ENCODE_LICENCE}",
            },
            {
                "accession": "GEO GSM9224207 (GSE306708) GSM9224207_Fresh_K562.h5ad",
                "version": f"downloaded {time.strftime('%Y-%m-%d')}; "
                + NON_ENCODE_SOURCES["K562"]["licence"],
            },
            {
                "accession": (
                    "CELLxGENE Discover dataset ec6ef004-5f30-4a40-b282-980cca20a561, Tabula Sapiens - Testis"
                ),
                "version": f"downloaded {time.strftime('%Y-%m-%d')}; "
                + NON_ENCODE_SOURCES["testis"]["licence"],
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": (
            "GRCh38 on every side; the GM12878 single-cell tables are GENCODE V24 and everything "
            "else is V29, and genes are matched on the UNVERSIONED Ensembl id or on the symbol, so "
            "the annotation difference costs genes and never mismatches one"
        ),
        "coordinates": "n/a: no interval is read or written. This result counts genes.",
        "parameters": parameters,
        "exclusions": [
            "every gene 7h did not admit, and every gene M2h's split puts on the TRAINING side",
            "a gene_id suffixed _PAR_Y (lane-pseudobulk's clause, carried)",
            "a gene with more than one row, or no row, in any contributing file on either side",
            "the PFA-fixed K562 file of the same GEO sample: a treatment of the cells",
            "every normalised, decontaminated, imputed or scaled layer of either h5ad",
            "every measured value: this script reads gene axes only",
        ],
        "code_revision": None,
        "entry": ENTRY,
        "own_code": list(OWN_CODE),
    }


def main() -> int:
    started = time.time()
    cells = pr.qualifying_cells()
    population = pr.held_out_population()
    keys = {
        pr.ENSEMBL: [g["gene_id"] for g in population],
        pr.SYMBOL: [g["symbol"] for g in population],
    }

    bulk = bulk_sources(cells)
    single = single_cell_sources(cells)

    bulk_counts = {c: [pr.row_counts_rsem_tsv(ROOT / f["path"]) for f in bulk[c]] for c in cells}
    single_counts = {c: [gene_axis_counts(s) for s in single[c]] for c in cells}

    per_cell: list[dict[str, Any]] = []
    resolved_sets: dict[str, set[int]] = {}
    for cell in cells:
        identifiers = {s["identifier"] for s in single[cell]}
        if len(identifiers) != 1:
            raise ValueError(f"{cell}: contributing files disagree on the identifier: {identifiers}")
        identifier = identifiers.pop()
        bulk_ok: set[int] = set()
        single_ok: set[int] = set()
        both: set[int] = set()
        for i, gene in enumerate(population):
            if pr.resolves(bulk_counts[cell], gene["gene_id"]):
                bulk_ok.add(i)
            if pr.resolves(single_counts[cell], keys[identifier][i]):
                single_ok.add(i)
            if i in bulk_ok and i in single_ok:
                both.add(i)
        resolved_sets[cell] = both
        per_cell.append(
            {
                "cell": cell,
                "identifier_matched_on": identifier,
                "modality": sorted({s["modality"] for s in single[cell]}),
                "contributing_single_cell_files": len(single[cell]),
                "bulk_files": len(bulk[cell]),
                "bulk_rows_resolving": len(bulk_ok),
                "single_cell_rows_resolving": len(single_ok),
                "usable_on_both_sides": len(both),
                "single_cell_gene_axis_rows": [len(c) for c in single_counts[cell]],
                "single_cell_sources": single[cell],
                "bulk_sources": bulk[cell],
            }
        )

    shared = set.intersection(*(resolved_sets[c] for c in cells)) if cells else set()
    shared_genes = sorted(population[i]["symbol"] for i in shared)
    answer = pr.verdict(len(cells), len(shared))

    record: dict[str, Any] = {
        "result": "pbrun_gene_count",
        "date": time.strftime("%Y-%m-%d"),
        "lane": "lane-pbrun",
        "stage": (
            "step 1 of 3: the 1,000-gene clause counted BEFORE the registration and before any "
            "comparison value is read"
        ),
        "question": (
            "How many of M2h's held-out genes of the multiplicative form have a usable value on BOTH "
            "sides -- bulk RNA-seq and pseudobulked single-cell or single-nucleus RNA -- in EVERY "
            "one of lane-scsweep's four qualifying cell types?"
        ),
        "gene_clause_counted_verbatim": pr.GENE_CLAUSE,
        "minimum": pr.MINIMUM,
        "minimum_is_lane_pseudobulks": pr.MINIMUM is pb.MINIMUM,
        "answer": answer,
        "answer_means": (
            "QUALIFIES: both of lane-pseudobulk's numbers are met, so the registration and the run "
            "follow. NO_SOURCE_QUALIFIES: the lane ENDS HERE, the modality swap stays unlicensed, "
            "the bile-duct out-of-domain test stays uninterpretable, and Albert's biosample question "
            "closes without a lift -- which lane-pseudobulk's rule declares an acceptable and "
            "COMPLETE answer."
        ),
        "qualifying_cell_count": len(cells),
        "qualifying_cells": cells,
        "cells_source": pr.CELLS_SOURCE,
        "form": pr.FORM,
        "form_reason": (
            "the one form of M2h whose reading is 'better' and whose scale is the +0.2146 this lane "
            "may use; the mean form is inconclusive and the promoter-competition form is a negative"
        ),
        "held_out_population": len(population),
        "qualifying_gene_count": len(shared),
        "qualifying_genes": shared_genes,
        "per_cell": per_cell,
        "sc_selection": pr.SC_SELECTION,
        "aggregation_to_be_registered": pr.AGGREGATION,
        "nothing_fitted": (
            "no value of any kind is read here, so nothing can be fitted. M2h's scale appears in "
            "this record only as the constant it is."
        ),
        "scale_constant": pr.SCALE,
        "no_value_read": (
            "the GENE AXIS ONLY of every file on both sides: a gene_id column, a STARsolo "
            "features.tsv, an h5ad var index. No TPM, no count and no matrix entry was opened."
        ),
        "wording_limit": pr.WORDING_LIMIT,
        "no_biosample_added": pb.SELECTION_RULE["no_biosample"],
        "stop_list_check": (
            "docs/ROADMAP.md section 5, 'The stop list of 2026-09-27, reinstated': NO biosample "
            "added -- only cell types the reader already covers -- no AlphaGenome request, no "
            "deletion request, nothing in area J, no node-model work, no BioForge calibration."
        ),
        "alphagenome_requests": 0,
        "money": "$0.00: no model request, no key, no data needing an application",
        "pinned_sha256": {
            "data/results/multrna_result.json": pr.sha256_of(pr.MULTRNA_RESULT),
            "data/results/magcal_result.json": pr.sha256_of(pr.MAGCAL_RESULT),
            "data/results/freshrna_result.json": pr.sha256_of(pr.FRESHRNA_RESULT),
            "data/results/freshrna_sources.json": pr.sha256_of(pr.FRESHRNA_SOURCES),
            "data/results/pseudobulk_inventory.json": pr.sha256_of(pr.PSEUDOBULK_INVENTORY),
            "data/results/scsweep_inventory.json": pr.sha256_of(pr.SCSWEEP_RESULT),
        },
    }
    record["seconds"] = round(time.time() - started, 2)
    record["result_manifest"] = result_manifest(
        mf.traced_reads(),
        {
            "sc_selection_version": pr.SC_SELECTION["version"],
            "form": pr.FORM,
            "minimum_cell_types": pr.MINIMUM["cell_types"],
            "minimum_genes": pr.MINIMUM["genes"],
            "held_out_population": len(population),
            "scale_constant": pr.SCALE,
        },
    )
    path = save_result("pbrun_gene_count", record)
    print(
        json.dumps(
            {
                "written": str(path),
                "answer": answer,
                "held_out_population": len(population),
                "qualifying_gene_count": len(shared),
                "per_cell": {
                    e["cell"]: {
                        "usable_on_both_sides": e["usable_on_both_sides"],
                        "modality": e["modality"],
                        "files": e["contributing_single_cell_files"],
                    }
                    for e in per_cell
                },
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
