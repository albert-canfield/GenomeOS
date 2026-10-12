#!/usr/bin/env python3
# SPDX-License-Identifier: AGPL-3.0-or-later
"""The registered run, executed ONCE: M2h's endpoint read twice, from bulk and from pseudobulk.

Reads `data/results/pbrun_registration.json` (`ccd5379`) and refuses on any disagreement with it.
The two arms differ in ONE thing, the measured side: the bulk arm's D comes from M2i's committed
ENCODE polyA plus RNA-seq files, the single-cell arm's from the pseudobulk the registration states
exactly. The prediction, the scale, the metric, the bootstrap and the control are identical between
them, and nothing is fitted anywhere.

The registered reading is written verbatim whatever it is, and the three outcomes -- agreement,
disagreement, inconclusive -- were fixed before any value here was opened.
"""

from __future__ import annotations

import io
import json
import sys
import tarfile
import time
from pathlib import Path
from typing import Any

import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import pbrun as pr  # noqa: E402
from genomeos.attribution import pseudobulk as pb  # noqa: E402
from genomeos.results import save_result  # noqa: E402
from scripts import pbrun_count as pc  # noqa: E402

ENTRY = "scripts/pbrun_run.py"
OWN_CODE = (
    "genomeos/attribution/pbrun.py",
    "scripts/pbrun_count.py",
    "scripts/pbrun_register.py",
    "scripts/pbrun_run.py",
    "tests/test_pbrun.py",
)


class HarnessDefect(SystemExit):
    """A pinned fact does not hold, so nothing is read and no reading is written."""


def cells_counted(source: dict[str, Any]) -> int:
    """How many cells or nuclei one contributing file supplies, read from the file itself."""
    path = ROOT / source["path"]
    if source["reader"] == "rsem_gene_id_column":
        return 1
    if source["reader"] == "starsolo_features_tsv":
        with tarfile.open(path, "r:gz") as tar:
            member = tar.extractfile(pr._one_member(tar.getnames(), "barcodes.tsv"))
            if member is None:
                raise HarnessDefect(f"{path}: barcodes.tsv is not readable")
            return sum(1 for line in io.TextIOWrapper(member, encoding="utf-8") if line.strip())
    if source["reader"] == "h5ad":
        with h5py.File(path, "r") as handle:
            return int(handle[source["counts_at"]].attrs["shape"][0])
    raise HarnessDefect(f"no cell counter for reader {source['reader']!r}")


def single_cell_counts(sources: list[dict[str, Any]]) -> tuple[dict[str, float], float]:
    """AGGREGATION step (1) for one cell type: the summed integer counts and the library size."""
    readers = {s["reader"] for s in sources}
    if len(readers) != 1:
        raise HarnessDefect(f"a cell type's files disagree on the reader: {readers}")
    reader = readers.pop()
    paths = [ROOT / s["path"] for s in sources]
    if reader == "rsem_gene_id_column":
        return pr.sum_counts_rsem(paths)
    if reader == "starsolo_features_tsv":
        return pr.sum_counts_starsolo(paths)
    if reader == "h5ad":
        if len(paths) != 1:
            raise HarnessDefect("an h5ad source must be a single file under SC_SELECTION")
        return pr.sum_counts_h5ad(paths[0], sources[0]["counts_at"], sources[0]["gene_axis"])
    raise HarnessDefect(f"no reader named {reader!r}")


def reading_of(arm: str, delta: float, ci: list[float], control: float, pairs: int, genes: int) -> str:
    """One arm's reading, in M2h's words under CARRIED['which_governs_m2h']."""
    if ci[0] > 0.0:
        verdict = (
            "better; on these genes the reader-gated run scaled by M2h's fitted constant +0.2146 "
            f"lowered the squared error of {arm} below the no-rule baseline's, so it is closer to "
            f"{arm} than the no-rule baseline"
        )
    elif ci[1] < 0.0:
        verdict = (
            f"worse; the scaled run's squared error is ABOVE the no-rule baseline's on {arm}, which "
            "is this lane's falsifier"
        )
    else:
        verdict = (
            "inconclusive; the interval straddles 0, and because the scale is a fixed constant there "
            "is no scale interval to split the straddle, so this measurement does not resolve "
            f"whether the scaled run is closer to or further from {arm} than the no-rule baseline"
        )
    return (
        f"{verdict} (Delta MSE {delta:+.6g} [98.33% CI {ci[0]:+.6g}, {ci[1]:+.6g}] over {pairs:,} "
        f"gene-cell pairs, {genes:,} genes; control p = {control:.6g}; the scale is the constant "
        "+0.2146 taken from ac5580b and NOTHING IS FITTED IN THIS LANE)"
    )


def result_manifest(read_now: list[str], parameters: dict[str, Any]) -> dict[str, Any]:
    inputs = [
        mf.input_entry(p, partition="M2h's HELD-OUT side, as the registration pins it")
        for p in sorted(set(read_now) | set(pc.UNTRACED_READS))
        if "/data/" in p.replace("\\", "/") or p.startswith("data/")
    ]
    return {
        "sources": [
            {
                "accession": (
                    "data/results/pbrun_registration.json and data/results/pbrun_gene_count.json "
                    "(this lane's committed registration and count), "
                    "data/results/multrna_result.json (7h's prediction rows), "
                    "data/results/magcal_result.json (M2h's constant and the margin), "
                    "data/results/freshrna_sources.json (the bulk files)"
                ),
                "version": "as committed in this checkout; pinned by sha256 in pinned_sha256",
            },
            {
                "accession": (
                    "ENCODE polyA plus RNA-seq gene quantifications of GM12878, K562, ovary and "
                    "testis (M2i's selection), and ENCODE scRNA-seq and snRNA-seq processed files "
                    "of GM12878 and ovary"
                ),
                "version": pb.ENCODE_LICENCE,
            },
            {
                "accession": "GEO GSM9224207 (GSE306708) GSM9224207_Fresh_K562.h5ad",
                "version": pc.NON_ENCODE_SOURCES["K562"]["licence"],
            },
            {
                "accession": (
                    "CELLxGENE Discover dataset ec6ef004-5f30-4a40-b282-980cca20a561, Tabula Sapiens - Testis"
                ),
                "version": pc.NON_ENCODE_SOURCES["testis"]["licence"],
            },
        ],
        "inputs": inputs,
        "input_count": len(inputs),
        "assembly": (
            "GRCh38 on every side; GM12878's single-cell tables are GENCODE V24 and everything else "
            "V29, and genes are matched on the unversioned Ensembl id or on the symbol"
        ),
        "coordinates": "n/a: no interval is read or written. The endpoint is a per-gene level.",
        "parameters": parameters,
        "exclusions": [
            "every gene outside the 8,051 the committed count record pins",
            "every cell type outside the registered four",
            "every normalised, decontaminated, imputed or scaled layer of either h5ad",
            "nothing else: no expression filter, no cell filter, no gene dropped for its value",
        ],
        "partitions": {
            "read": "M2h's HELD-OUT side of the multiplicative form, the registered 8,051 genes",
            "not_read": "M2h's TRAINING genes: the scale is the constant already fitted on them",
        },
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE),
        "entry": ENTRY,
        "own_code": list(OWN_CODE),
    }


def main() -> int:
    started = time.time()
    registration = json.loads(pr.REGISTRATION_PATH.read_text())
    count = json.loads(pr.GENE_COUNT_PATH.read_text())
    harness: dict[str, Any] = {}

    # (1) every pinned input is the file the registration pinned
    for path, digest in registration["pinned_sha256"].items():
        if pr.sha256_of(ROOT / path) != digest:
            raise HarnessDefect(f"HARNESS_DEFECT: {path} is not the file the registration pinned")
    harness["pinned_inputs_match"] = True

    # (2) the population is the registered one and clears both floors
    population: list[str] = count["qualifying_genes"]
    if registration["population_genes"] != len(population):
        raise HarnessDefect("HARNESS_DEFECT: the registration and the count disagree on the genes")
    if len(population) < max(pr.MINIMUM["genes"], pr.CARRIED["min_genes"]):
        raise HarnessDefect("HARNESS_DEFECT: the population is below a registered floor")
    harness["population_genes"] = len(population)

    # (3) the cell types are the registered four
    cells: list[str] = registration["cells"]
    if len(cells) != pr.MINIMUM["cell_types"] or cells != count["qualifying_cells"]:
        raise HarnessDefect("HARNESS_DEFECT: the cell types are not the registered four")
    harness["cell_types"] = len(cells)

    # (4) SCALE is M2h's committed constant
    committed = json.loads(pr.MAGCAL_RESULT.read_text())["per_form"][pr.FORM]
    if committed["scale"] != pr.SCALE or registration["scale"] != pr.SCALE:
        raise HarnessDefect("HARNESS_DEFECT: SCALE is not M2h's committed constant")
    harness["scale_is_m2hs_constant"] = True
    margin = pr.margin()
    if registration["margin"] != margin:
        raise HarnessDefect("HARNESS_DEFECT: the margin is not the registered one")

    # ---- the prediction and the gene identifiers: 7h's committed rows ----------------------------
    table = json.loads(pr.MULTRNA_RESULT.read_text())
    columns = table["genes_table"]["columns"]
    i_symbol, i_gene_id, i_p = (
        columns.index("symbol"),
        columns.index("gene_id"),
        columns.index("P"),
    )
    thirteen = table["cells"]
    column_of = [thirteen.index(c) for c in cells]
    rows = {str(r[i_symbol]): r for r in table["genes_table"]["rows"]}
    absent = [g for g in population if g not in rows]
    if absent:
        raise HarnessDefect(f"HARNESS_DEFECT: {len(absent)} registered genes are not in 7h's table")
    gene_id_of = {g: pr.unversioned(str(rows[g][i_gene_id])) for g in population}
    p = pr.centred(np.array([[rows[g][i_p][j] for j in column_of] for g in population], dtype=np.float64))

    # ---- the bulk arm: M2i's committed files, M2h's endpoint -------------------------------------
    bulk_sources = pc.bulk_sources(cells)
    bulk_y = np.zeros((len(population), len(cells)), dtype=np.float64)
    for j, cell in enumerate(cells):
        tpm = pr.bulk_tpm([ROOT / f["path"] for f in bulk_sources[cell]])
        gone = [g for g in population if gene_id_of[g] not in tpm]
        if gone:
            raise HarnessDefect(
                f"HARNESS_DEFECT: {len(gone)} registered genes do not resolve in {cell}'s bulk "
                "tables, which the committed count says they do"
            )
        bulk_y[:, j] = pr.log2p1(np.array([tpm[gene_id_of[g]] for g in population], dtype=np.float64))
    d_bulk = pr.centred(bulk_y)

    # ---- the single-cell arm: AGGREGATION, exactly as registered ---------------------------------
    single_sources = pc.single_cell_sources(cells)
    sc_y = np.zeros((len(population), len(cells)), dtype=np.float64)
    per_cell_source: list[dict[str, Any]] = []
    for j, cell in enumerate(cells):
        sources = single_sources[cell]
        counts, library = single_cell_counts(sources)
        identifier = sources[0]["identifier"]
        key_of = (lambda g: gene_id_of[g]) if identifier == pr.ENSEMBL else (lambda g: g)
        gone = [g for g in population if key_of(g) not in counts]
        if gone:
            raise HarnessDefect(
                f"HARNESS_DEFECT: {len(gone)} registered genes have no row in {cell}'s single-cell "
                "matrix, which the committed count says they do"
            )
        summed = np.array([counts[key_of(g)] for g in population], dtype=np.float64)
        sc_y[:, j] = pr.log2p1(1e6 * summed / library)
        per_cell_source.append(
            {
                "cell": cell,
                "modality": sorted({s["modality"] for s in sources}),
                "files": len(sources),
                "cells_or_nuclei": sum(cells_counted(s) for s in sources),
                "library_size_counts": library,
                "identifier_matched_on": identifier,
                "accessions": [s["accession"] for s in sources],
            }
        )
    d_sc = pr.centred(sc_y)

    # (5), (6) both arms are centred the same way and every gene resolved on both sides
    for name, d in (("bulk", d_bulk), ("single_cell", d_sc)):
        if float(np.max(np.abs(d.mean(axis=1)))) > pr.TOLERANCE:
            raise HarnessDefect(f"HARNESS_DEFECT: the {name} arm is not centred within every gene")
    harness["both_arms_centred_within_every_gene"] = True
    harness["every_registered_gene_resolved_on_both_sides"] = True

    # ---- the endpoint, at the fixed constant -----------------------------------------------------
    delta_sc = pr.delta_mse(d_sc, p, pr.SCALE)
    delta_bulk = pr.delta_mse(d_bulk, p, pr.SCALE)
    difference = delta_sc - delta_bulk

    # ---- the gene bootstrap, PAIRED on the same resample -----------------------------------------
    rng = np.random.default_rng(pr.BOOTSTRAP["seed"])
    n = len(population)
    boot_sc = np.empty(pr.BOOTSTRAP["resamples"], dtype=np.float64)
    boot_bulk = np.empty(pr.BOOTSTRAP["resamples"], dtype=np.float64)
    for b in range(pr.BOOTSTRAP["resamples"]):
        draw = rng.integers(0, n, size=n)
        boot_sc[b] = pr.delta_mse(d_sc[draw], p[draw], pr.SCALE)
        boot_bulk[b] = pr.delta_mse(d_bulk[draw], p[draw], pr.SCALE)
    ci_sc = pr.percentile_interval(boot_sc)
    ci_bulk = pr.percentile_interval(boot_bulk)
    ci_difference = pr.percentile_interval(boot_sc - boot_bulk)

    # ---- the control, on both arms ---------------------------------------------------------------
    permuted_sc = np.empty(pr.PERMUTATIONS["count"], dtype=np.float64)
    permuted_bulk = np.empty(pr.PERMUTATIONS["count"], dtype=np.float64)
    for i, permuted in enumerate(
        pr.permute_within_gene(p, pr.PERMUTATIONS["seed"], pr.PERMUTATIONS["count"])
    ):
        permuted_sc[i] = pr.delta_mse(d_sc, permuted, pr.SCALE)
        permuted_bulk[i] = pr.delta_mse(d_bulk, permuted, pr.SCALE)
    control_sc = pr.control_p(delta_sc, permuted_sc)
    control_bulk = pr.control_p(delta_bulk, permuted_bulk)

    # ---- the decision, exactly as registered -----------------------------------------------------
    sc_reads = pr.arm_reads(ci_sc)
    inside = -margin <= ci_difference[0] and ci_difference[1] <= margin
    wholly_outside = ci_difference[0] > margin or ci_difference[1] < -margin
    outcome = pr.decide(ci_sc, ci_difference, margin)

    pairs = n * len(cells)
    reading = (
        "pseudobulked single-cell and single-nucleus RNA against bulk RNA-seq on M2h's endpoint, "
        f"over the same {n:,} genes and the same four cell types at the fixed constant +0.2146. THE "
        "SINGLE-CELL ARM: "
        + reading_of("measured pseudobulk log2(CPM+1)", delta_sc, ci_sc, control_sc, pairs, n)
        + ". THE BULK ARM, read on the same genes, the same cell types and the same prediction: "
        + reading_of("measured bulk log2(TPM+1)", delta_bulk, ci_bulk, control_bulk, pairs, n)
        + f". THE PAIRED DIFFERENCE, single-cell minus bulk, on the same gene resamples: "
        f"{difference:+.6g} [98.33% CI {ci_difference[0]:+.6g}, {ci_difference[1]:+.6g}] against a "
        f"registered equivalence margin of +-{margin:.6g}, which is half of M2h's own committed "
        f"held-out gain. THE REGISTERED OUTCOME IS {outcome}. "
        + {
            "AGREEMENT": (
                "Agreement licenses THE MODALITY SWAP FOR THIS ENDPOINT ONLY and makes the bile-duct "
                "out-of-domain test interpretable, so the supervisor's recommendation to Albert can "
                "change. It is NEVER 'single-nucleus RNA equals bulk RNA', and it is never a "
                "statement about any other endpoint, any other gene set or any cell type not "
                "counted here."
            ),
            "DISAGREEMENT": (
                "Disagreement means pseudobulked single-cell RNA does not reproduce bulk on this "
                "endpoint: the modality swap stays UNLICENSED, the bile-duct out-of-domain test "
                "stays uninterpretable, and Albert's biosample question CLOSES WITHOUT A LIFT."
            ),
            "INCONCLUSIVE": (
                "Inconclusive LICENSES NOTHING: the bile-duct out-of-domain test stays "
                "uninterpretable and no recommendation changes. It is not a weak agreement and the "
                "registration says so in those words."
            ),
        }[outcome]
        + " CARRIED WITH EVERY READING: the bulk side is NOT a fresh measurement. It is RNA already "
        "read by items 7, 7d, 7h, M2h and M2i -- precisely, M2i's committed selection of ENCODE "
        "polyA plus RNA-seq of these same four biosamples under this same endpoint -- so the "
        "single-cell side is the only side of this comparison the project has never read. ONLY ONE "
        "OF THE FOUR CELL TYPES IS SINGLE-NUCLEUS (the ovary); GM12878, K562 and testis are "
        "single-CELL, and the bile-duct source is single-NUCLEUS, so the nucleus-specific part of "
        "anything this reading licenses rests on one of the four cell types."
    )

    per_cell_readings = []
    for j, cell in enumerate(cells):
        per_cell_readings.append(
            {
                "cell": cell,
                "single_cell_delta_mse": pr.delta_mse(d_sc[:, [j]], p[:, [j]], pr.SCALE),
                "bulk_delta_mse": pr.delta_mse(d_bulk[:, [j]], p[:, [j]], pr.SCALE),
                "correlation_of_the_two_arms_D": float(np.corrcoef(d_sc[:, j], d_bulk[:, j])[0, 1]),
            }
        )

    record: dict[str, Any] = {
        "result": "pbrun_result",
        "date": time.strftime("%Y-%m-%d"),
        "lane": "lane-pbrun",
        "registration": "data/results/pbrun_registration.json",
        "registration_sha256": pr.sha256_of(pr.REGISTRATION_PATH),
        "gene_count_record": "data/results/pbrun_gene_count.json",
        "question": registration["question"],
        "what_it_gates": registration["what_it_gates"],
        "outcome": outcome,
        "reading": reading,
        "cells": cells,
        "modality_per_cell": registration["modality_per_cell"],
        "genes": n,
        "pairs": pairs,
        "scale": pr.SCALE,
        "nothing_is_fitted": registration["nothing_is_fitted"],
        "single_cell_arm": {
            "delta_mse": delta_sc,
            "ci": ci_sc,
            "control_p": control_sc,
            "reads": sc_reads,
            "delta_mse_at_scale_1": pr.delta_mse(d_sc, p, 1.0),
        },
        "bulk_arm": {
            "delta_mse": delta_bulk,
            "ci": ci_bulk,
            "control_p": control_bulk,
            "reads": pr.arm_reads(ci_bulk),
            "delta_mse_at_scale_1": pr.delta_mse(d_bulk, p, 1.0),
        },
        "paired_difference": {
            "single_cell_minus_bulk": difference,
            "ci": ci_difference,
            "margin": margin,
            "margin_rule": pr.MARGIN_RULE,
            "interval_inside_the_margin": bool(inside),
            "interval_wholly_outside_the_margin": bool(wholly_outside),
        },
        "which_governs": registration["which_governs"],
        "decision": registration["decision"],
        "falsifier": registration["falsifier"],
        "harness_check": harness,
        "descriptives_not_verdicts": {
            "per_cell": per_cell_readings,
            "correlation_of_the_two_arms_D_over_all_pairs": float(
                np.corrcoef(d_sc.ravel(), d_bulk.ravel())[0, 1]
            ),
            "sources_per_cell": per_cell_source,
            "control_delta_mse_median_single_cell": float(np.median(permuted_sc)),
            "control_delta_mse_median_bulk": float(np.median(permuted_bulk)),
        },
        "aggregation": registration["aggregation"],
        "deviation_from_m2h": registration["deviation_from_m2h"],
        "prior_exposure": registration["prior_exposure"],
        "wording_limit": pr.WORDING_LIMIT,
        "limitations": registration["limitations"],
        "no_biosample_added": pb.SELECTION_RULE["no_biosample"],
        "stop_list_check": registration["stop_list_check"],
        "alphagenome_requests": 0,
        "money": "$0.00: no model request, no key, no data needing an application",
        "pinned_sha256": registration["pinned_sha256"],
    }
    record["seconds"] = round(time.time() - started, 2)
    record["result_manifest"] = result_manifest(
        mf.traced_reads(),
        {
            "scale": pr.SCALE,
            "margin": margin,
            "resamples": pr.BOOTSTRAP["resamples"],
            "bootstrap_seed": pr.BOOTSTRAP["seed"],
            "permutations": pr.PERMUTATIONS,
            "genes": n,
            "cell_count": len(cells),
            "sc_selection_version": pr.SC_SELECTION["version"],
            "tolerance": pr.TOLERANCE,
        },
    )
    path = save_result("pbrun_result", record)
    print(
        json.dumps(
            {
                "written": str(path),
                "outcome": outcome,
                "single_cell": [delta_sc, ci_sc, control_sc],
                "bulk": [delta_bulk, ci_bulk, control_bulk],
                "difference": [difference, ci_difference, margin],
                "reading": reading,
            },
            indent=1,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
