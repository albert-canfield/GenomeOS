# SPDX-License-Identifier: AGPL-3.0-or-later
"""Simpler genomes, the fly: inventory, census and registration, then the registered run.

    uv run --no-project --with pyreadr --with pandas --with numpy python scripts/simpler_fly_fetch.py
    uv run python scripts/simpler_fly_run.py inventory   # data/results/simpler_fly_inventory.json
    uv run python scripts/simpler_fly_run.py register    # census + registration, no RNA value read
    uv run python scripts/simpler_fly_run.py score       # refuses unless the registration is committed

No model request is sent and no key is read. docs/SIMPLER-GENOMES.md section 8 holds the reasoning.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import simpler  # noqa: E402
from genomeos.attribution import simpler_fly as sf  # noqa: E402
from genomeos.results import save_result  # noqa: E402

INVENTORY = "simpler_fly_inventory"
REGISTRATION = "simpler_fly_registration"
RESULT = "simpler_fly_result"
REGISTRATION_FILE = "data/results/simpler_fly_registration.json"
ENTRY = "scripts/simpler_fly_run.py"
OWN_CODE = (
    "genomeos/attribution/simpler.py",
    "genomeos/attribution/simpler_fly.py",
    "scripts/simpler_fly_fetch.py",
    "scripts/simpler_fly_run.py",
    "tests/test_simpler_fly.py",
)

#: Licences as the sources state them, checked on 2026-10-10. None is inferred.
GEO_TERMS = "GEO public record; no licence stated by the submitters"
LICENCES = {
    "SupplementaryTables.xlsx": "not stated on the authors' site (DEAP); the article is Science 377:eabn5800",
    "peaks.bed": "not stated on the authors' site (DEAP)",
    "rna_meta.rds": "not stated on the authors' site (DEAP)",
    "GSE190147_scirnaseq_gene_matrix.rows.txt.gz": GEO_TERMS,
    "GSE190147_scirnaseq_gene_matrix.columns.txt.gz": GEO_TERMS,
    "GSE190147_scirnaseq_gene_matrix.mtx.gz": GEO_TERMS,
    "Drosophila_melanogaster.BDGP6.32.104.gff3.gz": (
        "Ensembl states no restrictions on use of its data; the gene models are FlyBase's (CC BY 4.0)"
    ),
}
#: Bytes moved while scouting, outside the fetch record: the DEAP page, a GEO range read of the
#: matrix's names and header, and the GFF3 of Ensembl 99, 104 and 113 to find the release whose
#: protein-coding gene names match the RNA rows (104: all 13,968; 99: 13,837; 113: 13,594).
SCOUTING_BYTES = 862_284 + 88_224 + 233_003 + 5_056_879 + 5_713_526 + 5_814_950

NOT_USED = [
    {
        "source": (
            "Calderon 2022 ATAC peak-by-cell matrices (DEAP ATAC/revision/peak_matrices; GEO GSE190130)"
        ),
        "bytes": "119 MB to 2.1 GB per window and experiment, about 15 GB in all",
        "why_not": "over the 1.5 GB budget; Table S4 carries the source's own per-cluster call",
    },
    {
        "source": "Calderon 2022 ATAC bigWigs per window and cluster (DEAP ATAC/revision/bigwigs)",
        "bytes": "2 MB to 577 MB each, more than 100 files",
        "why_not": "a per-element read of every cluster track would move most of each file",
    },
    {
        "source": "Measured insulator or TAD boundary sets (dCTCF, BEAF-32, CP190; embryo Hi-C)",
        "bytes": "not checked",
        "why_not": (
            "the node decision (docs/SIMPLER-GENOMES.md section 8): this run replicates the worm's rule"
        ),
    },
]


def _manifest(
    inputs: list[tuple[str, str]], parameters: dict, exclusions: list[str], partitions: Any
) -> dict:
    return {
        "sources": [
            {"accession": sf.SOURCE, "version": "authors' site files as served on 2026-10-10; GEO GSE190147"},
            {"accession": sf.GENE_MODELS, "version": "Ensembl release 104"},
        ],
        "inputs": [mf.input_entry(path, partition=part) for path, part in inputs],
        "assembly": (
            "dm6 (BDGP6.32): the peaks, the source's RNA alignment and the gene models; nothing is lifted"
        ),
        "coordinates": {
            "base": 0,
            "interval": "half-open",
            "note": (
                "the peaks are BED (0-based, half-open) and the GFF3 is 1-based closed; the element "
                "midpoint (start + end) // 2 moves by at most 1 bp between the two conventions"
            ),
        },
        "parameters": parameters,
        "exclusions": exclusions,
        "partitions": partitions,
        "code_cleanliness": mf.code_cleanliness(ENTRY, OWN_CODE, ROOT),
    }


PARAMETERS = {
    "tissues": list(sf.TISSUES),
    "promoter_bp": sf.PROMOTER_BP,
    "distal_bp": sf.DISTAL_BP,
    "detection_cpm": sf.DETECTION_CPM,
    "marker_p_adj": sf.MARKER_P_ADJ,
    "distance_bins_bp": list(sf.DISTANCE_BINS),
    "control_draws": sf.CONTROL_DRAWS,
    "bootstrap": sf.BOOTSTRAP,
    "cluster_bp": sf.CLUSTER_BP,
    "seed": sf.SEED,
    "min_discordant": sf.MIN_DISCORDANT,
    "alpha": sf.ALPHA,
}
EXCLUSIONS = [
    "Table S4 rows with avg_logFC <= 0 (less accessible in the cluster)",
    "peaks with any positive marker call in a cluster whose label is excluded: "
    + ", ".join(sorted(sf.ATAC_EXCLUDED)),
    "peaks called in clusters of two or more tissues",
    "peaks neither distal (> 1,000 bp from every transcript 5' end) nor promoter (<= 500 bp from a "
    "protein-coding gene's 5' end)",
    "RNA cells whose label is excluded: " + ", ".join(sorted(sf.RNA_EXCLUDED)),
    "genes that are not protein_coding, and chromosomes other than 2L, 2R, 3L, 3R, 4, X, Y",
    "supplementary sheets never read: " + ", ".join(sf.NEVER_READ),
]


def _record() -> dict:
    return json.loads(sf.FETCH_RECORD.read_text())


def inventory() -> None:
    mf.trace_begin()
    rec = _record()
    fetched = []
    for name, f in rec["files"].items():
        fetched.append({**f, "name": name, "licence": LICENCES[name]})
    d = rec["distillation"]
    out = {
        "question": "the Calderon 2022 tables for the fly test: what was fetched, its size, hash and licence",
        "source": sf.SOURCE,
        "fetched": fetched,
        "distilled": {
            k: d[k] for k in ("what", "file", "bytes", "sha256", "matrix_shape", "matrix_nonzero_entries")
        }
        | {
            "cells_in_matrix": d["cells_in_matrix"],
            "cells_with_an_annotation": d["cells_with_an_annotation"],
            "annotations": d["annotations"],
            "cells_per_annotation": d["cells_per_annotation"],
            "method_follows": d["method_follows"],
        },
        "fetched_bytes": rec["download_bytes"],
        "scouting_bytes": SCOUTING_BYTES,
        "lane_download_bytes": rec["download_bytes"] + SCOUTING_BYTES,
        "kept_bytes": sum(f["bytes"] for f in rec["files"].values() if f["kept"]) + d["bytes"],
        "budget_bytes": 1_500_000_000,
        "not_used": NOT_USED,
        "money": "none; no model request",
    }
    out["result_manifest"] = _manifest(
        [(str(sf.FETCH_RECORD), "the fetch record written by scripts/simpler_fly_fetch.py")],
        {"download_budget_bytes": 1_500_000_000},
        ["the RNA matrix (697,773,156 bytes) was summed by annotation and deleted"],
        "n/a: an inventory has no evaluation partition",
    )
    print(save_result(INVENTORY, out))


def _rows() -> tuple[dict, list[dict[str, Any]], dict[str, Any], dict[str, str]]:
    chroms, tss, name_to_id = sf.read_annotation(sf.GFF3)
    calls, funnel = sf.tissue_calls(sf.read_markers(sf.SUPP), sf.read_cluster_labels(sf.SUPP))
    els, counts = sf.elements(calls, sf.read_peaks(sf.PEAKS), chroms, tss)
    return chroms, simpler.attribute(els, chroms), {"marker_calls": funnel, "classes": counts}, name_to_id


def _cells_per_tissue() -> dict[str, int]:
    per = dict.fromkeys(sf.TISSUES, 0)
    for lab, n in _record()["distillation"]["cells_per_annotation"].items():
        if lab in sf.RNA_TISSUE:
            per[sf.RNA_TISSUE[lab]] += n
    return per


def register() -> None:
    mf.trace_begin()
    chroms, rows, funnel, _ = _rows()
    cen = sf.census(rows)
    rec = _record()
    out = {
        "written_before_any_rna_value_was_read": True,
        "registered_inputs": {
            "supplementary_tables": {"path": str(sf.SUPP), "sha256": simpler.sha256(sf.SUPP)},
            "peaks": {"path": str(sf.PEAKS), "sha256": simpler.sha256(sf.PEAKS)},
            "gff3": {"path": str(sf.GFF3), "sha256": simpler.sha256(sf.GFF3)},
            "rna_by_annotation": {
                "path": str(sf.RNA),
                "sha256": rec["distillation"]["sha256"],
                "note": "distilled by scripts/simpler_fly_fetch.py; its values are first read by score",
            },
            "protein_coding_genes": sum(len(c.genes) for c in chroms.values()),
        },
        "rna_cells_per_tissue": _cells_per_tissue(),
        "registration": sf.REGISTRATION,
        "population_funnel": funnel,
        "census": cen,
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (str(sf.SUPP), "Table S2 and Table S4 only (ATAC cluster labels and marker calls)"),
            (str(sf.PEAKS), "the peak set"),
            (str(sf.GFF3), "gene and transcript positions and strands"),
            (str(sf.FETCH_RECORD), "cells per RNA annotation (labels only) and the distilled file's sha256"),
        ],
        PARAMETERS,
        EXCLUSIONS,
        {"population": "fixed here, before any RNA value was read"},
    )
    print(save_result(REGISTRATION, out))
    print(json.dumps({"funnel": funnel, "census": cen}, indent=1))


def _committed(rel: str) -> bool:
    logged = subprocess.run(
        ["git", "log", "-1", "--format=%H", "--", rel], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    clean = subprocess.run(["git", "diff", "--quiet", "HEAD", "--", rel], cwd=ROOT).returncode == 0
    return bool(logged) and clean


def score() -> None:
    for rel in (REGISTRATION_FILE, *OWN_CODE):
        if not _committed(rel):
            sys.exit(f"refused: {rel} must be committed and unchanged before any RNA value is read")
    mf.trace_begin()
    reg = json.loads(Path(REGISTRATION_FILE).read_text())
    if reg["registration"] != json.loads(json.dumps(sf.REGISTRATION)):
        sys.exit("refused: simpler_fly.REGISTRATION differs from the committed registration")
    for key in ("supplementary_tables", "peaks", "gff3", "rna_by_annotation"):
        entry = reg["registered_inputs"][key]
        if simpler.sha256(Path(entry["path"])) != entry["sha256"]:
            sys.exit(f"refused: {key} differs from the registered input")
    chroms, rows, funnel, name_to_id = _rows()
    if sf.census(rows) != reg["census"] or funnel != reg["population_funnel"]:
        sys.exit("refused: the census does not recompute to the registered one")
    expr, rna_summary = sf.read_expression(sf.RNA, name_to_id)
    res = sf.score(chroms, rows, expr)
    reg_sha = subprocess.run(
        ["git", "log", "-1", "--format=%h", "--", REGISTRATION_FILE], cwd=ROOT, capture_output=True, text=True
    ).stdout.strip()
    out = {
        "registration_commit": reg_sha,
        "registered_inputs": reg["registered_inputs"],
        "rna_summary": rna_summary,
        **res,
        "reading_is_the_registered_wording": (
            "the reading is REGISTRATION['outcome_wording'] with its numbers filled in "
            "and no other word changed"
        ),
        "alphagenome_requests": 0,
        "money": "none",
    }
    out["result_manifest"] = _manifest(
        [
            (REGISTRATION_FILE, "the committed registration"),
            (str(sf.SUPP), "Table S2 and Table S4 only"),
            (str(sf.PEAKS), "the peak set"),
            (str(sf.GFF3), "gene and transcript positions and strands"),
            (str(sf.RNA), "the endpoint: UMI sums by annotation"),
        ],
        PARAMETERS,
        EXCLUSIONS,
        {"population": "pinned by the committed registration's census", "positive_control": "promoter"},
    )
    print(save_result(RESULT, out))
    print(out["outcome"])
    print(out.get("reading", ""))


if __name__ == "__main__":
    os.chdir(ROOT)
    {"inventory": inventory, "register": register, "score": score}[sys.argv[1]]()
