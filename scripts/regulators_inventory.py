# SPDX-License-Identifier: AGPL-3.0-or-later
"""M1b inventory: what evidence the project holds or can fetch for WHICH factor acts at a compiled element.

    uv run python scripts/regulators_inventory.py [--no-save] [--refresh]

Metadata only. It reads the identity of each cached source (size, sha256), the factor and cell names
of the ENCODE peak clusters, and, from the ENCODE portal's JSON API, the identities of the published
transcription-factor knockdown RNA-seq experiments in the test cell types with their control
experiments and file listings (accessions, md5, sizes). It opens no expression value: no
quantification file is downloaded here. The experiment metadata is cached, git-ignored, at
data/knowledge/regulators/knockdown/experiments.json so the registration pins what the inventory saw.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomeos import manifest as mf  # noqa: E402
from genomeos.attribution import regulators as rg  # noqa: E402
from genomeos.results import save_result  # noqa: E402

RESULT = "regulators_inventory"
PORTAL = "https://www.encodeproject.org"
META = rg.KD_DIR / "experiments.json"
OWN_CODE = (
    "genomeos/attribution/regulators.py",
    "scripts/regulators_inventory.py",
    "scripts/regulators_census.py",
    "scripts/regulators_register.py",
    "scripts/regulators_fetch.py",
    "scripts/regulators_run.py",
    "tests/test_regulators.py",
)


def get_json(path: str, params: dict | None = None) -> dict:
    url = PORTAL + path + ("?" + urllib.parse.urlencode(params, doseq=True) if params else "")
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return json.load(r)
        except urllib.error.HTTPError as e:
            if e.code == 404 or attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))
        except Exception:
            if attempt == 3:
                raise
            time.sleep(5 * (attempt + 1))
    raise RuntimeError("unreachable")


def search(params: dict) -> list[dict]:
    """The portal answers a search with no hit by HTTP 404, which is an empty list, not an error."""
    base = {"type": "Experiment", "status": "released", "format": "json", "limit": "all"}
    try:
        return get_json("/search/", {**base, **params}).get("@graph", [])
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return []
        raise


def count(params: dict) -> int:
    base = {"type": "Experiment", "status": "released", "format": "json", "limit": "0"}
    return int(get_json("/search/", {**base, **params}).get("total", 0))


def fetch_metadata() -> dict:
    """Every TF knockdown experiment of the registered assays in the test cells, and its controls."""
    exps: dict[str, dict] = {}
    listed = []
    for cell in rg.TEST_CELLS:
        for assay in rg.KD_ASSAYS:
            hits = search(
                {
                    "assay_title": assay,
                    "biosample_ontology.term_name": cell,
                    "target.investigated_as": rg.TARGET_CLASS,
                }
            )
            for h in hits:
                listed.append({"accession": h["accession"], "cell": cell, "assay": assay})
    for row in listed:
        if row["accession"] not in exps:
            exp = get_json(f"/experiments/{row['accession']}/", {"format": "json"})
            exps[row["accession"]] = rg.slim_experiment(exp)
            time.sleep(0.2)
    for acc in sorted({c for e in list(exps.values()) for c in e["possible_controls"]}):
        if acc not in exps:
            exps[acc] = rg.slim_experiment(get_json(f"/experiments/{acc}/", {"format": "json"}))
            time.sleep(0.2)
    return {"fetched": time.strftime("%Y-%m-%d"), "listed": listed, "experiments": exps}


def files_of(exps: dict, acc: str) -> dict | None:
    e = exps[acc]
    return rg.choose_files(e, [exps[c] for c in e["possible_controls"] if c in exps])


def main(argv: list[str] | None = None) -> int:
    from genomeos.genome import motifs as mo

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    ap.add_argument("--refresh", action="store_true")
    args = ap.parse_args(argv)
    t0 = time.time()
    if args.refresh or not META.exists():
        META.parent.mkdir(parents=True, exist_ok=True)
        META.write_text(json.dumps(fetch_metadata(), indent=1, sort_keys=True))
    meta = json.loads(META.read_text())
    cl = rg.load_clusters()
    chip = cl.factors_in()
    per_cell_chip = {c: len(chip.get(c, ())) for c in rg.TEST_CELLS}
    exps = meta["experiments"]
    kd_rows = []
    for row in meta["listed"]:
        e = exps[row["accession"]]
        files = files_of(exps, row["accession"])
        kd_rows.append(
            {
                "accession": e["accession"],
                "cell": row["cell"],
                "assay": row["assay"],
                "target": e["target"],
                "lab": e["lab"],
                "controls": e["possible_controls"],
                "target_has_chip_in_cell": e["target"] in chip.get(rg.context_label(row["cell"]), set()),
                "files_annotation": files["annotation"] if files else None,
                "knockdown_files": len(files["knockdown"]) if files else 0,
                "control_files": len(files["control"]) if files else 0,
            }
        )
    usable = [
        r for r in kd_rows if r["target_has_chip_in_cell"] and r["knockdown_files"] and r["control_files"]
    ]
    units = sorted({(r["target"], r["cell"]) for r in usable})
    to_fetch: dict[str, int] = {}
    for r in usable:
        f = files_of(exps, r["accession"])
        for x in f["knockdown"] + f["control"]:
            to_fetch[x["accession"]] = x["file_size"] or 0
    encode4_chip = {
        c: count({"assay_title": "TF ChIP-seq", "biosample_ontology.term_name": c}) for c in rg.TEST_CELLS
    }
    out = {
        "lane": rg.LANE,
        "question": "what evidence the project holds or can fetch publicly for WHICH factor acts at each "
        "compiled element, each with its evidence kind; metadata only, no outcome read",
        "sources": {
            "measured_occupancy": {
                **rg.CLUSTER_SOURCE,
                "local_path": str(rg.CLUSTERS),
                "present": rg.CLUSTERS.exists(),
                "factors": len(cl.factors),
                "cell_labels": len(chip),
                "factors_per_test_cell": per_cell_chip,
                "committed": False,
            },
            "predicted_motifs": {
                "name": "JASPAR 2026 CORE vertebrates, non-redundant (area J's operator source)",
                "url": mo.JASPAR_URL,
                "local_path": str(mo.JASPAR_PATH),
                "licence": "CC BY 4.0",
                "threshold": f"relative score {mo.RELATIVE_SCORE} (genomeos/genome/motifs.py, area J)",
                "kind": "predicted: a motif is not occupancy",
                "identity": mf.input_entry(mo.JASPAR_PATH, partition=None),
            },
            "direction": {**rg.GO_SOURCE, "identity": mf.input_entry(rg.GOA, partition=None)},
            "knockdown_rna_seq": {
                "name": "ENCODE transcription-factor knockdown RNA-seq, with each experiment's "
                "non-targeting control experiments",
                "portal": PORTAL,
                "assays": list(rg.KD_ASSAYS),
                "cells": list(rg.TEST_CELLS),
                "target_class": rg.TARGET_CLASS,
                "licence": "ENCODE data, no restrictions on use (ENCODE data use policy)",
                "metadata_cache": str(META),
                "metadata_fetched": meta["fetched"],
                "experiments_listed": len(kd_rows),
                "by_cell_and_assay": dict(Counter(f"{r['cell']} {r['assay']}" for r in kd_rows)),
                "with_target_chip_in_same_cell_and_files": len(usable),
                "factor_cell_units_before_eligibility": len(units),
                "units": [f"{f} {c}" for f, c in units],
                "files_to_fetch": len(to_fetch),
                "bytes_to_fetch": sum(to_fetch.values()),
                "project_exposure": "none: no file, script, result or document of this project read these "
                "experiments before this inventory (searched for the assay names on 2026-10-10)",
                "experiments": kd_rows,
            },
            "larger_occupancy_not_fetched": {
                "name": "ENCODE4 TF ChIP-seq experiments per test cell (portal search count, released)",
                "counts": encode4_chip,
                "why_not_used": "the UCSC clusters are one 125 MB file with cell types per peak; the "
                "ENCODE4 set is several hundred files per cell type. It would add factors, not a new kind",
            },
            "held_not_used": {
                "name": "Replogle et al. 2022 K562 genome-wide Perturb-seq pseudobulk (cached by N1)",
                "local_path": "data/cache/n1/K562_gwps_raw_bulk_01.h5ad",
                "why_not_used": "reading its perturbed rows is N1b step 2, which ROADMAP section 5 item 16 "
                "lists as the owner's decision; this lane reads none of it",
            },
        },
        "what_no_source_gives": [
            "the direction of a factor's effect at a given element: only factor-level curated direction "
            "(GO) exists here, and it is not element-level",
            "a strength for a factor -> gene rule: no source here measures how much a factor moves a gene "
            "through one element",
            "occupancy in the compiled contexts ENCODE did not profile for factors",
        ],
        "alphagenome_requests": 0,
        "money": "none: no request is sent and no key is read",
        "seconds": round(time.time() - t0, 1),
    }
    man = {
        "sources": [
            {"accession": rg.CLUSTER_SOURCE["url"], "version": rg.CLUSTER_SOURCE["last_modified"]},
            {"accession": PORTAL + " JSON API", "version": "metadata fetched " + meta["fetched"]},
        ],
        "inputs": [
            mf.input_entry(rg.CLUSTERS, partition=None),
            mf.input_entry(META, partition=None),
            mf.input_entry(mo.JASPAR_PATH, partition=None),
            mf.input_entry(rg.GOA, partition=None),
        ],
        "input_count": 4,
        "assembly": "GRCh38",
        "coordinates": "n/a: identities and counts",
        "parameters": {"test_cells": list(rg.TEST_CELLS), "kd_assays": list(rg.KD_ASSAYS)},
        "exclusions": [],
        "partitions": "n/a: identities and counts, no evaluation partition",
        "code_cleanliness": mf.code_cleanliness(__file__, OWN_CODE),
    }
    if not args.no_save:
        print("saved", save_result(RESULT, out, manifest=man))
    ko = out["sources"]["knockdown_rna_seq"]
    keys = (
        "experiments_listed",
        "by_cell_and_assay",
        "with_target_chip_in_same_cell_and_files",
        "factor_cell_units_before_eligibility",
        "units",
        "files_to_fetch",
        "bytes_to_fetch",
    )
    print(json.dumps({k: ko[k] for k in keys}, indent=1))
    print(per_cell_chip, encode4_chip)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
