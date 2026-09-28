# SPDX-License-Identifier: AGPL-3.0-or-later
"""Census of the `confidence` field on predicted enhancer-to-gene links, before review R4 retires it.

`predict/enhancer_target.py` wrote `confidence = round(min(0.7, |log2 fold change|), 3)` on every
predicted link (`predicted` and `predicted_coding`): a confidence derived from the size of the
effect, which is what R4 retires. Before the write changes, this script records who reads the field
and how many stored links carry it, so that no reader silently loses a number.

Two parts:

- READERS: every code site that reads the field, found by grep on 2026-09-28 (files that mention a
  predicted link and `confidence`), with what it does and whether it reads both the old and the new
  form. Fixture dictionaries in tests that merely construct the field are listed apart.
- stored: per chromosome, the committed sampled runs (data/results/enhancer_targets_<chrom>.json)
  and the local all-element tables (data/knowledge/alphagenome/all_elements/<chrom>.json, not
  committed): links with a gene, links carrying `confidence`, links whose confidence differs from
  |log2 fold change|, links where the 0.7 cap bit, links already carrying a `certainty` record.

Reads stored files only: no model request, no per-element cache.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.results import RESULTS_DIR, save_result

CHROMS = [f"chr{c}" for c in [*range(1, 23), "X", "Y"]]
ARCHIVE = Path("data/knowledge/alphagenome/all_elements")
LINKS = ("predicted", "predicted_coding")

READERS: list[dict[str, str]] = [
    {
        "site": "genomeos/attribution/compile.py element_certainty",
        "reads": "pc.get('confidence') as the Certainty model_score",
        "both_forms": "old yes; new: model_score None unless it reads certainty.model_score",
        "owner": "lane-ontology (reported, not edited here)",
    },
    {
        "site": "genomeos/attribution/confidence_calibration.py stated_confidence",
        "reads": "pc['confidence'], falling back to |log2 fold change|, for the retired formula",
        "both_forms": "yes (the fallback is the same number before the 0.7/0.05 clip it applies)",
        "owner": "lane-evid2",
    },
    {
        "site": "genomeos/attribution/confidence_calibration.py model_score",
        "reads": "log2_fold_change first, confidence only when no fold change was kept",
        "both_forms": "yes",
        "owner": "lane-evid2",
    },
    {
        "site": "genomeos/decompile.py _elements and render",
        "reads": "pc.get('confidence'), printed as '[predicted: AlphaGenome, <confidence>]'",
        "both_forms": "old yes; new would print None",
        "owner": "lane-evid2",
    },
    {
        "site": "genomeos/cli.py cmd enhancer-target (one element)",
        "reads": "p['confidence'] printed beside the fold change",
        "both_forms": "old yes; new KeyError",
        "owner": "lane-evid2",
    },
    {
        "site": "scripts/node_containment_audit.py census",
        "reads": "counts links with a confidence (has_confidence, predicted_coding_has_confidence)",
        "both_forms": "old yes; new counts 0",
        "owner": "lane-evid2",
    },
]
FIXTURES = [
    "tests/test_enhancer_target.py (asserts the written confidence 0.45)",
    "tests/test_decompile.py (link dict with confidence 0.4)",
    "tests/test_confidence_calibration.py (links with confidence 0.4321 and 0.99)",
    "tests/test_compile_certainty.py, tests/test_grn_bridge.py, tests/test_ontology.py (confidence=|lfc|)",
]
NOT_READERS = [
    "genomeos/web/static/index.html and genomeos/web/server.py: no predicted-link confidence read",
    "genomeos/attribution/target_calibration.py raw_confidence: recomputes min(1, |lfc|) from features",
    "genomeos/attribution/measured.py, candidates.py, organise.py: their own confidences, not the link's",
    "scripts/enhancer_targets*.py: copy the link dict whole (compact), read no confidence",
]
RELATED_LEFTOVERS = [
    "genomeos/predict/alphagenome_adapter.py VariantEffect.confidence: magnitude-based, capped at 0.7",
    "genomeos/genome/missense.py: confidence = round(min(0.7, score), 3) on missense predictions",
]


def count_links(elements: list[dict[str, Any]]) -> dict[str, int]:
    out = {
        "links": 0,
        "with_confidence": 0,
        "confidence_differs_from_abs_lfc": 0,
        "cap_applied": 0,
        "with_certainty": 0,
        "predicted_coding_links": 0,
        "predicted_coding_differs_from_abs_lfc": 0,
    }
    for e in elements:
        for key in LINKS:
            p = e.get(key) or {}
            if not p.get("gene"):
                continue
            out["links"] += 1
            out["predicted_coding_links"] += key == "predicted_coding"
            if p.get("certainty") is not None:
                out["with_certainty"] += 1
            c = p.get("confidence")
            if c is None:
                continue
            out["with_confidence"] += 1
            lfc = p.get("log2_fold_change")
            if lfc is not None:
                out["confidence_differs_from_abs_lfc"] += c != abs(lfc)
                if key == "predicted_coding":
                    out["predicted_coding_differs_from_abs_lfc"] += c != abs(lfc)
                out["cap_applied"] += abs(lfc) > 0.7
    return out


def census(chroms: list[str], results_dir: Path = RESULTS_DIR, archive: Path = ARCHIVE) -> dict[str, Any]:
    committed: dict[str, dict[str, int]] = {}
    local: dict[str, dict[str, int]] = {}
    for chrom in chroms:
        p = results_dir / f"enhancer_targets_{chrom}.json"
        if p.exists():
            committed[chrom] = count_links(json.loads(p.read_text()).get("elements", []))
        a = archive / f"{chrom}.json"
        if a.exists():
            local[chrom] = count_links(json.loads(a.read_text()))

    def total(d: dict[str, dict[str, int]]) -> dict[str, int]:
        keys = next(iter(d.values()), {}).keys()
        return {k: sum(v[k] for v in d.values()) for k in keys}

    return {
        "readers": READERS,
        "fixtures": FIXTURES,
        "not_readers": NOT_READERS,
        "related_leftovers": RELATED_LEFTOVERS,
        "committed_sampled_runs": {"per_chromosome": committed, "total": total(committed)},
        "local_all_element_tables": {"per_chromosome": local, "total": total(local)},
        "reading": (
            "stored files keep their values: the write changes for new runs only, and every reader "
            "listed accepts a link with the old `confidence` or with the new `certainty` record"
        ),
    }


def manifest(chroms: list[str], results_dir: Path = RESULTS_DIR, archive: Path = ARCHIVE) -> dict[str, Any]:
    inputs = [
        mf.input_entry(p, partition=None)
        for chrom in chroms
        for p in (results_dir / f"enhancer_targets_{chrom}.json", archive / f"{chrom}.json")
        if p.exists()
    ]
    return {
        "sources": [
            {
                "accession": "this repository, enhancer_targets_<chrom> runs and the all-element tables",
                "version": "pinned by sha256",
            }
        ],
        "inputs": inputs,
        "assembly": "GRCh38",
        "coordinates": "n/a: counts of fields, no interval is read",
        "parameters": {"cap": 0.7, "links": list(LINKS)},
        "exclusions": ["links with no gene are not counted"],
        "partitions": "n/a: no evaluation partition",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--chroms", default="")
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    chroms = args.chroms.split(",") if args.chroms else CHROMS
    out = census(chroms)
    for part in ("committed_sampled_runs", "local_all_element_tables"):
        print(part, out[part]["total"])
    if not args.no_save:
        print(f"wrote {save_result('predicted_confidence_census', out, manifest=manifest(chroms))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
