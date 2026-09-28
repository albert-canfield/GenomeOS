# SPDX-License-Identifier: AGPL-3.0-or-later
"""Census of the last two magnitude-derived confidences, before review R4 retires them (R4f).

- `predict/alphagenome_adapter.py` `PredictedEffect.confidence = min(0.7, 0.2 + 0.25 * |log2 fold
  change|)`: a variant's predicted expression change turned into a confidence by its size.
- `genome/missense.py` writes `predicted.confidence = round(min(0.7, score), 3)`: the AlphaMissense
  pathogenicity score capped and called a confidence.

Before either write changes, this script records every code site that reads the two fields and how
many stored values carry them, so that no reader silently loses a number. READERS were found by grep
on 2026-09-28 (`e.confidence`, `PredictedEffect`, `to_module(` of the adapter, `missense_by_effect`,
`["predicted"]`). Stored: the committed `alphagenome_rs12740374` result (the CLI's --json output for
one variant) and the per-person AlphaMissense caches under data/individuals/<name>/ (git-ignored raw
scores; the confidence is computed when a person is read, never stored there).

Reads stored files only: no model request, no streaming.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from genomeos import manifest as mf
from genomeos.results import RESULTS_DIR, save_result

PEOPLE = Path("data/individuals")
STORED_VARIANT = "alphagenome_rs12740374"
CAP = 0.7

READERS: list[dict[str, str]] = [
    {
        "site": "genomeos/predict/alphagenome_adapter.py AlphaGenomeAdapter.to_module",
        "reads": "PredictedEffect.confidence as each predicted Rule's confidence",
        "callers": "tests/test_alphagenome.py only; no CLI, server or script builds this module",
    },
    {
        "site": "genomeos/cli.py cmd_predict (variant form)",
        "reads": "round(e.confidence, 3) into --json effects[], and 'confidence 0.xx' per printed line",
        "callers": "`genomeos predict chr:pos REF>ALT`",
    },
    {
        "site": "genomeos/web/server.py predict",
        "reads": "round(e.confidence, 3) into effects[] of /api/predict",
        "callers": "the /api/predict endpoint; no page of web/static/index.html calls it",
    },
    {
        "site": "genomeos/predict/alphagenome_adapter.py LICENCE_NOTE",
        "reads": "states predictions enter 'capped at confidence 0.7' (text, shown by status)",
        "callers": "`genomeos predict --status`, /api/predict/status",
    },
    {
        "site": "tests/test_alphagenome.py test_adapter_emits_predicted_rules_only",
        "reads": "pins rule confidence <= 0.7 and > 0.35 for |log2FC| 0.8: a pin on the retired formula",
        "callers": "test",
    },
    {
        "site": "genomeos/genome/missense.py annotate",
        "reads": "writes predicted.confidence = round(min(0.7, score), 3) on each scored missense row",
        "callers": "individuals.coding_variants (predict=True) -> CLI coding table, report, web panel",
    },
    {
        "site": "genomeos/cli.py, genomeos/genome/individuals.py, genomeos/web/static/index.html",
        "reads": "missense_by_effect rows: score, class, protein_variant; none reads predicted.confidence",
        "callers": "coding-variant views",
    },
    {
        "site": "tests/test_missense.py test_annotate",
        "reads": "pins predicted.confidence == 0.7 for a score of 0.95",
        "callers": "test",
    },
    {
        "site": "genomeos/attribution/compile.py element_certainty",
        "reads": "pc.get('confidence') as the model score of an enhancer link (R4e leftover)",
        "callers": "tests/test_compile_certainty.py; the compiled text does not emit the per-link score",
    },
]


def stored_variant(results_dir: Path = RESULTS_DIR) -> dict[str, Any]:
    p = results_dir / f"{STORED_VARIANT}.json"
    if not p.exists():
        return {"present": False}
    d = json.loads(p.read_text())
    effects = d.get("top_effects", [])
    conf = [e for e in effects if e.get("confidence") is not None]
    formula = [
        e for e in conf if round(min(CAP, 0.2 + 0.25 * abs(e["log2_fold_change"])), 3) == e["confidence"]
    ]
    return {
        "present": True,
        "effects": len(effects),
        "with_confidence": len(conf),
        "equal_to_the_retired_formula": len(formula),
        "cap_applied": sum(0.2 + 0.25 * abs(e["log2_fold_change"]) > CAP for e in conf),
        "confidence_range": [min(e["confidence"] for e in conf), max(e["confidence"] for e in conf)]
        if conf
        else None,
    }


def people_caches(root: Path = PEOPLE) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for p in sorted(root.glob("*/alphamissense.json")):
        d = json.loads(p.read_text())
        entries = [e for v in d.get("scores", {}).values() for e in v]
        out[p.parent.name] = {
            "variants_with_a_score": len(d.get("scores", {})),
            "transcript_entries": len(entries),
            "entries_above_cap": sum(e["score"] > CAP for e in entries),
            "stores_confidence": sum("confidence" in e for e in entries),
        }
    return out


def census(results_dir: Path = RESULTS_DIR, people: Path = PEOPLE) -> dict[str, Any]:
    return {
        "readers": READERS,
        "stored_variant_result": stored_variant(results_dir),
        "people_alphamissense_caches": people_caches(people),
        "reading": (
            "the variant-effect confidence is stored once, in the committed alphagenome_rs12740374 "
            "result, which keeps its values; the missense confidence is never stored, only computed "
            "when a person's coding variants are read, so retiring it moves no stored figure"
        ),
    }


def manifest(results_dir: Path = RESULTS_DIR, people: Path = PEOPLE) -> dict[str, Any]:
    paths = [results_dir / f"{STORED_VARIANT}.json", *sorted(people.glob("*/alphamissense.json"))]
    return {
        "sources": [
            {"accession": "this repository and the local per-person caches", "version": "pinned by sha256"}
        ],
        "inputs": [mf.input_entry(p, partition=None) for p in paths if p.exists()],
        "assembly": "GRCh38",
        "coordinates": "n/a: counts of fields, no interval is read",
        "parameters": {"cap": CAP, "adapter_formula": "min(0.7, 0.2 + 0.25*|log2FC|)"},
        "exclusions": ["none"],
        "partitions": "n/a: no evaluation partition",
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-save", action="store_true")
    args = ap.parse_args(argv)
    out = census()
    print("stored variant result", out["stored_variant_result"])
    print("people caches", out["people_alphamissense_caches"])
    if not args.no_save:
        print(f"wrote {save_result('variant_confidence_census', out, manifest=manifest())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
