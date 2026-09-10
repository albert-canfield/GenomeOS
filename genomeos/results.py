"""Results registry: the distilled outcome of every real-data run.

The storage principle of the project is *stream, distil, discard*: raw data
(genomes, VCFs, methylation matrices, read files) is either streamed or
downloaded once, the test or analysis is run, and only a small JSON summary
is kept under data/results/. Those summaries are committed, feed the tests
when the raw inputs are absent, and are what eventually gets embedded in the
software as extracted rules and calibrations.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

RESULTS_DIR = Path("data/results")


def save_result(name: str, payload: dict[str, Any], results_dir: Path = RESULTS_DIR) -> Path:
    results_dir.mkdir(parents=True, exist_ok=True)
    payload = {"result": name, "date": time.strftime("%Y-%m-%d"), **payload}
    p = results_dir / f"{name}.json"
    p.write_text(json.dumps(payload, indent=2, default=str))
    return p


def load_result(name: str, results_dir: Path = RESULTS_DIR) -> dict[str, Any] | None:
    p = results_dir / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def list_results(results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    out = []
    for p in sorted(results_dir.glob("*.json")) if results_dir.is_dir() else []:
        try:
            d = json.loads(p.read_text())
        except json.JSONDecodeError:
            continue
        out.append({"name": p.stem, "date": d.get("date", ""), "size": p.stat().st_size, "keys": list(d)[:8]})
    return out
