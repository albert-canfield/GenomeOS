# SPDX-License-Identifier: AGPL-3.0-or-later
"""The attributed elements of a chromosome, from whichever deletion runs exist.

Three runs name a target gene per element: the whole-chromosome scoring
(`enhancer_targets_all_<chrom>`, every element in a node; its committed result is
a summary whose `elements_where` points at the local table), the constrained
sample and the uniform sample. The whole-chromosome run is read first; the samples
add only elements it does not hold. The closure, the organiser and the compiler
all read through here, so none of them silently falls back to a sample when the
whole table exists.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from genomeos.results import RESULTS_DIR, load_result

RUNS = ("enhancer_targets_all", "constrained_targets", "enhancer_targets")
ORIGIN = {"enhancer_targets_all": "all", "constrained_targets": "constrained", "enhancer_targets": "uniform"}


def run_elements(name: str, chrom: str, results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    """One run's elements: inline, or from the local table a summary points at. A summary whose
    table is missing raises, rather than reading as a run with no elements."""
    r = load_result(f"{name}_{chrom}", results_dir) or {}
    if "elements" in r:
        return r["elements"]
    where = r.get("elements_where")
    if not where:
        return []
    p = Path(where)
    if not p.is_absolute() and not p.exists():
        p = results_dir.parent.parent / where  # a relative path is relative to the project root
    if not p.exists():
        raise FileNotFoundError(
            f"{name}_{chrom} points at {where}, which is not on this machine; rerun the job"
        )
    t = json.loads(p.read_text())
    return t if isinstance(t, list) else t.get("elements", [])


def attributed(chrom: str, results_dir: Path = RESULTS_DIR, coding: bool = True) -> list[dict[str, Any]]:
    """Every element with a named (coding) target across the runs, deduplicated by id, each with the
    run it came from under `origin`."""
    seen: set[str] = set()
    out: list[dict[str, Any]] = []
    for name in RUNS:
        for e in run_elements(name, chrom, results_dir):
            pc = e.get("predicted_coding" if coding else "predicted") or {}
            if not pc.get("gene") or e["id"] in seen:
                continue
            seen.add(e["id"])
            out.append({**e, "origin": ORIGIN[name]})
    return out


def runs_present(chrom: str, results_dir: Path = RESULTS_DIR) -> dict[str, int]:
    return {name: len(run_elements(name, chrom, results_dir)) for name in RUNS}
