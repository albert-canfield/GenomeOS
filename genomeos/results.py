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
import warnings
from pathlib import Path
from typing import Any

from genomeos import manifest as mf

RESULTS_DIR = Path("data/results")


class ManifestWarning(UserWarning):
    """A rewritten historical result that does not yet carry the full provenance contract."""


def save_result(
    name: str,
    payload: dict[str, Any],
    results_dir: Path = RESULTS_DIR,
    manifest: dict[str, Any] | None = None,
    strict: bool | None = None,
) -> Path:
    """Write a result with its manifest (genomeos/manifest.py, review item R9).

    Every write is stamped: the code revision is always filled in, the caller's manifest (the
    `manifest` argument, or a `result_manifest` key in the payload) is checked against the contract,
    and the verdict is kept as `result_manifest.complete` with `result_manifest.problems`.

    Enforcement is for new results. A result name that does not yet exist in the project's registry
    must meet the contract (`strict`, which defaults to exactly that case): the file is still
    written, so hours of compute are never lost, and then ManifestError names what is missing and
    the run fails. Rewriting a historical result, or writing outside the registry, warns instead,
    because the 955 results written before the contract must keep regenerating.
    """
    results_dir.mkdir(parents=True, exist_ok=True)
    p = results_dir / f"{name}.json"
    body = dict(payload)
    given = manifest if manifest is not None else body.pop(mf.KEY, None)
    body.pop(mf.KEY, None)
    if strict is None:
        strict = not p.exists() and _is_registry(results_dir)
    stamped = mf.stamp(given)
    out = {"result": name, "date": time.strftime("%Y-%m-%d"), **body, mf.KEY: stamped}
    p.write_text(json.dumps(out, indent=2, default=str))
    if not stamped["complete"]:
        msg = f"{p}: manifest incomplete: {'; '.join(stamped['problems'])}"
        if strict:
            raise mf.ManifestError(msg + " (written, but a new result must carry the contract)")
        warnings.warn(msg, ManifestWarning, stacklevel=2)
    return p


def _is_registry(results_dir: Path) -> bool:
    try:
        return results_dir.resolve() == RESULTS_DIR.resolve()
    except OSError:
        return False


def load_result(name: str, results_dir: Path = RESULTS_DIR) -> dict[str, Any] | None:
    p = results_dir / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else None


def load_manifest(name: str, results_dir: Path = RESULTS_DIR) -> dict[str, Any] | None:
    """What a result says about its provenance, read tolerantly (older results carry no manifest)."""
    r = load_result(name, results_dir)
    return None if r is None else mf.read(r)


def list_results(results_dir: Path = RESULTS_DIR) -> list[dict[str, Any]]:
    out = []
    for p in sorted(results_dir.glob("*.json")) if results_dir.is_dir() else []:
        try:
            d = json.loads(p.read_text())
        except json.JSONDecodeError:
            continue
        out.append({"name": p.stem, "date": d.get("date", ""), "size": p.stat().st_size, "keys": list(d)[:8]})
    return out
