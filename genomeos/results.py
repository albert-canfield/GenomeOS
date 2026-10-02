"""Results registry: the distilled outcome of every real-data run.

The storage principle of the project is *stream, distil, discard*: raw data
(genomes, VCFs, methylation matrices, read files) is either streamed or
downloaded once, the test or analysis is run, and only a small JSON summary
is kept under data/results/. Those summaries are committed, feed the tests
when the raw inputs are absent, and are what eventually gets embedded in the
software as extracted rules and calibrations.
"""

from __future__ import annotations

import contextlib
import json
import time
import warnings
from pathlib import Path
from typing import Any

from genomeos import manifest as mf

RESULTS_DIR = Path("data/results")
#: Item 12 S6 (docs/DATA.md, "The result registry"): a failed new result is written to
#: <results_dir's parent>/QUARANTINE/<results_dir's name>, never into the registry.
QUARANTINE = "quarantine"
#: The explicit list of the result names written before the contract (fe0880a), generated once from
#: git history by scripts/manifest_legacy.py and never extended by code.
LEGACY_ALLOWLIST = Path("data/results_legacy.txt")
LEGACY_COUNT = 955


class ManifestWarning(UserWarning):
    """A rewritten historical result that does not yet carry the full provenance contract."""


def save_result(
    name: str,
    payload: dict[str, Any],
    results_dir: Path = RESULTS_DIR,
    manifest: dict[str, Any] | None = None,
    strict: bool | None = None,
    compact: bool = False,
) -> Path:
    """Write a result with its manifest (genomeos/manifest.py, review item R9), or quarantine it.

    Every write is stamped: the code revision is always filled in, the caller's manifest (the
    `manifest` argument, or a `result_manifest` key in the payload) is checked against the contract,
    and the verdict is kept as `result_manifest.complete` with `result_manifest.problems`.

    Since item 12 S6 (2026-09-28, docs/DATA.md "The result registry") a result enters the registry
    only when its manifest is complete or its name is on the legacy allowlist (LEGACY_ALLOWLIST: the
    names written before the contract, a committed list). A rewritten legacy result warns, so the
    historical writers keep regenerating. Any other incomplete result in the registry is written to
    the quarantine (`quarantine_dir`) with its reason, never to `results_dir`, and ManifestError names
    where it went, so hours of compute are kept and the run still fails. Whether the file already
    exists does not matter: a retry of a failed result is quarantined again.

    Since 2026-10-02 (the supervisor's review order, `save_result` cleanliness) a name that is not on
    the allowlist must also carry `result_manifest.code_cleanliness` with every key of
    `manifest.CLEANLINESS_KEYS`, which is what the one shared `manifest.code_cleanliness` returns. A
    block built by hand is refused on the same path as any other incomplete manifest: quarantined with
    its reason, ManifestError naming where it went, nothing written to `results_dir`.

    `strict=True` enforces in any directory; `strict=False` cannot admit a name that is not on the
    allowlist into the registry. Outside the registry (tests, scratch) an incomplete result warns
    unless `strict=True`.

    `compact=True` writes the JSON without whitespace (`separators=(",", ":")`), for the per-cell and
    per-site tables that were written compact before they came under the contract (item 12 S6
    follow-up, lane-contract); the values are the same either way.
    """
    p = results_dir / f"{name}.json"
    body = dict(payload)
    given = manifest if manifest is not None else body.pop(mf.KEY, None)
    body.pop(mf.KEY, None)
    registry = _is_registry(results_dir)
    legacy = registry and name in legacy_names()
    enforce = True if registry and not legacy else bool(strict)
    stamped = mf.stamp(given)
    # A new name entering the registry must also carry the cleanliness block, from the one shared
    # function (genomeos/manifest.py: code_cleanliness). Every published result's honesty about which
    # code was uncommitted when it was written rests on that block; it was copy-pasted into ten files
    # whose answers had drifted apart, and nothing checked that they agreed. A legacy name is exempt,
    # as it is from the rest of the contract, so the historical writers keep regenerating.
    if registry and not legacy:
        unclean = mf.cleanliness_problems(stamped)
        if unclean:
            stamped["problems"] = [*stamped.get("problems", []), *unclean]
            stamped["complete"] = False
    out = {"result": name, "date": time.strftime("%Y-%m-%d"), **body, mf.KEY: stamped}
    q = quarantine_dir(results_dir) / f"{name}.json"
    if not stamped["complete"] and enforce:
        _quarantine(q, out, p, stamped["problems"], registry)
        why = "a name not on the legacy allowlist must carry the contract" if registry else "strict"
        raise mf.ManifestError(
            f"{p}: manifest incomplete: {'; '.join(stamped['problems'])} (nothing written to "
            f"{results_dir}; quarantined at {q}; {why})"
        )
    results_dir.mkdir(parents=True, exist_ok=True)
    fmt: dict[str, Any] = {"separators": (",", ":")} if compact else {"indent": 2}
    p.write_text(json.dumps(out, default=str, **fmt))
    if stamped["complete"]:
        with contextlib.suppress(OSError):
            q.unlink(missing_ok=True)  # an earlier failed attempt at this name is superseded
    else:
        warnings.warn(
            f"{p}: manifest incomplete: {'; '.join(stamped['problems'])}", ManifestWarning, stacklevel=2
        )
    return p


def legacy_names(path: Path | None = None) -> frozenset[str]:
    """The result names on the legacy allowlist (first column; `#` lines are comments). A list that
    cannot be read counts as empty, so enforcement fails closed."""
    try:
        text = (path or LEGACY_ALLOWLIST).read_text()
    except OSError:
        return frozenset()
    return frozenset(
        line.split("\t", 1)[0].strip()
        for line in text.splitlines()
        if line.strip() and not line.startswith("#")
    )


def quarantine_dir(results_dir: Path) -> Path:
    """Where a failed result meant for `results_dir` goes: data/quarantine/results for the registry."""
    return results_dir.parent / QUARANTINE / results_dir.name


_QUARANTINE_IGNORE = (
    "# Written by genomeos/results.py (item 12 S6): results that failed their manifest contract.\n"
    "# Never tracked: committing one would publish unvalidated numbers. See docs/DATA.md.\n*\n"
)


def _quarantine(q: Path, out: dict[str, Any], meant_for: Path, problems: list[str], registry: bool) -> None:
    root = q.parent.parent
    q.parent.mkdir(parents=True, exist_ok=True)
    ignore = root / ".gitignore"
    if not ignore.exists():
        ignore.write_text(_QUARANTINE_IGNORE)
    out = {
        **out,
        "quarantine": {
            "reason": "manifest incomplete, and the name is not on the legacy allowlist"
            if registry
            else "manifest incomplete, and strict was set",
            "problems": problems,
            "meant_for": str(meant_for),
            "written_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        },
    }
    q.write_text(json.dumps(out, indent=2, default=str))


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
