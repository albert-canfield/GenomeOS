# SPDX-License-Identifier: AGPL-3.0-or-later
"""Result manifests: what a result was made from, so a second environment can make it again.

Review item R9 (2026-09-28) found that `genomeos/results.py` wrote every result without a provenance
contract: nothing recorded which release of which source, which bytes, which assembly, which
coordinate convention, which revision of the code, which parameters, what was left out and which
evaluation partition each input belongs to. This module is that contract.

A manifest is a dict kept under the result's `result_manifest` key (26 older results already use a
`manifest` key for something else: a file name or a path) with eight required fields:

    sources      [{accession, version, ...}]   the public records the inputs come from
    inputs       [{path, sha256, bytes, partition}]   the bytes read, with R5's partition where known
    assembly     "GRCh38", or "n/a: <why>"
    coordinates  {base: 0 | 1, interval: "half-open" | "closed"}, or "n/a: <why>"
    code         {git_sha, dirty, dirty_code_paths, argv}   filled in by the writer, never by the caller
                 (and dirty_result_paths: results that differ from the commit, which are output, not code)
    parameters   {name: value}                  every knob that changes the numbers
    exclusions   [str | dict]                   what was dropped and why; [] says nothing was
    partitions   {name: description}, or "n/a: <why>"   the evaluation partitions the result reports

`validate` names what a manifest lacks. `read` is tolerant: it accepts the 955 results written before
the contract, reports which fields they carry, and never raises. The historical files are not
rewritten; `scripts/manifest_census.py` counts what they already hold.
"""

from __future__ import annotations

import contextlib
import hashlib
import subprocess
import sys
from pathlib import Path
from typing import Any

KEY = "result_manifest"
REQUIRED = ("sources", "inputs", "assembly", "coordinates", "code", "parameters", "exclusions", "partitions")
INTERVALS = ("half-open", "closed")
NOT_APPLICABLE = "n/a"

# Top-level keys that older results used for part of what a manifest now states. Presence of one is
# weak evidence that the field is covered in some form; the census reports it as "legacy", never as
# a declared manifest field.
LEGACY_KEYS: dict[str, tuple[str, ...]] = {
    "sources": ("sources", "source", "accession", "accessions", "release", "data_source", "gse", "encode"),
    "inputs": ("sha256", "md5", "checksum", "checksums", "inputs", "input", "files"),
    "assembly": ("assembly", "genome", "build", "reference", "genome_build"),
    "coordinates": ("coordinates", "coordinate_system", "zero_based", "half_open", "coordinate_convention"),
    "code": ("git_sha", "commit", "code_revision", "revision", "git"),
    "parameters": ("parameters", "params", "thresholds", "threshold", "config", "settings", "method"),
    "exclusions": ("exclusions", "excluded", "exclude", "dropped", "skipped", "filtered_out"),
    "partitions": ("partitions", "partition", "split", "folds", "heldout", "held_out", "training", "splits"),
}


class ManifestError(ValueError):
    """A new result that does not carry the provenance contract."""


def _git(root: Path, *args: str) -> str | None:
    try:
        out = subprocess.run(
            ["git", "-C", str(root), *args], capture_output=True, text=True, timeout=30, check=True
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout


#: tracked paths that hold results, not code: a result an earlier writer rewrote in the same checkout
#: is output, so it is recorded under `dirty_result_paths` and does not make the code dirty.
RESULT_PATHS = ("data/results/",)


def code_revision(root: Path | None = None) -> dict[str, Any]:
    """The revision of the code that wrote a result. `dirty` counts tracked files that differ from
    the commit, because in a shared checkout the sha alone does not say what ran."""
    # Since 2026-09-28 a modified file under data/results/ is not counted: a chain of writers in one
    # checkout would otherwise mark every result after the first dirty with its predecessor's output.
    # Those paths are recorded apart, under dirty_result_paths.
    root = root or Path.cwd()
    sha = _git(root, "rev-parse", "HEAD")
    status = _git(root, "status", "--porcelain", "--untracked-files=no")
    dirty_paths = sorted(line[3:] for line in (status or "").splitlines() if line.strip())
    result_paths = [p for p in dirty_paths if p.startswith(RESULT_PATHS)]
    dirty_paths = [p for p in dirty_paths if not p.startswith(RESULT_PATHS)]
    code_paths = [p for p in dirty_paths if p.startswith(("genomeos/", "scripts/"))]
    return {
        "git_sha": sha.strip() if sha else None,
        "dirty": bool(dirty_paths) if status is not None else None,
        "dirty_code_paths": code_paths,
        "dirty_result_paths": result_paths,
        "argv": _argv(root),
        "python": sys.version.split()[0],
    }


def _argv(root: Path) -> list[str]:
    """The command line, the script as a path relative to the repository so another checkout can run it."""
    argv = list(sys.argv)
    if argv:
        with contextlib.suppress(ValueError, OSError):
            argv[0] = str(Path(argv[0]).resolve().relative_to(root.resolve()))
    return argv


def sha256_of(path: Path) -> tuple[str, int, int]:
    """sha256, byte count and file count. A directory hashes its sorted files' relative paths and
    contents, so the same tree gives the same digest wherever it sits."""
    h = hashlib.sha256()
    files = sorted(p for p in path.rglob("*") if p.is_file()) if path.is_dir() else [path]
    total = 0
    for p in files:
        if path.is_dir():
            h.update(str(p.relative_to(path)).encode() + b"\0")
        with p.open("rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
                total += len(chunk)
    return h.hexdigest(), total, len(files)


def input_entry(path: str | Path, partition: str | None = None, **extra: Any) -> dict[str, Any]:
    """One input as the manifest records it. `partition` is R5's field: which evaluation partition
    the input came from (training, heldout, ...), None when not known."""
    p = Path(path)
    digest, size, count = sha256_of(p)
    entry: dict[str, Any] = {"path": str(path), "sha256": digest, "bytes": size, "partition": partition}
    if p.is_dir():
        entry["files"] = count
    return {**entry, **extra}


def _not_applicable(v: Any) -> bool:
    return isinstance(v, str) and v.startswith(NOT_APPLICABLE) and len(v) > len(NOT_APPLICABLE) + 2


def validate(manifest: Any) -> list[str]:
    """What a manifest lacks or gets wrong, as sentences; empty when it meets the contract."""
    if not isinstance(manifest, dict):
        return ["no manifest"]
    problems = [f"missing {f}" for f in REQUIRED if f not in manifest]
    src = manifest.get("sources")
    if "sources" in manifest:
        if not isinstance(src, list) or not src:
            problems.append("sources must be a non-empty list")
        else:
            for s in src:
                if not isinstance(s, dict) or not s.get("accession") or not s.get("version"):
                    problems.append(f"a source needs accession and version: {s!r}")
    inp = manifest.get("inputs")
    if "inputs" in manifest:
        if not isinstance(inp, list) or not inp:
            problems.append("inputs must be a non-empty list")
        else:
            for i in inp:
                if (
                    not isinstance(i, dict)
                    or not i.get("path")
                    or not i.get("sha256")
                    or "partition" not in i
                ):
                    problems.append(f"an input needs path, sha256 and partition (None if unknown): {i!r}")
    a = manifest.get("assembly")
    if "assembly" in manifest and not (isinstance(a, str) and a.strip()):
        problems.append("assembly must name a build, or say n/a: why")
    c = manifest.get("coordinates")
    well_formed = isinstance(c, dict) and c.get("base") in (0, 1) and c.get("interval") in INTERVALS
    if "coordinates" in manifest and not _not_applicable(c) and not well_formed:
        problems.append("coordinates must be {base: 0|1, interval: half-open|closed}, or n/a: why")
    code = manifest.get("code")
    if "code" in manifest and not (isinstance(code, dict) and code.get("git_sha") and "dirty" in code):
        problems.append("code must carry git_sha and dirty")
    if "parameters" in manifest and not isinstance(manifest["parameters"], dict):
        problems.append("parameters must be a dict")
    if "exclusions" in manifest and not isinstance(manifest["exclusions"], list):
        problems.append("exclusions must be a list ([] when nothing was excluded)")
    part = manifest.get("partitions")
    if "partitions" in manifest and not (isinstance(part, dict) or _not_applicable(part)):
        problems.append("partitions must be a dict, or n/a: why")
    return problems


def stamp(manifest: dict[str, Any] | None, root: Path | None = None) -> dict[str, Any]:
    """The caller's manifest with the code revision the writer fills in and the verdict attached."""
    m = {k: v for k, v in (manifest or {}).items() if k not in ("complete", "problems")}
    m["code"] = code_revision(root)
    problems = validate(m)
    m["complete"] = not problems
    if problems:
        m["problems"] = problems
    return m


def read(payload: Any) -> dict[str, Any]:
    """Tolerant reader: what a result says about its own provenance, whenever it was written.

    Never raises. `declared` is True when the result carries a `result_manifest` dict; `fields` gives, per
    required field, "declared", "legacy" (an older top-level key that covers part of it) or None.
    """
    if not isinstance(payload, dict):
        return {"declared": False, "complete": False, "fields": dict.fromkeys(REQUIRED), "manifest": None}
    m = payload.get(KEY)
    declared = isinstance(m, dict)
    fields: dict[str, str | None] = {}
    for f in REQUIRED:
        if declared and f in m:
            fields[f] = "declared"
        elif any(k in payload for k in LEGACY_KEYS[f]):
            fields[f] = "legacy"
        else:
            fields[f] = None
    return {
        "declared": declared,
        "complete": declared and not validate({k: v for k, v in m.items() if k in REQUIRED}),
        "fields": fields,
        "manifest": m if declared else None,
    }
