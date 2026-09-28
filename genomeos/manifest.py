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

#: fe0880a, the commit that introduced this contract (review R9). A result name is historical when it
#: was in the registry before it (item 12 S6, docs/DATA.md "The result registry").
CONTRACT_COMMIT = "fe0880ab7423b064a38da130e326fa57d4efeb28"
#: Item 12 S6: what counts as code in the revision stamp, tracked or untracked. tests/ because a writer
#: imports from it (scripts/grn_clamp_census.py); *.bio anywhere because writers execute BioLang
#: programs under data/organisms and data/demo. Untracked data files do not count (docs/DATA.md).
CODE_ROOTS = ("genomeos/", "scripts/", "tests/")
CODE_SUFFIXES = (".bio",)


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
    if "model_dependencies" in manifest:  # optional: stated where a result rests on a served model
        deps = manifest["model_dependencies"]
        if not isinstance(deps, list):
            problems.append("model dependencies must be a list, one model dependency per served model")
        else:
            for d in deps:
                if not _pins_or_says_why(d):
                    problems.append(
                        f"a model dependency needs a name and a model_version or unpinned: why: {d!r}"[:300]
                    )
    return problems


def _pins_or_says_why(d: Any) -> bool:
    if not isinstance(d, dict) or not d.get("name"):
        return False
    return bool(d.get("model_version")) or str(d.get("unpinned", "")).startswith("unpinned: ")


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


# --- model dependencies (R9 follow-up, 2026-09-28) ------------------------------------------------
#
# A result that rests on a served model cannot be rebuilt from its inputs' sha256 alone once those
# inputs are gone: a second machine must ask the model again, and gets the model the server serves
# that day. `model_dependencies` (optional in a manifest, checked by `validate` when present) says
# which model a result's inputs were made with, as far as the disk establishes it, and "unpinned: why"
# for each part it cannot. Everything below was read from files on this machine; no request was made.

LOCK = Path("uv.lock")


def lock_version(package: str, lock: Path = LOCK) -> str | None:
    """The version of `package` in uv.lock, or None when the lock or the package is not there."""
    try:
        text = lock.read_text()
    except OSError:
        return None
    for block in text.split("[[package]]"):
        lines = [x.strip() for x in block.strip().splitlines()]
        if f'name = "{package}"' in lines:
            for x in lines:
                if x.startswith("version = "):
                    return x.split("=", 1)[1].strip().strip('"')
    return None


def track_fingerprint(answers: Any) -> dict[str, Any]:
    """What cached deletion answers say about the tracks behind them, when the track table itself was
    never stored: how many tracks each answer read (its modal count) and the sorted set of tissue names
    its genes' largest moves were read on, with that set's sha256. A fingerprint of names observed in
    the answers, never a checksum of the model's output_metadata."""
    counts: dict[Any, int] = {}
    names: set[str] = set()
    n = 0
    for a in answers:
        n += 1
        counts[a.get("tracks")] = counts.get(a.get("tracks"), 0) + 1
        for g in a.get("genes") or []:
            for k in ("max_drop_tissue", "max_rise_tissue"):
                if g.get(k):
                    names.add(str(g[k]))
    modal = max(counts.items(), key=lambda kv: (kv[1], str(kv[0])))[0] if counts else None
    return {
        "elements_scanned": n,
        "tracks_per_element_modal": modal,
        "elements_at_modal": counts.get(modal, 0),
        "distinct_track_counts": len(counts),
        "tissue_names": len(names),
        "tissue_names_sha256": hashlib.sha256("\n".join(sorted(names)).encode()).hexdigest(),
    }


#: The AlphaGenome all-element deletion sweep (scripts/enhancer_targets_all.py, chained by
#: scripts/enhancer_targets_all_chain.py) and its per-element response cache, as established from disk
#: on 2026-09-28. Every value names its evidence; nothing was asked of the model.
ALPHAGENOME_SWEEP: dict[str, Any] = {
    "name": "AlphaGenome",
    "what": "the 2026-09 all-element deletion sweep (scripts/enhancer_targets_all.py) and its per-element "
    "response cache (data/knowledge/alphagenome/elements/chr*.json.gz, all_elements/chr*.json)",
    "client": {
        "package": "alphagenome",
        "version": "0.9.0",
        "wheel_sha256": "a4f35884341ae85b5d2cf088dfe0304961de7ae4a590d4538653069673de32f4",
        "evidence": "uv.lock [[package]] alphagenome, version 0.9.0 with this wheel hash since df1a184 "
        "(2026-09-10), the only commit that ever changed that entry, so every sweep run resolved it; "
        "the installed client says __version__ 0.9.0",
    },
    "api": {
        "service": "google.gdm.gdmscience.alphagenome.v1main.DnaModelService",
        "method": "ScoreVariant",
        "address": "dns:///gdmscience.googleapis.com:443",
        "evidence": "alphagenome 0.9.0 protos/dna_model_service.proto (package ...alphagenome.v1main) and "
        "dna_client.create's default address; neither alphagenome_adapter._live_scorer nor "
        "enhancer_targets_all.worker_scorer passes an address",
    },
    "scorer": {
        "name": "variant_scorers.RECOMMENDED_VARIANT_SCORERS['RNA_SEQ']",
        "repr": "GeneMaskLFCScorer(requested_output=RNA_SEQ)",
        "window": "dna_client.SEQUENCE_LENGTH_1MB around the deleted element",
        "evidence": "alphagenome_adapter._live_scorer, unchanged from d1d5652 (2026-09-11) through the "
        "sweep; the repr is the installed 0.9.0 client's",
    },
    "model_version": None,
    "unpinned": "unpinned: no model version was requested. dna_client.create(api_key) and "
    "create(api_key, timeout=300) leave model_version None, so every ScoreVariantRequest carried an "
    "empty model_version and the server chose. No response message carries a model version "
    "(model_version is a request field only in dna_model_service.proto) and no cache file, result or "
    "job log records one, so which model answered cannot be established from disk",
    "documented_default": "ALL_FOLDS, the distilled all-folds model, per the comment on "
    "alphagenome 0.9.0 dna_model.ModelVersion; a client document, not verified against the server",
    "run_dates": {
        "first": "2026-09-12",
        "last": "2026-09-16",
        "evidence": "answers are kept only with per-cell fields, added in 257dadd (2026-09-12 12:27; "
        "has_cells deletes older answers); chr21 committed 5701321 (2026-09-13 02:29); "
        "data/jobs/enhancer_targets_all_chain.log: started 2026-09-13 00:48, 'chain done: every "
        "chromosome scored' 2026-09-16 15:14; the 23 archives' mtimes 2026-09-13 10:23 to 2026-09-16 "
        "15:10; result files dated 2026-09-13 to 2026-09-16",
        "later_answers": "3,209 per-element files outside the archives, written 2026-09-17 to "
        "2026-09-21 by later deletion runs through the same unpinned client (file mtimes)",
    },
    "track_metadata": {
        "stored": False,
        "sha256": None,
        "unpinned": "unpinned: the sweep never read output_metadata and kept per gene only its mean, "
        "largest drop, largest rise and four cell lines' values, so the track table it averaged over is "
        "not on disk. The one later read (2026-09-27, genomeos/attribution/crispri.py) kept three cell "
        "lines' track counts, not the table",
        "observed": {
            "elements_scanned": 966615,
            "tracks_per_element_modal": 371,
            "elements_at_modal": 919248,
            "distinct_track_counts": 422,
            "tissue_names": 316,
            "tissue_names_sha256": "bdf63a526775766a713d2e197a5b216876388c7db6f9be7a34435b97c14a4f62",
            "evidence": "track_fingerprint over every answer in elements/chr*.json.gz (963,406) and the "
            "3,209 loose element files, read 2026-09-28",
        },
    },
}

MODEL_DEPENDENCIES: dict[str, dict[str, Any]] = {"alphagenome": ALPHAGENOME_SWEEP}


def model_dependency(name: str) -> dict[str, Any]:
    """The recorded block for one served model (KeyError for a model with none), as a fresh copy, with
    the client version uv.lock holds now beside the one the runs were made with."""
    import copy

    block = copy.deepcopy(MODEL_DEPENDENCIES[name.lower()])
    block["client"]["version_in_lock_now"] = lock_version(block["client"]["package"])
    return block


#: The label of an answer made with no model version requested (review R9 follow-up, 2026-09-28).
MODEL_VERSION_UNREQUESTED = "unrequested"


def answers_model_dependency(name: str, counts: dict[str, int]) -> dict[str, Any]:
    """The model dependency of a result built from cached answers, `counts` being answers per requested
    model version (MODEL_VERSION_UNREQUESTED for none). Answers asked of no version and answers asked of
    a named one are not guaranteed to come from the same model, so a result holding both says so in
    `mixed`; one holding only named answers of one version pins it."""
    block = model_dependency(name)
    counts = {str(k): int(v) for k, v in counts.items() if v}
    block["answers_by_model_version"] = counts
    named = sorted(k for k in counts if k != MODEL_VERSION_UNREQUESTED)
    unrequested = counts.get(MODEL_VERSION_UNREQUESTED, 0)
    block["mixed"] = None
    if len(named) + bool(unrequested) > 1:
        parts = [f"{counts[k]:,} asked for {k}" for k in named]
        if unrequested:
            parts.insert(0, f"{unrequested:,} asked for no model version (the server chose)")
        block["mixed"] = (
            "mixed: of the answers this result reads, "
            + "; ".join(parts)
            + ". They are not guaranteed to come from the same model, and nothing on disk says "
            "whether they did"
        )
    elif named:
        block["model_version"] = named[0]
        block["unpinned"] = None
        block["model_version_evidence"] = (
            f"every answer this result reads carries a run record that requested {named[0]} "
            "(dna_client.create(api_key, model_version=...)); the client and api fields above describe "
            "the earlier unpinned sweep and the same client"
        )
    return block


def with_model_dependencies(manifest: dict[str, Any], *names: str) -> dict[str, Any]:
    """A manifest with `model_dependencies` added for the named models; the caller's dict is untouched."""
    return {**manifest, "model_dependencies": [model_dependency(n) for n in names]}


def depends_on_models(*names: str):
    """Decorator for a writer's manifest function: its manifests carry `model_dependencies` for the named
    models (with_model_dependencies), so a writer states its model in one added line."""

    def wrap(fn):
        import functools

        @functools.wraps(fn)
        def inner(*a: Any, **k: Any) -> dict[str, Any]:
            return with_model_dependencies(fn(*a, **k), *names)

        return inner

    return wrap
