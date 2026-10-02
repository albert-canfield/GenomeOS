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
                 (and dirty_result_paths: results that differ from the commit, which are output, not code;
                 untracked_code_paths: code the commit does not hold, which makes it dirty; item 12 S6)
    parameters   {name: value}                  every knob that changes the numbers
    exclusions   [str | dict]                   what was dropped and why; [] says nothing was
    partitions   {name: description}, or "n/a: <why>"   the evaluation partitions the result reports

`validate` names what a manifest lacks. `read` is tolerant: it accepts the 955 results written before
the contract, reports which fields they carry, and never raises. The historical files are not
rewritten; `scripts/manifest_census.py` counts what they already hold.

Item 12 S6 (2026-09-28, docs/DATA.md "The result registry"): those 955 names are an explicit committed
list (data/results_legacy.txt, genomeos/results.py), a new result that fails `validate` is quarantined
rather than written, and the revision stamp counts untracked code.
"""

from __future__ import annotations

import ast
import contextlib
import hashlib
import subprocess
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path, PurePosixPath
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


def is_code(path: str) -> bool:
    """Whether a repository path counts as code in the revision stamp (CODE_ROOTS, CODE_SUFFIXES)."""
    return path.startswith(CODE_ROOTS) or path.endswith(CODE_SUFFIXES)


def _status_paths(status: str) -> tuple[list[str], list[str]]:
    """`git status --porcelain -z` as (tracked paths that differ from the commit, untracked paths).
    A rename or copy entry is followed by its original path, which is skipped."""
    tracked: list[str] = []
    untracked: list[str] = []
    entries = status.split("\0")
    i = 0
    while i < len(entries):
        entry = entries[i]
        i += 1
        if len(entry) < 4:
            continue
        xy, path = entry[:2], entry[3:]
        if xy == "??":
            untracked.append(path)
        elif xy != "!!":
            tracked.append(path)
            if "R" in xy or "C" in xy:
                i += 1
    return sorted(tracked), sorted(untracked)


def code_revision(root: Path | None = None) -> dict[str, Any]:
    """The revision of the code that wrote a result. `dirty` counts tracked files that differ from
    the commit and untracked code, because in a shared checkout the sha alone does not say what ran."""
    # Since 2026-09-28 a modified file under data/results/ is not counted: a chain of writers in one
    # checkout would otherwise mark every result after the first dirty with its predecessor's output.
    # Those paths are recorded apart, under dirty_result_paths.
    # Since item 12 S6 (2026-09-28) untracked files are read too (--untracked-files=all): one that is
    # code (is_code) makes the stamp dirty and is listed under untracked_code_paths and
    # dirty_code_paths, since the commit cannot reproduce a result an untracked script wrote. An
    # untracked data file does not count and is not listed: the manifest's inputs pin data by sha256.
    root = root or Path.cwd()
    sha = _git(root, "rev-parse", "HEAD")
    status = _git(root, "status", "--porcelain", "-z", "--untracked-files=all")
    tracked, untracked = _status_paths(status or "")
    result_paths = sorted(p for p in tracked + untracked if p.startswith(RESULT_PATHS))
    tracked = [p for p in tracked if not p.startswith(RESULT_PATHS)]
    untracked_code = [p for p in untracked if is_code(p) and not p.startswith(RESULT_PATHS)]
    return {
        "git_sha": sha.strip() if sha else None,
        "dirty": bool(tracked or untracked_code) if status is not None else None,
        "dirty_code_paths": sorted({p for p in tracked if is_code(p)} | set(untracked_code)),
        "untracked_code_paths": untracked_code,
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


def group_digest(names: list[str], root: Path | None = None) -> tuple[str, int, list[dict[str, Any]]]:
    """The digest of a group of files read together, the total bytes, and each member on its own.

    The group digest is sha256 over each member's name and bytes, in the order given, the way
    `sha256_of` hashes a directory, so the same files give the same digest wherever the checkout sits.
    `names` are already as the repository names them and already in the order they are to be hashed;
    `root` is the checkout to read them from, the working directory when it is None.

    The one implementation of this digest. `files_entry` records a group with it and
    `scripts/manifest_rebuild.py` checks a group with it, so the writer and the checker cannot answer
    differently about which bytes a group's digest covers.
    """
    h = hashlib.sha256()
    base = Path(root) if root is not None else Path()
    total = 0
    members: list[dict[str, Any]] = []
    for name in names:
        h.update(name.encode() + b"\0")
        member = hashlib.sha256()
        size = 0
        with (base / name).open("rb") as f:
            for chunk in iter(lambda f=f: f.read(1 << 20), b""):
                h.update(chunk)
                member.update(chunk)
                size += len(chunk)
        members.append({"path": name, "sha256": member.hexdigest(), "bytes": size})
        total += size
    return h.hexdigest(), total, members


def files_entry(label: str, paths: Any, partition: str | None = None, **extra: Any) -> dict[str, Any]:
    """One input made of several files read together (a person's per-chromosome calls, the tracked
    programs of a census): the group digest of `group_digest`, under a label rather than a path. Added
    by the item 12 S6 follow-up (lane-contract) for the writers it brought under the contract.

    Since 2026-10-02 the entry also names every file in the group, each with its own sha256 and byte
    count, and marks itself `group`. Until then it recorded only the label, the digest and a count, so
    nothing downstream could open the files it stood for: `scripts/manifest_rebuild.py` resolved an
    input by its `path`, found no file at the label, and reported an input absent that was on disk all
    along -- and a group was indistinguishable from `input_entry` on a directory, which also carries
    `files`. The group digest is unchanged, so a digest recorded before this date still verifies.
    """
    names = sorted({_repo_relative(x) for x in paths})
    digest, total, members = group_digest(names)
    return {
        "path": label,
        "group": True,
        "members": members,
        "sha256": digest,
        "bytes": total,
        "files": len(names),
        "partition": partition,
        **extra,
    }


def _repo_relative(path: str | Path) -> str:
    """A path as the repository names it (relative to the working directory) when it lies inside it, so
    the digest of a set of files does not depend on where the checkout sits."""
    import os

    with contextlib.suppress(ValueError, OSError):
        return str(Path(os.path.abspath(path)).relative_to(os.getcwd()))
    return str(path)


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


# --- the counting path and the cleanliness block -------------------------------------------------------
# One implementation each, imported by every writer. Before 2026-10-02 `code_cleanliness` was
# copy-pasted into ten files and the closure behind it into eight, and the copies had drifted into
# three algorithms that gave three different answers about the same tree: the same entry script came
# to 41, 68 or 44 paths depending on which copy ran. `scripts/cell2_eligibility.py` had never taken
# the cd263bc fix and published `genomeos/genome/Genome.py`, `genomeos/genome/Annotation.py` and
# `genomeos/knowledge/Reactome.py`, which are class names that no file spells, admitted only because
# APFS is case-insensitive. Nothing checked that the copies agreed.
#
# The closure here is the broadest of the three, because a counting path is a claim about what COULD
# have entered a number and the superset is the honest one: it follows every import that resolves to a
# file in this repository rather than only the package, it counts the `__init__.py` of each package
# above a module, and it resolves bare imports against the directories the source literally puts on
# `sys.path` as well as the repository root. That last part is what makes
# `scripts/response_map_coverage.py` appear on `scripts/placement_cause_198.py`'s path: that script
# reaches it through `sys.path.insert(0, str(ROOT / "scripts"))`, its code does enter the numbers, and
# every earlier copy reported it as off the path and therefore harmless.

#: Paths under these roots are the repository's own code for the closure (manifest.CODE_ROOTS covers
#: what makes a stamp dirty; this is what a closure may follow into).
_CLOSURE_ROOTS = ("genomeos/", "scripts/", "tests/")


def _is_file_exactly(root: Path, parts: Sequence[str]) -> bool:
    """A file at `parts` below `root`, spelled as the directories spell it (carried from cd263bc).

    `Path.is_file()` alone is not enough: on a case-insensitive filesystem it answers True for
    `genomeos/genome/Genome.py`, so a class name enters the closure as a file and the result cannot be
    reproduced where the filesystem is case-sensitive.
    """
    node = root
    for part in parts:
        if not part or part in (".", ".."):
            return False
        try:
            if part not in {p.name for p in node.iterdir()}:
                return False
        except OSError:
            return False
        node = node / part
    return node.is_file()


def _literal_path(node: ast.AST, names: dict[str, PurePosixPath]) -> PurePosixPath | None:
    """A path expression's value as a repository-relative path, or None when it is not literal.

    Reads the expression rather than guessing from its text, so `Path(__file__).resolve().parent` and
    `Path(__file__).resolve().parents[1]` are told apart, as are `ROOT` and `ROOT / "scripts"`.
    """
    if isinstance(node, ast.Constant):
        return PurePosixPath(node.value) if isinstance(node.value, str) else None
    if isinstance(node, ast.Name):
        return names.get(node.id)
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        left = _literal_path(node.left, names)
        right = _literal_path(node.right, names)
        return None if left is None or right is None else left / right
    if isinstance(node, ast.Call):
        func = node.func
        name = func.attr if isinstance(func, ast.Attribute) else getattr(func, "id", "")
        if name in ("resolve", "absolute", "expanduser") and isinstance(func, ast.Attribute):
            return _literal_path(func.value, names)
        if name in ("str", "Path", "PurePath", "PurePosixPath", "fspath") and node.args:
            return _literal_path(node.args[0], names)
        return None
    if isinstance(node, ast.Attribute):
        if node.attr == "parent":
            base = _literal_path(node.value, names)
            return None if base is None else base.parent
        return None
    if isinstance(node, ast.Subscript):  # X.parents[n]
        value = node.value
        if (
            isinstance(value, ast.Attribute)
            and value.attr == "parents"
            and isinstance(node.slice, ast.Constant)
            and isinstance(node.slice.value, int)
        ):
            base = _literal_path(value.value, names)
            if base is None:
                return None
            for _ in range(node.slice.value + 1):  # parents[0] is the directory holding the file
                base = base.parent
            return base
        return None
    return None


def _module_level_names(tree: ast.AST, file_rel: str) -> dict[str, PurePosixPath]:
    """`__file__` and every module-level `NAME = <path expression>` this file binds, as paths."""
    names: dict[str, PurePosixPath] = {"__file__": PurePosixPath(file_rel)}
    for node in getattr(tree, "body", []):
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        if not isinstance(target, ast.Name):
            continue
        value = _literal_path(node.value, names)
        if value is not None:
            names[target.id] = value
    return names


def _search_dirs(tree: ast.AST, file_rel: str) -> list[str]:
    """The repository directories this file puts on `sys.path`, read from its own insert/append calls.

    A target that is not literal (`sys.argv[1]`) is skipped: it names no directory at write time.
    """
    names = _module_level_names(tree, file_rel)
    out: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call) or not isinstance(node.func, ast.Attribute):
            continue
        if node.func.attr not in ("insert", "append"):
            continue
        holder = node.func.value
        if not (
            isinstance(holder, ast.Attribute)
            and holder.attr == "path"
            and isinstance(holder.value, ast.Name)
            and holder.value.id == "sys"
        ):
            continue
        args = [a for a in node.args if not isinstance(a, ast.Constant) or not isinstance(a.value, int)]
        if not args:
            continue
        value = _literal_path(args[-1], names)
        if value is None:
            continue
        rel = value.as_posix()
        out.append("" if rel in (".", "") else rel)
    return out


def _imported_modules(tree: ast.AST, file_rel: str) -> list[str]:
    """Every dotted module name one file imports, with relative imports resolved against its package.

    A `from X import a, b` contributes `X` and also `X.a` and `X.b`, because the name after `import`
    may itself be a submodule; the resolver keeps only the ones that are files in this repository.
    """
    package = list(PurePosixPath(file_rel).parts[:-1])
    out: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out.extend(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package[: len(package) - (node.level - 1)]
                prefix = ".".join([*base, *([node.module] if node.module else [])])
            else:
                prefix = node.module or ""
            if prefix:
                out.append(prefix)
                out.extend(f"{prefix}.{a.name}" for a in node.names)
    return out


def _module_files(module: str, root: Path, search: str) -> list[str]:
    """The repository files a dotted module name reads when it is imported from `search`.

    The module itself and the `__init__.py` of every package above it. Empty when the name resolves to
    nothing in this repository, which is how the standard library and the dependencies stay off a
    counting path.
    """
    base = list(PurePosixPath(search).parts) if search else []
    parts = [*base, *module.split(".")]
    found = [
        "/".join([*parts[: len(base) + i], "__init__.py"])
        for i in range(1, len(parts) - len(base))
        if _is_file_exactly(root, [*parts[: len(base) + i], "__init__.py"])
    ]
    for candidate in ([*parts[:-1], f"{parts[-1]}.py"], [*parts, "__init__.py"]):
        if _is_file_exactly(root, candidate):
            found.append("/".join(candidate))
            break
    return [f for f in found if f.startswith(_CLOSURE_ROOTS) or "/" not in f]


def counting_path(entry: str | Path, root: Path | None = None) -> list[str]:
    """Every file in this repository that `entry` can read by import, directly or at any remove.

    The one implementation: `scripts/` and `genomeos/` must not define their own
    (tests/test_code_cleanliness_shared.py fails if one does), because eight copies of this closure had
    drifted into three algorithms that disagreed about the same tree.

    Found by parsing the import statements, never by importing, so the answer is the same in any
    checkout of the same revision -- which is what lets a rebuild in a clean worktree reproduce it.
    Resolution follows the repository root and every directory the files literally put on `sys.path`,
    counts the `__init__.py` of each package above a module, and matches every path component against
    what its directory actually lists, so a class name cannot enter as a file.
    """
    root = Path(root) if root is not None else Path.cwd()
    first = _repo_relative_to(entry, root)
    seen = {first}
    dirs = {""}
    modules: set[str] = set()
    trees: dict[str, ast.AST | None] = {}

    def tree_of(rel: str) -> ast.AST | None:
        if rel not in trees:
            try:
                trees[rel] = ast.parse((root / rel).read_text())
            except (OSError, SyntaxError, ValueError):
                trees[rel] = None
        return trees[rel]

    changed = True
    while changed:
        changed = False
        for rel in sorted(seen):
            tree = tree_of(rel)
            if tree is None:
                continue
            for d in _search_dirs(tree, rel):
                if d not in dirs:
                    dirs.add(d)
                    changed = True
            for m in _imported_modules(tree, rel):
                if m not in modules:
                    modules.add(m)
                    changed = True
        for m in sorted(modules):
            for d in sorted(dirs):
                for f in _module_files(m, root, d):
                    if f not in seen:
                        seen.add(f)
                        changed = True
    return sorted(seen)


def _repo_relative_to(path: str | Path, root: Path) -> str:
    """`path` as a repository-relative posix string, whether it arrives absolute or already relative."""
    p = Path(path)
    if not p.is_absolute():
        return PurePosixPath(p.as_posix()).as_posix()
    try:
        return p.resolve().relative_to(root.resolve()).as_posix()
    except (OSError, ValueError):
        return p.as_posix()


def code_cleanliness(entry: str | Path, own_code: Iterable[str], root: Path | None = None) -> dict[str, Any]:
    """Which uncommitted code this shared checkout held when a result was written, split into the
    writing lane's own and other lanes', and whether any of it is on the counting path.

    The one implementation (CLEANLINESS_KEYS is its key set, `genomeos/results.py` refuses a registry
    write without it). Both of the things a caller may legitimately vary are arguments and neither is
    guessed: `entry` is the script whose import closure is the counting path, and `own_code` is the
    paths the writing lane is answerable for. Read from git and from the imports, never asserted.

    The field names are the convention `scripts/cell2_eligibility.py` set and
    `scripts/manifest_rebuild.py` matches, so `own_code_is_committed` and
    `foreign_uncommitted_code_on_the_counting_path` must hold on both sides of a rebuild while the
    lists that merely describe the tree a run happened in may take their clean-worktree values.
    """
    own_code = tuple(own_code)
    rev = code_revision(root)
    path = counting_path(entry, root)
    dirty = list(rev.get("dirty_code_paths") or [])
    own = [p for p in dirty if p in own_code]
    foreign = [p for p in dirty if p not in own_code]
    return {
        "git_sha": rev.get("git_sha"),
        "dirty": rev.get("dirty"),
        "own_uncommitted_code": own,
        "own_code_is_committed": not own,
        "foreign_uncommitted_code": foreign,
        "foreign_uncommitted_code_on_the_counting_path": [p for p in foreign if p in path],
        "counting_path": path,
        "counting_path_count": len(path),
        "counting_path_is_computed": (
            "the transitive import closure of this script, computed from the files' import statements "
            "at write time (genomeos.manifest.counting_path), never by importing: it follows every "
            "import that resolves to a file in this repository rather than only the package, counts "
            "the __init__.py of each package above a module, resolves bare imports against the "
            "directories the source literally puts on sys.path as well as the repository root, and "
            "matches every path component against what its directory lists, so a class name such as "
            "genomeos.genome.Genome cannot enter it as a file; not a hand-written list"
        ),
        "note": (
            "several lanes work in this one checkout. A file under foreign_uncommitted_code belongs to "
            "another lane; this lane did not write it and did not commit it. The counting path is the "
            "computed closure above, so a foreign file outside it cannot have entered a number here, "
            "and foreign_uncommitted_code_on_the_counting_path names any that could"
        ),
    }


#: The keys `code_cleanliness` returns. A result entering the registry under a name that is not on the
#: legacy allowlist must carry all of them (genomeos/results.py): the block is what says which code was
#: uncommitted when the numbers were written, and before 2026-10-02 it was copy-pasted into ten files
#: whose answers had drifted apart, so a hand-built block is refused rather than trusted.
CLEANLINESS_KEYS = (
    "git_sha",
    "dirty",
    "own_uncommitted_code",
    "own_code_is_committed",
    "foreign_uncommitted_code",
    "foreign_uncommitted_code_on_the_counting_path",
    "counting_path",
    "counting_path_count",
    "counting_path_is_computed",
    "note",
)


def cleanliness_problems(manifest: Any) -> list[str]:
    """What a manifest's `code_cleanliness` lacks, as sentences; empty when `code_cleanliness` wrote it.

    Checked by key set, because that is what can be checked: a block carrying every key of
    CLEANLINESS_KEYS is the shape the shared function returns, and one that lacks a key was built by
    hand and cannot be relied on to have counted the same things.
    """
    if not isinstance(manifest, dict) or "code_cleanliness" not in manifest:
        return ["missing code_cleanliness (call genomeos.manifest.code_cleanliness)"]
    block = manifest["code_cleanliness"]
    if not isinstance(block, dict):
        return ["code_cleanliness must be the dict genomeos.manifest.code_cleanliness returns"]
    missing = [k for k in CLEANLINESS_KEYS if k not in block]
    if missing:
        return [
            f"code_cleanliness lacks {', '.join(missing)}, so it did not come from "
            f"genomeos.manifest.code_cleanliness"
        ]
    return []


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
