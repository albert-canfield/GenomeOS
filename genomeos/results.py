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

    Since 2026-10-02 (lane-entrypoint) EVERY write into the registry, legacy name or not, is refused
    when the ENTRY SCRIPT is not committed: `sys.argv[0]` outside the repository, untracked, modified,
    `-c`, stdin, a test runner, or any form that cannot be classified (genomeos/manifest.py:
    ENTRY_FORMS, entry_script_problems). A cleanliness block carrying no `entry_script` at all is
    refused on the same footing, because the 24 `organised_chr*.json` results written by a script in
    /private/tmp all lacked the field while asserting `own_code_is_committed: True`
    (working tree at d4d0449, 19:41; superseded by e1def2e, 20:41), and nothing in them
    names the code that made the bytes. Outside the registry the fact is RECORDED and not refused,
    even under `strict=True`: no published result is strict-but-not-registry, and refusing there would
    cost every writer in the project its success-path test, since under pytest `sys.argv[0]` is the
    runner. `own_code_is_committed` is false for such a write either way, so nothing is certified.

    Since 2026-10-02 (lane-tracer) the declared inputs are also reconciled against what the writer
    actually opened: an audit hook records every file read under data/ (genomeos/manifest.py, "inputs
    are traced, not declared"), and a name that is not on the allowlist is refused when its traced
    reads include a file its `inputs` neither name nor contain. The refusal names every undeclared
    path. It exists because a declared input could be a POINTER at the bytes a run really read -- the
    24 `enhancer_targets_all_chr*` summaries declare 39,683 bytes of pointer files and read about
    1.4 GB of tables named by a field inside a result -- and because a rebuild from a clean checkout
    does not notice: the linked tables are found, the run completes, and the rebuild reports success
    without having hashed them. A legacy name warns and still writes, as it does for the rest of the
    contract. Whether a tracer was watching at all is recorded in `result_manifest.traced_inputs`, so
    an unwatched result is distinguishable from a clean one rather than reading the same.

    Since 2026-10-02 every write also records BOTH revision stamps and a flag when they disagree, as
    `result_manifest.revision_stamps` (genomeos/manifest.py: revision_stamps). The writer's cleanliness
    block reads HEAD while it computes and `stamp` reads it again here, so in a shared checkout a peer's
    commit between the two leaves a result naming two revisions. The disagreement is recorded, never
    silenced and never a reason to refuse the write: `scripts/manifest_rebuild.py` names the race and
    decides it with `git diff` over the counting path.

    `compact=True` writes the JSON without whitespace (`separators=(",", ":")`), for the per-cell and
    per-site tables that were written compact before they came under the contract (item 12 S6
    follow-up, lane-contract); the values are the same either way.
    """
    # The window closes before this function reads anything of its own, so what is reconciled is the
    # run that made this result and nothing else.
    held = mf.trace_close()
    p = results_dir / f"{name}.json"
    body = dict(payload)
    given = manifest if manifest is not None else body.pop(mf.KEY, None)
    body.pop(mf.KEY, None)
    registry = _is_registry(results_dir)
    legacy = registry and name in legacy_names()
    # The next result's window starts here, not above, because reading the legacy allowlist is this
    # function's own read and not the writer's: charged to the next result it would refuse a chain
    # writer's second result for a file the enforcement machinery opened. Nothing under data/ is read
    # after this point -- the revision stamp reads git and the cleanliness block reads code -- so the
    # window the next result is answerable for holds that result's reads and nothing else. It is opened
    # here rather than in a `finally` so that an unexpected failure below cannot leave it closed.
    mf.trace_begin()
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
    # --- the entry script must be committed (2026-10-02, lane-entrypoint) -------------------------
    # Measured defect: `data/results/organised_chr21.json` and 23 siblings were written at 19:41-19:42
    # by `/private/tmp/.../scratchpad/write.py`, a script OUTSIDE this repository, and every one of
    # them asserted `own_code_is_committed: True`
    # (working tree at d4d0449, 19:41; superseded by e1def2e, 20:41): the bytes are in no commit either
    # side, the cause is recorded in `2ebe9df`, and the scratchpad
    # script is the live artefact). The cleanliness block compared git's dirty list
    # against the paths the caller DECLARED and never looked at the script that was running, so
    # nothing in the manifest described the code that produced the bytes and no rebuild could run.
    # `mf.entry_script_problems` reads the `entry_script` sub-block the shared cleanliness function
    # now fills from `sys.argv[0]`, and the ABSENCE of that sub-block is refused too: the 24 results
    # all lack it, and reading its absence as a pass would be the same defect with an extra step.
    # The refusal takes the path that was already here: quarantined with its reason, ManifestError
    # naming where it went, nothing written to results_dir.
    #
    # WHERE IT FIRES, and this is a deliberate departure from the brief's wording ("a strict
    # save_result refuses"), measured rather than reasoned about. The condition is `registry`, NOT
    # `enforce`:
    #
    #   * EVERY write into the result registry is checked, legacy name or not. That is STRONGER than
    #     `enforce`, which for a legacy registry name is `bool(strict)` -- and `organised_chr*`, the
    #     24 names this check exists for, are legacy. A registry write is what publishing is.
    #   * A write outside the registry -- a test's `tmp_path`, a scratch directory -- records the fact
    #     and does not refuse, even under `strict=True`. No published result is ever
    #     strict-but-not-registry, so nothing that gets published escapes; what it buys is that a
    #     writer can still be driven from a test. Measured when `enforce` was the condition: nine
    #     tests across tests/test_traced_inputs.py, tests/test_organise_inputs.py and
    #     tests/test_attribution.py went red, all of them positive controls that run a real writer
    #     into a temporary directory, because under pytest `sys.argv[0]` is the runner. Refusing
    #     there costs every writer in the project its success-path test and protects nothing.
    #
    # NOTHING IS CERTIFIED EITHER WAY: `code_cleanliness` records `entry_script.form` as `test_runner`
    # and `own_code_is_committed` as FALSE for such a write, so a scratch result cannot be read as
    # having come from committed code, and `scripts/check_staged.py` refuses to let one be committed.
    # The exemption has its counterfactual, not an argument: tests/test_entry_script_counted.py
    # plants the same write twice, into the registry and into a scratch directory, and the first must
    # refuse while the second must still record the false verdict.
    entry_refusal = False
    if registry:
        unentered = mf.entry_script_problems(stamped)
        if unentered:
            stamped["problems"] = [*stamped.get("problems", []), *unentered]
            stamped["complete"] = False
            # A legacy name only WARNS about the rest of the contract; this one refuses, because the
            # names written by a script outside the repository are legacy names. Carried in its own
            # flag rather than by raising `enforce`, so that a legacy name's traced-input mismatch
            # below keeps warning exactly as it did -- one rule's strictness must not quietly become
            # another's.
            entry_refusal = True
    # --- end the entry script must be committed ---------------------------------------------------
    # --- both revision stamps, with a flag when they disagree (2026-10-02) ------------------------
    # A result carries HEAD twice: once from the writer's cleanliness block, once from the stamp above.
    # In this shared checkout a peer can commit between the two, and then they name two different
    # commits -- measured on astroreg2_astrocyte_activity (53f3b33 against 675e54a). Recorded here so
    # the reader is told, and so a rebuild can name the race and run the decisive diff instead of
    # inheriting an unexplained difference. Additive, and never a reason to refuse the write.
    stamped[mf.REVISION_STAMPS] = mf.revision_stamps(stamped)
    # --- traced-input reconciliation (2026-10-02, lane-tracer) ------------------------------------
    # A declared input is what the writer said it read; `traced_inputs` is what it opened. A new name
    # whose reads exceed its declared inputs is refused on the path that was already here: quarantined
    # with its reason, ManifestError naming where it went, nothing written to results_dir. A legacy
    # name warns and still writes. The block is attached either way, including when the reads were all
    # declared, so a result that was never traced does not read like one that reconciled cleanly.
    stamped["traced_inputs"] = mf.traced_inputs(stamped, held)
    mismatch = mf.reconciliation_problems(stamped["traced_inputs"])
    if mismatch and enforce:
        stamped["problems"] = [*stamped.get("problems", []), *mismatch]
        stamped["complete"] = False
    elif mismatch and legacy:
        warnings.warn(f"{p}: {'; '.join(mismatch)}", ManifestWarning, stacklevel=2)
    # --- end traced-input reconciliation ----------------------------------------------------------
    out = {"result": name, "date": time.strftime("%Y-%m-%d"), **body, mf.KEY: stamped}
    q = quarantine_dir(results_dir) / f"{name}.json"
    if not stamped["complete"] and (enforce or entry_refusal):
        _quarantine(q, out, p, stamped["problems"], registry)
        why = (
            f"the entry script must be committed for any write into the registry ({results_dir})"
            if entry_refusal
            else "a name not on the legacy allowlist must carry the contract"
            if registry
            else "strict"
        )
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
    """Whether `results_dir` is the result registry. A path that cannot be classified counts AS the
    registry, which is the strict answer and since 2026-10-02 the only safe one.

    Until then this returned False on an OSError, and the flag only chose between severities of
    warning, so failing open cost little. It now gates FOUR decisions in `save_result`, all of which
    go off together on that one exception: `legacy` (line 110), `enforce` (118, where a registry write
    stops being enforced at all and becomes a function of whatever the caller passed as `strict` -- a
    guard keyed to something the caller supplies is not a guard), the requirement to carry a
    cleanliness block at all (125), and the entry-script refusal (165). So an unclassifiable path is
    treated as the registry: a false-loud refusal gets read, a false-silent pass does not, and the
    refusal names the path so the reader can see which one could not be resolved.

    The try/except stays, because `resolve` can raise -- ELOOP on a symlink loop, ENAMETOOLONG -- and
    crashing every writer is worse than refusing one. tests/test_entry_script_counted.py plants it:
    `resolve` made to raise on a registry write must refuse, and a scratch directory that resolves
    normally must still write while recording the false verdict, so the fix cannot quietly have broken
    every test writer instead.
    """
    try:
        return results_dir.resolve() == RESULTS_DIR.resolve()
    except OSError:
        return True


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
