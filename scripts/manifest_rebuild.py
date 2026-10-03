# SPDX-License-Identifier: AGPL-3.0-or-later
"""Rebuild a result from its manifest in a clean second checkout, and say what differs or is missing.

    uv run python scripts/manifest_rebuild.py RESULT.json --where DIR [--venv fresh|shared] [--keep]

Review item R9's acceptance: a second environment reconstructs a representative result from its
manifest, or names exactly which dependency is unavailable. The steps, each of which can stop the
rebuild with the reason:

1. read `result_manifest` from RESULT.json (a result without one cannot be rebuilt: said so);
2. `git worktree add --detach` the recorded `code.git_sha` under DIR; uncommitted code at write time
   (`code.dirty_code_paths`) is named, since the commit alone did not run;
3. link the local data stores (data/reference, data/knowledge, data/cache) read-only, as the
   pre-push hook does, link the machine-local inputs under data/results read-only as well (see
   `link_machine_local_inputs`), and check every input's sha256 against the manifest, naming any input
   that is absent or has different bytes. Every declared input is accounted for: a group of files read
   together is opened and hashed member by member, and an input the tool cannot resolve stops the
   rebuild by name. `inputs_declared`, `inputs_checked` and `inputs_unchecked` carry the denominator,
   so an input that was not opened cannot leave the list and read as one that matched;
   `inputs_satisfied_from_machine_local_paths` says how much of the rebuild rested on this machine.
   Since 2026-10-02 the three symlinked stores are also protected by an audit hook in this process AND in
   the rebuilt command's own (scripts/rebuild_write_guard.py, `write_guard` in the report): a symlink
   carries its target's mode, so until then a rebuild could write to the only copy of an input store on
   this machine, which was measured and not supposed. The hook is a Python-level guard and the report
   says so; a clone at 0444 is the kernel-enforced one. A rebuild whose command did not arm it reports no
   verdict. `copied_not_cloned` names any machine-local input whose `cp -c` fell back to a real copy.
   Since 2026-10-03 an input declared as an ABSOLUTE path that RESOLVES INSIDE the repository root is
   read by its repository-relative form and checked like any other (`repo_relative_declaration`,
   `relativised_absolute_declarations`): `genomeos.manifest.input_entry` records `str(path)` where its
   sibling `files_entry` relativises, which left 9 of the 212 committed results with inputs declared
   absolutely and unrebuildable. An absolute path that resolves OUTSIDE the root is refused exactly as
   before, and the refusal now names the root it was tested against;
4. run `code.argv` in the worktree, in a fresh environment from the committed uv.lock (`--venv
   fresh`, offline) or the checkout's own (`--venv shared`);
5. compare the rebuilt result with RESULT.json field by field, ignoring only `date` and the
   `code` block of the manifest (which records the run, not the result).
   Since 2026-09-28 wall-clock timing keys are ignored too, at any depth (`seconds`, `*_seconds`,
   `seconds_*`, `*_per_second`); the report lists their paths under `timing_fields_ignored`.
   Since 2026-10-02 the run-resource readings named in `RESOURCE_KEYS` are ignored too, by exact key,
   and listed under `resource_fields_ignored`.
   Since 2026-10-02 the tracer's own `opens` count is ignored, by the exact paths in `RUN_FIELDS` and
   never by key, suffix or substring, and listed under `run_fields_ignored`; `cwd_at_open` and
   `cwd_at_close` are compared, and set aside as environment only where the recorded value is an
   absolute path (`CWD_FIELDS`, `cwd_fields_recorded_absolute`).
   Since 2026-10-02 neither name test applies to a key the DATA supplied rather than the writer: a
   timing or resource name reached through a key tally (`is_key_tally`: a `keys` or `vocabulary` token
   in the ancestry) is COMPARED, and those leaves are listed under
   `names_not_exempted_because_the_data_supplied_them`. The matcher had been setting aside three
   key-FREQUENCY counts (963,406, 3,209 and 705 under `schema_keys/seconds`) as though they were
   durations. Every leaf a name test does set aside is printed with its VALUE under
   `set_aside_by_a_name_test`.
   Whatever was set aside is printed with the reason (`ignored_because`), and `comparison_reading` is
   the one sentence that says whether anything was: "0 differences" and "0 differences after setting
   aside N leaves" are different claims and the second is never printable as the first.

Two holes in step 4's read reconciliation were measured by `lane-silencerebuild` on 2026-10-03 and are
closed in `reads_against_declared`: a RELATIVE read was resolved against the TOOL'S working directory
and then dropped as not under `data/` (26 of silenceragree's 3,340 recorded reads, silently, two of
which the manifest declared and which were then reported never opened), and a run that read back ITS OWN
OUTPUT was refused a verdict before a single leaf was compared, though the worktree is checked to have
held nothing at that path before the command ran.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from genomeos import manifest as mf

#: The write guard, loaded by path rather than imported: this script is itself loaded by path in the
#: tests, so a bare `import rebuild_write_guard` would not resolve, and the worktree's own `genomeos`
#: package is at an older revision and may predate the guard entirely.
_GUARD_FILE = Path(__file__).resolve().parent / "rebuild_write_guard.py"
_guard_spec = importlib.util.spec_from_file_location("rebuild_write_guard", _GUARD_FILE)
guard = importlib.util.module_from_spec(_guard_spec)
_guard_spec.loader.exec_module(guard)

STORES = ("reference", "knowledge", "cache")
IGNORED = ("date",)

#: Where a result's machine-local inputs live. data/results is git-ignored as a rule (docs/DATA.md), so a
#: result computed from other results declares inputs that no clean worktree can hold: reader v1's cached
#: DNase peak sets, 13 biosamples by 23 chromosomes, are on this machine and in no commit. Until 2026-10-02
#: the rebuild linked only the three data stores, so those inputs were reported absent one by one and the
#: rebuild stopped before any comparison: response_map_increment2 reached 111 of 410 inputs opened, and
#: context_evidence hit the same wall. They are linked here instead, read-only, and hashed like any other.
MACHINE_LOCAL_DIR = "data/results"

#: What an input's `sha256` says when nothing was hashed to produce it. Such an input cannot be
#: checked, which is a reason to stop and say so, not a reason to pass.
NOT_HASHED = "n/a"

#: The fields a rebuild in a clean worktree is allowed to differ on, by **exact path**, because they
#: describe the tree the run happened in and not anything the run computed. A result that honestly records
#: the other lanes' outstanding files in a shared checkout cannot reproduce those lists in a clean worktree,
#: where there are none; under the unamended rule such a result could only ever fail, while one that omitted
#: the record would pass. Adopted 2026-10-02 after the coordinator raised it and the supervisor ruled.
ENVIRONMENT_FIELDS = tuple(
    f"/{container}{field}"
    for container in ("code_cleanliness/", "result_manifest/code_cleanliness/", "result_manifest/code/")
    for field in (
        "dirty",
        "foreign_uncommitted_code",
        "dirty_code_paths",
        "untracked_code_paths",
        "dirty_result_paths",
    )
)

#: These must hold on **both** sides, so the exemption above can never excuse a result whose own code was
#: uncommitted, or one with a foreign uncommitted file on the path that produced its numbers.
MUST_HOLD: tuple[tuple[str, Any], ...] = (
    ("own_code_is_committed", True),
    ("foreign_uncommitted_code_on_the_counting_path", []),
)


# --- the tracer's own fields (2026-10-02, lane-tracercwd) -----------------------------------------
#
# Two separate lanes finished a rebuild that reconciled at every leaf and then had to inspect and argue
# away the same three differences by hand before they could report. lane-finemap's: 1,010 of 1,010 inputs
# checked, 0 absent, 0 differing bytes, 0 must_hold failures, 4,353 of 4,370 leaves reconciling -- and
# `opens` 1,684 against 1,855, plus two absolute working directories that a second worktree cannot ever
# hold. A difference list that is never empty teaches its reader to skim it, which is the opposite of
# what it is for.
#
# The two halves are not the same kind of fix, and conflating them was the first draft's error. A cwd
# recorded absolutely carries NO information outside its own checkout, so making it repository-relative
# loses nothing. `opens` DOES carry information -- see RUN_FIELDS, where what exempting it costs is
# written down -- so it is set aside as a known loss with the loss named in the report, and not as noise.
#
# `cwd_at_open` and `cwd_at_close` are NOT exempted. `genomeos.manifest.cwd_as_recorded` now writes them
# relative to the repository root, so the field keeps the purpose manifest.py states for it -- which
# directory the recorded relative paths were taken against, and whether it moved inside the window --
# and two checkouts at different absolute paths record the same value. Only the one case that cannot be
# made relative, a cwd outside any repository, is recorded absolute, and only an absolute value is set
# aside here. A result written BEFORE that change carries an absolute path and is covered by the same
# clause: its bytes are not rewritten, because a committed result's bytes are a pin.

#: The two leaves that say which directory a run's relative paths were taken against, by **exact path**.
#: `traced_inputs` is a dict at one place in a manifest, not a list, so these paths carry no index.
CWD_FIELDS = (
    "/result_manifest/traced_inputs/cwd_at_open",
    "/result_manifest/traced_inputs/cwd_at_close",
)

#: Exempt from comparison by **exact path**. This exemption is KNOWN TO BE LOSSY and the report says so:
#: it is not the case that `opens` is noise.
#:
#: What is measured, and it is one case: in lane-finemap's rebuild `opens` moved 1,684 to 1,855 because
#: THE WRITER'S OWN BEHAVIOUR CHANGED between the two runs -- its carry-forward had begun reading a git
#: blob in place of a file on disk. That is a real change in what the code opened. No claim is made here
#: that the count varies with the interpreter or the installed libraries; that was an unverified premise,
#: it was withdrawn by the lane that raised it, and nothing in this file rests on it. What is claimed is
#: only that the count is too coarse to diagnose anything -- it is one integer over every `open` the
#: process made, naming no path -- and that a rebuild cannot act on it.
#:
#: So what is given up by setting it aside, stated rather than left to be inferred: `canonical` returns
#: None for any path not under data/, and the hook then increments `opens` and records NOTHING ELSE. So
#: `opens` is the only leaf in a manifest with any sensitivity at all to what a writer opened OUTSIDE
#: data/, and the measured case is exactly that -- a git blob is not under data/. Exempting it makes the
#: comparison blind to a change in a writer's reads outside data/.
#:
#: What is NOT given up, and what makes the trade defensible: everything the tracer records about reads
#: UNDER data/ is still compared, leaf by leaf -- `opens_under_data`, `files_read`, `files_written`,
#: `active`, and the two path lists `undeclared` and `declared_not_read`. A writer that becomes an
#: undeclared input of itself, or that stops reading a file it declares, moves one of those.
#:
#: Why not compare the SET of traced paths instead of the count, which would be strictly better: the set
#: is not in a manifest to compare. `manifest.traced_inputs` keeps the counts and the two difference
#: lists and drops `trace_close()["reads"]`, so there is no set in any of the ~1,100 existing results;
#: and the set that WOULD be recorded is canonicalised through data/, so it would not have held the git
#: blob either and would not have caught the measured case. Recording something reproducible about reads
#: outside data/ is recorder work and a separate lane, not a comparison change.
#:
#: Exact paths, never a key, a suffix or a substring match, and this is not hypothetical:
#: `data/results/response_map_increment3.json` records `/per_element_response_cache/opens: 0`, which is
#: an ASSERTION THE RUN COMPUTED -- "the run opened no file under the per-element cache", recorded by a
#: patched `builtins.open` precisely because no reading of the source could establish it. A suffix match
#: on `opens` would set that leaf aside, and a rebuild in which the run did open the cache would pass.
RUN_FIELDS = ("/result_manifest/traced_inputs/opens",)


def _at_path(payload: Any, path: str) -> Any:
    """The value a `/a/b` path names, or None when nothing is there. Keys only: the paths that use this
    carry no list index, and a path that needed one would not be found rather than silently matched."""
    cur = payload
    for part in path.strip("/").split("/"):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def cwd_fields_recorded_absolute(original: Any, rebuilt: Any) -> list[str]:
    """Those of CWD_FIELDS whose recorded value is an ABSOLUTE path on either side, and only those.

    This is the one cwd case a rebuild cannot be asked to reproduce, since a second worktree is at a
    different absolute path by construction. Conditioned on the VALUE and not on the path, so the same
    two leaves are compared like any other when both sides record a repository-relative directory: `.`
    against `scripts` is a real difference and fails.
    """
    out = []
    for path in CWD_FIELDS:
        values = [_at_path(x, path) for x in (original, rebuilt)]
        if any(isinstance(v, str) and os.path.isabs(v) for v in values):
            out.append(path)
    return out


def comparison_reading(differences: list[str], ignored: dict[str, list[str]]) -> str:
    """The one sentence that says what was compared and what was set aside.

    "0 differences" and "0 differences after setting aside 3 leaves" are different claims, and this is
    where the second is kept from being printable as the first.
    """
    set_aside = {k: len(v) for k, v in ignored.items() if v}
    total = sum(set_aside.values())
    named = ", ".join(
        f"{n} {k.removesuffix('_ignored').replace('_', ' ')}" for k, n in sorted(set_aside.items())
    )
    # the one set-aside class that is a known loss rather than a free one has to reach the one sentence a
    # reader reads, or the sentence invites exactly the inference the loss has to be protected from.
    lossy = (
        " One class set aside is KNOWN-LOSSY: run fields (the tracer's total `opens`), whose exemption"
        " will hide a genuine change in the set of files this writer opens outside data/."
        if ignored.get("run_fields_ignored")
        else ""
    )
    if differences:
        tail = f"; {total} further leaves were set aside ({named})" if total else "; nothing was set aside"
        return f"{len(differences)} differences{tail}.{lossy}"
    if total:
        return (
            f"0 differences AFTER SETTING ASIDE {total} leaves ({named}); every set-aside path is listed "
            f"above with its reason, and this is not the same claim as 0 differences.{lossy}"
        )
    return "0 differences, and nothing was set aside: every leaf compared reconciles"


def no_verdict_reason(declared: int, checked: int) -> str | None:
    """Why no verdict may be reported over these inputs, or None when one may.

    The one rule the tool must never be able to break: a clean verdict over inputs it never opened. The
    project's lesson of 2026-10-01 is that a verification tool can fail in the direction that flatters, and
    this is where that would happen -- a report that compared the fields of two files and said nothing about
    the bytes behind them reads as a reproduction. So the decision lives in one function, is asked for both
    before the run and again before the verdict is written, and says yes only when at least one input was
    declared and every declared input was opened and hashed.
    """
    if declared <= 0:
        return (
            "the manifest declares no inputs, so no bytes were opened: a comparison of fields would say "
            "nothing about what the run read, and no verdict is reported"
        )
    if checked != declared:
        return f"{checked} of {declared} declared inputs were opened and hashed: no verdict is reported"
    return None


def leaves(x: Any, where: str = "") -> int:
    """How many leaf values a comparison covers. A count of top-level keys can hide a nested difference:
    `cell2_eligibility.json` has 14 top-level keys and several hundred leaves."""
    if isinstance(x, dict):
        return sum(leaves(v, f"{where}/{k}") for k, v in x.items())
    if isinstance(x, list):
        return sum(leaves(v, f"{where}[{i}]") for i, v in enumerate(x))
    return 1


def leaf_reconciliation(original: dict[str, Any]) -> dict[str, Any]:
    """Where every leaf of the original went, so a count carries its denominator.

    `compared` plus `not_compared` equals `total`: a comparison reported without the leaves it dropped
    invites the reader to assume it covered the file. The dropped leaves are the ones `comparable()`
    removes before any diff -- the date, the manifest's own `code` block, `model_dependencies` and every
    timing key -- and they are named here by cause rather than merely counted.
    """
    total = leaves(original)
    compared = leaves(comparable(original))
    return {
        "total": total,
        "compared": compared,
        "not_compared": total - compared,
        "not_compared_because": {
            "ignored_keys": list(IGNORED),
            "manifest_code_block": "result_manifest.code: the stamp of the run, not a value it computed",
            "model_dependencies": "recorded after the fact for results made before 2026-09-28",
            "timing_keys": (
                "every key is_timing() matches OUTSIDE a key tally, listed under timing_fields_ignored "
                "and printed with its value under set_aside_by_a_name_test; a timing name the data "
                "supplied is compared, and named under "
                "names_not_exempted_because_the_data_supplied_them"
            ),
            "resource_keys": (
                "the exact keys in RESOURCE_KEYS, a run-resource reading rather than a value the run "
                "computed, listed under resource_fields_ignored"
            ),
        },
        "reconciles": total == compared + (total - compared),
    }


def _base_path(difference: str) -> str:
    """The exact path a `diff()` line is about, with any list index dropped from its last step."""
    return difference.split(":")[0].split("[")[0]


def environment_differences(
    found: list[str], also_environment: tuple[str, ...] | list[str] = ()
) -> tuple[list[str], list[str]]:
    """(real differences, environment differences). A path counts as environmental only if it is exactly
    one of ENVIRONMENT_FIELDS or of `also_environment`, or an element of one of those lists; no pattern
    matching. `also_environment` is for an exemption that depends on the recorded VALUE rather than on
    the path alone -- `cwd_fields_recorded_absolute` is the only one -- so the caller decides per rebuild
    and the decision is still an exact path when it reaches here."""
    exempt = (*ENVIRONMENT_FIELDS, *also_environment)
    real, env = [], []
    for d in found:
        (env if _base_path(d) in exempt else real).append(d)
    return real, env


def run_differences(found: list[str]) -> tuple[list[str], list[str]]:
    """(everything else, the differences at a RUN_FIELDS path). Exact path, never a key or a suffix: an
    `opens` at any other path is a quantity some run computed and is compared like any other leaf."""
    rest, run = [], []
    for d in found:
        (run if _base_path(d) in RUN_FIELDS else rest).append(d)
    return rest, run


# --- the revision race, named and decided (2026-10-02, lane-raceandlocal) -------------------------
#
# A result records HEAD twice: `code_cleanliness.git_sha`, sampled while the writer computes, and
# `code.git_sha`, stamped when `save_result` writes. In the checkout several lanes share, HEAD moves
# between the two -- a peer commits -- and the result then names two commits. This rebuild builds its
# worktree at `code.git_sha`, so the cleanliness sha it recomputes CANNOT equal the recorded one, and the
# difference arrived in every such report with no explanation.
#
# An exemption for that leaf was proposed and REFUSED. Setting it aside would hide the one thing it
# carries that nothing else does -- that HEAD moved during the write, which makes the result's provenance
# genuinely ambiguous -- and the standing rule is that where a classification is uncertain the tool errs
# toward COMPARING (docs/LESSONS.md, 2026-10-02). So the leaf stays compared, and the race is made
# DECIDABLE instead: named, then settled by `git diff A B` restricted to the counting path.
#
# Benign is not silent. An empty diff is printed and named exactly like a non-empty one; the difference
# is the verdict, and the verdict is computed by git rather than argued in a report.

#: The leaves a revision race moves, by **exact path**. NOT an exemption list: none of these is set
#: aside. They are named so a report can say which of its differences the race accounts for while still
#: reporting every one of them, which is the whole difference between explaining a leaf and excusing it.
#: Both spellings of the cleanliness container occur, since some results carry the block at the top level
#: as well as inside the manifest. `revision_stamps` is here because its four moving leaves restate the
#: one fact: a rebuild's worktree is clean and at `code.git_sha`, so it recomputes agreement every time.
RACE_FIELDS = (
    "/code_cleanliness/git_sha",
    "/result_manifest/code_cleanliness/git_sha",
    "/result_manifest/revision_stamps/code_cleanliness_git_sha",
    "/result_manifest/revision_stamps/agree",
    "/result_manifest/revision_stamps/revision_race",
    "/result_manifest/revision_stamps/reading",
)

#: How many characters of a sha a report's prose uses. The full shas are always printed beside it.
SHA_IN_PROSE = 7


# --- the tool's own cleanliness (2026-10-02, the dirty-tool blind spot) ---------------------------
#
# Measured by a peer: the write guard above sat UNCOMMITTED in the working copy, 232 insertions against
# HEAD, and a rebuild run in that state is judged by verification code nobody can reproduce. It is the
# instrument-side twin of the stale code stamp found in a result the same morning.
#
# MUST_HOLD cannot catch it, and that is the point of a separate mechanism:
# `foreign_uncommitted_code_on_the_counting_path` is read off the REBUILT result's manifest, which was
# composed inside the worktree and knows nothing about the main checkout the tool itself ran from.
# Nothing in the report said the tool was dirty.
#
# So the tool reports its own sha and the cleanliness of its OWN import closure, and refuses a verdict
# when any file of that closure is uncommitted. A refusal and not a field, for the reason the unarmed
# guard is a refusal: a verdict nobody can reproduce is worse than no verdict.

#: The entry points whose import closures make up the verification tool. `rebuild_write_guard.py` is
#: named EXPLICITLY because `counting_path` follows import STATEMENTS and this script loads the guard
#: through `spec_from_file_location`, so the computed closure would not contain it -- the recursion the
#: ruling warned about, and it really was missed by the import walk.
#:
#: The `sitecustomize.py` the parent generates for the child is NOT in the closure, and the choice is
#: deliberate: it is not a source file of this tool but an OUTPUT of one, and its every byte is
#: `guard.SITECUSTOMIZE`, which is in the closure through the guard. Pinning the generated copy as well
#: would pin the same bytes twice and would make a temporary directory part of a committed closure.
TOOL_ENTRIES = ("scripts/manifest_rebuild.py", "scripts/rebuild_write_guard.py")


def tool_cleanliness(root: Path) -> dict[str, Any]:
    """The verification tool's own revision and whether any file it runs on is uncommitted.

    Computed from the SHARED `genomeos.manifest.code_cleanliness`, the one implementation, with
    `own_code=()` so that every uncommitted file on the closure is reported rather than split into mine
    and someone else's: for the tool it does not matter whose edit it is, only that the code which
    produced the verdict is not in any commit.
    """
    closure: set[str] = set()
    for entry in TOOL_ENTRIES:
        closure.add(entry)
        closure.update(mf.counting_path(entry, root))
    block = mf.code_cleanliness(TOOL_ENTRIES[0], (), root)
    uncommitted = sorted(p for p in (mf.code_revision(root).get("dirty_code_paths") or []) if p in closure)
    return {
        "git_sha": block.get("git_sha"),
        "entries": list(TOOL_ENTRIES),
        "closure_count": len(closure),
        "closure_is_computed": (
            "the union of the import closures of TOOL_ENTRIES (genomeos.manifest.counting_path) plus the "
            "entries themselves. rebuild_write_guard.py is named rather than found, because it is loaded "
            "by path and no import statement points at it"
        ),
        "uncommitted_code_on_the_tool_s_own_closure": uncommitted,
        "tool_is_committed": not uncommitted,
        "generated_sitecustomize_excluded_because": (
            "it is an output of this tool, not a source of it; its bytes are guard.SITECUSTOMIZE, which "
            "is on the closure through the guard"
        ),
    }


# --- the reads the rebuild actually made (2026-10-02) ----------------------------------------------
#
# `manifest_rebuild` was BLIND BY CONSTRUCTION to any file a manifest did not name. Demonstrated end to
# end on a live result: `constrained_unknown_targets` rebuilt at 193 of 193 inputs declared and checked,
# 193 opened and hashed, 3,858 of 3,896 leaves compared, no must_hold failures, "0 differences AFTER
# SETTING ASIDE 1 leaves" -- while a tracer on the same run saw a 194TH read,
# `data/results/unknown_chr21.json`, pinned by no declared sha256. Two instruments disagreed about one
# run and the flattering one was the clean bill of health.
#
# So the rebuilt command's own reads are RECORDED (scripts/rebuild_write_guard.py, by the one audit hook
# that already refuses its writes) and reconciled against the declared inputs here.

#: Reads every post-fe0880a writer makes because the FRAMEWORK makes them, not because the result needs
#: them, exempt by **exact path** on the RUN_FIELDS pattern -- never a prefix, a suffix or a directory.
#:
#: `data/results_legacy.txt` is `genomeos/results.py: LEGACY_ALLOWLIST`, read by `save_result` itself on
#: every registry write. Under the unamended rule it would flag some 1,100 results at once.
#:
#: THE EXEMPTION IS CONDITIONAL AND THE CONDITION IS PROVEN, not assumed: it holds only because that
#: file can never change a result's VALUES, only whether `save_result` ADMITS it.
#: `tests/test_legacy_allowlist_exemption.py` changes the allowlist in a temporary registry and re-runs
#: a writer: the bytes come back identical, or the write is quarantined and nothing is written. If a
#: value could move it would be an input and would have to be declared.
#:
#: WHAT IS GIVEN UP, stated rather than left to be inferred: a change to the legacy allowlist cannot
#: alter a value, only admission. A rebuild is therefore blind to a result that was admitted because a
#: name was added to that list, which is a change in what the registry ACCEPTS and not in what it says.
FRAMEWORK_READS = ("data/results_legacy.txt",)


def reads_against_declared(
    worktree: Path,
    manifest: Any,
    recorded: list[str],
    stores: list[str],
    output: str | None = None,
    output_held_bytes_before_the_run: bool | None = None,
) -> dict[str, Any]:
    """Which recorded reads the manifest declares, and which it does not.

    Only reads under the rebuild's own `data/` are considered: the interpreter, the virtual environment
    and the standard library are not a result's inputs. A read through a linked store resolves OUTSIDE
    the worktree, so the store roots are passed in and matched too -- without them every store read
    would be dropped as "not under data/" and the reconciliation would be vacuous.

    A RELATIVE read is resolved against the directory the CHILD was launched in, not this tool's own.
    Until 2026-10-03 every recorded path went through `os.path.realpath` here, which resolves a relative
    path against the CWD OF THE TOOL, so a read the rebuilt command made as `data/results/x.json`
    resolved into whatever checkout the tool happened to be run from and was then dropped as "not under
    data/" -- silently, with nothing in the report saying a path had been discarded. Measured on
    silenceragree's rebuild: 26 of 3,340 recorded reads were relative, all 26 were missing from
    `reads_under_data`, and the remaining 3,283 absolute data reads equalled `reads_under_data_count`
    exactly. Two of those 26 were cCRE members the manifest declares, which the reconciliation then
    reported as DECLARED BUT NEVER OPENED -- a false report about the record, from the same class of
    hole as the 194th-file blindness this reconciliation was built to close.

    WHAT THE RESOLUTION STILL CANNOT DO, stated rather than left to be inferred: the child's working
    directory is the one this tool launched it in. A command that `chdir`s and then opens a relative
    path resolves somewhere this tool does not know, and such a path lands in
    `relative_reads_that_are_not_there` instead of being dropped, so a wrong resolution is a named
    figure in the report rather than silence.

    `output` is the result path this run WRITES. A read of it is not a read of unpinned bytes when the
    path held nothing before the command ran, which is why the exemption is conditional on
    `output_held_bytes_before_the_run` being measured False by the caller and is never taken on None.
    The exemption is ONE EXACT PATH: not a prefix, not `data/results/`, not the name with another
    suffix. `link_machine_local_inputs` already refuses to link that same path over
    (`not_linked_because_it_is_the_output`); this is its counterpart on the read side.
    """
    declared_files, declared_dirs = mf.declared_inputs(manifest)
    child_cwd = str(worktree.resolve())
    data = str((worktree / "data").resolve())
    under: list[str] = []
    relative: list[str] = []
    relative_under: list[str] = []
    relative_missing: list[str] = []
    for raw in recorded:
        is_relative = not os.path.isabs(raw)
        if is_relative:
            relative.append(raw)
            real = os.path.realpath(os.path.join(child_cwd, raw))
        else:
            real = os.path.realpath(raw)
        rel = None
        if real == data or real.startswith(data + os.sep):
            rel = os.path.relpath(real, child_cwd)
        else:
            for store in stores:
                if real == store or real.startswith(store + os.sep):
                    rel = "data/" + os.path.relpath(real, os.path.dirname(store))
                    break
        if rel is not None:
            rel = rel.replace(os.sep, "/")
            under.append(rel)
            if is_relative:
                relative_under.append(rel)
        elif is_relative and not os.path.exists(real):
            relative_missing.append(raw)
    under = sorted(set(under))
    framework = [p for p in under if p in FRAMEWORK_READS]
    rest = [p for p in under if p not in FRAMEWORK_READS]
    undeclared_including_the_output = sorted(
        p
        for p in rest
        if p not in declared_files
        and not any(p == d or p.startswith(d.rstrip("/") + "/") for d in declared_dirs)
    )
    own_output = [
        p
        for p in undeclared_including_the_output
        if output is not None and p == output and output_held_bytes_before_the_run is False
    ]
    undeclared = [p for p in undeclared_including_the_output if p not in own_output]
    return {
        "reads_under_data": under,
        "reads_under_data_count": len(under),
        "declared": sorted(declared_files),
        "declared_directories": sorted(declared_dirs),
        "relative_reads_recorded": sorted(set(relative)),
        "relative_reads_recorded_count": len(set(relative)),
        "relative_reads_resolved_under_data": sorted(set(relative_under)),
        "relative_reads_resolved_under_data_count": len(set(relative_under)),
        "relative_reads_resolved_against": child_cwd,
        "relative_reads_that_are_not_there": sorted(set(relative_missing)),
        "relative_reads_because": (
            "a relative read is resolved against the directory the child was launched in. Before "
            "2026-10-03 it was resolved against THIS TOOL'S working directory and then dropped as not "
            "under data/, which lost 26 of silenceragree's 3,340 recorded reads without saying so. A "
            "path in relative_reads_that_are_not_there resolved to nothing, which is what a command "
            "that chdir'd would leave behind: the loss is named, not silent"
        ),
        "framework_reads_ignored": framework,
        "framework_reads_ignored_because": (
            "read by save_result itself on every registry write, exempt by exact path. The loss, "
            "measured and not assumed: a change to the legacy allowlist cannot alter a value, only "
            "admission (tests/test_legacy_allowlist_exemption.py)"
        )
        if framework
        else "nothing was exempted",
        "read_but_not_declared": undeclared,
        "read_but_not_declared_because_it_is_this_run_s_own_output": own_output,
        "own_output_exempt_because": (
            f"{own_output[0]} is the file this run writes and the worktree held nothing at that path "
            "before the command ran, so those bytes are this run's product and not an input it was "
            "computed from. One exact path, and only on a measured absence: if the path had existed "
            "beforehand the read would be of bytes no sha256 pins and the verdict would be withheld"
        )
        if own_output
        else "nothing was exempted as this run's own output",
        "every_read_is_declared": not undeclared_including_the_output,
        "every_read_is_declared_or_is_this_run_s_own_output": not undeclared,
    }


def guard_window(root: Path) -> dict[str, Any]:
    """From which commit a rebuild verdict is covered by the write guard, and what is not covered.

    The commit is read from git rather than written here, so the window cannot drift from the file.

    The caveat is on every report on purpose. A verdict taken before the guard existed rests on the
    stores not having been written by a rebuild, which the audit in `scripts/rebuild_write_guard.py`
    could not establish and no audit of mtimes can: several of those directories are caches the results'
    OWN writers populate, and a rebuild runs the result's own command, so a cache write by the original
    run and one by a rebuild of that run are the same writer at the same path. The guard makes the
    assumption evidenced from now on, not retrospectively.
    """
    sha = _git_output(root, "log", "--diff-filter=A", "--format=%H", "-1", "--", str(_GUARD_FILE))
    return {
        "covers_verdicts_from": sha or "not committed yet (this guard is uncommitted in this checkout)",
        "earlier_verdicts": (
            "a rebuild verdict taken before that commit is NOT covered: until then the three data stores "
            "were symlinked into the worktree and writable through the link, and whether any rebuild "
            "wrote to one cannot be established from mtimes -- a cache write by the original run and one "
            "by a rebuild of that run are the same writer at the same path. 'Every modification is "
            "explained by a known writer' was established by the audit; 'no rebuild ever wrote to an "
            "input store' was not"
        ),
        "what_the_audit_could_not_see": [
            "an mtime can be rewritten",
            "a write that restored identical bytes leaves no trace in mtime or in content",
            "attribution by path shows that a writer exists, not that it made that particular file",
            "everything outside the three stores, data/results and data/organisms included",
        ],
    }


def _git_output(root: Path, *args: str) -> str:
    r = subprocess.run(["git", "-C", str(root), *args], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else ""


def race_name(before: str, after: str) -> str:
    """The name a race is reported under: `revision race (A -> B)`, oldest stamp first.

    Before is the cleanliness sha (sampled first, while the writer computed) and after is the write
    stamp, so the arrow runs in the direction HEAD moved.
    """
    return f"revision race ({before[:SHA_IN_PROSE]} → {after[:SHA_IN_PROSE]})"


def counting_path_diff(
    root: Path, before: str, after: str, paths: Sequence[str]
) -> tuple[list[str], list[str], str | None]:
    """(argv, files that changed, why it could not be asked). The decisive check of a revision race.

    `git diff --name-only A B -- <counting path>`: the files that differ between the two revisions, out of
    the ones the result's own counting path names. Nothing else is consulted, and no verdict is reached by
    reasoning about what a commit "probably" touched.

    A missing counting path is a refusal and not a pass: a race whose counting path is unknown cannot be
    called benign, so the reason comes back and the caller reports it instead of a verdict.
    """
    if not paths:
        return (
            [],
            [],
            (
                "the result records no code_cleanliness.counting_path, so there is nothing to diff: a race "
                "over an unknown counting path is not called benign"
            ),
        )
    argv = ["git", "-C", str(root), "diff", "--name-only", before, after, "--", *paths]
    r = subprocess.run(argv, capture_output=True, text=True)
    if r.returncode:
        return (
            argv,
            [],
            (
                f"git diff {before[:SHA_IN_PROSE]} {after[:SHA_IN_PROSE]} could not be run, so the race is "
                f"undecided: {r.stderr.strip()[-300:]}"
            ),
        )
    return argv, sorted(line for line in r.stdout.splitlines() if line.strip()), None


def revision_race(original: Any, root: Path) -> dict[str, Any]:
    """Whether this result's two revision stamps disagree, and if they do, whether it mattered.

    Always reported, in all three states, because each is a different claim:

    * no race -- both stamps record one revision;
    * a race whose `git diff` over the counting path is EMPTY -- "benign race: the counting path is
      identical at both revisions."; still named, still printed;
    * a race whose diff is NON-EMPTY -- a real difference, with the changed files named. The code that
      produced this result is not the code at the revision the rebuild runs.

    A state that cannot be established says so and is never folded into "no race".
    """
    m = original.get(mf.KEY) if isinstance(original, dict) else None
    stamps = mf.revision_stamps(m if isinstance(m, dict) else {})
    after, before = stamps["code_git_sha"], stamps["code_cleanliness_git_sha"]
    clean = (m or {}).get("code_cleanliness")
    paths = list(clean.get("counting_path") or []) if isinstance(clean, dict) else []
    out: dict[str, Any] = {
        "stamps": stamps,
        "counting_path_count": len(paths),
        "recorded_revision_stamps_block": (m or {}).get(mf.REVISION_STAMPS) is not None,
    }
    if stamps["revision_race"] is None:
        return {**out, "found": None, "verdict": stamps["reading"]}
    if stamps["revision_race"] is False:
        return {
            **out,
            "found": False,
            "verdict": (
                f"no revision race: both revision stamps record {after[:SHA_IN_PROSE]}, so the worktree "
                f"this rebuild runs in is the revision the cleanliness block was sampled at"
            ),
        }
    out["found"] = True
    out["race"] = race_name(before, after)
    argv, changed, why = counting_path_diff(root, before, after, paths)
    out["git_diff_argv"] = argv
    out["counting_path_files_changed"] = changed
    if why:
        out["verdict"] = f"{out['race']}: undecided. {why}"
        return out
    if not changed:
        out["verdict"] = (
            f"{out['race']}: benign race: the counting path is identical at both revisions. "
            f"All {len(paths)} counting-path files are byte-identical between "
            f"{before[:SHA_IN_PROSE]} and {after[:SHA_IN_PROSE]}, so HEAD moved on code that did not "
            f"produce these numbers"
        )
        return out
    out["verdict"] = (
        f"{out['race']}: a REAL difference: {len(changed)} of the {len(paths)} counting-path files "
        f"changed between the revisions: {', '.join(changed)}. The code that wrote this result is not "
        f"the code at the revision this rebuild runs"
    )
    return out


def must_hold_failures(original: dict[str, Any], rebuilt: dict[str, Any]) -> list[str]:
    """Every place either side breaks MUST_HOLD. Checked at any depth, on both sides."""
    out = []

    def walk(x: Any, side: str, where: str = "") -> None:
        if isinstance(x, dict):
            for k, v in x.items():
                for field, wanted in MUST_HOLD:
                    if k == field and v != wanted:
                        out.append(f"{side}{where}/{k}: {v!r}, must be {wanted!r}")
                walk(v, side, f"{where}/{k}")
        elif isinstance(x, list):
            for i, v in enumerate(x):
                walk(v, side, f"{where}[{i}]")

    walk(original, "original")
    walk(rebuilt, "rebuilt")
    return out


def is_timing(key: Any) -> bool:
    """A wall-clock key: `seconds`, or a snake_case name with a `seconds` token or ending `per_second`.
    `second` alone (an ordinal, as in second_endpoint) and `duration` (often biological) are compared."""
    if not isinstance(key, str):
        return False
    tokens = key.lower().split("_")
    return "seconds" in tokens or tokens[-2:] == ["per", "second"]


#: Run-resource readings, exempt from comparison by **exact key**. Added 2026-10-02 on the supervisor's
#: ruling, and deliberately not through ENVIRONMENT_FIELDS, which is the code-cleanliness list and stays
#: that way. A peak memory figure is of the same class as a wall-clock second: what the machine spent to
#: produce the result, not anything the result asserts. Two otherwise identical rebuilds differed in this
#: one leaf and nothing else -- `/cost/peak_memory_mb` 1142.5 against 1148.8 in placement_audit, and
#: `/compute/peak_rss_mb` 1300.5 against 1298.2 in context_contrast_feasibility. Exact keys only: a pattern
#: such as `*_mb` would one day swallow a leaf a run did compute, and the point of a short list is that
#: adding to it is a decision somebody makes on the record.
RESOURCE_KEYS = ("peak_rss_mb", "peak_memory_mb")


def is_resource(key: Any) -> bool:
    """A run-resource reading, matched by exact key and never by pattern, so a name that merely looks like
    one -- `peak_signal_mb`, `peak_rss_mb_per_cell`, `memory_mb_budget` -- is compared like any other leaf."""
    return isinstance(key, str) and key in RESOURCE_KEYS


# --- a tallied key name is not a field name (2026-10-02, lane-tracercwd) ---------------------------
#
# `is_timing` and `is_resource` match a key NAME at any depth, which is the rule this file's own
# RUN_FIELDS comment rejects for `opens` -- and it was already firing. In lane-generow's result the
# timing matcher set aside FOUR leaves and only one was a timing:
#
#   /key_vocabulary/per_element_archives/schema_keys/seconds   963406
#   /key_vocabulary/loose_per_element_files/schema_keys/seconds   3209
#   /key_vocabulary/the_newer_hct116_answers/schema_keys/seconds    705
#   /seconds                                                       336.6   <- the only timing
#
# The first three are KEY-FREQUENCY COUNTS: how many cached records carry a key of that name. They sit
# beside `id`, `chrom`, `start`, `end` and `genes` at the same 963,406, they are quantities the run
# computed, and the result's credibility rests on `differences: []` meaning that they were compared.
#
# Why the fix is not an exact-path list, as it was for `opens`: across the 1,100 results in this
# checkout a timing name appears at 536 DISTINCT paths, and every new result invents more. An explicit
# list is the right shape for one leaf and the wrong shape for 536. (Resource readings appear at only 3
# distinct paths, so those COULD be listed; they are narrowed the same way here instead, because the
# defect is in the matching rule and not in either list, and one rule is easier to reason about.)
#
# So the name test stays and is narrowed STRUCTURALLY: it does not apply to a key reached through a
# container whose own name says its keys came from the data. The direction of error is the point. A
# container wrongly treated as a tally means a genuine timing leaf gets COMPARED, which produces a
# loud difference somebody reads; a container wrongly treated as ordinary means a computed count is
# SILENTLY set aside. Only the second is the failure this tool must not have, so the container test is
# deliberately generous -- any `keys` or `vocabulary` token, anywhere in the ancestry -- and the report
# prints both what was set aside WITH ITS VALUE and what the narrowing kept, so neither is invisible.


def is_key_tally(container_key: Any) -> bool:
    """Whether a dict's keys are names found in the DATA rather than field names the writer chose.

    Nothing in a JSON document says which it is, so the containers whose keys are data are recognised by
    their own name: a `keys` or `vocabulary` token. That catches `schema_keys`, `all_keys`,
    `key_vocabulary` and `low_frequency_keys_that_are_not_element_ids`. Matched generously and on
    purpose: see the note above on which direction the error has to run.
    """
    if not isinstance(container_key, str):
        return False
    tokens = container_key.lower().split("_")
    return "keys" in tokens or "vocabulary" in tokens


def _ignored_key(key: Any, inside_key_tally: bool = False) -> bool:
    """A key comparable() drops: a wall-clock reading or a run-resource one. Both are reported by path
    and by value, under timing_fields_ignored and resource_fields_ignored, so a reader sees what was set
    aside and can see whether it looks like a reading at all.

    `inside_key_tally` says the key is a name the DATA supplied, not a field name the writer chose, in
    which case no name test applies to it and the leaf is compared like any other.
    """
    if inside_key_tally:
        return False
    return is_timing(key) or is_resource(key)


#: Printed beside each list of set-aside paths, so an ignored field is visible and auditable with its
#: reason attached, and never merely absent from the differences.
IGNORED_BECAUSE = {
    "environment_fields_ignored": (
        "the tree the run happened in, not anything the run computed: the code-cleanliness lists in "
        f"ENVIRONMENT_FIELDS, and a working directory recorded as an absolute path "
        f"({', '.join(CWD_FIELDS)}), "
        "which is either a cwd outside any repository or a result written before cwd_as_recorded; a "
        "difference between two REPOSITORY-RELATIVE working directories is a real difference and is not "
        "set aside"
    ),
    "run_fields_ignored": (
        f"KNOWN-LOSSY EXEMPTION, not noise: {', '.join(RUN_FIELDS)} is the audit hook's count of every "
        "open the process made, one integer naming no path. Setting it aside WILL HIDE A GENUINE CHANGE "
        "IN THE SET OF FILES THIS WRITER OPENS, and in the one case measured it moved for a real reason "
        "(1,684 to 1,855, a carry-forward that had begun reading a git blob in place of a file on disk). "
        "It is set aside because the count is too coarse to diagnose anything and a rebuild cannot act on "
        "it. It is the only leaf sensitive to reads OUTSIDE data/, so that is what the comparison is now "
        "blind to; every leaf about reads under data/ is still compared -- opens_under_data, files_read, "
        "files_written, active, and the undeclared and declared_not_read path lists. Exempt by exact path "
        "only, so an `opens` anywhere else -- /per_element_response_cache/opens is one, and the run "
        "computed it -- is compared like any other leaf"
    ),
    "timing_fields_ignored": (
        "a wall-clock reading, matched by key (is_timing) at any depth: `seconds`, a name with a "
        "`seconds` token, or one ending `per_second`. How long a machine took is not a quantity the "
        "result asserts; `second` alone (an ordinal, as in second_endpoint) and `duration` (often "
        "biological) are compared like any other leaf. A NAME TEST, not an exact path, because a timing "
        "name appears at 536 distinct paths across this checkout's results -- so it is narrowed by "
        "is_key_tally: no name test applies to a key the DATA supplied, and the leaves that rule kept "
        "comparable are printed under names_not_exempted_because_the_data_supplied_them. Every leaf set "
        "aside is printed WITH ITS VALUE under set_aside_by_a_name_test, so a count that is plainly not "
        "a reading is visible rather than merely absent"
    ),
    "resource_fields_ignored": (
        f"a run-resource reading, by exact key at any depth: {', '.join(RESOURCE_KEYS)} -- what the "
        "machine spent, not what the result asserts. Narrowed by is_key_tally on the same rule as the "
        "timing test, and printed with its value under set_aside_by_a_name_test"
    ),
}


def _strip_timing(x: Any, inside_key_tally: bool = False) -> Any:
    """`inside_key_tally` is sticky: once the ancestry has passed through a key tally every key below it
    came from the data, so no name test applies anywhere beneath."""
    if isinstance(x, dict):
        return {
            k: _strip_timing(v, inside_key_tally or is_key_tally(k))
            for k, v in x.items()
            if not _ignored_key(k, inside_key_tally)
        }
    if isinstance(x, list):
        return [_strip_timing(v, inside_key_tally) for v in x]
    return x


def _paths_where(x: Any, pred: Any, where: str = "", inside_key_tally: bool = False) -> list[str]:
    """Every path whose key satisfies `pred(key, inside_key_tally)`, with the tally ancestry carried down
    so a predicate can decide differently for a key the data supplied."""
    if isinstance(x, dict):
        out = []
        for k, v in x.items():
            below = inside_key_tally or is_key_tally(k)
            if pred(k, inside_key_tally):
                out.append(f"{where}/{k}")
            else:
                out.extend(_paths_where(v, pred, f"{where}/{k}", below))
        return out
    if isinstance(x, list):
        return [p for i, v in enumerate(x) for p in _paths_where(v, pred, f"{where}[{i}]", inside_key_tally)]
    return []


def timing_paths(x: Any, where: str = "") -> list[str]:
    """Every path at which comparable() drops a timing key. A timing NAME inside a key tally is not one
    of these, because it is not dropped: see `names_rescued_from_a_key_tally`."""
    return _paths_where(x, lambda k, tally: not tally and is_timing(k), where)


def resource_paths(x: Any, where: str = "") -> list[str]:
    """Every path at which comparable() drops a run-resource key."""
    return _paths_where(x, lambda k, tally: not tally and is_resource(k), where)


def names_rescued_from_a_key_tally(x: Any, where: str = "") -> list[str]:
    """Every path where a key WOULD have been set aside by a name test but is compared, because the name
    came from the data. The narrowing's own record: without it this rule would be invisible, and a
    reader could not tell a tool that compared these leaves from one that quietly dropped them."""
    return _paths_where(x, lambda k, tally: tally and (is_timing(k) or is_resource(k)), where)


def value_at(payload: Any, path: str) -> Any:
    """The value a `diff()`-style path names, list indices included, or None when nothing is there. Used
    to print a set-aside leaf's VALUE beside its path, so a count that is plainly not a reading -- a
    key-frequency 963,406 under `seconds` -- is visible to a reader rather than only to a test."""
    cur = payload
    for step in path.strip("/").split("/"):
        name, _, rest = step.partition("[")
        if name:
            if not isinstance(cur, dict) or name not in cur:
                return None
            cur = cur[name]
        for index in (i for i in rest.rstrip("]").split("][") if i):
            if not isinstance(cur, list) or not index.isdigit() or int(index) >= len(cur):
                return None
            cur = cur[int(index)]
    return cur


def diff(a: Any, b: Any, where: str = "") -> list[str]:
    """Every path at which two JSON values differ."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for k in sorted(set(a) | set(b), key=str):
            if k not in a or k not in b:
                out.append(f"{where}/{k}: only in {'rebuilt' if k in b else 'original'}")
            else:
                out.extend(diff(a[k], b[k], f"{where}/{k}"))
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return [f"{where}: {len(a)} items vs {len(b)}"]
        return [d for i, (x, y) in enumerate(zip(a, b, strict=True)) for d in diff(x, y, f"{where}[{i}]")]
    return [] if a == b else [f"{where}: {a!r} vs {b!r}"[:300]]


def comparable(payload: dict[str, Any]) -> dict[str, Any]:
    out = {k: v for k, v in payload.items() if k not in IGNORED}
    m = dict(out.get(mf.KEY) or {})
    m.pop("code", None)
    # recorded after the fact for results made before 2026-09-28 (R9 follow-up), so a rerun at their
    # revision cannot write it: it states the inputs' provenance, not a value the run computes
    m.pop("model_dependencies", None)
    if m:
        out[mf.KEY] = m
    out = _strip_timing(out)
    return out


def is_group(entry: dict[str, Any]) -> bool:
    """Whether a declared input stands for several files read together rather than for one path.

    `files_entry` marks a group since 2026-10-02. Before that a group carried `files` and so was
    indistinguishable from `input_entry` on a directory, which carries `files` too; such an entry is
    not treated as a group here, and falls to the unresolvable branch, which says so by name.
    """
    return entry.get("group") is True or isinstance(entry.get("members"), list)


def declared_paths(inputs: Any) -> list[str]:
    """Every repository-relative path a manifest's inputs name, a group's members included.

    A label is not a path, so a group contributes its members and not its own `path`; a group recorded
    before its members were named contributes nothing, which is why such a group remains a named failure
    in `_check_group` rather than something this could quietly satisfy.
    """
    out: list[str] = []
    for i in inputs if isinstance(inputs, list) else []:
        if not isinstance(i, dict):
            continue
        if is_group(i):
            out += [
                m["path"]
                for m in i.get("members") or []
                if isinstance(m, dict) and isinstance(m.get("path"), str)
            ]
        elif isinstance(i.get("path"), str):
            out.append(i["path"])
    return out


def repo_relative_declaration(path: str, root: Path) -> str | None:
    """The repository-relative form of a declared input path, or None when it is outside `root`.

    `genomeos.manifest.input_entry` records `str(path)` verbatim where its sibling `files_entry` puts
    every member through `_repo_relative`, so 9 of the 212 committed results carrying
    `result_manifest.inputs` declare an input as an ABSOLUTE path -- `scripts/silenceragree_run.py`
    passes `ROOT / REGISTRATION`. That file is inside the paid study's frozen import closure and cannot
    be edited, so the declaration is absorbed HERE instead: an absolute path that resolves inside the
    repository root names a file this checkout holds, and naming it relatively is all the writer failed
    to do.

    The test is RESOLUTION, not a string prefix, and the difference is not cosmetic: `root` itself is a
    symlinked path on this machine (/tmp, /var and the per-session scratch directories all resolve
    elsewhere), a declaration can carry `..`, and a worktree under /private/var is not the prefix its
    unresolved spelling suggests. Both sides go through `os.path.realpath` and the boundary is a path
    SEPARATOR, so a sibling directory whose name merely starts with the root's (a second checkout at
    `GenomeOS-old`) is outside and stays outside.

    A relative declaration is returned unchanged, except one that climbs out of the repository with
    `..`, which is outside by the same rule.
    """
    if not os.path.isabs(path):
        return None if ".." in Path(path).parts else path
    real = os.path.realpath(path)
    real_root = os.path.realpath(str(root))
    if real == real_root:
        return "."
    if real.startswith(real_root + os.sep):
        return os.path.relpath(real, real_root).replace(os.sep, "/")
    return None


def git_ignored(root: Path, paths: list[str]) -> set[str]:
    """Which of `paths` git ignores in `root`: the machine-local ones, answered for all of them at once.

    A path git does not ignore is never linked from this machine, whatever is sitting at it: for a tracked
    path the committed tree at the rebuild's sha is what the rebuild must read, and a path that is neither
    tracked nor ignored is not a machine-local store but a stray file. Both are left to be named absent.
    """
    if not paths:
        return set()
    r = subprocess.run(
        ["git", "-C", str(root), "check-ignore", "--stdin"],
        input="\n".join(paths),
        capture_output=True,
        text=True,
    )
    return {line for line in r.stdout.splitlines() if line}


def link_read_only(src: Path, dst: Path) -> bool:
    """Link one machine-local file into the worktree so the rebuild can read it and cannot write it.

    `cp -c` asks APFS for a clone: no bytes are copied, both sides share them until one is written, and the
    clone is a separate inode -- so the 0o444 that follows takes the write bits off the worktree's side only
    and never off the machine's one and only copy. That is the reason for a clone rather than a symlink or a
    hard link: both of those carry the machine's own mode, so `open(path, "w")` through either would truncate
    the machine's file, and the only way to refuse it would be to chmod the original, which other sessions
    share. Through this link such an open fails with PermissionError before a byte is written, and the
    rebuild cannot write a result into the directory it is reading its inputs from. On a filesystem with no
    clone support the fallback is a plain copy, read-only in the same way; a worktree copy that is still
    writable after the chmod is an error, not a warning.

    Returns whether the clone succeeded. False means the fallback ran and the bytes really were copied,
    which is what happens across a volume boundary; it is reported (`copied_not_cloned`) rather than left
    silent, because a clone that costs nothing and a copy that costs the file's size read the same in a
    report otherwise, and a rebuild that filled a disk would have had no way to say why.
    """
    dst.parent.mkdir(parents=True, exist_ok=True)
    cloned = subprocess.run(["cp", "-c", str(src), str(dst)], capture_output=True).returncode == 0
    if not cloned:
        shutil.copyfile(src, dst)
    dst.chmod(0o444)
    if dst.stat().st_mode & 0o222:
        raise OSError(f"{dst} is writable after chmod 0444: the rebuild could overwrite the input it reads")
    return cloned


def link_machine_local_inputs(
    worktree: Path, root: Path, inputs: Any, output: str | None = None
) -> dict[str, Any]:
    """Link the git-ignored data/results inputs this manifest declares into the worktree, read-only.

    Only a path the manifest declares is linked, only under MACHINE_LOCAL_DIR, only when git ignores it in
    this checkout, and only when the worktree does not already hold that path -- a file the commit carries
    is read from the commit. Nothing here decides whether an input matches: the bytes linked are hashed
    afterwards against the sha256 the manifest declares, so a machine-local input whose bytes have since
    changed is a reported difference and not a silent pass. `output` is the result the rebuild will write,
    which is this run's product rather than its input and is never linked over.

    An input declared as an ABSOLUTE path that resolves inside this repository is linked by its
    repository-relative form (`repo_relative_declaration`), because that is the file it names; one that
    resolves outside is not a path this checkout holds and is left for `_check_one_path` to refuse.
    """
    out: dict[str, Any] = {
        "linked": [],
        "absent_on_this_machine": [],
        "not_linked_because_git_does_not_ignore_it": [],
        "not_linked_because_it_is_the_output": [],
        # a clone shares its bytes and costs nothing; the fallback copies them. Named so a rebuild that
        # moved real gigabytes across a volume boundary can say so instead of looking like a clone.
        "copied_not_cloned": [],
        "failures": [],
    }
    want = sorted(
        {
            rel
            for rel in (repo_relative_declaration(p, root) for p in declared_paths(inputs))
            if rel is not None
            and (rel == MACHINE_LOCAL_DIR or rel.startswith(MACHINE_LOCAL_DIR + "/"))
            and not (worktree / rel).exists()
        }
    )
    ignored = git_ignored(root, want)
    for p in want:
        if output is not None and p == output:
            out["not_linked_because_it_is_the_output"].append(p)
            continue
        if p not in ignored:
            out["not_linked_because_git_does_not_ignore_it"].append(p)
            continue
        src = root / p
        if not src.exists():
            out["absent_on_this_machine"].append(p)
            continue
        try:
            if src.is_dir():
                for f in sorted(x for x in src.rglob("*") if x.is_file()):
                    if not link_read_only(f, worktree / p / f.relative_to(src)):
                        out["copied_not_cloned"].append(
                            {"path": str(f.relative_to(src)), "under": p, "bytes": f.stat().st_size}
                        )
            else:
                if not link_read_only(src, worktree / p):
                    out["copied_not_cloned"].append({"path": p, "bytes": src.stat().st_size})
        except OSError as e:
            out["failures"].append(f"{p}: {e}")
            continue
        out["linked"].append(p)
    return out


def unlock_links(worktree: Path, paths: list[str]) -> None:
    """Put the write bits back on the linked copies, so removing the worktree cannot be refused by them.
    Each is a clone of its own, so this touches nothing the machine keeps."""
    for p in paths:
        target = worktree / p
        for f in sorted(x for x in target.rglob("*") if x.is_file()) if target.is_dir() else [target]:
            if f.is_file():
                os.chmod(f, 0o644)


def check_inputs(
    worktree: Path, inputs: Any, root: Path, machine_local: set[str] | None = None
) -> dict[str, Any]:
    """Open and hash every declared input, a group member by member, and account for all of them.

    The contract the caller rests on: `unavailable` is empty only when at least one input was declared
    and every declared input was opened and hashed, and `entries` holds one record per declared input,
    so nothing can leave the list silently. An input this cannot resolve is a named failure.

    Written 2026-10-02 to replace a loop that resolved every input by `entry["path"]`. Because
    `files_entry` records a label there and not a path, each group was reported "absent" with no file
    opened, and its entry never reached `inputs`: a list with no denominator, which read as though
    every input had matched. The direction of that error flattered the result -- a result whose inputs
    were grouped was excused as resting on an unavailable dependency, when its files were on disk and
    could have been compared.
    """
    out: dict[str, Any] = {
        "declared": 0,
        "checked": 0,
        "files_opened": 0,
        # how many declared inputs were opened and hashed from a git-ignored, machine-local path rather
        # than from the committed tree: the share of the rebuild that rests on this machine
        "machine_local_opened": 0,
        "entries": [],
        "unavailable": [],
        "checked_outside_the_worktree": [],
        # an absolute declaration read by its repository-relative form: named here, never implied, so a
        # reader can see which inputs were checked through a path the writer did not record
        "relativised_absolute_declarations": [],
    }
    local = machine_local or set()
    if not isinstance(inputs, list) or not inputs:
        out["unavailable"].append(
            "the manifest declares no inputs, so nothing was opened or hashed: a comparison of this "
            "result's fields would say nothing about the bytes the run read"
        )
        return out
    out["declared"] = len(inputs)
    linked = [str((root / "data" / d).resolve()) for d in STORES]
    for i in inputs:
        if not isinstance(i, dict) or not isinstance(i.get("path"), str) or not i["path"]:
            out["entries"].append({"path": None, "sha256_matches": False, "problem": "no path"})
            out["unavailable"].append(f"a declared input carries no path: {i!r}"[:200])
            continue
        path = i["path"]
        if is_group(i):
            _check_group(worktree, i, out, local)
        elif i.get("sha256") == NOT_HASHED:
            out["entries"].append(
                {"path": path, "sha256_matches": False, "problem": f"sha256 is {NOT_HASHED!r}"}
            )
            out["unavailable"].append(
                f"input {path} records sha256 {NOT_HASHED!r}: nothing was hashed when the result was "
                "written, so there is nothing for a rebuild to check these bytes against"
            )
        else:
            _check_one_path(worktree, i, linked, root, out, local)
    return out


def _check_group(
    worktree: Path, i: dict[str, Any], out: dict[str, Any], local: set[str] | None = None
) -> None:
    """One grouped input: every member opened and hashed on its own, then the group's own digest."""
    path, members = i["path"], i.get("members")
    named = [m for m in members or [] if isinstance(m, dict) and isinstance(m.get("path"), str)]
    if not members or len(named) != len(members):
        out["entries"].append(
            {
                "path": path,
                "group": True,
                "members_declared": i.get("files"),
                "members_opened_and_hashed": 0,
                "sha256_matches": False,
                "problem": "the group names no members",
            }
        )
        out["unavailable"].append(
            f"input {path} is a group of {i.get('files', 'an unstated number of')} files recorded under "
            "a label whose members are not named, so not one of them could be opened: it was written "
            "before files_entry named its members (2026-10-02) and cannot be checked from this manifest"
        )
        return
    names = sorted(m["path"] for m in named)
    absent = [n for n in names if not (worktree / n).exists()]
    if absent:
        out["entries"].append(
            {
                "path": path,
                "group": True,
                "members_declared": len(names),
                "members_opened_and_hashed": 0,
                "absent_members": absent[:20],
                "sha256_matches": False,
            }
        )
        out["unavailable"].append(
            f"input {path}: {len(absent)} of {len(names)} files in the group are absent, first {absent[:3]}"
        )
        return
    try:
        digest, size, seen = mf.group_digest(names, root=worktree)
    except OSError as e:
        out["entries"].append(
            {
                "path": path,
                "group": True,
                "members_declared": len(names),
                "members_opened_and_hashed": 0,
                "sha256_matches": False,
                "problem": f"a member cannot be read: {e}",
            }
        )
        out["unavailable"].append(f"input {path}: a file in the group cannot be read: {e}")
        return
    out["files_opened"] += len(seen)
    recorded = {m["path"]: m.get("sha256") for m in named}
    differing = [s["path"] for s in seen if recorded.get(s["path"]) not in (None, s["sha256"])]
    ok = digest == i.get("sha256")
    from_machine = sorted(set(names) & (local or set()))
    out["entries"].append(
        {
            "path": path,
            "group": True,
            "members_declared": len(names),
            "members_opened_and_hashed": len(seen),
            "members_with_different_bytes": differing,
            "members_from_machine_local_paths": len(from_machine),
            "sha256_matches": ok and not differing,
            "bytes": size,
        }
    )
    if from_machine:
        # a group counts once, as one declared input, whichever of its members came from this machine
        out["machine_local_opened"] += 1
    if differing:
        out["unavailable"].append(
            f"input {path}: {len(differing)} of {len(names)} files in the group have different bytes, "
            f"first {differing[:3]}"
        )
        return
    if not ok:
        out["unavailable"].append(
            f"input {path}: the group's own digest differs from the manifest's (sha256 {digest[:12]})"
        )
        return
    out["checked"] += 1


def _check_one_path(
    worktree: Path,
    i: dict[str, Any],
    linked: list[str],
    root: Path,
    out: dict[str, Any],
    local: set[str] | None = None,
) -> None:
    """One input recorded as a single path, which may be a file or a directory."""
    path = i["path"]
    outside: dict[str, Any] | None = None
    read_as = path
    if os.path.isabs(path):
        # `worktree / path` leaves the worktree when path is absolute, so the bytes hashed would not be
        # the bytes the rebuild reads -- except under data/reference, data/knowledge and data/cache,
        # which are linked into the worktree and so are the same file.
        real = str(Path(path).resolve())
        store = next((s for s in linked if real == s or real.startswith(s + os.sep)), None)
        if store is not None:
            outside = {"path": path, "linked_store": store}
        else:
            # Since 2026-10-03: an absolute declaration that RESOLVES INSIDE the repository root names a
            # file this checkout does hold, so it is read by its repository-relative form instead of
            # being refused. The refusal's own reason -- "the bytes there are not the bytes this worktree
            # reads" -- stopped applying to that case the moment the path was known to be in-repo. One
            # resolved in-repo check and nothing wider: a path outside the root is still refused below,
            # and the root tested against is named in the refusal, because the same declaration is
            # in-repo when the tool runs from the checkout and out-of-repo when it runs from a worktree.
            rel = repo_relative_declaration(path, root)
            if rel is None:
                out["entries"].append(
                    {"path": path, "sha256_matches": False, "problem": "absolute, outside the linked stores"}
                )
                out["unavailable"].append(
                    f"input {path} is an absolute path outside the linked data stores and outside the "
                    f"repository root {os.path.realpath(str(root))}: the bytes there are not the bytes "
                    "this worktree reads, so it cannot be checked against this manifest"
                )
                return
            read_as = rel
            out["relativised_absolute_declarations"].append(
                {
                    "declared": path,
                    "read_as": rel,
                    "because": (
                        "it resolves inside the repository root, so it names a file this checkout holds; "
                        "genomeos.manifest.input_entry records str(path) where files_entry relativises"
                    ),
                }
            )
    p = worktree / read_as
    if not p.exists():
        # `files` on an entry with no member list reads as a group recorded under a label before
        # files_entry named its members; saying so is the difference between an input that is missing
        # and an input this tool cannot resolve.
        reads_as_group = isinstance(i.get("files"), int) and i["files"] > 1
        why = (
            f"is absent, and reads as a group of {i['files']} files recorded under a label before "
            "files_entry named its members (2026-10-02), so no file could be opened"
            if reads_as_group
            else "is absent"
        )
        out["entries"].append(
            {
                "path": path,
                "sha256_matches": False,
                "problem": "unresolved group" if reads_as_group else "absent",
            }
        )
        out["unavailable"].append(f"input {path} {why}")
        return
    try:
        digest, size, count = mf.sha256_of(p)
    except OSError as e:
        out["entries"].append({"path": path, "sha256_matches": False, "problem": f"cannot be read: {e}"})
        out["unavailable"].append(f"input {path} cannot be read: {e}")
        return
    if outside is not None:
        out["checked_outside_the_worktree"].append(outside)
    ok = digest == i.get("sha256")
    out["files_opened"] += count
    entry = {"path": path, "sha256_matches": ok, "bytes": size, "files_opened": count}
    if read_as != path:
        entry["read_as"] = read_as
    if read_as in (local or set()):
        # opened and hashed, and counted here, whether or not it matched: a machine-local input whose
        # bytes differ is still an input this rebuild took from this machine, and the difference is
        # reported below. The count says what the rebuild rested on, not what passed.
        entry["machine_local"] = True
        out["machine_local_opened"] += 1
    out["entries"].append(entry)
    if not ok:
        out["unavailable"].append(f"input {path} has different bytes (sha256 {digest[:12]})")
        return
    out["checked"] += 1


def rebuild(
    result: Path, where: Path, venv: str = "fresh", keep: bool = False, inputs_only: bool = False
) -> dict[str, Any]:
    root = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], text=True).strip())
    original = json.loads(result.read_text())
    m = original.get(mf.KEY)
    if not isinstance(m, dict):
        return {"rebuilt": False, "unavailable": [f"{result} carries no {mf.KEY}: nothing to rebuild from"]}
    code = m.get("code") or {}
    sha, argv = code.get("git_sha"), code.get("argv")
    report: dict[str, Any] = {
        "result": str(result),
        "git_sha": sha,
        "argv": argv,
        "dirty_code_at_write": code.get("dirty_code_paths", []),
        # named and decided before anything is built: the check is `git diff` over the two recorded
        # revisions, so it costs nothing and is reported even when the rebuild stops below
        "revision_race": revision_race(original, root),
        # and the tool's own cleanliness, likewise before anything is built
        "verification_tool": tool_cleanliness(root),
        "unavailable": [],
    }
    if not report["verification_tool"]["tool_is_committed"]:
        report["unavailable"].append(
            "no verdict: the verification tool is uncommitted. "
            + ", ".join(report["verification_tool"]["uncommitted_code_on_the_tool_s_own_closure"])
            + " is on this tool's own import closure and is not in any commit, so a verdict here would "
            "rest on verification code nobody can reproduce. Run this from a worktree at the committed "
            "revision instead; nothing was built and no command was run"
        )
        return {**report, "rebuilt": False}
    if not sha or not argv:
        report["unavailable"].append("manifest has no git_sha or argv")
        return {**report, "rebuilt": False}
    wt = where / f"rebuild-{sha[:10]}"
    # Outside the worktree on purpose: a file inside it would show up as an untracked path in the rebuilt
    # result's own cleanliness block, which this tool would then have to explain away.
    guard_dir = where / f"guard-{sha[:10]}"
    # The stores this rebuild must not write to, by real path. data/results is protected by its clones at
    # mode 0444 (link_read_only); these three are symlinked, and a symlink carries the target's mode, so
    # until 2026-10-02 a write through one reached the only copy on this machine. Measured: an
    # open-for-update of data/reference/HG002_chr1.vcf.gz through the link succeeded.
    stores = [str((root / "data" / d).resolve()) for d in STORES if (root / "data" / d).is_dir()]
    guard.arm(stores)
    subprocess.run(["git", "-C", str(root), "worktree", "add", "-q", "--detach", str(wt), sha], check=True)
    links: dict[str, Any] = {"linked": []}
    try:
        for d in STORES:
            if (root / "data" / d).is_dir() and not (wt / "data" / d).exists():
                (wt / "data" / d).symlink_to(root / "data" / d)
        report["write_guard"] = {
            "protected_stores": stores,
            "armed_in_this_process": list(guard.armed()),
            "events_watched": ["open (writing mode or flag)", *sorted(guard.WATCHED)],
            "how": (
                "an audit hook raises PermissionError on a write whose real path resolves into one of "
                "the stores, in this process and in the rebuilt command's own process"
            ),
            "window": guard_window(root),
            "limit": (
                "a Python-level guard only: a C library that opens a file itself (htslib, so pysam) is "
                "invisible to an audit hook, as is any write by a subprocess the guard was not put into. "
                "The clones at mode 0444 under data/results are kernel-enforced; this is not"
            ),
        }
        links = link_machine_local_inputs(wt, root, m.get("inputs"), output=f"data/results/{result.name}")
        found = check_inputs(wt, m.get("inputs"), root, machine_local=set(links["linked"]))
        report["inputs"] = found["entries"]
        report["inputs_declared"] = found["declared"]
        report["inputs_checked"] = found["checked"]
        report["inputs_unchecked"] = found["declared"] - found["checked"]
        report["files_opened_and_hashed"] = found["files_opened"]
        report["inputs_satisfied_from_machine_local_paths"] = found["machine_local_opened"]
        report["machine_local_inputs"] = {
            "count": found["machine_local_opened"],
            "of_declared": found["declared"],
            "linked_read_only_from": str(root / MACHINE_LOCAL_DIR),
            "how": "a clone at mode 0o444: the rebuild reads these bytes and a writer gets PermissionError",
            "paths": links["linked"],
            **{k: v for k, v in links.items() if k != "linked" and v},
        }
        report["unavailable"] += [
            f"a machine-local input could not be linked read-only: {f}" for f in links["failures"]
        ]
        if found["checked_outside_the_worktree"]:
            report["checked_outside_the_worktree"] = found["checked_outside_the_worktree"]
        if found["relativised_absolute_declarations"]:
            report["relativised_absolute_declarations"] = found["relativised_absolute_declarations"]
        report["unavailable"] += found["unavailable"]
        if report["unavailable"]:
            return {**report, "rebuilt": False}
        if inputs_only:
            # Steps 1 to 3 only. Item 10's capacity rule allows two lanes at once to read the per-element
            # cache, whose reader peaks near 2.9 GB, so a rebuild whose command reads it cannot be run by a
            # third lane. The inputs can still be linked and hashed, which costs a megabyte at a time, and
            # what that buys is the input figure and nothing else: no command ran, so nothing was compared.
            report["unavailable"].append(
                "--inputs-only: the command was not run, so no field was compared and no verdict is "
                "reported; the input figures above are all this says"
            )
            return {**report, "rebuilt": False}
        env = {k: v for k, v in os.environ.items() if k not in ("VIRTUAL_ENV", "UV_PROJECT_ENVIRONMENT")}
        if venv == "shared":
            env.update(UV_NO_SYNC="1", UV_PROJECT_ENVIRONMENT=str(root / ".venv"))
        else:
            sync = subprocess.run(
                ["uv", "sync", "--frozen", "--offline", "--quiet"],
                cwd=wt,
                env=env,
                capture_output=True,
                text=True,
            )
            report["venv"] = "fresh, from the committed uv.lock, offline"
            if sync.returncode:
                report["unavailable"].append(
                    f"fresh environment: uv sync --offline failed: {sync.stderr[-400:]}"
                )
                return {**report, "rebuilt": False}
        # sitecustomize is written beside the worktree and put FIRST on the child's path, so CPython's
        # `site` arms the guard before the result's writer runs its first line. An audit hook in this
        # process does not reach a subprocess, and the subprocess is where the risk is: this parent only
        # reads and hashes, while the rebuilt command is arbitrary committed code.
        env.update(guard.write_sitecustomize(guard_dir, _GUARD_FILE))
        env[guard.ENV_ROOTS] = os.pathsep.join(stores)
        env["PYTHONPATH"] = os.pathsep.join([str(guard_dir), str(wt)])
        # The committed uv.lock is used as it is. Without UV_FROZEN `uv run` re-locks (the lock at
        # fe0880a still names genomeos 0.9.0 against pyproject's 1.0.0) and the rewritten lock makes
        # the second checkout dirty, which the rebuilt manifest then records as a one-byte difference.
        env.update(UV_FROZEN="1", UV_OFFLINE="1")
        # Measured HERE, before the command runs, because the command is about to write this path and
        # after that no observation can tell what was there first. A run that reads back the result it
        # writes -- `save_result` does, on the registration results -- was refused a verdict outright
        # until 2026-10-03, before a single leaf was compared. The read is exempt only on this measured
        # absence: if the worktree already held bytes at that path, the run read bytes no sha256 pins
        # and the refusal stands.
        name = original.get("result") or result.stem
        output_rel = f"data/results/{name}.json"
        output_held_bytes_before_the_run = (wt / output_rel).exists()
        report["output_of_this_run"] = {
            "path": output_rel,
            "held_bytes_in_the_worktree_before_the_command_ran": output_held_bytes_before_the_run,
            "measured": "by os.path.exists in the worktree immediately before the command was started",
            "what_it_decides": (
                "whether a recorded read of that one exact path is this run's own product (exempt) or a "
                "read of pre-existing bytes pinned by no sha256 (no verdict)"
            ),
        }
        run = subprocess.run(["uv", "run", "python", *argv], cwd=wt, env=env, capture_output=True, text=True)
        report["exit"] = run.returncode
        report["stdout_tail"] = run.stdout[-1500:]
        # evidence, not a hope: the child records its pid once the hook is armed. An unarmed child means
        # the command ran with the stores writable, which is a reason to report no verdict rather than a
        # detail -- the inputs the comparison rests on were not protected while the comparison was made.
        armed_pids = guard.marker_pids(guard_dir)
        report["write_guard"]["armed_in_the_run"] = bool(armed_pids)
        report["write_guard"]["processes_that_armed_it"] = armed_pids
        # what the run READ, from the same audit hook that refused its writes
        recorded = guard.reads_recorded(guard_dir)
        report["reads_recorded"] = recorded
        report["reads_reconciled"] = reads_against_declared(
            wt,
            m,
            recorded["paths"],
            stores,
            output=output_rel,
            output_held_bytes_before_the_run=output_held_bytes_before_the_run,
        )
        if not armed_pids:
            report["unavailable"].append(
                "the write guard did not arm in the rebuilt command's process, so the command ran with "
                f"the data stores writable through the worktree's symlinks ({guard.MARKER} is absent "
                f"from {guard_dir}): no verdict is reported"
            )
            return {**report, "rebuilt": False}
        # No verdict unless reads were SEEN and the guard was IN FORCE. An unwritten record is not
        # "no undeclared reads": that is exactly how a pool child that read nothing produced a result
        # with 0 paths opened, satisfying "every read is declared" vacuously.
        if not recorded["recorded"]:
            report["unavailable"].append(
                "no verdict: no read was recorded for this run, so nothing is known about what it "
                f"opened. An empty read record is not a record of no reads (looked in {guard_dir} for "
                f"{guard.READS_PREFIX}*.tsv)"
            )
            return {**report, "rebuilt": False}
        if report["reads_reconciled"]["read_but_not_declared"]:
            report["unavailable"].append(
                "no verdict: the rebuilt command read "
                + ", ".join(report["reads_reconciled"]["read_but_not_declared"][:20])
                + " under data/, which this manifest's inputs neither name nor contain, so those bytes "
                "are pinned by no sha256 and a comparison of fields says nothing about them. This is "
                "the blindness measured on constrained_unknown_targets: 193 of 193 declared inputs "
                "hashed and a 194th file read unpinned"
            )
            return {**report, "rebuilt": False}
        if run.returncode:
            report["unavailable"].append(f"the command failed: {run.stderr[-800:]}")
            return {**report, "rebuilt": False}
        name = original.get("result") or result.stem
        written = wt / "data" / "results" / f"{name}.json"
        if written.name != result.name:  # the manifest was read from a copy under another name
            shutil.copyfile(written, wt / "data" / "results" / result.name)
        rebuilt = json.loads((wt / "data" / "results" / result.name).read_text())
        report["rebuilt_sha256"] = mf.sha256_of(written)[0]
        found = diff(comparable(original), comparable(rebuilt))
        rest, run = run_differences(found)
        real, env = environment_differences(rest, cwd_fields_recorded_absolute(original, rebuilt))
        report["differences"] = real
        report["revision_race"]["differences_it_accounts_for"] = {
            "paths": [d for d in real if _base_path(d) in RACE_FIELDS],
            "set_aside": False,
            "why_not": (
                "the race explains this leaf, it does not excuse it: an exemption here would hide that "
                "HEAD moved during the write, so the leaf is reported as a difference and the verdict "
                "above says whether it mattered"
            ),
        }
        report["environment_fields_ignored"] = env
        report["run_fields_ignored"] = run
        report["must_hold_failures"] = must_hold_failures(original, rebuilt)
        report["differences"] += report["must_hold_failures"]
        report["timing_fields_ignored"] = sorted(set(timing_paths(original)) | set(timing_paths(rebuilt)))
        report["resource_fields_ignored"] = sorted(
            set(resource_paths(original)) | set(resource_paths(rebuilt))
        )
        # the VALUE beside the path, because a path alone cannot show that a leaf set aside as a reading
        # is not one: `/key_vocabulary/.../schema_keys/seconds` reads as a timing until its 963,406 is
        # printed next to it. An over-exemption has to be visible to a reader, not only to a test.
        report["set_aside_by_a_name_test"] = {
            path: value_at(original, path)
            for path in report["timing_fields_ignored"] + report["resource_fields_ignored"]
        }
        # and the narrowing's own record, so the rule is auditable rather than invisible
        report["names_not_exempted_because_the_data_supplied_them"] = {
            path: value_at(original, path)
            for path in sorted(
                set(names_rescued_from_a_key_tally(original)) | set(names_rescued_from_a_key_tally(rebuilt))
            )
        }
        report["ignored_because"] = {k: IGNORED_BECAUSE[k] for k in IGNORED_BECAUSE if report.get(k)} or {
            "nothing": "no leaf was set aside by any of the four rules"
        }
        report["fields_compared"] = len(comparable(original))
        report["leaves"] = leaf_reconciliation(original)
        report["leaves_compared"] = report["leaves"]["compared"]
        # despite its name this key compares fields (date and the code block ignored); the next is bytes
        # (timing keys are ignored as well, and listed in timing_fields_ignored)
        why = no_verdict_reason(report["inputs_declared"], report["inputs_checked"])
        if why:
            # unreachable through the branch above; asked again so that no later edit can write a verdict
            # for a result whose inputs were not all opened, the failure this tool had on 2026-10-02
            report["unavailable"].append(why)
            return {**report, "rebuilt": False}
        report["identical_bytes_except_date_and_run"] = not report["differences"]
        report["identical_bytes"] = written.read_bytes() == result.read_bytes()
        set_aside = {k: report[k] for k in IGNORED_BECAUSE if report.get(k)}
        report["leaves_set_aside"] = sum(len(v) for v in set_aside.values())
        report["no_differences_and_nothing_set_aside"] = (
            not report["differences"] and not report["leaves_set_aside"]
        )
        report["comparison_reading"] = comparison_reading(report["differences"], set_aside)
        return {**report, "rebuilt": True}
    finally:
        # the hook cannot be uninstalled, so its roots are emptied: it stays in place and inert, which is
        # what the rest of this process (a test session above all) must be left with
        guard.disarm()
        if not keep:
            unlock_links(wt, links["linked"])
            subprocess.run(["git", "-C", str(root), "worktree", "remove", "--force", str(wt)], check=False)
            shutil.rmtree(guard_dir, ignore_errors=True)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("result", type=Path)
    ap.add_argument("--where", type=Path, required=True, help="directory for the second checkout")
    ap.add_argument("--venv", choices=("fresh", "shared"), default="fresh")
    ap.add_argument("--keep", action="store_true", help="leave the worktree in place")
    ap.add_argument(
        "--inputs-only",
        action="store_true",
        help="link and hash every declared input, then stop: no command is run and no verdict is reported",
    )
    args = ap.parse_args(argv)
    r = rebuild(args.result.resolve(), args.where.resolve(), args.venv, args.keep, args.inputs_only)
    json.dump({k: v for k, v in r.items() if k != "stdout_tail"}, sys.stdout, indent=1)
    print()
    return 0 if r.get("rebuilt") and not r.get("differences") else 1


if __name__ == "__main__":
    raise SystemExit(main())
