# SPDX-License-Identifier: AGPL-3.0-or-later
# Part of the GenomeOS application; see LICENSING.md.
"""The tracer's own three leaves, which made every rebuild report three differences that carried nothing.

Two separate lanes finished a rebuild that reconciled at every leaf and then had to inspect and argue the
same three away by hand before they could report. A comparison that always reports three differences
teaches its readers to skim the difference list, which is the opposite of what it is for.

The two halves are NOT the same kind of fix and the tests keep them apart:

* `cwd_at_open` and `cwd_at_close` are made COMPARABLE, not exempt. An absolute path carries no
  information outside the checkout that wrote it -- a second worktree is at a different absolute path by
  construction -- so recording it relative to the repository root loses nothing and keeps the purpose
  manifest.py states for the field. Only a cwd outside any repository, which cannot be made relative, is
  set aside, and only on the value being absolute.
* `opens` is EXEMPT, and the exemption is a known loss rather than a free one. In the one case measured
  the count moved for a real reason: a writer's carry-forward had begun reading a git blob in place of a
  file on disk. The tests below require the report to say that out loud, and require that the exemption
  is by exact path, because a real result asserts a computed `opens` at another path.

Every exemption here is locked by a counterfactual: a real change is planted and the comparison is shown
to still fail, and for the exact-path rule the wrong implementation is written out and shown to swallow a
quantity a run computed. An exemption with no counterfactual is a hole.
"""

import importlib.util
import os
import subprocess
from pathlib import Path

import pytest

from genomeos import manifest as mf

spec = importlib.util.spec_from_file_location("mr", Path("scripts/manifest_rebuild.py"))
mr = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mr)

CWD_OPEN = "/result_manifest/traced_inputs/cwd_at_open"
CWD_CLOSE = "/result_manifest/traced_inputs/cwd_at_close"
OPENS = "/result_manifest/traced_inputs/opens"


def _a_repo(where: Path) -> Path:
    """A git repository at `where`, which is all `_repo_root` looks for."""
    where.mkdir(parents=True, exist_ok=True)
    subprocess.run(["git", "init", "-q", str(where)], check=True)
    return where


def _traced(**fields):
    """A result carrying a traced_inputs block, as a writer records it."""
    return {"result": "x", "result_manifest": {"traced_inputs": dict(fields)}}


# --------------------------------------------------------------------------------------------------
# 1. The recorder: a working directory as a second checkout can compare it.
# --------------------------------------------------------------------------------------------------


def test_the_same_directory_in_two_checkouts_records_one_value(tmp_path, monkeypatch):
    """The whole defect, and the whole fix. Two checkouts at different absolute paths; the cwd recorded
    for each must be the same string, or no rebuild can ever match the field."""
    one, two = _a_repo(tmp_path / "checkout-one"), _a_repo(tmp_path / "a" / "b" / "checkout-two")
    monkeypatch.chdir(one)
    first = mf.cwd_as_recorded()
    monkeypatch.chdir(two)
    second = mf.cwd_as_recorded()
    assert first == second == mf.CWD_AT_ROOT
    # the two checkouts must be at different absolute paths for this to mean anything, and the old
    # recording -- the absolute path -- is what every rebuild compared until now
    assert str(one.resolve()) != str(two.resolve())


def test_a_cwd_inside_the_repository_records_its_relative_path(tmp_path, monkeypatch):
    root = _a_repo(tmp_path / "checkout")
    (root / "scripts" / "deep").mkdir(parents=True)
    monkeypatch.chdir(root / "scripts" / "deep")
    assert mf.cwd_as_recorded() == "scripts/deep"


def test_a_cwd_outside_any_repository_is_recorded_absolute(tmp_path, monkeypatch):
    """The one case that cannot be made relative to a repository root, because there is none."""
    outside = tmp_path / "no-repo-here"
    outside.mkdir()
    monkeypatch.chdir(outside)
    recorded = mf.cwd_as_recorded()
    assert os.path.isabs(recorded), recorded
    assert mf._repo_root(outside) is None


def test_a_linked_worktree_is_its_own_root(tmp_path, monkeypatch):
    """A worktree's `.git` is a FILE, a clone's is a directory. A run inside `.claude/worktrees/<name>`
    is at the root of its own checkout, and `.` is the honest answer there: that is the directory its
    recorded relative paths were taken against. Existence is the test, never `is_dir`."""
    root = _a_repo(tmp_path / "checkout")
    worktree = root / ".claude" / "worktrees" / "agent-abc"
    worktree.mkdir(parents=True)
    (worktree / ".git").write_text("gitdir: ../../../.git/worktrees/agent-abc\n")
    assert (worktree / ".git").is_file()
    monkeypatch.chdir(worktree)
    assert mf.cwd_as_recorded() == mf.CWD_AT_ROOT


def test_the_window_records_the_repository_relative_cwd_at_both_ends(tmp_path, monkeypatch):
    root = _a_repo(tmp_path / "checkout")
    monkeypatch.chdir(root)
    mf.trace_begin()
    held = mf.trace_close()
    assert held["cwd_at_open"] == held["cwd_at_close"] == mf.CWD_AT_ROOT


def test_a_cwd_that_moved_inside_the_window_still_records_differently(tmp_path, monkeypatch):
    """The purpose the field exists for, preserved. manifest.py records these two "rather than leaving it
    to be assumed" that the working directory held still; a fix that made every cwd record the same value
    would have destroyed that, which is why `.` is the root and not every directory."""
    root = _a_repo(tmp_path / "checkout")
    (root / "elsewhere").mkdir()
    monkeypatch.chdir(root)
    mf.trace_begin()
    monkeypatch.chdir(root / "elsewhere")
    held = mf.trace_close()
    assert held["cwd_at_open"] == mf.CWD_AT_ROOT
    assert held["cwd_at_close"] == "elsewhere"
    assert held["cwd_at_open"] != held["cwd_at_close"], "a cwd that moved must still be visible"


# --------------------------------------------------------------------------------------------------
# 2. The comparison: a relative cwd difference is real; only an absolute value is environment.
# --------------------------------------------------------------------------------------------------


def test_an_absolute_working_directory_is_set_aside_as_environment():
    """Either a cwd outside any repository, or a result written before cwd_as_recorded. A committed
    result's bytes are a pin and are not rewritten, so the second case is served here or nowhere."""
    original = _traced(cwd_at_open="/Users/someone/GenomeOS", cwd_at_close="/Users/someone/GenomeOS")
    rebuilt = _traced(cwd_at_open=".", cwd_at_close=".")
    found = mr.diff(mr.comparable(original), mr.comparable(rebuilt))
    real, env = mr.environment_differences(found, mr.cwd_fields_recorded_absolute(original, rebuilt))
    assert real == [], real
    assert len(env) == 2


def test_two_relative_working_directories_that_differ_are_a_real_difference():
    """THE COUNTERFACTUAL for the cwd clause. The run was in another directory, so every relative path it
    recorded means something else than the paths it is being compared against. That must fail."""
    original = _traced(cwd_at_open=".", cwd_at_close=".")
    rebuilt = _traced(cwd_at_open="scripts", cwd_at_close="scripts")
    found = mr.diff(mr.comparable(original), mr.comparable(rebuilt))
    assert mr.cwd_fields_recorded_absolute(original, rebuilt) == [], "neither value is absolute"
    real, env = mr.environment_differences(found, mr.cwd_fields_recorded_absolute(original, rebuilt))
    assert len(real) == 2, real
    assert env == []


def test_the_absolute_clause_is_decided_per_field_and_not_for_the_block():
    """One leaf absolute must not carry the other. A run that began outside the repository and ended
    inside it is the case: the open is set aside, the close is compared."""
    original = _traced(cwd_at_open="/tmp/elsewhere", cwd_at_close=".")
    rebuilt = _traced(cwd_at_open=".", cwd_at_close="scripts")
    assert mr.cwd_fields_recorded_absolute(original, rebuilt) == [CWD_OPEN]
    found = mr.diff(mr.comparable(original), mr.comparable(rebuilt))
    real, env = mr.environment_differences(found, mr.cwd_fields_recorded_absolute(original, rebuilt))
    assert [d.split(":")[0] for d in env] == [CWD_OPEN]
    assert [d.split(":")[0] for d in real] == [CWD_CLOSE]


def test_a_cwd_field_at_any_other_path_is_never_set_aside():
    """Exact paths. A result that computed something it called cwd_at_open elsewhere is compared."""
    original = {"survey": {"cwd_at_open": "/a"}, "result_manifest": {}}
    rebuilt = {"survey": {"cwd_at_open": "/b"}, "result_manifest": {}}
    assert mr.cwd_fields_recorded_absolute(original, rebuilt) == []
    real, env = mr.environment_differences(
        mr.diff(mr.comparable(original), mr.comparable(rebuilt)),
        mr.cwd_fields_recorded_absolute(original, rebuilt),
    )
    assert real == ["/survey/cwd_at_open: '/a' vs '/b'"]
    assert env == []


def test_the_cwd_fields_are_not_in_the_unconditional_exempt_list():
    """They are compared leaves. Putting them in ENVIRONMENT_FIELDS would exempt them whatever they hold,
    which is the fix this one was chosen over."""
    assert CWD_OPEN not in mr.ENVIRONMENT_FIELDS
    assert CWD_CLOSE not in mr.ENVIRONMENT_FIELDS
    assert CWD_OPEN not in mr.RUN_FIELDS


# --------------------------------------------------------------------------------------------------
# 3. `opens`: exempt by exact path, and the exact path is load-bearing.
# --------------------------------------------------------------------------------------------------


def test_the_tracer_opens_count_is_exempt_by_its_exact_path():
    original = _traced(opens=1684, opens_under_data=3)
    rebuilt = _traced(opens=1855, opens_under_data=3)
    rest, run = mr.run_differences(mr.diff(mr.comparable(original), mr.comparable(rebuilt)))
    assert rest == []
    assert [d.split(":")[0] for d in run] == [OPENS]


def test_the_exemption_is_one_exact_path_and_not_a_key_a_suffix_or_a_substring():
    assert mr.RUN_FIELDS == (OPENS,)
    assert all(f.startswith("/") for f in mr.RUN_FIELDS), "a bare key would match at every depth"


def test_a_computed_opens_at_another_path_is_compared():
    """THE CASE THAT MAKES THE EXACT PATH LOAD-BEARING, and it is a real one, not a hypothetical:
    data/results/response_map_increment3.json records

        "per_element_response_cache": {"opens": 0, "files": [], "how_this_is_known": "recorded, not
         inferred: `builtins.open` was patched for the whole run ..."}

    That zero is an assertion the run COMPUTED -- the run opened no file under the per-element cache --
    recorded that way precisely because no reading of the source could establish it. A rebuild in which
    the run did open the cache must fail.
    """
    original = {"per_element_response_cache": {"opens": 0, "files": []}, "result_manifest": {}}
    rebuilt = {"per_element_response_cache": {"opens": 7, "files": ["chr1.json"]}, "result_manifest": {}}
    rest, run = mr.run_differences(mr.diff(mr.comparable(original), mr.comparable(rebuilt)))
    assert run == [], "nothing here is a run field"
    assert any(d.startswith("/per_element_response_cache/opens") for d in rest), rest
    real, _ = mr.environment_differences(rest, mr.cwd_fields_recorded_absolute(original, rebuilt))
    assert len(real) == 2, real


# -------------------------------------------------- the counterfactual, by writing the wrong rule out


def _run_differences_by_suffix(found):
    """`run_differences` as it would be with a suffix match in place of an exact path, and nothing else
    changed. This is the obvious implementation and the one the exact-path rule was chosen over."""
    rest, run = [], []
    for d in found:
        (run if d.split(":")[0].split("[")[0].endswith("opens") else rest).append(d)
    return rest, run


def test_a_suffix_match_would_swallow_the_computed_opens_and_the_exact_path_does_not():
    """The exemption, credited only with what it is shown not to do. Same two results, same difference:
    under a suffix match the computed assertion is set aside and the rebuild passes; under the exact path
    it is a difference and the rebuild fails. The harm happens when the rule is loosened, so this is a
    protection and not a diagnosis."""
    original = {"per_element_response_cache": {"opens": 0}, "result_manifest": {}}
    rebuilt = {"per_element_response_cache": {"opens": 7}, "result_manifest": {}}
    found = mr.diff(mr.comparable(original), mr.comparable(rebuilt))

    swallowed_rest, swallowed_run = _run_differences_by_suffix(found)
    assert swallowed_rest == [], "the wrong rule leaves nothing to fail on"
    assert len(swallowed_run) == 1

    kept_rest, kept_run = mr.run_differences(found)
    assert kept_run == []
    assert len(kept_rest) == 1, "the exact-path rule must leave the computed difference to fail on"


def test_the_run_exemption_cannot_hide_a_real_change_in_a_computed_quantity():
    """THE COUNTERFACTUAL for the run-field exemption: a real change planted beside the exempt one. The
    exempt leaf is set aside and the planted one still fails, so the exemption cannot turn a defect into
    a pass by being present."""
    original = {
        "verdict": {"pooled_independent_loci_genome_wide": 22},
        "result_manifest": {"traced_inputs": {"opens": 1684, "opens_under_data": 3, "files_read": 1}},
    }
    rebuilt = {
        "verdict": {"pooled_independent_loci_genome_wide": 21},
        "result_manifest": {"traced_inputs": {"opens": 1855, "opens_under_data": 3, "files_read": 1}},
    }
    rest, run = mr.run_differences(mr.diff(mr.comparable(original), mr.comparable(rebuilt)))
    real, _ = mr.environment_differences(rest, mr.cwd_fields_recorded_absolute(original, rebuilt))
    assert [d.split(":")[0] for d in run] == [OPENS]
    assert real == ["/verdict/pooled_independent_loci_genome_wide: 22 vs 21"]


@pytest.mark.parametrize(
    "field,before,after",
    [
        ("opens_under_data", 3, 4),
        ("files_read", 1, 2),
        ("files_written", 0, 1),
        ("active", True, False),
        ("undeclared", [], ["data/knowledge/alphagenome/all_elements/chr1.json"]),
        ("declared_not_read", [], ["data/results/carried_forward.json"]),
    ],
)
def test_every_other_tracer_leaf_is_still_compared(field, before, after):
    """What is NOT given up. The exemption covers the total count and nothing else, so everything the
    tracer records about reads UNDER data/ still has to reproduce. A writer that becomes an undeclared
    input of itself moves `undeclared`; one that stops reading a file it declares moves
    `declared_not_read`; one that reads a different number of files moves `files_read`."""
    original = _traced(opens=1684, **{field: before})
    rebuilt = _traced(opens=1855, **{field: after})
    rest, run = mr.run_differences(mr.diff(mr.comparable(original), mr.comparable(rebuilt)))
    real, _ = mr.environment_differences(rest, mr.cwd_fields_recorded_absolute(original, rebuilt))
    assert [d.split(":")[0] for d in run] == [OPENS]
    assert [d.split(":")[0] for d in real] == [f"/result_manifest/traced_inputs/{field}"], real


# --------------------------------------------------------------------------------------------------
# 4. The report: what was set aside is printed, with its reason, and the two claims stay apart.
# --------------------------------------------------------------------------------------------------


def test_zero_differences_after_setting_aside_is_not_printable_as_zero_differences():
    """The claim this is all for. Those are different claims and the second must say so."""
    clean = mr.comparison_reading([], {})
    after = mr.comparison_reading([], {"run_fields_ignored": [f"{OPENS}: 1684 vs 1855"]})
    assert clean != after
    assert "nothing was set aside" in clean
    assert "SETTING ASIDE" in after and "not the same claim" in after
    assert "1 run fields" in after, after


def test_a_reading_that_set_nothing_aside_says_exactly_that():
    assert mr.comparison_reading([], {}) == (
        "0 differences, and nothing was set aside: every leaf compared reconciles"
    )


def test_differences_and_set_aside_leaves_are_both_counted():
    reading = mr.comparison_reading(
        ["/verdict/x: 1 vs 2"], {"environment_fields_ignored": ["/code_cleanliness/dirty: True vs False"]}
    )
    assert reading.startswith("1 differences")
    assert "1 further leaves were set aside" in reading


def test_the_run_field_exemption_is_declared_lossy_wherever_it_is_reported():
    """The correction that produced this wording: `opens` is NOT noise. In the one case measured it moved
    because the writer's own behaviour changed. A report that leaves a reader to infer the exemption is
    free is the defect, so the loss is named in the one sentence a reader reads and in the reason."""
    reading = mr.comparison_reading([], {"run_fields_ignored": [f"{OPENS}: 1684 vs 1855"]})
    assert "KNOWN-LOSSY" in reading
    assert "hide a genuine change" in reading
    why = mr.IGNORED_BECAUSE["run_fields_ignored"]
    assert "KNOWN-LOSSY" in why
    assert "HIDE A GENUINE CHANGE" in why
    assert "git blob" in why, "the measured case must be named, not generalised away"


def test_no_claim_is_made_that_the_opens_count_is_interpreter_noise():
    """An unverified premise, withdrawn by the lane that raised it. Nothing here may rest on it, and a
    later edit that reintroduces it as the justification must fail this test."""
    why = mr.IGNORED_BECAUSE["run_fields_ignored"]
    source = Path("scripts/manifest_rebuild.py").read_text()
    for claim in ("interpreter incidental", "library incidental", "moves with the Python build"):
        assert claim not in why
        assert claim not in source, f"{claim!r} is the withdrawn premise"
    assert "too coarse" in why, "what IS claimed: the count is too coarse to act on"


def test_a_lossy_exemption_is_announced_even_when_the_comparison_also_failed():
    """A reader who sees differences must still be told what was set aside, or the failure list reads as
    the whole of what was looked at."""
    reading = mr.comparison_reading(["/verdict/x: 1 vs 2"], {"run_fields_ignored": [f"{OPENS}: 1 vs 2"]})
    assert "KNOWN-LOSSY" in reading


def test_every_class_of_set_aside_leaf_carries_a_printed_reason():
    """An ignored field must be visible and auditable, never invisible."""
    for key in (
        "environment_fields_ignored",
        "run_fields_ignored",
        "timing_fields_ignored",
        "resource_fields_ignored",
    ):
        assert key in mr.IGNORED_BECAUSE
        assert len(mr.IGNORED_BECAUSE[key]) > 60, f"{key} needs a reason, not a label"
    assert mr.CWD_FIELDS[0] in mr.IGNORED_BECAUSE["environment_fields_ignored"]
    assert mr.RUN_FIELDS[0] in mr.IGNORED_BECAUSE["run_fields_ignored"]


def test_the_environment_reason_says_a_relative_cwd_difference_is_real():
    why = mr.IGNORED_BECAUSE["environment_fields_ignored"]
    assert "REPOSITORY-RELATIVE" in why and "real difference" in why


def test_the_existing_environment_exemption_still_works_unchanged():
    """`environment_differences` gained a second argument; its old single-argument behaviour must not
    have moved, because every other caller and every other test relies on it."""
    real, env = mr.environment_differences(
        ["/code_cleanliness/dirty: True vs False", "/verdict/positives: 12 vs 11"]
    )
    assert real == ["/verdict/positives: 12 vs 11"]
    assert env == ["/code_cleanliness/dirty: True vs False"]
