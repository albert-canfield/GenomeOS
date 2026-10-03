#!/usr/bin/env bash
# The check every commit to dev must pass, in the order CI runs it:
# lint, format, the test suite, and the BioLang programs testing themselves.
# Usage: scripts/check.sh            (whole project, what CI runs)
#        scripts/check.sh FILE...    (lint and format only the given .py files, then the full tests)
#
# ONE OF THESE RUNS AT A TIME ON THE MACHINE, AND BOTH FORMS COUNT. `check.sh FILE...` narrows the
# LINT ONLY -- the pytest leg below is the whole suite either way -- so it is not a light run and it
# is not a lane's acceptance. On 2026-10-02 three of these ran at once because the acceptance rule
# then in force asked every lane for a project-wide verdict: the machine paged (swapouts rising
# 90,000 to 114,000 per sample, 14 GB of swap on disk, free disk under the 10 GB floor), the capacity
# gate began refusing the heavy jobs real work needed, and a hygiene rule ended up holding a release.
# So runs QUEUE here, by the same lock scripts/pre-push.sh uses for one push's check at a time, and
# CONTRIBUTING.md now asks a lane for its TARGETED tests instead. See scripts/suite_lock.sh.
#
# A NON-PYTHON ARGUMENT IS REFUSED BY NAME AND DOES NOT KILL THE RUN. ruff is a Python linter and it
# reads whatever path it is handed as Python: on 2026-10-02 a lane passed
# data/results/manifest_headlines.json, ruff reported 96 errors in it, `set -e` stopped the run at the
# ruff-check leg with no tests run at all, and a RED status file was written for the tree. Twice, from
# nothing but an argument -- status-b17ef4d0b2800f21f64e34c4dbee62ea707cc180 and
# status-f8a9cf884bdf2f665f9dde009400d718749dc560. Both named a .json among their scope_files and
# neither was about the code. So the arguments are triaged first: .py and .pyi go to ruff, everything
# else is named in a refusal line and in the verdict's note, and the tests run either way. The
# refusal is not a leg and cannot be red: an argument is the caller's mistake, not the tree's.
#
# THE VERDICT IS A FILE, NOT AN EXIT CODE. Every run ends by writing a status file that names the
# TREE IT JUDGED, its exit code and its pass/fail/error/skip counts; genomeos/verdict.py holds the
# reasons at length and the rule for which tree hash is the right one to record. The short version:
# four times in one night an exit code lied -- through a pipe, through `tail`, through a background
# launch -- and a fifth time an honest exit code was about a tree that had moved underneath the run.
# scripts/pre-push.sh and scripts/push_own.sh therefore read the file and compare its tree against
# the tree they are about to push, and refuse when it is a different tree or when there is no file
# at all. A missing verdict is not a pass.
#
# Nothing about what this script RUNS has changed, and nothing about the codes it exits with: the
# legs are the same four in the same order, `set -e` still stops at the first red, and the exit code
# is still that leg's. The status file is written from an EXIT trap so that a red run -- which never
# reaches the last line -- still records its verdict instead of leaving the previous run's in place
# to be mistaken for this one's.
set -euo pipefail
cd "$(dirname "$0")/.."

# Taken HERE, before started_at and tree_begin are measured: a wait of half an hour must not be
# reported as this run's duration, and must not widen the window in which the tree can move
# underneath the run. The lock is never a red -- every path in suite_lock_acquire returns 0 and an
# unobtainable lock becomes a loud line plus a note in the verdict -- because a check that exits
# non-zero over a lock writes a red status file about nothing in the tree.
# shellcheck source=scripts/suite_lock.sh
. scripts/suite_lock.sh
suite_lock_acquire "${*:-whole project}"
# The real trap is installed further down, once finish() exists; until then this one is here so a
# run killed in the next few lines does not leave its lock for the next run to reclaim.
trap suite_lock_release EXIT

started_at=$(date +%s)
scope=project
scope_args=()
lint_args=()
refused_args=()
if [ "$#" -gt 0 ]; then
  scope=files
  for f in "$@"; do
    scope_args+=(--scope-file "$f")
    case "$f" in
      *.py | *.pyi) lint_args+=("$f") ;;
      *) refused_args+=("$f") ;;
    esac
  done
fi

# Said in a line of its own AND carried into the verdict, for the same reason the shellcheck note is:
# a check that quietly drops part of what it was asked to do makes its own green mean less than the
# reader thinks. The caller is told what ruff was NOT given and why, by name.
refusal_note=""
if [ "${#refused_args[@]}" -gt 0 ]; then
  refusal_note="not linted, ruff reads a path as Python and these are not: ${refused_args[*]}"
  echo "check: REFUSED these lint arguments by name: ${refused_args[*]}" >&2
  echo "check: ruff is a Python linter; a .json handed to it is reported as broken Python, not as" >&2
  echo "check: broken JSON. For a JSON file: python3 -c 'import json,sys; json.load(open(sys.argv[1]))' FILE" >&2
  echo "check: the run continues and the tests below are unaffected; this is not a verdict on the tree" >&2
fi
if [ "$#" -gt 0 ] && [ "${#lint_args[@]}" -eq 0 ]; then
  echo "check: no .py argument was given, so NOTHING was linted in this run" >&2
fi

# Measured before the legs and measured again by `write` afterwards. Two readings, because one
# cannot tell a run that judged one tree from a run the tree moved underneath.
tree_begin=$(python3 -m genomeos.verdict tree)
status_dir=$(python3 -m genomeos.verdict status-dir)
mkdir -p "$status_dir"
pytest_log="$status_dir/pytest-$$.log"
leg=startup

# Every command here is kept from failing the trap: `set -e` is in force inside a trap too, so a
# write that goes wrong would otherwise replace the real exit code with its own and the caller
# would be told the wrong thing by the very code meant to stop that happening.
# The verdict's note carries every reason this run is less than it looks: a run that was not
# serialised against the other full suites on the machine, a shellcheck that could not run, and a
# lint argument refused by name. Joined here rather than in the `write` call, because an unset
# variable inside a nested expansion inside a trap is how a trap starts replacing real exit codes
# with its own.
combined_note() {
  local note=""
  local extra
  for extra in "${SUITE_LOCK_NOTE:-}" "${shellcheck_note:-}" "${refusal_note:-}" "${skip_report_note:-}"; do
    [ -n "$extra" ] || continue
    if [ -n "$note" ]; then note="$note; $extra"; else note="$extra"; fi
  done
  printf '%s' "$note"
}

finish() {
  local code=$1
  local failed_leg=""
  [ "$code" -eq 0 ] || failed_leg=$leg
  python3 -m genomeos.verdict write \
    --exit-code "$code" \
    --tree-begin "$tree_begin" \
    --scope "$scope" \
    ${scope_args[@]+"${scope_args[@]}"} \
    --started-at "$started_at" \
    --pytest-log "$pytest_log" \
    --failed-leg "$failed_leg" \
    --note "$(combined_note)" ||
    echo "check: WARNING the verdict file could not be written; treat this run as having no verdict" >&2
  # After the verdict is written, never before: the next run in the queue should find the status
  # file this one earned. Only this process's own lock is removed; suite_lock_release decides that
  # from the pid file, so an interrupted peer's lock is left where it is.
  suite_lock_release
  return 0
}
# Replaces the lock-only trap set above; bash keeps one EXIT trap, and this one does both.
trap 'finish "$?"' EXIT

# `ruff check` with no paths lints the whole project, so an empty argument list must SKIP the legs
# rather than fall through to them: a file-scoped run that quietly became a project lint is the same
# class of lie as a green that linted nothing.
leg=ruff-check
if [ "$#" -eq 0 ]; then
  uv run ruff check .
  leg=ruff-format
  uv run ruff format --check .
elif [ "${#lint_args[@]}" -gt 0 ]; then
  uv run ruff check "${lint_args[@]}"
  leg=ruff-format
  uv run ruff format --check "${lint_args[@]}"
fi
# The shell scripts carry the commit and push guards, and until 2026-10-02 nothing linted them:
# ruff is given the file list, and a .sh path named there is read as Python and fails as one. Two
# routes, because shellcheck is not installed on every machine and installing it is not this
# script's business. The pinned PyPI wheel is transient and deliberately NOT in pyproject:
# its licence is GPL-3.0 and this project takes permissive licences only, so it stays a dev tool
# that is fetched when needed and never shipped or locked. (A comment must not BEGIN with the
# tool's name: it reads one as a directive, which is how this very leg failed on its first run.)
#
# What matters more than either route: when neither runs, this says so in a line of its own AND in
# the verdict. A check that quietly skips is worse than no check, because the green line then means
# less than the reader thinks it does, and that is the exact failure this project keeps writing up.
leg=shellcheck
shell_scripts=(scripts/*.sh)
shellcheck_note=""
if command -v shellcheck >/dev/null 2>&1; then
  shellcheck "${shell_scripts[@]}"
elif uv run --quiet --with shellcheck-py==0.11.0.1 --no-project -- shellcheck "${shell_scripts[@]}"; then
  :
else
  status=$?
  # 1 is shellcheck's "I found something"; anything else means it never ran (offline, no index)
  if [ $status -eq 1 ]; then exit 1; fi
  shellcheck_note="shellcheck not installed and the pinned wheel could not be fetched: ${#shell_scripts[@]} shell scripts not linted"
  echo "check: $shellcheck_note" >&2
fi

# Teed, not redirected: the counts in the verdict are read back from this log, and the run is still
# watchable while it happens. PIPESTATUS, because `tee`'s exit code is not pytest's -- the whole
# family of faults this verdict file exists for begins with somebody reading the wrong member of a
# pipeline. (pipefail is on as well, so `set -e` stops here on red either way.)
#
# PYTHONUNBUFFERED, because a pipe makes the interpreter's stdout block-buffered and the run then
# shows nothing for minutes at a time. A check that looks stalled is what gets backgrounded, and a
# backgrounded check is where the exit code stopped meaning anything. Buffering only: the same
# command, the same arguments, the same codes.
# A COUNT OF SKIPS IS NOT A VERDICT ON A TEST, and until 2026-10-03 this leg ran bare `pytest -q`,
# so the log a status file points at recorded "12 skipped" and could not say WHICH twelve. push19's
# verdict for tree c78041b1 read "5,053 passed, 12 skipped, 8 xfailed", and no requirement of the
# form "the signed tree must show test X PASSED, not skipped" could be settled from it: somebody
# re-ran 260 tests by hand that night to establish what the kept file should already have carried.
# CI has passed `-rs` all along (.github/workflows/ci.yml), so this also ends an asymmetry in which
# the throwaway CI log said more about the tree than the verdict file this project signs from.
#
# `-rs` ALONE IS NOT ENOUGH, which is the part worth knowing: pytest FOLDS its skip summary by
# (file, line, reason) and prints `SKIPPED [3] tests/test_x.py:12: reason` -- a count and a location,
# no test id, and three parametrised cases collapsed into one line. `--no-fold-skipped` (pytest 8.3
# and later) prints one line per test, `SKIPPED tests/test_x.py::test_y - Skipped: reason`, which is
# the form a reader can match against a requirement that names a test.
#
# PROBED, NOT ASSUMED. The lock pins pytest 9.1.1 but the dev pin is `pytest>=8`, so a checkout can
# hold an 8.0-8.2 that does not know the flag -- and pytest answers an unknown flag with a USAGE
# ERROR, exit 4, NO TEST RUN AT ALL. That is the exact shape of the two false reds of 2026-10-02,
# where an argument and not the tree made a status file red. So support is read off `pytest --help`
# and the flag is dropped where it is absent, with the loss said in a line of its own AND in the
# verdict's note: the shellcheck leg's rule, for its reason -- a check that quietly carries less than
# its reader thinks is worse than one that says so.
#
# NOTHING ELSE HERE MOVES. These flags change what the log PRINTS ABOUT SKIPS and nothing else: the
# summary line genomeos/verdict.py reads is byte-identical, `parse_pytest_counts` is untouched, and
# the exit code is still pytest's own. tests/test_check_skip_reasons.py plants a real skipping test,
# runs this script, and reads the id out of the log it wrote.
leg=pytest
skip_report_args=(-rs)
skip_report_note=""
# Captured into a variable rather than piped into grep: with pipefail on, a `grep -q` that closes the
# pipe early can leave the pipeline reporting uv's SIGPIPE, and the probe would then answer
# "unsupported" about a pytest that supports the flag perfectly well.
pytest_help=$(uv run pytest --help 2>/dev/null || true)
case "$pytest_help" in
*--no-fold-skipped*)
  skip_report_args+=(--no-fold-skipped)
  ;;
*)
  skip_report_note="skips reported FOLDED: this pytest has no --no-fold-skipped, so the log names each skip's file, line and reason but NOT its test id"
  echo "check: $skip_report_note" >&2
  ;;
esac
set +e
PYTHONUNBUFFERED=1 uv run pytest -q "${skip_report_args[@]}" 2>&1 | tee "$pytest_log"
pytest_code=${PIPESTATUS[0]}
set -e
[ "$pytest_code" -eq 0 ] || exit "$pytest_code"

leg=biolang
uv run bio test data/demo genomeos/std data/organisms
leg=report
report_note=$(combined_note)
if [ -n "$report_note" ]; then
  echo "check: green EXCEPT $report_note"
else
  echo "check: green"
fi
