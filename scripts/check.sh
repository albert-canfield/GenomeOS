#!/usr/bin/env bash
# The check every commit to dev must pass, in the order CI runs it:
# lint, format, the test suite, and the BioLang programs testing themselves.
# Usage: scripts/check.sh            (whole project, what CI runs)
#        scripts/check.sh FILE...    (lint and format only the given files, then the full tests)
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

started_at=$(date +%s)
scope=project
scope_args=()
if [ "$#" -gt 0 ]; then
  scope=files
  for f in "$@"; do scope_args+=(--scope-file "$f"); done
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
    --note "${shellcheck_note:-}" ||
    echo "check: WARNING the verdict file could not be written; treat this run as having no verdict" >&2
  return 0
}
trap 'finish "$?"' EXIT

leg=ruff-check
if [ "$#" -gt 0 ]; then
  uv run ruff check "$@"
  leg=ruff-format
  uv run ruff format --check "$@"
else
  uv run ruff check .
  leg=ruff-format
  uv run ruff format --check .
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
leg=pytest
set +e
PYTHONUNBUFFERED=1 uv run pytest -q 2>&1 | tee "$pytest_log"
pytest_code=${PIPESTATUS[0]}
set -e
[ "$pytest_code" -eq 0 ] || exit "$pytest_code"

leg=biolang
uv run bio test data/demo genomeos/std data/organisms
leg=report
if [ -n "$shellcheck_note" ]; then
  echo "check: green EXCEPT $shellcheck_note"
else
  echo "check: green"
fi
