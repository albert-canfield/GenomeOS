#!/usr/bin/env bash
# The check every commit to dev must pass, in the order CI runs it:
# lint, format, the test suite, and the BioLang programs testing themselves.
# Usage: scripts/check.sh            (whole project, what CI runs)
#        scripts/check.sh FILE...    (lint and format only the given files, then the full tests)
set -euo pipefail
cd "$(dirname "$0")/.."
if [ "$#" -gt 0 ]; then
  uv run ruff check "$@"
  uv run ruff format --check "$@"
else
  uv run ruff check .
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

uv run pytest -q
uv run bio test data/demo genomeos/std data/organisms
if [ -n "$shellcheck_note" ]; then
  echo "check: green EXCEPT $shellcheck_note"
else
  echo "check: green"
fi
