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
uv run pytest -q
uv run bio test data/demo genomeos/std data/organisms
echo "check: green"
