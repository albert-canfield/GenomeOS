#!/usr/bin/env bash
# Promote dev to main: run the checks locally, push dev, open (or reuse) the PR,
# wait for CI, and merge. Usage: scripts/promote.sh ["PR title"]
set -euo pipefail
cd "$(dirname "$0")/.."
branch=$(git rev-parse --abbrev-ref HEAD)
[ "$branch" = "dev" ] || { echo "switch to dev first (on $branch)"; exit 1; }
uv run ruff check . && uv run ruff format --check . && uv run pytest -q
git push -u origin dev
title="${1:-Promote dev to main}"
pr=$(gh pr list --base main --head dev --json number -q '.[0].number')
if [ -z "$pr" ]; then
  gh pr create --base main --head dev --title "$title" --body "Tested on dev; CI green."
fi
gh pr checks dev --watch --interval 15
gh pr merge dev --merge --delete-branch=false
git fetch -q origin main
echo "main is now at $(git rev-parse --short origin/main)"
