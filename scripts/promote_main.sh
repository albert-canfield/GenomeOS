#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Fast-forward main to a dev sha, only if CI's `test` job passed on that exact sha.
#
# Usage: scripts/promote_main.sh [--dry-run | --push] [--remote NAME] SHA
#        scripts/promote_main.sh [--dry-run | --push] [--remote NAME] --latest-green
#
# A dry run is the default: it makes every check and prints the push it would make. --push makes
# it, and is for the coordinator on the owner's go, at most once a day (CONTRIBUTING.md, "Release
# checks"). --latest-green picks the newest sha between main and the tip of dev that has a green
# CI run; the daily scheduled run on dev is what usually provides one, and
# `gh workflow run ci.yml --ref dev` asks for one on dev's current tip.
#
# Refuses, changing nothing, when the sha is not on the remote's dev, when main is not its
# ancestor (not a fast-forward), or when the most recent completed `test` check run that GitHub
# Actions recorded on that sha is anything but success: failure, cancelled, skipped, or none at
# all. It reads the same check that main's branch rule names as required, so what it lets through
# is what the rule would let through for someone the rule binds.
set -euo pipefail

usage() { sed -n '6,7p' "$0" | sed 's/^# //' >&2; }
refuse() {
  echo "promote_main: refused: $*" >&2
  exit 1
}

mode=dry
remote=origin
target=""
latest=0
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) mode=dry ;;
    --push) mode=push ;;
    --remote)
      [ $# -ge 2 ] || { usage; exit 2; }
      remote="$2"
      shift
      ;;
    --latest-green) latest=1 ;;
    -h | --help)
      usage
      exit 0
      ;;
    -*)
      echo "promote_main: unknown option $1" >&2
      usage
      exit 2
      ;;
    *)
      [ -z "$target" ] || { echo "promote_main: one sha only" >&2; exit 2; }
      target="$1"
      ;;
  esac
  shift
done
if [ -z "$target" ] && [ "$latest" = 0 ]; then usage; exit 2; fi
if [ -n "$target" ] && [ "$latest" = 1 ]; then
  echo "promote_main: a sha or --latest-green, not both" >&2
  exit 2
fi

# remote-tracking refs only; the working tree and the shared index are not touched
git fetch -q "$remote" dev main || refuse "cannot fetch dev and main from $remote"
main_sha=$(git rev-parse --verify "refs/remotes/$remote/main^{commit}")
dev_sha=$(git rev-parse --verify "refs/remotes/$remote/dev^{commit}")
repo="${PROMOTE_REPO:-$(gh repo view --json nameWithOwner -q .nameWithOwner)}" ||
  refuse "cannot tell which GitHub repository $remote is"

between() { # true when $1 is on dev and main is its ancestor
  git merge-base --is-ancestor "$1" "$dev_sha" && git merge-base --is-ancestor "$main_sha" "$1"
}

if [ "$latest" = 1 ]; then
  runs=$(gh api "repos/$repo/actions/workflows/ci.yml/runs?branch=dev&status=success&per_page=50") ||
    refuse "cannot list CI runs of $repo"
  target=""
  for s in $(python3 -c 'import json,sys; print(*[r["head_sha"] for r in json.load(sys.stdin)["workflow_runs"]])' <<<"$runs"); do
    if git cat-file -e "$s^{commit}" 2>/dev/null && between "$s"; then
      target="$s"
      break
    fi
  done
  [ -n "$target" ] || refuse "no sha between $remote/main and $remote/dev has a green CI run; ask for one with: gh workflow run ci.yml --ref dev"
fi

sha=$(git rev-parse --verify --quiet "$target^{commit}") || refuse "$target is not a commit in this clone"
short=$(git rev-parse --short "$sha")
if [ "$sha" = "$main_sha" ]; then
  echo "promote_main: $remote/main is already at $short; nothing to do"
  exit 0
fi
git merge-base --is-ancestor "$sha" "$dev_sha" || refuse "$short is not on $remote/dev"
git merge-base --is-ancestor "$main_sha" "$sha" ||
  refuse "$short does not contain $remote/main ($(git rev-parse --short "$main_sha")); moving main there is not a fast-forward"

checks=$(gh api "repos/$repo/commits/$sha/check-runs?check_name=test&filter=all&per_page=100") ||
  refuse "cannot read the check runs of $short from $repo"
verdict=$(python3 -c '
import json, sys
runs = [r for r in json.load(sys.stdin).get("check_runs", [])
        if (r.get("app") or {}).get("slug") == "github-actions" and r.get("status") == "completed"]
runs.sort(key=lambda r: r.get("completed_at") or "")
if not runs:
    print("none")
else:
    r = runs[-1]
    print(r.get("conclusion"), r.get("completed_at"), r.get("html_url") or "", len(runs))
' <<<"$checks") || refuse "cannot read the check runs of $short: GitHub's answer did not parse"
read -r conclusion completed url count <<<"$verdict"
[ "$conclusion" != "none" ] || refuse "$short has no completed \`test\` run; ask for one with: gh workflow run ci.yml --ref dev (it tests dev's tip at that moment)"
[ "$conclusion" = "success" ] || refuse "the latest \`test\` run on $short ($completed, $count on record) concluded $conclusion: $url"

ahead=$(git rev-list --count "$main_sha..$sha")
echo "promote_main: $remote/main $(git rev-parse --short "$main_sha") -> $short, $ahead commits; test green at $completed ($url)"
if [ "$mode" = dry ]; then
  echo "promote_main: dry run; would run: git push $remote $sha:refs/heads/main"
  exit 0
fi
git push "$remote" "$sha:refs/heads/main"
git fetch -q "$remote" main
[ "$(git rev-parse "refs/remotes/$remote/main")" = "$sha" ] || refuse "pushed, but $remote/main is not at $short"
echo "promote_main: $remote/main is at $short"
