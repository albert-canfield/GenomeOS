#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Fast-forward main to a dev sha, only if CI's `test` job passed on that exact sha.
#
# Usage: scripts/promote_main.sh [--dry-run | --push] [--remote NAME] SHA
#        scripts/promote_main.sh [--dry-run | --push] [--remote NAME] --latest-green
# Usage: scripts/promote_main.sh [--dry-run | --prepare] [--remote NAME] SHA
#        scripts/promote_main.sh [--dry-run] [--remote NAME] --latest-green
#        scripts/promote_main.sh --verify [--remote NAME] SHA
# (2026-09-29: the three lines above replace the two before them, whose --push is retired and refused.)
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
#
# 2026-09-29, beside the original: the sentence above is false. With enforce_admins on, `--push 9a59faf`
# passed every check here and GitHub refused it (GH006, "Required status check \"test\" is expected"),
# although GitHub Actions had recorded a successful `test` check run on that exact sha: a manually
# triggered (workflow_dispatch) run on dev. A successful manually triggered check is not enough to
# establish that GitHub will accept the push, so passing this gate is necessary, not sufficient. Every
# earlier promotion went through as an admin, whom the rule did not bind. Promotion now goes through a
# pull request from a branch frozen at the chosen sha, merged by the owner.
#
# 2026-09-29, later, beside both notes above (lane-gate): the gate now reads the event of the `test`
# run it relies on and keeps two findings apart. "CI pre-check passed": that run concluded success.
# "Eligible for protected promotion": that run was also triggered by `pull_request` or `push`. GitHub's
# documentation ("Troubleshooting required status checks", section "Checks from some workflow jobs are
# not evaluated") evaluates checks for pull requests and rulesets only from runs triggered by push,
# pull_request, pull_request_review, pull_request_target, deployment or deployment_status; that the
# same rule decided the refused direct push of 9a59faf is an inference, not established. A green
# `workflow_dispatch` or `schedule` run on dev is therefore a pre-check only; eligibility comes from
# the pull request's own run.
#
# --push is retired, and refused before the arguments are read: sessions never push to main. Its case
# in the argument loop, the first usage() definition, usage lines 6-7 and the push block at the end
# stay, unreachable or unused, only because they fall inside the checkout guard's two-day window for
# b50c873; they can be removed after 2026-09-30 23:03. In its place, --prepare SHA makes the same
# checks, then pushes one branch, promote-<first 7 characters of SHA>, at exactly that sha (refusing if
# the branch already exists at another sha), and prints the compare URL for the owner to open and
# merge the pull request, and the --verify command for after the merge. It never opens a pull request,
# never merges and never touches main. --verify SHA fetches main and checks that its file tree is
# SHA's and that SHA is its ancestor; it asks nothing of GitHub's API.
set -euo pipefail

usage() { sed -n '6,7p' "$0" | sed 's/^# //' >&2; }
# 2026-09-29: this definition replaces the one above, which prints the retired --push usage;
# bash keeps the last definition.
usage() { sed -n '8,10p' "$0" | sed 's/^# //' >&2; }
refuse() {
  echo "promote_main: refused: $*" >&2
  exit 1
}

mode=dry
remote=origin
target=""
latest=0
# 2026-09-29: --push is refused here, before the loop below can read it; and one mode flag at most
nmodes=0
for arg in "$@"; do
  case "$arg" in
    --push)
      echo "promote_main: --push is retired: GitHub refused a direct push under the bound rule; use --prepare SHA and open a pull request" >&2
      exit 2
      ;;
    --dry-run | --prepare | --verify) nmodes=$((nmodes + 1)) ;;
  esac
done
if [ "$nmodes" -gt 1 ]; then
  echo "promote_main: one of --dry-run, --prepare, --verify" >&2
  exit 2
fi
while [ $# -gt 0 ]; do
  case "$1" in
    --dry-run) mode=dry ;;
    # unreachable since 2026-09-29 (refused above); remove after 2026-09-30 23:03
    --push) mode=push ;;
    --remote)
      [ $# -ge 2 ] || { usage; exit 2; }
      remote="$2"
      shift
      ;;
    --prepare) mode=prepare ;;
    --verify) mode=verify ;;
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
if [ "$mode" != dry ] && [ "$latest" = 1 ]; then
  echo "promote_main: --$mode takes the sha the owner named, not --latest-green" >&2
  exit 2
fi

if [ "$mode" = verify ]; then
  # after the owner's merge: main carries exactly the chosen sha's files, and the sha itself
  git fetch -q "$remote" main || refuse "cannot fetch main from $remote"
  main_sha=$(git rev-parse --verify "refs/remotes/$remote/main^{commit}")
  sha=$(git rev-parse --verify --quiet "$target^{commit}") || refuse "$target is not a commit in this clone"
  short=$(git rev-parse --short "$sha")
  main_short=$(git rev-parse --short "$main_sha")
  if [ "$(git rev-parse "$main_sha^{tree}")" != "$(git rev-parse "$sha^{tree}")" ]; then
    echo "promote_main: not verified: $remote/main ($main_short) does not have $short's file tree" >&2
    exit 1
  fi
  if ! git merge-base --is-ancestor "$sha" "$main_sha"; then
    echo "promote_main: not verified: $short is not an ancestor of $remote/main ($main_short)" >&2
    exit 1
  fi
  echo "promote_main: verified: $remote/main ($main_short) has $short's file tree, and $short is its ancestor"
  exit 0
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

# 2026-09-29: the event that triggered that same run (same selection as above, joined on its check suite)
wruns=$(gh api "repos/$repo/actions/runs?head_sha=$sha&per_page=100") ||
  refuse "cannot read the workflow runs of $short from $repo"
event=$(CHECKS="$checks" python3 -c '
import json, os, sys
runs = [r for r in json.loads(os.environ["CHECKS"]).get("check_runs", [])
        if (r.get("app") or {}).get("slug") == "github-actions" and r.get("status") == "completed"]
runs.sort(key=lambda r: r.get("completed_at") or "")
suite = (runs[-1].get("check_suite") or {}).get("id")
events = {w.get("check_suite_id"): w.get("event") for w in json.load(sys.stdin).get("workflow_runs", [])}
print((suite is not None and events.get(suite)) or "unknown")
' <<<"$wruns") || refuse "cannot read the workflow runs of $short: GitHub's answer did not parse"

ahead=$(git rev-list --count "$main_sha..$sha")
echo "promote_main: $remote/main $(git rev-parse --short "$main_sha") -> $short, $ahead commits; test green at $completed ($url)"
case "$event" in
  pull_request | push)
    echo "promote_main: eligible for protected promotion: the \`test\` run relied on was triggered by $event (GitHub's documentation evaluates checks for pull requests and rulesets from runs of this event)"
    ;;
  *)
    echo "promote_main: CI pre-check passed: the \`test\` run relied on was triggered by $event, a pre-check only; it does not make the sha eligible for protected promotion"
    echo "promote_main: eligibility comes from the pull request's own run. GitHub's documentation evaluates checks for pull requests and rulesets only from runs triggered by push, pull_request, pull_request_review, pull_request_target, deployment or deployment_status; that the same rule refused the direct push of 9a59faf is an inference, not established"
    ;;
esac
if [ "$mode" = dry ]; then
  echo "promote_main: dry run; would run: git push $remote $sha:refs/heads/main"
  echo "promote_main: that push is no longer made: --push is retired (2026-09-29). Instead, on the owner's go, scripts/promote_main.sh --prepare $sha pushes branch promote-${sha:0:7} at that sha and nothing else, for a pull request the owner opens and merges"
  exit 0
fi
if [ "$mode" = prepare ]; then
  # one branch frozen at the chosen sha; never a pull request, never a merge, never main
  branch="promote-${sha:0:7}"
  at_remote() { git ls-remote "$remote" "refs/heads/$branch" | awk -v ref="refs/heads/$branch" '$2 == ref { print $1 }'; }
  existing=$(at_remote) || refuse "cannot read the branches of $remote"
  if [ -z "$existing" ]; then
    git push -q "$remote" "$sha:refs/heads/$branch" || refuse "the push of $branch to $remote failed"
  elif [ "$existing" = "$sha" ]; then
    echo "promote_main: $remote already has $branch at $short; nothing pushed"
  else
    refuse "$remote already has $branch at ${existing:0:7}, not at $short; nothing pushed"
  fi
  [ "$(at_remote)" = "$sha" ] || refuse "$remote's $branch is not at $short"
  echo "promote_main: $remote/$branch is at $short; main was not touched and no pull request was opened"
  echo "promote_main: for the owner, to open and merge the pull request: https://github.com/$repo/compare/main...$branch?expand=1"
  echo "promote_main: after the merge: scripts/promote_main.sh --verify $sha   (main's file tree is $short's, and $short is an ancestor of main)"
  exit 0
fi

# Unreachable since 2026-09-29: every mode has exited above, and --push is refused before the
# arguments are read. The four lines after this refusal stay only because they fall inside the
# checkout guard's window for b50c873; remove them, and this refusal, after 2026-09-30 23:03.
refuse "no mode pushes to main"
git push "$remote" "$sha:refs/heads/main"
git fetch -q "$remote" main
[ "$(git rev-parse "refs/remotes/$remote/main")" = "$sha" ] || refuse "pushed, but $remote/main is not at $short"
echo "promote_main: $remote/main is at $short"
