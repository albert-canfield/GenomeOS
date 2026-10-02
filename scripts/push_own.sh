#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Push the sha that was checked, and recover from the one race this checkout keeps losing.
#
#     scripts/push_own.sh                 # pushes HEAD to origin/dev
#     scripts/push_own.sh --remote NAME --sha SHA
#
# Two problems, both seen on 2026-10-02 within an hour of each other.
#
# FIRST: `git push` names a BRANCH, and the pre-push check takes nine minutes. In a checkout several
# sessions commit to, local dev can move between the check starting and the push happening, so the
# verdict is about one tree and the push carries another. Pushing an explicit sha fixes that: the
# content that reaches the remote is exactly the content the check passed.
#
# SECOND: git sends the remote-tracking ref as the expected old value, and a peer pushing during
# those nine minutes makes it stale, so the remote refuses with "cannot lock ref ... is at X but
# expected Y". Both a lane and the coordinator hit this, and both diagnosed it by hand.
#
# The recovery: fetch, and retry ONLY if the remote's new tip is an ancestor of the sha being pushed,
# which makes the retry a plain fast-forward carrying exactly the checked tree.
#
# What that ancestor test does and does NOT do, because the first version of this comment claimed
# more than the test supports. It is NOT what stops a peer's commits being discarded: the retry is an
# ordinary `git push` with no --force, and git refuses a non-fast-forward by itself. Removing the
# ancestor test and running the diverged case proves it -- the push still fails and the peer's commit
# is still on the remote. What the test actually buys is a diagnosis instead of a second confusing
# rejection, a refusal before another round trip, and a written rule at the place where someone would
# otherwise be tempted to add --force to make the retry "work". That is worth having; it is not a
# safety property, and calling it one would have been a guard credited with someone else's work.
#
# The retry skips the check, and that is safe for exactly one reason: it is the SAME sha this script
# has just watched pass, in this invocation. The skip is never set from a flag, never read from the
# environment, and never applied to a sha this run did not check -- so it cannot become a way to push
# something unchecked. If the first attempt's check went red, there is no retry at all.
#
# "Watched pass" USED TO MEAN an exit code and a grep of captured output, and both of those have been
# caught lying in one night: a guard piped into `tail` exited 0 over a REFUSED, and a backgrounded
# check returned 0 while its log said `19 failed`. So the retry no longer rests on having watched
# anything. It requires the status file scripts/check.sh writes, and requires it to name THE TREE OF
# THE SHA BEING PUSHED -- which is the only thing that makes "the same sha" a fact rather than an
# assumption in a checkout where local dev moves during a nine-minute check. The same requirement is
# applied BEFORE the first push whenever GENOMEOS_SKIP_CHECK is already set in the environment, since
# that is the one way to reach the remote with no check having run at all.
set -euo pipefail

remote=origin
branch=dev
sha=""
while [ $# -gt 0 ]; do
  case "$1" in
    --remote) remote="$2"; shift 2 ;;
    --sha) sha="$2"; shift 2 ;;
    --branch) branch="$2"; shift 2 ;;
    -*) echo "unknown flag: $1" >&2; exit 64 ;;
    *) echo "unexpected argument: $1" >&2; exit 64 ;;
  esac
done
[ "$branch" != "main" ] || { echo "push_own: refusing to push main; promotion is scripts/promote_main.sh" >&2; exit 65; }
[ -n "$sha" ] || sha=$(git rev-parse HEAD)
sha=$(git rev-parse "$sha")
short=$(git rev-parse --short "$sha")
tree=$(git rev-parse "$sha^{tree}")

# The verdict, read from the file scripts/check.sh writes and never from an exit code. No pipe:
# `v=$(cmd)` keeps the command's own status, and both the status and the printed token are required,
# so a reader that consults only one of them is still safe. It fails closed -- if the reader cannot
# run, there is no verdict, and a missing verdict is not a pass.
require_verdict() {
  local why=$1 v
  if v=$(python3 -m genomeos.verdict require --tree "$tree" --scope project 2>&1); then
    case "$v" in
      "GENOMEOS_VERDICT_OK $tree "*)
        echo "push_own: $why: the status file's verdict is for this exact tree: $v"
        return 0
        ;;
    esac
    echo "push_own: REFUSED ($why): the verdict reader exited 0 without naming tree $tree:" >&2
    echo "$v" >&2
    return 1
  fi
  echo "$v" >&2
  echo "push_own: REFUSED ($why): no green project-wide verdict exists for $short's tree $tree." >&2
  echo "  Run scripts/check.sh on that tree; a verdict for another tree does not answer for this one." >&2
  return 1
}

if [ "${GENOMEOS_SKIP_CHECK:-0}" = "1" ]; then
  echo "push_own: GENOMEOS_SKIP_CHECK is set in the environment, so scripts/pre-push.sh will not check"
  echo "  anything. The verdict must then already exist for this tree, in writing."
  require_verdict "the check is being skipped" || exit 1
fi

attempt() {
  # unredirected on the terminal AND captured, so the check's output is visible while it runs and
  # still readable afterwards: a guard whose reason is thrown away has cost this project a commit
  set +e
  "$@" 2>&1 | tee "$out"
  local st=${PIPESTATUS[0]}
  set -e
  return "$st"
}

out=$(mktemp -t genomeos-push)
trap 'rm -f "$out"' EXIT

echo "push_own: pushing $short to $remote/$branch (the sha, not the branch name)"
if attempt git push "$remote" "$sha:refs/heads/$branch"; then
  echo "push_own: $short is on $remote/$branch"
  exit 0
fi

if ! grep -q 'cannot lock ref' "$out"; then
  echo "push_own: the push failed for a reason that is not the ref-lock race; nothing retried." >&2
  echo "  Read the output above. A red check is a red check: fix it and push again." >&2
  exit 1
fi
if grep -q 'pre-push: red' "$out"; then
  echo "push_own: the check was red as well as the ref being stale; nothing retried." >&2
  exit 1
fi

echo "push_own: the ref was stale (a peer pushed during the check). Fetching to see whether this is a fast-forward."
git fetch -q "$remote" "$branch"
tip=$(git rev-parse "refs/remotes/$remote/$branch")
if ! git merge-base --is-ancestor "$tip" "$sha"; then
  echo "push_own: REFUSED to retry: $remote/$branch is at $(git rev-parse --short "$tip"), which is NOT an" >&2
  echo "  ancestor of $short. The histories have diverged, so a retry would not be a fast-forward." >&2
  echo "  git would refuse it anyway; this stops first so you get a reason instead of a second" >&2
  echo "  rejection. Rebase your work onto $remote/$branch and run the check again." >&2
  echo "  Never --force here: forcing is the one thing that WOULD discard their commits." >&2
  exit 1
fi
echo "push_own: $remote/$branch is at $(git rev-parse --short "$tip"), an ancestor of $short, so the retry is a"
echo "  fast-forward carrying exactly the tree the check passed. Retrying without re-running the check."
# The retry skips the check, so the claim in the line above is now required in writing rather than
# inferred from the first attempt's exit code and output.
require_verdict "the retry skips the check" || exit 1
GENOMEOS_SKIP_CHECK=1 git push "$remote" "$sha:refs/heads/$branch"
echo "push_own: $short is on $remote/$branch"
