#!/usr/bin/env bash
# Shared checkout: nothing reaches origin/dev without the checks CI runs.
# The commit being pushed is checked out into a temporary worktree, so another
# session's half-edited files in this checkout neither block nor pass the push.
# Skip once with GENOMEOS_SKIP_CHECK=1 only for a documentation-only push.
#
# One push's check at a time. On 2026-10-02 three pushes failed with a missing
# temporary directory: several sessions push the same checkout, each adding and
# removing a worktree under one shared TMPDIR, and the checkout still carries
# prunable registrations from runs that were interrupted mid-check. The lock
# below lets the pushes queue instead of interleaving, and the trap removes this
# push's own worktree and its own lock and nothing else, so an interrupted push
# leaves neither behind and no push can take a directory from another.
set -u
[ "${GENOMEOS_SKIP_CHECK:-0}" = "1" ] && exit 0
root=$(git rev-parse --show-toplevel)
lock="${TMPDIR:-/tmp}/genomeos-push.lock"
tmp=""
held=""

cleanup() {
  if [ -n "$tmp" ]; then
    git -C "$root" worktree remove --force "$tmp" 2>/dev/null
    rm -rf "$tmp"
  fi
  # only this push's lock: the pid file says whose it is
  if [ "$(cat "$lock/pid" 2>/dev/null || echo)" = "$$" ]; then rm -rf "$lock"; fi
}
trap cleanup EXIT INT TERM

waited=0
while ! mkdir "$lock" 2>/dev/null; do
  held=$(cat "$lock/pid" 2>/dev/null || echo)
  if [ -n "$held" ] && ! kill -0 "$held" 2>/dev/null; then
    echo "pre-push: the lock was left by process $held, which is gone; taking it"
    rm -rf "$lock"
    continue
  fi
  if [ -z "$held" ] && [ "$waited" -ge 60 ]; then
    echo "pre-push: the lock names no process after ${waited}s; taking it"
    rm -rf "$lock"
    continue
  fi
  if [ "$waited" -ge 3600 ]; then
    echo "pre-push: process ${held:-unknown} has held the push lock for an hour; not waiting longer"
    exit 1
  fi
  [ "$waited" -eq 0 ] && echo "pre-push: waiting for the check running in process ${held:-unknown}"
  sleep 5
  waited=$((waited + 5))
done
echo "$$" > "$lock/pid"

status=0
while read -r local_ref local_sha remote_ref remote_sha; do
  case "$remote_ref" in refs/heads/dev|refs/heads/main) ;; *) continue ;; esac
  [ "$local_sha" = "0000000000000000000000000000000000000000" ] && continue
  tmp=$(mktemp -d "${TMPDIR:-/tmp}/genomeos-push.XXXXXX")
  echo "pre-push: checking $(git rev-parse --short "$local_sha") in a clean worktree"
  if git -C "$root" worktree add -q --detach "$tmp" "$local_sha"; then
    # the local caches are shared read-only, so real-data tests run as they do here
    for d in reference knowledge cache; do
      [ -d "$root/data/$d" ] && [ ! -e "$tmp/data/$d" ] && ln -s "$root/data/$d" "$tmp/data/$d"
    done
    (cd "$tmp" && UV_NO_SYNC=1 UV_PROJECT_ENVIRONMENT="$root/.venv" PYTHONPATH="$tmp" scripts/check.sh) || status=1
    git -C "$root" worktree remove --force "$tmp"
  else
    status=1
  fi
  tmp=""
done
[ "$status" -eq 0 ] || echo "pre-push: red; fix and push again (see CONTRIBUTING.md)"
exit "$status"
