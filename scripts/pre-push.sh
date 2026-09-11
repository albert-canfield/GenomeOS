#!/usr/bin/env bash
# Shared checkout: nothing reaches origin/dev without the checks CI runs.
# The commit being pushed is checked out into a temporary worktree, so another
# session's half-edited files in this checkout neither block nor pass the push.
# Skip once with GENOMEOS_SKIP_CHECK=1 only for a documentation-only push.
set -u
[ "${GENOMEOS_SKIP_CHECK:-0}" = "1" ] && exit 0
root=$(git rev-parse --show-toplevel)
status=0
while read -r local_ref local_sha remote_ref remote_sha; do
  case "$remote_ref" in refs/heads/dev|refs/heads/main) ;; *) continue ;; esac
  [ "$local_sha" = "0000000000000000000000000000000000000000" ] && continue
  tmp=$(mktemp -d "${TMPDIR:-/tmp}/genomeos-push.XXXXXX")
  echo "pre-push: checking $(git rev-parse --short "$local_sha") in a clean worktree"
  if git -C "$root" worktree add -q --detach "$tmp" "$local_sha"; then
    (cd "$tmp" && UV_NO_SYNC=1 UV_PROJECT_ENVIRONMENT="$root/.venv" PYTHONPATH="$tmp" scripts/check.sh) || status=1
    git -C "$root" worktree remove --force "$tmp"
  else
    status=1
  fi
done
[ "$status" -eq 0 ] || echo "pre-push: red; fix and push again (see CONTRIBUTING.md)"
exit "$status"
