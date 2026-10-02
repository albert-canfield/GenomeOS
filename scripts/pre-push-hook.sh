#!/usr/bin/env bash
# THE INSTALLED HOOK, and it is a FIXED WRAPPER that never needs re-copying.
#
# WHY IT EXISTS. `scripts/install-hooks.sh` used to COPY scripts/pre-push.sh into .git/hooks/pre-push,
# so the hook was not the script: an edit to the script changed nothing the push ran until somebody
# remembered to reinstall. On 2026-10-02 that cost a day. The push's verification worktree got its
# git-ignored data stores as SYMLINKS, because the installed copy was taken at 04:34 and still did
# `ln -s "$root/data/$d" "$tmp/data/$d"`, while scripts/pre-push.sh had moved at 17:52 to read-only
# APFS clones. `git check-ignore` cannot answer past a symbolic link, 20 tests ERRORED at setup, and
# two RED status files were written for trees that were nothing of the kind -- one of them refusing a
# push and leaving origin/dev thirty commits behind. The clone switch was correct for its own reason
# (kernel-mode protection of the inputs, which no audit hook can give) and could not have fixed this
# and did not: it never reached what the push ran.
#
# WHAT IT DOES, in the order it matters.
#
# 1. It RUNS THE COMMITTED SCRIPT -- `git show HEAD:scripts/pre-push.sh` -- and not the working-tree
#    copy. In a checkout several sessions edit at once, the working copy is somebody's draft about
#    half the time; the commit being pushed is the thing under verification, so the verifier comes
#    from the commit too.
# 2. It REFUSES BY NAME when the working-tree copy differs from HEAD's, because an uncommitted
#    verification script gives no verdict. That is the same rule scripts/manifest_rebuild.py applies
#    to itself (`tool_is_committed`), for the same reason: a tool nobody can name cannot be audited,
#    and a green it produced cannot be reproduced. So an edit to pre-push.sh takes effect WHEN IT IS
#    COMMITTED, and never silently not at all.
# 3. Steps 1 and 2 are redundant on purpose. Once they agree, HEAD's bytes and the working copy's are
#    the same bytes; extracting from HEAD anyway closes the gap between the comparison and the run,
#    in which a peer in this shared checkout really can save the file.
#
# WHAT IT DOES NOT CLOSE, stated rather than left to be discovered. A COMMIT CAN STILL WEAKEN ITS OWN
# PRE-PUSH CHECK: the hook runs the pre-push.sh of the commit being pushed, so a commit that gutted it
# would be verified by the gutted version. That case is covered by the commit review, which reads the
# diff -- it is a visible change in a tracked file with a message attached. What this wrapper closes is
# the SILENT case, where nobody changed anything and the push quietly ran last week's checker. That is
# the one that bit us, and it is the one no review could have caught.
#
# This file's own bytes are checked: tests/test_pre_push_hook_is_the_wrapper.py asserts that the
# installed hook IS this file, and names `scripts/install-hooks.sh` when it is not.
set -euo pipefail

root=$(git rev-parse --show-toplevel)
cd "$root"

script=scripts/pre-push.sh
committed=$(git rev-parse --quiet --verify "HEAD:$script" 2>/dev/null || true)
if [ -z "$committed" ]; then
  echo "pre-push: REFUSED: $script is not in HEAD, so there is no committed checker to run." >&2
  echo "pre-push:   commit it, or push with GENOMEOS_SKIP_CHECK=1 only for a documentation-only push." >&2
  exit 1
fi

if [ ! -f "$script" ]; then
  echo "pre-push: REFUSED: $script is in HEAD but absent from the working tree." >&2
  exit 1
fi
working=$(git hash-object -- "$script")
if [ "$committed" != "$working" ]; then
  echo "pre-push: REFUSED: $script in the working tree differs from HEAD's copy." >&2
  echo "pre-push:   HEAD $committed, working tree $working" >&2
  echo "pre-push:   An uncommitted verification script gives no verdict: commit the change (or stage" >&2
  echo "pre-push:   nothing of it and restore the file) and push again. This is the same rule" >&2
  echo "pre-push:   scripts/manifest_rebuild.py applies to itself." >&2
  echo "pre-push:   In this shared checkout the edit may be a PEER'S -- ask before touching it." >&2
  exit 1
fi

# Extracted to a file rather than piped into bash, because the pre-push protocol sends its four
# fields per ref ON STDIN and the checker reads them: a `git show ... | bash` would eat that stdin.
tmp=$(mktemp "${TMPDIR:-/tmp}/genomeos-prepush-checker.XXXXXX")
trap 'rm -f "$tmp"' EXIT
git show "HEAD:$script" > "$tmp"
bash "$tmp" "$@"
