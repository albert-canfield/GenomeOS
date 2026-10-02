#!/usr/bin/env bash
# Install the pre-push check into this clone (see CONTRIBUTING.md).
#
# WHAT IS INSTALLED IS A FIXED WRAPPER, scripts/pre-push-hook.sh, and NOT the checker itself. The
# wrapper runs `git show HEAD:scripts/pre-push.sh` and refuses when the working copy differs from
# HEAD's, so an edit to the checker takes effect when it is committed and this script never has to be
# run again. Until 2026-10-02 this copied scripts/pre-push.sh directly, and the copy in .git/hooks
# went stale: a push ran a checker from 04:34 that symlinked the data stores into its verification
# worktree, `git check-ignore` could not answer past the symlink, and two RED status files were
# written for trees that were nothing of the kind. scripts/pre-push-hook.sh says that at length.
set -euo pipefail
cd "$(dirname "$0")/.."
cp scripts/pre-push-hook.sh .git/hooks/pre-push
chmod +x .git/hooks/pre-push
echo "pre-push hook installed: a fixed wrapper that runs HEAD:scripts/pre-push.sh"
