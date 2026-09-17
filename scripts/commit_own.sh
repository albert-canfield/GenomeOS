#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Commit your own files from a private index, with the guard unskippable.
#
#     scripts/commit_own.sh -F msg.txt genomeos/compare.py tests/test_compare.py
#     scripts/commit_own.sh -F msg.txt --force docs/ROADMAP.md
#
# Several sessions share this checkout, so the documented sequence is: private index, read-tree HEAD,
# add explicit paths, run the stale-base guard, commit-tree, update-ref pinning the old value. Every
# step of that is in LESSONS.md and in scripts/check_staged.py, and on 2026-09-17 the session that
# wrote both of those documents broke two of the steps in one day:
#
#   - it piped the guard into `tail`, so a REFUSED exited 0 and the `&&` meant to stop the commit ran
#     it instead. Third occurrence, same session, warning in bold in the tool it wrote;
#   - it called `git update-ref refs/heads/dev "$NEW"` with no old value, every time. A peer commit
#     landing between `rev-parse HEAD` and `update-ref` would have been dropped without a word. That
#     never happened, which is luck: the window is seconds wide and there were peers committing.
#
# So the sequence lives here instead of in anyone's fingers. The guard is run unpiped and
# unredirected, its exit status is checked directly, and update-ref always pins the parent it read.
set -euo pipefail

force=""
msg_file=""
paths=()
while [ $# -gt 0 ]; do
  case "$1" in
    -F) msg_file="$2"; shift 2 ;;
    --force) force="--force"; shift ;;
    -*) echo "unknown flag: $1" >&2; exit 64 ;;
    *) paths+=("$1"); shift ;;
  esac
done

[ -n "$msg_file" ] || { echo "a commit message file is required: -F msg.txt" >&2; exit 64; }
[ -f "$msg_file" ] || { echo "no such message file: $msg_file" >&2; exit 64; }
[ ${#paths[@]} -gt 0 ] || { echo "name the paths to commit; this never adds -A" >&2; exit 64; }

branch=$(git rev-parse --abbrev-ref HEAD)
[ "$branch" != "main" ] || { echo "refusing to commit to main; work on dev" >&2; exit 65; }

parent=$(git rev-parse HEAD)
index=$(mktemp -t "genomeos-index")
trap 'rm -f "$index"' EXIT
export GIT_INDEX_FILE="$index"

git read-tree "$parent"
git add -- "${paths[@]}"

# Unpiped and unredirected on purpose: a pipeline's status is its last command's, and a redirect
# throws away the reason. Both of those have already cost this project a commit.
set +e
python3 scripts/check_staged.py $force
guard=$?
set -e
if [ $guard -ne 0 ]; then
  echo "commit_own: the guard refused (exit $guard). Look at what it named, then re-run with --force" >&2
  exit "$guard"
fi

tree=$(git write-tree)
commit=$(git commit-tree "$tree" -p "$parent" -F "$msg_file")
# the old value is pinned: if a peer moved the branch since `rev-parse HEAD`, this fails rather than
# dropping their commit. Re-run after rebasing your paths onto the new tip.
git update-ref "refs/heads/$branch" "$commit" "$parent"
echo "commit_own: $(git log --oneline -1 "$branch")"
