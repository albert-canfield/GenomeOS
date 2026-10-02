#!/usr/bin/env bash
# Shared checkout: nothing reaches origin/dev without the checks CI runs.
# The commit being pushed is checked out into a temporary worktree, so another
# session's half-edited files in this checkout neither block nor pass the push.
# Skip once with GENOMEOS_SKIP_CHECK=1 only for a documentation-only push.
#
# THE VERDICT IS READ FROM THE FILE THE CHECK WROTE, NOT FROM ITS EXIT CODE. The exit
# code is still required to be 0 -- that refusal is untouched -- but on its own it is not
# evidence: on 2026-10-02 a `scripts/check.sh ... | tail` reported 0 over four real
# failures, and a backgrounded check returned 0 while its log said `19 failed`, because a
# pipeline's `$?` belongs to the last command and a background launch's to the launcher.
# Worse, an exit code says nothing about WHICH TREE it judged, and this check takes nine
# to twelve minutes in a checkout several sessions commit to. So below, the status file
# written by scripts/check.sh must exist, must name THE TREE OF THE COMMIT BEING PUSHED in
# both its before and after readings, must be project-scoped, and must say green. No file
# means no verdict, and a missing verdict is not a pass.
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
announced=""

# shellcheck disable=SC2329  # invoked by the trap below, which shellcheck cannot see
cleanup() {
  if [ -n "$tmp" ]; then
    # the linked stores are at 0444 and their directories at 0555, so `rm -rf` is refused without this.
    # A shell chmod rather than the linking script: a trap must not depend on uv being runnable, and an
    # interrupted push would otherwise leave a read-only tree nothing could take away. Every entry under
    # there is a clone of its own, so this touches nothing the machine keeps.
    chmod -R u+w "$tmp" 2>/dev/null || true
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
  if [ "$waited" -eq 0 ] || [ "$held" != "$announced" ]; then
    echo "pre-push: waiting for the check running in process ${held:-unknown}"
    announced="$held"
  fi
  sleep 5
  waited=$((waited + 5))
done
echo "$$" > "$lock/pid"

status=0
# shellcheck disable=SC2034  # git's pre-push protocol sends four fields per ref on stdin and they
# are positional: local_ref and remote_sha are read because they are there, not because we use them
while read -r local_ref local_sha remote_ref remote_sha; do
  case "$remote_ref" in refs/heads/dev|refs/heads/main) ;; *) continue ;; esac
  [ "$local_sha" = "0000000000000000000000000000000000000000" ] && continue
  tmp=$(mktemp -d "${TMPDIR:-/tmp}/genomeos-push.XXXXXX")
  # The tree of the commit being pushed, which is what the verdict must name. A clean worktree
  # of this commit hashes to exactly this, so the comparison is an identity and not a tolerance:
  # if the link lines below ever materialise something .gitignore does not cover, the hashes stop
  # matching and this refuses, naming the paths the two trees differ in, rather than quietly
  # allowing a difference nobody looked at.
  expect_tree=$(git -C "$root" rev-parse "$local_sha^{tree}")
  echo "pre-push: checking $(git rev-parse --short "$local_sha") in a clean worktree (tree $expect_tree)"
  if git -C "$root" worktree add -q --detach "$tmp" "$local_sha"; then
    # The git-ignored stores are put here READ-ONLY -- APFS clones at mode 0444, not symlinks. Until
    # 2026-10-02 these three lines were `ln -s`, and a symlink carries its TARGET's mode: the stores'
    # directories are 0755, so every push verdict this project has ever taken was produced in a worktree
    # whose input stores were WRITABLE. Measured, not reasoned:
    # `open("$tmp/data/reference/HG002_chr1.vcf.gz", "r+b")` through the link succeeded, and a probe file
    # created through such a link appeared in the real store. A verification run that can modify the
    # inputs it verifies against can produce a clean verdict by changing what it compared to.
    #
    # A clone is a separate inode, so taking the write bits off takes them off this worktree's side only
    # and never off the machine's one copy. Measured cost: 11.2 GB of stores, 47,716 entries, about 9
    # seconds and 17 MB of disk -- a clone that silently became a real copy would show as minutes and
    # gigabytes, which is why the script prints the figures rather than assuming them.
    #
    # Tree identity is preserved and CHECKED rather than assumed: .gitignore names the three stores
    # WITHOUT a trailing slash, which matches a real directory as well as a symlink, so the worktree
    # still hashes to the committed tree. `--check-tree` prints the hash before and after, and
    # tests/test_link_stores_read_only.py asserts they are equal.
    if ! (cd "$root" && UV_NO_SYNC=1 UV_PROJECT_ENVIRONMENT="$root/.venv" \
          uv run python scripts/link_stores_read_only.py "$tmp" --check-tree); then
      echo "pre-push: REFUSED: the data stores could not be linked read-only into $tmp" >&2
      status=1
    fi
    (cd "$tmp" && UV_NO_SYNC=1 UV_PROJECT_ENVIRONMENT="$root/.venv" PYTHONPATH="$tmp" scripts/check.sh) || status=1
    # And now the verdict itself, from the file. Run from the clean worktree so the reader is the
    # committed one and not whatever a peer has half-edited in the shared checkout; the status
    # directory is the common .git either way, so both see the same files. No pipe anywhere here:
    # `v=$(cmd)` keeps the command's own exit status, which `cmd | tail` does not, and BOTH the
    # status and the printed token are required, so a shell that reads only one is still safe.
    # If the reader cannot run at all it exits non-zero and this refuses: it fails closed.
    if v=$(cd "$tmp" && python3 -m genomeos.verdict require --tree "$expect_tree" --scope project 2>&1); then
      case "$v" in
        "GENOMEOS_VERDICT_OK $expect_tree "*) echo "pre-push: verdict read from the status file: $v" ;;
        *)
          echo "pre-push: REFUSED: the verdict reader exited 0 without naming this tree:" >&2
          echo "$v" >&2
          status=1
          ;;
      esac
    else
      echo "$v" >&2
      echo "pre-push: REFUSED: no green verdict exists for the tree being pushed ($expect_tree)." >&2
      status=1
    fi
    # Wait for anything still running IN the worktree before taking the directory away. On 2026-10-02
    # a peer's `pkill -f "scripts/check.sh"` killed the check here; this line removed the worktree at
    # once, and `bio test` -- a child of the killed check, still running -- lost its working directory.
    # The push then reported red with pytest already green at 3,485 passed, on an absent TMPDIR path
    # and then an absent RELATIVE path to a tracked, present file. The kill was the fault; removing the
    # directory under a live child is what made it a puzzle instead of a message.
    waited_for_children=0
    while pgrep -f "^bash $tmp/scripts/check.sh" >/dev/null 2>&1 || pgrep -f "PYTHONPATH=$tmp" >/dev/null 2>&1; do
      [ "$waited_for_children" -lt 60 ] || { echo "pre-push: something is still running in $tmp after ${waited_for_children}s; removing anyway" >&2; break; }
      sleep 2
      waited_for_children=$((waited_for_children + 2))
    done
    chmod -R u+w "$tmp" 2>/dev/null || true  # the linked stores are read-only; see cleanup()
    git -C "$root" worktree remove --force "$tmp"
  else
    status=1
  fi
  tmp=""
done
[ "$status" -eq 0 ] || echo "pre-push: red; fix and push again (see CONTRIBUTING.md)"
exit "$status"
