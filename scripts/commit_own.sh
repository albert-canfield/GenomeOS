#!/usr/bin/env bash
# SPDX-License-Identifier: AGPL-3.0-or-later
#
# Commit your own files from a private index, with the guard unskippable.
#
#     scripts/commit_own.sh -F msg-<lane>-<purpose>-$(date +%s).txt genomeos/compare.py tests/...
#     scripts/commit_own.sh -F msg-coord-roadmap-1759400000.txt --force docs/ROADMAP.md
#
# The message file's name is checked, not decoration: it must carry the committing lane and the
# second it was written (-L <lane>, or GENOMEOS_LANE, says which lane you are). Five commits have
# gone in under another lane's message from this shared scratchpad; see the three checks below.
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
lane="${GENOMEOS_LANE:-}"
paths=()
while [ $# -gt 0 ]; do
  case "$1" in
    -F) msg_file="$2"; shift 2 ;;
    -L) lane="$2"; shift 2 ;;
    --force) force="--force"; shift ;;
    -*) echo "unknown flag: $1" >&2; exit 64 ;;
    *) paths+=("$1"); shift ;;
  esac
done

[ -n "$msg_file" ] || { echo "a commit message file is required: -F msg.txt" >&2; exit 64; }
[ -f "$msg_file" ] || { echo "no such message file: $msg_file" >&2; exit 64; }
[ ${#paths[@]} -gt 0 ] || { echo "name the paths to commit; this never adds -A" >&2; exit 64; }

# A message file left in the shared scratchpad by another lane is silently reusable, and on
# 2026-09-22 one was: a lane's own heredoc had been refused by the guard inside a compound call,
# so `-F msg1.txt` picked up the previous lane's text and 45e4f92 went in describing work it does
# not contain. Nothing in git noticed, because a commit message is never wrong to git. The check
# that went in then compared the whole message to HEAD's only, and it caught none of the four that
# followed (ec8536d, dd5a223, a3f1c04, d0bc88a): the reused text was older than HEAD every time.
#
# Three checks replace it. Each refuses on its own, with its own exit code, and each is overridable
# with --force once you have opened the file and seen that the repeat is deliberate.

# 1. The name. A file any lane can write is a file any lane can read by accident, so the name
#    carries the lane that is committing and the second it was written. The epoch is not
#    decoration and it is the strongest part of this guard: write the file as
#    msg-<lane>-<purpose>-$(date +%s).txt in the SAME call that commits, and a heredoc the hook
#    refuses can no longer be papered over, because the name the retry expands to does not exist
#    and the script stops at "no such message file" instead of reading whatever was there. Every
#    one of the five bad commits came from a name that still existed after the write was lost.
base=$(basename "$msg_file")
case "$base" in
  msg-*-[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]*.txt) ;;
  *)
    echo "REFUSED: $base is not a lane-unique message name." >&2
    echo "  Use msg-<lane>-<purpose>-\$(date +%s).txt. The scratchpad is shared between lanes, and" >&2
    echo "  a name any of them would pick is how four commits so far carried another lane's text." >&2
    [ -n "$force" ] || exit 67
    ;;
esac
if [ -n "$lane" ]; then
  case "$base" in
    *"$lane"*) ;;
    *)
      echo "REFUSED: $base does not carry the lane name '$lane'." >&2
      echo "  Name the file for the lane that is committing, or name the lane you mean with -L." >&2
      [ -n "$force" ] || exit 67
      ;;
  esac
fi

# 2. When it was written. A message written before the work it describes is a file from an earlier
#    round: d0bc88a's predated its own files by five hours. The comparison is against the OLDEST
#    uncommitted change among the named paths, and it allows the message to be half an hour older
#    than that, because writing the message and then rerunning the generator that rewrites a result
#    is a normal order of work. Thirty minutes is this script's own tolerance, not a measured one,
#    and it has a known hole: the stale file behind one of the five was only 15 minutes old, well
#    inside the slack, and check 1 is what would have stopped that one. This check catches the
#    hours-old case, nothing finer. Tightening it below about a quarter of an hour would start
#    refusing the legitimate order of work above, which is why it has not been tightened.
MSG_MAY_PREDATE_WORK_BY=1800
case "$(uname)" in
  Darwin|*BSD) mtime() { stat -f %m "$1"; }; when() { date -r "$1" '+%F %T'; } ;;
  *) mtime() { stat -c %Y "$1"; }; when() { date -d "@$1" '+%F %T'; } ;;
esac
msg_mtime=$(mtime "$msg_file")
oldest=""
while IFS= read -r changed; do
  # Written as an explicit `if` rather than `[ A ] && [ B ] || continue`. The two are equivalent
  # here, but shellcheck 0.9.0 -- which is what the ubuntu-24.04 runner ships -- flags the second as
  # SC2015, and in CI that is a fatal lint, not a note. Said positively: skip a line that is empty
  # or is not a file that exists.
  if [ -z "$changed" ] || [ ! -f "$changed" ]; then continue; fi
  t=$(mtime "$changed")
  if [ -z "$oldest" ] || [ "$t" -lt "$oldest" ]; then oldest="$t"; fi
done <<EOF
$( { git diff --name-only HEAD -- "${paths[@]}"; git ls-files --others --exclude-standard -- "${paths[@]}"; } | sort -u )
EOF
if [ -n "$oldest" ] && [ "$msg_mtime" -lt $(( oldest - MSG_MAY_PREDATE_WORK_BY )) ]; then
  echo "REFUSED: $base was written $(( (oldest - msg_mtime) / 60 )) minutes before the oldest change it describes." >&2
  echo "  message written $(when "$msg_mtime"); oldest uncommitted change $(when "$oldest")." >&2
  echo "  A message older than the work is a stale scratchpad file. Write this commit's own message." >&2
  [ -n "$force" ] || exit 68
fi

# 3. What it says. A first line already in recent history is the signature of a reused file. The
#    count is taken without a pipe into a short-circuiting reader on purpose: `grep -q` closing the
#    pipe early makes the pipeline's status 141 under pipefail, which would read as "no match" and
#    silently retire this check. Piping a guard into something that swallows its status has cost
#    this project a commit once already.
first=$(head -1 "$msg_file")
recent=$(git log -200 --pretty=%s)
dups=0
if [ -n "$first" ]; then
  dups=$(printf '%s\n' "$recent" | grep -Fxc -- "$first" || true)
fi
if [ "$dups" -gt 0 ]; then
  echo "REFUSED: this message's first line is already in the last 200 commits:" >&2
  hits=$(printf '%s\n' "$recent" | grep -Fxn -- "$first" | cut -d: -f1 | head -3 || true)
  shas=$(git log -200 --pretty=%h)
  for n in $hits; do
    echo "  $(printf '%s\n' "$shas" | sed -n "${n}p") already begins with this line" >&2
  done
  echo "  $msg_file is almost certainly another lane's file, or your own from an earlier commit:" >&2
  echo "  the scratchpad is shared between lanes. Write this commit's message to a lane-unique" >&2
  echo "  name (msg-<lane>-<purpose>-\$(date +%s).txt) and run again; --force if the repeat is" >&2
  echo "  deliberate." >&2
  [ -n "$force" ] || exit 66
fi

branch=$(git rev-parse --abbrev-ref HEAD)
[ "$branch" != "main" ] || { echo "refusing to commit to main; work on dev" >&2; exit 65; }

parent=$(git rev-parse HEAD)
# Full path with six X's, not `mktemp -t <prefix>`. GNU mktemp -- which is what CI runs --
# REFUSES a -t template with fewer than three trailing X's: `mktemp: too few X's in template
# 'genomeos-index'`, which is verbatim what the run of 2026-10-02 12:08 UTC printed. BSD mktemp, which
# is what every lane has locally, accepts a bare prefix and invents its own suffix, so this
# script worked on each developer machine and could not work on Linux at all. This form is the
# one scripts/pre-push.sh already uses, it is what GNU documents, and BSD takes it too (probed:
# it appends a further suffix of its own, which is harmless since the name is never parsed).
index=$(mktemp "${TMPDIR:-/tmp}/genomeos-index.XXXXXX")
trap 'rm -f "$index"' EXIT
export GIT_INDEX_FILE="$index"

git read-tree "$parent"
git add -- "${paths[@]}"

# Unpiped and unredirected on purpose: a pipeline's status is its last command's, and a redirect
# throws away the reason. Both of those have already cost this project a commit.
set +e
# --lane and --message-name are what this script already knows and the guard cannot find out: the
# guard refuses a staged path that a DIFFERENT live lane holds on the work board, so it has to be
# told which lane is committing. Both are optional there; without them a lane's own held file would
# read as a peer's. The message name is accepted as a second form because check 1 above already
# requires it to carry the lane, so the test that says the message is yours says the entry is yours.
python3 scripts/check_staged.py $force --lane "$lane" --message-name "$base"
guard=$?
set -e
if [ $guard -ne 0 ]; then
  echo "commit_own: the guard refused (exit $guard). Look at what it named, then re-run with --force" >&2
  # Every commit in this checkout is authored by Albert, so the author the guard prints does NOT
  # say whose lane wrote it. Read the named commit's MESSAGE to tell your own earlier work from a
  # peer's: your own is the common case (a fixture or a list you are editing again) and is safe to
  # force past; a peer's is not, and is what the check exists for.
  echo "  The author it names is Albert for every commit here: read the commit's message, not its" >&2
  echo "  author, to tell your own earlier commit from a peer's." >&2
  exit "$guard"
fi

tree=$(git write-tree)
commit=$(git commit-tree "$tree" -p "$parent" -F "$msg_file")
# the old value is pinned: if a peer moved the branch since `rev-parse HEAD`, this fails rather than
# dropping their commit. Re-run after rebasing your paths onto the new tip.
git update-ref "refs/heads/$branch" "$commit" "$parent"

# Clean up the trap this script would otherwise leave. Committing through a private index leaves the
# SHARED index still describing the old state, so `git status` shows the paths just committed as
# staged deletions and untracked files at once. They are neither — the working copies match HEAD —
# but a peer running `git add -A` would stage those deletions for real. Resetting only the paths this
# run touched fixes it without disturbing anything anyone else has staged.
unset GIT_INDEX_FILE
git reset -q HEAD -- "${paths[@]}"
echo "commit_own: $(git log --oneline -1 "$branch")"
