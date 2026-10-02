#!/usr/bin/env bash
# ONE FULL TEST SUITE AT A TIME, MACHINE-WIDE. Sourced by scripts/check.sh.
#
# WHY. On 2026-10-02 the acceptance rule in force was "a status-file verdict of the committed
# tree", and every lane satisfied it by running the whole ~4,500-test suite in a worktree of its
# own commit. Three ran at once. The machine began paging -- swapouts rising 90,000 to 114,000 per
# sample, 14 GB of swap on disk, free disk under the 10 GB floor -- and the capacity gate then
# refused the heavy jobs that real work needed, which is how a hygiene rule came to hold a release.
# Two things came out of that: the acceptance rule itself was narrowed (CONTRIBUTING.md, "A lane's
# acceptance is its targeted tests"), and the runs that remain QUEUE instead of interleaving.
#
# THIS IS THE PUSH LOCK'S ALGORITHM, DELIBERATELY. scripts/pre-push.sh already serialises one
# push's check at a time at "${TMPDIR:-/tmp}/genomeos-push.lock", and tests/test_push_lock.py pins
# its two cases. A second lock in the same project that behaved differently would be worse than no
# second lock, so the wording, the takeover rule and the re-announcement below are that lock's:
#   - a waiter RE-PRINTS THE HOLDER WHEN IT CHANGES, so somebody watching learns who they are
#     waiting for instead of reading one line and then silence for half an hour;
#   - a DEAD holder is taken over, because a crashed run must not block the machine for an hour;
#   - a LIVE holder is NEVER stolen. On the same day a stale log line naming a holder that had
#     since exited was read as this branch failing, and the proposed remedy -- deleting the lock by
#     hand -- would have taken a live peer's lock.
# tests/test_suite_lock.py plants both, through the real scripts/check.sh, with its own TMPDIR.
#
# THE LOCK IS NEVER A RED. Every path here returns 0. A check that exits non-zero because it could
# not get a lock writes a RED STATUS FILE ABOUT NOTHING IN THE TREE, which is the exact fault this
# project spent 2026-10-02 undoing: five reds that day were tooling faults and each read like a test
# failure. So when the lock cannot be had -- an unusable TMPDIR, or a holder that has sat there past
# the cap -- this says so loudly on stderr AND hands check.sh a note for the verdict's `note` field,
# and the run proceeds unserialised. A run that quietly skipped its own serialisation would make the
# verdict mean less than its reader thinks, and that is the failure this project keeps writing up.
#
# WHAT IT DOES NOT FIX. The status file is keyed by TREE ALONE, so two checks of the SAME tree still
# collide and the later writer still destroys the earlier verdict. Serialising makes that RARER --
# the two runs no longer overlap -- but not impossible: the second run still writes the same path
# the moment the first has released. Nothing here should be read as closing it.

# shellcheck shell=bash

#: The lock this process holds, empty when it holds none. Read by suite_lock_release.
SUITE_LOCK=""

#: Set when the run is going ahead WITHOUT the lock, for the verdict's note. Consumed by
#: scripts/check.sh's combined_note; assigned here, which is why shellcheck cannot see its use.
# shellcheck disable=SC2034
SUITE_LOCK_NOTE=""

suite_lock_acquire() {
  local label="${1:-a check}"
  local lock parent held announced waited cap
  lock="${TMPDIR:-/tmp}/genomeos-suite.lock"
  cap="${GENOMEOS_SUITE_LOCK_WAIT:-5400}"

  # RE-ENTRANCY, WHICH IS NOT A BYPASS. tests/test_check_sh_arguments.py runs the real
  # scripts/check.sh inside the suite, so a project-wide run reaches this script from inside a
  # process that already holds the lock. Waiting there would deadlock until the cap. The outer run
  # already serialises everything it starts, so a nested one needs no lock of its own. The marker
  # must name a LIVE process: a variable left behind in somebody's environment cannot stand in for
  # the lock, and a dead pid falls through to the real acquisition below.
  held="${GENOMEOS_SUITE_LOCK_HELD:-}"
  if [ -n "$held" ] && kill -0 "$held" 2>/dev/null; then
    echo "suite lock: process $held holds it and started this run; not waiting for its own lock"
    SUITE_LOCK=""
    return 0
  fi

  parent=$(dirname "$lock")
  if [ ! -d "$parent" ] && ! mkdir -p "$parent" 2>/dev/null; then
    SUITE_LOCK=""
    SUITE_LOCK_NOTE="not serialised: $parent does not exist and could not be created, so this run did not queue behind another full suite"
    echo "check: WARNING $SUITE_LOCK_NOTE" >&2
    return 0
  fi

  announced=""
  waited=0
  while ! mkdir "$lock" 2>/dev/null; do
    # mkdir failed and there is nothing there: the path itself is unusable (a file in the way, a
    # read-only parent). Looping would spin to the cap and then say the wrong thing about a holder.
    if [ ! -d "$lock" ]; then
      SUITE_LOCK=""
      SUITE_LOCK_NOTE="not serialised: $lock could not be created, so this run did not queue behind another full suite"
      echo "check: WARNING $SUITE_LOCK_NOTE" >&2
      return 0
    fi
    held=$(cat "$lock/pid" 2>/dev/null || echo)
    if [ -n "$held" ] && ! kill -0 "$held" 2>/dev/null; then
      echo "suite lock: the lock was left by process $held, which is gone; taking it"
      rm -rf "$lock"
      continue
    fi
    # A directory with no pid file is a lock caught between mkdir and the write below. Sixty
    # seconds is far longer than that window and far shorter than a suite.
    if [ -z "$held" ] && [ "$waited" -ge 60 ]; then
      echo "suite lock: the lock names no process after ${waited}s; taking it"
      rm -rf "$lock"
      continue
    fi
    if [ "$waited" -ge "$cap" ]; then
      SUITE_LOCK=""
      SUITE_LOCK_NOTE="not serialised: process ${held:-unknown} held the suite lock for ${waited}s, past the ${cap}s cap, and this run went ahead beside it"
      echo "check: WARNING $SUITE_LOCK_NOTE" >&2
      return 0
    fi
    if [ "$waited" -eq 0 ] || [ "$held" != "$announced" ]; then
      echo "suite lock: waiting for the suite running in process ${held:-unknown} ($label)"
      announced="$held"
    fi
    sleep 5
    waited=$((waited + 5))
  done
  printf '%s\n' "$$" > "$lock/pid"
  SUITE_LOCK="$lock"
  # Inherited by anything this run starts, including a nested scripts/check.sh; see RE-ENTRANCY.
  export GENOMEOS_SUITE_LOCK_HELD="$$"
  [ "$waited" -eq 0 ] || echo "suite lock: taken after ${waited}s"
  return 0
}

# Only this run's lock, decided by the pid file and never by the path, for the reason
# scripts/pre-push.sh gives: a release that removes "the lock" removes a live peer's.
suite_lock_release() {
  [ -n "${SUITE_LOCK:-}" ] || return 0
  if [ "$(cat "$SUITE_LOCK/pid" 2>/dev/null || echo)" = "$$" ]; then rm -rf "$SUITE_LOCK"; fi
  SUITE_LOCK=""
  return 0
}
