# SPDX-License-Identifier: AGPL-3.0-or-later
"""Proving a resource guard judges the figure it was HANDED and not the live machine.

WHY THIS EXISTS (2026-10-03, lane-rssceiling). Four guards in this project compared a ceiling
against `resource.getrusage(resource.RUSAGE_SELF).ru_maxrss`, read live, from inside a test. That
figure is the PROCESS high-water mark and it is MONOTONIC: measured on this machine, 15.6 MB before
allocating 300 MB, 330.2 MB after, and still 330.2 MB after freeing it and collecting. A pytest
session is one process, so by the time a late test reads it the figure is every earlier test's peak
and none of the late test's own. In a full-suite run under the suite lock the process had already
peaked at 4,406 MB, and ceilings of 900 MB and 4 GiB were breached by tests that allocated none of
it: 4 failed, 5,101 passed, and all 58 tests in the two files passed when run alone.

THE FIX IS INJECTION, which is the pattern a disk guard in this tree already uses for the same
class of defect (a money guard that read the live volume and reddened the suite whenever free disk
was low). The guard takes the figure as a parameter; a real run passes nothing and gets the live
reading, a test passes the figure it means to test.

NOT A DELTA. A before-and-after difference of `RUSAGE_SELF` is NOT a fix and is not offered here: if
an earlier test peaked higher, the high-water does not move, the delta reads 0 and the guard passes
vacuously. That is a guard that cannot fire, which is worse than the false red it replaces.

WHAT THIS MODULE ADDS over asserting the two obvious plants by hand. A guard that accepts the
parameter and then reads the live figure anyway passes both obvious plants on most machines: a 5 GB
injection "refuses" because the live figure happens to be over the ceiling too, and a 100 MB
injection "passes" because the live figure happens to be under it. The third plant below is the one
that separates wired from assumed-wired, and on the disk lane it cost exactly one line to get
wrong. `wiring_problems` runs all three and is itself checked against a deliberately unwired copy.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

#: Every exception type a ceiling guard in this project raises when its ceiling is breached. Caught
#: together on purpose: a plant must never pass because the guard refused with a type it did not
#: expect, which would read as "no refusal" and hide an unwired guard.
REFUSALS = (MemoryError, RuntimeError, AssertionError)

#: A guard adapted to one signature: `(ceiling, figure) -> the figure judged`, raising on a breach.
#: The caller writes the adapter, so a guard taking a positional `where` or a differently named
#: keyword is still testable without this module knowing anything about it.
Guard = Callable[[int, "int | None"], Any]


def refusal(guard: Guard, ceiling: int, figure: int | None) -> str | None:
    """The text of the guard's refusal, or None when it returned a figure instead of refusing."""
    try:
        guard(ceiling, figure)
    except REFUSALS as e:
        return str(e)
    return None


def answer(guard: Guard, ceiling: int, figure: int | None) -> Any:
    """What the guard returned, or None when it refused."""
    try:
        return guard(ceiling, figure)
    except REFUSALS:
        return None


def wiring_problems(guard: Guard, *, ceiling: int, over: int, under: int, live: int) -> list[str]:
    """The three plants, as a list of what failed. Empty means the guard judges what it was handed.

    `ceiling` is the guard's own registered ceiling, `over` a figure above it, `under` a figure
    below it, and `live` the guard's own live reading of the machine right now.

    1. OVER THE CEILING REFUSES, and the refusal NAMES the injected figure. A guard that ignores the
       injection refuses with the machine's figure, so reading the number back out of the message is
       what distinguishes the two; "it raised" on its own does not.
    2. UNDER THE CEILING PASSES, and returns exactly the injected figure.
    3. THE LIVE READING IS NOT CONSULTED. The ceiling is set one unit under the live figure and 0 is
       injected. A wired guard judges 0 and passes; a guard that reads the machine sees `live` over
       `live - 1` and refuses. This plant is decided by the live figure alone, so it discriminates
       identically whether the suite has peaked high or the file is run on its own.
    """
    problems: list[str] = []
    if over <= ceiling or under >= ceiling:
        raise ValueError(f"the plants need over > {ceiling} > under, got over={over} under={under}")
    if live < 2:
        raise ValueError(f"a live reading of {live} cannot decide plant 3; it is never this small")

    said = refusal(guard, ceiling, over)
    if said is None:
        problems.append(f"plant 1: an injected {over} over the ceiling {ceiling} did not refuse")
    elif str(over) not in said:
        problems.append(f"plant 1: the refusal does not name the injected {over}: {said!r}")

    got = answer(guard, ceiling, under)
    if got != under:
        problems.append(f"plant 2: an injected {under} under the ceiling {ceiling} returned {got!r}")

    unconsulted = answer(guard, live - 1, 0)
    if unconsulted != 0:
        problems.append(
            f"plant 3: with the ceiling at {live - 1}, one under the live reading {live}, an "
            f"injected 0 returned {unconsulted!r}; the guard is reading the machine, not the figure"
        )
    return problems


def unwired(live_reading: Callable[[], int]) -> Guard:
    """The defect, as a guard: it accepts the figure, and then reads `live_reading()` anyway.

    This is the one-line mistake `wiring_problems` exists to catch, built here so each caller can
    show its own plants FAIL against it. A plant suite that cannot fail proves nothing about the
    guard it passes for.
    """

    def guard(ceiling: int, figure: int | None) -> int:
        peak = live_reading()
        if peak >= ceiling:
            raise RuntimeError(f"peak {peak} reached the ceiling of {ceiling}")
        return peak

    return guard
