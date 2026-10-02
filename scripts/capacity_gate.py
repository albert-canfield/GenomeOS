# SPDX-License-Identifier: AGPL-3.0-or-later
"""Refuse a heavy job while the machine is paging, and print the figures that said so.

Several sessions share this checkout. A heavy job -- a writer re-run, a full nine-minute suite, a
reader of the per-element cache -- can make every other lane slow rather than fail, which is worse,
because nothing reports it.

Choosing the signal took three wrong answers on 2026-10-02, all in the same shape: a number read
without knowing what it measures.

1. `vm_stat` "Pages free" alone. The coordinator stopped a useful rebuild on 99 MB free. macOS keeps
   free low on purpose and that number is close to meaningless on its own.
2. free plus inactive. Better -- inactive pages are reclaimable, and it stood at 6.7 GB when free
   said 99 MB -- but it says nothing about whether the machine is paging.
3. the compressor. "Pages occupied by compressor" is in-kernel compression, not a process, and
   5.6 GB of it was once read as a job that could be stopped. Nothing can be stopped there.

What actually hurt was SWAP TRAFFIC: 833,625 swapouts against 389,240 swapins since boot, about 13 GB
out, with swap 5,034 MB used of 6,144. The machine was not short of memory at any instant; it had been
paging for hours, and that is what slows a suite. So the gate reads the RATE, not a level: two samples
a few seconds apart, and a job starts only if swapouts did not rise between them.

Disk is in the same gate because macOS grows swap files on disk, so a paging night eats both. The data
volume was at 97% with 17 GB free while swap had 1.1 GB left.

The thresholds here are this tool's own tolerances, set from one night's evidence. None of them is a
measured limit, and the refusal prints every figure so a reader can disagree with the number rather
than with the decision.

What this gate is, stated plainly because a guard credited with more than it does is its own defect:
it is a PRECONDITION, not a protection. Removing it does not let anything be corrupted or lost; it
only lets a heavy job start while the machine is paging, which makes every other lane slower. There is
no harm here that something else was already preventing, and equally no harm it prevents beyond that.

An unreadable signal is never a yes. If swapouts cannot be read, or free disk cannot be read, or
memory_pressure gives neither a level nor a free percentage, the gate refuses and says which reading
failed. This machine's memory_pressure prints a free percentage and no level line, so an absent level
alone is normal and is not treated as a failure.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import time
from dataclasses import dataclass

SAMPLE_SECONDS = 10
MIN_FREE_DISK_GB = 10.0
PRESSURE_REFUSED = ("warn", "critical")


@dataclass(frozen=True)
class Reading:
    swapouts: int | None
    pressure_free_percent: int | None
    pressure_level: str | None
    free_disk_gb: float | None


def vm_swapouts(text: str) -> int | None:
    m = re.search(r"^Swapouts:\s+(\d+)", text, re.M)
    return int(m.group(1)) if m else None


def read_vm_stat() -> str:
    return subprocess.run(["vm_stat"], capture_output=True, text=True, check=False).stdout


def read_pressure() -> tuple[int | None, str | None]:
    out = subprocess.run(["memory_pressure"], capture_output=True, text=True, check=False).stdout
    pct = re.search(r"free percentage:\s*(\d+)", out, re.I)
    level = re.search(r"pressure level:\s*(\w+)", out, re.I)
    return (int(pct.group(1)) if pct else None, level.group(1).lower() if level else None)


def read_free_disk_gb(path: str = ".") -> float | None:
    out = subprocess.run(["df", "-k", path], capture_output=True, text=True, check=False).stdout
    rows = out.strip().splitlines()
    if len(rows) < 2:
        return None
    parts = rows[-1].split()
    try:
        return int(parts[3]) / (1024 * 1024)
    except (IndexError, ValueError):
        return None


def decide(before: int | None, after: int | None, reading: Reading) -> tuple[bool, list[str]]:
    """Return (may_start, reasons). An unreadable signal never silently permits."""
    reasons: list[str] = []
    if before is None or after is None:
        reasons.append("swapouts could not be read from vm_stat, so paging cannot be ruled out")
    elif after > before:
        reasons.append(
            f"the machine is paging: swapouts rose {before} -> {after} ({after - before}) while sampling"
        )
    if reading.pressure_level in PRESSURE_REFUSED:
        reasons.append(f"memory_pressure reports {reading.pressure_level}")
    elif reading.pressure_level is None and reading.pressure_free_percent is None:
        # this machine's memory_pressure prints a free percentage and no level line, so an absent
        # level is normal here; both absent means the tool told us nothing and must not pass as a yes
        reasons.append("memory_pressure reported neither a level nor a free percentage")
    if reading.free_disk_gb is None:
        reasons.append("free disk could not be read, so the swap file's room cannot be ruled out")
    elif reading.free_disk_gb < MIN_FREE_DISK_GB:
        reasons.append(
            f"free disk {reading.free_disk_gb:.1f} GB is below this tool's floor of "
            f"{MIN_FREE_DISK_GB:.0f} GB, and macOS grows its swap files on disk"
        )
    return (not reasons, reasons)


def sample(seconds: int = SAMPLE_SECONDS) -> tuple[int | None, int | None, Reading]:
    before = vm_swapouts(read_vm_stat())
    time.sleep(seconds)
    after = vm_swapouts(read_vm_stat())
    pct, level = read_pressure()
    return before, after, Reading(after, pct, level, read_free_disk_gb())


def figures(before: int | None, after: int | None, r: Reading) -> str:
    rose = "unknown" if before is None or after is None else str(after - before)
    disk = "unknown" if r.free_disk_gb is None else f"{r.free_disk_gb:.1f} GB"
    return (
        f"swapouts {before} -> {after} (rose {rose}), "
        f"memory_pressure level {r.pressure_level or 'unknown'} "
        f"free {r.pressure_free_percent or 'unknown'}%, free disk {disk}"
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--seconds", type=int, default=SAMPLE_SECONDS, help="gap between the two samples")
    ap.add_argument("--job", default="a heavy job", help="what is being gated, for the message")
    ap.add_argument(
        "--override",
        metavar="REASON",
        help="start anyway, naming why; the reason is printed, never blank, beside the figures",
    )
    a = ap.parse_args(argv)
    before, after, r = sample(a.seconds)
    ok, reasons = decide(before, after, r)
    line = figures(before, after, r)
    if ok:
        print(f"capacity_gate: {a.job} may start. {line}")
        return 0
    print(f"capacity_gate: REFUSED {a.job}.", file=sys.stderr)
    for why in reasons:
        print(f"  - {why}", file=sys.stderr)
    print(f"  {line}", file=sys.stderr)
    if a.override:
        if not a.override.strip():
            print("capacity_gate: --override needs a reason", file=sys.stderr)
            return 2
        print(f"capacity_gate: starting anyway, because: {a.override}", file=sys.stderr)
        return 0
    print("  Wait, or pass --override with a reason that will be read later.", file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
