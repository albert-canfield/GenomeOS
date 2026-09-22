# SPDX-License-Identifier: AGPL-3.0-or-later
"""Lateral inhibition as noise plus selection: the reference behaviour a contact-signalling runtime must meet.

Collier et al. 1996 (J Theor Biol 183:429): each cell's Notch activity rises with its neighbours' Delta, and
its own Delta falls with its Notch activity. Two equivalent cells started exactly equal stay equal forever:
the symmetric state is a steady state, and a deterministic rule reading the same inputs gives the same answer.
The symmetric state is unstable, so any difference is amplified until one cell is high-Delta (the one that
keeps the primary fate) and the other is high-Notch. Which cell wins is decided by the noise, not by the
genome or the position: in the worm's AC/VU decision Z1.ppp and Z4.aaa each become the anchor cell about
half the time (Kimble & Hirsh 1979, Dev Biol 70:396; Seydoux & Greenwald 1989, Cell 57:1237), while the
reference lineage records one outcome.

    dn_i/dt = d_j^k / (a + d_j^k) - n_i        (Notch activity from the neighbours' Delta)
    dd_i/dt = 1 / (1 + b n_i^h) - d_i          (Delta repressed by the cell's own Notch)

with Collier's a = 0.01, b = 100, k = h = 2, integrated by Euler steps with optional Gaussian noise on d.
"""

from __future__ import annotations

import random
from dataclasses import dataclass


@dataclass(slots=True)
class LateralResult:
    delta: list[float]
    notch: list[float]
    diverged: bool
    winner: int | None  # index of the high-Delta cell, None when the cells stayed equivalent


def two_cells(
    noise: float = 0.0,
    seed: int | None = None,
    start: tuple[float, float] = (0.5, 0.5),
    steps: int = 4000,
    dt: float = 0.01,
    a: float = 0.01,
    b: float = 100.0,
    k: float = 2.0,
    h: float = 2.0,
    gap: float = 0.5,
) -> LateralResult:
    rng = random.Random(seed)
    d = list(start)
    n = [0.5, 0.5]
    for _ in range(steps):
        dn = [(d[1 - i] ** k) / (a + d[1 - i] ** k) - n[i] for i in range(2)]
        dd = [1.0 / (1.0 + b * n[i] ** h) - d[i] for i in range(2)]
        n = [n[i] + dt * dn[i] for i in range(2)]
        d = [
            max(0.0, d[i] + dt * dd[i] + (rng.gauss(0.0, noise) * dt**0.5 if noise else 0.0))
            for i in range(2)
        ]
    diverged = abs(d[0] - d[1]) >= gap
    return LateralResult(d, n, diverged, (0 if d[0] > d[1] else 1) if diverged else None)


def selection_statistics(noise: float = 0.02, runs: int = 200, seed: int = 0) -> dict:
    """Over `runs` pairs of equivalent cells: how often they diverge and how often the first cell wins."""
    wins = diverged = 0
    for r in range(runs):
        res = two_cells(noise=noise, seed=seed + r)
        diverged += res.diverged
        wins += res.winner == 0
    deterministic = two_cells(noise=0.0)
    return {
        "model": "Collier et al. 1996 two-cell lateral inhibition",
        "deterministic_equal_start_diverges": deterministic.diverged,
        "noise_sd": noise,
        "runs": runs,
        "diverged": diverged,
        "first_cell_wins": wins,
        "first_cell_win_share": round(wins / diverged, 3) if diverged else None,
    }
