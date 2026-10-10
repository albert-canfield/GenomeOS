# SPDX-License-Identifier: AGPL-3.0-or-later
"""Reference features for curated BioModels, and the score of an SBML engine against them.

ROADMAP §5 item 13 (area C): the libRoadRunner adapter is scored first on Schoeberl2002, then on the
models the in-house stand-in already ran. The reference for each model is BioModels' own curation
figure (`curation_image` in the entry's files, with the simulator the curators named in
`curation_notes`): the figure that certified the SBML reproduces the paper. The numbers below were
read from those images by pixel colour, on 2026-10-10, BEFORE either engine was run on these models
for this score; the reading, its axis calibration and the pixel scale of each axis are recorded per
feature so that the reading error can be told from the engine error.

The tolerance is fixed here and pinned literally in tests/test_roadrunner.py:

    a level or a time   |sim - ref| <= max(10% of |ref|, 3 pixels of that figure axis)
    a period            |sim - ref| <= 3% of ref

Three pixels covers the line width and the calibration of a figure read by eye-free pixel picking;
10% covers what the pixel bound cannot (a curve the curators' solver sampled coarsely, a peak a few
pixels wide). A model PASSES on an engine when every one of its features passes; nothing is averaged.

McClean2007 has no time-course reference (its curation figure is a parameter scan), so it is run on
both engines and compared, never scored.
"""

from __future__ import annotations

import math
from typing import Any

TOLERANCE = {"relative": 0.10, "pixels": 3, "period_relative": 0.03}
#: two engines agree on a model when no species, at any recorded time, differs by more than this
#: fraction of its range over the run (the larger of the two engines' ranges)
AGREEMENT = 0.02

# Feature tuples, by kind:
#   ("peak", quantity, (t0, t1), level, time)    the largest level in the window, and when
#   ("at", quantity, t, level)                   the level at time t
#   ("max" | "min", quantity, (t0, t1), level)   the largest / smallest level in the window
#   ("period", quantity, (t0, t1), period)       mean spacing of the window's peaks above mid-range
#   ("first_peak", quantity, t0, time)           the first such peak after t0
# `axes` gives the figure's units per pixel: "t" for time, and one per quantity for its level.

SCHOEBERL_FIGURE = (
    "BioModels BIOMD0000000019 curation_image.png: 'time course of active species as in fig. 2 c,d,e,f'; "
    "integration by SBMLODESolver (20080507); 4 panels, 0-60 min; red 50, blue 0.5, green 0.125 ng/ml EGF"
)


def _schoeberl(erk, mek, raf, ras) -> list[tuple]:
    """Features of one EGF dose: each total's peak (level, minute), and ERK-PP at 20, 40 and 59.6 min."""
    return [
        ("peak", "ERK_PP", (0.0, 60.0), *erk[:2]),
        ("at", "ERK_PP", 20.0, erk[2]),
        ("at", "ERK_PP", 40.0, erk[3]),
        ("at", "ERK_PP", 59.6, erk[4]),
        ("peak", "MEK_PP", (0.0, 60.0), *mek),
        ("peak", "Raf_act", (0.0, 60.0), *raf),
        ("peak", "Ras_GTP", (0.0, 60.0), *ras),
    ]


REFERENCES: dict[str, dict[str, Any]] = {
    "BIOMD0000000019": {
        "name": "Schoeberl2002 - EGF MAPK",
        "figure": SCHOEBERL_FIGURE,
        "time_unit": "min (the SBML's own)",
        "duration": 60.0,
        "standin_dt": 0.01,  # the committed run's dt (data/results/kinetic_BIOMD0000000019.json)
        "record_dt": 0.01,
        "axes": {
            "t": 60 / 277,
            "ERK_PP": 1.2e7 / 188,
            "MEK_PP": 2e5 / 188,
            "Raf_act": 2e4 / 195,
            "Ras_GTP": 2e4 / 195,
        },
        # the figure's totals are the SBML's own assignment rules (ERK_PP = x59 + x83, and so on)
        "doses": {
            "50 ng/ml (x1 = 4962, the SBML default)": {
                "initial": {},
                "features": _schoeberl(
                    (1.034e7, 3.68, 6.83e6, 3.96e6, 1.72e6),
                    (1.319e5, 3.90),
                    (1.651e4, 1.30),
                    (1.631e4, 1.08),
                ),
            },
            "0.5 ng/ml (x1 = 49.62)": {
                "initial": {"x1": 49.62},
                "features": _schoeberl(
                    (9.26e6, 6.72, 5.94e6, 3.32e6, 1.34e6),
                    (3.72e4, 6.72),
                    (1.056e4, 3.25),
                    (7.90e3, 2.82),
                ),
            },
            "0.125 ng/ml (x1 = 12.405)": {
                "initial": {"x1": 12.405},
                "features": _schoeberl(
                    (6.26e6, 11.26, 4.53e6, 1.91e6, 5.1e5),
                    (1.28e4, 10.83),
                    (6.56e3, 4.77),
                    (4.41e3, 4.77),
                ),
            },
        },
    },
    "BIOMD0000000010": {
        "name": "Kholodenko2000 - MAPK oscillations",
        "figure": (
            "BioModels BIOMD0000000010 curation_image.png: Figure 2A reproduced with Copasi 4.10 (Build 55); "
            "MAPK (green) and MAPK-PP (blue), nM, 0-150 min; x calibrated on the tick labels, 3.653 px/min"
        ),
        "time_unit": "s (SBML default); the figure's minutes x 60",
        "duration": 9000.0,
        "standin_dt": 0.05,
        "record_dt": 1.0,
        "axes": {"t": 60 / 3.653, "MAPK": 300 / 348, "MAPK_PP": 300 / 348},
        "doses": {
            "as curated": {
                "initial": {},
                "features": [
                    ("period", "MAPK", (900.0, 9000.0), 21.63 * 60),
                    ("first_peak", "MAPK", 900.0, 26.96 * 60),
                    ("max", "MAPK", (900.0, 9000.0), 204.3),
                    ("min", "MAPK", (900.0, 9000.0), 0.9),
                    ("max", "MAPK_PP", (1800.0, 9000.0), 296.6),
                    ("min", "MAPK_PP", (900.0, 9000.0), 37.1),
                    ("max", "MAPK_PP", (0.0, 900.0), 298.3),  # a plateau from 3.4 to ~15 min: no peak time
                ],
            }
        },
    },
    "BIOMD0000000012": {
        "name": "Elowitz2000 - Repressilator",
        "figure": (
            "BioModels BIOMD0000000012 curation_image.png: protein numbers as in figure 1c, deterministic "
            "run with SBML ODESolver; the notes say it departs from the paper before 400 min (initial "
            "conditions), so only 400-1000 min is read; red LacI (PX), orange TetR (PY), blue cI (PZ)"
        ),
        "time_unit": "min (the SBML's own)",
        "duration": 1000.0,
        "standin_dt": 0.05,  # as tests/test_sbml.py runs it
        "record_dt": 0.5,
        "axes": {"t": 1000 / 543, "PX": 2500 / 400, "PY": 2500 / 400, "PZ": 2500 / 400},
        "doses": {
            "as curated": {
                "initial": {},
                "features": [
                    ("period", "PX", (400.0, 1000.0), 125.7),
                    ("period", "PY", (400.0, 1000.0), 125.2),
                    ("period", "PZ", (400.0, 1000.0), 125.2),
                    # from 380 min: each protein's previous peak is before 367 min, its next one the reference
                    ("first_peak", "PX", 380.0, 489.9),
                    ("first_peak", "PY", 380.0, 449.4),
                    ("first_peak", "PZ", 380.0, 407.0),
                    ("max", "PX", (400.0, 1000.0), 2369.0),
                    ("max", "PY", (400.0, 1000.0), 2369.0),
                    ("max", "PZ", (400.0, 1000.0), 2369.0),
                    ("min", "PX", (400.0, 1000.0), 50.0),
                    ("min", "PY", (400.0, 1000.0), 50.0),
                    ("min", "PZ", (400.0, 1000.0), 50.0),
                ],
            }
        },
    },
    "BIOMD0000000084": {
        "name": "Hornberg2005 - ERK cascade",
        "figure": (
            "BioModels BIOMD0000000084 curation_image.jpeg: the plot of Fig 2, simulated with MathSBML; "
            "top R over 0-2, bottom x1p (green), x2p (blue), x3p (red) over 0-50; axes calibrated on the "
            "tick labels (top 139.5 px per unit time, 362.5 px per unit level; bottom 5.8 and 320)"
        ),
        "time_unit": "the SBML's own",
        "duration": 50.0,
        "standin_dt": 0.01,
        "record_dt": 0.01,
        "axes": {
            "t": 1 / 5.8,
            "t_R": 1 / 139.5,
            "R": 1 / 362.5,
            "x1p": 1 / 320,
            "x2p": 1 / 320,
            "x3p": 1 / 320,
        },
        "doses": {
            "as curated": {
                "initial": {},
                "features": [
                    ("at", "R", 0.25, 0.301),
                    ("at", "R", 0.5, 0.132),
                    ("at", "R", 0.75, 0.030),
                    ("peak", "x1p", (0.0, 50.0), 0.1359, 0.52),
                    ("peak", "x2p", (0.0, 50.0), 0.2016, 3.97),
                    ("peak", "x3p", (0.0, 50.0), 0.5609, 8.28),
                    ("at", "x2p", 10.0, 0.1078),
                    ("at", "x2p", 20.0, 0.0203),
                    ("at", "x3p", 20.0, 0.1797),
                    ("at", "x3p", 30.0, 0.0422),
                    ("at", "x3p", 40.0, 0.0266),
                    ("at", "x3p", 49.0, 0.0234),
                ],
            }
        },
    },
    "BIOMD0000000116": {
        "name": "McClean2007 - CrossTalk",
        "figure": (
            "BioModels BIOMD0000000116 curation_image.png is a parameter scan (Supplementary Figure 2, "
            "COPASI 4.0.20), not a time course: no feature is read, the engines are compared only"
        ),
        "time_unit": "the SBML's own",
        "duration": 100.0,
        "standin_dt": 0.01,
        "record_dt": 0.1,
        "axes": {},
        "doses": {"as curated": {"initial": {}, "features": []}},
    },
}


# ---- reading a trajectory -------------------------------------------------------------------------


def series(model, traj, quantity: str, initial: dict[str, float] | None = None) -> list[float]:
    """A quantity over the trajectory: a recorded species, or an assignment-rule variable evaluated on
    each recorded state with the model's own rules (the stand-in's compiled MathML, the same for both
    engines, so the totals are read identically whichever engine produced the states)."""
    if quantity in traj.levels:
        return list(traj.levels[quantity])
    rules = dict(model.assignments)
    if quantity not in rules:
        raise KeyError(f"{quantity!r} is neither a recorded species nor an assignment rule")
    fixed = {**model.parameters, **model.compartments, **model.species, **(initial or {})}
    out = []
    for i in range(len(traj.times)):
        env = dict(fixed)
        env.update({s: traj.levels[s][i] for s in traj.species})
        for var, fn in model.assignments:
            env[var] = fn(env)
        out.append(env[quantity])
    return out


def _window(times, xs, t0, t1):
    return [(t, x) for t, x in zip(times, xs, strict=True) if t0 - 1e-9 <= t <= t1 + 1e-9]


def _peaks(points) -> list[tuple[float, float]]:
    """Local maxima above the window's mid-range: an oscillation's peaks, not its numerical ripples."""
    if len(points) < 3:
        return []
    lo, hi = min(x for _, x in points), max(x for _, x in points)
    mid = (lo + hi) / 2
    return [
        points[i]
        for i in range(1, len(points) - 1)
        if points[i - 1][1] < points[i][1] >= points[i + 1][1] and points[i][1] > mid
    ]


def _bound(ref: float, units_per_pixel: float | None) -> float:
    pixel = TOLERANCE["pixels"] * units_per_pixel if units_per_pixel else 0.0
    return max(TOLERANCE["relative"] * abs(ref), pixel)


def score_feature(feature: tuple, times: list[float], xs: list[float], axes: dict[str, float]) -> list[dict]:
    """One feature against one simulated series: a check per number the feature names."""
    kind, quantity = feature[0], feature[1]
    t_axis = axes.get("t_R") if quantity == "R" else axes.get("t")
    v_axis = axes.get(quantity)
    checks: list[tuple[str, float | None, float, float]] = []  # (what, sim, ref, bound)
    if kind == "peak":
        (t0, t1), level, when = feature[2], feature[3], feature[4]
        pts = _window(times, xs, t0, t1)
        t_sim, x_sim = max(pts, key=lambda p: p[1]) if pts else (None, None)
        checks += [
            ("level", x_sim, level, _bound(level, v_axis)),
            ("time", t_sim, when, _bound(when, t_axis)),
        ]
    elif kind == "at":
        t, level = feature[2], feature[3]
        i = min(range(len(times)), key=lambda j: abs(times[j] - t))
        checks.append(("level", xs[i], level, _bound(level, v_axis)))
    elif kind in ("max", "min"):
        (t0, t1), level = feature[2], feature[3]
        pts = [x for _, x in _window(times, xs, t0, t1)]
        checks.append(
            ("level", (max if kind == "max" else min)(pts) if pts else None, level, _bound(level, v_axis))
        )
    elif kind == "period":
        (t0, t1), period = feature[2], feature[3]
        pk = _peaks(_window(times, xs, t0, t1))
        sim = (pk[-1][0] - pk[0][0]) / (len(pk) - 1) if len(pk) >= 2 else None
        checks.append(("period", sim, period, TOLERANCE["period_relative"] * period))
    elif kind == "first_peak":
        t0, when = feature[2], feature[3]
        pk = _peaks(_window(times, xs, t0, times[-1]))
        checks.append(("time", pk[0][0] if pk else None, when, _bound(when, t_axis)))
    else:
        raise ValueError(f"unknown feature kind {kind!r}")
    out = []
    for what, sim, ref, bound in checks:
        ok = sim is not None and math.isfinite(sim) and abs(sim - ref) <= bound
        out.append(
            {
                "feature": kind,
                "quantity": quantity,
                "where": feature[2],
                "check": what,
                "reference": ref,
                "simulated": None if sim is None else round(sim, 6),
                "bound": round(bound, 6),
                "pass": bool(ok),
            }
        )
    return out


def agreement(a, b) -> dict[str, Any]:
    """The largest difference between two trajectories on the same grid, per species as a fraction of
    its range (the larger of the two), and the species where it is largest."""
    if a.species != b.species or len(a.times) != len(b.times):
        raise ValueError("the trajectories are not on the same species and grid")
    worst, where = 0.0, None
    for s in a.species:
        xa, xb = a.levels[s], b.levels[s]
        span = max(max(xa) - min(xa), max(xb) - min(xb))
        if span <= 0:
            continue
        d = max(abs(p - q) for p, q in zip(xa, xb, strict=True)) / span
        if not math.isfinite(d):
            d = math.inf
        if d > worst:
            worst, where = d, s
    return {"max_relative_difference": worst, "species": where, "agree": worst <= AGREEMENT}
