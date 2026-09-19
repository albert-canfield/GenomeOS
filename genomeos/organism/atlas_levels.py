# SPDX-License-Identifier: AGPL-3.0-or-later
"""The Ma et al. 2021 atlas as levels over time, not only presence.

`tf_atlas` keeps one presence call per factor per cell (maximum adjusted expression at or above 20% of the
factor's maximum). Two questions of area E need more than that: whether a cell reads a factor as an
instantaneous threshold or as exposure integrated over time, and how a factor's level differs between
sister cells (POP-1 asymmetry, Notch targets). This module streams the same archive (88 MB, Zenodo 4737593,
CC BY 4.0) and keeps per factor per cell:

    peak       maximum adjusted expression over the cell's frames, as a fraction of the factor's maximum
    mean       mean over the frames the cell was tracked, same scale
    exposure   sum over frames x 1.25 min, same scale: minutes at the factor's maximum level
    first_on   first frame at or above the presence fraction (-1 if never)
    frames     frames the cell was tracked in the strain used

When a factor has several reporter strains, the strain with the higher peak in that cell is used, and each
strain's peak is also kept for the cells both tracked (`strain_peaks`: a measurement noise floor). Values
below 2% of the factor's maximum are dropped (sparse). The table stays local under data/knowledge; the
analyses that read it commit their summaries.
"""

from __future__ import annotations

import csv
import io
import json
import zipfile
from collections import defaultdict
from pathlib import Path

from .tf_atlas import FRAME_MIN, PRESENCE_FRACTION

LEVELS_FILE = Path("data/knowledge/celegans/atlas_levels.json")
ARCHIVE_FILE = Path("data/knowledge/celegans/ma2021_atlas.zip")
FLOOR = 0.02


def distil_levels(archive: bytes, fraction: float = PRESENCE_FRACTION, floor: float = FLOOR) -> dict:
    z = zipfile.ZipFile(io.BytesIO(archive))
    # factor -> strain -> cell -> [peak, total, frames, first frame at the threshold is decided later]
    series: dict[str, dict[str, dict[str, list[tuple[int, float]]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )
    for name in z.namelist():
        if not name.endswith(".csv"):
            continue
        stem = Path(name).stem
        tf = stem.split("_")[0]
        for row in csv.DictReader(io.StringIO(z.read(name).decode("utf-8", "replace"))):
            try:
                v = float(row["Adjustment-expression"])
                t = int(row["Time"])
            except ValueError:
                continue
            series[tf][stem][row["Cell-name"]].append((t, v))
    out: dict[str, dict] = {}
    per_strain: dict[str, dict[str, list[float]]] = {}
    for tf, strains in series.items():
        top = max(
            (max(v for _, v in pts) for cells in strains.values() for pts in cells.values()), default=0.0
        )
        if top <= 0:
            out[tf] = {"max": 0.0, "cells": {}}
            continue
        if len(strains) > 1:  # the same factor seen by two reporters: a per-cell measurement noise floor
            names = sorted(strains)
            both = set.intersection(*(set(strains[n]) for n in names))
            tops = {n: max(max(v for _, v in pts) for pts in strains[n].values()) for n in names}
            per_strain[tf] = {
                c: [
                    round(max(v for _, v in strains[n][c]) / tops[n], 3) if tops[n] > 0 else 0.0
                    for n in names
                ]
                for c in sorted(both)
            }
        best: dict[str, list[tuple[int, float]]] = {}
        for cells in strains.values():
            for cell, pts in cells.items():
                if cell not in best or max(v for _, v in pts) > max(v for _, v in best[cell]):
                    best[cell] = pts
        cells_out = {}
        for cell, pts in best.items():
            pts.sort()
            peak = max(v for _, v in pts) / top
            if peak < floor:
                continue
            total = sum(v for _, v in pts) / top
            on = next((t for t, v in pts if v >= fraction * top), -1)
            cells_out[cell] = [
                round(peak, 3),
                round(total / len(pts), 3),
                round(total * FRAME_MIN, 2),
                on,
                len(pts),
            ]
        out[tf] = {"max": round(top, 2), "cells": cells_out}
    return {
        "source": "Ma et al. 2021, Nat Methods 18:893 (Zenodo 4737593)",
        "fields": ["peak", "mean", "exposure_min", "first_on_frame", "frames"],
        "presence_fraction_of_max": fraction,
        "frame_min": FRAME_MIN,
        "factors": out,
        "strain_peaks": per_strain,
    }


def build(archive_path: Path = ARCHIVE_FILE, path: Path = LEVELS_FILE) -> Path:
    table = distil_levels(Path(archive_path).read_bytes())
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(table, separators=(",", ":")))
    return path


def load_strain_peaks(path: Path = LEVELS_FILE) -> dict[str, dict[str, list[float]]]:
    """factor -> cell -> peak per reporter strain (fraction of that strain's own max)."""
    return json.loads(Path(path).read_text()).get("strain_peaks", {})


def load_levels(path: Path = LEVELS_FILE) -> dict[str, dict[str, list]]:
    """factor -> cell -> [peak, mean, exposure_min, first_on_frame, frames]."""
    data = json.loads(Path(path).read_text())
    return {tf: rec["cells"] for tf, rec in data["factors"].items()}
