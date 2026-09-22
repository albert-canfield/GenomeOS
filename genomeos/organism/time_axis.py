# SPDX-License-Identifier: AGPL-3.0-or-later
"""Is the reference lineage's minute a real minute? The WormWeb tree gives every division a time, but
its absolute scale was never re-verified. Two independent clocks exist for the same cells: the Ma
2021 atlas tracks each cell from its first to its last imaging frame (1.25 min per frame), and Packer
2019 gives every sequenced cell an estimated embryo time. Fitting both against the WormWeb times of
the same named cells gives the scale of the reference axis; the result is stored, and the reference
itself is left as the source gave it, so that topology, fates and deaths stay exact and timing
statements carry the conversion.
"""

from __future__ import annotations

import csv
import gzip
import io
import json
import statistics
import urllib.request
from collections import defaultdict
from pathlib import Path

from .packer import GEO_URL
from .reference import ReferenceLineage
from .tf_atlas import CELLS_FILE, FRAME_MIN


def _fit(xs: list[float], ys: list[float]) -> dict:
    n = len(xs)
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    slope = sum((x - mx) * (y - my) for x, y in zip(xs, ys, strict=True)) / sxx if sxx else 0.0
    icpt = my - slope * mx
    res = [y - (slope * x + icpt) for x, y in zip(xs, ys, strict=True)]
    return {
        "n": n,
        "slope": round(slope, 4),
        "intercept_min": round(icpt, 1),
        "residual_sd_min": round(statistics.pstdev(res), 1),
    }


def atlas_lifetimes(path: Path = CELLS_FILE) -> dict[str, tuple[float, float]]:
    data = json.loads(Path(path).read_text())
    return {c: (f * FRAME_MIN, last * FRAME_MIN) for c, (f, last) in data.get("lifetimes", {}).items()}


def packer_times(url: str = GEO_URL, timeout: float = 120.0) -> dict[str, list[float]]:
    """exact lineage id -> embryo times (minutes) of the cells assigned to it."""
    with urllib.request.urlopen(url, timeout=timeout) as r:  # noqa: S310 (fixed public URL)
        text = gzip.decompress(r.read()).decode()
    out: dict[str, list[float]] = defaultdict(list)
    for row in csv.DictReader(io.StringIO(text)):
        lin, t = row.get("lineage", ""), row.get("embryo.time", "")
        if lin in ("", "NA") or t in ("", "NA") or "x" in lin or "/" in lin:
            continue
        out[lin].append(float(t))
    return out


def compare(
    ref: ReferenceLineage, lifetimes: dict[str, tuple[float, float]], packer: dict[str, list[float]]
) -> dict:
    born_x, born_y, div_x, div_y = [], [], [], []
    end_of_imaging = max((last for _, last in lifetimes.values()), default=0.0)
    for cell, (first, last) in lifetimes.items():
        c = ref.cells.get(cell)
        if c is None or c.divides is None or c.born <= 0:
            continue  # the first frame of the root and the cells alive at the start do not mark births
        born_x.append(c.born)
        born_y.append(first)
        if last < end_of_imaging - 2 * FRAME_MIN:  # a last frame at the end of imaging is not a division
            div_x.append(c.divides)
            div_y.append(last)
    mid_x, mid_y, inside = [], [], 0
    for lin, ts in packer.items():
        c = ref.cells.get(lin)
        if c is None or c.divides is None or len(ts) < 5:
            continue
        med = statistics.median(ts)
        mid_x.append((c.born + c.divides) / 2)
        mid_y.append(med)
        inside += c.born <= med <= c.divides
    atlas_born = _fit(born_x, born_y) if len(born_x) > 10 else None
    atlas_div = _fit(div_x, div_y) if len(div_x) > 10 else None
    pack = _fit(mid_x, mid_y) if len(mid_x) > 10 else None
    slopes = [f["slope"] for f in (atlas_born, atlas_div, pack) if f]
    scale = round(statistics.median(slopes), 3) if slopes else None
    factor = round(1 / scale, 2) if scale else None
    return {
        "reference_axis": "WormWeb minutes from first cleavage (source scale unverified until now)",
        "atlas_births": atlas_born,
        "atlas_divisions": atlas_div,
        "packer_midlife": pack,
        "packer_median_inside_lifetime": {"inside": inside, "of": len(mid_x)},
        "minutes_per_reference_minute": scale,
        "reading": (
            f"one reference minute is about {scale} minutes on the two independent clocks (Ma 2021 imaging "
            f"at {FRAME_MIN} min per frame; Packer 2019 estimated embryo times), so the reference axis runs "
            f"slow by a factor of about {factor}; topology, fates and deaths are unaffected, and a timing "
            "statement in reference minutes converts by this factor"
        )
        if scale
        else "not enough shared cells to fit",
    }


def distil() -> dict:
    return compare(ReferenceLineage.load(), atlas_lifetimes(), packer_times())
