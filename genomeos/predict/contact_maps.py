# SPDX-License-Identifier: AGPL-3.0-or-later
"""Node boundaries from a predicted contact map (AlphaGenome), held against the CTCF-only inference.

`genomeos domains` draws a node between CTCF-only elements and calls it inferred at 0.4, because the
project holds no Hi-C. AlphaGenome predicts a contact map at 2 kb over a 1 Mb window for 28 cell types
(4DN Micro-C and Hi-C tracks). The insulation score of each bin, the mean contact across the diagonal
in a 20 kb window, dips where a domain boundary is; its local minima are predicted boundaries. Where a
CTCF-only boundary sits on such a dip the node has a second, independent kind of evidence for the same
edge; where it does not, the inference is on its own. Cached per window under
data/knowledge/alphagenome/contact (local); predicted evidence, capped at 0.7.
"""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

CACHE = Path("data/knowledge/alphagenome/contact")
WINDOW = 1_048_576
BIN = 2_048
INSULATION_BINS = 10  # 20 kb on each side of the bin
MIN_DEPTH = 0.05  # a local minimum must sit this far below the profile's median; swept on chr21 (0.02 to 0.3)
TOLERANCE = 20_000  # a CTCF-only boundary within this of an insulation minimum agrees with it


def window_starts(length: int, window: int = WINDOW) -> list[int]:
    if length <= window:
        return [0]
    starts = list(range(0, length - window, window))
    starts.append(length - window)
    return starts


def insulation(values, w: int = INSULATION_BINS) -> list[float]:
    """Per bin: mean contact between the w bins upstream and the w bins downstream, averaged over the
    cell types; NaN where the window would leave the map."""
    import numpy as np

    v = np.nan_to_num(np.asarray(values, dtype=np.float32))
    n = v.shape[0]
    prof = np.full(n, np.nan, dtype=np.float32)
    mean_map = v.mean(axis=2) if v.ndim == 3 else v
    for i in range(w, n - w):
        prof[i] = float(mean_map[i - w : i, i : i + w].mean())
    return [float(x) for x in prof]


def fetch_window(client, chrom: str, start: int) -> dict[str, Any]:
    from alphagenome.data import genome as ag  # type: ignore[import-not-found]
    from alphagenome.models import dna_client  # type: ignore[import-not-found]

    iv = ag.Interval(chromosome=chrom, start=start, end=start + WINDOW)
    out = client.predict_interval(
        interval=iv, requested_outputs=[dna_client.OutputType.CONTACT_MAPS], ontology_terms=None
    )
    cm = out.contact_maps
    return {
        "start": start,
        "bin": int(getattr(cm, "resolution", BIN)),
        "tracks": int(cm.values.shape[2]) if cm.values.ndim == 3 else 1,
        "insulation": [round(x, 4) if x == x else None for x in insulation(cm.values)],
    }


def local_minima(profile: list[float | None], min_depth: float = MIN_DEPTH, spacing: int = 5) -> list[int]:
    """Bin indices of local minima that sit at least `min_depth` below the profile's median."""
    import statistics

    vals = [x for x in profile if x is not None]
    if not vals:
        return []
    med = statistics.median(vals)
    out = []
    n = len(profile)
    for i in range(spacing, n - spacing):
        x = profile[i]
        if x is None or med - x < min_depth:
            continue
        neigh = [profile[j] for j in range(i - spacing, i + spacing + 1) if j != i and profile[j] is not None]
        if neigh and all(x <= y for y in neigh):
            out.append(i)
    return out


class Insulation:
    def __init__(self, chrom: str, length: int, cache: Path = CACHE) -> None:
        self.chrom = chrom
        self.length = length
        self.cache = cache
        self.windows: list[dict[str, Any]] = []
        self.requests = 0

    def load(self, client_factory, progress=None) -> Insulation:
        client = None
        for start in window_starts(self.length):
            p = self.cache / self.chrom / f"{start}.json"
            if p.exists():
                d = json.loads(p.read_text())
            else:
                if client is None:
                    client = client_factory()
                failures = 0
                while True:
                    try:
                        d = fetch_window(client, self.chrom, start)
                        break
                    except Exception as ex:  # noqa: BLE001 - wait and retry
                        failures += 1
                        if failures > 8:
                            raise
                        if progress:
                            progress(f"{self.chrom}:{start:,}: {type(ex).__name__}, retry {failures}")
                        time.sleep(min(120, 5 * failures))
                        client = client_factory()
                self.requests += 1
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(json.dumps(d))
                if progress:
                    progress(f"{self.chrom}:{start:,}: contact map predicted ({self.requests} requests)")
            self.windows.append(d)
        return self

    def boundaries(self, min_depth: float = MIN_DEPTH) -> list[dict[str, Any]]:
        """Predicted boundaries in genomic coordinates (the bin's centre), deduplicated across windows."""
        out: dict[int, float] = {}
        for w in self.windows:
            prof = w["insulation"]
            for i in local_minima(prof, min_depth):
                pos = w["start"] + i * w["bin"] + w["bin"] // 2
                key = pos // 10_000
                depth = prof[i]
                if key not in out or depth < out[key]:
                    out[key] = depth
        return [{"pos": k * 10_000 + 5_000, "insulation": round(v, 4)} for k, v in sorted(out.items())]


def compare(inferred: list[int], predicted: list[int], tolerance: int = TOLERANCE) -> dict[str, Any]:
    """How the two boundary sets sit on each other: each inferred boundary's nearest predicted one and
    the reverse, and the fractions within tolerance."""
    import bisect

    def nearest(x: int, xs: list[int]) -> int | None:
        if not xs:
            return None
        i = bisect.bisect_left(xs, x)
        cands = [xs[j] for j in (i - 1, i) if 0 <= j < len(xs)]
        return min(cands, key=lambda c: abs(c - x)) if cands else None

    inf, pred = sorted(inferred), sorted(predicted)
    d_inf = [abs(x - nearest(x, pred)) for x in inf if nearest(x, pred) is not None]
    d_pred = [abs(x - nearest(x, inf)) for x in pred if nearest(x, inf) is not None]
    return {
        "inferred": len(inf),
        "predicted": len(pred),
        "tolerance": tolerance,
        "inferred_on_a_predicted_boundary": sum(1 for d in d_inf if d <= tolerance),
        "fraction_inferred_supported": round(sum(1 for d in d_inf if d <= tolerance) / len(inf), 3)
        if inf
        else None,
        "predicted_on_an_inferred_boundary": sum(1 for d in d_pred if d <= tolerance),
        "fraction_predicted_matched": round(sum(1 for d in d_pred if d <= tolerance) / len(pred), 3)
        if pred
        else None,
        "median_distance_inferred_to_predicted": sorted(d_inf)[len(d_inf) // 2] if d_inf else None,
    }


def random_control(
    inferred: list[int], predicted: list[int], length: int, seed: int = 1, n: int = 5
) -> float:
    """The supported fraction that as many boundaries placed uniformly at random would get: what the
    agreement is worth above chance, averaged over `n` draws."""
    import random

    rng = random.Random(seed)
    fr = []
    for _ in range(n):
        ctrl = sorted(rng.randint(1, length - 1) for _ in inferred)
        fr.append(compare(ctrl, predicted)["fraction_inferred_supported"] or 0.0)
    return round(sum(fr) / len(fr), 3)
