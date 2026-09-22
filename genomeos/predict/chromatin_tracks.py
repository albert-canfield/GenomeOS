# SPDX-License-Identifier: AGPL-3.0-or-later
"""Predicted chromatin per cell line over elements (AlphaGenome use 5: the reader where it is not measured).

AlphaGenome predicts a DNase-seq track at 1 bp for a named biosample over a 1 Mb window in about
a second. For a list of elements this module groups them by window, asks once per window for the
requested biosamples, keeps the mean predicted signal over each element per biosample, and caches
those means under data/knowledge/alphagenome/dnase (local). The 1 Mb arrays are never stored.
Everything that comes out is `predicted` evidence.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

CACHE = Path("data/knowledge/alphagenome/dnase")
WINDOW = 1_048_576  # dna_client.SEQUENCE_LENGTH_1MB
BIOSAMPLES = {"K562": "EFO:0002067", "HepG2": "EFO:0001187"}


def window_for(start: int, end: int, length: int, window: int = WINDOW) -> int:
    """The window start that holds the element whole, pulled back at the chromosome end."""
    ws = max(0, min(start - (start % window), length - window))
    if end > ws + window:  # element straddles a window edge: centre a window on it
        ws = max(0, min(start - window // 2, length - window))
    return ws


def fetch_means(
    client, chrom: str, ws: int, intervals: list[tuple[int, int, str]], biosamples: dict[str, str]
) -> dict[str, dict[str, float]]:
    """One live request; the mean predicted DNase over each interval per biosample."""
    from alphagenome.data import genome as ag  # type: ignore[import-not-found]
    from alphagenome.models import dna_client  # type: ignore[import-not-found]

    iv = ag.Interval(chromosome=chrom, start=ws, end=ws + WINDOW)
    out = client.predict_interval(
        interval=iv,
        requested_outputs=[dna_client.OutputType.DNASE],
        ontology_terms=list(biosamples.values()),
    )
    values = out.dnase.values
    col = {r["ontology_curie"]: k for k, r in out.dnase.metadata.iterrows()}
    means: dict[str, dict[str, float]] = {}
    for s, e, key in intervals:
        a, b = max(0, s - ws), min(WINDOW, e - ws)
        if b <= a:
            continue
        means[key] = {
            cell: round(float(values[a:b, col[curie]].mean()), 5)
            for cell, curie in biosamples.items()
            if curie in col
        }
    return means


class PredictedDnase:
    """Mean predicted DNase per element per biosample, one request per 1 Mb window, cached."""

    def __init__(
        self, chrom: str, length: int, biosamples: dict[str, str] | None = None, cache: Path = CACHE
    ) -> None:
        self.chrom = chrom
        self.length = length
        self.biosamples = biosamples or BIOSAMPLES
        self.cache = cache
        self.requests = 0
        self.seconds = 0.0

    def _path(self, ws: int) -> Path:
        tag = "_".join(sorted(self.biosamples))
        return self.cache / self.chrom / f"{ws}_{tag}.json"

    def means(
        self, intervals: list[tuple[int, int, str]], client_factory, progress=None
    ) -> dict[str, dict[str, float]]:
        by_window: dict[int, list[tuple[int, int, str]]] = {}
        for s, e, key in intervals:
            by_window.setdefault(window_for(s, e, self.length), []).append((s, e, key))
        out: dict[str, dict[str, float]] = {}
        client = None
        for n, (ws, items) in enumerate(sorted(by_window.items()), 1):
            p = self._path(ws)
            cached = json.loads(p.read_text()) if p.exists() else {}
            missing = [it for it in items if it[2] not in cached]
            if missing:
                if client is None:
                    client = client_factory()
                t0 = time.time()
                failures = 0
                while True:
                    try:
                        fresh = fetch_means(client, self.chrom, ws, missing, self.biosamples)
                        break
                    except Exception as ex:  # noqa: BLE001 - the API drops calls; wait and retry
                        failures += 1
                        if failures > 8:
                            raise
                        if progress:
                            progress(f"{self.chrom}:{ws:,}: {type(ex).__name__}, retry {failures}")
                        time.sleep(min(120, 5 * failures))
                        client = client_factory()
                self.requests += 1
                self.seconds += time.time() - t0
                cached.update(fresh)
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(json.dumps(cached))
                if progress:
                    progress(f"{self.chrom}: window {n}/{len(by_window)} predicted ({len(items)} elements)")
            for _, _, key in items:
                if key in cached:
                    out[key] = cached[key]
        return out
