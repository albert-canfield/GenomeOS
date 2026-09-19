# SPDX-License-Identifier: AGPL-3.0-or-later
"""Splice sites predicted at 1 bp (AlphaGenome feature c), as a second opinion for the segment parser.

The learned donor and acceptor matrices reach about 90% recall at seven false
hits per kilobase, which is why the grammar parser asserts nothing. AlphaGenome
predicts splice sites at base resolution over a 1 Mb window in under a second:
four tracks, donor and acceptor on each strand, each a probability per base.
Calibrated against GENCODE on chr21 the tracks mark the exonic boundary base in
genomic coordinates: the donor track peaks on the last exonic base before GT,
the acceptor track on the first exonic base after AG (mean probability 0.91 to
1.00 on canonical exon boundaries in three 1 Mb windows, both strands).

Answers are cached per 1 Mb window (positions with probability ≥ 0.02 only) under
data/knowledge/alphagenome/splice, local and never committed. Everything that
comes out is `predicted` evidence; the grammar stays the fallback where the key
is missing.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any

CACHE = Path("data/knowledge/alphagenome/splice")
WINDOW = 1_048_576  # dna_client.SEQUENCE_LENGTH_1MB
KEEP_ABOVE = 0.02  # probabilities stored in the cache
MIN_PROB = 0.05  # probabilities offered to the parser as candidates
# bits added to log2(p/(1-p)) so that a confident site lands on the scale of the learned matrices
# (max about 12 bits); 10.0 beat 13.3 on chr21 exons (40.1%/46.2% against 39.0%/41.9%)
BACKGROUND_LOG_ODDS = 10.0
TRACKS = ("donor+", "acceptor+", "donor-", "acceptor-")


def score_of(prob: float) -> float:
    """A probability as bits over background, on the scale of the learned matrices (max about 12)."""
    p = min(max(prob, 1e-6), 1 - 1e-6)
    return math.log2(p / (1 - p)) + BACKGROUND_LOG_ODDS


def window_starts(length: int, window: int = WINDOW) -> list[int]:
    """Window starts covering 0..length; the last window is pulled back so it ends at the chromosome end."""
    if length <= window:
        return [0]
    starts = list(range(0, length - window, window))
    starts.append(length - window)
    return starts


def cache_path(chrom: str, start: int, cache: Path = CACHE) -> Path:
    return cache / chrom / f"{start}.json"


def fetch_window(client, chrom: str, start: int, window: int = WINDOW) -> dict[str, list[list[float]]]:
    """One live request: sparse (position, probability) per track, positions genomic 0-based."""
    from alphagenome.data import genome as ag  # type: ignore[import-not-found]
    from alphagenome.models import dna_client  # type: ignore[import-not-found]

    iv = ag.Interval(chromosome=chrom, start=start, end=start + window)
    out = client.predict_interval(
        interval=iv, requested_outputs=[dna_client.OutputType.SPLICE_SITES], ontology_terms=None
    )
    values = out.splice_sites.values  # (window, 4): donor+, acceptor+, donor-, acceptor-
    names = [f"{r['name']}{r['strand']}" for _, r in out.splice_sites.metadata.iterrows()]
    sparse: dict[str, list[list[float]]] = {}
    for k, name in enumerate(names):
        col = values[:, k]
        idx = (col >= KEEP_ABOVE).nonzero()[0]
        sparse[name] = [[int(start + i), round(float(col[i]), 4)] for i in idx]
    return sparse


class SpliceSites:
    """Predicted donors and acceptors for one chromosome, served to the parser in local coordinates."""

    def __init__(self, chrom: str, length: int, cache: Path = CACHE) -> None:
        self.chrom = chrom
        self.length = length
        self.cache = cache
        self.tracks: dict[str, dict[int, float]] = {t: {} for t in TRACKS}
        self.requests = 0
        self.seconds = 0.0

    def load(self, client_factory, progress=None, min_prob: float = MIN_PROB) -> SpliceSites:
        """Every window from the cache, fetching what is missing; overlaps keep the larger value."""
        client = None
        for start in window_starts(self.length):
            p = cache_path(self.chrom, start, self.cache)
            if p.exists():
                sparse = json.loads(p.read_text())
            else:
                if client is None:
                    client = client_factory()
                t0 = time.time()
                failures = 0
                while True:
                    try:
                        sparse = fetch_window(client, self.chrom, start)
                        break
                    except Exception as ex:  # noqa: BLE001 - the API drops calls; wait and retry
                        failures += 1
                        if failures > 8:
                            raise
                        if progress:
                            progress(f"{self.chrom}:{start:,}: {type(ex).__name__}, retry {failures}")
                        time.sleep(min(120, 5 * failures))
                        client = client_factory()
                self.requests += 1
                self.seconds += time.time() - t0
                p.parent.mkdir(parents=True, exist_ok=True)
                p.write_text(json.dumps(sparse))
                if progress:
                    progress(f"{self.chrom}:{start:,}: splice sites predicted ({self.requests} requests)")
            for name, rows in sparse.items():
                d = self.tracks.setdefault(name, {})
                for pos, prob in rows:
                    if prob >= min_prob and prob > d.get(pos, 0.0):
                        d[pos] = prob
        return self

    def counts(self, min_prob: float = 0.5) -> dict[str, int]:
        return {t: sum(1 for v in d.values() if v >= min_prob) for t, d in self.tracks.items()}

    def local(
        self, strand: str, offset: int, n: int
    ) -> tuple[list[tuple[int, float]], list[tuple[int, float]]]:
        """Donor and acceptor candidates for the parser's string: genomic [offset, offset+n) read on
        `strand` ('+' as is, '-' as its reverse complement). Positions follow the parser's conventions:
        a donor is the first intronic base (the G of GT), an acceptor the first exonic base after AG."""
        donors: list[tuple[int, float]] = []
        acceptors: list[tuple[int, float]] = []
        lo, hi = offset, offset + n
        if strand == "+":
            for pos, prob in self.tracks["donor+"].items():
                if lo <= pos < hi - 1:
                    donors.append((pos - offset + 1, score_of(prob)))  # mark = last exonic base
            for pos, prob in self.tracks["acceptor+"].items():
                if lo <= pos < hi:
                    acceptors.append((pos - offset, score_of(prob)))  # mark = first exonic base
        else:
            for pos, prob in self.tracks["donor-"].items():
                if lo < pos < hi:
                    x = n - 1 - (pos - offset)  # local index of the marked (last exonic) base
                    donors.append((x + 1, score_of(prob)))
            for pos, prob in self.tracks["acceptor-"].items():
                if lo <= pos < hi:
                    acceptors.append((n - 1 - (pos - offset), score_of(prob)))
        donors.sort()
        acceptors.sort()
        return donors, acceptors

    def summary(self) -> dict[str, Any]:
        return {
            "chrom": self.chrom,
            "windows": len(window_starts(self.length)),
            "requests_this_run": self.requests,
            "seconds_this_run": round(self.seconds, 1),
            "candidates_at_or_above_0.05": {t: len(d) for t, d in self.tracks.items()},
            "confident_at_or_above_0.5": self.counts(0.5),
        }


def client_factory():
    """A live AlphaGenome client from the key in the environment or .env (raises without one)."""
    from alphagenome.models import dna_client  # type: ignore[import-not-found]

    from genomeos.predict.alphagenome_adapter import AlphaGenomeAdapter

    key = AlphaGenomeAdapter().api_key
    if not key:
        raise RuntimeError("ALPHAGENOME_API_KEY is not set")
    return dna_client.create(key)
