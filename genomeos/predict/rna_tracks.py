# SPDX-License-Identifier: AGPL-3.0-or-later
"""Predicted RNA over the genome (AlphaGenome RNA-seq tracks), as evidence of transcription for the parser.

Chromatin at a start says a promoter could be used; RNA over the exons says the gene is made. The
model's RNA-seq tracks give predicted coverage per base and strand for named tissues; this module
fetches a fixed panel of tissues per 1 Mb window, keeps the runs of bases whose maximal predicted
coverage clears a threshold (sparse, cached under data/knowledge/alphagenome/rna, local), and answers
"what fraction of this interval on this strand is predicted to be transcribed?" so a candidate gene
can be held to it. Predicted evidence; the panel is a handful of tissues, so a gene expressed only
elsewhere reads as silent here, and that is stated with every result.
"""

from __future__ import annotations

import bisect
import json
import time
from pathlib import Path
from typing import Any

CACHE = Path("data/knowledge/alphagenome/rna")
WINDOW = 1_048_576
COVERAGE = 0.5  # predicted coverage a base must reach in some track of the panel
PANEL = {  # UBERON terms the model has RNA-seq tracks for; a spread, not the whole body
    "UBERON:0002107": "liver",
    "UBERON:0000955": "brain",
    "UBERON:0002048": "lung",
    "UBERON:0001157": "transverse colon",
    "UBERON:0002113": "kidney",
    "UBERON:0000948": "heart",
    "UBERON:0002367": "prostate gland",
    "UBERON:0000992": "ovary",
}


def window_starts(length: int, window: int = WINDOW) -> list[int]:
    if length <= window:
        return [0]
    starts = list(range(0, length - window, window))
    starts.append(length - window)
    return starts


def fetch_window(client, chrom: str, start: int, coverage: float = COVERAGE) -> dict[str, Any]:
    """One request: per strand, runs [a, b) (genomic, 0-based) where max predicted coverage ≥ threshold."""
    import numpy as np
    from alphagenome.data import genome as ag  # type: ignore[import-not-found]
    from alphagenome.models import dna_client  # type: ignore[import-not-found]

    iv = ag.Interval(chromosome=chrom, start=start, end=start + WINDOW)
    out = client.predict_interval(
        interval=iv, requested_outputs=[dna_client.OutputType.RNA_SEQ], ontology_terms=list(PANEL)
    )
    md = out.rna_seq.metadata
    values = out.rna_seq.values
    runs: dict[str, list[list[int]]] = {}
    for strand in ("+", "-"):
        cols = [i for i, s in enumerate(md["strand"]) if s == strand]
        if not cols:
            runs[strand] = []
            continue
        mx = values[:, cols].max(axis=1) >= coverage
        edges = np.flatnonzero(np.diff(np.concatenate(([0], mx.astype(np.int8), [0]))))
        runs[strand] = [
            [int(start + a), int(start + b)] for a, b in zip(edges[::2], edges[1::2], strict=False)
        ]
    return {"tracks": int(values.shape[1]), "runs": runs}


class RnaCoverage:
    def __init__(self, chrom: str, length: int, cache: Path = CACHE) -> None:
        self.chrom = chrom
        self.length = length
        self.cache = cache
        self.runs: dict[str, list[tuple[int, int]]] = {"+": [], "-": []}
        self._starts: dict[str, list[int]] = {"+": [], "-": []}
        self.requests = 0
        self.seconds = 0.0

    def load(self, client_factory, progress=None) -> RnaCoverage:
        client = None
        for start in window_starts(self.length):
            p = self.cache / self.chrom / f"{start}.json"
            if p.exists():
                d = json.loads(p.read_text())
            else:
                if client is None:
                    client = client_factory()
                t0 = time.time()
                failures = 0
                while True:
                    try:
                        d = fetch_window(client, self.chrom, start)
                        break
                    except Exception as ex:  # noqa: BLE001 - wait and retry, the API drops calls
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
                p.write_text(json.dumps(d))
                if progress:
                    progress(f"{self.chrom}:{start:,}: RNA predicted ({self.requests} requests)")
            for strand, rs in d["runs"].items():
                self.runs[strand].extend((a, b) for a, b in rs)
        for strand in self.runs:
            merged: list[tuple[int, int]] = []
            for a, b in sorted(self.runs[strand]):
                if merged and a <= merged[-1][1]:
                    merged[-1] = (merged[-1][0], max(merged[-1][1], b))
                else:
                    merged.append((a, b))
            self.runs[strand] = merged
            self._starts[strand] = [a for a, _ in merged]
        return self

    def covered_fraction(self, strand: str, start: int, end: int) -> float:
        """Fraction of [start, end) inside a predicted-transcribed run on that strand."""
        if end <= start:
            return 0.0
        runs, starts = self.runs[strand], self._starts[strand]
        i = max(0, bisect.bisect_right(starts, start) - 1)
        covered = 0
        while i < len(runs) and runs[i][0] < end:
            a, b = runs[i]
            covered += max(0, min(b, end) - max(a, start))
            i += 1
        return covered / (end - start)

    def summary(self) -> dict[str, Any]:
        return {
            "chrom": self.chrom,
            "panel": PANEL,
            "coverage_threshold": COVERAGE,
            "windows": len(window_starts(self.length)),
            "requests_this_run": self.requests,
            "seconds_this_run": round(self.seconds, 1),
            "transcribed_bases": {s: sum(b - a for a, b in rs) for s, rs in self.runs.items()},
        }


def filter_predictions(preds: list, coverage: RnaCoverage, min_fraction: float = 0.5) -> list:
    """Keep the parser's candidates whose exons are predicted transcribed (mean covered fraction)."""
    kept = []
    for p in preds:
        strand = p.strand.value
        fr = [coverage.covered_fraction(strand, a, b) for a, b in p.exons]
        if fr and sum(fr) / len(fr) >= min_fraction:
            kept.append(p)
    return kept
