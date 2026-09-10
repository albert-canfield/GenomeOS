"""Methylation at clock CpGs from nanopore bedMethyl files, streamed (no disk).

modkit/wf-human-variation bedMethyl rows: chrom, start, end, mod code (m = 5mC,
h = 5hmC, a = 6mA), score, strand, ..., column 10 = valid coverage, column 11 =
percent modified. Only the 5mC rows at the clock probe positions are kept; both
strands of a CpG (C at pos and pos+1) are pooled by coverage, as are haplotype
files. The result is a beta vector the clocks can score directly.
"""

from __future__ import annotations

import gzip
import io
import json
import urllib.request
from collections.abc import Iterable
from dataclasses import dataclass, field
from pathlib import Path

PACKAGED_PROBES = Path(__file__).resolve().parent / "data" / "clock_probes_hg38.json"


def load_probes(path: str | Path = PACKAGED_PROBES) -> dict[str, tuple[str, int]]:
    d = json.loads(Path(path).read_text())
    return {pid: (chrom, int(beg)) for pid, (chrom, beg) in d["probes"].items()}


@dataclass(slots=True)
class MethylationSummary:
    probes_total: int
    probes_covered: int
    betas: dict[str, float] = field(default_factory=dict)
    coverage: dict[str, int] = field(default_factory=dict)
    rows_scanned: int = 0

    @property
    def fraction_covered(self) -> float:
        return self.probes_covered / self.probes_total if self.probes_total else 0.0


def iter_bedmethyl(source: str, chunk: int = 1 << 20):
    """Yield (chrom, start, code, coverage, percent) for 5mC rows; source is a path or http(s) URL."""
    if source.startswith(("http://", "https://")):
        req = urllib.request.Request(source, headers={"User-Agent": "GenomeOS/0.1 (stream)"})
        resp = urllib.request.urlopen(req, timeout=300)  # noqa: S310, SIM115
        raw = gzip.GzipFile(fileobj=io.BufferedReader(resp, chunk)) if source.endswith(".gz") else resp
    else:
        raw = gzip.open(source, "rb") if source.endswith(".gz") else open(source, "rb")  # noqa: SIM115
    with raw:
        for line in io.TextIOWrapper(raw, encoding="ascii", errors="replace"):
            f = line.split("\t", 11)
            if len(f) < 11 or f[3] != "m":
                continue
            yield f[0], int(f[1]), float(f[9]), float(f[10])


def summarise(
    sources: Iterable[str], probes: dict[str, tuple[str, int]], min_coverage: int = 5, progress=None
) -> MethylationSummary:
    """Scan bedMethyl sources (files or URLs) and pool methylation at the probe CpGs."""
    by_pos: dict[tuple[str, int], str] = {}
    for pid, (chrom, beg) in probes.items():
        by_pos[(chrom, beg)] = pid
        by_pos[(chrom, beg + 1)] = pid  # the C on the opposite strand
    meth: dict[str, float] = dict.fromkeys(probes, 0.0)  # modified-read counts
    cov: dict[str, float] = dict.fromkeys(probes, 0.0)
    n = 0
    for src in sources:
        for chrom, start, coverage, percent in iter_bedmethyl(src):
            n += 1
            if progress and n % 5_000_000 == 0:
                progress(n)
            pid = by_pos.get((chrom, start))
            if pid is None:
                continue
            meth[pid] += coverage * percent / 100.0
            cov[pid] += coverage
    betas = {pid: meth[pid] / cov[pid] for pid in probes if cov[pid] >= min_coverage}
    return MethylationSummary(len(probes), len(betas), betas, {p: int(c) for p, c in cov.items() if c > 0}, n)
