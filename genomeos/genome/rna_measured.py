# SPDX-License-Identifier: AGPL-3.0-or-later
"""Measured RNA over the genome: ENCODE RNA-seq signal read from bigWigs, as evidence for the parser.

The predicted RNA filter (genomeos/predict/rna_tracks.py) asked a model; this asks an experiment. ENCODE
holds strand-specific total RNA-seq signal tracks (bigWig, GRCh38) for its cell lines; the bigWig
reader (genomeos/attribution/bigwig.py) reads only the sections that cover the intervals asked for,
over HTTP ranges, so a chromosome's worth of candidate exons costs a few tens of MB and no download.
A candidate gene stays when its exons carry signal on its strand.

Experimental evidence, one cell line at a time: a gene the line does not express is invisible here,
which is the same limit the reader has, stated with every result.
"""

from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

ENCODE = "https://www.encodeproject.org"
KNOWLEDGE = Path("data/knowledge/rna")
SIGNAL = 0.05  # signal a base must reach to count as transcribed; swept on chr21 (0.05, 0.2, 0.5, 1.0)
EVIDENCE = "experimental: ENCODE total RNA-seq, strand-specific signal of unique reads (GRCh38, released)"


def find_tracks(cell_type: str, timeout: int = 60) -> dict[str, Any]:
    """The plus- and minus-strand signal bigWigs of one released total RNA-seq experiment of a cell type."""
    KNOWLEDGE.mkdir(parents=True, exist_ok=True)
    p = KNOWLEDGE / f"tracks_{cell_type.replace(' ', '_')}.json"
    if p.exists():
        return json.loads(p.read_text())
    found: dict[str, dict[str, Any]] = {}
    for strand, out_type in (
        ("+", "plus strand signal of unique reads"),
        ("-", "minus strand signal of unique reads"),
    ):
        q = (
            f"{ENCODE}/search/?type=File&assay_title=total+RNA-seq&file_format=bigWig&assembly=GRCh38"
            f"&status=released&biosample_ontology.term_name={urllib.parse.quote(cell_type)}"
            f"&output_type={urllib.parse.quote(out_type)}&limit=20&format=json"
        )
        req = urllib.request.Request(q, headers={"Accept": "application/json", "User-Agent": "GenomeOS/0.1"})
        with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
            d = json.loads(r.read())
        for f in d.get("@graph", []):
            found.setdefault(f["dataset"], {})[strand] = {
                "accession": f["accession"],
                "href": ENCODE + f["href"],
                "size": f.get("file_size"),
            }
    pairs = [(ds, t) for ds, t in found.items() if "+" in t and "-" in t]
    if not pairs:
        raise LookupError(
            f"no released GRCh38 strand-specific total RNA-seq bigWigs for {cell_type!r} on ENCODE"
        )
    pairs.sort(key=lambda x: (x[1]["+"]["size"] or 0) + (x[1]["-"]["size"] or 0))  # the smallest pair
    ds, tracks = pairs[0]
    out = {"cell_type": cell_type, "experiment": ENCODE + ds, "tracks": tracks, "evidence": EVIDENCE}
    p.write_text(json.dumps(out))
    return out


class MeasuredRna:
    """Signal over the exons of the parser's candidates, one strand's bigWig at a time; with several cell
    types (a comma-separated panel) an exon's covered fraction is the best over the panel, so a gene made
    in any of them counts."""

    def __init__(
        self, cell_type: str, chrom: str, signal: float = SIGNAL, tracks: dict | None = None
    ) -> None:
        self.cell_type = cell_type
        self.chrom = chrom
        self.signal = signal
        self.panel = [c.strip() for c in cell_type.split(",") if c.strip()]
        self.missing: list[str] = []
        if tracks is not None:
            self.tracks_by_cell = {self.panel[0]: tracks}
        else:
            self.tracks_by_cell = {}
            for c in self.panel:
                try:
                    self.tracks_by_cell[c] = find_tracks(c)
                except (LookupError, OSError) as ex:
                    self.missing.append(f"{c}: {str(ex)[:60]}")
        self.tracks = next(iter(self.tracks_by_cell.values()), {"tracks": {}})
        self.fractions: dict[tuple[str, int, int], float] = {}
        self.bytes_fetched = 0

    def prepare(self, exons_by_strand: dict[str, list[tuple[int, int]]], progress=None) -> None:
        """One pass per cell type and strand over every exon interval (non-overlapping), keeping the best
        covered fraction seen."""
        from genomeos.attribution.bigwig import BigWig

        for tracks in self.tracks_by_cell.values():
            for strand, exons in exons_by_strand.items():
                if not exons or strand not in tracks["tracks"]:
                    continue
                bw = BigWig(tracks["tracks"][strand]["href"])
                try:
                    stats = bw.summarise(self.chrom, exons, self.signal, progress=progress)
                    for (a, b), st in zip(exons, stats, strict=False):
                        fr = min(1.0, st.above / (b - a)) if b > a else 0.0
                        key = (strand, a, b)
                        if fr > self.fractions.get(key, 0.0):
                            self.fractions[key] = fr
                    self.bytes_fetched += getattr(bw.src, "bytes_fetched", 0)
                finally:
                    bw.close()

    def covered_fraction(self, strand: str, start: int, end: int) -> float:
        return self.fractions.get((strand, start, end), 0.0)

    def filter(self, preds: list, min_fraction: float = 0.3, progress=None) -> list:
        """Keep the candidates whose exons carry signal on their strand (mean covered fraction)."""
        by_strand: dict[str, list[tuple[int, int]]] = {"+": [], "-": []}
        for p in preds:
            by_strand[p.strand.value].extend(p.exons)
        for s in by_strand:
            by_strand[s] = sorted(set(by_strand[s]))
        self.prepare(by_strand, progress)
        kept = []
        for p in preds:
            fr = [self.covered_fraction(p.strand.value, a, b) for a, b in p.exons]
            if fr and sum(fr) / len(fr) >= min_fraction:
                kept.append(p)
        return kept

    def summary(self) -> dict[str, Any]:
        return {
            "cell_type": self.cell_type,
            "panel": list(self.tracks_by_cell),
            "missing": self.missing,
            "experiments": {c: tr.get("experiment") for c, tr in self.tracks_by_cell.items()},
            "tracks": {
                c: {s: t["accession"] for s, t in tr["tracks"].items()}
                for c, tr in self.tracks_by_cell.items()
            },
            "signal_threshold": self.signal,
            "exons_measured": len(self.fractions),
            "bytes_fetched": self.bytes_fetched,
            "evidence": EVIDENCE,
        }
