"""Reader v1: which nodes a cell type reads.

The genome is the same in every cell; what differs is which parts are open.
ENCODE's DNase-seq peaks for a cell type say where the chromatin is open.
Laid over the nodes (domains between CTCF boundaries), the promoters and the
enhancers, they give a first reader: per node, how much is open; per gene,
whether its promoter is open (read) or closed (silent); per enhancer,
whether it is active in that cell type. Two cell types compared give the
genes one reads and the other does not.

Data: one narrowPeak file per cell type from the ENCODE portal (1-2 MB),
streamed, rows for the requested chromosomes kept under data/results.
Evidence: experimental (DNase-seq) for the peaks; inferred for "read": an
open promoter is necessary for transcription, not proof of it.
"""

from __future__ import annotations

import bisect
import gzip
import io
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

from genomeos.coords import Strand

ENCODE = "https://www.encodeproject.org"
RESULTS = Path("data/results")
PROMOTER_WINDOW = 1_000
EVIDENCE = "experimental: ENCODE DNase-seq narrowPeak (GRCh38, released)"


def slug(cell_type: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", cell_type).strip("_")


def find_dnase_file(cell_type: str, timeout: int = 60) -> dict[str, Any]:
    """The largest released GRCh38 DNase narrowPeak file for a biosample term (K562, HepG2, liver…)."""
    q = (
        f"{ENCODE}/search/?type=File&assay_title=DNase-seq&file_format=bed&file_type=bed+narrowPeak"
        f"&assembly=GRCh38&status=released&output_type=peaks&limit=50&format=json"
        f"&biosample_ontology.term_name={urllib.parse.quote(cell_type)}"
    )
    req = urllib.request.Request(q, headers={"Accept": "application/json", "User-Agent": "GenomeOS/0.1"})
    with urllib.request.urlopen(req, timeout=timeout) as r:  # noqa: S310
        d = json.load(r)
    files = d.get("@graph", [])
    if not files:
        raise LookupError(f"no released GRCh38 DNase peaks for biosample {cell_type!r} on the ENCODE portal")
    best = max(files, key=lambda f: f.get("file_size", 0))
    return {
        "accession": best["accession"],
        "href": ENCODE + best["href"],
        "size": best.get("file_size"),
        "cell_type": cell_type,
        "candidates": len(files),
    }


def peaks_path(cell_type: str, chrom: str) -> Path:
    return RESULTS / f"dnase_{slug(cell_type)}_{chrom}.bed.gz"


def fetch_peaks(cell_type: str, chroms: set[str], timeout: int = 300) -> dict[str, Any]:
    """Stream the peak file once and keep the rows of the chromosomes asked for."""
    info = find_dnase_file(cell_type)
    req = urllib.request.Request(info["href"], headers={"User-Agent": "GenomeOS/0.1 (stream)"})
    kept: dict[str, list[tuple[int, int, float]]] = {c: [] for c in chroms}
    with (
        urllib.request.urlopen(req, timeout=timeout) as resp,  # noqa: S310
        gzip.open(io.BufferedReader(resp, 1 << 20), "rt") as fh,
    ):
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7 or f[0] not in kept:
                continue
            kept[f[0]].append((int(f[1]), int(f[2]), float(f[6])))
    RESULTS.mkdir(parents=True, exist_ok=True)
    for c, rows in kept.items():
        rows.sort()
        with gzip.open(peaks_path(cell_type, c), "wt") as out:
            out.write(f"# {info['accession']} {cell_type} DNase-seq narrowPeak, {c}, signalValue\n")
            for s, e, v in rows:
                out.write(f"{s}\t{e}\t{v:.2f}\n")
    info["kept"] = {c: len(v) for c, v in kept.items()}
    return info


def load_peaks(cell_type: str, chrom: str) -> list[tuple[int, int, float]]:
    p = peaks_path(cell_type, chrom)
    if not p.exists():
        return []
    out = []
    with gzip.open(p, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            s, e, v = line.rstrip("\n").split("\t")
            out.append((int(s), int(e), float(v)))
    return out


class PeakIndex:
    def __init__(self, peaks: list[tuple[int, int, float]]) -> None:
        self.peaks = sorted(peaks)
        self.starts = [p[0] for p in self.peaks]
        self.max_len = max((e - s for s, e, _ in self.peaks), default=0)

    def overlapping(self, start: int, end: int) -> list[tuple[int, int, float]]:
        i = bisect.bisect_left(self.starts, start - self.max_len)
        out = []
        for s, e, v in self.peaks[i:]:
            if s >= end:
                break
            if e > start:
                out.append((s, e, v))
        return out

    def covered_bp(self, start: int, end: int) -> int:
        return sum(min(e, end) - max(s, start) for s, e, _ in self.overlapping(start, end))


def read_chromosome(cell_type: str, chrom: str, annotation, domains: list, ccres: list) -> dict[str, Any]:
    """What this cell type reads on one chromosome: open nodes, read genes, active enhancers."""
    idx = PeakIndex(load_peaks(cell_type, chrom))
    if not idx.peaks:
        raise FileNotFoundError(f"no peaks for {cell_type} on {chrom}; run fetch_peaks first")
    genes = [g for g in annotation.genes.values() if g.locus.chrom == chrom]
    read, silent = [], []
    for g in genes:
        if g.type != "protein_coding":
            continue
        tss = g.locus.end - 1 if g.locus.strand is Strand.MINUS else g.locus.start
        hits = idx.overlapping(tss - PROMOTER_WINDOW, tss + PROMOTER_WINDOW)
        (read if hits else silent).append(
            {"gene": g.symbol, "signal": round(max((v for _, _, v in hits), default=0.0), 2)}
        )
    read.sort(key=lambda r: -r["signal"])
    nodes = []
    for d in domains:
        cov = idx.covered_bp(d.start, d.end)
        n_peaks = len(idx.overlapping(d.start, d.end))
        nodes.append(
            {
                "id": d.id,
                "start": d.start,
                "end": d.end,
                "peaks": n_peaks,
                "open_fraction": round(cov / max(1, d.length), 4),
                "peaks_per_100kb": round(n_peaks / max(1, d.length) * 100_000, 2),
                "coding_genes": d.coding_genes,
            }
        )
    enh = [c for c in ccres if c.cls in ("pELS", "dELS")]
    active = sum(1 for c in enh if idx.overlapping(c.start, c.end))
    median_density = sorted(n["peaks_per_100kb"] for n in nodes)[len(nodes) // 2] if nodes else 0
    open_nodes = [n for n in nodes if n["peaks_per_100kb"] >= max(1.0, median_density)]
    silent_nodes = [n for n in nodes if n["peaks"] == 0]
    return {
        "cell_type": cell_type,
        "chrom": chrom,
        "peaks": len(idx.peaks),
        "coding_genes": len(read) + len(silent),
        "genes_read": len(read),
        "genes_silent": len(silent),
        "read_fraction": round(len(read) / max(1, len(read) + len(silent)), 4),
        "top_read": read[:25],
        "_read_all": [r["gene"] for r in read],
        "silent_genes": [s["gene"] for s in silent][:200],
        "enhancers": len(enh),
        "enhancers_active": active,
        "enhancers_active_fraction": round(active / max(1, len(enh)), 4),
        "nodes": len(nodes),
        "nodes_open": len(open_nodes),
        "nodes_silent": len(silent_nodes),
        "silent_node_ids": [n["id"] for n in silent_nodes][:100],
        "node_table": nodes,
        "evidence": {
            "peaks": EVIDENCE,
            "read": "inferred: promoter (TSS ± 1 kb) overlaps a DNase peak; "
            "open is necessary for transcription, not proof of it",
            "nodes": "inferred: CTCF domains",
        },
    }


def compare(a: dict[str, Any], b: dict[str, Any], annotation=None) -> dict[str, Any]:
    """Genes read in one cell type and silent in the other."""
    ra = {r["gene"] for r in a["top_read"]} | ({g for g in _all_read(a)})
    rb = {r["gene"] for r in b["top_read"]} | ({g for g in _all_read(b)})
    sa, sb = set(a["silent_genes"]), set(b["silent_genes"])
    return {
        "a": a["cell_type"],
        "b": b["cell_type"],
        "chrom": a["chrom"],
        "read_in_a_only": sorted(ra & sb)[:60],
        "read_in_b_only": sorted(rb & sa)[:60],
        "read_in_both": len(ra & rb),
        "evidence": "inferred from DNase peaks at promoters in each cell type",
    }


def _all_read(r: dict[str, Any]) -> set[str]:
    return set(r.get("_read_all", []))
